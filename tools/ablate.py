#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
 multiphaseEulerFoam 算例消融诊断器（ablation probe runner）
================================================================================
 在本机 blueCFD(OpenFOAM 8) 环境里跑「单因素对照」探针，用于把
 "求解发散/崩溃" 归因到具体的一个配置因素，而不是靠猜。

 核心思路：每个探针 = base 算例的独立副本 + 只改一个因素 + 跑同一个求解器，
 汇总成一张表，看哪个因素改变了结果。

 用法：
   # 列出所有内置探针
   python tools/ablate.py --list

   # 跑内置探针
   python tools/ablate.py conv          # 收敛基线（含残差轨迹）
   python tools/ablate.py --custom NAME  # 跑自定义探针（见 CUSTOM 段）

 依赖：
   * 工程内置Python（含 numpy）：C:/Users/LELE/.workbuddy/binaries/python/envs/default/
   * blueCFD-Core 2020-1：默认 D:/APP/openfoam（用 OF8_HOME 覆盖）
================================================================================
"""
from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from of8 import run as of_run  # noqa: E402

BASE = os.path.join(ROOT, 'case')
PROBE_ROOT = os.path.join(ROOT, '_probes')
SOLVER = 'multiphaseEulerFoam'


# ======================================================================
#  探针基础设施
# ======================================================================
def fresh(name, endtime='0.5', writeint='0.25'):
    """从 base 复制一份探针算例；用时间戳避免 rmtree（沙箱禁止删除）。"""
    stamp = time.strftime('%H%M%S')
    d = os.path.join(PROBE_ROOT, f'{name}_{stamp}')
    shutil.copytree(BASE, d)
    cd = os.path.join(d, 'system', 'controlDict')
    t = open(cd, encoding='utf-8').read()
    t = re.sub(r'endTime\s+[\d.eE+-]+;', f'endTime         {endtime};', t)
    t = re.sub(r'writeInterval   [\d.eE+-]+;', f'writeInterval   {writeint};', t)
    open(cd, 'w', encoding='utf-8').write(t)
    return d


def _sub(path, pattern, repl, flags=0):
    t = open(path, encoding='utf-8').read()
    t2 = re.sub(pattern, repl, t, flags=flags)
    if t2 == t:
        raise RuntimeError(f'未匹配到内容，请检查模式: {pattern[:60]}  ({path})')
    open(path, 'w', encoding='utf-8').write(t2)


def set_internal(d, field, val):
    _sub(os.path.join(d, '0', field),
         r'internalField\s+uniform\s+[\d.eE+-]+;',
         f'internalField   uniform {val};')


def set_dt(d, dt, adjust=False):
    cd = os.path.join(d, 'system', 'controlDict')
    _sub(cd, r'deltaT\s+[\d.eE+-]+;', f'deltaT          {dt};')
    if not adjust:
        _sub(cd, r'adjustTimeStep\s+yes;', 'adjustTimeStep  no;')


def set_wall_T(d, val, field='T.gas'):
    _sub(os.path.join(d, '0', field),
         r'(coldWall\s*\{[^}]*?)value\s+uniform\s+[\d.]+;',
         rf'\g<1>value           uniform {val};', flags=re.S)


def set_wall_adiabatic(d, field='T.gas'):
    _sub(os.path.join(d, '0', field),
         r'(coldWall\s*\{[^}]*?)type\s+fixedValue;\s*value\s+uniform\s+[\d.]+;',
         r'\g<1>type            zeroGradient;', flags=re.S)


def add_pimple(d, n_outer=5, n_corr=1):
    p = os.path.join(d, 'system', 'fvSolution')
    t = open(p, encoding='utf-8').read()
    if 'PIMPLE' in t:
        return
    open(p, 'w', encoding='utf-8').write(
        t + f'\nPIMPLE\n{{\n    nOuterCorrectors {n_outer};\n'
            f'    nCorrectors      {n_corr};\n'
            f'    nNonOrthogonalCorrectors 0;\n}}\n')


def append_fos(d, extra_fos: str):
    """把 functionObjects 段插到 controlDict 的 functions{ 之后。"""
    p = os.path.join(d, 'system', 'controlDict')
    _sub(p, r'functions\s*\n\{', 'functions\n{\n' + extra_fos)


# ======================================================================
#  求解 + 结果解析
# ======================================================================
def solve(d, solver=SOLVER, timeout=1200):
    log = os.path.join(d, f'{solver.replace(" ", "_")}.log')
    try:
        rc, out = of_run(d, solver, timeout=timeout, log_file=log, echo=False)
    except subprocess.TimeoutExpired:
        rc, out = -9, 'TIMEOUT'
    return rc, out, log


def summarize(rc, out):
    times = re.findall(r'Time = ([\d.eE+-]+)', out)
    h2o = re.findall(
        r'Solving for H2O\.\w+, Initial residual = ([\d.eE+-]+), Final residual = ([\d.eE+-]+)', out)
    eg = re.findall(
        r'Solving for e\.gas, Initial residual = ([\d.eE+-]+), Final residual = ([\d.eE+-]+)', out)
    Ts = [float(x) for x in re.findall(r'T = ([\d.]+)', out)]
    # 残差统一转 float（正则捕获的是字符串，直接格式化会报 'g' 不支持 str）
    h2o = [(float(a), float(b)) for a, b in h2o]
    eg = [(float(a), float(b)) for a, b in eg]
    return dict(
        rc=rc,
        steps=len(times),
        last_time=times[-1] if times else None,
        janaf=out.count('out of temperature range'),
        minT=min(Ts) if Ts else None,
        h2o=h2o,
        egas=eg,
        fatal=(re.search(r'--> FOAM FATAL[^\n]*\n\s*(.{0,160})', out) or [None, None])[1],
    )


def print_row(name, s):
    h = f"{s['h2o'][-1][1]:.3g}" if s['h2o'] else '  --  '
    e = f"{s['egas'][-1][1]:.3g}" if s['egas'] else '  --  '
    print(f"{name:22s} rc={s['rc']:4d} steps={s['steps']:3d} "
          f"janafWarn={s['janaf']:4d} H2O={h:>10s} e.gas={e:>10s} minT={s['minT']}",
          flush=True)
    if s['fatal']:
        print(f"{'':22s} └─ FATAL: {s['fatal'].strip()[:120]}")


# ======================================================================
#  内置探针
# ======================================================================
def probe_baseline(d):
    pass


def probe_wall_adiabatic(d):
    set_wall_adiabatic(d, 'T.gas')
    set_wall_adiabatic(d, 'T.liquid')


def probe_liquid_dense(d):
    set_internal(d, 'alpha.liquid', '1e-3')


def probe_small_dt(d):
    set_dt(d, '1e-5')


def probe_wall_warm(d):
    set_wall_T(d, '300', 'T.gas')
    set_wall_T(d, '300', 'T.liquid')


def probe_pimple(d):
    add_pimple(d, 5)


def probe_every_step_fo(d):
    """每步都算 FO 并落盘 —— 崩溃前拿到唯一可靠的数据源。"""
    p = os.path.join(d, 'system', 'controlDict')
    _sub(p, r'writeControl\s+adjustableRunTime;', 'writeControl    timeStep;')
    _sub(p, r'writeInterval   [\d.eE+-]+;', 'writeInterval   1;')
    _sub(p, r'writeControl    writeTime;', 'writeControl    timeStep;')
    append_fos(d, """
    statAlphaT
    {
        type            fieldMinMax;
        libs            ("libfieldFunctionObjects.so");
        writeControl    timeStep;
        writeFields     no;
        regionType      all;
        fields          ( alpha.liquid alpha.gas T.gas T.liquid );
    }
""")


BUILTIN = {
    'baseline':    (probe_baseline,      '原样基线（用于对照）'),
    'wall_adiab':  (probe_wall_adiabatic, '冷壁改绝热 —— 判别是否热冲击'),
    'liquid_dense': (probe_liquid_dense, 'alpha.liquid 1e-6→1e-3 —— 判别是否液相稀疏'),
    'small_dt':    (probe_small_dt,     'dt 压到 1e-5 且不自适应 —— 判别是否时间步'),
    'wall_warm':   (probe_wall_warm,     '冷壁 278→300 K —— 判别是否热冲击幅度'),
    'pimple':      (probe_pimple,       '补 PIMPLE 段(nOuterCorrectors 5)'),
    'every_step':  (probe_every_step_fo, '每步落盘 FO —— 崩溃前抓场数据'),
}


def dump_fo(d):
    pp = os.path.join(d, 'postProcessing')
    if not os.path.isdir(pp):
        print('  (无 postProcessing)')
        return
    for fo in sorted(os.listdir(pp)):
        for root, _, files in os.walk(os.path.join(pp, fo)):
            for fn in sorted(files):
                if fn.endswith('.dat'):
                    print(f'\n--- {fo}/{os.path.relpath(root, pp)}/{fn}')
                    txt = open(os.path.join(root, fn), encoding='utf-8',
                               errors='replace').read()
                    print('\n'.join(txt.splitlines()[-8:]))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('probes', nargs='*', help='要跑的探针名')
    ap.add_argument('--list', action='store_true', help='列出内置探针')
    ap.add_argument('--endtime', default='0.5')
    ap.add_argument('--timeout', type=int, default=1200)
    args = ap.parse_args()

    if args.list or not args.probes:
        print('内置探针：')
        for k, (_, desc) in BUILTIN.items():
            print(f'  {k:14s} {desc}')
        print('\n用法：python tools/ablate.py baseline pimple every_step')
        return 0

    if not os.path.isdir(os.path.join(BASE, 'constant', 'polyMesh')):
        raise SystemExit(f'base 算例不存在或无网格: {BASE}')

    os.makedirs(PROBE_ROOT, exist_ok=True)
    print(f'{"探针":22s} {"rc":>4s} {"steps":>5s} {"janaf":>6s}  '
          f'{"H2O末残差":>10s} {"e.gas末残差":>10s}  minT')
    print('-' * 88)
    for name in args.probes:
        if name not in BUILTIN:
            print(f'未知探针 {name}（--list 查看）')
            continue
        fn, _ = BUILTIN[name]
        d = fresh(name, endtime=args.endtime)
        fn(d)
        rc, out, log = solve(d, timeout=args.timeout)
        s = summarize(rc, out)
        print_row(name, s)
        if name == 'every_step':
            dump_fo(d)
        print(f'{"":22s}    log: {log}\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())