#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
隔离探针第 3 轮：
  tut_nums : 教程对齐配置 + 整套教程 fvSchemes/fvSolution（判别是否数值设置问题）
  fo_minmax: 教程对齐配置 + fieldMinMax 每步落盘（抓 alpha/T 越界的单元位置，
             与重复面涉及的单元 20474/22767 对照）
用法:
  python _num_probe.py tut_nums fo_minmax --endtime 0.05
"""
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'tools'))
import _liq_probe as LP  # noqa: E402
from ablate import fresh, solve, summarize, print_row, dump_fo, append_fos  # noqa: E402

TUT = ('D:/APP/openfoam/OpenFOAM-8/tutorials/multiphase/multiphaseEulerFoam/'
       'laminar/bubbleColumnEvaporatingDissolving')


def aligned_config(d):
    """教程对齐的相配置（液相多组分 + Henry + 气相 eConst）。"""
    LP.make_liq_multicomponent(d)
    open(os.path.join(d, 'constant', 'thermophysicalProperties.gas'),
         'w', encoding='utf-8').write(LP.GAS_ECONST)


def probe_tut_nums(d):
    aligned_config(d)
    for fn in ('fvSchemes', 'fvSolution'):
        shutil.copy(os.path.join(TUT, 'system', fn),
                    os.path.join(d, 'system', fn))


def probe_fo_minmax(d):
    aligned_config(d)
    p = os.path.join(d, 'system', 'controlDict')
    # 每步落盘 + fieldMinMax 抓 min/max 位置
    t = open(p, encoding='utf-8').read()
    t = t.replace('writeControl    adjustableRunTime;',
                  'writeControl    timeStep;')
    t = t.replace('writeInterval   0.5;', 'writeInterval   1;')
    open(p, 'w', encoding='utf-8').write(t)
    append_fos(d, """
    statAlphaT
    {
        type            fieldMinMax;
        libs            ("libfieldFunctionObjects.so");
        writeControl    timeStep;
        writeFields     no;
        regionType      all;
        fields          ( alpha.liquid alpha.gas T.gas T.liquid U.gas U.liquid );
    }
""")


PROBES = {
    'tut_nums': probe_tut_nums,
    'fo_minmax': probe_fo_minmax,
}


def main():
    argv = sys.argv[1:]
    names = [a for a in argv if not a.startswith('--')]
    endtime = '0.05'
    for i, a in enumerate(argv):
        if a.startswith('--endtime='):
            endtime = a.split('=', 1)[1]
        elif a == '--endtime' and i + 1 < len(argv):
            endtime = argv[i + 1]
            names = [n for n in names if n != argv[i + 1]]
    if not names:
        print('用法: python _num_probe.py tut_nums fo_minmax [--endtime=0.05]')
        return 1
    for name in names:
        fn = PROBES[name]
        d = fresh(name, endtime=endtime)
        fn(d)
        rc, out, log = solve(d, timeout=900)
        s = summarize(rc, out)
        print_row(name, s)
        if name == 'fo_minmax':
            dump_fo(d)
        print(f'{"":22s}    log: {log}\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
