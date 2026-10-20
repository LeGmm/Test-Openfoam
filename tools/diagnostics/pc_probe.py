#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
隔离探针第 2 轮：在「教程对齐配置」（liq_mc_gasconst）基础上定位 dmdtf 爆炸。
  no_pc        : 关闭相变（interfaceComposition 两侧置空）—— 判别是否传质源项
  no_pc_janaf  : no_pc + 恢复原气相 janaf 热物性 —— 判别 janaf 本身是否有问题
  nalpha1      : 保留相变，nAlphaCorr 2→1 —— 判别是否第 2 次 MULES 修正越界
用法:
  python _pc_probe.py no_pc no_pc_janaf nalpha1 --endtime 0.05
"""
import os
import re
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'tools'))
import _liq_probe as LP  # noqa: E402  (复用教程对齐配置)
from ablate import fresh, solve, summarize, print_row  # noqa: E402


def disable_phase_change(d):
    pp = os.path.join(d, 'constant', 'phaseProperties')
    t = open(pp, encoding='utf-8').read()
    # 置空 gas 侧饱和模型（注意：不能用 .*? 到第一个 ');'，会撞上 'species ( H2O );'
    # 必须匹配到行首的 '\n);' 才是块结尾）
    t2 = re.sub(
        r'interfaceComposition\.gas\s*\n\(.*?\n\);',
        'interfaceComposition.gas\n();',
        t, flags=re.S)
    assert t2 != t, 'interfaceComposition.gas not matched'
    t = t2
    # 置空 liquid 侧 Henry 模型
    t2 = re.sub(
        r'interfaceComposition\.liquid\s*\n\(.*?\n\);',
        'interfaceComposition.liquid\n();',
        t, flags=re.S)
    assert t2 != t, 'interfaceComposition.liquid not matched'
    open(pp, 'w', encoding='utf-8').write(t2)


def restore_gas_janaf(d):
    src = os.path.join(HERE, 'case', 'constant', 'thermophysicalProperties.gas')
    shutil.copy(src, os.path.join(d, 'constant', 'thermophysicalProperties.gas'))


def set_nalpha(d, n):
    p = os.path.join(d, 'system', 'fvSolution')
    t = open(p, encoding='utf-8').read()
    t2 = t.replace('nAlphaCorr      2;', f'nAlphaCorr      {n};')
    assert t2 != t
    open(p, 'w', encoding='utf-8').write(t2)


def probe_no_pc(d):
    LP.make_liq_multicomponent(d)
    open(os.path.join(d, 'constant', 'thermophysicalProperties.gas'),
         'w', encoding='utf-8').write(LP.GAS_ECONST)
    disable_phase_change(d)


def probe_no_pc_janaf(d):
    LP.make_liq_multicomponent(d)
    disable_phase_change(d)
    restore_gas_janaf(d)


def probe_nalpha1(d):
    LP.make_liq_multicomponent(d)
    open(os.path.join(d, 'constant', 'thermophysicalProperties.gas'),
         'w', encoding='utf-8').write(LP.GAS_ECONST)
    set_nalpha(d, 1)


PROBES = {
    'no_pc': probe_no_pc,
    'no_pc_janaf': probe_no_pc_janaf,
    'nalpha1': probe_nalpha1,
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
        print('用法: python _pc_probe.py no_pc no_pc_janaf nalpha1 [--endtime=0.05]')
        return 1
    for name in names:
        fn = PROBES[name]
        d = fresh(name, endtime=endtime)
        fn(d)
        rc, out, log = solve(d, timeout=900)
        s = summarize(rc, out)
        print_row(name, s)
        print(f'{"":22s}    log: {log}\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
