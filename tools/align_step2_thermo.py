#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
 把 Step2（multiphaseEulerFoam）算例的热物性/相模型对齐成"能跑"的那一套
================================================================================
 为什么需要它：generate_case.py 生成的是 janaf 气相 + purePhaseModel 液相，
 在相变启动时会 FATAL（thermoI.H: Maximum number of iterations exceeded: 100）。
 实测有效的配方（与官方 bubbleColumnEvaporatingDissolving 教程一致）是：
   * 液相 multiComponentPhaseModel + Henry 溶解 + 多组分热物性
   * 气相 eConst 热物性（避开 janaf 反解温度不收敛）
 这套改动原本散落在 tools/diagnostics/liq_probe.py 里（探针脚本），本工具把它变成
 一条可重复执行的流水线步骤，避免每次重新生成算例后又踩回同一个坑。

 用法：
   python tools/align_step2_thermo.py case_qishi
 幂等：已经对齐过的算例会直接跳过，不会重复修改。
================================================================================
"""
from __future__ import annotations

import importlib
import os
import sys

TOOLS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TOOLS)


def main() -> int:
    if len(sys.argv) < 2:
        print(__doc__)
        return 1
    case = sys.argv[1].rstrip("\\/")

    phase = os.path.join(case, "constant", "phaseProperties")
    gas = os.path.join(case, "constant", "thermophysicalProperties.gas")
    liq = os.path.join(case, "constant", "thermophysicalProperties.liquid")
    for p in (phase, gas, liq):
        if not os.path.isfile(p):
            print(f"缺文件：{p}\n（先用 tools/generate_case.py 生成算例）")
            return 2

    phase_txt = open(phase, encoding="utf-8").read()
    gas_txt = open(gas, encoding="utf-8").read()
    if "multiComponentPhaseModel" in phase_txt and "thermo          eConst" in gas_txt:
        print("已经对齐过（液相 multiComponent + 气相 eConst），跳过。")
        return 0

    try:
        lp = importlib.import_module("diagnostics.liq_probe")
        lp.probe_liq_mc_gasconst(case)
    except Exception as e:                                   # 归档脚本被挪动/改动时的兜底
        print(f"调用 diagnostics/liq_probe.py 失败：{type(e).__name__}: {e}")
        return 3

    print(f"已对齐 {case}：")
    print("  液相 phaseProperties → multiComponentPhaseModel + Henry 溶解（air）")
    print("  thermophysicalProperties.liquid → 多组分")
    print("  thermophysicalProperties.gas    → eConst（不用 janaf）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
