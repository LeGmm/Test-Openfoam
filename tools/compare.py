#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
 跨结构「冷凝效率」对比：从各算例的 postProcessing 里取指标，排名并输出表格
================================================================================
 用法（装好 OpenFOAM 并跑完算例后）：
   python tools/compare.py                    # 对比 cases/ 下全部结构
   python tools/compare.py v8_1 v8_8          # 只对比指定结构

 指标（对应 Fluent Report Definition）：
   mdot_cond  出口凝液质量流率 kg/s   ← 核心指标①（越大越好）
   Y_H2O_out  出口水蒸气质量分数       ← 核心指标②（越小=去湿越彻底）
   Q_wall     冷壁换热量 W            ← 定壁温下热流是输出（越大越好）
   dp         压降 Pa                 ← 结构强化的代价（越小越好）
   mdot_per_dp 单位压降凝液量          ← 综合评价（越大越好）
================================================================================
"""
from __future__ import annotations

import csv
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CASES = os.path.join(ROOT, "cases")

# ---- 与 generate_case.py 顶部参数保持一致 ----
MDOT   = 2.0e-5    # 入口总质量流量 kg/s
Y_H2O  = 0.0289    # 入口水蒸气质量分数
T_IN   = 308.15    # 入口温度 K (35℃)
CP_AIR = 1006.43   # J/(kg·K)
H_FG   = 2.442e6   # 汽化潜热 J/kg

# 入口水蒸气质量流率（凝结量不可能超过它）——用于合理性校核
MDOT_H2O_IN = MDOT * Y_H2O


def read_latest(func_dir):
    """读 postProcessing/<func>/<time>/<file> 里最新时间的最后一个数值。"""
    if not os.path.isdir(func_dir):
        return None
    times = [t for t in os.listdir(func_dir)
             if os.path.isdir(os.path.join(func_dir, t))]
    if not times:
        return None
    # 时间目录按数值排序
    def tval(t):
        try:
            return float(t)
        except ValueError:
            return -1.0
    times.sort(key=tval)
    latest = os.path.join(func_dir, times[-1])
    files = [f for f in os.listdir(latest) if f.endswith(".dat")]
    if not files:
        return None
    path = os.path.join(latest, files[0])
    val = None
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            try:
                val = float(parts[-1])
            except ValueError:
                continue
    return val


def collect(name):
    case = os.path.join(CASES, name)
    pp = os.path.join(case, "postProcessing")
    if not os.path.isdir(pp):
        return None
    # mdot_cond 已经是质量流率（alphaRhoPhi = alpha*rho*phi 是质量通量），不要再乘密度
    mdot = read_latest(os.path.join(pp, "mdot_cond"))
    y = read_latest(os.path.join(pp, "YH2O_out"))
    tout = read_latest(os.path.join(pp, "T_out"))
    pin = read_latest(os.path.join(pp, "p_in"))
    pout = read_latest(os.path.join(pp, "p_out"))

    dp = (pin - pout) if (pin is not None and pout is not None) else None
    # Q_wall 由能量平衡反算（SOP §8.1）：显热降温 + 凝结潜热
    q = None
    if tout is not None:
        q = MDOT * CP_AIR * (T_IN - tout) + (mdot or 0.0) * H_FG
    ratio = (mdot / dp * 1e6) if (mdot not in (None, 0) and dp not in (None, 0)) else None

    # 合理性校核：凝液量不可能超过入口带进来的水蒸气量
    warn = ""
    if mdot is not None and mdot > MDOT_H2O_IN * 1.05:
        warn = "⚠凝液>入口水汽(不合理)"
    elif mdot is not None and mdot < 0:
        warn = "⚠凝液为负"

    return {"struct": name, "mdot_cond": mdot, "Y_H2O_out": y,
            "T_out": tout, "Q_wall": q, "dp": dp,
            "mdot_per_dp": ratio, "check": warn}


def fmt(v, spec="{:.6g}"):
    return spec.format(v) if isinstance(v, float) else "  --  "


def main():
    names = sys.argv[1:] if len(sys.argv) > 1 else sorted(
        d for d in os.listdir(CASES) if os.path.isdir(os.path.join(CASES, d)))
    rows = [r for r in (collect(n) for n in names) if r]

    if not rows:
        print("没有任何结构产出指标。请先跑完算例（postProcessing/ 下有数据）。")
        print(f"检查目录：{CASES}")
        return 1

    # 按凝液质量流率降序排名
    rows.sort(key=lambda r: (r["mdot_cond"] is not None, r["mdot_cond"] or 0), reverse=True)

    print("\n" + "=" * 96)
    print("跨结构冷凝效率对比（按凝液质量流率 ṁ_cond 降序）")
    print("=" * 96)
    print(f"校核基准：入口水蒸气质量流率 ṁ_H2O,in = {MDOT_H2O_IN:.3e} kg/s"
          f"（凝液量不应超过它）")
    print()
    hdr = (f"{'排名':<4}{'结构':<8}{'ṁ_cond(kg/s)':>15}{'Y_H2O_out':>13}"
           f"{'T_out(K)':>11}{'Q_wall(W)':>12}{'Δp(Pa)':>11}{'ṁ/Δp×1e6':>13}{'校核':>10}")
    print(hdr)
    print("-" * 100)
    for i, r in enumerate(rows, 1):
        print(f"{i:<4}{r['struct']:<8}{fmt(r['mdot_cond']):>15}"
              f"{fmt(r['Y_H2O_out']):>13}{fmt(r['T_out']):>11}"
              f"{fmt(r['Q_wall']):>12}{fmt(r['dp']):>11}"
              f"{fmt(r['mdot_per_dp']):>13}{r.get('check',''):>10}")
    print("=" * 100)
    print("说明：Q_wall 由能量平衡反算 = ṁ·cp·(T_in−T_out) + ṁ_cond·h_fg"
          "（SOP §8.1 口径，优于直接读壁面热流）")

    out_csv = os.path.join(ROOT, "comparison.csv")
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"已写出 CSV：{out_csv}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
