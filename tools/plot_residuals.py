#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
  OpenFOAM 求解器日志残差解析 + 残差-迭代曲线绘制
================================================================================
  用法：
    python tools/plot_residuals.py <输出.png> "<图题>" <日志1> [日志2 ...]

  日志按给定顺序拼接解析（同一算例的续跑日志传两个即可）。
  每个 "Time = N" 迭代块内只取每条方程【第一次】出现的 Initial residual，
  避免 PIMPLE 外迭代重复行把曲线打乱。

  方程归组（chtMultiRegionFoam 双区域）：
    流体动量  U    = max(Ux, Uy, Uz)
    流体能量  h
    流体压力  p_rgh
    固体能量  e
  容差参考线取 system/fluid|solid/fvSolution 的 residualControl：
    U 1e-4 / h 1e-6 / p_rgh 1e-5 / e 1e-6（与 caseCHT_v88、caseCHT_qishi 一致）。
================================================================================
"""
from __future__ import annotations

import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

BG, PANEL, LINE, TX, TX2 = "#0f1216", "#161b22", "#2a3340", "#e6edf3", "#9aa7b4"

RE_TIME = re.compile(r"^Time = (\S+)")
RE_SOLVE = re.compile(
    r"Solving for (\w+), Initial residual = (\S+), Final residual = (\S+)")

# 关注方程 → 图上分组/颜色
GROUPS = {
    "流体动量 U": (["Ux", "Uy", "Uz"], "#f0883e", 1e-4),
    "流体能量 h": (["h"], "#58a6ff", 1e-6),
    "流体压力 p_rgh": (["p_rgh"], "#3fb950", 1e-5),
    "固体能量 e": (["e"], "#bc8cff", 1e-6),
}


def parse_logs(paths):
    """返回 (iters, {方程名: [初残差]})。多日志按顺序拼接，迭代号连续累加。"""
    data = {}
    iters = []
    t0 = 0  # 续跑日志的 Time 起点（仅用于占位，实际迭代号用累计计数）
    for path in paths:
        with open(path, encoding="utf-8", errors="replace") as f:
            cur = None          # 当前 Time 号
            seen = set()        # 本迭代块已记录的方程
            for line in f:
                m = RE_TIME.match(line)
                if m:
                    try:
                        cur = float(m.group(1))
                    except ValueError:
                        cur = None
                    seen = set()
                    continue
                if cur is None:
                    continue
                m = RE_SOLVE.search(line)
                if not m:
                    continue
                eq, res = m.group(1), m.group(2)
                if eq in seen:
                    continue
                seen.add(eq)
                data.setdefault(eq, []).append(float(res))
                iters.append(len(iters) + 1 if not iters else iters[-1] + 1)
    # iters 以"方程出现次数"累计太粗，改为按迭代块推进：重新扫一遍更稳
    return data


def parse_logs2(paths):
    """按迭代块推进迭代号：每遇到新 Time 块，迭代号 +1。"""
    data = {}
    iters = []
    step = 0
    for path in paths:
        with open(path, encoding="utf-8", errors="replace") as f:
            cur = None
            seen = set()
            for line in f:
                m = RE_TIME.match(line)
                if m:
                    step += 1
                    cur = True
                    seen = set()
                    continue
                if cur is None:
                    continue
                m = RE_SOLVE.search(line)
                if not m:
                    continue
                eq, res = m.group(1), m.group(2)
                if eq in seen:
                    continue
                seen.add(eq)
                data.setdefault(eq, []).append(float(res))
                iters.append(step)
    # 同一 step 内多条方程共享同一步号 → 去重排序
    return iters, data


def style(ax, title, xl="外迭代步", yl="初始残差"):
    ax.set_facecolor(PANEL)
    for s in ax.spines.values():
        s.set_color(LINE)
    ax.tick_params(colors=TX2, labelsize=9)
    ax.set_title(title, color=TX, fontsize=11)
    ax.set_xlabel(xl, color=TX2)
    ax.set_ylabel(yl, color=TX2)
    ax.set_yscale("log")
    ax.grid(True, which="both", color=LINE, lw=0.4, alpha=0.5)


def main():
    out, title = sys.argv[1], sys.argv[2]
    logs = sys.argv[3:]

    iters, data = parse_logs2(logs)
    n_iter = max(iters) if iters else 0
    x = np.array(sorted(set(iters)))

    fig, axes = plt.subplots(1, 2, figsize=(13.6, 5.4), facecolor=BG)

    # --- 左：流体三方程
    ax = axes[0]
    for gname, (eqs, color, tol) in GROUPS.items():
        if gname == "固体能量 e":
            continue
        if not any(e in data for e in eqs):
            continue
        y = None
        for e in eqs:
            if e not in data:
                continue
            v = np.array(data[e])
            y = v if y is None else np.maximum(y, v)
        ax.plot(x[:len(y)], y, color=color, lw=1.1, label=gname)
        ax.axhline(tol, color=color, lw=0.8, ls=":", alpha=0.7)
        ax.annotate(f"容差 {tol:g}", xy=(x[0], tol), xytext=(2, 3),
                    textcoords="offset points", color=color, fontsize=7.5)
        ax.annotate(f"末值 {y[-1]:.2e}", xy=(x[len(y)-1], y[-1]),
                    xytext=(-70, 6), textcoords="offset points",
                    color=color, fontsize=8)
    style(ax, f"{title} · 流体区残差")
    lg = ax.legend(fontsize=8.5, facecolor=PANEL, edgecolor=LINE, loc="upper right")
    [t.set_color("#c9d6e4") for t in lg.get_texts()]

    # --- 右：固体能量 e
    ax = axes[1]
    eqs, color, tol = GROUPS["固体能量 e"]
    if "e" in data:
        y = np.array(data["e"])
        ax.plot(x[:len(y)], y, color=color, lw=1.1, label="固体能量 e")
        ax.axhline(tol, color=color, lw=0.8, ls=":", alpha=0.7)
        ax.annotate(f"容差 {tol:g}", xy=(x[0], tol), xytext=(2, 3),
                    textcoords="offset points", color=color, fontsize=7.5)
        ax.annotate(f"末值 {y[-1]:.2e}", xy=(x[len(y)-1], y[-1]),
                    xytext=(-70, 6), textcoords="offset points",
                    color=color, fontsize=8)
        style(ax, f"{title} · 固体区残差")
        lg = ax.legend(fontsize=8.5, facecolor=PANEL, edgecolor=LINE, loc="upper right")
        [t.set_color("#c9d6e4") for t in lg.get_texts()]
    else:
        ax.text(0.5, 0.5, "无固体残差记录", ha="center", va="center",
                color=TX2, transform=ax.transAxes)
        style(ax, f"{title} · 固体区残差")

    # 首残差概览写进总标题
    firsts = []
    for gname, (eqs, _, _) in GROUPS.items():
        for e in eqs:
            if e in data and data[e]:
                firsts.append(f"{e}:{data[e][0]:.1e}")
    fig.suptitle(f"{title}　残差-迭代曲线（共 {n_iter} 个外迭代步）",
                 color=TX, fontsize=13)
    fig.text(0.99, 0.01, "  ".join(firsts[:6]), ha="right", color=TX2, fontsize=8)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.savefig(out, dpi=150, facecolor=BG)
    print(f"saved {out}  ({n_iter} steps)")

    # 摘要
    for gname, (eqs, _, _) in GROUPS.items():
        for e in eqs:
            if e in data and data[e]:
                print(f"  {e:8s} 首 {data[e][0]:.3e}  末 {data[e][-1]:.3e}  n={len(data[e])}")


if __name__ == "__main__":
    main()
