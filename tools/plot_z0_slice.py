#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
 z=0（xy）平面的温度 / 流场图：读 OpenFOAM surfaces FO 切出的 VTK
================================================================================
 先在算例目录里切面：
   postProcess -region fluid -func surfaces -latestTime   → postProcessing/surfaces/<t>/{T,U,p_rgh}_z0.vtk
   postProcess -region solid -func surfaces -latestTime
 ⚠️ 多区域算例的 postProcessing 路径**不带 region 前缀**，两个区域互相覆盖，
    跑完一个必须立刻拷走改名。

 用法：
   python tools/plot_z0_slice.py fluid_T.vtk fluid_U.vtk out.png [solid_T.vtk]

 实现说明：cuttingPlane 输出的是 4~11 边多边形，扇形三角化会产生退化三角形
 （matplotlib 的 TriFinder 会报 "Triangulation is invalid"），所以这里改用
 scipy griddata 把散点插到规则网格再画。
================================================================================
"""
from __future__ import annotations

import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

BG, PANEL, LINE, TX, TX2 = "#0f1216", "#161b22", "#2a3340", "#e6edf3", "#9aa7b4"


def read_vtk(path):
    """OpenFOAM surfaces FO 写 legacy VTK 的 FIELD 格式：
         FIELD attributes 1
         U 3 4479 float
    """
    tok = open(path, encoding="utf-8", errors="replace").read().split()
    i, n = 0, len(tok)
    pts, polys, fields = None, None, {}
    while i < n:
        k = tok[i]
        if k == "POINTS":
            m = int(tok[i + 1])
            pts = np.array([float(x) for x in tok[i + 3: i + 3 + 3 * m]]).reshape(m, 3)
            i += 3 + 3 * m
        elif k == "POLYGONS":
            size = int(tok[i + 2])
            raw = [int(x) for x in tok[i + 3: i + 3 + size]]
            polys, j = [], 0
            while j < len(raw):
                cnt = raw[j]
                polys.append(raw[j + 1: j + 1 + cnt])
                j += 1 + cnt
            i += 3 + size
        elif k == "FIELD":
            nattr = int(tok[i + 2])
            j = i + 3
            for _ in range(nattr):
                name, ncomp, ntuple = tok[j], int(tok[j + 1]), int(tok[j + 2])
                j += 4
                vals = np.array([float(x) for x in tok[j: j + ncomp * ntuple]])
                fields[name] = vals.reshape(ntuple, ncomp) if ncomp > 1 else vals
                j += ncomp * ntuple
            i = j
        else:
            i += 1
    return pts, polys, fields


def style(ax, title, xl="x (mm)", yl="y (mm)"):
    ax.set_facecolor(PANEL)
    for s in ax.spines.values():
        s.set_color(LINE)
    ax.tick_params(colors=TX2, labelsize=9)
    ax.set_title(title, color=TX, fontsize=11)
    ax.set_xlabel(xl, color=TX2)
    ax.set_ylabel(yl, color=TX2)
    ax.set_aspect("equal")


def to_grid(X, Y, values, nx=320, ny=340):
    gx = np.linspace(X.min(), X.max(), nx)
    gy = np.linspace(Y.min(), Y.max(), ny)
    GX, GY = np.meshgrid(gx, gy)
    V = griddata((X, Y), values, (GX, GY), method="linear")
    Vn = griddata((X, Y), values, (GX, GY), method="nearest")
    return gx, gy, GX, GY, np.where(np.isnan(V), Vn, V)


def main():
    fT, fU, out = sys.argv[1:4]
    sT = sys.argv[4] if len(sys.argv) > 4 else None

    pt, _, fld = read_vtk(fT)
    X, Y = pt[:, 0] * 1e3, pt[:, 1] * 1e3
    T = fld["T"]

    _, _, fldU = read_vtk(fU)
    U = fldU["U"]
    mag = np.linalg.norm(U, axis=1)

    panels = 3 if sT else 2
    fig, axes = plt.subplots(1, panels, figsize=(6.4 * panels, 5.6), facecolor=BG)
    axes = np.atleast_1d(axes)

    # --- 流体温度
    ax = axes[0]
    gx, gy, GX, GY, GT = to_grid(X, Y, T)
    lv = np.linspace(np.nanpercentile(GT, 0.5), np.nanpercentile(GT, 99.5), 28)
    cf = ax.contourf(GX, GY, GT, levels=lv, cmap="coolwarm", extend="both")
    cb = fig.colorbar(cf, ax=ax, shrink=0.86, pad=0.02)
    cb.set_label("T (K)", color=TX2)
    cb.ax.tick_params(colors=TX2, labelsize=8)
    style(ax, "流体温度 T（z=0）")
    ax.plot([12], [29], marker="v", color="#3fb950", ms=9, ls="none", label="inlet x=+12")
    ax.plot([-17.5], [29], marker="^", color="#f0883e", ms=9, ls="none", label="outlet x=−17.5")
    lg = ax.legend(fontsize=8, facecolor=PANEL, edgecolor=LINE, loc="lower left")
    [t.set_color("#c9d6e4") for t in lg.get_texts()]

    # --- 速度幅值 + 流线
    ax = axes[1]
    gx, gy, GX, GY, GM = to_grid(X, Y, mag)
    cf = ax.contourf(GX, GY, GM, levels=np.linspace(0, np.nanpercentile(GM, 99), 28), cmap="viridis")
    cb = fig.colorbar(cf, ax=ax, shrink=0.86, pad=0.02)
    cb.set_label("|U| (m/s)", color=TX2)
    cb.ax.tick_params(colors=TX2, labelsize=8)
    iU = griddata((X, Y), U[:, 0], (GX, GY), method="linear")
    iV = griddata((X, Y), U[:, 1], (GX, GY), method="linear")
    UU, VV = np.nan_to_num(iU), np.nan_to_num(iV)
    ax.streamplot(gx, gy, UU, VV, density=2.0, color="#e6edf3", linewidth=0.65,
                  arrowsize=0.7, maxlength=5.0)
    style(ax, "流体速度 |U| 与流线（z=0）")

    if sT:
        sp, _, fldS = read_vtk(sT)
        sX, sY = sp[:, 0] * 1e3, sp[:, 1] * 1e3
        gx, gy, GX, GY, GS = to_grid(sX, sY, fldS["T"])
        ax = axes[2]
        cf = ax.contourf(GX, GY, GS, levels=np.linspace(np.nanpercentile(GS, 0.5),
                                                        np.nanpercentile(GS, 99.5), 28),
                         cmap="coolwarm", extend="both")
        cb = fig.colorbar(cf, ax=ax, shrink=0.86, pad=0.02)
        cb.set_label("T (K)", color=TX2)
        cb.ax.tick_params(colors=TX2, labelsize=8)
        style(ax, "固体温度 T（z=0）")

    plt.tight_layout()
    plt.savefig(out, dpi=150, facecolor=BG)
    print(f"saved {out}")
    print(f"流体切片 {len(X)} 点；T {np.nanmin(T):.2f}~{np.nanmax(T):.2f} K；"
          f"|U| max {mag.max():.4g} / 中位 {np.median(mag):.4g} m/s")


if __name__ == "__main__":
    main()
