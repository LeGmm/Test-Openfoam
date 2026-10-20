#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""结构冷凝对比的量化指标提取（稳态口径，纯后处理，不再跑瞬态）。

──────────────────────────────────────────────────────────────────────────────
为什么不用瞬态多相（Step2）
──────────────────────────────────────────────────────────────────────────────
Step2（multiphaseEulerFoam + 界面组分相变）实测单步 0.949 s / 8 进程、dt≈2.93e-5，
跑满 1τ(2.97 s) 需 ≈27 h/结构，且当前在 t=0.0467 s 因极小 α_l 单元液相温度发散而崩。
本机负荷与工期都不允许，因此结构对比改用【稳态 + 传质控制冷凝模型】：

  ① 稳态共轭传热（chtMultiRegionFoam，ddt=steadyState）给出
     界面温度分布 T_w(x)、壁面热流 q''(x)、流场与进出口状态；
  ② 冷凝速率用 Colburn–Hougen / Chilton–Colburn 传质类比求：
     壁面被液膜润湿 → 壁面气相水蒸气质量分数 = 饱和值 Y_sat(T_w)；
     水蒸气穿过含不凝气体（空气）的边界层扩散到壁面 → Stefan 流修正：

         ṁ'' = ρ · h_m · ln[(1 − Y_sat) / (1 − Y_b)]
         h_m = h / (ρ · cp · Le^{2/3})        （Chilton–Colburn，Le = α/D ≈ 0.857）
         h   = −q'' / (T_ref − T_w)           （由 CFD 的壁面热流反算）

     出口含湿量由干空气质量守恒 + 冷凝量闭合：
         ṁ_cond = ṁ_da · (ω_in − ω_out)，  ω = Y/(1−Y) [kg 汽 / kg 干空气]
     用对数平均传质推动力 ΔY_lm 迭代求解 ω_out（固定点迭代，收敛快）。

适用条件（不满足即结论不成立，见报告「不确定性来源」）：
  * 壁面完全润湿、液膜热阻可忽略 → 给出的是【传质控制上限】；
  * 不凝气体（空气）扩散阻力主导（本工况 Y_H2O≈0.029，成立）；
  * 忽略潜热对固体温度场的反馈（Step1 是单相干空气，单向耦合）；
  * 网格尺度须在同一量级，否则 h 与 ṁ 不可直接横比（本工程 v8_1/v8_8 差 7 倍，已披露）。

用法：
    python tools/condensation_metrics.py caseCHT_qishi 15000
    python tools/condensation_metrics.py caseCHT_v88 9000 --json _m_v88.json
    python tools/condensation_metrics.py A tA B tB --compare     # 同口径对比
"""
from __future__ import annotations

import argparse
import json
import math
import os
import re
import sys

import numpy as np

try:
    from ofmesh import read_mesh, face_area_centre
except ImportError:                                   # 允许直接 python tools/xxx.py 调用
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from ofmesh import read_mesh, face_area_centre

# ── 物性与工况常数（与 tools/generate_case.py 顶部保持一致）────────────────
P_OP = 101325.0        # Pa
Y_H2O_IN = 0.0289      # 入口水蒸气质量分数（35℃ / 80%RH）
LE = 0.857             # Le = α/D（水蒸气-空气），Chilton–Colburn  exponent 2/3
MW_RATIO = 0.62198     # M_H2O / M_air


def p_sat(T: np.ndarray) -> np.ndarray:
    """Arden Buck 饱和水蒸气压 [Pa]，T in K。"""
    Tc = T - 273.15
    return 611.21 * np.exp((18.678 - Tc / 234.5) * Tc / (257.14 + Tc))


def y_sat(T: np.ndarray, p: float = P_OP) -> np.ndarray:
    """饱和水蒸气质量分数。"""
    ps = p_sat(np.asarray(T, dtype=float))
    return MW_RATIO * ps / np.maximum(p - ps, 1.0)


def t_dew(Y: float, p: float = P_OP) -> float:
    """由质量分数反算露点温度 [K]（二分）。"""
    if Y <= 0:
        return 200.0
    lo, hi = 230.0, 360.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if y_sat(mid, p) < Y:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# ── OpenFOAM 场读取 ────────────────────────────────────────────────────────
_NUM = r"([-\d.eE+]+)"


def _boundary_block(txt: str, patch: str) -> str:
    """截取 boundaryField 里某个补丁块（到下一个同名缩进的补丁块为止）。"""
    i = txt.find("boundaryField")
    if i < 0:
        raise ValueError("场文件里找不到 boundaryField")
    txt = txt[i:]
    m = re.search(r"^\s{4}" + re.escape(patch) + r"\s*\n\s*\{", txt, re.M)
    if not m:
        raise KeyError(f"补丁 {patch} 不在 boundaryField 里")
    start = m.end()
    depth = 1
    j = start
    while depth and j < len(txt):
        if txt[j] == "{":
            depth += 1
        elif txt[j] == "}":
            depth -= 1
        j += 1
    return txt[start:j]


def read_patch_scalar(path: str, patch: str) -> np.ndarray:
    """读 volScalarField / surfaceScalarField 在 patch 上的值（uniform 或非均匀）。"""
    txt = open(path, encoding="utf-8", errors="replace").read()
    blk = _boundary_block(txt, patch)
    m = re.search(r"value\s+nonuniform\s+List<scalar>\s*\n\s*(\d+)\s*\n\((.*?)\n\s*\)",
                  blk, re.S)
    if m:
        n = int(m.group(1))
        v = np.fromstring(m.group(2).replace("\n", " "), sep=" ")
        if len(v) != n:
            raise ValueError(f"{path}:{patch} 声明 {n} 个值，实读 {len(v)}")
        return v
    m = re.search(r"value\s+uniform\s+" + _NUM, blk)
    if m:
        return np.full(1, float(m.group(1)))
    m = re.search(r"uniform\s+" + _NUM, blk)          # 只有 uniform 的简写
    if m:
        return np.full(1, float(m.group(1)))
    raise ValueError(f"{path}:{patch} 读不到值")


def read_internal(path: str) -> np.ndarray:
    txt = open(path, encoding="utf-8", errors="replace").read()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*\n\s*(\d+)\s*\n\((.*?)\n\s*\)",
                  txt, re.S)
    if m:
        return np.fromstring(m.group(2).replace("\n", " "), sep=" ")
    m = re.search(r"internalField\s+uniform\s+" + _NUM, txt)
    if m:
        return None, float(m.group(1))
    raise ValueError(f"{path} 读不到 internalField")


def read_thermo(case: str, region: str = "fluid") -> dict:
    """从 constant/<region>/thermophysicalProperties 取 rho / Cp / Pr。"""
    p = os.path.join(case, "constant", region, "thermophysicalProperties")
    txt = open(p, encoding="utf-8", errors="replace").read()
    out = {}
    for key, pat in (("rho", r"rho\s+" + _NUM), ("Cp", r"Cp\s+" + _NUM),
                     ("Pr", r"Pr\s+" + _NUM), ("mu", r"mu\s+" + _NUM)):
        m = re.search(pat, txt)
        if m:
            out[key] = float(m.group(1))
    return out


# ── 主计算 ─────────────────────────────────────────────────────────────────
def compute(case: str, time: str, patch: str = "interface",
            y_in: float = Y_H2O_IN, verbose: bool = True) -> dict:
    case = case.rstrip("\\/")
    mesh = os.path.join(case, "constant", "fluid", "polyMesh")
    tdir = os.path.join(case, str(time), "fluid")

    pts, faces, own, nei, patches = read_mesh(mesh)
    if patch not in patches:
        raise KeyError(f"{case} 的 fluid 区没有补丁 {patch}；可选：{list(patches)}")
    b = patches[patch]
    i0, nf = b["startFace"], b["nFaces"]
    idx = list(range(i0, i0 + nf))
    area = np.array([face_area_centre(pts, faces[i])[0] for i in idx])
    owner = own[i0:i0 + nf]                     # 壁面相邻单元

    q = read_patch_scalar(os.path.join(tdir, "wallHeatFlux"), patch)   # W/m²
    Tw = read_patch_scalar(os.path.join(tdir, "T"), patch)             # K
    if len(Tw) == 1:
        Tw = np.full(nf, Tw[0])
    if len(q) == 1:
        q = np.full(nf, q[0])

    Tcell_raw = read_internal(os.path.join(tdir, "T"))
    if isinstance(Tcell_raw, tuple):
        Tcell = np.full(int(nei.max()) + 1, Tcell_raw[1])
    else:
        Tcell = Tcell_raw
    Tc = Tcell[owner]

    rho_f = read_internal(os.path.join(tdir, "rho"))
    rho = float(np.mean(rho_f)) if not isinstance(rho_f, tuple) else rho_f[1]
    th = read_thermo(case)
    cp = th.get("Cp", 1006.0)

    # 体积加权主流温度（作为整体参考温度）
    Tbulk = float(np.mean(Tcell))

    # ── 换热系数：局部以壁面相邻单元温度为参考，整体以体积平均温度为参考
    dT_loc = Tc - Tw
    dT_glb = Tbulk - Tw
    with np.errstate(divide="ignore", invalid="ignore"):
        h_loc = np.where(dT_loc > 1e-3, -q / np.maximum(dT_loc, 1e-9), np.nan)
        # 只在「壁面确实在冷却气体」的面给 h：q<0 且局部温差 >0.5 K。
        # 否则壁面比主流还热（回流/远离冷板区）会给出 ±1e4 的无意义值。
        ok = (q < 0) & (dT_glb > 0.5)
        h_glb = np.where(ok, -q / np.maximum(dT_glb, 1e-9), np.nan)
    Q = float(np.sum(q * area))                     # W（负 = 流体放热）

    # ── Colburn–Hougen 冷凝速率（固定点迭代求出口含湿量）
    A = float(area.sum())
    Ys = y_sat(Tw, P_OP)
    w_in = y_in / (1.0 - y_in)
    m_da = None            # 干空气质量流率，由入口总流量×工况给出，见下方 m_total
    m_total = None
    # 入口总质量流率：优先读 postProcessing/m_in
    mi = os.path.join(case, "postProcessing", "fluid", "m_in")
    if os.path.isdir(mi):
        latest = sorted(os.listdir(mi), key=lambda s: int(s) if s.isdigit() else 0)[-1]
        f = os.path.join(mi, latest, "surfaceFieldValue.dat")
        if os.path.isfile(f):
            lines = [l for l in open(f, encoding="utf-8", errors="replace") if l.strip()]
            m_total = abs(float(lines[-1].split()[-1]))
    if m_total is None:
        m_total = 2.0e-5                             # generate_case.py 默认 MDOT
    m_da = m_total * (1.0 - y_in)

    # 出口质量加权温度（能量平衡判据与 UA 都要用）
    Tout_mw = None
    fo = os.path.join(case, "postProcessing", "fluid", "T_out_mw")
    if os.path.isdir(fo):
        latest = sorted(os.listdir(fo), key=lambda s: int(s) if s.isdigit() else 0)[-1]
        f = os.path.join(fo, latest, "surfaceFieldValue.dat")
        if os.path.isfile(f):
            lines = [l for l in open(f, encoding="utf-8", errors="replace") if l.strip()]
            Tout_mw = float(lines[-1].split()[-1])

    h_bar = float(np.nansum(h_glb * area) / np.nansum(area[np.isfinite(h_glb)]))
    hm_bar = h_bar / (rho * cp * LE ** (2.0 / 3.0))  # m/s

    # ── 总热导 UA：由 CFD 的 Q 与对数平均温差反算（结构差异的直接量化量）
    Tw_bar = float(np.average(Tw, weights=area))
    T_out_dry = Tout_mw if Tout_mw is not None else float(np.mean(Tcell))
    d1, d2 = 308.15 - Tw_bar, T_out_dry - Tw_bar
    dT_lm = (d1 - d2) / math.log(d1 / d2) if (d1 > 0 and d2 > 0 and abs(d1 - d2) > 1e-9) \
        else max(0.5 * (d1 + d2), 1e-6)
    # 焓平衡热负荷（守恒、不依赖近壁梯度）：Q_flux = ṁ·cp·(T_in − T_out,通量加权)
    Q_flux = m_total * cp * (308.15 - T_out_dry) if Tout_mw is not None else abs(Q)
    UA = Q_flux / dT_lm                       # W/K（用守恒口径，不用近壁梯度积分）
    Q_vol = m_total / rho                     # m³/s
    hm_A = hm_bar * A                         # m³/s（传质"通风量"）
    NTU_m = hm_A / Q_vol                      # 传质传递单元数

    # ── 湿工况耦合求解：UA 由 CFD 给出；出口按平衡（NTU_m>>1）取饱和
    #    能量：UA·ΔT_lm = ṁ_da[h(T_in,ω_in) − h(T_out,ω_out)] − ṁ_cond·cp_l·T_w
    CP_A, CP_V, CP_L, H_FG0 = 1006.0, 1860.0, 4218.0, 2.501e6

    def h_moist(T_C: float, w: float) -> float:
        return CP_A * T_C + w * (H_FG0 + CP_V * T_C)

    def resid(T_out: float) -> float:
        w_s = float(y_sat(T_out, P_OP))
        w_s = w_s / (1.0 - w_s)
        w_s = min(w_s, w_in)
        m_c = m_da * (w_in - w_s)
        dd1, dd2 = (308.15 - 273.15) - (Tw_bar - 273.15), (T_out - 273.15) - (Tw_bar - 273.15)
        if dd1 <= 0 or dd2 <= 1e-6:
            return 1e9
        lm = (dd1 - dd2) / math.log(dd1 / dd2) if abs(dd1 - dd2) > 1e-9 \
            else 0.5 * (dd1 + dd2)
        q_wall = UA * lm
        q_gas = m_da * (h_moist(308.15 - 273.15, w_in) - h_moist(T_out - 273.15, w_s)) \
            - m_c * CP_L * (Tw_bar - 273.15)
        return q_wall - q_gas

    lo, hi = Tw_bar + 1e-3, 308.15
    flo, fhi = resid(lo), resid(hi)
    if flo * fhi > 0:                          # 无根：退化为纯显热解
        T_out_wet = T_out_dry
    else:
        for _ in range(120):
            mid = 0.5 * (lo + hi)
            if resid(mid) * flo > 0:
                lo, flo = mid, resid(mid)
            else:
                hi = mid
        T_out_wet = 0.5 * (lo + hi)

    y_out_sat = float(y_sat(T_out_wet, P_OP))
    w_out = y_out_sat / (1.0 - y_out_sat)
    w_out = min(w_out, w_in)
    m_cond = max(m_da * (w_in - w_out), 0.0)
    Y_out = w_out / (1.0 + w_out)
    m_vap_in = m_da * w_in
    eta = m_cond / m_vap_in if m_vap_in > 0 else float("nan")

    # 传质能力交叉校核：在平衡出口状态下，传质方程能提供的最大速率
    drive_in = max(y_in - float(np.mean(Ys)), 1e-12)
    drive_out = max(Y_out - float(np.mean(Ys)), 1e-12)
    dY_lm = (drive_in - drive_out) / math.log(drive_in / drive_out) \
        if abs(drive_in - drive_out) > 1e-12 else drive_in
    m_cond_transfer = rho * hm_A * dY_lm
    transfer_limited = bool(m_cond_transfer < m_cond)
    m_cond_max = max(m_da * (w_in - y_sat(T_out_dry, P_OP) / (1 - y_sat(T_out_dry, P_OP))), 0.0)
    T_out = T_out_wet

    # ── 能量平衡（Step1 口径）：Q vs ṁ·cp·ΔT
    ratio = None
    if Tout_mw is not None:
        ratio = abs(Q) / (m_total * cp * (308.15 - Tout_mw))

    def stats(a, w=None):
        a = np.asarray(a, dtype=float)
        m = np.isfinite(a)
        if not m.any():
            return dict(mean=float("nan"), min=float("nan"), max=float("nan"),
                        std=float("nan"), p05=float("nan"), p95=float("nan"))
        v = a[m]
        ww = None if w is None else np.asarray(w)[m]
        mean = float(np.average(v, weights=ww)) if ww is not None else float(v.mean())
        return dict(mean=mean, min=float(v.min()), max=float(v.max()), std=float(v.std()),
                    median=float(np.median(v)),
                    p05=float(np.percentile(v, 5)), p95=float(np.percentile(v, 95)))

    res = {
        "case": case, "time": str(time), "patch": patch,
        "mesh": {"fluid_cells": int(len(nei) and (own.max() + 1)) or int(own.max() + 1),
                 "patch_faces": nf},
        "area_m2": A,
        "T_wall": stats(Tw, area),
        "q_flux": stats(q, area),
        "Q_W": Q,                       # wallHeatFlux 积分（近壁梯度，网格敏感）
        "Q_flux_W": Q_flux,             # ṁ·cp·ΔT 焓平衡（守恒，对比主口径）
        "h_local": stats(h_loc, area),
        "h_global": stats(h_glb, area),
        "h_bar": h_bar, "h_m": hm_bar,
        "h_UA": UA / A if A > 0 else float("nan"),   # 可横比的整体换热系数 = UA/A
        "h_valid_area_m2": float(np.nansum(area[np.isfinite(h_glb)])),
        "T_wall_bar": Tw_bar,
        "UA_W_per_K": UA,
        "NTU_m": NTU_m, "hm_A_m3_s": hm_A, "Q_vol_m3_s": Q_vol,
        "transfer_limited": transfer_limited,
        "m_cond_transfer_kg_s": m_cond_transfer,
        "T_out_wet": T_out_wet,
        "T_bulk": Tbulk, "T_out_mw": Tout_mw, "energy_ratio": ratio,
        "rho": rho, "cp": cp,
        "m_total": m_total, "m_dry_air": m_da,
        "Y_in": y_in, "w_in": w_in,
        "Y_sat_wall": stats(Ys, area),
        "m_cond_kg_s": m_cond,
        "m_cond_mg_h": m_cond * 3.6e6,
        "m_cond_max_kg_s": m_cond_max,
        "eta_cond": eta,
        "Y_out": Y_out, "w_out": w_out, "w_out_g_per_kg": w_out * 1000.0,
        "T_dew_in": t_dew(y_in), "T_dew_out": t_dew(Y_out),
        "condensation_possible": bool(np.mean(Ys) < y_in),
    }
    if verbose:
        print(json.dumps(res, ensure_ascii=False, indent=2, default=float))
    return res


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("case")
    ap.add_argument("time")
    ap.add_argument("--patch", default="interface")
    ap.add_argument("--y-in", type=float, default=Y_H2O_IN)
    ap.add_argument("--json", default=None)
    args = ap.parse_args()
    r = compute(args.case, args.time, args.patch, args.y_in, verbose=args.json is None)
    if args.json:
        with open(args.json, "w", encoding="utf-8") as f:
            json.dump(r, f, ensure_ascii=False, indent=2, default=float)
        print(f"-> {args.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
