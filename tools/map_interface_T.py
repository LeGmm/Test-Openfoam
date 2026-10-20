#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
 map_interface_T.py —— 两步共轭方案的 Step1→Step2 温度映射工具
================================================================================
 功能：
   extract  从 chtMultiRegionFoam 写出的 T 场文件里提取某耦合补丁的
            （面心坐标, T）列表 → CSV。这是 Step 1 的界面温度分布 T_int。
   map      把 CSV 里的 T_int 按最近邻映射到目标算例某补丁的每个面上，
            并把目标 0 文件里该补丁的 BC 改写为 fixedValue + nonuniform。
            （--out-only 时只写出数值列表文件，不改目标文件）
   verify   同网格一致性断言：确认 CSV 与目标补丁面数一致、面心逐一一对应
            （容差 --tol，默认 1e-6 m）。正式管线在 map 前必须先 verify 通过。

 用法：
   extract  <region polyMesh 目录> <T 场文件> <补丁名> <out.csv>
   map      <T_int.csv> <目标 polyMesh 目录> <目标补丁名> <目标 0 文件>
            [--out-only <值列表文件>] [--clamp <min> <max>]
   verify   <T_int.csv> <目标 polyMesh 目录> <目标补丁名> [--tol <米>]

 说明：
   * written T 场的 boundaryField 里，耦合补丁的 value 是 nonuniform 列表
     （turbulentTemperatureCoupledBaffleMixed 每步更新），直接可用。
   * 同一几何的两套网格（caseCHT fluid ↔ case/）面心完全重合，映射距离≈0；
     不同几何之间用最近邻（用于管线验证）。
================================================================================
"""
from __future__ import annotations

import os
import re
import sys

import numpy as np


# ---------------------------------------------------------------- 基础读取
def read_list(p):
    txt = open(p, encoding="utf-8").read()
    body = txt.split("\n(", 1)[1].rsplit(")", 1)[0]
    return [int(x) for x in body.split()]


def read_faces(p):
    raw = open(p, "rb").read()
    head = raw[:1500].decode("latin-1")
    if "format      binary" in head or "format          binary" in head:
        # 二进制 faceList：N\n(\n 后为 int32 流（每面：n + n 个索引）
        i = raw.index(b"\n(") + 2          # 列表括号后紧跟二进制数据
        n_faces = int(raw[:i].rsplit(b"\n", 2)[-2])
        stream = np.frombuffer(raw[i:i + ((len(raw) - i) // 4) * 4], dtype="<i4")
        faces, pos = [], 0
        for _ in range(n_faces):
            n = int(stream[pos]); pos += 1
            faces.append([int(x) for x in stream[pos:pos + n]])
            pos += n
        return faces
    txt = raw.decode("utf-8", "replace")
    body = txt.split("\n(", 1)[1].rsplit(")", 1)[0]
    out = []
    for line in body.strip().splitlines():
        line = line.strip()
        if not line:
            continue
        n, rest = line.split("(", 1)
        out.append([int(x) for x in rest.rstrip(")").split()])
    return out


def read_points(p):
    txt = open(p, encoding="utf-8").read()
    body = txt.split("\n(", 1)[1].rsplit(")", 1)[0]
    vals = [float(x) for x in body.replace("(", " ").replace(")", " ").split()]
    return np.array(vals).reshape(-1, 3)


def read_boundary(p):
    txt = open(p, encoding="utf-8").read()
    body = txt.split("\n(", 1)[1].rsplit("\n)", 1)[0]
    patches = {}
    for m in re.finditer(r"\n\s*(\w+)\s*\n\s*\{([^}]*)\}", body, re.S):
        blk = m.group(2)
        patches[m.group(1)] = (int(re.search(r"startFace\s+(\d+)", blk).group(1)),
                               int(re.search(r"nFaces\s+(\d+)", blk).group(1)))
    return patches


def face_centers(pmesh, patch):
    pts = read_points(os.path.join(pmesh, "points"))
    faces = read_faces(os.path.join(pmesh, "faces"))
    sf, nf = read_boundary(os.path.join(pmesh, "boundary"))[patch]
    cen = np.zeros((nf, 3))
    for k, f in enumerate(range(sf, sf + nf)):
        cen[k] = pts[faces[f]].mean(axis=0)
    return cen


# ------------------------------------------------- 场文件补丁值提取
def extract_patch_values(field_file, patch):
    """从写出的场文件 boundaryField 里取补丁 value 的数值列表。"""
    txt = open(field_file, encoding="utf-8").read()
    i = txt.index("boundaryField")
    seg = txt[i:]
    # 顶层逐块扫描（花括号配平）
    j = seg.index("{")
    pos = j + 1
    depth = 1
    entries = []   # (name_start, name, block_start, block_end)
    while depth > 0 and pos < len(seg):
        ch = seg[pos]
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
        elif ch in " \t\n" and depth == 1:
            m = re.match(r"\s*([\w.]+)\s*\n?\s*\{", seg[pos:])
            if m:
                name = m.group(1)
                b0 = pos + m.group(0).index("{") + pos - pos + pos  # placeholder
                b0 = pos + seg[pos:].index("{")
                # 找配平的块尾
                d = 1
                k = b0 + 1
                while d > 0:
                    if seg[k] == "{":
                        d += 1
                    elif seg[k] == "}":
                        d -= 1
                    k += 1
                entries.append((name, b0 + 1, k - 1))
                pos = k
                continue
        pos += 1

    for name, b0, b1 in entries:
        if name != patch:
            continue
        blk = seg[b0:b1]
        m = re.search(r"value\s+nonuniform\s+List<scalar>\s*(\d+)\s*\((.*?)\)", blk, re.S)
        if not m:
            raise SystemExit(f"补丁 {patch} 的 value 不是 nonuniform 列表")
        n = int(m.group(1))
        vals = [float(x) for x in m.group(2).split()]
        assert len(vals) == n, f"数量不符 {len(vals)} vs {n}"
        return np.array(vals)
    raise SystemExit(f"boundaryField 里找不到补丁 {patch}")


# ------------------------------------------------- extract 子命令
def do_extract(pmesh, field_file, patch, out_csv):
    cen = face_centers(pmesh, patch)
    T = extract_patch_values(field_file, patch)
    assert len(cen) == len(T)
    with open(out_csv, "w", encoding="utf-8") as f:
        f.write("x,y,z,T\n")
        for c, t in zip(cen, T):
            f.write(f"{c[0]:.8g},{c[1]:.8g},{c[2]:.8g},{t:.8g}\n")
    print(f"提取 {patch}: {len(T)} 面 → {out_csv}")
    print(f"  T 范围: {T.min():.3f} ~ {T.max():.3f} K, 面心范围 "
          f"x[{cen[:,0].min()*1000:.1f},{cen[:,0].max()*1000:.1f}] "
          f"y[{cen[:,1].min()*1000:.1f},{cen[:,1].max()*1000:.1f}] "
          f"z[{cen[:,2].min()*1000:.1f},{cen[:,2].max()*1000:.1f}] mm")




def do_extract_vtk(vtk_file, out_csv):
    """从 foamToVTK 的 ascii POLYDATA VTK 里提取（面心, T）。"""
    txt = open(vtk_file, encoding="utf-8").read()
    npts = int(re.search(r"POINTS\s+(\d+)", txt).group(1))
    seg = txt.split("POINTS", 1)[1].split("POLYGONS", 1)[0]
    vals = [float(x) for x in seg.split()[2:]]      # 去掉 "15488 float"
    P = np.array(vals).reshape(-1, 3)[:npts]
    nf = int(re.search(r"POLYGONS\s+(\d+)", txt).group(1))
    seg = txt.split("POLYGONS", 1)[1].split("CELL_DATA", 1)[0]
    poly = [int(x) for x in seg.split()[2:]]        # 去掉 "15360 <size>"
    assert len(poly) >= nf
    # 解析面（每面：n + n 个点索引）
    faces, pos = [], 0
    for _ in range(nf):
        n = poly[pos]; pos += 1
        faces.append(poly[pos:pos + n]); pos += n
    cen = np.array([P[ids].mean(axis=0) for ids in faces])
    # CELL_DATA 段：FIELD attributes 里逐字段找 T（行头 "T 1 N float"）
    cd = txt.split("CELL_DATA", 1)[1]
    lines = cd.splitlines()
    T = None
    k = 0
    while k < len(lines):
        toks = lines[k].split()
        if len(toks) == 4 and toks[0] == "T" and toks[1] == "1":
            n = int(toks[2])
            vals = []
            k += 1
            while len(vals) < n and k < len(lines):
                vals += [float(x) for x in lines[k].split()]
                k += 1
            T = np.array(vals)
            break
        k += 1
    assert T is not None and len(T) == nf, "VTK 里没找到面温度字段 T"
    print(f"（VTK 字段名 T，{nf} 个面值）")
    with open(out_csv, "w", encoding="utf-8") as f:
        f.write("x,y,z,T\n")
        for c, t in zip(cen, T):
            f.write(f"{c[0]:.8g},{c[1]:.8g},{c[2]:.8g},{t:.8g}\n")
    print(f"提取(VTK) {len(T)} 面 → {out_csv}; T 范围 {T.min():.3f} ~ {T.max():.3f} K")

# ------------------------------------------------- map 子命令
def do_map(csv, pmesh, patch, target_file, out_only=None, clamp=None):
    data = np.genfromtxt(csv, delimiter=",", names=True)
    src_xyz = np.stack([data["x"], data["y"], data["z"]], axis=1)
    src_T = data["T"]
    cen = face_centers(pmesh, patch)
    # 最近邻（分块防内存爆）
    Tmap = np.zeros(len(cen))
    dmax = 0.0
    for i in range(0, len(cen), 512):
        d = np.linalg.norm(cen[i:i+512, None, :] - src_xyz[None, :, :], axis=2)
        j = d.argmin(axis=1)
        Tmap[i:i+512] = src_T[j]
        dmax = max(dmax, d.min(axis=1).max())
    if clamp:
        Tmap = np.clip(Tmap, clamp[0], clamp[1])
    print(f"映射 {patch}: {len(cen)} 面, 最近邻最大距离 = {dmax*1000:.3f} mm, "
          f"T 范围 {Tmap.min():.3f} ~ {Tmap.max():.3f} K")

    values = "\n".join(f"{t:.8g}" for t in Tmap)
    new_entry = (f"        type            fixedValue;\n"
                 f"        value           nonuniform List<scalar>\n"
                 f"{len(Tmap)}\n"
                 f"(\n{values}\n)\n        ;")
    if out_only:
        with open(out_only, "w", encoding="utf-8") as f:
            f.write(new_entry)
        print(f"值列表 → {out_only}")
        return 0

    # 改写目标 0 文件的补丁块
    txt = open(target_file, encoding="utf-8").read()
    i = txt.index("boundaryField")
    seg = txt[i:]
    j = seg.index("{")
    pos = j + 1
    depth = 1
    while depth > 0 and pos < len(seg):
        ch = seg[pos]
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
        elif depth == 1:
            m = re.match(r"\s*([\w.]+)\s*\n?\s*\{", seg[pos:])
            if m:
                name = m.group(1)
                b0 = pos + seg[pos:].index("{")
                d = 1
                k = b0 + 1
                while d > 0:
                    if seg[k] == "{":
                        d += 1
                    elif seg[k] == "}":
                        d -= 1
                    k += 1
                if name == patch:
                    new_txt = txt[:i + b0 + 1] + "\n" + new_entry + "\n            " + txt[i + k - 1:]
                    open(target_file, "w", encoding="utf-8").write(new_txt)
                    print(f"已改写 {target_file} 的 {patch} 补丁 → fixedValue nonuniform")
                    return 0
                pos = k
                continue
        pos += 1
    raise SystemExit(f"目标文件里找不到补丁 {patch}")


# ------------------------------------------------- verify 子命令
def do_verify(csv, pmesh, patch, tol=1e-6):
    """同网格一致性断言：映射前确认源 CSV 与目标补丁是同一几何。
    判据：每个目标面心都能在源里找到 ≤ tol 的精确对应（默认 1e-6 m），
    且对应关系是一一映射（防面序错位/重复点）。失败返回 1，供管线 gate。"""
    data = np.genfromtxt(csv, delimiter=",", names=True)
    src_xyz = np.stack([data["x"], data["y"], data["z"]], axis=1)
    cen = face_centers(pmesh, patch)
    print(f"verify: 源 {len(src_xyz)} 面 vs 目标 {patch} {len(cen)} 面, 容差 {tol:g} m")

    if len(src_xyz) != len(cen):
        print(f"FAIL: 面数不一致（{len(src_xyz)} ≠ {len(cen)}）—— 不是同一网格界面")
        return 1

    dmin_all = np.zeros(len(cen))
    jmin_all = np.zeros(len(cen), dtype=int)
    for i in range(0, len(cen), 512):
        d = np.linalg.norm(cen[i:i+512, None, :] - src_xyz[None, :, :], axis=2)
        dmin_all[i:i+512] = d.min(axis=1)
        jmin_all[i:i+512] = d.argmin(axis=1)
    maxd = float(dmin_all.max())
    print(f"  最近邻最大距离 = {maxd*1000:.6f} mm")

    if maxd > tol:
        n_bad = int((dmin_all > tol).sum())
        worst = int(dmin_all.argmax())
        print(f"FAIL: {n_bad}/{len(cen)} 个面心无精确对应（最差面 #{worst}，"
              f"偏差 {dmin_all[worst]*1000:.6f} mm）")
        print("      → 源/目标不是同一几何，请改用 map（最近邻）并人工确认合理性")
        return 1

    n_uniq = len(set(jmin_all.tolist()))
    if n_uniq != len(cen):
        print(f"FAIL: 对应关系不是一一映射（{len(cen)} 个目标面只用了 {n_uniq} 个源点）")
        return 1
    print("PASS: 面数一致、逐面精确对应、一一映射 —— 同网格映射安全")
    return 0


def main():
    argv = sys.argv[1:]
    mode = argv[0]
    if mode == "extract":
        do_extract(argv[1], argv[2], argv[3], argv[4])
    elif mode == "extract-vtk":
        do_extract_vtk(argv[1], argv[2])
    elif mode == "verify":
        tol = 1e-6
        if "--tol" in argv:
            tol = float(argv[argv.index("--tol") + 1])
        sys.exit(do_verify(argv[1], argv[2], argv[3], tol))
    elif mode == "map":
        csv, pmesh, patch, target = argv[1:5]
        out_only = None
        clamp = None
        rest = argv[5:]
        for i, a in enumerate(rest):
            if a == "--out-only" and i + 1 < len(rest):
                out_only = rest[i + 1]
            if a == "--clamp" and i + 2 < len(rest):
                clamp = (float(rest[i + 1]), float(rest[i + 2]))
        do_map(csv, pmesh, patch, target, out_only, clamp)
    else:
        print(__doc__)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
