#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
 polyMesh 修复：给零体积单元注入真实体积（节点偏移法）
================================================================================
 背景：Fluent 网格存在「零体积枢纽薄片」——一个单元的全部面几乎共面（体积 ~1e-26 m³），
 却通过十几个微小面连接周围一圈单元。后果：动量矩阵对角元 Ap=V/dt≈0 →
 压力方程系数 1/Ap 发散 → 全场压力解损坏（速度爆炸、T 飞掉），
 且与入口 BC/重力/线性求解器类型无关（对照实验已排除）。
 修复：对每个近零体积单元，取其最大面，把该面的一个节点沿面法向偏移 delta，
 使单元获得真实体积。delta 默认 2e-5 m（0.02mm，远小于网格尺度 ~0.4mm），
 对共享该节点的邻居单元的扰动可忽略。
 用法：
   python tools/fix_zero_volume_cells.py <polyMesh目录> [--apply] [--delta 2e-5] [--vmin 1e-15]
 不带 --apply 只报告。
================================================================================
"""
from __future__ import annotations

import os
import sys

import numpy as np


def header(cls, obj):
    return (
        "/*--------------------------------*- C++ -*----------------------------------*\\\n"
        "| =========                 |                                                 |\n"
        "|  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n"
        "|   \\\\    /   O peration    | Version:  8                                     |\n"
        "|    \\\\  /    A nd           | Website:  https://openfoam.org                  |\n"
        "|    \\\\/     M anipulation  |                                                 |\n"
        "\\*---------------------------------------------------------------------------*/\n"
        "FoamFile\n"
        "{\n"
        "    version     2.0;\n"
        "    format      ascii;\n"
        f"    class       {cls};\n"
        "    location    \"constant/polyMesh\";\n"
        f"    object      {obj};\n"
        "}\n"
        "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n"
    )


def read_list(p):
    txt = open(p, encoding="utf-8").read()
    body = txt.split("(", 1)[1].rsplit(")", 1)[0]
    return [int(x) for x in body.split()]


def read_faces(p):
    txt = open(p, encoding="utf-8").read()
    body = txt.split("(", 1)[1].rsplit(")", 1)[0]
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
    body = txt.split("(", 1)[1].rsplit(")", 1)[0]
    vals = [float(x) for x in body.replace("(", " ").replace(")", " ").split()]
    return np.array(vals).reshape(-1, 3)


def write_points(p, pts):
    with open(p, "w", encoding="utf-8") as f:
        f.write(header("vectorField", "points"))
        f.write(f"\n{len(pts)}\n(\n")
        for x in pts:
            f.write(f"({x[0]:.10g} {x[1]:.10g} {x[2]:.10g})\n")
        f.write(")\n\n// ************************************************************************* //\n")


def face_area_vec(pts, ids):
    v = np.zeros(3)
    for j in range(len(ids)):
        v += np.cross(pts[ids[j]], pts[ids[(j + 1) % len(ids)]]) / 2
    return v


def cell_volume(pts, faces, c, cell_faces):
    """OpenFOAM 同款：面金字塔对单元形心求和。"""
    fs = cell_faces.get(c, [])
    nodes = set()
    for f in fs:
        nodes.update(faces[f])
    cen = pts[list(nodes)].mean(axis=0)
    vol = 0.0
    for f in fs:
        ids = faces[f]
        a = face_area_vec(pts, ids)
        vol += abs(a @ (pts[ids].mean(axis=0) - cen)) / 3
    return vol, fs, cen


def main():
    mesh = sys.argv[1]
    apply = "--apply" in sys.argv
    delta = 2e-5
    vmin = 1e-15
    for i, a in enumerate(sys.argv):
        if a == "--delta" and i + 1 < len(sys.argv):
            delta = float(sys.argv[i + 1])
        if a == "--vmin" and i + 1 < len(sys.argv):
            vmin = float(sys.argv[i + 1])

    pts = read_points(os.path.join(mesh, "points"))
    faces = read_faces(os.path.join(mesh, "faces"))
    owner = read_list(os.path.join(mesh, "owner"))
    n_cells = max(owner) + 1

    neigh = read_list(os.path.join(mesh, "neighbour"))
    n_cells = max(owner + neigh) + 1
    n_int = len(neigh)

    cell_faces = {}
    for f, o in enumerate(owner):
        cell_faces.setdefault(o, []).append(f)
    for f in range(n_int):
        cell_faces.setdefault(neigh[f], []).append(f)

    # 扫描零体积单元
    bad = []
    for c in range(n_cells):
        vol, fs, cen = cell_volume(pts, faces, c, cell_faces)
        if vol < vmin:
            bad.append((c, vol, fs, cen))
    print(f"扫描 {n_cells} 单元：体积 < {vmin:.1e} m³ 的单元 {len(bad)} 个")
    for c, vol, fs, cen in bad:
        print(f"  cell {c}: V={vol:.3e} m³, 面数={len(fs)}, "
              f"质心=({cen[0]*1000:+.2f},{cen[1]*1000:+.2f},{cen[2]*1000:+.2f}) mm")
    if not bad:
        return 0
    if not apply:
        print("(干跑，未修改；加 --apply 执行修复)")
        return 0

    # 逐个修复：最大面的一个节点沿法向偏移
    for c, vol, fs, cen in bad:
        # 最大面
        areas = [(np.linalg.norm(face_area_vec(pts, faces[f])), f) for f in fs]
        areas.sort(reverse=True)
        a0, f0 = areas[0]
        n = face_area_vec(pts, faces[f0])
        n = n / (np.linalg.norm(n) + 1e-300)
        # 选该面上「只被少量单元共享」的节点优先（减小扰动面）——这里直接取第一个节点
        node = faces[f0][0]
        best = (vol, None)
        for sgn in (+1, -1):
            trial = pts.copy()
            trial[node] = pts[node] + sgn * delta * n
            v2, _, _ = cell_volume(trial, faces, c, cell_faces)
            if v2 > best[0]:
                best = (v2, trial)
        if best[1] is None:
            print(f"  cell {c}: 偏移后仍无体积，跳过（需人工处理）")
            continue
        pts = best[1]
        print(f"  cell {c}: V {vol:.2e} → {best[0]:.2e} m³（节点 {node} 偏移 {delta:.1e} m）")

    write_points(os.path.join(mesh, "points"), pts)
    print("points 已写回")
    return 0


if __name__ == "__main__":
    sys.exit(main())
