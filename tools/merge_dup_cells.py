#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
 polyMesh 手术：消除「同一对单元被 ≥2 个内部面连接」的非流形缺陷
================================================================================
 背景：Fluent 网格在 v8_1 存在退化薄片单元（流体区 20474 与 5 个邻居各有 2 面相连，
 固体区 41281 同样），导致通量双重计入 + ldu 寻址异常 → 压力-速度耦合畸形
 （Courant 爆炸 ~2.5e4）→ 能量场飞掉。caseCHT 与 case（multiphaseEulerFoam）共用。
 手术（并查集融合，迭代直到干净）：
   对每个多面连接单元对 (a,b)：删除这些连接面，a、b 并入同一簇；
   簇内自连接面全部丢弃；最后单元重编号、内部面 owner<neighbour +
   (owner,neighbour) 字典序重排、边界面保持补丁分组（startFace 平移）。
 用法：
   python tools/merge_dup_cells.py <polyMesh目录> [--apply]
 不带 --apply 只报告（干跑）。
================================================================================

  ⚠️ 已弃用（2026-10-20）：迭代并单元在 polyhedral 网格上会级联误判，
  实测把 64729 单元并成 28378（网格被毁）。判定/修复请用 tools/mesh_topology_report.py。
"""
from __future__ import annotations

import os
import re
import sys

from collections import defaultdict


def header(cls, obj):
    return (
        "/*--------------------------------*- C++ -*----------------------------------*\\\n"
        "| =========                 |                                                 |\n"
        "|  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n"
        "|   \\\\    /   O peration    | Version:  8                                     |\n"
        "|    \\\\  /    A nd           | Website:  https://openfoam.org                  |\n"
        "|     \\\\/     M anipulation  |                                                 |\n"
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


def read_boundary(p):
    txt = open(p, encoding="utf-8").read()
    body = txt.split("(", 1)[1].rsplit(")\n", 1)[0]
    patches = []
    for m in re.finditer(r"\n\s*(\w+)\s*\n\s*\{(.*?)\}", body, re.S):
        blk = m.group(2)
        nf = re.search(r"nFaces\s+(\d+)", blk)
        patches.append((m.group(1), {"lines": blk.strip(), "nFaces": int(nf.group(1))}))
    return patches


def write_faces(p, faces):
    with open(p, "w", encoding="utf-8") as f:
        f.write(header("faceList", "faces"))
        f.write(f"\n{len(faces)}\n(\n")
        for ids in faces:
            f.write(f"{len(ids)}(" + " ".join(str(i) for i in ids) + ")\n")
        f.write(")\n\n// ************************************************************************* //\n")


def write_labels(p, vals):
    with open(p, "w", encoding="utf-8") as f:
        f.write(header("labelList", os.path.basename(p)))
        f.write(f"\n{len(vals)}\n(\n")
        for v in vals:
            f.write(f"{v}\n")
        f.write(")\n\n// ************************************************************************* //\n")


def write_boundary(p, patches):
    with open(p, "w", encoding="utf-8") as f:
        f.write(header("polyBoundaryMesh", "boundary"))
        f.write(f"\n{len(patches)}\n(\n")
        for name, d in patches:
            f.write(f"    {name}\n    {{\n{d['lines']}\n    }}\n")
        f.write(")\n\n// ************************************************************************* //\n")


def find_multi(owner, neigh):
    pairs = defaultdict(list)
    for i, (o, n) in enumerate(zip(owner, neigh)):
        pairs[(min(o, n), max(o, n))].append(i)
    return {k: v for k, v in pairs.items() if len(v) > 1}


def main():
    mesh = sys.argv[1]
    apply = "--apply" in sys.argv
    faces = read_faces(os.path.join(mesh, "faces"))
    owner = read_list(os.path.join(mesh, "owner"))
    neigh = read_list(os.path.join(mesh, "neighbour"))
    boundary = read_boundary(os.path.join(mesh, "boundary"))
    n_int0 = len(neigh)
    print(f"总面 {len(faces)}  内部 {n_int0}  边界 {len(faces) - n_int0}")

    multi = find_multi(owner[:n_int0], neigh[:n_int0])
    print(f"多面连接单元对: {len(multi)} 个")
    for k, v in multi.items():
        print(f"  {k}: {len(v)} 面 {v[:8]}")
    if not multi:
        return 0
    if not apply:
        print("(干跑，未修改；加 --apply 执行手术)")
        return 0

    # ---- 并查集融合 ----
    parent = {}

    def find(c):
        parent.setdefault(c, c)
        while parent[c] != c:
            parent[c] = parent[parent[c]]
            c = parent[c]
        return c

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            lo, hi = min(ra, rb), max(ra, rb)
            parent[hi] = lo

    removed = set()
    for it in range(50):
        cur_o = [find(c) for c in owner[:n_int0]]
        cur_n = [find(c) for c in neigh[:n_int0]]
        multi = find_multi(cur_o, cur_n)
        multi = {k: [i for i in v if i not in removed] for k, v in multi.items()}
        multi = {k: v for k, v in multi.items() if len(v) > 1}
        if not multi:
            print(f"第 {it} 轮后无多面连接，结束")
            break
        print(f"第 {it} 轮: 融合 {len(multi)} 对")
        for (a, b), fids in sorted(multi.items()):
            removed.update(fids)
            union(a, b)

    # ---- 重建面/owner/neighbour ----
    owner2, neigh2, faces2 = [], [], []
    for i in range(len(faces)):
        if i in removed:
            continue
        o = find(owner[i])
        n = find(neigh[i]) if i < n_int0 else -1
        if i < n_int0 and o == n:
            continue          # 簇内自连接面，丢弃
        owner2.append(o)
        neigh2.append(n)
        faces2.append(faces[i])

    used = sorted(set(owner2) | {c for c in neigh2 if c >= 0})
    renum = {c: i for i, c in enumerate(used)}
    owner2 = [renum[c] for c in owner2]
    neigh2 = [renum[c] if c >= 0 else -1 for c in neigh2]
    print(f"单元数 → {len(used)}；面数 → {len(faces2)}")

    int_idx = [i for i, c in enumerate(neigh2) if c >= 0]
    for i in int_idx:
        if owner2[i] > neigh2[i]:
            owner2[i], neigh2[i] = neigh2[i], owner2[i]
            faces2[i] = faces2[i][::-1]
    int_idx.sort(key=lambda i: (owner2[i], neigh2[i]))
    bnd_idx = [i for i, c in enumerate(neigh2) if c < 0]   # 原相对顺序（=补丁分组顺序）

    new_faces = [faces2[i] for i in int_idx] + [faces2[i] for i in bnd_idx]
    new_owner = [owner2[i] for i in int_idx] + [owner2[i] for i in bnd_idx]
    new_neigh = [neigh2[i] for i in int_idx]

    sf = len(int_idx)
    new_boundary = []
    for name, d in boundary:
        lines = re.sub(r"startFace\s+\d+", f"startFace       {sf}", d["lines"])
        new_boundary.append((name, {"lines": lines, "nFaces": d["nFaces"]}))
        sf += d["nFaces"]

    write_faces(os.path.join(mesh, "faces"), new_faces)
    write_labels(os.path.join(mesh, "owner"), new_owner)
    write_labels(os.path.join(mesh, "neighbour"), new_neigh)
    write_boundary(os.path.join(mesh, "boundary"), new_boundary)
    print(f"手术完成：面 {len(faces)} → {len(new_faces)}  内部 {n_int0} → {len(int_idx)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
