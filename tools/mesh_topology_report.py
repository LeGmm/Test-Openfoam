#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
 诊断「一对单元之间有多张面」到底是哪种情况，并只对可修的那种动手
================================================================================
 背景：checkMesh 的 `Found N neighbouring cells with multiple in between faces`
       只报告现象，不区分原因。实际上有两种完全不同的情况：

   ① 真重复（coincident duplicate）：两张面**几何重合**（面心几乎同点、面积几乎相等）。
      这是 Fluent 导出/拼接留下的废面 —— 删掉多余的一张，两个单元仍然闭合。
      → 程序可修，本脚本 `--fix` 直接修。

   ② 铰链（pinch / edge-connected）：两张面**不重合**，各自面积、位置都不同，
      而且共用一条棱（至少 2 个公共节点）。它们是同一对单元绕一条棱"折"起来的
      两个面瓣，本来就是两个单元各自的闭合边界之一。
      → 删任何一张都会让单元开口（checkMesh 报 open cells）；
        把两个单元并成一个只会把同一现象扩散给邻居（实测可级联到半个网格）。
        程序**不能**正确修，正解是回 CAD/网格器让两体变成面接触或留间隙；
        若只是求解，可以原样保留（单元闭合，有限体积逐面独立积分，不影响守恒）。

 用法：
   python tools/mesh_topology_report.py <polyMesh>          # 只诊断
   python tools/mesh_topology_report.py <polyMesh> --fix    # 仅删真重复面
================================================================================
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ofmesh import classify_multi_face, duplicate_pairs, read_mesh  # noqa: E402
from fluent_cas_to_foam_cht import (  # noqa: E402
    write_boundary, write_faces, write_labels, write_points,
)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("mesh")
    ap.add_argument("--fix", action="store_true")
    ap.add_argument("--tol", type=float, default=1e-9, help="面心重合容差（米）")
    a = ap.parse_args()
    mesh = a.mesh.rstrip("\\/")

    pts, faces, owner, neigh, patches = read_mesh(mesh)
    n_int = len(neigh)
    print(f"网格 {mesh}\n  节点 {len(pts)}  面 {len(faces)}  内部面 {n_int}")

    groups = duplicate_pairs(owner, neigh)

    if not groups:
        print("\n没有多面连接的单元对 —— 无需任何处理。")
        return 0

    dup_fixable, hinge = [], []
    print(f"\n发现 {len(groups)} 对单元之间存在多张面：\n")
    for key, ids in groups.items():
        info = classify_multi_face(pts, faces, ids, a.tol)
        areas, centres, dmax = info["areas"], info["centres"], info["dmax"]

        print(f"  单元对 {key}：{len(ids)} 张面")
        for k, i in enumerate(ids):
            c = centres[k]
            print(f"    面 {i}: 面积 {areas[k]*1e6:.4f} mm²，"
                  f"中心 ({c[0]*1e3:.3f},{c[1]*1e3:.3f},{c[2]*1e3:.3f}) mm")
        print(f"    面心最大间距 {dmax*1e6:.4f} µm，共用节点 {info['shared_nodes']} 个")
        if info["verdict"] == "duplicate":
            keep = int(np.argmax(areas))
            drop = [i for i in ids if i != ids[keep]]
            print(f"    → 判定【真重复】：可删 {drop}，保留面积最大的一张\n")
            dup_fixable.append((key, ids, drop))
        else:
            print("    → 判定【铰链 pinch】：合法拓扑，删面会让单元开口，合并会让缺陷扩散")
            print("      建议：原样保留（不影响守恒）；要彻底消掉需回 CAD/网格器留间隙\n")
            hinge.append((key, ids))

    if not a.fix:
        print(f"汇总：真重复 {len(dup_fixable)} 对（可用 --fix 修）；铰链 {len(hinge)} 对（不建议程序修）")
        return 0

    if not dup_fixable:
        print("没有可自动修复的真重复面，未改动网格。")
        return 0

    drop = set()
    for _, _, d in dup_fixable:
        drop.update(d)
    keep_idx = [i for i in range(len(faces)) if i not in drop]
    n_drop_int = sum(1 for i in drop if i < n_int)

    new_faces = [faces[i] for i in keep_idx]
    new_owner = [int(owner[i]) for i in keep_idx]
    new_neigh = [int(neigh[i]) for i in keep_idx if i < n_int]

    used = np.zeros(len(pts), dtype=bool)
    for ids in new_faces:
        used[ids] = True
    pmap = np.full(len(pts), -1, dtype=np.int64)
    pmap[used] = np.arange(int(used.sum()), dtype=np.int64)
    new_faces = [[int(pmap[i]) for i in ids] for ids in new_faces]

    for p in patches.values():
        p["startFace"] -= n_drop_int

    write_points(os.path.join(mesh, "points"), pts[used])
    write_faces(os.path.join(mesh, "faces"), new_faces)
    write_labels(os.path.join(mesh, "owner"), new_owner)
    write_labels(os.path.join(mesh, "neighbour"), new_neigh)
    write_boundary(os.path.join(mesh, "boundary"), patches)
    print(f"\n已删除 {len(drop)} 张真重复面并重写 polyMesh（节点 {len(pts)} → {int(used.sum())}）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
