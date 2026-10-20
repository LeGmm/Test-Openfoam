#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
 在 OpenFOAM 侧消除「捏点」：把被 2 张面连接的那一对单元并成一个
================================================================================
症状：两个单元之间挂了 **2 张内部面**（Fluent 里两个实体点/线接触，网格器在接触处
      生成面簇）。checkMesh 只告警（`Found 1 neighbouring cells with multiple in
      between faces`）仍打印 Mesh OK；Fluent 自己的网格质量检查也不报（它只查单元
      几何质量，不含拓扑多重连接）。

为什么不能简单删掉一张面：每张面都是两侧单元的闭合面，删掉一张会让两个单元
「开口」（实测 checkMesh → `open cells found, max cell openness 0.024`，Failed 1）。

正确做法：把这一对单元 **并成一个**。两张共享面互为反向（对 A 出、对 B 入），
并集后相互抵消，合并后的单元仍然闭合。

注意：不要用 tools/merge_dup_cells.py 迭代融合——那份脚本在 polyhedral 网格上会
级联误判，实测把 64729 单元并成 28378（网格被毁）。本脚本只处理检测到的单元对，
不迭代、不扩散。

用法：
  python tools/fix_pinched_cells.py <polyMesh 目录> [-n]
================================================================================

  ⚠️ 已弃用（2026-10-20）：铰链 pinch 合并单对会让缺陷扩散（1 对 → 7 对 → 36352 单元）。
  默认已加 --force 硬守护拒绝执行；判定请用 tools/mesh_topology_report.py。
"""
from __future__ import annotations

import argparse
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fluent_cas_to_foam_cht import (  # noqa: E402
    write_boundary, write_faces, write_labels, write_points,
)


def read_labels(path):
    txt = open(path, encoding="utf-8").read()
    body = txt.split("(", 1)[1].rsplit(")", 1)[0]
    return np.array([int(x) for x in body.split()], dtype=np.int64)


def read_points(path):
    pts = []
    for line in open(path, encoding="utf-8"):
        s = line.strip()
        if s.startswith("(") and s.endswith(")") and s.count(" ") >= 2:
            t = s.strip("()").split()
            if len(t) == 3:
                pts.append([float(x) for x in t])
    return np.array(pts, dtype=np.float64)


def read_faces(path):
    faces = []
    for line in open(path, encoding="utf-8"):
        s = line.strip()
        if "(" in s and s.endswith(")") and s[0].isdigit():
            n, rest = s.split("(", 1)
            ids = [int(x) for x in rest[:-1].split()]
            if len(ids) == int(n):
                faces.append(ids)
    return faces


def read_boundary(path):
    txt = open(path, encoding="utf-8").read()
    patches = {}
    for m in re.finditer(r"^\s{4}(\S+)\s*\n\s*\{(.*?)^\s{4}\}", txt, re.S | re.M):
        name, body = m.group(1), m.group(2)
        nf = re.search(r"nFaces\s+(\d+)", body)
        sf = re.search(r"startFace\s+(\d+)", body)
        ty = re.search(r"type\s+(\S+);", body)
        if not (nf and sf and ty):
            continue
        extra = []
        for key in ("sampleMode", "sampleRegion", "samplePatch"):
            mm = re.search(key + r"\s+(\S+);", body)
            if mm:
                extra.append((key, mm.group(1)))
        patches[name] = {"type": ty.group(1), "nFaces": int(nf.group(1)),
                         "startFace": int(sf.group(1)), "extra": extra,
                         "mapped": ty.group(1) == "mappedWall"}
    return patches


class DSU:
    def __init__(self, n):
        self.p = list(range(n))

    def find(self, a):
        while self.p[a] != a:
            self.p[a] = self.p[self.p[a]]
            a = self.p[a]
        return a

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[rb] = ra


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("mesh")
    ap.add_argument("--force", action="store_true",
                    help="确实要在铰链拓扑上强行并单元才加（会把缺陷扩散，通常不该用）")
    ap.add_argument("-n", "--dry-run", action="store_true")
    a = ap.parse_args()
    mesh = a.mesh.rstrip("\\/")

    pts = read_points(os.path.join(mesh, "points"))
    faces = read_faces(os.path.join(mesh, "faces"))
    owner = read_labels(os.path.join(mesh, "owner"))
    neigh = read_labels(os.path.join(mesh, "neighbour"))
    n_int = len(neigh)
    n_cells = int(max(owner.max(), neigh.max())) + 1
    print(f"节点 {len(pts)}  面 {len(faces)}  内部面 {n_int}  单元 {n_cells}")

    pairs = np.stack([owner[:n_int], neigh[:n_int]], axis=1)
    order = np.lexsort((pairs[:, 1], pairs[:, 0]))
    sp = pairs[order]
    same = np.all(sp[1:] == sp[:-1], axis=1)
    groups = {}
    for r in np.where(same)[0]:
        key = (int(sp[r][0]), int(sp[r][1]))
        groups.setdefault(key, set()).update([int(order[r]), int(order[r + 1])])

    if not groups:
        print("没有重复内部面，无需处理。")
        return 0

    dsu = DSU(n_cells)
    drop = set()
    for key, ids in groups.items():
        print(f"  单元对 {key} 由面 {sorted(ids)} 双连接 → 合并为一个单元")
        dsu.union(*key)
        drop |= ids

    if a.dry_run:
        print(f"[dry-run] 将合并 {len(groups)} 对、删除 {len(drop)} 张面")
        return 0

    if not a.force:
        print("\n已停止：这个操作用在铰链 pinch 上会把缺陷扩散给邻居（实测 1 对 → 7 对 → 级联\n"
              "到 36352 单元），网格会被毁。先跑 tools/mesh_topology_report.py 判断类型；\n"
              "若确实是【真重复面】，用 mesh_topology_report.py --fix。\n"
              "确实要强行合并请加 --force。")
        return 2

    # 单元重编号（保留代表元顺序）
    reps = {}
    for c in range(n_cells):
        r = dsu.find(c)
        reps.setdefault(r, len(reps))
    cell_map = np.array([reps[dsu.find(c)] for c in range(n_cells)], dtype=np.int64)
    n_new = len(reps)

    keep = [i for i in range(len(faces)) if i not in drop]
    n_bnd = len(faces) - n_int
    ints, bnds = [], []
    for i in keep:
        o = int(cell_map[owner[i]])
        if i < n_int:
            n = int(cell_map[neigh[i]])
            if o == n:          # 合并后落入同一单元的面，丢弃
                continue
            ints.append((o, n, i))
        else:
            bnds.append((o, int(owner[i]), i))   # (new owner, old owner, face)

    ints.sort(key=lambda t: (t[0], t[1]))
    new_faces, new_owner, new_neigh = [], [], []
    for o, n, i in ints:
        new_faces.append(faces[i])
        new_owner.append(o)
        new_neigh.append(n)
    n_int_new = len(new_faces)

    for o, _, i in bnds:
        new_faces.append(faces[i])
        new_owner.append(o)

    # 补丁 startFace 重算（边界面按原补丁顺序连续排列）
    patches = read_boundary(os.path.join(mesh, "boundary"))
    cursor = n_int_new
    for name, p in patches.items():
        p["startFace"] = cursor
        cursor += p["nFaces"]
    assert cursor == len(new_faces), f"边界面数不匹配 {cursor} vs {len(new_faces)}"

    used = np.zeros(len(pts), dtype=bool)
    for ids in new_faces:
        used[ids] = True
    pmap = np.full(len(pts), -1, dtype=np.int64)
    pmap[used] = np.arange(int(used.sum()), dtype=np.int64)
    new_faces = [[int(pmap[i]) for i in ids] for ids in new_faces]

    write_points(os.path.join(mesh, "points"), pts[used])
    write_faces(os.path.join(mesh, "faces"), new_faces)
    write_labels(os.path.join(mesh, "owner"), new_owner)
    write_labels(os.path.join(mesh, "neighbour"), new_neigh)
    write_boundary(os.path.join(mesh, "boundary"), patches)

    print(f"完成：单元 {n_cells} → {n_new}；内部面 {n_int} → {n_int_new}；"
          f"节点 {len(pts)} → {int(used.sum())}")
    print("注意：单元数变了，时间目录/单元场（如 V）必须重生成，0/ 下的 uniform 场可直接用。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
