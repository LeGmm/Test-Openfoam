#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
 在 OpenFOAM 侧清掉 Fluent 网格带来的「捏点」（非流形双连接面）
================================================================================
症状：两个单元之间挂了**两张内部面**（Fluent 里两个实体点/线接触，网格器在
      接触处生成面簇）。checkMesh 只告警：`Found 1 neighbouring cells with
      multiple in between faces`，仍打印 Mesh OK；Fluent 自己的网格质量检查
      也不报——它的口径是单元几何质量，不含这种拓扑多重连接。

处理：同一对 (owner, neighbour) 出现在多张内部面时，保留面积最大的一张，
      其余删掉；同步重排 owner/neighbour、修 boundary 的 startFace、清掉
      因删面而不再被引用的节点。删完两个单元仍共用一张面，网格保持闭合。

用法：
  python tools/drop_duplicate_faces.py <case>/constant/<region>/polyMesh
  python tools/drop_duplicate_faces.py caseCHT_qishi/constant/fluid/polyMesh -n   # 只看不改
================================================================================

  ⚠️ 已弃用（2026-10-20）：删一张面会让铰链单元开口（checkMesh: open cells 2, Failed 1）。
  只有"真重复面"才能删，判定请用 tools/mesh_topology_report.py --fix。
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fluent_cas_to_foam_cht import write_boundary, write_faces, write_labels, write_points  # noqa: E402


def read_labels(path: str) -> np.ndarray:
    txt = open(path, encoding="utf-8").read()
    body = txt.split("(", 1)[1].rsplit(")", 1)[0]
    return np.array([int(x) for x in body.split()], dtype=np.int64)


def read_points(path: str) -> np.ndarray:
    pts = []
    for line in open(path, encoding="utf-8"):
        s = line.strip()
        if s.startswith("(") and s.endswith(")") and s.count(" ") >= 2:
            t = s.strip("()").split()
            if len(t) == 3:
                pts.append([float(x) for x in t])
    return np.array(pts, dtype=np.float64)


def read_faces(path: str) -> list:
    faces = []
    for line in open(path, encoding="utf-8"):
        s = line.strip()
        if "(" in s and s.endswith(")") and s[0].isdigit():
            n, rest = s.split("(", 1)
            ids = [int(x) for x in rest[:-1].split()]
            if len(ids) == int(n):
                faces.append(ids)
    return faces


def read_boundary(path: str):
    import re
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


def face_area(pts, ids):
    p = pts[ids]
    a = np.zeros(3)
    for j in range(len(p)):
        a += np.cross(p[j], p[(j + 1) % len(p)])
    return float(np.linalg.norm(a * 0.5))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("mesh")
    ap.add_argument("-n", "--dry-run", action="store_true")
    a = ap.parse_args()
    mesh = a.mesh.rstrip("\\/")

    pts = read_points(os.path.join(mesh, "points"))
    faces = read_faces(os.path.join(mesh, "faces"))
    owner = read_labels(os.path.join(mesh, "owner"))
    neigh = read_labels(os.path.join(mesh, "neighbour"))
    n_int = len(neigh)
    print(f"节点 {len(pts)}  面 {len(faces)}  内部面 {n_int}  边界面 {len(faces) - n_int}")

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

    drop = set()
    for key, ids in groups.items():
        ids = sorted(ids)
        areas = {i: face_area(pts, faces[i]) for i in ids}
        keep = max(areas, key=areas.get)
        print(f"  (owner {key[0]}, neighbour {key[1]}): 面 {ids} → 保留 {keep}"
              f"（面积 {areas[keep]:.4e} m²），删除 "
              f"{ {i: f'{areas[i]:.4e}' for i in ids if i != keep} }")
        drop.update(i for i in ids if i != keep)

    if a.dry_run:
        print(f"[dry-run] 将删除 {len(drop)} 张面：{sorted(drop)}")
        return 0

    keep_idx = [i for i in range(len(faces)) if i not in drop]
    remap = {old: new for new, old in enumerate(keep_idx)}
    n_drop_int = sum(1 for i in drop if i < n_int)
    assert len(drop) == n_drop_int, "本次只支持删内部面"

    new_faces = [faces[i] for i in keep_idx]
    new_owner = [int(owner[i]) for i in keep_idx]
    new_neigh = [int(neigh[i]) for i in keep_idx if i < n_int]

    # 未使用节点清理
    used = np.zeros(len(pts), dtype=bool)
    for ids in new_faces:
        used[ids] = True
    pmap = np.full(len(pts), -1, dtype=np.int64)
    pmap[used] = np.arange(int(used.sum()), dtype=np.int64)
    new_faces = [[int(pmap[i]) for i in ids] for ids in new_faces]
    new_pts = pts[used]
    print(f"删面 {len(drop)}，内部面 {n_int} → {len(new_neigh)}；"
          f"节点 {len(pts)} → {int(used.sum())}")

    patches = read_boundary(os.path.join(mesh, "boundary"))
    for name, p in patches.items():
        p["startFace"] -= n_drop_int

    write_points(os.path.join(mesh, "points"), new_pts)
    write_faces(os.path.join(mesh, "faces"), new_faces)
    write_labels(os.path.join(mesh, "owner"), new_owner)
    write_labels(os.path.join(mesh, "neighbour"), new_neigh)
    write_boundary(os.path.join(mesh, "boundary"), patches)
    print("已重写 polyMesh（points/faces/owner/neighbour/boundary）")
    return 0


if __name__ == "__main__":
    sys.exit(main())
