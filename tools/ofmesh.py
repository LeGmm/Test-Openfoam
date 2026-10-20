#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
 polyMesh 读写与几何小工具（多个脚本共用的那部分）
================================================================================
 由来：read_points / read_faces / read_labels / read_boundary / face_area 这几段
 逻辑原本在 check_face_orientation、drop_duplicate_faces、fix_pinched_cells、
 fix_zero_volume_cells、map_interface_T、merge_dup_cells、mesh_stats 里各抄了一份。
 此处集中实现，新脚本一律从这里 import；老脚本按需逐步迁移（避免动已验证的
 转换/映射管线）。

 约定（与 fluent_cas_to_foam.py 保持一致）：
   * points 每行 "(x y z)"，按出现顺序编号，0-based
   * faces  每行 "n(i0 i1 ...)"，节点 0-based
   * owner/neighbour 为 labelList：先计数再 "( ... )"
   * boundary 每个补丁块 "    name { type ...; nFaces ..; startFace ..; }"
     补丁顺序即补丁索引顺序，startFace 必须随之递增
================================================================================
"""
from __future__ import annotations

import os
import re

import numpy as np

PT = re.compile(r"^\(([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\)$")
FACE = re.compile(r"^(\d+)\((.*)\)$")
BLOCK = re.compile(r"^\s{4}(\S+)\s*\n\s*\{(.*?)^\s{4}\}", re.S | re.M)


def read_points(path: str) -> np.ndarray:
    """读 constant/<region>/polyMesh/points → (N,3) float64。"""
    pts = []
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            m = PT.match(line.strip())
            if m:
                pts.append([float(x) for x in m.groups()])
    if not pts:
        raise ValueError(f"没读到任何点：{path}")
    return np.asarray(pts, dtype=np.float64)


def read_faces(path: str) -> list[list[int]]:
    """读 faces → 每张面的节点编号列表（0-based）。"""
    faces = []
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.strip()
            if not s or not s[0].isdigit():
                continue
            m = FACE.match(s)
            if not m:
                continue
            ids = [int(x) for x in m.group(2).split()]
            if len(ids) != int(m.group(1)):
                raise ValueError(f"面节点数与声明不符：{s[:60]}")
            faces.append(ids)
    if not faces:
        raise ValueError(f"没读到任何面：{path}")
    return faces


def read_labels(path: str) -> np.ndarray:
    """读 owner / neighbour 这类 labelList。"""
    txt = open(path, encoding="utf-8", errors="replace").read()
    head, sep, body = txt.partition("(")
    if not sep:
        raise ValueError(f"不是 labelList 格式：{path}")
    body = body.rsplit(")", 1)[0]
    return np.asarray([int(x) for x in body.split()], dtype=np.int64)


def read_boundary(path: str) -> dict:
    """读 boundary → {补丁名: {type, nFaces, startFace, extra, mapped}}（保持文件顺序）。"""
    txt = open(path, encoding="utf-8", errors="replace").read()
    patches: dict[str, dict] = {}
    for m in BLOCK.finditer(txt):
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
        patches[name] = {
            "type": ty.group(1),
            "nFaces": int(nf.group(1)),
            "startFace": int(sf.group(1)),
            "extra": extra,
            "mapped": ty.group(1) == "mappedWall",
        }
    if not patches:
        raise ValueError(f"没读到任何补丁：{path}")
    return patches


def read_mesh(mesh: str) -> tuple:
    """一次读齐，返回 (points, faces, owner, neighbour, patches)。"""
    mesh = mesh.rstrip("\\/")
    pts = read_points(os.path.join(mesh, "points"))
    faces = read_faces(os.path.join(mesh, "faces"))
    own = read_labels(os.path.join(mesh, "owner"))
    nei = read_labels(os.path.join(mesh, "neighbour"))
    patches = read_boundary(os.path.join(mesh, "boundary"))
    return pts, faces, own, nei, patches


def face_area_centre(pts: np.ndarray, ids) -> tuple[float, np.ndarray]:
    """Newell 法：向量面积 = Σ cross/2，面积 = |v|（不要再除以 2）。"""
    p = pts[list(ids)]
    a = np.zeros(3)
    for j in range(len(p)):
        a += np.cross(p[j], p[(j + 1) % len(p)])
    a *= 0.5
    return float(np.linalg.norm(a)), p.mean(axis=0)


def duplicate_pairs(owner: np.ndarray, neighbour: np.ndarray) -> dict:
    """找出 (owner, neighbour) 相同的内部面组 → {(o, n): [face idx, ...]}。"""
    n_int = len(neighbour)
    pairs = np.stack([owner[:n_int], neighbour[:n_int]], axis=1)
    order = np.lexsort((pairs[:, 1], pairs[:, 0]))
    sp = pairs[order]
    same = np.all(sp[1:] == sp[:-1], axis=1)
    groups: dict[tuple, list] = {}
    for r in np.where(same)[0]:
        key = (int(sp[r][0]), int(sp[r][1]))
        groups.setdefault(key, set()).update([int(order[r]), int(order[r + 1])])
    return {k: sorted(v) for k, v in groups.items()}


def classify_multi_face(pts, faces, ids, tol: float = 1e-9) -> dict:
    """判定一对单元间的多张面属于【真重复】还是【铰链】。

    * 真重复：所有面心在 tol 内重合且面积基本相等 → 可删多余面；
    * 铰链  ：面不重合（通常共用一条棱）→ 删面会让单元开口，不可程序修。
    """
    geo = {i: face_area_centre(pts, faces[i]) for i in ids}
    areas = [geo[i][0] for i in ids]
    centres = [geo[i][1] for i in ids]
    dmax = max((np.linalg.norm(a - b) for k, a in enumerate(centres)
                for b in centres[k + 1:]), default=0.0)
    shared = set(faces[ids[0]])
    for i in ids[1:]:
        shared &= set(faces[i])
    coincident = dmax < tol and max(areas) / max(min(areas), 1e-30) < 1.05
    return {
        "ids": list(ids),
        "areas": areas,
        "centres": centres,
        "dmax": float(dmax),
        "shared_nodes": len(shared),
        "verdict": "duplicate" if coincident else "hinge",
    }
