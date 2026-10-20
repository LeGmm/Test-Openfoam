#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""统计各结构算例的流体单元数 + coldWall（冷凝面）面数与面积。"""
from __future__ import annotations

import os
import re
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CASES = os.path.join(ROOT, "cases")


def read_points(path):
    pts = []
    with open(path, encoding="utf-8") as f:
        txt = f.read()
    blk = txt[txt.index("(\n") + 2: txt.rindex(")")]
    for line in blk.splitlines():
        line = line.strip()
        if not line:
            continue
        t = line.strip("()").split()
        if len(t) == 3:
            pts.append(t)
    return np.array(pts, dtype=np.float64)


def read_faces(path):
    """返回 faces(list[list[int]]) 与 nInternal（从 owner/neighbour 长度推）。"""
    faces = []
    with open(path, encoding="utf-8") as f:
        txt = f.read()
    blk = txt[txt.index("(\n") + 2: txt.rindex(")")]
    for line in blk.splitlines():
        line = line.strip()
        if not line:
            continue
        head, ids = line.split("(")
        faces.append([int(x) for x in ids.rstrip(")").split()])
    return faces


def read_owner_ncells(path):
    """owner 文件里最大值+1 = 单元数。
    注意：必须只取括号内的值——文件头的『条数行』是面数，比任何单元索引都大，
    若把它算进 max 会得到面数而非单元数。"""
    with open(path, encoding="utf-8") as f:
        txt = f.read()
    blk = txt[txt.index("(\n") + 2: txt.rindex(")")]
    mx = -1
    for line in blk.splitlines():
        line = line.strip()
        if line.lstrip("-").isdigit():
            v = int(line)
            if v > mx:
                mx = v
    return mx + 1 if mx >= 0 else 0


def read_patch(path, name):
    with open(path, encoding="utf-8") as f:
        txt = f.read()
    m = re.search(r"\b" + re.escape(name) + r"\s*\{(.*?)\}", txt, re.S)
    if not m:
        return None
    body = m.group(1)
    s = int(re.search(r"startFace\s+(\d+)", body).group(1))
    n = int(re.search(r"nFaces\s+(\d+)", body).group(1))
    return s, n


def area_of(poly, pts):
    a = np.zeros(3)
    nn = len(poly)
    for j in range(nn):
        a += np.cross(pts[poly[j]], pts[poly[(j + 1) % nn]])
    return 0.5 * np.linalg.norm(a)


def main():
    names = sys.argv[1:] if len(sys.argv) > 1 else sorted(os.listdir(CASES))
    rows = []
    for n in names:
        d = os.path.join(CASES, n)
        pm = os.path.join(d, "constant", "polyMesh")
        if not os.path.isfile(os.path.join(pm, "points")):
            continue
        pts = read_points(os.path.join(pm, "points"))
        faces = read_faces(os.path.join(pm, "faces"))
        ncell = read_owner_ncells(os.path.join(pm, "owner"))
        sp = read_patch(os.path.join(pm, "boundary"), "coldWall")
        if sp is None:
            print(f"{n}: 无 coldWall")
            continue
        s, cnt = sp
        area = sum(area_of(faces[s + i], pts) for i in range(cnt))
        rows.append((n, ncell, cnt, area))
        print(f"{n}: fluid={ncell} coldWall={cnt}面  A={area*1e4:.2f} cm²")

    if rows:
        print("\n=== 冷凝面面积排序（面积越大，潜在冷凝量越大）===")
        for n, nc, cnt, a in sorted(rows, key=lambda r: -r[3]):
            print(f"  {n:8s} 流体单元={nc:7d}  coldWall={cnt:6d} 面   A={a*1e4:8.2f} cm²")


if __name__ == "__main__":
    main()
