#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Check a converted OpenFOAM polyMesh for internal face orientation / ordering conformance.

Checks:
  A. owner < neighbour for every internal face (OpenFOAM lower-face convention)
  B. owner list non-decreasing over internal faces (required face ordering)
  C. boundary faces ordered after internal faces, startFace consistent
  D. geometric: face area vector points owner -> neighbour (internal)
     and outward from owner cell (boundary)
"""
import sys
import numpy as np


def read_foam_list(path, dtype=float):
    with open(path) as f:
        txt = f.read()
    i = txt.find("(")
    j = txt.rfind(")")
    body = txt[i + 1 : j].strip()
    if dtype is float:
        out = [tuple(map(float, ln.strip("()").split())) for ln in body.splitlines() if ln.strip()]
        return np.array(out)
    return np.array([int(ln.strip()) for ln in body.splitlines() if ln.strip()], dtype=np.int64)


def read_faces(path):
    with open(path) as f:
        txt = f.read()
    i = txt.find("(")
    j = txt.rfind(")")
    faces = []
    for ln in txt[i + 1 : j].splitlines():
        ln = ln.strip()
        if not ln:
            continue
        k = ln.index("(")
        n = int(ln[:k])
        ids = np.array(list(map(int, ln[k + 1 : -1].split())), dtype=np.int64)
        assert len(ids) == n, f"face count mismatch: {ln[:40]}"
        faces.append(ids)
    return faces


def read_boundary(path):
    patches = {}
    with open(path) as f:
        txt = f.read()
    import re
    for m in re.finditer(r"(\w+)\s*\{([^}]*)\}", txt):
        name, body = m.group(1), m.group(2)
        m2 = re.search(r"nFaces\s+(\d+)", body)
        if not m2:
            continue  # skip FoamFile header block
        nf = int(m2.group(1))
        sf = int(re.search(r"startFace\s+(\d+)", body).group(1))
        patches[name] = (sf, nf)
    return patches


def main(pm_dir):
    pts = read_foam_list(pm_dir + "/points")
    faces = read_faces(pm_dir + "/faces")
    owner = read_foam_list(pm_dir + "/owner", int)
    neigh = read_foam_list(pm_dir + "/neighbour", int)
    patches = read_boundary(pm_dir + "/boundary")

    nf = len(faces)
    n_int = len(neigh)
    print(f"points={len(pts)} faces={nf} internal={n_int} boundary={nf - n_int}")
    assert len(owner) == nf

    ok = True

    # A. owner < neighbour
    bad_a = np.where(owner[:n_int] >= neigh)[0]
    print(f"A. owner<neighbour violated: {len(bad_a)}")
    if len(bad_a):
        ok = False
        print("   examples (face, owner, neighbour):",
              [(int(f), int(owner[f]), int(neigh[f])) for f in bad_a[:5]])

    # B. internal faces sorted by (owner, neighbour) lexicographically
    # OpenFOAM 要求的不只是 owner 单调不减：同一 owner 的多个面之间 neighbour 也必须递增，
    # 否则 checkMesh 会报 "Faces not in upper triangular order"。
    bad_b = np.where(
        (np.diff(owner[:n_int]) < 0)
        | ((owner[1:n_int] == owner[: n_int - 1]) & (np.diff(neigh) < 0))
    )[0] + 1
    print(f"B. internal (owner,neighbour) ordering violations: {len(bad_b)}")
    if len(bad_b):
        ok = False
        print("   at face index:", [int(x) for x in bad_b[:5]])

    # C. boundary block consistency
    prev_end = n_int
    for name, (sf, nb) in patches.items():
        good = sf == prev_end and all(owner[sf + k] >= 0 for k in range(nb))
        print(f"C. patch {name:14s} startFace={sf} nFaces={nb} {'OK' if good else 'MISMATCH'}")
        if not good:
            ok = False
        prev_end = sf + nb
    if prev_end != nf:
        ok = False
        print(f"C. boundary block end {prev_end} != nFaces {nf}")

    # cell centers (face-center average proxy, same as converter)
    n_cells = max(owner.max(), neigh.max() if n_int else 0) + 1
    fc = np.array([pts[ids].mean(axis=0) for ids in faces])
    cc_sum = np.zeros((n_cells, 3))
    cc_cnt = np.zeros(n_cells)
    for f in range(nf):
        cc_sum[owner[f]] += fc[f]
        cc_cnt[owner[f]] += 1
    for f in range(n_int):
        cc_sum[neigh[f]] += fc[f]
        cc_cnt[neigh[f]] += 1
    cc = cc_sum / np.maximum(cc_cnt, 1)[:, None]

    # D. geometric orientation
    bad_di = []
    for f in range(n_int):
        ids = faces[f]
        a = np.zeros(3)
        p = pts[ids]
        for k in range(len(p)):
            a += np.cross(p[k], p[(k + 1) % len(p)])
        a *= 0.5
        if float(np.dot(a, cc[neigh[f]] - cc[owner[f]])) < 0:
            bad_di.append(f)
    print(f"D1. internal faces normal not owner->neighbour: {len(bad_di)}")
    if bad_di:
        ok = False

    bad_db = []
    for f in range(n_int, nf):
        ids = faces[f]
        a = np.zeros(3)
        p = pts[ids]
        for k in range(len(p)):
            a += np.cross(p[k], p[(k + 1) % len(p)])
        a *= 0.5
        if float(np.dot(a, fc[f] - cc[owner[f]])) < 0:
            bad_db.append(f)
    print(f"D2. boundary faces normal not outward: {len(bad_db)}")
    if bad_db:
        ok = False

    # E. unused points (checkMesh 的 "Unused points found" 失败项)
    used = np.zeros(len(pts), dtype=bool)
    for ids in faces:
        used[ids] = True
    n_unused = int((~used).sum())
    print(f"E. unused points: {n_unused} / {len(pts)}")
    if n_unused:
        ok = False

    print("\nRESULT:", "CONFORMS" if ok else "VIOLATIONS FOUND")


if __name__ == "__main__":
    main(sys.argv[1])
