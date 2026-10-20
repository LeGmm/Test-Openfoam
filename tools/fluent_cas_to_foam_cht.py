#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
 Fluent .cas.h5 → chtMultiRegionFoam 双区域网格（fluid + solid）
================================================================================
 与 fluent_cas_to_foam.py（只取流体区）不同，本脚本把【流体和固体两个 cell zone】
 分别导出为两个区域网格，供两步共轭方案的 Step 1（chtMultiRegionFoam）使用：
   <out>/constant/fluid/polyMesh   流体区
        inlet / outlet / fluid_to_solid(mappedWall) / wall
   <out>/constant/solid/polyMesh   固体区
        solid_to_fluid(mappedWall) / coldWall(wall_cool, 5℃ 面) / outer(其余外表面)

 格式要点（同 fluent_cas_to_foam.py，均已实测）：
   * nodes/coords  节点坐标，node id = index+1
   * faces/nodes/1/{nnodes,nodes}  面-节点连通，fnodes 1-based
   * faces/c0/1 (Nf,)、faces/c1/1 (Ni,)  边界面 fid ∈ [1, nBoundary] 在前、
     内部面 fid ∈ [nBoundary+1, Nf] 在后（c1 索引 = fid-1-nBoundary）
   * faces/zoneTopology  name 列是分号连接串（1 行），其余列逐区行
   * 每个面的归属区 = 其 zone 行的 c0（cell zone id）；耦合面两侧各有一条
     独立面带（shadow≠0），正好成为两个区域各自的边界面
 用法：
   python tools/fluent_cas_to_foam_cht.py <case.cas.h5> <out_case_dir>
================================================================================
"""
from __future__ import annotations

import argparse
import os

import h5py
import numpy as np


def header(cls: str, loc: str, obj: str) -> str:
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
        f"    location    \"{loc}\";\n"
        f"    object      {obj};\n"
        "}\n"
        "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n"
    )


def write_labels(path: str, values) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write(header("labelList", "constant/polyMesh", os.path.basename(path)))
        f.write(f"\n{len(values)}\n(\n")
        for v in values:
            f.write(f"{v}\n")
        f.write(")\n\n// ************************************************************************* //\n")


def write_faces(path: str, faces) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write(header("faceList", "constant/polyMesh", "faces"))
        f.write(f"\n{len(faces)}\n(\n")
        for ids in faces:
            f.write(f"{len(ids)}(" + " ".join(str(i) for i in ids) + ")\n")
        f.write(")\n\n// ************************************************************************* //\n")


def write_points(path: str, pts) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write(header("vectorField", "constant/polyMesh", "points"))
        f.write(f"\n{len(pts)}\n(\n")
        for p in pts:
            f.write(f"({p[0]:.10g} {p[1]:.10g} {p[2]:.10g})\n")
        f.write(")\n\n// ************************************************************************* //\n")


def write_boundary(path: str, patches) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write(header("polyBoundaryMesh", "constant/polyMesh", "boundary"))
        f.write(f"\n{len(patches)}\n(\n")
        for name, p in patches.items():
            f.write(f"    {name}\n    {{\n")
            f.write(f"        type            {p['type']};\n")
            if p.get("mapped"):
                f.write(f"        inGroups        1 (wall);\n")
            f.write(f"        nFaces          {p['nFaces']};\n")
            f.write(f"        startFace       {p['startFace']};\n")
            for k, v in p.get("extra", []):
                f.write(f"        {k:16s}{v};\n")
            f.write("    }\n")
        f.write(")\n\n// ************************************************************************* //\n")


def orient_faces(coords, faces, owner, neigh, internal_mask):
    """几何兜底：内部面法向 owner→neighbour，边界面朝外。
    单元中心由面形心累计得到（与 fluent_cas_to_foam.py 已验证实现一致）。"""
    nf = len(faces)
    area = np.zeros((nf, 3))
    fc = np.zeros((nf, 3))
    for i, ids in enumerate(faces):
        p = coords[ids]
        a = np.zeros(3)
        for j in range(len(p)):
            a += np.cross(p[j], p[(j + 1) % len(p)])
        area[i] = a * 0.5
        fc[i] = p.mean(axis=0)

    n_cells = int(max(owner.max(), neigh.max() if len(neigh) else 0)) + 1
    cc_sum = np.zeros((n_cells, 3))
    cc_cnt = np.zeros(n_cells, dtype=np.int64)
    for i in range(nf):
        cc_sum[owner[i]] += fc[i]
        cc_cnt[owner[i]] += 1
        if internal_mask[i]:
            cc_sum[neigh[i]] += fc[i]
            cc_cnt[neigh[i]] += 1
    cc_cnt[cc_cnt == 0] = 1
    cc = cc_sum / cc_cnt[:, None]

    faces = [list(ids) for ids in faces]
    for i in range(nf):
        d = (cc[neigh[i]] - cc[owner[i]]) if internal_mask[i] else (fc[i] - cc[owner[i]])
        if float(np.dot(area[i], d)) < 0.0:
            faces[i] = faces[i][::-1]
    return faces


def convert(cas_path: str, out_dir: str) -> None:
    h = h5py.File(cas_path, "r")
    m = h["meshes/1"]

    ckey = list(m["nodes/coords"].keys())[0]
    coords = np.asarray(m["nodes/coords"][ckey][:], dtype=np.float64)
    nnodes = np.asarray(m["faces/nodes/1/nnodes"][:], dtype=np.int64)
    fnodes = np.asarray(m["faces/nodes/1/nodes"][:], dtype=np.int64)
    offs = np.zeros(len(nnodes) + 1, dtype=np.int64)
    np.cumsum(nnodes, out=offs[1:])
    c0 = np.asarray(m["faces/c0/1"][:], dtype=np.int64)
    c1 = np.asarray(m["faces/c1/1"][:], dtype=np.int64)
    nf = len(c0)
    n_boundary = nf - len(c1)
    has_c1 = np.zeros(nf, dtype=bool)
    has_c1[n_boundary:] = True

    fz = m["faces/zoneTopology"]
    fz_id = np.asarray(fz["id"][:], dtype=np.int64)
    fz_min = np.asarray(fz["minId"][:], dtype=np.int64)
    fz_max = np.asarray(fz["maxId"][:], dtype=np.int64)
    fz_c0 = np.asarray(fz["c0"][:], dtype=np.int64)
    fz_shadow = np.asarray(fz["shadowZoneId"][:], dtype=np.int64)
    fz_name = fz["name"][0].decode("latin-1").split(";")

    cz = m["cells/zoneTopology"]
    cz_id = [int(z) for z in cz["id"][:]]
    cz_min = np.asarray(cz["minId"][:], dtype=np.int64)
    cz_max = np.asarray(cz["maxId"][:], dtype=np.int64)
    cz_name = cz["name"][0].decode("latin-1").split(";")

    # ---- 流体/固体 cell zone ----
    fluid_zone = None
    for i in range(len(fz_id)):
        if fz_name[i] == "inlet":
            fluid_zone = int(fz_c0[i])
            break
    if fluid_zone is None:
        raise SystemExit("未找到 inlet")
    others = [z for z in cz_id if z != fluid_zone]
    if len(others) != 1:
        raise SystemExit(f"预期恰好 2 个 cell zone，实际 {len(cz_id)}")
    solid_zone = others[0]
    zones = {"fluid": fluid_zone, "solid": solid_zone}
    print(f"cell zones: fluid={fluid_zone} ({cz_name[cz_id.index(fluid_zone)]}), "
          f"solid={solid_zone} ({cz_name[cz_id.index(solid_zone)]})")

    # ---- face → zone 行归属（各行不相交且铺满全部面，已实测）----
    zone_row_of_face = np.full(nf, -1, dtype=np.int64)
    for i in range(len(fz_id)):
        zone_row_of_face[fz_min[i] - 1: fz_max[i]] = i

    def patch_of(row: int):
        nm = fz_name[row]
        c0z = int(fz_c0[row])
        sh = int(fz_shadow[row])
        if c0z == fluid_zone:
            if nm == "inlet":
                return "fluid", "inlet", "patch"
            if nm == "outlet":
                return "fluid", "outlet", "patch"
            if sh != 0:
                return "fluid", "fluid_to_solid", "mappedWall"
            return "fluid", "wall", "wall"
        if sh != 0:
            return "solid", "solid_to_fluid", "mappedWall"
        if "cool" in nm.lower():
            return "solid", "coldWall", "wall"
        return "solid", "outer", "wall"

    for region in ("fluid", "solid"):
        zid = zones[region]
        zi = cz_id.index(zid)
        lo, hi = int(cz_min[zi]), int(cz_max[zi])
        n_cells = hi - lo + 1
        cell_map = np.full(nf + 1, -1, dtype=np.int64)   # c0/c1 是 1-based cell id
        cell_map[np.arange(lo, hi + 1, dtype=np.int64)] = np.arange(n_cells, dtype=np.int64)

        rows = [i for i in range(len(fz_id)) if int(fz_c0[i]) == zid]
        idx = np.concatenate([np.arange(fz_min[i] - 1, fz_max[i]) for i in rows])
        kept = np.sort(idx)
        owner = cell_map[c0[kept]]
        neigh = np.where(has_c1[kept], cell_map[c1[kept - n_boundary]], -1)
        internal = neigh >= 0

        swap = internal & (owner > neigh)
        owner, neigh = (np.where(swap, neigh, owner),
                        np.where(swap, owner, neigh))

        faces = []
        for k, fid in enumerate(kept):
            ids = fnodes[offs[fid]: offs[fid + 1]] - 1
            faces.append(ids[::-1] if swap[k] else ids)
        print(f"[{region}] 单元 {n_cells}  面数 {len(kept)}  内部 {int(internal.sum())}  "
              f"边界 {int((~internal).sum())}  交换 {int(swap.sum())}")

        # 内部面 (owner,neighbour) 字典序；边界面按补丁分组排序 → 每补丁连续
        int_pos = np.where(internal)[0]
        int_pos = int_pos[np.lexsort((neigh[int_pos], owner[int_pos]))]
        bnd_pos = np.where(~internal)[0]
        patch_order = (["inlet", "outlet", "fluid_to_solid", "wall"] if region == "fluid"
                       else ["solid_to_fluid", "coldWall", "outer"])
        patch_of_bnd = [patch_of(int(zone_row_of_face[kept[k]])) for k in bnd_pos]
        bnd_sorted = sorted(range(len(bnd_pos)),
                            key=lambda j: (patch_order.index(patch_of_bnd[j][1]), j))
        bnd_pos = bnd_pos[np.array(bnd_sorted, dtype=np.int64)]

        order = np.concatenate([int_pos, bnd_pos])
        n_int = len(int_pos)
        faces_o = [faces[k] for k in order]
        owner_o = owner[order]
        neigh_o = neigh[order]
        mask = np.zeros(len(order), dtype=bool)
        mask[:n_int] = True
        faces_o = orient_faces(coords, faces_o, owner_o, neigh_o, mask)
        faces_o = [list(map(int, ids)) for ids in faces_o]
        owner_o = owner_o.tolist()

        used = np.zeros(coords.shape[0], dtype=bool)
        for ids in faces_o:
            used[ids] = True
        remap = np.full(coords.shape[0], -1, dtype=np.int64)
        remap[used] = np.arange(int(used.sum()), dtype=np.int64)
        pts = coords[used]
        faces_o = [[int(remap[i]) for i in ids] for ids in faces_o]
        print(f"[{region}] 节点 {coords.shape[0]} → 使用 {int(used.sum())}")

        # 补丁统计（边界段在 order 中已是连续分组）
        patches = {}
        start = n_int
        for pname in patch_order:
            members = [j for j in range(n_int, len(order))
                       if patch_of(int(zone_row_of_face[kept[order[j]]]))[1] == pname]
            if not members:
                continue
            ptype = patch_of(int(zone_row_of_face[kept[order[members[0]]]]))[2]
            assert members == list(range(members[0], members[0] + len(members))), \
                f"补丁 {pname} 不连续"
            entry = {"type": ptype, "nFaces": len(members), "startFace": members[0]}
            if ptype == "mappedWall":
                entry["mapped"] = True
                if pname == "fluid_to_solid":
                    entry["extra"] = [("sampleMode", "nearestPatchFace"),
                                      ("sampleRegion", "solid"),
                                      ("samplePatch", "solid_to_fluid")]
                else:
                    entry["extra"] = [("sampleMode", "nearestPatchFace"),
                                      ("sampleRegion", "fluid"),
                                      ("samplePatch", "fluid_to_solid")]
            patches[pname] = entry
            start += len(members)

        pd = os.path.join(out_dir, "constant", region, "polyMesh")
        os.makedirs(pd, exist_ok=True)
        write_points(os.path.join(pd, "points"), pts)
        write_faces(os.path.join(pd, "faces"), faces_o)
        write_labels(os.path.join(pd, "owner"), owner_o)
        write_labels(os.path.join(pd, "neighbour"), [int(x) for x in neigh_o[:n_int]])
        write_boundary(os.path.join(pd, "boundary"), patches)
        print(f"[{region}] 写出 {pd}  patches={list(patches.keys())}")

    h.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cas")
    ap.add_argument("out")
    a = ap.parse_args()
    convert(a.cas, a.out)
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
