#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
 Fluent .cas.h5 网格 → OpenFOAM polyMesh 转换器（流体区提取版）
================================================================================
 把 ANSYS Fluent 2023 R2 的 .cas.h5 网格转成 OpenFOAM 的 constant/polyMesh。
 针对本工程「共轭传热（固体铝 + 流体湿空气）」结构，默认【只提取流体区】，
 并把「固/流耦合界面（带 shadow 的面区）」作为冷壁 coldWall 导出
 （对应 wall_cool = T_c 定壁温的简化：忽略铝结构导热，界面即冷面）。

 用法（用工程内置 Python，勿用系统 anaconda）：
   C:/Users/LELE/.workbuddy/binaries/python/envs/default/Scripts/python.exe \
       tools/fluent_cas_to_foam.py <case.cas.h5> <output_case_dir>

 输出：<output_case_dir>/constant/polyMesh/{points,faces,owner,neighbour,boundary}

 已实测的格式要点（Fluent 2023 R2 HDF5）：
   * nodes/coords/<key>  (Nn,3) float64，node id = index+1
   * faces/nodes/1/nnodes (Nf,)；faces/nodes/1/nodes 扁平 1-based
   * faces/c0/1 (Nf,)、faces/c1/1 (Ni,) —— c1 为紧凑数组：
     边界面 fid ∈ [1, nBoundary] 在前、内部面 fid ∈ [nBoundary+1, Nf] 在后，
     c1 索引 = fid-1-nBoundary
   * faces/zoneTopology：id/minId/maxId/zoneType/shadowZoneId/c0(所属cell zone)/name(';'分隔)
================================================================================
"""

from __future__ import annotations

import argparse
import os
from collections import OrderedDict

import h5py
import numpy as np

ZT_INTERIOR = 2
ZT_WALL = 3
ZT_PRESSURE_OUTLET = 5
ZT_MASS_FLOW_INLET = 20


def fo_header(cls: str, obj: str, loc: str = "constant/polyMesh") -> str:
    return (
        "/*--------------------------------*- C++ -*----------------------------------*\\\n"
        "| =========                 |                                                 |\n"
        "| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n"
        "|  \\\\    /   O peration     | Version:  v10                                   |\n"
        "|   \\\\  /    A nd           | Website:  www.openfoam.com                      |\n"
        "|    \\\\/     M anipulation  |                                                 |\n"
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


def convert(cas_path: str, out_dir: str) -> None:
    h = h5py.File(cas_path, "r")
    try:
        m = h["meshes/1"]

        # ---------- 读网格 ----------
        ckey = list(m["nodes/coords"].keys())[0]
        coords = np.asarray(m["nodes/coords"][ckey][:], dtype=np.float64)   # (Nn,3)

        nnodes = np.asarray(m["faces/nodes/1/nnodes"][:], dtype=np.int64)   # (Nf,)
        fnodes = np.asarray(m["faces/nodes/1/nodes"][:], dtype=np.int64)    # flat, 1-based
        offs = np.zeros(len(nnodes) + 1, dtype=np.int64)
        np.cumsum(nnodes, out=offs[1:])

        c0 = np.asarray(m["faces/c0/1"][:], dtype=np.int64)   # (Nf,) 1-based
        c1 = np.asarray(m["faces/c1/1"][:], dtype=np.int64)   # (Ni,) 1-based
        nf = len(c0)
        n_boundary = nf - len(c1)

        fz = m["faces/zoneTopology"]
        fz_id = np.asarray(fz["id"][:], dtype=np.int64)
        fz_min = np.asarray(fz["minId"][:], dtype=np.int64)
        fz_max = np.asarray(fz["maxId"][:], dtype=np.int64)
        fz_type = np.asarray(fz["zoneType"][:], dtype=np.int64)
        fz_shadow = np.asarray(fz["shadowZoneId"][:], dtype=np.int64)
        fz_c0 = np.asarray(fz["c0"][:], dtype=np.int64)
        fz_name = fz["name"][0].decode("latin-1").split(";")

        cz = m["cells/zoneTopology"]
        cz_id = np.asarray(cz["id"][:], dtype=np.int64)
        cz_min = np.asarray(cz["minId"][:], dtype=np.int64)
        cz_max = np.asarray(cz["maxId"][:], dtype=np.int64)
        cz_name = cz["name"][0].decode("latin-1").split(";")

        # ---------- 识别流体 cell zone ----------
        fluid_zone = None
        for i, t in enumerate(fz_type):
            if t in (ZT_MASS_FLOW_INLET, ZT_PRESSURE_OUTLET):
                fluid_zone = int(fz_c0[i])
                break
        if fluid_zone is None:
            raise SystemExit("未找到 inlet/outlet，无法识别流体 cell zone")

        zi = int(np.where(cz_id == fluid_zone)[0][0])
        fluid_lo, fluid_hi = int(cz_min[zi]), int(cz_max[zi])
        n_cells_total = int(cz_max.max())
        print(f"流体 cell zone id={fluid_zone}  单元 {fluid_lo}..{fluid_hi} "
              f"({fluid_hi - fluid_lo + 1} 个)")
        print(f"单元分区：{list(zip(cz_id.tolist(), cz_name))}")
        print("面分区：")
        for i in range(len(fz_id)):
            print(f"  {fz_name[i]:30s} type={fz_type[i]:2d} c0zone={fz_c0[i]:3d} "
                  f"shadow={fz_shadow[i]:3d} fid {fz_min[i]}..{fz_max[i]}")

        # ---------- cell id 重映射（只保留流体单元）----------
        cell_map = np.full(n_cells_total + 1, -1, dtype=np.int64)   # 1-based id -> 0-based
        cell_map[np.arange(fluid_lo, fluid_hi + 1, dtype=np.int64)] = \
            np.arange(fluid_hi - fluid_lo + 1, dtype=np.int64)

        # ---------- 完整 c1（内部面才有值）----------
        c1_full = np.zeros(nf, dtype=np.int64)
        c1_full[n_boundary:] = c1        # 内部面 fid = n_boundary+1 .. nf  (0-based 索引 = fid-1)

        # ---------- 面分类 ----------
        keep = cell_map[c0] >= 0                 # c0 是流体单元
        is_internal = keep & (cell_map[c1_full] >= 0)   # c1 也是流体单元 -> 流体内部面
        # 注：c0=流体、c1=固体的面（固/流界面）会自动落为边界面（is_internal=False）

        kept_idx = np.where(keep)[0]
        int_idx = np.where(is_internal[kept_idx])[0]     # 相对 kept 列表的内部面位置
        bnd_idx = np.where(~is_internal[kept_idx])[0]
        n_int = len(int_idx)

        # ---------- owner / neighbour ----------
        owner = cell_map[c0[kept_idx]]
        neigh_full = cell_map[c1_full[kept_idx]]   # 全长，边界面 = -1（受 internal_mask 保护）

        # ---------- 面-节点（0-based）----------
        faces = []
        for idx in kept_idx:
            faces.append(fnodes[offs[idx]:offs[idx + 1]] - 1)

        internal_kept = is_internal[kept_idx]

        # ---------- 标号约定（OpenFOAM：内部面 owner < neighbour）----------
        # Fluent 只保证面法向 c0->c1，不保证 c0 编号小于 c1。
        # 对 owner > neighbour 的内部面：交换 owner/neighbour，并反转节点序，
        # 使法向仍从（新）owner 指向（新）neighbour。
        swap = internal_kept & (owner > neigh_full)
        owner, neigh_full = np.where(swap, neigh_full, owner), np.where(swap, owner, neigh_full)
        for i in np.where(swap)[0]:
            faces[i] = faces[i][::-1]
        print(f"owner/neighbour 交换的内部面：{int(swap.sum())}")

        # ---------- 重排：内部面在前，边界面在后 ----------
        # OpenFOAM 要求内部面按 (owner, neighbour) 字典序排列，而不只是 owner 单调不减。
        # 只按 owner 排序会让「同一 owner 的多个面之间 neighbour 乱序」，
        # checkMesh 会报 "Faces not in upper triangular order"（实测 v8_1 达 127895 个）。
        # 因此这里用 lexsort：主键 owner、次键 neighbour。
        int_pos = np.where(internal_kept)[0]
        bnd_pos = np.where(~internal_kept)[0]
        int_pos = int_pos[np.lexsort((neigh_full[int_pos], owner[int_pos]))]
        order = np.concatenate([int_pos, bnd_pos])
        n_int = len(int_pos)
        faces_o = [faces[i] for i in order]
        owner_o = owner[order]
        neighbour_o = neigh_full[order[:n_int]]

        # ---------- 朝向校正（OpenFOAM：内部面法向 owner->neighbour；边界面朝外）----------
        # 上面交换后内部面朝向应已正确；此处为几何兜底校验/修正。
        orient_mask = np.zeros(len(order), dtype=bool)
        orient_mask[:n_int] = True
        faces_o = orient_faces(coords, faces_o, owner_o, neighbour_o, orient_mask)
        faces_o = [list(map(int, ids)) for ids in faces_o]
        owner_o = owner_o.tolist()
        neighbour_o = neighbour_o.tolist()

        # ---------- 边界面 patch 命名 ----------
        orig_fid = kept_idx + 1                       # 原 1-based face id
        patch_of_face = np.array(["wall"] * nf, dtype=object)
        for i in range(len(fz_id)):
            lo, hi = fz_min[i], fz_max[i]
            nm = fz_name[i]
            if nm == "inlet":
                p = "inlet"
            elif nm == "outlet":
                p = "outlet"
            elif fz_c0[i] == fluid_zone and fz_type[i] == ZT_WALL:
                # 固/流耦合界面：流体侧、且带 shadow -> 冷壁
                p = "coldWall" if fz_shadow[i] != 0 else "wall"
            else:
                p = "__drop__"                          # 非流体区，不导出
            patch_of_face[lo - 1:hi] = p

        # 边界面的 patch（按原顺序，保持连续）
        bnd_patches = OrderedDict()
        for b in bnd_idx:
            nm = patch_of_face[orig_fid[b] - 1]
            bnd_patches.setdefault(nm, []).append(b)

        # 校验：所有导出的边界面都要有明确 patch
        assert "__drop__" not in bnd_patches, "存在不应导出的边界面"

        # ---------- 写文件 ----------
        d = os.path.join(out_dir, "constant", "polyMesh")
        os.makedirs(d, exist_ok=True)

        # ---------- 清理未使用的节点 ----------
        # Fluent 的节点表包含固体区节点，只提取流体区后大部分节点不再被任何面引用。
        # OpenFOAM 的 checkMesh 会报 "Unused points found"（实测 v8_1 有 316832 个，
        # 占全部节点的一半以上），既拖慢求解又污染网格体积统计，故在此压缩重编号。
        used = np.zeros(coords.shape[0], dtype=bool)
        for ids in faces_o:
            used[ids] = True
        if not used.all():
            remap = np.full(coords.shape[0], -1, dtype=np.int64)
            remap[used] = np.arange(int(used.sum()), dtype=np.int64)
            print(f"清理未使用节点：{coords.shape[0]} -> {int(used.sum())} "
                  f"(丢弃 {coords.shape[0] - int(used.sum())})")
            coords = coords[used]
            faces_o = [remap[ids].tolist() for ids in faces_o]

        with open(os.path.join(d, "points"), "w") as f:
            f.write(fo_header("vectorField", "points"))
            f.write(f"{coords.shape[0]}\n(\n")
            for p in coords:
                f.write(f"({p[0]:.10g} {p[1]:.10g} {p[2]:.10g})\n")
            f.write(")\n")

        with open(os.path.join(d, "faces"), "w") as f:
            f.write(fo_header("faceList", "faces"))
            f.write(f"{len(faces_o)}\n(\n")
            for ids in faces_o:
                f.write(f"{len(ids)}(" + " ".join(str(x) for x in ids) + ")\n")
            f.write(")\n")

        with open(os.path.join(d, "owner"), "w") as f:
            f.write(fo_header("labelList", "owner"))
            f.write(f"{len(owner_o)}\n(\n")
            for o in owner_o:
                f.write(f"{o}\n")
            f.write(")\n")

        with open(os.path.join(d, "neighbour"), "w") as f:
            f.write(fo_header("labelList", "neighbour"))
            f.write(f"{len(neighbour_o)}\n(\n")
            for nb in neighbour_o:
                f.write(f"{nb}\n")
            f.write(")\n")

        start = n_int
        lines = []
        for nm, blist in bnd_patches.items():
            ty = "patch" if nm in ("inlet", "outlet") else "wall"
            lines.append((nm, ty, start, len(blist)))
            start += len(blist)
        with open(os.path.join(d, "boundary"), "w") as f:
            f.write(fo_header("polyBoundaryMesh", "boundary"))
            f.write(f"{len(lines)}\n(\n")
            for nm, ty, s, n in lines:
                f.write(f"    {nm}\n    {{\n        type            {ty};\n"
                        f"        nFaces          {n};\n"
                        f"        startFace       {s};\n    }}\n")
            f.write(")\n")

        print("\n===== 转换完成 =====")
        print(f"节点：{coords.shape[0]}  单元：{int(max(owner_o)) + 1}")
        print(f"面：{len(faces_o)}（内部 {n_int} + 边界 {len(bnd_idx)}）")
        print("边界分区：")
        for nm, ty, s, n in lines:
            print(f"  {nm:14s} {ty:6s} nFaces={n:6d} startFace={s}")
        print(f"输出：{d}")
    finally:
        h.close()


def orient_faces(coords, faces, owner, neighbour, internal_mask):
    """按 OpenFOAM 约定校正朝向：内部面法向 owner->neighbour；边界面法向朝外。"""
    nf = len(faces)
    area = np.zeros((nf, 3))
    fc = np.zeros((nf, 3))
    for i, ids in enumerate(faces):
        p = coords[ids]
        a = np.zeros(3)
        for j in range(len(p)):
            a += np.cross(p[j], p[(j + 1) % len(p)])
        a *= 0.5
        area[i] = a
        fc[i] = p.mean(axis=0)

    n_cells = int(owner.max()) + 1
    cc_sum = np.zeros((n_cells, 3))
    cc_cnt = np.zeros(n_cells, dtype=np.int64)
    for i in range(nf):
        cc_sum[owner[i]] += fc[i]
        cc_cnt[owner[i]] += 1
        if internal_mask[i]:
            cc_sum[neighbour[i]] += fc[i]
            cc_cnt[neighbour[i]] += 1
    cc_cnt[cc_cnt == 0] = 1
    cc = cc_sum / cc_cnt[:, None]

    for i in range(nf):
        d = (cc[neighbour[i]] - cc[owner[i]]) if internal_mask[i] else (fc[i] - cc[owner[i]])
        if float(np.dot(area[i], d)) < 0.0:
            faces[i] = faces[i][::-1]
    return faces


def main():
    ap = argparse.ArgumentParser(description="Fluent .cas.h5 -> OpenFOAM polyMesh")
    ap.add_argument("cas", help=".cas.h5 路径")
    ap.add_argument("out", help="输出算例目录（写入 <out>/constant/polyMesh）")
    args = ap.parse_args()
    if not os.path.isfile(args.cas):
        raise SystemExit(f"文件不存在：{args.cas}")
    convert(args.cas, args.out)


if __name__ == "__main__":
    main()
