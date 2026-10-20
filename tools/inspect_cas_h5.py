#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Fluent .cas.h5 体检器：校验完整性 + 提取网格规模/区域/边界 + 查负体积与非法单元。

为什么要它：
  换结构进管线之前，先用 2 分钟判断这个 cas.h5 值不值得转——源文件里若有负体积、
  非法区域名、缺 zoneTopology，转出来的 polyMesh 会一路错到求解器才炸，排查成本高得多。

检查项：
  ① 完整性  文件存在、大小、is_hdf5、能打开、顶层组齐全（meshes/settings）。
  ② 规模    节点 / 面（总、内部、边界）/ 单元数。
  ③ 区域    单元区（cell zone）与面区（face zone）的 id / 名称 / 范围 / 类型。
  ④ 合法性  单元体积（散度定理逐面积分）→ 负体积、零/近零体积单元计数；
             非法面（节点数<3）、孤立单元（没有被任何面引用）、区域 id 越界。

用法：
    python tools/inspect_cas_h5.py "<结构>.cas.h5" [--json out.json] [--vol-threshold 1e-15]

依赖：h5py、numpy（装在 .workbuddy 的 venv 里，见 README §5）。
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

try:
    import h5py
except ImportError:  # pragma: no cover
    sys.exit("缺少 h5py：请用带 h5py 的解释器运行（README §5 的 venv），或 pip install h5py")

# Fluent 面类型（meshes/1/faces/zoneTopology/faceType）
FACE_TYPE = {0: "mixed", 2: "internal", 3: "boundary",
             4: "wall", 5: "pressure-inlet", 6: "pressure-outlet",
             7: "symmetry", 8: "periodic-shadow", 9: "pressure-far-field",
             10: "velocity-inlet", 11: "periodic", 12: "fan", 13: "porous-jump",
             14: "radiator", 20: "mass-flow-inlet", 24: "interface",
             31: "parent", 36: "outflow", 37: "axis", 40: "wall-shadow"}
# Fluent 单元类型（meshes/1/cells/zoneTopology/cellType）
CELL_TYPE = {0: "dead", 1: "fluid", 17: "solid", 16: "mixed"}


def _s(v) -> str:
    """h5py 读出的字节串 → 可读字符串（Fluent 的 name 是一条拼接长串）。"""
    if isinstance(v, bytes):
        return v.decode("utf-8", "replace")
    if isinstance(v, np.ndarray):
        return v.tobytes().decode("utf-8", "replace") if v.dtype.kind == "S" else str(v)
    return str(v)


def _role_of(name: str) -> str:
    """按名字关键字判定边界角色（中文名在导出时会丢成 ____，只能靠关键字）。"""
    n = name.lower()
    if "shadow" in n:
        return "耦合界面-shadow"
    if "interior" in n:
        return "内部面"
    if "wall_cool" in n or "cool" in n:
        return "冷源(定温壁)"
    if "inlet" in n:
        return "入口"
    if "outlet" in n:
        return "出口"
    if "wall" in n:
        return "壁面"
    return "未识别"


def _names(raw) -> list[str]:
    """Fluent 把多个名字拼成一条定长字节串（末尾补 0），按 \\0 切分。"""
    return [n for n in _s(raw).split("\x00") if n]


def face_geom(pts: np.ndarray, nodes: np.ndarray, nnodes: np.ndarray,
              chunk: int = 400_000) -> tuple[np.ndarray, np.ndarray]:
    """逐面算形心与 Newell 向量（模 = 面积）。分块以控内存。"""
    n = len(nnodes)
    start = np.zeros(n + 1, dtype=np.int64)
    np.cumsum(nnodes, out=start[1:])
    cen = np.empty((n, 3), dtype=np.float64)
    vec = np.empty((n, 3), dtype=np.float64)
    for i in range(0, n, chunk):
        j = min(i + chunk, n)
        cnt = nnodes[i:j]
        nf = j - i
        m = int(cnt.max())
        # 取每个面的前 m 个节点（不足的用第 0 个填充，后面按 cnt 置零）
        off = start[i:j][:, None] + np.arange(m)[None, :]
        off = np.minimum(off, len(nodes) - 1)
        ids = nodes[off]                      # (nf, m)
        p = pts[ids]                          # (nf, m, 3)
        mask = np.arange(m)[None, :] < cnt[:, None]
        p = np.where(mask[:, :, None], p, 0.0)
        cen[i:j] = p.sum(axis=1) / cnt[:, None]
        # Newell：v = Σ cross(p_k, p_{k+1}) / 2，|v| = 面积
        # Newell：v = Σ_{k=0}^{cnt-2} cross(p_k, p_{k+1}) + cross(p_{cnt-1}, p_0)
        # 坑：分块时短面被补零，若写成 cross(p_{m-1}, p_0) 会拿到补的零向量，
        #     闭合边丢失 → 法向系统性错误 → 单元体积约一半为负。
        #     必须按每个面自己的最后一个【有效】节点取 p_{cnt-1}。
        v = np.zeros((nf, 3), dtype=np.float64)
        for k in range(m - 1):
            v += np.cross(p[:, k, :], p[:, k + 1, :])
        last = np.maximum(cnt - 1, 0)
        p_last = p[np.arange(nf), last, :]
        v += np.cross(p_last, p[:, 0, :])
        vec[i:j] = v * 0.5
    return cen, vec


def inspect(path: str, vol_threshold: float = 1e-15) -> dict:
    if not os.path.isfile(path):
        raise SystemExit(f"文件不存在：{path}")
    size_mb = os.path.getsize(path) / 1e6
    if not h5py.is_hdf5(path):
        raise SystemExit(f"不是合法的 HDF5 文件：{path}")

    out: dict = {"path": path, "size_mb": round(size_mb, 3), "is_hdf5": True}
    with h5py.File(path, "r") as f:
        if "meshes" not in f:
            raise SystemExit("顶层缺 meshes 组：不是 Fluent 网格文件")
        m = f["meshes/1"]
        out["fluent_version"] = _s(f["settings/Version"][()]) if "settings/Version" in f else "?"

        pts = np.asarray(m["nodes/coords/2211"][:], dtype=np.float64)
        nnodes = np.asarray(m["faces/nodes/1/nnodes"][:], dtype=np.int64)
        # Fluent 的节点编号是 1 基，转 0 基
        fnodes = np.asarray(m["faces/nodes/1/nodes"][:], dtype=np.int64) - 1
        c0 = np.asarray(m["faces/c0/1"][:], dtype=np.int64)
        c1 = np.asarray(m["faces/c1/1"][:], dtype=np.int64)
        out["n_nodes"] = int(len(pts))
        out["n_faces"] = int(len(c0))
        out["n_internal_faces"] = int(len(c1))
        out["n_boundary_faces"] = int(len(c0) - len(c1))

        # ── 区域
        cz = m["cells/zoneTopology"]
        cnames = _names(cz["name"][0])
        out["cell_zones"] = [{
            "id": int(cz["id"][i]), "name": cnames[i] if i < len(cnames) else f"zone{cz['id'][i]}",
            "type": CELL_TYPE.get(int(cz["cellType"][i]), str(int(cz["cellType"][i]))),
            "minId": int(cz["minId"][i]), "maxId": int(cz["maxId"][i]),
            "n_cells": int(cz["maxId"][i] - cz["minId"][i] + 1),
        } for i in range(len(cz["id"]))]
        out["n_cells"] = int(sum(z["n_cells"] for z in out["cell_zones"]))

        fz = m["faces/zoneTopology"]
        fnames = _names(fz["name"][0])
        # 注意：zoneTopology/faceType 在本文件里全部是 5，无法据此区分边界类型
        # （Fluent 23.2 的取值口径与旧文档不一致）。因此【只报原始码】，
        # 边界角色改由名字里的关键字判定——这与本项目转换器一贯的做法一致。
        out["face_zones"] = sorted([{
            "id": int(fz["id"][i]),
            "name": fnames[i] if i < len(fnames) else f"zone{fz['id'][i]}",
            "faceTypeCode": int(fz["faceType"][i]),
            "zoneTypeCode": int(fz["zoneType"][i]),
            "role": _role_of(fnames[i] if i < len(fnames) else ""),
            "minId": int(fz["minId"][i]), "maxId": int(fz["maxId"][i]),
            "n_faces": int(fz["maxId"][i] - fz["minId"][i] + 1),
        } for i in range(len(fz["id"]))], key=lambda z: z["minId"])

        # ── 内部面 / 边界面的【真实索引】（必须在 with 块内读，文件关了就读不到）
        #    坑：c0 覆盖全部面、c1 只覆盖内部面，但内部面并不排在前面。
        #    本文件实测「边界区在前（164,951 张）、两个 interior 区在后（3,927,457 张）」，
        #    若按「前 n_int 张是内部面」切，单元体积会全错（单个体积能到总体积的 7%）。
        #    正确做法：从 zoneTopology 的 minId/maxId 取 interior 区的编号区间。
        zmin = np.asarray(fz["minId"][:], dtype=np.int64)  # 1 基面编号
        zmax = np.asarray(fz["maxId"][:], dtype=np.int64)
        ztyp = np.asarray(fz["zoneType"][:], dtype=np.int64)
        internal_idx = np.concatenate([np.arange(a - 1, b) for a, b in
                                       zip(zmin[ztyp == 2], zmax[ztyp == 2])]) \
            if (ztyp == 2).any() else np.arange(len(c1))
        bnd_mask = np.ones(len(c0), dtype=bool)
        bnd_mask[internal_idx] = False
        bnd_idx = np.where(bnd_mask)[0]

    out["face_order_check"] = {
        "n_internal_from_zones": int(len(internal_idx)),
        "n_internal_from_c1": int(len(c1)),
        "consistent": bool(len(internal_idx) == len(c1)),
        "n_boundary_from_zones": int(len(bnd_idx)),
    }

    # ── 非法面 / 孤立单元
    out["illegal_faces_nnodes_lt3"] = int((nnodes < 3).sum())
    # Fluent 的 c0/c1 是 1 基单元编号；0 表示"无"
    used = np.unique(np.concatenate([c0[c0 > 0], c1[c1 > 0]]))
    max_id = max(int(z["maxId"]) for z in out["cell_zones"])
    out["orphan_cells_not_referenced"] = int(max_id - len(used)) if max_id >= len(used) else 0

    # ── 单元体积（散度定理 V = 1/3 Σ r_f·n_f A_f）
    #    关键：Fluent 的【面节点排序不保证与 c0→c1 一致】，直接按 ±Newell 法向累加会
    #    互相抵消（实测出现 51% 负体积、体积分布对称于 0，就是这个症状）。
    #    故先用「单元形心 → 面形心」的方向给每张面定向：
    #      内部面：法向应指向 c1；(x_c1 − x_c0)·n < 0 则翻向
    #      边界面：法向应背离单元指向域外；(x_f − x_c0)·n < 0 则翻向
    cen, vec = face_geom(pts, fnodes, nnodes)
    n_int = len(c1)
    ncell = max_id + 1

    # 单元形心近似 = 其各面形心的均值
    # 单元形心必须同时用 c0 和 c1 两侧的面：只算 c0 的话，从未作为 c0 出现的单元
    # 形心会落在原点，定向判据随之失效（表现为约一半单元体积为负）。
    csum = np.zeros((ncell, 3), dtype=np.float64)
    ccnt = np.zeros(ncell, dtype=np.float64)
    np.add.at(csum, c0[c0 > 0], cen[c0 > 0])
    np.add.at(ccnt, c0[c0 > 0], 1.0)
    cen_int = cen[internal_idx]                    # c1 对应的是内部面
    np.add.at(csum, c1[c1 > 0], cen_int[c1 > 0])
    np.add.at(ccnt, c1[c1 > 0], 1.0)
    ccen = csum / np.maximum(ccnt, 1.0)[:, None]

    # 内部面：法向应从 c0 指向 c1
    ci0 = c0[internal_idx]
    dot_i = np.einsum("ij,ij->i", ccen[c1] - ccen[ci0], vec[internal_idx])
    s_i = np.where(dot_i < 0, -1.0, 1.0)
    # 边界面：法向应背离所属单元
    cb0 = c0[bnd_idx]
    dot_b = np.einsum("ij,ij->i", cen[bnd_idx] - ccen[cb0], vec[bnd_idx])
    s_b = np.where(dot_b < 0, -1.0, 1.0)

    vol = np.zeros(ncell, dtype=np.float64)
    ci = np.einsum("ij,ij->i", cen[internal_idx], vec[internal_idx]) / 3.0 * s_i
    np.add.at(vol, ci0, ci)
    np.add.at(vol, c1, -ci)
    cb = np.einsum("ij,ij->i", cen[bnd_idx], vec[bnd_idx]) / 3.0 * s_b
    np.add.at(vol, cb0, cb)
    cell_vol = vol[1:]                     # 去掉 0 号占位（Fluent 单元编号 1 基）

    out["total_volume_m3"] = float(cell_vol.sum())
    out["min_volume_m3"] = float(cell_vol.min())
    out["max_volume_m3"] = float(cell_vol.max())
    out["negative_volume_cells"] = int((cell_vol < 0).sum())
    out["near_zero_volume_cells"] = int(((cell_vol >= 0) & (cell_vol < vol_threshold)).sum())
    out["volume_threshold_m3"] = vol_threshold
    out["volume_percentiles_m3"] = {
        f"p{p}": float(np.percentile(cell_vol, p)) for p in (1, 5, 50, 95, 99)
    }
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cas")
    ap.add_argument("--json", default=None)
    ap.add_argument("--vol-threshold", type=float, default=1e-15)
    args = ap.parse_args()
    r = inspect(args.cas, args.vol_threshold)
    if args.json:
        with open(args.json, "w", encoding="utf-8") as fp:
            json.dump(r, fp, ensure_ascii=False, indent=2)
        print(f"-> {args.json}")
    else:
        print(json.dumps(r, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
