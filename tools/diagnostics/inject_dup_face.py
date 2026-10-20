#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
inject_dup_face.py —— 向 polyMesh 注入「同对单元双面连接」缺陷（根因判定实验用）

用法:
    python inject_dup_face.py <polyMesh目录>

操作: 复制排序在最后的一个内部面（face[n_int-1] 的 owner/neighbour/节点序），
     插到内部面块末尾（index n_int），使该对单元之间出现 2 个面。
     由于插在排序末尾，(owner,neighbour) 字典序保持非严格递增，
     只需把所有 boundary 补丁的 startFace +1，并更新 faces/neighbour 的数量头。
"""
import os
import re
import sys


def header(cls, obj):
    return (
        "/*--------------------------------*- C++ -*----------------------------------*\\\n"
        "| =========                 |                                                 |\n"
        "|  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n"
        "|   \\\\    /   O peration    | Version:  8                                     |\n"
        "|    \\\\  /    A nd           | Website:  https://openfoam.org                  |\n"
        "|     \\\\ /     M anipulation  |                                                 |\n"
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


def read_ascii_list(p, cast=int):
    txt = open(p, encoding="utf-8").read()
    body = txt.split("\n(", 1)[1].rsplit(")\n", 1)[0]
    return [cast(x) for x in body.split()]


def read_ascii_faces(p):
    txt = open(p, encoding="utf-8").read()
    body = txt.split("\n(", 1)[1].rsplit(")\n", 1)[0]
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
    body = txt.split("\n(", 1)[1].rsplit("\n)", 1)[0]
    patches = []
    for m in re.finditer(r"\n\s*(\w+)\s*\n\s*\{([^}]*)\}", body, re.S):
        blk = m.group(2)
        patches.append({
            "name": m.group(1),
            "startFace": int(re.search(r"startFace\s+(\d+)", blk).group(1)),
            "nFaces": int(re.search(r"nFaces\s+(\d+)\s*;", blk).group(1)),
            "type": re.search(r"type\s+(\w+)\s*;", blk).group(1),
            "extra": re.sub(r"(startFace|nFaces)\s+\d+\s*;", "", blk).strip(),
        })
    n_total = int(re.search(r"\n\s*(\d+)\s*\n\s*\(", txt).group(1))
    return patches, n_total


def write_list(p, cls, obj, vals, one_per_line=False):
    with open(p, "w", encoding="utf-8", newline="\n") as f:
        f.write(header(cls, obj))
        f.write(f"\n{len(vals)}\n(\n")
        for v in vals:
            f.write(f"{v}\n" if one_per_line else f"{v} ")
        if not one_per_line:
            f.write("\n")
        f.write(")\n\n// ************************************************************************* //\n")


def main():
    mesh = sys.argv[1]
    owner = read_ascii_list(os.path.join(mesh, "owner"))
    neigh = read_ascii_list(os.path.join(mesh, "neighbour"))
    faces = read_ascii_faces(os.path.join(mesh, "faces"))
    patches, n_faces_declared = read_boundary(os.path.join(mesh, "boundary"))

    n_int = len(neigh)
    assert len(owner) == len(faces), \
        f"owner {len(owner)} != faces {len(faces)}"
    assert n_int > 0
    for pt in patches:
        assert pt["startFace"] + pt["nFaces"] <= len(faces), \
            f"补丁 {pt['name']} 越界: {pt['startFace']}+{pt['nFaces']} > {len(faces)}"
    if n_int < 1:
        raise SystemExit("没有内部面，无法注入")

    # 复制最后一个内部面
    src = n_int - 1
    o, n = owner[src], neigh[src]
    assert o < n, "内部面 owner 应小于 neighbour"
    faces.insert(n_int, list(faces[src]))
    owner.insert(n_int, o)
    neigh.append(n)

    # boundary 全部 startFace +1（数量不变）
    for pt in patches:
        pt["startFace"] += 1

    # 写回
    write_list(os.path.join(mesh, "faces"), "faceCompactList", "faces", [], False)  # placeholder replaced below
    with open(os.path.join(mesh, "faces"), "w", encoding="utf-8", newline="\n") as f:
        f.write(header("faceList", "faces"))
        f.write(f"\n{len(faces)}\n(\n")
        for ids in faces:
            f.write(f"{len(ids)}(" + " ".join(str(i) for i in ids) + ")\n")
        f.write(")\n\n// ************************************************************************* //\n")
    write_list(os.path.join(mesh, "owner"), "labelList", "owner", owner)
    write_list(os.path.join(mesh, "neighbour"), "labelList", "neighbour", neigh)

    with open(os.path.join(mesh, "boundary"), "w", encoding="utf-8", newline="\n") as f:
        f.write(header("polyBoundaryMesh", "boundary"))
        f.write(f"\n{len(patches)}\n(\n")
        for pt in patches:
            f.write(f"    {pt['name']}\n    {{\n")
            f.write(f"        type            {pt['type']};\n")
            if pt["extra"]:
                f.write(f"        {pt['extra']}\n")
            f.write(f"        nFaces          {pt['nFaces']};\n")
            f.write(f"        startFace       {pt['startFace']};\n")
            f.write("    }\n")
        f.write(")\n\n// ************************************************************************* //\n")

    print(f"注入完成: 复制内部面 #{src} (owner={o}, neighbour={n}) 到 index {n_int}")
    print(f"  faces {len(faces)-1} -> {len(faces)}; neighbour {n_int} -> {n_int+1}; "
          f"boundary startFace 全部 +1")
    print(f"  验证: 重跑 merge_dup_cells.py 应报告单元对 ({o},{n}) 为 2 面")


if __name__ == "__main__":
    main()
