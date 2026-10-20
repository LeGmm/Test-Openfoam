#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
 把 chtMultiRegionFoam 的双区域界面接上「真正的共轭传热」
================================================================================
问题：fluent_cas_to_foam_cht.py 给两侧界面起了两个不同的名字
       （fluid_to_solid / solid_to_fluid），而 OF8 的共轭边界
       compressible::turbulentTemperatureCoupledBaffleMixed 是**按 patch 名配对**的，
       加上 0/T 里两侧都写成 zeroGradient → 结果 Q_interface ≡ 0，一分热都不传。

本脚本做三件事：
  1. 两侧界面补丁统一改名（默认 interface）。补丁类型必须保持 mappedWall
     （coupledBaffleMixed 内部会 cast 成 mappedPatchBase，写 wall 会 FATAL：
      "'...' not type 'mappedPatchBase'"），samplePatch 也同步成新名字；
  2. 0/<fluid>/T 与 0/<solid>/T 的界面补丁写 coupledBaffleMixed
     （OF8 里 kappa 由各区域 thermo 自动取，不需要 kappaMethod 条目）；
  3. 0/<fluid>/U 的界面补丁写 noSlip；system/controlDict 里 FO 的 patches 名同步。

用法：
  python tools/setup_cht_coupling.py caseCHT_qishi [interface]
================================================================================
"""
from __future__ import annotations

import os
import re
import sys

FLUID_OLD = "fluid_to_solid"
SOLID_OLD = "solid_to_fluid"
OTHER = {"fluid": "solid", "solid": "fluid"}

T_COUPLED = """    {name}
    {{
        type            compressible::turbulentTemperatureCoupledBaffleMixed;
        neighbourFieldName T;
        Tnbr            T;
        value           $internalField;
    }}"""

MAPPED_BLOCK = """    {name}
    {{
        type            mappedWall;
        inGroups        1 (wall);
        nFaces          {nfaces};
        startFace       {startface};
        sampleMode      nearestPatchFace;
        sampleRegion    {other};
        samplePatch     {name};
    }}"""

BLOCK = re.compile(r"(\s{4})(\S+)(\s*\n\s*\{)(.*?)(\n\s{4}\})", re.S)


def renames(text: str, new: str) -> str:
    return text.replace(FLUID_OLD, new).replace(SOLID_OLD, new)


def read_blocks(path: str):
    txt = open(path, encoding="utf-8").read()
    return txt, [(m.group(2), m.group(4), m.start(), m.end()) for m in BLOCK.finditer(txt)]


def fix_boundary(path: str, region: str, new: str) -> None:
    txt, blocks = read_blocks(path)
    txt = renames(txt, new)          # 先把 patch 名和 samplePatch 都改掉
    open(path, "w", encoding="utf-8").write(txt)

    txt, blocks = read_blocks(path)  # 再按新名字重写补丁块
    for name, body, s, e in blocks:
        if name != new:
            continue
        nf = re.search(r"nFaces\s+(\d+)", body).group(1)
        sf = re.search(r"startFace\s+(\d+)", body).group(1)
        newblock = MAPPED_BLOCK.format(name=new, nfaces=nf, startface=sf,
                                       other=OTHER[region])
        # 用新块替换（保留行首缩进）
        line_start = txt.rfind("\n", 0, s) + 1
        txt = txt[:line_start] + newblock + txt[e:]
        break
    open(path, "w", encoding="utf-8").write(txt)
    print(f"  boundary: {region} 侧界面 → {new} (mappedWall, 指向 {OTHER[region]})")


def replace_patch_block(path: str, new: str, template: str) -> None:
    txt = open(path, encoding="utf-8").read()
    pat = re.compile(r"(\s{4}" + re.escape(new) + r"\s*\n\s*\{)(.*?)(\n\s{4}\})", re.S)
    m = pat.search(txt)
    if not m:
        print(f"  !! {path} 里没找到补丁 {new}，跳过")
        return
    txt = txt[:m.start()] + "\n" + template + txt[m.end():]
    open(path, "w", encoding="utf-8").write(txt)
    print(f"  场文件已改: {path}")


def main() -> int:
    if len(sys.argv) < 2:
        print(__doc__)
        return 1
    case = sys.argv[1]
    new = sys.argv[2] if len(sys.argv) > 2 else "interface"

    for region, old in (("fluid", FLUID_OLD), ("solid", SOLID_OLD)):
        b = os.path.join(case, "constant", region, "polyMesh", "boundary")
        if not os.path.isfile(b):
            raise SystemExit(f"缺文件: {b}")
        if old not in open(b, encoding="utf-8").read():
            print(f"  {region}: 已是 {new}，跳过 boundary")
            continue
        fix_boundary(b, region, new)

    # 0/ 场文件：patch 名同步 + BC 定制
    for rel in (("0", "fluid", "T"), ("0", "solid", "T"), ("0", "fluid", "U")):
        p = os.path.join(case, *rel)
        if not os.path.isfile(p):
            continue
        txt = open(p, encoding="utf-8").read()
        if FLUID_OLD in txt or SOLID_OLD in txt:
            open(p, "w", encoding="utf-8").write(renames(txt, new))
            print(f"  场文件 patch 名同步: {p}")

    replace_patch_block(os.path.join(case, "0", "fluid", "T"), new,
                        T_COUPLED.format(name=new))
    replace_patch_block(os.path.join(case, "0", "solid", "T"), new,
                        T_COUPLED.format(name=new))
    up = os.path.join(case, "0", "fluid", "U")
    if os.path.isfile(up):
        replace_patch_block(up, new,
                            f"    {new}\n    {{\n        type            noSlip;\n    }}")

    cd = os.path.join(case, "system", "controlDict")
    if os.path.isfile(cd):
        txt = open(cd, encoding="utf-8").read()
        if FLUID_OLD in txt or SOLID_OLD in txt:
            open(cd, "w", encoding="utf-8").write(renames(txt, new))
            print(f"  controlDict 里 FO 的 patches 名同步")

    print(f"\n完成。两侧界面 patch 名统一为 {new}，T 用 coupledBaffleMixed。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
