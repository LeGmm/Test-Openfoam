#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
批量审计 8 个结构的 Fluent 数据源（FLTG*-Setup-Output.cas.h5）：
模型 / 重力 / 操作压力 / 材料 / 边界条件 是否正确且统一。
"""
import os
import re
import sys

import h5py

ROOT = "F:/B_Simulation/fangzhen-3"
CASES = {
    "v8_1":  "QiShi_v8_1_20261009/qishi_v8_1_files/dp0/FLTG/Fluent/FLTG-Setup-Output.cas.h5",
    "v8_8":  "Qishi_v8_8_20261012/qishi_v8_8_files/dp0/FLTG/Fluent/FLTG.2-Setup-Output.cas.h5",
    "v8_9":  "Qishi_v8_9_20261012/qishi_v8_9_files/dp0/FLTG/Fluent/FLTG.1-Setup-Output.cas.h5",
    "v8_10": "Qishi_v8_10_20261014/qishi-v8-10_files/dp0/FLTG/Fluent/FLTG-Setup-Output.cas.h5",
    "v8_11": "Qishi_v8_11_20261014/qishi-v8-11_files/dp0/FLTG/Fluent/FLTG-Setup-Output.cas.h5",
    "v8_12": "Qishi_v8_12_20261014/qishi-v8-12_files/dp0/FLTG/Fluent/FLTG-Setup-Output.cas.h5",
    "v8_13": "Qishi_v8_13_20261014/qishi-v8-13_files/dp0/FLTG/Fluent/FLTG-Setup-Output.cas.h5",
    "v8_14": "Qishi_v8_14_20261014/qishi-v8-14_files/dp0/FLTG/Fluent/FLTG-Setup-Output.cas.h5",
}


def load_text(path):
    parts = []
    with h5py.File(path) as f:
        for k in f["settings"]:
            parts.append(f[f"settings/{k}"][()].tobytes().decode("utf-8", "replace"))
        # 单元区拓扑
        zones = []
        zt = f.get("meshes/1/cells/zoneTopology")
        if zt is not None:
            for key in zt:
                d = zt[key][()]
                try:
                    zones.append(d.tobytes().decode("utf-8", "replace").strip())
                except Exception:
                    zones.append(str(d))
        ncells = None
        try:
            ncells = int(f["meshes/1/cells/coords"].shape[0]) if "meshes/1/cells/coords" in f else None
        except Exception:
            pass
    return "\n".join(parts), zones


def rx(pat, txt, cast=str, default="?"):
    m = re.search(pat, txt)
    if not m:
        return default
    return cast(m.group(1))


def flag(pat, txt):
    m = re.search(pat, txt)
    if not m:
        return "?"
    return "Y" if m.group(1) in ("#t", "yes") else "N"


def bc_block(txt, name):
    """取 (nn (zone type name ...)... 的线程块（到下一个 (nn ( 为止）。"""
    m = re.search(r"\(39 \(\d+ \w+ [\w\-.:]+ \d+\)\(", txt)
    blocks = re.findall(r"\(39 \((\d+) (\S+) ([\w\-.:]+) (\d+)\)\((.*?)(?=\(39 \(|\Z)", txt, re.S)
    for zid, ztype, zname, _ in blocks:
        if zname.lower() == name.lower():
            return ztype, zname
    return None, None


def zone_blocks(txt):
    return re.findall(r"\(39 \((\d+) (\S+) ([\w\-.:]+) (\d+)\)\((.*?)(?=\(39 \(|\Z)", txt, re.S)


def audit(name, rel):
    path = os.path.join(ROOT, rel)
    txt, zones = load_text(path)
    r = {"name": name}
    r["turb"] = ("laminar" if flag(r"\(rp-lam\? (#t|#f)", txt) == "Y"
                 else "kw" if flag(r"\(rp-kw\? (#t|#f)", txt) == "Y"
                 else "ke" if flag(r"\(rp-ke\? (#t|#f)", txt) == "Y" else "?")
    r["sst"] = flag(r"\(rp-trans-sst\? (#t|#f)", txt)
    r["energy"] = flag(r"\(energy\? (#t|#f)", txt)
    r["mp/lee"] = flag(r"\(mp/evaporation-condensation-model\? (#t|#f)", txt)
    r["species"] = "Y" if "(h2o" in txt or "species (h2o" in txt else "N"
    r["water-liquid"] = "Y" if "water-liquid" in txt else "N"
    r["humid-air"] = "Y" if "humid-air" in txt else "N"
    r["gravity_y"] = rx(r"\(gravity/y ([\d.eE+-]+)", txt, float)
    r["op_pres"] = rx(r"\(operating-pressure ([\d.eE+-]+)", txt, float)
    r["op_dens"] = rx(r"\(operating-density ([\d.eE+-]+)", txt, float)

    # 材料
    mats = re.findall(r"\(([\w\-.]+) (fluid|solid) \(chemical-formula", txt)
    r["materials"] = ",".join(f"{n}:{t}" for n, t in mats)

    # 边界条件（线程块）
    blocks = zone_blocks(txt)
    inlet = outlet = cool = None
    walls = []
    for zid, ztype, zname, _n, body in blocks:
        b = dict(ztype=ztype, name=zname)
        b["mf"] = rx(r"\(mass-flow \(constant\s*(?:\.\s*)?([\d.eE+-]+)", body, float, None)
        b["t0"] = rx(r"\(t0 \(constant\s*(?:\.\s*)?([\d.eE+-]+)", body, float, None)
        b["p"] = rx(r"\(p \(constant\s*(?:\.\s*)?([\d.eE+-]+)", body, float, None)
        b["tbc"] = rx(r"\(thermal-bc \.? (\d+)", body, int, None)
        b["t"] = rx(r"\(t \(constant\s*(?:\.\s*)?([\d.eE+-]+)", body, float, None)
        b["q"] = rx(r"\(q \(constant\s*(?:\.\s*)?([\d.eE+-]+)", body, float, None)
        n = zname.lower()
        if n == "inlet":
            inlet = b
        elif n == "outlet":
            outlet = b
        elif n == "wall_cool":
            cool = b
        elif ztype == "wall":
            walls.append(b)
    r["inlet"] = f"{inlet['ztype']}, mdot={inlet['mf']}, T0={inlet['t0']}" if inlet else "未找到"
    r["outlet"] = f"p={outlet['p']} (gauge)" if outlet else "未找到"
    r["wall_cool"] = (f"tbc={cool['tbc']}({'定温' if cool['tbc']==1 else '热流' if cool['tbc']==0 else cool['tbc']})"
                      f" T={cool['t']} q={cool['q']}") if cool else "未找到"
    other = [f"{w['name']}:tbc={w['tbc']},T={w['t']},q={w['q']}" for w in walls]
    r["other_walls"] = " | ".join(other) if other else "-"
    r["cell_zones"] = ";".join(zones)
    return r


def main():
    rows = [audit(n, rel) for n, rel in CASES.items()]
    keys = ["name", "turb", "sst", "energy", "mp/lee", "species", "water-liquid",
            "humid-air", "gravity_y", "op_pres", "op_dens", "inlet", "outlet",
            "wall_cool", "other_walls", "materials", "cell_zones"]
    for k in keys:
        print(f"\n== {k} ==")
        for r in rows:
            v = str(r[k])
            print(f"  {r['name']:6s} {v[:150]}")


if __name__ == "__main__":
    sys.exit(main())
