#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""把 Step1（chtMultiRegionFoam）收敛的体场热启动进 Step2（multiphaseEulerFoam）的 0/ 目录。

目的：消除 Step2 启动冲击——0/ 原是均匀 308.15K 的初场，第一步就要被 278K 冷壁
按 ArdenBuck 指数拉相变源项，薄片单元随即爆掉。用 Step1 已收敛的 T/U 作初场，
Step2 是从物理上合理的状态出发，启动刚性大幅下降。

前提（本工具不做校验，调用前自行确认）：两侧 polyMesh 必须同一单元序
（同一 cas.h5 同一转换器产物；用单元质心比对核验，见 2026-10-20 日志）。

用法:
    python tools/warmstart_step2.py <cht_case> <time> <step2_case>

映射（按 cell 一一对应）:
    <time>/fluid/T     -> 0/T.gas, 0/T.liquid 的 internalField
    <time>/fluid/U     -> 0/U.gas, 0/U.liquid 的 internalField
    <time>/fluid/p_rgh -> 0/p_rgh 的 internalField
附带（可 --no-yinit 关闭）: 按热启动 T 场对 H2O.gas/air.gas 做饱和度自洽初始化
    Y0 = min(Y_inlet, Y_sat(T))，air = 1 - Y0（ArdenBuck，与相变模型同式）。
    不做这一步时，冷 T 场 + 均匀 Y 会在第一步产生人为过饱和爆发（实测 t=4e-6 即爆）。
边界条件（含 map_interface_T 写入的 coldWall 非均匀值）原样保留。
"""
import math
import re
import sys
from pathlib import Path

IF_RE = re.compile(r"internalField\s+(.*?);", re.S)

# Step1（chtMultiRegionFoam，rhoConst）的 p_rgh 是零参考的"表压"；
# Step2（multiphaseEulerFoam，perfectGas）的 p_rgh 必须给"绝压"。
# 依据 OF8 源码 pU/pEqn.H:391 —— 每步无条件执行 p = max(p_rgh + rho*gh, pMin)，
# 若 p_rgh 按表压给 ~0，p 会被 pMin=10000 Pa 夹住 → 密度只剩 1/10 → 速度爆炸 SIGFPE。
P_OP = 101325.0


def y_sat(T: float, p_total: float = 101325.0) -> float:
    """ArdenBuck 饱和水汽质量分数（与 phaseProperties 的 saturationPressure 同式）。"""
    tc = T - 273.15
    psat_kpa = 0.61121 * math.exp((18.678 - tc / 234.5) * (tc / (257.14 + tc)))
    psat = psat_kpa * 1000.0
    if psat >= p_total * 0.99:
        return 0.99
    return 0.622 * psat / (p_total - psat)


def read_internal(path: Path):
    """返回 (kind, values, raw)：kind ∈ {uniform, nonuniform}。"""
    txt = path.read_text(encoding="utf-8", errors="replace")
    m = IF_RE.search(txt)
    if not m:
        raise ValueError(f"{path} 里没有 internalField")
    body = m.group(1).strip()
    if body.startswith("uniform"):
        tok = body.split()[1]
        val = float(tok) if not tok.startswith("(") else tok
        return "uniform", val, body
    mm = re.match(r"nonuniform\s+List<(scalar|vector)>\s*\n\s*(\d+)\s*\n\((.*)\)\s*$",
                  body, re.S)
    if not mm:
        raise ValueError(f"{path} 的 internalField 格式无法解析：{body[:80]}")
    n = int(mm.group(2))
    inner = mm.group(3)
    if mm.group(1) == "scalar":
        vals = [float(x) for x in inner.split()]
    else:
        vals = re.findall(r"\(([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\)", inner)
        vals = [tuple(float(x) for x in v) for v in vals]
    if len(vals) != n:
        raise ValueError(f"{path} 声明 {n} 个值但解析出 {len(vals)} 个")
    return "nonuniform", vals, body


def render_internal(kind: str, values, elem: str) -> str:
    if kind == "uniform":
        return f"internalField   uniform {values};"
    lines = [f"internalField   nonuniform List<{elem}> ", f"{len(values)}", "("]
    if elem == "scalar":
        lines += [repr(float(v)) if False else f"{v:.10g}" for v in values]
    else:
        lines += ["({0:.10g} {1:.10g} {2:.10g})".format(*v) for v in values]
    lines += [");"]
    return "\n".join(lines)


def set_internal(path: Path, kind: str, values, elem: str):
    txt = path.read_text(encoding="utf-8", errors="replace")
    new_if = render_internal(kind, values, elem)
    out, n = IF_RE.subn(lambda _m: new_if, txt, count=1)
    if n != 1:
        raise ValueError(f"{path} 替换 internalField 失败")
    path.write_text(out, encoding="utf-8", newline="\n")


def main():
    args = [a for a in sys.argv[1:] if a != "--no-yinit"]
    do_yinit = "--no-yinit" not in sys.argv
    if len(args) != 3:
        print(__doc__)
        sys.exit(2)
    cht, time, s2 = Path(args[0]), args[1], Path(args[2])
    src = cht / time / "fluid"
    if not src.is_dir():
        sys.exit(f"源目录不存在：{src}")

    plan = [
        ("T", ["T.gas", "T.liquid"], "scalar"),
        ("U", ["U.gas", "U.liquid"], "vector"),
        ("p_rgh", ["p_rgh"], "scalar"),
    ]
    t_vals = None
    for field, targets, elem in plan:
        sp = src / field
        if not sp.is_file():
            print(f"[SKIP] 源场缺失：{sp}")
            continue
        kind, vals, _ = read_internal(sp)
        if field == "T" and kind == "nonuniform":
            t_vals = vals
        if field == "p_rgh":
            # 表压 → 绝压（否则 p 被 pMin 夹住，密度塌陷，第一步即速度爆炸）
            vals = [v + P_OP for v in vals] if kind == "nonuniform" else vals + P_OP
            print(f"[SHIFT] p_rgh 零参考 → 绝压 (+{P_OP:g} Pa)")
        for t in targets:
            tp = s2 / "0" / t
            if not tp.is_file():
                print(f"[SKIP] 目标缺失：{tp}")
                continue
            set_internal(tp, kind, vals, elem)
            n = len(vals) if kind == "nonuniform" else 1
            print(f"[OK]   {field} -> {t}  ({kind}, {n} 值)")

    if do_yinit and t_vals is not None:
        h2o = s2 / "0" / "H2O.gas"
        air = s2 / "0" / "air.gas"
        if h2o.is_file() and air.is_file():
            kind, yvals, _ = read_internal(h2o)
            y_in = max(yvals) if kind == "nonuniform" else yvals
            y0 = [min(y_in, y_sat(t)) for t in t_vals]
            set_internal(h2o, "nonuniform", y0, "scalar")
            set_internal(air, "nonuniform", [1.0 - y for y in y0], "scalar")
            print(f"[OK]   组分饱和度自洽初始化: Y_in={y_in:.4g}, "
                  f"Y0 范围 {min(y0):.4g}~{max(y0):.4g}, air=1-Y0")
        else:
            print("[SKIP] 未找到 H2O.gas / air.gas，跳过组分初始化")
    print("完成。边界条件未动。")


if __name__ == "__main__":
    main()
