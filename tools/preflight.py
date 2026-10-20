#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
 preflight —— 开跑前的门禁：把"会让人白跑一遍"的坑提前挡掉
================================================================================
 由来：2026-10-20 审计发现的三处严重问题，两处是"配置残留/默认值不一致"，
 症状都是**跑完才发现白跑**（startFrom=latestTime 直接空转；janaf 必 FATAL；
 界面补丁不同名导致 Q≡0 却看起来正常）。这些都可以在开跑前 2 秒内查出来。

 检查项（FAIL 会让退出码变 1，WARN 只提示）：
   1 python    tools/*.py 语法编译（含 SyntaxWarning）
   2 deps      numpy/h5py（必需）、scipy/matplotlib（画图可选）
   3 openfoam  blueCFD 的 setvars_OF8.bat（of8.py 能否工作）
   4 layout    单区域 / 多区域（regionProperties），constant 与 0 目录是否配套
   5 polymesh  点数/面数/owner/neighbour 计数自洽，补丁 startFace 连续，单元编号无洞
   6 topology  一对单元多张面：真重复→FAIL（可修）；铰链→WARN（按既定决策保留）
   7 patches   0/<region>/ 各场的补丁名与 polyMesh/boundary 一致
   8 controld  startFrom 空转、endTime、adjustTimeStep/maxCo
   9 cht       多区域：界面两侧同名 + mappedWall + 0/T 用 coupledBaffleMixed（否则 FATAL 或 Q≡0）
  10 step2     单区域多相：janaf + purePhaseModel 组合必 FATAL，提示跑 align_step2_thermo.py
  11 pressure  单区域多相：p_rgh 必须是绝压（表压口径会让 p 被 pMin 夹住 → 第一步 SIGFPE）

 用法：
   python tools/preflight.py caseCHT_qishi
   python tools/preflight.py caseCHT_qishi --checkmesh      # 追加 checkMesh（慢，约 10 s）
   python tools/preflight.py case_qishi --json              # 机器可读输出
================================================================================
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import re
import sys
import warnings

TOOLS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, TOOLS)

PASS, WARN, FAIL, SKIP = "PASS", "WARN", "FAIL", "SKIP"
_res: list[tuple[str, str, str]] = []


def add(status: str, name: str, msg: str) -> None:
    _res.append((status, name, msg))


def _read(p: str) -> str:
    return open(p, encoding="utf-8", errors="replace").read()


def _kv(text: str, key: str):
    m = re.search(rf"^\s*{key}\s+([^;]+);", text, re.M)
    return m.group(1).strip() if m else None


def time_dirs(case: str) -> list[str]:
    out = []
    for d in os.listdir(case):
        if not os.path.isdir(os.path.join(case, d)) or d == "0":
            continue
        try:
            float(d)
        except ValueError:
            continue
        out.append(d)
    return out


# ------------------------------------------------------------------ 1 语法
def check_tools() -> None:
    bad = []
    for f in sorted(glob.glob(os.path.join(TOOLS, "*.py"))):
        with warnings.catch_warnings(record=True) as ws:
            warnings.simplefilter("always")
            try:
                compile(_read(f), f, "exec")
            except SyntaxError as e:
                bad.append(f"{os.path.basename(f)}: {e.msg} (line {e.lineno})")
                continue
        for w in ws:
            if issubclass(w.category, SyntaxWarning):
                bad.append(f"{os.path.basename(f)}: {w.message}")
    if bad:
        add(FAIL, "python", "语法/告警问题：\n      " + "\n      ".join(bad))
    else:
        n = len(glob.glob(os.path.join(TOOLS, "*.py")))
        add(PASS, "python", f"{n} 个脚本编译通过，无 SyntaxWarning")


# ------------------------------------------------------------------ 2 依赖
def check_deps() -> None:
    missing = []
    for mod in ("numpy", "h5py"):
        try:
            __import__(mod)
        except ImportError:
            missing.append(mod)
    if missing:
        add(FAIL, "deps", f"缺少必需模块 {missing}（h5py 只有 cas.h5 转换需要）")
        return
    opt = []
    for mod in ("scipy", "matplotlib"):
        try:
            __import__(mod)
        except ImportError:
            opt.append(mod)
    add(PASS, "deps", "numpy / h5py 就绪" + (f"；缺 {opt}（仅出图需要）" if opt else ""))


# ------------------------------------------------------------- 3 OpenFOAM
def check_of() -> None:
    try:
        import of8
    except Exception as e:                                   # 理论上不会发生
        add(WARN, "openfoam", f"无法导入 of8.py：{e}")
        return
    if os.path.isfile(of8.SETVARS):
        add(PASS, "openfoam", f"blueCFD 环境就绪（{of8.BLUECFD}）")
    else:
        add(FAIL, "openfoam", f"找不到 {of8.SETVARS}；用环境变量 OF8_HOME 指向安装根目录")


# ---------------------------------------------------------------- 4 结构
def detect_regions(case: str) -> list[str]:
    const = os.path.join(case, "constant")
    if not os.path.isdir(const):
        return []
    rp = os.path.join(const, "regionProperties")
    if os.path.isfile(rp):
        txt = _read(rp)
        m = re.search(r"regions\s*\(([^;]*)\);", txt, re.S)
        if m:
            names = []
            for kind, body in re.findall(r"(fluid|solid)\s*\(([^)]*)\)", m.group(1)):
                names += body.split()
            names = [n for n in names if os.path.isdir(os.path.join(const, n))]
            if names:
                return names
    if os.path.isdir(os.path.join(const, "polyMesh")):
        return [""]                                   # 单区域，region 名为空
    return []


def check_layout(case: str, regions: list[str]) -> None:
    if not regions:
        add(FAIL, "layout", f"{case} 里没有 constant/polyMesh，也不是有效的多区域算例")
        return
    kind = "单区域" if regions == [""] else f"多区域（{', '.join(regions)}）"
    problems = []
    for r in regions:
        pm = os.path.join(case, "constant", r, "polyMesh") if r else os.path.join(case, "constant", "polyMesh")
        f0 = os.path.join(case, "0", r) if r else os.path.join(case, "0")
        for req in ("points", "faces", "owner", "neighbour", "boundary"):
            if not os.path.isfile(os.path.join(pm, req)):
                problems.append(f"缺 {os.path.relpath(os.path.join(pm, req), case)}")
        if not os.path.isdir(f0):
            problems.append(f"缺 0/{r or ''} 目录")
    if problems:
        add(FAIL, "layout", f"{kind}；" + "；".join(problems[:4]))
    else:
        add(PASS, "layout", kind + "，polyMesh 与 0/ 目录齐备")


# --------------------------------------------------------------- 5 polyMesh
def mesh_dir(case: str, region: str) -> str:
    return os.path.join(case, "constant", region, "polyMesh") if region \
        else os.path.join(case, "constant", "polyMesh")


def check_polymesh(case: str, regions: list[str]) -> None:
    from ofmesh import read_boundary, read_labels
    for r in regions:
        md = mesh_dir(case, r)
        tag = r or "case"
        own = read_labels(os.path.join(md, "owner"))
        nei = read_labels(os.path.join(md, "neighbour"))
        nf = len(own)
        n_int = len(nei)
        patches = read_boundary(os.path.join(md, "boundary"))
        errs = []
        n_bound = sum(p["nFaces"] for p in patches.values())
        if n_int + n_bound != nf:
            errs.append(f"内部面 {n_int} + 边界面 {n_bound} ≠ 总面数 {nf}")
        cursor = n_int
        for name, p in patches.items():
            if p["startFace"] != cursor:
                errs.append(f"补丁 {name} startFace={p['startFace']} 应为 {cursor}")
            cursor += p["nFaces"]
        used = set(own.tolist()) | set(nei.tolist())
        if used and max(used) + 1 != len(used):
            errs.append(f"单元编号有洞：max={max(used)} 但只出现 {len(used)} 个")
        if errs:
            add(FAIL, f"polymesh[{tag}]", "；".join(errs))
        else:
            add(PASS, f"polymesh[{tag}]", f"{nf} 面（内部 {n_int}）/ {len(used)} 单元 / "
                                          f"{len(patches)} 个补丁，计数与 startFace 自洽")


# --------------------------------------------------------------- 6 拓扑
def check_topology(case: str, regions: list[str]) -> None:
    from ofmesh import classify_multi_face, duplicate_pairs, read_labels, read_mesh
    for r in regions:
        md = mesh_dir(case, r)
        tag = r or "case"
        pts, faces, own, nei, _ = read_mesh(md)
        groups = duplicate_pairs(own, nei)
        if not groups:
            add(PASS, f"topology[{tag}]", "没有多面连接的单元对")
            continue
        dup, hinge = [], []
        for key, ids in groups.items():
            info = classify_multi_face(pts, faces, ids)
            (dup if info["verdict"] == "duplicate" else hinge).append((key, info))
        if dup:
            msg = "；".join(f"{k} 面心距 {i['dmax']*1e6:.1f}µm（真重复）" for k, i in dup[:3])
            add(FAIL, f"topology[{tag}]",
                f"{len(dup)} 对真重复面：{msg}\n"
                f"      → python tools/mesh_topology_report.py {md} --fix")
        if hinge:
            msg = "；".join(f"{k} 面心距 {i['dmax']*1e6:.1f}µm 共 {i['shared_nodes']} 节点"
                           for k, i in hinge[:3])
            add(WARN, f"topology[{tag}]",
                f"{len(hinge)} 对铰链（合法拓扑，按既定决策保留）：{msg}")


# --------------------------------------------------------------- 7 补丁名
def patch_names_in_field(path: str):
    """返回 (补丁名集合, 是否用了引号正则键)。正则键（如 ".*"）视为覆盖全部补丁。"""
    txt = _read(path)
    if "boundaryField" not in txt:
        return None
    bf = txt.split("boundaryField", 1)[1]
    names = re.findall(r"^\s*(\S+)\s*\n\s*\{", bf, re.M)
    has_regex = any(n.startswith('"') for n in names)
    return {n for n in names if not n.startswith(("#", '"'))}, has_regex


def check_patch_names(case: str, regions: list[str]) -> None:
    from ofmesh import read_boundary
    bad = []
    for r in regions:
        pset = set(read_boundary(os.path.join(mesh_dir(case, r), "boundary")))
        d = os.path.join(case, "0", r) if r else os.path.join(case, "0")
        for f in sorted(glob.glob(os.path.join(d, "*"))):
            if os.path.isdir(f):
                continue
            got = patch_names_in_field(f)
            if got is None:
                continue
            ns, has_regex = got
            if has_regex:
                continue
            missing = pset - ns
            if missing:
                bad.append(f"{os.path.relpath(f, case)} 缺补丁 {sorted(missing)}")
    if bad:
        add(FAIL, "patches", "；".join(bad[:4]))
    else:
        add(PASS, "patches", "0/ 各场的补丁名与 polyMesh 一致（正则补丁 \".*\" 视为全覆盖）")


# ------------------------------------------------------------ 8 controlDict
def check_control_dict(case: str) -> None:
    cd = os.path.join(case, "system", "controlDict")
    if not os.path.isfile(cd):
        add(FAIL, "controld", f"缺 {cd}")
        return
    txt = _read(cd)
    msgs, status = [], PASS
    start_from = _kv(txt, "startFrom")
    end = _kv(txt, "endTime")
    dt = _kv(txt, "deltaT")
    adj = _kv(txt, "adjustTimeStep")
    co = _kv(txt, "maxCo")

    if start_from == "latestTime":
        tds = time_dirs(case)
        if tds:
            latest = max(tds, key=float)
            status = WARN
            msgs.append(f"startFrom=latestTime 且已有时间目录（最新 {latest}）→ 重跑会从这里续算，"
                        f"可能直接结束；确认无误后改 startFrom startTime")
    if end is None or dt is None:
        status = FAIL
        msgs.append("缺 endTime 或 deltaT")
    else:
        try:
            if float(end) <= 0:
                status = FAIL
                msgs.append(f"endTime={end} 不会推进")
        except ValueError:
            msgs.append(f"endTime 非数值（{end}），按迭代次数解释")
    if (adj or "no").lower().startswith("y") and co is None:
        status = WARN if status != FAIL else status
        msgs.append("adjustTimeStep=yes 但没有 maxCo")
    add(status, "controld", "；".join(msgs) if msgs else
        f"startFrom={start_from} endTime={end} deltaT={dt} adjust={adj} maxCo={co}")


# -------------------------------------------------------------- 9 cht 界面
def check_cht_interface(case: str, regions: list[str]) -> None:
    if len(regions) < 2:
        add(SKIP, "cht", "单区域算例，跳过界面耦合检查")
        return
    from ofmesh import read_boundary
    maps = {r: {n: p for n, p in read_boundary(
        os.path.join(mesh_dir(case, r), "boundary")).items() if p["type"] == "mappedWall"}
        for r in regions}
    if not any(maps.values()):
        add(FAIL, "cht", "没有任何 mappedWall 补丁，两区域不会发生共轭传热")
        return
    problems, info = [], []
    for r, ps in maps.items():
        others = [x for x in regions if x != r]
        for name, p in ps.items():
            tgt = dict(p.get("extra", [])).get("sampleRegion", others[0])
            if name not in maps.get(tgt, {}):
                problems.append(f"{r}:{name} 在 {tgt} 侧没有同名补丁（OF8 按名配对，会配不上）")
                continue
            tfile = os.path.join(case, "0", r, "T")
            if os.path.isfile(tfile):
                body = _read(tfile)
                blk = re.search(r"^\s{4}" + re.escape(name) + r"\s*\n\s*\{(.*?)\n\s{4}\}", body, re.S | re.M)
                if not blk or "CoupledBaffleMixed" not in blk.group(1):
                    problems.append(f"0/{r}/T 的 {name} 不是 coupledBaffleMixed（会 Q≡0 不传热）")
                else:
                    info.append(f"{r}:{name}")
    if problems:
        add(FAIL, "cht", "；".join(problems[:4]) + "\n"
                         "      → python tools/setup_cht_coupling.py " + case)
    else:
        add(PASS, "cht", "界面双侧同名、mappedWall、T=coupledBaffleMixed：" + ", ".join(info))


# --------------------------------------------------------- 11 p_rgh 压力口径
def check_p_ref(case: str, regions: list[str]) -> None:
    """multiphaseEulerFoam 的 p_rgh 必须是【绝压】。

    OF8 源码 pU/pEqn.H:391 每步无条件执行 p = max(p_rgh + rho*gh, pMin)。
    气相是 perfectGas，密度由 p 算。若 p_rgh 按表压给 ~0，p 会被 pMin=10000 Pa
    夹住 → 密度只剩常压的 1/10 → 出口与内部差一个大气压 → 速度一步冲到 3000 m/s
    → FOAM_SIGFPE。这条门禁就是为了不再"跑完才发现白跑"。
    """
    if regions != [""]:
        add(SKIP, "pressure", "多区域算例（Step1）p_rgh 为表压口径，跳过")
        return
    pr = os.path.join(case, "0", "p_rgh")
    pp = os.path.join(case, "constant", "phaseProperties")
    if not os.path.isfile(pr) or not os.path.isfile(pp):
        add(SKIP, "pressure", "无 0/p_rgh 或 phaseProperties，跳过")
        return
    txt = _read(pr)
    m = re.search(r"internalField\s+(?:uniform\s+([-\d.eE+]+)|nonuniform\s+List<scalar>\s*\d+\s*\(\s*([-\d.eE+]+))",
                  txt)
    if not m:
        add(WARN, "pressure", "0/p_rgh internalField 解析不了，请人工确认是绝压")
        return
    val = float(m.group(1) or m.group(2))
    pm = re.search(r"pMin\s+([-\d.eE+]+)", _read(pp))
    pmin = float(pm.group(1)) if pm else 10000.0
    msg = f"0/p_rgh internalField ≈ {val:.0f} Pa（pMin={pmin:.0f}）"
    if val < 10 * pmin:
        add(FAIL, "pressure", msg + "\n"
                         "      → 这是表压口径。multiphaseEulerFoam 要求绝压，否则 p 被 pMin 夹住、"
                         "密度塌 10 倍、第一步即 SIGFPE。\n"
                         f"      → 改成 uniform 101325，出口用 prghPressure p=101325；"
                         f"或重跑 python tools/generate_case.py {case}")
        return
    out = re.search(r"outlet\s*\{(.*?)\}", txt, re.S)
    if out and re.search(r"type\s+fixedValue", out.group(1)) and \
            re.search(r"value\s+uniform\s+0\s*;", out.group(1)):
        add(FAIL, "pressure", "出口 p_rgh 用 fixedValue 0（表压口径）\n"
                              "      → 应改 prghPressure p=101325")
        return
    add(PASS, "pressure", msg + "，绝压口径正确")


# ------------------------------------------------------------- 10 Step2 热物性
def check_step2_thermo(case: str, regions: list[str]) -> None:
    if regions != [""]:
        add(SKIP, "step2", "多区域算例，跳过相变热物性检查")
        return
    pp = os.path.join(case, "constant", "phaseProperties")
    gas = os.path.join(case, "constant", "thermophysicalProperties.gas")
    if not os.path.isfile(pp):
        add(SKIP, "step2", "没有 phaseProperties（单相算例），跳过")
        return
    ptxt = _read(pp)
    gtxt = _read(gas) if os.path.isfile(gas) else ""
    if "janaf" in gtxt and "purePhaseModel" in ptxt:
        add(FAIL, "step2", "气相 janaf + 液相 purePhaseModel：相变启动必 FATAL\n"
                           "      （thermoI.H: Maximum number of iterations exceeded: 100）\n"
                           f"      → python tools/align_step2_thermo.py {case}")
        return
    phases = [os.path.basename(p).split(".")[-1]
              for p in glob.glob(os.path.join(case, "constant", "thermophysicalProperties.*"))]
    miss = [p for p in phases if not os.path.isfile(os.path.join(case, "0", f"T.{p}"))]
    if miss:
        add(WARN, "step2", f"缺 0/T.{miss}（对应相场的初始温度）")
    else:
        add(PASS, "step2", f"热物性已对齐（{', '.join(phases)}），相场 0/ 文件齐备")


# ------------------------------------------------------------------ checkMesh
def check_checkmesh(case: str, regions: list[str]) -> None:
    import subprocess
    for r in regions:
        cmd = [sys.executable, os.path.join(TOOLS, "of8.py"), case,
               f"checkMesh -region {r} -constant" if r else "checkMesh -constant"]
        p = subprocess.run(cmd, capture_output=True, text=True, errors="replace")
        out = (p.stdout or "") + (p.stderr or "")
        if "Mesh OK" in out:
            add(PASS, f"checkMesh[{r or 'case'}]", "Mesh OK")
        elif "Failed" in out:
            fails = [l.strip() for l in out.splitlines() if "Failed" in l]
            add(FAIL, f"checkMesh[{r or 'case'}]", "；".join(fails[:2]))
        else:
            add(WARN, f"checkMesh[{r or 'case'}]", "无法判定结果，请手动看输出")


def main() -> int:
    ap = argparse.ArgumentParser(description="开跑前门禁")
    ap.add_argument("case")
    ap.add_argument("--checkmesh", action="store_true", help="追加 checkMesh（慢）")
    ap.add_argument("--json", action="store_true", help="机器可读输出")
    a = ap.parse_args()
    case = a.case.rstrip("\\/")
    if not os.path.isdir(case):
        print(f"算例目录不存在：{case}")
        return 2

    check_tools()
    check_deps()
    check_of()
    regions = detect_regions(case)
    check_layout(case, regions)
    if regions:
        try:
            check_polymesh(case, regions)
            check_topology(case, regions)
            check_patch_names(case, regions)
        except Exception as e:                              # 网格读不动时不要连带其它检查一起崩
            add(FAIL, "mesh", f"读取 polyMesh 失败：{type(e).__name__}: {e}")
        check_control_dict(case)
        check_cht_interface(case, regions)
        check_step2_thermo(case, regions)
        check_p_ref(case, regions)
        if a.checkmesh:
            check_checkmesh(case, regions)

    n_fail = sum(1 for s, _, _ in _res if s == FAIL)
    n_warn = sum(1 for s, _, _ in _res if s == WARN)
    if a.json:
        print(json.dumps({"case": case, "fail": n_fail, "warn": n_warn,
                          "results": [{"status": s, "name": n, "msg": m} for s, n, m in _res]},
                         ensure_ascii=False, indent=2))
    else:
        print(f"\n=== preflight: {case} ===\n")
        for s, n, m in _res:
            mark = {PASS: "[PASS]", WARN: "[WARN]", FAIL: "[FAIL]", SKIP: "[SKIP]"}[s]
            print(f"{mark} {n:<18} {m}")
        verdict = "可以开跑" if n_fail == 0 else f"有 {n_fail} 项必须先处理"
        print(f"\n结论：{verdict}" + (f"（{n_warn} 项提示）" if n_warn else ""))
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
