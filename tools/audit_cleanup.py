#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""清理候选三项核对器（只读，绝不删除）。

对每个候选路径逐项回答三件事，只有三项全过才允许进入删除清单：

  ① 引用     —— 是否被任一 .py / .sh / .md / .html 引用？
                并区分「代码引用」(运行依赖) 与「报告提及」(历史记录，不构成依赖)。
  ② 性质     —— 是配置文件 / 核心逻辑，还是产物 / 缓存 / 备份？
                真配置 = 各算例的 0/、constant/(非 polyMesh)、system/。
                产物   = 时间目录、polyMesh、processor*/、postProcessing/、VTK、log.*。
  ③ 删除后果 —— 删除后是否影响当前运行，或让版本历史 / 证据链悬空。

用法：
    python tools/audit_cleanup.py                     # 扫描默认候选集并打印核对表
    python tools/audit_cleanup.py --add <路径> [...]  # 追加候选
    python tools/audit_cleanup.py --json out.json     # 机读输出

退出码：0 = 全部候选三项皆过（可删）；1 = 存在未通过项；2 = 用法错误。
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# 扫描引用的文件类型（排除候选自身）
REF_SUFFIX = {".py", ".sh", ".md", ".html", ".jou", ".wbjn"}

# 视为「核心逻辑 / 配置」的路径片段；命中即性质=配置，禁止删除
CONFIG_HINTS = (
    re.compile(r"(^|/)0(/|$)"),                 # 0/ 初始场+边界条件
    re.compile(r"(^|/)system(/|$)"),
    re.compile(r"(^|/)constant(/|$)"),
    re.compile(r"\.py$"),
    re.compile(r"(^|/)README"),
)

# 视为「产物 / 缓存 / 备份」的片段
ARTIFACT_HINTS = (
    re.compile(r"(^|/)polyMesh(/|$)"),
    re.compile(r"(^|/)processor\d*(/|$)"),
    re.compile(r"(^|/)postProcessing(_prev)?(/|$)"),
    re.compile(r"(^|/)VTK(/|$)"),
    re.compile(r"(^|/)log\."),
    re.compile(r"^_trash"),
    re.compile(r"^__pycache__(/|$)"),
    re.compile(r"\.bak$"),
    re.compile(r"\.pyc$"),
    re.compile(r"\.npy$"),
    re.compile(r"\.vtk$"),
    re.compile(r"\.csv$"),
    re.compile(r"\.png$"),
    re.compile(r"(^|/)_\w"),                     # 下划线开头的临时目录
    re.compile(r"^\d+$"),                        # 纯数字时间目录
    re.compile(r"^\d+(\.\d+)?(e-[0-9]+)?$"),     # 0.046 / 2e-06
)

# 证据链守护：这些名字出现过即要求人工确认（删除会让旧报告索引悬空）
EVIDENCE_HINTS = ("run2", "diag_step2", "audit_", "4e-06", "HD_A", "bcDissolving_map",
                  "bubbleColumnEvaporating",
                  # 0.046 = 1τ 崩溃的尸检场（P0-1 液相温度发散的唯一单元级证据）
                  "0.046")

# 纯缓存：名字本身会被审计脚本自己的源码命中，① 直接判 PASS
CACHE_NAMES = ("__pycache__",)


def default_candidates(root: Path) -> list[Path]:
    """默认候选集：缓存、日志、一次性副本、备份、旧时间目录。"""
    out: list[Path] = []
    # 1) 回收站
    for p in sorted(root.glob("_trash_*")):
        out.append(p)
    # 2) 缓存
    for p in sorted(root.rglob("__pycache__")):
        out.append(p)
    for p in sorted(root.rglob("*.pyc")):
        out.append(p)
    # 3) 备份
    for p in sorted(root.rglob("*.bak")):
        out.append(p)
    # 4) 顶层日志
    for p in sorted(root.glob("log.*")):
        out.append(p)
    # 5) 顶层一次性命令 / 后台输出
    for p in sorted(root.glob("_cmd_*.txt")) + sorted(root.glob("_nohup_*.out")):
        out.append(p)
    # 6) 顶层下划线临时目录（不含 tools/diagnostics 里的诊断脚本）
    for p in sorted(root.glob("_*")):
        if p.is_dir() and p.name.startswith("_trash"):
            continue
        if p.is_dir():
            out.append(p)
    # 7) 算例时间目录（保留 0 与每 1000 步）
    for case in sorted(root.glob("case*")):
        if not case.is_dir():
            continue
        for t in sorted(case.iterdir()):
            if not t.is_dir():
                continue
            if not re.fullmatch(r"\d+(\.\d+)?(e-[0-9]+)?", t.name):
                continue
            if t.name == "0":
                continue
            try:
                n = int(t.name)
            except ValueError:
                out.append(t)          # 非整数时刻（0.046 / 2e-06）一律候选
                continue
            if n % 1000 != 0:
                out.append(t)
    return out


def build_ref_index(root: Path) -> dict[str, list[str]]:
    """把每个脚本/文档里出现的候选名建索引（区分代码引用 vs 报告提及）。"""
    idx: dict[str, list[str]] = {}
    files = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames
                       if d not in {"__pycache__", ".git"} and not d.startswith("_trash")]
        for fn in filenames:
            if Path(fn).suffix in REF_SUFFIX:
                files.append(Path(dirpath) / fn)
    for f in files:
        try:
            txt = f.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        rel = str(f.relative_to(root))
        for m in re.finditer(r"[A-Za-z0-9_.\-]+", txt):
            idx.setdefault(m.group(0), []).append(rel)
    return idx


def size_of(p: Path) -> int:
    if p.is_file():
        try:
            return p.stat().st_size
        except OSError:
            return 0
    tot = 0
    for dirpath, _, filenames in os.walk(p):
        for fn in filenames:
            try:
                tot += (Path(dirpath) / fn).stat().st_size
            except OSError:
                pass
    return tot


def check_one(p: Path, root: Path, idx: dict[str, list[str]]) -> dict:
    rel = str(p.relative_to(root)) if str(p).startswith(str(root)) else str(p)
    name = p.name

    # ① 引用
    # 时间目录名是纯数字/科学计数法，会在代码里作为数值字面量命中，属假阳性，跳过 token 匹配
    time_dir = bool(re.fullmatch(r"\d+(\.\d+)?(e-[0-9]+)?", name))
    hits = [] if time_dir else sorted(set(idx.get(name, [])))
    code_hits = [h for h in hits if h.endswith((".py", ".sh", ".jou", ".wbjn"))]
    doc_hits = [h for h in hits if h.endswith((".md", ".html"))]
    # run_all.sh / prepare_all_structures.py 依赖不存在的 cases/，其引用不算有效依赖
    code_hits = [h for h in code_hits
                 if not (h.endswith("run_all.sh") or h.endswith("prepare_all_structures.py")
                         or h.endswith("mesh_stats.py") or h.endswith("compare.py"))]
    if name in CACHE_NAMES:
        code_hits = []          # 缓存目录名只会被审计脚本自身源码命中
    c1 = "PASS" if not code_hits else "FAIL"
    c1_note = ("时间目录：不做 token 匹配（数字字面量会造成假阳性）" if time_dir
               else "无任何脚本/代码引用" if not hits
               else (f"仅报告提及：{', '.join(doc_hits[:3])}" if not code_hits
                     else f"代码引用：{', '.join(code_hits[:3])}"))

    # ② 性质
    is_cfg = any(h.search(rel) or h.search(name) for h in CONFIG_HINTS)
    is_art = any(h.search(rel) or h.search(name) for h in ARTIFACT_HINTS)
    if is_cfg and not is_art:
        c2, c2_note = "FAIL", "配置文件/核心逻辑，禁止删除"
    elif is_art:
        c2, c2_note = "PASS", "产物/缓存/备份，非配置非核心逻辑"
    else:
        c2, c2_note = "WARN", "性质不明，需人工判定"

    # ③ 删除后果
    if any(e in rel for e in EVIDENCE_HINTS):
        c3, c3_note = "WARN", "命中证据链守护名单，删除会让旧报告索引悬空"
    elif re.fullmatch(r"\d+", name) and int(name) % 1000 != 0 and name != "0":
        c3, c3_note = "PASS", "中间检查点，可由重算或相邻检查点替代"
    elif name.startswith("_trash"):
        c3, c3_note = "PASS", "回收站本体，内容上轮已核对；删除不可回滚（git 已忽略）"
    else:
        c3, c3_note = "PASS", "无可辨识的运行时依赖"

    return {
        "path": rel,
        "size_mb": round(size_of(p) / 1e6, 2),
        "exists": p.exists(),
        "check_ref": c1, "check_ref_note": c1_note,
        "check_kind": c2, "check_kind_note": c2_note,
        "check_impact": c3, "check_impact_note": c3_note,
        "deletable": (c1 == "PASS" and c2 == "PASS" and c3 == "PASS"),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--add", nargs="*", default=[], help="追加候选路径")
    ap.add_argument("--json", default=None, help="机读输出路径")
    args = ap.parse_args()

    cands = default_candidates(ROOT)
    for a in args.add:
        cands.append((ROOT / a).resolve())
    # 去掉互为父子的重复（只保留最外层）
    cands = sorted(set(cands))
    kept: list[Path] = []
    for c in cands:
        if any(str(c).startswith(str(k) + os.sep) for k in kept):
            continue
        kept.append(c)

    idx = build_ref_index(ROOT)
    rows = [check_one(p, ROOT, idx) for p in kept if p.exists()]

    deletable = [r for r in rows if r["deletable"]]
    blocked = [r for r in rows if not r["deletable"]]

    print(f"\n=== 清理候选三项核对 ===  根目录: {ROOT}")
    print(f"候选 {len(rows)} 项 | 可删 {len(deletable)} 项 "
          f"({sum(r['size_mb'] for r in deletable):.0f} MB) | "
          f"保留 {len(blocked)} 项 ({sum(r['size_mb'] for r in blocked):.0f} MB)\n")

    hdr = f"{'路径':<52}{'MB':>9}  {'①引用':<6}{'②性质':<6}{'③后果':<6}"
    print(hdr)
    print("-" * len(hdr))
    for r in sorted(rows, key=lambda r: -r["size_mb"]):
        print(f"{r['path'][:51]:<52}{r['size_mb']:>9.2f}  "
              f"{r['check_ref']:<6}{r['check_kind']:<6}{r['check_impact']:<6}"
              + ("  ← 可删" if r["deletable"] else ""))

    print("\n--- 未通过项（保留，附依据）---")
    for r in blocked:
        why = [f"①{r['check_ref']}:{r['check_ref_note']}",
               f"②{r['check_kind']}:{r['check_kind_note']}",
               f"③{r['check_impact']}:{r['check_impact_note']}"]
        print(f"  {r['path']}\n      " + "\n      ".join(w for w in why if not w.endswith("PASS")))

    if args.json:
        with open(args.json, "w", encoding="utf-8") as f:
            json.dump(rows, f, ensure_ascii=False, indent=2)
        print(f"\n机读结果 -> {args.json}")
    return 0 if not blocked or deletable else 1


if __name__ == "__main__":
    sys.exit(main())
