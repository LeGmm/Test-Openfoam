#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""按 tools/audit_cleanup.py 的机读结果执行删除（只删三项全过的项）。

一次性执行脚本，用完即归档：真正的判据在 audit_cleanup.py，本脚本不重新判定。
用法：
    python tools/audit_cleanup.py --json _cleanup_audit.json
    python tools/_do_cleanup.py _cleanup_audit.json
"""
import json
import os
import shutil
import sys


def main() -> int:
    src = sys.argv[1] if len(sys.argv) > 1 else "_cleanup_audit.json"
    rows = json.load(open(src, encoding="utf-8"))
    tgt = [r for r in rows if r["deletable"]]
    print(f"待删除 {len(tgt)} 项，合计 {sum(r['size_mb'] for r in tgt):.0f} MB")

    ok, fail = [], []
    for r in tgt:
        p = r["path"]
        try:
            if os.path.isdir(p):
                shutil.rmtree(p)
            elif os.path.exists(p):
                os.remove(p)
            ok.append(r)
        except Exception as e:      # noqa: BLE001
            fail.append((p, str(e)[:140]))
    print(f"成功 {len(ok)}  失败 {len(fail)}")
    for p, e in fail[:20]:
        print("  FAIL", p, e)

    by_top: dict[str, list] = {}
    for r in ok:
        top = r["path"].replace("\\", "/").split("/")[0]
        by_top.setdefault(top, [0, 0.0])
        by_top[top][0] += 1
        by_top[top][1] += r["size_mb"]

    lines = ["# 清理清单（2026-10-21）", "",
             "回滚基线：`git tag baseline-2026-10-21`（提交 a2a88e4）。",
             "判据工具：`tools/audit_cleanup.py`（只读三项核对：引用 / 性质 / 删除后果）。", "",
             f"候选 {len(rows)} 项；通过 {len(ok)} 项 / {sum(r['size_mb'] for r in ok):.0f} MB；"
             f"保留 {len(rows) - len(ok)} 项。", "",
             "## 已删除（三项全过）", ""]
    for k in sorted(by_top, key=lambda k: -by_top[k][1]):
        lines.append(f"- `{k}` —— {by_top[k][0]} 项 / {by_top[k][1]:.0f} MB")
    lines += ["", "## 保留（未通过核对，附依据）", ""]
    for r in rows:
        if r["deletable"]:
            continue
        why = []
        if r["check_ref"] != "PASS":
            why.append(f"①{r['check_ref']} {r['check_ref_note']}")
        if r["check_kind"] != "PASS":
            why.append(f"②{r['check_kind']} {r['check_kind_note']}")
        if r["check_impact"] != "PASS":
            why.append(f"③{r['check_impact']} {r['check_impact_note']}")
        lines.append(f"- `{r['path']}`（{r['size_mb']:.1f} MB）：" + "；".join(why))
    with open("清理记录_2026-10-21.md", "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print("清单 -> 清理记录_2026-10-21.md")
    return 1 if fail else 0


if __name__ == "__main__":
    sys.exit(main())
