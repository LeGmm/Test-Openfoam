#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
批量：把 8 个结构版本的 Fluent 网格 + 算例一次性全部准备好。

 对每个结构：
   1) fluent_cas_to_foam.py  .cas.h5 -> constant/polyMesh
   2) generate_case.py       生成 0/ constant/ system/ 全套字典（OpenFOAM 8 语法）

 输出到 OpenFOAM/cases/<结构名>/
================================================================================
"""
from __future__ import annotations

import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)                      # OpenFOAM/
CASES = os.path.join(ROOT, "cases")
PY = r"C:/Users/LELE/.workbuddy/binaries/python/envs/default/Scripts/python.exe"

# 结构名 -> .cas.h5 相对路径（相对工程根 f:/B_Simulation/fangzhen-3）
STRUCTURES = {
    "v8_1":  "QiShi_v8_1_20261009/qishi_v8_1_files/dp0/FLTG/Fluent/FLTG-Setup-Output.cas.h5",
    "v8_8":  "Qishi_v8_8_20261012/qishi_v8_8_files/dp0/FLTG/Fluent/FLTG.2-Setup-Output.cas.h5",
    "v8_9":  "Qishi_v8_9_20261012/qishi_v8_9_files/dp0/FLTG/Fluent/FLTG.1-Setup-Output.cas.h5",
    "v8_10": "Qishi_v8_10_20261014/qishi-v8-10_files/dp0/FLTG/Fluent/FLTG-Setup-Output.cas.h5",
    "v8_11": "Qishi_v8_11_20261014/qishi-v8-11_files/dp0/FLTG/Fluent/FLTG-Setup-Output.cas.h5",
    "v8_12": "Qishi_v8_12_20261014/qishi-v8-12_files/dp0/FLTG/Fluent/FLTG-Setup-Output.cas.h5",
    "v8_13": "Qishi_v8_13_20261014/qishi-v8-13_files/dp0/FLTG/Fluent/FLTG-Setup-Output.cas.h5",
    "v8_14": "Qishi_v8_14_20261014/qishi-v8-14_files/dp0/FLTG/Fluent/FLTG-Setup-Output.cas.h5",
}

PROJ = r"F:/B_Simulation/fangzhen-3"


def run(cmd):
    r = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8", errors="replace")
    return r.returncode, (r.stdout or "") + (r.stderr or "")


def main():
    only = sys.argv[1:] if len(sys.argv) > 1 else list(STRUCTURES)
    os.makedirs(CASES, exist_ok=True)
    summary = []
    for name in only:
        cas = STRUCTURES[name]
        cas_abs = os.path.join(PROJ, cas)
        out = os.path.join(CASES, name)
        print(f"\n{'='*70}\n### {name}  ({os.path.basename(cas)})\n{'='*70}")
        if not os.path.isfile(cas_abs):
            print(f"  !! 算例文件不存在，跳过：{cas_abs}")
            summary.append((name, "缺失", ""))
            continue
        t0 = time.time()
        # 1) 网格
        rc, out1 = run([PY, os.path.join(HERE, "fluent_cas_to_foam.py"), cas_abs, out])
        print(out1[-1500:] if len(out1) > 1500 else out1)
        if rc != 0:
            print(f"  !! 网格转换失败 rc={rc}")
            summary.append((name, "网格失败", ""))
            continue
        # 2) 算例字典
        rc, out2 = run([PY, os.path.join(HERE, "generate_case.py"), out])
        if rc != 0:
            print(f"  !! 字典生成失败 rc={rc}\n{out2[-800:]}")
            summary.append((name, "字典失败", ""))
            continue
        dt = time.time() - t0
        # 提取单元数
        ncells = ""
        for line in out1.splitlines():
            if line.startswith("单元："):
                ncells = line.strip()
        print(f"  >> {name} 完成，用时 {dt:.0f}s   {ncells}")
        summary.append((name, "就绪", f"{dt:.0f}s"))

    print("\n\n===== 汇总 =====")
    for n, s, d in summary:
        print(f"  {n:8s} {s:10s} {d}")


if __name__ == "__main__":
    main()
