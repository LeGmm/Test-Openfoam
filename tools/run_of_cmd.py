#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""在 blueCFD 环境里执行 OF 命令，并把输出【实时】写进日志文件。

与 of8.py 的区别：of8.py 用 capture_output 收集，命令结束才落盘，长跑无法监控。
本脚本用 Popen 直接把 stdout/stderr 重定向到文件，可 tail -f 跟踪，适合挂长任务。

用法：
    python tools/run_of_cmd.py <case 目录> "<命令>" <日志文件>
示例：
    python tools/run_of_cmd.py case_v81_1tau "mpiexec -localonly -np 8 multiphaseEulerFoam -parallel" log.1tau
"""
import os
import subprocess
import sys

BLUECFD = os.environ.get("OF8_HOME", r"D:\APP\openfoam")
SETVARS = os.path.join(BLUECFD, "setvars_OF8.bat")
_UNSAFE = set('"&|<>^%')


def main() -> int:
    if len(sys.argv) < 4:
        print(__doc__)
        return 2
    case = os.path.abspath(sys.argv[1])
    command = sys.argv[2]
    log = sys.argv[3]
    # 支持 @命令文件：命令里含空格/引号（如 mpiexec 全路径）时可避开 shell 引号地狱
    if command.startswith("@"):
        command = open(command[1:], encoding="utf-8").read().strip()
    if not os.path.isdir(case):
        raise SystemExit(f"算例目录不存在: {case}")
    if not os.path.isfile(SETVARS):
        raise SystemExit(f"找不到 {SETVARS}（用 OF8_HOME 覆盖）")
    bad = _UNSAFE & set(case)
    if bad:
        raise SystemExit(f"路径含 cmd 元字符 {''.join(sorted(bad))}: {case}")

    full = f'call "{SETVARS}" > nul 2>&1 && cd /d "{case}" && {command}'
    with open(log, "w", encoding="utf-8", errors="replace", buffering=1) as f:
        f.write(f"# cmd: {command}\n# case: {case}\n# started: {__import__('datetime').datetime.now()}\n")
        f.flush()
        p = subprocess.Popen(full, shell=True, stdout=f, stderr=subprocess.STDOUT,
                             cwd=case, text=True, errors="replace")
        rc = p.wait()
        f.write(f"\n# finished: {__import__('datetime').datetime.now()}  rc={rc}\n")
    print(f"[run_of_cmd] rc={rc}  log={log}")
    return rc


if __name__ == "__main__":
    sys.exit(main())
