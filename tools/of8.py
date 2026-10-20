#!/usr/bin/env python3
"""在 blueCFD-Core 2020 (OpenFOAM 8) 环境中运行任意 OpenFOAM 命令。

为什么需要它：
  blueCFD 的环境必须由 setvars_OF8.bat 建立（WM_OPTIONS/WM_TARGET_ARCH/PATH/DLL 目录），
  且该 bat 会 `cd` 到安装根目录，因此不能直接在本机 bash 里 source etc/bashrc。
  这里用 cmd 调用 setvars，再切回目标算例目录执行命令。

环境变量：
  OF8_HOME    blueCFD 安装根目录；默认 D:\\APP\\openfoam
  （与 ablate.py 文档里写的覆盖方式保持一致）

用法：
  python of8.py <case 目录> "<OpenFOAM 命令>" [超时秒数]
示例：
  python of8.py ../cases/v8_1 "checkMesh"
  python of8.py ../cases/v8_1 "multiphaseEulerFoam" 600
"""
import os
import subprocess
import sys

BLUECFD = os.environ.get("OF8_HOME", r"D:\APP\openfoam")
SETVARS = os.path.join(BLUECFD, "setvars_OF8.bat")

# cmd 元字符：shell=True 的字符串拼接，路径里出现这些会被解释甚至截断
_UNSAFE = set('"&|<>^%')


def check_env() -> None:
    if not os.path.isfile(SETVARS):
        raise SystemExit(
            f"找不到 blueCFD 环境脚本：{SETVARS}\n"
            f"  请确认 blueCFD-Core 已安装，或用环境变量 OF8_HOME 指向安装根目录，例如：\n"
            f"  set OF8_HOME=D:\\APP\\openfoam")


def run(case_dir, command, timeout=None, log_file=None, echo=True):
    """在 blueCFD 环境中执行 command；返回 (returncode, output)。"""
    case_dir = os.path.abspath(case_dir)
    if not os.path.isdir(case_dir):
        raise SystemExit(f"算例目录不存在: {case_dir}")
    bad = _UNSAFE & set(case_dir)
    if bad:
        raise SystemExit(f"算例路径含 cmd 元字符 {''.join(sorted(bad))}，请改名后重试: {case_dir}")

    # 用 && 串联：先建环境，再回到算例目录，最后执行
    full = f'call "{SETVARS}" && cd /d "{case_dir}" && {command}'

    p = subprocess.run(full, shell=True, capture_output=True, text=True,
                       errors="replace", timeout=timeout)
    out = (p.stdout or "") + (p.stderr or "")

    if log_file:
        with open(log_file, "w", encoding="utf-8", errors="replace") as f:
            f.write(out)
    if echo:
        print(out)
    return p.returncode, out


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 1
    check_env()
    case_dir, command = sys.argv[1], sys.argv[2]
    timeout = int(sys.argv[3]) if len(sys.argv) > 3 else None
    rc, _ = run(case_dir, command, timeout=timeout)
    print(f"\n[of8] 退出码 = {rc}")
    return rc


if __name__ == "__main__":
    sys.exit(main())
