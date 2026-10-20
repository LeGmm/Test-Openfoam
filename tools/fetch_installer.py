#!/usr/bin/env python3
"""多源轮换 + 断点续传 + SHA1 校验的安装包下载器（对抗 GitHub 限速/断连）。

针对本项目实测环境：
  - 管道上限约 1MB/s（阿里云镜像实测 986KB/s）
  - GitHub 各通道被限到 40-170KB/s，且并发会被 RemoteDisconnected 打断
  → 策略：多镜像轮换（一个源失败换下一个，不空等）+ 分块续传 + 少量并发

用法：
  python fetch_installer.py [输出文件] [线程数]
特性：
  - 若目标文件已存在且大小正确，直接校验 SHA1 后退出（便于用户手动下载后接管）
  - 每块独立 .part，中断后重跑自动续传
"""
import hashlib
import os
import shutil
import sys
import threading
import time
import urllib.error
import urllib.request

UA = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 Chrome/120.0 Safari/537.36"
STAMP = "blueCFD-Core-2020-1-win64-setup.exe"
REPO = f"blueCFD/Core/releases/download/blueCFD-Core-2020-1/{STAMP}"
GH = f"https://github.com/{REPO}"

# 实测可用性排序（会逐个探测，取第一个报告出大小的作为主源，其余作轮换源）
SOURCES = [
    f"https://gh-proxy.com/{GH}",
    f"https://gh.llkk.cc/{GH}",
    f"https://ghproxy.net/{GH}",
    f"https://ghfast.top/{GH}",
    GH,
]

EXPECT_SHA1 = "31c78a935309ece053035e4d2a81ebac5f766c43"
DEFAULT_OUT = r"D:\APP\openfoam\_installer\blueCFD-Core-2020-1-win64-setup.exe"

CHUNK = 4 * 1024 * 1024
BLOCK = 131072
_lock = threading.Lock()


def log(msg):
    with _lock:
        print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def probe(url, timeout=20):
    req = urllib.request.Request(url, headers={"Range": "bytes=0-0", "User-Agent": UA})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        cr = r.headers.get("Content-Range")
        cl = r.headers.get("Content-Length")
        if cr and "/" in cr:
            return int(cr.split("/")[-1]), True
        if cl:
            return int(cl), r.status == 206
    return None, False


def alive_sources():
    out = []
    for u in SOURCES:
        try:
            size, ranged = probe(u)
            if size and ranged:
                out.append((u, size))
                log(f"可用源: {u[:52]}… size={size/1e6:.1f}MB")
        except Exception as e:
            log(f"跳过源: {u[:52]}… ({type(e).__name__})")
    if not out:
        raise SystemExit("没有可用下载源")
    return out


def fetch_chunk(idx, start, end, sources, parts, counter):
    """下载一块；在 sources 间轮换，每轮每源尝试一次。"""
    part = os.path.join(parts, f"{idx:05d}.part")
    need = end - start + 1

    for rnd in range(40):
        for url, _ in sources:
            have = os.path.getsize(part) if os.path.exists(part) else 0
            if have >= need:
                return True
            s = start + have
            req = urllib.request.Request(
                url, headers={"Range": f"bytes={s}-{end}", "User-Agent": UA})
            try:
                with urllib.request.urlopen(req, timeout=45) as r:
                    if r.status != 206 and s > 0:
                        skip = s
                        while skip > 0:
                            b = r.read(min(BLOCK, skip))
                            if not b:
                                break
                            skip -= len(b)
                    with open(part, "ab") as f:
                        while have < need:
                            b = r.read(BLOCK)
                            if not b:
                                break
                            if have + len(b) > need:
                                b = b[:need - have]
                            f.write(b)
                            have += len(b)
                            counter[0] += len(b)
                if have >= need:
                    return True
            except Exception:
                continue
        # 一整轮都没成功，退避
        time.sleep(min(3 + rnd, 15))

    return False


def sha1_of(path):
    h = hashlib.sha1()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(4 * 1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_OUT
    nthr = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    os.makedirs(os.path.dirname(out), exist_ok=True)

    # 已有人工下载的成品？直接校验
    if os.path.exists(out) and os.path.getsize(out) > 700_000_000:
        log(f"发现已存在文件 {os.path.getsize(out)/1e6:.1f}MB，校验 SHA1…")
        got = sha1_of(out)
        if got == EXPECT_SHA1:
            log("SHA1 通过，安装包就绪")
            return 0
        log(f"SHA1 不符（{got[:12]}…），将重新下载缺失部分")

    sources = alive_sources()
    size = max(s for _, s in sources)
    parts = out + ".parts"
    os.makedirs(parts, exist_ok=True)

    chunks = [(i, s, min(s + CHUNK - 1, size - 1))
              for i, s in enumerate(range(0, size, CHUNK))]
    counter = [0]
    todo = []
    for idx, s, e in chunks:
        p = os.path.join(parts, f"{idx:05d}.part")
        need = e - s + 1
        have = os.path.getsize(p) if os.path.exists(p) else 0
        if have >= need:
            counter[0] += have
        else:
            todo.append((idx, s, e))

    log(f"总 {len(chunks)} 块 / 待下 {len(todo)} 块 / "
        f"已有 {counter[0]/1e6:.1f}MB / 目标 {size/1e6:.1f}MB / {nthr} 线程")

    stop = threading.Event()

    def reporter():
        last, lt = counter[0], time.time()
        while not stop.wait(30):
            cur, now = counter[0], time.time()
            spd = (cur - last) / max(now - last, 1e-6)
            last, lt = cur, now
            eta = (size - cur) / spd / 60 if spd > 1024 else -1
            log(f"进度 {cur/1e6:.1f}/{size/1e6:.1f}MB ({100*cur/size:.1f}%)  "
                f"{spd/1024:.0f} KB/s  剩余 {'—' if eta < 0 else f'{eta:.0f} 分钟'}")

    threading.Thread(target=reporter, daemon=True).start()

    q, qlock, results = list(todo), threading.Lock(), []

    def worker():
        while True:
            with qlock:
                if not q:
                    return
                i, s, e = q.pop(0)
            results.append(fetch_chunk(i, s, e, sources, parts, counter))

    ths = [threading.Thread(target=worker, daemon=True) for _ in range(nthr)]
    [t.start() for t in ths]
    [t.join() for t in ths]
    stop.set()

    if not all(results) or len(results) != len(todo):
        log("仍有块未完成，重跑本脚本可续传")
        return 2

    log("拼接中…")
    tmp = out + ".merging"
    with open(tmp, "wb") as w:
        for i, _, _ in chunks:
            with open(os.path.join(parts, f"{i:05d}.part"), "rb") as r:
                shutil.copyfileobj(r, w, 4 * 1024 * 1024)
    if os.path.getsize(tmp) != size:
        log("拼接大小不符，保留 parts")
        return 3
    os.replace(tmp, out)

    log("校验 SHA1…")
    got = sha1_of(out)
    if got == EXPECT_SHA1:
        log(f"SHA1 通过 → {out}")
        shutil.rmtree(parts, ignore_errors=True)
        return 0
    log(f"SHA1 不符：期望 {EXPECT_SHA1} 实际 {got}（保留 parts）")
    return 4


if __name__ == "__main__":
    sys.exit(main())
