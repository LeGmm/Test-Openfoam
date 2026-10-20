#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""找出 polyMesh 内部面中 (owner, neighbour) 成对重复的面。"""
import os
import sys

import numpy as np

mesh = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    os.path.dirname(os.path.abspath(__file__)), 'case', 'constant', 'polyMesh')


def read_scalar_list(path):
    with open(path, encoding='utf-8', errors='replace') as f:
        txt = f.read()
    body = txt.split('(', 1)[1].rsplit(')', 1)[0]
    return np.array([int(x) for x in body.split()], dtype=np.int64)


owner = read_scalar_list(os.path.join(mesh, 'owner'))
neigh = read_scalar_list(os.path.join(mesh, 'neighbour'))
n_internal = len(neigh)

pairs = np.stack([owner[:n_internal], neigh[:n_internal]], axis=1)
# 排序后找重复（同 owner 同 neighbour 出现 >= 2 次的内部面）
order = np.lexsort((pairs[:, 1], pairs[:, 0]))
sp = pairs[order]
same_prev = np.all(sp[1:] == sp[:-1], axis=1)
dup_rows = np.where(same_prev)[0]
dup_face_set = set()
for r in dup_rows:
    dup_face_set.add(order[r])
    dup_face_set.add(order[r + 1])
dups = sorted(dup_face_set)

print(f'内部面总数: {n_internal}')
print(f'(owner,neighbour) 重复的内部面数: {len(dups)}')
for fidx in dups:
    o, n = pairs[fidx]
    same = [i for i in range(n_internal)
            if pairs[i, 0] == o and pairs[i, 1] == n]
    print(f'  面 {fidx}: owner={o} neighbour={n}  重复组面号={same}')

# 顺便查反向重复（owner/neighbour 互换，即两 cell 间还有一条反向面）
rev = np.stack([pairs[:, 1], pairs[:, 0]], axis=1)
allpairs = np.vstack([pairs, rev])
o2 = np.lexsort((allpairs[:, 1], allpairs[:, 0]))
sp2 = allpairs[o2]
rd = np.all(sp2[1:] == sp2[:-1], axis=1)
print(f'owner/neighbour 互换意义下的重复对数: {rd.sum()} (含正常重复)')
