# -*- coding: utf-8 -*-
"""诊断：读某时刻的 internalField，按单元统计最大量并定位单元。"""
import sys, re, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ofmesh import read_mesh

def read_internal(path, ncomp):
    txt = open(path, encoding='utf-8', errors='replace').read()
    m = re.search(r'internalField\s+nonuniform\s+List<(?:vector|scalar)>\s*(\d+)\s*\((.*?)\)\s*;', txt, re.S)
    if m:
        n = int(m.group(1))
        v = [float(x) for x in m.group(2).replace('(',' ').replace(')',' ').split()]
        return np.array(v).reshape(-1, ncomp), n
    m = re.search(r'internalField\s+uniform\s+([^;]+);', txt)
    return None, None

case, t = sys.argv[1], sys.argv[2]
pts, faces, own, nei, patches = read_mesh(os.path.join(case, 'constant', 'polyMesh'))
nf = len(faces); nint = len(nei)
avec = np.zeros((nf,3)); fc = np.zeros((nf,3)); area=np.zeros(nf)
for i, ids in enumerate(faces):
    p = pts[ids]; a = np.zeros(3)
    for j in range(len(p)): a += np.cross(p[j], p[(j+1)%len(p)])
    avec[i] = a*0.5; fc[i] = p.mean(0); area[i]=np.linalg.norm(a)*0.5
term = np.einsum('ij,ij->i', avec, fc)
nc = int(max(own.max(), nei.max()))+1
v = np.bincount(own, weights=term, minlength=nc)
v = np.bincount(nei, weights=-term[:nint], minlength=nc) + v
v = np.abs(v)/3.
amax = np.zeros(nc); np.maximum.at(amax, own, area); np.maximum.at(amax, nei, area[:nint])
csum = np.zeros((nc,3)); np.add.at(csum, own, fc*area[:,None]); np.add.at(csum, nei, fc[:nint]*area[:nint,None])
asum = np.zeros(nc); np.add.at(asum, own, area); np.add.at(asum, nei, area[:nint])
cc = csum/np.maximum(asum,1e-30)[:,None]

def stat(name, ncomp, log=True):
    f = os.path.join(case, t, name)
    if not os.path.isfile(f):
        print(f'  {name:14s} (无此文件)'); return None
    arr, n = read_internal(f, ncomp)
    if arr is None:
        print(f'  {name:14s} uniform'); return None
    if ncomp == 3:
        mag = np.linalg.norm(arr, axis=1)
        k = int(mag.argmax())
        print(f'  {name:14s} max|.|={mag[k]:12.4e} @cell {k:6d} '
              f'xyz=({cc[k,0]*1e3:6.2f},{cc[k,1]*1e3:6.2f},{cc[k,2]*1e3:6.2f})mm V={v[k]:.2e}')
        print(f'                 p99.9={np.percentile(mag,99.9):.4e}  中位={np.median(mag):.4e}')
        return mag, v, cc, k
    else:
        a = arr.ravel(); k = int(np.argmax(np.abs(a)))
        print(f'  {name:14s} min={a.min():12.4e} max={a.max():12.4e} @cell {k:6d} '
              f'xyz=({cc[k,0]*1e3:6.2f},{cc[k,1]*1e3:6.2f},{cc[k,2]*1e3:6.2f})mm V={v[k]:.2e}')
        return a, v, cc, k

print(f'=== {case} @ t={t}  单元 {nc} ===')
for nm, ncp in (('U.gas',3),('U.liquid',3),('p_rgh',1),('p',1),('alpha.liquid',1),('alpha.gas',1),
                ('T.gas',1),('T.liquid',1),('H2O.gas',1)):
    stat(nm, ncp)

print('\n--- 体积最小的 12 个单元 ---')
idx = np.argsort(v)[:12]
for i in idx:
    print(f'  cell {i:6d}  V={v[i]:.3e}  xyz=({cc[i,0]*1e3:6.2f},{cc[i,1]*1e3:6.2f},{cc[i,2]*1e3:6.2f})mm')
