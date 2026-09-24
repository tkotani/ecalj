#!/usr/bin/env python3
"""LiTi2O4 6^3 MLO-QSGW: t2g bands along Gamma-X, one row per iteration.
usage: mlo_chain_rows.py DIR OUT.png [label ...]"""
import sys, os, glob, numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
D, OUT = sys.argv[1], sys.argv[2]
def rd(f):
    xs=[];es=[];cx=[];ce=[];ef=None
    for ln in open(f):
        if ln.startswith('#'):
            if 'base=eferm' in ln: ef=float(ln.split()[1])
            continue
        t=ln.split()
        if len(t)<3:
            if cx: xs.append(cx);es.append(ce);cx=[];ce=[]
            continue
        cx.append(float(t[1])); ce.append(float(t[2]))
    if cx: xs.append(cx);es.append(ce)
    return np.array(xs),np.array(es),ef
files=[('LDA',f'{D}/bnd_lda.dat')]
for f in sorted(glob.glob(f'{D}/bnd_iter*.dat'), key=lambda s:int(''.join(c for c in os.path.basename(s) if c.isdigit()))):
    files.append((f"iter {''.join(c for c in os.path.basename(f) if c.isdigit())}", f))
files=[(l,f) for l,f in files if os.path.exists(f)]
MESH=[0.0,1/3,2/3,1.0]; IDX=[0,70,140,210]
def chord(y,x):
    o=[]
    for a,b in zip(IDX[:-1],IDX[1:]):
        c=np.interp(x[a:b+1],[x[a],x[b]],[y[a],y[b]]); o.append(np.abs(y[a:b+1]-c).max())
    return np.array(o)*1e3
n=len(files)
fig,ax=plt.subplots(n,1,figsize=(7.5,3.1*n), squeeze=False)
for r,(lab,f) in enumerate(files):
    x_,E,ef = rd(f); x=x_[0]; a=ax[r][0]
    lo,hi = (32,44) if E.shape[0]>44 else (0,min(12,E.shape[0]))
    for ib in range(lo,hi): a.plot(x,E[ib],'-',lw=1.2)
    t=np.array([chord(E[ib],x) for ib in range(lo,hi)])
    a.set_ylim(-0.8,1.4); a.set_xlim(0,1); a.axhline(0,color='k',lw=0.9)
    for m in MESH: a.axvline(m,color='k',ls=':',lw=0.8)
    a.grid(axis='y',which='major',color='0.75',lw=0.6); a.grid(axis='y',which='minor',color='0.9',lw=0.4)
    a.yaxis.set_major_locator(MultipleLocator(0.5)); a.yaxis.set_minor_locator(MultipleLocator(0.1))
    a.set_ylabel('$E-E_F$ [eV]'); a.set_title(f'{lab}  —  chord dev: all {t.mean():.1f}, middle {t[:,1].mean():.1f} meV', fontsize=10)
ax[-1][0].set_xlabel('$\\Gamma \\to X$')
fig.suptitle('LiTi$_2$O$_4$ 6$^3$ MLO-QSGW (76 orb, pwmode=11) — t$_{2g}$ bands; dotted = $\\Sigma$ mesh points', fontsize=11)
plt.tight_layout(rect=[0,0,1,0.97]); plt.savefig(OUT,dpi=130)
print('wrote',OUT)
