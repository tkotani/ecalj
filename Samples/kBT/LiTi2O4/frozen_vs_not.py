#!/usr/bin/env python3
"""LiTi2O4 6^3 iteration 1: LDA / conventional / MLO without the chi~ freeze / MLO with it."""
import sys, os, numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.sans-serif'] = ['IPAGothic', 'Droid Sans Fallback', 'DejaVu Sans']
matplotlib.rcParams['axes.unicode_minus'] = False
from matplotlib.ticker import MultipleLocator
ROOT, OUT = sys.argv[1], sys.argv[2]
def rd(f):
    xs, es, cx, ce = [], [], [], []
    for ln in open(f):
        if ln.startswith('#'): continue
        t = ln.split()
        if len(t) < 3:
            if cx: xs.append(cx); es.append(ce); cx = []; ce = []
            continue
        cx.append(float(t[1])); ce.append(float(t[2]))
    if cx: xs.append(cx); es.append(ce)
    return np.array(xs[0]), np.array(es)
MESH=[0.0,1/3,2/3,1.0]; IDX=[0,70,140,210]; LO,HI=32,44
def chord(y,x):
    return np.array([np.abs(y[a:b+1]-np.interp(x[a:b+1],[x[a],x[b]],[y[a],y[b]])).max()
                     for a,b in zip(IDX[:-1],IDX[1:])])*1e3
C=[('LDA', f'{ROOT}/liti_mlo/bnd_lda.dat','0.35'),
   ('従来 QSGW (MTO) iter 1', f'{ROOT}/liti_ref/bnd_iter1.dat','tab:blue'),
   ('MLO iter 1 — $\\tilde\\chi$ 凍結 **なし**（バグ）', f'{ROOT}/liti_mlo/bnd_iter1.dat','tab:red'),
   ('MLO iter 1 — $\\tilde\\chi$ 凍結 **あり**（修正後）', f'{ROOT}/frozen/bnd_iter1.dat','tab:green')]
C=[c for c in C if os.path.exists(c[1])]
fig,AX=plt.subplots(1,len(C),figsize=(4.0*len(C),4.6),sharey=True)
for ax,(lab,f,col) in zip(np.atleast_1d(AX),C):
    x,E=rd(f); S=np.sort(E[LO:HI],axis=0)
    for b in S: ax.plot(x,b,'-',lw=1.2,color=col)
    w=np.mean([b.max()-b.min() for b in S])*1e3
    c=np.array([chord(b,x) for b in S]).mean()
    ax.set_ylim(-0.9,1.5); ax.set_xlim(0,1); ax.axhline(0,color='k',lw=0.9)
    for q in MESH: ax.axvline(q,color='k',ls=':',lw=0.9)
    ax.yaxis.set_major_locator(MultipleLocator(0.5)); ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax.grid(axis='y',which='major',color='0.78',lw=0.6); ax.grid(axis='y',which='minor',color='0.92',lw=0.4)
    ax.set_title(f'{lab}\n帯幅 {w:.0f} meV / 弦ずれ {c:.1f} = {c/w*100:.1f} %',fontsize=10)
    ax.set_xlabel('$\\Gamma \\to X$')
np.atleast_1d(AX)[0].set_ylabel('$E-E_F$ [eV]')
fig.suptitle('LiTi$_2$O$_4$ 6$^3$ iter 1: $\\tilde\\chi$ の凍結が入っていないと MLO の $\\Sigma$ は行き過ぎる '
             '— 点線 = $\\Sigma$ の q メッシュ点', fontsize=12)
plt.tight_layout(rect=[0,0,1,0.92]); plt.savefig(OUT,dpi=130)
print('wrote',OUT)
