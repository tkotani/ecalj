#!/usr/bin/env python3
"""LiTi2O4 6^3, iterations 1 and 2: conventional / MLO without the chi~ freeze / MLO with it."""
import sys, os, numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.sans-serif'] = ['IPAGothic','Droid Sans Fallback','DejaVu Sans']
matplotlib.rcParams['axes.unicode_minus'] = False
from matplotlib.ticker import MultipleLocator
ROOT, OUT = sys.argv[1], sys.argv[2]
def rd(f):
    xs,es,cx,ce=[],[],[],[]
    for ln in open(f):
        if ln.startswith('#'): continue
        t=ln.split()
        if len(t)<3:
            if cx: xs.append(cx); es.append(ce); cx=[];ce=[]
            continue
        cx.append(float(t[1])); ce.append(float(t[2]))
    if cx: xs.append(cx); es.append(ce)
    return np.array(xs[0]), np.array(es)
MESH=[0.0,1/3,2/3,1.0]; IDX=[0,70,140,210]; LO,HI=32,44
def chord(y,x):
    return np.array([np.abs(y[a:b+1]-np.interp(x[a:b+1],[x[a],x[b]],[y[a],y[b]])).max()
                     for a,b in zip(IDX[:-1],IDX[1:])])*1e3
ROWS=[('従来 QSGW (MTO)', [f'{ROOT}/liti_ref/bnd_iter1.dat', f'{ROOT}/refpair/bnd_iter2.dat'],'tab:blue'),
      ('MLO — $\\tilde\\chi$ 凍結なし（バグ）',[f'{ROOT}/liti_mlo/bnd_iter1.dat',f'{ROOT}/liti_mlo/bnd_iter2.dat'],'tab:red'),
      ('MLO — $\\tilde\\chi$ 凍結あり（修正後）',[f'{ROOT}/frozen/bnd_iter1.dat',f'{ROOT}/frozen/bnd_iter2.dat'],'tab:green')]
x0,Elda = rd(f'{ROOT}/liti_mlo/bnd_lda.dat'); Slda=np.sort(Elda[LO:HI],axis=0)
wl=np.mean([b.max()-b.min() for b in Slda])*1e3; cl=np.array([chord(b,x0) for b in Slda]).mean()
fig,AX=plt.subplots(3,3,figsize=(14.5,11.2),sharex=True,sharey=True)
for r,(lab,fs,col) in enumerate(ROWS):
    ax=AX[r][0]
    for b in Slda: ax.plot(x0,b,'-',lw=1.0,color='0.45')
    ax.set_title(f'LDA（内挿なし）\n帯幅 {wl:.0f} / 弦 {cl:.1f} = {cl/wl*100:.1f} %',fontsize=9.5)
    for c,f in enumerate(fs,1):
        ax=AX[r][c]
        if not os.path.exists(f): ax.axis('off'); continue
        x,E=rd(f); S=np.sort(E[LO:HI],axis=0)
        for b in S: ax.plot(x,b,'-',lw=1.15,color=col)
        w=np.mean([b.max()-b.min() for b in S])*1e3
        ch=np.array([chord(b,x) for b in S]).mean()
        ax.set_title(f'{lab}  iter {c}\n帯幅 {w:.0f} / 弦 {ch:.1f} = {ch/w*100:.1f} %',fontsize=9.5)
    for c in range(3):
        a=AX[r][c]; a.set_ylim(-0.9,1.5); a.set_xlim(0,1); a.axhline(0,color='k',lw=0.9)
        for q in MESH: a.axvline(q,color='k',ls=':',lw=0.9)
        a.yaxis.set_major_locator(MultipleLocator(0.5)); a.yaxis.set_minor_locator(MultipleLocator(0.1))
        a.grid(axis='y',which='major',color='0.78',lw=0.6); a.grid(axis='y',which='minor',color='0.92',lw=0.4)
    AX[r][0].set_ylabel('$E-E_F$ [eV]')
for a in AX[-1]: a.set_xlabel('$\\Gamma \\to X$')
fig.suptitle('LiTi$_2$O$_4$ 6$^3$ t$_{2g}$: 「弦ずれ/帯幅」が小さいほどメッシュ点間が滑らか。'
             'LDA の 7.9 % が下限（内挿なし）', fontsize=12.5)
plt.tight_layout(rect=[0,0,1,0.955]); plt.savefig(OUT,dpi=120)
print('wrote',OUT)
