#!/usr/bin/env python3
"""NiO: conventional QSGW vs MLO-QSGW, LDA / iter 1 / iter 2 in rows.
usage: mloqsgw_rows.py REFDIR MLODIR OUT.png [nk]
  each DIR holds bnd_lda.dat, bnd_iter1.dat, bnd_iter2.dat"""
import sys, os, numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
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
def gap(E):
    """Fundamental gap: take the HIGHEST band whose maximum is still at/below E_F,
    then the distance to the next band's minimum. (E_F sits on the VBM here, so a
    sign test on E alone fails, and taking the largest gap finds the semicore one.)"""
    hi = E.max(axis=1); lo = E.min(axis=1)
    occ = np.where(hi <= 0.05)[0]
    if occ.size == 0 or occ[-1]+1 >= E.shape[0]: return (np.nan,)*3
    n = occ[-1]
    return (lo[n+1]-hi[n], hi[n], lo[n+1])

REF,MLO,OUT = sys.argv[1],sys.argv[2],sys.argv[3]
nk = sys.argv[4] if len(sys.argv)>4 else '3'
ROWS = [('LDA','bnd_lda.dat'), ('QSGW iter 1','bnd_iter1.dat'), ('QSGW iter 2','bnd_iter2.dat')]
fig,ax = plt.subplots(3,2,figsize=(11.5,12.5))
for r,(label,fn) in enumerate(ROWS):
    fr, fm = os.path.join(REF,fn), os.path.join(MLO,fn)
    if not (os.path.exists(fr) and os.path.exists(fm)):
        for c in (0,1): ax[r][c].text(0.5,0.5,f'{label}: missing',ha='center',transform=ax[r][c].transAxes)
        continue
    xr,er,efr = rd(fr); xm,em,efm = rd(fm)
    x = xr[0]; NB = min(24, er.shape[0], em.shape[0])
    ER = er; EM = em + efm - efr           # common zero = the conventional E_F of this row
    gr = gap(er); gm = gap(em)
    d = EM[:NB]-ER[:NB]
    for c,(ylo,yhi,ttl) in enumerate([(-9,12,'full range'),(-4,3,'zoom $-4$ to $+3$ eV')]):
        a = ax[r][c]
        for ib in range(NB):
            a.plot(x, ER[ib], '-',  lw=1.2, color='tab:blue')
            a.plot(x, EM[ib], '--', lw=1.2, color='tab:red')
        a.set_ylim(ylo,yhi); a.set_xlim(0,1); a.axhline(0,color='k',lw=0.9)
        a.grid(axis='y', which='major', color='0.75', lw=0.6)
        a.grid(axis='y', which='minor', color='0.9', lw=0.4)
        a.yaxis.set_major_locator(MultipleLocator(2.0 if c==0 else 1.0))
        a.yaxis.set_minor_locator(MultipleLocator(0.5 if c==0 else 0.25))
        a.set_ylabel(r'$E-E_F^{\rm conv}$ [eV]'); a.set_xlabel(r'$\Gamma \to X$')
        a.set_title(f'{label} — {ttl}', fontsize=10)
    ax[r][0].text(0.02,0.04,
        f'gap: conv {gr[0]:.2f} eV, MLO {gm[0]:.2f} eV   |   |diff| mean {np.abs(d).mean()*1e3:.0f}, max {np.abs(d).max()*1e3:.0f} meV',
        transform=ax[r][0].transAxes, fontsize=8)
ax[0][0].plot([],[],'-', color='tab:blue',label='conventional ($\\Sigma^{\\rm MTO}$ interpolated)')
ax[0][0].plot([],[],'--',color='tab:red', label='MLO ($\\Sigma^{\\rm MLO}$, 50 orb)')
ax[0][0].legend(fontsize=8, loc='lower right')
fig.suptitle(f'NiO {nk}$^3$ (AF, $n_k=n_1n_2n_3={nk}^3$), spin 1 — rows: LDA / QSGW iter 1 / iter 2', fontsize=11)
plt.tight_layout(rect=[0,0,1,0.975]); plt.savefig(OUT, dpi=130)
print('wrote', OUT)
