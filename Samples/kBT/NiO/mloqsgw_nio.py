#!/usr/bin/env python3
"""NiO 2^3 (AF, nspin=2): conventional QSGW vs MLO-QSGW.
usage: mloqsgw_nio.py REF.dat MLO.dat OUT.png"""
import sys, numpy as np, matplotlib; matplotlib.use('Agg')
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
xr,er,efr = rd(sys.argv[1]); xm,em,efm = rd(sys.argv[2])
x = xr[0]; NB = min(24, er.shape[0], em.shape[0])
# plot BOTH on the conventional E_F, otherwise the E_F difference shows up as a fake rigid shift
ER = er + efr - efr; EM = em + efm - efr
fig,ax = plt.subplots(1,2,figsize=(11.5,4.6))
for ib in range(NB):
    ax[0].plot(x, ER[ib], '-',  lw=1.2, color='tab:blue')
    ax[0].plot(x, EM[ib], '--', lw=1.2, color='tab:red')
ax[0].plot([],[],'-', color='tab:blue', label='conventional QSGW ($\\Sigma^{\\rm MTO}$ interpolated)')
ax[0].plot([],[],'--',color='tab:red',  label='MLO-QSGW ($\\Sigma^{\\rm MLO}$, 50 orb)')
ax[0].set_ylabel('$E-E_F^{\\rm conv}$ [eV]  (common zero)'); ax[0].legend(fontsize=8, loc='lower right')
ax[0].set_title('NiO (AF), self-consistent bands, spin 1 (2 QSGW iterations)', fontsize=10)
d = EM[:NB]-ER[:NB]
for ib in range(NB):
    ax[1].plot(x, ER[ib], '-',  lw=1.3, color='tab:blue')
    ax[1].plot(x, EM[ib], '--', lw=1.3, color='tab:red')
ax[1].set_ylabel(r'$E-E_F^{\rm conv}$ [eV]')
ax[1].set_ylim(-4, 3)
ax[1].set_title('zoom: $-4$ to $+3$ eV', fontsize=10)
ax[1].text(0.02,0.04,f'|diff| over b1-b{NB}: mean {np.abs(d).mean()*1e3:.0f}, max {np.abs(d).max()*1e3:.0f} meV',
           transform=ax[1].transAxes, fontsize=8)
for a in ax:
    for m in (0.0,1.0): a.axvline(m,color='k',ls=':',lw=0.8)
    a.set_xlabel('$\\Gamma \\to X$'); a.set_xlim(0,1)
    a.grid(axis='y', which='major', color='0.75', lw=0.6)
    a.grid(axis='y', which='minor', color='0.9',  lw=0.4)
ax[0].yaxis.set_major_locator(MultipleLocator(2.0))   # eV
ax[0].set_ylim(-9, 12)
ax[0].yaxis.set_minor_locator(MultipleLocator(0.5))
ax[0].axhline(0.0, color='k', lw=0.9)                 # E_F
ax[1].yaxis.set_major_locator(MultipleLocator(1.0))   # eV
ax[1].yaxis.set_minor_locator(MultipleLocator(0.25))
ax[1].axhline(0.0, color='k', lw=0.9)
nk = sys.argv[5] if len(sys.argv)>5 else '2'
fig.suptitle(f'NiO {nk}$^3$ ($n_k=n_1n_2n_3={nk}^3$, AF) — dotted = $\\Sigma$ mesh points', fontsize=10)
plt.tight_layout(rect=[0,0,1,0.94]); plt.savefig(sys.argv[3], dpi=140)
print('wrote', sys.argv[3])
