#!/usr/bin/env python3
"""Si 2^3: conventional QSGW vs MLO-QSGW (Sigma interpolated in the MLO space).
usage: mloqsgw_si.py REF.dat MLO.dat OUT.png"""
import sys, numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
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
x = xr[0]; NB = 8
fig,ax = plt.subplots(1,2,figsize=(11.5,4.6))
for ib in range(NB):
    ax[0].plot(x, er[ib], '-',  lw=1.2, color='tab:blue')
    ax[0].plot(x, em[ib], '--', lw=1.2, color='tab:red')
ax[0].plot([],[],'-', color='tab:blue', label='conventional QSGW ($\\Sigma^{\\rm MTO}$ interpolated)')
ax[0].plot([],[],'--',color='tab:red',  label='MLO-QSGW ($\\Sigma^{\\rm MLO}$, 18 orb)')
ax[0].set_ylabel('$E-E_F$ [eV]'); ax[0].legend(fontsize=8, loc='lower right')
ax[0].set_title('Si, self-consistent bands (2 QSGW iterations)', fontsize=10)
for ib in range(NB):
    ax[1].plot(x, ((em[ib]+efm)-(er[ib]+efr))*1e3, '-', lw=1.0, label=f'b{ib+1}')
ax[1].set_ylabel('MLO-QSGW $-$ conventional [meV]')
ax[1].set_title('difference of the self-consistent solutions', fontsize=10)
ax[1].legend(fontsize=7, ncol=4)
d = (em[:NB]+efm)-(er[:NB]+efr)
ax[1].text(0.02,0.05,f'|diff| mean {np.abs(d).mean()*1e3:.1f}, max {np.abs(d).max()*1e3:.1f} meV',
           transform=ax[1].transAxes, fontsize=8)
for a in ax:
    for m in (0.0,1.0): a.axvline(m,color='k',ls=':',lw=0.8)
    a.set_xlabel('$\\Gamma \\to X$'); a.set_xlim(0,1)
fig.suptitle('Si 2$^3$ ($n_{k}=n_{1}n_{2}n_{3}=2^3$) — dotted = $\\Sigma$ mesh points ($\\Gamma$, X)', fontsize=10)
plt.tight_layout(rect=[0,0,1,0.94]); plt.savefig(sys.argv[3], dpi=140)
print('wrote', sys.argv[3])
