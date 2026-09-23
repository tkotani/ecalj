#!/usr/bin/env python3
"""Occupied t2g bands 33/34 along Gamma-X: LDA vs 6^3 vs 9^3, plus the residual after a smooth (5th-order) fit."""
import glob, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from collections import defaultdict
plt.rcParams.update({'font.size':12})
def seg1(d):
    f=sorted(glob.glob(d+'/bnd0*.spin1')) [0]
    b=defaultdict(list)
    for line in open(f):
        if line.startswith('#') or not line.strip(): continue
        t=line.split(); b[int(t[0])].append((float(t[1]),float(t[2])))
    return b
runs=(('LDA','GXfineLDA',0,'0.45'),('6^3 P iter 8','GXfine666',6,'C0'),('9^3 P iter 15','GXfine999',9,'C3'))
fig,axs=plt.subplots(2,2,figsize=(15,10))
for col,ib in enumerate((33,34)):
    for lab,d,N,c in runs:
        b=seg1(d); X=np.array([x for x,_ in b[ib]]); X=X/X.max(); E=np.array([e for _,e in b[ib]])
        axs[0,col].plot(X,E,'-',lw=1.8,color=c,label=lab)
        m=(X>0.02)&(X<0.98)
        res=E[m]-np.polyval(np.polyfit(X[m],E[m],5),X[m])
        axs[1,col].plot(X[m],res*1e3,'-',lw=1.8,color=c,label=f'{lab}  (rms {res.std()*1e3:.1f} meV)')
        if N:
            mesh=[2*mm/N for mm in range(N//2+1) if 2*mm/N<=1.0]
            idx=[int(round(t*(len(X)-1))) for t in mesh]
            axs[0,col].plot(X[idx],E[idx],'o',ms=9,mfc='none',mew=2,color=c)
            idm=[int(round(t*(m.sum()-1))) for t in mesh]
            axs[1,col].plot(X[m][idm],res[idm]*1e3,'o',ms=9,mfc='none',mew=2,color=c)
    axs[0,col].axhline(0,color='0.4',lw=0.7,ls='--')
    axs[0,col].set_title(f'band {ib}  (occupied t2g, crosses E_F)'); axs[0,col].set_ylabel('E − E_F (eV)')
    axs[1,col].axhline(0,color='0.4',lw=0.7,ls='--')
    axs[1,col].set_title(f'band {ib}: residual after a smooth 5th-order fit'); axs[1,col].set_ylabel('residual (meV)')
    for r in (0,1):
        axs[r,col].set_xlim(0,1); axs[r,col].grid(color='0.92'); axs[r,col].set_xlabel('Γ → X  (k / k_X)')
        axs[r,col].legend(fontsize=10)
fig.suptitle('Occupied t2g along Γ–X (211 k points); open circles = Σ q-mesh points of each run',fontsize=13)
fig.tight_layout(rect=(0,0,1,0.96)); fig.savefig('t2g_occ.png',dpi=105)
