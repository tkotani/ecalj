#!/usr/bin/env python3
"""Fine Gamma-X with band indices; k points where a band is 2-fold degenerate (gap < 5 meV
with a neighbour) are overlaid with thick grey, so a doublet can be followed through crossings.
usage: zoom_gx_fine2.py OUT.png "title" YLO YHI BLO BHI DIR LABEL [DIR LABEL ...]"""
import sys, glob, os
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
plt.rcParams.update({'font.size':12})
def seg1(d):
    f=sorted(glob.glob(os.path.join(d,'bnd0*.spin1')))[0]
    b=defaultdict(list)
    for line in open(f):
        if line.startswith('#') or not line.strip(): continue
        t=line.split(); b[int(t[0])].append((float(t[1]),float(t[2])))
    return b
out,title=sys.argv[1],sys.argv[2]; ylo,yhi=float(sys.argv[3]),float(sys.argv[4])
lo,hi=int(sys.argv[5]),int(sys.argv[6])
runs=[(sys.argv[i],sys.argv[i+1]) for i in range(7,len(sys.argv),2)]
fig,axs=plt.subplots(1,len(runs),figsize=(7.6*len(runs),9),squeeze=False)
cmap=plt.get_cmap('tab10')
for c,(d,lab) in enumerate(runs):
    ax=axs[0][c]; b=seg1(d)
    X=np.array([x for x,_ in b[lo]]); X=X/X.max()
    E={ib:np.array([e for _,e in b[ib]]) for ib in range(lo,hi+1) if ib in b}
    for ib in sorted(E):
        deg=np.zeros(len(X),bool)
        for jb in (ib-1,ib+1):
            if jb in E: deg |= np.abs(E[ib]-E[jb])<0.005
        ax.plot(X,np.where(deg,E[ib],np.nan),'-',lw=5,color='0.8',zorder=1)
        ax.plot(X,E[ib],'-',lw=1.5,color=cmap((ib-lo)%10),zorder=2)
        for frac in (0.02,0.5,0.98):
            k=int(frac*(len(X)-1))
            if ylo<E[ib][k]<yhi:
                ax.annotate(str(ib),(X[k],E[ib][k]),textcoords='offset points',xytext=(0,5),
                            color=cmap((ib-lo)%10),fontsize=10,fontweight='bold',ha='center',zorder=3)
    ax.axhline(0,color='0.4',lw=0.6,ls='--')
    ax.set_xlim(0,1); ax.set_ylim(ylo,yhi); ax.set_xlabel('Γ → X  (k / k_X)'); ax.set_title(lab)
    ax.grid(color='0.92',lw=0.5)
    if c==0: ax.set_ylabel('E − E_F (eV)')
fig.suptitle(title+'   [thick grey = 2-fold degenerate within 5 meV]',fontsize=13)
fig.tight_layout(rect=(0,0,1,0.96)); fig.savefig(out,dpi=110)
