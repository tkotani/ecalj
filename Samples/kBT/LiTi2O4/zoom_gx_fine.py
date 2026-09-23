#!/usr/bin/env python3
"""High-resolution Gamma-X zoom of the t2g/eg manifold, bands labelled.
usage: zoom_gx_fine.py OUT.png "title" YLO YHI DIR LABEL [DIR LABEL ...]"""
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
runs=[(sys.argv[i],sys.argv[i+1]) for i in range(5,len(sys.argv),2)]
lo,hi=37,46
fig,axs=plt.subplots(1,len(runs),figsize=(8.5*len(runs),9),squeeze=False)
cmap=plt.get_cmap('tab10')
for c,(d,lab) in enumerate(runs):
    ax=axs[0][c]; b=seg1(d)
    for ib in range(lo,hi+1):
        if ib not in b: continue
        X=np.array([x for x,_ in b[ib]]); E=np.array([e for _,e in b[ib]])
        col=cmap((ib-lo)%10)
        ax.plot(X/X.max(),E,'-',lw=1.5,color=col)
        m=(E>ylo)&(E<yhi)
        if m.any():
            for frac in (0.02,0.5,0.98):
                k=int(frac*(len(X)-1))
                if ylo<E[k]<yhi:
                    ax.annotate(str(ib),(X[k]/X.max(),E[k]),textcoords='offset points',xytext=(0,5),
                                color=col,fontsize=10,fontweight='bold',ha='center')
    ax.set_xlim(0,1); ax.set_ylim(ylo,yhi); ax.set_xlabel('Γ → X  (k / k_X)'); ax.set_title(lab)
    ax.grid(color='0.9',lw=0.5)
    if c==0: ax.set_ylabel('E − E_F (eV)')
fig.suptitle(title,fontsize=14); fig.tight_layout(rect=(0,0,1,0.96)); fig.savefig(out,dpi=110)
