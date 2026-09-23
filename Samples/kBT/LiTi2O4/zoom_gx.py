#!/usr/bin/env python3
"""Zoom of the Gamma-X segment with band indices labelled on each branch.
usage: zoom_gx.py OUT.png "title" DIR1 LABEL1 [DIR2 LABEL2 ...]  (uses bnd001.spin1 = first symline segment)"""
import sys, glob, os
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
plt.rcParams.update({'font.size':11})
def seg1(d):
    f=sorted(glob.glob(os.path.join(d,'bnd0*.spin1')))[0]
    b=defaultdict(list)
    for line in open(f):
        if line.startswith('#') or not line.strip(): continue
        t=line.split(); b[int(t[0])].append((float(t[1]),float(t[2])))
    return b
out,title=sys.argv[1],sys.argv[2]
runs=[(sys.argv[i],sys.argv[i+1]) for i in range(3,len(sys.argv),2)]
lo,hi=33,46
fig,axs=plt.subplots(1,len(runs),figsize=(7.5*len(runs),8.5),squeeze=False)
cmap=plt.get_cmap('tab20')
for c,(d,lab) in enumerate(runs):
    ax=axs[0][c]; b=seg1(d)
    for ib in range(lo,hi+1):
        if ib not in b: continue
        X=np.array([x for x,_ in b[ib]]); E=np.array([e for _,e in b[ib]])
        col=cmap((ib-lo)%20)
        ax.plot(X,E,'-o',lw=1.4,ms=3.2,color=col)
        # label at both ends and at the point of largest |2nd diff|
        ax.annotate(str(ib),(X[0],E[0]),textcoords='offset points',xytext=(-16,-4),color=col,fontsize=9,fontweight='bold')
        ax.annotate(str(ib),(X[-1],E[-1]),textcoords='offset points',xytext=(5,-4),color=col,fontsize=9,fontweight='bold')
        if len(E)>4:
            k=int(np.abs(np.diff(E,2)).argmax())+1
            if abs(np.diff(E,2)[k-1])>0.03:
                ax.annotate(str(ib),(X[k],E[k]),textcoords='offset points',xytext=(0,7),color=col,fontsize=9,fontweight='bold',ha='center')
    ax.axhline(0,color='0.4',lw=0.6,ls='--')
    ax.set_xlim(-0.06,1.06); ax.set_ylim(-0.6,1.65); ax.set_xlabel('Γ → X'); ax.set_title(lab)
    ax.grid(axis='y',color='0.9',lw=0.5)
    if c==0: ax.set_ylabel('E − E_F (eV)')
fig.suptitle(title,fontsize=13); fig.tight_layout(rect=(0,0,1,0.96)); fig.savefig(out,dpi=110)
