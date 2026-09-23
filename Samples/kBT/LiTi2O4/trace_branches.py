#!/usr/bin/env python3
"""Follow band branches through crossings on a fine k line by linear extrapolation matching,
then plot the traced branches (colour = branch, not band index).
usage: trace_branches.py OUT.png "title" YLO YHI DIR LABEL [DIR LABEL ...]"""
import sys, glob, os
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
from scipy.optimize import linear_sum_assignment
plt.rcParams.update({'font.size':12})
def seg1(d):
    f=sorted(glob.glob(os.path.join(d,'bnd0*.spin1')))[0]
    b=defaultdict(list)
    for line in open(f):
        if line.startswith('#') or not line.strip(): continue
        t=line.split(); b[int(t[0])].append((float(t[1]),float(t[2])))
    return b
def trace(b, lo, hi):
    ibs=[ib for ib in range(lo,hi+1) if ib in b]
    X=np.array([x for x,_ in b[ibs[0]]]); X=X/X.max()
    E=np.array([[b[ib][i][1] for ib in ibs] for i in range(len(X))])  # (nk, nb), sorted by index
    n=E.shape[1]
    tr=np.zeros_like(E); tr[0]=E[0]; tr[1]=E[1]
    order=np.arange(n)
    for i in range(2,len(X)):
        pred=2*tr[i-1]-tr[i-2]                      # linear extrapolation of each branch
        C=np.abs(pred[:,None]-E[i][None,:])         # cost branch x eigenvalue
        r,c=linear_sum_assignment(C)
        tr[i]=E[i][c]
    return X, tr, ibs
out,title=sys.argv[1],sys.argv[2]; ylo,yhi=float(sys.argv[3]),float(sys.argv[4])
runs=[(sys.argv[i],sys.argv[i+1]) for i in range(5,len(sys.argv),2)]
fig,axs=plt.subplots(1,len(runs),figsize=(7.6*len(runs),8.5),squeeze=False)
cmap=plt.get_cmap('tab10')
for c,(d,lab) in enumerate(runs):
    ax=axs[0][c]; b=seg1(d)
    X,tr,ibs=trace(b,33,48)
    for j in range(tr.shape[1]):
        E=tr[:,j]
        if E.max()<ylo or E.min()>yhi: continue
        col=cmap(j%10)
        ax.plot(X,E,'-',lw=1.8,color=col)
        k=int(0.5*(len(X)-1))
        if ylo<E[k]<yhi: ax.annotate(f'br{j}',(X[k],E[k]),textcoords='offset points',xytext=(0,5),
                                     color=col,fontsize=10,fontweight='bold',ha='center')
    ax.set_xlim(0,1); ax.set_ylim(ylo,yhi); ax.set_xlabel('Γ → X  (k / k_X)'); ax.set_title(lab)
    ax.grid(color='0.92',lw=0.5)
    if c==0: ax.set_ylabel('E − E_F (eV)')
fig.suptitle(title,fontsize=13); fig.tight_layout(rect=(0,0,1,0.96)); fig.savefig(out,dpi=110)
# report
for d,lab in runs:
    X,tr,ibs=trace(seg1(d),33,48)
    print(f'--- {lab}')
    for j in range(tr.shape[1]):
        E=tr[:,j]
        m=(X>0.05)&(X<0.95)
        if E[m].mean()<0.78 or E[m].mean()>1.05: continue
        d1=np.diff(E); sg=np.sign(d1); idx=np.where(np.diff(sg)!=0)[0]+1
        xs=[X[i] for i in idx if 0.05<X[i]<0.95]
        amp=(E[m].max()-E[m].min())*1e3
        per=2*np.mean(np.diff(xs)) if len(xs)>1 else float('nan')
        print(f'  branch {j}: mean {E[m].mean():.3f} eV  peak-to-peak {amp:5.1f} meV  extrema {len(xs)}  period {per:.3f}')
