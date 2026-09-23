#!/usr/bin/env python3
"""br5 along Gamma-X with the sigma q-mesh points marked: the oscillation lives only BETWEEN mesh points."""
import glob, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from collections import defaultdict
from scipy.optimize import linear_sum_assignment
plt.rcParams.update({'font.size':12})
def seg1(d):
    f=sorted(glob.glob(d+'/bnd0*.spin1'))[0]
    b=defaultdict(list)
    for line in open(f):
        if line.startswith('#') or not line.strip(): continue
        t=line.split(); b[int(t[0])].append((float(t[1]),float(t[2])))
    return b
def trace(b,lo=33,hi=48):
    ibs=[ib for ib in range(lo,hi+1) if ib in b]
    X=np.array([x for x,_ in b[ibs[0]]]); X=X/X.max()
    E=np.array([[b[ib][i][1] for ib in ibs] for i in range(len(X))])
    tr=np.zeros_like(E); tr[0]=E[0]; tr[1]=E[1]
    for i in range(2,len(X)):
        pred=2*tr[i-1]-tr[i-2]
        r,c=linear_sum_assignment(np.abs(pred[:,None]-E[i][None,:])); tr[i]=E[i][c]
    k0=int(0.10*(len(X)-1)); return X, tr[:,np.argsort(tr[k0])]
fig,ax=plt.subplots(figsize=(11,7))
for lab,d,N,col,off in (('LDA (br5, shifted by -0.53 eV)','GXfineLDA',0,'0.45',-0.53),
                        ('6^3 P iter 8 (br5)','GXfine666',6,'C0',0.0),
                        ('9^3 P iter 15 (br5)','GXfine999',9,'C3',0.0)):
    X,tr=trace(seg1(d)); y=tr[:,5]+off
    ax.plot(X,y,'-',lw=2,color=col,label=lab)
    if N:
        mesh=[2*m/N for m in range(0,N//2+1) if 2*m/N<=1.0]
        idx=[int(round(t*(len(X)-1))) for t in mesh]
        ax.plot(X[idx],y[idx],'o',ms=11,mfc='none',mew=2.5,color=col)
        ax.plot(X[idx],y[idx],':',lw=1.6,color=col,alpha=0.8)
ax.set_xlabel('Γ → X  (k / k_X)'); ax.set_ylabel('E − E_F (eV)')
ax.set_title('br5 along Γ–X: open circles = points of the Σ q-mesh (x = 2m/N), dotted = line through them')
ax.grid(color='0.92'); ax.legend(loc='upper left')
fig.tight_layout(); fig.savefig('br5_mesh.png',dpi=110)
