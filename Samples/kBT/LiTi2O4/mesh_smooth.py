#!/usr/bin/env python3
"""Interpolated band curves (211 k points) vs the sigma q-mesh points, for several cases.
Rows = cases, columns = occupied t2g / unoccupied t2g-eg."""
import glob, sys, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from collections import defaultdict
plt.rcParams.update({'font.size':11.5})
def seg1(d):
    b=defaultdict(list)
    for line in open(sorted(glob.glob(d+'/bnd0*.spin1'))[0]):
        if line.startswith('#') or not line.strip(): continue
        t=line.split(); b[int(t[0])].append((float(t[1]),float(t[2])))
    return b
# (label, dir, N, line label)
cases=[('9^3 P iter 15  (09-23 10:20)','GXfine999',9,'Γ → X'),
       ('6^3 P iter 8   (09-22 22:51)','GXfine666',6,'Γ → X')]
import os
if os.path.isdir('gl_fine_999'): cases.append(('9^3 P iter 15  (09-23 10:20)','gl_fine_999',9,'Γ → L'))
if os.path.isdir('gl_fine_666'): cases.append(('6^3 P iter 8   (09-22 22:51)','gl_fine_666',6,'Γ → L'))
cmap=plt.get_cmap('tab10')
fig,axs=plt.subplots(len(cases),2,figsize=(15.5,5.6*len(cases)),squeeze=False)
for r,(lab,d,N,line) in enumerate(cases):
    b=seg1(d); X=np.array([x for x,_ in b[33]]); X=X/X.max()
    mesh=[2*m/N for m in range(N//2+1) if 2*m/N<=1.0]
    idx=[int(round(t*(len(X)-1))) for t in mesh]
    for c,(lo,hi,ylo,yhi,tt) in enumerate(((33,34,-0.62,0.32,'occupied t2g (33, 34)'),
                                           (35,44,0.28,1.15,'unoccupied t2g / eg (35-44)'))):
        ax=axs[r,c]
        for j,ib in enumerate(range(lo,hi+1)):
            if ib not in b: continue
            E=np.array([e for _,e in b[ib]]); col=cmap(j%10)
            ax.plot(X,E,'-',lw=1.5,color=col)
            ax.plot(X[idx],E[idx],'o',ms=9,mfc='none',mew=2.2,color=col)
            ax.plot(X[idx],E[idx],':',lw=1.3,color=col,alpha=0.85)
            if ylo<E[-1]<yhi: ax.annotate(str(ib),(1.0,E[-1]),textcoords='offset points',xytext=(5,-3),
                                          fontsize=9,color=col,annotation_clip=False)
        for t in mesh: ax.axvline(t,color='0.86',lw=1.0,zorder=0)
        ax.set_xlim(-0.01,1.05); ax.set_ylim(ylo,yhi); ax.set_xlabel(f'{line}   (k / k_end)')
        ax.set_title(f'{lab}   {line}   {tt}',fontsize=11)
        ax.grid(color='0.94')
        if c==0:
            ax.set_ylabel('E − E_F (eV)'); ax.axhline(0,color='0.4',lw=0.7,ls='--')
fig.suptitle('Interpolated bands (solid, 211 k points) vs the Σ q-mesh points (circles, x = 2m/N; dotted = line through them)',fontsize=13)
fig.tight_layout(rect=(0,0,1,1-0.025/len(cases)*2)); fig.savefig('mesh_smooth2.png',dpi=100)
# numbers
print('case                          band | mesh-pt residual about a cubic (max) | full-curve residual (rms)')
for lab,d,N,line in cases:
    b=seg1(d); X=np.array([x for x,_ in b[33]]); X=X/X.max()
    mesh=np.array([2*m/N for m in range(N//2+1) if 2*m/N<=1.0]); idx=[int(round(t*(len(X)-1))) for t in mesh]
    for ib in (33,34,41,44):
        if ib not in b: continue
        E=np.array([e for _,e in b[ib]]); ym=E[idx]
        deg=min(3,len(mesh)-2)
        rm=ym-np.polyval(np.polyfit(mesh,ym,deg),mesh)
        m=(X>0.02)&(X<0.98); rf=E[m]-np.polyval(np.polyfit(X[m],E[m],5),X[m])
        print(f'{lab[:28]:28s} {line} {ib:3d} |  {np.abs(rm).max()*1e3:7.2f} meV                     |  {rf.std()*1e3:7.2f} meV')
