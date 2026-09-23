#!/usr/bin/env python3
"""Full band structure along the whole symmetry line: LDA vs 6^3 vs 9^3."""
import glob, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from collections import defaultdict
plt.rcParams.update({'font.size':12})
def read(d):
    segs=[]
    for f in sorted(glob.glob(d+'/bnd0*.spin1')):
        b=defaultdict(list)
        for line in open(f):
            if line.startswith('#') or not line.strip(): continue
            t=line.split(); b[int(t[0])].append((float(t[1]),float(t[2])))
        segs.append(b)
    return segs
labels=['Γ','X','U|K','Γ','L','W','X']
runs=(('LDA','Pchain/lda'),('6^3 P iter 8','Pchain/iter8'),('9^3 P iter 15','P999/iter15'))
fig,axs=plt.subplots(2,3,figsize=(21,11))
for c,(lab,d) in enumerate(runs):
    segs=read(d)
    for r,(lo,hi,ylo,yhi,tag) in enumerate(((5,33,-11,-3,'O 2p'),(20,60,-1.0,1.6,'near E_F'))):
        ax=axs[r,c]; ticks=[0]
        for b in segs:
            for ib in range(lo,hi):
                if ib in b:
                    xs=[x for x,_ in b[ib]]; ys=[e for _,e in b[ib]]
                    cc=('C3' if 9<=ib<=12 else '0.3') if tag=='O 2p' else ('C2' if ib<=52 else '0.6')
                    ax.plot(xs,ys,'-',lw=1.0,color=cc)
            ticks.append(max(x for x,_ in b[1]))
        for t in ticks: ax.axvline(t,color='0.75',lw=0.6)
        if tag!='O 2p': ax.axhline(0,color='0.4',lw=0.8,ls='--')
        ax.set_xticks(ticks); ax.set_xticklabels(labels[:len(ticks)])
        ax.set_xlim(ticks[0],ticks[-1]); ax.set_ylim(ylo,yhi)
        ax.set_title(f'{lab}   ({tag})')
        if c==0: ax.set_ylabel('E − E_F (eV)')
fig.suptitle('LiTi2O4 band structure: LDA vs 6^3 QSGW (converged, iter 8) vs 9^3 QSGW (iter 15), P setting',fontsize=14)
fig.tight_layout(rect=(0,0,1,0.965)); fig.savefig('full_bands.png',dpi=100)
