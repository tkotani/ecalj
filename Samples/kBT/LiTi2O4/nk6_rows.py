#!/usr/bin/env python3
"""Row layout for the nk6 series: rows = LDA, iter 1, iter 2, ...  columns = O 2p | near E_F.
usage: nk6_rows.py OUT.png TITLE DIR "0,1,2,..."   (0 = LDA -> DIR/lda, n -> DIR/iterN)"""
import sys, glob, os
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from collections import defaultdict
plt.rcParams.update({'font.size':12})
def read(d):
    segs=[]
    for f in sorted(glob.glob(os.path.join(d,'bnd0*.spin1'))):
        b=defaultdict(list)
        for line in open(f):
            if line.startswith('#') or not line.strip(): continue
            t=line.split(); b[int(t[0])].append((float(t[1]),float(t[2])))
        segs.append(b)
    return segs
labels=['Γ','X','U|K','Γ','L','W','X']
out,title,d0=sys.argv[1],sys.argv[2],sys.argv[3]; its=[int(x) for x in sys.argv[4].split(',')]
rows=[('lda','LDA') if i==0 else (f'iter{i}',f'iter {i}') for i in its]
fig,axs=plt.subplots(len(rows),2,figsize=(15,4.6*len(rows)),squeeze=False)
for r,(sub,lab) in enumerate(rows):
    path=os.path.join(d0,sub)
    if not os.path.isdir(path):
        axs[r,0].set_visible(False); axs[r,1].set_visible(False); continue
    segs=read(path)
    for c,(lo,hi,ylo,yhi,tag) in enumerate(((5,33,-11,-3,'O 2p'),(20,60,-1.0,1.6,'near E_F'))):
        ax=axs[r,c]; ticks=[0]
        for b in segs:
            for ib in range(lo,hi):
                if ib in b:
                    xs=[x for x,_ in b[ib]]; ys=[e for _,e in b[ib]]
                    cc=('C3' if 9<=ib<=12 else '0.3') if tag=='O 2p' else ('C2' if ib<=52 else '0.6')
                    ax.plot(xs,ys,'-',lw=1.1,color=cc)
            ticks.append(max(x for x,_ in b[1]))
        for t in ticks: ax.axvline(t,color='0.75',lw=0.6)
        if tag!='O 2p': ax.axhline(0,color='0.4',lw=0.8,ls='--')
        ax.set_xticks(ticks); ax.set_xticklabels(labels[:len(ticks)])
        ax.set_xlim(ticks[0],ticks[-1]); ax.set_ylim(ylo,yhi); ax.set_title(f'{lab}   ({tag})')
    axs[r,0].set_ylabel('E − E_F (eV)')
fig.suptitle(title,fontsize=13); fig.tight_layout(rect=(0,0,1,1-0.03/len(rows)))
fig.savefig(out,dpi=100)
