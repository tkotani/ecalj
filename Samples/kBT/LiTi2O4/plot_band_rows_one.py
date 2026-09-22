#!/usr/bin/env python3
"""Rows: iterations of ONE run. Columns: O 2p | near E_F.
usage: plot_rows9on.py OUT.png TITLE DIR "9,10,11,..."   (DIR/iterN)"""
import sys, glob, os
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':11})
from collections import defaultdict
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
out,title,d9=sys.argv[1],sys.argv[2],sys.argv[3]; its=[int(x) for x in sys.argv[4].split(',')]
rows=[(f'iter{i}',f'iter {i}') for i in its]
fig,axs=plt.subplots(len(rows),2,figsize=(18,5.2*len(rows)),squeeze=False)
def panel(ax,segs,lo,hi,ylo,yhi,tag,ttl):
    ticks=[0]
    for b in segs:
        for ib in range(lo,hi):
            if ib in b:
                xs=[x for x,_ in b[ib]]; ys=[e for _,e in b[ib]]
                cc=('C3' if 9<=ib<=12 else '0.3') if tag=='O 2p' else ('C2' if ib<=52 else '0.6')
                ax.plot(xs,ys,'-',lw=0.9,color=cc)
        ticks.append(max(x for x,_ in b[1]))
    for t in ticks: ax.axvline(t,color='0.75',lw=0.5)
    if tag!='O 2p': ax.axhline(0,color='0.4',lw=0.6,ls='--')
    ax.set_xticks(ticks); ax.set_xticklabels(labels[:len(ticks)]); ax.set_xlim(ticks[0],ticks[-1]); ax.set_ylim(ylo,yhi); ax.set_title(ttl,fontsize=11)
for r,(sub,lab) in enumerate(rows):
    path=os.path.join(d9,sub)
    if not os.path.isdir(path):
        axs[r,0].set_visible(False); axs[r,1].set_visible(False); continue
    segs=read(path)
    panel(axs[r,0],segs,5,22,-11,-3,'O 2p',f'9^3 {lab}  (O 2p)')
    panel(axs[r,1],segs,20,60,-2.5,3.0,'near E_F',f'9^3 {lab}  (near E_F)')
    axs[r,0].set_ylabel('E − E_F (eV)')
fig.suptitle(title,fontsize=13); fig.tight_layout(rect=(0,0,1,0.985)); fig.savefig(out,dpi=90)
