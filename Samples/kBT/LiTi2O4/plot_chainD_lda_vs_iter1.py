import glob, os
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
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
runs=[('lda','LDA (starting point)','0.35'),('iter1','QSGW iter 1: chi0 1000 K + chi0_filterw=[3,0.2] Drude filtered, t_sigmaw 1000, wcsmear','C3')]
fig,axs=plt.subplots(2,3,figsize=(20,10))
for col,(d,title,c) in enumerate(runs):
    segs=read(d)
    for row,(lo,hi,ylo,yhi,tag) in enumerate(((5,22,-11,-3,'O 2p'),(20,60,-2.5,3.0,'near E_F'))):
        ax=axs[row,col]; ticks=[0]
        for b in segs:
            for ib in range(lo,hi):
                if ib in b:
                    xs=[x for x,_ in b[ib]]; ys=[e for _,e in b[ib]]
                    cc=('C3' if 9<=ib<=12 else '0.3') if tag=='O 2p' else ('C2' if ib<=52 else '0.6')
                    ax.plot(xs,ys,'-',lw=1.0,color=cc)
            ticks.append(max(x for x,_ in b[1]))
        for t in ticks: ax.axvline(t,color='0.75',lw=0.5)
        if tag!='O 2p': ax.axhline(0,color='0.4',lw=0.6,ls='--')
        ax.set_xticks(ticks); ax.set_xticklabels(labels[:len(ticks)],fontsize=9)
        ax.set_xlim(ticks[0],ticks[-1]); ax.set_ylim(ylo,yhi)
        ax.set_title(f'{title}  ({tag})',fontsize=9)
for row,(lo,hi,ylo,yhi,tag) in enumerate(((5,22,-11,-3,'O 2p'),(20,60,-2.5,3.0,'near E_F'))):
    ax=axs[row,2]
    for d,title,c in runs:
        segs=read(d); ticks=[0]
        for b in segs:
            for ib in range(lo,hi):
                if ib in b:
                    xs=[x for x,_ in b[ib]]; ys=[e for _,e in b[ib]]
                    ax.plot(xs,ys,'-',lw=0.9,color=c,alpha=0.85)
            ticks.append(max(x for x,_ in b[1]))
    for t in ticks: ax.axvline(t,color='0.75',lw=0.5)
    if tag!='O 2p': ax.axhline(0,color='0.4',lw=0.6,ls='--')
    ax.set_xticks(ticks); ax.set_xticklabels(labels[:len(ticks)],fontsize=9)
    ax.set_xlim(ticks[0],ticks[-1]); ax.set_ylim(ylo,yhi)
    ax.set_title(f'overlay: gray LDA / red QSGW iter 1  ({tag})',fontsize=10)
for a in axs[:,0]: a.set_ylabel('E − E_F (eV)')
fig.suptitle('LiTi2O4 6^3: LDA vs one QSGW iteration with all smoothing (chain D, 2026-09-21 10:27). O 2p bands 9-12 red; near-E_F t2g bands 21-52 green',fontsize=12)
fig.tight_layout(rect=(0,0,1,0.96)); fig.savefig('chainD_lda_vs_iter1.png',dpi=75)
