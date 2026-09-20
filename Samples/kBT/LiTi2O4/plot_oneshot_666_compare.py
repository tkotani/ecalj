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
runs=[('REF','REF: no wcsmear (09-18)','0.45'),('OMP01','yesterday: wcsmear (09-19; pole FD, imag-axis Gaussian)','C0'),('NEW20260920','today: wcsmear (09-20; FD on both contour parts)','C3')]
fig,axs=plt.subplots(2,4,figsize=(22,10))
for col,(d,title,c) in enumerate(runs):
    segs=read(d)
    for row,(lo,hi,ylo,yhi,tag) in enumerate(((7,20,-10.5,-4,'O 2p'),(24,60,-2.5,2.5,'near E_F'))):
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
        ax.set_title(f'{title}  ({tag})',fontsize=10)
# overlay column
for row,(lo,hi,ylo,yhi,tag) in enumerate(((7,20,-10.5,-4,'O 2p'),(24,60,-2.5,2.5,'near E_F'))):
    ax=axs[row,3]; ticks=[0]
    for d,title,c in runs:
        segs=read(d); ticks=[0]
        for b in segs:
            for ib in range(lo,hi):
                if ib in b:
                    xs=[x for x,_ in b[ib]]; ys=[e for _,e in b[ib]]
                    ax.plot(xs,ys,'-',lw=0.9 if d!='REF' else 0.7,color=c,alpha=0.9 if d!='REF' else 0.7,label=title if ib==lo else None)
            ticks.append(max(x for x,_ in b[1]))
    for t in ticks: ax.axvline(t,color='0.75',lw=0.5)
    if tag!='O 2p': ax.axhline(0,color='0.4',lw=0.6,ls='--')
    ax.set_xticks(ticks); ax.set_xticklabels(labels[:len(ticks)],fontsize=9)
    ax.set_xlim(ticks[0],ticks[-1]); ax.set_ylim(ylo,yhi)
    ax.set_title(f'overlay: gray REF / blue yesterday / red today  ({tag})',fontsize=10)
for a in axs[:,0]: a.set_ylabel('E − E_F (eV)')
fig.suptitle('LiTi2O4 6^3 dq=0.1, 1000 K, one-shot QSGW from LDA: O 2p bottom (bands 9-12 red) and near-E_F t2g (bands 25-52 green). Today\'s consistent code (NEW) reproduces yesterday\'s wcsmear run; the REF needle on Gamma-L is absent in both',fontsize=12)
fig.tight_layout(rect=(0,0,1,0.96)); fig.savefig('oneshot_666_REF_OMP01_NEW_20260920.png',dpi=80)
