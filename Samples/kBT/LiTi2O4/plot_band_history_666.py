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
panels=[('june',it) for it in [1,2,3,4,6,7,8,10,16,17,18,20,25,30]]
fig,axs=plt.subplots(4,7,figsize=(28,18))
for i,(src,it) in enumerate(panels):
    row=(i//7)*2; col=i%7
    segs=read(f'{src}/iter{it}')
    for ax,lo,hi,ylo,yhi,tag in ((axs[row,col],7,20,-10.5,-4,'O 2p'),(axs[row+1,col],24,60,-2.5,2.5,'near E_F')):
        ticks=[0]
        for b in segs:
            for ib in range(lo,hi):
                if ib in b:
                    xs=[x for x,_ in b[ib]]; ys=[e for _,e in b[ib]]
                    c=('C3' if 9<=ib<=12 else 'C0') if tag=='O 2p' else ('C2' if ib<=52 else '0.5')
                    ax.plot(xs,ys,'-',lw=1.0,color=c)
            ticks.append(max(x for x,_ in b[1]))
        for t in ticks: ax.axvline(t,color='0.75',lw=0.5)
        if tag!='O 2p': ax.axhline(0,color='0.4',lw=0.6,ls='--')
        ax.set_xticks(ticks); ax.set_xticklabels(labels[:len(ticks)],fontsize=8)
        ax.set_xlim(ticks[0],ticks[-1]); ax.set_ylim(ylo,yhi)
        ax.set_title(('June ' if src=='june' else 'TODAY ')+f'iter {it}  ({tag})',fontsize=10,color='k' if src=='june' else 'C2')
for a in axs[:,0]: a.set_ylabel('E − E_F (eV)')
fig.suptitle('LiTi2O4 6³ dq=0.1, 1000 K (tetrakbt + t_sigmakbt, June): O 2p bottom (bands 9–12 red) and near-E_F Ti t2g (green, bands 25–52) per QSGW iteration',fontsize=13)
fig.tight_layout(rect=(0,0,1,0.97)); fig.savefig('bands_history_666_T1000.png',dpi=80)
