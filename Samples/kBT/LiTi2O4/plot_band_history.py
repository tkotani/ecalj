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
        segs.append((0.0,b))
    return segs
labels=['Γ','X','U|K','Γ','L','W','X']
panels=[('june',it) for it in [1,2,3,4,6,7,10,17,18,20,30,45]]+[('today',it) for it in [1,2]]
fig,axs=plt.subplots(2,7,figsize=(28,9),sharey=True)
for ax,(src,it) in zip(axs.flat,panels):
    segs=read(f'{src}/iter{it}'); ticks=[0]
    for off,b in segs:
        for ib in range(7,20):
            if ib in b:
                xs=[off+x for x,_ in b[ib]]; ys=[e for _,e in b[ib]]
                ax.plot(xs,ys,'-',lw=1.2,color='C3' if 9<=ib<=12 else 'C0')
        ticks.append(max(x for x,_ in b[1]))
    for t in ticks: ax.axvline(t,color='0.7',lw=0.5)
    ax.set_xticks(ticks); ax.set_xticklabels(labels[:len(ticks)],fontsize=8)
    ax.set_xlim(ticks[0],ticks[-1]); ax.set_ylim(-10.5,-4)
    ax.set_title(('June ' if src=='june' else 'TODAY ')+f'iter {it}',fontsize=10,color='k' if src=='june' else 'C2')
axs[0,0].set_ylabel('E − E_F (eV)'); axs[1,0].set_ylabel('E − E_F (eV)')
fig.suptitle('LiTi2O4 9³ dq=0.1, 1000 K (tetrakbt + t_sigmakbt): O 2p bottom (bands 9–12 red) per QSGW iteration',fontsize=12)
fig.tight_layout(rect=(0,0,1,0.96)); fig.savefig('bands_history_999_T1000.png',dpi=90)
