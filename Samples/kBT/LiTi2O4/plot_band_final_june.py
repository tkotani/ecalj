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
panels=[('n666_dq0.1_T1000K_sigmakbt1000','6³ 1000 K, χ₀+Σ, iter 30 (final)'),
        ('n666_dq0.1_T2000K_sigmakbt2000','6³ 2000 K, χ₀+Σ, iter 11 (final)'),
        ('n666_dq0.1_T3000K_sigmakbt3000','6³ 3000 K, χ₀+Σ, iter 29 (final)'),
        ('n666_dq0.1_T3000K_baseline','6³ 3000 K, χ₀ only, iter 60 (final)'),
        ('n666_dq0.1_T3000K_esmr019','6³ 3000 K, χ₀ only, esmr=0.019, iter 15'),
        ('n666_dq0.1_T5000K_scf10','6³ 5000 K, χ₀ only, mixbeta 1, iter 10'),
        ('n999_dq0.1_T1000K_sigmakbt1000','9³ 1000 K, χ₀+Σ, iter 45 (final)'),
        ('n999_dq0.1_T3000K_sigmakbt3000','9³ 3000 K, χ₀+Σ, iter 2 (last done)'),
        ('n666_dq0.1_T1000K_scf10','6³ 1000 K, χ₀ only, iter 1 (one-shot)'),
        ('n999_dq0.1_T3000K_scf10','9³ 3000 K, χ₀ only, iter 1 (one-shot)')]
fig,axs=plt.subplots(2,5,figsize=(26,11),sharey=True)
for ax,(d,title) in zip(axs.flat,panels):
    segs=read(d); ticks=[0]
    for b in segs:
        for ib in range(7,60):
            if ib in b:
                xs=[x for x,_ in b[ib]]; ys=[e for _,e in b[ib]]
                ax.plot(xs,ys,'-',lw=1.0,color='C3' if 9<=ib<=12 else ('C0' if ib<25 else '0.4'))
        ticks.append(max(x for x,_ in b[1]))
    for t in ticks: ax.axvline(t,color='0.75',lw=0.5)
    ax.axhline(0,color='0.5',lw=0.6,ls='--')
    ax.set_xticks(ticks); ax.set_xticklabels(labels[:len(ticks)],fontsize=9)
    ax.set_xlim(ticks[0],ticks[-1]); ax.set_ylim(-10.5,3.5)
    ax.set_title(title,fontsize=11)
for a in axs[:,0]: a.set_ylabel('E − E_F (eV)')
fig.suptitle('LiTi2O4 dq=0.1: final-iteration QSGW bands of the June finite-T runs (red: O 2p bottom, bands 9–12; blue: 13–24; grey: Ti t2g/eg)',fontsize=13)
fig.tight_layout(rect=(0,0,1,0.96)); fig.savefig('bands_final_june_runs.png',dpi=90)
