import os, numpy as np, datetime, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
def rd(f):
    xs, es, cx, ce = [], [], [], []
    for ln in open(f):
        if ln.startswith('#'): continue
        t = ln.split()
        if len(t) < 3:
            if cx: xs.append(cx); es.append(ce); cx=[]; ce=[]
            continue
        cx.append(float(t[1])); ce.append(float(t[2]))
    if cx: xs.append(cx); es.append(ce)
    return np.sort(np.array(es)[32:44], axis=0)
def chain(d):
    o={}
    for i in range(0,31):
        f=f'{d}/bnd_lda.dat' if i==0 else f'{d}/bnd_iter{i}.dat'
        if os.path.exists(f): o[i]=rd(f)
    return o
R=chain('ref'); M=chain('v6')
fig,ax=plt.subplots(figsize=(7.2,4.2))
for lab,C,col in [('conventional QSGW (MTO)',R,'tab:blue'),('MLO-QSGW (new scheme, unmixed)',M,'tab:green')]:
    ks=sorted(C); xs=[];ys=[]
    for a,b in zip(ks,ks[1:]):
        if b!=a+1: continue
        xs.append(b); ys.append(1000*np.abs(C[b]-C[a]).max())
    ax.plot(xs,ys,'o-',color=col,label=lab,lw=1.6,ms=6)
    for x,y in zip(xs,ys): ax.annotate(f'{y:.0f}',(x,y),textcoords='offset points',xytext=(0,7),ha='center',fontsize=8,color=col)
ax.set_yscale('log'); ax.set_xlabel('QSGW iteration N  (change from N-1 to N)')
ax.set_ylabel(r'max$_{k,b}\;|E_N-E_{N-1}|$  [meV]')
ax.set_title('LiTi$_2$O$_4$ 6$^3$: how much the t$_{2g}$ bands still move each iteration\n'
             'conventional $\\Sigma$ is Anderson-mixed ($\\beta$=0.5); $\\Sigma^{MLO}$ is not - and is kept that way,\n'
             'since the overshoot at 3 damps out by itself\n'
             f'generated {datetime.datetime.now():%Y-%m-%d %H:%M}',fontsize=10)
ax.grid(alpha=.3,which='both'); ax.legend(fontsize=9); ax.set_xticks(range(1,11))
plt.tight_layout(); plt.savefig('/home/takao/ecalj/Samples/kBT/LiTi2O4/mlo_conv.png',dpi=120)
print('ref :',{b:round(1000*np.abs(R[b]-R[b-1]).max(),1) for b in sorted(R) if b>0 and b-1 in R})
print('mlo :',{b:round(1000*np.abs(M[b]-M[b-1]).max(),1) for b in sorted(M) if b>0 and b-1 in M})
