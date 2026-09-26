import numpy as np, datetime, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.optimize import linear_sum_assignment
NK, NB = 211, 88
cols = ['s','py','pz','px','dxy','dyz','dz2','dxz','dx2y2']
NB = 60   # only the lowest 60 bands are needed (t2g are 33..44); the count per k varies above that
def load(a):
    blocks=[]; cur=[]; lastx=None
    for l in open(f'bandweight_atom{a}.spin1'):
        if l.startswith('#'): continue
        t=l.split()
        if len(t)<=3: continue
        xv=float(t[0])
        if lastx is not None and xv!=lastx: blocks.append(cur); cur=[]
        cur.append([float(v) for v in t[:11]]); lastx=xv
    blocks.append(cur)
    assert len(blocks)==NK, (a, len(blocks))
    arr=np.array([b[:NB] for b in blocks])                 # (NK, NB, 11)
    return arr[:,:,0], arr[:,:,1], arr[:,:,2:]
x, E, _ = load(3); x=x[:,0]
W={a:load(a)[2] for a in range(1,15)}
Li=[1,2]; Ti=[3,4,5,6]; O=list(range(7,15))
def ch(k,b):   # character of sorted band b (0-based among all 88) at k
    g=lambda at,names: sum(W[a][k,b,cols.index(n)] for a in at for n in names)
    return dict(Ti_xy=g(Ti,['dxy']), Ti_yz=g(Ti,['dyz']), Ti_zx=g(Ti,['dxz']), Ti_eg=g(Ti,['dz2','dx2y2']),
                O_p=g(O,['px','py','pz']), Li=g(Li,cols), rest=g(Ti,['s','px','py','pz'])+g(O,['s']))
# t2g block, bands 33..44 -> sorted indices 32..43
S=E[:,32:44].T                                  # (12, NK)
T=np.empty_like(S); idx=np.empty((12,NK),int); T[:,0]=S[:,0]; T[:,1]=S[:,1]; idx[:,0]=idx[:,1]=np.arange(12)
for i in range(2,NK):
    pred=2*T[:,i-1]-T[:,i-2]; r,c=linear_sum_assignment(np.abs(pred[:,None]-S[None,:,i])); T[r,i]=S[c,i]; idx[r,i]=c
mid=(x>0.40)&(x<0.70); best,bb=0.05,None
for b in range(12):
    bulge=T[b,mid].max()-max(np.interp(0.30,x,T[b]),np.interp(0.80,x,T[b]))
    if bulge>best: best,bb=bulge,b
print(f'blue = tracked t2g band {bb} (bulge {1000*best:.0f} meV)')
red_sorted=np.zeros(NK,int); blue_sorted=idx[bb]
keys=['Ti_xy','Ti_yz','Ti_zx','Ti_eg','O_p','Li','rest']
C={'red':np.array([[ch(k,32+red_sorted[k])[q] for q in keys] for k in range(NK)]),
   'blue':np.array([[ch(k,32+blue_sorted[k])[q] for q in keys] for k in range(NK)])}
for lab in ('red','blue'):
    print(f'\n{lab}: x    E[eV]  ' + '  '.join(f'{q:>6s}' for q in keys) + '   sum')
    for xi in (0.0,0.15,0.30,0.40,0.45,0.50,0.55,0.60,0.70,0.80,0.90,1.0):
        k=int(round(xi*210)); e=(T[0,k] if lab=='red' else T[bb,k]); c=C[lab][k]
        print(f'      {x[k]:.3f} {e:+.3f}  ' + '  '.join(f'{v:6.3f}' for v in c) + f'  {c.sum():.3f}')
# figure
fig,AX=plt.subplots(3,1,figsize=(8.2,10.5),sharex=True,gridspec_kw=dict(height_ratios=[1.3,1,1]))
ax=AX[0]
for b in range(12): ax.plot(x,T[b],color='0.75',lw=0.9)
ax.plot(x,S[0],color='tab:red',lw=2); ax.plot(x,T[bb],color='tab:blue',lw=2)
for q in [0,2/9,4/9,6/9,8/9]: ax.axvline(q,color='k',ls=':',lw=0.8)
ax.axhline(0,color='k',lw=0.8); ax.set_ylim(-0.9,1.5); ax.set_ylabel('$E-E_F$ [eV]')
ax.set_title('t$_{2g}$ bands, PMT band drawn via sigm (as in the left column)')
colr=dict(Ti_xy='#1f77b4',Ti_yz='#17becf',Ti_zx='#2ca02c',Ti_eg='#9467bd',O_p='#d62728',Li='#8c564b',rest='#bbbbbb')
lbl=dict(Ti_xy='Ti d$_{xy}$',Ti_yz='Ti d$_{yz}$',Ti_zx='Ti d$_{zx}$',Ti_eg='Ti e$_g$',O_p='O p',Li='Li',rest='Ti s,p + O s')
for ax,lab in ((AX[1],'red'),(AX[2],'blue')):
    ax.stackplot(x,*[C[lab][:,i] for i in range(len(keys))],colors=[colr[q] for q in keys],labels=[lbl[q] for q in keys])
    for q in [0,2/9,4/9,6/9,8/9]: ax.axvline(q,color='k',ls=':',lw=0.8)
    ax.set_ylim(0,1); ax.set_ylabel('weight'); ax.set_title(f'character of the {"lowest (red)" if lab=="red" else "bulging (blue)"} band')
AX[1].legend(loc='upper left',fontsize=8,ncol=4)
AX[2].set_xlabel('$\\Gamma\\to X$  (path along $-y$)')
fig.suptitle(f'LiTi$_2$O$_4$ 9$^3$ MLO-QSGW iter 4: orbital character (fat band, atom/lm-projected; the gap to 1 is the unprojected part)\ngenerated {datetime.datetime.now():%Y-%m-%d %H:%M}',fontsize=10)
plt.tight_layout(rect=[0,0,1,0.95]); plt.savefig('/home/takao/ecalj/Samples/kBT/LiTi2O4/fatchar_k9_iter4.png',dpi=115)
