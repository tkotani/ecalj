import numpy as np, sys, os
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':12})
cases=[c for c in sys.argv[1:] if os.path.isdir(c)]
qlab={1:'iq=1 Gamma cell (offset)',2:'iq=2 (-1/6,1/6,1/6) 1st shell',3:'iq=3 (-1/3,1/3,1/3) 2nd shell',5:'iq=5 (0,0,1/3) 1st shell'}
fig,axs=plt.subplots(3,4,figsize=(22,13))
for col,iq in enumerate([1,2,3,5]):
    for c,ls in zip(cases,['-','--',':']):
        h=np.loadtxt(f'{c}/wc_head_iq{iq}.dat'); w=h[:,1]
        m=(w>0.2)&(w<6)
        axs[0,col].plot(w[m],h[m,2],ls,label=f'{c} Re',lw=1.6)
        axs[1,col].plot(w[m],h[m,3],ls,label=f'{c} Im',lw=1.6)
        d=np.loadtxt(f'{c}/wc_diag_iq{iq}.dat'); wd=d[:,0]; md=(wd>0.2)&(wd<6)
        axs[2,col].plot(wd[md],-d[md,1:].sum(axis=1),ls,label=f'{c}',lw=1.6)
    axs[0,col].set_title(qlab[iq]+'\nRe W_c(1,1)'); axs[1,col].set_title('Im W_c(1,1)'); axs[2,col].set_title('-Im tr W_c (diag 1..30) = plasmon spectrum')
    for r in range(3): axs[r,col].axhline(0,color='0.5',lw=.5); axs[r,col].set_xlabel('omega (eV)')
axs[0,0].legend(fontsize=9); axs[2,0].legend(fontsize=9)
fig.suptitle('LiTi2O4 6^3, chi0 1000 K, W_c(q, omega) head from --dumpW (2026-09-21 23:40): '+' vs '.join(cases),fontsize=14)
fig.tight_layout(rect=(0,0,1,0.95)); fig.savefig('wc_head_compare.png',dpi=80)
