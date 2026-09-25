#!/usr/bin/env python3
"""LiTi2O4 6^3: is either QSGW chain converging?
Left  = max over (k,band): sensitive to a single crossing, jumpy.
Right = rms over (k,band): what says whether the whole band is still moving.
usage: mlo_conv.py [ROOT]   (ROOT holds ref/ and v6/)"""
import os, sys, datetime, numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT = sys.argv[1] if len(sys.argv) > 1 else '.'
def rd(f):
    es, ce = [], []
    for ln in open(f):
        if ln.startswith('#'): continue
        t = ln.split()
        if len(t) < 3:
            if ce: es.append(ce); ce = []
            continue
        ce.append(float(t[2]))
    if ce: es.append(ce)
    return np.sort(np.array(es)[32:44], axis=0)
def chain(d):
    o = {}
    for i in range(0, 31):
        f = f'{d}/bnd_lda.dat' if i == 0 else f'{d}/bnd_iter{i}.dat'
        if os.path.exists(f): o[i] = rd(f)
    return o
C = [('conventional MTO, sigm mixed $\\beta$=0.5',        chain(f'{ROOT}/ref'), 'tab:blue'),
     ('MLO 126, $\\Sigma^{MLO}$ unmixed ($\\beta$=1)',      chain(f'{ROOT}/v7'),  'tab:red'),
     ('MLO 126, $\\Sigma^{MLO}$ mixed $\\beta$=0.5 - drawn via sigm', chain(f'{ROOT}/v9'),  'tab:purple'),
     ('MLO 126, $\\Sigma^{MLO}$ mixed $\\beta$=0.5 - MLO band', chain(f'{ROOT}/v9mlo'), 'tab:orange')]
C = [c for c in C if len(c[1]) > 1]
fig, AX = plt.subplots(1, 2, figsize=(12.4, 4.6), sharex=True)
for ax, how, lab in [(AX[0], 'max', r'max$_{k,b}\;|E_N-E_{N-1}|$  [meV]'),
                     (AX[1], 'rms', r'rms$_{k,b}\;|E_N-E_{N-1}|$  [meV]')]:
    for name, ch, col in C:
        ks = sorted(ch); xs, ys = [], []
        for a, b in zip(ks, ks[1:]):
            if b != a + 1: continue
            d = 1000 * (ch[b] - ch[a])
            xs.append(b); ys.append(np.abs(d).max() if how == 'max' else np.sqrt((d**2).mean()))
        ax.plot(xs, ys, 'o-', color=col, label=name, lw=1.6, ms=6)
        for x, y in zip(xs, ys):
            ax.annotate(f'{y:.0f}' if y >= 10 else f'{y:.1f}', (x, y), textcoords='offset points',
                        xytext=(0, 7), ha='center', fontsize=8, color=col)
    ax.set_yscale('log'); ax.set_xlabel('QSGW iteration N  (change from N-1 to N)')
    ax.set_ylabel(lab); ax.grid(alpha=.3, which='both'); ax.set_xticks(range(1, 11))
AX[0].legend(fontsize=9)
AX[0].set_title('max - one crossing can dominate it', fontsize=10)
AX[1].set_title('rms - whether the whole band is still moving', fontsize=10)
fig.suptitle('LiTi$_2$O$_4$ 6$^3$: how much the t$_{2g}$ bands still move each iteration\n'
             'the conventional chain flattens from N=6 on and stops falling - that is a floor, not convergence\n'
             'all three are nmlo 126 or conventional MTO - only the damping differs between the last two\n'
             f'generated {datetime.datetime.now():%Y-%m-%d %H:%M}', fontsize=10.5)
plt.tight_layout(rect=[0, 0, 1, 0.88])
plt.savefig('/home/takao/ecalj/Samples/kBT/LiTi2O4/mlo_conv.png', dpi=118)
print('wrote mlo_conv.png')
