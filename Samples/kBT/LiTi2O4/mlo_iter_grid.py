#!/usr/bin/env python3
"""LiTi2O4 6^3 MLO-QSGW: every iteration 1..10 of the t2g bands along Gamma-X.
usage: mlo_iter_grid.py DIR OUT.png"""
import sys, os, numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.sans-serif'] = ['IPAGothic', 'Droid Sans Fallback', 'DejaVu Sans']
matplotlib.rcParams['axes.unicode_minus'] = False
from matplotlib.ticker import MultipleLocator

D, OUT = sys.argv[1], sys.argv[2]
def rd(f):
    xs, es, cx, ce = [], [], [], []
    for ln in open(f):
        if ln.startswith('#'): continue
        t = ln.split()
        if len(t) < 3:
            if cx: xs.append(cx); es.append(ce); cx = []; ce = []
            continue
        cx.append(float(t[1])); ce.append(float(t[2]))
    if cx: xs.append(cx); es.append(ce)
    return np.array(xs[0]), np.array(es)

MESH = [0.0, 1/3, 2/3, 1.0]; IDX = [0, 70, 140, 210]; LO, HI = 32, 44
def chord(y, x):
    return np.array([np.abs(y[a:b+1] - np.interp(x[a:b+1], [x[a], x[b]], [y[a], y[b]])).max()
                     for a, b in zip(IDX[:-1], IDX[1:])]) * 1e3

cases = [('LDA', f'{D}/bnd_lda.dat')] + [(f'iter {i}', f'{D}/bnd_iter{i}.dat') for i in range(1, 11)]
cases = [(l, f) for l, f in cases if os.path.exists(f)]
cmap = plt.get_cmap('viridis')

fig, AX = plt.subplots(3, 4, figsize=(16, 10.5), sharex=True, sharey=True)
for k, (lab, f) in enumerate(cases):
    ax = AX.flat[k]; x, E = rd(f)
    col = '0.35' if lab == 'LDA' else cmap(0.08 + 0.85 * (k - 1) / max(1, len(cases) - 2))
    for ib in range(LO, HI): ax.plot(x, E[ib], '-', lw=1.1, color=col)
    t = np.array([chord(E[ib], x) for ib in range(LO, HI)])
    w = np.mean([E[ib].max() - E[ib].min() for ib in range(LO, HI)])
    ax.set_ylim(-0.9, 1.6); ax.set_xlim(0, 1); ax.axhline(0, color='k', lw=0.9)
    for q in MESH: ax.axvline(q, color='k', ls=':', lw=0.9)
    ax.yaxis.set_major_locator(MultipleLocator(0.5)); ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax.grid(axis='y', which='major', color='0.78', lw=0.6)
    ax.grid(axis='y', which='minor', color='0.92', lw=0.4)
    ax.set_title(f'{lab}   弦ずれ {t.mean():.0f} meV / 帯幅 {w*1e3:.0f} meV = {t.mean()/(w*1e3)*100:.1f} %',
                 fontsize=10)
for k in range(len(cases), AX.size): AX.flat[k].axis('off')
for ax in AX[:, 0]: ax.set_ylabel('$E-E_F$ [eV]')
for ax in AX[-1, :]: ax.set_xlabel('$\\Gamma \\to X$')
fig.suptitle('LiTi$_2$O$_4$ 6$^3$ MLO-QSGW (76 orb, pwmode=11): t$_{2g}$ バンドの反復変化  '
             '— 点線 = $\\Sigma$ の q メッシュ点', fontsize=13)
plt.tight_layout(rect=[0, 0, 1, 0.96]); plt.savefig(OUT, dpi=120)
print('wrote', OUT)
