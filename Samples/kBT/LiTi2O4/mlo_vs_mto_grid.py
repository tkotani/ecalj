#!/usr/bin/env python3
"""LiTi2O4 6^3: MLO-QSGW vs conventional (MTO) QSGW, t2g bands along Gamma-X.
Top row  : one-shot comparison (LDA / MTO iter1 / MLO 76 orb iter1 / MLO 154 orb iter1)
Mid row  : MLO 76-orbital chain, iterations 2 / 4 / 7 / 10
Bottom   : chord deviation (departure from a straight line between Sigma mesh points)
usage: mlo_vs_mto_grid.py ROOT OUT.png"""
import sys, os, numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.sans-serif'] = ['IPAGothic', 'Droid Sans Fallback', 'DejaVu Sans']
matplotlib.rcParams['axes.unicode_minus'] = False
from matplotlib.ticker import MultipleLocator
import matplotlib.gridspec as gs

ROOT, OUT = sys.argv[1], sys.argv[2]

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

MESH = [0.0, 1/3, 2/3, 1.0]; IDX = [0, 70, 140, 210]
LO, HI = 32, 44                      # t2g manifold (0-based band index)

def chord(y, x):
    o = []
    for a, b in zip(IDX[:-1], IDX[1:]):
        c = np.interp(x[a:b+1], [x[a], x[b]], [y[a], y[b]])
        o.append(np.abs(y[a:b+1] - c).max())
    return np.array(o) * 1e3          # meV

def dev(f):
    x, E = rd(f)
    t = np.array([chord(E[ib], x) for ib in range(LO, HI)])
    return x, E, t.mean(), t.max()

TOP = [('LDA',                 f'{ROOT}/liti_mlo/bnd_lda.dat',      '0.35'),
       ('従来 QSGW (MTO) iter 1', f'{ROOT}/liti_ref/bnd_iter1.dat',    'tab:blue'),
       ('MLO-QSGW 76 orb iter 1', f'{ROOT}/liti_mlo/bnd_iter1.dat',    'tab:red'),
       ('MLO-QSGW 154 orb iter 1',f'{ROOT}/liti_mlo154/bnd_iter1.dat', 'tab:purple')]
MID = [(f'MLO 76 orb iter {i}', f'{ROOT}/liti_mlo/bnd_iter{i}.dat', 'tab:red') for i in (2, 4, 7, 10)]

fig = plt.figure(figsize=(15.5, 11.0))
G = gs.GridSpec(3, 4, height_ratios=[1, 1, 0.75], hspace=0.42, wspace=0.20)

def panel(ax, lab, f, col):
    x, E, m, mx = dev(f)
    for ib in range(LO, HI):
        ax.plot(x, E[ib], '-', lw=1.1, color=col, alpha=0.85)
    ax.set_ylim(-0.9, 1.5); ax.set_xlim(0, 1); ax.axhline(0, color='k', lw=0.9)
    for q in MESH: ax.axvline(q, color='k', ls=':', lw=0.9)
    ax.yaxis.set_major_locator(MultipleLocator(0.5)); ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax.grid(axis='y', which='major', color='0.78', lw=0.6)
    ax.grid(axis='y', which='minor', color='0.92', lw=0.4)
    ax.set_title(f'{lab}\nchord dev  mean {m:.1f} / max {mx:.1f} meV', fontsize=9.5)
    return m, mx

for c, (lab, f, col) in enumerate(TOP):
    if os.path.exists(f): panel(fig.add_subplot(G[0, c]), lab, f, col)
for c, (lab, f, col) in enumerate(MID):
    if os.path.exists(f): panel(fig.add_subplot(G[1, c]), lab, f, col)

# bottom: chord deviation vs iteration
ax = fig.add_subplot(G[2, :])
it, mm, xx = [], [], []
for i in range(1, 11):
    f = f'{ROOT}/liti_mlo/bnd_iter{i}.dat'
    if os.path.exists(f):
        _, _, m, mx = dev(f); it.append(i); mm.append(m); xx.append(mx)
ax.plot(it, mm, 'o-', color='tab:red', lw=2, ms=6, label='MLO 76 orb  mean')
ax.plot(it, xx, 's--', color='tab:red', lw=1.4, ms=5, alpha=0.6, label='MLO 76 orb  max')
for lab, f, col, mk in [('従来 (MTO) iter 1', f'{ROOT}/liti_ref/bnd_iter1.dat', 'tab:blue', 'o'),
                        ('MLO 154 orb iter 1', f'{ROOT}/liti_mlo154/bnd_iter1.dat', 'tab:purple', 'o')]:
    if os.path.exists(f):
        _, _, m, mx = dev(f)
        ax.plot([1], [m], mk, color=col, ms=9, label=f'{lab}  mean {m:.1f}')
        ax.plot([1], [mx], 's', color=col, ms=7, alpha=0.55)
_, _, m0, _ = dev(f'{ROOT}/liti_mlo/bnd_lda.dat')
ax.axhline(m0, color='0.45', ls='-.', lw=1.2, label=f'LDA (内挿なし) mean {m0:.1f} meV')
ax.set_xlabel('QSGW iteration'); ax.set_ylabel('chord deviation [meV]')
ax.set_xticks(range(1, 11)); ax.grid(color='0.85', lw=0.6); ax.legend(fontsize=9, ncol=3)
ax.set_title('$\\Sigma$ メッシュ点の間での折れ（直線からのずれ）: MLO 版は反復とともに増大する', fontsize=10.5)

fig.suptitle('LiTi$_2$O$_4$  6$^3$  $\\Gamma\\to X$  t$_{2g}$ バンド  —  点線 = $\\Sigma$ の q メッシュ点 '
             '(pwmode=11, nkabc=6$^3$ = GW メッシュ)', fontsize=12.5)
plt.savefig(OUT, dpi=125, bbox_inches='tight')
print('wrote', OUT)
