#!/usr/bin/env python3
"""Where in k does each MLO-QSGW iteration change the bands?
Plots dE(k) between successive iterations for the t2g manifold.
usage: mlo_iter_diff.py DIR OUT.png"""
import sys, os, numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.sans-serif'] = ['IPAGothic', 'Droid Sans Fallback', 'DejaVu Sans']
matplotlib.rcParams['axes.unicode_minus'] = False
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
def load(tag):
    f = f'{D}/bnd_lda.dat' if tag == 'LDA' else f'{D}/bnd_iter{tag}.dat'
    x, E = rd(f); return x, np.sort(E[LO:HI], axis=0)   # sort within the manifold at each k
pairs = [('LDA', 1), (1, 2), (2, 3), (3, 4), (4, 5), (9, 10)]
fig, AX = plt.subplots(2, 3, figsize=(15, 7.6), sharex=True)
for k, (a, b) in enumerate(pairs):
    ax = AX.flat[k]
    x, Ea = load(a); _, Eb = load(b)
    d = (Eb - Ea) * 1e3
    for r in range(d.shape[0]): ax.plot(x, d[r], '-', lw=1.0)
    # how much of the change sits at the mesh points vs between them
    at = np.abs(d[:, IDX]).mean()
    mid = np.abs(d[:, [35, 105, 175]]).mean()
    ax.set_xlim(0, 1); ax.axhline(0, color='k', lw=0.9)
    for q in MESH: ax.axvline(q, color='k', ls=':', lw=0.9)
    ax.grid(color='0.88', lw=0.5)
    ax.set_title(f'$\\Delta E$ = iter {b} − {a}\nメッシュ点上 {at:.0f} meV / 中点 {mid:.0f} meV'
                 f'  (比 {mid/max(at,1e-9):.2f})', fontsize=10)
    if k % 3 == 0: ax.set_ylabel('$\\Delta E$ [meV]')
    if k >= 3: ax.set_xlabel('$\\Gamma \\to X$')
fig.suptitle('LiTi$_2$O$_4$ 6$^3$ MLO-QSGW: 反復ごとの t$_{2g}$ バンド変化  '
             '— 点線 = $\\Sigma$ の q メッシュ点（そこでは内挿は厳密）', fontsize=12.5)
plt.tight_layout(rect=[0, 0, 1, 0.94]); plt.savefig(OUT, dpi=125)
print('wrote', OUT)
