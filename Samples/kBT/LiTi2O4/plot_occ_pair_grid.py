#!/usr/bin/env python3
"""The two lowest t2g bands (occupied) near Gamma, one small panel per iteration (2026-09-28).

  plot_occ_pair_grid.py <out.png> <title> <N>:<label>:<file pattern with {i}>:<i1>-<i2> [...]

Each run is one row (or two, if more than 5 iterations): N = Sigma mesh (dotted lines at x = 2n/N), the file pattern gives the
MLO band of iteration i (e.g. liti_mlo_k9/mloband/bnd_iter{i}.dat).  Panel title: the dip of the lowest band below its Gamma value
on x < 0.3 (meV; 0 = rises monotonically) and the roughness of the pair on x < 0.5 (mean |second difference|, meV).
"""
import sys, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

out, title = sys.argv[1], sys.argv[2]
specs = []
for a in sys.argv[3:]:
    N, lab, pat, rng = a.split(':')
    i1, i2 = map(int, rng.split('-'))
    specs.append((int(N), lab, pat, list(range(i1, i2 + 1))))


def bands(f):
    blocks, cur = [], []
    for l in open(f):
        if l.startswith('#'):
            continue
        w = l.split()
        if len(w) < 3:
            if cur:
                blocks.append(cur)
                cur = []
            continue
        cur.append((float(w[1]), float(w[2])))
    if cur:
        blocks.append(cur)
    nk = min(len(b) for b in blocks)
    return np.array([p[0] for p in blocks[0][:nk]]), np.array([[p[1] for p in b[:nk]] for b in blocks])


rows = []
for N, lab, pat, its in specs:
    for k in range(0, len(its), 5):
        rows.append((N, lab, pat, its[k:k + 5]))
fig, axs = plt.subplots(len(rows), 5, figsize=(16, 3.1 * len(rows)), squeeze=False)
for r, (N, lab, pat, its) in enumerate(rows):
    for c in range(5):
        a = axs[r][c]
        if c >= len(its):
            a.axis('off')
            continue
        i = its[c]
        f = pat.format(i=i)
        if not os.path.exists(f):
            a.axis('off')
            continue
        x, E = bands(f)
        B = np.sort(E[32:44], axis=0)
        dip = (B[0][0] - B[0][x < 0.3].min()) * 1e3
        h = x < 0.5
        rough = np.mean([np.abs(b[h][:-2] - 2 * b[h][1:-1] + b[h][2:]).mean() for b in B[:2]]) * 1e3
        for b, col in zip(B[:2], ('C0', 'C3')):
            a.plot(x, b, color=col, lw=0.8)
        for n in range(1, N):
            if 2 * n / N < 0.5:
                a.axvline(2 * n / N, color='0.5', ls=':', lw=0.7)
        a.set_xlim(0, 0.5); a.set_ylim(-0.58, -0.15)
        a.set_xticks([0, 0.5]); a.set_xticklabels([r'$\Gamma$', '0.5'])
        a.set_title(f'{lab} iter {i}: dip {dip:.1f}, rough {rough:.3f} meV', fontsize=8)
        if c == 0:
            a.set_ylabel(r'$E-E_F$ (eV)')
fig.suptitle(title, fontsize=11)
fig.tight_layout()
fig.savefig(out, dpi=110)
print('->', out)
