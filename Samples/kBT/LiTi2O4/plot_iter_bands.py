#!/usr/bin/env python3
"""Bands of successive iterations on top of each other near E_F, and how rough they are (2026-09-28).

  plot_iter_bands.py <N> <out.png> <title> <label1>=<band file 1> <label2>=<band file 2> ...

N is the Sigma mesh (6 or 9): the dotted lines mark its q points on Gamma-X (x = 2n/N).  Band files are the 'i x E-EF' blocks of
the MLO band (bnd_mlo_*.dat) or of job_band (bndPMT_*.dat).  Upper panel: the t2g bands near E_F; lower panel: the two lowest
t2g bands (occupied, below E_F) on the first half of Gamma-X, where 9^3 shows a wiggle and 6^3 is smooth (user, 2026-09-28).
One thin line per iteration, light to dark.  Printed per iteration (bands sorted per k first): the roughness of the two
lowest t2g bands on x < 0.5 (mean |E(i-1) - 2E(i) + E(i+1)|, meV), the smallest gap between them on 0.02 < x < 0.45
(near 0 = they touch or cross), how far the lowest band sinks below its Gamma value on x < 0.3 (0 when it rises
monotonically, as on 6^3), and the roughness of all t2g bands 33-44.
"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

N, out, title = int(sys.argv[1]), sys.argv[2], sys.argv[3]
runs = [a.split('=', 1) for a in sys.argv[4:]]


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


mesh = [2 * n / N for n in range(N) if 2 * n / N <= 1 + 1e-9]
if abs(mesh[-1] - 1) > 1e-9 and N % 2 == 0:
    mesh.append(1.0)


def metrics(x, E):
    B = np.sort(E[32:44], axis=0)
    rough_all = np.mean([np.abs(b[:-2] - 2 * b[1:-1] + b[2:]).mean() for b in B]) * 1e3
    h = x < 0.5
    rough_occ = np.mean([np.abs(b[h][:-2] - 2 * b[h][1:-1] + b[h][2:]).mean() for b in B[:2]]) * 1e3
    g = (x > 0.02) & (x < 0.45)
    gap = (B[1][g] - B[0][g]).min() * 1e3
    dip = (B[0][0] - B[0][x < 0.3].min()) * 1e3     # how far the lowest band sinks below its Gamma value (0 = monotonic)
    return rough_occ, gap, dip, rough_all, B


fig, ax = plt.subplots(2, 1, figsize=(9, 11))
cmap = plt.get_cmap('viridis')
print('iteration  occupied-2 roughness x<0.5 (meV)  min gap 0.02<x<0.45 (meV)  dip below Gamma (meV)  t2g roughness (meV)')
for i, (lab, f) in enumerate(runs):
    x, E = bands(f)
    r, gmin, dip, ra, B = metrics(x, E)
    print(f'{lab:>9s}  {r:10.3f}  {gmin:10.1f}  {dip:10.1f}  {ra:10.3f}')
    col = cmap(0.15 + 0.8 * i / max(1, len(runs) - 1))
    for a, lo, hi, bb in ((ax[0], -0.8, 1.3, B), (ax[1], -0.6, 0.0, B[:2])):
        first = True
        for b in bb:
            a.plot(x, b, color=col, lw=0.6, label=f'{lab}: occupied rough {r:.3f} meV, dip {dip:.1f} meV' if first else None)
            first = False
for a, lo, hi in ((ax[0], -0.8, 1.3), (ax[1], -0.6, 0.0)):
    for q in mesh[1:-1] if abs(mesh[-1] - 1) < 1e-9 else mesh[1:]:
        a.axvline(q, color='0.5', ls=':', lw=0.7)
    a.axhline(0, color='0.6', lw=0.4)
    a.set_xlim(0, 1); a.set_ylim(lo, hi)
    a.set_xticks([0, 1]); a.set_xticklabels([r'$\Gamma$', 'X']); a.set_ylabel(r'$E-E_F$ (eV)')
ax[0].legend(fontsize=7, loc='lower right')
ax[0].set_title(title, fontsize=10)
ax[1].set_xlim(0, 0.5)
ax[1].set_xticks([0, 0.5]); ax[1].set_xticklabels([r'$\Gamma$', 'x = 0.5'])
ax[1].set_title(r'the two lowest t2g bands (occupied) near $\Gamma$', fontsize=9)
fig.tight_layout()
fig.savefig(out, dpi=130)
print('->', out)
