#!/usr/bin/env python3
"""Runs side by side, not on top of each other (2026-09-28, user: overlays are hard to read).

  plot_bands_side.py <out.png> <title> <N>:<label>:<band file> [...]

One column per run (N = its Sigma mesh; dotted lines at x = 2n/N).  Upper row: the t2g bands near E_F (-0.8..1.3 eV) on Gamma-X;
lower row: the two lowest t2g bands (occupied) on the first half, where 9^3 wiggles.  Black, 0.7 pt.  The lower panels give the
dip of the lowest band below its Gamma value (x < 0.3) and the splitting of the pair at the 9^3 mesh point 2/9 (meV).
Band files: the 'i x E-EF' blocks that draw_mloband.sh writes.
"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

out, title = sys.argv[1], sys.argv[2]
runs = [a.split(':', 2) for a in sys.argv[3:]]


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


fig, axs = plt.subplots(2, len(runs), figsize=(5.2 * len(runs), 9.5), squeeze=False, sharey='row')
for c, (N, lab, f) in enumerate(runs):
    N = int(N)
    x, E = bands(f)
    B = np.sort(E[32:44], axis=0)
    dip = (B[0][0] - B[0][x < 0.3].min()) * 1e3
    s29 = (np.interp(2 / 9, x, B[1]) - np.interp(2 / 9, x, B[0])) * 1e3
    for r, (lo, hi, xmax, bb) in enumerate(((-0.8, 1.3, 1.0, B), (-0.56, -0.18, 0.5, B[:2]))):
        a = axs[r][c]
        for b in bb:
            a.plot(x, b, color='k', lw=0.7)
        for n in range(1, N):
            if 2 * n / N < xmax:
                a.axvline(2 * n / N, color='0.6', ls=':', lw=0.7)
        a.axhline(0, color='0.6', lw=0.4)
        a.set_xlim(0, xmax); a.set_ylim(lo, hi)
        a.set_xticks([0, xmax]); a.set_xticklabels([r'$\Gamma$', 'X' if xmax == 1 else 'x = 0.5'])
        if c == 0:
            a.set_ylabel(r'$E-E_F$ (eV)')
    axs[0][c].set_title(lab, fontsize=10)
    axs[1][c].set_title(f'occupied pair: dip {dip:.1f} meV, split at 2/9 {s29:.1f} meV', fontsize=9)
fig.suptitle(title, fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig(out, dpi=110)
print('->', out)
