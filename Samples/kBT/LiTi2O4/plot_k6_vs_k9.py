#!/usr/bin/env python3
"""6^3 and 9^3 bands on top of each other near E_F, with the Sigma mesh points marked (2026-09-28).

  plot_k6_vs_k9.py <6^3 MLO band> <6^3 sigm band> <9^3 MLO band> <9^3 sigm band> <out.png> [title [label 6^3 [label 9^3]]]

Columns: MLO band, sigm band; rows: all bands from -9 to 7 eV, the t2g bands near E_F, and the O 2p bands
(ROWS='lo:hi:label,...' in the environment chooses other windows, e.g. ROWS='-0.7:4.6:Ti 3d (t2g + eg)').  Thin lines (6^3 blue, 9^3 red).  Dotted
vertical lines are the q points of the Sigma mesh on Gamma-X: x = 2n/N (6^3: 0, 1/3, 2/3, 1; 9^3: 0, 2/9, 4/9, 6/9, 8/9).
The roughness printed in each title is the mean |E(i-1) - 2E(i) + E(i+1)| over the 211 path points of the bands in the panel (meV).
"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

f6m, f6s, f9m, f9s, out = sys.argv[1:6]
title = sys.argv[6] if len(sys.argv) > 6 else ''
lab6 = sys.argv[7] if len(sys.argv) > 7 else '6$^3$'
lab9 = sys.argv[8] if len(sys.argv) > 8 else '9$^3$'


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


def rough(E, lo, hi):
    sel = [b for b in E if b.max() > lo and b.min() < hi]
    return np.mean([np.abs(b[:-2] - 2 * b[1:-1] + b[2:]).mean() for b in sel]) * 1e3 if sel else float('nan')


import os
ROWS = [(-9.0, 7.0, 'all'), (-0.8, 1.3, 't2g near $E_F$'), (-8.9, -3.9, 'O 2p')]
if os.environ.get('ROWS'):      # e.g. ROWS='-0.7:4.6:Ti 3d (t2g + eg)'
    ROWS = [(float(a), float(b), c) for a, b, c in (r.split(':', 2) for r in os.environ['ROWS'].split(','))]
fig, axs = plt.subplots(len(ROWS), 2, figsize=(13, 5.3 * len(ROWS) + 0.5), squeeze=False)
for col, (f6, f9, name) in enumerate([(f6m, f9m, 'MLO band'), (f6s, f9s, 'sigm band')]):
    x6, e6 = bands(f6)
    x9, e9 = bands(f9)
    for row, (lo, hi, lab) in enumerate(ROWS):
        a = axs[row][col]
        for x, E, c, l in ((x6, e6, 'C0', lab6), (x9, e9, 'C3', lab9)):
            first = True
            for b in E:
                if b.max() < lo - 0.5 or b.min() > hi + 0.5:
                    continue
                a.plot(x, b, color=c, lw=0.6, label=l if first else None)
                first = False
        for n in range(1, 3):
            a.axvline(2 * n / 6, color='C0', ls=':', lw=0.7)
        for n in range(1, 5):
            a.axvline(2 * n / 9, color='C3', ls=':', lw=0.7)
        a.axhline(0, color='0.6', lw=0.4)
        a.set_xlim(0, 1); a.set_ylim(lo, hi)
        a.set_xticks([0, 1]); a.set_xticklabels([r'$\Gamma$', 'X']); a.set_ylabel(r'$E-E_F$ (eV)')
        a.set_title(f'{name}, {lab}: roughness 6$^3$ {rough(e6, lo, hi):.3f}, 9$^3$ {rough(e9, lo, hi):.3f} meV', fontsize=9)
        a.legend(loc='lower right', fontsize=8)
if title:
    fig.suptitle(title, fontsize=11)
fig.tight_layout()
fig.savefig(out, dpi=130)
print('->', out)
