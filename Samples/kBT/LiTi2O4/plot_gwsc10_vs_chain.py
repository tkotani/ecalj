#!/usr/bin/env python3
"""Final MLO band of a gwsc 10 run (QMLO_* code) on top of the gwsc 1 x 10 chain it reproduces.

  plot_gwsc10_vs_chain.py <new run dir> <old chain dir> <out.png> [title]

Left: both bands along the path (old thick grey, new thin red), E - E_F in eV.  Right: new - old per band
index (meV) for the bands within +-3 eV; spikes at band crossings are index swaps, not shifts.  Band files: new bnd_mlo_final.dat, old mloband/bnd_iter10.dat
(the 'i x E-EF' blocks written by draw_mloband.sh).
"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

new, old, out = sys.argv[1], sys.argv[2], sys.argv[3]
title = sys.argv[4] if len(sys.argv) > 4 else ''

def bands(f):
    blocks, cur = [], []
    for l in open(f):
        if l.startswith('#'): continue
        w = l.split()
        if not w:
            if cur: blocks.append(cur); cur = []
            continue
        cur.append((float(w[1]), float(w[2])))
    if cur: blocks.append(cur)
    nk = min(len(b) for b in blocks)
    x = np.array([p[0] for p in blocks[0][:nk]])
    return x, np.array([[p[1] for p in b[:nk]] for b in blocks])

xn, en = bands(f'{new}/bnd_mlo_final.dat')
xo, eo = bands(f'{old}/mloband/bnd_iter10.dat')
nb, nk = min(len(en), len(eo)), min(en.shape[1], eo.shape[1])
en, eo, x = en[:nb, :nk], eo[:nb, :nk], xn[:nk]
win = 3.0
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 5.2), gridspec_kw={'width_ratios': [1.3, 1]})
first = True
for b in range(nb):
    if np.all(np.abs(eo[b]) > win + 1): continue
    a1.plot(x, eo[b], color='0.55', lw=2.2, ls='-', label='gwsc 1 x 10 (old names)' if first else None)
    a1.plot(x, en[b], color='C3', lw=0.9, label='gwsc 10 (QMLO_*)' if first else None)
    first = False
a1.axhline(0, color='k', lw=0.5)
a1.set_xlim(x[0], x[-1]); a1.set_ylim(-win, win)
a1.set_xticks([x[0], x[-1]]); a1.set_xticklabels([r'$\Gamma$', 'X'])
a1.set_ylabel(r'$E-E_F$ (eV)'); a1.set_title('MLO band after iteration 10')
h, l = a1.get_legend_handles_labels(); a1.legend(h[:2], l[:2], loc='lower right', fontsize=8)
d = (en - eo) * 1000
for b in range(nb):
    m = np.abs(eo[b]) < win
    if m.any(): a2.plot(x[m], d[b][m], lw=0.8)
mask = np.abs(eo) < win
rms, mx = np.sqrt(np.mean(d[mask]**2)), np.max(np.abs(d[mask]))
a2.axhline(0, color='k', lw=0.5)
a2.set_xlim(x[0], x[-1]); a2.set_xticks([x[0], x[-1]]); a2.set_xticklabels([r'$\Gamma$', 'X'])
a2.set_ylabel('new - old (meV)'); a2.set_title(f'bands within $\\pm${win:g} eV: rms {rms:.1f}, max {mx:.1f} meV')
if title: fig.suptitle(title)
fig.tight_layout(); fig.savefig(out, dpi=130)
print(f'rms {rms:.2f} meV  max {mx:.2f} meV  -> {out}')
