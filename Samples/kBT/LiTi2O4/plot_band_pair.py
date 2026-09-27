#!/usr/bin/env python3
"""Two runs on top of each other along the path: MLO band, conventional sigm band, and new - old (2026-09-28).

  plot_band_pair.py <new MLO band> <new sigm band> <old MLO band> <old sigm band> <out.png> [title [label new [label old]]]

Band files are the 'i x E-EF' blocks (one block per band) that run_gwsc10.sh (bnd_mlo_final.dat, bndPMT_final.dat) and the
older chains (mloband/bnd_iter10.dat, bndPMT_iter10.dat) write.  Left: MLO band, middle: sigm band (old solid black, new dashed
red, both 0.6 pt), with a lower row zoomed on the t2g bands near E_F; right: new - old per band index (meV) for the bands within
+-3 eV (spikes at crossings are index swaps).
The MLO band files have 1e-5 Ry (0.136 meV) steps, hence the staircase.
"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fmn, fsn, fmo, fso, out = sys.argv[1:6]
title = sys.argv[6] if len(sys.argv) > 6 else ''
labn = sys.argv[7] if len(sys.argv) > 7 else 'new'
labo = sys.argv[8] if len(sys.argv) > 8 else 'old'
win = 3.0


def bands(f):
    blocks, cur = [], []
    for l in open(f):
        if l.startswith('#'):
            continue
        w = l.split()
        if not w:
            if cur:
                blocks.append(cur)
                cur = []
            continue
        cur.append((float(w[1]), float(w[2])))
    if cur:
        blocks.append(cur)
    nk = min(len(b) for b in blocks)
    return np.array([p[0] for p in blocks[0][:nk]]), np.array([[p[1] for p in b[:nk]] for b in blocks])


fig, axs = plt.subplots(2, 3, figsize=(15, 10), gridspec_kw={'width_ratios': [1, 1, 1.1]})
ax = axs[0]
for a, az, fn, fo, name, col in [(axs[0][0], axs[1][0], fmn, fmo, 'MLO band', 'C0'),
                                 (axs[0][1], axs[1][1], fsn, fso, 'sigm band', 'C1')]:
    x, en = bands(fn)
    _, eo = bands(fo)
    nb, nk = min(len(en), len(eo)), min(en.shape[1], eo.shape[1])
    en, eo, x = en[:nb, :nk], eo[:nb, :nk], x[:nk]
    for aa, ylim in ((a, (-9, 7)), (az, (-0.8, 1.3))):   # full range, and the t2g bands near E_F
        first = True
        for b in range(nb):
            if np.all((eo[b] < ylim[0] - 1) | (eo[b] > ylim[1] + 1)):
                continue
            aa.plot(x, eo[b], color='k', lw=0.6, label=labo if first else None)
            aa.plot(x, en[b], color='C3', lw=0.6, ls=(0, (4, 3)), label=labn if first else None)
            first = False
        aa.axhline(0, color='0.5', lw=0.4)
        aa.set_xlim(x[0], x[-1]); aa.set_ylim(*ylim)
        aa.set_xticks([x[0], x[-1]]); aa.set_xticklabels([r'$\Gamma$', 'X']); aa.set_ylabel(r'$E-E_F$ (eV)')
        aa.legend(loc='lower right', fontsize=8)
    az.set_title(f'{name} near $E_F$ (t2g)', fontsize=9)
    m = np.abs(eo) < win
    d = (en - eo) * 1000
    a.set_title(f'{name}: rms {np.sqrt(np.mean(d[m]**2)):.1f}, max {np.max(np.abs(d[m])):.1f} meV (|E|<{win:g} eV)', fontsize=9)
    first = True
    for b in range(nb):
        mb = np.abs(eo[b]) < win
        if mb.any():
            ax[2].plot(x[mb], d[b][mb], lw=0.7, color=col, label=name if first else None)
            first = False
ax[2].axhline(0, color='k', lw=0.5)
ax[2].set_xlim(x[0], x[-1]); ax[2].set_xticks([x[0], x[-1]]); ax[2].set_xticklabels([r'$\Gamma$', 'X'])
ax[2].set_ylabel(f'{labn} - {labo} (meV)'); ax[2].legend(fontsize=8)
ax[2].set_title(f'bands within $\\pm${win:g} eV (spikes at crossings = index swaps)', fontsize=9)
axs[1][2].axis('off')
axs[1][2].text(0.02, 0.95, f'black solid: {labo}\nred dashed: {labn}\nlines 0.6 pt; lower row: the t2g bands near $E_F$',
               va='top', fontsize=10, transform=axs[1][2].transAxes)
if title:
    fig.suptitle(title, fontsize=11)
fig.tight_layout()
fig.savefig(out, dpi=130)
print('->', out)
