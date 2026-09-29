#!/usr/bin/env python3
'''Figure of the mass fit: |E(k)-E(Gamma)| against k for the three lines, with the fitted curves.
usage: python3 massplot.py [title]     (reads massfit.npz written by massfit.py, writes massfit.png)
The plotted numbers are those of massfit.npz (k_<line>_<ib>, e_<line>_<ib>, ab_<line>_<ib>).
'''
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

d = np.load('massfit.npz')
emin, emax = d['window']
lines = ['111', '100', '110']
bands = {'so': 1, 'lh': 3, 'hh': 5, 'e': 7}          # first band of each pair
ymax = 0.1
labelat = {'e': (0.092, 'right'), 'lh': (0.080, 'left'), 'so': (0.062, 'left'), 'hh': (0.090, 'right')}   # energy (eV), side
fig, axes = plt.subplots(1, 3, figsize=(10, 3.6), sharey=True)
for ax, ln in zip(axes, lines):
    for name, ib in bands.items():
        k, e = d[f'k_{ln}_{ib}'], d[f'e_{ln}_{ib}']
        a, b = d[f'ab_{ln}_{ib}']
        sel = e <= ymax
        ax.plot(k[sel], e[sel], 'o', ms=2.5, mfc='none', mec='0.35', mew=0.5)
        ee = np.linspace(0, ymax, 200)
        k2 = a*ee + b*ee**2
        ok = k2 >= 0
        ax.plot(np.sqrt(k2[ok]), ee[ok], '-', color='k', lw=0.6)
        elab, side = labelat[name]                       # direct label on the fitted curve
        klab = np.sqrt(max(a*elab + b*elab**2, 0.))
        ax.annotate(name, (klab, elab), xytext=(-4 if side == 'right' else 4, 0), textcoords='offset points',
                    fontsize=8, va='center', ha=side)
    ax.axhspan(emin, emax, color='0.92', zorder=0)       # fitting window
    ax.set_xlabel('k (1/bohr)')
    ax.set_title(f'[{ln}]', fontsize=9)
    ax.set_ylim(0, ymax)
    ax.set_xlim(left=0)
    ax.tick_params(labelsize=8)
axes[0].set_ylabel('|E(k) - E(Gamma)| (eV)')
title = sys.argv[1] if len(sys.argv) > 1 else 'GaAs, bands near Gamma with SOC'
fig.suptitle(f'{title}: circles = lmf, lines = fit of eq. (1), gray = fitting window ({d["made"]})', fontsize=9)
fig.tight_layout()
fig.savefig('massfit.png', dpi=150)
print('massfit.png written; data: massfit.npz from', d['source'])
