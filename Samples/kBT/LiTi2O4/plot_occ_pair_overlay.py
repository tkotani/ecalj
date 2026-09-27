#!/usr/bin/env python3
"""The two lowest t2g bands (occupied) near Gamma of two runs on top of each other, one pair of runs per row (2026-09-28).

  plot_occ_pair_overlay.py <out.png> <title> '<row title>|<blue file>|<red file>|<blue label>|<red label>' [...]

For "does 9^3 match 6^3 below E_F after 10 iterations, and is tf32 the cause" (user, 2026-09-28): e.g. the 09-26 fp32 chains
(the runs of the ecaljdoc mlo_gwsc.md figure), the tf32 runs of now, and 6^3 tf32 against 6^3 fp32.  Blue solid, red dashed,
0.7 pt; bands sorted per k (33 and 34 of the MLO band files); dotted lines at the Sigma mesh points of 6^3 (1/3) and 9^3 (2/9, 4/9).
Each row says red - blue on the pair for x < 0.5 (rms and max, meV).
"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

out, title = sys.argv[1], sys.argv[2]
rows = [a.split('|') for a in sys.argv[3:]]


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


fig, axs = plt.subplots(len(rows), 1, figsize=(9, 3.6 * len(rows) + 0.5), squeeze=False)
for a, (rt, fb, fr, lb, lr) in zip(axs[:, 0], rows):
    xb, Bb = bands(fb)
    xr, Br = bands(fr)
    Bb, Br = np.sort(Bb[32:44], axis=0)[:2], np.sort(Br[32:44], axis=0)[:2]
    h = xr < 0.5
    d = np.array([Br[i][h] - np.interp(xr[h], xb, Bb[i]) for i in (0, 1)]) * 1e3
    for i in (0, 1):
        a.plot(xb, Bb[i], color='C0', lw=0.7, label=lb if i == 0 else None)
        a.plot(xr, Br[i], color='C3', lw=0.7, ls=(0, (4, 3)), label=lr if i == 0 else None)
    for q, s in ((1 / 3, '6^3'), (2 / 9, '9^3'), (4 / 9, '9^3')):
        a.axvline(q, color='0.6', ls=':', lw=0.7)
        a.text(q, 0.02, s, transform=a.get_xaxis_transform(), fontsize=7, color='0.4', ha='center')
    a.set_xlim(0, 0.5); a.set_ylim(-0.56, -0.18)
    a.set_xticks([0, 0.5]); a.set_xticklabels([r'$\Gamma$', 'x = 0.5 (half way to X)']); a.set_ylabel(r'$E-E_F$ (eV)')
    a.set_title(f'{rt}: red - blue on x < 0.5: rms {np.sqrt((d ** 2).mean()):.1f}, max {np.abs(d).max():.1f} meV', fontsize=9)
    a.legend(fontsize=8, loc='lower right')
fig.suptitle(title, fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 1 - 0.4 / (3.6 * len(rows) + 0.5)))
fig.savefig(out, dpi=110)
print('->', out)
