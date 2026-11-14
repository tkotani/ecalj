#!/usr/bin/env python3
"""Schematic of the muffin-tin division of space (2026-11-13).

    python3 pmt_mtdivision.py [--en]      # writes pmt_mtdivision.png and pmt_mtdivision.dat next to this script

(a) a 2D cut of a rock-salt-like plane: spheres of two species, radii 0.97 x the touching radii (the default --tratio=0.97
of ctrlgenToml.py), and the interstitial region shaded. (b) a model potential along the diagonal through the large spheres (dashed in (a)):
-Z/r-like inside the spheres, smooth in the interstitial. Illustration only.
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle

here = os.path.dirname(os.path.abspath(__file__))
# --en draws the English version (<name>_en.png, for ecaljdoc/en/theory/pmt.md); the default is Japanese (2026-11-13).
EN = '--en' in sys.argv
T = lambda ja, en: en if EN else ja
SUF = '_en' if EN else ''
plt.rcParams['font.family'] = ['Noto Sans CJK JP', 'DejaVu Sans']
teal, ink, inter = '#1f6f78', '#333', '#dcebed'

d = 1.0                                  # nearest-neighbour distance A-B
ratio = 0.97                             # default touching ratio
RA, RB = 0.58 * d * ratio, 0.42 * d * ratio
fig, (ax, bx) = plt.subplots(1, 2, figsize=(10, 4.0), gridspec_kw={'width_ratios': [1, 1.25]})

# (a) MT spheres in a plane
n = 4
ax.add_patch(Rectangle((-0.5, -0.5), n, n, color=inter, lw=0))
rows = []
for i in range(n):
    for j in range(n):
        big = (i + j) % 2 == 0
        R = RA if big else RB
        ax.add_patch(Circle((i, j), R, facecolor='white', edgecolor=teal, lw=1.2))
        ax.plot(i, j, 'o', ms=3.5 if big else 2.5, color=ink)
        rows.append((i, j, R))
ax.annotate('', xy=(2 + RA, 2), xytext=(2, 2), arrowprops=dict(arrowstyle='->', color='#c0392b', lw=1.2))
ax.text(2 + RA / 2, 2.06, '$R_a$', color='#c0392b', ha='center', va='bottom', fontsize=10)
ax.text(1.5, 1.5, T('格子間', 'inter-\nstitial'), ha='center', va='center', fontsize=9 if not EN else 7, color=ink, linespacing=0.9)
ax.plot([-0.5, n - 0.5], [-0.5, n - 0.5], color='#c0392b', lw=1.0, ls='--')
ax.set_xlim(-0.5, n - 0.5); ax.set_ylim(-0.5, n - 0.5); ax.set_aspect('equal'); ax.axis('off')
ax.set_title(T('(a) MT 球（白）と格子間（色）', '(a) MT spheres (white) and interstitial (shaded)'), fontsize=11)

# (b) model potential along the diagonal (t, t) through the large spheres A at t = 0, 1, 2, 3 (spacing sqrt(2) d)
x = np.linspace(-0.5, n - 0.5, 3001) * np.sqrt(2)
V = np.zeros_like(x)
for i in range(-20, n + 20):           # images outside the window too, so that both ends look periodic (2026-11-14)
    V -= 0.35 / np.sqrt((x - i * np.sqrt(2))**2 + 0.004)
V -= V.max()
V = np.maximum(V, -6.0)
for i in range(n):
    c = i * np.sqrt(2)
    bx.axvspan(c - RA, c + RA, color='white', lw=0, zorder=0)
    bx.axvline(c - RA, color=teal, lw=0.8, ls=':'); bx.axvline(c + RA, color=teal, lw=0.8, ls=':')
bx.set_facecolor(inter)
bx.plot(x, V, color=ink, lw=1.6)
bx.text(1.5 * np.sqrt(2), 0.35, T('格子間: なだらか', 'interstitial: smooth'), ha='center', va='bottom', fontsize=9)
bx.annotate(T('球の中: 深くほぼ球対称', 'in the sphere: deep, nearly spherical'), xy=(np.sqrt(2) + 0.04, -5.0), xytext=(np.sqrt(2) + 0.35, -4.2), fontsize=9, color=teal,
            ha='left', va='center', bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='none'),
            arrowprops=dict(arrowstyle='->', color=teal, lw=1.0))
bx.set_xlim(x[0], x[-1]); bx.set_ylim(-6.2, 0.6)
bx.set_xticks([]); bx.set_yticks([]); bx.set_ylabel('$V$', rotation=0, labelpad=10)
for s in ('top', 'right'):
    bx.spines[s].set_visible(False)
bx.set_title(T('(b) (a) の破線上のポテンシャル（模式）', '(b) potential along the dashed line in (a) (schematic)'), fontsize=11)

fig.tight_layout()
fig.savefig(os.path.join(here, 'pmt_mtdivision' + SUF + '.png'), dpi=150)
np.savetxt(os.path.join(here, 'pmt_mtdivision.dat'), np.array(rows), fmt='%.4f',
           header='pmt_mtdivision.py (2026-11-13), schematic: x y R of the spheres in (a); A-B distance 1, R = 0.97 x touching')
