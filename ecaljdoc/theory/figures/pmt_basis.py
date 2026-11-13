#!/usr/bin/env python3
"""Schematic of the PMT basis (APW and MTO) in the 3-component formalism, a 1D cut through two MT spheres (2026-11-13).

    python3 pmt_basis.py [--en]      # writes pmt_basis.png and pmt_basis.dat next to this script

Illustration only, not from an ecalj run: the 0th component is a plane wave (APW) or a smooth Hankel h0 centred on atom A (MTO);
the 2nd component is set equal to the 0th within the MT (condition A'); the 1st component is the 2nd plus a bump that has zero
value and slope at the MT boundary (condition B) and a sharp peak at the nucleus, standing in for the true partial waves.
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.special import erfc

here = os.path.dirname(os.path.abspath(__file__))
# --en draws the English version (<name>_en.png, for ecaljdoc/en/theory/pmt.md); the default is Japanese (2026-11-13).
EN = '--en' in sys.argv
T = lambda ja, en: en if EN else ja
SUF = '_en' if EN else ''
plt.rcParams['font.family'] = ['Noto Sans CJK JP', 'DejaVu Sans']

xa, xb, R = -2.0, 2.0, 1.0           # atoms A, B and the MT radius
x = np.linspace(-4.2, 4.2, 2001)


def smhankel(r, rsm=0.5, kappa=1.0):
    """Smooth Hankel h0 (Eq. 25 of the page): finite at r = 0, -> exp(-kappa r)/r for r >> rsm."""
    r = np.maximum(r, 1e-6)
    a = rsm * kappa / 2
    return (np.exp(-kappa * r) * erfc(a - r / rsm) - np.exp(kappa * r) * erfc(a + r / rsm)) / (2 * r)


def bump(x, xc, amp):
    """Zero value and slope at |x - xc| = R; peaked at the nucleus like a true partial wave."""
    d = np.abs(x - xc) / R
    w = np.where(d < 1, (1 - d**2)**2, 0.0)
    return amp * w * (np.exp(-d / 0.12) - 0.35 * np.cos(2.5 * np.pi * d))


cases = {
    'APW': np.cos(1.4 * x + 0.4),
    'MTO': smhankel(np.abs(x - xa)) / smhankel(np.array([0.0]))[0],
}
amps = {'APW': (0.9, -0.45), 'MTO': (1.1, 0.25)}

fig, axs = plt.subplots(1, 2, figsize=(10, 3.6), sharey=True)
cols = [x]
for ax, (name, f0) in zip(axs, cases.items()):
    f1 = f0 + bump(x, xa, amps[name][0]) + bump(x, xb, amps[name][1])
    inmt = (np.abs(x - xa) < R) | (np.abs(x - xb) < R)
    aug = np.where(inmt, f1, f0)                      # A[F] = F0 + F1 - F2, with F2 = F0 inside the MT
    for xc, lab in ((xa, 'A'), (xb, 'B')):
        ax.axvspan(xc - R, xc + R, color='#1f6f78', alpha=0.10, lw=0)
        ax.text(xc, -1.78, f'MT {lab}', ha='center', va='bottom', fontsize=9, color='#1f6f78')
        ax.plot([xc], [-1.52], marker='o', ms=5, color='#1f6f78')  # nucleus
    ax.plot(x, f0, color='#c4c9cc', lw=6, solid_capstyle='butt', label=T('第 0 成分 $F_0$（envelope、全空間）', '0th component $F_0$ (envelope, all space)'))
    ax.plot(x, np.where(inmt, f0, np.nan), color='#2e8b57', lw=1.4, label=T('第 2 成分 $F_2$（MT 内、$F_0$ に重なる）', '2nd component $F_2$ (in MT, on top of $F_0$)'))
    ax.plot(x, np.where(inmt, f1, np.nan), color='#c0392b', lw=1.8, label=T('第 1 成分 $F_1$（MT 内の本物）', '1st component $F_1$ (the true function in MT)'))
    ax.set_title({'APW': T('(a) APW: 平面波の envelope', '(a) APW: plane-wave envelope'), 'MTO': T('(b) MTO: smooth Hankel の envelope（A 中心）', '(b) MTO: smooth-Hankel envelope (centred on A)')}[name], fontsize=11)
    ax.axhline(0, color='#bbb', lw=0.6)
    ax.set_xlim(x[0], x[-1]); ax.set_ylim(-1.8, 1.75)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ('top', 'right', 'left'):
        ax.spines[s].set_visible(False)
    cols += [f0, f1, aug]
axs[0].legend(loc='upper left', fontsize=8.5, frameon=False, ncol=1)
fig.tight_layout()
fig.savefig(os.path.join(here, 'pmt_basis' + SUF + '.png'), dpi=150)
np.savetxt(os.path.join(here, 'pmt_basis.dat'), np.array(cols).T, fmt='%.6f',
           header='pmt_basis.py (2026-11-13), schematic, MT A at x=-2, B at x=2, R=1\n'
                  'x  APW:F0 F1 A[F]  MTO:F0 F1 A[F]   (F2 = F0 inside the MT)')
