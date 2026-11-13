#!/usr/bin/env python3
"""Schematic of the frozen core with the smooth-Hankel tail n^c_sH,a (2026-11-13).

    python3 pmt_core.py [--en]      # writes pmt_core.png and pmt_core.dat next to this script

Illustration only, not from an ecalj run. A model core density n^c(r) (two exponentials, like 1s and an outer shell)
is fitted for r > R by one smooth Hankel function h0 (Eq. 25 of the page, l=0), as lmfa and corprm do for n^c_sH,a. (a) radial plot on a log scale; (b) the three components of n^Zc (Eq. 14, without the
nucleus) along a line through atoms A and B: the 0th component is the sum of the n^c_sH of all sites, so inside MT A it
also holds the tail from B, while the 2nd component of A is only its own n^c_sH.
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
R = 2.0                                   # MT radius (a.u.)
rsm = 1.0                                 # smoothing radius of the smooth Hankel (rfoca)


def smhankel(r, rsm, kappa):
    """Smooth Hankel h0: finite at r = 0, -> exp(-kappa r)/r for r >> rsm."""
    r = np.maximum(r, 1e-6)
    a = rsm * kappa / 2
    return (np.exp(-kappa * r) * erfc(a - r / rsm) - np.exp(kappa * r) * erfc(a + r / rsm)) / (2 * r)


def ncore(r):
    """Model core density: a deep 1s-like part and an outer shell that spills out of the MT."""
    return 60.0 * np.exp(-12.0 * r) + 0.4 * (1 + r) * np.exp(-2.8 * r)


# fit n^c for R <= r <= 3 R by one smooth Hankel c*h0(kappa) (l=0), as corprm/lmfa do (one function, ctail and etail);
# kappa is scanned and c is the least-squares coefficient on the relative error
rf = np.linspace(R, 3.0 * R, 200)
w = 1 / ncore(rf)
best = None
for kappa in np.linspace(1.5, 4.0, 251):
    h = smhankel(rf, rsm, kappa) * w
    c = max(np.dot(h, ncore(rf) * w) / np.dot(h, h), 0.0)
    res = np.sum((c * h - ncore(rf) * w)**2)
    if best is None or res < best[0]:
        best = (res, kappa, c)
_, kappa, coef = best


def nsh(r):
    return coef * smhankel(r, rsm, kappa)


fig, (ax, bx) = plt.subplots(1, 2, figsize=(10, 4.4), gridspec_kw={'width_ratios': [1, 1.35]})

# (a) one site, radial
r = np.linspace(0.005, 4.5, 1500)
ax.axvspan(0, R, color='#1f6f78', alpha=0.10, lw=0)
ax.text(R / 2, 3e-4, T('MT の中', 'inside MT'), ha='center', fontsize=9, color='#1f6f78')
ax.axvline(R, color='#1f6f78', lw=0.8, ls=':')
ax.text(R + 0.05, 2e-4, '$R_a$', color='#1f6f78', fontsize=10)
ax.semilogy(r, ncore(r), color='#c0392b', lw=1.8, label=T('$n^{\\rm c}_a$（lmfa のコア密度）', '$n^{\\rm c}_a$ (core density from lmfa)'))
ax.semilogy(r, nsh(r), color='#2e8b57', lw=1.8, ls='--', label='$n^{\\rm c}_{{\\rm sH},a}$（smooth Hankel）')
ax.set_ylim(1e-4, 1e2); ax.set_xlim(0, 4.5)
ax.set_xlabel('$r$ (a.u.)'); ax.set_ylabel(T('密度（対数）', 'density (log)'))
ax.legend(loc='upper right', fontsize=8.5, frameon=False)
ax.set_title(T('(a) 1 サイト: MT の外で $n^{\\rm c}_a$ に合わせる', '(a) one site: fitted to $n^{\\rm c}_a$ outside the MT'), fontsize=11)

# (b) two sites along a line, the three components (nucleus omitted)
xa, xb = -2.6, 2.6
x = np.linspace(-5.6, 5.6, 3001)
ra, rb = np.abs(x - xa), np.abs(x - xb)
c0 = nsh(ra) + nsh(rb)                                   # 0th: sum over sites
inA, inB = ra < R, rb < R
c1 = np.where(inA, ncore(ra), np.where(inB, ncore(rb), np.nan))
c2 = np.where(inA, nsh(ra), np.where(inB, nsh(rb), np.nan))
for xc, lab in ((xa, 'A'), (xb, 'B')):
    bx.axvspan(xc - R, xc + R, color='#1f6f78', alpha=0.10, lw=0)
    bx.text(xc - R + 0.1, 2e-6, f'MT {lab}', ha='left', fontsize=9, color='#1f6f78')
bx.semilogy(x, c0, color='#c4c9cc', lw=6, solid_capstyle='butt', label=T('第 0 成分: $\\sum n^{\\rm c}_{\\rm sH}$（全サイト）', '0th component: $\\sum n^{\\rm c}_{\\rm sH}$ (all sites)'))
bx.semilogy(x, c2, color='#2e8b57', lw=1.4, label=T('第 2 成分: そのサイトの $n^{\\rm c}_{{\\rm sH},a}$', '2nd component: $n^{\\rm c}_{{\\rm sH},a}$ of the site'))
bx.semilogy(x, c1, color='#c0392b', lw=1.6, label=T('第 1 成分: $n^{\\rm c}_a$', '1st component: $n^{\\rm c}_a$'))
bx.semilogy(x, nsh(rb), color='#7a6fb0', lw=1.0, ls=':', label=T('B の裾 $n^{\\rm c}_{{\\rm sH},B}$', 'tail of B, $n^{\\rm c}_{{\\rm sH},B}$'))
bx.annotate(T('A の中の B の裾:\n第 0 成分にあり、第 2 成分にない', 'tail of B inside A:\nin the 0th component, not in the 2nd'), xy=(xa + 0.9, nsh(np.array([xb - xa - 0.9]))[0]),
            xytext=(-5.5, 2e2), fontsize=8.5, color='#7a6fb0',
            arrowprops=dict(arrowstyle='->', color='#7a6fb0', lw=1.0))
bx.set_ylim(1e-6, 3e3); bx.set_xlim(x[0], x[-1])
bx.set_xticks([]); bx.set_yticks([])
for s in ('top', 'right'):
    bx.spines[s].set_visible(False)
bx.legend(loc='lower center', fontsize=8, frameon=False, ncol=2, bbox_to_anchor=(0.5, -0.32))
bx.set_title(T('(b) 2 サイト: $n^{\\rm Zc}$ の 3 成分（原子核は省略）', '(b) two sites: 3 components of $n^{\\rm Zc}$ (nucleus omitted)'), fontsize=11)

fig.tight_layout()
fig.savefig(os.path.join(here, 'pmt_core' + SUF + '.png'), dpi=150)
np.savetxt(os.path.join(here, 'pmt_core.dat'), np.array([r, ncore(r), nsh(r)]).T, fmt='%.6e',
           header='pmt_core.py (2026-11-13), schematic, R=%.2f, rsm=%.2f, kappa=%.3f, coef=%.6f\nr  n^c_a  n^c_sH,a'
                  % (R, rsm, kappa, coef))
print('kappa', kappa, 'coef', coef)
