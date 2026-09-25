#!/usr/bin/env python3
"""Relative real-space decay |Sigma(R)|/|Sigma(0)| of the conventional (MTO) and the
MLO representation.  The absolute scales are NOT comparable: sigm carries O^-1 on both
sides (on-site 73.8 eV) while Sigma^MLO = c^dag Sigma^psi c keeps O^-1 for getsenex.
Only the decay relative to the on-site block is meaningful."""
import sys, collections, numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.sans-serif'] = ['IPAGothic', 'Droid Sans Fallback', 'DejaVu Sans']
matplotlib.rcParams['axes.unicode_minus'] = False
DAT, OUT = sys.argv[1], sys.argv[2]
b = collections.defaultdict(lambda: [0., 0., 0])
for ln in open(DAT):
    if ln.startswith('#'): continue
    t = ln.split(); r = round(float(t[3]), 4)
    a = b[r]; a[0] = max(a[0], float(t[4])); a[1] += float(t[5]) * int(t[8]); a[2] += int(t[8])
mto = {r: (v[0], v[1]/v[2]) for r, v in b.items()}
mlo1 = {0.0:(187.05,3.629),0.7071:(50.21,0.548),1.0:(2.34,0.168),1.2247:(4.85,0.100),
        1.4142:(2.54,0.065),1.5811:(2.54,0.045),1.7321:(1.54,0.036),1.8708:(1.54,0.030),
        2.0:(2.54,0.031),2.1213:(2.54,0.027),2.2361:(1.24,0.024),2.3452:(13.04,0.026),
        2.4495:(2.54,0.023),2.5495:(3.70,0.032),2.7386:(2.54,0.032),2.9155:(1.54,0.034),
        3.0:(2.34,0.028),3.0822:(0.67,0.019),3.1623:(1.24,0.023),3.3166:(0.25,0.019),
        3.4641:(0.13,0.019),3.5355:(0.05,0.013),3.6742:(0.13,0.020)}
mlo10 = {0.0:(155.52,3.697),0.7071:(43.30,0.568),1.0:(3.96,0.200),1.2247:(9.21,0.129),
         1.4142:(9.21,0.086),1.5811:(4.82,0.064),1.7321:(2.03,0.048),1.8708:(9.21,0.047),
         2.0:(4.82,0.051),2.1213:(4.82,0.048),2.2361:(2.49,0.042),2.3452:(13.80,0.041),
         2.4495:(9.21,0.038),2.5495:(9.21,0.050),2.7386:(4.82,0.041),2.9155:(1.65,0.034),
         3.0:(2.05,0.056),3.0822:(0.51,0.017),3.1623:(0.60,0.019),3.3166:(0.15,0.015),
         3.4641:(0.12,0.015),3.5355:(0.03,0.011),3.6742:(0.10,0.015)}
sets = [('従来 QSGW (MTO) iter 1', mto, 'tab:blue', 'o-'),
        ('MLO 76 orb iter 1 (LDA 由来)', mlo1, 'tab:red', 'o-'),
        ('MLO 76 orb iter 10', mlo10, 'tab:orange', 's--')]
fig, AX = plt.subplots(1, 2, figsize=(13, 5))
for k, (col, ttl) in enumerate([(1, 'ブロックごとの平均'), (0, 'ブロックごとの最大')]):
    ax = AX[k]
    for lab, d, c, m in sets:
        r = sorted(d); y = [d[x][col]/d[0.0][col] for x in r]
        ax.semilogy(r, y, m, color=c, lw=1.8, ms=5, label=lab)
    ax.set_xlabel('$|R|$ [alat]'); ax.set_ylabel('$|\\Sigma(R)|\\,/\\,|\\Sigma(0)|$')
    ax.grid(which='both', color='0.88', lw=0.5); ax.set_title(ttl, fontsize=11)
    ax.legend(fontsize=9)
AX[0].set_ylim(1e-3, 1.5)
fig.suptitle('LiTi$_2$O$_4$ 6$^3$: $\\Sigma(R)$ のオンサイト規格化した減衰  '
             '— 絶対値は正規化が違うので比較できない', fontsize=12)
plt.tight_layout(rect=[0, 0, 1, 0.94]); plt.savefig(OUT, dpi=130)
print('wrote', OUT)
