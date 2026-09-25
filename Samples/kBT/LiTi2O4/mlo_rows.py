#!/usr/bin/env python3
"""LiTi2O4 6^3: conventional (MTO) QSGW beside MLO-QSGW (chi~ frozen).
One ROW per QSGW iteration, growing downward as the chains proceed.
usage: mlo_rows.py ROOT OUT.png
ROOT holds  ref/bnd_lda.dat, ref/bnd_iter<N>.dat  and  frozen/bnd_lda.dat, frozen/bnd_iter<N>.dat"""
import sys, os, datetime, numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

ROOT, OUT = sys.argv[1], sys.argv[2]
# argv[3] and argv[4] are comma-separated: one MLO chain per column, left to right
MLODIR = (sys.argv[3] if len(sys.argv) > 3 else 'frozen').split(',')
MLOLAB = (sys.argv[4] if len(sys.argv) > 4 else 'MLO-QSGW').split(',')
MLOCOL = (sys.argv[5].split(',') if len(sys.argv) > 5
          else ['tab:green', 'tab:red', 'tab:purple', 'tab:orange'])
NOREF = len(sys.argv) > 6 and sys.argv[6] == 'noref'   # no conventional chain for this mesh
COLS = ([] if NOREF else [('conventional QSGW (MTO), pwmode=1', f'{ROOT}/ref', 'tab:blue')]) + \
       [(MLOLAB[i] if i < len(MLOLAB) else d, f'{ROOT}/{d}', MLOCOL[i % len(MLOCOL)])
        for i, d in enumerate(MLODIR)]

def rd(f):
    xs, es, cx, ce = [], [], [], []
    for ln in open(f):
        if ln.startswith('#'): continue
        t = ln.split()
        if len(t) < 3:
            if cx: xs.append(cx); es.append(ce); cx = []; ce = []
            continue
        cx.append(float(t[1])); ce.append(float(t[2]))
    if cx: xs.append(cx); es.append(ce)
    return np.array(xs[0]), np.array(es)

## Sigma q-mesh points on Gamma->X: (0,0,x) = (n/N)(b1+b2) with (0,0,1) = (b1+b2)/2, so x = 2n/N.
## N=6: 0, 1/3, 2/3, 1.  N=9: 0, 2/9, 4/9, 6/9, 8/9 -- X itself is not on an odd mesh.
_N = int(os.environ.get('MESH', '6'))
QMESH = [2*n/_N for n in range(_N) if 2*n/_N <= 1 + 1e-9]; LO, HI = 32, 44
def path(d, it): return f'{d}/bnd_lda.dat' if it == 0 else f'{d}/bnd_iter{it}.dat'

# rows are driven by the MLO chain (last column); the MTO column is shown alongside
iters = [0] + [i for i in range(1, 31) if any(os.path.exists(path(c[1], i)) for c in (COLS if NOREF else COLS[1:]))]
n = len(iters)
fig, AX = plt.subplots(n, len(COLS), figsize=(4.9 * len(COLS) + 0.4, 2.5 * n + 1.0),
                       sharex=True, sharey=True, squeeze=False)
store = {}
for r, it in enumerate(iters):
    for c, (lab, d, col) in enumerate(COLS):
        ax = AX[r][c]; f = path(d, it)
        if not os.path.exists(f):
            ax.axis('off'); continue
        x, E = rd(f); S = np.sort(E[LO:HI], axis=0); store[(c, it)] = S
        for b in S: ax.plot(x, b, '-', lw=1.2, color='0.45' if it == 0 else col)
        ax.set_title(f'{lab}   {"LDA" if it == 0 else f"iter {it}"}', fontsize=10)
        ax.set_ylim(-0.9, 1.5); ax.set_xlim(0, 1); ax.axhline(0, color='k', lw=0.9)
        for q in QMESH: ax.axvline(q, color='k', ls=':', lw=0.9)
        ax.yaxis.set_major_locator(MultipleLocator(0.5)); ax.yaxis.set_minor_locator(MultipleLocator(0.1))
        ax.grid(axis='y', which='major', color='0.78', lw=0.6)
        ax.grid(axis='y', which='minor', color='0.92', lw=0.4)
    AX[r][0].set_ylabel('$E-E_F$ [eV]')
for a in AX[-1]: a.set_xlabel('$\\Gamma \\to X$')
stamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M')
MESH = os.environ.get('MESH', '6')
fig.suptitle(f'LiTi$_2$O$_4$  {MESH}$^3$ (nkabc = n1n2n3 = {MESH}$^3$)   t$_{{2g}}$ (b33-44) along $\\Gamma\\to X$, 211 points\n'
             'dotted = $\\Sigma$ q-mesh points (interpolation is exact there)\n'
             + ('' if NOREF else
                'CAUTION: different APW cutoffs - MTO chain pwmode=1 (|G|), MLO chains pwmode=11 (|q+G|).\n'
                'Each column states its own damping; the MTO chain always Anderson-mixes sigm at [gw] mixbeta\n') +
             
             f'generated {stamp}', fontsize=10.5)
plt.tight_layout(rect=[0, 0, 1, 1 - 0.95/(2.5*n + 1.0)])
plt.savefig(OUT, dpi=115)
print('wrote', OUT, f'({n} rows x {len(COLS)} cols)')
