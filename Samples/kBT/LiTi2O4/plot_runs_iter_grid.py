#!/usr/bin/env python3
"""Runs side by side, iterations down the rows, in an energy window of the MLO band (2026-09-28).

  plot_runs_iter_grid.py <out.png> <title> <lo>:<hi> <i1>-<i2> <N>:<label>:<pattern with {i}>[:<reference file>] [...]
                        [+<row label>:<file>,<file>,...]

One column per run (N = its Sigma mesh, dotted lines at x = 2n/N), rows i1..i2 read top to bottom (user, 2026-09-28:
the 09-26 fp32 chains iteration by iteration).  An argument starting with '+' adds one more row at the bottom with one
band file per column (e.g. the tf32 runs after 10 iterations under the chains' iteration 10); an empty entry leaves the
panel blank.  A missing file (an iteration that was not kept) leaves the panel blank with a note.  A reference file
(4th field, no {i}) is drawn in blue under the column's bands, e.g. 6^3 under 9^3.  Each panel says how far the bands moved from the row above (max and rms over the points inside the window,
meV, bands compared by index) and the dip of the lowest t2g band (33) below its Gamma value on x < 0.3 (meV).
Band files: the 'i x E-EF' blocks that draw_mloband.sh writes.
"""
import sys, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

out, title, win, rng = sys.argv[1:5]
lo, hi = map(float, win.split(':'))
i1, i2 = map(int, rng.split('-'))
its = list(range(i1, i2 + 1))
cols, extra = [], None
for a in sys.argv[5:]:
    if a.startswith('+'):
        lab, fl = a[1:].split(':', 1)
        extra = (lab, fl.split(','))
    else:
        w = a.split(':')
        cols.append((int(w[0]), w[1], w[2], w[3] if len(w) > 3 else ''))


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


def dev(E, E0):
    nb, nk = min(len(E), len(E0)), min(E.shape[1], E0.shape[1])
    d = (E[:nb, :nk] - E0[:nb, :nk]) * 1e3
    m = (E0[:nb, :nk] > lo) & (E0[:nb, :nk] < hi)
    return np.abs(d[m]).max(), np.sqrt(np.mean(d[m] ** 2))


def dip(x, E):
    b = np.sort(E[32:44], axis=0)[0]
    return (b[0] - b[x < 0.3].min()) * 1e3


rows = [(f'iteration {i}', [c[2].format(i=i) for c in cols]) for i in its]
if extra:
    rows.append(extra)
nr, nc = len(rows), len(cols)
fig, axs = plt.subplots(nr, nc, figsize=(6 * nc, 2.6 * nr + 0.6), squeeze=False)
prev = [None] * nc
for r, (rlab, files) in enumerate(rows):
    for c, (N, lab, _, ref) in enumerate(cols):
        a = axs[r][c]
        f = files[c] if c < len(files) else ''
        if not f or not os.path.isfile(f):
            a.axis('off')
            if f:
                a.text(0.5, 0.5, f'{lab}, {rlab}: not kept', transform=a.transAxes, ha='center', fontsize=8, color='0.5')
            continue
        x, E = bands(f)
        col = 'C3' if N == 9 else 'C0'
        if ref:
            xr, Er = bands(ref)
            for b in Er:
                if b.max() > lo - 0.5 and b.min() < hi + 0.5:
                    a.plot(xr, b, color='C0', lw=0.6)
        for b in E:
            if b.max() > lo - 0.5 and b.min() < hi + 0.5:
                a.plot(x, b, color=col, lw=0.6)
        for n in range(1, N):
            if 2 * n / N < 1:
                a.axvline(2 * n / N, color='0.6', ls=':', lw=0.6)
        a.axhline(0, color='0.6', lw=0.4)
        a.set_xlim(0, 1); a.set_ylim(lo, hi)
        a.set_xticks([0, 1]); a.set_xticklabels([r'$\Gamma$', 'X'] if r == nr - 1 else ['', ''])
        a.tick_params(labelsize=7)
        s = f'{lab}, {rlab}: dip {dip(x, E):.1f} meV'
        if prev[c] is not None:
            mx, rms = dev(E, prev[c])
            s += f'; from the row above max {mx:.0f}, rms {rms:.1f} meV'
        a.text(0.01, 0.97, s, transform=a.transAxes, va='top', fontsize=7.5, backgroundcolor='white')
        prev[c] = E
        if c == 0:
            a.set_ylabel(r'$E-E_F$ (eV)', fontsize=8)
fig.suptitle(title, fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 1 - 0.4 / (2.6 * nr + 0.6)))
fig.savefig(out, dpi=90)
print('->', out)
