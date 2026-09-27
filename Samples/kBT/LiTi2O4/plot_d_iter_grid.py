#!/usr/bin/env python3
"""The whole Ti 3d manifold iteration by iteration, one row per iteration (2026-09-28, user: arrange vertically).

  plot_d_iter_grid.py <out.png> <title> <ref MLO band> <ref sigm band> <MLO pattern with {i}> <sigm pattern with {i}> <i1>-<i2> [N]

Left column: MLO band, right column: conventional sigm band.  Each row draws iteration i (red) on the reference (blue, e.g. 6^3
after 10 iterations), from -0.7 to 4.6 eV (t2g and eg).  Dotted lines: the Sigma mesh points of the run (x = 2n/N, N default 9).
Each panel says how far the bands moved from the previous iteration (max and rms over the points within the window, meV;
bands compared by index, i.e. sorted per k) and how far they are from the reference; the same numbers are printed as a table.
"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

out, title, rm, rs, pm, ps, rng = sys.argv[1:8]
N = int(sys.argv[8]) if len(sys.argv) > 8 else 9
i1, i2 = map(int, rng.split('-'))
its = list(range(i1, i2 + 1))
lo, hi = -0.7, 4.6


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


ref = [bands(rm), bands(rs)]
prev = [None, None]
print('iteration   MLO band: change max/rms, vs ref max/rms (meV)   sigm band: change max/rms, vs ref max/rms (meV)')
fig, axs = plt.subplots(len(its), 2, figsize=(12, 2.9 * len(its) + 0.6), squeeze=False)
for r, i in enumerate(its):
    row = f'{i:9d}'
    for c, (pat, name) in enumerate(((pm, 'MLO band'), (ps, 'sigm band'))):
        a = axs[r][c]
        x, E = bands(pat.format(i=i))
        xr, Er = ref[c]
        rmx, rrms = dev(E, Er)
        if prev[c] is None:
            lab = f'{name}, iteration {i}: vs ref max {rmx:.0f}, rms {rrms:.1f} meV'
            row += f'   {"-":>7s} {"-":>6s}  {rmx:7.1f} {rrms:6.1f}          '
        else:
            cmx, crms = dev(E, prev[c])
            lab = f'{name}, iteration {i}: change from {i - 1} max {cmx:.1f}, rms {crms:.1f} meV; vs ref rms {rrms:.1f} meV'
            row += f'   {cmx:7.1f} {crms:6.1f}  {rmx:7.1f} {rrms:6.1f}          '
        prev[c] = E
        for b in Er:
            if b.max() > lo - 0.5 and b.min() < hi + 0.5:
                a.plot(xr, b, color='C0', lw=0.6)
        for b in E:
            if b.max() > lo - 0.5 and b.min() < hi + 0.5:
                a.plot(x, b, color='C3', lw=0.6)
        for n in range(1, N):
            if 2 * n / N < 1:
                a.axvline(2 * n / N, color='0.6', ls=':', lw=0.6)
        a.axhline(0, color='0.6', lw=0.4)
        a.set_xlim(0, 1); a.set_ylim(lo, hi)
        a.set_xticks([0, 1]); a.set_xticklabels([r'$\Gamma$', 'X'] if r == len(its) - 1 else ['', ''])
        a.tick_params(labelsize=7)
        a.text(0.01, 0.97, lab, transform=a.transAxes, va='top', fontsize=8, backgroundcolor='white')
        if c == 0:
            a.set_ylabel(r'$E-E_F$ (eV)', fontsize=8)
    print(row)
fig.suptitle(title, fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 1 - 0.4 / (2.9 * len(its) + 0.6)))
fig.savefig(out, dpi=90)
print('->', out)
