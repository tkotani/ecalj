#!/usr/bin/env python3
"""How much the MLO bands move from one iteration to the next, for the six patterns (2026-09-28, user: "MLO バンドの反復ごとの動きで見て").

  plot_mlo_conv_six.py <root> <out.png> <out.npz>

<root> is the panel directory of update_six_rows.sh: <root>/<column>/bnd_iter<N>.dat (links) with .src (the band file on kt1)
and .label (a stand-in from another run).  For every column and N the change E(N) - E(N-1) of the MLO band (energy-sorted at
each k) is taken only when both come from the same run (the same directory in .src), so a column that switches from an older
run to its own run gets no point at the seam.  Three panels, max |dE| (meV, log scale) against N:
  t2g (b33-44), the occupied pair (the two lowest t2g on x < 0.5), eg (the bands at 2.5-4.8 eV).
Solid (9^3) / dashed (6^3): the current code (the column's own run or the earlier tf32 gwsc 10 run); dotted: an older code
standing in (09-26 / 09-27).
The numbers go to <out.npz>: <column>__N, __t2g_max, __t2g_rms, __occ_max, __eg_max, __standin (1 = a stand-in).
"""
import sys, os, glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root, out, outnpz = sys.argv[1:4]
COLS = [('k9tf32', '9$^3$ tf32', 'C0', '-'), ('k9tf32es', '9$^3$ tf32 ES', 'C4', '-'),   # ES: empty spheres (2026-09-28)
        ('k9fp32', '9$^3$ fp32', 'C1', '-'), ('k9fp64', '9$^3$ fp64', 'C2', '-'),
        ('k6tf32', '6$^3$ tf32', 'C0', '--'), ('k6fp32', '6$^3$ fp32', 'C1', '--'), ('k6fp64', '6$^3$ fp64', 'C2', '--')]


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


def side(f, ext):
    return open(f + '.' + ext).read().strip() if os.path.exists(f + '.' + ext) else ''


fig, ax = plt.subplots(1, 3, figsize=(17, 5.2))
arrs = {}
for col, lab, color, ls in COLS:
    files = {int(os.path.basename(f)[8:-4]): f for f in glob.glob(f'{root}/{col}/bnd_iter*.dat')}
    rows = []
    for n in sorted(files):
        if n - 1 not in files:
            continue
        f1, f0 = files[n], files[n - 1]
        if os.path.dirname(side(f1, 'src')) != os.path.dirname(side(f0, 'src')):
            continue                                            # different runs: no change across the seam
        x, E1 = bands(f1); _, E0 = bands(f0)
        nb, nk = min(len(E1), len(E0)), min(E1.shape[1], E0.shape[1])
        E1, E0, x = np.sort(E1[:nb, :nk], axis=0), np.sort(E0[:nb, :nk], axis=0), x[:nk]
        d = (E1 - E0) * 1e3
        t = d[32:44]
        occ = np.abs(d[32:34][:, x < 0.5]).max()
        m = (E0 > 2.5) & (E0 < 4.8)
        older = 'code' in side(f1, 'label')                     # [09-26 code] / [09-27 code]; the tf32 gwsc 10 run is current code
        rows.append((n, np.abs(t).max(), np.sqrt((t ** 2).mean()), occ, np.abs(d[m]).max(), 1 if older else 0))
    if not rows:
        continue
    R = np.array(rows)
    for k, name in enumerate(('N', 't2g_max', 't2g_rms', 'occ_max', 'eg_max', 'standin')):
        arrs[f'{col}__{name}'] = R[:, k]
    for a, j in zip(ax, (1, 3, 4)):
        own = R[:, 5] == 0
        # own run solid, stand-ins dotted; points joined only within a run (the seam has no point anyway)
        for sel, style, alpha in ((own, ls, 1.0), (~own, ':', 0.6)):
            if sel.any():
                a.plot(R[sel, 0], R[sel, j], style, marker='o', ms=3.5, lw=1.2, color=color, alpha=alpha,
                       label=lab + ('' if style != ':' else ' (stand-in)'))
for a, t in zip(ax, ('t2g (b33-44): max |$\\Delta E$| from the previous iteration',
                     'the occupied pair (two lowest t2g, x < 0.5)', 'eg (2.5-4.8 eV)')):
    a.set_yscale('log'); a.set_xlabel('iteration N'); a.set_ylabel('max |E(N) - E(N-1)| (meV)')
    a.set_title(t, fontsize=10); a.grid(which='both', color='0.9', lw=0.5)
    a.set_xticks(range(0, 23, 2))
ax[0].legend(fontsize=8, ncol=2)
fig.suptitle('LiTi$_2$O$_4$ MLO band, change per iteration (solid: the current code, 9$^3$ solid / 6$^3$ dashed line; '
             'dotted: stand-ins from older runs)', fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.94))
fig.savefig(out, dpi=110)
np.savez_compressed(outnpz, **arrs)
print('wrote', out, outnpz)
