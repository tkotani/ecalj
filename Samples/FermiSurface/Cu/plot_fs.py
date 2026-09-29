#!/usr/bin/env python3
"""Cross sections of the Fermi surface of fcc Cu from fermiup.bxsf (written by job_fermisurface).

  python3 plot_fs.py [fermiup.bxsf]

Writes fs_section.png and fs_section.npz (the plotted contour lines, the source file and the date).
The bxsf file holds, for each band, the eigenvalues (Ry) on the grid k = (i1/N1) b1 + (i2/N2) b2 + (i3/N3) b3,
i = 0..N (both ends included, the last index runs fastest); b1, b2, b3 are in units of 2 pi/alat.
Here the eigenvalues are interpolated linearly to two planes through Gamma: (001) and (1-10).
"""
import datetime
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def read_bxsf(fname):
    lines = open(fname).read().split('\n')
    ef = [float(l.split(':')[1]) for l in lines if 'Fermi Energy' in l][0]
    i0 = [i for i, l in enumerate(lines) if 'BEGIN_BANDGRID_3D' in l][0]
    nband = int(lines[i0 + 1])
    n = [int(x) for x in lines[i0 + 2].split()]
    qlat = np.array([[float(x) for x in lines[i0 + 4 + j].split()] for j in range(3)])   # rows b1, b2, b3
    bands, cur = {}, None
    for l in lines[i0 + 7:]:
        if 'BAND:' in l:
            cur = int(l.split()[1])
            bands[cur] = []
        elif 'END_BANDGRID' in l:
            break
        elif cur is not None:
            bands[cur] += [float(x) for x in l.split()]
    assert len(bands) == nband
    return ef, n, qlat, {b: np.array(v).reshape(n) for b, v in bands.items()}


def interpolate(e, qlat, k):
    """Trilinear interpolation of the periodic grid e[i1,i2,i3] at the Cartesian points k[...,3]."""
    n = np.array(e.shape) - 1
    f = k @ np.linalg.inv(qlat)              # k = f1 b1 + f2 b2 + f3 b3
    f = (f - np.floor(f)) * n
    i = np.minimum(np.floor(f).astype(int), n - 1)
    d = f - i
    out = 0.0
    for a in (0, 1):
        for b in (0, 1):
            for c in (0, 1):
                w = (d[..., 0] if a else 1 - d[..., 0]) * (d[..., 1] if b else 1 - d[..., 1]) * (d[..., 2] if c else 1 - d[..., 2])
                out = out + w * e[i[..., 0] + a, i[..., 1] + b, i[..., 2] + c]
    return out


def bzfun(k):
    """1 on the boundary of the first Brillouin zone of fcc (units 2 pi/alat), < 1 inside."""
    a = abs(k)
    return np.maximum(a.max(axis=-1), a.sum(axis=-1) / 1.5)


def crossing(e, qlat, ef, k0, d, tmax):
    """Smallest t in (0, tmax) with E(k0 + t d) = Ef (d is normalized here); nan when there is none."""
    t = np.linspace(0, tmax, 2001)
    v = interpolate(e, qlat, np.array(k0) + t[:, None] * np.array(d) / np.linalg.norm(d)) - ef
    i = np.nonzero(v[:-1] * v[1:] < 0)[0]
    return t[i[0]] - v[i[0]] * (t[i[0] + 1] - t[i[0]]) / (v[i[0] + 1] - v[i[0]]) if len(i) else np.nan


fname = sys.argv[1] if len(sys.argv) > 1 else 'fermiup.bxsf'
ef, n, qlat, bands = read_bxsf(fname)
cross = [b for b, e in bands.items() if e.min() < ef < e.max()]
print(f'{fname}: grid {n}, Ef = {ef:.5f} Ry, bands in the file {sorted(bands)}, bands crossing Ef {cross}')
for b in cross:
    k100 = crossing(bands[b], qlat, ef, [0, 0, 0], [1, 0, 0], 1.0)
    k110 = crossing(bands[b], qlat, ef, [0, 0, 0], [1, 1, 0], 1.0)
    neck = crossing(bands[b], qlat, ef, [0.5, 0.5, 0.5], [1, 1, -2], 0.5)
    print(f'band {b}: k_F along [100] = {k100:.3f}, along [110] = {k110:.3f}, radius of the neck at L (along [11-2]) = {neck:.3f}  (2 pi/a)')
now = datetime.datetime.now().strftime('%Y-%m-%d %H:%M')
s = np.linspace(-1.3, 1.3, 521)
x, y = np.meshgrid(s, s, indexing='ij')
planes = {'(001) plane through Gamma, X, K': (np.stack([x, y, 0 * x], axis=-1), '$k_x$', '$k_y$'),
          '(1-10) plane through Gamma, X, L, K': (np.stack([x / np.sqrt(2), x / np.sqrt(2), y], axis=-1), '$k_{[110]}$', '$k_z$')}
fig, axes = plt.subplots(1, 2, figsize=(8.4, 4.4))
save = {'date': now, 'source': os.path.abspath(fname), 'ef_ry': ef, 'grid': n,
        'kf100_kf110_neck': np.array([k100, k110, neck])}
for ip, (ax, (title, (k, xl, yl))) in enumerate(zip(axes, planes.items())):
    bz = bzfun(k)
    ax.contour(x, y, bz, levels=[1.0], colors='0.5', linewidths=0.5)
    for b in cross:
        e = np.where(bz <= 1.0, interpolate(bands[b], qlat, k) - ef, np.nan)
        cs = ax.contour(x, y, e, levels=[0.0], colors='k', linewidths=0.6)
        # the contour lines E = Ef in the coordinates of the plot; the pieces are separated by a row of nan
        save[f'plane{ip}_band{b}_contour'] = np.concatenate([np.vstack([seg, [[np.nan, np.nan]]]) for seg in cs.allsegs[0]])
    save[f'plane{ip}_title'] = title
    ax.set_aspect('equal')
    ax.set_xlabel(xl + r' ($2\pi/a$)')
    ax.set_ylabel(yl + r' ($2\pi/a$)')
    ax.set_title(title, fontsize=9)
    ax.plot(0, 0, 'k.', ms=3)
    ax.annotate(r'$\Gamma$', (0, 0), xytext=(3, 3), textcoords='offset points', fontsize=8)
axes[0].annotate('X', (1, 0), xytext=(3, 0), textcoords='offset points', fontsize=8)
axes[0].annotate('K', (0.75, 0.75), xytext=(3, 3), textcoords='offset points', fontsize=8)
axes[1].annotate('X', (0, 1), xytext=(0, 3), textcoords='offset points', fontsize=8)
axes[1].annotate('L', (np.sqrt(0.5), 0.5), xytext=(3, 3), textcoords='offset points', fontsize=8)
axes[1].annotate('K', (0.75 * np.sqrt(2), 0), xytext=(3, 0), textcoords='offset points', fontsize=8)
fig.suptitle(f'Cu Fermi surface (band {cross}), Ef = {ef:.5f} Ry, grid {n[0]}x{n[1]}x{n[2]}   [{now}]', fontsize=9)
fig.tight_layout()
fig.savefig('fs_section.png', dpi=150)
np.savez_compressed('fs_section.npz', **save)
print('fs_section.png and fs_section.npz written')
