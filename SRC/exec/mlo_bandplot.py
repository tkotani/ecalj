#!/usr/bin/env python3
"""MLO vs first-principles band plot.

Grey lines are the PMT first-principles bands (bnd00*.spin*), red markers the
MLO model (band_MLO_spin*.dat). Both are put on the same zero, E_F:

  bnd00*.spin*        column 3, already (E - basel) in eV, basel = eferm
  band_MLO_spin*.dat  column 2, RAW Ry -- we subtract ef and convert

ef comes from bandplot_MLO.isp1.glt ("ef=") when present, else from the first
line of qplist.dat. Both are written by the same writeband call, so they agree.

Symmetry-line labels are read from bandplot.isp1.glt's `set xtics (...)`.

Usage:
    mlo_bandplot.py <sampledir> [-o out.png] [--emin -6] [--emax 8] [--full]

--full draws the whole range the MLO model spans (useful for checking that the
truncation is smooth rather than ringing); otherwise the window around E_F.

2026-09-16 T.K. + Claude
"""
import argparse
import glob
import os
import re
import sys

RY = 13.6057

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def read_ef(d):
    """Fermi level in Ry, as used for the MLO file's zero."""
    g = os.path.join(d, 'bandplot_MLO.isp1.glt')
    if os.path.exists(g) and os.path.getsize(g) > 0:
        m = re.search(r'ef=\s*([-\d.eEdD+]+)', open(g).read())
        if m:
            return float(m.group(1).replace('D', 'E').replace('d', 'e'))
    q = os.path.join(d, 'qplist.dat')
    if os.path.exists(q):
        return float(open(q).readline().split()[0])
    raise SystemExit(f'{d}: no ef (need bandplot_MLO.isp1.glt or qplist.dat)')


def read_xtics(d):
    """[(x, label), ...] from the gnuplot file job_band wrote."""
    for name in ('bandplot.isp1.glt', 'bandplot.isp2.glt'):
        p = os.path.join(d, name)
        if not os.path.exists(p):
            continue
        txt = open(p).read()
        m = re.search(r'set xtics \((.*?)\)\n', txt, re.S)
        if not m:
            continue
        out = []
        for lab, x in re.findall(r"'([^']*)'\s+([-\d.]+)", m.group(1)):
            out.append((float(x), lab.replace('GAMMA', r'$\Gamma$')))
        return out
    return []


def load_dft(d, isp):
    """[(x, E)] in eV relative to E_F, split into per-band segments."""
    segs = []
    for f in sorted(glob.glob(os.path.join(d, f'bnd0*.spin{isp}'))):
        cur = {}
        for line in open(f):
            t = line.split()
            if not t or line.startswith('#'):
                continue
            try:
                ib, x, e = int(t[0]), float(t[1]), float(t[2])
            except ValueError:
                continue
            cur.setdefault(ib, []).append((x, e))
        for ib in sorted(cur):
            segs.append(np.array(cur[ib]))
    return segs


def load_mlo(d, isp, ef):
    """(x, E) points in eV relative to E_F."""
    f = os.path.join(d, f'band_MLO_spin{isp}.dat')
    if not os.path.exists(f):
        return None
    pts = []
    for line in open(f):
        t = line.split()
        if not t or line.startswith('#'):
            continue
        try:
            pts.append((float(t[0]), (float(t[1]) - ef) * RY))
        except ValueError:
            continue
    return np.array(pts) if pts else None


def panel(ax, d, isp, ef, emin, emax, title):
    for s in load_dft(d, isp):
        ax.plot(s[:, 0], s[:, 1], '-', color='0.55', lw=0.8, zorder=1)
    m = load_mlo(d, isp, ef)
    if m is not None:
        ax.plot(m[:, 0], m[:, 1], 'x', color='#c0392b', ms=3.2, mew=0.8, zorder=3)
    ax.axhline(0.0, color='#2c3e50', lw=0.7, ls='--', zorder=2)
    tics = read_xtics(d)
    if tics:
        ax.set_xticks([x for x, _ in tics])
        ax.set_xticklabels([l for _, l in tics])
        for x, _ in tics[1:-1]:
            ax.axvline(x, color='0.85', lw=0.6, zorder=0)
        ax.set_xlim(tics[0][0], tics[-1][0])
    ax.set_ylim(emin, emax)
    ax.set_ylabel(r'$E-E_F$  (eV)')
    ax.set_title(title, fontsize=10)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('dir')
    ap.add_argument('-o', '--out', default=None)
    ap.add_argument('--emin', type=float, default=None)
    ap.add_argument('--emax', type=float, default=None)
    ap.add_argument('--full', action='store_true',
                    help='span the whole range the MLO model reaches')
    ap.add_argument('--label', default=None)
    a = ap.parse_args()

    d = a.dir.rstrip('/')
    name = a.label or os.path.basename(d).replace('_work', '')
    ef = read_ef(d)

    spins = [s for s in (1, 2)
             if os.path.exists(os.path.join(d, f'band_MLO_spin{s}.dat'))]
    if not spins:
        sys.exit(f'{d}: no band_MLO_spin*.dat')

    if a.emin is None or a.emax is None:
        allm = [load_mlo(d, s, ef) for s in spins]
        allm = [m for m in allm if m is not None]
        lo = min(m[:, 1].min() for m in allm)
        hi = max(m[:, 1].max() for m in allm)
        if a.full:
            pad = 0.05 * (hi - lo)
            emin, emax = lo - pad, hi + pad
        else:
            emin, emax = max(lo - 1.0, -8.0), min(hi + 1.0, 10.0)
    else:
        emin, emax = a.emin, a.emax

    fig, axes = plt.subplots(1, len(spins), figsize=(5.2 * len(spins), 4.2),
                             squeeze=False)
    for ax, s in zip(axes[0], spins):
        tag = f'{name}' if len(spins) == 1 else f'{name}  spin{s}'
        panel(ax, d, s, ef, emin, emax, tag)
    fig.tight_layout()
    out = a.out or f'{name}.png'
    fig.savefig(out, dpi=130)
    print(f'{out}  ({len(spins)} spin, E in [{emin:.2f},{emax:.2f}] eV)')


if __name__ == '__main__':
    main()
