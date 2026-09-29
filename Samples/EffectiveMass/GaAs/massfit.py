#!/usr/bin/env python3
'''Effective masses of a zincblende semiconductor from the mass-mode band files of lmf (numpy fit; no gnuplot).
usage: python3 massfit.py [--emin 0.01] [--emax 0.05] [--dir .]

lmf --band with a syml file in the mass mode writes Band<ib>Syml<is>Spin1.dat for the bands near the gap:
  column 5: |k| (1/bohr) measured from the left end of the line (here Gamma)
  column 7: E(k) - E(left end) (eV)
With SOC, the eight files along each line are, from below (pairs are Kramers or nearly degenerate pairs):
  1,2 split-off (so), 3,4 light hole (lh), 5,6 heavy hole (hh), 7,8 conduction (e).
For each band, the points with emin <= |E| <= emax are fitted by linear least squares to
  k^2 = a |E| + b E^2 ,  mass m*/m = a*Ry,  E0 = a/b           (Ry = 13.605 eV)
that is  hbar^2 k^2/(2 m*) = |E| (1 + |E|/E0)  (non-parabolicity E0 in eV).
Output: mass.txt (masses, gap and spin-orbit splitting at Gamma), mass_detail.txt (with E0 and the number of points)
        and massfit.npz (the fitted points and coefficients, for massplot.py).
'''
import argparse, datetime, os, sys
import numpy as np
RY = 13.605
LINES = {1: '111', 2: '100', 3: '110'}                 # order of the lines in syml.<sname>
BANDS = {'so': (1, 2), 'lh': (3, 4), 'hh': (5, 6), 'e': (7, 8)}

def readband(fname):
    'k (1/bohr), |E - E(Gamma)| (eV) of one band file, and E(Gamma) relative to the valence band maximum (eV)'
    k, e, eg = [], [], None
    for line in open(fname):
        if line.startswith('#') or not line.strip(): continue
        w = line.split()
        if eg is None: eg = float(w[5])
        k.append(float(w[4])); e.append(abs(float(w[6])))
    return np.array(k), np.array(e), eg

def fit(k, e, emin, emax):
    'least squares for k^2 = a e + b e^2 in emin <= e <= emax; returns mass, E0, number of points'
    sel = (e >= emin) & (e <= emax)
    if sel.sum() < 3: return np.nan, np.nan, int(sel.sum()), np.nan, np.nan
    A = np.vstack([e[sel], e[sel]**2]).T
    (a, b), *_ = np.linalg.lstsq(A, k[sel]**2, rcond=None)
    return a*RY, a/b, int(sel.sum()), a, b

if __name__ == '__main__':
    p = argparse.ArgumentParser(description='effective mass by a fit of k^2 = a|E| + b E^2')
    p.add_argument('--emin', type=float, default=0.01, help='lower end of the fitting window (eV)')
    p.add_argument('--emax', type=float, default=0.05, help='upper end of the fitting window (eV)')
    p.add_argument('--dir', default='.', help='directory with Band*Syml*Spin1.dat')
    args = p.parse_args()
    head = ['# effective mass m*/m from hbar^2 k^2/(2 m*) = |E|(1+|E|/E0), E measured from the band edge at Gamma',
            f'# fitting window {args.emin} <= |E| <= {args.emax} eV']
    out = head + ['# line band  file_ib     mass']                                  # mass.txt (compared by test.py)
    det = head + ['# npts = number of points in the window; E0 is ill-determined when the band is nearly parabolic',
                  '# line band  file_ib     mass        E0(eV)   npts']              # mass_detail.txt
    store = {}
    for isyml, dname in LINES.items():
        for bname, ibs in BANDS.items():
            for ib in ibs:
                fname = os.path.join(args.dir, f'Band{ib:03d}Syml{isyml:03d}Spin1.dat')
                if not os.path.exists(fname): sys.exit(f'massfit.py: {fname} not found')
                k, e, eg = readband(fname)
                mass, e0, n, a, b = fit(k, e, args.emin, args.emax)
                out.append(f'  {dname}  {bname:3s}  {ib:4d}   {mass:10.4f}')
                det.append(f'  {dname}  {bname:3s}  {ib:4d}   {mass:10.4f}  {e0:12.3f}  {n:5d}')
                store[f'k_{dname}_{ib}'] = k; store[f'e_{dname}_{ib}'] = e
                store[f'ab_{dname}_{ib}'] = np.array([a, b])
                if isyml == 1 and ib in (1, 7): store[f'edge_{ib}'] = eg
    # band edges at Gamma relative to the valence band maximum (column 6 of the first point)
    out += ['# energies at Gamma relative to the valence band maximum (eV): gap and spin-orbit splitting',
            f'  gap_Gamma      {store["edge_7"]:10.4f}',
            f'  split_off      {store["edge_1"]:10.4f}']
    det += out[-3:]
    store['window'] = np.array([args.emin, args.emax])
    store['source'] = os.path.abspath(args.dir)                      # where the band files were
    store['made'] = datetime.datetime.now().strftime('%Y-%m-%d %H:%M')
    open('mass.txt', 'w').write('\n'.join(out) + '\n')
    open('mass_detail.txt', 'w').write('\n'.join(det) + '\n')
    np.savez('massfit.npz', **store)
    print('\n'.join(det))
