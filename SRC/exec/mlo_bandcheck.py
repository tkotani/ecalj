#!/usr/bin/env python3
"""MLO model bands against the first-principles bands, along the band path.

    mlo_bandcheck.py <dir> [<dir> ...] [--json out.json] [--tsv out.tsv]

<dir> is a directory where job_band and job_mlo (or job_mlo_soc) have run:
  DFT  bnd*.spin{1,2}           columns: index, x, E - E_F (eV)
  MLO  band_MLO_spin{1,2}.dat   columns: x, E (Ry, raw; E_F is the ef= of bandplot_MLO.isp1.glt, else qplist.dat)
  efermi.lmf                    E_F, top of the valence and bottom of the conduction band of the SCF mesh (Ry)

Metal or insulator is decided by efermi.lmf: an insulator has E_F on the top of the valence and the bottom of the
conduction band above it (gap_mesh > 0.05 eV); a metal has Top < E_F < Bottom (the eigenvalues next to E_F). The path
alone misjudges metals whose bands jump over E_F between two k points.

Insulator: VBM_D = max{E_D <= s}, CBM_D = min{E_D > s} over both spins along the path, s = gap_mesh/2 (a coarse mesh can
miss the VBM of the path); the MLO edges are taken on either side of the DFT mid-gap of the path.  gapD, gapM, dVBM = VBM_M - VBM_D, dCBM = CBM_M - CBM_D.
Accuracy, with the nearest band at the same k (|dx| < 0.01):
  rms_m2d, max_m2d : each MLO point in [VBM_D - 8, CBM_D + 3] eV -> nearest DFT band   (a wrong band shows up here)
  rms_d2m, max_d2m : each DFT point in [VBM_D - 8, CBM_D + 1] eV -> nearest MLO band   (a missing band shows up here)
Metals: the same with VBM_D = CBM_D = E_F = 0.

2026-10-01: rewritten for any directory (the test of Samples/MATERIALS, MD/research_log.md 2026-10-01).
The earlier version (v6.2, effective-mass ratios, fixed to Samples/MLOsamples/*__m3_work) is in the git history.
"""
import argparse, glob, json, os, re, sys
import numpy as np

RY = 13.6057


def read_ef(d):
    p = os.path.join(d, 'bandplot_MLO.isp1.glt')
    if os.path.exists(p):
        m = re.search(r'ef=\s*([-+0-9.eEdD]+)', open(p).read())
        if m: return float(m.group(1).replace('D', 'E'))
    return float(open(os.path.join(d, 'qplist.dat')).readline().split()[0])


def load(files, conv, bnd=False):
    P = {}
    for f in files:
        for l in open(f):
            if l.lstrip().startswith('#'): continue
            t = l.split()
            if len(t) < 2: continue
            try: x, e = (float(t[1]), conv(float(t[2]))) if bnd else (float(t[0]), conv(float(t[1])))
            except (ValueError, IndexError): continue
            P.setdefault(round(x, 5), []).append(e)
    return P


def nearest(P, Q, lo, hi):
    """for each point of P in [lo,hi]: |E - nearest E of Q at the same x|"""
    xs = np.array(sorted(Q)); dev = []
    for x, es in P.items():
        i = np.searchsorted(xs, x)
        cand = [xs[k] for k in (i - 1, i) if 0 <= k < len(xs)]
        if not cand: continue
        xq = min(cand, key=lambda c: abs(c - x))
        if abs(xq - x) > 0.01: continue
        qe = np.array(Q[xq])
        dev += [np.min(np.abs(qe - e)) for e in es if lo <= e <= hi]
    return np.array(dev)


def check(d):
    ef = read_ef(d)
    spins = [s for s in (1, 2) if glob.glob(os.path.join(d, f'bnd*.spin{s}'))]
    D, M = {}, {}
    for s in spins:
        for x, es in load(sorted(glob.glob(os.path.join(d, f'bnd*.spin{s}'))), lambda e: e, bnd=True).items():
            D.setdefault(x, []).extend(es)
        mf = os.path.join(d, f'band_MLO_spin{s}.dat')
        if os.path.exists(mf) and os.path.getsize(mf) > 0:
            for x, es in load([mf], lambda e: (e - ef) * RY).items(): M.setdefault(x, []).extend(es)
    if not D or not M: return None
    eD = np.concatenate([np.array(v) for v in D.values()]); eM = np.concatenate([np.array(v) for v in M.values()])
    ln = [float(l.split()[0].replace('D', 'E')) for l in open(os.path.join(d, 'efermi.lmf')).readlines()[:3]]
    gmesh = (ln[2] - ln[1]) * RY
    # the band edges along the path are separated at the mid-gap of the mesh: with a coarse mesh the path goes through
    # k points above the mesh VBM (LaGaO3, 3x2x3: E - E_F up to +0.1 eV in the valence band)
    sep = 0.5 * gmesh if (ln[0] - ln[1]) < 1e-6 and gmesh > 0.05 else 0.02
    vbD = eD[eD <= sep].max(); cbD = eD[eD > sep].min()
    ins = bool((ln[0] - ln[1]) < 1e-6 and gmesh > 0.05 and cbD - vbD > 0.05)
    r = dict(insulator=ins, gap_mesh=(ln[2] - ln[1]) * RY if ins else 0.0, nspin=len(spins))
    if ins:
        mid = 0.5 * (vbD + cbD)
        vbM = eM[eM <= mid].max() if (eM <= mid).any() else float('nan')
        cbM = eM[eM > mid].min() if (eM > mid).any() else float('nan')
        r.update(gapD=float(cbD - vbD), gapM=float(cbM - vbM), dVBM=float(vbM - vbD), dCBM=float(cbM - cbD))
        top = cbD
    else:
        vbD, top = 0.0, 0.0
    m2d = nearest(M, D, vbD - 8, top + 3); d2m = nearest(D, M, vbD - 8, top + 1)
    for k, a in (('m2d', m2d), ('d2m', d2m)):
        r['rms_' + k] = float(np.sqrt((a ** 2).mean())) if a.size else float('nan')
        r['max_' + k] = float(a.max()) if a.size else float('nan')
    return r


KEYS = ['insulator', 'gap_mesh', 'gapD', 'gapM', 'dVBM', 'dCBM', 'rms_m2d', 'max_m2d', 'rms_d2m', 'max_d2m', 'nspin']

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('dirs', nargs='+')
    ap.add_argument('--json'); ap.add_argument('--tsv')
    a = ap.parse_args()
    res = {}
    for d in a.dirs:
        try: v = check(d)
        except Exception as e: v = dict(error=str(e))
        if v is not None: res[os.path.basename(os.path.normpath(d))] = v
    if a.json: json.dump(res, open(a.json, 'w'), indent=1)
    if a.tsv:
        with open(a.tsv, 'w') as f:
            f.write('name\t' + '\t'.join(KEYS) + '\n')
            for k, v in res.items():
                f.write(k + '\t' + '\t'.join(('%.4f' % v[x]) if isinstance(v.get(x), float) else str(v.get(x, '')) for x in KEYS) + '\n')
    print(f"{'name':14s} {'gap_mesh':>8s} {'gapD':>7s} {'gapM':>7s} {'dVBM':>7s} {'dCBM':>7s} | {'rms_m2d':>7s} {'max':>5s} {'rms_d2m':>7s} {'max':>5s}   (eV)")
    for k, v in res.items():
        if 'error' in v: print(f'{k:14s} ERROR {v["error"]}'); continue
        g = (f"{v['gap_mesh']:8.3f} {v['gapD']:7.3f} {v['gapM']:7.3f} {v['dVBM']:+7.3f} {v['dCBM']:+7.3f}" if v['insulator']
             else f"{'metal':>8s} {'':7s} {'':7s} {'':7s} {'':7s}")
        print(f"{k:14s} {g} | {v['rms_m2d']:7.3f} {v['max_m2d']:5.2f} {v['rms_d2m']:7.3f} {v['max_d2m']:5.2f}")
