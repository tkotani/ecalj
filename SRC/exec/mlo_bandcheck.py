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
Accuracy, with the nearest band at the same k (|dx| < 0.01), in the window [VBM_D - 8, CBM_D + mlo_delta] eV:
  rms_m2d, max_m2d : each MLO point in the window -> nearest DFT band   (a wrong band shows up here)
  rms_d2m, max_d2m : each DFT point in the window -> nearest MLO band   (a missing band shows up here)
Metals: the same with VBM_D = CBM_D = E_F = 0.
mlo_delta is that of [mlo] in ctrlg.*.toml of the directory (2 eV when absent, the default of m_GWinput): it is how far
above the band edge the model is asked to be right, so the model is judged in that same window. Until 2026-10-01 17:54
the upper edge was CBM + 3 (MLO -> DFT) and CBM + 1 (DFT -> MLO), whatever mlo_delta was.

SOC: the MLO bands of job_mlo_soc (2 N_MLO bands per k in band_MLO_spin1.dat, N_MLO from lmlo) must be compared with DFT bands
with SOC; job_mlo_soc draws them itself (step 1b, the same flags; HAM_SO of llmf_band). A mismatch is
reported as 'warning' (2026-10-01 20:1x: Samples/MLOsamples/GaAsSoc looked 0.1 eV off, it was the non-SOC DFT, Delta_SO/3).

Check of a broken model (2026-10-01 21:1x, user: the model can break at a few k points, as Cu and Ni with EH2 on every
atom did, and the rms hides it): FAIL when
  (1) the largest single deviation in the window (both directions) exceeds --fail-max (0.1 eV);
  (2) the deviation of a DFT band (bnd* band index, D -> M) jumps between neighbouring k of a line by more than --fail-jump
      (0.1 eV);
  (3) the MLOs lose linear independence: the smallest eigenvalue of the MLO overlap matrix along the path (MLO_ovlpmin.dat,
      written by mlo) is <= 0 or below 1/100 of its median along the path. (Cu and Ni with EH2: it went negative where the
      bands collapsed; baseline-1 models stay at 3e-3..1e-2.)
The verdict is 'check' (PASS/FAIL) with the reasons in 'fail'.

2026-10-01: rewritten for any directory (the test of Samples/MATERIALS, MD/research_log.md 2026-10-01).
The earlier version (v6.2, effective-mass ratios, fixed to Samples/MLOsamples/*__m3_work) is in the git history.
"""
import argparse, glob, json, os, re, sys
import numpy as np
try:
    import tomllib
except ImportError:
    tomllib = None

RY = 13.6057


def read_ef(d):
    p = os.path.join(d, 'bandplot_MLO.isp1.glt')
    if os.path.exists(p):
        m = re.search(r'ef=\s*([-+0-9.eEdD]+)', open(p).read())
        if m: return float(m.group(1).replace('D', 'E'))
    return float(open(os.path.join(d, 'qplist.dat')).readline().split()[0])


def read_delta(d):
    """mlo_delta (eV) of the ctrlg of the directory; 2.0 (the default of m_GWinput) when absent"""
    for f in glob.glob(os.path.join(d, 'ctrlg.*.toml')):
        if tomllib is None: break
        try: return float(tomllib.load(open(f, 'rb')).get('mlo', {}).get('mlo_delta', 2.0))
        except Exception: pass
    return 2.0


def soc_of(d):
    """(MLO with SOC?, DFT with SOC?), None when it cannot be told. MLO: bands per k of band_MLO_spin1.dat against
    N_MLO of lmlo (2 N with job_mlo_soc); DFT: HAM_SO in llmf_band (the log of job_band), else so of [ham] in the ctrlg."""
    mso = dso = None
    try:
        n = int(re.findall(r'HamRsMTO=\s*(\d+)', open(os.path.join(d, 'lmlo'), errors='replace').read())[-1])
        xs = [l.split()[0] for l in open(os.path.join(d, 'band_MLO_spin1.dat')) if l.strip() and not l.lstrip().startswith('#')]
        mso = xs.count(xs[0]) == 2 * n
    except (OSError, IndexError, ValueError):
        pass
    try:
        m = re.search(r'HAM_SO\s+\S+\s+n=\s*\d+\s+val=\s*([-0-9.]+)', open(os.path.join(d, 'llmf_band'), errors='replace').read())
        if m: dso = float(m.group(1)) > 0
    except OSError:
        pass
    if dso is None and tomllib is not None:          # no log of job_band: the so of the ctrlg it read by default
        for f in glob.glob(os.path.join(d, 'ctrlg.*.toml')):
            try: dso = int(tomllib.load(open(f, 'rb')).get('ham', {}).get('so', 0)) > 0
            except Exception: pass
    return mso, dso


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


def worst_point(P, Q, lo, hi):
    """(deviation, x, E) of the point of P in [lo,hi] farthest from the nearest band of Q at the same x"""
    xs = np.array(sorted(Q)); best = (0.0, None, None)
    for x, es in P.items():
        i = np.searchsorted(xs, x)
        cand = [xs[k] for k in (i - 1, i) if 0 <= k < len(xs)]
        if not cand: continue
        xq = min(cand, key=lambda c: abs(c - x))
        if abs(xq - x) > 0.01: continue
        qe = np.array(Q[xq])
        for e in es:
            if lo <= e <= hi:
                dv = float(np.min(np.abs(qe - e)))
                if dv > best[0]: best = (dv, x, e)
    return best


def band_jump(d, spins, M, lo, hi):
    """largest change of the deviation of one DFT band (index of bnd*) between neighbouring k of a line, both points in
    [lo,hi]: (jump, x, E)"""
    xs = np.array(sorted(M)); best = (0.0, None, None)
    def dev(x, e):
        i = np.searchsorted(xs, x); cand = [xs[k] for k in (i - 1, i) if 0 <= k < len(xs)]
        if not cand: return None
        xq = min(cand, key=lambda c: abs(c - x))
        return float(np.min(np.abs(np.array(M[xq]) - e))) if abs(xq - x) <= 0.01 else None
    for s in spins:
        for f in sorted(glob.glob(os.path.join(d, f'bnd0*.spin{s}'))):
            B = {}
            for l in open(f):
                t = l.split()
                if len(t) < 3 or l.lstrip().startswith('#'): continue
                try: B.setdefault(int(t[0]), []).append((round(float(t[1]), 5), float(t[2])))
                except ValueError: pass
            for pts in B.values():
                prev = None
                for x, e in pts:
                    cur = dev(x, e) if lo <= e <= hi else None
                    if prev is not None and cur is not None and abs(cur - prev) > best[0]: best = (abs(cur - prev), x, e)
                    prev = cur
    return best


def overlap_min(d, soc):
    """(min, x, median) of the smallest eigenvalue of the MLO overlap along the path, from MLO_ovlpmin.dat"""
    f = os.path.join(d, 'MLO_ovlpmin.dat')
    if not os.path.exists(f): return None
    a = np.loadtxt(f, ndmin=2)
    v = a[:, 1] if soc else a[:, 1:].min(axis=1)
    return float(v.min()), float(a[v.argmin(), 0]), float(np.median(v))


def check(d, fail_max=0.1, fail_jump=0.1):
    ef = read_ef(d)
    spins = [s for s in (1, 2) if glob.glob(os.path.join(d, f'bnd*.spin{s}'))]
    mso, dso = soc_of(d)
    if mso: spins = [1]   # SOC: every band is in spin1, a spin2 file is left over from a run without SOC (2026-10-01 20:1x)
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
    r = dict(insulator=ins, gap_mesh=(ln[2] - ln[1]) * RY if ins else 0.0, nspin=len(spins), soc=bool(mso))
    if mso is not None and dso is not None and mso != dso:
        r['warning'] = ('MLO with SOC against DFT without SOC: rerun job_mlo_soc (its step 1b draws the DFT bands with SOC)'
                        if mso else 'MLO without SOC against DFT with SOC: rerun job_band (job_mlo_soc left its SOC bands here)')
    if ins:
        mid = 0.5 * (vbD + cbD)
        vbM = eM[eM <= mid].max() if (eM <= mid).any() else float('nan')
        cbM = eM[eM > mid].min() if (eM > mid).any() else float('nan')
        r.update(gapD=float(cbD - vbD), gapM=float(cbM - vbM), dVBM=float(vbM - vbD), dCBM=float(cbM - cbD))
        top = cbD
    else:
        vbD, top = 0.0, 0.0
    delta = read_delta(d)
    r.update(VBM=float(vbD), CBM=float(top), emin=float(vbD - 8), emax=float(top + delta), mlo_delta=delta)
    m2d = nearest(M, D, vbD - 8, top + delta); d2m = nearest(D, M, vbD - 8, top + delta)
    for k, a in (('m2d', m2d), ('d2m', d2m)):
        r['rms_' + k] = float(np.sqrt((a ** 2).mean())) if a.size else float('nan')
        r['max_' + k] = float(a.max()) if a.size else float('nan')
    # check of a broken model
    lo, hi = vbD - 8, top + delta
    w1 = worst_point(M, D, lo, hi); w2 = worst_point(D, M, lo, hi)
    w = max(w1, w2, key=lambda p: p[0])
    j = band_jump(d, spins, M, lo, hi)
    ov = overlap_min(d, bool(mso))
    r.update(max=w[0], max_x=w[1], max_E=w[2], jump=j[0], jump_x=j[1], jump_E=j[2])
    fail = []
    if w[0] > fail_max: fail.append(f'max {w[0]:.3f} eV at x={w[1]:.3f} E={w[2]:+.2f} (> {fail_max})')
    if j[0] > fail_jump: fail.append(f'jump {j[0]:.3f} eV at x={j[1]:.3f} E={j[2]:+.2f} (> {fail_jump})')
    if ov is not None:
        r.update(ovlp_min=ov[0], ovlp_min_x=ov[1], ovlp_median=ov[2])
        if ov[0] <= 0 or ov[0] < 1e-2 * ov[2]:
            fail.append(f'MLO overlap min {ov[0]:.2e} at x={ov[1]:.3f} (median {ov[2]:.2e}): linear dependence')
    r['check'] = 'FAIL' if fail else 'PASS'
    r['fail'] = fail
    return r


KEYS = ['insulator', 'gap_mesh', 'gapD', 'gapM', 'dVBM', 'dCBM', 'VBM', 'CBM', 'emin', 'emax', 'rms_m2d', 'max_m2d', 'rms_d2m', 'max_d2m', 'nspin',
        'max', 'jump', 'ovlp_min', 'check']

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('dirs', nargs='+')
    ap.add_argument('--json'); ap.add_argument('--tsv')
    ap.add_argument('--fail-max', type=float, default=0.1, help='largest single deviation allowed (eV)')
    ap.add_argument('--fail-jump', type=float, default=0.1, help='largest jump of the deviation of a band between neighbouring k (eV)')
    a = ap.parse_args()
    res = {}
    for d in a.dirs:
        try: v = check(d, a.fail_max, a.fail_jump)
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
        if 'warning' in v: print(f'{k:14s} WARNING {v["warning"]}')
        g = (f"{v['gap_mesh']:8.3f} {v['gapD']:7.3f} {v['gapM']:7.3f} {v['dVBM']:+7.3f} {v['dCBM']:+7.3f}" if v['insulator']
             else f"{'metal':>8s} {'':7s} {'':7s} {'':7s} {'':7s}")
        print(f"{k:14s} {g} | {v['rms_m2d']:7.3f} {v['max_m2d']:5.2f} {v['rms_d2m']:7.3f} {v['max_d2m']:5.2f}")
        ovs = f", overlap min {v['ovlp_min']:.1e}" if 'ovlp_min' in v else ', no MLO_ovlpmin.dat'
        print(f"{'':14s} CHECK {v['check']}: max {v['max']:.3f}, jump {v['jump']:.3f} eV{ovs}" + ''.join('\n' + ' ' * 21 + '- ' + s for s in v['fail']))
