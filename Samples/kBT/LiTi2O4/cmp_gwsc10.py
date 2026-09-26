#!/usr/bin/env python3
"""Compare a run_gwsc10.sh run (one `gwsc 10`, QMLO_* code) with a run_snap.sh chain (`gwsc 1` x 10).

  cmp_gwsc10.py <new run dir> <old chain dir> [window_eV]

Per iteration: ehf of the lmf that closes the iteration (new: llmf.<N>run, old: steps.log) and the
largest change of the QP energies in QPU.<N>run (states within the window around E_F).
Bands: THE MLO band and the conventional sigm band after iteration 10
(new: bnd_mlo_final.dat / bndPMT_final.dat, old: mloband/bnd_iter10.dat / bndPMT_iter10.dat),
rms and max over the bands inside the window along the path.
"""
import sys, re, os
import numpy as np

new, old = sys.argv[1], sys.argv[2]
win = float(sys.argv[3]) if len(sys.argv) > 3 else 3.0

def ehf_new(n):
    f = f'{new}/llmf.{n}run'
    if not os.path.isfile(f): return None
    v = None
    for l in open(f, errors='replace'):
        m = re.search(r'ehf\(eV\)=\s*(-?[\d.]+)', l)
        if m: v = float(m.group(1))
    return v

def ehf_old():
    d = {}
    for l in open(f'{old}/steps.log', errors='replace'):
        m = re.search(r' iter (\d+) .*ehf\(eV\)=\s*(-?[\d.]+)', l)
        if m: d[int(m.group(1))] = float(m.group(2))
    return d

def qpu(f):
    """(q, state) -> (e, dSE) from a QPU file: e = column 'eQP(starting by lmf)', the eigenvalue of the
    lmf that started this iteration (eV from E_F); dSE = dSEnoZ, the QP shift (main_hqpe.sc.f90:231)."""
    rows = {}
    if not os.path.isfile(f): return rows
    for l in open(f, errors='replace'):
        w = l.split()
        if len(w) < 15: continue
        try: v = [float(x) for x in w[:15]]
        except ValueError: continue
        key = (round(v[0], 5), round(v[1], 5), round(v[2], 5), int(v[3]))
        rows[key] = (v[10], v[9])
    return rows

def bands(f):
    """blocks of 'i x E' separated by blank lines -> array (nband, nk) of E and x."""
    blocks, cur = [], []
    for l in open(f):
        if l.startswith('#'): continue
        w = l.split()
        if not w:
            if cur: blocks.append(cur); cur = []
            continue
        cur.append((float(w[1]), float(w[2])))
    if cur: blocks.append(cur)
    nk = min(len(b) for b in blocks)
    x = np.array([p[0] for p in blocks[0][:nk]])
    e = np.array([[p[1] for p in b[:nk]] for b in blocks])
    return x, e

eo = ehf_old()
print(f'# new: {new}\n# old: {old}\n# window for QP energies and bands: |E-EF| < {win} eV')
print(f'{"iter":>4} {"ehf new (eV)":>16} {"ehf old (eV)":>16} {"new-old (meV)":>14} '
      f'{"max|de| (meV)":>14} {"max|dSE| (meV)":>15}   (e, dSE from QPU.<N>run, |e|<window)')
for n in range(1, 11):
    en, eo_n = ehf_new(n), eo.get(n)
    qn, qo = qpu(f'{new}/QPU.{n}run'), qpu(f'{old}/QPU.{n}run')
    common = [k for k in qn if k in qo and abs(qo[k][0]) < win]
    de = max((abs(qn[k][0] - qo[k][0]) for k in common), default=float('nan')) * 1000
    ds = max((abs(qn[k][1] - qo[k][1]) for k in common), default=float('nan')) * 1000
    d = (en - eo_n) * 1000 if (en is not None and eo_n is not None) else float('nan')
    print(f'{n:4d} {en if en is not None else float("nan"):16.6f} {eo_n if eo_n is not None else float("nan"):16.6f} '
          f'{d:14.1f} {de:14.1f} {ds:15.1f}')

for name, fn, fo in [('MLO band', 'bnd_mlo_final.dat', 'mloband/bnd_iter10.dat'),
                     ('sigm band', 'bndPMT_final.dat', 'bndPMT_iter10.dat')]:
    pn, po = f'{new}/{fn}', f'{old}/{fo}'
    if not (os.path.isfile(pn) and os.path.isfile(po)):
        print(f'{name}: missing {pn if not os.path.isfile(pn) else po}'); continue
    xn, en = bands(pn); xo, eo_b = bands(po)
    nb = min(len(en), len(eo_b)); nk = min(en.shape[1], eo_b.shape[1])
    en, eo_b = en[:nb, :nk], eo_b[:nb, :nk]
    mask = (np.abs(eo_b) < win)
    d = (en - eo_b)[mask]
    print(f'{name}: {mask.sum()} points in the window; rms {1000*np.sqrt(np.mean(d**2)):.1f} meV, '
          f'max {1000*np.max(np.abs(d)):.1f} meV')
