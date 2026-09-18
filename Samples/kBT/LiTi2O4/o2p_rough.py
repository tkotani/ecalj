#!/usr/bin/env python3
"""Roughness of bands along the k path from job_band's bnd*.spin1 files.

rough(ib) = mean over interior path points of |E(i-1) - 2 E(i) + E(i+1)|  [meV],
segment by segment (bnd001..bnd006 are the syml segments).  Bands 9-12 are the
O 2p bottom of LiTi2O4 (about -8 eV); bands 1-8 are Li 1s / O 2s / semicore.
usage: o2p_rough.py DIR [DIR ...]   (each DIR holds bnd00?.spin1)
"""
import sys, glob, os
from collections import defaultdict

def read(dirn):
    bands = defaultdict(list)          # ib -> list of segments, each list of E
    for f in sorted(glob.glob(os.path.join(dirn, 'bnd0*.spin1'))):
        seg = defaultdict(list)
        for line in open(f):
            if line.startswith('#') or not line.strip():
                continue
            t = line.split()
            seg[int(t[0])].append(float(t[2]))
        for ib, e in seg.items():
            bands[ib].append(e)
    return bands

def rough(bands, ib):
    tot, n = 0.0, 0
    for e in bands[ib]:
        for i in range(1, len(e) - 1):
            tot += abs(e[i-1] - 2*e[i] + e[i+1]); n += 1
    return 1e3*tot/n if n else float('nan')

if __name__ == '__main__':
    lo, hi = 9, 12
    print('# dir  rough(meV) bands %d-%d : mean  max   [per band]   Emin(b9)  Emax(b12)' % (lo, hi))
    for d in sys.argv[1:]:
        b = read(d)
        if not b: print(d, 'no bands'); continue
        r = [rough(b, ib) for ib in range(lo, hi+1)]
        emin = min(min(s) for s in b[lo]); emax = max(max(s) for s in b[hi])
        print('%-40s %7.1f %7.1f   [%s]  %8.3f %8.3f' % (
            os.path.basename(d.rstrip('/')) or d, sum(r)/len(r), max(r),
            ' '.join('%.1f' % x for x in r), emin, emax))
