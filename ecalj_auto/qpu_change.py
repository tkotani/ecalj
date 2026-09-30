#!/usr/bin/env python3
# 2026-09-30: change of the QSGW eigenvalues between iterations, for materials whose gap cannot be used (metals).
# Reads QPU.<n>run in each directory: the column "eQP(starting by lmf)" (index 10) is the lmf eigenvalue of that
# iteration, relative to E_F (eV). For each pair of iterations n-1 -> n prints max and rms |change| over the states
# with |e| < WIN eV (both spins, all q of QPU).
import sys, os, re, glob
WIN = 5.0

def read_qpu(p):
    d = {}; isp = 0
    for l in open(p):
        if 'isp=' in l:
            isp = int(l.split('isp=')[1].split()[0]); continue
        f = l.split()
        if len(f) >= 11 and re.match(r'^-?\d+\.\d+$', f[0]) and f[3].isdigit():
            try: d[(isp, f[0], f[1], f[2], int(f[3]))] = float(f[10])
            except ValueError: pass
    return d

for dd in sys.argv[1:]:
    runs = sorted(int(re.search(r'QPU\.(\d+)run', p).group(1)) for p in glob.glob(os.path.join(dd, 'QPU.*run')))
    out = []
    prev = None
    for n in runs:
        cur = read_qpu(os.path.join(dd, f'QPU.{n}run'))
        if prev is not None:
            ch = [abs(cur[k] - prev[k]) for k in cur if k in prev and abs(cur[k]) < WIN]
            if ch:
                rms = (sum(c * c for c in ch) / len(ch)) ** 0.5
                out.append(f'{n-1}->{n}: max {max(ch):.3f} rms {rms:.3f}')
        prev = cur
    print(os.path.basename(dd.rstrip('/')), f'iters {len(runs)} |', ' ; '.join(out))
