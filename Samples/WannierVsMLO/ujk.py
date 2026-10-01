#!/usr/bin/env python3
"""U, U', J at omega=0, R=0 from Coulomb_v.UP, Screening_W-v.UP, Screening_W-v_crpa.UP (hwmatK_MPI). eV.
Element (i1 i2|i3 i4): U = mean (ii|ii), U' = mean_{i!=j} (ii|jj), J = mean_{i!=j} (ij|ji).  (2026-10-01)"""
import sys, os, numpy as np
def read(f, nre):
    d = {}
    for l in open(f):
        w = l.split()
        if not w or w[0] != 'Wannier': continue
        R = tuple(float(x) for x in w[3:6])
        if any(abs(x) > 1e-6 for x in R): continue
        i = tuple(int(x) for x in w[7:11])
        if nre == 4 and abs(float(w[11])) > 1e-8: continue     # omega = 0 only
        if i not in d: d[i] = float(w[-2])
    return d
def ujj(d):
    n = max(k[0] for k in d)
    U = np.mean([d[(i,i,i,i)] for i in range(1,n+1)])
    Up = np.mean([d[(i,i,j,j)] for i in range(1,n+1) for j in range(1,n+1) if i!=j])
    J = np.mean([d[(i,j,j,i)] for i in range(1,n+1) for j in range(1,n+1) if i!=j])
    return n, U, Up, J
for dd in sys.argv[1:]:
    v = read(os.path.join(dd,'Coulomb_v.UP'), 2)
    out = [('v', ujj(v))]
    for tag, f in (('W(RPA)','Screening_W-v.UP'), ('W(cRPA)','Screening_W-v_crpa.UP')):
        p = os.path.join(dd, f)
        if os.path.exists(p):
            w = read(p, 4); out.append((tag, ujj({k: v[k]+w[k] for k in w})))
    for tag, (n, U, Up, J) in out:
        print(f'{dd:12s} {tag:8s} n={n}  U={U:8.3f}  Up={Up:8.3f}  J={J:7.3f}')
