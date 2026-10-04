#!/usr/bin/env python3
"""Is a structure a molecular crystal? (2026-10-05, GW1500; user: CO2, CO are molecular, they want empty spheres)

    gw1500_molecular.py <POSCAR> [...]     prints "<file> molecular|extended <largest cluster size> <atoms per cell>"

Two atoms are bonded when their distance is below 1.2 x the sum of their covalent radii (Cordero et al. 2008). The bonded
atoms are followed through the periodic images; the structure is molecular when no cluster reaches its own periodic image
(every cluster is finite) and every atom belongs to a cluster of two or more atoms: an isolated atom means an ionic
crystal whose ions are farther apart than the covalent radii (TlCl, PbF2, BiF3 came out molecular without this).
"""
import sys, itertools
import numpy as np
R = dict(H=0.31, He=0.28, Li=1.28, Be=0.96, B=0.84, C=0.76, N=0.71, O=0.66, F=0.57, Ne=0.58, Na=1.66, Mg=1.41, Al=1.21, Si=1.11,
         P=1.07, S=1.05, Cl=1.02, Ar=1.06, K=2.03, Ca=1.76, Sc=1.70, Ti=1.60, V=1.53, Cr=1.39, Mn=1.39, Fe=1.32, Co=1.26, Ni=1.24,
         Cu=1.32, Zn=1.22, Ga=1.22, Ge=1.20, As=1.19, Se=1.20, Br=1.20, Kr=1.16, Rb=2.20, Sr=1.95, Y=1.90, Zr=1.75, Nb=1.64,
         Mo=1.54, Tc=1.47, Ru=1.46, Rh=1.42, Pd=1.39, Ag=1.45, Cd=1.44, In=1.42, Sn=1.39, Sb=1.39, Te=1.38, I=1.39, Xe=1.40,
         Cs=2.44, Ba=2.15, La=2.07, Hf=1.75, Ta=1.70, W=1.62, Re=1.51, Os=1.44, Ir=1.41, Pt=1.36, Au=1.36, Hg=1.32, Tl=1.45,
         Pb=1.46, Bi=1.48)
def read_poscar(f):
    L = open(f).read().split('\n'); s = float(L[1].split()[0])
    A = np.array([[float(x) for x in L[i].split()[:3]] for i in (2, 3, 4)]) * s
    el = L[5].split(); n = [int(x) for x in L[6].split()]
    k = 7; 
    if L[k].strip()[0] in 'sS': k += 1
    direct = L[k].strip()[0] in 'dD'
    P = np.array([[float(x) for x in L[k + 1 + i].split()[:3]] for i in range(sum(n))])
    if not direct: P = P @ np.linalg.inv(A)
    sp = [e for e, m in zip(el, n) for _ in range(m)]
    return A, P % 1.0, sp
def molecular(f):
    A, P, sp = read_poscar(f); N = len(P)
    bonds = []
    for i in range(N):
        for j in range(N):
            for t in itertools.product((-1, 0, 1), repeat=3):
                if i == j and t == (0, 0, 0): continue
                d = np.linalg.norm((P[j] + t - P[i]) @ A)
                if d < 1.2 * (R.get(sp[i], 1.5) + R.get(sp[j], 1.5)): bonds.append((i, j, np.array(t)))
    # follow the bonds; an atom met again with another cell shift means an infinite cluster
    seen = {}; big = 0; small = 10**9; infinite = False
    for s in range(N):
        if s in seen: continue
        seen[s] = np.zeros(3, int); stack = [s]; size = 1
        while stack:
            i = stack.pop()
            for a, j, t in bonds:
                if a != i: continue
                sh = seen[i] + t
                if j not in seen:
                    seen[j] = sh; stack.append(j); size += 1
                elif not (seen[j] == sh).all():
                    infinite = True
        big = max(big, size); small = min(small, size)
    return ('extended' if infinite or small < 2 else 'molecular'), big, N
for f in sys.argv[1:]:
    try:
        k, b, n = molecular(f); print(f, k, b, n)
    except Exception as e:
        print(f, 'error', type(e).__name__)
