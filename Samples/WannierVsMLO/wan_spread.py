#!/usr/bin/env python3
"""One-off (2026-10-02): spreads of the Wannier functions of genMLWFx with the formula of mlo_spread.py,
M^W(k,b) = dnk(k)^+ UU(k,b) dnk(k+b), from UUU (band basis, window iko) and MLWU (dnk). Also the MV form, to check
against hmaxloc. Usage: wan_spread.py <dir> ; reads BBVEC, UUU, MLWU, PlatQlat.chk there."""
import sys, os, numpy as np
from scipy.io import FortranFile
d = sys.argv[1]; os.chdir(d)
L = open('BBVEC').read().split('\n')
nbb, nqbz = map(int, L[1].split()[:2])
bb = np.array([[float(x) for x in L[2+i].split()[:3]] for i in range(nbb)]); wb = np.array([float(L[2+i].split()[3]) for i in range(nbb)])
alat = float(open('PlatQlat.chk').read().split('\n')[4].split()[0]); u = alat/(2*np.pi)
f = FortranFile('MLWU','r'); nq, nwf, i0, i1 = f.read_ints(np.int32); nb = i1-i0+1
dnk = np.zeros((nq, nb, nwf), complex)
for iq in range(nq):
    f.read_record(np.int32, np.float64) if False else f.read_record('i4,3f8')
    dnk[iq] = f.read_record(np.complex128).reshape((nb, nwf), order='F')
f.close()
f = FortranFile('UUU','r'); f.read_record('S40,i4') if False else f.read_record(np.uint8); hdr = f.read_ints(np.int32)
UU = np.zeros((nqbz, nbb, nb, nb), complex); ikb = np.zeros((nqbz, nbb), int); pair = {}
while True:
    try: tag = f.read_ints(np.int32)
    except Exception: break
    if len(tag) > 1: continue
    tag = tag[0]
    if tag == -10:
        iq, ib, ik = f.read_ints(np.int32); UU[iq-1, ib-1] = f.read_record(np.complex128).reshape((nb, nb), order='F'); ikb[iq-1, ib-1] = ik-1
    elif tag == -20:
        iq, ib, iqt, ibt = f.read_ints(np.int32); pair[(iq-1, ib-1)] = (iqt-1, ibt-1)
for (iq, ib), (iqt, ibt) in pair.items():
    UU[iq, ib] = UU[iqt, ibt].conj().T; ikb[iq, ib] = iqt
M = np.zeros((nqbz, nbb, nwf), complex)
for iq in range(nqbz):
    for ib in range(nbb):
        M[iq, ib] = np.einsum('ni,nm,mi->i', dnk[iq].conj(), UU[iq, ib], dnk[ikb[iq, ib]])
r = -np.einsum('b,ba,kbi->ia', wb, bb, M.imag)/nqbz*u
r2s = 2*np.einsum('b,kbi->i', wb, 1-M.real)/nqbz*u*u
rmv = -np.einsum('b,ba,kbi->ia', wb, bb, np.angle(M))/nqbz*u
r2mv = np.einsum('b,kbi->i', wb, 1-np.abs(M)**2 + np.angle(M)**2)/nqbz*u*u
print(f'{d}: Wannier, {nwf} functions, bands {i0}..{i1}, {nbb} b')
for i in range(nwf):
    print(f'  {i+1}: Omega simple {r2s[i]-np.sum(r[i]**2):8.4f}   Omega MV {r2mv[i]-np.sum(rmv[i]**2):8.4f}  bohr^2')
