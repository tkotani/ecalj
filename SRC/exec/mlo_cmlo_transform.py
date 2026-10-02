#!/usr/bin/env python3
"""mlo_cmlo_transform.py: change the basis of the MLO subspace in __cmlo.data (prototype, 2026-10-02).

    mlo_cmlo_transform.py lowdin                      C(k) -> C(k) O(k)^-1/2,           O = C^+ C
    mlo_cmlo_transform.py mv  U_isp1.npz [U_isp2.npz]  C(k) -> C(k) O(k)^-1/2 U(k)     (U of mlo_maxloc.py --save-u)

Run it in a COPY of a directory where the MLO model and its W were made (job_mloW, job_mlo_magnon): it rewrites
__cmlo.data in place (the old one is kept as __cmlo.data.orig). The columns of C are the MLOs on the bands
(|F_j(k)> = sum_n |psi_kn> C_nj); the new columns are orthonormal (Loewdin), and with mv also maximally localized
within the subspace (mlo_maxloc.py, --sym keeps the point group). The bands of the model do not change; quantities taken
on the MLOs (hwmatK: v, W, their on-site and off-site elements; mlo_magnon) then refer to the new functions.
Only the irreducible q of __cmlo.data are stored; the others are rotated from them by the symmetry operations, which is
consistent with U only if U(gk) = X(g) U(k) X(g)^+ (mlo_maxloc.py --sym) and with the Loewdin step always.
__cmlo.info: (nmlo, nqbz, nqirr, nMTO, mrecbb), (ix, qplistgw); __cmlo.data: direct access, record isp + nspx (iq - 1)
holding C(nbandmx, nmlo) (m_mlo_wfs.f90).
"""
import sys, os, shutil
import numpy as np
from scipy.io import FortranFile


def read_info():
    f = FortranFile('__cmlo.info', 'r')
    nmlo, nqbz, nqirr, nmto, mrecbb = f.read_ints(np.int32)
    ix, q = f.read_record(np.dtype(('i4', (nmlo,))), np.dtype(('f8', (nqirr, 3))))
    f.close()
    return int(nmlo), int(nqirr), int(mrecbb), np.array(q).reshape(nqirr, 3)


def invsqrt(O):
    w, v = np.linalg.eigh(O)
    return (v / np.sqrt(w)) @ v.conj().T


def main():
    if len(sys.argv) < 2 or sys.argv[1] not in ('lowdin', 'mv'):
        sys.exit(__doc__)
    mode = sys.argv[1]
    nmlo, nqirr, mrecbb, qirr = read_info()
    nbmx = mrecbb // (16 * nmlo)
    raw = np.fromfile('__cmlo.data' if not os.path.exists('__cmlo.data.orig') else '__cmlo.data.orig', dtype=np.complex128)
    nrec = raw.size // (nbmx * nmlo)
    nspx = nrec // nqirr
    C = raw.reshape(nrec, nmlo, nbmx).transpose(0, 2, 1)           # C[irec, band, mlo] (Fortran order in the record)
    if not os.path.exists('__cmlo.data.orig'):
        shutil.copy('__cmlo.data', '__cmlo.data.orig')
    U = {}
    if mode == 'mv':
        L = open('PlatQlat.chk').read().split('\n')
        qlat = np.array([[float(x) for x in L[i].split()[3:6]] for i in range(3)]).T
        qinv = np.linalg.inv(qlat)
        key = lambda q: tuple(np.round((qinv @ q) % 1.0, 6) % 1.0)
        for isp, fn in enumerate(sys.argv[2:4], 1):
            d = np.load(fn)
            if int(d['isp']) != isp or len(d['orb']) != nmlo:
                sys.exit(f'{fn}: isp or the MLOs do not match')
            kidx = {key(k): i for i, k in enumerate(d['qbz'])}
            U[isp] = (d['U'], d['O'], kidx, key)
    out = C.copy()
    dev = 0.0
    near = 0
    for iq in range(nqirr):
        for isp in range(1, nspx + 1):
            irec = isp - 1 + nspx * iq
            c = C[irec]
            O = c.conj().T @ c
            X = invsqrt(O)
            if mode == 'mv':
                Uk, Ok, kidx, key = U[isp if isp in U else 1]
                ik = kidx.get(key(qirr[iq]))
                if ik is None:   # the small q0 near Gamma of the GW q list: U of the nearest mesh point (U is smooth in k)
                    fr = np.array(key(qirr[iq])); fr = fr - np.rint(fr)
                    ks = np.array([np.array(kk) - np.rint(kk) for kk in kidx])
                    ik = list(kidx.values())[int(np.argmin(np.linalg.norm(ks - fr, axis=1)))]
                    near += 1
                dev = max(dev, np.abs(Ok[ik] - O).max())                # the O of huumat (UUU, b=0) against C^+ C
                X = X @ Uk[ik]
            out[irec] = c @ X
            assert np.abs(out[irec].conj().T @ out[irec] - np.eye(nmlo)).max() < 1e-8
    out.transpose(0, 2, 1).astype(np.complex128).tofile('__cmlo.data')
    msg = f' (max |O_huumat - C^+C| = {dev:.1e}; {near} q off the mesh took U of the nearest mesh point)' if mode == 'mv' else ''
    print(f'mlo_cmlo_transform: {mode}, {nqirr} irreducible q x {nspx} spins, {nmlo} MLOs, {nbmx} bands per record{msg}')


if __name__ == '__main__':
    main()
