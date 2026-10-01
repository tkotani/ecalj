#!/usr/bin/env python3
"""mlo_maxloc.py: maximal localization within the MLO subspace (prototype, 2026-10-02; MD/TODOandQuestion.md, user: "MLO を作った後に
Marzari の方法で最大局在化させることはありえる。バンドは変えない。Sakuma の方法のように対称性をレスペクトしないといけない").

    mlo_maxloc.py <sname> [--niter 400] [--isp 1|2] [--orb 5,6,7,8,9]

Run it where mlo_spread.py has run (BBVEC, UUU.* / UUD.*, QBZ.chk, PlatQlat.chk, HamRsMLO). It does not change any file of the
MLO model; it reports the spreads of
  (a) the MLOs as they are (not orthonormal; eq. (2) of mlo_spread.py),
  (b) the Loewdin orthonormalized MLOs  F~(k) = F(k) O(k)^{-1/2},  M~(k,b) = O(k)^{-1/2} M(k,b) O(k+b)^{-1/2},
  (c) (b) followed by the steepest descent of Marzari and Vanderbilt (PRB 56, 12847) over unitary U(k):
        M(k,b) -> U(k)^+ M(k,b) U(k+b),   dOmega/dW(k) = 4 sum_b w_b (A[R] - S[T])                          (MV eq. 52)
        R_mn = M_mn M_nn^*,  T_mn = (M_mn / M_nn) q_n,  q_n = Im ln M_nn + b.r_n,  A[X] = (X - X^+)/2,  S[X] = (X + X^+)/2i
      with the spread  Omega = sum_n [ (1/N) sum_kb w_b (1 - |M_nn|^2 + (Im ln M_nn)^2) - r_n^2 ],
      r_n = -(1/N) sum_kb w_b b Im ln M_nn                                                                   (MV eqs. 31, 32)
The unitary mixing stays in the MLO subspace at each k, so the interpolated bands and the cRPA weights p_kn do not change.
No symmetry constraint yet (Sakuma, PRB 87, 235109: U(gk) = D(g) U(k) d(g)^-1 is the next step); the result may break the
site symmetry slightly, which the report shows as the spread of Omega within a shell (t2g, eg).
--orb restricts (b),(c) to a subset of MLOs (1-based): their own subspace, e.g. the d of one atom. Without it all MLOs.
"""
import argparse, sys
from pathlib import Path
import numpy as np
from scipy.linalg import expm
sys.path.insert(0, str(Path(__file__).parent))
from mlo_spread import read_lattice, read_qbz, read_uu
from scipy.io import FortranFile


def read_bbvec():
    """b vectors (2pi/alat), weights, and the index of k+b for each (k,b), from BBVEC of mlo_spread.py (the last b is 0)."""
    L = open('BBVEC').read().split('\n')
    nbb, nqbz = map(int, L[1].split())
    bb = np.array([[float(x) for x in L[2 + i].split()[:3]] for i in range(nbb)])
    wb = np.array([float(L[2 + i].split()[3]) for i in range(nbb)])
    kb = np.zeros((nqbz, nbb), int)
    p = 2 + nbb
    for iq in range(nqbz):
        p += 1
        for ib in range(nbb):
            kb[iq, ib] = int(L[p].split()[2]) - 1
            p += 1
    return bb, wb, kb


def invsqrt(O):
    w, v = np.linalg.eigh(O)
    return (v / np.sqrt(w)) @ v.conj().T


def spread_mv(M, bb, wb, u):
    """MV spread of orthonormal functions: Omega per function (bohr^2) and centres (bohr). M[k,b,m,n] without b=0."""
    Mnn = np.einsum('kbnn->kbn', M)
    ph = np.angle(Mnn)                                           # Im ln M_nn
    N = M.shape[0]
    r = -np.einsum('b,ba,kbn->na', wb, bb, ph) / N               # 2pi/alat^-1 = alat/2pi units
    r2 = np.einsum('b,kbn->n', wb, 1 - np.abs(Mnn) ** 2 + ph ** 2) / N
    om = r2 - np.sum(r * r, axis=1)
    return om * u * u, r * u


def gradient(M, bb, wb, rn):
    Mnn = np.einsum('kbnn->kbn', M)
    q = np.angle(Mnn) + np.einsum('ba,na->bn', bb, rn)[None]
    R = M * Mnn.conj()[:, :, None, :]
    T = M / Mnn[:, :, None, :] * q[:, :, None, :]
    A = lambda X: (X - np.conj(np.swapaxes(X, -1, -2))) / 2
    S = lambda X: (X + np.conj(np.swapaxes(X, -1, -2))) / 2j
    return 4 * np.einsum('b,kbmn->kmn', wb, A(R) - S(T))


def rotate(M0, U, kb):
    """U(k)^+ M0(k,b) U(k+b)."""
    return np.einsum('kmi,kbmn,kbnj->kbij', U.conj(), M0, U[kb])


def main():
    p = argparse.ArgumentParser(prog='mlo_maxloc.py', description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('sname')
    p.add_argument('--niter', type=int, default=400)
    p.add_argument('--isp', type=int, default=1)
    p.add_argument('--orb', default='', help='subset of MLOs (1-based, comma separated)')
    a = p.parse_args()
    plat, qlat, alat = read_lattice()
    qbz = read_qbz()
    bb, wb, kb = read_bbvec()
    f = FortranFile('HamRsMLO', 'r'); nbasis, _, nspx = f.read_ints(np.int32); f.close()
    M = read_uu(len(qbz), len(bb), nbasis, ['UUU.', 'UUD.'][a.isp - 1])
    u = alat / (2 * np.pi)
    O = M[:, -1]                                                  # b = 0: the overlap O(k)
    Mb, bbx, wbx, kbx = M[:, :-1], bb[:-1], wb[:-1], kb[:, :-1]
    sel = [int(x) - 1 for x in a.orb.split(',')] if a.orb else list(range(nbasis))
    # (a) raw MLOs, mlo_spread.py eq. (2)
    Od = np.real(np.einsum('kii->ki', O))
    Md = np.einsum('kbii->kbi', Mb)
    r_a = -np.einsum('b,ba,kbi->ia', wbx, bbx, Md.imag) / len(qbz) * u
    om_a = 2 * np.einsum('b,kbi->i', wbx, Od[:, None, :] - Md.real) / len(qbz) * u * u - np.sum(r_a * r_a, axis=1)
    # (b) Loewdin within the subset
    Os = O[:, sel][:, :, sel]
    X = np.array([invsqrt(Os[k]) for k in range(len(qbz))])
    M0 = np.einsum('kmi,kbmn,kbnj->kbij', X.conj(), Mb[:, :, sel][:, :, :, sel], X[kbx])
    om_b, r_b = spread_mv(M0, bbx, wbx, u)
    # (c) MV steepest descent
    nw = len(sel)
    U = np.tile(np.eye(nw, dtype=complex), (len(qbz), 1, 1))
    Mc = M0.copy()
    om, r = spread_mv(Mc, bbx, wbx, u)
    tot = om.sum()
    step = 0.25 / (4 * wbx.sum())
    G = gradient(Mc, bbx, wbx, r / u)
    sgn = 1.0                                   # the direction that lowers Omega (conventions differ between papers)
    for s in (1.0, -1.0):
        Ut = np.array([U[k] @ expm(s * 1e-3 * step * G[k]) for k in range(len(qbz))])
        if spread_mv(rotate(M0, Ut, kbx), bbx, wbx, u)[0].sum() < tot:
            sgn = s
            break
    for it in range(a.niter):
        G = sgn * gradient(Mc, bbx, wbx, r / u)
        while True:
            Ut = np.array([U[k] @ expm(step * G[k]) for k in range(len(qbz))])
            Mt = rotate(M0, Ut, kbx)
            omt, rt = spread_mv(Mt, bbx, wbx, u)
            if omt.sum() < tot or step < 1e-8:
                break
            step /= 2
        if omt.sum() >= tot:
            break
        dt = tot - omt.sum()
        U, Mc, om, r, tot = Ut, Mt, omt, rt, omt.sum()
        step *= 1.2
        if dt < 1e-7:
            break
    print(f'mlo_maxloc: {a.sname} isp={a.isp}, {len(qbz)} k, {len(bbx)} b, MLOs {[s + 1 for s in sel]}; MV steps {it + 1}')
    print('   MLO   Omega(bohr^2): raw MLO   Loewdin   Loewdin+MV      centre after MV (bohr)')
    for i, s in enumerate(sel):
        print(f'  {s + 1:4d}          {om_a[s]:10.4f} {om_b[i]:9.4f} {om[i]:10.4f}     {r[i, 0]:8.4f} {r[i, 1]:8.4f} {r[i, 2]:8.4f}')
    print(f'   sum          {om_a[sel].sum():10.4f} {om_b.sum():9.4f} {om.sum():10.4f}')
    # gauge-invariant part (MV eq. 34): Omega_I = (1/N) sum_kb w_b (n_w - sum_mn |M_mn|^2); no unitary mixing lowers it
    omI = np.einsum('b,kb->', wbx, nw - np.sum(np.abs(M0) ** 2, axis=(2, 3))) / len(qbz) * u * u
    print(f'   Omega_I (gauge invariant) = {omI:.4f} bohr^2;  Omega - Omega_I: Loewdin {om_b.sum() - omI:.4f}, Loewdin+MV {om.sum() - omI:.4f}')


if __name__ == '__main__':
    main()
