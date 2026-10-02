#!/usr/bin/env python3
"""mlo_maxloc.py: maximal localization within the MLO subspace (prototype, 2026-10-02; MD/TODOandQuestion.md, user: "MLO を作った後に
Marzari の方法で最大局在化させることはありえる。バンドは変えない。Sakuma の方法のように対称性をレスペクトしないといけない").

    mlo_maxloc.py <sname> [--niter 400] [--isp 1|2] [--orb 5,6,7,8,9] [--sym [--lblocks 0,1,2]]

Run it where mlo_spread.py has run (BBVEC, UUU.* / UUD.*, QBZ.chk, PlatQlat.chk, HamRsMLO). It does not change any file of the
MLO model; it reports the spreads of
  (a) the MLOs as they are (not orthonormal; eq. (2) of mlo_spread.py),
  (b) the Loewdin orthonormalized MLOs  F~(k) = F(k) O(k)^{-1/2},  M~(k,b) = O(k)^{-1/2} M(k,b) O(k+b)^{-1/2},
  (c) (b) followed by the steepest descent of Marzari and Vanderbilt (PRB 56, 12847) over unitary U(k):
        M(k,b) -> U(k)^+ M(k,b) U(k+b),   dOmega/dW(k) = 4 sum_b w_b (A[R] - S[T])                          (MV eq. 52)
        R_mn = M_mn M_nn^*,  T_mn = (M_mn / M_nn) q_n,  q_n = Im ln M_nn + b.r_n,  A[X] = (X - X^+)/2,  S[X] = (X + X^+)/2i
      with the spread  Omega = sum_n [ (1/N) sum_kb w_b (1 - |M_nn|^2 + (Im ln M_nn)^2) - r_n^2 ],
      r_n = -(1/N) sum_kb w_b b Im ln M_nn                                                                   (MV eqs. 31, 32)
(2026-10-02, later the same day: the MLOs are now the Loewdin orthonormalized ones, and so are the M(k,b) of mlo_spread.py;
(a) and (b) then coincide. The spread of the raw MLOs needs a model made with job_mlo --mlo_raw.)
The unitary mixing stays in the MLO subspace at each k, so the interpolated bands and the cRPA weights p_kn do not change.
Without --sym there is no symmetry constraint; the result may break the site symmetry (Fe spd: s, p and eg mix into
hybrids off the atom). --sym (06:30) keeps it, as Sakuma (PRB 87, 235109) does: the Loewdin MLOs transform under the point
group like the orbitals, M~(gk,gb) = X(g) M~(k,b) X(g)^+ with X(g) the rotation of the real harmonics (rotdlmm of libecaljF),
and the gauge is kept in the form U(gk) = X(g) U(k) X(g)^+ by averaging the gradient over the group,
    G_s(k) = (1/N_g) sum_g X(g)^+ G(gk) X(g).
The relation of M~ is checked first (it fixes whether X is D or D^T and catches a wrong k map). --sym handles one atom at the
origin and a symmorphic group (Ni, Fe); the MLOs ordered by l (--lblocks, default from their number: 5 d, 9 s p d).
--orb restricts (b),(c) to a subset of MLOs (1-based): their own subspace, e.g. the d of one atom. Without it all MLOs.
"""
import argparse, sys
from pathlib import Path
import numpy as np
from scipy.linalg import expm
EXEC = Path(__file__).resolve().parents[2] / 'SRC' / 'exec'   # mlo_spread.py, symfind.py (2026-10-02: moved to TOOLS/gadget)
sys.path.insert(0, str(EXEC))
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


def symmetry_X(sname, plat, nbasis, lblocks):
    """X(g) (nbasis x nbasis) for the point-group operations of the crystal: blocks of the real-harmonic rotation matrices
    D^l(R) of ecalj (rotdlmm), one block per l in lblocks. One atom at the origin, no translations."""
    import ctypes, importlib.util, os
    spec = importlib.util.spec_from_file_location('symfind', str(EXEC / 'symfind.py'))
    sf = importlib.util.module_from_spec(spec); spec.loader.exec_module(sf)
    alat_c, plat_c, names, frac, af = sf.read_ctrlg(sname)
    if len(names) != 1:
        sys.exit('mlo_maxloc --sym: one atom per cell only (prototype)')
    _, get = sf.spglib_dataset(alat_c, plat_c, names, frac, 1e-5)
    R, T = np.array(get('rotations')), np.array(get('translations'))
    if np.abs(T - np.rint(T)).max() > 1e-8:
        sys.exit('mlo_maxloc --sym: a non-symmorphic group (prototype)')
    Rc = np.array([plat_c @ r @ np.linalg.inv(plat_c) for r in R])        # Cartesian rotations
    if sum(2 * l + 1 for l in lblocks) != nbasis:
        sys.exit(f'mlo_maxloc --sym: lblocks {lblocks} do not make {nbasis} MLOs')
    lib = None
    import shutil
    lmf = shutil.which('lmf')   # the bindir (~/bin) holds lmf and a link to libecaljF.so
    for d in ([Path(lmf).parent] if lmf else []) + [EXEC.parent / 'build_gfortran']:
        if (d / 'libecaljF.so').exists() and (d / 'lmf').exists():
            # the libraries lmf is linked with (MKL's FFTW and LAPACK, MPI) first, global: libecaljF.so leaves them to the executable
            import subprocess
            for line in subprocess.run(['ldd', str(d / 'lmf')], capture_output=True, text=True).stdout.split('\n'):
                w = line.split()
                if len(w) >= 3 and w[1] == '=>' and Path(w[2]).exists() and 'libecaljF' not in w[2]:
                    ctypes.CDLL(w[2], mode=ctypes.RTLD_GLOBAL)
            lib = ctypes.CDLL(str(d / 'libecaljF.so')); break
    if lib is None:
        sys.exit('mlo_maxloc --sym: libecaljF.so and lmf not found (lmf on PATH, or SRC/build_gfortran)')
    ng, nl = len(Rc), max(lblocks) + 1
    sym = np.asfortranarray(np.transpose(Rc, (1, 2, 0)))                   # symops(3,3,ng)
    dl = np.zeros((2 * nl - 1, 2 * nl - 1, nl, ng), order='F')
    lib.__m_symderive_MOD_rotdlmm(sym.ctypes.data_as(ctypes.c_void_p), ctypes.byref(ctypes.c_int(ng)),
                                  ctypes.byref(ctypes.c_int(nl)), dl.ctypes.data_as(ctypes.c_void_p))
    X = np.zeros((ng, nbasis, nbasis))
    for ig in range(ng):
        o = 0
        for l in lblocks:
            X[ig, o:o + 2 * l + 1, o:o + 2 * l + 1] = dl[nl - 1 - l:nl + l, nl - 1 - l:nl + l, l, ig]
            o += 2 * l + 1
    return Rc, X


def kb_maps(Rc, qbz, qlat, bb):
    """index of g k on the mesh and of g b in the b list, for each operation."""
    qinv = np.linalg.inv(qlat)
    key = lambda q: tuple(np.round((qinv @ q) % 1.0, 6) % 1.0)
    kidx = {key(k): i for i, k in enumerate(qbz)}
    gk = np.array([[kidx[key(r @ k)] for k in qbz] for r in Rc])
    gb = np.array([[int(np.argmin(np.linalg.norm(bb - r @ b, axis=1))) for b in bb] for r in Rc])
    for r, row in zip(Rc, gb):
        if np.abs(bb[row] - (r @ bb.T).T).max() > 1e-6:
            sys.exit('mlo_maxloc --sym: the b shells are not closed under the group')
    return gk, gb


def main():
    p = argparse.ArgumentParser(prog='mlo_maxloc.py', description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('sname')
    p.add_argument('--niter', type=int, default=400)
    p.add_argument('--isp', type=int, default=1)
    p.add_argument('--orb', default='', help='subset of MLOs (1-based, comma separated)')
    p.add_argument('--sym', action='store_true', help='keep the point-group symmetry (Sakuma)')
    p.add_argument('--lblocks', default='', help='l of the MLO blocks in order, e.g. 0,1,2 (with --sym)')
    p.add_argument('--save-u', default='', help='write qbz, the overlap O(k) and U(k) to this npz (2026-10-02 10:00; mlo_cmlo_transform.py applies them)')
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
    if a.sym:
        if a.orb:
            sys.exit('mlo_maxloc: --sym with all MLOs only (prototype)')
        lbl = [int(x) for x in a.lblocks.split(',')] if a.lblocks else {5: [2], 9: [0, 1, 2], 4: [0, 1], 1: [0]}[nbasis]
        Rc, Xd = symmetry_X(a.sname, plat, nbasis, lbl)
        gk, gb = kb_maps(Rc, qbz, qlat, bbx)
        best = None
        for name, Xg in (('D', Xd), ('D^T', np.transpose(Xd, (0, 2, 1)))):
            err = max(np.abs(M0[gk[ig]][:, gb[ig]] - np.einsum('ij,kbjl,ml->kbim', Xg[ig], M0, Xg[ig])).max() for ig in range(len(Rc)))
            if best is None or err < best[1]:
                best = (name, err, Xg)
        print(f'mlo_maxloc --sym: {len(Rc)} operations, X = {best[0]}, max |M(gk,gb) - X M(k,b) X^+| = {best[1]:.2e}')
        if best[1] > 1e-4:
            sys.exit('mlo_maxloc --sym: the MLOs do not transform as assumed (lblocks, origin, or the k map)')
        Xs = best[2]
        def symgrad(G):
            return np.mean([np.einsum('ji,kjl,lm->kim', Xs[ig], G[gk[ig]], Xs[ig]) for ig in range(len(Rc))], axis=0)
    else:
        symgrad = lambda G: G
    U = np.tile(np.eye(nw, dtype=complex), (len(qbz), 1, 1))
    Mc = M0.copy()
    om, r = spread_mv(Mc, bbx, wbx, u)
    tot = om.sum()
    step = 0.25 / (4 * wbx.sum())
    G = symgrad(gradient(Mc, bbx, wbx, r / u))
    sgn = 1.0                                   # the direction that lowers Omega (conventions differ between papers)
    for s in (1.0, -1.0):
        Ut = np.array([U[k] @ expm(s * 1e-3 * step * G[k]) for k in range(len(qbz))])
        if spread_mv(rotate(M0, Ut, kbx), bbx, wbx, u)[0].sum() < tot:
            sgn = s
            break
    for it in range(a.niter):
        G = sgn * symgrad(gradient(Mc, bbx, wbx, r / u))
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
    if a.save_u:
        np.savez(a.save_u, qbz=qbz, O=Os, U=U, isp=a.isp, orb=np.array(sel), sym=a.sym)
        print(f'   wrote {a.save_u}: qbz, O(k), U(k) (U from the Loewdin MLOs; C -> C O^-1/2 U)')


if __name__ == '__main__':
    main()
