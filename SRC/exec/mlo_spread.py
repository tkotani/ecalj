#!/usr/bin/env python3
"""mlo_spread.py: the spread of the MLOs in real space, <r>, <r^2> and Omega = <r^2> - <r>^2 of each MLO.

    mlo_spread.py <sname> [-np N] [--norun]

Run it after job_mloW (it needs __cmlo.data/.info of `mlo --mlo`, the eigenfunctions of `lmf --jobgw=1 --mlo`,
HamRsMLO, QBZ.chk and PlatQlat.chk). It
  1. writes BBVEC: the vectors b that join each k of the GW mesh to its neighbours, with weights w_b such that
     sum_b w_b b_a b_b = delta_ab (the shells nearest first, as Wannier90 does; Marzari-Vanderbilt, PRB 56, 12847,
     eq. B1), and b = 0 with weight 0 for <u_k|u_k>;
  2. runs `huumat --dwnb=mlo --job=2`: M_ij(k,b) = <u_ik|u_j,k+b> of the MLO Bloch functions (files UUU.*, UUD.*);
  3. prints, for each MLO i,
        <r>_i   = -(1/N_k) sum_k sum_b w_b b Im M_ii(k,b)                                           (1)
        <r^2>_i =  (2/N_k) sum_k sum_b w_b [ O_ii(k) - Re M_ii(k,b) ],   O_ii(k) = M_ii(k,0)          (2)
        Omega_i = <r^2>_i - |<r>_i|^2
     in bohr. Eq. (2) is |u(k+b) - u(k)|^2 summed; it does not assume <u_ik|u_ik> = 1, which the MLOs, normalized
     in real space (one constant per orbital), do not satisfy at each k. It equals the second moment of the
     real-space MLO to O(b^2): the finer the GW mesh, the closer. (2026-10-02)
"""
import argparse, os, shutil, subprocess, sys, tomllib
from pathlib import Path
import numpy as np
from scipy.io import FortranFile

EXEC_DIR = Path(__file__).parent


def read_lattice():
    """plat, qlat (columns, alat and 2pi/alat units) and alat (bohr) from PlatQlat.chk."""
    L = open('PlatQlat.chk').read().split('\n')
    plat = np.array([[float(x) for x in L[i].split()[0:3]] for i in range(3)]).T
    qlat = np.array([[float(x) for x in L[i].split()[3:6]] for i in range(3)]).T
    alat = float(L[4].split()[0])
    return plat, qlat, alat


def read_qbz():
    L = open('QBZ.chk').read().replace('D', 'E').split('\n')
    n = int(L[0].split()[0])
    return np.array([[float(x) for x in L[1 + i].split()[:3]] for i in range(n)])


def bshells(qlat, n123, tol=1e-6):
    """b vectors and weights: add the shells of mesh neighbours, nearest first, until sum_b w_b b b^T = 1."""
    cand = []
    for i in range(-3, 4):
        for j in range(-3, 4):
            for k in range(-3, 4):
                if (i, j, k) == (0, 0, 0):
                    continue
                cand.append(qlat @ np.array([i / n123[0], j / n123[1], k / n123[2]]))
    cand = sorted(cand, key=lambda b: np.linalg.norm(b))
    shells, cur, r0 = [], [], None
    for b in cand:
        r = np.linalg.norm(b)
        if r0 is not None and abs(r - r0) > tol * max(1, r0):
            shells.append(cur); cur = []
        cur.append(b); r0 = r
    shells.append(cur)
    target = np.array([1, 0, 0, 1, 0, 1], dtype=float)          # xx xy xz yy yz zz of the identity
    use = []
    for s in shells:
        trial = use + [s]
        A = np.array([[sum(b[a] * b[c] for b in sh) for sh in trial] for a, c in ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))])
        w, *_ = np.linalg.lstsq(A, target, rcond=None)
        res = np.linalg.norm(A @ w - target)
        if np.linalg.matrix_rank(A, tol=1e-8) > len(use):         # the shell adds a new direction
            use = trial
            if res < 1e-8:
                bb, wb = [], []
                for sh, ws in zip(use, w):
                    for b in sh:
                        bb.append(b); wb.append(ws)
                return np.array(bb), np.array(wb)
    sys.exit('mlo_spread: no set of b shells satisfies sum_b w_b b b^T = 1')


def kindex(qbz, qlat):
    """index of k + b on the mesh, by the fractional coordinates modulo 1."""
    qinv = np.linalg.inv(qlat)
    frac = {tuple(np.round((qinv @ k) % 1.0, 6) % 1.0): i for i, k in enumerate(qbz)}
    def find(q):
        f = tuple(np.round((qinv @ q) % 1.0, 6) % 1.0)
        if f not in frac:
            sys.exit(f'mlo_spread: k+b = {q} is not on the mesh')
        return frac[f]
    return find


def write_bbvec(qbz, bb, nbasis, nspx, find):
    nbb = len(bb) + 1                                            # the last one is b = 0
    with open('BBVEC', 'w') as f:
        f.write(' nbb,nqbz (mlo_spread.py: the last b is 0, weight 0)\n')
        f.write(f' {nbb} {len(qbz)}\n')
        for b, w in zip(bb, wbb_global):
            f.write(f' {b[0]:.16f} {b[1]:.16f} {b[2]:.16f} {w:.16f}\n')
        f.write(' 0.0 0.0 0.0 0.0\n')
        for iq, k in enumerate(qbz):
            f.write(f' {iq + 1} {k[0]:.16f} {k[1]:.16f} {k[2]:.16f}\n')
            for ib, b in enumerate(list(bb) + [np.zeros(3)]):
                f.write(f' {iq + 1} {ib + 1} {find(k + b) + 1} {k[0] + b[0]:.16f} {k[1] + b[1]:.16f} {k[2] + b[2]:.16f}\n')
        f.write(' nspx\n')
        f.write(f' {nspx}\n')
        for _ in range(nspx):
            f.write(f' 1 {nbasis} 0\n')


def read_uu(nqbz, nbb, nbasis, head):
    """M[k, b, i, j] = <u_ik|u_j,k+b> from the files of huumat --job=2."""
    M = np.zeros((nqbz, nbb, nbasis, nbasis), dtype=complex)
    pair = {}
    for iq in range(nqbz):
        f = FortranFile(f'{head}{iq + 1:04d}', 'r')
        for _ in range(nbb):
            tag = f.read_ints(np.int32)[0]
            if tag == -10:
                iqx, ibx, ikb = f.read_ints(np.int32)
                M[iqx - 1, ibx - 1] = f.read_record(np.complex128).reshape((nbasis, nbasis), order='F')
            elif tag == -20:
                iqx, ibx, iqt, ibt = f.read_ints(np.int32)
                pair[(iqx - 1, ibx - 1)] = (iqt - 1, ibt - 1)
            else:
                sys.exit(f'mlo_spread: unexpected record {tag} in {head}{iq + 1:04d}')
        f.close()
    for (iq, ib), (iqt, ibt) in pair.items():                   # M(k,b) = M(k+b,-b)^+
        M[iq, ib] = M[iqt, ibt].conj().T
    return M


def main():
    global wbb_global
    p = argparse.ArgumentParser(prog='mlo_spread.py', description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('sname')
    p.add_argument('-np', type=int, default=4, help='MPI ranks for huumat')
    p.add_argument('--norun', action='store_true', help='do not run huumat (read the UUU.* there are)')
    p.add_argument('--dwnb', default='mlo', help=argparse.SUPPRESS)   # wan: the Wannier functions, for the comparison of 2026-10-02
    p.add_argument('--nbasis', type=int, default=0, help=argparse.SUPPRESS)
    args = p.parse_args()
    s = args.sname
    gw = tomllib.load(open(f'ctrlg.{s}.toml', 'rb')).get('gw', {})
    n123 = gw.get('n1n2n3')
    plat, qlat, alat = read_lattice()
    qbz = read_qbz()
    if n123 is None:
        sys.exit('mlo_spread: [gw] n1n2n3 not in ctrlg')
    if np.prod(n123) != len(qbz):
        sys.exit(f'mlo_spread: n1n2n3 {n123} does not match QBZ.chk ({len(qbz)} k)')
    if args.dwnb == 'mlo':
        f = FortranFile('HamRsMLO', 'r'); nbasis, _, nspx = f.read_ints(np.int32); f.close()
    else:
        nbasis, nspx = args.nbasis, tomllib.load(open(f'ctrlg.{s}.toml', 'rb')).get('ham', {}).get('nspin', 1)
    bb, wbb_global = bshells(qlat, n123)
    find = kindex(qbz, qlat)
    write_bbvec(qbz, bb, nbasis, nspx, find)
    print(f'mlo_spread: {len(bb)} b vectors (+ b=0), |b| = {sorted(set(np.round(np.linalg.norm(bb, axis=1), 6)))} 2pi/alat, '
          f'{nbasis} MLOs, nspx={nspx}, mesh {n123}')
    nbb = len(bb) + 1
    if not args.norun:
        exe = EXEC_DIR / 'huumat_MPI'                             # next to this script in the bin directory
        if not exe.exists():
            exe = shutil.which('huumat_MPI') or sys.exit('mlo_spread: huumat_MPI not found')
        cmd = ['mpirun', '-np', str(args.np), str(exe), s, f'--dwnb={args.dwnb}', '--job=2']
        with open('lhuumat_spread', 'w') as lf:
            r = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)
        if r.returncode != 0:
            sys.exit('mlo_spread: huumat failed, see lhuumat_spread')
    u = alat / (2 * np.pi)                                       # 2pi/alat -> bohr^-1
    for isp, head in enumerate(['UUU.', 'UUD.'][:nspx]):
        M = read_uu(len(qbz), nbb, nbasis, head)
        O = np.real(np.einsum('kii->ki', M[:, -1]))              # b = 0
        Mb = np.einsum('kbii->kbi', M[:, :-1])
        r = -np.einsum('b,ba,kbi->ia', wbb_global, bb, Mb.imag) / len(qbz) * u
        r2 = 2 * np.einsum('b,kbi->i', wbb_global, O[:, None, :] - Mb.real) / len(qbz) * u * u
        om = r2 - np.sum(r * r, axis=1)
        print(f'isp={isp + 1}  MLO  <O(k)>     <r> (bohr)                      <r^2> (bohr^2)  Omega (bohr^2)  sqrt(Omega) (bohr)')
        for i in range(nbasis):
            print(f'      {i + 1:4d} {O[:, i].mean():8.4f}  {r[i, 0]:9.4f} {r[i, 1]:9.4f} {r[i, 2]:9.4f}  {r2[i]:12.4f} {om[i]:14.4f} {np.sqrt(max(om[i], 0)):12.4f}')
        print(f'      sum Omega = {om.sum():.4f} bohr^2')


if __name__ == '__main__':
    main()
