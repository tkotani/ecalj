#!/usr/bin/env python3
"""symfind.py: the space group of ctrlg.<sname>.toml by spglib, written to symmetry.json (step S0 of MD/symmetry_spglib.md).

    symfind.py <sname> [--symprec 1e-5] [--out symmetry.json] [--check [llmchk]]

symmetry.json holds the operations in the standard notation of spglib: x' = R x + t with R an integer matrix and t a translation,
both in the fractional coordinates of plat (the columns of plat are the lattice vectors), together with the space group number,
its international symbol and Hall number, and the structure the operations belong to (a normalized string and its SHA-256), so
that a reader can tell whether symmetry.json still matches the ctrlg. Units are named in the file; alat is not used.
The atoms are told apart by the species name of [[site]] atom (so spin-split species such as Niup/Nidn are different).

--check compares with the operations ecalj found (lmchk; the first 'gensym: ig group ops' table of its output, i.e. the
group without the AF operations): the number of operations and whether each ecalj operation is one of spglib's (mod lattice).
Without a file name it runs lmchk itself. Fortran is not changed in S0; ecalj does not read symmetry.json yet. (2026-10-02)
"""
import argparse, hashlib, json, os, re, subprocess, sys, tomllib, datetime
import numpy as np

BOHR_ANGSTROM = 0.529177210903


def read_ctrlg(sname):
    c = tomllib.load(open(f'ctrlg.{sname}.toml', 'rb'))
    st = c['struc']
    alat = float(st['alat'])                                   # bohr
    plat = np.array(st['plat'], float).T                        # columns = lattice vectors, alat units
    names, frac = [], []
    pinv = np.linalg.inv(plat)
    for s in c['site']:
        names.append(s['atom'])
        if 'xpos' in s:
            frac.append(np.array(s['xpos'], float))
        else:                                                   # pos: Cartesian, alat units
            frac.append(pinv @ np.array(s['pos'], float))
    return alat, plat, names, np.array(frac)


def structure_key(alat, plat, names, frac, symprec):
    """The normalized string the hash is taken of: alat (bohr), plat (alat), the sites (species, fractional mod 1), symprec."""
    f = np.mod(frac, 1.0)
    f[np.abs(f - 1.0) < 1e-10] = 0.0
    lines = [f'alat {alat:.10f}', 'plat ' + ' '.join(f'{x:.10f}' for x in plat.T.ravel())]
    lines += [f'site {n} ' + ' '.join(f'{x:.10f}' for x in p) for n, p in zip(names, f)]
    lines += [f'symprec {symprec:.3e}']
    return '\n'.join(lines)


def spglib_dataset(alat, plat, names, frac, symprec):
    import spglib
    species = sorted(set(names), key=names.index)
    numbers = [species.index(n) + 1 for n in names]
    lattice = (alat * BOHR_ANGSTROM * plat).T                   # rows = lattice vectors in angstrom
    ds = spglib.get_symmetry_dataset((lattice, frac, numbers), symprec=symprec)
    if ds is None:
        sys.exit('symfind: spglib found no symmetry (check the structure or symprec)')
    get = (lambda k: getattr(ds, k)) if not isinstance(ds, dict) else (lambda k: ds[k])
    return spglib, get


def ecalj_ops(llmchk):
    """(R_cart, t_cart) of the first 'gensym: ig group ops' table of a lmchk output (symbols as asymop prints them)."""
    lines = open(llmchk, errors='replace').read().split('\n')
    i0 = next((i for i, l in enumerate(lines) if 'gensym: ig group ops' in l), None)
    if i0 is None:
        sys.exit(f'symfind: no gensym table in {llmchk}')
    ops = []
    for l in lines[i0 + 1:]:
        m = re.match(r'\s*(\d+)\s+(\S+)\s*$', l)
        if not m:
            break
        ops.append(parse_op(m.group(2)))
    return ops


def axis(s):
    return {'x': [1, 0, 0], 'y': [0, 1, 0], 'z': [0, 0, 1], 'd': [1, 1, 1]}.get(s) or [float(x) for x in s.strip('()').split(',')]


def rot(n, v):
    """parsop of m_mksym_util: a right-handed rotation by 2 pi / n about v."""
    v = np.array(v, float); v /= np.linalg.norm(v)
    c, s = np.cos(2 * np.pi / n), np.sin(2 * np.pi / n)
    k = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return c * np.eye(3) + (1 - c) * np.outer(v, v) + s * k


def parse_op(sym):
    """e, i, rN<axis>, m<axis>, products with '*', and ':(tx,ty,tz)' (Cartesian, alat) for the translation."""
    t = np.zeros(3)
    if ':' in sym:
        sym, tr = sym.split(':', 1)
        t = np.array([float(x) for x in tr.strip('()').split(',')])
    g = np.eye(3)
    for p in sym.split('*'):
        if p in ('e', 'E'):
            a = np.eye(3)
        elif p in ('i', 'I'):
            a = -np.eye(3)
        elif p[0] in 'rR':
            a = rot(int(p[1]), axis(p[2:]))
        elif p[0] in 'mM':
            v = np.array(axis(p[1:]), float)
            a = np.eye(3) - 2 * np.outer(v, v) / v.dot(v)
        else:
            sys.exit(f'symfind: cannot read the operation {sym}')
        g = g @ a
    return g, t


def main():
    p = argparse.ArgumentParser(prog='symfind.py', description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('sname')
    p.add_argument('--symprec', type=float, default=1e-5, help='spglib tolerance in angstrom (default 1e-5)')
    p.add_argument('--out', default='symmetry.json')
    p.add_argument('--check', nargs='?', const='', default=None, help='compare with the lmchk output (a file, or run lmchk)')
    a = p.parse_args()
    alat, plat, names, frac = read_ctrlg(a.sname)
    spg, get = spglib_dataset(alat, plat, names, frac, a.symprec)
    R, T = np.array(get('rotations')), np.array(get('translations'))
    key = structure_key(alat, plat, names, frac, a.symprec)
    npure = int(sum(1 for r in R if (r == np.eye(3, dtype=int)).all()))
    out = {
        'format': 'ecalj-symmetry-1',
        'made_by': f'symfind.py (spglib {spg.__version__})', 'date': datetime.datetime.now().isoformat(timespec='minutes'),
        'ctrlg': f'ctrlg.{a.sname}.toml',
        'structure_key': key, 'structure_sha256': hashlib.sha256(key.encode()).hexdigest(),
        'symprec': a.symprec, 'symprec_unit': 'angstrom',
        'spacegroup': {'number': int(get('number')), 'international': str(get('international')), 'hall_number': int(get('hall_number'))},
        'n_operations': int(len(R)), 'n_pure_translations': npure,
        'rotation_unit': "integer matrix in the basis of plat: x' = R x + t, x fractional (columns of plat)",
        'translation_unit': 'fractional coordinates of plat',
        'operations': [{'rotation': r.tolist(), 'translation': [round(float(x), 10) for x in t], 'time_reversal': False} for r, t in zip(R, T)],
        'equivalent_atoms': [int(x) for x in get('equivalent_atoms')],
    }
    json.dump(out, open(a.out, 'w'), indent=1)
    print(f'symfind: {a.sname}: {out["spacegroup"]["international"]} ({out["spacegroup"]["number"]}), {len(R)} operations '
          f'({npure} pure translations incl. the identity), symprec {a.symprec} angstrom -> {a.out}')
    if a.check is not None:
        f = a.check
        if f == '':
            f = 'llmchk_symfind'
            with open(f, 'w') as fo:
                subprocess.run(['mpirun', '-np', '1', 'lmchk', a.sname], stdout=fo, stderr=subprocess.STDOUT)
        E = ecalj_ops(f)
        pinv = np.linalg.inv(plat)
        miss = 0
        for g, t in E:
            W = pinv @ g @ plat
            w = pinv @ t
            Wi = np.rint(W)
            if np.abs(W - Wi).max() > 1e-6:
                miss += 1; continue
            ok = any((Wi == r).all() and np.abs((w - tt) - np.rint(w - tt)).max() < 1e-4 for r, tt in zip(R, T))
            miss += 0 if ok else 1
        rots_e = {tuple(np.rint(pinv @ g @ plat).astype(int).ravel()) for g, _ in E}
        rots_s = {tuple(r.ravel()) for r in R}
        status = ('SAME' if (miss == 0 and len(E) == len(R)) else
                  'SUBGROUP' if miss == 0 else 'DIFFERENT')
        print(f'symcheck: {a.sname}: ecalj {len(E)} ops ({len(rots_e)} rotations), spglib {len(R)} ops ({len(rots_s)} rotations), '
              f'ecalj ops not in spglib: {miss} -> {status}')


if __name__ == '__main__':
    main()
