#!/usr/bin/env python3
"""symfind.py: the space group of ctrlg.<sname>.toml by spglib (Python), written to symmetry.<sname>.json.

lmf and lmchk make the same file themselves with the vendored spglib (C) when it is missing or its structure differs
(2026-10-02 08:33, MD/symmetry_spglib.md 4.7g); this tool makes it outside them.

    symfind.py <sname> [--symprec 1e-5] [--out symmetry.<sname>.json] [--ctrlg:<path>=<value> ...]

symmetry.<sname>.json holds the operations in the standard notation of spglib: x' = R x + t with R an integer matrix and t a translation,
both in the fractional coordinates of plat (the columns of plat are the lattice vectors), together with the space group number,
its international symbol and Hall number, and the structure the operations belong to (a normalized string and its SHA-256), so
that a reader can tell whether the file still matches the ctrlg. The structure is also written as numbers ('structure'):
Fortran (m_symfind) compares them with the ctrlg it reads and stops when they differ (2026-10-02 05:33: comparing the
normalized string would need the same formatting and rounding in Fortran). Units are named in the file.
The file name carries <sname>: one directory can hold two ctrlg (Samples/LDAU/ReN).
The atoms are told apart by the species name of [[site]] atom (so spin-split species such as Niup/Nidn are different).

--check (comparison with the old search of ecalj) was removed on 2026-10-02 22:43 with that search.
The file is written to a temporary name and renamed, so that ranks starting at once (a test wrapper) do not read half a file.
"""
import argparse, hashlib, json, os, sys, tomllib, datetime
import numpy as np

BOHR_ANGSTROM = 0.529177210903


def apply_overrides(c, overrides):
    """--ctrlg:<dotted.path>=<value> as lmf takes them (m_toml_override.f90): a TOML value; an integer in the path is a
    1-based index into an array of tables (site.2.pos). So that symmetry.<sname>.json is made for the structure lmf will
    read when it is run with such overrides (2026-10-02 06:05: Samples/AtomDimer sets the positions this way)."""
    for o in overrides:
        path, val = o[len('--ctrlg:'):].split('=', 1)
        v = tomllib.loads(f'v = {val}')['v']
        keys = path.split('.')
        d = c
        for k in keys[:-1]:
            d = d[int(k) - 1] if k.isdigit() else d.setdefault(k, {})
        k = keys[-1]
        if k.isdigit():
            d[int(k) - 1] = v
        else:
            d[k] = v
    return c


def read_ctrlg(sname, overrides=()):
    c = apply_overrides(tomllib.load(open(f'ctrlg.{sname}.toml', 'rb')), overrides)
    st = c['struc']
    alat = float(st['alat'])                                   # bohr
    plat = np.array(st['plat'], float).T                        # columns = lattice vectors, alat units
    names, frac, af = [], [], []
    pinv = np.linalg.inv(plat)
    for s in c['site']:
        names.append(s['atom'])
        af.append(int(s.get('af', 0)))                          # AF pair label: +k and -k are exchanged by the AF operations
        if 'xpos' in s:
            frac.append(np.array(s['xpos'], float))
        else:                                                   # pos: Cartesian, alat units
            frac.append(pinv @ np.array(s['pos'], float))
    return alat, plat, names, np.array(frac), af


def structure_key(alat, plat, names, frac, symprec):
    """The normalized string the hash is taken of: alat (bohr), plat (alat), the sites (species, fractional mod 1), symprec."""
    f = np.mod(frac, 1.0)
    f[np.abs(f - 1.0) < 1e-10] = 0.0
    lines = [f'alat {alat:.10f}', 'plat ' + ' '.join(f'{x:.10f}' for x in plat.T.ravel())]
    lines += [f'site {n} ' + ' '.join(f'{x:.10f}' for x in p) for n, p in zip(names, f)]
    lines += [f'symprec {symprec:.3e}']
    return '\n'.join(lines)


def magnetic_operations(alat, plat, names, frac, af, symprec):
    """AF (step S5, 2026-10-02 06:14): the sites +k and -k of an AF pair are one type with the magnetic moments +1 and -1
    (other sites 0), as ecalj merges their species for the lattice+AF group (ipsAF of m_mksym). Returns (R, T, time_reversal)
    of the magnetic space group: the operations without time reversal keep the spins (the group of the crystal), those with
    it exchange up and down (the AF operations that SYMGRPAF generated)."""
    import spglib
    rep = {}
    for i, a in enumerate(af):
        if a > 0:
            rep[a] = names[i]
    tnames = []
    for i, a in enumerate(af):
        if a < 0:
            if -a not in rep:
                sys.exit(f'symfind: site {i+1} has af={a} but no site has af={-a}')
            tnames.append(rep[-a])
        else:
            tnames.append(names[i])
    species = sorted(set(tnames), key=tnames.index)
    numbers = [species.index(n) + 1 for n in tnames]
    magmoms = [float(np.sign(a)) for a in af]
    lattice = (alat * BOHR_ANGSTROM * plat).T
    ds = spglib.get_magnetic_symmetry_dataset((lattice, frac, numbers, magmoms), symprec=symprec)
    if ds is None:
        sys.exit('symfind: spglib found no magnetic symmetry')
    get = (lambda k: getattr(ds, k)) if not isinstance(ds, dict) else (lambda k: ds[k])
    return np.array(get('rotations')), np.array(get('translations')), np.array(get('time_reversals'), bool), get


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


def main():
    p = argparse.ArgumentParser(prog='symfind.py', description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('sname')
    p.add_argument('--symprec', type=float, default=1e-5, help='spglib tolerance in angstrom (default 1e-5)')
    p.add_argument('--out', default=None, help='default symmetry.<sname>.json')
    a, rest = p.parse_known_args()       # other arguments of lmf (--quit=band, ...) are ignored; --ctrlg: overrides are applied
    overrides = [x for x in rest if x.startswith('--ctrlg:')]
    if a.out is None:
        a.out = f'symmetry.{a.sname}.json'
    alat, plat, names, frac, af = read_ctrlg(a.sname, overrides)
    spg, get = spglib_dataset(alat, plat, names, frac, a.symprec)
    R, T = np.array(get('rotations')), np.array(get('translations'))
    TR = np.zeros(len(R), bool)
    magnetic = None
    if any(af):
        Rm, Tm, TRm, gm = magnetic_operations(alat, plat, names, frac, af, a.symprec)
        # the operations without time reversal must be those of the crystal with up and down told apart
        same = sorted(map(lambda x: (tuple(x[0].ravel()), tuple(np.round(x[1] % 1.0, 6) % 1.0)), zip(Rm[~TRm], Tm[~TRm]))) == \
               sorted(map(lambda x: (tuple(x[0].ravel()), tuple(np.round(x[1] % 1.0, 6) % 1.0)), zip(R, T)))
        if not same:
            sys.exit('symfind: the magnetic group without time reversal differs from the group of the crystal (check af= and the species)')
        R, T, TR = Rm, Tm, TRm
        magnetic = {'uni_number': int(gm('uni_number')), 'msg_type': int(gm('msg_type')),
                    'n_time_reversed': int(TRm.sum())}
    T = np.where(np.abs(T - np.rint(T)) < 1e-10, np.rint(T), T) % 1.0      # translations in [0,1)
    i0 = next(i for i, (r, tt, tr) in enumerate(zip(R, T, TR))
              if (r == np.eye(3, dtype=int)).all() and np.abs(tt).max() < 1e-10 and not tr)
    order = [i0] + [i for i in range(len(R)) if i != i0 and not TR[i]] + [i for i in range(len(R)) if TR[i]]
    R, T, TR = R[order], T[order], TR[order]          # the identity first (ecalj assumes it), then the rest of the crystal group, then AF
    key = structure_key(alat, plat, names, frac, a.symprec)
    npure = int(sum(1 for r, tr in zip(R, TR) if (r == np.eye(3, dtype=int)).all() and not tr))
    out = {
        'format': 'ecalj-symmetry-1',
        'made_by': f'symfind.py (spglib {spg.__version__})', 'date': datetime.datetime.now().isoformat(timespec='minutes'),
        'ctrlg': f'ctrlg.{a.sname}.toml',
        'structure_key': key, 'structure_sha256': hashlib.sha256(key.encode()).hexdigest(),
        'structure': {'alat_bohr': alat, 'plat_alat': plat.T.tolist(), 'plat_unit': 'rows = lattice vectors in alat (as [struc] plat)',
                      'species': names, 'af': af, 'frac': [[float(x) for x in p] for p in frac],
                      'frac_unit': 'fractional coordinates of plat (not reduced mod 1)'},
        'symprec': a.symprec, 'symprec_unit': 'angstrom',
        'spacegroup': {'number': int(get('number')), 'international': str(get('international')), 'hall_number': int(get('hall_number'))},
        'n_operations': int(len(R)), 'n_pure_translations': npure,
        'rotation_unit': "integer matrix in the basis of plat: x' = R x + t, x fractional (columns of plat)",
        'translation_unit': 'fractional coordinates of plat',
        'magnetic': magnetic,
        'operations': [{'rotation': r.tolist(), 'translation': [round(float(x), 10) for x in t], 'time_reversal': bool(tr)}
                       for r, t, tr in zip(R, T, TR)],
        'equivalent_atoms': [int(x) for x in get('equivalent_atoms')],
    }
    tmp = f'{a.out}.{os.getpid()}.tmp'
    json.dump(out, open(tmp, 'w'), indent=1)
    os.replace(tmp, a.out)
    print(f'symfind: {a.sname}: {out["spacegroup"]["international"]} ({out["spacegroup"]["number"]}), {len(R)} operations '
          f'({npure} pure translations incl. the identity), symprec {a.symprec} angstrom -> {a.out}')

if __name__ == '__main__':
    main()
