#!/usr/bin/env python3
"""Put empty spheres (ES) at the large voids of a structure, for the MLO model (and for the basis).

    ctrlg_addes.py <sname> [--rmin 3.0] [--scale 0.9] [--rmax 4.0] [--dry-run]

Reads ctrlg.<sname>.toml. The "void radius" of a point is its distance to the nearest MT sphere surface
(|r - R_a| - r_a, over all atoms a). Every local maximum of it above --rmin (a.u.) gets an ES:
  - its position is snapped to a multiple of 1/48 in the plat basis when that is within 0.02, to keep the symmetry;
  - whole Wyckoff orbits (2026-10-05, wyckoff_voids: spglib; user: "ES at the Wyckoff positions"; an ES on part of an orbit
    lowers the symmetry, which makes GW heavier and spoils the band plot); the space-group operations with and without ES are printed;
  - its radius is --scale x its own void radius, at most --rmax (2026-10-05: 0.9 and 4.0 a.u.: an ES is expanded up to l = 2,
    too little for a sphere of 6-15 a.u. in a molecular crystal; the APWs of PMT take the rest of a wide void; before, 0.6 x the smallest void radius for all ES, which
    left ES of 1.8 a.u. in voids of 3-6 a.u.); the largest voids are taken first, a smaller one only when its ES does
    not overlap the kept ones;
  - the ES are added as the LAST sites (the indices of the atoms, and the mlo_lm rows, do not change), with one
    [[spec]] per radius, "E" (or E1, E2, ... for several radii; z = 0, lmxa = 2, rsmh = r/2, eh = -0.3; the basis s only
    (lmx = 0) below --rsp = 3.0 a.u., s,p (lmx = 1) from there; 2026-10-05, user: "keep the Wyckoff positions, fill by number
    when they get small, s only when small" -- the number of orbitals sets the time; s,p,d for every ES before), and a row
    "<i> E 1" or "<i> E 1 2 3 4" in [mlo] mlo_lm.
The basis changes, so the calculation starts again from lmfa. For GW, run gwinit again (the per-atom tables of
[product_basis] do not list the ES).

Why 3.0 a.u. (2026-10-01, Samples/MATERIALS): the void radius is 4.30 a.u. for SiO2 cristobalite; zinc blende and
wurtzite compounds have 2.0-2.9 (EH2 s,p in mlo_lm2 is enough for them), rock salt and perovskite 1.2-1.4, metals below
1.1. The van der Waals gap of Bi2Te3 is 2.61 a.u. at its centre (1/2,1/2,1/2), 2.68 at the best point (2026-10-02 06:26, with
enough images; 3.56 written here before was wrong), below 3.0: no ES there. SiO2 with ES at its two voids: MLO bands within
0.001 eV of DFT (the conduction band was missing from the model without them).
"""
import argparse, itertools, re, sys
import numpy as np
try:
    import tomllib
except ImportError:
    sys.exit('ctrlg_addes.py: python 3.11 or later is needed (tomllib)')


def structure(d):
    alat = d['struc']['alat']; plat = np.array(d['struc']['plat'], float)
    r = {s['atom']: s['r'] for s in d['spec']}
    pos = []
    for s in d['site']:
        p = np.array(s['pos'], float) if 'pos' in s else np.array(s['xpos'], float) @ plat
        pos.append((p * alat if 'pos' in s else p * alat, r[s['atom']]))
    return alat, plat * alat, pos


def clearance_fn(plat, pos, dcut=8.0):
    # The images: every lattice translation that can bring an atom within dcut (a.u., > void radius + MT radius) of a point
    # of the cell, from the spacing of the lattice planes. (2026-10-02 06:26: the images -1..1 along each vector missed near
    # atoms in the long rhombohedral cell of Bi2Te3, |a_i| = 19.8 a.u.; the void came out larger than it is, and the local
    # refinement below could walk away from the atoms without end, where the clearance grows without bound.)
    q = np.linalg.inv(plat).T                                    # rows: 1/|q_i| is the spacing of the lattice planes i
    nimg = [int(np.ceil(dcut * np.linalg.norm(qi))) + 1 for qi in q]
    shifts = [np.array(t) @ plat for t in itertools.product(*[range(-k, k + 1) for k in nimg])]
    A = np.array([p + sh for p, _ in pos for sh in shifts]); R = np.array([rr for _, rr in pos for sh in shifts])
    def c(x):
        x = np.atleast_2d(x)
        return (np.linalg.norm(x[:, None, :] - A[None, :, :], axis=2) - R[None, :]).min(axis=1)
    return c


def grid(plat):
    # about 0.6 a.u. between the points along each axis (a long axis, e.g. c of Bi2Te3, needs more than a fixed 24;
    # 2026-10-01: 24 missed the van der Waals gap of Bi2Te3)
    return [max(16, int(np.ceil(np.linalg.norm(a) / 0.6))) for a in plat]


def min_image(df, plat):
    """The fractional vector equal to df modulo the lattice that is shortest in Cartesian space. Bug fixed 2026-10-06 12:05:
    the distances were taken after wrapping each fractional component into [-1/2, 1/2), which misses the nearest image in an
    oblique cell (hexagonal, the primitive cell of a centred lattice); the shortest distance in a Wyckoff orbit came out too
    long and the ES overlapped one another (Ba2CuClO2 4.9 %, h-BN 22 %, NiO2 142 %). Those ES runs of GW1500 all stopped:
    a NaN in hvccfp0 (the overlap matrix of the interstitial plane waves is not positive definite) or CholQR in lmf."""
    w = (np.asarray(df, float) + 0.5) % 1.0 - 0.5
    c = w + _SHIFTS
    return c[np.argmin(np.linalg.norm(c @ plat, axis=1))]


_SHIFTS = np.array(list(itertools.product((-2, -1, 0, 1, 2), repeat=3)), float)


def find_voids(plat, pos, rmin):
    c = clearance_fn(plat, pos)
    n = grid(plat)
    frac = np.array(list(itertools.product(range(n[0]), range(n[1]), range(n[2]))), float) / np.array(n)
    val = np.concatenate([c(frac[i:i + 20000] @ plat) for i in range(0, len(frac), 20000)]).reshape(n)
    peaks = []
    for i, j, k in zip(*np.where(val > rmin - 1.0)):   # refined below; a peak may rise above rmin
        v = val[i, j, k]
        nb = [val[(i + a) % n[0], (j + b) % n[1], (k + e) % n[2]] for a, b, e in itertools.product((-1, 0, 1), repeat=3) if (a, b, e) != (0, 0, 0)]
        if v >= max(nb) - 1e-12: peaks.append(np.array([i, j, k], float) / np.array(n))
    out = []
    inv = np.linalg.inv(plat)
    for f in peaks:
        x = f @ plat; best = c(x)[0]
        for step in (0.25, 0.12, 0.06, 0.03, 0.015):   # local refinement (a.u.)
            moved = True
            while moved:
                moved = False
                for t in itertools.product((-1, 0, 1), repeat=3):
                    y = x + np.array(t) * step; v = c(y)[0]
                    if v > best + 1e-9: best, x, moved = v, y, True
        fr = (x @ inv) % 1.0
        snap = np.round(fr * 48) / 48
        if np.all(np.abs(snap - fr) < 0.02): fr = snap % 1.0
        x = fr @ plat; best = c(x)[0]
        if best <= rmin: continue
        if any(np.linalg.norm(min_image(fr - g, plat) @ plat) < 0.5 for g, _ in out): continue
        out.append((fr, best))
    return out


def select_voids(voids, plat, scale, rmax=4.0):
    """Keep the largest voids first; a smaller one only when its ES would not overlap the kept ones, i.e. the centres are at
    least r_i + r_j apart, r = min(scale v, rmax) the ES radius. Bug fixed 2026-10-05: all local maxima above rmin were kept, so a wide
    void (a molecular crystal, a slab-like structure) got up to 24 ES whose common radius, half the distance to the next ES,
    fell to 0.25 a.u. (53 of the 209 GW1500 structures with a void above 3.0 a.u. got ES below 1.5 a.u.)."""
    kept = []
    for fr, v in sorted(voids, key=lambda a: -a[1]):
        if all(np.linalg.norm(min_image(fr - g, plat) @ plat) >= min(scale * v, rmax) + min(scale * w, rmax) for g, w in kept):
            kept.append((fr, v))
    return kept



def wyckoff_voids(voids, plat, d, scale, rmax, clear, rfloor=1.5):
    """ES on whole Wyckoff orbits (2026-10-05, user: "ES at the Wyckoff positions"; a part of an orbit lowers the symmetry,
    and the GW of the ES runs of GW1500 got 1.4-2.7 times the irreducible q points). For each void, from the largest: the images
    by the space group (spglib), the centre moved to the special position (the mean over its site symmetry), the orbit taken
    whole when it overlaps no ES taken before. ES radius: min(scale x void radius, rmax, half the shortest distance in the
    orbit); an orbit whose radius would fall below rfloor (a.u.) is left out. Returns [(fractional position, void radius,
    ES radius)]."""
    import spglib
    lab = [s['atom'] for s in d['site']]; ids = {a: i + 1 for i, a in enumerate(dict.fromkeys(lab))}
    alat = d['struc']['alat']; inv = np.linalg.inv(plat)
    frac = [(np.array(s['pos'], float) * alat) @ inv if 'pos' in s else np.array(s['xpos'], float) for s in d['site']]
    sym = spglib.get_symmetry((plat, frac, [ids[a] for a in lab]), symprec=1e-3)
    ops = list(zip(sym['rotations'], sym['translations']))
    dist = lambda f, g: np.linalg.norm(min_image(f - g, plat) @ plat)
    lmin = min(np.linalg.norm(n @ plat) for n in _SHIFTS if np.any(n))   # the shortest lattice vector
    kept, seen = [], []
    for fr, v in sorted(voids, key=lambda a: -a[1]):
        if any(dist(fr, g) < 0.3 for g in seen):
            continue
        x = fr
        for _ in range(4):   # images closer than an ES diameter are merged into their centre, a position of higher symmetry
            stab = [R @ x + t for R, t in ops if dist(R @ x + t, x) < 2 * rfloor]
            x = (x + np.mean([min_image(y - x, plat) for y in stab], axis=0)) % 1.0     # the special position
            orb = []
            for R, t in ops:
                y = (R @ x + t) % 1.0
                if all(dist(y, g) > 1e-3 for g in orb):
                    orb.append(y)
            dmin = min([dist(a, b) for i, a in enumerate(orb) for b in orb[i + 1:]], default=1e9)
            if dmin >= 2 * rfloor:
                break
        v = float(clear(x @ plat)[0])                                                              # the void radius there
        seen += orb + [fr]
        r = min(scale * v, rmax, 0.5 * min(dmin, lmin) - 0.01)   # lmin: an ES and its own lattice images (2026-10-06)
        if r < rfloor:
            continue
        if any(dist(y, g) < r + rg for y in orb for g, _, rg in kept):
            continue
        kept += [(y, v, r) for y in orb]
    return kept, len(ops)


def nops_with_es(plat, d, es):
    """number of space-group operations of the structure with the ES (each ES radius a species of its own)"""
    import spglib
    lab = [s['atom'] for s in d['site']]; alat = d['struc']['alat']; inv = np.linalg.inv(plat)
    frac = [(np.array(s['pos'], float) * alat) @ inv if 'pos' in s else np.array(s['xpos'], float) for s in d['site']]
    ids = {a: i + 1 for i, a in enumerate(dict.fromkeys(lab))}; n0 = len(ids)
    rs = sorted({round(r, 2) for _, _, r in es})
    num = [ids[a] for a in lab] + [n0 + 1 + rs.index(round(r, 2)) for _, _, r in es]
    return len(spglib.get_symmetry((plat, frac + [f for f, _, _ in es], num), symprec=1e-3)['rotations'])

def main():
    ap = argparse.ArgumentParser(description='put empty spheres at the voids larger than --rmin (a.u.)')
    ap.add_argument('sname'); ap.add_argument('--rmin', type=float, default=3.0); ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--scale', type=float, default=0.9, help='ES radius / void radius')
    ap.add_argument('--rmax', type=float, default=4.0, help='largest ES radius (a.u.)')
    ap.add_argument('--rsp', type=float, default=3.0, help='an ES of this radius (a.u.) or larger gets s,p; a smaller one s only')
    a = ap.parse_args()
    f = f'ctrlg.{a.sname}.toml'; t = open(f).read(); d = tomllib.loads(t)
    if any(s['z'] == 0 for s in d['spec']): sys.exit(f'{f} has an empty sphere already (z = 0); nothing done')
    alat, plat, pos = structure(d)
    es, nops = wyckoff_voids(find_voids(plat, pos, a.rmin), plat, d, a.scale, a.rmax, clearance_fn(plat, pos))
    voids = [(f, v) for f, v, _ in es]
    if not voids:
        print(f'no void larger than {a.rmin} a.u.; nothing done')
        return
    cart = [fr @ plat for fr, _ in voids]
    rad = [np.floor(r / 0.05) * 0.05 for _, _, r in es]          # on a 0.05 a.u. grid: one [[spec]] per radius
    n2 = nops_with_es(plat, d, es)
    print(f'space-group operations: {nops} without ES, {n2} with ES' + ('' if n2 == nops else '  (LOWER: check)'))
    radii = sorted(set(round(float(r), 2) for r in rad), reverse=True)
    name = {r: ('E' if len(radii) == 1 else f'E{k + 1}') for k, r in enumerate(radii)}
    nat = len(d['site'])
    for i, ((fr, v), x, r) in enumerate(zip(voids, cart, rad), nat + 1):
        print(f'ES site {i}: plat fraction {np.round(fr, 4).tolist()}  pos/alat {np.round(x / alat, 5).tolist()}  void radius {v:.2f} a.u.'
              f'  ES {name[round(r, 2)]} radius {r:.2f} a.u.')
    print(f'ES radius {min(radii):.2f} a.u. (smallest; {len(voids)} spheres, radii {radii})')
    if a.dry_run: return
    sites = ''.join(f'[[site]]   # empty sphere at a void, void radius {v:.2f} a.u. (ctrlg_addes.py)\natom = "{name[round(r, 2)]}"\n'
                    f'pos  = [{x[0]/alat:.8f}, {x[1]/alat:.8f}, {x[2]/alat:.8f}]\n'
                    '# ----------------------------------------------------------------\n' for (fr, v), x, r in zip(voids, cart, rad))
    # after the last [[site]] table
    last = [m.start() for m in re.finditer(r'(?m)^\[\[site\]\]', t)][-1]
    m = re.search(r'(?m)^(\[|# ===)', t[last + 8:]); end = last + 8 + (m.start() if m else len(t) - last - 8)
    t = t[:end] + sites + t[end:]
    LSP = lambda r: 1 if r >= a.rsp else 0     # s,p for a large ES, s only for a small one (2026-10-05, user)
    spec = ''.join(f'[[spec]]   # empty sphere (ctrlg_addes.py; {"s,p" if LSP(r) else "s only"}, 2026-10-05)\natom   = "{name[r]}"\nz      = 0\nr      = {r:.2f}\n'
                   f'lmx    = {LSP(r)}\nlmxa   = 2\nrsmh   = [{", ".join([f"{r/2:.3f}"] * (LSP(r) + 1))}]\neh     = [{", ".join(["-0.3"] * (LSP(r) + 1))}]\n'
                   '# ----------------------------------------------------------------\n' for r in radii)
    last = [m.start() for m in re.finditer(r'(?m)^\[\[spec\]\]', t)][-1]
    m = re.search(r'(?m)^(\[|# ===)', t[last + 8:]); end = last + 8 + (m.start() if m else len(t) - last - 8)
    t = t[:end] + spec + t[end:]
    rows = ''.join(f'{i} {name[round(r, 2)]}    {"1 2 3 4" if LSP(round(r, 2)) else "1"}\n' for i, r in zip(range(nat + 1, nat + 1 + len(voids)), rad))
    if re.search(r'(?m)^mlo_lm = """\n', t):
        t = re.sub(r'(?ms)^(mlo_lm = """\n.*?)^"""', lambda mm: mm.group(1) + rows + '"""', t, count=1)
    open(f + '.bak_addes', 'w').write(open(f).read())
    open(f, 'w').write(t)
    tomllib.loads(t)
    print(f'{f} rewritten (the old one is {f}.bak_addes). Start again from lmfa; for GW run gwinit again.')


if __name__ == '__main__':
    main()
