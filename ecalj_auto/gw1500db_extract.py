#!/usr/bin/env python3
"""Extract the band lines and the total DOS of GW1500 runs into one small npz per (material, source).  2026-10-02.

    gw1500db_extract.py <out dir> < list        list lines: <mpid> <tag> <PlotBand dir>
A tag ending in 'mlo' (2026-10-04) reads the MLO bands of gw1500_mlo.sh in that directory (band_MLO_spin{1,2}.dat, E_F the
ef= of bandplot_MLO.isp1.glt, the same reference as the bnd files) with the check bandcheck.json and mlo_version.txt.
Writes <out dir>/<mpid>.<tag>.npz with
    x (nx,), E (nb, nx) eV from E_F (bands with any point within [-25, 25] eV), seg (nx,) segment index,
    lab, labx (k-point labels and positions), ef_ry, dosE (eV from E_F), dos (states/eV, both spins), nspin,
    and the gap along the path: vbm, cbm (eV), gap_path, gap_direct_path, xv, xc (where VBM and CBM sit), metal (bool).
The bnd files of job_band: one file per segment of the path; columns: band index, x, E - E_F (eV), dE/dx, q.
The header '#   97   <E_F in Ry>  QPE(ev)'.  dos.tot.<mpid>.chk: E (Ry, absolute), DOS (per Ry), per spin.
"""
import glob, json, os, re, sys
import numpy as np

RY = 13.605693122994


def read_bnd(d, spin):
    fs = sorted(glob.glob(f'{d}/bnd[0-9][0-9][0-9].spin{spin}'))
    if not fs:
        return None
    xs, Es, segs, ef = [], [], [], None
    for iseg, f in enumerate(fs):
        lines = open(f).read().split('\n')
        if ef is None:
            ef = float(lines[0].split()[2])      # '#', nband, E_F (Ry)
        a = np.array([[float(v) for v in l.split()[:3]] for l in lines if l.strip() and not l.startswith('#')])
        ib = a[:, 0].astype(int)
        nb = ib.max()
        nk = int((ib == ib[0]).sum())
        E = a[:, 2].reshape(nb, nk) if len(a) == nb * nk else None
        if E is None:
            return None
        xs.append(a[:nk, 1]); Es.append(E); segs.append(np.full(nk, iseg))
    nbmin = min(e.shape[0] for e in Es)
    return np.concatenate(xs), np.concatenate([e[:nbmin] for e in Es], axis=1), np.concatenate(segs), ef


def read_labels(d):
    g = os.path.join(d, 'bandplot.isp1.glt')
    if not os.path.exists(g):
        return [], []
    s = open(g).read()
    m = re.search(r'set xtics \((.*?)\)\s*\n', s, re.S)
    if not m:
        return [], []
    lab, labx = [], []
    for name, x in re.findall(r"'([^']*)'\s+([-0-9.Ee+]+)", m.group(1)):
        lab.append(name.replace('GAMMA', 'Γ')); labx.append(float(x))
    return lab, labx


def read_dos(d, mpid, ef_ry):
    for f in (f'{d}/dos.tot.{mpid}.chk', f'{d}/../dos.tot.{mpid}.chk'):
        if os.path.exists(f):
            a = np.loadtxt(f)
            if a.ndim != 2:
                return None, None, 0
            nsp = a.shape[1] - 1
            return ((a[:, 0] - ef_ry) * RY).astype(np.float32), (a[:, 1:].sum(axis=1) / RY).astype(np.float32), nsp
    return None, None, 0


def gaps(x, E):
    """The band gap along the path: the band count n where max(E[:n]) < min(E[n:]) with E_F (0) inside or within 0.05 eV
    of that window (the VBM on the path can sit a few meV above the E_F fixed on the k mesh). No such n: a metal."""
    E = np.sort(E, axis=0)
    top = np.maximum.accumulate(E.max(axis=1)); bot = np.minimum.accumulate(E.min(axis=1)[::-1])[::-1]
    for n in range(1, len(E)):
        if top[n - 1] < bot[n] and top[n - 1] <= 0.05 and bot[n] >= -0.05:
            v = E[:n].max(axis=0); c = E[n:].min(axis=0)
            iv, ic = int(v.argmax()), int(c.argmin())
            return dict(metal=False, vbm=float(v[iv]), cbm=float(c[ic]), gap_path=float(c[ic] - v[iv]),
                        gap_direct_path=float((c - v).min()), xv=float(x[iv]), xc=float(x[ic]), nocc=n)
    return dict(metal=True, vbm=np.nan, cbm=np.nan, gap_path=0.0, gap_direct_path=0.0, xv=np.nan, xc=np.nan, nocc=0)


def read_mlo(d):
    """MLO bands of job_mlo: rows x, E (Ry, raw), spin, band index; one block of N_MLO rows per k point.
    A segment boundary repeats its x."""
    g = open(os.path.join(d, 'bandplot_MLO.isp1.glt')).read()
    ef = float(re.search(r'ef=\s*([-+0-9.eEdD]+)', g).group(1).replace('D', 'E'))
    Es = []
    for spin in (1, 2):
        f = os.path.join(d, f'band_MLO_spin{spin}.dat')
        if not os.path.exists(f):
            continue
        a = np.loadtxt(f, ndmin=2)
        nb = int(a[:, 3].max()); nk = len(a) // nb
        if nk * nb != len(a) or not (a[:nb, 3] == np.arange(1, nb + 1)).all():
            return None
        E = ((a[:, 1] - ef) * RY).reshape(nk, nb).T
        x = a[::nb, 0]
        Es.append(E)
    if not Es:
        return None
    nk = min(e.shape[1] for e in Es)
    E = np.concatenate([e[:, :nk] for e in Es], axis=0); x = x[:nk]
    seg = np.concatenate([[0], np.cumsum(np.diff(x) <= 1e-9)])
    return x, E, seg, len(Es)


def mlo_npz(out, mpid, tag, d):
    r = read_mlo(d)
    if r is None:
        print(f'{mpid} {tag} NOBAND'); return
    x, E, seg, nspin = r
    c = {}
    jf = os.path.join(d, 'bandcheck.json')
    if os.path.exists(jf):
        c = next(iter(json.load(open(jf)).values()), {})
    vf = os.path.join(d, 'mlo_version.txt')
    ver = open(vf).read().split()[1] if os.path.exists(vf) else ''
    keep = (np.abs(E) < 25).any(axis=1)
    kw = {k: float(c[k]) for k in ('gap_mesh', 'gapD', 'gapM', 'dVBM', 'dCBM', 'rms_m2d', 'max_m2d', 'rms_d2m', 'max_d2m', 'max',
                                    'jump', 'spike', 'emin', 'emax', 'mlo_delta', 'ovlp_min') if c.get(k) is not None}
    np.savez_compressed(f'{out}/{mpid}.{tag}.npz', x=x.astype(np.float32), E=E[keep].astype(np.float32), seg=seg.astype(np.int16),
                        nmlo=int(E.shape[0] // nspin), nspin=nspin, check=c.get('check', ''), fail=' / '.join(c.get('fail', [])),
                        version=ver, **kw)
    print(f'{mpid} {tag} OK nmlo={E.shape[0] // nspin} check={c.get("check", "")}')


def main():
    out = sys.argv[1]; os.makedirs(out, exist_ok=True)
    for line in sys.stdin:
        p = line.split()
        if len(p) < 3:
            continue
        mpid, tag, d = p[:3]
        try:
            if tag.endswith('mlo'):
                mlo_npz(out, mpid, tag, d); continue
            r1 = read_bnd(d, 1)
            if r1 is None:
                print(f'{mpid} {tag} NOBAND'); continue
            x, E, seg, ef = r1
            r2 = read_bnd(d, 2)
            nspin = 2 if r2 is not None else 1
            if r2 is not None:
                E = np.concatenate([E, r2[1][:, :E.shape[1]]], axis=0)
            keep = (np.abs(E) < 25).any(axis=1)
            lab, labx = read_labels(d)
            dE, dos, _ = read_dos(d, mpid, ef)
            g = gaps(x, E)
            np.savez_compressed(f'{out}/{mpid}.{tag}.npz', x=x.astype(np.float32), E=E[keep].astype(np.float32),
                                seg=seg.astype(np.int16), lab=np.array(lab), labx=np.array(labx, np.float32), ef_ry=ef,
                                dosE=dE if dE is not None else np.zeros(0, np.float32),
                                dos=dos if dos is not None else np.zeros(0, np.float32), nspin=nspin, **g)
            print(f'{mpid} {tag} OK gap_path={g["gap_path"]:.3f} direct={g["gap_direct_path"]:.3f} metal={g["metal"]}')
        except Exception as e:
            print(f'{mpid} {tag} ERROR {type(e).__name__}: {e}')


if __name__ == '__main__':
    main()
