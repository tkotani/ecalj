#!/usr/bin/env python3
"""Build the GW1500 QSGW80 database (markdown + figures + tsv) from the extracted npz and the status tables.  2026-10-02.

    gw1500db_build.py <db dir> [--figs] [--only mp-1,mp-2]

Inputs (paths below): <db dir>/npz/<mpid>.<tag>.npz (gw1500db_extract.py; tags may_qsgw, may_lda, rr_<run>, db_qsgw, db_lda),
<db dir>/logs/*.log (rerun.log of the database run on kt1 and kr7), ecalj_auto tables, MP json, POSCARs, the 2025 table of
tkotani/DOSnpSupplement (dosnp_README_2025.md).
Outputs: <db dir>/README.md (conditions, summary, how to read), table.md and gw1500db.tsv (one row per material),
bands_<n>atoms.md (figures), fig/<mpid>.png.
Condition sets: M (2026-04/05 production), R (reruns 2026-09-30..10-02), N (database run 2026-10-02 night).
The value adopted: N converged > R converged > M GOOD.
"""
import csv, glob, json, os, re, sys
from collections import Counter, defaultdict
import numpy as np

ECALJ = os.path.expanduser('~/ecalj/ecalj_auto')
WORK = os.path.expanduser('~/work/gw1500db')
AGREE, LARGE = 0.05, 0.2      # eV: |difference| <= AGREE is written as one number; > LARGE is marked to be checked

COND = {
    'M': ('2026-04/05', 'production of April-May 2026: ecalj of that time (~/bin2), input ctrl+GWinput (most) or the first TOML; '
          'GPU products all TF32 (--mp); chi0 by the T=0 tetrahedron method; Sigma levels Gaussian-smeared (esmr = 0.003 Ry); '
          'QSGW80 (scaledsigma 0.8), k 8x8x8 for lmf, 4x4x4 for GW; 5 iterations, then gwscconv until the gap of the last 3 '
          'iterations stays within 0.1 eV (at most 10)'),
    'R': ('2026-09-30/10-02', 'reruns from scratch: ecalj 989a18637 (2026-09-29); POSCAR -> vasp2ctrl -> ctrlgenToml.py --ssig=0.8; '
          '--prec=fp32; [gw] t_tetrakbt = +300 K (finite-T tetrahedron), t_sigmaw = 300 K; same k meshes; gwscconv: gap of the '
          'last 3 iterations within 0.1 eV twice (metals: eigenvalues within 5 eV of E_F change < 0.03 eV twice), at most 10'),
    'N': ('2026-10-02 night', 'database run: as R, but --prec=tf32 and t_tetrakbt = -300 K (chi0 by the T=0 tetrahedron method, '
          'Im chi0 smoothed by a Gaussian of the width of a 300 K Fermi-Dirac), t_sigmaw = 300 K; ecalj 989a18637 (kt1) or '
          'b695fa65c (2026-09-30, kr7) with the same gwscconv'),
}


def fnum(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def load():
    S = {r['mpid']: r for r in csv.DictReader(open(f'{ECALJ}/gw1500_status_20260930.tsv'), delimiter='\t')}
    N = {r['mpid']: r for r in csv.DictReader(open(f'{ECALJ}/gw1500_notes_20261001.tsv'), delimiter='\t')}
    MP = json.load(open(f'{ECALJ}/gw1500_mp_20261001.json'))
    old = {}
    for l in open(f'{WORK}/dosnp_README_2025.md'):
        if re.match(r'\|\s*\d+\s*\|', l):
            c = [x.strip() for x in l.strip().strip('|').split('|')]
            old['mp-' + c[0]] = dict(sgnum=c[3], sgsym=c[4], lda=fnum(c[5]), shot1=fnum(c[6]), shot2=fnum(c[7]), di=c[8])
    return S, N, MP, old


def load_db_logs(db):
    """rerun.log lines of the database run: <date> <time> <worker> <mpid> <verdict> iter=.. gapLDA=.. gap=.. <s>s dqp=.."""
    D = {}
    for f in glob.glob(f'{db}/logs/*.log'):
        host = os.path.basename(f).split('.')[0]
        for l in open(f):
            p = l.split()
            if len(p) < 6 or not p[3].startswith('mp-'):
                continue
            kv = dict(x.split('=', 1) for x in p[5:] if '=' in x)
            sec = next((x[:-1] for x in p[5:] if re.fullmatch(r'\d+s', x)), None)
            D[p[3]] = dict(verdict=p[4], iter=kv.get('iter'), gaplda=fnum(kv.get('gapLDA')), gap=fnum(kv.get('gap')),
                           sec=fnum(sec), host=host, when=f'{p[0]} {p[1]}')
    return D


def spacegroups():
    import spglib
    out = {}
    for f in glob.glob(f'{ECALJ}/INPUT/gw1500/POSCARALL/POSCAR.mp-*'):
        m = f.split('POSCAR.')[1]
        try:
            L = open(f).read().split('\n')
            sc = float(L[1]); lat = np.array([[float(x) for x in L[i].split()[:3]] for i in (2, 3, 4)]) * sc
            names = L[5].split(); counts = [int(x) for x in L[6].split()]
            mode = L[7].strip()[0].lower(); n = sum(counts)
            pos = np.array([[float(x) for x in L[8 + i].split()[:3]] for i in range(n)])
            if mode in 'ck':
                pos = pos @ np.linalg.inv(lat)
            typ = sum(([i] * c for i, c in enumerate(counts)), [])
            ds = spglib.get_symmetry_dataset((lat, pos, typ), symprec=1e-3)
            out[m] = (ds.number, ds.international) if ds is not None else ('', '')
        except Exception:
            out[m] = ('', '')
    return out


def npz(db, m, tag):
    f = f'{db}/npz/{m}.{tag}.npz'
    return np.load(f) if os.path.exists(f) else None


def rr_tag(db, m):
    for r in ('run3x', 'run3', 'run2m', 'run1', 'try'):
        if os.path.exists(f'{db}/npz/{m}.rr_{r}.npz'):
            return f'rr_{r}'
    return None


def decide(m, S, N, D):
    """values per condition set and the one adopted"""
    s, n = S[m], N[m]
    v = {}
    if s['final_state'] == 'GOOD' and fnum(s['gap_QSGW80_last_llmf_eV']) is not None:
        v['M'] = dict(gap=fnum(s['gap_QSGW80_last_llmf_eV']), lda=fnum(s['gap_LDA_eV']), state='GOOD', iter=s['qsgw_iter_final'])
    if n['rerun_verdict'] in ('CONVERGED', 'CONVERGED_METAL'):
        g = 0.0 if n['rerun_verdict'] == 'CONVERGED_METAL' else fnum(n['rerun_gap_fp32_eV'])
        v['R'] = dict(gap=g, lda=fnum(n['rerun_gapLDA_eV']), state=n['rerun_verdict'], iter=n['rerun_iter'])
    d = D.get(m)
    if d and d['verdict'] in ('CONVERGED', 'CONVERGED_METAL'):
        g = 0.0 if d['verdict'] == 'CONVERGED_METAL' else d['gap']
        v['N'] = dict(gap=g, lda=d['gaplda'], state=d['verdict'], iter=d['iter'], host=d['host'])
    adopt = next((k for k in ('N', 'R', 'M') if k in v and v[k]['gap'] is not None), None)
    return v, adopt


def fmt(x, nd=2):
    return '' if x is None else f'{x:.{nd}f}'


def figure(db, m, title, curves, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, (ax, ad) = plt.subplots(1, 2, figsize=(7.2, 4.2), gridspec_kw=dict(width_ratios=[4, 1], wspace=0.04), sharey=True)
    ref = None
    for (d, color, ls, lw, label) in curves:
        if d is None:
            continue
        x, E = d['x'], d['E']
        sh = 0.0 if bool(d['metal']) else float(d['vbm'])   # VBM at 0 for an insulator, E_F for a metal
        seg = d['seg']
        for b in E:
            for s in np.unique(seg):
                k = seg == s
                ax.plot(x[k], b[k] - sh, color=color, ls=ls, lw=lw)
        ax.plot([], [], color=color, ls=ls, lw=1.2, label=label)
        if len(d['dos']):
            ad.plot(d['dos'], d['dosE'] - sh, color=color, ls=ls, lw=lw)
        if ref is None:
            ref = d
    if ref is not None:
        for lx, lb in zip(ref['labx'], ref['lab']):
            ax.axvline(float(lx), color='0.75', lw=0.5)
        ax.set_xticks([float(v) for v in ref['labx']]); ax.set_xticklabels([str(v) for v in ref['lab']], fontsize=8)
        ax.set_xlim(float(ref['x'].min()), float(ref['x'].max()))
    ax.axhline(0, color='0.5', lw=0.5, ls=':'); ad.axhline(0, color='0.5', lw=0.5, ls=':')
    ax.set_ylim(-8, 10); ax.set_ylabel('E (eV), VBM (or E_F) = 0')
    ad.set_xlabel('DOS'); ad.set_xticks([]); ad.set_xlim(left=0)
    ax.set_title(title, fontsize=9, loc='left')
    ax.legend(fontsize=7, loc='upper right', framealpha=0.8)
    fig.savefig(path, dpi=80, bbox_inches='tight'); plt.close(fig)


def main():
    db = sys.argv[1]
    figs = '--figs' in sys.argv
    only = None
    if '--only' in sys.argv:
        only = set(sys.argv[sys.argv.index('--only') + 1].split(','))
    S, N, MP, OLD = load()
    D = load_db_logs(db)
    SG = spacegroups()
    os.makedirs(f'{db}/fig', exist_ok=True)
    rows = []; jobs = []
    for m in sorted(S, key=lambda k: (int(S[k]['natom']), S[k]['formula'], k)):
        s, n = S[m], N[m]
        v, adopt = decide(m, S, N, D)
        cat = n['category']
        lda = next((v[k]['lda'] for k in ('N', 'R', 'M') if k in v and v[k]['lda'] is not None), fnum(s['gap_LDA_eV']))
        gap = v[adopt]['gap'] if adopt else None
        others = {k: v[k]['gap'] for k in v if k != adopt and v[k]['gap'] is not None}
        diff = max((abs(g - gap) for g in others.values()), default=0.0) if gap is not None else 0.0
        if adopt is None:
            flag = 'no result'
        elif diff > LARGE:
            flag = 'CHECK'
        elif diff > AGREE:
            flag = 'differs'
        else:
            flag = ''
        tagq = {'N': 'db_qsgw', 'R': rr_tag(db, m), 'M': 'may_qsgw'}
        dq = npz(db, m, tagq[adopt]) if adopt and tagq[adopt] else None
        dl = npz(db, m, 'db_lda') or npz(db, m, 'may_lda')
        di = ''
        if dq is not None and not bool(dq['metal']):
            di = 'D' if abs(float(dq['gap_path']) - float(dq['gap_direct_path'])) < 0.01 else 'I'
        elif dq is not None:
            di = 'metal'
        o = OLD.get(m, {})
        mp = MP.get(m, {})
        rows.append(dict(mpid=m, formula=s['formula'], natom=int(s['natom']), sg=f'{SG.get(m, ("", ""))[1]} ({SG.get(m, ("", ""))[0]})',
                         category=cat, lda=lda, gap=gap, adopt=adopt or '', others=others, flag=flag, di=di,
                         iters={k: v[k].get('iter') for k in v}, mp_pbe=mp.get('band_gap'), ehull=mp.get('energy_above_hull'),
                         icsd=bool((mp.get('database_IDs') or {}).get('icsd')), old2=o.get('shot2'), old1=o.get('shot1'),
                         note=n['note'], host=v.get('N', {}).get('host', '')))
        if figs and (only is None or m in only):
            spec = [('db_lda' if os.path.exists(f'{db}/npz/{m}.db_lda.npz') else 'may_lda', '0.6', '-', 0.7, 'LDA')]
            spec.append((tagq[adopt] if adopt else None, 'tab:red', '-', 0.9, f'QSGW80 ({adopt})'))
            for k, g in others.items():
                if gap is not None and abs(g - gap) > 0.1 and tagq.get(k):
                    spec.append((tagq[k], 'tab:blue', '--', 0.7, f'QSGW80 ({k}) {g:.2f} eV'))
            title = f'{m}  {s["formula"]}  {SG.get(m, ("", ""))[1]}   LDA {fmt(lda)} eV,  QSGW80 {fmt(gap)} eV ({adopt}) {di}'
            jobs.append((db, m, title, spec))
    if jobs:
        from multiprocessing import Pool
        with Pool(8) as p:
            p.map(figjob, jobs, chunksize=8)
    write_pages(db, rows, D)


def figjob(job):
    db, m, title, spec = job
    curves = [(npz(db, m, t) if t else None, c, ls, lw, lab) for (t, c, ls, lw, lab) in spec]
    if any(c[0] is not None for c in curves):
        try:
            figure(db, m, title, curves, f'{db}/fig/{m}.png')
        except Exception as e:
            print(m, 'figure failed', e)


def gapcell(r):
    if r['gap'] is None:
        return '—'
    s = f"**{r['gap']:.2f}** {r['adopt']}"
    for k, g in sorted(r['others'].items()):
        s += f" / {k} {g:.2f}" if abs(g - r['gap']) > AGREE else f" / {k} ≈"
    return s


def write_pages(db, rows, D):
    cnt = Counter(r['adopt'] or 'none' for r in rows)
    fl = Counter(r['flag'] for r in rows)
    with open(f'{db}/gw1500db.tsv', 'w') as f:
        keys = ['mpid', 'formula', 'natom', 'sg', 'category', 'lda', 'gap', 'adopt', 'di', 'flag', 'mp_pbe', 'ehull', 'icsd', 'old1', 'old2', 'host']
        f.write('\t'.join(keys + ['others']) + '\n')
        for r in rows:
            f.write('\t'.join('' if r[k] is None else (f'{r[k]:.4f}' if isinstance(r[k], float) else str(r[k])) for k in keys)
                    + '\t' + ';'.join(f'{k}={g:.4f}' for k, g in sorted(r['others'].items())) + '\n')
    pages = sorted({r['natom'] for r in rows})
    with open(f'{db}/table.md', 'w') as f:
        f.write('# GW1500 QSGW80 — table\n\n[README](README.md) for the conditions M, R, N and how to read this table. '
                'Bands: ' + ', '.join(f'[{n} atoms](bands_{n}atoms.md)' for n in pages) + '.\n\n')
        f.write('| mpid | formula | atoms | space group | LDA (eV) | QSGW80 (eV) | D/I | check | MP PBE | 2025 2shot | category | fig |\n')
        f.write('| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |\n')
        for r in rows:
            link = f"[fig](bands_{r['natom']}atoms.md#{r['mpid']})" if os.path.exists(f"{db}/fig/{r['mpid']}.png") else ''
            mpl = f"[{r['mpid']}](https://next-gen.materialsproject.org/materials/{r['mpid']})"
            f.write(f"| {mpl} | {r['formula']} | {r['natom']} | {r['sg']} | {fmt(r['lda'])} | {gapcell(r)} | {r['di']} | {r['flag']} | "
                    f"{fmt(r['mp_pbe'])} | {fmt(r['old2'])} | {r['category']} | {link} |\n")
    for n in pages:
        with open(f'{db}/bands_{n}atoms.md', 'w') as f:
            f.write(f'# GW1500 QSGW80 — bands and DOS, {n} atoms per cell\n\n[README](README.md) · [table](table.md)\n\n'
                    'Gray: LDA. Red: QSGW80 adopted (condition in parentheses). Blue dashed: another QSGW80 result when its gap differs '
                    'by more than 0.1 eV. Energies from the VBM (insulators) or E_F (metals). Right: total DOS.\n\n')
            for r in rows:
                if r['natom'] == n and os.path.exists(f"{db}/fig/{r['mpid']}.png"):
                    f.write(f'<a id="{r["mpid"]}"></a>\n**{r["mpid"]}** {r["formula"]} — QSGW80 {gapcell(r)} eV {r["di"]} {r["flag"]}\n\n'
                            f'<img src="fig/{r["mpid"]}.png" width="70%">\n\n')
    with open(f'{db}/summary_counts.json', 'w') as f:
        json.dump(dict(adopt=cnt, flag=fl, n=len(rows), db_logged=len(D)), f, indent=1)
    readme(db, rows, D, cnt, fl, pages)


def pairstats(rows, a, b):
    d = []
    for r in rows:
        va = r['gap'] if r['adopt'] == a else r['others'].get(a)
        vb = r['gap'] if r['adopt'] == b else r['others'].get(b)
        if va is not None and vb is not None:
            d.append(va - vb)
    if not d:
        return None
    d = np.array(d)
    return dict(n=len(d), med=float(np.median(np.abs(d))), mean=float(d.mean()), n005=int((np.abs(d) > AGREE).sum()),
                n02=int((np.abs(d) > LARGE).sum()), mx=float(np.abs(d).max()))


def readme(db, rows, D, cnt, fl, pages):
    import datetime
    now = datetime.datetime.now().strftime('%Y-%m-%d %H:%M')
    st = Counter(v['verdict'] for v in D.values())
    with open(f'{db}/README.md', 'w') as f:
        f.write(f'''# GW1500 — QSGW80 band gaps, bands and DOS of 1546 materials (ecalj)

Made {now} by `gw1500db_build.py` (ecalj, `ecalj_auto`). An update of the 2025 database
[tkotani/DOSnpSupplement](https://github.com/tkotani/DOSnpSupplement) (1shot/2shot QSGW of 1516 materials, the supplement of arXiv:2507.19189):
here every material is QSGW80 (scaledsigma = 0.8) iterated to convergence, the conditions of each value are written,
and results of different conditions are compared.

- [table.md](table.md): one row per material (gaps, conditions, direct/indirect, checks). Machine-readable: [gw1500db.tsv](gw1500db.tsv)
- Bands and DOS: ''' + ', '.join(f'[{n} atoms](bands_{n}atoms.md)' for n in pages) + f'''
- Materials: Materials Project, paramagnetic, no lanthanides or actinides, PBE gap > 0, at most 8 atoms per cell
  (structures `ecalj_auto/INPUT/gw1500/POSCARALL`). The structures of MP are not always the experimental ones.

## Conditions

Every gap carries the letter of its conditions. **Bold** is the value adopted: N if the database run converged, else R, else M.

| | when | conditions |
| --- | --- | --- |
''')
        for k, (w, t) in COND.items():
            f.write(f'| **{k}** | {w} | {t} |\n')
        f.write(f'''
Common to all: LDA (VWN) as the starting point and for the LDA column, PMT basis (APW + MTO) made by `ctrlgenToml.py`
(M: its predecessors), no spin polarization, no spin-orbit coupling.

## How to read the QSGW80 column

`**2.96** N / M ≈` : 2.96 eV from N; M agrees within {AGREE} eV. `**3.04** R / M 2.75` : M differs, its value written.
Column *check*: `differs` = some other result differs by {AGREE}–{LARGE} eV; `CHECK` = by more than {LARGE} eV (to be looked at).
D/I: direct or indirect gap along the band path of the figure (the gap value itself is from the k mesh of lmf).

## Status ({now})

| adopted from | materials |
| --- | --- |
''')
        for k in ('N', 'R', 'M', 'none'):
            f.write(f'| {k} | {cnt.get(k, 0)} |\n')
        f.write(f'''
Database run (N) so far: {len(D)} materials logged ({', '.join(f'{k} {v}' for k, v in st.most_common())}).
Order of the run: the materials with trouble in May first (SUSPECT_GOOD, MAY_WRONG, DRIFT_GOOD, UNKNOWN), then the
FAILED/NOTCONV of May (reruns R exist) alternating with a random sample of the May GOOD (by number of atoms), and the
two-atom GOOD on the second machine. The 16 materials with invalid or suspect structures were not run.

**Agreement between the condition sets** (|difference of gaps|, eV)

| pair | materials | median | mean (first − second) | > {AGREE} | > {LARGE} | max |
| --- | --- | --- | --- | --- | --- | --- |
''')
        for a, b in (('N', 'M'), ('N', 'R'), ('R', 'M')):
            p = pairstats(rows, a, b)
            if p:
                f.write(f"| {a} − {b} | {p['n']} | {p['med']:.3f} | {p['mean']:+.3f} | {p['n005']} | {p['n02']} | {p['mx']:.2f} |\n")
        chk = [r for r in rows if r['flag'] == 'CHECK']
        f.write(f'''
## Materials to check ({len(chk)}): results differ by more than {LARGE} eV

| mpid | formula | QSGW80 | category | note |
| --- | --- | --- | --- | --- |
''')
        for r in chk:
            f.write(f"| [{r['mpid']}](bands_{r['natom']}atoms.md#{r['mpid']}) | {r['formula']} | {gapcell(r)} | {r['category']} | {r['note'][:160].replace('|', '/')} |\n")
        f.write('''
## Categories (from the May production, `ecalj_auto/GW1500_status.md`)

GOOD: converged in May. DRIFT_GOOD: converged but the gap drifted over the iterations. SUSPECT_GOOD: converged, but the
iterations oscillated. NOTCONV_MAY: not converged within 10 iterations in May. FAILED_MAY: crashed in May (mostly the
all-TF32 precision of that time). UNKNOWN_MAY: no gap (metals, semimetals). MAY_WRONG: the May value was wrong.
INVALID_STRUCTURE / SUSPECT_STRUCTURE: the MP structure is a fragment cut out of another compound, or doubtful
(listed with their values, which have no physical meaning).
''')


if __name__ == '__main__':
    main()
