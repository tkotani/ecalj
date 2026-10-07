#!/usr/bin/env python3
"""A dated snapshot of the GW1500 database for publishing (2026-10-07; user: a database of today's date with its own cover, and the data
needed to reproduce it: the input ctrlg, the band data, the gaps and the conditions).

    gw1500db_snapshot.py <db dir> <inputs dir> <out dir> <label>

<db dir>      made by gw1500db_build.py <db> --figs (README.md, table.md, bands_*.md, fig/, *.tsv, npz/)
<inputs dir>  the inputs of the adopted runs, collected from the machines: {N,E,R}/<mpid>/ctrlg.<mpid>.toml and
              {N,E,R}/<mpid>/PlotBand/syml.<mpid> (R/<mpid>/run.txt names the rerun), M/<mpid>/ (the May input, legacy ctrl)
<out dir>     written: README.md (the cover of this snapshot: the database README with a section on the files), table.md,
              bands_*.md, gw1500db.tsv, history.tsv, later_list.tsv, summary_counts.json, fig/*.png,
              inputs/<mpid>/ (ctrlg.<mpid>.toml, syml.<mpid>, source.txt), bands/<mpid>.npz (QSGW80, LDA and the MLO model)
<label>       e.g. 2026-10-07
"""
import sys, os, csv, shutil, glob
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gw1500db_extract import gaps

db, inp, out, label = sys.argv[1:5]
os.makedirs(out, exist_ok=True)
for f in ('table.md', 'gw1500db.tsv', 'history.tsv', 'later_list.tsv', 'summary_counts.json'):
    if os.path.exists(f'{db}/{f}'): shutil.copy2(f'{db}/{f}', out)
for f in glob.glob(f'{db}/bands_*.md'): shutil.copy2(f, out)
os.makedirs(f'{out}/fig', exist_ok=True)
for f in glob.glob(f'{db}/fig/*.png'): shutil.copy2(f, f'{out}/fig/')

rows = list(csv.DictReader(open(f'{db}/gw1500db.tsv'), delimiter='\t'))
hist = [h for h in csv.DictReader((l for l in open(f'{db}/history.tsv') if not l.startswith('#')), delimiter='\t')]
eh = {h['mpid']: h for h in hist if h['set'] == 'E' and h['verdict'].startswith('CONVERGED')}

def rr_tag(m):
    for r in ('run3x', 'run3', 'run2m', 'run1', 'try'):
        if os.path.exists(f'{db}/npz/{m}.rr_{r}.npz'): return f'rr_{r}'

def load(m, tag):
    f = f'{db}/npz/{m}.{tag}.npz'
    return np.load(f) if tag and os.path.exists(f) else None

nin = nb = 0
for r in rows:
    m, a = r['mpid'], r['adopt']
    if not a: continue
    # inputs
    src = {'N': 'N', 'E': 'E', 'R': 'R', 'M': 'M'}[a]
    d = f'{inp}/{src}/{m}'
    if os.path.isdir(d):
        o = f'{out}/inputs/{m}'; os.makedirs(o, exist_ok=True)
        for f in glob.glob(f'{d}/ctrlg.{m}.toml') + glob.glob(f'{d}/PlotBand/syml.{m}') + glob.glob(f'{d}/ctrl.{m}') + glob.glob(f'{d}/ctrlG.{m}.toml'):
            shutil.copy2(f, o)
        how = {'N': 'the database run N (2026-10-02..07)', 'E': 'the run with empty spheres E (2026-10-05..07)',
               'R': f'the rerun R ({open(d + "/run.txt").read().strip() if os.path.exists(d + "/run.txt") else ""}, 2026-09-30..10-02)',
               'M': 'the May 2026 production M (legacy ctrl of before the TOML input)'}[a]
        if a == 'E' and m in eh and eh[m]['es'] == '0':
            how += ' (the rule found no void: no ES, the same input as N run again with the later ecalj)'
        open(f'{o}/source.txt', 'w').write(f'{m}: input of {how}; ecalj {r["version"]}; host {eh[m]["host"] if a == "E" and m in eh else r["host"]}\n'
                                           'ctrlg.<mpid>.toml is the input as it was at the end of the run (gwsc may have set keys such as readp);\n'
                                           'syml.<mpid> is the band path of the figure. See README.md, "Files for reproducing".\n')
        nin += 1
    # bands: QSGW80 of the adopted run, LDA, the MLO model
    tq = {'E': 'es_qsgw', 'N': 'db_qsgw', 'R': rr_tag(m), 'M': 'may_qsgw_fixed' if os.path.exists(f'{db}/npz/{m}.may_qsgw_fixed.npz') else 'may_qsgw'}[a]
    tl = 'es_lda' if a == 'E' else ('db_lda' if os.path.exists(f'{db}/npz/{m}.db_lda.npz') else 'may_lda')
    parts = {'qsgw80': load(m, tq), 'lda': load(m, tl), 'mlo': load(m, 'db_mlo')}
    pack = {}
    def fnum(s):
        try: return float(s)
        except (TypeError, ValueError): return None
    gq, gl = fnum(r['gap']), fnum(r['lda'])
    for k, v in parts.items():
        if v is None: continue
        for key in v.files: pack[f'{k}_{key}'] = v[key]
        if k in ('qsgw80', 'lda') and 'x' in v.files:
            # the gap keys again by the rule of the table (gw1500db_build.pathgap): npz files extracted before 2026-10-07 19:09
            # have metal=True for 61 insulators (gaps of gw1500db_extract)
            g = gq if k == 'qsgw80' else gl
            for key, val in gaps(v['x'], v['E'], insul=g is not None and g > 0.05).items(): pack[f'{k}_{key}'] = np.array(val)
    if pack:
        os.makedirs(f'{out}/bands', exist_ok=True)
        pack['source'] = np.array(f'{m}: qsgw80 from {tq}, lda from {tl}, mlo from db_mlo (adopted set {a}); energies in eV from the VBM '
                                  '(or E_F for metals); E[band, k] along x with the segments seg and labels lab at labx; dos on dosE')
        np.savez_compressed(f'{out}/bands/{m}.npz', **pack); nb += 1

# the cover of this snapshot: the README of the database with a section on the files
t = open(f'{db}/README.md').read()
sec = f'''
## Files for reproducing (this snapshot, {label})

- `inputs/<mpid>/ctrlg.<mpid>.toml`: the input of the run whose QSGW80 value is adopted (the bold value of the table), as it was at the
  end of that run; `syml.<mpid>`: its band path; `source.txt`: which run (set N, E, R or M), the ecalj version and the machine.
  Run it with the ecalj of that version: `lmfa`, then `gwsc`/`gwscconv` with the conditions of the README (`--ssig` 0.8 is in the ctrlg).
  The May value M (one material) has the legacy `ctrl.<mpid>` of before the TOML input. {nin} materials.
- `bands/<mpid>.npz`: the band data of the figure: `qsgw80_*` (the adopted run), `lda_*`, `mlo_*` (the MLO model and its check),
  each with `x`, `E` (bands × k, eV from the VBM or E_F), `seg`, `lab`/`labx` (symmetry points), and the DOS `dosE`/`dos`. {nb} materials.
  Read with `numpy.load`.
- `history.tsv`: every run of every material; `later_list.tsv`: the materials left for later (and runs still RUNNING).
- The raw run directories (wave functions, self-energies) are not published; they are kept by the maintainers.
'''
i = t.find('\n## ')
t = (t[:i] + '\n' + sec + t[i:]) if i > 0 else t + sec
t = t.replace('# GW1500', f'# GW1500 ({label})', 1) if '# GW1500' in t[:200] else f'<!-- snapshot {label} -->\n' + t
open(f'{out}/README.md', 'w').write(t)
print(f'snapshot {label}: {out}: inputs {nin}, bands {nb}, figures {len(glob.glob(out + "/fig/*.png"))}')
