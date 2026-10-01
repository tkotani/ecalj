# 2026-10-01: MLO test of Samples/MATERIALS at the DFT level.
# Copies ctrlg.<sname>.toml of each material into <work>/<Material>/ and rewrites its [mlo] section uniformly:
#   mlo_lm : per site, Z<=10 -> s,p (1..4); Z>=11 -> s,p,d (1..9); lanthanides Z=57..71 -> s,p,d,f (1..16)
#   mlo_method=4, mlo_delta=2.0, mlo_w=2.0 (the defaults), mlo_nkabc = [bz] nkabc
#   ham.rdsig -> 0 (DFT level, no sigm)
import os, re, glob, shutil, tomllib, json
REPO = os.path.expanduser('~/ecalj/Samples/MATERIALS')
WORK = os.path.dirname(os.path.abspath(__file__))
src = {os.path.basename(os.path.dirname(f)): f for f in glob.glob(f'{REPO}/*/ctrlg.*.toml')}   # (the 62 were under Database/ until 2026-10-01 15:4x)
src['La2CuO4'] = glob.glob(f'{REPO}/La2CuO4/ctrlg.*.toml')[0]
src['InAsGaSb_n4'] = glob.glob(f'{REPO}/InAsGaSb/n4/ctrlg.*.toml')[0]
src['InAsGaSb_n10'] = glob.glob(f'{REPO}/InAsGaSb/n10/ctrlg.*.toml')[0]
src['BaTiO3'] = glob.glob(f'{REPO}/BaTiO3/QSGW5run/ctrlg.*.toml')[0]

def lms(z, pf=None):
    # f is listed when the atom has f: the lanthanides, and an atom whose valence f is the 4f (P_f = 4.x, Hf in HfO2;
    # added by hand to HfO2 on 2026-10-01 07:15 after the first run with spd only, which was 0.057 eV off at the VBM)
    if 57 <= z <= 71 or (pf is not None and int(pf) == 4): return range(1, 17)
    return range(1, 5) if z <= 10 else range(1, 10)

info = {}
for mat, f in sorted(src.items()):
    sname = os.path.basename(f)[6:-5]
    t = open(f).read()
    d = tomllib.loads(t)
    zof = {s['atom']: s['z'] for s in d['spec']}
    pfof = {s['atom']: (s.get('p', [None]*4) + [None]*4)[3] for s in d['spec']}
    rows = [f"{i} {s['atom']:4s} " + ' '.join(map(str, lms(zof[s['atom']], pfof[s['atom']]))) for i, s in enumerate(d['site'], 1)]
    nk = d['bz']['nkabc']
    nk = nk if isinstance(nk, list) else [nk] * 3
    mlo = ('[mlo]\n'
           '# 2026-10-01 mlocheck: rewritten by prep.py (Z<=10 sp, Z>=11 spd, lanthanides spdf)\n'
           'mlo_method = 4\nmlo_delta  = 2.0\nmlo_w      = 2.0\n'
           f'mlo_nkabc  = [{nk[0]}, {nk[1]}, {nk[2]}]\n'
           'mlo_lm = """\n' + '\n'.join(rows) + '\n"""\n\n')
    t2, n = re.subn(r'(?ms)^\[mlo\]\n.*?(?=^# === |^\[(?!mlo\])[a-z_]+\]\s*$)', lambda m: mlo, t, count=1)
    if n != 1: raise SystemExit(f'{mat}: [mlo] section not found')
    t2 = re.sub(r'(?m)^(rdsig\s*=\s*)\d+', r'\g<1>0', t2)
    wd = os.path.join(WORK, mat)
    os.makedirs(wd, exist_ok=True)
    open(os.path.join(wd, f'ctrlg.{sname}.toml'), 'w').write(t2)
    tomllib.loads(t2)  # still valid TOML
    info[mat] = dict(sname=sname, src=f, nsite=len(d['site']), nspin=d['ham'].get('nspin', 1), so=d['ham'].get('so', 0),
                     nkabc=nk, f=[s['atom'] for s in d['spec'] if 57 <= s['z'] <= 71])
json.dump(info, open(os.path.join(WORK, 'materials.json'), 'w'), indent=1)
for k, v in sorted(info.items(), key=lambda x: x[1]['nsite']):
    print(f"{k:14s} {v['sname']:14s} nsite={v['nsite']:3d} nspin={v['nspin']} so={v['so']} nk={v['nkabc']} f={v['f']}")
