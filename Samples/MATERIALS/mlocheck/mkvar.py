# 2026-10-01 MLO variants: copy a finished material and add mlo_lm2 / mlo_lm3 blocks, then job_mlo only.
import sys, os, shutil, tomllib, subprocess
base, name, kind = sys.argv[1], sys.argv[2], sys.argv[3]   # kind: lm2sp (EH2 s,p on all sites) | lm3d (LO d on sites with a semicore d LO)
src = f'../mlocheck_20261001/{base}'
if os.path.exists(name): shutil.rmtree(name)
shutil.copytree(src, name, symlinks=True)
f = [x for x in os.listdir(name) if x.startswith('ctrlg.') and x.endswith('.toml')][0]
t = open(f'{name}/{f}').read(); d = tomllib.loads(t)
sites = [s['atom'] for s in d['site']]; spec = {s['atom']: s for s in d['spec']}
if kind == 'lm2sp':
    rows = ''.join(f'{i} {a:4s} 1 2 3 4\n' for i, a in enumerate(sites, 1)); key = 'mlo_lm2'
elif kind == 'lm3semi':   # every semicore LO channel (pz shell below the valence shell), s, p or d
    import math
    lmr = {0: '1', 1: '2 3 4', 2: '5 6 7 8 9'}
    rows = ''
    for i, a in enumerate(sites, 1):
        pz = spec[a].get('pz', []); p = spec[a].get('p', [])
        ls = [l for l in range(min(len(pz), 3)) if pz[l] > 0 and int(pz[l] % 10) < (int(p[l]) if l < len(p) and p[l] else 99)]
        if ls: rows += f'{i} {a:4s} ' + ' '.join(lmr[l] for l in ls) + '\n'
    key = 'mlo_lm3'
elif kind == 'lm3d':
    rows = ''.join(f'{i} {a:4s} 5 6 7 8 9\n' for i, a in enumerate(sites, 1) if len(spec[a].get('pz', [])) > 2 and spec[a]['pz'][2] > 0); key = 'mlo_lm3'
if kind == 'both':   # lm2sp and lm3semi together
    import math
    lmr = {0: '1', 1: '2 3 4', 2: '5 6 7 8 9'}
    rows3 = ''
    for i, a in enumerate(sites, 1):
        pz = spec[a].get('pz', []); p = spec[a].get('p', [])
        ls = [l for l in range(min(len(pz), 3)) if pz[l] > 0 and int(pz[l] % 10) < (int(p[l]) if l < len(p) and p[l] else 99)]
        if ls: rows3 += f'{i} {a:4s} ' + ' '.join(lmr[l] for l in ls) + '\n'
    rows = ''.join(f'{i} {a:4s} 1 2 3 4\n' for i, a in enumerate(sites, 1)); key = 'mlo_lm2'
    if rows3: t = t.replace('mlo_lm = """\n', f'mlo_lm3 = """\n{rows3}"""\nmlo_lm = """\n', 1)
t = t.replace('mlo_lm = """\n', f'{key} = """\n{rows}"""\nmlo_lm = """\n', 1)
open(f'{name}/{f}', 'w').write(t)
sname = f[6:-5]
for x in os.listdir(name):   # no stale model of the base run
    if x.startswith('band_MLO') or x in ('lmlo', 'HamRsMLO', 'ljob_mlo', 'lwriteham') or x.endswith('.png'): os.remove(os.path.join(name, x))
env = dict(os.environ, OMP_NUM_THREADS='1')
subprocess.run(f'{os.path.expanduser("~/bin/job_mlo")} {sname} -np 2 --nognuplot > ljob_mlo 2>&1', shell=True, cwd=name, env=env)
print(name, key, rows.replace('\n', '; '))
