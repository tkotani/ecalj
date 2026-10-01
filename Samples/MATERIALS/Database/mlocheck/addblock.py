# 2026-10-01 14:3x: copy a finished material and insert one block (mlo_lm2 or mlo_lm3) with the given rows, then job_mlo only.
# usage: addblock.py <base> <name> <key> "<row>;<row>..."
import sys, os, shutil, subprocess
base, name, key, rows = sys.argv[1:5]
src = f'../mlocheck_20261001/{base}'
if os.path.exists(name): shutil.rmtree(name)
shutil.copytree(src, name, symlinks=True)
for x in os.listdir(name):
    if x.startswith('band_MLO') or x in ('lmlo', 'HamRsMLO', 'ljob_mlo', 'lwriteham') or x.endswith('.png'): os.remove(os.path.join(name, x))
f = [x for x in os.listdir(name) if x.startswith('ctrlg.') and x.endswith('.toml')][0]
t = open(f'{name}/{f}').read()
block = '\n'.join(r.strip() for r in rows.split(';')) + '\n'
t = t.replace('mlo_lm = """\n', f'{key} = """\n{block}"""\nmlo_lm = """\n', 1)
open(f'{name}/{f}', 'w').write(t)
subprocess.run(f'{os.path.expanduser("~/bin/job_mlo")} {f[6:-5]} -np 4 --nognuplot > ljob_mlo 2>&1', shell=True, cwd=name, env=dict(os.environ, OMP_NUM_THREADS='1'))
print(name, key, block.replace('\n', '; '))
