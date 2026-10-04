#!/usr/bin/env python3
# 2026-10-04 GW1500 MLO: rewrite the [mlo] section of <dir>/ctrlg.<m>.toml with the gwinit rules of the ecalj in use
# (ctrlgenToml --addgw on a copy; the other sections are kept). For runs whose ctrlg was written before the present rules
# (the semicore local orbitals of 2026-10-01, mlo_lm2). The old file is kept as ctrlg.<m>.toml.orig. Usage: gw1500_mlo_regen.py <dir> <m>
import os, sys, shutil, subprocess
d, m = sys.argv[1], sys.argv[2]
BIN = os.environ.get('ECALJ_BIN', os.path.expanduser('~/bin'))
def region(t):
    a = t.index('\n[mlo]\n') + 1; b = t.index('\n# === BLOCKS', a) + 1; return a, b
f = os.path.join(d, f'ctrlg.{m}.toml'); t = open(f).read(); a, b = region(t)
w = os.path.join(d, '_regen'); shutil.rmtree(w, ignore_errors=True); os.makedirs(w)
open(os.path.join(w, os.path.basename(f)), 'w').write(t[:t.index('\n# === GW') + 1])
p = subprocess.run([f'{BIN}/ctrlgenToml.py', m, '--addgw'], cwd=w, capture_output=True, text=True, env={**os.environ, 'OMP_NUM_THREADS': '1'})
try:
    tn = open(os.path.join(w, os.path.basename(f))).read(); an, bn = region(tn)
except Exception:
    print('FAILED', p.returncode, (p.stdout + p.stderr)[-400:]); sys.exit(1)
if not os.path.exists(f + '.orig'): shutil.copy(f, f + '.orig')
open(f, 'w').write(t[:a] + tn[an:bn] + t[b:]); shutil.rmtree(w)
print('ok')
