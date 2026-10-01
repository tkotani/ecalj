#!/usr/bin/env python3
# 2026-10-01 21:1x: regenerate the [mlo] section of every ctrlg of Samples/MATERIALS with the present gwinit
# (through ctrlgenToml --addgw, for the same layout), without touching the other sections. Writes <rel>.new and <rel>.diff
# under OUT; applying them is a separate step.
import os, re, glob, shutil, subprocess, difflib, sys
REPO = '/home/takao/ecalj/Samples/MATERIALS'; OUT = '/home/takao/work/regen_mlo/out'
BIN = os.path.expanduser('~/bin')
def mlo_region(t):
    a = t.index('\n[mlo]\n') + 1
    b = t.index('\n# === BLOCKS', a) + 1
    return a, b
files = sorted(f for f in glob.glob(f'{REPO}/**/ctrlg.*.toml', recursive=True) if '/mlocheck/' not in f)
for f in files:
    rel = os.path.relpath(f, REPO); d = os.path.join(OUT, os.path.dirname(rel)); os.makedirs(d, exist_ok=True)
    sname = os.path.basename(f)[6:-5]
    t = open(f).read()
    try: a, b = mlo_region(t)
    except ValueError: print(rel, 'NO [mlo] REGION'); continue
    w = os.path.join(d, 'work'); shutil.rmtree(w, ignore_errors=True); os.makedirs(w)
    i = t.index('\n# === GW') + 1
    open(os.path.join(w, os.path.basename(f)), 'w').write(t[:i])
    p = subprocess.run([f'{BIN}/ctrlgenToml.py', sname, '--addgw'], cwd=w, capture_output=True, text=True, env={**os.environ, 'OMP_NUM_THREADS': '1'})
    nf = os.path.join(w, os.path.basename(f))
    try:
        tn = open(nf).read(); an, bn = mlo_region(tn)
    except Exception as e:
        print(rel, 'FAILED', p.returncode, (p.stdout + p.stderr)[-300:]); continue
    new = t[:a] + tn[an:bn] + t[b:]
    open(os.path.join(OUT, rel + '.new'), 'w').write(new)
    diff = ''.join(difflib.unified_diff(t[a:b].splitlines(True), tn[an:bn].splitlines(True), 'old', 'new', n=0))
    open(os.path.join(OUT, rel + '.diff'), 'w').write(diff)
    print(rel, 'ok', len(diff.splitlines()), 'diff lines', flush=True)
