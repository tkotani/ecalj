# 2026-10-01 14:4x: add, as mlo_lm3, every semicore local orbital whose band top (lmlo of the default run) is above EF - 17 eV
import sys, re, tomllib, glob, subprocess
m = sys.argv[1]; thr = -17.0
base = f'../mlocheck_20261001/{m}'
lm = {0: '1', 1: '2 3 4', 2: '5 6 7 8 9'}
d = tomllib.load(open(glob.glob(f'{base}/ctrlg.*.toml')[0], 'rb')); sites = [s['atom'] for s in d['site']]
rows = {}
for l in open(f'{base}/lmlo', errors='replace'):
    g = re.search(r'local orbital atom\s+(\d+)\s+l=\s*(\d)\s+top of its band EF\s+(-?[0-9.]+)', l)
    if g and float(g.group(3)) > thr:
        ia, ll = int(g.group(1)), int(g.group(2)); rows.setdefault(ia, set()).add(ll)
blk = ';'.join(f'{ia} {sites[ia-1]} ' + ' '.join(lm[ll] for ll in sorted(ls)) for ia, ls in sorted(rows.items()))
if not blk: print(m, 'no LO above', thr); sys.exit()
subprocess.run(['python3', 'addblock.py', m, f'{m}_lo17', 'mlo_lm3', blk])
