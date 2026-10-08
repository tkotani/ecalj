#!/usr/bin/env python3
"""ES rule es-d for a slab (a vacuum layer in the cell), 2026-10-08 (user: the new rule; PbS mp-561320, ecaljdoc mlo.md table M11).

    gw1500_es_slab_layers.py <in ctrlg> <out ctrlg> [--r 2.6] [--nlayer 1] [--nolayer] [--nop]

<in ctrlg> has the ES of a slab placed by hand (the es_custom file of gw1500_addes.py): species "E" on the next layer outside both
surfaces and "E2" at the centre of the vacuum. Written: --nlayer more ES layers "E3" (radius --r a.u.) on each side, the surface ES moved
along the normal of the slab evenly towards the centre ES (one layer: half way), so that the vacuum is filled (the empty gap between the layers about 2 a.u.
or less), and the basis s,p (lmx 1, rsmh and eh for l = 0, 1) on every ES. The overlap of all spheres (nearest images) is checked.
"""
import argparse, itertools, re, sys
import numpy as np
import tomllib

ap = argparse.ArgumentParser()
ap.add_argument('src'); ap.add_argument('out')
ap.add_argument('--r', type=float, default=2.6, help='radius of the added ES layer (a.u.)')
ap.add_argument('--nlayer', type=int, default=1, help='ES layers added on each side, evenly spaced between the surface ES and the centre')
ap.add_argument('--nolayer', action='store_true', help='no added layer (s,p only)')
ap.add_argument('--nop', action='store_true', help='keep the basis of the ES (s only)')
a = ap.parse_args()

src = open(a.src).read(); d = tomllib.loads(src)
al = d['struc']['alat']; P = np.array(d['struc']['plat'], float) * al
n = np.cross(P[0], P[1]); n /= np.linalg.norm(n); c = abs(P[2] @ n)        # normal of the slab: the plane of the first two vectors
sites = [(s['atom'], np.array(s['pos'], float) * al if 'pos' in s else np.array(s['xpos'], float) @ P) for s in d['site']]
h = lambda x: (x @ n) % c
cen = [x for s, x in sites if s == 'E2']
if not cen or not any(s == 'E' for s, _ in sites):
    sys.exit('gw1500_es_slab_layers: needs ES "E" outside the surfaces and "E2" at the vacuum centre')
hc = h(cen[0])
new = []
for s, x in sites:
    if s == 'E':
        dz = hc - h(x)
        if abs(dz) > c / 2: dz -= np.sign(dz) * c
        new += [x + n * dz * i / (a.nlayer + 1) for i in range(1, a.nlayer + 1)]
t = src
if not a.nolayer:
    blk = ''.join(f'[[site]]   # ES layer between the surface ES and the vacuum centre (rule es-d, gw1500_es_slab_layers.py)\natom = "E3"\n'
                  f'pos  = [{x[0]/al:.8f}, {x[1]/al:.8f}, {x[2]/al:.8f}]\n# ----------------------------------------------------------------\n' for x in new)
    i = t.index('# === SPEC'); t = t[:i] + blk + '\n' + t[i:]
    j = t.index('# === BZ')
    t = t[:j] + (f'[[spec]]   # ES layer (rule es-d)\natom   = "E3"\nz      = 0\nr      = {a.r:.2f}\nlmx    = 0\nlmxa   = 2\n'
                 f'rsmh   = [{a.r/2:.3f}]\neh     = [-0.3]\n\n') + t[j:]
if not a.nop:
    for e in ('E', 'E2') + (() if a.nolayer else ('E3',)):
        pat = re.compile(r'(atom\s*=\s*"%s"\s*\n(?:(?!\[\[).*\n)*?)lmx\s*=\s*0\s*\n((?:(?!\[\[).*\n)*?)rsmh\s*=\s*\[([0-9.]+)\]\s*\n'
                         r'((?:(?!\[\[).*\n)*?)eh\s*=\s*\[([-0-9.]+)\]' % e)
        t, k = pat.subn(lambda m: f'{m.group(1)}lmx    = 1\n{m.group(2)}rsmh   = [{m.group(3)}, {m.group(3)}]\n{m.group(4)}'
                        f'eh     = [{m.group(5)}, {m.group(5)}]', t)
        if k != 1: sys.exit(f'gw1500_es_slab_layers: the species {e} with lmx = 0 not found once ({k})')
o = tomllib.loads(t)
S = [(s['atom'], np.array(s['pos'], float) * al if 'pos' in s else np.array(s['xpos'], float) @ P) for s in o['site']]
R = {s['atom']: s['r'] for s in o['spec']}
L = np.array(list(itertools.product(range(-2, 3), repeat=3))) @ P
worst = max(((R[p] + R[q]) / np.min(np.linalg.norm(y - x + L, axis=1)[np.linalg.norm(y - x + L, axis=1) > 1e-6]) - 1, p, q)
            for i, (p, x) in enumerate(S) for j, (q, y) in enumerate(S) if j >= i)
if worst[0] > 0:
    sys.exit(f'gw1500_es_slab_layers: spheres overlap {100*worst[0]:.1f} % ({worst[1]}, {worst[2]}); give a smaller --r')
open(a.out, 'w').write(t)
hs = sorted(round(h(x) / 1.8897, 2) for s, x in S if s.startswith('E'))
print(f'{a.out}: {sum(1 for s, _ in S if s.startswith("E"))} ES, heights along the normal (A) {hs}; largest overlap {100*worst[0]:.1f} %')
