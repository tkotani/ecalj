#!/usr/bin/env python3
"""Empty spheres (ES) at the large voids of a material, before its LDA and QSGW (2026-10-04, GW1500; user: bands tied to the
vacuum level, nearly free states in the voids, need ES, and ES are a more natural choice than more basis on the atoms).

    gw1500_addes.py <sname> <rmin> <bindir>        (in the directory of ctrlg.<sname>.toml)

<rmin> = auto (2026-10-05, user: CO2, CO are molecular): 2.0 a.u. for a molecular crystal (gw1500_molecular.py on POSCAR of the
directory), else 3.0 a.u.

1. ctrlg_addes.py <sname> --rmin <rmin>: ES at every local maximum of the void radius above rmin (a.u.), [[spec]] E and a row
   "<i> E 1 2 3 4" in [mlo] mlo_lm
2. the GW part of the ctrlg ([gw], [mlo], [blocks], [product_basis]) written again from the ctrl part with the ES by
   <bindir>/ctrlgenToml.py --addgw (gwinit), so that the per-atom tables of [product_basis] list the ES; [ham] scaledsigma is kept
Prints "ES <n>" (n = 0: no void above rmin, nothing changed). ctrlg_addes.py is taken from $ADDES, this directory, ../SRC/exec.
"""
import os, re, shutil, subprocess, sys
m, rmin, B = sys.argv[1], sys.argv[2], sys.argv[3]
here = os.path.dirname(os.path.abspath(__file__))
if rmin == 'auto':
    r = subprocess.run([sys.executable, f'{here}/gw1500_molecular.py', 'POSCAR'], capture_output=True, text=True)
    mol = r.stdout.split()[1:2] == ['molecular']
    rmin = '2.0' if mol else '3.0'
    print(f'{"molecular" if mol else "extended"} crystal: rmin {rmin} a.u.')
custom = os.path.join(os.environ.get('ES_CUSTOM_DIR', f'{here}/es_custom'), f'ctrlg.{m}.toml')
if os.path.exists(custom):
    # 2026-10-07 03:11: ES placed by hand for this material (the E* tables of [[site]] and [[spec]] in es_custom/ctrlg.<sname>.toml),
    # e.g. just outside the surface of a slab (PbS, SnS: ES at the centre of the vacuum missed the state the MLO needs) or below
    # the void threshold (C3N4, a void of 2.96 a.u.). The atoms of the custom file must be those of this ctrlg.
    import tomllib
    f = f'ctrlg.{m}.toml'; t = open(f).read(); tc = open(custom).read()
    at = [(x['atom'], x.get('pos')) for x in tomllib.loads(t)['site']]
    atc = [(x['atom'], x.get('pos')) for x in tomllib.loads(tc)['site'] if not x['atom'].startswith('E')]
    if len(at) != len(atc) or any(a != b or (p is not None and max(abs(u - v) for u, v in zip(p, q)) > 1e-4)
                                  for (a, p), (b, q) in zip(at, atc)):
        print(f'{custom}: its atoms differ from {f}'); sys.exit(1)
    blocks = lambda head, key: ''.join(re.findall(r'(?ms)^' + re.escape(head) + r'[^\n]*\n' + key + r' *= "E[0-9]*"\n.*?(?=^\[|^# ===|\Z)', tc))
    sites, specs = blocks('[[site]]', 'atom'), blocks('[[spec]]', 'atom')

    def after_last(t, head, add):
        i = [mm.start() for mm in re.finditer(r'(?m)^' + re.escape(head), t)][-1]
        mm = re.search(r'(?m)^(\[|# ===)', t[i + len(head):])
        e = i + len(head) + (mm.start() if mm else len(t) - i - len(head))
        return t[:e] + add + t[e:]
    open(f + '.bak_addes', 'w').write(t)
    t = after_last(after_last(t, '[[site]]', sites), '[[spec]]', specs)
    open(f, 'w').write(t); tomllib.loads(t)
    n = sites.count('[[site]]')
    print(f'ES from {custom}: {n} sites')
else:
    addes = next(p for p in (os.environ.get('ADDES', ''), f'{here}/ctrlg_addes.py', f'{here}/../SRC/exec/ctrlg_addes.py', f'{B}/ctrlg_addes.py')
                 if p and os.path.exists(p))
    r = subprocess.run([sys.executable, addes, m, '--rmin', rmin], capture_output=True, text=True)
    print(r.stdout.strip())
    if r.returncode != 0:
        print(r.stderr[-500:]); sys.exit(1)
    n = len(re.findall(r'(?m)^ES site', r.stdout))
if n == 0:
    print('ES 0'); sys.exit(0)
f = f'ctrlg.{m}.toml'; t = open(f).read()
w = '_addes_gw'; shutil.rmtree(w, ignore_errors=True); os.makedirs(w)
open(f'{w}/{f}', 'w').write(t[:t.index('\n# === GW') + 1])
p = subprocess.run([f'{B}/ctrlgenToml.py', m, '--addgw'], cwd=w, capture_output=True, text=True, env={**os.environ, 'OMP_NUM_THREADS': '1'})
tn = open(f'{w}/{f}').read()
if '[product_basis]' not in tn or tn.count('\n[[site]]') != t.count('\n[[site]]'):
    print('GW part not written:', (p.stdout + p.stderr)[-500:]); sys.exit(1)
open(f, 'w').write(tn); shutil.rmtree(w)
print(f'ES {n}')
