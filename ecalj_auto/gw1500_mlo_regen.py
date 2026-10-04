#!/usr/bin/env python3
# 2026-10-04 GW1500 MLO: rewrite the [mlo] section of <dir>/ctrlg.<m>.toml with the gwinit rules of the ecalj in use
# (ctrlgenToml --addgw on a copy; the other sections are kept). For runs whose ctrlg was written before the present rules
# (the semicore local orbitals of 2026-10-01, mlo_lm2). The old file is kept as ctrlg.<m>.toml.orig. Usage: gw1500_mlo_regen.py <dir> <m> [--lm2 | --lm2all]
# --lm2 (2026-10-04): baseline 2 of ecaljdoc mlo section 9, the mlo_lm2 rows written by gwinit (EH2 s,p of the cations) taken
# --lm2all (2026-10-04, a test; user: a full model for every material, by adding basis): EH2 s,p of every atom that is not a
#   transition metal Sc-Cu, Y-Ag, La-Au or Ac- (anions and H too; NaN3: max 0.47 -> 0.05 eV with N added)
import os, re, sys, shutil, subprocess
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
mlo = tn[an:bn]
if '--lm2all' in sys.argv:
    EL = ('H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru '
          'Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn '
          'Fr Ra Ac Th Pa U').split()
    i = mlo.index('mlo_lm = """'); j = mlo.index('"""', i + 13)
    rows = []
    for l in mlo[i + 13:j].split('\n'):
        c = l.split()
        if len(c) < 2 or c[0].startswith('!'):
            continue
        el = re.match(r'[A-Z][a-z]?', c[1]).group(0)
        z = EL.index(el) + 1 if el in EL else 0
        if z == 0 or 21 <= z <= 29 or 39 <= z <= 47 or 57 <= z <= 79 or z >= 89:
            continue
        rows.append(f'{c[0]} {c[1]}   1 2 3 4')
    body = 'mlo_lm2 = """\n' + '\n'.join(rows) + '\n"""\n'
    if 'mlo_lm2 = """' in mlo:
        a2 = mlo.index('mlo_lm2 = """'); b2 = mlo.index('"""', a2 + 13) + 4
        mlo = mlo[:a2] + body + mlo[b2:]
    else:
        mlo = mlo.rstrip('\n') + '\n' + body
elif '--lm2' in sys.argv and 'mlo_lm2 = """' in mlo:
    i = mlo.index('mlo_lm2 = """'); j = mlo.index('"""', i + 14)
    mlo = mlo[:i] + '\n'.join(l[2:] if l.startswith('! ') else l for l in mlo[i:j].split('\n')) + mlo[j:]
open(f, 'w').write(t[:a] + mlo + t[b:]); shutil.rmtree(w)
print('ok')
