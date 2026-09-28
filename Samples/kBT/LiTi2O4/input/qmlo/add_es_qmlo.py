#!/usr/bin/env python3
"""Add empty spheres to the LiTi2O4 MLO-QSGW input (2026-09-28, user: 9^3 tf32 with empty spheres from LDA).

  add_es_qmlo.py merge <ctrlg in> <ctrlg with ES, GW sections from gwinit> <ctrlg out>
  add_es_qmlo.py sites <ctrlg in> <ctrlg out>

The empty sites of the spinel (Fd-3m, origin 1: Li 8a, Ti 16d, O 32e) in the primitive cell, as in the 2026-09-23 trial
(kt1 mlo_es_lda): E1 at 16c (1/8,1/8,1/8)-type, 4 sites, r = 1.30; E2 at 8b (1/2,1/2,1/2) and (3/4,3/4,3/4), r = 1.45.
Each gets s and p MTOs (lmx = 1) and s + p MLOs (mlo_lm 1 2 3 4), so the MLO model grows from 126 to 150 orbitals.

  sites: the [[site]] blocks after the last site, the [[spec]] blocks before '# === BZ', and the ES lines at the end of
         mlo_lm; nbas / nspec, if [struc] still has them, are removed (the counts are those of the tables; lmf stops on
         them since 2026-09-28, and before that a stale count silently dropped the added sites).
         The other sections (the [gw] settings of the run, [mlo], [product_basis] of the 14 atoms) are kept as they are.
  merge: the per-atom product-basis rows of the new atoms (pb_lcutmx entries, nlx / valence / core rows with iatom > 14)
         taken from a copy that gwinit has regenerated (lmfa -> lmf --jobgw=0 -> gwinit) are appended to the tables of
         <ctrlg in> (the output of 'sites'), so that the 14 atoms keep exactly the tables of the production input.
"""
import re
import sys

SITES = [("E1", (0.125, 0.125, 0.125)), ("E1", (-0.125, -0.125, 0.125)), ("E1", (-0.125, 0.125, -0.125)),
         ("E1", (0.125, -0.125, -0.125)), ("E2", (0.5, 0.5, 0.5)), ("E2", (0.75, 0.75, 0.75))]
SPECS = '''[[spec]]   # @4  empty sphere at 16c (2026-09-28, add_es_qmlo.py)
atom   = "E1"
z      = 0
r      = 1.30
lmx    = 1
lmxa   = 2
rsmh   = [0.65, 0.65]
eh     = [-0.3, -0.3]

[[spec]]   # @5  empty sphere at 8b
atom   = "E2"
z      = 0
r      = 1.45
lmx    = 1
lmxa   = 2
rsmh   = [0.72, 0.72]
eh     = [-0.3, -0.3]

'''


def add_sites(s):
    nat = len(re.findall(r'^\[\[site\]\]', s, re.M))
    assert nat == 14, f'expected the 14 atoms of LiTi2O4, found {nat}'
    blk = ''.join(f'[[site]]   # @{nat + 1 + k}  empty sphere\natom = "{a}"\npos  = [{q[0]}, {q[1]}, {q[2]}]\n\n'
                  for k, (a, q) in enumerate(SITES))
    i = s.index('# === SPEC')
    s = s[:i] + blk + s[i:]
    j = s.index('# === BZ')
    s = s[:j] + SPECS + s[j:]
    s = re.sub(r'^\s*(nspec|nbas)\s*=.*\n', '', s, flags=re.M)   # [struc] must not carry them (lmf stops; 2026-09-28):
                                                                   # the counts are those of the tables
    lm = ''.join(f'{nat + 1 + k} {a}   1 2 3 4\n' for k, (a, q) in enumerate(SITES))
    m = re.search(r'^mlo_lm = """\n(.*?)^"""', s, re.M | re.S)
    s = s[:m.end(1)] + lm + s[m.end(1):]
    return s


def table(s, key):
    """[start, end) of the list 'key = [ ... ]' in the [product_basis] section (rows one per line)."""
    m = re.search(r'^' + key + r' = \[\n', s, re.M)
    e = s.index('\n]', m.end())
    return m.end(), e


def merge(a, b):
    pa = a.index('[product_basis]'); pb = b.index('[product_basis]')
    A, B = a[pa:], b[pb:]
    # pb_lcutmx: keep the 14 entries of A, append the ES entries of B
    la = re.search(r'^pb_lcutmx\s*=\s*\[([^\]]*)\]', A, re.M); lb = re.search(r'^pb_lcutmx\s*=\s*\[([^\]]*)\]', B, re.M)
    va = [v.strip() for v in la.group(1).split(',')]; vb = [v.strip() for v in lb.group(1).split(',')]
    assert len(va) == 14 and len(vb) == 20, (len(va), len(vb))
    A = A[:la.start(1)] + ', '.join(va + vb[14:]) + A[la.end(1):]
    for key in ('nlx', 'valence', 'core'):
        s0, e0 = table(B, key)
        rows = [r for r in B[s0:e0].split('\n') if re.match(r'\s*\[\s*(\d+)', r) and int(re.match(r'\s*\[\s*(\d+)', r).group(1)) > 14]
        s1, e1 = table(A, key)
        body = A[s1:e1].rstrip('\n')
        if rows:
            if body.strip() and not body.rstrip().endswith(','):
                body += ','
            body += '\n' + '\n'.join(r if r.rstrip().endswith(',') else r.rstrip() + ',' for r in rows)
        A = A[:s1] + body + A[e1:]
        print(f'{key}: {len(rows)} rows for the empty spheres')
    return a[:pa] + A


if __name__ == '__main__':
    if sys.argv[1] == 'sites':
        open(sys.argv[3], 'w').write(add_sites(open(sys.argv[2]).read()))
    elif sys.argv[1] == 'merge':
        open(sys.argv[4], 'w').write(merge(open(sys.argv[2]).read(), open(sys.argv[3]).read()))
    else:
        sys.exit(__doc__)
