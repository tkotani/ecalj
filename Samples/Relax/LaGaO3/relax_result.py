#!/usr/bin/env python3
'''Collect the results of a relaxation run of lmf: energies and forces at the end of the self-consistency of
each geometry.
usage: python3 relax_result.py <sname> [log of lmf]     (default log: llmf)
Reads  save.<sname> : lines h, i (iterations) and c, x (end of the self-consistency of one geometry; x = the
                      number of iterations iter.nit was reached before conv)
       the log      : tables "Forces, with eigenvalue correction" (mRy/bohr) and "Maximum Harris force"
Writes relax_energy.txt : geometry, number of iterations, ehf and ehk (eV), c or x
       relax_force.txt  : geometry, maximum force, and the total force on every site (mRy/bohr)
The positions of every step are in AtomPos.<sname> (written by lmf, in units of alat); the number of sites is read from it.
'''
import re, sys
m = sys.argv[1]
log = sys.argv[2] if len(sys.argv) > 2 else 'llmf'
# --- energies from save.<sname>
geo, nit = [], 0
for line in open(f'save.{m}'):
    r = re.match(r'([hicx]) .*ehf\(eV\)=\s*(\S+)\s+ehk\(eV\)=\s*(\S+)', line)
    if not r: continue
    nit += 1
    if r.group(1) in 'cx':
        geo.append((nit, float(r.group(2)), float(r.group(3)), r.group(1))); nit = 0
with open('relax_energy.txt', 'w') as f:
    f.write(f'# {m}: end of the self-consistency of each geometry, from save.{m}\n')
    f.write('# geometry iterations       ehf(eV)            ehk(eV)     c=converged x=iter.nit reached\n')
    for i, (n, ehf, ehk, tag) in enumerate(geo):
        f.write(f'  {i+1:4d} {n:8d} {ehf:18.6f} {ehk:18.6f}   {tag}\n')
# --- forces from the log: the last table before each relaxation step
L = open(log).read().split('\n')
nbas = int(next(l for l in open(f'AtomPos.{m}') if '!nbas' in l).split()[0])
res, table, fmax = [], None, None
for i, l in enumerate(L):
    if 'Forces, with eigenvalue correction' in l:
        rows, j = {}, i + 1
        while j < len(L) and len(rows) < nbas:     # lines of other ranks may come in between the rows
            r = re.match(r'\s*(\d+)((?:\s*-?\d+\.\d\d){9})\s*$', L[j])
            if r: rows[int(r.group(1))] = [float(x) for x in re.findall(r'-?\d+\.\d\d', r.group(2))][-3:]
            j += 1
        table = rows
    r = re.search(r'Maximum Harris force =\s*(\S+) mRy/au \(site\s*(\d+)', l)
    if r: fmax = (float(r.group(1)), int(r.group(2)))
    if 'Updated atom positions' in l and table is not None:
        res.append((fmax, table)); table = None
if table is not None: res.append((fmax, table))
with open('relax_force.txt', 'w') as f:
    f.write(f'# {m}: forces at the end of the self-consistency of each geometry (mRy/bohr), from {log}\n')
    for g, (fm, tab) in enumerate(res):
        # the site of the largest force goes to a comment line (test.py does not compare it): sites that are equal by
        # symmetry have the same force, and which of them lmf names depends on the last digits (2026-09-30: 12 with
        # gfortran, 11 with nvfortran)
        f.write(f'# geometry {g+1}: site, force x y z.  The largest force is on site {fm[1]}\n')
        f.write(f'maxforce {g+1} {fm[0]:12.4f}\n')
        for ib in sorted(tab):
            f.write(f'  {g+1} {ib:4d} {tab[ib][0]:10.2f} {tab[ib][1]:10.2f} {tab[ib][2]:10.2f}\n')
print(open('relax_energy.txt').read(), end='')
print(open('relax_force.txt').read(), end='')
