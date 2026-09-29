#!/usr/bin/env python3
'''Collect the results of an LDA+U run of lmf into result.<sname>.txt (small, for reading and for test.py).
usage: python3 ldau_result.py <sname>
Reads  save.<sname>      : first iteration (h) and the last line of the self-consistency (c, or x when not converged)
       mmom.<sname>.chk  : charge and spin moment in the MT spheres
       orbitalmom.chk    : orbital moments (written when so is not 0). The name carries no <sname>, so the file is
                           copied to orbitalmom.<sname>.chk at the first call: run this script right after lmf
       dmats.<sname>     : density matrix of the LDA+U shell (first blocks: spherical harmonics, spin 1 and 2)
'''
import os, re, shutil, sys
m = sys.argv[1]
orbfile = f'orbitalmom.{m}.chk'
if not os.path.exists(orbfile): shutil.copy('orbitalmom.chk', orbfile)
out = [f'# {m}: written by ldau_result.py from save.{m}, mmom.{m}.chk, {orbfile}, dmats.{m}',
       '# total spin moment mmom (Bohr magneton) and total energy ehk (eV): first iteration (h), end of the self-consistency (c)']
first = last = None
for line in open(f'save.{m}'):
    r = re.match(r'([hicx]) mmom=\s*(\S+).*ehk\(eV\)=\s*(\S+)', line)
    if not r: continue
    if r.group(1) == 'h': first = r.groups()
    if r.group(1) in 'cx': last = r.groups()
if first is None or last is None: sys.exit(f'ldau_result.py: no h or c/x line in save.{m}')
for tag, mm, e in (first, last):
    out.append(f'{tag}  {float(mm):9.4f}  {float(e):16.6f}')
if last[0] == 'x': out.append('# WARNING: x = not converged within iter.nit')
# site moments
spin = [l.split() for l in open(f'mmom.{m}.chk') if l.strip() and not l.startswith('#')]
orb = [float(l.split(':')[1]) for l in open(orbfile) if 'total orbital moment' in l]
out.append('# site, spin moment in the MT sphere and orbital moment (Bohr magneton)')
for s, o in zip(spin, orb):
    out.append(f'site {s[0]} {s[4]:3s} {float(s[2]):10.4f} {o:10.4f}')
# occupations: diagonal of the real part of the first two blocks of dmats
blocks, cur = [], None
for line in open(f'dmats.{m}'):
    if line.startswith('#'): break                    # the rest of the file is not read by lmf either
    if line.startswith('%'):
        cur = []; blocks.append((line, cur)); continue
    if line.strip() and cur is not None: cur.append([float(x) for x in line.split()])
out.append('# occupation of the LDA+U orbitals: diagonal of the density matrix (spherical harmonics, m=-l..l)')
for head, rows in blocks:
    n = len(rows)//2                                   # real part, then imaginary part
    r = re.search(r'l=\s*(\d+)\s+site\s+(\d+)\s+spin\s+(\d+)', head)
    occ = ' '.join(f'{rows[i][i]:8.4f}' for i in range(n))
    out.append(f'occ l={r.group(1)} site={r.group(2)} spin={r.group(3)} {occ}   sum {sum(rows[i][i] for i in range(n)):8.4f}')
open(f'result.{m}.txt', 'w').write('\n'.join(out) + '\n')
print('\n'.join(out))
