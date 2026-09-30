#!/usr/bin/env python3
"""Table of the anomalous Hall conductivity sigma_xy from the files of hahc / hx0ahc.py.

Reads, in the current directory,
  ahc_tet.isp11.dat, ahc_tet.isp22.dat : tetrahedron weights, spin channels 1 and 2
  ahc_sp.isp11.dat,  ahc_sp.isp22.dat  : direct sum over the k points
(columns: 1 dEf(eV), 2-10 sigma_xx, xy, xz, yx, yy, yz, zx, zy, zz in Ohm^-1 cm^-1;
 isp22 exists for two spin channels, ham.so = 0 or 2; for ham.so = 1 there is isp11 only)
and save.<sname> of lmf.  Writes
  ahc.txt : dEf, sigma_xy of the tetrahedron method (spin 1, spin 2, sum) and of the direct sum (sum of the spins)
  scf.txt : magnetic moment and total energy of the first converged SCF (line c of save.<sname>)
Usage: python3 ahc_table.py [sname]      (sname = fe by default)
"""
import os, sys, re, datetime

def read_sigma(fname):
    """{dEf: [sigma_xx, sigma_xy, ..., sigma_zz]} of one file; lines starting with # are headers."""
    out = {}
    if not os.path.exists(fname): return out
    for line in open(fname):
        if not line.strip() or line.lstrip().startswith('#'): continue
        v = [float(x) for x in line.split()]
        out[round(v[0], 7)] = v[1:10]
    return out

sname = sys.argv[1] if len(sys.argv) > 1 else 'fe'
IXY = 1   # index of sigma_xy in xx, xy, xz, yx, ...
tet = [read_sigma(f'ahc_tet.isp{i}{i}.dat') for i in (1, 2)]
sp  = [read_sigma(f'ahc_sp.isp{i}{i}.dat')  for i in (1, 2)]
if not tet[0]: sys.exit('ahc_table.py: ahc_tet.isp11.dat is missing or empty')

with open('ahc.txt', 'w') as f:
    f.write('# Anomalous Hall conductivity sigma_xy (Ohm^-1 cm^-1) against the shift dEf (eV) of the Fermi energy\n')
    f.write('# made by ahc_table.py on %s in %s\n' % (datetime.datetime.now().strftime('%Y-%m-%d %H:%M'), os.getcwd()))
    f.write('# tet: tetrahedron weights, direct: sum over the k points; spin2 = 0 when there is one channel (ham.so=1)\n')
    f.write('#  dEf(eV)    tet_spin1    tet_spin2      tet_sum   direct_sum\n')
    for e in sorted(tet[0]):
        t1 = tet[0][e][IXY]
        t2 = tet[1][e][IXY] if e in tet[1] else 0.
        d  = sp[0][e][IXY] + (sp[1][e][IXY] if e in sp[1] else 0.)
        f.write('%9.4f %12.4f %12.4f %12.4f %12.4f\n' % (e, t1, t2, t1 + t2, d))

# first converged SCF: the later c lines of save.<sname> are restarts by job_band and job_AHC
for line in open(f'save.{sname}'):
    if line.startswith('c '):
        mmom = float(re.search(r'mmom=\s*(\S+)', line).group(1))
        ehf  = float(re.search(r'ehf\(eV\)=\s*(\S+)', line).group(1))
        ehk  = float(re.search(r'ehk\(eV\)=\s*(\S+)', line).group(1))
        with open('scf.txt', 'w') as f:
            f.write('# first converged SCF of save.%s: magnetic moment (mu_B), Harris-Foulkes and Kohn-Sham energies (eV)\n' % sname)
            f.write('%8.4f %16.6f %16.6f\n' % (mmom, ehf, ehk))
        break
else:
    sys.exit(f'ahc_table.py: no converged SCF (line c) in save.{sname}')
print(open('ahc.txt').read(), end='')
print(open('scf.txt').read(), end='')
