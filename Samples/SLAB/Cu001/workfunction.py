#!/usr/bin/env python3
"""Vacuum level and work function of the slab from estaticpot.dat and efermi.lmf (written by lmf).

  python3 workfunction.py cu001          -> workfunction.cu001.txt

estaticpot.dat holds the smooth part of the electrostatic potential (eV) on lines along the third
lattice vector; outside the muffin-tin spheres it is the electrostatic potential itself. Its zero is
that of the eigenvalues and of the Fermi energy in efermi.lmf. The vacuum level is the value at the
middle of the vacuum (the middle of the widest gap between the atomic layers), eq. (1) of README.md.
"""
import sys
import tomllib
import numpy as np

RYDBERG = 13.605693       # eV, as in the programs
sname = sys.argv[1] if len(sys.argv) > 1 else 'cu001'
with open(f'ctrlg.{sname}.toml', 'rb') as f:
    ctrl = tomllib.load(f)
c = ctrl['struc']['plat'][2][2]                      # length of the third lattice vector (units of alat)
z = np.sort(np.array([s['pos'][2] for s in ctrl['site']]) % c)
gap = np.diff(np.append(z, z[0] + c))                # distances between the layers, the last one across the cell boundary
i = np.argmax(gap)
zvac = (z[i] + gap[i] / 2) % c                       # middle of the vacuum

rows = np.array([[float(x) for x in l.split()[:4]] for l in open('estaticpot.dat') if l.strip() and l[0] != '#'])
n3 = int(rows[:, 2].max())
line = rows[(rows[:, 0] == 1) & (rows[:, 1] == 1)]   # the line through the origin
zz = (line[:, 2] - 1) / n3 * c
iv = np.argmin(abs(zz - zvac))
vvac = line[iv, 3]
sel = abs(zz - zvac) < 0.1 * gap[i]                  # flatness of the potential around the middle of the vacuum
flat = line[sel, 3].max() - line[sel, 3].min()
ef = float(open('efermi.lmf').readline().split()[0].replace('D', 'e')) * RYDBERG
unit = ctrl['struc']['alat'] * 0.529177              # alat in angstrom
with open(f'workfunction.{sname}.txt', 'w') as f:
    f.write('# vacuum level, Fermi energy and work function (eV); flatness = max - min of the potential\n'
            f'# within +-10% of the vacuum width around its middle, z = {zvac * unit:.4f} angstrom\n'
            f'vacuum_level   {vvac:12.5f}\n'
            f'fermi_energy   {ef:12.5f}\n'
            f'work_function  {vvac - ef:12.5f}\n'
            f'flatness       {flat:12.5f}\n')
print(open(f'workfunction.{sname}.txt').read(), end='')
