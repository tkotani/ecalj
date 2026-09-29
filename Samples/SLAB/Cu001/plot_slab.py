#!/usr/bin/env python3
"""Bands, total DOS and electrostatic potential of the Cu(001) slab.

  python3 plot_slab.py [cu001]       (after lmf, workfunction.py, job_band and job_tdos)

Writes slab_<sname>.png and slab_<sname>.npz (the plotted numbers, the source directory and the date).
  bnd00?.spin1      : band index, x along the path, E - Ef (eV), ...     (job_band)
  dos.tot.<sname>   : E - Ef (Ry), DOS (states/Ry/cell)                  (job_tdos)
  estaticpot.dat    : ix iy iz, smooth part of the electrostatic potential (eV)   (lmf)
"""
import datetime
import glob
import os
import re
import sys
import tomllib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

RYDBERG = 13.605693
EMIN, EMAX = -10.0, 5.0            # eV from Ef: window of the bands and of the DOS
sname = sys.argv[1] if len(sys.argv) > 1 else 'cu001'
now = datetime.datetime.now().strftime('%Y-%m-%d %H:%M')
save = {'date': now, 'source': os.getcwd()}
fig, (axb, axd, axv) = plt.subplots(1, 3, figsize=(11, 4.2), gridspec_kw={'width_ratios': [2.2, 1, 2]})

# (a) bands: one block for each band in each file
xb, eb = [], []
for f in sorted(glob.glob('bnd[0-9][0-9][0-9].spin1')):
    d = np.array([[float(x) for x in l.split()[:3]] for l in open(f) if l.strip() and l[0] != '#'])
    for ib in np.unique(d[:, 0]):
        b = d[d[:, 0] == ib]
        if b[:, 2].max() < EMIN or b[:, 2].min() > EMAX:
            continue
        axb.plot(b[:, 1], b[:, 2], 'k-', lw=0.5)
        xb.append(np.append(b[:, 1], np.nan))
        eb.append(np.append(b[:, 2], np.nan))
save['band_x'], save['band_e_minus_ef_ev'] = np.concatenate(xb), np.concatenate(eb)
glt = open('bandplot.isp1.glt').read()                      # labels and positions of the symmetry points
ticks = re.findall(r"'(\w+)'\s+([0-9.]+)", glt)
for name, x in ticks:
    axb.axvline(float(x), color='0.6', lw=0.4)
axb.set_xticks([float(x) for _, x in ticks], [r'$\bar\Gamma$' if n == 'Gamma' else r'$\bar{\rm %s}$' % n for n, _ in ticks])
axb.axhline(0, color='0.6', lw=0.4)
axb.set_xlim(0, float(ticks[-1][1]))
axb.set_ylim(EMIN, EMAX)
axb.set_ylabel(r'$E - E_{\rm F}$ (eV)')
axb.set_title('(a) bands', fontsize=9)
save['band_ticks_x'] = np.array([float(x) for _, x in ticks])
save['band_ticks_label'] = np.array([n for n, _ in ticks])

# (b) total DOS
d = np.loadtxt(f'dos.tot.{sname}')
e, dos = d[:, 0] * RYDBERG, d[:, 1] / RYDBERG
sel = (e >= EMIN) & (e <= EMAX)
axd.plot(dos[sel], e[sel], 'k-', lw=0.5)
axd.axhline(0, color='0.6', lw=0.4)
axd.set_ylim(EMIN, EMAX)
axd.set_xlim(0, None)
axd.set_xlabel('DOS (states/eV/cell)')
axd.set_yticklabels([])
axd.set_title('(b) total DOS', fontsize=9)
save['dos_e_minus_ef_ev'], save['dos_states_per_ev_cell'] = e[sel], dos[sel]

# (c) electrostatic potential on the line through the origin, measured from Ef
with open(f'ctrlg.{sname}.toml', 'rb') as f:
    ctrl = tomllib.load(f)
unit = ctrl['struc']['alat'] * 0.529177                     # alat in angstrom
c = ctrl['struc']['plat'][2][2] * unit
zlayer = np.array([s['pos'][2] for s in ctrl['site']]) * unit
wf = dict(l.split() for l in open(f'workfunction.{sname}.txt') if l[0] != '#')
ef, vvac = float(wf['fermi_energy']), float(wf['vacuum_level'])
rows = np.array([[float(x) for x in l.split()[:4]] for l in open('estaticpot.dat') if l.strip() and l[0] != '#'])
line = rows[(rows[:, 0] == 1) & (rows[:, 1] == 1)]
z = (line[:, 2] - 1) / rows[:, 2].max() * c
axv.plot(z, line[:, 3] - ef, 'k-', lw=0.5)
axv.axhline(0, color='0.6', lw=0.4)
axv.axhline(vvac - ef, color='k', lw=0.4, ls='--')
for zl in zlayer:
    axv.plot([zl, zl], [-19.5, -18], 'k-', lw=0.8)
axv.annotate('atomic layers', (zlayer.mean(), -17.5), ha='center', fontsize=8)
axv.annotate(f'vacuum level: work function {vvac - ef:.2f} eV', (c, vvac - ef), xytext=(0, 4),
             textcoords='offset points', ha='right', fontsize=8)
axv.annotate(r'$E_{\rm F}$', (c, 0), xytext=(0, 3), textcoords='offset points', ha='right', fontsize=8)
axv.set_xlim(0, c)
axv.set_ylim(-20, 8)
axv.set_xlabel(r'$z$ ($\rm\AA$)')
axv.set_ylabel(r'$V_{\rm es} - E_{\rm F}$ (eV)')
axv.set_title('(c) smooth electrostatic potential on the line x = y = 0', fontsize=9)
save['pot_z_angstrom'], save['pot_minus_ef_ev'] = z, line[:, 3] - ef
save['ef_ev'], save['vacuum_level_ev'], save['layers_z_angstrom'] = ef, vvac, zlayer

fig.suptitle(f'Cu(001) slab, 4 layers, LDA ({sname})   [{now}]', fontsize=9)
fig.tight_layout()
fig.savefig(f'slab_{sname}.png', dpi=150)
np.savez_compressed(f'slab_{sname}.npz', **save)
print(f'slab_{sname}.png and slab_{sname}.npz written')
