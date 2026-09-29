#!/usr/bin/env python3
'''Tables and a figure of the total DOS (job_tdos) and the partial DOS (job_pdos).
usage: python3 dos_table.py <sname> [--emin -15] [--emax 15] [--title "ZnS"] [--sites Zn,S] [--dir .]

Input (written by lmf; energies in Ry from the Fermi energy, or from the top of the valence band for an insulator):
  dos.tot.<sname>              : energy, total DOS of spin 1 (and spin 2)      (states/Ry/cell)
  dosi.tot.<sname>             : energy, number of states below the energy
  dos.isp<s>.site<iii>.<sname> : energy, DOS in the sphere of the site iii for the 25 real harmonics
                                 lm = 1 (s), 2-4 (p), 5-9 (d), 10-16 (f), 17-25 (g)   (states/Ry)
Output (energies in eV, DOS in states/eV, only emin <= E <= emax):
  tdos.txt   : energy, total DOS for each spin, number of states for each spin
  pdos_l.txt : energy, then s p d f g (sums over m) for each site and spin
  dos.png    : total DOS and partial DOS of each site; spin 2 is drawn downwards
The numbers of the figure are those of tdos.txt and pdos_l.txt.
'''
import argparse, datetime, glob, os, re, sys
import numpy as np
RY = 13.605693
LNAME = ['s', 'p', 'd', 'f', 'g']
LRANGE = [(1, 2), (2, 5), (5, 10), (10, 17), (17, 26)]      # columns of the pdos files for each l (column 0 is the energy)

p = argparse.ArgumentParser(description='tables of the total and partial DOS in eV')
p.add_argument('sname')
p.add_argument('--emin', type=float, default=-15., help='lower end of the tables (eV)')
p.add_argument('--emax', type=float, default=15., help='upper end of the tables (eV)')
p.add_argument('--title', default=None, help='title of the figure')
p.add_argument('--sites', default=None, help='names of the sites for the figure, separated by commas (e.g. Zn,S)')
p.add_argument('--dir', default='.', help='directory with the dos files')
args = p.parse_args()
M = args.sname
src = os.path.abspath(args.dir)
made = datetime.datetime.now().strftime('%Y-%m-%d %H:%M')
head = '# made by dos_table.py on %s from %s\n' % (made, src)

def window(a):
    e = a[:, 0]*RY
    return a[(e >= args.emin) & (e <= args.emax)]

# --- total DOS
t = window(np.loadtxt(os.path.join(src, f'dos.tot.{M}')))
ti = window(np.loadtxt(os.path.join(src, f'dosi.tot.{M}')))
nsp = t.shape[1] - 1
with open('tdos.txt', 'w') as f:
    f.write('# total DOS (states/eV/cell) and number of states below E, for each spin; E (eV) from the Fermi energy\n')
    f.write('# (from the top of the valence band for an insulator)\n' + head)
    f.write('#   E(eV)  ' + ''.join('   dos_spin%d' % (i+1) for i in range(nsp)) + ''.join('  nstates_spin%d' % (i+1) for i in range(nsp)) + '\n')
    for r, ri in zip(t, ti):
        f.write('%10.4f' % (r[0]*RY) + ''.join('%12.5f' % (x/RY) for x in r[1:]) + ''.join('%15.5f' % x for x in ri[1:]) + '\n')

# --- partial DOS summed over m
files = sorted(glob.glob(os.path.join(src, f'dos.isp1.site*.{M}')))
if not files: sys.exit(f'dos_table.py: no dos.isp1.site*.{M} in {src} (run job_pdos)')
sites = [int(re.search(r'site(\d+)', os.path.basename(x)).group(1)) for x in files]
pdos = {}                                                   # (site, spin) -> array (npoint, 5)
for ib in sites:
    for isp in range(1, nsp+1):
        a = window(np.loadtxt(os.path.join(src, 'dos.isp%d.site%03d.%s' % (isp, ib, M)), comments='#'))
        epdos = a[:, 0]*RY
        pdos[ib, isp] = np.array([a[:, i:j].sum(axis=1) for i, j in LRANGE]).T/RY
with open('pdos_l.txt', 'w') as f:
    f.write('# partial DOS in the spheres (states/eV), summed over m; E (eV) as in tdos.txt\n' + head)
    f.write('# columns: E, then s p d f g of site 1 spin 1, (site 1 spin 2,) site 2 spin 1, ...\n')
    f.write('#   E(eV)  ' + ''.join(' %s_site%d_sp%d' % (l, ib, isp) for ib in sites for isp in range(1, nsp+1) for l in LNAME) + '\n')
    for i, e in enumerate(epdos):
        f.write('%10.4f' % e + ''.join('%12.5f' % x for ib in sites for isp in range(1, nsp+1) for x in pdos[ib, isp][i]) + '\n')
print('tdos.txt: %d energies, pdos_l.txt: %d energies, %d site(s), %d spin(s)' % (len(t), len(epdos), len(sites), nsp))

# --- figure
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except ImportError:
    sys.exit('dos_table.py: matplotlib not found; dos.png is not written')
INK, MUTED, SURFACE = '#0b0b0b', '#52514e', '#fcfcfb'
COLOR = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100']       # s, p, d, f in this order for every site
sign = [1., -1.]                                           # spin 2 downwards
sitenames = args.sites.split(',') if args.sites else None
fig, axes = plt.subplots(len(sites)+1, 1, figsize=(6.4, 2.2*(len(sites)+1)), sharex=True, facecolor=SURFACE)
ax = axes[0]
for isp in range(nsp):
    ax.plot(t[:, 0]*RY, sign[isp]*t[:, 1+isp]/RY, '-', color=INK, lw=0.8)
ax.annotate('total', (0.01, 0.92), xycoords='axes fraction', fontsize=8, color=INK, va='top')
for ax, ib in zip(axes[1:], sites):
    for il in range(4):                                    # g is small: in pdos_l.txt only
        for isp in range(1, nsp+1):
            y = sign[isp-1]*pdos[ib, isp][:, il]
            ax.plot(epdos, y, '-', color=COLOR[il], lw=0.8, label=LNAME[il] if isp == 1 else None)
        y1 = pdos[ib, 1][:, il]
        if y1.max() > 0.15*max(pdos[ib, 1].max(), 1e-9):   # direct label at the highest peak of the larger channels
            i = int(y1.argmax())
            ax.annotate(LNAME[il], (epdos[i], y1[i]), xytext=(4, 0), textcoords='offset points', fontsize=8, color=INK, va='center')
    name = 'site %d' % ib + (' (%s)' % sitenames[ib-1] if sitenames and ib <= len(sitenames) else '')
    ax.annotate(name, (0.01, 0.92), xycoords='axes fraction', fontsize=8, color=INK, va='top')
    lo, hi = ax.get_ylim()
    ax.set_ylim(lo, hi + 0.1*(hi - lo))                     # room for the labels at the peaks
    ax.legend(fontsize=7, frameon=False, loc='upper right', ncol=4, labelcolor=INK)
for ax in axes:
    ax.set_facecolor(SURFACE)
    ax.axvline(0., color=MUTED, lw=0.5, ls='--')
    ax.axhline(0., color='0.6', lw=0.5)
    ax.grid(True, color='0.90', lw=0.5)
    ax.set_axisbelow(True)
    ax.set_xlim(args.emin, args.emax)
    ax.set_ylabel('states/eV', fontsize=8, color=INK)
    ax.tick_params(labelsize=8, colors=MUTED, width=0.6)
    for s in ax.spines.values(): s.set_color('0.6'); s.set_linewidth(0.6)
axes[-1].set_xlabel('energy (eV)', color=INK, fontsize=9)
axes[0].set_title('%s: total and partial DOS%s (%s)' % (args.title or M, ', spin 2 downwards' if nsp == 2 else '', made), fontsize=9, color=INK)
fig.tight_layout()
fig.savefig('dos.png', dpi=150, facecolor=SURFACE)
print('dos.png written; data: tdos.txt, pdos_l.txt')
