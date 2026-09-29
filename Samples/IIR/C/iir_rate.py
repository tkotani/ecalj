#!/usr/bin/env python3
'''Impact ionization rate from QPU_life of gw_lmfh (hqpe).
usage: python3 iir_rate.py [--dir .] [--title "diamond, 6x6x6"]

QPU_life has one line for each state (q, band):
  columns 1-3: q (2 pi/alat),  4: band index,  5: energy of the state (eV, from the Fermi energy of heftet),
  column 6: FWHM = 2 Z Im Sigma (eV) at the energy of the state; > 0 for occupied states, < 0 for unoccupied states.
The rate is  1/tau = |FWHM| / hbar.
Output:
  iir.txt : q, band, energy from the valence band maximum of the mesh (eV), FWHM (eV), rate (1/s)
  iir.png : rate against the energy (log scale); states with |FWHM| < 1e-7 eV (below the threshold) are not drawn
The numbers of the figure are those of iir.txt.
'''
import argparse, datetime, os, sys
HBAR = 0.6582119569      # eV fs

p = argparse.ArgumentParser(description='impact ionization rate from QPU_life')
p.add_argument('--dir', default='.', help='directory with QPU_life')
p.add_argument('--title', default='diamond', help='title of the figure')
args = p.parse_args()
src = os.path.join(os.path.abspath(args.dir), 'QPU_life')
rows = []
for line in open(src):
    w = line.split()
    if len(w) < 6 or line.lstrip().startswith('#'): continue
    rows.append((float(w[0]), float(w[1]), float(w[2]), int(w[3]), float(w[4]), float(w[5])))
if not rows: sys.exit('iir_rate.py: no data in ' + src)
# the Fermi energy of an insulator lies in the gap: band edges among the states of the mesh
evbm = max(r[4] for r in rows if r[4] < 0.)
ecbm = min(r[4] for r in rows if r[4] > 0.)
gap = ecbm - evbm
made = datetime.datetime.now().strftime('%Y-%m-%d %H:%M')
with open('iir.txt', 'w') as f:
    f.write('# impact ionization rate 1/tau = |FWHM|/hbar, FWHM = 2 Z Im Sigma; made by iir_rate.py on %s\n' % made)
    f.write('# from %s\n' % src)
    f.write('# band edges on this mesh, from the Fermi energy (eV): VBM %.3f CBM %.3f ; gap %.3f\n' % (evbm, ecbm, gap))
    f.write('#   q1       q2       q3    band  E-E_VBM(eV)    FWHM(eV)      rate(1/s)\n')
    for q1, q2, q3, ib, e, w in rows:
        f.write('%8.5f %8.5f %8.5f %4d %10.3f %14.10f %14.6e\n' % (q1, q2, q3, ib, e - evbm, w, abs(w)/HBAR*1e15))
print('iir.txt written: %d states, gap on the mesh = %.3f eV' % (len(rows), gap))

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except ImportError:
    sys.exit('iir_rate.py: matplotlib not found; iir.png is not written')
INK, MUTED, SURFACE = '#0b0b0b', '#52514e', '#fcfcfb'
series = [('holes (occupied states)',      '#2a78d6', 'o', lambda e: e <= 0.),
          ('electrons (unoccupied states)', '#eb6834', 's', lambda e: e > 0.)]
fig, ax = plt.subplots(figsize=(6.4, 4.0), facecolor=SURFACE)
ax.set_facecolor(SURFACE)
for label, color, marker, sel in series:
    pts = [(e - evbm, abs(w)/HBAR*1e15) for _, _, _, _, e, w in rows if sel(e - evbm) and abs(w) > 1e-7]
    ax.plot([x for x, _ in pts], [y for _, y in pts], marker, ms=4.5, mfc=color, mec=SURFACE, mew=0.6, ls='none', label=label)
ax.set_yscale('log')
ymin, ymax = ax.get_ylim()
ax.axvspan(0., gap, color='0.90', zorder=0, lw=0)                       # gap: no states
for x in (-gap, 2.*gap):                                                # thresholds of energy conservation
    ax.axvline(x, color=MUTED, lw=0.6, ls='--', zorder=1)
ax.annotate('gap', (gap/2., ymin), xytext=(0, 4), textcoords='offset points', ha='center', va='bottom', fontsize=8, color=MUTED)
ax.annotate('$-E_g$', (-gap, ymax), xytext=(-3, -3), textcoords='offset points', ha='right', va='top', fontsize=8, color=MUTED)
ax.annotate('$2E_g$', (2.*gap, ymax), xytext=(3, -3), textcoords='offset points', ha='left', va='top', fontsize=8, color=MUTED)
ax.set_ylim(ymin, ymax)
ax.set_xlabel('energy from the valence band maximum (eV)', color=INK)
ax.set_ylabel('impact ionization rate (1/s)', color=INK)
ax.grid(True, which='major', color='0.88', lw=0.5)
ax.set_axisbelow(True)
for s in ax.spines.values(): s.set_color('0.6'); s.set_linewidth(0.6)
ax.tick_params(labelsize=8, colors=MUTED, width=0.6)
ax.legend(fontsize=8, frameon=False, loc='lower left', labelcolor=INK)
ax.set_title('%s: gap %.2f eV on the mesh (%s)' % (args.title, gap, made), fontsize=9, color=INK)
fig.tight_layout()
fig.savefig('iir.png', dpi=150, facecolor=SURFACE)
print('iir.png written; data: iir.txt')
