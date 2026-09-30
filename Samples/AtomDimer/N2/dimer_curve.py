#!/usr/bin/env python3
"""Binding curve of a homonuclear dimer in a box from the runs of test.py.
usage: dimer_curve.py <sname of the dimer> <sname of the atom> <d1> <d2> ...   (d in A; the runs are in d_<d>/ and atom/)
Reads the last converged line (c) of save.<sname> in each directory: ehk(eV).
Writes dimer.txt: E(d) of the dimer, E of the atom, and the parabola through the three points nearest to the minimum:
  r_e (A), E_min (eV), the binding energy D_e = 2 E_atom - E_min (eV); dimer.png (the curve and the parabola), dimer.npz.
(2026-09-30, Samples/AtomDimer)"""
import sys, re, numpy as np
mol, atom = sys.argv[1], sys.argv[2]; ds = [float(x) for x in sys.argv[3:]]
def ehk(f):
    last = None
    for l in open(f):
        m = re.match(r'[cx] .*ehk\(eV\)=\s*(\S+)', l)
        if m: last = float(m.group(1))
    if last is None: sys.exit(f'{f}: no converged line')
    return last
E = np.array([ehk(f'd_{d:.3f}/save.{mol}') for d in ds]); Ea = ehk(f'atom/save.{atom}')
i = int(np.argmin(E)); sel = slice(max(0, i-1), min(len(ds), i+2))
c = np.polyfit(np.array(ds)[sel], E[sel], 2)            # E = c0 d^2 + c1 d + c2
re_, emin = -c[1]/(2*c[0]), np.polyval(c, -c[1]/(2*c[0]))
de = 2*Ea - emin
with open('dimer.txt', 'w') as f:
    f.write(f'# {mol}: total energy ehk (eV) of the dimer at the bond length d (A), and of the atom\n')
    for d, e in zip(ds, E): f.write(f'  {d:8.4f} {e:16.6f}\n')
    f.write(f'# atom {atom}\n  {Ea:16.6f}\n# parabola through the points nearest to the minimum: r_e (A), E_min (eV), D_e = 2 E_atom - E_min (eV)\n')
    f.write(f'  {re_:8.4f} {emin:16.6f} {de:10.4f}\n')
np.savez('dimer.npz', d=np.array(ds), E=E, E_atom=Ea, coef=c, r_e=re_, E_min=emin, D_e=de)
print(open('dimer.txt').read(), end='')
try:
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    x = np.linspace(min(ds)-0.02, max(ds)+0.02, 100)
    fig, ax = plt.subplots(figsize=(4.5, 3.5))
    ax.plot(ds, E-2*Ea, 'ko', label=f'{mol}: E(d) - 2 E({atom})'); ax.plot(x, np.polyval(c, x)-2*Ea, 'k-', lw=1, label='parabola')
    ax.axvline(re_, color='gray', lw=0.7, ls='--'); ax.set_xlabel('d (A)'); ax.set_ylabel('E - 2 E_atom (eV)')
    ax.set_title(f'r_e = {re_:.3f} A, D_e = {de:.2f} eV', fontsize=10); ax.legend(fontsize=8); fig.tight_layout(); fig.savefig('dimer.png', dpi=130)
except Exception as e: print('no figure:', e)
