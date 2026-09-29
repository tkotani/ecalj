#!/usr/bin/env python3
"""How much the lowest t2g MLO band of LiTi2O4 wiggles between the Sigma mesh points on Gamma-X (README.md table 3).
usage: band_wiggle.py [iterations...]      (default 10 15 20 25 30 35 40)
Reads liti2o4_k6tf32_10_40_data.npz and liti2o4_k9tf32_10_40_data.npz (written by update_six_rows.sh / mlo_rows.py).

  D(x) = E_N(x) - E_LDA(x)     the QSGW correction of the lowest band (energy-sorted at each x) after iteration N, x on Gamma-X
  fit    D(x) ~ sum_{n=0..3} c_n cos(n pi x)     (even about Gamma and X, as a band along Gamma-X is)
  wiggle = max |D - fit| and its rms; a smooth correction gives a small residual
  dip    = E(Gamma) - min_{x<0.3} E(x): how far the band goes below its value at Gamma
A rough indicator only: the fit order is a choice, and D has real structure where the two lowest bands approach each other."""
import sys, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
its = [int(a) for a in sys.argv[1:]] or [10, 15, 20, 25, 30, 35, 40]
def wiggle(x, d, n=3):
    A = np.stack([np.cos(k * np.pi * x) for k in range(n + 1)], 1)
    c = np.linalg.lstsq(A, d, rcond=None)[0]
    r = d - A @ c
    return 1e3 * abs(r).max(), 1e3 * np.sqrt((r * r).mean())
print('# mesh iter  wiggle_max[meV] wiggle_rms[meV] dip[meV]')
for m in (6, 9):
    d = np.load(f'{HERE}/liti2o4_k{m}tf32_10_40_data.npz')
    x, lda = d[f'k{m}tf32_lda__x'], d[f'k{m}tf32_lda__E'][0]
    for it in its:
        e = d[f'k{m}tf32_iter{it}__E'][0]
        w = wiggle(x, e - lda)
        print(f'{m}^3 {it:4d}  {w[0]:6.1f} {w[1]:6.1f} {1e3 * (e[0] - e[x < 0.3].min()):6.1f}')
