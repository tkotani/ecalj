#!/usr/bin/env python3
"""LiTi2O4 6^3 MLO-QSGW after gwsc 10 from LDA: the tf32 run on top of the fp32 run (Gamma-X), 2026-09-27.
  plot_tf32_fp32.py <dir holding the run dirs qmlo_k6_gwsc10 (fp32) and qmlo_k6_tf32n (tf32) of run_gwsc10.sh>
Left: MLO band, middle: conventional sigm band (job_band), right: tf32 - fp32 per band index (meV), |E-EF|<3 eV.
The MLO band file has 1e-5 Ry steps (0.136 meV), hence the staircase."""
import sys, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
d = sys.argv[1]
def bands(f):
    blocks, cur = [], []
    for l in open(f):
        if l.startswith('#'): continue
        w = l.split()
        if not w:
            if cur: blocks.append(cur); cur = []
            continue
        cur.append((float(w[1]), float(w[2])))
    if cur: blocks.append(cur)
    nk = min(len(b) for b in blocks)
    return np.array([p[0] for p in blocks[0][:nk]]), np.array([[p[1] for p in b[:nk]] for b in blocks])
fig, ax = plt.subplots(1, 3, figsize=(14, 5.4), gridspec_kw={'width_ratios': [1, 1, 1.1]})
win = 3.0
for a, kind, title in [(ax[0], 'bnd_mlo_final.dat', 'MLO band'), (ax[1], 'bndPMT_final.dat', 'sigm band (job_band)')]:
    x, e32 = bands(f'{d}/qmlo_k6_gwsc10/{kind}')
    _, etf = bands(f'{d}/qmlo_k6_tf32n/{kind}')
    nb, nk = min(len(e32), len(etf)), min(e32.shape[1], etf.shape[1])
    e32, etf, x = e32[:nb, :nk], etf[:nb, :nk], x[:nk]
    for b in range(nb):
        if np.all((e32[b] < -9) | (e32[b] > 7)): continue
        a.plot(x, e32[b], color='0.6', lw=2.4, label='fp32' if b == 0 or 'fp32' not in [t.get_label() for t in a.lines] else None)
        a.plot(x, etf[b], color='C3', lw=0.8, label='tf32' if 'tf32' not in [t.get_label() for t in a.lines] else None)
    a.axhline(0, color='k', lw=0.5); a.set_xlim(x[0], x[-1]); a.set_ylim(-9, 7)
    a.set_xticks([x[0], x[-1]]); a.set_xticklabels([r'$\Gamma$', 'X']); a.set_ylabel(r'$E-E_F$ (eV)')
    m = np.abs(e32) < win; dd = (etf - e32) * 1000
    a.set_title(f'{title}: rms {np.sqrt(np.mean(dd[m]**2)):.1f}, max {np.max(np.abs(dd[m])):.1f} meV (|E|<{win:g} eV)', fontsize=9)
    a.legend(loc='lower right', fontsize=8)
    for b in range(nb):
        mb = np.abs(e32[b]) < win
        if mb.any(): ax[2].plot(x[mb], dd[b][mb], lw=0.7, color='C0' if kind.startswith('bnd_mlo') else 'C1',
                                label=(title if b == np.argmax(np.any(np.abs(e32) < win, axis=1)) else None))
ax[2].axhline(0, color='k', lw=0.5); ax[2].set_xlim(x[0], x[-1]); ax[2].set_xticks([x[0], x[-1]])
ax[2].set_xticklabels([r'$\Gamma$', 'X']); ax[2].set_ylabel('tf32 - fp32 (meV)'); ax[2].legend(fontsize=8)
ax[2].set_title('bands within $\\pm$3 eV (spikes at crossings = index swaps)', fontsize=9)
fig.suptitle(r'LiTi$_2$O$_4$ $6^3$ MLO-QSGW, gwsc 10 from LDA (kt1, 2026-09-27): tf32 (red) on fp32 (grey)', fontsize=11)
fig.tight_layout(); fig.savefig(f'liti2o4_k6_tf32_vs_fp32_bands.png', dpi=130)
print('ok')
