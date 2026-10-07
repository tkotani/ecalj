#!/usr/bin/env python3
"""Figures and data of ecaljdoc manual/mtradius.md (2026-10-07). LDA bands of Rb2Sc2O4 (mp-7650, P6_3/mmc) for Rb MT radii 2.3-2.8 a.u.
(the scan of 2026-10-07, ecalj b81da2342 built on kt1, made on kt1 in /mnt/data1/rscan and on t14 in ~/work/mlo_vac/Rr<r>, L_mp-7650,
newcap/mp-7650: the band path of getsyml with at least 5 points per segment), and the levels on A-L without a degenerate partner.
   make_figs.py <dir of Rr2.3 ... Rr2.7, L_mp-7650 (Rb 2.8), newcap/mp-7650>     -> fig1.png, fig2.png, bands.npz"""
import sys, glob, os, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
top = sys.argv[1]; here = os.path.dirname(os.path.abspath(__file__))
RUNS = [('2.3', 'Rr2.3'), ('2.4', 'Rr2.4'), ('2.5', 'Rr2.5'), ('2.6', 'Rr2.6'), ('2.7', 'Rr2.7'), ('2.8', 'L_mp-7650'), ('2.4 (rsmh 1.2)', 'newcap/mp-7650')]

def load(d):
    X, E, S = [], [], []
    for f in sorted(glob.glob(f'{d}/bnd0*.spin1')):
        for l in open(f):
            if l.startswith('#'): continue
            w = l.split()
            if len(w) >= 3: X.append(float(w[1])); E.append(float(w[2])); S.append(int(os.path.basename(f)[3:6]))
    return np.array(X), np.array(E), np.array(S)

def unpaired(X, E, S, seg=5, lo=3, hi=25, tol=0.02):
    """levels on segment seg (A-L) with no partner within tol eV (P6_3/mmc: all doubly degenerate there)"""
    out = []
    for x in np.unique(X[S == seg]):
        e = np.sort(E[(S == seg) & (X == x)]); e = e[(e > lo) & (e < hi)]; i = 0
        while i < len(e):
            if i + 1 < len(e) and e[i + 1] - e[i] < tol: i += 2
            else: out.append((x, e[i])); i += 1
    return np.array(out).reshape(-1, 2)

def ticks(d):
    segs = [l.split() for l in open(glob.glob(f'{d}/syml.*')[0]) if l.split() and l.split()[0] != '0']
    xs, labs = [], []
    for f, s in zip(sorted(glob.glob(f'{d}/bnd0*.spin1')), segs):
        x = sorted(float(l.split()[1]) for l in open(f) if not l.startswith('#') and len(l.split()) >= 3)
        for v, lab in ((x[0], s[7]), (x[-1], s[8])):
            lab = 'Γ' if lab == 'GAMMA' else lab
            if xs and abs(xs[-1] - v) < 1e-6:
                if labs[-1] != lab: labs[-1] += '|' + lab
            else: xs.append(v); labs.append(lab)
    return xs, labs

data = {}
for r, d in RUNS:
    X, E, S = load(f'{top}/{d}'); U = unpaired(X, E, S)
    data[r] = (X, E, U)
np.savez_compressed(os.path.join(here, 'bands.npz'), **{f'x_{r}': v[0] for r, v in data.items()}, **{f'e_{r}': v[1] for r, v in data.items()},
                    **{f'unpaired_{r}': v[2] for r, v in data.items()},
                    note=np.array('LDA bands of Rb2Sc2O4 mp-7650 (eV from the VBM) along the getsyml path, per Rb MT radius (a.u.); '
                                  'unpaired_*: (x, E) of levels on A-L without a partner within 0.02 eV. Made by make_figs.py on 2026-10-07 from '
                                  + top))
# Fig. 1: Rb 2.8 (production before 2026-10-07) and 2.4 with the present rule (rsmh = r/2)
xs, labs = ticks(f'{top}/L_mp-7650')
fig, axs = plt.subplots(1, 2, figsize=(14, 5.6), sharey=True)
for ax, r, t in ((axs[0], '2.8', 'Rb r = 2.8 a.u. (before 2026-10-07)'), (axs[1], '2.4 (rsmh 1.2)', 'Rb r = 2.4 a.u., rsmh = 1.2 (from 2026-10-07)')):
    X, E, U = data[r]
    ax.plot(X, E, 'o', color='0.25', ms=1.8)
    if len(U): ax.plot(U[:, 0], U[:, 1], 'o', mfc='none', mec='darkorange', ms=11, mew=2, label=f'no degenerate partner on A-L ({len(U)} up to 25 eV, {int((U[:, 1] < 10).sum())} in this window)'); ax.legend(loc='upper right', fontsize=9)
    for v in xs: ax.axvline(v, color='k', lw=0.5)
    ax.set_xticks(xs); ax.set_xticklabels(labs, fontsize=9); ax.set_xlim(xs[0], xs[-1]); ax.set_ylim(-4, 10); ax.axhline(0, lw=0.3, color='k')
    ax.set_title(f'Rb2Sc2O4 mp-7650, LDA, {t}', fontsize=10)
axs[0].set_ylabel('E (eV) from the VBM'); plt.tight_layout(); plt.savefig(os.path.join(here, 'fig1.png'), dpi=100)
# Fig. 2: radii 2.3 ... 2.8 (rsmh 1.4), -4 ... 25 eV
fig, axs = plt.subplots(2, 3, figsize=(18, 12), sharey=True)
for ax, r in zip(axs.flat, ('2.3', '2.4', '2.5', '2.6', '2.7', '2.8')):
    X, E, U = data[r]
    ax.plot(X, E, 'o', color='0.25', ms=1.3)
    if len(U): ax.plot(U[:, 0], U[:, 1], 'o', mfc='none', mec='darkorange', ms=8, mew=1.5)
    for v in xs: ax.axvline(v, color='k', lw=0.5)
    ax.set_xticks(xs); ax.set_xticklabels(labs, fontsize=8); ax.set_xlim(xs[0], xs[-1]); ax.set_ylim(-4, 25); ax.axhline(0, lw=0.3, color='k')
    ax.set_title(f'Rb r = {r} a.u. (rsmh 1.4): unpaired on A-L {len(U)}, lowest {U[:, 1].min():.1f} eV' if len(U) else f'Rb r = {r} a.u.: none', fontsize=10)
axs[0, 0].set_ylabel('E (eV) from the VBM'); axs[1, 0].set_ylabel('E (eV) from the VBM'); plt.tight_layout(); plt.savefig(os.path.join(here, 'fig2.png'), dpi=90)
for r, v in data.items(): print(r, 'unpaired', len(v[2]), 'lowest', v[2][:, 1].min() if len(v[2]) else '-')
