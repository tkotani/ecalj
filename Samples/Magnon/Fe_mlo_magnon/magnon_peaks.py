#!/usr/bin/env python3
"""Positions of the magnon peaks (the maximum of |Im R(q,omega)| at each q, 0 < omega < 1.6 eV) of one or more
MagSuscep.syml001 (job_mlo_magnon, Im R in column 9) or TrRpm.syml001 (job_magnon, Im Tr R in column 7) files.
usage: magnon_peaks.py [--title T] [--xlabel X] <label>=<file> [<label>=<file> ...]  -> magnon_peaks.png, magnon_peaks.npz, a table
(2026-09-30: the figure of README.md)
The peak at each q is the largest |Im R|, except near Gamma (the peak of the previous q below 0.05 eV), where it is the
lowest-frequency local maximum that reaches 30 % of the largest: the acoustic magnon. (2026-10-02 09:09: with two magnetic atoms
(FeCo) Tr R holds the optical branch too; at small q the acoustic peak falls between the points of the frequency mesh and
the optical one near 1.2 eV was taken. Away from Gamma the largest is kept: in the Stoner continuum (Fe, high q) the
lowest maximum is a weak side bump.)"""
import sys, os, collections, numpy as np
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
def peaks(f):
    col = 8 if os.path.basename(f).startswith('MagSuscep') else 6
    d = collections.OrderedDict()
    for l in open(f):
        t = l.split()
        if len(t) < 9 or l.lstrip().startswith('#'): continue
        d.setdefault(tuple(round(float(x), 3) for x in t[:3]), []).append((float(t[4]), float(t[col])))
    q, w = [], []
    for k, v in d.items():
        a = np.array([(x, y) for x, y in v if 1e-6 < x < 1.6]); y = np.abs(a[:, 1])
        i = int(np.argmax(y))
        if not w or w[-1] < 0.05:      # near Gamma: the acoustic magnon, the lowest strong local maximum
            loc = [j for j in range(len(y)) if (j == 0 or y[j] >= y[j - 1]) and (j == len(y) - 1 or y[j] >= y[j + 1]) and y[j] >= 0.3 * y.max()]
            if loc: i = loc[0]
        q.append(np.linalg.norm(k)); w.append(a[i, 0])
    return np.array(q), np.array(w)
title, xlabel, args = 'bcc Fe: magnon energy along Gamma-H', '|q| (2pi/a), Gamma to H', sys.argv[1:]
while args and args[0] in ('--title', '--xlabel'):
    if args[0] == '--title': title = args[1]
    else: xlabel = args[1]
    args = args[2:]
data = {}
for arg in args:
    lab, f = arg.split('=', 1); data[lab] = (peaks(f), os.path.abspath(f))
fig, ax = plt.subplots(figsize=(5, 4))
for lab, ((q, w), f) in data.items(): ax.plot(q, w, 'o-', ms=4, lw=1, label=lab)
ax.set_xlabel(xlabel); ax.set_ylabel('peak of Im R (eV)'); ax.legend(fontsize=8); ax.grid(alpha=.3)
ax.set_title(title, fontsize=10); fig.tight_layout(); fig.savefig('magnon_peaks.png', dpi=130)
np.savez('magnon_peaks.npz', **{f'{lab}_q': v[0][0] for lab, v in data.items()}, **{f'{lab}_w': v[0][1] for lab, v in data.items()},
         files=np.array([f'{lab}: {v[1]}' for lab, v in data.items()]))
print('| q |', ' | '.join(data), '|'); print('|---|', '|'.join('---' for _ in data), '|')
q0 = next(iter(data.values()))[0][0]
for i, q in enumerate(q0): print(f'| {q:.1f} |', ' | '.join(f'{v[0][1][i]:.3f}' for v in data.values()), '|')
