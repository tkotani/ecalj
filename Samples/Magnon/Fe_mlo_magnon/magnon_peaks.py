#!/usr/bin/env python3
"""Positions of the magnon peaks (the maximum of |Im R(q,omega)| at each q, 0 < omega < 1.6 eV) of one or more
MagSuscep.syml001 (job_mlo_magnon, Im R in column 9) or TrRpm.syml001 (job_magnon, Im Tr R in column 7) files.
usage: magnon_peaks.py <label>=<file> [<label>=<file> ...]     -> magnon_peaks.png, magnon_peaks.npz, a table on stdout
(2026-09-30: the figure of README.md)"""
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
        a = np.array([(x, y) for x, y in v if 1e-6 < x < 1.6]); i = np.argmax(np.abs(a[:, 1]))
        q.append(np.linalg.norm(k)); w.append(a[i, 0])
    return np.array(q), np.array(w)
data = {}
for arg in sys.argv[1:]:
    lab, f = arg.split('=', 1); data[lab] = (peaks(f), os.path.abspath(f))
fig, ax = plt.subplots(figsize=(5, 4))
for lab, ((q, w), f) in data.items(): ax.plot(q, w, 'o-', ms=4, lw=1, label=lab)
ax.set_xlabel('|q| (2pi/a), Gamma to H'); ax.set_ylabel('peak of Im R (eV)'); ax.legend(fontsize=8); ax.grid(alpha=.3)
ax.set_title('bcc Fe: magnon energy along Gamma-H', fontsize=10); fig.tight_layout(); fig.savefig('magnon_peaks.png', dpi=130)
np.savez('magnon_peaks.npz', **{f'{lab}_q': v[0][0] for lab, v in data.items()}, **{f'{lab}_w': v[0][1] for lab, v in data.items()},
         files=np.array([f'{lab}: {v[1]}' for lab, v in data.items()]))
print('| q |', ' | '.join(data), '|'); print('|---|', '|'.join('---' for _ in data), '|')
q0 = next(iter(data.values()))[0][0]
for i, q in enumerate(q0): print(f'| {q:.1f} |', ' | '.join(f'{v[0][1][i]:.3f}' for v in data.values()), '|')
