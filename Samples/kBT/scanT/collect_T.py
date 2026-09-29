#!/usr/bin/env python3
"""Collect a temperature scan made by scan_T.sh: a table, the numbers (npz) and two figures.

usage: collect_T.py <work dir> <sname> <out prefix> [--emin=-14] [--emax=8] [--title='Si 4x4x4']

Reads <work dir>/scan_T.log (the labels and the two keys) and, in each <work dir>/<label>/:
  llmf.<N>run   the SCF after QSGW iteration N: the last 'gap' line (insulators) -> gap per iteration
  EFERMI        E_F of the T=0 tetrahedron (heftet);  EFERMI_kbt  E_F at t_tetrakbt > 0 (ecaljdoc kBT, equation (6))
  bnd*.spin*    the bands of the last iteration along syml.<sname> (job_band): energy in eV from the zero job_band takes
                (the valence band maximum for an insulator, E_F for a metal); both spins when there are two
  save.<sname>  the magnetic moment of the last SCF (mmom=), when spin polarized
Writes <out prefix>_table.md, <out prefix>_data.npz, <out prefix>_bands.png (one panel per setting, side by side) and,
when there is a gap, <out prefix>_gap.png (gap against T).  The energies at the ends of the path segments are in the table
for the bands nearest to the zero."""
import sys, os, re, glob, datetime
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

args = [a for a in sys.argv[1:] if not a.startswith('--')]
opt = dict(a[2:].split('=', 1) for a in sys.argv[1:] if a.startswith('--'))
W, S, OUT = args[0], args[1], args[2]
EMIN, EMAX = float(opt.get('emin', -14)), float(opt.get('emax', 8))
RY = 13.605693

def settings():
    out = []
    for ln in open(f'{W}/scan_T.log'):
        m = re.match(r'(\S+) t_tetrakbt=(\S+) t_sigmaw=(\S+) gwsc rc=(\d+) job_band rc=(\d+) (\d+) s', ln)
        if m: out.append(dict(label=m[1], tt=float(m[2]), ts=float(m[3]), rc=int(m[4]), rcb=int(m[5]), secs=int(m[6])))
    return out

def gaps(d):
    """gap (eV), VBM and CBM (Ry) after each iteration; empty for a metal"""
    g = []
    for f in sorted(glob.glob(f'{d}/llmf.*run'), key=lambda f: int(re.search(r'llmf\.(\d+)run', f)[1])):
        ls = [l for l in open(f, errors='replace') if 'gap =' in l and 'VBmax' in l]
        if ls:
            t = ls[-1].replace('=', ' ').split()
            g.append((float(t[7]), float(t[1]), float(t[3])))     # gap eV, VBmax Ry, CBmin Ry
    return np.array(g)

def first_number(f):
    return float(open(f).readline().split()[0].replace('D', 'E')) if os.path.exists(f) else np.nan

def bands(d, isp=1):
    """x (path coordinate), E[band, point] (eV), the x of the segment ends and their labels"""
    segs = sorted(glob.glob(f'{d}/bnd[0-9][0-9][0-9].spin{isp}'))
    xs, es, ends = [], [], []
    for f in segs:
        a = np.array([[float(v) for v in l.split()[:3]] for l in open(f) if l.strip() and not l.startswith('#')])
        nb = int(a[:, 0].max()); npnt = len(a) // nb
        x = a[:npnt, 1]; e = a[:, 2].reshape(nb, npnt)
        xs.append(x); es.append(e); ends.append((x[0], x[-1]))
    nb = min(e.shape[0] for e in es)
    lab = [l.split()[-2:] for l in open(f'{d}/syml.{S}') if l.strip() and not l.startswith('0') and len(l.split()) >= 9]
    return np.concatenate(xs), np.concatenate([e[:nb] for e in es], axis=1), ends, lab

def mmom(d):
    ls = [l for l in open(f'{d}/save.{S}', errors='replace') if 'mmom=' in l] if os.path.exists(f'{d}/save.{S}') else []
    return float(ls[-1].split('mmom=')[1].split()[0]) if ls else np.nan

sets = settings()
rows, data = [], {}
NSP = 2 if glob.glob(f"{W}/{sets[0]['label']}/bnd001.spin2") else 1
for s in sets:
    d = f"{W}/{s['label']}"
    g = gaps(d)
    x, e, ends, lab = bands(d)
    s.update(gap=g, x=x, e=e, ends=ends, lab=lab, ef=first_number(f'{d}/EFERMI'), efk=first_number(f'{d}/EFERMI_kbt') if s['tt'] > 0 else np.nan,
             mmom=mmom(d), e2=bands(d, 2)[1] if NSP == 2 else None)
    k = s['label']
    data.update({f'{k}__x': x, f'{k}__E': e, f'{k}__gap_per_iteration': g, f'{k}__t_tetrakbt': s['tt'], f'{k}__t_sigmaw': s['ts'],
                 f'{k}__EFERMI': s['ef'], f'{k}__EFERMI_kbt': s['efk'], f'{k}__dir': os.path.abspath(d), f'{k}__mmom': s['mmom']})
    if NSP == 2: data[f'{k}__E_spin2'] = s['e2']
np.savez_compressed(f'{OUT}_data.npz', **data)

insul = all(len(s['gap']) for s in sets)
ref = sets[0]
with open(f'{OUT}_table.md', 'w') as o:
    o.write(f'<!-- written by collect_T.py {datetime.datetime.now():%Y-%m-%d %H:%M} from {os.path.abspath(W)} -->\n')
    if insul:
        o.write('| label | t_tetrakbt (K) | t_sigmaw (K) | iterations | gap (eV) | change from the first row (meV) | last change of the gap (meV) | E_F(T) − E_F(0) (eV) | time (s) |\n')
        o.write('| --- | --- | --- | --- | --- | --- | --- | --- | --- |\n')
        for s in sets:
            g = s['gap'][:, 0]
            o.write(f"| {s['label']} | {s['tt']:g} | {s['ts']:g} | {len(g)} | {g[-1]:.3f} | {1e3 * (g[-1] - ref['gap'][-1, 0]):+.0f} | "
                    f"{1e3 * (g[-1] - g[-2]) if len(g) > 1 else float('nan'):+.1f} | "
                    f"{(s['efk'] - s['ef']) * RY if s['tt'] > 0 else 0:+.3f} | {s['secs']} |\n")
    else:
        o.write('| label | t_tetrakbt (K) | t_sigmaw (K) | E_F(T) − E_F(0) (eV) | magnetic moment (mu_B) | time (s) |\n| --- | --- | --- | --- | --- | --- |\n')
        for s in sets:
            o.write(f"| {s['label']} | {s['tt']:g} | {s['ts']:g} | {(s['efk'] - s['ef']) * RY if s['tt'] > 0 else 0:+.3f} | {s['mmom']:.3f} | {s['secs']} |\n")
    # energies at the segment ends, for the bands within [emin, emax] of the zero at that point
    o.write('\n| point | spin | band | ' + ' | '.join(s['label'] for s in sets) + ' |\n| --- | --- | --- | ' + ' | '.join('---' for s in sets) + ' |\n')
    for isp, key in ((1, 'e'), (2, 'e2'))[:NSP]:
        seen = set()
        for (x0, x1), (l0, l1) in zip(ref['ends'], ref['lab']):
            for xx, ll in ((x0, l0), (x1, l1)):
                if ll in seen: continue
                seen.add(ll)
                i = int(np.argmin(abs(ref['x'] - xx)))
                for b in range(ref[key].shape[0]):
                    if float(opt.get('tmin', -6)) <= ref[key][b, i] <= float(opt.get('tmax', 6)):
                        o.write(f'| {ll} | {isp} | {b + 1} | ' + ' | '.join(f"{s[key][b, int(np.argmin(abs(s['x'] - xx)))]:.3f}" for s in sets) + ' |\n')

n = len(sets)
fig, axs = plt.subplots(NSP, n, figsize=(3.2 * n, 5.2 * NSP), sharey=True, squeeze=False)
for isp in range(1, NSP): # the second spin in the second row, same layout
    for a, s in zip(axs[isp], sets):
        for b in s['e2']: a.plot(s['x'], b, '-', color='k', lw=0.6)
        for (x0, x1) in s['ends']: a.axvline(x1, color='0.6', lw=0.5)
        a.axhline(0, color='0.5', lw=0.5, ls='--'); a.set_xlim(s['x'].min(), s['x'].max()); a.set_ylim(EMIN, EMAX)
        a.set_xticks([s['ends'][0][0]] + [e[1] for e in s['ends']]); a.set_xticklabels([s['lab'][0][0]] + [l[1] for l in s['lab']], fontsize=7)
        a.set_title(f"{s['label']}  spin 2", fontsize=8)
    axs[isp][0].set_ylabel('E (eV) from E_F, spin 2')
ax = axs
for a, s in zip(ax[0], sets):
    for b in s['e']:
        a.plot(s['x'], b, '-', color='k', lw=0.6)
    for (x0, x1) in s['ends']: a.axvline(x1, color='0.6', lw=0.5)
    a.axhline(0, color='0.5', lw=0.5, ls='--')
    a.set_xlim(s['x'].min(), s['x'].max()); a.set_ylim(EMIN, EMAX)
    a.set_xticks([s['ends'][0][0]] + [e[1] for e in s['ends']]); a.set_xticklabels([s['lab'][0][0]] + [l[1] for l in s['lab']], fontsize=7)
    tl = f"{s['label']}\nt_tetrakbt {s['tt']:g} K, t_sigmaw {s['ts']:g} K"
    if insul: tl += f"\ngap {s['gap'][-1, 0]:.3f} eV"
    if NSP == 2: tl += f"\nspin 1, moment {s['mmom']:.3f}"
    a.set_title(tl, fontsize=8)
ax[0][0].set_ylabel('E (eV) from ' + ('the valence band maximum' if insul else 'E_F') + (', spin 1' if NSP == 2 else ''))
fig.suptitle(f"{opt.get('title', S)}: QSGW bands after {len(glob.glob(W + '/' + ref['label'] + '/llmf.*run'))} iterations   "
             f"(collect_T.py {datetime.datetime.now():%Y-%m-%d %H:%M})", fontsize=9)
plt.tight_layout(rect=[0, 0, 1, 0.95]); plt.savefig(f'{OUT}_bands.png', dpi=110); plt.close()

if insul:
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
    T = [s['ts'] if s['tt'] == 0 else abs(s['tt']) for s in sets]
    for s, t in zip(sets, T):
        ax[1].plot(range(1, len(s['gap']) + 1), s['gap'][:, 0], 'o-', ms=3, lw=0.7, label=s['label'])
    ax[0].plot([t for s, t in zip(sets, T) if s['tt'] > 0], [s['gap'][-1, 0] for s in sets if s['tt'] > 0], 'ko-', ms=4, lw=0.7, label='t_tetrakbt = t_sigmaw = T')
    for s, t in zip(sets, T):
        if s['tt'] <= 0: ax[0].plot([t], [s['gap'][-1, 0]], 's', ms=5, mfc='none', label=f"{s['label']} (t_tetrakbt {s['tt']:g})")
    ax[0].set_xlabel('T (K)'); ax[0].set_ylabel('gap (eV)'); ax[0].legend(fontsize=7); ax[0].grid(lw=0.3)
    ax[1].set_xlabel('QSGW iteration'); ax[1].set_ylabel('gap (eV)'); ax[1].legend(fontsize=7); ax[1].grid(lw=0.3)
    fig.suptitle(f"{opt.get('title', S)}: band gap against the electron temperature   (collect_T.py {datetime.datetime.now():%Y-%m-%d %H:%M})", fontsize=9)
    plt.tight_layout(rect=[0, 0, 1, 0.94]); plt.savefig(f'{OUT}_gap.png', dpi=110); plt.close()
print(open(f'{OUT}_table.md').read())
