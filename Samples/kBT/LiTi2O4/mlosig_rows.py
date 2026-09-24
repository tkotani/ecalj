#!/usr/bin/env python3
"""Sigma interpolation: conventional (Sigma^MTO, full PMT) vs MLO (154 orbitals).
LiTi2O4 6^3, Gamma-X 211 points.  usage: mlosig_rows.py DIR OUT.png"""
import sys, numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
RY = 13.605693122994
D, OUT = sys.argv[1], sys.argv[2]
def rd_bnd(f):
    xs=[];es=[];cx=[];ce=[];ef=None
    for ln in open(f):
        if ln.startswith('#'):
            if 'base=eferm' in ln: ef=float(ln.split()[1])
            continue
        t=ln.split()
        if len(t)<3:
            if cx: xs.append(cx);es.append(ce);cx=[];ce=[]
            continue
        cx.append(float(t[1])); ce.append(float(t[2]))
    if cx: xs.append(cx);es.append(ce)
    return np.array(xs),np.array(es),ef
def rd_mlo(f):
    d={}
    for ln in open(f):
        t=ln.split()
        if len(t)==4: d.setdefault(int(t[3]),[]).append((float(t[0]),float(t[1])))
    nb=max(d)
    return np.array([p[0] for p in d[1]]), np.array([[p[1] for p in d[ib]] for ib in range(1,nb+1)])
xq,eq,efq = rd_bnd(f'{D}/bnd_qsgw.dat')
Xm,Es     = rd_mlo(f'{D}/mlo_withSig.dat')
x = xq[0]; EM = Es*RY - efq                      # MLO, eV rel. E_F(QSGW)
MESH = [0.0, 1/3, 2/3, 1.0]                      # 6^3 mesh points on Gamma-X
def resid(y, w=5):
    r = np.full(len(y), np.nan)
    for i in range(len(y)):
        a, b = max(0, i-w), min(len(y), i+w+1)
        if b-a >= 5: r[i] = y[i] - np.polyval(np.polyfit(x[a:b], y[a:b], 2), x[i])
    return r
fig, ax = plt.subplots(2, 2, figsize=(11.5, 7.2),
                       gridspec_kw=dict(height_ratios=[2, 1]))
sets = [('conventional: $\\Sigma^{\\rm MTO}$ interpolated (230 ch, full PMT)', eq[32:44], 'tab:blue'),
        ('new: $\\Sigma^{\\rm MLO}$ interpolated (154 orb)',                   EM[32:44], 'tab:red')]
for col, (title, E, c) in enumerate(sets):
    for b in E: ax[0][col].plot(x, b, '-', lw=1.0, color=c)
    ax[0][col].set_title(title, fontsize=10)
    ax[0][col].set_ylabel('$E-E_F$ [eV]'); ax[0][col].set_ylim(-0.7, 1.3)
    w = np.array([resid(b) for b in E])*1e3
    for r in w: ax[1][col].plot(x, r, '-', lw=0.9, color=c)
    ax[1][col].set_ylabel('residual [meV]'); ax[1][col].set_ylim(-7, 7)
    ax[1][col].set_xlabel('$\\Gamma \\to X$')
    mx = np.nanmax(np.abs(w))
    ax[1][col].text(0.02, 0.06, f'max |resid| = {mx:.2f} meV', transform=ax[1][col].transAxes, fontsize=9)
    for a in (ax[0][col], ax[1][col]):
        for m in MESH: a.axvline(m, color='k', ls=':', lw=0.7)
        a.set_xlim(0, 1)
fig.suptitle('LiTi$_2$O$_4$ 6$^3$ (from n666_nk6_from_lda iter 10) — t$_{2g}$ bands 33–44 along $\\Gamma$–X, '
             '211 pts; dotted = $\\Sigma$ mesh points', fontsize=10)
plt.tight_layout(rect=[0,0,1,0.96]); plt.savefig(OUT, dpi=140)
print('wrote', OUT)
