#!/usr/bin/env python3
"""Compare x0homo.dat of hhomogas with the analytic Lindhard function (README.md, eqs. (1)-(4)).

  python3 plot_lindhard.py [label=file ...]        default: 10x10x10=x0homo.dat

Writes lindhard.png and lindhard.npz (the plotted numbers, the source files and the date).
x0homo.dat: one block for each q (iq = 1..10); columns
  iq iw qx qy qz(2pi/alat) |q|(1/bohr) omega(Ry) Re Im of chi0*Omega (1/Hartree; Omega = cell volume)
The Fermi energy and the q points are fixed in SRC/subroutines/main_hhomogas.f90.
"""
import datetime
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ALAT = 6.35                      # bohr, [struc] alat of ctrlg.es.toml (bcc: Omega = alat^3/2)
OMEGA = ALAT**3 / 2
EF_RY = 0.237356                 # Ry, "efermi_egas::" in the output of hhomogas
KF = np.sqrt(EF_RY)              # 1/bohr (Rydberg units: E = k^2)
NF = KF / np.pi**2               # states/(Hartree bohr^3), both spins: eq. (2)
IQ = 10                          # the block that is plotted: q = 0.01 (2pi/alat) along x


def blocks(fname):
    out, cur = [], []
    for line in open(fname):
        if line.strip():
            cur.append([float(x) for x in line.split()])
        elif cur:
            out.append(np.array(cur))
            cur = []
    if cur:
        out.append(np.array(cur))
    return out


def lindhard(u, z):
    """chi0/N_F at u = omega/(q v_F), z = q/(2 k_F): eqs. (3) and (4)."""
    a, b = z - u, z + u
    with np.errstate(divide='ignore', invalid='ignore'):
        re = -(0.5 + (1 - a * a) / (8 * z) * np.log(abs((a + 1) / (a - 1)))
               + (1 - b * b) / (8 * z) * np.log(abs((b + 1) / (b - 1))))
    im = np.where(u < 1 - z, -np.pi / 2 * u, np.where(u < 1 + z, -np.pi / (8 * z) * (1 - a * a), 0.0))
    return re, im


args = sys.argv[1:] or ['10x10x10=x0homo.dat']
now = datetime.datetime.now().strftime('%Y-%m-%d %H:%M')
fig, axes = plt.subplots(1, 2, figsize=(9, 3.8), sharex=True)
markers = ['o', 'x', '+', 's']
save = {'date': now, 'kf_per_bohr': KF, 'nf': NF, 'omega_bohr3': OMEGA, 'iq': IQ}
for n, arg in enumerate(args):
    label, fname = arg.split('=', 1)
    b = blocks(fname)[IQ - 1]
    q = b[0, 5]
    z = q / (2 * KF)
    u = b[:, 6] / 2 / (q * KF)               # omega(Hartree)/(q v_F), v_F = k_F in Hartree units
    sel = u <= 2.0
    u = u[sel]
    re = b[sel, 7] / (NF * OMEGA)
    im = b[sel, 8] / (NF * OMEGA)
    for ax, y in zip(axes, (re, im)):
        ax.plot(u, y, markers[n % 4], ms=3.5, mfc='none', mec='k', mew=0.6, label='tetrahedron ' + label)
    save[f'u_{n}'], save[f're_{n}'], save[f'im_{n}'] = u, re, im
    save[f'label_{n}'], save[f'source_{n}'] = label, os.path.abspath(fname)
    i0 = np.argmin(abs(u))
    print(f'{label}: q = {q:.5f} /bohr, q/2kF = {z:.5f}, chi0/N_F(omega=0) = {re[i0]:.4f} (exact {lindhard(np.array([0.0]), z)[0][0]:.4f})')
ue = np.linspace(0, 2, 2001)
ree, ime = lindhard(ue, z)
for ax, y, name in zip(axes, (ree, ime), ('Re', 'Im')):
    ax.plot(ue, y, 'k-', lw=0.6, label='Lindhard function')
    ax.axhline(0, color='0.6', lw=0.4)
    ax.set_xlabel(r'$\omega/(q v_{\rm F})$')
    ax.set_ylabel(name + r' $\chi_0/N_{\rm F}$')
    ax.set_xlim(0, 2)
    ax.legend(frameon=False, fontsize=8)
axes[0].set_ylim(-1.5, 2.5)
save['u_exact'], save['re_exact'], save['im_exact'] = ue, ree, ime
fig.suptitle(f'Electron gas, q = {q:.5f}/bohr (q/2k_F = {z:.4f}), k_F = {KF:.5f}/bohr   [{now}]', fontsize=9)
fig.tight_layout()
fig.savefig('lindhard.png', dpi=150)
np.savez('lindhard.npz', **save)
print('lindhard.png and lindhard.npz written')
