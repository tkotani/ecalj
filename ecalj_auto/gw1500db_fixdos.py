#!/usr/bin/env python3
"""Correct the DOS energies of the npz files of gw1500db_extract.py made before 2026-10-07 19:14 (in place).

    gw1500db_fixdos.py <npz dir> ...

dos.tot.<mpid>.chk holds E - E_F (Ry), but the extraction subtracted E_F once more: dosE was (E - E_F - ef_ry) * Ry.
dosE += ef_ry * Ry undoes it; dos_ref='E_F' marks a corrected (or newly extracted) file, so running this twice does nothing.
"""
import sys, glob, os
import numpy as np
RY = 13.605693122994
n = 0
for d in sys.argv[1:]:
    for f in sorted(glob.glob(f'{d}/*.npz')):
        z = dict(np.load(f))
        if 'dosE' not in z or 'dos_ref' in z or 'ef_ry' not in z:
            continue
        if z['dosE'].size:
            z['dosE'] = (z['dosE'].astype(np.float64) + float(z['ef_ry']) * RY).astype(np.float32)
        z['dos_ref'] = np.array('E_F')
        tmp = f[:-4] + '.tmp.npz'
        np.savez_compressed(tmp, **z); os.replace(tmp, f); n += 1
print(f'gw1500db_fixdos: {n} files corrected')
