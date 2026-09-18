#!/usr/bin/env python3
"""Read SEX2U / SEC2U (hsfp0_sc, --ntqxx): <i|Sigma|j> in the band basis, Hartree.
usage: read_se2.py DIR   -> prints, for chosen q, the largest off-diagonal |Sigma_ij| (eV) of rows near E_F"""
import sys, struct, numpy as np
HA = 27.211386
def readf(fn):
    b = open(fn, 'rb').read(); p = 0
    def rec():
        nonlocal p
        n = struct.unpack('<i', b[p:p+4])[0]; p += 4
        d = b[p:p+n]; p += n + 4
        return d
    h = struct.unpack('<8i', rec()); nspin, nq, ntq = h[0], h[1], h[2]
    out = {}
    for ip in range(nq):
        d = rec()
        isx = struct.unpack('<i', d[:4])[0]; q = struct.unpack('<3d', d[4:28])
        z = np.frombuffer(d[28:], dtype=np.complex128, count=ntq*ntq).reshape(ntq, ntq, order='F')
        out[tuple(round(x, 5) for x in q)] = z
    return out, ntq
d = sys.argv[1]
sx, ntq = readf(d + '/SEBK/SEX2U'); sc, _ = readf(d + '/SEBK/SEC2U')
try: sxc, _ = readf(d + '/SEBK/SEXcore2U')
except Exception: sxc = {k: 0*v for k, v in sx.items()}
rows = [int(x) for x in sys.argv[2].split(',')] if len(sys.argv) > 2 else [33, 34]
for q in sorted(sx):
    S = (sx[q] + sxc[q] + 0.5*(sc[q] + sc[q].conj().T)) * HA   # eV, hermitian part
    nx = int(np.argmax(np.abs(np.diag(sx[q])) == 0)) or ntq
    for i in rows:
        r = np.abs(S[i-1, :nx]).copy(); r[i-1] = 0
        j = np.argsort(r)[::-1][:5]
        print('q=(%7.4f,%7.4f,%7.4f) nx=%3d row %2d diag %7.3f | top off-diag |S_ij| (eV): %s' % (
            *q, nx, i, S[i-1, i-1].real, ' '.join('j%d:%.2f' % (jj+1, r[jj]) for jj in j)))
