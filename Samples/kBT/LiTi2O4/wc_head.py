#!/usr/bin/env python3
"""W_c(q,omega) from --dumpW __WVR.<iq> (direct access; records nblochpmx^2 complex(4) with --fp32).
usage: wc_head.py DIR iq nblochpmx [kp] [ndiag]
writes DIR/wc_head_iq<iq>.dat: iw omega(eV) Re Wc11 Im Wc11 Re Wc22 Re Wc33 max|diag| -Im tr(diag 1..ndiag)
and DIR/wc_diag_iq<iq>.dat: omega then Im Wc(i,i) for i=1..ndiag (plasmon spectrum per component)"""
import sys, numpy as np
d, iq, nb = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]); kp = int(sys.argv[4]) if len(sys.argv) > 4 else 4
nd = int(sys.argv[5]) if len(sys.argv) > 5 else 30
fr = [float(l.split()[0].replace('D','E')) for l in open(d+'/freq_r').readlines()[1:]]
dt = np.complex64 if kp == 4 else np.complex128
rec = nb*nb*2*kp
f = open(d+'/__WVR.%d' % iq, 'rb')
h = open(d+'/wc_head_iq%d.dat' % iq, 'w'); g = open(d+'/wc_diag_iq%d.dat' % iq, 'w')
h.write('# iw omega(eV) Re Wc(1,1) Im Wc(1,1)  Re Wc(2,2) Re Wc(3,3) max|diag(1:%d)|  -Im tr diag(1:%d)\n' % (nd, nd))
g.write('# omega(eV)  Im Wc(i,i) i=1..%d\n' % nd)
for iw, w in enumerate(fr):
    b = f.read(rec)
    if len(b) < rec: break
    m = np.frombuffer(b, dtype=dt).reshape(nb, nb, order='F')
    dg = np.diag(m)[:nd]
    h.write('%d %.4f %.4e %.4e %.4e %.4e %.4e %.4e\n' % (iw, w*27.211386, m[0,0].real, m[0,0].imag, m[1,1].real, m[2,2].real, np.abs(dg).max(), -dg.imag.sum()))
    g.write('%.4f ' % (w*27.211386) + ' '.join('%.3e' % x for x in dg.imag) + '\n')
