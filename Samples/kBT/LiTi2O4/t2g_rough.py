#!/usr/bin/env python3
"""Roughness (2nd difference along k, meV) of the near-E_F t2g bands 33-40 (occupied bottom .. +1.5 eV).
usage: t2g_rough.py DIR [DIR...]"""
import sys
sys.path.insert(0, "/home/takao/LiTi2O4/kbt/runs")
from o2p_rough import read, rough
lo, hi = 33, 40
print("# dir  t2g rough(meV) bands %d-%d: mean max [per band]" % (lo, hi))
for d in sys.argv[1:]:
    b = read(d)
    r = [rough(b, ib) for ib in range(lo, hi+1)]
    print("%-50s %7.1f %7.1f  [%s]" % (d, sum(r)/len(r), max(r), " ".join("%.1f" % x for x in r)))
