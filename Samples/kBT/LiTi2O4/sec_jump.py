#!/usr/bin/env python3
"""Iteration-to-iteration jumps of SEc in QPU.<n>run files.
For each iteration n>1: count of (k, state) with state in [lo,hi] whose SEc moved by more
than THR eV since iteration n-1, and the largest |dSEc|.  usage: sec_jump.py DIR [lo hi thr]"""
import sys, glob, re, os
d = sys.argv[1]; lo = int(sys.argv[2]) if len(sys.argv) > 2 else 9
hi = int(sys.argv[3]) if len(sys.argv) > 3 else 24; thr = float(sys.argv[4]) if len(sys.argv) > 4 else 1.5
def read(f):
    out = {}
    for line in open(f):
        t = line.split()
        if len(t) < 11: continue
        try: st = int(t[3]); q = tuple(round(float(x), 5) for x in t[:3]); sec = float(t[6])
        except ValueError: continue
        if lo <= st <= hi: out[(q, st)] = sec
    return out
files = sorted(glob.glob(os.path.join(d, 'QPU.*run')), key=lambda f: int(re.search(r'QPU\.(\d+)run', f).group(1)))
prev = None; print('# %s  states %d-%d  thr %.1f eV' % (d, lo, hi, thr))
for f in files:
    it = int(re.search(r'QPU\.(\d+)run', f).group(1)); cur = read(f)
    if prev:
        dd = [(abs(cur[k] - prev[k]), k) for k in cur if k in prev]
        n = sum(1 for x, _ in dd if x > thr); mx, km = max(dd)
        print('iter%-3d jumps=%3d/%d  max|dSEc|=%6.2f eV at q=%s st%d' % (it, n, len(dd), mx, km[0], km[1]))
    prev = cur
