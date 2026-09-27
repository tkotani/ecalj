#!/usr/bin/env python3
# cmpef.py ref res [res ...] : Re Sigma_c (SECU col 7) and Sigma_x (SEXU col 7) differences from ref, eV,
# over all states and over |E - E_F| < 1 eV and < 3 eV (col 6 = E - E_F in eV).
import sys
def rows(f):
    out = []
    for l in open(f, errors="replace"):
        w = l.split()
        if len(w) < 8 or not w[0].isdigit():
            continue
        try:
            out.append([float(x.replace("D", "E")) for x in w])
        except ValueError:
            pass
    return out
ref = sys.argv[1]
for r in sys.argv[2:]:
    out = [r]
    for f in ("SECU", "SEXU"):
        a, b = rows(f"{ref}/{f}"), rows(f"{r}/{f}")
        n = min(len(a), len(b))
        for w in (1e9, 3.0, 1.0):
            d = [abs(a[i][7] - b[i][7]) for i in range(n) if abs(a[i][6]) < w]
            lab = "all" if w > 1e8 else f"|E-EF|<{w:g}"
            out.append(f"{f[:3]} {lab}: max {max(d)*1e3:.3f} rms {((sum(x*x for x in d)/len(d))**0.5)*1e3:.3f} meV (n={len(d)})")
    print("\n   ".join(out))
