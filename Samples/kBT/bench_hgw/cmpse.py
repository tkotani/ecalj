#!/usr/bin/env python3
# cmpse.py res_A res_B : max |difference| of the numbers in SECU / SEXU (text), per column.
import sys
def rows(f):
    out=[]
    for l in open(f, errors="replace"):
        w=l.split()
        try: v=[float(x.replace("D","E")) for x in w]
        except ValueError: continue
        if len(v)>=6: out.append(v)
    return out
a,b=sys.argv[1],sys.argv[2]
for f in ["SEXU","SECU"]:
    ra,rb=rows(a+"/"+f),rows(b+"/"+f)
    n=min(len(ra),len(rb)); ncol=min(len(ra[0]),len(rb[0]))
    dmax=[max(abs(ra[i][c]-rb[i][c]) for i in range(n)) for c in range(ncol)]
    amax=[max(abs(ra[i][c]) for i in range(n)) for c in range(ncol)]
    print(f"{f}: rows {len(ra)} {len(rb)}  max|diff| per column: "+" ".join(f"{d:.2e}" for d in dmax))
    print(f"{f}:                 max|value| per column: "+" ".join(f"{x:.2e}" for x in amax))
