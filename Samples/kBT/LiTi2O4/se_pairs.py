#!/usr/bin/env python3
"""Near-E_F pairs: Sigma_c_ij evaluated at e_i (row i) vs at e_j (row j, conj): the omega dependence of the
off-diagonal seen by the static QSGW average. usage: se_pairs.py DIR [ilo ihi dmax_eV]"""
import sys, struct, re, numpy as np
HA=27.211386
def readf(fn):
    b=open(fn,"rb").read(); p=0
    def rec():
        nonlocal p
        n=struct.unpack("<i",b[p:p+4])[0]; p+=4; d=b[p:p+n]; p+=n+4; return d
    h=struct.unpack("<8i",rec()); nq,ntq=h[1],h[2]; out={}
    for ip in range(nq):
        d=rec(); q=struct.unpack("<3d",d[4:28])
        out[tuple(round(x,5) for x in q)]=np.frombuffer(d[28:],dtype=np.complex128,count=ntq*ntq).reshape(ntq,ntq,order="F")
    return out
d=sys.argv[1]; lo=int(sys.argv[2]) if len(sys.argv)>2 else 30; hi=int(sys.argv[3]) if len(sys.argv)>3 else 44; dmax=float(sys.argv[4]) if len(sys.argv)>4 else 1.5
sc=readf(d+"/SEBK/SEC2U")
e={}
for l in open(d+"/QPU"):
    t=l.split()
    if len(t)>=13 and re.match(r"^-?\d",t[0]):
        try: e[(tuple(round(float(x),5) for x in t[:3]),int(t[3]))]=float(t[10])
        except: pass
print("# q   i  j   e_i    e_j   Delta |  Sc_ij(e_i)  Sc_ji(e_j)*  (complex, eV) | |diff|  |diff|/Delta  |avg|   avg^2/Delta")
rows=[]
for q in sorted(sc):
    S=sc[q]*HA
    for i in range(lo-1,hi):
        for j in range(i+1,hi):
            ei=e.get((q,i+1)); ej=e.get((q,j+1))
            if ei is None or ej is None: continue
            D=abs(ei-ej)
            if D>dmax or D<0.02: continue
            a=S[i,j]; b=np.conj(S[j,i]); dd=abs(a-b); av=abs(a+b)/2
            rows.append((dd/D,q,i+1,j+1,ei,ej,D,a,b,dd,av))
rows.sort(reverse=True)
for r in rows[:14]:
    print("(%5.2f,%5.2f,%5.2f) %2d %2d %6.2f %6.2f %5.2f | %6.2f%+6.2fi %6.2f%+6.2fi | %5.2f  %5.2f  %5.2f  %5.2f"%(*r[1],r[2],r[3],r[4],r[5],r[6],r[7].real,r[7].imag,r[8].real,r[8].imag,r[9],r[0],r[10],r[10]**2/r[6]))
print("# pairs %d; median |diff|/Delta = %.2f, median |avg| = %.2f eV"%(len(rows),np.median([r[0] for r in rows]),np.median([r[10] for r in rows])))
