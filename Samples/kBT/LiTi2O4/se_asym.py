#!/usr/bin/env python3
"""Anti-hermitian part of <i|Sigma_c(e_i)|j> vs <j|Sigma_c(e_j)|i>* in the band basis: the seesaw diagnostic.
usage: se_asym.py DIR [ilo ihi]  (bands 25..52 = t2g region by default). Prints, per q, the largest |A_ij| with e_i,e_j."""
import sys, struct, numpy as np
HA=27.211386
def readf(fn):
    b=open(fn,"rb").read(); p=0
    def rec():
        nonlocal p
        n=struct.unpack("<i",b[p:p+4])[0]; p+=4; d=b[p:p+n]; p+=n+4; return d
    h=struct.unpack("<8i",rec()); nq,ntq=h[1],h[2]; out={}
    for ip in range(nq):
        d=rec(); q=struct.unpack("<3d",d[4:28])
        z=np.frombuffer(d[28:],dtype=np.complex128,count=ntq*ntq).reshape(ntq,ntq,order="F")
        out[tuple(round(x,5) for x in q)]=z
    return out,ntq
d=sys.argv[1]; lo=int(sys.argv[2]) if len(sys.argv)>2 else 25; hi=int(sys.argv[3]) if len(sys.argv)>3 else 52
sc,ntq=readf(d+"/SEBK/SEC2U")
# eigenvalues from QPU (elda column) 
import re
e={}
for l in open(d+"/QPU"):
    t=l.split()
    if len(t)>=13 and re.match(r"^-?\d",t[0]):
        try: e[(tuple(round(float(x),5) for x in t[:3]),int(t[3]))]=float(t[10])
        except: pass
print("# q  i j  e_i e_j (eV)  |S_ij(e_i)| |S_ji(e_j)*|  |A_ij|=|S_ij(e_i)-S_ji(e_j)*|/2   (eV)")
tot=[]
for q in sorted(sc):
    S=sc[q]*HA; A=0.5*(S-S.conj().T)
    best=[]
    for i in range(lo-1,hi):
        for j in range(i+1,hi):
            best.append((abs(A[i,j]),i+1,j+1,abs(S[i,j]),abs(S[j,i])))
    best.sort(reverse=True)
    for a,i,j,s1,s2 in best[:3]:
        ei=e.get((q,i),float("nan")); ej=e.get((q,j),float("nan"))
        print("q=(%7.4f,%7.4f,%7.4f) %3d %3d  %7.3f %7.3f   %6.2f %6.2f   %6.2f"%(*q,i,j,ei,ej,s1,s2,a))
    tot.append(best[0][0])
print("# max over q of max|A_ij|: %.2f eV, mean %.2f"%(max(tot),sum(tot)/len(tot)))
