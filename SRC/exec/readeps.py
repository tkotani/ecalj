#!/usr/bin/env python3
'''
usage:
>readeps.py 
gives Im(eps(omega)),Re(eps) from EPS*.nlfc.dat.interbandonly files. 
This is needed for error cancellation.
'''
import sys,os,re,glob
import numpy as np

def cleanl(iline):
    iout= [float(re.sub('D','e',x)) for x in re.split('\s+',iline) if(len(x)>0)]
    return iout

def fd0(x,f,dig):
    c=np.polyfit(x[0:dig+1],f[0:dig+1],dig)
    return c[dig-1]
    
files=glob.glob('EPS*.nlfc.dat.interbandonly')
files.sort()
epsfile=[open(fff,'rt').read().split('\n') for fff in files]
dlines=['']*len(epsfile)
q2=['']*len(epsfile)

nfile=len(epsfile)
for idf,iff in enumerate(range(nfile)):
    dlines[idf]=[cleanl(iline) for iline in epsfile[iff][1:] if(len(cleanl(iline)))>0]
for idf in range(nfile):   # idf is file ID for EPS*
    q= dlines[idf][0][0:4] # q vector
    q2[idf]=np.dot(q,q)    # q**2
    print('# ',idf,q2[idf],dlines[idf][0][0:4],files[idf])
for idat in range(len(dlines[0])):
    idf=1
    omega=  dlines[0][idat][3]     #eps2  = np.array([        dlines[idf][idat][5] for idf in range(nfile)])
    q2eps1= [q2[idf]*dlines[idf][idat][4] for idf in range(nfile)] #list of q**2*eps1 in EPS*
    q2eps2= [q2[idf]*dlines[idf][idat][5] for idf in range(nfile)] #list of q**2*eps2 in EPS*
    eps1out = fd0(q2, q2eps1,2)
    eps2out = fd0(q2, q2eps2,2)
    print(idat,omega, \
          eps2out,fd0(q2[1:],q2eps2[1:],2),fd0(q2[2:],q2eps2[2:],2),'  ', \
          eps1out,fd0(q2[1:],q2eps1[1:],2),fd0(q2[2:],q2eps1[2:],2))
    #,fd0(q2[3:],q2eps2[3:],2),fd0(q2[4:],q2eps2[4:],2),\#,fd0(q2[3:],q2eps1[3:],2),fd0(q2[4:],q2eps1[4:],2))
#gnuplot -p -e 'set xran[0:30];plot "epsinter.dat" using ($2)*13.605:($3) w l,"" using ($2)*13.605:($4) w l,"" using ($2)*13.605:($5) w l,"" using ($2)*13.605:($6) w l,"" using ($2)*13.605:($7) w l,"" using ($2)*13.605:($8) w l'
sys.exit()
