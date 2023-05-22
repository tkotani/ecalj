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
    return (f[1]-f[0])/(x[1]-x[0])
    #return c[dig-1]
### interband part ###############################    
files=glob.glob('EPS*.nlfc.dat.interbandonly')
files.sort()
epsfile=[open(fff,'rt').read().split('\n') for fff in files]
dlines=['']*len(epsfile)
q2=['']*len(epsfile)
nfile=len(epsfile)
for idf,iff in enumerate(range(nfile)): # idf is file ID for EPS*
    dlines[idf]=[cleanl(iline) for iline in epsfile[iff][1:] if(len(cleanl(iline)))>0] #data line for the file idf
    q= dlines[idf][0][0:4] # q vector
    q2[idf]=np.dot(q,q)    # q**2
    print('# ',idf,dlines[idf][0][0:4],files[idf])
dinter=[]
dig=1
for idat in range(len(dlines[0])):
    omega=  dlines[0][idat][3]     #Note: bare eps2= np.array([ dlines[idf][idat][5] for idf in range(nfile)])
    q2eps1= [q2[idf]*dlines[idf][idat][4] for idf in range(nfile)] #list of q**2*eps1 in EPS*
    q2eps2= [q2[idf]*dlines[idf][idat][5] for idf in range(nfile)] #list of q**2*eps2 in EPS*
    eps1=[fd0(q2[idf:],q2eps1[idf:],dig) for idf in range(nfile-dig)]
    eps2=[fd0(q2[idf:],q2eps2[idf:],dig) for idf in range(nfile-dig)]
    dinter.append([idat,omega,*eps1,*eps2])
if(nfile>0):
    f=open('epsinter.dat','w')
    for idf in range(nfile): print('# ',idf,dlines[idf][0][0:4],files[idf],file=f)
    print('# recid omega eps1.... eps2... where num of eps is the number of EPSfiles',file=f) 
    for idat in dinter:
        print(*idat,file=f)
    f.close()
    f=open('epsinter.glt','w') 
    print()
    fdata=''
    for iddat in range(3,nfile+3-dig):
        fdata=fdata+' '+'"epsinter.dat" using ($2)*13.605:($'+str(iddat)+') w l title "RealPart:'+ files[iddat-3]+'"'
        fdata=fdata+',\\\n'
    for iddat in range(3,nfile+3-dig):
        fdata=fdata+' '+'"epsinter.dat" using ($2)*13.605:($'+str(iddat+nfile-dig)+') w l title   "ImagPart:'+ files[iddat-3]+'",\\\n'
    aaa='set title "IntraBand part of Epsilon(omega(eV))"\nset xlabel "(eV)"\nset datafile fortran\nset xran[0:30]\nset yran[-30:30]\nplot\\\n' + fdata
    print(aaa,file=f)
    f.close()
### intraband part ###############################    
files=glob.glob('EPS*.nlfc.dat.intrabandonly')
files.sort()
epsfile=[open(fff,'rt').read().split('\n') for fff in files]
dlines=['']*len(epsfile)
q2=['']*len(epsfile)
nfile=len(epsfile)
for idf,iff in enumerate(range(nfile)): # idf is file ID for EPS*
    dlines[idf]=[cleanl(iline) for iline in epsfile[iff][1:] if(len(cleanl(iline)))>0] #data line for the file idf
    print('# ',idf,dlines[idf][0][0:3],files[idf])
dintra=[]
for idat in range(len(dlines[0])):
    omega= dlines[0][idat][3]     
    eps1= [dlines[idf][idat][4] for idf in range(nfile)] 
    eps2= [dlines[idf][idat][5] for idf in range(nfile)] 
    dintra.append([idat,omega,*eps1,*eps2])
f=open('epsintra.dat','w')
for idf in range(nfile): print('# ',idf,dlines[idf][0][0:4],files[idf],file=f)
print('# recid omega eps1.... eps2... where num of eps is the number of EPSfiles',file=f) 
for idat in dintra:
    print(*idat,file=f)
f.close()    
print()
f=open('epsintra.glt','w')
fdata=''
for iddat in range(3,nfile+3):
    fdata=fdata+' '+'"epsintra.dat" using ($2)*13.605:($'+str(iddat)+') w l title "RealPart:'+ files[iddat-3]+'",\\\n'
for iddat in range(3,nfile+3):
    fdata=fdata+' '+'"epsintra.dat" using ($2)*13.605:($'+str(iddat+nfile)+') w l title  "ImagPart:'+ files[iddat-3]+'",\\\n'
aaa='set title "IntraBand part of Epsilon(omega(eV))"\nset xlabel "(eV)"\nset datafile fortran\nset xran[0:30]\nset yran[-30:30]\n'+\
    'plot\\\n' + fdata
print(aaa,file=f)
f.close()
#### epsall.data = epsinter+epsintra. Only with 1st data ####    
f=open('epsall.dat','w')
for idat in range(len(dintra)):
    eps1= dinter[idat][2]
    eps2= dinter[idat][2+nfile-2]
    id,omega=dintra[idat][0:2]
    if(len(dintra)!=0):
        eps1+=dintra[idat][2]
        eps2+=dintra[idat][2+nfile]
    print(id,omega, eps1-1, eps2,file=f)
f.close()    
f=open('epsall.glt','w')
aaa='set title "Epsilon(omega(eV))"\nset xlabel "(eV)"\nset datafile fortran\nset xran[0:30]\nset yran[-30:30]\n'+\
    'plot\\\n' + '"epsall.dat" using ($2)*13.605:($3) w l title "RealPart","" using ($2*13.605):($4) w l title "ImagPart"'
print(aaa,file=f)
f.close()
print('OK! We have generated epsinter.dat,epsinter.glt, epsintra.*(if metal) epsall.*')
