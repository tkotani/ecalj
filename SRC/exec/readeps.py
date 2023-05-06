#!/usr/bin/env python3
'''
usage:
>readeps.py 
gives Im(eps(omega)) for EPS*.nlfc.dat.interbandonly files. 
This is needed for error cancellation.
'''
import sys,os,re,glob
import numpy as np
files=glob.glob('EPS*.nlfc.dat.interbandonly')
files.sort()
#print(files)
f1 = open(files[0],'rt').read().split('\n')

for iff,fff in enumerate(files):
    #print(iff,fff)
    if(iff==0): continue
    f2 = open(fff,'rt').read().split('\n')
    for ifnum,ix in enumerate(f1):
        if(ifnum==0): continue
        if(len(ix)==0): continue
        iline1= [float(re.sub('D','e',x)) for x in re.split('\s+',f1[ifnum]) if(len(x)>0)]
        iline2= [float(re.sub('D','e',x)) for x in re.split('\s+',f2[ifnum]) if(len(x)>0)]
        q1=iline1[0:3]
        q2=iline2[0:3]
        q1a=np.dot(q1,q1)**.5
        q2a=np.dot(q2,q2)**.5
        #print(iline1,iline2)
        ilineout=[]
        for ix,_ in enumerate(iline2):
            if(ix<4):
                ilineout.append(iline2[ix])
            else:
                val= ( q2a**2*iline2[ix] - q1a**2*iline1[ix])/(q2a**2-q1a**2)
                ilineout.append(val)
        print(ilineout[3],ilineout[4],ilineout[5])
    print()
    print()
    
