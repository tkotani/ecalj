#!/usr/bin/env python
##!/usr/bin/env python2
# takao kotani jan2017 for symmetry line based on https://github.com/giovannipizzi/seekpath
#
import numpy as np
import sys,string,re,os,math,seekpath
from math import *
import getpaths
np.set_printoptions(precision=16)

if len(sys.argv)!=2:
    print ' Usage: >getsyml.py nio !!!'
    print '   for ctrl.nio. We run lmchk internally in getsyml !!!'
    print 'this is based on https://github.com/giovannipizzi/seekpath'
    print 'we have to cite Y. Hinuma, G. Pizzi, Y. Kumagai, F. Oba, I. Tanaka, Band structure diagram paths based on crystallography, Comp. Mat. Sci. 128, 140 (2017) (JOURNAL LINK, arXiv link).'
    sys.exit(-1)


### We have to improve lmchk, so that plinfo and siteinfo are directry written.
os.system('lmchk '+sys.argv[1]+'> outlmchk')
#os.system('grep Plat -A3 outlmchk > plinfo')
plfile = open('PlatQlat.lmchk','r').read().split('\n')
i=0
plat=['']*3
for iline in plfile:
    if(i==3): break
    i=i+1
    plat[i-1] = [float(re.split('\s+',iline)[ix]) for ix in range(1,4)]
cell = plat
print 'primitive cell=',cell

### qlat = inverse plat
qlat=[]
vol=np.dot(cell[0],np.cross(cell[1],cell[2]))
qlat.append(np.cross(cell[1],cell[2])/vol)
qlat.append(np.cross(cell[2],cell[0])/vol)
qlat.append(np.cross(cell[0],cell[1])/vol)
for i in range(0,3):
    for j in range(0,3):
        if(i!=j):
            if np.dot(qlat[i],cell[j])>1e-5:
                print 'qlat i/=j'
                sys.exit(-1)
        else:
            if abs(np.dot(qlat[i],cell[j])-1e0)>1e-5:
                print 'qlat i=j'
                sys.exit(-1)

###
aaa='"Site     Spec            Rmax"'
#os.system('grep -A10000 '+aaa+ ' outlmchk > siteinfo')
sitefile = open('SiteInfo.lmchk','r').read().split('\n')

i=0
pos=['']*3
positions=[]
numbers=[]
#print sitefile
for iline in sitefile:
    i=i+1
    #if i==1: continue
    if len(iline)==0: 
        nsite=i-1
        break
    vec=[re.split('\s+',iline)[ix] for ix in range(1,8)]
#    print iline
    pos = [float(vec[ix]) for ix in range(4,7)]
    positions.append([np.dot(pos,qlat[ix]) for ix in range(0,3)])
    numbers.append(int(vec[1]))
    print i,' pos=',pos,' id=',int(vec[1])
print 'nsite=',nsite#,positions


###############################################
structure = (cell, positions, numbers)
with_time_reversal=True
a = getpaths.get_path(structure,with_time_reversal)
#for i in a.keys():
#    print i,a[i]
symlfile='syml.'+sys.argv[1]
sfile = open(symlfile,'w')
ppp = ['']*3


###############################################
# print '=== qlat==='
# print np.array([qlat[i].tolist() for i in range(0,3)])
# print '=== qreci ==='
# print np.array([qrr[i].tolist() for i in range(0,3)])
# print '=== <qlat|qlat> matrix ==='

### qlat is given primitive lattice.
qqmat=np.array([[0e0]*3]*3)
for i in range(0,3):
    for j in range(0,3):
        qqmat[i][j]=np.dot(qlat[i],qlat[j])
#print qqmat

### qrr is standerized primitive lattice
#print '=== <qreci|qreci> matrix ==='
qrr = np.array(a['reciprocal_primitive_lattice'])/2e0/math.pi
qcmat=np.array([[0e0]*3]*3)
for i in range(0,3):
    for j in range(0,3):
        qcmat[i][j]=np.dot(qrr[i],qrr[j])
#print qcmat

###########################################################
### Among folloing points3d, we look for Q1,Q2,Q3, which gives <Q_i|Q_j>= qcmat(i,j)
supercell_size = 4 # Is this enough?
points3d = []
for i in range(-supercell_size, supercell_size+1):
    for j in range(-supercell_size, supercell_size+1):
        for k in range(-supercell_size, supercell_size+1):
            points3d.append( i*np.array(qlat[0]) + j*np.array(qlat[1]) + k*np.array(qlat[2]) )
#print 'len(possible point3d)=',len(points3d)
ipx=0
candidate=[[],[],[]] # candidates for Q1,Q2,Q3
for ip in points3d:
    ipx=ipx+1
    norm=np.dot(ip,ip)
    for i in range(0,3):
        if(np.abs(qcmat[i][i]-np.dot(ip,ip))<1e-5):
            candidate[i].append(ip)
#            print ipx,'->',i
# ix=0
# for i in candidate:
#     ix=ix+1
#     print '--- candidate for i=',ix
#     for j in i:
#         print j
# print '---------------'
qmat=np.array([[1e9]*3]*3)
qqq =np.array([[1e9]*3]*3)
ifound=False
for q1 in candidate[0]:         #candidates of Q1 as linear combination of qlat
    if ifound:break
    for q2 in candidate[1]:     #candidates of Q2 as linear combination of qlat
        if ifound:break
        for q3 in candidate[2]: #candidates of Q3 as linear combination of qlat
            if ifound:break
            qq=[q1,q2,q3]        
            vol=np.dot(q1,np.cross(q2,q3))
            if(abs(vol)<1e-5): continue
            ix=0
            for i in qq:
                ix=ix+1
                iy=0
                for j in qq:
                    iy=iy+1
                    qmat[ix-1,iy-1]=np.dot(i,j)
            diff = sum(sum(abs(qcmat-qmat)))
            if(diff<1e-5):
                #print 'qcmat-qmat=',diff
                qqq[0]=q1
                qqq[1]=q2
                qqq[2]=q3
                print 'OK! We found qqq vectors, which are equilanent to the reciprocal_primitive_lattice.'
                print 'We can use qqq1,qqq2,qqq3 instead of reciprocal_primitive_lattice'
                print 'Following qqq1,qqq2,qqq3 are linear combinations of primitive vectors qlat.'
                print '  qqq1=',qqq[0] # We can use qqq1,qqq2,qqq3 instead of reciprocal_primitive_lattice
                print '  qqq2=',qqq[1] # Here qqq1,qqq2,qqq3 are linear cobination of qlat.
                print '  qqq3=',qqq[2]
                ifound=True
                break
if not ifound:
    print 'can not find qqq equivalent to reciprocal_primitive_lattice'
    sys.exit()

distot = 0e0
for ipath in a['path']:
    ii=ipath[0]
    ppp= a['point_coords'][ii]
    ppp= ppp[0]*qqq[0]+ppp[1]*qqq[1] +ppp[2]*qqq[2]
    ii1,ii2,ii3 = [ float(ppp[ixx]) for ixx in range(0,3)]
    ee=ipath[1]
    ppp= a['point_coords'][ee]
    ppp= ppp[0]*qqq[0]+ppp[1]*qqq[1] +ppp[2]*qqq[2]
    ee1,ee2,ee3 = [ float(ppp[ix]) for ix in range(0,3)]
    dis= ((ee1-ii1)**2+(ee2-ii2)**2+(ee3-ii3)**2)**.5
    distot = distot+dis

###############################################
totalbandpoint=100
for ipath in a['path']:
    ii=ipath[0]
    ppp= a['point_coords'][ii]
    #print 'pppii=',ii,ppp
    ppp= ppp[0]*qqq[0]+ppp[1]*qqq[1] +ppp[2]*qqq[2]
    ii1,ii2,ii3 = [ float(ppp[ixx]) for ixx in range(0,3)]
    ee=ipath[1]
    ppp= a['point_coords'][ee]
    ppp= ppp[0]*qqq[0]+ppp[1]*qqq[1] +ppp[2]*qqq[2]
    ee1,ee2,ee3 = [ float(ppp[ix]) for ix in range(0,3)]
    dis= ((ee1-ii1)**2+(ee2-ii2)**2+(ee3-ii3)**2)**.5
    ndiv=int(totalbandpoint*dis/distot) #this controls number of divisions
    print ii1,ii2,ii3,ii
    print >> sfile, ndiv,ii1,ii2,ii3,' ',ee1,ee2,ee3,'  ', ii,ee 
print >>sfile, 0
print
sfile.close()
print 'OK! Check ',symlfile,' file!----------'
print 
from brillouinzone import brillouinzone_takao
brillouinzone_takao.plotws(qlat[0],qlat[1],qlat[2])

