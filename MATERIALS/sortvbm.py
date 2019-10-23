#!/usr/bin/env python2
import re,sys,string,os,signal,time
vbm  = open("VBM.dat",'rt').read().split('\n')
vbmref= open("VBMref.dat",'rt').read().split('\n')

vbmdat={}
for line in vbm:
    if line=='' : continue
    if line[0]=='#': continue
    if line[0]=='=': continue
    line= re.sub('###','',line)
    liner=line.split('/')
    matid=liner[0].strip()
    #print matid,liner[1]
    vdat=liner[1].split('VBM=')[1]
    vdat=vdat.strip().split(' ')
    vbmdat[matid]=vdat[0]
print vbmdat
print vbmdat.keys()

########################################
for line in vbmref:
    if line=='' : continue
    if line[0]=='#': continue
    if line[0]=='=': continue
    dat=re.split(r'\s+',line)
    mat = dat[0]
    print mat,vbmdat[mat],dat[1],dat[2],dat[3],dat[4]

sys.exit()
