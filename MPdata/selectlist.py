#!/usr/bin/env python3
import sys,os

args = sys.argv
ffile=args[1]

fff=open(ffile).read()
ff=fff.split("\n")
count=0
for line in ff:
    if(len(line)==0): continue
    bg=float(line.split('}')[1].split()[0])
    if(bg>0): count=count+1 #band gap
length=count
print(length)

ndiv=10
nsize=length//ndiv+1 
#print(ndiv,nsize)
count=0
i=0
fnum=0
for line in ff:
    if(len(line)==0): continue
    #print(line)
    bg=float(line.split('}')[1].split()[0])
    #print(bg)
    if(bg>0):
        count=count+1
        #print(count, line)
        if count % nsize==1:
            fnum=fnum+1
            if(fnum!=1): ffw.close()
            ffw=open(ffile+'.'+str(fnum),'w',)
        print(line,file=ffw)
        
sys.exit()
