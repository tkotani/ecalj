#!/usr/bin/env python
import sys,os
sys.path.insert(0,'../TOOLS/f2py/fparser')
from api import parse
#for dirpath,dirnames,filenames in os.walk('.'):
#    print dirpath,dirnames,filenames
#sys.exit()

srcdir='./FPLOTdir/'
files =os.listdir(srcdir)
print files
src=[]
for file in files:
    if(os.path.splitext(file)[1]=='.F'): src.append(srcdir+file)
print src

#./FPLOTdir/fplot.F  ./FPLOTdir/fpsub.F ./FPLOTdir/plbnds.F ./FPLOTdir/pldos.F  ./FPLOTdir/plsub.F'
#fff='./bndfp.F'
#fff='./x.F'
#/home/takao/ecal//lm-7.0betaK001/fp/bndfp.F'
#print fff
#src=[srcdir+'fpsub.F']
for file in src:
    print '--------------',file,'------------------'
    tree = parse(file,isfree=False,isstrict=False,ignore_comments=False)
    print tree.content
#print tree




