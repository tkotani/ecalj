#!/usr/bin/env python
import sys,os
thisdir= os.path.dirname(os.path.abspath(__file__))
sys.path.append(thisdir+'/f2py/')
from fparser.api import parse

nargv = len(sys.argv) -1
argset= sys.argv[1:]
print argset

#fff='./xxx.F'
#fff='./bndfp.F'
#fff='./x.F'
#/home/takao/ecal//lm-7.0betaK001/fp/bndfp.F'
#print fff
#for ffile in argset:
#    print ffile
#
#sys.exit()

for ffile in argset:
    print '@@@@@ '+ffile+' @@@@@'
    tree = parse(ffile,isfree=False,isstrict=False,ignore_comments=False,analyze=True)
    #print dir(tree)
    #print tree.content
    print tree.torepr()
    #print tree.torepr(3)
    #tree
    #print tree.item


sys.exit()

#print tree
