#!/usr/bin/env python
# Program unit analyzer.
import sys,os
thisdir= os.path.dirname(os.path.abspath(__file__))
sys.path.append(thisdir+'/f2py/fparser')
from api import parse,walk
from readfortran import *
from parsefortran import FortranParser
from inspect import *
from base_classes import classes
import block_statements

nargv = len(sys.argv) -1
argset= sys.argv[1:]

for ffile in argset:
#    print '---- '+ffile+' start -----'
    reader = FortranFileReader(ffile)
    reader.set_mode(isfree=False,isstrict=False)
    parser=FortranParser(reader,ignore_comments=False)
    parser.parse()
    parser.analyze()
    sstack=[]
    deptho=-1
    inso=None
    for c in walk(parser.block):
        ins=c[0]
        depth=c[1]
        if(isinstance(ins, classes.Comment)): continue
        if(depth>deptho):         sstack.append(inso)
        while len(sstack)>depth:  sstack.pop()
        if( "item" in dir(ins)) :
            item = ins.item
        else:
            item = ins
        print depth,item.span,item.strline,'--->',item.strlinemap
        #print len(sstack),depth, type(ins), item.span,  item.strline, item.label
        deptho=depth
        inso=ins
        continue
