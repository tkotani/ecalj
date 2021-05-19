#!/usr/bin/env python3
import sys,os,re
thisdir= os.path.dirname(os.path.abspath(__file__))
#sys.path.append(thisdir+'/FparserTools/fparser/')
#from api import parse,walk
#from readfortran import *
#from parsefortran import FortranParser
#from inspect import *
#from base_classes import classes
#import block_statements

#nargv = len(sys.argv) -1
argset= sys.argv[1:]
#print argset

#print classes.typedecl_statements.Type
#print classes.block_statements.Type
#
#print typedecl_statements.Type
#sys.exit()
#for i in classes:
#    print i


functions=[] # save defined functions in array.

#print '<graphviz>'
#print 'digraph G {'


def replacer(modz,lineb):
    modx=modz[0]
    mody=modz[1]
    val =modz[2]
    zzz=None
    if re.match(r'.*use '+modx,lineb):
        print('yyy',ixb,lineb)
        zzz=re.sub(val+'[\s]*?,','',lineb)
        zzz=zzz+'      use '+mody+',only: '+val
    return zzz

modz=['m_lmfinit','m_lattic','lat_nkd']


for ffile in argset:
    print('=========== '+ffile+' start ===========')
    f=open(ffile)
    lines = f.readlines()
    f.close()
    aaa=''
    lineb=''
    ixb=-1
    zzz=''
    for ix,line in enumerate(lines):
        if len(line)<6 or re.match('\S',line[0]) or \
           (not re.match(r'\S',line[5])):
            zzzr = replacer(modz,lineb)
            #print(ixb,lineb)
            if(zzzr):
                print(zzzr)
                zzz=zzz+zzzr+'\n'
                print('yyy',zzzr)
            else:
                zzz=zzz+lineb
            lineb=line
            ixb=ix+1
            if ix==len(lines)-1 : zzz=zzz+line
        else: #get a block or lines
            lineb=lineb+line
    #print(aaa)        
    fout=open(ffile,mode='w')
    fout.write(zzz)
    fout.close()
