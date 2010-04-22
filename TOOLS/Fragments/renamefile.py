#!/usr/bin/env python
import sys,string,re
#sys.path.append('/home/takao/ecal/TOOLS/f2py')
#from api import parse
#fff='./suham.F'
nargv = len(sys.argv) -1
argset= sys.argv[1:]
print argset
fff='./xxx.F'
#fff='./bndfp.F'
#fff='./x.F'
#/home/takao/ecal//lm-7.0betaK001/fp/bndfp.F'
#print fff
re_head = re.compile( '^\#', re.IGNORECASE)
for ffile in argset:
    fin = open(ffile,'rt').read()
    fout=ffile.split('.o')[0]+'.F'
    oxx= string.split(fin,'\n') 
    fileF=''
    for ix in range( len(oxx)):
        fline=oxx[ix]
        if(re_head.search(fline) ): 
            pass
        else:
            fileF=fileF+fline+'\n'
    #tree = parse(fileF,isfree=False,isstrict=False,ignore_comments=False)
    #print tree.content
    f=open(fout,'w')
    f.write(fileF)
sys.exit()

#print tree


