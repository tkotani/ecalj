#!/usr/bin/env python
import sys,os
thisdir= os.path.dirname(os.path.abspath(__file__))
sys.path.append(thisdir+'/f2py/')
from fparser.api import parse
from fparser.readfortran import *
#nargv = len(sys.argv) -1
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

#for ffile in argset:
#    print '@@@@@ '+ffile+' @@@@@ start -----'
#    reader=FortranFileReader(ffile)
    #reader.isfree=False
    #reader.isstrict=False
    #for item in reader:

# def get_item(reader):
#     try:
#         item = reader.next(ignore_comments = False)
#         return item
#     except StopIteration:
#         return 'eof'
#     return



for ffile in argset:
    print '@@@@@ '+ffile+' @@@@@ start -----'
    reader = FortranFileReader(ffile)
    reader.set_mode(isfree=False,isstrict=False)
    #for item in reader:
    for item in reader:
         print 'xxxxxx: ',item.__class__.__name__
         print 
         #for i in dir(item):
         #    print '   ', i
#         print '   ', dir(item.__class__)
#         print '   ', dir(item.__class__.__name__)
#     while 1:
#         item = get_item(reader) #reader.next(ignore_comments = False)
         #if(item=='eof'): break
         if(item.__class__.__name__=='Line'): 
             print item.span, item.__class__.__name__,item.line
             print '   reader    : ', item.reader
             print '   label     : ', item.label
             print '   srtline   : ', item.strline
             print '   strlinemap: ', item.strlinemap
         elif(item.__class__.__name__=='Comment'):
             print item.span, item.comment
             
    print '@@@@@ '+ffile+' @@@@@ end -------'
    print 
    #print 'RRRRRR@ ', i,aaa
    #tree = parse(ffile,isfree=False,isstrict=False,ignore_comments=False,analyze=True)
    #print tree.content
    #print tree.torepr(6)
    #tree
    #print tree.item

sys.exit()


