#!/usr/bin/env python
import sys,os
thisdir= os.path.dirname(os.path.abspath(__file__))
#sys.path.append(thisdir+'/fparser')
sys.path.append(thisdir+'/f2py/fparser')
#from api import parse
from readfortran import *
#from parsefortran import FortranParser
#nargv = len(sys.argv) -1
argset= sys.argv[1:]
print argset
from inspect import *


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
#   print dir(reader)
#   sys.exit()

    for item in reader:
        classname = item.__class__.__name__
        if(classname=='Line'):
            print item.span,item
        else:
            print ' xxx',item.span,item
#            sys.exit() 


#    for item in reader:
#        print item
#        print item.span()
    sys.exit()









#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#from api import parse
#from parsefortran import FortranParser


### test1 #############################################################
#    for i in getmembers(reader):
#        print i

### test2 #####################
    reader.set_mode(isfree=False,isstrict=False)
    parser=FortranParser(reader,ignore_comments=False)
    parser.parse()
    parser.analyze()
#    print parser.block.torepr(3)
#    print parser.block
    for i in dir(parser.block):
        print i
    print '----------------------'
    print 'aaa ',parser.block.name
    print 'aaa ',parser.block.blocktype

    ix=0
    for i in parser.block.content:
        ix=ix+1
        try:
            print 'zzzzzzzzz', ix,'  ',i.blocktype #,i.name,i.item #,i.a
            print ' a =',i.a
            iy=0
            for j in i.content:
               iy=iy+1
               try:
                 print 'bbb ',iy,'  ',j.__class__.__name__,
                 print ' blockname=',j.name,
                 print ' item=',j.item.strline,' map=', j.item.strlinemap,' span=',j.item.span
                 print ' item=',j.item.line
                 print ' a =',dir(j.a)
                 print ' content=',j.content,
                 print ' blocktype=',j.blocktype
               except:
                 print 
                 print 'eee ',iy,'  ',j.__class__.__name__

        except:
            print 'eeeeeeee', ix,i.__class__.__name__,'\n',i.item.comment,i.item.span

#    for i in parser.block.classes:
#        print i



    sys.exit()


###########################################################
    ix=0
    print '-------dir parser.block  ------------'
    for i in dir(parser.block):
        ix=ix+1
        print ' zzzzzzz ix=',ix,' ',i
    print '-------------------'
    print 
    print '-------dir parser.block.a ------------'
    ix=0
    for i in dir(parser.block.a):
        ix=ix+1
        print ' vvvvvv  ix',ix,' ',i
    print '-------------------'
    print 'class name of parser.block.a=',parser.block.a.__class__.__name__
    print parser.block.a #this is given by AttributeHolder.__repr__ base_classes.py L83
    #print parser.block.a.torepr.external_subprogram #a is Attribute class
################################################################



###########################################
#    print parser.block.__doc__
    sys.exit()

###########################################
#    ix=0
#    for i in getmembers(parser.block):
#        ix=ix+1
#        print 'zzzzzzz ix=',ix,' ',i
#    sys.exit()
###################


    print 
    #print parser.block.__dict__
    aaa=parser.block.__dict__
    for k in aaa.keys():
        print 'key=', k
    #print aaa['a']
    sys.exit()

    sys.exit()

################################################################



    reader.set_mode(isfree=False,isstrict=False)
    #for item in reader:
    parser=FortranParser(reader,ignore_comments=False)
    parser.parse()
    parser.analyze()
#    print parser.block.torepr(depth=2) 

    print parser.block.torepr(depth=2) 
          #this is in FortranParser->BeginSource->Beginstatement
sys.exit()





