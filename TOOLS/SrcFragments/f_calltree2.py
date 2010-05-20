#!/usr/bin/env python
import sys,os
thisdir= os.path.dirname(os.path.abspath(__file__))
sys.path.append(thisdir+'/f2py/fparser')
from api import parse,walk
from readfortran import *
from parsefortran import FortranParser
from inspect import *
from base_classes import classes
import block_statements

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

print '<graphviz>'
print 'diagram G {'


for ffile in argset:
    #print '---- '+ffile+' start -----'
    reader = FortranFileReader(ffile)
    reader.set_mode(isfree=False,isstrict=False)
    parser=FortranParser(reader,ignore_comments=False)
    parser.parse()
    #parser.analyze()

    neg='(dpzero|defcc|ll|defrc|lgunit|iprint|rx|ival|isum|' \
    + 'nglob|cmdopt|rlse|defrr|zsum|pack|upack|dpadd|awrit|dcopy|daxpy|dscal|def|daxpy|' \
    + 'r8tos8|dvset|dpscop|dmpy|fexit|sp2cls|redfrr|upot|ubz|dmscop|getpr|fclr|usite|uspec|fclose|' \
    + 'errmsg|ulat|pshpr|isanrg|tcn|tcx|zmpy|nwordg|poppr|dpcopy|dinv33|fsanrg|dpdamp|query' \
    + ')|mpi*|word*|spack*|info*'
    p=re.compile(neg,re.I)
    sstack=[]
    deptho=-1
    inso=None
    ranktag=["Program", "Subroutine", "Function", "Module","Type"]

    for c in walk(parser.block):
        ins=c[0]
        depth=c[1]
        if(isinstance(ins, classes.Comment)): continue
        if(depth>deptho):         sstack.append(inso)
        while len(sstack)>depth:  sstack.pop()
#        print len(sstack),depth, type(ins), ins.item.span,  ins.item.line, ins.item.label #, ins.item.name
#        deptho=depth
#        inso=ins
#        continue
        frame2x= [x.__class__.__name__+':'+x.name for x in sstack[1:] if x.__class__.__name__ in ranktag]
        frame2='@'.join(frame2x)
#        if(isinstance(ins,classes.Type)):
#            print ins.__class__.__name__ #,ins.item.line
#        continue
        loc1 =  ' in '+frame2+' '
        loc  =  ' '+ffile+(":L%d" % ins.item.span[0])
######
#grep smvxc2 callcaller.dat |egrep -v '(dpzero|defcc|ll|defrc|lgunit|iprint|rx|ival|isum|nglob|cmdopt|rlse|defrr|zsum|pack|dpadd|awrit|dcopy|daxpy|dscal)' 
        if(isinstance(ins, classes.Call)):
            if(p.match(ins.designator)): continue
            if(p.match(x.name)): continue
#            print sstack[len(sstack)-1].name+"->"+ins.designator #ins.item.line
            frame3= [x.name for x in sstack[1:] if x.__class__.__name__ in ranktag]
            print frame3[-1]+"->"+ins.designator+";" #ins.item.line
#
#        if(isinstance(ins, classes.Call)):      print "@cal Subr:"+ins.designator+loc1+loc #ins.item.line
#########
        deptho=depth
        inso=ins
        continue
#########

        if(isinstance(ins, classes.Type)):      print "@dcl Type:"+ins.name+loc1+loc
        if(isinstance(ins, classes.Use)):       print "@use Modu:"+ins.name
        #label=ins.item.label
        if(isinstance(ins, classes.Program)):   print "@def Prog:"+ins.name+loc1+loc #ins.item.name
        if(isinstance(ins, classes.Interface)) :print "@def ifc :"+ins.name+loc1+loc #,ins.a #,ins.a.variable_names
        if(isinstance(ins, classes.Subroutine)):print "@def Subr:"+ins.name+loc1+loc #,ins.a #,ins.a.variable_names
        if(isinstance(ins, classes.Function)):  
            print "@def Func:"+ins.name+loc1+loc,type(ins) #,ins.a
            functions.append(ins.name)
        if(isinstance(ins, classes.Entry)):         print "@def Entr:"+ins.name+loc1+loc
        if(isinstance(ins, classes.Module)):        print "@def Modu:"+ins.name+loc1+loc
        if(isinstance(ins, block_statements.Type)): print "@def Type:"+ins.name+loc1+loc

        deptho=depth
        inso=ins
print '}'
print '</graphviz>'

    #print '---- '+ffile+'  end   -----'
    #print 
sys.exit()

print '--- Looking for function calls=',functions
# Find function calls
for ffile in argset:
    print '---- '+ffile+' start -----'
    reader = FortranFileReader(ffile)
    reader.set_mode(isfree=False,isstrict=False)
    parser=FortranParser(reader,ignore_comments=False)
    parser.parse()

    sstack=[]
    deptho=-1
    inso=None
    ranktag=["Program", "Subroutine", "Function", "Module","Type"]

    for c in walk(parser.block):
        ins=c[0]
        depth=c[1]
        if(isinstance(ins, classes.Comment)): continue
        if(depth>deptho):         sstack.append(inso)
        while len(sstack)>depth:  sstack.pop()

        frame2= [x.__class__.__name__+':'+x.name for x in sstack[1:] if x.__class__.__name__ in ranktag]
        frame2='@'.join(frame2)
        loc1 =  ' in '+frame2+' '
        loc  =  ' '+ffile+(":L%d" % ins.item.span[0])

        line= ins.item.line
        #print 'ppppp',line
        for fun in functions:
            if(fun in line) :
                print "@cal Func:"+fun+loc1+loc,type(ins) #ins.item.line
        deptho=depth
        inso=ins
    print '---- '+ffile+' end   -----'
    print 
