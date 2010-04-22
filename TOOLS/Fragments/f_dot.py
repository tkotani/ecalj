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
subs=[] # save defined functions in array.
#print '<graphviz>'
print 'digraph G {'

# find defined functions and subroutines.
neg='\Wll\W|defrc|lgunit|iprint|rx|ival|isum|tkadd|nswadd|toksw|lgors|\Wnm\W|' \
    + 'nglob|cmdopt|zsum|pack|upack|dpadd|awrit|defd*|defc*|defi*' \
    + 'r8tos8|dvset|dpscop|dmpy|fexit|sp2cls|redfrr|upot|ubz|ugw|utb|ustr|iget|dget|dmscop|getpr|fclr|usite|uspec|fclose|' \
    + 'errmsg|ulat|pshpr|isanrg|tcn|tcx|zmpy|nwordg|poppr|dpcopy|dinv33|fsanrg|dpdamp|query' \
    + '|mpi*|word*|spack*|info*|struc_*'

for ffile in argset:
#    print '---- '+ffile+' start -----'
    reader = FortranFileReader(ffile)
    reader.set_mode(isfree=False,isstrict=False)
    parser=FortranParser(reader,ignore_comments=False)
    parser.parse()
    pn=re.compile(neg,re.I)
    sstack=[]
    deptho=-1
    inso=None
    ranktag=["Subroutine", "Function"]
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
#        if(isinstance(ins, classes.Type)):      print "@dcl Type:"+ins.name+loc1+loc
#        if(isinstance(ins, classes.Use)):       print "@use Modu:"+ins.name
         #label=ins.item.label
#        if(isinstance(ins, classes.Program)):   print "@def Prog:"+ins.name+loc1+loc #ins.item.name
#        if(isinstance(ins, classes.Interface)) :print "@def ifc :"+ins.name+loc1+loc #,ins.a #,ins.a.variable_names
        if(isinstance(ins, classes.Subroutine)):
            #print "@def Subr:"+ins.name+loc1+loc #,ins.a #,ins.a.variable_names
            subs.append(ins.name)
        if(isinstance(ins, classes.Function)):  
            #print "@def Func:"+ins.name+loc1+loc,type(ins) #,ins.a
            functions.append(ins.name)
# now all entry is classified to function (not correct but safer treatment)
        if(isinstance(ins, classes.Entry)):
            #print "@def Entr:"+ins.name+loc1+loc
            functions.append(ins.name)
#        if(isinstance(ins, classes.Module)):        print "@def Modu:"+ins.name+loc1+loc
#        if(isinstance(ins, block_statements.Type)): print "@def Type:"+ins.name+loc1+loc
        deptho=depth
        inso=ins

targetf = "("+'\Z|'.join(functions)+"\Z)"
if('|'.join(functions)==''): targetf='xx xx xx xx xx'
pf=re.compile(targetf,re.I)

targets = "("+'\Z|'.join(subs)+"\Z)"
if('|'.join(subs)==''): targets='xx xx xx xx xx'
ps=re.compile(targets,re.I)

targetsf = "("+'\Z|'.join(subs)+'\Z|'.join(functions)+"\Z)"
if('|'.join(subs)+'|'.join(functions)==''): targetsf='xx xx xx xx xx'
psf=re.compile(targetsf,re.I)

#print targets
#print '}'
#print 'targetf=',functions

#####################################################
for ffile in argset:
    #print '---- '+ffile+' start -----'
    reader = FortranFileReader(ffile)
    reader.set_mode(isfree=False,isstrict=False)
    parser=FortranParser(reader,ignore_comments=False)
    parser.parse()
    #parser.analyze()

    #p=re.compile(neg,re.I)
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
        frame2x= [x.__class__.__name__+':'+x.name for x in sstack[1:] if x.__class__.__name__ in ranktag]
        frame2='@'.join(frame2x)
        loc1 =  ' in '+frame2+' '
        loc  =  ' '+ffile+(":L%d" % ins.item.span[0])
        if(isinstance(ins, classes.Call)):
            frame3= [x.name for x in sstack[1:] if x.__class__.__name__ in ranktag]
            mother= frame3[-1]
            child = ins.designator
            if(psf.match(child)):
            #print 'pf.mother=',pf.match(mother)
            #print 'pf.child =',pf.match(child)
                lll=(mother+"->"+child+";").lower()
                if(pn.search(lll)): continue
                print lll #ins.item.line
        deptho=depth
        inso=ins

# Find function calls
sys.stderr.write('--- Looking for calls for='+'|'.join(functions))
#print 
#print 'xxxxxxxxxxxxxxxxxxxxxxxx'
#print targetf
rr=re.compile("\'.*?\'")

for ffile in argset:
    #print '---- '+ffile+' start -----'
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
        if(depth>deptho):
            #print inso
            sstack.append(inso)
        while len(sstack)>depth:  sstack.pop()

        frame2= [x.__class__.__name__+':'+x.name for x in sstack[1:] if x.__class__.__name__ in ranktag]
        frame2='@'.join(frame2)
        loc1 =  ' in '+frame2+' '
        loc  =  ' '+ffile+(":L%d" % ins.item.span[0])
        line= ins.item.line
        line= rr.sub("''",line)
#
#        print 'ppppp',depth,line,[x.name for x in sstack[1:]]
#        print 'pppp2',depth,rr.sub("''",line)
#        deptho=depth
#        inso=ins
#        continue

###################################3
        if(isinstance(ins, classes.Function) or isinstance(ins,classes.EndStatement) ):  
            deptho=depth
            inso=ins
            continue
        a=pf.search(line)
        if(a):
            #print 'aaa',a.group().lower(),'bbb'
            frame3= [x.name for x in sstack[1:] if x.__class__.__name__ in ranktag]
            #print frame3
            mother= frame3[-1].lower()
            child=a.group().lower()
            if(mother != child): 
                lll=(mother+"->"+child+";").lower()
                if(pn.search(lll)): continue
                print lll #ins.item.line
        deptho=depth
        inso=ins
#    print '---- '+ffile+' end   -----'
#    print 
print '}'
#print '</graphviz>'
sys.exit()
