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

dotdata = open("callcaller.dotdata",'wt')
moddep = open("Make.mod.dependency",'wt')

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

rr=re.compile("\'.*?\'|\".*?\"")
functions=[] # save defined functions in array.
subs=[]
ranktag=["Program", "Subroutine", "Function", "Module","Type","Interface"]

moddic={}
### find subroutine and function definitions ############################
for ffile in argset:
#    print '---- '+ffile+' start -----'
    reader = FortranFileReader(ffile)
    reader.set_mode(isfree=False,isstrict=False)
    parser=FortranParser(reader,ignore_comments=False)
    parser.parse()
    sstack=[]
    deptho=-1
    inso=None
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
#        if(isinstance(ins, classes.Call)):      print "@cal Subr:"+ins.designator+loc1+loc #ins.item.line
        if(isinstance(ins, classes.Type)):      print "@dcl Type:"+ins.name+loc1+loc
#        if(isinstance(ins, classes.Use)):       print "@use Modu:"+ins.name
         #label=ins.item.label
        if(isinstance(ins, classes.Program)):   print "@def Prog:"+ins.name+loc1+loc #ins.item.name
        if(isinstance(ins, classes.Interface)) :print "@def ifc :"+ins.name+loc1+loc #,ins.a #,ins.a.variable_names
        if(isinstance(ins, classes.Subroutine)):
            print "@def Subr:"+ins.name+loc1+loc #,ins.a #,ins.a.variable_names
            subs.append(ins.name)
        if(isinstance(ins, classes.Function)):  
            print "@def Func:"+ins.name+loc1+loc #,type(ins) #,ins.a
            functions.append(ins.name)
# now all entry is classified to function (not correct but safer treatment)
        if(isinstance(ins, classes.Entry)):
            mother=sstack[-1]
            if(mother.__class__.__name__=='Function')  : 
                aaa="@def Entr_Func:"
                functions.append(ins.name)
            if(mother.__class__.__name__=='Subroutine'): 
                aaa="@def Enty_Subr:"
                subs.append(ins.name)
            print aaa+ins.name+loc1+loc
        if(isinstance(ins, classes.Module)):        
            print "@def Modu:"+ins.name+loc1+loc
            ffileo = re.sub('subs/','$(subs_obj_path)/',ffile)
            ffileo = re.sub('fp/'  ,'$(fp_obj_path)/',ffileo)
            ffileo = re.sub('slatsm/','$(slatsm_obj_path)/',ffileo)
            ffileo = re.sub('.F','.o',ffileo)
            moddep.write("# $(moddir)/"+ins.name+".mod :"+ffile+"\n")
            moddic[ins.name]=ffileo
        if(isinstance(ins, block_statements.Type)): print "@def Type:"+ins.name+loc1+loc
        deptho=depth
        inso=ins

############################
print 'moddic=',moddic


targetf = "("+'\Z|'.join(functions)+"\Z)"
if('|'.join(functions)==''): targetf='xx xx xx xx xx' #this is when targetf is empty. you know better procedure?
pf=re.compile(targetf,re.I)

targets = "("+'\Z|'.join(subs)+"\Z)"
if('|'.join(subs)==''): targets='xx xx xx xx xx'
ps=re.compile(targets,re.I)

targetsf = "("+'\Z|'.join(subs)+'\Z|'.join(functions)+"\Z)"
if('|'.join(subs)+'|'.join(functions)==''): targetsf='xx xx xx xx xx'
psf=re.compile(targetsf,re.I)

targets0 = "(\W"+'\W|\W'.join(subs)+"\W)"
if('|'.join(subs)==''): targets0='xx xx xx xx xx'
ps0=re.compile(targets0,re.I)
#print targets0
#sys.exit()


### Find subroutine calls ########################################
for ffile in argset:
#    print '---- '+ffile+' start -----'
    modd=[]
    reader = FortranFileReader(ffile)
    reader.set_mode(isfree=False,isstrict=False)
    parser=FortranParser(reader,ignore_comments=False)
    parser.parse()
    #parser.analyze()
    sstack=[]
    deptho=-1
    inso=None
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
        line= rr.sub("''",ins.item.line)
        #print line
        if(ps0.search(line)) :
            if(isinstance(ins, classes.Call)):      
                print "@cal Subr:"+ins.designator+loc1+loc #ins.item.line
                frame3= [x.name for x in sstack[1:] if x.__class__.__name__ in ranktag]
                mother= frame3[-1]
                child = ins.designator
                dotdata.write( (mother+"->"+child+";\n").lower() )
        if(isinstance(ins, classes.Use)):       
            print "@use Modu:"+ins.name+loc1+loc
            modd.append(moddic[ins.name])
            #modd.append(' $(moddir)/'+ins.name+'.mod')
            
        deptho=depth
        inso=ins
    #print modd
    moddx   = ' '.join(list(set(modd)))
    #print moddx
    #print 'xxxxxxx',ffile
    ffileo = re.sub('subs/','$(subs_obj_path)/',ffile)
    ffileo = re.sub('fp/'  ,'$(fp_obj_path)/',ffileo)
    ffileo = re.sub('slatsm/','$(slatsm_obj_path)/',ffileo)
    ffileo = re.sub('.F','.o',ffileo)
    if(moddx != '') :
        moddep.write( ffileo+' : '+ffile+' '+moddx+'\n')

#    print '---- '+ffile+'  end   -----'
#    print 

#print '--- Looking for function calls ---'
#print functions
#targetf = "("+'\Z|'.join(functions)+"\Z)"
#if('|'.join(functions)==''): targetf='xx xx xx xx xx' #dummy when targetf is empty. Do you know better method?
#pf=re.compile(targetf,re.I)


#### Find function calls ########################################
for ffile in argset:
#    print '---- '+ffile+' start -----'
    reader = FortranFileReader(ffile)
    reader.set_mode(isfree=False,isstrict=False)
    parser=FortranParser(reader,ignore_comments=False)
    parser.parse()

    sstack=[]
    deptho=-1
    inso=None
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
        line= rr.sub("''",line)
#        print 'ppppp',depth,line,[x.name for x in sstack[1:]]
#        print 'pppp2',depth,rr.sub("''",line)
#        deptho=depth
#        inso=ins
#        continue

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
            if(mother != child):  #mother = child means a case where a function name is in the function.
                #lll=(mother+"->"+child+";").lower()
                #if(pn.search(lll)): continue
                print "@cal Func:"+child+loc1+loc #,type(ins)
                dotdata.write( (mother+"->"+child+";\n").lower() )
        deptho=depth
        inso=ins
#print 'ppppp',line
#         for fun in functions:
#             if(fun in line) :
#                 print "@cal Func:"+fun+loc1+loc,type(ins) #ins.item.line
#         deptho=depth
#         inso=ins
#    print '---- '+ffile+' end   -----'
#    print 
