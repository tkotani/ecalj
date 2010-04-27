#!/usr/bin/env python
# Program unit analyzer.
def addfsdist(fsdict,key,data):
    #print 'vvv',key,data
    if(key in fsdict):
        fsdict[key]= fsdict[key]+' '+data
    else:
        fsdict[key]=data
    return


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
moddep  = open("Make.mod.dependency",'wt')
mover   = open("MoveUnUsed",'wt')
notused = open("notused",'wt')

nargv = len(sys.argv) -1
if(nargv==0 or sys.argv[1]=='--help'):
    print ' === A Fortran code analyzer based on f2py/fparser r41 http://code.google.com/p/f2py/ ==='
    print ' usage: f_calltree.py  file1.F file2.F ... '
    print '        files can be specified by wild cards as "subs/*.F fp/*.F"'
    print '  Standard out put   : text data base. Easy to understand. Use grep to read it'
    print '  Make.mod.dependency: module dependency, which can be inlucded in Makefile'
    print '  callcaller.dotdata : in graphviz format. '
    print '                       sort callcaller.dotdata|uniq >dot.dat'
    print '                       Then add "digraph G{" at top of dot.dat, and "}" at end of dot.dat.'  
    print '                       dot -Tps dot.dat -o dot.ps'
    print ' NOTE: T.Kotani slightly modified f2py/fparser r41 so that it works standalone. No numpy.'
    print '       This works with python 2.6 or so.'
    print '       T.Kotani highly acknowledge to Dr.pearu.peterson, who mainly develop the fparser.'
    print ' --> An unique feature of this f_calltree.py is "Easy custmization". Look into this code.'
    print ' --> This code is a litle specialized for ecal package, but it should be not so difficult to fit it to your codes.'
    sys.exit()

argset= sys.argv[1:]

#print argset
#sys.exit()
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
mods=[]
srcfiles=set(argset)
ranktag=["Program", "Subroutine", "Function", "Module","Type","Interface"]

moddict={}
fsdict={}
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
        if(isinstance(ins, classes.Type)):      print "@dcl Type:"+ins.name+loc1+loc + ' ---> '+ins.item.line
#        if(isinstance(ins, classes.Use)):       print "@use Modu:"+ins.name
         #label=ins.item.label
        if(isinstance(ins, classes.Program)):   print "@def Prog:"+ins.name+loc1+loc + ' ---> '+ins.item.line
        if(isinstance(ins, classes.Interface)) :print "@def ifc :"+ins.name+loc1+loc + ' ---> '+ins.item.line #,ins.a #,ins.a.variable_names
        if(isinstance(ins, classes.Subroutine)):
            print "@def Subr:"+ins.name+loc1+loc + ' ---> '+ins.item.line #,ins.a #,ins.a.variable_names
            subs.append(ins.name)
            addfsdist(fsdict,ins.name,ffile)
        if(isinstance(ins, classes.Function)):  
            print "@def Func:"+ins.name+loc1+loc + ' ---> '+ins.item.line #,type(ins) #,ins.a
            functions.append(ins.name)
            addfsdist(fsdict,ins.name,ffile)
# now all entry is classified to function (not correct but safer treatment)
        if(isinstance(ins, classes.Entry)):
            mother=sstack[-1]
            if(mother.__class__.__name__=='Function')  : 
                aaa="@def Entr_Func:"
                functions.append(ins.name)
                addfsdist(fsdict,ins.name,ffile)
            if(mother.__class__.__name__=='Subroutine'): 
                aaa="@def Enty_Subr:"
                subs.append(ins.name)
                addfsdist(fsdict,ins.name,ffile)
            print aaa+ins.name+loc1+loc + ' ---> '+ins.item.line
        if(isinstance(ins, classes.Module)):        
            mods.append(ins.name)
            print "@def Modu:"+ins.name+loc1+loc + ' ---> '+ins.item.line
            addfsdist(fsdict,ins.name,ffile)
            ffileo = re.sub('subs/','$(subs_obj_path)/',ffile)
            ffileo = re.sub('fp/'  ,'$(fp_obj_path)/',ffileo)
            ffileo = re.sub('slatsm/','$(slatsm_obj_path)/',ffileo)
            ffileo = re.sub('.F','.o',ffileo)
            moddep.write("# $(moddir)/"+ins.name+".mod :"+ffile+"\n")
            moddict[ins.name]=ffileo
        if(isinstance(ins, block_statements.Type)): print "@def Type:"+ins.name+loc1+loc + ' ---> '+ins.item.line
        deptho=depth
        inso=ins

############################
print >>sys.stderr, 'moddict=',moddict 

targetf = "(\W"+'\s*\(|\W'.join(functions)+"\s*\()"
if('|'.join(functions)==''): targetf='xx xx xx xx xx' #this is when targetf is empty. you know better procedure?
pf=re.compile(targetf,re.I)

targets = "(\W"+'\s*\(|\W'.join(subs)+"\s*\()"
if('|'.join(subs)==''): targets='xx xx xx xx xx'
ps=re.compile(targets,re.I)

targetsf = "(\W"+'\s*\(|\W'.join(subs)+'\s*\(|\W'.join(functions)+"\s*\()"
if('|'.join(subs)+'|'.join(functions)==''): targetsf='xx xx xx xx xx'
psf=re.compile(targetsf,re.I)

targetsfm = "(\W"+'\s*\(|\W'.join(mods)+"\s*\()"
psm=re.compile(targetsfm,re.I)


targets0 = "(\W"+'\s*\(|\W'.join(subs)+"\s*\()"
if('|'.join(subs)==''): targets0='xx xx xx xx xx'
ps0=re.compile(targets0,re.I)

print >> sys.stderr, 'pf=',targetf


##################
sets=set(subs)
### check names of subs are not duplicated with names of functions.
for i in functions:
    if i in sets:
        print i,' is in subs'

#setf=set(functions)
#setsf= sets.union(setf)
#print >> sys.stderr,'set sub number=',len(sets)
#print >> sys.stderr,'set fun number=',len(setf)
#print >> sys.stderr,'set sub+function number=',len(setsf),'\n',setsf
#sys.exit()

usedsf=set([])
##############
#>  f_calltree.py subs/m_struc_*.F subs/struc_main.F subs/struc_sub.F
#ttt="(\Wstruc_eval_io_r8\s*\(|\Wstruc_eval_io_r8v\s*\(|\Wstruc_eval_io_i8\s*\(|\Wstruc_eval_io_i8v\s*\(|\Wstruc_strtok\s*\(|\Wstruc_eval_io_r8_realbody\s*\(|\Wstruc_eval_io_i8_realbody\s*\(|\Wstruc_checkclass\s*\(|\Wstruc_spackv_iv\s*\(|\Wstruc_spackv_r8v\s*\(|\Wspackv\s*\(|\Wspacks\s*\(|\Wsp2cls\s*\(|\Wshstru\s*\(|\Wstruc_packupack_val1\s*\(|\Wstruc_uarray_io\s*\(|\Wstruc_ubz_io\s*\(|\Wstruc_uctrl_io\s*\(|\Wstruc_ugw_io\s*\(|\Wstruc_uham_io\s*\(|\Wstruc_ulat_io\s*\(|\Wstruc_umix_io\s*\(|\Wstruc_umove_io\s*\(|\Wstruc_uoptic_io\s*\(|\Wstruc_uordn_io\s*\(|\Wstruc_upot_io\s*\(|\Wstruc_usite_io\s*\(|\Wstruc_uspec_io\s*\(|\Wstruc_ustr_io\s*\(|\Wstruc_utb_io\s*\(|\Wuarray_init\s*\(|\Wuarray_show\s*\(|\Wuarray\s*\(|\Wubz_init\s*\(|\Wubz_show\s*\(|\Wubz\s*\(|\Wuctrl_init\s*\(|\Wuctrl_show\s*\(|\Wuctrl\s*\(|\Wugw_init\s*\(|\Wugw_show\s*\(|\Wugw\s*\(|\Wuham_init\s*\(|\Wuham_show\s*\(|\Wuham\s*\(|\Wulat_init\s*\(|\Wulat_show\s*\(|\Wulat\s*\(|\Wumix_init\s*\(|\Wumix_show\s*\(|\Wumix\s*\(|\Wumove_init\s*\(|\Wumove_show\s*\(|\Wumove\s*\(|\Wuoptic_init\s*\(|\Wuoptic_show\s*\(|\Wuoptic\s*\(|\Wuordn_init\s*\(|\Wuordn_show\s*\(|\Wuordn\s*\(|\Wupot_init\s*\(|\Wupot_show\s*\(|\Wupot\s*\(|\Wusite_init\s*\(|\Wusite_show\s*\(|\Wusite\s*\(|\Wuspec_init\s*\(|\Wuspec_show\s*\(|\Wuspec\s*\(|\Wustr_init\s*\(|\Wustr_show\s*\(|\Wustr\s*\(|\Wutb_init\s*\(|\Wutb_show\s*\(|\Wutbuarray_size\s*\(|\Wubz_size\s*\(|\Wuctrl_size\s*\(|\Wugw_size\s*\(|\Wuham_size\s*\(|\Wulat_size\s*\(|\Wumix_size\s*\(|\Wumove_size\s*\(|\Wuoptic_size\s*\(|\Wuordn_size\s*\(|\Wupot_size\s*\(|\Wusite_size\s*\(|\Wuspec_size\s*\(|\Wustr_size\s*\(|\Wutb_size\s*\()"
#ps=re.compile(ttt,re.I)
#psf=ps
#pf=ps
#ps0=ps
#print targetsf
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
    ffileo = re.sub('subs/','$(subs_obj_path)/',ffile)
    ffileo = re.sub('fp/'  ,'$(fp_obj_path)/',ffileo)
    ffileo = re.sub('slatsm/','$(slatsm_obj_path)/',ffileo)
    ffileo = re.sub('.F','.o',ffileo)
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
                print "@cal Subr:"+ins.designator+loc1+loc+ ' ---> '+ins.item.line
                frame3= [x.name for x in sstack[1:] if x.__class__.__name__ in ranktag]
                mother= frame3[-1]
                child = ins.designator
                dotdata.write( (mother+"->"+child+";\n").lower() )
                #setsf.discard(child)
                #srcfiles.discard(ffile)
                usedsf.add(child)
        if(isinstance(ins, classes.Use)):
            print "@use Modu:"+ins.name+loc1+loc + ' ---> '+ins.item.line
            usedsf.add(ins.name)
## For test, you may need to comment out a following line because it can cause key error 
## if all used modules are not contained in input files.
            if(moddict[ins.name] != ffileo):
                #print 'checkwrite qqqqq', moddict[ins.name],ffileo
                modd.append(moddict[ins.name])
        deptho=depth
        inso=ins
    #print modd
    moddx   = ' '.join(list(set(modd)))
    #print moddx
    #print 'xxxxxxx',ffile
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
#            print 'aaa',a.group().lower(),'bbb'
            frame3= [x.name for x in sstack[1:] if x.__class__.__name__ in ranktag]
#            print frame3
            mother= frame3[-1].lower()
            b=re.match('(\w|_)*', a.group().lower()[1:])
            child=b.group().lower()
            #setsf.discard(child)
            #srcfiles.discard(ffile)
            usedsf.add(child)
            if(mother != child):  #mother = child means a case where a function name is in the function.
                #lll=(mother+"->"+child+";").lower()
                #if(pn.search(lll)): continue
                print "@cal Func:"+child+loc1+loc +' ---> '+ins.item.line  #,type(ins)
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


###############################################################
print >>sys.stderr, '------------ dictionary ------------------'
print >>sys.stderr, fsdict

print >>sys.stderr, '------------ Used sub/fun ----------------'
print >>sys.stderr, usedsf

print >>sys.stderr, '------------ Used files ------------------'
using=[]
for i in usedsf:
    if i in fsdict:
        fff = fsdict[i]
        fff = re.split(' ',fff)
        print >>sys.stderr,fff
        using=using+ fff
print >>sys.stderr, set(using)

print >>sys.stderr, '------------ UnUsed files ------------\n'
setun=srcfiles.difference(set(using))
print >>sys.stderr, setun
for i in setun:
    mover.write("mv "+i+" SRCbank/"+i+"\n")
    notused.write(i+" ")
