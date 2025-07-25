#!/usr/bin/env python3
# Convert ctrl file to ctrlp, stdin to stdout
# example: ctrl can contain
#    %const a=2.4 b=2.3 c=a*b R=c
#  , where >lmf -va=3.8 replace value of a.
# Except %const, such as %if, are not supported. T.K. 2023feb
import sys,re
from math import sin,cos,tan,log,exp,sqrt,pi
instr=sys.stdin.read()
instrl=instr.split('\n')

aaa=''
labels=[]
T,t,F,f=1,1,0,0
pat=r'(\r?\n)|(\r\n?)'
aft=r'\n'
for linei in instrl: # Values from %const section !line cannot contain python keyword
    line=re.sub(pat,aft,linei) #for DOS compatible
    if(len(line)==0): continue
    if(line[0]!='%'): continue
    if(line[0]=='#'): continue
    line=re.sub(r'=\s+','=',line) #remove space after =
    constdata0=(line.split('const')[1]).split('#')[0].split(' ')
    for ix in constdata0:
        if(ix!=''):
            ix=ix.replace('^','**') #math ^ ==> ** .replace('d','e')
            exec(ix)
            labels.append(ix.split('=')[0])
constrep={i:str(eval(i)) for i in labels} #replacement dic
for i in constrep.keys():
    exec('del '+str(i)) #delete const variables after we get numerical constrep 
outfile0=''
for iarg in sys.argv[1:]: # -vfoobar replacemebt by args
    try:
        vin=iarg.split('-v')[1]
        label,val=vin.split('=')
    except:
        continue
    try:
        outfile0=outfile0+' '+str(label)+' '+str(eval(val))+' ! -vfoobar\n'
    except:
        outfile0=outfile0+' '+str(label)+' '+str(val)+' ! -vfoobar char\n'
    constrep[label]=str(val)
#print(sys.argv[1:])
#print(outfile0.split('\n'))
outfile0=outfile0[0:-1] #remove final \n
#print(outfile0.split('\n'))

midfile=instr
for i,irep in constrep.items():
    ix='{'+i+'}'
    midfile=midfile.replace(ix,irep) #print('ix rep=',ix,irep)

### Pure math section. # we replaced {foobar} with numerical values.
outfile=''
for ilinex in midfile.split('\n'): #line by line, for pure mathematical operations.
    iline=re.sub(pat,aft,ilinex) #for DOS compatible
    if(len(iline)==0): continue
    if(iline[0]=='%'): continue 
    iii= iline.split('#')[0].split('!')[0].split('%')[0]  #print('iii ',iii)
    iii= re.sub(r'^\s+$','',iii)
    if(len(iii)==0): continue    #print('line=',iii,'=========')
    mmm=re.split('([=, ])',iii) #(...) means separators in lists.
    if('SYMGR' in mmm[0]): #insereted at 2023-5-9 for NiS
        iii=re.sub(r'\(','( ',iii)
        iii=re.sub(r'\)',' )',iii)
        mmm=re.split('([=, ])',iii) #(...) means separators in lists.
    nnn=[]
    for ix in mmm:
        eout=ix
        #print(eout)
        try:
            eout=eval(ix) # math
        except:
            pass
        nnn.append(str(eout).replace(',',' '))
        #print('input:',ix,' output:',eout)
        #print('input :'+''.join(mmm))
        #print('output:'+''.join(nnn))
    outfile=outfile+''.join(nnn)+'\n'
outfile=re.sub(r'\t',' ',outfile)
lll=''
init=False
ix=0
for line in outfile.split('\n'):
     if(len(line)==0): continue
     if(line[0]!=' '):
         ix=ix+1
         ladd=line
         if(len(ladd)!=0): lll= lll+'\n'+ladd
         init=True
     elif(init==True):
         lll=lll+line.rstrip(' ')
lll=lll[1:] #remove initial \n
lll=re.sub(r'=\s+','=',lll) 
lmax=0
for line in lll.split('\n'):
    #print(len(line),line)
    if(len(line)>lmax): lmax=len(line)
#print(len(outfile0.split('\n')))
#print('lxxxxxxx',lll)
catok=''
nbas=0
nspec=0
for iline in lll.split('\n'):
    line=iline
    cat=iline.split(' ')[0]
    if cat=='ITER': line=line.replace(',',' ')
    
    line=[i for i in line.split(' ') if i!='']
    tok=''
    tokk=[]
    cat=line[0].split(' ')[0]
    id=0
    idx=''
    line.append('EOL') #for re.match satisfied at the end 
    for ix in line[1:]:
        i=ix
        if cat=='SITE':
            tokx=i
            if(tokx[0:4]=='ATOM'):
                id=id+1
                nbas=nbas+1
        if cat=='SPEC':
            tokx=i
            if(tokx[0:4]=='ATOM'):
                id=id+1
                nspec=nspec+1
        if re.match('[a-zA-Z]',i[0]):
            ic='_'
            #if(cat=='SYMGRP'): ic=' '
            dat=cat+ic+tok
            dat=dat.replace('=',idx+' ',1)
            #if(id>0  and (not 'ATOM' in dat) ):
            #    dat=dat.replace('_','_')
            tokk.append(dat) #tok write
            idx=''
            if(id>0):
                idx='@'+str(id)
            tok=i
        else:
            tok=tok +' '+ i
    #print('tokk=',tokk) 
    for itk in tokk:
        if(itk[-1]!='_' and itk[0:6]!='SYMGRP'):  #Skip SYMGRP. We don't need process SYMGRP SYMGRPAF
            #print('itk=',itk)
            catok=catok+itk+'\n'
catok=catok+'STRUC_NSPEC '+str(nspec)+'\n' # less prior if STRUC_NSPEC already in ctrl
catok=catok+'STRUC_NBAS '+str(nbas)+'\n'   # less prior
llx=len(outfile0.split('\n'))+len(lll.split('\n'))

for iline in lll.split('\n'):
    if('SYMGRP' in iline): print(iline)
#print(catok)
for iline in catok.split('\n'):
    ilinex=iline
    if(re.search('_ATOM@',iline)): ilinex=re.sub(r' 0',' F',iline)
    print(ilinex)
