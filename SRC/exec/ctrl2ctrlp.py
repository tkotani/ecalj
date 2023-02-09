#!/usr/bin/env python3
# Convert ctrl file to ctrlp, stdin to stdout
# example: ctrl can contain
#     %const a=2.4 b=2.3 c=a*b R=c
#  , where >lmf-MPIK -va=3.8 replace value of a.
# Except %const, such as %if, are not supported. T.K. 2023feb
import sys,re
from math import *
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
    #print(line)
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
    except:
        continue
    label,val=vin.split('=')
    try:
        outfile0=outfile0+' '+str(label)+' '+str(eval(val))+' ! -vfoobar\n'
    except:
        outfile0=outfile0+' '+str(label)+' '+str(val)+' ! -vfoobar char\n'
    constrep[label]=str(val)
#print(sys.argv[1:])
#print(outfile0.split('\n'))
outfile0=outfile0[0:-1] #remove final \n
#print(outfile0.split('\n'))
#sys.exit()

midfile=instr
for i,irep in constrep.items():
    ix='{'+i+'}'
    midfile=midfile.replace(ix,irep) #print('ix rep=',ix,irep)

#Pure math section. # we replaced {foobar} with numerical values.
outfile=''
for ilinex in midfile.split('\n'): #line by line, for pure mathematical operations.
    iline=re.sub(pat,aft,ilinex) #for DOS compatible
    if(len(iline)==0): continue
    if(iline[0]=='%'): continue 
    iii= iline.split('#')[0].split('!')[0].split('%')[0]  #print('iii ',iii)
    iii= re.sub('^\s+$','',iii)
    if(len(iii)==0): continue    #print('line=',iii,'=========')
    mmm=re.split('([=, ])',iii) #(...) means separators in lists.
    nnn=[]
    for ix in mmm:
        eout=ix
        try:
            eout=eval(ix) # math
        except:
            pass
        nnn.append(str(eout)) #.replace(',',' ').replace('=','= '))
    #print('input :'+''.join(mmm))  #print('output:'+''.join(nnn))
    outfile=outfile+''.join(nnn)+'\n'
#llx=len(outfile.split('\n'))
#print(outfile.split('\n'))
#print('outfile=###'+outfile+'###')

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
#print(lll)
lll=lll[1:] #remove initial \n

lmax=0
for line in lll.split('\n'):
    #print(len(line),line)
    if(len(line)>lmax): lmax=len(line)
#print(len(outfile0.split('\n')))
#print(lll.split('\n'))
llx=len(outfile0.split('\n'))+len(lll.split('\n'))
print(llx,lmax,'# of line, # of reclen Category per line\n'+outfile0+'\n'+lll)

#print('=====following is for human readable ===========')
#print(outfile)
