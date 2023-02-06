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
for line in instrl: # Values from %const section !line cannot contain python keyword
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
#print(constrep)
outfile=''
for iarg in sys.argv[1:]: # -vfoobar replacemebt by args
    try:
        vin=iarg.split('-v')[1]
    except:
        continue
    label,val=vin.split('=')
    try:
        outfile=outfile+' '+str(label)+' '+str(eval(val))+'\n'
    except:
        outfile=outfile+' '+str(label)+' '+str(val)+'\n'
    constrep[label]=str(val)
outfile=outfile+'=== end of -vfoobar ===\n'
#print(constrep)

midfile=instr
for i,irep in constrep.items():
    ix='{'+i+'}'
    midfile=midfile.replace(ix,irep) #print('ix rep=',ix,irep)

#Pure math section. # we replaced {foobar} with numerical values.
for iline in midfile.split('\n'): #line by line, for pure mathematical operations.
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
        nnn.append(str(eout))
    #print('input :'+''.join(mmm))  #print('output:'+''.join(nnn))
    outfile=outfile+''.join(nnn)+'\n'
ll=len(outfile.split('\n'))
print(ll,' !line number')
print(outfile)
