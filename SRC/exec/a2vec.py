#!/usr/bin/env python3
#pure math conversion
import sys,re 
from math import *
mmm=re.split('[=, ]',sys.stdin.read())
mmm=[i for i in mmm if len(i)!=0]
aaa='\n'
nx=0
for ix in mmm:
    eout=ix
    try:
        eout=eval(ix) # math
    except:
        break
    nx=nx+1
    aaa=aaa+' '+str(eout)
print(nx, aaa)
    
