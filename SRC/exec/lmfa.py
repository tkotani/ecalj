#!/usr/bin/env python3 
#>mpirun -np 4 pysample
from setcomm import callF,setcommF,getlib
from mpi4py import MPI
import sys,os
import numpy as np
arglist=' '.join(sys.argv[1:])
scriptpath = os.path.dirname(os.path.realpath(__file__))+'/'
ecaljF = getlib(scriptpath+'/libecaljF.so')

# world
commw = MPI.COMM_WORLD
rankw = commw.Get_rank()
# group1 #only for single core
grp1=[0]
commF1 = setcommF(grp=grp1)

if(rankw in grp1): 
	callF(ecaljF.setcmdpathc,[scriptpath])  # Set path for ctrl2ctrlp.py at m_setcmdpath
	callF(ecaljF.m_setargsc,[arglist])      # Set args at m_args
	callF(ecaljF.lmfa,[commF1]) 
	#print(f'rrrrrrrrr rank={rankw} is doing something')
else:
        #print(f'rrrrrrrrr rank={rankw} is not doing anything')
        pass
