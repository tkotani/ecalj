#!/usr/bin/env python3 
#>mpirun -np 4 pysample
from setcomm import callF,setcommF,getlibF
from mpi4py import MPI
import sys,os
import numpy as np
arglist=' '.join(sys.argv[1:])
scriptpath = os.path.dirname(os.path.realpath(__file__))+'/'
ecaljF = getlibF(scriptpath+'/libecaljF.so')

# world
commw = MPI.COMM_WORLD
rankw = commw.Get_rank()
master_mpi= rankw==0

# group1 #only for single core
grp1=[0,1,2,3]
commF1 = setcommF(grp=grp1)
if(master_mpi): print('gpr1=',grp1)

if(rankw in grp1): 
	callF(ecaljF.setcmdpathc,[scriptpath])  # Set path for ctrl2ctrlp.py at m_setcmdpath
	callF(ecaljF.m_setargsc,[arglist])      # Set args at m_args
#	callF(ecaljF.sopen) 
	callF(ecaljF.lmf,[commF1]) 
#	callF(ecaljF.sclose) 

if(rankw==0): print('end eeeeee grp1=',grp1)

