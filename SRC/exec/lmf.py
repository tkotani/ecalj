#!/usr/bin/env python3 
from setcomm import callF,setcommF,getlibF
from mpi4py import MPI
import sys,os
import numpy as np
arglist=' '.join(sys.argv[1:])
scriptpath = os.path.dirname(os.path.realpath(__file__))+'/'
flib = getlibF(scriptpath+'/libecaljF.so') #fortran library

# world
commw = MPI.COMM_WORLD
rankw = commw.Get_rank()
master_mpi= rankw==0

# group1 #only for single core
grp1=[0,1,2,3]
comm1 = setcommF(grp=grp1) #communicator for grp

if(rankw in grp1): 
	callF(flib.setcmdpathc,[scriptpath])  # Set path for ctrl2ctrlp.py at m_setcmdpath
	callF(flib.m_setargsc, [arglist])     # Set args at m_args
#	callF(flib.sopen) 
	callF(flib.lmf,        [comm1]) 
#	callF(flib.sclose) 

if(rankw==0): print('end eeeeee grp1=',grp1)
