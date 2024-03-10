#!/usr/bin/env python3 
from setcomm import callF,setcommF,getlibF
from mpi4py import MPI
import sys,os
#import numpy as np
arglist=' '.join(sys.argv[1:])
scriptpath = os.path.dirname(os.path.realpath(__file__))+'/'
ecaljF = getlibF(scriptpath+'/libecaljF.so')

# world
commw = MPI.COMM_WORLD
rankw = commw.Get_rank()

#initialization
callF(ecaljF.setcmdpathc,[scriptpath])  # Set path for ctrl2ctrlp.py at m_setcmdpath
callF(ecaljF.m_setargsc,[arglist])      # Set args at m_args

# group1 #only for single core
grp1=[0]
comm1 = setcommF(grp=grp1)
if(rankw in grp1): 
	callF(ecaljF.lmfa,[comm1]) 
