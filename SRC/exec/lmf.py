#!/usr/bin/env python3 
from setcomm import callF,setcommF,getlibF
from mpi4py import MPI
import sys,os,glob
#import numpy as np
arglist=' '.join(sys.argv[1:])
scriptpath = os.path.dirname(os.path.realpath(__file__))+'/'

# MPIworld
commw = MPI.COMM_WORLD
sizew = commw.Get_size()
rankw = commw.Get_rank()
master_mpi= rankw==0

# lmfa   #files=glob.glob(os.getcwd()+'/atmpnu.*'); #if(len(files)==0):
group=[0] #only rank=0 is used
stdout='llmfa'
if(master_mpi): print('=== Run lmfa by ranks=',group, ' See stdout=',stdout)
flib = getlibF(scriptpath+'/libecaljF.so',prt=master_mpi) #load dynamic library
comm = setcommF(grp=group)
if(rankw in group): 
    callF(flib.setcmdpathc,[scriptpath,master_mpi])  # Set path for ctrl2ctrlp.py at m_setcmdpath
    callF(flib.m_setargsc, [arglist,master_mpi])     # Set args at m_args
    callF(flib.sopen,[stdout]) #standard output
    callF(flib.lmfa,[comm])
    callF(flib.sclose)
flib.dlclose(flib._handle) #close library to delete all allocations and for initialization.
if(master_mpi): print('=== end of lmfa ===')

# lmf
group=[i for i in range(sizew)] #used ranks
stdout='llmf'
if(master_mpi): print('=== Run lmf by ranks=',group, ' See stdout=',stdout)
flib = getlibF(scriptpath+'/libecaljF.so',prt=master_mpi) #load dynamic library
comm = setcommF(grp=group) #communicator for group
if(rankw in group): 
    callF(flib.setcmdpathc,[scriptpath,master_mpi])  # Set path for ctrl2ctrlp.py at m_setcmdpath
    callF(flib.m_setargsc, [arglist,master_mpi])     # Set args at m_args
    callF(flib.sopen,[stdout]) #standard output
    callF(flib.lmf,  [comm]) 
    callF(flib.sclose) 
flib.dlclose(flib._handle) #close library
if(master_mpi): print('=== end of lmf ===')
