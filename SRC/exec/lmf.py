#!/usr/bin/env python3
def unload_library(library):
    import ctypes
    try:
        dlclose = ctypes.cdll.LoadLibrary("libdl.so.1").dlclose
        dlclose.restype = ctypes.c_int
        dlclose.argtypes = [ctypes.c_void_p]
        dlclose(library._handle)
    except :
        pass
    try:
        dlclose = ctypes.cdll.LoadLibrary("libdl.so.2").dlclose
        dlclose.restype = ctypes.c_int
        dlclose.argtypes = [ctypes.c_void_p]
        dlclose(library._handle)
    except Exception as e:
        print(f"Error during dlclose: {e}")
    del library

from setcomm import callF,setcommF,getlibF
from mpi4py import MPI
import sys,os,glob
import functools
print = functools.partial(print, flush=True)
arglist=' '.join(sys.argv[1:])
scriptpath = os.path.dirname(os.path.realpath(__file__))+'/'

# MPIworld
commw = MPI.COMM_WORLD
sizew = commw.Get_size()
rankw = commw.Get_rank()
master_mpi= rankw==0

# lmfa -------------------
group=[0] #only rank=0 is used
stdout='llmfa'
if(master_mpi): print('=== Run lmfa by ranks=',group, ' See stdout=',stdout)
flib = getlibF(scriptpath+'/libecaljF.so',prt=master_mpi) #load dynamic library
comm = setcommF(grp=group)
print('rankw=',rankw,' comm=',comm,' group=',group)
if(rankw in group): 
    callF(flib.setcmdpathc,[scriptpath,master_mpi])  # Set path for ctrl2ctrlp.py at m_setcmdpath
    callF(flib.m_setargsc, [arglist,master_mpi])     # Set args at m_args
    callF(flib.sopen,[stdout]) #standard output
    callF(flib.lmfa,[comm])
    callF(flib.sclose)
if(master_mpi): print('=== end of lmfa ===')
unload_library(flib)

# lmf ---------------------
if(master_mpi): print('=== Run lmf by ranks=',group, ' See stdout=',stdout)
print('rankw=',rankw,' comm=',comm,' group=',group)
group=[i for i in range(sizew)] #used ranks
stdout='llmf'
flib = getlibF(scriptpath+'/libecaljF.so',prt=master_mpi) #load dynamic library
comm = setcommF(grp=group) #communicator for group
if(rankw in group): 
    callF(flib.setcmdpathc,[scriptpath,master_mpi])  # Set path for ctrl2ctrlp.py at m_setcmdpath
    callF(flib.m_setargsc, [arglist,   master_mpi])  # Set args at m_args
    callF(flib.sopen,[stdout]) #standard output
    callF(flib.lmf,  [comm])   #main part
    callF(flib.sclose) 
unload_library(flib)
if(master_mpi): print('=== end of lmf ===')
