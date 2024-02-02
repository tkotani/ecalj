#!/usr/bin/env python3 
#>mpirun -np 4 pysample
from setcomm import callF,setcomm,mklpath
from mpi4py import MPI
import sys,os

# load ecalj and MKL libraries. Return ecalj.foobar and communicator
scriptpath = os.path.dirname(os.path.realpath(__file__))+'/'
ecaljF,rank,size,comm= setcomm(scriptpath+'/libecaljF.so', mklpath) #We use libecaljF.so in the same directory as this script.
#print('scriptpath=',scriptpath)

callF(ecaljF.setcmdpathc,[scriptpath,len(scriptpath)]) # Set path for ctrl2ctrlp.py in m_setcmdpath
arglist=' '.join(sys.argv[1:])
callF(ecaljF.m_setargsc,[arglist,len(arglist)])        # Set args in m_args
print('rank=',rank,'size=',size)
callF(ecaljF.lmfa) # Initialize MPI
print("OK! end of lmfaxxxx",rank)
#prgnam='LMFA'
#callF(ecaljF.set_prgnamc,[prgnam]) # Initialize MPI