#!/usr/bin/env python3 
#>mpirun -np 4 pysample
from setcomm import callF,ecaljF,rank,size,comm
from mpi4py import MPI
import sys,os
arglist=' '.join(sys.argv[1:])
scriptpath = os.path.dirname(os.path.realpath(__file__))+'/'
callF(ecaljF.setcmdpathc,[scriptpath])  # Set path for ctrl2ctrlp.py in m_setcmdpath
callF(ecaljF.m_setargsc,[arglist])      # Set args in m_args
callF(ecaljF.lmfa) 
