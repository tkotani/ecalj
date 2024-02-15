#!/usr/bin/env python3 
#>mpirun -np 4 pysample
from setcomm import callF,setcomm,ecaljF,rank,size,comm
from mpi4py import MPI
import sys,os
scriptpath = os.path.dirname(os.path.realpath(__file__))+'/'
arglist=' '.join(sys.argv[1:])
callF(ecaljF.setcmdpathc,[scriptpath]) # Set path for ctrl2ctrlp.py in m_setcmdpath
callF(ecaljF.m_setargsc,[arglist])        # Set args in m_args
callF(ecaljF.lmf) 
	