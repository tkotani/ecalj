#!/usr/bin/env python3
# sample mpirun -np 4 pysample
from setcomm import callF,setcomm
import ctypes
from mpi4py import MPI
#'''An example to run MPI-fortran codes successively.
#   To get ecaljF.so, run >mpif90 -j -shared -fPIC -o ecaljfortran.so *.f90'''

# 1. Set bind(C) for main_foobar
# 2. Remove mpi_init
# 3. arguments, m_ext_init, and cmdopt for getarg

mkl = '/usr/lib/x86_64-linux-gnu/libmkl_rt.so'
ctypes.CDLL(mkl, mode=ctypes.RTLD_GLOBAL)     #load mkl

comm = MPI.COMM_WORLD # MPI
comm = comm.py2f()
ecaljF = setcomm("gfortran/lib/libecaljF.so",comm) #load ecalj and send comm to fortran

callF( ecaljF.hello   )
callF( ecaljF.hello2, [3429] )
callF( ecaljF.hello3, [True, 34, 8.2, 9.3] )
