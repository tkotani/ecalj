# This version allow only integer:: real(8):: logical:: can be passed to fortran. No return.
# This should be easily replaced by fortran code.
#mpif90 -shared -fPIC -o ecaljF.so *.f90
#mpiexec -n 4 python3 ./hello.py | sort
from ctypes import (CDLL, POINTER, c_int32, c_double, c_bool, c_char, ARRAY, byref, create_string_buffer)
#import numpy as np
from mpi4py import MPI
import ctypes,sys,os

def setcommF(grp):
    ''' Pass communicator to module m_comm.f90 
        Setcomm requires a shared library libfoobar.so, generaged by
        >mpif90 -shared -fPIC -o ecaljF.so m_comm.f90 fmath.f90'''
    commw = MPI.COMM_WORLD
    group = commw.Get_group()
    new_group = group.Incl(grp) #if(rankw==0) else group.Incl([1,2,3])
    comm = commw.Create(new_group)
    if( not (commw.Get_rank() in grp) ):
        return None,None
    commF = comm.py2f()
    return commF

def getlibF(fortranso,prt):
    ''' Get fortran library flib callable from python. mkl is taken from mklloc.txt '''
    scriptpath = os.path.dirname(os.path.realpath(__file__))+'/'
    mkl=[]
    with open(scriptpath+'./mklloc.txt') as f:
        for line in f:
            if 'mkl' in line:
                mmm=line.split('=>')[1].strip().split(' ')[0].strip()             #print(f'mmm=',mmm)
                mkl.append(mmm)
    for mkll in mkl:
        CDLL(mkll, mode=ctypes.RTLD_GLOBAL)     #flib = np.ctypeslib.load_library(fortranso,".")
        if(prt): print(f' Load mkl=',mkll)
    flib = CDLL(fortranso, mode=ctypes.RTLD_GLOBAL)
    if(prt): print(f' Load fortran lib=',fortranso)
    return flib

class callF:
    import sys
    def __init__(self,foobar,arguments=[]):
        '''Equivalent to 'call foobar(a,b,c,...)' in fortran, where we supply arguments=[a,b,c,...]. 
        a,b,c,... are integer,logical, or real(8) in this version of callF.'''
        MAXSTRLEN=512
        c_char_array = ARRAY(c_char,MAXSTRLEN)
        argtypess=[]
        data=[]
        for ii in arguments:
            if(type(ii)==type(1)):  
                argtypess.append( POINTER(c_int32) ) # data type
                data.append(c_int32(ii))             # data 
            elif(type(ii)==type(1.0)):  
                argtypess.append( POINTER(c_double) )
                data.append(c_double(ii))
            elif(type(ii)==type(True)):
                argtypess.append( POINTER(c_bool) )
                data.append(c_bool(ii))
            elif(type(ii)==type('a')): #Not working well, because bind(C) in fortran allows only char(1)::aaa(:)
                argtypess.append( POINTER(c_char_array) )
                data.append(byref(create_string_buffer(ii.encode(),MAXSTRLEN)))
            #print(f' outputargs=',argtypess)
            #print(f' data=',data)
            foobar.argtypes= argtypess
        if(len(argtypess)==0):
            foobar()
        else:
            foobar(*data)
        return
