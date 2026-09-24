from setcomm import callF,setcommF,getlibF
import sys,os,glob
import functools
print = functools.partial(print, flush=True)

#import numpy as np
arglist=' '.join(sys.argv[1:])
scriptpath = os.path.dirname(os.path.realpath(__file__))+'/'

# lmf
#group=[i for i in range(sizew)] #used ranks
stdout='llmf'
flib = getlibF(scriptpath+'/libecaljF.so',True)#,prt=master_mpi) #load dynamic library
#comm = setcommF(grp=group) #communicator for group
#callF(flib.setcmdpathc,[scriptpath,master_mpi])  # Set path for ctrl2ctrlp.py at m_setcmdpath
#callF(flib.m_setargsc, [arglist,   master_mpi])  # Set args at m_args
#callF(flib.sopen,[stdout]) #standard output
callF(flib.lmf) #,  [comm])   #main part
#callF(flib.sclose) 
#flib.dlclose(flib._handle) #close library
#if(master_mpi): print('=== end of lmf ===')
