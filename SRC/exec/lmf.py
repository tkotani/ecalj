from setcomm import callF,setcommF,getlibF,unload_library
from mpi4py import MPI
import sys,os,glob
import functools

#------------------------------------------------
print = functools.partial(print, flush=True)
arglist=' '.join(sys.argv[1:])
scriptpath = os.path.dirname(os.path.realpath(__file__))+'/'
# MPIworld
commw = MPI.COMM_WORLD
sizew = commw.Get_size()
rankw = commw.Get_rank()
master_mpi= rankw==0
if(master_mpi): print()

# lmfa -------------------
group=[0] #only rank=0 is used
stdout='llmfa'
if(master_mpi): print('=== Run lmfa by ranks=',group, ' See stdout=',stdout)

comm = setcommF(grp=group)
if(rankw in group): 
    if(master_mpi): print(' rankw=',rankw,' group=',group)
    flib = getlibF(scriptpath+'/libecaljF.so',prt=master_mpi) #load dynamic library
    callF(flib. setcmdpathc,[scriptpath,master_mpi])  # Set path for ctrl2ctrlp.py at m_setcmdpath
    callF(flib. m_setargsc, [arglist,master_mpi])     # Set args at m_args
    callF(flib. m_ext_init, [])
    callF(flib. sopen,[stdout]) #standard output
    callF(flib. convertctrl2ctrlpbypython,[])
    callF(flib. lmfa,[comm])
    callF(flib. sclose)
    if(master_mpi): print('=== end of lmfa ===')
    unload_library(flib)
commw.Barrier()

# lmf ---------------------
group=[i for i in range(sizew)] #used ranks
stdout='llmf'
comm = setcommF(grp=group) #communicator for group
if(rankw in group): 
    if(master_mpi): print('\n=== start lmf by ranks=',group, ' See stdout=',stdout)
    flib = getlibF(scriptpath+'/libecaljF.so',prt=master_mpi) #load dynamic library
    callF(flib. setcmdpathc,[scriptpath,master_mpi])  # Set path for ctrl2ctrlp.py at m_setcmdpath
    callF(flib. m_setargsc, [arglist,   master_mpi])  # Set args at m_args
    callF(flib. m_ext_init, [])
    callF(flib. sopen, [stdout]) #standard output
    callF(flib. convertctrl2ctrlpbypython,[])
    callF(flib. lmf,  [comm])    #main part
    callF(flib. sclose) 
    unload_library(flib)
commw.Barrier()
if(master_mpi): print('=== end of lmf ===')

if(master_mpi): print()
