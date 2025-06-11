# >mpirun -np 4 python lmf.py si
from setcomm import callF,setcommF,getlibF,unload_library,readflib
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
    fl = getlibF(scriptpath+'/libecaljF.so',prt=master_mpi) #load dynamic library
    callF(fl. setcmdpathc, scriptpath,master_mpi)  #Set path for ctrl2ctrlp.py at m_setcmdpath
    callF(fl. m_setargsc,  arglist,master_mpi)     #Set args at m_args
    callF(fl. m_ext_init )
    callF(fl. sopen,  stdout) #standard output
    callF(fl. convertctrl2ctrlpbypython)
    callF(fl. lmfa, comm)
    callF(fl. sclose )
    if(master_mpi): print('=== end of lmfa ===')
    unload_library(fl)
commw.Barrier()

# lmf ---------------------
group=[i for i in range(sizew)] #used ranks
stdout='llmf'
comm = setcommF(grp=group) #communicator for group
if(rankw in group): 
    if(master_mpi): print('\n=== start lmf by ranks=',group, ' See stdout=',stdout)
    fl = getlibF(scriptpath+'/libecaljF.so',prt=master_mpi) #load dynamic library
    callF(fl. setcmdpathc, scriptpath,master_mpi)  # Set path for ctrl2ctrlp.py at m_setcmdpath
    callF(fl. m_setargsc,  arglist,   master_mpi)  # Set args at m_args
    callF(fl. m_ext_init )
    callF(fl. sopen,  stdout) #standard output
    callF(fl. convertctrl2ctrlpbypython)
    callF(fl. lmf,  comm)    #main part
    callF(fl. sclose) 
    unload_library(fl)
commw.Barrier()
if(master_mpi): print('=== end of lmf ===')

if(master_mpi): print()
