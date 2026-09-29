import sys,os,time,pathlib,shlex
import numpy as np
# Two ways to run (2026-09-30):
#  python3 hx0ahc.py ...             one process; hahc runs through "mpirun -np 1", one point after another.
#                                    mpi4py is not loaded (it is not needed, and need not be installed).
#  mpirun -np N python3 hx0ahc.py ...  the points are divided among the N ranks (mpi4py built with the MPI of this
#                                    mpirun); each rank starts hahc itself, without mpirun.
under_launcher = any(k in os.environ for k in ('OMPI_COMM_WORLD_SIZE','PMI_SIZE','PMIX_RANK','PMI_RANK'))
MPI = None
if under_launcher:
    try:
        from mpi4py import MPI
    except ImportError:
        MPI = None

start = time.perf_counter()
usage = """ USAGE: mpirun -np 4 python hx0ahc.py -4. 4. 101 [options for hahc, e.g. --ctrlg:ham.so=1] """
epath=os.path.dirname(os.path.abspath(__file__))
args = sys.argv
if (len(args) < 4):
    print(usage)
    sys.exit()
# Options for hahc as one string, each quoted for the shell (--ctrlg:gw.n1n2n3=[4,4,4] has brackets).
# Bug fixed 2026-09-30: with exactly one option the list itself was formatted into the command, as ['--opt'].
options = " ".join(shlex.quote(op) for op in args[4:])
efs = float(args[1])
eff = float(args[2])
nd = int(args[3])
efshift = np.linspace(efs,eff,nd)

if MPI is None:
    class _Serial:        # stands for the communicator of one process
        def Get_size(self): return 1
        def Get_rank(self): return 0
        def barrier(self): pass
    comm = _Serial()
else:
    comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
# a program started without mpirun stops in MPI_Init with the OpenMPI of the NVIDIA HPC SDK: one rank through mpirun
hahc = f"{epath}/hahc" if MPI is not None else f"mpirun -np 1 {epath}/hahc"

if rank == 0:
    pahc = pathlib.Path('ahc_tet.isp11.dat')
    psum = pathlib.Path('sum.isp11.dat')
    if pahc.exists():
        os.system('rm ahc*')
    if psum.exists():
        os.system("rm sum*")

# devide jobs
njob = nd//size
residue = [i for i in range(nd-njob*size)]
comm.barrier()
for i in range(njob):
    print('rank=',rank,rank*njob+i,efshift[rank*njob+i],flush=True)
    print(f"{hahc} --job=202 --ahc --interbandonly --EfermiShifteV={efshift[rank*njob+i]} {options} > lahc.{rank}")
    os.system("{exe} --job=202 --ahc --interbandonly --EfermiShifteV={ef} {op} > lahc.{rank}"
              .format(exe=hahc,ef=efshift[rank*njob+i],op=options,rank=rank))
if rank in residue:
    print('rank=',rank,size*njob+rank,efshift[size*njob+rank],flush=True)
    os.system("{exe} --job=202 --ahc --interbandonly --EfermiShifteV={ef} {op} > lahc.{rank}"
              .format(exe=hahc,ef=efshift[size*njob+rank],op=options,rank=rank))

comm.barrier()
if rank == 0:
    pahc = pathlib.Path('ahc_tet.isp22.dat')
    files = ["ahc_tet","ahc_sp","sum"]
    for f in files:
        os.system("grep -v '#' {head}.isp11.dat | sort --sort=numeric > {head}_converted.isp11.dat".format(head=f))
        os.system("mv {head}_converted.isp11.dat {head}.isp11.dat".format(head=f))
        if pahc.exists():
            os.system("grep -v '#' {head}.isp22.dat | sort --sort=numeric > {head}_converted.isp22.dat".format(head=f))
            os.system("mv {head}_converted.isp22.dat {head}.isp22.dat".format(head=f))
    end=time.perf_counter()
    print("Computation time: ", end-start)
