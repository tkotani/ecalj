import sys,os,time,pathlib
import numpy as np
from mpi4py import MPI

start = time.perf_counter()
usage = """ USAGE: mpirun -np 4 python hx0ahc.py -4. 4. 101 """
epath=os.path.dirname(os.path.abspath(__file__))
options=""
args = sys.argv
if (len(args) < 4):
    print(usage)
    sys.exit()
elif (len(args) > 4):
    options = args[4:]
if (not len(options)==1):
    op0 = ""
    for op in options:
        op0 += op+" "
    options = op0
efs = float(args[1])
eff = float(args[2])
nd = int(args[3])
efshift = np.linspace(efs,eff,nd)

comm = MPI.COMM_WORLD 
size = comm.Get_size()
rank = comm.Get_rank()

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
    print(f"{epath}/hahc --job=202 --ahc --interbandonly -EfermiShifteV={efshift[rank*njob+i]} {options} > lahc.{rank}")
    os.system("{exe}/hahc --job=202 --ahc --interbandonly -EfermiShifteV={ef} {op} > lahc.{rank}"
              .format(exe=epath,ef=efshift[rank*njob+i],op=options,rank=rank))
if rank in residue:
    print('rank=',rank,size*njob+rank,efshift[size*njob+rank],flush=True)
    os.system("{exe}/hahc --job=202 --ahc --interbandonly -EfermiShifteV={ef} {op} > lahc.{rank}"
              .format(exe=epath,ef=efshift[size*njob+rank],op=options,rank=rank))

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
