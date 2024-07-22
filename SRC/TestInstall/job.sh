#!/bin/sh
#PBS -q i1accs
##PBS -q i2cpu
#PBS -l select=1:ncpus=16:mpiprocs=4:ompthreads=1

# ulimit -s unlimited
# module purge
# module use --append /home/k0413/k041300/local/nvidia/hpc_sdk/modulefiles
# module load nvhpc-nompi/24.1 openmpi_nvhpc intel

BIN=~/ecaljdeveloper/bin
export PATH=$BIN:$PATH

# ./testecalj.py -np 8 gas_epsPP_lmfh
./testecalj.py -np 4 gwall
