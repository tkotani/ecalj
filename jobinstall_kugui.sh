#!/bin/sh
#PBS -q i1accs
#PBS -l select=1:ncpus=64:mpiprocs=64:ompthreads=1
#PBS -N installi1accs
##PBS -q i2cpu
##PBS -l select=1:ncpus=16:mpiprocs=8:ompthreads=1
##PBS -N install 

ulimit -s unlimited
module purge
module load nvhpc-nompi/24.7 openmpi_nvhpc compiler-rt tbb mkl

export PATH=$HOME/bin:$PATH

cp SRC/ISSPkugui/run_arg.py SRC/ISSPkugui/job_tdos SRC/exec/
rm -rf SRC/TestInstall/bin
FC=nvfortran ./InstallAll --gpu
