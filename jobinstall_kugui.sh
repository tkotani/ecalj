#!/bin/sh
#PBS -q i2cpu
#PBS -l select=1:ncpus=16:mpiprocs=8:ompthreads=1
#PBS -N install 

ulimit -s unlimited
module purge
module load nvhpc-nompi/24.7 openmpi_nvhpc compiler-rt tbb mkl

export PATH=$HOME/bin:$PATH

rm -rf SRC/TestInstall/bin
FC=nvfortran ./InstallAll --gpu
