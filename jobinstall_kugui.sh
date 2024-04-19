#!/bin/sh
#PBS -q i2cpu
#PBS -l select=1:ncpus=128:mpiprocs=4:ompthreads=6

#module load nvhpc-nompi openmpi_nvhpc intel
module purge
module load intel intel-mpi
ulimit -s unlimited

ulimit -a
rm -rf SRC/TestInstall/bin
FC=ifort ./InstallAll
