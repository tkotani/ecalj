#!/bin/bash

BINDIR=${HOME}/bin

### This is an example for ubuntu
# For each machine, we have to prepare
#  lm7k/MAKEINC/Make.inc.gfortran
#  lm7k/MAKEINC/Make.inc.gfortran_mpik
#  fpgw/exec/Make.inc.gfortran_mpik
#  fpgw/Wannier/

mkdir ${BINDIR}
echo Going to install required binaries and scripts to ${BINDIR} !

### Use lm7K/Makefile
cd lm7K/
make PLATFORM=gfortran 
make PLATFORM=gfortran      BINDIR=$BINDIR install
make PLATFORM=gfortran_mpik 
make PLATFORM=gfortran_mpik BINDIR=$BINDIR install

### Use fpgw/exec/makefile
cd ../fpgw/exec/
make PLATFORM=gfortran
make PLATFORM=gfortran BINDIR=$BINDIR  install
make PLATFORM=gfortran BINDIR=$BINDIR install2

### Use fpgw/Wannier/Makefile
cd ../Wannier/
make PLATFORM=gfortran
make BINDIR=$BINDIR PLATFORM=gfortran install

cd ../../TestInstall/
make mpi_size=4 all
