#!/bin/bash

BINDIR=${HOME}/bin

### This is an example for hakozaki (ifort14) ###
# For each machine, we have to prepare
#  lm7k/MAKEINC/Make.inc.ifort
#  lm7k/MAKEINC/Make.inc.ifort_mpik
#  fpgw/exec/Make.inc.ifort_mpik
#  fpgw/Wannier/

mkdir ${BINDIR}
echo Going to install required binaries and scripts to ${BINDIR} !

### Use lm7K/Makefile
cd lm7K/
make PLATFORM=ifort 
make PLATFORM=ifort      BINDIR=$BINDIR install
make PLATFORM=ifort_mpik 
make PLATFORM=ifort_mpik BINDIR=$BINDIR install

### Use fpgw/exec/makefile
cd ../fpgw/exec/
make PLATFORM=ifort
make PLATFORM=ifort BINDIR=$BINDIR  install
make PLATFORM=ifort BINDIR=$BINDIR install2

### Use fpgw/Wannier/Makefile
cd ../Wannier/
make PLATFORM=ifort
make BINDIR=$BINDIR PLATFORM=ifort install

cd ../../TestInstall/
make mpi_size=4 all
