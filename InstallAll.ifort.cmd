#!/bin/bash

BINDIR=${HOME}/bin
MATH='-mkl'

### This is an example for hakozaki (ifort14) ###
# For each machine, we have to prepare
#  lm7k/MAKEINC/Make.inc.ifort
#  lm7k/MAKEINC/Make.inc.ifort_mpik.cmd
#  fpgw/exec/Make.inc.ifort_mpik.cmd
#  fpgw/Wannier/

mkdir ${BINDIR}
echo Going to install required binaries and scripts to ${BINDIR} !

### viewvesta
cd StructureTool/
./makelink

### Use lm7K/Makefile
cd lm7K/
make PLATFORM=ifort.cmd LIBMATH="$MATH"
make PLATFORM=ifort.cmd      BINDIR=$BINDIR install
make PLATFORM=ifort_mpik.cmd LIBMATH="$MATH"
make PLATFORM=ifort_mpik.cmd BINDIR=$BINDIR install

### Use fpgw/exec/makefile
cd ../fpgw/exec/
make PLATFORM=ifort.cmd LIBMATH="$MATH"
make PLATFORM=ifort.cmd BINDIR=$BINDIR  install
make PLATFORM=ifort.cmd BINDIR=$BINDIR install2

### Use fpgw/Wannier/Makefile
cd ../Wannier/
make PLATFORM=ifort.cmd LIBMATH="$MATH"
make BINDIR=$BINDIR PLATFORM=ifort.cmd install

cd ../../TestInstall/
make mpi_size=6 all
