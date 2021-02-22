#!/bin/bash

BINDIR=${HOME}/bin
MATH='-mkl'

function my_error() {
    exit 1
}

###
mkdir ${BINDIR}
echo Going to install required binaries and scripts to ${BINDIR} !

### make a link to getsyml
rm -f ${BINDIR}/getsyml
ln -s ${PWD}/GetSyml/getsyml.py ${BINDIR}/getsyml

### viewvesta
pushd .
cd StructureTool/
./makelink $BINDIR
popd

### 

pushd .
cd SRC/exec
make PLATFORM=ifort LIBMATH="$MATH" FC=mpifort LK=mpifort ||my_error
make PLATFORM=ifort BINDIR=$BINDIR  FC=mpifort LK=mpifort install ||my_error
popd

echo '=== goto test ==='
cd SRC/TestInstall/
make mpi_size=4 all

