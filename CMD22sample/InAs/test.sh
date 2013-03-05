#!/bin/csh
#$ -N InAs
#$ -q c6.q
#$ -pe mpi 1

cd $HOME/cmdsampl1/InAs
gwsc 5 -np 1 inas >& log_temp
