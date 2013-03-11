#!/bin/csh
#$ -N NiO
#$ -q c6.q
#$ -pe mpi 4

cd $HOME/cmdsampl1/NiO
gwsc 5 -np 4 nio >& log_temp
