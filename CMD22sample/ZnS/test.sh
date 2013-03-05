#!/bin/csh
#$ -N ZnS
#$ -q c6.q
#$ -pe mpi 2

cd $HOME/cmdsampl1/ZnS
gwsc 5 -np 2 zns >& log_temp
