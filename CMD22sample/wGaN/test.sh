#!/bin/csh
#$ -N wGaN
#$ -q c6.q
#$ -pe mpi 4

cd $HOME/cmdsampl1/wGaN
gwsc 5 -np 4 wgan >& log_temp
