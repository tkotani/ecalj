#!/bin/sh
#$ -N La2CuO4swjAGAIN
#$ -pe smp 32
#$ -cwd
#$ -M takaokotani@gmail.com
#$ -m be
#$ -V
#$ -S /usr/bin/bash

export OMP_NUM_THREADS=1
EXEPATH=/home/takao/bin/
$EXEPATH/gwsc -np 32 10 la2cuo4swj -vssig=1.0 >& out2


