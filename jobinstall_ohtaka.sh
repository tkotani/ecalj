#!/bin/sh
#SBATCH -p i8cpu
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -c 1
#SBATCH --job-name=install
#SBATCH --ntasks-per-node=8

cp SRC/ISSPohtaka/run_arg.py SRC/ISSPohtaka/job_tdos SRC/exec/
cp SRC/ISSPohtaka/creplot.py SRC/ISSPohtaka/jobsubmit.py ecalj_auto/auto/
rm -rf SRC/TestInstall/bin
FC=ifort ./InstallAll
