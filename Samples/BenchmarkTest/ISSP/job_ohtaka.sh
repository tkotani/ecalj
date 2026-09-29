#!/bin/sh
#SBATCH -p i8cpu
#SBATCH -N 8
#SBATCH -n 512
#SBATCH -c 2
#SBATCH --job-name=inas2gasb2
#SBATCH --ntasks-per-node=64
#
# ISSP system B (ohtaka), CPU: one QSGW iteration of Samples/BenchmarkTest/inas2gasb2 (or inas4gasb4: set id).
# Copy the sample directory to your work directory, put this script into it, and sbatch job_ohtaka.sh there.
# (2026-09-30: brought to the present input ctrlg.<id>.toml and options; not run on ohtaka since.)
id=inas2gasb2
gwsc 1 -np 512 $id > lgwsc
# compare QPU.1run with the one of the sample directory
exit

### total DOS, partial DOS and bands with spin-orbit coupling, after the QSGW cycle
#job_tdos $id -np 128 --ctrlg:ham.nspin=2 --ctrlg:ham.so=1 --ctrlg:bz.nkabc=[12,12,4] > ljob_tdos
#job_pdos $id -np 128 --ctrlg:ham.nspin=2 --ctrlg:ham.so=1 --ctrlg:bz.nkabc=[12,12,4] > ljob_pdos
#job_band $id -np 128 --ctrlg:ham.nspin=2 --ctrlg:ham.so=1 > ljob_band
