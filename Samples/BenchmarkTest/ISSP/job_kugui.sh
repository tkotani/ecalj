#!/bin/sh
#PBS -q i1accs
#PBS -l select=1:ncpus=64:mpiprocs=64:ompthreads=1
#PBS -N inas2gasb2
#
# ISSP system C (kugui), GPU node: one QSGW iteration of Samples/BenchmarkTest/inas2gasb2 (or inas4gasb4: set id).
# Copy the sample directory to your work directory, put this script into it, and qsub job_kugui.sh there.
# (2026-09-30: brought to the present input ctrlg.<id>.toml and options; not run on kugui since.)
export NV_ACC_TIME=1
id=inas2gasb2
gwsc 1 -np 64 -np2 4 --gpu $id > lgwsc
# compare QPU.1run with the one of the sample directory
exit

### total DOS, partial DOS and bands with spin-orbit coupling, after the QSGW cycle
#job_tdos $id -np 64 --ctrlg:ham.nspin=2 --ctrlg:ham.so=1 --ctrlg:bz.nkabc=[12,12,4] > ljob_tdos
#job_pdos $id -np 64 --ctrlg:ham.nspin=2 --ctrlg:ham.so=1 --ctrlg:bz.nkabc=[12,12,4] > ljob_pdos
#job_band $id -np 64 --ctrlg:ham.nspin=2 --ctrlg:ham.so=1 > ljob_band
