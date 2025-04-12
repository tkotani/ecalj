#!/bin/sh
#PBS -q i1accs
#PBS -l select=1:ncpus=64:mpiprocs=64:ompthreads=1
#PBS -N inas2gasb2
#
export NV_ACC_TIME=1
id=inas2gasb2
## QSGW cycle: number of cycles set to 1
gwsc -np 64 -np2 4 --gpu 1 $id > lgwsc

exit

### total dos calculation with spin-orbit coupling
#mpirun lmf $id  > llmf
#job_tdos $id -np 128 -vnspin=2 -vso=1 -vnk1=12 -vnk2=12 -vnk3=4 > ljob_tdos
#
###  partial dos calculation with spin-orbit coupling
#mpirun lmf $id  > llmf
#job_pdos $id -np 128 -vnspin=2 -vso=1 -vnk1=12 -vnk2=12 -vnk3=4 > ljob_tdos
#
### band structure calculation with spin-orbit coupling
#mpirun lmf $id  > llmf
#job_band $id -np 128  -vnspin=2 -vso=1 >ljob_band
