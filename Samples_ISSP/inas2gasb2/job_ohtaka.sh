#!/bin/sh
#SBATCH -p i8cpu
#SBATCH -N 8
#SBATCH -n 512
#SBATCH -c 2
#SBATCH --job-name=inas2gasb2
#SBATCH --ntasks-per-node=64

id=inas2gasb2
## 
## QSGW cycle: number of cycles set to 1
gwsc -np 512 1 $id > lgwsc

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
