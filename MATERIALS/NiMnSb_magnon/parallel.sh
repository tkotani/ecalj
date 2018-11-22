#!/bin/sh
#$ -N nms_Lm
#$ -pe smp 10
#$ -cwd
#$ -V
#$ -S /bin/bash

export OMP_NUM_THREADS=1
MATERIAL=nimnsb
ITERATION=10
EXEPATH=/home/okumura/bin/

# ### 1. band calculation and create MLWFs
# # job_band
cp GWinput_for_MLWF GWinput
$EXEPATH/job_band $MATERIAL -np $NSLOTS &> job_band.log
genMLWF_vw $MATERIAL -np $NSLOTS &> genmlwf_vw.log

# ### 2. magnon calculation 
# cp GWinput_for_magnon GWinput
# epsPP_magnon_chipm_mpi -np $NSLOTS $MATERIAL
