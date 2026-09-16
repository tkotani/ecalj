#!/bin/sh
#$ -N SmP23READPSKIPF
#$ -pe smp 32
#$ -q all.q
#$ -cwd
#$ -V
#$ -S /bin/bash


# # # # # # # PARAMETERS # # # # # # # # 
export OMP_NUM_THREADS=1
MATERIAL=smp
# # # LMF # # #
NSLOTS_LMF=20
# # # QSGW # # #
ITERATION=3
NSLOTS_QSGW=32
# # # Band & DOS # # #
NSLOTS_BAND=32
NSLOTS_DOS=32
# # # ML Wannier Function fitting # # #
mpisize=3
# # # magnon calculation # # #
mpi_magnon=10
# # # PATH # # #
EXEPATH=/home/takao/bin

# # # # # # # # # # # # CALCULATION # # # # # # # # # # # # #
# # # LMF start # # #
# $EXEPATH/lmfa $MATERIAL >& llmfa
# mpirun -np $NSLOTS_LMF $EXEPATH/lmf-MPIK $MATERIAL >& llmf

# # # Band & Density of States calculation of LMF # # #
# $EXEPATH/job_band $MATERIAL -np $NSLOTS_BAND >& job_band.out
# $EXEPATH/job_tdos $MATERIAL -np $NSLOTS_DOS >& job_tdos.out
# $EXEPATH/job_pdos $MATERIAL -np $NSLOTS_DOS >& job_pdos.out

# # # QSGW calculation start # # #
$EXEPATH/gwsc $ITERATION -np $NSLOTS_QSGW $MATERIAL > gwsc.out

# # # Band & Density of States calculation of QSGW # # #
$EXEPATH/job_band $MATERIAL -np $NSLOTS_BAND >& job_band_gw.out
$EXEPATH/job_tdos $MATERIAL -np $NSLOTS_DOS >& job_tdos.out
$EXEPATH/job_pdos $MATERIAL -np $NSLOTS_DOS >& job_pdos.out

# # # Wannier fitting setup # # #
# $EXEPATH/gwsc 0 -np $NSLOTS $MATERIAL  >& out
# $EXEPATH/job_pdos $MATERIAL -np $NSLOTS &> job_pdos.log
# $EXEPATH/job_band $MATERIAL -np $NSLOTS >& job_band.log

# # # Maximum Localized Wannier Function fitting # # #
# cp GWinput_for_MLWF GWinput
# ~/ecalj2/bin/genMLWF_vw_hmaxloc2 $MATERIAL -np $mpisize >& genmlwf_vw.log 
# gnuplot fbplot.glt
# ~/ecalj_dipole/bin/genMLWFdipole $MATERIAL -np $mpisize >& genmlwf_dipole.log

# # # magnon calculation start # # # 
# cp GWinput_for_magnon GWinput
# ~/ecaljmagnon/fpgw/exec/epsPP_magnon_chipm_mpi -np $mpi_magnon $MATERIAL >& epsPP.out  ### magnon
# ~/ecaljmagnon/fpgw/exec/epsPP_magnon_chipm_mpi -np $mpi_magnon $MATERIAL >& epsPP.out  ### shift
# gnuplot mag3d.glt
# gnuplot wanplot.glt

# # # # # # # # # # # calculation is over # # # # # # # # # # # #
