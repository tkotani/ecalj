#!/bin/bash +x
# we need GWinput,ctrl,syml.fe, fbplot.glt, mag3d.glt, wanplot.glt
MATERIAL=ni
NSLOTS=8
lmfa fe >& llmfa
mpirun -np $NSLOTS lmf fe >&llmf
# ### 1. band calculation and create MLWFs
job_band $MATERIAL -np $NSLOTS NoGnuplot # &> job_band.log
genMLWF_vw $MATERIAL -np $NSLOTS # &> genmlwf_vw.log

### 2. magnon calculation
### If we implement more detailed k, change n1n2n3 in "GWinput_for_magnon".
epsPP_magnon_chipm_mpi -np $NSLOTS $MATERIAL

### 3. plotting Stoner and Magnon (eps output) 
# fat band plot (Wannier) --> "fbplot.eps"
# Stoner and Magnon --> "fe9_kmat.eps" and "fe9_rmat.eps"
# 3D plot for Magnon --> "magnon3d.eps"
gnuplot fbplot.glt 
gnuplot wanplot.glt
gnuplot mag3d.glt

echo " ======================================================"
echo " Magnon calculation finished"
echo " 'wan_ChiPMr.dat' <--- R(q,omega)"
echo " 'wan_ChiPMz.dat' <--- K(q,omega)"
echo " Check '*.eps' files"
echo " Compare the results to the prepared eps file in ./eps/"
echo " ======================================================"
