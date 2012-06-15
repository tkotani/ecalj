#!/bin/bash
#@$-lP 1
#@$-lp 1
#@$-lm 3.5gb
#@$-eo
#@$-oi

#cd ${QSUB_WORKDIR}


#LD_LIBRARY_PATH=/home/usr5/f70205a/local/lib:/home/usr5/f70205a/local/lib64:$LD_LIBRARY_PATH
#LD_RUN_PATH=/home/usr5/f70205a/local/lib:/home/usr5/f70205a/local/lib64:$LD_RUN_PATH
#PATH=${HOME}/bin:${HOME}/local/bin:$PATH

source atomlist.bash
source homodimerdistance.bash
source extra.bash
#source extra_nrel.bash
REPLACEHERE___
date
