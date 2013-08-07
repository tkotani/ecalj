#!/bin/bash -x 
IN=../gw_lmfh
IN2=gw_lmfh.ref
OUT=../gw_lmfh.new

IN=../eps_lmfh
IN2=eps_lmfh.ref
OUT=../eps_lmfh.new

IN=../gwsc
IN2=gwsc.ref
OUT=../gwsc.new

#IN=../epsPP_lmfh
#IN2=epsPP_lmfh.ref
#OUT=../epsPP_lmfh.new

#IN=../epsPP_lmfh_chipm
#IN2=epsPP_lmfh_chipm.ref
#OUT=../epsPP_lmfh_chipm.new

gawk -f makeref.awk $IN >$IN2
gawk -f ref2cmd.awk $IN2 >$OUT

chmod +x $OUT
echo $OUT is made.
