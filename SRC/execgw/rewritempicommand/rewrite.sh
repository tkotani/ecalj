#!/bin/bash 

olddir="../"
newdir="../"

i=1
IN[i]=${olddir}gw_lmfh
IN2[i]=gw_lmfh.ref
OUT[i]=${newdir}gw_lmfh

i=$i+1
IN[i]=${oldir}eps_lmfh
IN2[i]=eps_lmfh.ref
OUT[i]=${newdir}eps_lmfh

i=$i+1
IN[i]=${olddir}gwsc
IN2[i]=gwsc.ref
OUT[i]=${newdir}gwsc

i=$i+1
IN[i]=${olddir}epsPP_lmfh
IN2[i]=epsPP_lmfh.ref
OUT[i]=${newdir}epsPP_lmfh

i=$i+1
IN[i]=${olddir}epsPP_lmfh_chipm
IN2[i]=epsPP_lmfh_chipm.ref
OUT[i]=${newdir}epsPP_lmfh_chipm

n=$i

echo process $n items

for (( i = 1 ; i <= $n ; ++i ))
do
echo ${IN[$i]} "->" ${IN2[$i]} "->" ${OUT[$i]} 

#gawk -f makeref.awk ${IN[$i]} > ${IN2[$i]}

gawk -f ref2cmd.awk ${IN2[$i]} > ${OUT[$i]}
chmod +x ${OUT[$i]}

echo ${OUT[$i]} is made.

done
