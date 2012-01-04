dir=tmp
for file in fp/bndfp.F fp/chimedit.F fp/ioden.F fp/lmaux.F fp/lmfp.F \
fp/rsibl_ev.F fp/supot.F subs/lattic.F subs/mksym.F \
subs/rdctrl2.F subs/rdsigm2.F subs/suham.F
do echo $file; gawk -f $dir/sub_rv_p_opos.awk $file >x;mv x $file; done

