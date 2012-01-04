dir=tmp
for nn in subs/rdctrl2.F lmfav7.F lmv7.F lmv7util.F
do gawk -f $dir/sspec.awk $nn >x;mv x $nn; done
