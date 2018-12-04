set term postscript enhanced color eps
set output "fe9_rmat.eps"

set multiplot layout 2,2
set key right bottom
set grid
set ytics
# set xrange[0:0.6]
set xlabel "{/Symbol w} [meV]"
set ylabel "Im[R] [1/Htr]"

ratio=7.5
set xrange[-0.02:200]
set yrange[-300*ratio:0] 
plot "wan_ChiPMr.mat" every :::0::0 u (13600*$5):($7) lw 4 lt  2  ti " 0/12" w l,\
     "wan_ChiPMr.mat" every :::1::1 u (13600*$5):($7) lw 4 lt  1  ti " 1/12" w l,\
     "wan_ChiPMr.mat" every :::2::2 u (13600*$5):($7) lw 4 lt  3  ti " 2/12" w l
    		      	    	      		      	 
set xrange[-0.01:700]
set yrange[-17.5*ratio:0]
plot "wan_ChiPMr.mat" every :::3::3 u (13600*$5):($7) lw 4 lt  2  ti " 3/12" w l,\
     "wan_ChiPMr.mat" every :::4::4 u (13600*$5):($7) lw 4 lt  1  ti " 4/12" w l,\
     "wan_ChiPMr.mat" every :::5::5 u (13600*$5):($7) lw 4 lt  3  ti " 5/12" w l

set xrange[-0.0001:700]
set yrange[-7*ratio:0]
plot "wan_ChiPMr.mat" every :::6::6 u (13600*$5):($7) lw 4 lt  2  ti " 6/12" w l,\
     "wan_ChiPMr.mat" every :::7::7 u (13600*$5):($7) lw 4 lt  1  ti " 7/12" w l,\
     "wan_ChiPMr.mat" every :::8::8 u (13600*$5):($7) lw 4 lt  3  ti " 8/12" w l

set xrange[-0.0001:700]
set yrange[-7*ratio:0]
plot "wan_ChiPMr.mat" every ::: 9:: 9 u (13600*$5):($7) lw 4 lt 2  ti " 9/12" w l,\
     "wan_ChiPMr.mat" every :::10::10 u (13600*$5):($7) lw 4 lt 1  ti "10/12" w l,\
     "wan_ChiPMr.mat" every :::11::11 u (13600*$5):($7) lw 4 lt 3  ti "11/12" w l
