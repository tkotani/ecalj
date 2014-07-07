set terminal postscript enhanced color eps
set output 'BandCuMLWfig1b.eps'
set grid
set ylab "Energy(Ry*13.60500)"
set yra [   -12.0:    10.00000]
set xtics ( ""         0.0000000000,\
 ""         1.0000000000,\
 ""         1.5000000000,\
 ""         2.2071067812,\
 ""         3.0731321850,\
 ""         4.1337923568)
 plot \
 "bnd1.dat" u 2:3 lt 1 pt 1 not w l,\
 "bnd2.dat" u 2:3 lt 1 pt 1 not w l,\
 "bnd3.dat" u 2:3 lt 1 pt 1 not w l,\
 "bnd4.dat" u 2:3 lt 1 pt 1 not w l,\
 "bnd5.dat" u 2:3 lt 1 pt 1 not w l, \
 'bnds.maxloc.up' using ($5):($6) w p lt 3
 # pause -1 (instead, gnuplot -p ThisScript)
