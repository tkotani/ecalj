#!/usr/bin/gnuplot
set terminal postscript enhanced color eps
set output 'Sisigma.eps'
set grid
set title "Sigma (near Ef)"
set ylab "Energy(eV)"
plot 'SEComg.UP' u ($10):($11) w l,'' u ($10):($12) w l
