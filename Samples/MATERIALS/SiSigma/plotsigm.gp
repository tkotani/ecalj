#!/usr/bin/gnuplot
set terminal postscript enhanced color eps
set output 'Sisigma.eps'
set grid
set title "Sigma (near Ef)"
set ylab "Energy(eV)"
plot 'SEComg.UP' u ($9):($10) w l,'' u ($9):($11) w l
