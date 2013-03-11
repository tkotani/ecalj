#!/usr/bin/gnuplot -persist
set title "Total DOS"
set xlabel " eV (zero at the Fermi energy or top of valence)"
set ylabel " number of states/(cell Ry)"
set yrange [0:50]
plot 'dos.tot.wgan' u ($1*13.605):2 w l
