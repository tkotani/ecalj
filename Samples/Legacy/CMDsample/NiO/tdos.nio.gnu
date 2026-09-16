#!/usr/bin/gnuplot -persist
set title "Total DOS"
set xlabel " eV (zero at the Fermi energy or top of valence)"
set ylabel " number of states/(cell Ry)"
set yrange [-50:50]
plot 'dos.tot.nio' u ($1*13.605):2 w l, '' u ($1*13.605):(-$3) w l
