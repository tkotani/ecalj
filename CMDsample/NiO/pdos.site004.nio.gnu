#!/usr/bin/gnuplot -persist
set title "PDOS s,p,d,f,g division in MuffinTin"
set xlabel " eV (zero at the Fermi energy or top of valence)"
set ylabel " number of states/(cell Ry)"
set yrange [-50:50]
plot "dos.isp1.site004.nio" u ($1*13.605):($2) w l, '' u ($1*13.605):($3+$4+$5) w l, "" u ($1*13.605):($6+$7+$8+$9+$10) w l, "" u ($1*13.605):($11+$12+$13+$14+$15+$16+$17) w l, "" u ($1*13.605):($18+$19+$20+$21+$22+$23+$24+$25+$26) w l, "dos.isp2.site004.nio" u ($1*13.605):(-$2) w l, '' u ($1*13.605):(-$3-$4-$5) w l, "" u ($1*13.605):(-$6-$7-$8-$9-$10) w l, "" u ($1*13.605):(-$11-$12-$13-$14-$15-$16-$17) w l, "" u ($1*13.605):(-$18-$19-$20-$21-$22-$23-$24-$25-$26) w l
