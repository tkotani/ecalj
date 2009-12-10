#! /bin/sh

fplot -frme 0,1,0,0.5 -frmt 2 \
      -p0 -x -0.9,0.1 -y 0,6.5 -tmx 0.1:1 -tmy 1:1 \
      -xl 'Energy (Ry)' -yl 'Density of States\(states/Ry/atom)' \
      -font t14 -k -0.85,5.5 \
      -lt 1,bold=3,col=.2,.4,.6 -l 'Al' -nc=2 ./qdos.al \
      -lt 3,bold=5,col=.6,.4,.2 -l 'Free-electrons' \
      -ord '5.6026*(x+0.8575)^0.5' -tp -0.857:0.1:0.001 \
      -lt 2,bold=2,2,.3,.5,.3 -tp 2~-0.0237,6.5,-0.0237,0 \
      -font t18 \
      -lblm 10.8,17 rc 'density of states in aluminium'


