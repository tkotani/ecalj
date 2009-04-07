 fplot: ------ starting new frame ------
 pltsts:  set  logx=F  logy=F
 pltstp:  frame corners (GU): lb = (0 0)  rt = (1 .5)
 fplot file "examples/qdos.al"  pass 1
 fplot file "-tp"  pass 1
 fplot file "-tp"  pass 1
 fplot: set line type 1 bold 3 col= .2 .4 .6
 fplot file "examples/qdos.al"  pass 2
 mapdat: 1001 points
 FRME : axes at x=-.9 .1  y=0 6.5  bold=3
       xt1=-.9 tsx=.1 mtx=1  yt1=0 tsy=1 mty=1
 key:  x,y=-.85,5.5(UU)  lt 1  sym 0  len=.125  legend=Al
 plcrv: new curve, 2 pts lt 1 bold 3 col= .2 .4 .6
 plcrv: new curve, 1001 pts lt 1 bold 3 col= .2 .4 .6
 fplot: set line type 3 bold 5 col= .6 .4 .2
 expand:  958 points generated
 fplot file "-tp"  pass 2
 mapdat: 958 points
 key:  x,y=-.85,5.067(UU)  lt 3  sym 0  len=.125  legend=Free-electrons
 plcrv: new curve, 2 pts lt 3 bold 5 col= .6 .4 .2
 plcrv: new curve, 958 pts lt 3 bold 5 col= .6 .4 .2
 fplot: set line type 2 bold 2 pat ( 2 .3 .5 .3)
 expand:  4 points generated
 fplot file "-tp"  pass 2
 mapdat: 2 points
 plcrv: new curve, 2 pts lt 2 bold 2 len=(14.4 2.16 3.6 2.16)
 pstr: rc: blk=T  "density of states in aluminium" at -1.236,-4.169(u) 10.8,17.0
