 fplot: ------ starting new frame ------
 pltsts:  set  logx=F  logy=F
 pltstp:  frame corners (GU): lb = (0 0)  rt = (.707107 1)
 fplot file "examples/chgd.cr"  pass 1
 mapdat: 101 points
 fplot:ic  contour     line type
        1  .045        1  3  2  .5
        2  .055        1  3  2  .5
        3  .065        1  3  2  .5
        4  .075        1  3  2  .5
 fplot file "examples/chgd.cr"  pass 2
 mapdat: 101 points
 FRME : axes at x=0 1  y=0 1  bold=3
       xt1=0 tsx=1 mtx=2  yt1=0 tsy=1 mty=2
 pstr: cc: blk=T  "45" at .17,.556(u) 191.6,364.2(m)
 pstr: cc: blk=T  "55" at .28,.355(u) 222.4,284.6(m)
 pstr: rc: blk=T  "charge density in bcc Chromium" at -.4757,-.2263(u) 10.8,54.4
 pstr: rc: blk=T  "contours: 45,55,65,75 (10^{-3} a.u.)" at -.4757,-.302(u) 10.8
