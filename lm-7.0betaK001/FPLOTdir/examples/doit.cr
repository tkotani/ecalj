#! /bin/sh

fplot -frme 0,'sqrt(0.5)',0,1 -tmx '1;0' -tmy '1;0' -noxn -noyn \
      -x 0,1 -y 0,1 -con 0.045,0.055,0.065,0.075 -nc=101 chgd.cr \
      -font t18 \
      -font h14 \
      -lblu 0.17,0.556 cc '45' \
      -lblu 0.28,0.355 cc '55' \
      -lblm 10.8,54.4 rc 'charge density in bcc Chromium' \
      -lblm 10.8,24.4 rc 'contours: 45,55,65,75 (10^{-3} a.u.)' 

