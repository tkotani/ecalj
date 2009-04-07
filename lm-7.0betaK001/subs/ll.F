      integer function ll(ilm)
C- Returns l, given lm index
C     implicit none
C Passed parameters
      integer ilm
C Local parameters
      integer lmaxx
      parameter (lmaxx=17)
      integer lla(lmaxx**2)
      data lla/0,3*1,5*2,7*3,9*4,11*5,13*6,15*7,17*8,19*9,
     .   21*10,23*11,25*12,27*13,29*14,31*15,33*16 /
      ll = lla(ilm)
      end

