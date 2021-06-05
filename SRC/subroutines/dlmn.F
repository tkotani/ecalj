      subroutine dlmn(nbas,plat,bas,iax,wk)
C- Calculate direction cosines and distances for one atom pair
C ----------------------------------------------------------------
Ci Inputs
Ci   nbas,plat,bas,iax(1),iax(2),itr
Co Outputs
Co   wk(0): distance between pairs
Co   wk(1..3): direction cosines
Cr Remarks
Cr
C ----------------------------------------------------------------
C     implicit none
      integer iax(5),nbas
      double precision plat(3,3),bas(3,nbas),wk(0:3)
      integer i,j

      do  101  i = 1, 3
        wk(i) = bas(i,iax(2))  - bas(i,iax(1))
      do  100  j = 1, 3
        wk(i) = wk(i) + plat(i,j)*iax(2+j)
  100 continue
  101 continue
      wk(0) = dsqrt(wk(1)**2 + wk(2)**2 + wk(3)**2)
      if (wk(0) .ne. 0) then
        wk(1) = wk(1)/wk(0)
        wk(2) = wk(2)/wk(0)
        wk(3) = wk(3)/wk(0)
      endif
      end

