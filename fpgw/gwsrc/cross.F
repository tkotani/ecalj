c --- external product
      subroutine cross_x(a,b,c)
      implicit real*8 (a-h,p-z), integer(o)
      dimension a(3),b(3),c(3)
      c(1)=a(2)*b(3)-a(3)*b(2)
      c(2)=a(3)*b(1)-a(1)*b(3)
      c(3)=a(1)*b(2)-a(2)*b(1)
      return
      end
c
      subroutine dinv33x(plat,qlat)
C- This is a replacement of dinv33 of Ferdi's GW  => dinv33(plat,1,qlat,det) --------------
Cr THIS IS the SAME as the one of dinv33 in extens.f in ferdi/lmto/extens.f
      implicit none
      double precision plat(3,3),qlat(3,3),det
      call cross_x(plat(1,2),plat(1,3), qlat     )
      call cross_x(plat(1,3),plat     , qlat(1,2))
      call cross_x(plat     ,plat(1,2), qlat(1,3))
      det  = sum( plat(1:3,1)*qlat(1:3,1) )
      qlat = qlat/det
      end

      subroutine dinv33y(plat,qlat,det)
C- This is a replacement of dinv33 of Ferdi's GW  => dinv33(plat,1,qlat,det) --------------
Cr THIS IS the SAME as the one of dinv33 in extens.f in ferdi/lmto/extens.f
      implicit none
      double precision plat(3,3),qlat(3,3),det
      call cross_x(plat(1,2),plat(1,3), qlat     )
      call cross_x(plat(1,3),plat     , qlat(1,2))
      call cross_x(plat     ,plat(1,2), qlat(1,3))
      det  = sum( plat(1:3,1)*qlat(1:3,1) )
      qlat = qlat/det
      end

c-Taken from Ferdi's GW  -----------------------------------------------------
      subroutine dinv33(matrix,iopt,inverse,det)
C- Inverts 3X3 matrix
C ----------------------------------------------------------------
Ci Inputs
Ci   inverse: input matrix
Ci   iopt:  if 0, usual inverse
Ci             1, transpose of inverse
Co Outputs
Co   inverse, as modified according to iopt
Co   det:      determinant
C ----------------------------------------------------------------
      implicit none
      integer iopt,i,j
      double precision matrix(3,3),inverse(3,3),det,ddot
      if(iopt<0.or.iopt>1) stop 'dinv33:wrong iopt'
      call cross_x(matrix(1,2),matrix(1,3),inverse     )
      call cross_x(matrix(1,3),matrix     ,inverse(1,2))
      call cross_x(matrix     ,matrix(1,2),inverse(1,3))
      det = ddot(3,matrix,1,inverse,1)
      if (abs(det) ==0d0) stop 'dinv33: vanishing determinant'
      if (iopt == 0) inverse = transpose(inverse)
      inverse = inverse/det
c      double precision xx
c      if (iopt .ge. 2) det = det/(8*datan(1d0))
c      if (mod(iopt,2) == 0) then
c        do  10  i = 1, 3
c          do  10  j = i+1, 3
c          xx = inverse(i,j)
c          inverse(i,j) = inverse(j,i)
c          inverse(j,i) = xx
c   10   continue
c      endif
c     call dscal(9,1/det,inverse,1)
c      return
      end
