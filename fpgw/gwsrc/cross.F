      subroutine cross_x(a,b,c)
      implicit none
      real(8), intent(in) :: a(3), b(3)
      real(8), intent(out) :: c(3)
      c(1)=a(2)*b(3)-a(3)*b(2)
      c(2)=a(3)*b(1)-a(1)*b(3)
      c(3)=a(1)*b(2)-a(2)*b(1)
      return
      end

!> This is a replacement of dinv33 of Ferdi's GW  => dinv33(plat,1,qlat,det) 
!! the SAME as the one of dinv33 in extens.f in ferdi/lmto/extens.f
      subroutine dinv33x(plat,qlat)
      implicit none
      real(8),intent(in)::  plat(3,3)
      real(8),intent(out):: qlat(3,3)
      real(8):: det
      call cross_x(plat(1,2),plat(1,3), qlat     )
      call cross_x(plat(1,3),plat     , qlat(1,2))
      call cross_x(plat     ,plat(1,2), qlat(1,3))
      det  = sum( plat(1:3,1)*qlat(1:3,1) )
      qlat = qlat/det
      end

!> This is a replacement of dinv33 of Ferdi's GW  => dinv33(plat,1,qlat,det)
!! the SAME as the one of dinv33 in extens.f in ferdi/lmto/extens.f
      subroutine dinv33y(plat,qlat,det)
      implicit none
      real(8),intent(in)::  plat(3,3)
      real(8),intent(out):: qlat(3,3)
      real(8):: det
      call cross_x(plat(1,2),plat(1,3), qlat     )
      call cross_x(plat(1,3),plat     , qlat(1,2))
      call cross_x(plat     ,plat(1,2), qlat(1,3))
      det  = sum( plat(1:3,1)*qlat(1:3,1) )
      qlat = qlat/det
      end

!>- Inverts 3X3 matrix
      subroutine dinv33(matrix,iopt,inverse,det)
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
      integer iopt
      real(8), intent(in) :: matrix(3,3)
      real(8), intent(out) :: inverse(3,3),det
      real(8) :: ddot           ! BLAS external function
      call cross_x(matrix(1,2),matrix(1,3),inverse     )
      call cross_x(matrix(1,3),matrix     ,inverse(1,2))
      call cross_x(matrix     ,matrix(1,2),inverse(1,3))
      det = ddot(3,matrix,1,inverse,1)
      if (abs(det) ==0d0) call rx( 'dinv33: vanishing determinant')
      if (iopt == 0) inverse = transpose(inverse)
      inverse = inverse/det
      end

C   $Id: cross.f,v 1.2 1998/04/28 18:45:23 coa Exp $
c --- cross
      subroutine cross(a,b,c)
      implicit real*8 (a-h,p-z), integer(o)
      dimension a(3),b(3),c(3)
      c(1)=a(2)*b(3)-a(3)*b(2)
      c(2)=a(3)*b(1)-a(1)*b(3)
      c(3)=a(1)*b(2)-a(2)*b(1)
      return
      end

c$$$c --- tripl
c$$$      doubleprecision function tripl(a,b,c)
c$$$      implicit real*8 (a-h,p-z), integer(o)
c$$$c     dimension a(3),b(3),c(3)
c$$$c ... g77 needs this rewriting for opt > O1
c$$$      dimension a(0:2),b(0:2),c(0:2)
c$$$c     tripl=a(1)*b(2)*c(3)+a(2)*b(3)*c(1)+a(3)*b(1)*c(2)
c$$$c    .     -a(3)*b(2)*c(1)-a(2)*b(1)*c(3)-a(1)*b(3)*c(2)
c$$$c ... g77 needs this rewriting
c$$$      tmp = 0.d0
c$$$      do i = 0,2
c$$$        j = mod(i + 1, 3)
c$$$        k = mod(i + 2, 3)
c$$$        tmp = tmp + a(i) * (b(j)*c(k) - b(k)*c(j))
c$$$      enddo
c$$$      tripl = tmp
c$$$      end

c --- sprod
      doubleprecision function sprod(a,b)
      implicit real*8 (a-h,p-z), integer(o)
      dimension a(3),b(3)
      sprod=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
      return
      end
