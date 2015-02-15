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
