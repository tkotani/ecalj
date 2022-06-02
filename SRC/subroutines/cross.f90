real(8) function tripl(a,b,c)
  !! == tripl (determinant of 3x3 matrix) ==
  !     implicit none
  double precision :: a(3),b(3),c(3)
  tripl = a(1)*b(2)*c(3) + a(2)*b(3)*c(1) + a(3)*b(1)*c(2) &
       -a(3)*b(2)*c(1) - a(2)*b(1)*c(3) - a(1)*b(3)*c(2)
END function tripl

subroutine cross(a,b,c)
  implicit none
  intent(in)  ::   a,b
  intent(out) ::       c
  real(8):: a(3),b(3),c(3)
  c(1)=a(2)*b(3)-a(3)*b(2)
  c(2)=a(3)*b(1)-a(1)*b(3)
  c(3)=a(1)*b(2)-a(2)*b(1)
  return
end subroutine cross

!> This is a replacement of dinv33 of Ferdi's GW  => dinv33(plat,1,qlat,det)
!! the SAME as the one of dinv33 in extens.f in ferdi/lmto/extens.f
subroutine minv33tp(  plat,qlat)
  implicit none
  real(8),intent(in)::  plat(3,3)
  real(8),intent(out):: qlat(3,3)
  real(8):: det
  call cross(plat(1,2),plat(1,3), qlat     )
  call cross(plat(1,3),plat     , qlat(1,2))
  call cross(plat     ,plat(1,2), qlat(1,3))
  det  = sum( plat(1:3,1)*qlat(1:3,1) )
  qlat = qlat/det
end subroutine minv33tp

!>- Inverts 3X3 matrix
subroutine minv33(matrix,inverse)
  !o Outputs
  !o   inverse, as modified according to iopt
  !o   det:      determinant
  implicit none
  !      integer:: iopt=0
  real(8), intent(in) :: matrix(3,3)
  real(8), intent(out) :: inverse(3,3)
  real(8) :: det,ddot
  call cross(matrix(1,2),matrix(1,3),inverse     )
  call cross(matrix(1,3),matrix     ,inverse(1,2))
  call cross(matrix     ,matrix(1,2),inverse(1,3))
  det = ddot(3,matrix,1,inverse,1)
  if (abs(det) ==0d0) call rx( 'minv33: vanishing determinant')
  inverse = transpose(inverse)
  inverse = inverse/det
end subroutine minv33
