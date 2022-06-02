subroutine dinv33(matrix,iopt,invrse,det)
  !- Inverts 3x3 matrix
  ! ----------------------------------------------------------------
  !i Inputs
  !i   matrix:  matrix to be inverted
  !i   iopt:  if 0, usual inverse
  !i             1, transpose of inverse
  !i             2, 2*pi*inverse
  !i             3, 2*pi*transpose of inverse
  !o Outputs
  !o   invrse    see iopt
  !o   det:      determinant, or det/2*pi (sign ok ??)
  !r Remarks
  !r   To generate reciprocal lattice vectors, call dinv33(plat,3,rlat)
  ! ----------------------------------------------------------------
  !     implicit none
  integer :: iopt,i,j
  double precision :: matrix(3,3),invrse(3,3),det,ddot
  double precision :: xx
  call cross(matrix(1,2),matrix(1,3),invrse     )
  call cross(matrix(1,3),matrix     ,invrse(1,2))
  call cross(matrix     ,matrix(1,2),invrse(1,3))
  det = ddot(3,matrix,1,invrse,1)
  if (det == 0d0) call rx('INV33: vanishing determinant')
  if (iopt >= 2) det = det/(8*datan(1d0))
  if (mod(iopt,2) == 0) then
     do    i = 1, 3
        do  j = i+1, 3
           xx = invrse(i,j)
           invrse(i,j) = invrse(j,i)
           invrse(j,i) = xx
        enddo
     enddo
  endif
  call dscal(9,1/det,invrse,1)
end subroutine dinv33

