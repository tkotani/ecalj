subroutine latlim(pqlat,rmax,i1,i2,i3)
  !- Set limits in X Y Z direction
  ! ----------------------------------------------------------------
  !i Inputs
  !i   pqlat: primitive lattice vectors (real or reciprocal space)
  !i   rmax:  maximum length of connecting vector.
  !o Outputs
  !o   i1,i2,i3: all connecting vectors lie within these
  !o             multiples of the lattice vectors.
  !r Remarks
  !r Remarks
  !r   Define I_jk = \vec r_i \dot \vec r_j, where r is a lattice vector
  !r   and define a vector \vec v as some (non-integral) multiples
  !r   \vec \alpha of the \vec r_j, which has square length
  !r   v^2 = \sum_jk \alpha_j \alpha_k I_jk.  Then want to maximize the
  !r   \alpha_j subject to constraint that  v^2 = rmax^2.
  !r   The result turns out to be \alpha_j = rmax * (I^-1_jj)^1/2.
  !r
  !r   Bugs:
  !r   This routine only returns the integer part of
  !r   \alpha, when it should return next higher integer.
  ! ----------------------------------------------------------------
  !     implicit none
  ! Passed parameters
  double precision :: pqlat(3,3),rmax
  integer :: i1,i2,i3
  ! Local parameters
  double precision :: a(3,3),b(3,3),det
  integer :: makint,i
  makint(i) = int(rmax*dsqrt(b(i,i)) + .999d0)
  call dmpy(pqlat,1,3,pqlat,3,1,a,1,3,3,3,3)
  call dinv33(a,0,b,det)
  i1 = makint(1)
  i2 = makint(2)
  i3 = makint(3)
end subroutine latlim

