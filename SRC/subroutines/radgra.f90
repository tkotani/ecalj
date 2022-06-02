subroutine radgra(a,b,nr,r,v,gradv)

  !- radial gradient
  ! ------------------------------------------------------------------
  !i Inputs:
  !i    a,b,nr,r(nr),v(nr)
  !o Outputs:
  !o    gradv(nr)
  !r Remarks:
  !r    makes the derivative of the function v defined in a mesh
  !r    of this kind:
  !r                r(i)=B (exp A(i-1)-1)
  ! ------------------------------------------------------------------
  !     implicit none
  ! Passed parameters
  integer :: nr
  double precision :: a,b,r(nr),v(nr),gradv(nr)
  ! Local parameters
  integer :: nm2,i

  ! Forward diffs for first and second point (Handbook,25.3.9 with 25.1.1)
  gradv(1) = ((6d0*v(2)+20d0/3d0*v(4)+1.2d0*v(6)) &
       -(2.45d0*v(1)+7.5d0*v(3)+3.75d0*v(5)+1d0/6d0*v(7)))/a
  gradv(2) = ((6d0*v(3)+20d0/3d0*v(5)+1.2d0*v(7)) &
       -(2.45d0*v(2)+7.5d0*v(4)+3.75d0*v(6)+1d0/6d0*v(8)))/a

  ! Five points formula  (25.3.6)
  nm2 = nr-2
  do  1  i = 3, nm2
     gradv(i) = ((v(i-2)+8*v(i+1)) - (8*v(i-1)+v(i+2)))/12/a
1 enddo

  ! Five points formula  (25.3.6)
  gradv(nr-1) = (-1d0/12d0*v(nr-4)+0.5d0*v(nr-3)-1.5d0*v(nr-2) &
       +5d0/6d0*v(nr-1)+0.25d0*v(nr))/a
  gradv(nr) = (0.25d0*v(nr-4)-4d0/3d0*v(nr-3)+3d0*v(nr-2) &
       -4d0*v(nr-1)+25d0/12d0*v(nr))/a
  ! Three points formula  (25.3.4)
  !     gradv(nr-1)=(v(nr)-v(nr-2))/2d0/a
  !     gradv(nr)=(v(nr-2)/2d0-2d0*v(nr-1)+1.5d0*v(nr))/a

  do  2  i = 1, nr
     gradv(i) = gradv(i)/(r(i)+b)
2 enddo
end subroutine radgra


