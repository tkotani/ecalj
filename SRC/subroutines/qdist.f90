subroutine qdist(q,q1,ux,uy,uz,gam)
  !  (from msm) new version of 03.10.89. u=(ux,uy,uz) gives direction,
  !  gam is multiplier in real space along u, volume is conserved.
  !     implicit none
  double precision :: q(3),q1(3),ux,uy,uz,gam
  double precision :: g,u2,a,xxx,yyy
  g = 1d0/gam
  u2 = ux*ux+uy*uy+uz*uz
  a = (ux*q(1)+uy*q(2)+uz*q(3))/u2
  xxx = 1d0/dsqrt(dabs(g))
  if (g < 0d0) xxx = 1
  yyy = a*(g-xxx)
  q1(1) = xxx*q(1)+yyy*ux
  q1(2) = xxx*q(2)+yyy*uy
  q1(3) = xxx*q(3)+yyy*uz
end subroutine qdist
subroutine rdist(v,v1,ux,uy,uz,g)
  !     implicit none
  double precision :: v(3),v1(3),ux,uy,uz,g,u2,a,xxx,yyy
  u2 = ux*ux+uy*uy+uz*uz
  a = (ux*v(1)+uy*v(2)+uz*v(3))/u2
  xxx = 1d0/dsqrt(dabs(g))
  if (g < 0d0) xxx=1d0
  yyy = a*(g-xxx)
  v1(1) = xxx*v(1)+yyy*ux
  v1(2) = xxx*v(2)+yyy*uy
  v1(3) = xxx*v(3)+yyy*uz
end subroutine rdist
subroutine rdistn(a1,a2,na,gx,gy,gz,gt)
  !- distort na real-space vectors in array a1 into array a2
  !     implicit none
  integer :: na,ia
  double precision :: a1(3,na),a2(3,na),gx,gy,gz,gt
  do 10 ia = 1, na
     call rdist(a1(1,ia),a2(1,ia),gx,gy,gz,gt)
10 enddo
end subroutine rdistn
subroutine qdistn(a1,a2,na,gx,gy,gz,gt)
  !- distort na q-space vectors in array a1 into array a2
  !     implicit none
  integer :: na,ia
  double precision :: a1(3,na),a2(3,na),gx,gy,gz,gt
  do 10 ia = 1,na
     call qdist(a1(1,ia),a2(1,ia),gx,gy,gz,gt)
10 enddo
end subroutine qdistn

