! Simple mathematic functions

! ccccccccccccccccccccccccccccccccccc
! Factrization

integer function myfact(n)
  implicit none
  integer,intent(in) :: n
  integer :: i
  myfact=1
  do i=1,n
     myfact=myfact*i
  enddo
end function myfact
! ccccccccccccccccccccccccccccccccccc
! Permutation

integer function myperm(n,m)
  implicit none
  integer,intent(in) :: n,m
  integer :: i
  myperm=1
  do i=0,m-1
     myperm=myperm*(n-i)
  enddo
end function myperm
! ccccccccccccccccccccccccccccccccccc
! Combination

integer function mycomb(n,m)
  implicit none
  integer,intent(in) :: n,m
  integer :: i,nmm
  !  function
  integer :: myperm,myfact

  nmm=m
  if (n-m < m) nmm=n-m

  mycomb=myperm(n,nmm)/myfact(nmm)
end function mycomb
! cccccccccccccccccccccccccccccccccc
! Calculate m-th derivative of x^n
! (d/dx^m) x^n = n*(n-1)*...*(n-m+1)*x^(n-m)

double precision function mydif(x,n,m)
  implicit none
  double precision,intent(in) :: x
  integer,intent(in) :: n,m
  ! function
  integer :: myperm

  if (m > n) then
     mydif=0
     return
  endif

  mydif=dble(myperm(n,m))*x**(n-m)
end function mydif

! Simple matrix routines

! ccccccccccccccccccccccccc
! cross product V1 X V2 =>U
subroutine mycross(V1,V2,U)
  implicit none
  double precision :: V1(3),V2(3),U(3)
  U(1)=V1(2)*V2(3)-V1(3)*V2(2)
  U(2)=V1(3)*V2(1)-V1(1)*V2(3)
  U(3)=V1(1)*V2(2)-V1(2)*V2(1)
end subroutine mycross
! ccccccccccccccccccccccccc
! multiply matrix A(m,n) and vector V(n) => AV(m)
subroutine mymatvec(A,V,AV,m,n)
  implicit none
  integer :: m,n
  double precision :: A(m,n),V(n),AV(m)
  ! local
  integer :: i,j

  AV(1:m)=0.0d0
  do j=1,n
     do i=1,m
        AV(i)=AV(i)+A(i,j)*V(j)
     enddo
  enddo
end subroutine mymatvec
! ccccccccccccccccccccccccc
! multiply matrix A(l,m) and B(m,n) => C(l,n)
subroutine mymatmat(A,B,C,l,m,n)
  implicit none
  integer :: l,m,n
  double precision :: A(l,m),B(m,n),C(l,n)
  ! local
  double precision ::  At(m,l)
  integer :: i,j,k

  call mytranspose(A,At,l,m)
  C(1:l,1:n)=0.0d0

  do j=1,n
     do i=1,l
        do k=1,m
           C(i,j)=C(i,j)+At(k,i)*B(k,j)
        enddo
     enddo
  enddo
end subroutine mymatmat
! cccccccccccccccccccccccc
! determinant of 3x3 matrix A
subroutine mydet3(A,det3)
  implicit none
  double precision ::  A(3,3)
  double precision ::  det3

  det3=A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+ &
       A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+ &
       A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
end subroutine mydet3
! cccccccccccccccccccccccc
subroutine mytranspose(A,At,m,n)
  implicit none

  ! input
  integer :: m,n
  double precision ::  A(m,n)
  ! output
  double precision ::  At(n,m)
  ! local
  integer :: i,j

  do j=1,n
     do i=1,m
        At(j,i)=A(i,j)
     enddo
  enddo
end subroutine mytranspose
! cccccccccccccccccccccccc
subroutine myinv3(A,invA)
  implicit none
  ! input
  double precision ::  A(3,3)
  ! output
  double precision ::  invA(3,3)
  ! local
  double precision ::  eps
  parameter (eps=1.0d-20)
  double precision ::  detA,At(3,3) ! the transpose of A
  double precision ::  at1(3),at2(3),at3(3)

  call mydet3(A,detA)
  if (abs(detA) <= eps) then
     write(*,*) 'Error in myinv3: detA<eps'
     stop 'Error in myinv3: detA<eps'
  endif
  call mytranspose(A,At,3,3)

  at1(1:3)=At(1:3,1)
  at2(1:3)=At(1:3,2)
  at3(1:3)=At(1:3,3)

  call mycross(at2,at3,invA(1:3,1))
  call mycross(at3,at1,invA(1:3,2))
  call mycross(at1,at2,invA(1:3,3))

  invA(1:3,1:3)=invA(1:3,1:3)/detA

end subroutine myinv3
! ccccccccccccccccccccccccc

