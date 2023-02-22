!> interface for matrix multiplication  c=a*b -------------------------
subroutine matm(a,b,c,n1,n2,n3)
  integer, intent(in) :: n1,n2,n3
  complex(8), intent(in) :: a(n1,n2), b(n2,n3)
  complex(8), intent(out) :: c(n1,n3)
  c=matmul(a,b)
end subroutine matm
subroutine matmaw(a,b,c,n1,n2,n3,ww)
  integer, intent(in) :: n1,n2,n3
  complex(8), intent(in) :: a(n1,n2), b(n2,n3)
  real(8), intent(in) :: ww(n1)
  complex(8), intent(out) :: c(n1,n3)
  complex(8):: aa(n1,n2)
  integer:: ix
  do ix=1,n2
     aa(:,ix) = ww(:)*a(:,ix)
  enddo
  c= c+ matmul(aa,b)
end subroutine matmaw

function matcinvf(b) result(a)
  !!== Test routine for Inversion ==
  implicit none
  integer :: info,n,n2(2)
  integer,allocatable :: ipiv(:)
  complex(8):: b(:,:)
  complex(8),allocatable:: a(:,:)
  complex(8),allocatable:: work(:)
  n2= SHAPE(b)
  n=n2(1)
  allocate(a,source=b) !call zcopy(n,b,1,a,1)
  call zgetrf(n,n,a,n,ipiv,info)
  if(info/=0) then
     write(6,*)' matcinvf: zegtrf info=',info
     call rx( ' matcinvf: zegtrf ')
  endif
  allocate(work(n*n))
  call zgetri(n,a,n,ipiv,work,n*n,info)
  deallocate(work)
  if(info/=0) then
     write(6,*)'matcinv: zegtri info=',info
     call rx( 'matcinv: zegtri ')
  endif
end function matcinvf

subroutine matcinv(n,a)
  implicit none
  integer :: n, info, ipiv(n)
  complex(8):: a(n,n)
  complex(8),allocatable:: work(:)
  call zgetrf(n,n,a,n,ipiv,info)
  if(info/=0) then
     print *,' matcinv: zegtrf info=',info
     call rx( ' matcinv: zegtrf ')
  endif
  allocate(work(n*n))
  call zgetri(n,a,n,ipiv,work,n*n,info)
  deallocate(work)
  if(info/=0) then
     print *,'matcinv: zegtri info=',info
     call rx( 'matcinv: zegtri ')
  endif
end subroutine matcinv
subroutine matinv(n,a)
  implicit none
  integer :: n, info, ipiv(n)
  real(8):: a(n,n)
  real(8),allocatable:: work(:)
  call dgetrf(n,n,a,n,ipiv,info)
  if(info/=0) then
     print *,' matinv: degtrf info=',info
     call rx( ' matinv: degtrf ')
  endif
  allocate(work(n*n))
  call dgetri(n,a,n,ipiv,work,n*n,info)
  deallocate(work)
  if(info/=0) then
     print *,'matinv: degtri info=',info
     call rx( 'matinv: degtri ')
  endif
end subroutine matinv
subroutine matinv2(n,a,info)
  implicit none
  integer :: n, info, ipiv(n)
  real(8):: a(n,n)
  real(8),allocatable:: work(:)
  call dgetrf(n,n,a,n,ipiv,info)
  if(info/=0) then
     print *,' matinv: degtrf info=',info
     return
  endif
  allocate(work(n*n))
  call dgetri(n,a,n,ipiv,work,n*n,info)
  deallocate(work)
  if(info/=0) then
     print *,'matinv: degtri info=',info
     return
  endif
end subroutine matinv2

