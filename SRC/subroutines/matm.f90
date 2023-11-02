! routines for matrix multiplication and matrix inversion
subroutine matm(a,b,c,n1,n2,n3)
  integer, intent(in) :: n1,n2,n3
  complex(8), intent(in) :: a(n1,n2), b(n2,n3)
  complex(8), intent(out) :: c(n1,n3)
  c=matmul(a,b)
end subroutine matm
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
subroutine zgesvdnn(ngb,zzz, SS,UU,VT)! Sigular Value Decomp. zzz= matmul(UU,matmul(SS,VT))
  implicit none
  integer(4)::lwork,info,ngb,i
  complex(8):: zzz(ngb,ngb),UU(ngb,ngb),VT(ngb,ngb)
  real(8):: ss(ngb)
  real(8),allocatable:: rwork(:)
  complex(8),allocatable:: work(:),zw0bk(:,:),vtt(:,:)
  lwork=4*ngb
  allocate(zw0bk(ngb,ngb))
  allocate(work(LWORK),rwork(5*ngb)) !,VTT(ngb,ngb))
  zw0bk = zzz
  call zgesvd('A','A',ngb,ngb,zzz,ngb,SS,UU,ngb,VT,ngb,work,lwork,rwork,info)
  !      do i=1,ngb
  !         write(6,"(' i ss=',i4,' ', d13.5 )")i,SS(i)
  !    write(6,"(' i ss=',i4,'  ', d13.5,' ss0*ss=',d13.5 )")i,SS(i),ss(i)*ss0(ngb-i+1)
  !         vtt(i,:)=ss(i)*vt(i,:)
  !      enddo
  !      write(6,"('sumcheck zzz  zzz-uu*s*vt=',d13.5,d13.5)")
  !     &  sum(abs(zw0bk)), sum(abs(zw0bk - matmul(uu,vtt)))
  !      if(abs(sum(abs(zw0bk - matmul(uu,vtt))))>1d-8*sum(abs(zw0bk)))
  !     &  stop 'sumcheck zzz  zzz-uu*s*vt= error'
  !      deallocate(vtt)
end subroutine zgesvdnn

