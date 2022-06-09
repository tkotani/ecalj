subroutine diagcvh2(hh,ngb,eb)
  implicit none
  integer(4):: nmx,nev,i,ngb
  complex(8):: hh(ngb,ngb),hhx(ngb,ngb),oo(ngb,ngb),zz(ngb,ngb)
  real(8):: eb(ngb)
  hhx=hh
  nmx=ngb
  oo = 0d0
  do i=1,ngb
     oo(i,i) = 1d0
  enddo
  call diagcv(oo,hhx,zz,ngb, eb,nmx,1d99,nev)
  !      print *,' diagcvv: ngb,nev=',ngb,nev
  !      do i=1,nev
  !        write(6,'(i4,d23.16)')i, eb(i)
  !      enddo
end subroutine diagcvh2

subroutine diagcv(s,h,t,n,evl,nmx,emx,nev)
  !  diagonalizes and returns the nev lowest eigenstates
  !  eigenvecs are returned for i.lt.nmx. emx is dummy
  implicit none
  integer :: n,nmx,nev
  double precision :: emx,evl(n),abstol,dlamch
  complex*16 s(n,n),h(n,n),t(n,n)
  logical :: lx,lnv
  integer :: i,ipr,iprint,j,oww
  complex(8),allocatable:: work(:),s_(:,:),h_(:,:)
  integer(4),allocatable:: iwork(:),ifail(:)
  real(8),allocatable:: rwork(:)
  real(8):: vl,vu  !bugfix (whis was not decleared) jan2013
  integer:: lwork,info

  !      print *,'diagcv in diagcv2:'
  abstol= 2d0*DLAMCH('S')
  !      print *,' xxx abstol=',abstol

  !---find optimum work size
  allocate( work(1))
  LWORK = -1
  call ZHEGVX( 1, 'V', 'I', 'U', n, h, n, s, n, &
       VL, VU, 1, nmx, abstol, nev, evl, t, n, &
       work,LWORK, RWORK, IWORK, IFAIL, INFO )
  lwork = work(1) + 1
  deallocate(work)

  !      LWORK = max(1,2*n-1)    ! This caused a  problem---why?
  !      LWORK = max(1,2*n-1)+1  ! This caused no problem---why?

  allocate( work(LWORK),rwork(7*n),iwork(5*n),ifail(n))
  !      print *,' goto zhegvx'
  call ZHEGVX( 1, 'V', 'I', 'U', n, h, n, s, n, &
       VL, VU, 1, nmx, abstol, nev, evl, t, n, &
       work,LWORK, RWORK, IWORK, IFAIL, INFO )
  if(INFO/=0) then
     print *, 'See http://www.netlib.org/cgi-bin/netlibget.pl', &
          '/lapack/complex16/zhegvx.f'
     print *, 'ZHEGVX error info=',info
     print *, 'ZHEGVX error info=',ifail
  endif
  deallocate(work,rwork,iwork,ifail)
  return
end subroutine diagcv

