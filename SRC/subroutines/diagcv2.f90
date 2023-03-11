subroutine diagcvh2(hh,ngb,eb)!  diagonalizes and returns eigenstates
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
  abstol= 2d0*DLAMCH('S')
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
subroutine rss(nn,omat,eb,zz,ierr)!get eigenvalues and eigenvectors for real-symmetric funciton.
  integer:: ierr,lwork,nn
  real(8):: omat(nn,nn),omatin(nn,nn),eb(nn),zz(nn,nn)
  real(8),allocatable:: work(:)
  lwork=max(1,3*nn-1)
  allocate(work(max(1,lwork)))
  omatin=omat
  call dsyev('V','U',nn,omatin,nn,eb,work,lwork,ierr)
  zz=omatin
  deallocate(work)
end subroutine rss
!     For diagonalizing Non-Hermite matrix (2017/09/22, Okumura)
!     For Hermite matrix, use 'diagcv2.F'
subroutine diagcvuh3(uh,n,ev)
  ! uh:Non-Hermite matrix, n:dimenstion, ev:eigenvalue(COMPLEX)
  implicit none
  integer(4):: nmx,nev,i,n
  complex(8):: uh(n,n),hhx(n,n),oo(n,n),zz(n,n)
  complex(8):: ev(n) !complex eigenvalue
  hhx=uh
  nmx=n
  oo = 0d0
  do i=1,n
     oo(i,i) = 1d0
  enddo
  call diagcvz(oo,hhx,zz,n, ev,nmx,1d99,nev)
end subroutine diagcvuh3
subroutine diagcvz(s,uh,t,n,evl,nmx,emx,nev)
  !  diagonalizes and returns the nev lowest eigenstates
  !  eigenvecs are returned for i.lt.nmx. emx is dummy
  implicit none
  integer :: n,nmx,nev
  double precision :: emx,abstol,dlamch
  complex(8)::evl(n)
  complex*16 s(n,n),uh(n,n),t(n,n)
  complex(8),allocatable:: WORK(:)
  real(8),allocatable:: RWORK(:)
  !      complex*16:: vl,vr       !bugfix (whis was not decleared) jan2013
  complex*16:: vl(1),vr(1)
  integer:: LWORK,info
  abstol= 2d0*DLAMCH('S')
  !---find optimum work size
  allocate( WORK(1))
  lwork = -1
  call ZGEEV( 'N', 'N', n, uh, n, evl, vl, 1,  vr, 1, &
       WORK, LWORK, RWORK, INFO)
  LWORK = WORK(1) + 1
  deallocate(WORK)
  !      LWORK = max(1,2*n-1)    ! This caused a  problem---why?
  !      LWORK = max(1,2*n-1)+1  ! This caused no problem---why?
  allocate( WORK(LWORK),RWORK(7*n))
!!! not return eigenvector (vl and vr)
  call ZGEEV( 'N', 'N', n, uh, n, evl, vl, 1,  vr, 1, &
       WORK, LWORK, RWORK, INFO)
  if(INFO/=0) then
     print *, 'See http://www.netlib.org/cgi-bin/netlibget.pl', &
          '/lapack/complex16/zgeev.f'
     print *, 'ZGEEV error info=',info
  endif
  deallocate(WORK,RWORK)
  return
end subroutine diagcvz
!     For Hermite matrix, use 'diagcv2.F'
subroutine diagcvuh3_vec(uh,n,ev,vr)
  ! uh:Non-Hermite matrix, n:dimenstion, ev:eigenvalue(COMPLEX), vr:right-eigenvecter
  implicit none
  integer(4):: nmx,nev,i,n
  complex(8):: uh(n,n),hhx(n,n),oo(n,n),zz(n,n)
  complex(8):: ev(n) !complex eigenvalue
  complex(8):: vr(n,n)
  hhx=uh
  nmx=n
  oo = 0d0
  do i=1,n
     oo(i,i) = 1d0
  enddo
  call diagcvz_r(oo,hhx,zz,n, ev,nmx,1d99,nev,vr)
end subroutine diagcvuh3_vec
subroutine diagcvz_r(s,uh,t,n,evl,nmx,emx,nev,vr)
  !  diagonalizes and returns the nev lowest eigenstates
  !  eigenvecs are returned for i.lt.nmx. emx is dummy
  implicit none
  integer :: n,nmx,nev
  double precision :: emx,abstol,dlamch
  complex(8)::evl(n)
  complex*16 s(n,n),uh(n,n),t(n,n)
  complex(8),allocatable:: WORK(:)
  real(8),allocatable:: RWORK(:)
  !      complex*16:: vl,vr       !bugfix (whis was not decleared) jan2013
  complex*16:: vl(1) !vr(1)
  complex(8):: vr(n,n)
  integer:: LWORK,info
  abstol= 2d0*DLAMCH('S')
  !---find optimum work size
  allocate( WORK(1))
  lwork = -1
  call ZGEEV( 'N', 'V', n, uh, n, evl, vl, 1,  vr, n,    WORK, LWORK, RWORK, INFO)
  LWORK = WORK(1) + 1
  deallocate(WORK)
  !      LWORK = max(1,2*n-1)    ! This caused a  problem---why?
  !      LWORK = max(1,2*n-1)+1  ! This caused no problem---why?
  allocate( WORK(LWORK),RWORK(7*n))
!!! return eigenvector (vl and vr)
  call ZGEEV( 'N', 'V', n, uh, n, evl, vl, 1,  vr, n,     WORK, LWORK, RWORK, INFO)
  if(INFO/=0) then
     print *, 'See http://www.netlib.org/cgi-bin/netlibget.pl', &
          '/lapack/complex16/zgeev.f'
     print *, 'ZGEEV error info=',info
  endif
  deallocate(WORK,RWORK)
  return
end subroutine diagcvz_r
