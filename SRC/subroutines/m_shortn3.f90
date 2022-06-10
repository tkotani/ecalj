module m_shortvec !== Find shortest vector in modulo of rlat ===
  public:: shortvec,shortvecinitialize
contains
  ! NOTE: shortn3 is better than shorbz. we will have to replace shorbz with shortn3.
  ! NOTE: In advance, we will need to check speed and convenience.
  subroutine shortvec(pin,rlatp,xmx2,noutmx,nout,nlatout)
    !! To call shortn3 for given rlat, 
    !! we have to call shorn33initialize to set rlatp and xmx2
    !!
    ! i pin is the fractional coodinate on rlat.
    ! i rlatp,xmx2 are passed from shortvecinitialize
    !!   rlatp(i,j)= sum( rlat(:,i)*rlat(:,j) )
    !!   rlat(3,i) i-th vertor for modulo
    ! i noutmax: upper limit of nlatout
    ! o nout
    ! o nlatout
    !!  Shortest vectors are
    !!    pin+nlatout(:,ix), where ix=1:nout, is the shortest vectors.
    !!    We may have multiple nlatout (# is nout).
    !!
    !!==========================================================================
    implicit none
    integer:: noutmx
    integer:: nlatout(3,noutmx)
    integer:: nmax(3),nknknk,ik1,ik2,ik3,nout,nk,ik,i,j
    real(8):: rmax2,pin(3),eps=1d-8,rlat(3,3),xmax2(3),rr(3),rmin,nrmax(3)
    integer,allocatable:: nlat0(:,:)
    real(8),allocatable:: rnorm(:)
    real(8):: rlatp(3,3),xmx2(3)
    rmax2 = sum(pin*matmul(rlatp,pin)) + eps  ! eps is to make degeneracy safe.
    nrmax(:) =  sqrt(rmax2*xmx2(:))+abs(pin(:)) ! range of ix
    nmax =  nrmax
    ! we are looking for shortest vectors
    ik=0
    rmin=1d9
    nknknk= (2*nmax(1)+1)*(2*nmax(2)+1)*(2*nmax(3)+1)
    allocate( nlat0(3, nknknk), rnorm(nknknk) )
    do ik1=-nmax(1),nmax(1)
       do ik2=-nmax(2),nmax(2)
          do ik3=-nmax(3),nmax(3)
             ik=ik+1
             nlat0(:,ik) = (/ik1,ik2,ik3/)
             rr= pin + nlat0(:,ik)
             rnorm(ik) = sum(rr*matmul(rlatp,rr))
             if(rnorm(ik)<rmin) rmin=rnorm(ik)
          enddo
       enddo
    enddo
    nk=ik
    ! rint *,'nk rmin=',nk,rmin
    nout=0
    do ik=1,nk
       rr= pin + nlat0(:,ik)
       ! rint *,'ik rr   =',ik,rr
       ! rint *,'ik rnorm=',ik,rnorm(ik)
       if(rnorm(ik)<rmin+eps) then
          nout=nout+1
          if(nout>noutmx) stop 'shortn3: enlarge noutmx'
          nlatout(:,nout)=nlat0(:,ik)
          ! rint *,'ik nlat0',nlat0(:,ik)
       endif
    enddo
    ! rite(6,"('pin=',3f8.3,' nmax=',3i4,' nout=',i3)")pin, nmax(1:3),nout
    deallocate(rnorm,nlat0)
    return
  end subroutine shortvec
  !------------------------------------
  subroutine shortvecinitialize(rlat,rlatp,xmx2)
    !!== Set translation vactors rlat(:,i),i=1,3 ==
    ! i rlat
    ! o rlatp,xmx2: these are passed to shortn3
    !     !=============================================
    integer:: i,j
    real(8):: rlat(3,3)
    real(8):: rlatp(3,3),xmx2(3)
    do i=1,3
       do j=1,3
          rlatp(i,j) = sum(rlat(:,i)*rlat(:,j))
       enddo
    enddo
    call ellipsoidxmax(rlatp,xmx2)
  end subroutine shortvecinitialize
end module m_shortvec

module m_shortn3!gives shortest vector modulo of rlat
  ! We need to set rlat into this module in advance to call shortn3
  implicit none
  public:: shortn3_initialize, shortn3
  integer,parameter,private:: noutmx=48
  integer,public:: nout,nlatout(3,noutmx)
  real(8),private:: rlatp(3,3),xmx2(3)
contains
  subroutine shortn3_initialize(rlat)
    use m_shortvec,only: shortvecinitialize
    real(8):: rlat(3,3)
    call shortvecinitialize(rlat,rlatp,xmx2)
  end subroutine shortn3_initialize
  subroutine shortn3(pin) !pin is in fractional coordinate in rlat(:,i)
    use m_shortvec,only: shortvec
    real(8):: pin(3)
    call shortvec(pin,rlatp,xmx2,noutmx,nout,nlatout)
  end subroutine shortn3
end module m_shortn3

subroutine ellipsoidxmax(nn, xmx2)
  !!== Maximum value for x_i for ellipsoid ==
  !!  Ellipsoid is given as 1d0 = sum x_i nn(i,j) x_j.
  ! i nn(3,3)
  ! o xmx2(i)  Maximum of x_i**2
  !!==========================================
  implicit none
  real(8),target:: nn(3,3)
  real(8):: v2(2),ainv(2,2), rmax2, xmx2(3),det,fac,nv2(2)
  real(8),pointer::n11,n12,n13,n21,n22,n23,n31,n32,n33
  n11=>nn(1,1)
  n12=>nn(1,2)
  n13=>nn(1,3)
  n21=>nn(2,1)
  n22=>nn(2,2)
  n23=>nn(2,3)
  n31=>nn(3,1)
  n32=>nn(3,2)
  n33=>nn(3,3)
  det= n22*n33-n23*n32
  ainv(1,1)=  n33/det
  ainv(2,2)=  n22/det
  ainv(1,2)= -n23/det
  ainv(2,1)= -n32/det
  nv2  = (/n12,n13/)
  fac = n11-sum(nv2 *matmul(ainv,nv2))
  ! rint *,'ainv=',ainv
  ! rint *,'ainv*nv=',matmul(ainv,nv2)
  xmx2(1) = 1d0/fac

  det= n33*n11-n31*n13
  ainv(1,1)=  n11/det
  ainv(2,2)=  n33/det
  ainv(1,2)= -n31/det
  ainv(2,1)= -n13/det
  nv2  = (/n23,n21/)
  fac = n22-sum(nv2 *matmul(ainv,nv2))
  xmx2(2) = 1d0/fac

  det= n11*n22-n12*n21
  ainv(1,1)=  n22/det
  ainv(2,2)=  n11/det
  ainv(1,2)= -n12/det
  ainv(2,1)= -n21/det
  nv2  = (/n31,n32/)
  fac = n33-sum(nv2 *matmul(ainv,nv2))
  xmx2(3) = 1d0/fac
end subroutine ellipsoidxmax

! !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! module m_shortn4
!   implicit none
!   public shortn4,shortn4_initialize
!   integer,private,parameter:: noutmx=48
!   integer,public:: nout,nlatout(3,noutmx)
  
!   real(8),private:: rlatp(3,3),xmx2(3)
!   logical,private:: init=.true.
! contains
!   subroutine shortn4(pin)
!     use m_shortvec,only: shortvec
!     use m_lattic,only: qlat=>lat_qlat
!     real(8):: pin(3)
!     call shortvec(pin,rlatp,xmx2,noutmx,nout,nlatout)
!   end subroutine shortn4
!   subroutine shortn4_initialize(rlat)
!     use m_shortvec,only: shortvecinitialize
!     real(8):: rlat(3,3)
!     call shortvecinitialize(rlat,rlatp,xmx2)
!   end subroutine shortn4_initialize
! end module m_shortn4




