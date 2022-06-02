module m_shortn3
  !!== Find shortest vector in modulo of rlat ===
  public:: shortn3,shortn3_initialize
  real(8),private:: rlatp(3,3),xmx2(3)
contains
  ! NOTE: shortn3 is better than shorbz. we will have to replace shorbz with shortn3.
  ! NOTE: In advance, we will need to check speed and convenience.
  subroutine shortn3(pin,noutmx, nout,nlatout)
    !! To call shortn3 for given rlat,
    !! we have to call shorn3_initialize in advane to obtain rlatp and xmx2, which are passed to shortn3.
    !!
    ! i pin is on the rlat coodinate.
    ! i rlatp,xmx2 are passed from shortn3_initialize
    !!  rlatp(i,j)= sum( rlat(:,i)*rlat(:,j) )
    !!  rlat(3,i) i-th vertor for modulo
    ! o noutmax: upper limit of nlatout
    ! o nout
    ! o nlatout
    !!  Shortest vectors are
    !!  pin + nlatout \in  [ pin + any integer linear-combination of (rlat(:,1),rlat(:,2),rlat(:,3)) ].
    !!
    !!  length = sum_i sum_j pin(i)*rlatp(i,j)*pin(j)
    !!  pin+nlatout(:,ix), where ix=1:nout, is the shortest vectors. We may have multiple nlatout (# is nout).
    !!
    !! Takao think shorbz will be almost OK, but not perfect as an algolism.
    !! Takao think all shorbz should be replaced by shortn3 in future.
    !!==========================================================================
    implicit none
    integer:: nmax(3),nknknk,ik1,ik2,ik3,nout,nk,ik,i,j
    real(8):: rmax2,pin(3),eps=1d-8,rlat(3,3),xmax2(3),rr(3),rmin,nrmax(3)
    integer,allocatable:: nlat0(:,:)
    real(8),allocatable:: rnorm(:)
    integer:: noutmx
    integer:: nlatout(3,noutmx)
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
  end subroutine shortn3
  !------------------------------------
  subroutine shortn3_initialize(rlat)
    !!== Set translation vactors rlat(:,i),i=1,3 ==
    ! i rlat
    ! o rlatp,xmx2: these are passed to shortn3
    !     !=============================================
    integer:: i,j
    real(8):: rlat(3,3)
    do i=1,3
       do j=1,3
          rlatp(i,j) = sum(rlat(:,i)*rlat(:,j))
       enddo
    enddo
    call ellipsoidxmax(rlatp,xmx2)
  end subroutine shortn3_initialize
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
end module m_shortn3
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


subroutine shorbz(p,pout,plat,qlat)
  !- Shortens vector to equivalent in first Brillouin zone.
  ! ----------------------------------------------------------------
  !i Inputs:
  !i   plat,qlat lattice vectors and inverse
  !i   p         vector to shorten
  !o Outputs:
  !o   pout      shortened p.   pout and p can point to the same address.
  !r Remarks
  !r   Switch around plat,qlat to shorten reciprocal space vectors.
  !r   Jan 1997 Adapted from shorps to fix bug:  Example:
  !r   plat=  -0.5  0.5  1.7517  0.5  -0.5  1.7517  0.5  0.5  -1.7517
  !r   p= 0.0 -0.5 -1.26384
  !r   Should get pout -> 0.5 0.0 0.48786, not -0.5 1.0 0.48786.
  ! ----------------------------------------------------------------
  !     implicit none
  double precision :: p(3),pout(3),plat(3,3),qlat(3,3),x(3),x0,xx,a2,ap
  double precision :: tol
  parameter (tol=-1d-10)
  integer :: i,j,m,j2min,j3min,j1,j2,j3

  ! --- Reduce to unit cell centered at origin ---
  do  1  i = 1, 3
     ! ... x is projection of pin along plat(i), with multiples of p removed
     x0 = p(1)*qlat(1,i)+p(2)*qlat(2,i)+p(3)*qlat(3,i)
     xx = idnint(x0)
     x(i) = x0-xx
1 enddo
  ! ... pout is x rotated back to Cartesian coordinates
  do  2  m = 1, 3
     pout(m) = x(1)*plat(m,1)+x(2)*plat(m,2)+x(3)*plat(m,3)
2 enddo

  ! --- Try shortening by adding +/- basis vectors ---
15 continue
  do  j1 =  0, 1
     j2min = -1
     if (j1 == 0) j2min = 0
     do  j2 = j2min, 1
        j3min = -1
        if (j1 == 0 .AND. j2 == 0) j3min = 0
        do  j3 = j3min, 1

           !     ... (-1,0,1) (plat(1) + (-1,0,1) plat(2)) + (-1,0,1) plat(3))
           ! takao think this alglorism cause problems (=wrong) for very anisotropic cases.

           do  17  i = 1, 3
              x(i) = plat(i,1)*j1 + plat(i,2)*j2 + plat(i,3)*j3
17         enddo
           a2 = x(1)*x(1) + x(2)*x(2) + x(3)*x(3)
           ap = pout(1)*x(1) + pout(2)*x(2) + pout(3)*x(3)
           j = 0
           if (a2 + 2*ap < tol) j = 1
           if (a2 - 2*ap < tol) j = -1
           if (j /= 0) then
              pout(1) = pout(1) + j*x(1)
              pout(2) = pout(2) + j*x(2)
              pout(3) = pout(3) + j*x(3)
              goto 15
           endif
        enddo
     enddo
  enddo
end subroutine shorbz
!      subroutine fmain
!      implicit none
!      double precision plat(9),qlat(9),p(3),p1(3),xx
!      integer mode(3)

!      data plat /-0.5d0,0.5d0,1.7517d0,
!     .            0.5d0,-.5d0,1.7517d0,
!     .            0.5d0,0.5d0,-1.7517d0/

!      integer w(10000)
!      common /w/ w

!      data p /0.0d0,-0.5d0,-1.2638400000000001d0/

!      call wkinit(10000)

!C ... qlat = (plat^-1)^T so that qlat^T . plat = 1
!      call mkqlat(plat,qlat,xx)

!      call shorbz(p,p1,plat,qlat)
!      call prmx('p1 from shorbz',p1,1,1,3)

!      mode(1) = 2
!      mode(2) = 2
!      mode(3) = 3
!      call shorps(1,plat,mode,p,p1)
!      call prmx('p1 from shorps',p1,1,1,3)
!      end


subroutine shorbzm(flags,p,pout,plat,qlat)
  !- Shortens vector to equivalent in first Brillouin zone.
  ! ----------------------------------------------------------------
  !i Inputs:
  !i   plat,qlat lattice vectors and inverse
  !i   p         vector to shorten
  !o Outputs:
  !o   pout      shortened p.   pout and p can point to the same address.
  !r Remarks
  !r   Switch around plat,qlat to shorten reciprocal space vectors.
  !r   Jan 1997 Adapted from shorps to fix bug:  Example:
  !r   plat=  -0.5  0.5  1.7517  0.5  -0.5  1.7517  0.5  0.5  -1.7517
  !r   p= 0.0 -0.5 -1.26384
  !r   Should get pout -> 0.5 0.0 0.48786, not -0.5 1.0 0.48786.
  ! ----------------------------------------------------------------
  !     implicit none
  integer,intent(in)::flags(3)
  double precision :: p(3),pout(3),plat(3,3),qlat(3,3),x(3),x0,xx,a2,ap
  double precision :: tol
  parameter (tol=-1d-10)
  integer :: i,j,m,j2min,j3min,j1,j2,j3
  real(8),parameter:: eps=1.0d-10

  ! --- Reduce to unit cell centered at origin ---
  do  1  i = 1, 3
     ! ... x is projection of pin along plat(i), with multiples of p removed
     x0 = p(1)*qlat(1,i)+p(2)*qlat(2,i)+p(3)*qlat(3,i)
     if (flags(i) /= 0) then
        xx = idnint(x0-eps)
     else
        xx = idnint(x0)
     endif
     x(i) = x0-xx
1 enddo
  ! ... pout is x rotated back to Cartesian coordinates
  do  2  m = 1, 3
     pout(m) = x(1)*plat(m,1)+x(2)*plat(m,2)+x(3)*plat(m,3)
2 enddo

  ! --- Try shortening by adding +/- basis vectors ---
15 continue
  do   j1 =  0, 1
     j2min = -1
     if (j1 == 0) j2min = 0
     do   j2 = j2min, 1
        j3min = -1
        if (j1 == 0 .AND. j2 == 0) j3min = 0
        do   j3 = j3min, 1

           !     ... (-1,0,1) (plat(1) + (-1,0,1) plat(2)) + (-1,0,1) plat(3))
           ! takao think this alglorism cause problems (=wrong) for very anisotropic cases.

           do  17  i = 1, 3
              x(i) = plat(i,1)*j1 + plat(i,2)*j2 + plat(i,3)*j3
17         enddo
           a2 = x(1)*x(1) + x(2)*x(2) + x(3)*x(3)
           ap = pout(1)*x(1) + pout(2)*x(2) + pout(3)*x(3)
           j = 0
           if (a2 + 2*ap < tol) j = 1
           if (a2 - 2*ap < tol) j = -1
           if (j /= 0) then
              pout(1) = pout(1) + j*x(1)
              pout(2) = pout(2) + j*x(2)
              pout(3) = pout(3) + j*x(3)
              goto 15
           endif
        enddo
     enddo
  enddo
end subroutine shorbzm

