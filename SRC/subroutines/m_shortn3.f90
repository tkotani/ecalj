module m_shortvec ! core of m_shortn3: Find shortest vectors in modulo of rlat 
  public:: shortvec,shortvecinitialize
contains
  subroutine shortvec(pin,rlatp,xmx2,noutmx,nout,nlatout)
    !i pin is the fractional coodinate on rlat.
    !i rlatp,xmx2 are passed from shortvecinitialize
    !   rlatp(i,j)= sum( rlat(:,i)*rlat(:,j) )
    !   rlat(3,i) i-th vertor for modulo
    !i noutmax: upper limit of nlatout
    !o nout
    !o nlatout
    ! Shortest vectors are
    !    pin+nlatout(:,ix), where ix=1:nout, is the shortest vectors.
    !    We may have multiple nlatout (# is nout).
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
    nknknk= (2*nmax(1)+1)*(2*nmax(2)+1)*(2*nmax(3)+1)
    allocate( nlat0(3, nknknk), rnorm(nknknk) )
    do ik1=-nmax(1),nmax(1)
       do ik2=-nmax(2),nmax(2)
          do ik3=-nmax(3),nmax(3)
             ik=ik+1
             nlat0(:,ik) = [ik1,ik2,ik3]
             rr= pin + nlat0(:,ik)
             rnorm(ik) = sum(rr*matmul(rlatp,rr))
          enddo
       enddo
    enddo
    nk=ik
    rmin=minval(rnorm(1:nk))
    nout=0
    do ik=1,nk
       rr= pin + nlat0(:,ik)
       if(rnorm(ik)<rmin+eps) then
          nout=nout+1
          if(nout>noutmx) call rx('shortn3: enlarge noutmx')
          nlatout(:,nout)=nlat0(:,ik)
       endif
    enddo
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
module m_shortn3! shortest vectors. modulo of rlat(1:3,i),i=1,3)
  ! Set rlat at first, then call shortn3(pin)
  ! Shortest p= pin + matmul(rlat(:,:),nlatout(:,i)) for i=1,nout
  ! Only when pin is on the Volonoi boundaries, nout>1.
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
  implicit none
  !!== Maximum value for x_i for ellipsoid ==
  !!  Ellipsoid is given as 1d0 = sum x_i nn(i,j) x_j.
  !i nn(3,3)
  !o xmx2(i)  Maximum of x_i**2
  real(8):: nn(3,3),v2(2),ainv(2,2), rmax2, xmx2(3),det,fac1,fac2,fac3,nv2(2)
  associate(&
       n11=>nn(1,1),n12=>nn(1,2),n13=>nn(1,3), &
       n21=>nn(2,1),n22=>nn(2,2),n23=>nn(2,3), &
       n31=>nn(3,1),n32=>nn(3,2),n33=>nn(3,3))
    ainv(1,1:2)=  [ n33,-n23]
    ainv(2,1:2)=  [-n32, n22]
    fac1 = n11-sum([n12,n13]*matmul(ainv,[n12,n13]))/(n22*n33-n23*n32)
    ainv(1,1:2)=  [ n11,-n31]
    ainv(2,1:2)=  [-n13, n33]
    fac2 = n22-sum([n23,n21]*matmul(ainv,[n23,n21]))/(n33*n11-n31*n13)
    ainv(1,1:2)=  [ n22,-n12]
    ainv(2,1:2)=  [-n21, n11]
    fac3 = n33-sum([n31,n32] *matmul(ainv,[n31,n32]))/(n11*n22-n12*n21)
    xmx2 = [1d0/fac1,1d0/fac2,1d0/fac3]
  endassociate
end subroutine ellipsoidxmax
