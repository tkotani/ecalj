module m_lattic
  use m_xlgen,only:xlgen
  public m_lattic_init,setopos,lctoff
  real(8), allocatable,protected,public ::  rv_a_odlv (:)
  real(8), allocatable,protected,public ::  rv_a_oqlv (:)
  real(8), protected,public:: lat_plat(3,3),lat_qlat(3,3),lat_awald,lat_vol
  real(8), allocatable,protected,public :: rv_a_opos(:,:)
  integer,protected,public:: lat_nkd,lat_nkq
  logical,protected,public:: lattic_init=.false.
  ! real(8),protected:: lat_dist(3,3) !unused because no deformation of cell allowed currently.
  private
contains
  subroutine Setopos(posin) !called from lmfp to revise atomic position by reading rst file.
    use m_lmfinit,only: nbas
    real(8):: posin(:,:)
    rv_a_opos(:,1:nbas)= posin(:,1:nbas)
  end subroutine Setopos
  subroutine m_lattic_init() ! Sets up the real and rmeciprocal space lattice vectors !no shear now  ldist=0
    use m_lmfinit,only:nbas,lat_alat,lat_as,lat_tol,lat_rpad,lat_nkdmx,lat_nkqmx,lat_gam,&
         lat_platin, poss=>pos
    implicit none
    integer::  lmxst , nkd , nkdmx , nkq , nkqmx,i_data_size,ib
    real(8),allocatable:: rv_a_tmp(:)
    real(8):: alat,awald,awald0,gam(4),gx,gy,gz,gt,tol,vol, &
         xx1,xx2,dotprd,pi,rpad, plat0(3,3),plat(3,3),qlat(3,3) 
    equivalence (gam(1), gx), (gam(2), gy), (gam(3), gz), (gam(4), gt)
    call tcn('m_lattic_init')
    lattic_init=.true.
    alat=lat_alat
    awald0=lat_as
    tol=lat_tol
    rpad=lat_rpad
    nkdmx=lat_nkdmx
    nkqmx=lat_nkqmx
    gam = lat_gam
    alat = lat_alat
    plat0=lat_platin
    allocate(rv_a_opos(3,nbas))
    rv_a_opos(:,1:nbas)= poss(:,1:nbas)
    if(abs(gt-1d0)>1d-10) call rdistn( rv_a_opos , rv_a_opos , nbas , gx , gy , gz , gt )
    allocate(rv_a_odlv(abs(3*nkdmx)))
    allocate(rv_a_oqlv(abs(3*nkqmx)))
    lmxst = 6
    call lattc ( awald0 , tol , rpad , alat , alat , plat0 , gx , &
         gy , gz , gt , plat , qlat , lmxst , vol , awald , rv_a_odlv &
         , nkd , rv_a_oqlv , nkq , nkdmx , nkqmx )
    lat_vol  =vol
    lat_plat =plat
    lat_qlat =qlat
    !! reduce size. necessary?
    i_data_size=size(rv_a_oqlv)
    allocate(rv_a_tmp(i_data_size))
    rv_a_tmp=rv_a_oqlv
    deallocate(rv_a_oqlv)
    i_data_size=min(i_data_size,3*nkq)
    allocate(rv_a_oqlv(3*nkq))
    rv_a_oqlv(:i_data_size)=rv_a_tmp(:i_data_size)
    deallocate(rv_a_tmp)
    !! reduce size. necessary?
    i_data_size=size(rv_a_odlv)
    allocate(rv_a_tmp(i_data_size))
    rv_a_tmp=rv_a_odlv
    deallocate(rv_a_odlv)
    i_data_size=min(i_data_size,3*nkd)
    allocate(rv_a_odlv(3*nkd))
    rv_a_odlv(:i_data_size)=rv_a_tmp(:i_data_size)
    deallocate(rv_a_tmp)
    lat_awald=awald
    lat_nkd=nkd
    lat_nkq=nkq
    call tcx('m_lattic_init')
  end subroutine m_lattic_init
  subroutine lattc(as,tol,rpad,alat,alat0,platin,g1,g2,g3,gt,plat,qlat,lmax,vol,awald,dlat,nkd,glat,nkg,nkdmx,nkgmx)! Sets
    !up the real and reciprocal space lattice vectors for Ewald
    use m_ftox
    use m_lgunit,only:stdo
    ! ----------------------------------------------------------------
    !i Inputs
    !i   as    :dimensionless Ewald parameter (2 is suggested).
    !i         :Ewald parameter awald scales with the lattice as
    !i         :as/(vol)**(1/3)
    !i   tol   :tolerance for ewald sums
    !i   alat  :length scale of lattice and basis vectors, a.u.
    !i   platin:primitive lattice translation vectors, in units of alat
    !i   g1,g2,g3:x,y,z distortions of platin
    !i   gt:    multiplier of g1,g2,g3.  gt=1 => no distortion
    !i   lmax:  Ewald sums will be taken to maximum L.
    !i   nkdmx: maximum number of direct lattice vectors
    !i   nkgmx: maximum number of reciprocal lattice vectors
    !o Outputs
    !o    plat:   distorted lattice vectors
    !o    qlat:   distorted reciprocal vectors
    !o     vol:   cell volume
    !o   awald:   ewald parameter
    !o   dlat,nkd: direct lattice vectors and number
    !o   glat,nkg: reciprocal lattice vectors and number
    !r Remarks
    !r   awald is in (atomic units)^-1.
    !r   Tolerance tol is the estimated error for a lattice of unit volume.
    !r   plat(*,k) holds lattice vector k
    !r   plat can point to same address as platin; then platin => plat
    !r
    !r   Local variables r0 and q0 are the ranges for
    !r   real-space and reciprocal lattice vectors for a unit lattice.
    !u Updates
    !u   2 Mar 04 New rpad: truncate radius of lattice vectors to rpad*rmax
    !u            when list has to be padded in order to include at
    !r            least one lattice vector.
    ! ----------------------------------------------------------------
    implicit none
    integer :: lmax,nkd,nkg,nkdmx,nkgmx
    double precision :: as,tol,alat,g1,g2,g3,gt,vol,awald,alat0,rpad, &
         glat(3,nkgmx),dlat(3,nkdmx),platin(3,3),plat(3,3),qlat(3,3)
    integer :: k,iprint,m,i1mach,modeg(3),isw
    double precision :: qlat0(3,3),vol0,plat0(3,3),radd,qadd
    double precision :: qdist0,a0,rdist0,tol1,r0,q0,one(3,3),oned(3,3)
    integer:: ifp,i
    real(8):: rxx,qxx
    plat0=platin !call dcopy(9,platin,1,plat0,1)
    call dinv33(plat0,1,qlat0,vol0)
    vol0 = dabs(vol0)
    call rdistn(plat0,plat,3,g1,g2,g3,gt)
    call dinv33(plat,1,qlat,vol)
    vol = dabs(vol)*(alat**3)
    if (iprint() < 20) goto 20
    write(stdo,351)
351 format(/t17,'Plat',t55,'Qlat')
    write (stdo,350) ((plat0(m,k),m=1,3),(qlat0(m,k),m=1,3),k=1,3)
    open(newunit=ifp,file='PlatQlat.chk')
    write (ifp,350) ((plat0(m,k),m=1,3),(qlat0(m,k),m=1,3),k=1,3)
    write (ifp,"('             PLAT              and         QLAT    ')")
    close(ifp)
350 format(3f11.6,5x,3f11.6)
    if (dabs(gt-1d0) > 1.d-5) then
       write(stdo,ftox)' Distorted with gx,y,z=',ftof(g1),ftof(g2),ftof(g3),'gt=:',ftof(gt)
       write(stdo,350) ((plat(m,k),m=1,3),(qlat(m,k),m=1,3),k=1,3)
       one=0
       one(1,1) = 1
       one(2,2) = 1
       one(3,3) = 1
       call rdistn(one,oned,3,g1,g2,g3,gt)
       if (iprint() > 40) write(stdo,352) ((oned(m,k),m=1,3), k=1,3)
352    format(t14,'shear matrix'/(3f11.7))
    endif
    write(stdo,ftox)'  Cell vol= ',ftof(vol)
    if((dabs(vol-vol0*(alat**3)) > 1d-9))write(stdo,ftox)'(undistorted vol=',ftof(vol0*(alat**3))
20  continue
    ! --- Set up real and reciprocal vectors ---
    ! The errors are estimated making a continuum approximation to a
    ! discrete set of lattice sums.  Adding .7, slightly more than
    ! half the average spacing makes the continuum approximation err
    ! on the safe side.
    rdist0 = vol0**(1d0/3d0)
    qdist0 = 1d0/rdist0
    radd = 0.7d0*rdist0 
    qadd = 0.7d0*qdist0 
    a0 = as/rdist0
    awald = a0/alat
    tol1 = tol*alat0**(lmax+1)
    call lctoff(a0,vol0,lmax,tol1,r0,q0)
    modeg(1) = 2
    modeg(2) = 2
    modeg(3) = 2
    !  rxx=r0+radd
    !  qxx=q0+qadd
    rxx= maxval([r0+radd,(sum(plat0(:,i)**2)**.5*1.05,i=1,3)]) !2022-10-12 for very anisotropic cases safer.
    qxx= maxval([q0+qadd,(sum(qlat0(:,i)**2)**.5*1.05,i=1,3)]) !2022-10-12
    call xlgen(plat0,rxx,rpad*(r0+radd),nkdmx,11,modeg,nkd,dlat)
    call xlgen(qlat0,qxx,rpad*(q0+qadd),nkgmx,11,modeg,nkg,glat)
    call rdistn(dlat,dlat,nkd,g1,g2,g3,gt)
    call qdistn(glat,glat,nkg,g1,g2,g3,gt)
    ! --- Printout ---
    if (iprint() < 30) goto 60
    write (stdo,340) as,tol,alat,awald
    write (stdo,342) r0+radd,nkd,q0+qadd,nkg
340 format(/'LATTC:  as=',f6.3,'   tol=',1p,e9.2, &
         '   alat=',0p,f8.5,'   awald=',f6.3)
342 format(9x,'r1=',f7.3,'   nkd=',i4,'      q1=',f7.3,'   nkg=',i4)
60  if(dabs(alat0/alat-1d0) > 0.04d0) call rx('lattc: alat and alat0 deviate by more than 4 %')
  end subroutine lattc
  subroutine lctoff(a0,v0,lmax,tol,r0,q0) !- makes limits r0,q0 for sums in real and recip space for a lattice
    use m_lgunit,only:stdo
    !  with lattice constant 1.
    !u Updates
    !u   25 Jun 03 (Kino) bug fix in dimension of f and g
    !     implicit none
    integer :: lmax,i
    double precision :: a0,q0,r0,tol,v0
    double precision :: gq0,gq1,pi,q1,q2,r1,r2
    double precision :: f(0:lmax),g(0:lmax)
    pi = 4d0*datan(1d0)
    q1 = 0.001d0
    if (lmax > 2) q1 = dsqrt(.5d0*(lmax-2))*a0/pi
    gq1 = (2d0*pi*q1)**(lmax-2)*dexp(-(pi*q1/a0)**2)*4d0*pi/v0
    if (tol > gq1) write(stdo,*)' lctoff (warning): tol gt gq1'
    q2 = 50d0
    q0 = 5d0
    do  33  i = 1, 25
       gq0 = (2d0*pi*q0)**(lmax-2)*dexp(-(pi*q0/a0)**2)*4d0*pi/v0
       if(gq0 > tol) q1 = q0
       if(gq0 < tol) q2 = q0
       q0 = .5d0*(q1+q2)
33  enddo
    r1 = 0.1d0
    r2 = 50d0
    r0 = 5d0
    do  15  i = 1, 25
       call dlmtor(r0,a0,lmax,f,g)
       if(f(lmax) > tol) r1 = r0
       if(f(lmax) <= tol) r2 = r0
       r0 = .5d0*(r1+r2)
15  enddo
    !|    try = (2d0*pi*q0)**(lmax-2)*dexp(-(pi*q0/a0)**2)*4d0*pi/v0
    !|    write(6,957) q0,try,r0,f(lmax)
    ! 957 format(' lcut: q0=',f12.6,'   try=',f12.6,'   r0=',f12.6,
    !|   .  '   f=',f12.6)
  end subroutine lctoff
  subroutine dlmtor(r,a,lmax,f,fbar)  !- Radial part of damped lmtos f and fbar, l=0 to lmax
    !     implicit none
    integer :: l,lmax
    double precision :: a,f(0:lmax),fbar(0:lmax),r
    double precision :: derfc,emz2,erfc0,erfc1,erfc2,fbsrpi, flm2,g,ta2r,z
    fbsrpi = 0.564189835d0
    z = a*r
    emz2 = dexp(-z*z)
    erfc0 = derfc(z)
    erfc1 = -z*erfc0 + fbsrpi*emz2
    erfc2 = -0.5d0*z*erfc1 + 0.25d0*erfc0
    f(0) = erfc0/r
    fbar(0) = -erfc2/(a*a*r)
    ta2r = 2d0*a*a*r
    g = 2d0*a*emz2*fbsrpi/r
    flm2 = fbsrpi*emz2/z - erfc0
    do  10  l = 1, lmax
       f(l) = ((l+l-1)/r)*f(l-1) + g
       fbar(l) = ((l+l-1)/r)*fbar(l-1) - flm2
       flm2 = f(l-1)
       g = g*ta2r
10  enddo
  end subroutine dlmtor
  subroutine qdist(q,q1,ux,uy,uz,gam) !  (from msm) new version of 03.10.89. u=(ux,uy,uz) gives direction,
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
  subroutine rdistn(a1,a2,na,gx,gy,gz,gt)  !- distort na real-space vectors in array a1 into array a2
    !     implicit none
    integer :: na,ia
    double precision :: a1(3,na),a2(3,na),gx,gy,gz,gt
    do 10 ia = 1, na
       call rdist(a1(1,ia),a2(1,ia),gx,gy,gz,gt)
10  enddo
  end subroutine rdistn
  subroutine qdistn(a1,a2,na,gx,gy,gz,gt)  !- distort na q-space vectors in array a1 into array a2
    !     implicit none
    integer :: na,ia
    double precision :: a1(3,na),a2(3,na),gx,gy,gz,gt
    do 10 ia = 1,na
       call qdist(a1(1,ia),a2(1,ia),gx,gy,gz,gt)
10  enddo
  end subroutine qdistn

end module m_lattic

