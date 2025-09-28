module m_lattic ! Sets up the real and rmeciprocal space lattice vectors !no shear now  ldist=0
  use m_ftox
  use m_lgunit,only:stdo
  use m_xlgen,only:xlgen
  public m_lattic_init,lctoff ,setopos,readatompos
  real(8), allocatable,protected,public ::  rv_a_odlv (:,:)
  real(8), allocatable,protected,public ::  rv_a_oqlv (:,:)
  real(8), protected,public:: lat_plat(3,3),lat_qlat(3,3),lat_awald,lat_vol
  integer,protected,public:: lat_nkd,lat_nkq
  logical,protected,public:: lattic_init=.false.
  real(8), allocatable,protected,public :: rv_a_opos(:,:)
  private
  real(8),parameter:: pi = 4d0*datan(1d0)
contains
  
  subroutine readatompos(irpos) !readin atomic position from AtomPos file if available
    use m_ext,only:sname
    use m_lmfinit,only:nbas
    implicit none
    intent(out)::        irpos
    real(8):: pos(3,nbas)
    logical:: irpos
    block
      integer::i,nbaso,ifipos
      real(8):: p(3)
      irpos=.false.
      open(newunit=ifipos,file='AtomPos.'//trim(sname),status='old',err=1010)
      do 
         read(ifipos,*,end=1010)
         read(ifipos,*)
         read(ifipos,*) nbaso
         do i=1,nbaso
            read(ifipos,*) p
            if(i<=nbas) pos(:,i)=p    !write(stdo,ftox)i,ftof(p)
         enddo
         irpos=.true.
      enddo
1010  continue
      close(ifipos)
    endblock
    if(irpos) rv_a_opos=pos
  endsubroutine readatompos
  
  subroutine setopos(posin) !called from lmfp to revise atomic position
    real(8):: posin(:,:)
    rv_a_opos= posin
  end subroutine Setopos
  
  subroutine m_lattic_init() 
    use m_lmfinit,only:nbas,alat=>lat_alat,as=>lat_as,tol=>lat_tol,rpad=>lat_rpad,nkdmx=>lat_nkdmx,nkqmx=>lat_nkqmx,lat_platin,pos
    implicit none
    integer::  lmxst , nkd,nkq,ib
    real(8) :: dlv(3,nkdmx),qlv(3,nkqmx), awald,vol,xx1,xx2,dotprd, plat0(3,3),plat(3,3),qlat(3,3) 
    !i Inputs
    !i   as    :Ewald smoothing parameter. dimensionless Ewald parameter (2 is suggested).
    !i         :Ewald parameter awald scales with the lattice as as/(vol)**(1/3)
    !i   tol   :tolerance for ewald sums
    !i   alat  :length scale of lattice and basis vectors, a.u.
    !i   platin:primitive lattice translation vectors, in units of alat
    !i   lmax:  Ewald sums will be taken to maximum L.
    !i   nkdmx: maximum number of direct lattice vectors
    !i   nkqmx: maximum number of reciprocal lattice vectors
    !o Outputs
    !o    lat_plat:   distorted lattice vectors
    !o    lat_qlat:   distorted reciprocal vectors
    !o    lat_vol:   cell volume
    !o    lat_awald:   ewald parameter
    !o    dlv,nkd: direct lattice vectors and number for Ewald sum
    !o    qlv,nkq: reciprocal lattice vectors and number for Ewald sum
    !r Remarks
    !r   awald is in (atomic units)^-1.
    !r   Tolerance tol is the estimated error for a lattice of unit volume.
    !r   plat(*,k) holds lattice vector k
    !r   plat can point to same address as platin; then platin => plat
    !r
    !r   Local variables r0 and q0 are the ranges for
    !r   real-space and reciprocal lattice vectors for a unit lattice.
    !r  rpad..  truncate Ewald to rpad*rmax when lattice vector
    !r          list has to be padded in order to include at
    !r          least one lattice vector
    !r  tol     Ewald tolerance
    !u Updates
    !u   2 Mar 04 New rpad: truncate radius of lattice vectors to rpad*rmax when list has to be padded in order to include at least one lattice vector.
    call tcn('m_lattic_init')
    lattic_init=.true.
    plat = lat_platin
    lattc: block !subroutine lattc(as,tol,rpad,alat,plat,qlat,lmax,vol,awald,dlv,nkd,qlv,nkq,nkdmx,nkqmx)
      ! Sets up the real and reciprocal space lattice vectors for Ewald
      integer :: lmax=6 , k,iprint,m,modeg(3),ifp,i
      real(8) :: qlat0(3,3),vol0,plat0(3,3),radd,qadd,qdist0,a0,rdist0,tol1,r0,q0,one(3,3),oned(3,3),rxx,qxx
      real(8),external:: bohr
      plat0=plat 
      call dinv33(plat0,1,qlat0,vol0)
      vol0 = dabs(vol0)
      plat=plat0 
      call dinv33(plat,1,qlat,vol)
      print *,'aaaaaaa alat=',alat
      vol = dabs(vol)*(alat**3)
      if(iprint()>0) then
         write(stdo,"(/t17,'Plat',t55,'Qlat')")
         write(stdo,350) ((plat0(m,k),m=1,3),(qlat0(m,k),m=1,3),k=1,3)
         open(newunit=ifp,file='PlatQlat.chk')
         write(ifp,350) ((plat0(m,k),m=1,3),(qlat0(m,k),m=1,3),k=1,3)
         write(ifp,"('             PLAT              and         QLAT    ')")
         write(ifp,ftox) ftof(alat),ftof(alat*bohr()), ' ! Unit:P alat/bohr, alat/AA'
         write(ifp,ftox) ftof(2*pi/alat),ftof(2*pi/(alat*bohr())),' ! Unit:Q 2pi/alat /(bohr^-1), 2pi/alat/(AA^-1)'
         write(ifp,ftox) ftof(vol),ftof(vol*bohr()**3),' ! Cell vol (bohr**3) (AA**3)'
         close(ifp)
350      format(3f11.6,5x,3f11.6)
         write(stdo,ftox)'  Cell vol (borh**3)= ',ftof(vol)
         if((dabs(vol-vol0*(alat**3)) > 1d-9))write(stdo,ftox)'(undistorted vol=',ftof(vol0*(alat**3))
      endif
      ! --- Set up real and reciprocal vectors ---
      ! The errors are estimated making a continuum approximation to a
      ! discrete set of lattice sums.  Adding .7, slightly more than half the average spacing makes the continuum approximation err on the safe side.
      rdist0 = vol0**(1d0/3d0)
      qdist0 = 1d0/rdist0
      radd = .7d0*rdist0  
      qadd = .7d0*qdist0  != 1.2d0*qdist0 
      a0 = as/rdist0
      awald = a0/alat
      tol1 = tol*alat**(lmax+1)
      call lctoff(a0,vol0,lmax,tol1,r0,q0)
      modeg = 2
      rxx= maxval([r0+radd,(sum(plat0(:,i)**2)**.5*1.05,i=1,3)]) !2022-10-12 for very anisotropic cases safer.
      qxx= maxval([q0+qadd,(sum(qlat0(:,i)**2)**.5*1.05,i=1,3)]) !2022-10-12
      call xlgen(plat0,rxx,rpad*(r0+radd),nkdmx,11,modeg,nkd,dlv)
      call xlgen(qlat0,qxx,rpad*(q0+qadd),nkqmx,11,modeg,nkq,qlv)
      if(iprint()>0)write(stdo,"(/'m_lattic_init:  as=',f6.3,'  tol=',1p,e9.2,'  alat=',0p,f8.5,'   awald=',f6.3)")as,tol,alat,awald
      if(iprint()>0)write(stdo,"(9x,'r1=',f7.3,'   nkd=',i4,'      q1=',f7.3,'   nkq=',i4)") r0+radd,nkd,q0+qadd,nkq
    endblock lattc
    lat_vol  =vol
    lat_plat =plat
    lat_qlat =qlat
    lat_awald=awald
    lat_nkd=nkd
    lat_nkq=nkq
    allocate(rv_a_opos,source=pos)
    allocate(rv_a_oqlv, source=qlv(1:3,1:nkq))
    allocate(rv_a_odlv, source=dlv(1:3,1:nkd))
    call tcx('m_lattic_init')
  end subroutine m_lattic_init
  subroutine lctoff(a0,v0,lmax,tol,r0,q0) !- makes limits r0,q0 for sums in real and recip space for a lattice
    !  with lattice constant 1.    !u   25 Jun 03 (Kino) bug fix in dimension of f and g
    implicit none
    integer :: lmax,i
    real(8) :: a0,q0,r0,tol,v0, gq0,gq1,q1,q2,r1,r2, f(0:lmax),g(0:lmax)
!    pi = 4d0*datan(1d0)
    q1 = 0.001d0
    if (lmax > 2) q1 = dsqrt(.5d0*(lmax-2))*a0/pi
    gq1 = (2d0*pi*q1)**(lmax-2)*dexp(-(pi*q1/a0)**2)*4d0*pi/v0
    if (tol > gq1) write(stdo,*)' lctoff (warning): tol gt gq1'
    q2 = 50d0
    q0 = 5d0
    do  33  i = 1,25
       gq0 = (2d0*pi*q0)**(lmax-2)*dexp(-(pi*q0/a0)**2)*4d0*pi/v0 !!this can cause
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
  end subroutine lctoff
  subroutine dlmtor(r,a,lmax,f,fbar)  !- Radial part of damped lmtos f and fbar, l=0 to lmax
    implicit none
    integer :: l,lmax
    real(8):: a,f(0:lmax),fbar(0:lmax),r, derfc,emz2,erfc0,erfc1,erfc2,fbsrpi, flm2,g,ta2r,z
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
end module m_lattic
