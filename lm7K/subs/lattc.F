#define F90 1
      subroutine lattc(as,tol,rpad,alat,alat0,platin,g1,g2,g3,gt,plat,
     .  qlat,lmax,vol,awald,dlat,nkd,glat,nkg,nkdmx,nkgmx)
C- Sets up the real and reciprocal space lattice vectors for Ewald
C ----------------------------------------------------------------
Ci Inputs
Ci   as    :dimensionless Ewald parameter (2 is suggested).
Ci         :Ewald parameter awald scales with the lattice as
Ci         :as/(vol)**(1/3)
Ci   tol   :tolerance for ewald sums
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   platin:primitive lattice translation vectors, in units of alat
Ci   g1,g2,g3:x,y,z distortions of platin
Ci   gt:    multiplier of g1,g2,g3.  gt=1 => no distortion
Ci   lmax:  Ewald sums will be taken to maximum L.
Ci   nkdmx: maximum number of direct lattice vectors
Ci   nkgmx: maximum number of reciprocal lattice vectors
Co Outputs
Co    plat:   distorted lattice vectors
Co    qlat:   distorted reciprocal vectors
Co     vol:   cell volume
Co   awald:   ewald parameter
Co   dlat,nkd: direct lattice vectors and number
Co   glat,nkg: reciprocal lattice vectors and number
Cr Remarks
Cr   awald is in (atomic units)^-1.
Cr   Tolerance tol is the estimated error for a lattice of unit volume.
Cr   plat(*,k) holds lattice vector k
Cr   plat can point to same address as platin; then platin => plat
Cr
Cr   Local variables r0 and q0 are the ranges for
Cr   real-space and reciprocal lattice vectors for a unit lattice.
Cu Updates
Cu   2 Mar 04 New rpad: truncate radius of lattice vectors to rpad*rmax
Cu            when list has to be padded in order to include at
Cr            least one lattice vector.
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer lmax,nkd,nkg,nkdmx,nkgmx
      double precision as,tol,alat,g1,g2,g3,gt,vol,awald,alat0,rpad,
     . glat(3,nkgmx),dlat(3,nkdmx),platin(3,3),plat(3,3),qlat(3,3)
C Local parameters
      integer k,iprint,m,i1mach,modeg(3),isw,stdo,lgunit
      double precision qlat0(3,3),vol0,plat0(3,3),radd,qadd
      double precision qdist0,a0,rdist0,tol1,r0,q0,one(3,3),oned(3,3)

      stdo = lgunit(1)
      call dcopy(9,platin,1,plat0,1)
      call dinv33(plat0,1,qlat0,vol0)
      vol0 = dabs(vol0)
      call rdistn(plat0,plat,3,g1,g2,g3,gt)
C     qlat = (plat+)^-1
      call dinv33(plat,1,qlat,vol)
      vol = dabs(vol)*(alat**3)
      if (iprint() .lt. 20) goto 20
      write(stdo,351)
  351 format(/t17,'Plat',t55,'Qlat')
      write (stdo,350) ((plat0(m,k),m=1,3),(qlat0(m,k),m=1,3),k=1,3)
  350 format(3f11.6,5x,3f11.6)
      if (dabs(gt-1d0) .gt. 1.d-5) then
        call awrit4(' Distorted with gx,y,z= %d %d %d  gt=%d:',
     .    ' ',80,i1mach(2),g1,g2,g3,gt)
        write(stdo,350) ((plat(m,k),m=1,3),(qlat(m,k),m=1,3),k=1,3)
        call dpzero(one,9)
        one(1,1) = 1
        one(2,2) = 1
        one(3,3) = 1
        call rdistn(one,oned,3,g1,g2,g3,gt)
        if (iprint() .gt. 40) write(stdo,352) ((oned(m,k),m=1,3), k=1,3)
  352   format(t14,'shear matrix'/(3f11.7))
      endif
      call awrit3('  Cell vol= %1,6;6d'//
     .  '%?#n# (undistorted vol=%1,6;6d)#%j#',' ',80,i1mach(2),
     .  vol,isw(dabs(vol-vol0*(alat**3)).gt.1d-9),vol0*(alat**3))
   20 continue

C --- Set up real and reciprocal vectors ---
C The errors are estimated making a continuum approximation to a
C discrete set of lattice sums.  Adding .7, slightly more than
C half the average spacing makes the continuum approximation err
C on the safe side.
      rdist0 = vol0**(1d0/3d0)
      qdist0 = 1d0/rdist0
      radd = .7d0*rdist0
      qadd = .7d0*qdist0
      a0 = as/rdist0
      awald = a0/alat
      tol1 = tol*alat0**(lmax+1)
      call lctoff(a0,vol0,lmax,tol1,r0,q0)
      modeg(1) = 2
      modeg(2) = 2
      modeg(3) = 2
      call xlgen(plat0,r0+radd,rpad*(r0+radd),nkdmx,11,modeg,nkd,dlat)
      call xlgen(qlat0,q0+qadd,rpad*(q0+qadd),nkgmx,11,modeg,nkg,glat)
      call rdistn(dlat,dlat,nkd,g1,g2,g3,gt)
      call qdistn(glat,glat,nkg,g1,g2,g3,gt)

C --- Printout ---
      if (iprint() .lt. 30) goto 60
      write (stdo,340) as,tol,alat,awald
      write (stdo,342) r0+radd,nkd,q0+qadd,nkg
  340 format(/' LATTC:  as=',f6.3,'   tol=',1p,e8.2,
     .        '   alat=',0p,f8.5,'   awald=',f6.3)
  342 format(9x,'r1=',f7.3,'   nkd=',i4,'      q1=',f7.3,'   nkg=',i4)
   60 if(dabs(alat0/alat-1d0) .gt. 0.04d0)
     .  call rx('lattc: alat and alat0 deviate by more than 4 %')

      end
      subroutine lctoff(a0,v0,lmax,tol,r0,q0)
C- makes limits r0,q0 for sums in real and recip space for a lattice
c  with lattice constant 1.
Cu Updates
Cu   25 Jun 03 (Kino) bug fix in dimension of f and g
C     implicit none
      integer lmax,i
      double precision a0,q0,r0,tol,v0
      double precision gq0,gq1,pi,q1,q2,r1,r2
#if F90 | AUTO_ARRAY
      double precision f(0:lmax),g(0:lmax)
#else
      double precision f(0:40),g(0:40)
#endif

      pi = 4d0*datan(1d0)
      q1 = 0.001d0
      if (lmax .gt. 2) q1 = dsqrt(.5d0*(lmax-2))*a0/pi
      gq1 = (2d0*pi*q1)**(lmax-2)*dexp(-(pi*q1/a0)**2)*4d0*pi/v0
      if (tol.gt.gq1) call info0(10,0,0,' lctoff (warning): tol gt gq1')
      q2 = 50d0
      q0 = 5d0
      do  33  i = 1, 25
        gq0 = (2d0*pi*q0)**(lmax-2)*dexp(-(pi*q0/a0)**2)*4d0*pi/v0
        if(gq0 .gt. tol) q1 = q0
        if(gq0 .lt. tol) q2 = q0
        q0 = .5d0*(q1+q2)
   33 continue
      r1 = 0.1d0
      r2 = 50d0
      r0 = 5d0
      do  15  i = 1, 25
        call dlmtor(r0,a0,lmax,f,g)
        if(f(lmax) .gt. tol) r1 = r0
        if(f(lmax) .le. tol) r2 = r0
        r0 = .5d0*(r1+r2)
   15 continue
c|    try = (2d0*pi*q0)**(lmax-2)*dexp(-(pi*q0/a0)**2)*4d0*pi/v0
c|    write(6,957) q0,try,r0,f(lmax)
c|957 format(' lcut: q0=',f12.6,'   try=',f12.6,'   r0=',f12.6,
c|   .  '   f=',f12.6)
      end
      subroutine dlmtor(r,a,lmax,f,fbar)
C- Radial part of damped lmtos f and fbar, l=0 to lmax
C     implicit none
      integer l,lmax
      double precision a,f(0:lmax),fbar(0:lmax),r
      double precision derfc,emz2,erfc0,erfc1,erfc2,fbsrpi,
     .  flm2,g,ta2r,z

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
   10 continue
      end

