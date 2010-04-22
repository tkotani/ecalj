      subroutine soldhj(r,e,loka,lmax,hl,bl,cy)
C- Real solid hankel and bessel functions.
C ----------------------------------------------------------------
Ci Inputs
Ci   r     :radius (or radius/avw using OKA's conventions)
Ci   e     :energy (or energy*avw**2 using OKA's conventions)
Ci   loka  :conventions for Hankels, Bessels; see besslr
Ci   lmax  :maximum l for a given site
Ci   cy    :Normalization constants for spherical harmonics
Co Outputs
Co   HL,BL: Hankel and Bessel functions:  calculated to (lmax+1)**2
Cr Remarks
Cr   Generates bessel function * real spherical harmonic
Cr   MSM's standard defs, notes IV-43.
Cu Updates
Cu   19 May 04 Changed loka from logical to integer
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer lmax,loka
      double precision e
      double precision bl(1),cy(1),hl(1),r(3)
C Local parameters
      integer ilm,l,m,ll
      double precision rfac,xx,r2
      double precision phi(0:10),psi(0:10)
      external ll

      call sylm(r,hl,lmax,r2)
      call besslr(e*r2,loka,0,lmax,phi,psi)
      ilm = 0
      rfac = dsqrt(r2)
      if (r2 .lt. 1.d-10) r2 = 1.d0
      do  10  l = 0, lmax
C       rfac = 1/r**(2l+1), or 1/(r/w)**(2l+1) using OKA conventions
        rfac = rfac/r2
        do  10  m = -l, l
        ilm = ilm+1
        xx = cy(ilm)*hl(ilm)
        bl(ilm) = phi(l)*xx
        hl(ilm) = (rfac*psi(l))*xx
  10  continue
      end

