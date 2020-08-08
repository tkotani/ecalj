      subroutine hklos(rsm,e,kmax,h0k)
C- k,L-dependent smooth hankel functions at (0,0,0)
C ----------------------------------------------------------------------
Ci Inputs
Ci   rsm   :smoothing radius
Ci   e     :energy of smoothed Hankel
Ci   kmax  :polynomial cutoff
Co Outputs
Ci   h0k   :Smoothed Hankel at origin
Cr Remarks
Cr   Only the functions for l=0 are generated (remaining are zero)
Cr   Equivalent to hkl_ml for p=(0,0,0) and lmax=0.
Cu Updates
Cu   24 Apr 00 Adapted from nfp hkl_os.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer kmax
      double precision rsm,e,h0k(0:kmax)
C ... Local parameters
      integer k
      double precision akap,asm,cc,derfc,gam,gg,hh,pi,y0

      pi = 4d0*datan(1d0)
      y0 = 1d0/dsqrt(4*pi)
      gam = rsm*rsm/4d0
      asm = 0.5d0/dsqrt(gam)

C ... Make smooth Hankel at zero
      if (e .gt. 0d0) call rx('hklos: e is positive')
      akap = dsqrt(dabs(e))
      cc = 4d0*asm**3 * dexp(gam*e) / dsqrt(pi)
      hh = cc/(2*asm*asm)-akap*derfc(akap/(2*asm))
      h0k(0) = hh*y0

C ... Upward recursion for laplace**k h0
      gg = dexp(gam*e) * (asm*asm/pi)**1.5d0
      do k = 1,kmax
        hh = -e*hh - 4*pi*gg
        h0k(k) = hh*y0
        gg = -2*asm*asm*(2*k+1)*gg
      enddo

      end

