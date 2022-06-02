subroutine hklos(rsm,e,kmax,h0k)
  !- k,L-dependent smooth hankel functions at (0,0,0)
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   rsm   :smoothing radius
  !i   e     :energy of smoothed Hankel
  !i   kmax  :polynomial cutoff
  !o Outputs
  !i   h0k   :Smoothed Hankel at origin
  !r Remarks
  !r   Only the functions for l=0 are generated (remaining are zero)
  !r   Equivalent to hkl_ml for p=(0,0,0) and lmax=0.
  !u Updates
  !u   24 Apr 00 Adapted from nfp hkl_os.f
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: kmax
  double precision :: rsm,e,h0k(0:kmax)
  ! ... Local parameters
  integer :: k
  double precision :: akap,asm,cc,derfc,gam,gg,hh,pi,y0

  pi = 4d0*datan(1d0)
  y0 = 1d0/dsqrt(4*pi)
  gam = rsm*rsm/4d0
  asm = 0.5d0/dsqrt(gam)

  ! ... Make smooth Hankel at zero
  if (e > 0d0) call rx('hklos: e is positive')
  akap = dsqrt(dabs(e))
  cc = 4d0*asm**3 * dexp(gam*e) / dsqrt(pi)
  hh = cc/(2*asm*asm)-akap*derfc(akap/(2*asm))
  h0k(0) = hh*y0

  ! ... Upward recursion for laplace**k h0
  gg = dexp(gam*e) * (asm*asm/pi)**1.5d0
  do k = 1,kmax
     hh = -e*hh - 4*pi*gg
     h0k(k) = hh*y0
     gg = -2*asm*asm*(2*k+1)*gg
  enddo

end subroutine hklos

