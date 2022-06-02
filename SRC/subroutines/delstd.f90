subroutine delstd(n,x,d,s,e,ep)
  !- Returns generalised delta and step functions (Methfessel & Paxton)
  !-----------------------------------------------------------------------
  !i  Inputs
  !i    n  : order of approximant; see Remarks
  !i       : n>=0 returns Methfessel-Paxton broadening
  !i       : n<0  returns Fermi-Dirac broadening
  !i    x  : (efermi - e) / width
  !i       : width should be gaussian width (n>=0)
  !i       : or temperature (n<0)
  !o  Outputs
  !o    D_n (x): smeared delta-function
  !o    S_n (x): smeared heaviside function
  !o    e_n (x): entropy
  !o    ep     : de/dx (Fermi-Dirac case only)
  !r  Remarks
  !r    For Methfessel-Paxton (generalized gaussian) broadening
  !r    (see Phys Rev B40, 3616 (1989))
  !r      D_n (x) = exp -x^2 * sum_i=0^n A_i H_2i(x)
  !r      S_n (x) = (1 - erf x)/2 + exp -x^2 * sum_i=1^n A_i H_{2i-1}(x)
  !r      where H is a Hermite polynomial and
  !r      A_i = (-1)^i / ( i! 4^i sqrt(pi) )
  !r    For Fermi-Dirac broadening
  !r      s = 1/(exp(x)+1)   (fermi function)
  !r      d = ds/dx = exp(x)*s*s
  !r      e = -( s*log(s) + (1-s)*log(1-s) )
  !r     ep = log(s) - log(1-s)
  !u Updates
  !u   04 Aug 07 extended delstp to returns ep in Fermi-dirac case
  !u   23 May 00 extended to handle Fermi-Dirac broadening
  !-----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: n
  double precision :: x,d,s,e,ep
  ! ... Local parameters
  integer :: i,k
  double precision :: a,h1,h2,h3,s0,ex2,derfc,srpi
  !      intrinsic dsqrt,datan,dexp
  external derfc

  srpi = dsqrt(4d0*datan(1d0))

  ! ... Fermi-Dirac broadening
  if (n < 0) then
     if (x < -36d0) goto 91
     if (x >  36d0) goto 92
     s = 1d0/(dexp(x)+1d0)
     d = dexp(x)*s*s
     e = -( s*dlog(s) + (1-s)*dlog(1-s) )
     ep = (dlog(s) - dlog(1-s)) * s**2 * exp(x)
     return
  endif

  ! ... Methfessel-Paxton broadening
  if (x < -6d0) goto 91
  if (x >  6d0) goto 92
  ex2 = dexp(-x*x)
  s0 = .5d0 * derfc(x)
  a = 1d0/srpi
  k = 0
  h1 = 1d0
  h2 = 2d0 * x
  s = 0d0
  d = a
  do  1  i = 1, n
     a = -a / ( 4d0*i )
     k = k+1
     h3 = h1
     h1 = h2
     h2 = 2*x*h2 - 2*k*h3
     s = s + a*h1
     k = k+1
     h3 = h1
     h1 = h2
     h2 = 2*x*h2 - 2*k*h3
     d = d + a*h1
1 enddo
  d = d * ex2
  s = s0 + s*ex2
  e = 0.5d0*a*h1*ex2
  ep = 0
  return

  ! ... Branch for very small or very large x
91 s = 1d0
  e = 0d0
  d = 0d0
  ep = 0d0
  return

92 s = 0d0
  e = 0d0
  d = 0d0
  ep = 0d0
  return

end subroutine delstd

