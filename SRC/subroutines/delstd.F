      subroutine delstd(n,x,d,s,e,ep)
C- Returns generalised delta and step functions (Methfessel & Paxton)
C-----------------------------------------------------------------------
Ci  Inputs
Ci    n  : order of approximant; see Remarks
Ci       : n>=0 returns Methfessel-Paxton broadening
Ci       : n<0  returns Fermi-Dirac broadening
Ci    x  : (efermi - e) / width
Ci       : width should be gaussian width (n>=0)
Ci       : or temperature (n<0)
Co  Outputs
Co    D_n (x): smeared delta-function
Co    S_n (x): smeared heaviside function
Co    e_n (x): entropy
Co    ep     : de/dx (Fermi-Dirac case only)
Cr  Remarks
Cr    For Methfessel-Paxton (generalized gaussian) broadening
Cr    (see Phys Rev B40, 3616 (1989))
Cr      D_n (x) = exp -x^2 * sum_i=0^n A_i H_2i(x)
Cr      S_n (x) = (1 - erf x)/2 + exp -x^2 * sum_i=1^n A_i H_{2i-1}(x)
Cr      where H is a Hermite polynomial and
Cr      A_i = (-1)^i / ( i! 4^i sqrt(pi) )
Cr    For Fermi-Dirac broadening
Cr      s = 1/(exp(x)+1)   (fermi function)
Cr      d = ds/dx = exp(x)*s*s
Cr      e = -( s*log(s) + (1-s)*log(1-s) )
Cr     ep = log(s) - log(1-s)
Cu Updates
Cu   04 Aug 07 extended delstp to returns ep in Fermi-dirac case
Cu   23 May 00 extended to handle Fermi-Dirac broadening
C-----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer n
      double precision x,d,s,e,ep
C ... Local parameters
      integer i,k
      double precision a,h1,h2,h3,s0,ex2,derfc,srpi
c      intrinsic dsqrt,datan,dexp
      external derfc

      srpi = dsqrt(4d0*datan(1d0))

C ... Fermi-Dirac broadening
      if (n .lt. 0) then
        if (x .lt. -36d0) goto 91
        if (x .gt.  36d0) goto 92
        s = 1d0/(dexp(x)+1d0)
        d = dexp(x)*s*s
        e = -( s*dlog(s) + (1-s)*dlog(1-s) )
        ep = (dlog(s) - dlog(1-s)) * s**2 * exp(x)
        return
      endif

C ... Methfessel-Paxton broadening
      if (x .lt. -6d0) goto 91
      if (x .gt.  6d0) goto 92
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
    1 continue
      d = d * ex2
      s = s0 + s*ex2
      e = 0.5d0*a*h1*ex2
      ep = 0
      return

C ... Branch for very small or very large x
  91  s = 1d0
      e = 0d0
      d = 0d0
      ep = 0d0
      return

  92  s = 0d0
      e = 0d0
      d = 0d0
      ep = 0d0
      return

      end

