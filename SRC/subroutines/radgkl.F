      subroutine radgkl(r,rsm,kmax,lmax,n0,g)
C- Make radial parts g_kl of G_kL, divided by r**l.
C ----------------------------------------------------------------------
Ci Inputs
Ci   r     :radius
Ci   rsm   :smoothing radius
Ci   kmax  :make g_kl for k=0...kmax
Ci   lmax  :make g_kl for l=0...lmax
Ci   n0    :leading dimension of g
Co Outputs
Co   g     :radial part of generalized gaussians G_kL; see Remarks
Co         :radgkl produces g such that G_kL = g(k,l) r**l Y_L
Cr Remarks
Cr   Definition:  G_kL = (lap)**k YL(-grad) g0 = g_kl*r^l*Y_L
Cr   Energy factor beta = dexp(e/(4*a*a)) is not included here.
Cr   See J. Math. Phys. {\bf 39},3393 (1998), Section V:
Cr
Cr     G_kL  = phi_kl (a^2/pi)^(3/2) exp(-a^2 r^2) r^l Y_L
Cr           = phi_kl r^l Y_L g_00
Cr     G_0L  = (2a^2)^l r^l Y_L g_00
Cr
Cr   where (Eq. 12.3)
Cr
Cr    phi_kl = (-1)^k (2a^2)^(k+l) 2^k k! L_k^(l+1/2)(a^2 r^2)
Cr    L_k    = generalized Laguerre polynomial
Cr    phi_0l = (2a^2)^l
Cr
Cr   also
Cr
Cr                          a^l (2l+1)!!
Cr      p_kl = phi_kl ------------------------ r^l
Cr                    (2a^2)^(k+l) (2k+2l+1)!!
Cr
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer kmax,lmax,n0
      double precision r,rsm,g(0:n0,0:lmax)
C ... Local parameters
      integer l,k
      double precision pi,a,ta2,g0,x,y

      pi = 4d0*datan(1d0)
      a = 1d0/rsm
      ta2 = 2*a*a
      if (kmax.lt.0 .or. lmax.lt.0) return
      g0 = dexp(-a*a*r*r) * (a*a/pi)**1.5d0

C --- Do explicitly for k=0,1; See Eqs. 5.14 and 5.15 ---
      do  6  l = 0, lmax
        g(0,l) = ta2**l * g0
    6 continue
      if (kmax .eq. 0) return
      do  7  l = 0, lmax
        g(1,l) = ta2**(l+1)*(2d0*a*a*r*r-3-2*l) * g0
    7 continue

C --- Recursion for higher k; see Eq. 5.19 ---
      do    k = 2, kmax
      do    l = 0, lmax
        x = 2*(k-1)*(2*k+2*l-1)
        y = 4*k+2*l-1
        g(k,l) = ta2*(ta2*r*r*g(k-1,l) - y*g(k-1,l) - x*ta2*g(k-2,l))
      enddo
      enddo
      end
C      subroutine fmain
C      integer n0,kmax,lmax,k,l
C      parameter (n0=4,kmax=3,lmax=3)
C      double precision pkl(0:n0,0:lmax),gkl(0:n0,0:lmax)
C      double precision r,rsm,pi,a,g0
C
C      r = 1.7d0
C      rsm = 1d0
C      pi = 4d0*datan(1d0)
C      a = 1d0/rsm
C      g0 = dexp(-a*a*r*r) * (a*a/pi)**1.5d0
C
C      call radgkl(r,rsm,kmax,lmax,n0,gkl)
C      call radpkl(r,rsm,kmax,lmax,n0,pkl)
C
C      do  k = 0, kmax
C        print 333, (gkl(k,l), l=0,lmax)
C        print 333, (gkl(k,l)/g0, l=0,lmax)
C        print 333, (pkl(k,l), l=0,lmax)
C        print 333, (gkl(k,l)/pkl(k,l), l=0,lmax)
C        print *, ' '
C  333   format(5f14.8)
C      enddo
C      end

