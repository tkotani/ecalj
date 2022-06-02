subroutine radgkl(r,rsm,kmax,lmax,n0,g)
  !- Make radial parts g_kl of G_kL, divided by r**l.
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   r     :radius
  !i   rsm   :smoothing radius
  !i   kmax  :make g_kl for k=0...kmax
  !i   lmax  :make g_kl for l=0...lmax
  !i   n0    :leading dimension of g
  !o Outputs
  !o   g     :radial part of generalized gaussians G_kL; see Remarks
  !o         :radgkl produces g such that G_kL = g(k,l) r**l Y_L
  !r Remarks
  !r   Definition:  G_kL = (lap)**k YL(-grad) g0 = g_kl*r^l*Y_L
  !r   Energy factor beta = dexp(e/(4*a*a)) is not included here.
  !r   See J. Math. Phys. {\bf 39},3393 (1998), Section V:
  !r
  !r     G_kL  = phi_kl (a^2/pi)^(3/2) exp(-a^2 r^2) r^l Y_L
  !r           = phi_kl r^l Y_L g_00
  !r     G_0L  = (2a^2)^l r^l Y_L g_00
  !r
  !r   where (Eq. 12.3)
  !r
  !r    phi_kl = (-1)^k (2a^2)^(k+l) 2^k k! L_k^(l+1/2)(a^2 r^2)
  !r    L_k    = generalized Laguerre polynomial
  !r    phi_0l = (2a^2)^l
  !r
  !r   also
  !r
  !r                          a^l (2l+1)!!
  !r      p_kl = phi_kl ------------------------ r^l
  !r                    (2a^2)^(k+l) (2k+2l+1)!!
  !r
  !u Updates
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: kmax,lmax,n0
  double precision :: r,rsm,g(0:n0,0:lmax)
  ! ... Local parameters
  integer :: l,k
  double precision :: pi,a,ta2,g0,x,y

  pi = 4d0*datan(1d0)
  a = 1d0/rsm
  ta2 = 2*a*a
  if (kmax < 0 .OR. lmax < 0) return
  g0 = dexp(-a*a*r*r) * (a*a/pi)**1.5d0

  ! --- Do explicitly for k=0,1; See Eqs. 5.14 and 5.15 ---
  do  6  l = 0, lmax
     g(0,l) = ta2**l * g0
6 enddo
  if (kmax == 0) return
  do  7  l = 0, lmax
     g(1,l) = ta2**(l+1)*(2d0*a*a*r*r-3-2*l) * g0
7 enddo

  ! --- Recursion for higher k; see Eq. 5.19 ---
  do    k = 2, kmax
     do    l = 0, lmax
        x = 2*(k-1)*(2*k+2*l-1)
        y = 4*k+2*l-1
        g(k,l) = ta2*(ta2*r*r*g(k-1,l) - y*g(k-1,l) - x*ta2*g(k-2,l))
     enddo
  enddo
end subroutine radgkl
!      subroutine fmain
!      integer n0,kmax,lmax,k,l
!      parameter (n0=4,kmax=3,lmax=3)
!      double precision pkl(0:n0,0:lmax),gkl(0:n0,0:lmax)
!      double precision r,rsm,pi,a,g0

!      r = 1.7d0
!      rsm = 1d0
!      pi = 4d0*datan(1d0)
!      a = 1d0/rsm
!      g0 = dexp(-a*a*r*r) * (a*a/pi)**1.5d0

!      call radgkl(r,rsm,kmax,lmax,n0,gkl)
!      call radpkl(r,rsm,kmax,lmax,n0,pkl)

!      do  k = 0, kmax
!        print 333, (gkl(k,l), l=0,lmax)
!        print 333, (gkl(k,l)/g0, l=0,lmax)
!        print 333, (pkl(k,l), l=0,lmax)
!        print 333, (gkl(k,l)/pkl(k,l), l=0,lmax)
!        print *, ' '
!  333   format(5f14.8)
!      enddo
!      end

