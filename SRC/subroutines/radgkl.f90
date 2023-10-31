subroutine radgkl(r,rsm,kmax,lmax,n0,g)  ! Make radial parts g_kl of G_kL, divided by r**l.
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
  !r     G_kL  = phi_kl (a^2/pi)^(3/2) exp(-a^2 r^2) r^l Y_L = phi_kl r^l Y_L g_00
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
  implicit none
  integer :: kmax,lmax,n0,l,k
  real(8):: r,rsm,g(0:n0,0:lmax), a,ta2,g0,x,y
  real(8),parameter:: pi = 4d0*datan(1d0)
  a = 1d0/rsm
  ta2 = 2d0*a*a
  if (kmax < 0 .OR. lmax < 0) return
  g0 = dexp(-a*a*r*r) * (a*a/pi)**1.5d0 
  g(0,0:lmax) = [(ta2**l, l=0,lmax)]*g0 ! --- Do explicitly for k=0,1; See Eqs. 5.14 and 5.15 ---
  if (kmax == 0) return
  g(1,0:lmax) = [(ta2**(l+1)*(ta2*r*r-3-2*l),l=0,lmax)]*g0
  do   k = 2, kmax !Recursion for higher k; see Eq. 5.19 ---
     g(k,0:lmax) = [(ta2*(ta2*r*r*g(k-1,l) - (4*k+2*l-1)*g(k-1,l) - 2*(k-1)*(2*k+2*l-1)*ta2*g(k-2,l)),l=0,lmax)]
  enddo
end subroutine radgkl
