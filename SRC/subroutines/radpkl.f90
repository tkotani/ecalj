subroutine radpkl(r,rsm,kmax,lmax,n0,p)
  !- Radial parts of polynomials P_kL /  r**l
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   r     :radius
  !i   rsm   :smoothing radius
  !i   kmax  :make p for k=0..kmax
  !i   lmax  :make p for l=0..lmax
  !i   n0    :leading dimension of p
  !o Outputs
  !o   p     :radial part of spherical polynomials P_kL; see Remarks
  !o         :radpkl produces p such that P_kL = p_kl r**l Y_L
  !r Remarks
  !r   P_kL are polyonomials orthogonal in the following sense:
  !r                                          (4a^2)^k a^l k! (2l+1)!!
  !r    int P_kL G_k'L' = delta_kk'*delta_ll'  ----------------------
  !r                                                    4pi
  !r   and are defined in J. Math. Phys. 39, 3393 (1988).
  !r   Combining eqns 12.7 and 5.19 in that paper, we obtain
  !r    p_kl = a**l / (2a**2)^(k+l) (2l+1)!! / (2k+2l+1)!! phi_kl
  !r    p_0l = a**l
  !r    p_1l = a**l (2*(ar)**2/(2l+3) - 1)
  !r    p_kl = [(2*(ar)**2 - (4k+2l-1))p_k-1,l - 2(k-1)p_k-2,l]
  !r           / (2k+2l+1)
  ! ----------------------------------------------------------------------
  !     implicit none
  integer :: kmax,lmax,n0
  double precision :: r,rsm,p(0:n0,0:lmax)
  double precision :: tar2,a
  integer :: l,k

  if (kmax < 0 .OR. lmax < 0) return
  a = 1d0/rsm
  tar2 = 2d0*a*a*r*r

  ! --- Do explicitly for k=0,1 ---
  do  7  l = 0, lmax
     p(0,l) = a**l
     if (kmax >= 1) p(1,l) = (a**l)*(tar2/(2*l+3)-1d0)
7 enddo

  ! ---- Recursion for higher k ---
  do    k = 2, kmax
     do    l = 0, lmax
        p(k,l) = ((tar2-(4*k+2*l-1))*p(k-1,l)-2*(k-1)*p(k-2,l)) &
             /(2*k+2*l+1)
     enddo
  enddo
end subroutine radpkl

