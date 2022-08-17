subroutine vecpkl(r,rsm,nr,kmax,lmax,nrx,k0,wk,lrl,p,gp)! Vector of p_kl polynomials, or r^l p_kl
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   r     :vector of points
  !i   rsm   :smoothing radius
  !i   nr    :number of points
  !i   kmax  :make p for k=0..kmax
  !i   lmax  :make p for l=0..lmax
  !i   nrx   :leading dimension of p
  !i   k0    :second dimension of p
  !i   wk    :work array of length nr
  !i   lrl   :if 1s digit = 0, returns p_kl; otherwise returns r^l p_kl
  !i         :if 10s digit nonzero, returns gp; otherwise gp is not touched.
  !o Outputs
  !o   p     :radial part of spherical polynomials P_kL; see Remarks
  !o   gp    :radial derivative of p from l=0..lmax-1 (depending on lrl).
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
  !u Updates
  !u   22 Aug 01 bug fix for gp when kmax=0
  !u   25 Jan 00 veckl generates gp as well as p.
  ! ----------------------------------------------------------------------
  !     implicit none
  integer :: nr,kmax,lmax,nrx,k0,lrl
  double precision :: r(nrx),wk(nr),rsm,p(nrx,0:k0,0:*), &
       gp(nrx,0:k0,0:*)
  integer :: i,l,k
  double precision :: a,xx,xx2,xx3

  if (kmax < 0 .OR. lmax < 0) return
  if (kmax > k0) call rx('vecpkl: kmax gt k0')
  if (rsm <= 0) call rx('vecpkl: rsm <= 0')
  a = 1d0/rsm

  ! --- Set wk = 2*a**2*r**2 ---
  xx = 2*a*a
  do  6  i = 1, nr
     wk(i) = xx*r(i)**2
6 enddo

  ! --- Do explicitly for k=0,1 ---
  do    l = 0, lmax
     xx = a**l
     do    i = 1, nr
        p(i,0,l) = xx
     enddo
  enddo

  if (kmax > 0) then
     do   l = 0, lmax
        xx = a**l
        xx2 = 1/dble(2*l+3)
        do   i = 1, nr
           p(i,1,l)=xx*(wk(i)*xx2-1d0)
        enddo
     enddo
  endif

  ! --- Recursion for higher k ---
  do    k = 2, kmax
     xx3 = 2*(k-1)
     do    l = 0, lmax
        xx2 = (4*k+2*l-1)
        xx = 1/dble(2*k+2*l+1)
        do    i = 1, nr
           p(i,k,l) = xx*((wk(i)-xx2)*p(i,k-1,l) - xx3*p(i,k-2,l))
        enddo
     enddo
  enddo

  ! --- Radial derivative of p ---
  if (mod(lrl/10,10) /= 0) then

     !  ... Set wk = 2*a**2*r**2
     xx = 2*a*a
     do  16  i = 1, nr
        wk(i) = xx*r(i)
16   enddo

     do    k = 0, kmax
        do    l = 0, lmax-1
           xx2 = dble(2*k+2*l+3)/(a*(2*l+3))
           do    i = 1, nr
              gp(i,k,l) = wk(i)*(p(i,k,l) - xx2*p(i,k,l+1))
           enddo
        enddo
     enddo
  endif

  ! --- Scale by r^l if lrl nonzero ---
  if (mod(lrl,10) == 0) return
  do  50  i = 1, nr
     wk(i) = 1
50 enddo

  do  52  l = 1, lmax

     !   ... gP scales as  r*l gP +  l*r^l-1 P
     if (mod(lrl/10,10) /= 0 .AND. l < lmax) then
        do   k = 0, kmax
           do   i = 1, nr
              gp(i,k,l) = wk(i)*r(i)*gp(i,k,l) + l*wk(i)*p(i,k,l)
           enddo
        enddo
     endif

     do  54  i = 1, nr
        wk(i) = wk(i)*r(i)
54   enddo
     do   k = 0, kmax
        do   i = 1, nr
           p(i,k,l) = p(i,k,l)*wk(i)
        enddo
     enddo

52 enddo

end subroutine vecpkl

