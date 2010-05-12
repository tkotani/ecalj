      subroutine vecpkl(r,rsm,nr,kmax,lmax,nrx,k0,wk,lrl,p,gp)
C- Vector of p_kl polynomials, or r^l p_kl
C ----------------------------------------------------------------------
Ci Inputs
Ci   r     :vector of points
Ci   rsm   :smoothing radius
Ci   nr    :number of points
Ci   kmax  :make p for k=0..kmax
Ci   lmax  :make p for l=0..lmax
Ci   nrx   :leading dimension of p
Ci   k0    :second dimension of p
Ci   wk    :work array of length nr
Ci   lrl   :if 1s digit = 0, returns p_kl; otherwise returns r^l p_kl
Ci         :if 10s digit nonzero, returns gp; otherwise gp is not touched.
Co Outputs
Co   p     :radial part of spherical polynomials P_kL; see Remarks
Co   gp    :radial derivative of p from l=0..lmax-1 (depending on lrl).
Cr Remarks
Cr   P_kL are polyonomials orthogonal in the following sense:
Cr                                          (4a^2)^k a^l k! (2l+1)!!
Cr    int P_kL G_k'L' = delta_kk'*delta_ll'  ----------------------
Cr                                                    4pi
Cr   and are defined in J. Math. Phys. 39, 3393 (1988).
Cr   Combining eqns 12.7 and 5.19 in that paper, we obtain
Cr    p_kl = a**l / (2a**2)^(k+l) (2l+1)!! / (2k+2l+1)!! phi_kl
Cr    p_0l = a**l
Cr    p_1l = a**l (2*(ar)**2/(2l+3) - 1)
Cr    p_kl = [(2*(ar)**2 - (4k+2l-1))p_k-1,l - 2(k-1)p_k-2,l]
Cr           / (2k+2l+1)
Cu Updates
Cu   22 Aug 01 bug fix for gp when kmax=0
Cu   25 Jan 00 veckl generates gp as well as p.
C ----------------------------------------------------------------------
C     implicit none
      integer nr,kmax,lmax,nrx,k0,lrl
      double precision r(nrx),wk(nr),rsm,p(nrx,0:k0,0:1),
     .  gp(nrx,0:k0,0:1)
      integer i,l,k
      double precision a,xx,xx2,xx3

      if (kmax .lt. 0 .or. lmax .lt. 0) return
      if (kmax .gt. k0) call rx('vecpkl: kmax gt k0')
      if (rsm .le. 0) call rx('vecpkl: rsm <= 0')
      a = 1d0/rsm

C --- Set wk = 2*a**2*r**2 ---
      xx = 2*a*a
      do  6  i = 1, nr
        wk(i) = xx*r(i)**2
    6 continue

C --- Do explicitly for k=0,1 ---
      do  10  l = 0, lmax
        xx = a**l
      do  10  i = 1, nr
        p(i,0,l) = xx
   10 continue

      if (kmax .gt. 0) then
        do  20  l = 0, lmax
          xx = a**l
          xx2 = 1/dble(2*l+3)
        do  20  i = 1, nr
          p(i,1,l)=xx*(wk(i)*xx2-1d0)
   20   continue
      endif

C --- Recursion for higher k ---
      do  30  k = 2, kmax
        xx3 = 2*(k-1)
      do  30  l = 0, lmax
        xx2 = (4*k+2*l-1)
        xx = 1/dble(2*k+2*l+1)
      do  30  i = 1, nr
        p(i,k,l) = xx*((wk(i)-xx2)*p(i,k-1,l) - xx3*p(i,k-2,l))
   30 continue

C --- Radial derivative of p ---
      if (mod(lrl/10,10) .ne. 0) then

C  ... Set wk = 2*a**2*r**2
        xx = 2*a*a
        do  16  i = 1, nr
          wk(i) = xx*r(i)
   16   continue

        do  40  k = 0, kmax
        do  40  l = 0, lmax-1
          xx2 = dble(2*k+2*l+3)/(a*(2*l+3))
        do  40  i = 1, nr
          gp(i,k,l) = wk(i)*(p(i,k,l) - xx2*p(i,k,l+1))
   40   continue

      endif

C --- Scale by r^l if lrl nonzero ---
      if (mod(lrl,10) .eq. 0) return
      do  50  i = 1, nr
        wk(i) = 1
   50 continue

      do  52  l = 1, lmax

C   ... gP scales as  r*l gP +  l*r^l-1 P
        if (mod(lrl/10,10) .ne. 0 .and. l .lt. lmax) then
          do  53  k = 0, kmax
          do  53  i = 1, nr
            gp(i,k,l) = wk(i)*r(i)*gp(i,k,l) + l*wk(i)*p(i,k,l)
   53     continue
        endif

        do  54  i = 1, nr
          wk(i) = wk(i)*r(i)
   54   continue
        do  56  k = 0, kmax
        do  56  i = 1, nr
          p(i,k,l) = p(i,k,l)*wk(i)
   56   continue

   52 continue

      end

