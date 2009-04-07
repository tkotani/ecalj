      subroutine mstrx2(e,dr,nlma,nlmb,ndim,cg,indxcg,jcg,cy,iop,s,sdot)
C- Calculate one-center expansion coefficents to j of h (struc. consts)
C  and their energy derivative
C ----------------------------------------------------------------
Ci Inputs:
Ci   e:   energy of Hankel function
Ci   dr:  vector R_source - R_expansion
Ci   nlma:Generate coefficients S_R'L',RL for L' < nlma
Ci   nlmb:Generate coefficients S_R'L',RL for L  < nlmb
Ci   ndim:row dimension of structure constant matrix s
Ci   cy  :normalization constants for spherical harmonics (sylmnc.f)
Ci   cg,indxcg,jcg:Clebsch Gordon coefficients, and indices (scg.f)
Ci   iop :1's digit
Ci        1: calculate s only
Ci        2: calculate sdot only
Ci        any other number: calculate both s and sdot
Ci       10's digit: nonzero => use OKA's conventions
Ci        1 strux to satisfy Andersen's definitions, with factor (-1)
Co Outputs
Co   s   :structure constant matrix S_R'L',RL
Co   sdot:Energy derivative of s
Cr Remarks
Cr
Cr  Expansion Theorem:
Cr  (1)  H_L(r-DR) = sum_L' J_L' (r) S_L'L (DR)
Cr       or
Cr       H_L(r-R) = sum_L' J_L' (r-R') S_L'L (R-R')
Cr       where DR=R-R', R=R_source; R'=R_expansion and
Cr  (2)  S_MK (DR) = 4 pi sum_L C_KLM (-E)^((k+m-l)/2) H_L (DR) (-1)^k
Cr
Cr  Note: the expansion theorem can also be written
Cr  (3)  H_RL(E,r) = \sum_L' J_R'L'(E,r) S_R'L',RL
Cr  where  H_RL(r) = H_L(r), and S_R'L',RL = S_L'L (R-R')
Cr
Cr   MSM definitions of solid Hankel, Bessel:
Cr    H_L  =      k^(l+1) n_l (kr) Y_L   (E > 0)
Cr    H_L  =  (i k)^(l+1) n_l (kr) Y_L   (E < 0)
Cr    J_L  =      k^(-l)  j_l (kr) Y_L   (E > 0)
Cr    J_L  =  (i k)^(-l)  j_l (kr) Y_L   (E < 0)
Cr
Cr  where  h_0 = Re e^(ikr)/r and j = sin(kr)/kr, k = sqrt(e), Im k >=0.
Cr  Generate w/ H_L = Y_L(-grad) h(r);   J_L = E^(-l) Y_L (-grad) j(r)
Cr
Cr  These satisfy:
Cr    (1)  cases e .ne. 0 and e .eq. 0 same equations
Cr    (2)  functions are real
Cr    (3)  S_KL (dr)  =  S^*_LK(-dr)
Cr    (4)  no extra factors (2l-1)!! etc as in Andersen's definitions.
Cr
Cr  Indices for L: middle is m=0; (z axis for m=0)
Cr                 after is real part of spherical harmonics: m=1,2, ...
Cr                 before is imaginary part m=1,2, ...
Cr  For s,p,d,f the indices 1-16 correspond to:
Cr  index   l    m     polynomial
Cr     1    0    0        1
Cr     2    1   -1        y
Cr     3    1    0        z
Cr     4    1    1        x
Cr     5    2   -2        xy
Cr     6    2   -1        yz
Cr     7    2    0        3z^2-1
Cr     8    2    1        xz
Cr     9    2    2        x^2-y^2
Cr     10   3   -3        y(3x^2-y^2)
Cr     11   3   -2        xyz
Cr     12   3   -1        y(5z^2-1)
Cr     13   3    0        z(5z^2-3)
Cr     14   3    1        x(5z^2-1)
Cr     15   3    2        z(x^2-y^2)
Cr     16   3    3        x(x^2-3y^2)
Cr
Cr  Andersen's definitions:
Cr  HOKA_L(k^2,r)  =  H_L / (2l-1)!!
Cr  JOKA_L(k^2,r)  = -J_L * (2l-1)!!/2
Cr
Cr  Thus Andersen structure constants differ by the ratio
Cr    S_RL,R'L' (OKA) / S_RL,R'L' =  2/((2l-1)!!(2l'-1)!!)
Cr  NB: Andersen defines H = -J S; mstrx2 returns -S (OKA)
Cr
Cr  Note: the e-> zero limit of the ss(sigma) sdot diverges as
Cr  1/sqrt(e).  This singularity is removed for case e=0, making
Cr  a discontinuity in Sdot at e=0.
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer ndim,nlma,nlmb,iop
      integer indxcg(1),jcg(1)
      double precision e
      double precision s(ndim,nlmb),sdot(ndim,nlmb),cy(*),cg(1),dr(3)
C Local parameters
      logical loka
      integer lmxx
      parameter (lmxx=11)
      integer icg,icg1,icg2,ii,mlm,indx,ipow,klm,l,lk,llm,lm,lmax,
     .        m,ll,klm0,sig(0:lmxx),iiop
      double precision fpi,rfac,sum,r2,sumd,gdot,xx,
     .  edot(0:lmxx),efac(0:lmxx),fac2l(0:lmxx),
     .  psi(-1:lmxx),phi(-1:lmxx),
     .  yl((lmxx+1)**2),hl((lmxx+1)**2),hdot((lmxx+1)**2)
      external ll,sylm

C ... Setup
      lmax = ll(nlma) + ll(nlmb)
      fac2l(0) = 1
      do  5  l = 1, lmax
    5 fac2l(l) = fac2l(l-1) * (2*l-1)
      iiop = mod(iop,10)
      loka = mod(iop/10,10) .ne. 0
      fpi = 16.d0*datan(1.d0)
      lmax = ll(nlma) + ll(nlmb)
      if (nlma .gt. ndim)
     .  call rxi('MSTRX2: nlma gt ndim: need ndim=',nlma)
      if (lmax .gt. lmxx)
     .  call rxi('MSTRX2: lmax too big, need lmax=',lmax)
      call sylm(dr,yl,lmax,r2)
      if (r2 .lt. 1d-10) return
C ... Bessel functions, Andersen's definitions
c     call besslr(e*r2,1,-1,lmax,phi,psi)
C ... Bessel functions, MSM definitions
      call besslr(e*r2,0,-1,lmax,phi,psi)
      mlm = 0
      rfac = dsqrt(r2)
C --- Generate H_L(DR) = 1/r^(2l+1) Psi_l  Y_lm ---
C     and psidot, defined so that Hdot_L = H_L * psidot(l)
      do  10  l = 0, lmax
C       gdot = psi(l-1)*r2/(4*l-2)
        gdot = psi(l-1)*r2/2
C   ... rfac = r^-(1+2*l); note yl has a factor r**l.
        rfac = rfac/r2
        do  12  m = -l, l
          mlm = mlm+1
C         xx = rfac*cy(mlm)*yl(mlm)*fac2l(l)
          xx = rfac*cy(mlm)*yl(mlm)
C     ... hl,hdot same as standard MSM IV-43.
          hl(mlm)   = psi(l)*xx
          hdot(mlm) = gdot * xx
   12   continue
   10 continue

C ... efac(l)=(-e)**l; sig(l)=(-)**l * sign of iop
      efac(0) = 1
      edot(0) = 0
      sig(0)  = 1
      do  1  l = 1, lmax
        efac(l) = -e*efac(l-1)
        edot(l) = -l*efac(l-1)
        sig(l)  = -sig(l-1)
    1 continue

      if (iiop .ne. 2) then

C --- Generate S_MK ---
      do  11  mlm = 1, nlma
        lm = ll(mlm)
        klm0 = mlm
C   ... Handles subblock nlmb < mlm <= nlma, if there is one
        if (mlm .gt. nlmb) klm0 = 1
        do  14  klm = klm0, nlmb
          lk = ll(klm)
          sum = 0d0
          ii = max0(mlm,klm)
          indx = (ii*(ii-1))/2 + min0(mlm,klm)
          icg1 = indxcg(indx)
          icg2 = indxcg(indx+1) - 1
          do  15  icg = icg1, icg2
            llm  = jcg(icg)
            ipow = (lm+lk-ll(llm))/2
            sum  = sum + cg(icg)*efac(ipow)*hl(llm)
   15     continue
          s(mlm,klm) = fpi*sig(lk)*sum
          if (loka) s(mlm,klm) = s(mlm,klm)*2/(fac2l(lm)*fac2l(lk))
          if (klm .le. nlma .and. mlm .le. nlmb)
     .      s(klm,mlm) = s(mlm,klm)*sig(lk)*sig(lm)
   14 continue
   11 continue

      endif

      if (iiop .ne. 1) then
C --- Generate Sdot_MK ---
      do  31  mlm = 1, nlma
        lm = ll(mlm)
        klm0 = mlm
C   ... Handles subblock nlmb < mlm <= nlma, if there is one
        if (mlm .gt. nlmb) klm0 = 1
        do  31  klm = klm0, nlmb
        lk = ll(klm)
        sumd = 0d0
        ii = max0(mlm,klm)
        indx = (ii*(ii-1))/2 + min0(mlm,klm)
        icg1 = indxcg(indx)
        icg2 = indxcg(indx+1) - 1
        do  32  icg = icg1, icg2
          llm = jcg(icg)
          ipow = (lm+lk-ll(llm))/2
          sumd = sumd + cg(icg)*
     .                  (efac(ipow)*hdot(llm)+edot(ipow)*hl(llm))
   32   continue
        sdot(mlm,klm) = fpi*sumd*sig(lk)
        if (loka) sdot(mlm,klm) = sdot(mlm,klm)*2/(fac2l(lm)*fac2l(lk))
        if (klm .le. nlma .and. mlm .le. nlmb)
     .  sdot(klm,mlm) = sdot(mlm,klm)*sig(lk)*sig(lm)
   31 continue
      endif

      end

