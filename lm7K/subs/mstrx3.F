      subroutine mstrx3(kap2,dr,nlsqri,nlsqrj,nlsqr,cg,indxcg,jcg,cy,
     .                  iop,s,sdot)
C- Calculate one-center expansion coefficents to J of J (struc. consts)
C  and their energy derivative
C ----------------------------------------------------------------------
Ci Inputs:
Ci   kap2  :(kappa*average WS-radius)**2
Ci   dr    :vector connecting sphere at R' to that at R (units of avw)
Ci   nlsqri:Generate coefficients S_R'L',RL for L' < nlsqri
Ci   nlsqrj:Generate coefficients S_R'L',RL for L  < nlsqrj
Ci   nlsqr :row dimension of structure constant matrix s
Ci   cg    :Clebsch Gordon coefficients
Ci   indxcg:index for cg
Ci   jcg   :jcg(icg)=lm"
Ci   cy    :normalization constants for spherical harmonics (see SYLMNC)
Ci   iop   :1: calculate s
Ci          2: calculate sdot
Co Outputs: (depending on iop)
Co   s     :structure constant matrix
Co          multiplied with a factor -1
Co   sdot  :(kappa*avw)**2 derivative of structure constant matrix
Co          multiplied with a factor -1
Cr Remarks:
Cr  Conventions for structure constants: (in the following equations
Cr  the lengths and energies are given in atomic units, whereas
Cr  in the program r is scaled by w (=avw) and kap2 means {kw}^2.
Cr  The solutions of the Helmholtz equation are defined as (OKA)
Cr  (k=kappa):
Cr
Cr  H_L(k^2,r)  = - (kw)^(l+1)/(2l-1)!!   n_l(kr) Y_L(r)
Cr  J_L(k^2,r)  =   (kw)^(-l) *(2l-1)!!/2 j_l(kr) Y_L(r)
Cr
Cr  they depend on kappa only through kappa^2 and apply both for
Cr  positive and negative kappa^2.
Cr  n_0(kr) = - cos(kr)/kr, j_0(kr) = sin(kr)/kr, ... are the
Cr  usual Bessel and Neumann Functions of fractional order.
Cr
Cr  These satisfy:
Cr    (1)  cases kap2 .ne. 0 and kap2 .eq. 0 same equations
Cr    (2)  functions are real
Cr    (3)  S_{R'L',RL}  =  S^*_{RL,R'L'}
Cr
Cr  Expansion Theorem: J_RL(r) = J_L(r-R)
Cr
Cr  (1)  J_{RL}(k^2,r) = - J_{R'L'}(k^2,r) * S_{R'L',RL}
Cr  where
Cr  (2)  S_R'L',RL  = - Sum_l" 4 pi C_{LL'L"} i^{-l+l'-l"}  *
Cr
Cr                          -l+l'+l"
Cr                     2 (kw)       (2l-1)!!
Cr                  * -----------------------  J_L"(k^2,R-R')
Cr                      (2l'-1)!!  (2l"-1)!!
Cr
Cr  Indices for L: middle is m=0; (z axis for m=0)
Cr          after is real part of spherical harmonics: m=1,2, ...
Cr          before is imaginary part m=1,2, ... (* +/- sqrt(2))
Cr  For s,p,d the indices 1-9 correspond to:
Cr         1     2     3     4     5     6       7     8     9
Cr         s     y     z     x     xy    yz  3z^2-r^2  xz  x^2-y^2
Cr   the structure constants differ by the ratio
Cr   S_RL,R'L'(MSM) / S_RL,R'L'(OKA) = -(2l'-1)!!/(2l-1)!!
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer nlsqr,nlsqri,nlsqrj,iop,indxcg(*),jcg(*)
      double precision kap2,s(nlsqr,*),sdot(nlsqr,*),cy(*),cg(*),dr(3)
C Local variables:
      integer nlmax,lmaxx,nlsqrx
      parameter (nlmax=20,lmaxx=2*nlmax-2,nlsqrx=(lmaxx+1)*(lmaxx+1))
      integer icg,icg1,icg2,ii,ilm,indx,ipow,jlm,klm,
     .        l,li,lj,lk,ll,lmax,lgunit,mk,sig(0:lmaxx)
      double precision edot(0:lmaxx),efac(0:lmaxx),fac2l(0:lmaxx),
     .                 fi(0:lmaxx+1),fpi,fdot,gi(0:lmaxx+1),r2,
     .                 bdot(nlsqrx),bl(nlsqrx),sum,sumd,ylm(nlsqrx)
      parameter(fpi=12.5663706143591725d0)
C External calls:
      external errmsg,ll,lgunit,sylm,bessl2
C Intrinsic functions:
      intrinsic max0,min0

      if (nlsqri.gt.nlsqr) call rx('MSTRX3: nlsqri gt nlsqr')
      lmax = ll(nlsqri) + ll(nlsqrj)
      if (lmax .gt. lmaxx)
     .  call rxi('MSTRX3: increase lmaxx, need',lmax)

C --- A table of (2l-1)!!
c     data fac2l /1,1,3,15,105,945,10395,135135,2027025,34459425/
      fac2l(0) = 1.d0
      do  10  l = 1, lmax
   10 fac2l(l) = fac2l(l-1) * (2*l-1)

      call sylm(dr,ylm,lmax,r2)
      call bessl2(kap2*r2,0,lmax+1,fi(0),gi(0))
      klm = 0
      do  20  lk = 0, lmax
        fdot  = -fi(lk+1)*r2/(4*lk+2)
        do  22  mk = -lk, lk
          klm = klm+1
          bl  (klm) = fi(lk)*cy(klm)*ylm(klm)/fac2l(lk)
          bdot(klm) = fdot  *cy(klm)*ylm(klm)/fac2l(lk)
   22   continue
   20 continue

C --- efac(l)=(-kap2)^l
C --- sig(l)=(-)^l
      efac(0) = 1.d0
      edot(0) = 0.d0
      sig(0)  = 1.d0
      do  30  l = 1, lmax
        efac(l) = -kap2*efac(l-1)
        edot(l) = -l   *efac(l-1)
        sig(l)  = -sig(l-1)
   30 continue

      if (iop .ne. 2) then
C ----- Generate S_IJ
        do  40  ilm = 1, nlsqri
          li = ll(ilm)
          do  42  jlm = 1, nlsqrj
            lj = ll(jlm)
            ii = max0(ilm,jlm)
            indx = (ii*(ii-1))/2 + min0(ilm,jlm)
            icg1 = indxcg(indx)
            icg2 = indxcg(indx+1) - 1
            sum = 0.d0
            do  44  icg = icg1, icg2
              klm  = jcg(icg)
              lk   = ll(klm)
              ipow = (li-lj+lk)/2
              sum  = sum + cg(icg)*efac(ipow)*bl(klm)*sig(lk)
   44       continue
            s(ilm,jlm) = fpi*sum*2.d0*fac2l(lj)/fac2l(li)
   42     continue
   40   continue
      endif
      if (iop .ne. 1) then
C ----- Generate Sdot_IJ
        do  50  ilm = 1, nlsqri
          li = ll(ilm)
          do  52  jlm = 1, nlsqrj
            lj = ll(jlm)
            ii = max0(ilm,jlm)
            indx = (ii*(ii-1))/2 + min0(ilm,jlm)
            icg1 = indxcg(indx)
            icg2 = indxcg(indx+1) - 1
            sumd = 0d0
            do  54  icg = icg1, icg2
              klm  = jcg(icg)
              lk   = ll(klm)
              ipow = (li-lj+lk)/2
              sumd = sumd + cg(icg)*sig(lk)*
     .                      (efac(ipow)*bdot(klm)+edot(ipow)*bl(klm))
   54       continue
            sdot(ilm,jlm) = fpi*sumd*2.d0*fac2l(lj)/fac2l(li)
   52     continue
   50   continue
      endif
      end

