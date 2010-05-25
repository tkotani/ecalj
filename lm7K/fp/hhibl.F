      subroutine hhibl(mode,p1,p2,q,rsm1,rsm2,e1,e2,nlm1,nlm2,kmax,
     .ndim1,ndim2,cg,indxcg,jcg,cy,slat,s)
C- Integrals between smooth Hankels with k-th power of Laplace operator.
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1's digit (not implemented: always vectors)
Ci         :0 rsm1,rsm2,e1,e2 are scalars
Ci         :1 rsm1,rsm2,e1,e2 are l-dependent vectors
Ci         :10's digit
Ci         :1: do not calculate strux for any (l1,l2) pairs for
Ci         :   if rsm(l1) or rsm(l2) is zero.
Ci   p1    :first center
Ci   p2    :second center
Ci   q     :wave number for Bloch sum
Ci   rsm1  :smoothing radii of Hankels at p1 (l-dependent)
Ci   rsm2  :smoothing radii of Hankels at p2 (l-dependent)
Ci   e1    :energies of smooth Hankels at p1 (l-dependent)
Ci   e2    :energies of smooth Hankels at p2 (l-dependent)
Ci   nlm1  :L-cutoff for functions at p1
Ci   nlm2  :L-cutoff for functions at p2
Ci   kmax  :cutoff in power of Laplace operator
Ci   ndim1 :leading dimension of s
Ci   ndim2 :second dimension of s
Ci   cg    :Clebsch Gordon coefficients (scg.f)
Ci   indxcg:index for Clebsch Gordon coefficients
Ci   jcg   :L q.n. for the C.G. coefficients (scg.f)
Ci   cy    :Normalization constants for spherical harmonics
Ci   slat  :struct containing information about the lattice
Co Outputs
Co   s     :integrals; see Remarks
Cr Remarks
Cr   s(L,M) contains integral of H_L^*(r-p1) (laplace)^k H_M(r-p2)
Cr   Row L corresponds to p1 and col M corresponds to p2.
Cr   Strux s(L,M) is computed for dr=p1-p2
Cr   See JMP 39, 3393, Section 10
Cu Updates
Cu   8 Jun 00  Added 10s digit mode.  New arg list.
Cu   19 May 00 Made rsm1,e1,rsm2,e2 l-dependent.  Elements for which
Cu             rsm1 =0 or rsm2 = 0 are not computed.
Cu   18 May 00 Adapted from nfp hhi_bl.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,jcg(1),indxcg(1),nlm1,nlm2,kmax,ndim1,ndim2
      double precision p1(3),p2(3),q(3),cg(1),cy(1),slat(1),
     .rsm1(0:*),rsm2(0:*),e1(0:*),e2(0:*)
      double complex s(ndim1,ndim2,0:kmax)
C ... Local parameters
      integer m,lmx1,lmx2,ll,l1,l2,k,jlm,ilm,lm11,lm21,lm12,lm22,l1t,l2t
      double precision dr(3)

      if (nlm1 .eq. 0 .or. nlm2 .eq. 0) return
      do  1  m = 1, 3
        dr(m) = p1(m)-p2(m)
    1 continue
      lmx1 = ll(nlm1)
      lmx2 = ll(nlm2)

      do  3  k = 0, kmax
        do  3  jlm = 1, nlm2
          do  3  ilm = 1, nlm1
            s(ilm,jlm,k) = dcmplx(0d0,0d0)
    3 continue

      l1t = -1
      do  20  l1 = 0, lmx1
        if (l1 .le. l1t) goto 20
        call gtbsl2(l1,lmx1,e1,rsm1,l1t)
C       l1t = l1

        l2t = -1
        do  22  l2 = 0, lmx2
          if (l2 .le. l2t) goto 22
          call gtbsl2(l2,lmx2,e2,rsm2,l2t)
C         l2t = l2

          lm11 = l1**2+1
          lm12 = (l1t+1)**2
          lm21 = l2**2+1
          lm22 = (l2t+1)**2
          if (mode/10.eq.1 .and. rsm1(l1)*rsm2(l2) .eq. 0) goto 22
          call phhibl(dr,q,rsm1(l1),rsm2(l2),e1(l1),e2(l2),lm11,lm12,
     .    lm21,lm22,kmax,ndim1,ndim2,cg,indxcg,jcg,cy,slat,s)
   22   continue
   20 continue
      end
      subroutine phhibl(dr,q,rsm1,rsm2,e1,e2,mlm1,nlm1,mlm2,nlm2,kmax,
     .ndim1,ndim2,cg,indxcg,jcg,cy,slat,s)
C- Integrals between smooth Hankels with k-th power of Laplace operator.
C ----------------------------------------------------------------------
Ci Inputs
Ci   dr    :p1-p2
Ci   q     :wave number for Bloch sum
Ci   rsm1  :smoothing radius of Hankels at p1
Ci   rsm2  :smoothing radius of Hankels at p2
Ci   e1    :energy of smooth Hankels at p1
Ci   e2    :energy of smooth Hankels at p2
Ci   nlm1  :L-cutoff for functions at p1
Ci   nlm2  :L-cutoff for functions at p2
Ci   kmax  :cutoff in power of Laplace operator
Ci   ndim1 :leading dimension of s
Ci   ndim2 :second dimension of s
Ci   cg    :Clebsch Gordon coefficients (scg.f)
Ci   indxcg:index for Clebsch Gordon coefficients
Ci   jcg   :L q.n. for the C.G. coefficients (scg.f)
Ci   cy    :Normalization constants for spherical harmonics
Ci   slat  :struct containing information about the lattice
Co Outputs
Co   s     :integrals; see Remarks
Cr Remarks
Cr   s(L,M) contains integral of H_L^*(r-p1) (laplace)^k H_M(r-p2)
Cu Updates
Cu   18 May 00 Adapted from nfp hhi_bl.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer jcg(1),indxcg(1),mlm1,nlm1,mlm2,nlm2,kmax,ndim1,ndim2
      double precision dr(3),q(3),cg(1),cy(1),slat(1),rsm1,
     .rsm2,e1,e2
      double complex s(ndim1,ndim2,0:kmax)
C ... Local parameters
      integer nlm0,ktop0,icg,icg1,icg2,ii,ilm,ilm1,ilm2,indx,ip,k,
     .ktop,l1,l2,ll,lm,lmax1,lmax2,lmaxx,nlmx
      parameter( nlm0=100, ktop0=10 )
      double precision fpi,e,fac,fac1,fac2,gam1,gam2,gamx,rsmx
      double complex hkl1(0:ktop0,nlm0),hkl2(0:ktop0,nlm0),
     .hsm(nlm0),hsmp(nlm0)

      fpi = 16d0*datan(1.d0)
      gam1 = 0.25d0*rsm1*rsm1
      gam2 = 0.25d0*rsm2*rsm2
      gamx = gam1+gam2
      rsmx = 2d0*dsqrt(gamx)
      lmax1 = ll(nlm1)
      lmax2 = ll(nlm2)
      lmaxx = lmax1+lmax2
      nlmx = (lmaxx+1)**2
      ktop = max0(lmax1,lmax2)+kmax
      if (nlmx .gt. nlm0) call rxi('increase nlm0 in hhibl need',nlmx)
      if (ktop .gt. ktop0) call rx('hhibl: increase ktop0')

C ... Set up functions for connecting vector p1-p2
      if (dabs(e1-e2) .gt. 1d-5) then
        fac1 = dexp(gam2*(e2-e1))/(e1-e2)
        fac2 = dexp(gam1*(e1-e2))/(e2-e1)
        call hklbl(dr,rsmx,e1,q,ktop,nlmx,ktop0,cy,slat, hkl1)
        call hklbl(dr,rsmx,e2,q,ktop,nlmx,ktop0,cy,slat, hkl2)
        do  4  ilm = 1, nlmx
          do  5  k = 0, ktop
            hkl1(k,ilm) = fac1*hkl1(k,ilm) + fac2*hkl2(k,ilm)
    5     continue
    4   continue
      else
        e = .5d0*(e1+e2)
        call hklbl(dr,rsmx,e,q,ktop,nlmx,ktop0,cy,slat, hkl2)
        call hsmbl(dr,rsmx,e,q,lmaxx,cy,slat, hsm,hsmp)
        do  6  ilm = 1, nlmx
          hkl1(0,ilm) = hsmp(ilm) - gamx*hsm(ilm)
          do  7  k = 1, ktop
            hkl1(k,ilm) = -e*hkl1(k-1,ilm) - hkl2(k-1,ilm)
    7     continue
    6   continue

      endif

C ... Combine with Clebsch-Gordan coefficients
      do  11  ilm1 = mlm1, nlm1
        l1 = ll(ilm1)
        do  11  ilm2 = mlm2, nlm2
          l2 = ll(ilm2)
          ii = max0(ilm1,ilm2)
          indx = (ii*(ii-1))/2 + min0(ilm1,ilm2)
          icg1 = indxcg(indx)
          icg2 = indxcg(indx+1)-1
          do  11  icg = icg1, icg2
            ilm = jcg(icg)
            lm = ll(ilm)
            k = (l1+l2-lm)/2
            fac = fpi*(-1d0)**l1*cg(icg)
            do  12  ip = 0, kmax
              s(ilm1,ilm2,ip) = s(ilm1,ilm2,ip) + fac*hkl1(k+ip,ilm)
   12       continue
   11 continue

      end

