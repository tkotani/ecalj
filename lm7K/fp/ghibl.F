#define F90 1
      subroutine ghibl(pg,ph,q,rsmg,rsmh,eg,eh,nlmg,nlmh,kmax,
     .  ndim1,ndim2,cg,indxcg,jcg,cy,slat,s)
C- Block of integrals between smooth hankels and gaussians with some power
C  of the laplace operator.
C ----------------------------------------------------------------------
Ci Inputs
Ci   ph    :Function is centered at ph; see Remarks
Ci   pg    :Function it expansed is at pg; see Remarks
Ci   q     :wave number for Bloch sum
Ci   rsmg  :smoothing radius of gaussian
Ci   rsmh  :vector of l-dependent smoothing radii of smoothed hankel
Ci         :rsmh must be specified for 1..ll(nlmh)+1
Ci   eg    :gkL scaled by exp(e*rsm**2/4)
Ci   eh    :vector of l-dependent energies of smoothed Hankel
Ci         :eh must be specified for 1..ll(nlmh)+1
Ci   nlmg  :L-cutoff for P_kL expansion
Ci   nlmh  :L-cutoff for smoothed Hankel functions
Ci   kmax  :polynomial cutoff
Ci   ndim1 :leading dimension of s
Ci   ndim2 :second dimension of s
Ci   cg    :Clebsch Gordon coefficients, stored in condensed form (scg.f)
Ci   indxcg:index for Clebsch Gordon coefficients
Ci   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
Ci   cy    :Normalization constants for spherical harmonics
Ci   slat  :struct containing information about the lattice
Co Outputs
Co   s     :integrals of gaussian and Hankels; see Remarks
Cr Remarks
Cr   s(L,M,k) contains integral of G_L^*(r-pg) (lap)^k H_M(r-ph)
Cr   See J. Math. Phys. {\bf 39},3393 (1998), Eq. 8.4
Cu Updates
Cu   18 May 00 Made rsmh,eh l-dependent
Cu   24 Apr 00 Adapted from nfp ghi_bl.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nlmg,nlmh,kmax,ndim1,ndim2,jcg(1),indxcg(1)
      double precision rsmg,rsmh(1),eg,eh(1)
      double precision ph(3),pg(3),cg(1),cy(1),q(3),slat(1)
      double complex s(ndim1,ndim2,0:1)
C ... Local parameters
      integer nlm0,ktop0,icg,icg1,icg2,ii,ilg,ilh,ilm,indx,ip,jlm,k,
     .  ktop,lg,lh,ll,lm,lmaxg,lmaxh,lmaxx,m,nlmx,l1,l2,ilm1,ilm2
#if F90
      complex(8),allocatable:: hkl(:,:)
#else
      parameter (nlm0=144, ktop0=21)
      double complex hkl(0:ktop0,nlm0)
#endif
      double precision ee,fac,gamg,gamh,rsmx,dr(3),e,rsm

      if (nlmh .eq. 0 .or. nlmg .eq. 0) return

C ... rsmh- and eh- independent setup
      do  1  m = 1, 3
        dr(m) = pg(m)-ph(m)
    1 continue
      lmaxh = ll(nlmh)
      lmaxg = ll(nlmg)
      lmaxx = lmaxg+lmaxh
      nlmx = (lmaxx+1)**2
      ktop = max0(lmaxg,lmaxh)+kmax

#if F90
      ktop0 = ktop
      nlm0  = nlmx
      allocate( hkl(0:ktop0,nlm0))
#else
      if (nlmx .gt. nlm0)  call rxi('ghibl: increase nlm0 to',nlmx)
      if (ktop .gt. ktop0) call rxi('ghibl: increase ktop0 to',ktop)
#endif
      do  3  k = 0, kmax
      do  3  jlm = 1, nlmh
      do  3  ilm = 1, nlmg
        s(ilm,jlm,k) = dcmplx(0d0,0d0)
    3 continue

C --- Loop over sequences of l with a common rsm,e ---
      l2 = -1
      do  20  l1 = 0, lmaxh
        if (l1 .le. l2) goto 20
        call gtbsl2(l1,lmaxh,eh,rsmh,l2)
        rsm  = rsmh(l1+1)
        e    = eh(l1+1)
        if (rsm .le. 0 .or. e .gt. 0) goto 20
        ilm1 = l1**2+1
        ilm2 = (l2+1)**2
        lmaxx= lmaxg+l2
        nlmx = (lmaxx+1)**2
        gamh = 0.25d0*rsm*rsm
        gamg = 0.25d0*rsmg*rsmg
        rsmx = 2d0*dsqrt(gamg+gamh)
        ktop = max0(lmaxg,l2)+kmax
        call hklbl(dr,rsmx,e,q,ktop,nlmx,ktop0,cy,slat,hkl)

C   ... Combine with Clebsch-Gordan coefficients
        ee = dexp(gamg*(eg-e))
        do  11  ilg = 1, nlmg
          lg = ll(ilg)
        do  11  ilh = ilm1, ilm2
          lh = ll(ilh)
          ii = max0(ilg,ilh)
          indx = (ii*(ii-1))/2 + min0(ilg,ilh)
          icg1 = indxcg(indx)
          icg2 = indxcg(indx+1)-1
        do  11  icg = icg1, icg2
          ilm = jcg(icg)
          lm = ll(ilm)
          k = (lg+lh-lm)/2
          fac = ee*(-1d0)**lg*cg(icg)
          do  12  ip = 0, kmax
            s(ilg,ilh,ip) = s(ilg,ilh,ip) + fac*hkl(k+ip,ilm)
   12     continue
   11   continue
   20 continue

#if F90
      deallocate(hkl)
#endif

      end

