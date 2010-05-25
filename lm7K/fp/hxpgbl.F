#define F90 1
      subroutine hxpgbl(ph,pg,q,rsmh,rsmg,eh,kmax,nlmh,nlmg,k0,ndimh,
     .ndimg,cg,indxcg,jcg,cy,slat,c,dc)
C- Coefficients to expand smooth bloch hankels and grads centered at ph
C  into a sum of polynomials P_kL centered at pg.
C ----------------------------------------------------------------------
Ci Inputs
Ci   ph    :Function is centered at ph; see Remarks
Ci   pg    :Function it expansed is at pg; see Remarks
Ci   q     :wave number for Bloch sum
Ci   rsmh  :smoothing radius of smoothed Hankel (l-dependent)
Ci   rsmg  :smoothing radius of gaussian
Ci   eh    :energies of smoothed Hankel (l-dependent)
Ci   kmax  :polynomial cutoff
Ci   nlmh  :L-cutoff for smoothed Hankel functions being expanded
Ci   nlmg  :L-cutoff for P_kL expansion
Ci   k0    :leading dimension of coefficient arrays c,dc
Ci   ndimg :second dimension of coefficient arrays c,dc
Ci   ndimh :htird dimension of coefficient arrays c,dc
Ci   cg    :Clebsch Gordon coefficients, stored in condensed form (scg.f)
Ci   indxcg:index for Clebsch Gordon coefficients
Ci   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
Ci   cy    :Normalization constants for spherical harmonics
Ci   slat  :struct containing information about the lattice
Co Outputs
Co   c     :cofficients c(k,M,L); see Remarks
Co   dc    :cofficients dc(k,M,L,1..3) (gradients of c); see Remarks
Co         :here k=0..kmax, M=1..nlmg, L=1..nlmh
Cr Remarks
Cr   Expansion is:  H_L(r-ph) = sum(kM)  c(k,M,L) * P_kM(r-pg)
Cr             grad H_L(r-ph) = sum(kM) dc(k,M,L) * P_kM(r-pg)
Cr   As rsmg -> 0, expansion turns into a Taylor series of H_L.
Cr   As rsmg increases the error is spread out over a progressively
Cr   larger radius, thus reducing the error for larger r while
Cr   reducing accuracy at small r.
Cr
Cr   Grads are wrt to ph; take negative for grad wrt to pg.
Cu Updates
Cu   25 May 00 Made rsmh,eh l-dependent
Cu   25 May 00 Adapted from nfp hxpg_bl.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer k0,kmax,ndimg,ndimh,nlmg,nlmh,jcg(1),indxcg(1)
      double precision eh(1),rsmg,rsmh(1),ph(3),pg(3),cg(1),cy(1),q(3),
     .slat(1)
      double complex c(0:k0,ndimg,ndimh),dc(0:k0,ndimg,ndimh,3)
C ... Local parameters
      integer ndim1,ndim2,ktop0,ilmg,ilmh,k,l,ll,lmaxg,m,nm
      double precision a,dfact,eg,fac,factk,fpi
#if F90
      complex(8),allocatable:: s(:,:,:),ds(:,:,:,:)
#else
      parameter (ndim1=49, ndim2=36, ktop0=20)
      double complex s(ndim1,ndim2,0:ktop0),ds(ndim1,ndim2,0:ktop0,3)
#endif

      if (nlmg .eq. 0 .or. nlmh .eq. 0) return
      fpi = 16d0*datan(1d0)
C     Memory allocation
#if F90
      ndim1 = nlmg
      ndim2 = nlmh
      ktop0 = kmax
      allocate(s(ndim1,ndim2,0:ktop0),ds(ndim1,ndim2,0:ktop0,3))
#endif
      if (kmax .gt. ktop0) call rxi('hxpgbl: increase ktop0, need',kmax)
      if (nlmg .gt. ndim1) call rxi('hxpgbl: increase ndim1, need',nlmg)
      if (nlmh .gt. ndim2) call rxi('hxpgbl: increase ndim2, need',nlmh)

C ... Integrals of Hankels with Gaussians
      eg = 0d0
      call ghigbl(pg,ph,q,rsmg,rsmh,eg,eh,nlmg,nlmh,kmax,ndim1,ndim2,
     .ktop0,cg,indxcg,jcg,cy,slat,s,ds)

C ... Scale to get coefficients of the PkL
      a = 1d0/rsmg
      lmaxg = ll(nlmg)

      ilmg = 0
      dfact = 1d0
      do  1  l = 0, lmaxg
        nm = 2*l+1
        do  2  m = 1, nm
          ilmg = ilmg+1
          factk = 1d0
          do  3  k = 0, kmax
            fac = fpi/ ((4*a*a)**k * a**l * factk * dfact)
            do  5  ilmh = 1, nlmh
              c(k,ilmg,ilmh) = s(ilmg,ilmh,k)*fac
              dc(k,ilmg,ilmh,1) = -ds(ilmg,ilmh,k,1)*fac
              dc(k,ilmg,ilmh,2) = -ds(ilmg,ilmh,k,2)*fac
              dc(k,ilmg,ilmh,3) = -ds(ilmg,ilmh,k,3)*fac
  5         continue
            factk = factk*(k+1)
  3       continue
  2     continue
        dfact = dfact*(2*l+3)
  1   continue

#if F90
      deallocate(s,ds)
#endif
      end

