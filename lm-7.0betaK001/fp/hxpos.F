#define F90 1
      subroutine hxpos(rsmh,rsmg,eh,kmax,nlmh,k0,c)
C- Coefficients to expand smooth hankels at (0,0,0) into P_kl's.
C ----------------------------------------------------------------------
Ci Inputs
Ci   rsmh  :vector of l-dependent smoothing radii of smoothed hankel
Ci         :rsmh must be specified for 1..ll(nlmh)+1
Ci   rsmg  :smoothing radius of gaussian
Ci   eh    :vector of l-dependent energies of smoothed Hankel
Ci         :eh must be specified for 1..ll(nlmh)+1
Ci   kmax  :polynomial cutoff
Ci   nlmh  :L-cutoff for smoothed Hankel functions and P_kL
Ci   k0    :leading dimension of coefficient array c
Co Outputs
Co   c     :structure constants c(k,M,L); see Remarks
Co         :c(*,ilm) for rsmh(ll(ilm))<=0 are SET TO ZERO
Cr Remarks
Cr   Expansion is:  H_L = sum(k) c(k,L) * P_kL
Cr   As rsmg -> 0, expansion turns into a Taylor series of H_L.
Cr   As rsmg increases the error is spread out over a progressively
Cr   larger radius, thus reducing the error for larger r while
Cr   reducing accuracy at small r.
Cr
Cr    Only diagonal elements ilmh=ilmg are nonzero and returned.
Cr    Routine is equivalent to hxpml with ph=pg.
Cu Updates
Cu   18 May 00 Made rsmh,eh l-dependent
Cu   24 Apr 00 Adapted from nfp hxp_os.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer k0,kmax,nlmh
      double precision eh(1),rsmg,rsmh(1),c(0:k0,nlmh)
C ... Local parameters
      integer ndim,ktop0,ilm,k,l,ll,lmax,m,nm
      double precision a,dfact,eg,fac,factk,fpi,sig
#if F90
      real(8),allocatable:: s(:,:)
#else
      parameter ( ndim=25, ktop0=25 )
      double precision s(ndim,0:ktop0)
#endif

      if (nlmh .eq. 0) return
      fpi = 16d0*datan(1d0)

#if F90
      ktop0 = kmax
      ndim  = nlmh
      allocate(s(ndim,0:ktop0))
#else
      if (kmax .gt. ktop0) call rxi('hxpos: increase ktop0, need',kmax)
      if (nlmh .gt. ndim)   call rxi('hxpos: increase ndim, need', nlmh)
#endif

C ... Make integrals with gaussians G_kL
      eg = 0d0
      call ghios(rsmg,rsmh,eg,eh,nlmh,kmax,ndim,s)

C ... Scale integrals to get coefficients of the P_kL
      a = 1d0/rsmg
      lmax = ll(nlmh)
      ilm = 0
      dfact = 1d0
      do  l = 0, lmax
        nm = 2*l+1
        do  m = 1, nm
          ilm = ilm+1
          factk = 1d0
          if (rsmh(l+1) .gt. 0) then
            do  k = 0, kmax
              fac = (4*a*a)**k * a**l * factk * dfact
              sig = 1d0
              c(k,ilm) = s(ilm,k)*(fpi/fac)*sig
              factk = factk*(k+1)
            enddo
          else
            do  k = 0, kmax
              c(k,ilm) = 0
            enddo
          endif
        enddo
        dfact = dfact*(2*l+3)
      enddo

#if F90
      deallocate(s)
#endif

      end

