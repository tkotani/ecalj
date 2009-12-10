#define F90 1
      subroutine gklml(p,rsm,e,kmax,nlm,k0,gkl)
C- Bloch-sums of k,L-dependent gaussians
C ----------------------------------------------------------------------
Ci Inputs
Ci   p     :Function is centered at p
Ci   rsm   :smoothing radius
Ci   e     :gkL scaled by exp(e*rsm**2/4)
Ci   kmax  :polynomial cutoff
Ci   nlm   :L-cutoff for gkl
Ci   k0    :leading dimension of gkl
Co Outputs
Co   gkl   :Bloch-summed Gaussians
Cu Updates
Cu   24 Jan 07 Adapted from gklbl.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer k0,kmax,nlm
      double precision e,rsm,p(3)
      double precision gkl(0:k0,nlm)
C ... Local parameters
      integer ilm,k,ll,l,m,lmax
      double precision r1,r2,cfac
#if F90 | AUTO_ARRAY
      double precision wk(0:kmax,0:nlm-1),yl(nlm)
#else
      double precision wk(1,1),yl(1)
      call rx('gklml requires automatic arrays')
#endif

      if (nlm .eq. 0) return
      do  ilm = 1, nlm
        do  k = 0, kmax
          gkl(k,ilm) = 0d0
        enddo
      enddo

      lmax = ll(nlm)
      call ropyln(1,p(1),p(2),p(3),lmax,1,yl,r2)
      r1 = dsqrt(r2)
      call radgkl(r1,rsm,kmax,lmax,kmax,wk)
      if (e .ne. 0) then
        cfac = dexp(0.25d0*e*rsm*rsm)
        do  l = 0, lmax
          do  k = 0, kmax
            wk(k,l) = wk(k,l)*cfac
          enddo
        enddo
      endif
      ilm = 0
      do  l = 0, lmax
        do  m = -l, l
          ilm = ilm+1
          do  k = 0, kmax
            gkl(k,ilm) = gkl(k,ilm) + yl(ilm)*wk(k,l)
          enddo
        enddo
      enddo

      end

