! FCPP#define F90 1
subroutine hxpos(rsmh,rsmg,eh,kmax,nlmh,k0,c)
  !- Coefficients to expand smooth hankels at (0,0,0) into P_kl's.
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   rsmh  :vector of l-dependent smoothing radii of smoothed hankel
  !i         :rsmh must be specified for 1..ll(nlmh)+1
  !i   rsmg  :smoothing radius of gaussian
  !i   eh    :vector of l-dependent energies of smoothed Hankel
  !i         :eh must be specified for 1..ll(nlmh)+1
  !i   kmax  :polynomial cutoff
  !i   nlmh  :L-cutoff for smoothed Hankel functions and P_kL
  !i   k0    :leading dimension of coefficient array c
  !o Outputs
  !o   c     :structure constants c(k,M,L); see Remarks
  !o         :c(*,ilm) for rsmh(ll(ilm))<=0 are SET TO ZERO
  !r Remarks
  !r   Expansion is:  H_L = sum(k) c(k,L) * P_kL
  !r   As rsmg -> 0, expansion turns into a Taylor series of H_L.
  !r   As rsmg increases the error is spread out over a progressively
  !r   larger radius, thus reducing the error for larger r while
  !r   reducing accuracy at small r.
  !r
  !r    Only diagonal elements ilmh=ilmg are nonzero and returned.
  !r    Routine is equivalent to hxpml with ph=pg.
  !u Updates
  !u   18 May 00 Made rsmh,eh l-dependent
  !u   24 Apr 00 Adapted from nfp hxp_os.f
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: k0,kmax,nlmh
  double precision :: eh(1),rsmg,rsmh(1),c(0:k0,nlmh)
  ! ... Local parameters
  integer :: ndim,ktop0,ilm,k,l,ll,lmax,m,nm
  double precision :: a,dfact,eg,fac,factk,fpi,sig
  ! FCPP#if F90
  real(8),allocatable:: s(:,:)
  ! FCPP#else
  ! FCPP      parameter ( ndim=25, ktop0=25 )
  ! FCPP      double precision s(ndim,0:ktop0)
  ! FCPP#endif

  if (nlmh == 0) return
  fpi = 16d0*datan(1d0)

  ! FCPP#if F90
  ktop0 = kmax
  ndim  = nlmh
  allocate(s(ndim,0:ktop0))
  ! FCPP#else
  ! FCPP      if (kmax .gt. ktop0) call rxi('hxpos: increase ktop0, need',kmax)
  ! FCPP      if (nlmh .gt. ndim)   call rxi('hxpos: increase ndim, need', nlmh)
  ! FCPP#endif

  ! ... Make integrals with gaussians G_kL
  eg = 0d0
  call ghios(rsmg,rsmh,eg,eh,nlmh,kmax,ndim,s)

  ! ... Scale integrals to get coefficients of the P_kL
  a = 1d0/rsmg
  lmax = ll(nlmh)
  ilm = 0
  dfact = 1d0
  do  l = 0, lmax
     nm = 2*l+1
     do  m = 1, nm
        ilm = ilm+1
        factk = 1d0
        if (rsmh(l+1) > 0) then
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

  ! FCPP#if F90
  deallocate(s)
  ! FCPP#endif

end subroutine hxpos

