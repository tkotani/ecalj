! FCPP#define F90 1
subroutine hxpbl(ph,pg,q,rsmh,rsmg,eh,kmax,nlmh,nlmg,k0,ndim, &
     cg,indxcg,jcg,cy,c) !,slat
  !- Coefficients to expand smooth bloch hankels centered at ph
  !  into a sum of polynomials P_kL centered at pg.
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   ph    :Function is centered at ph; see Remarks
  !i   pg    :Function it expansed is at pg; see Remarks
  !i   q     :wave number for Bloch sum
  !i   rsmh  :vector of l-dependent smoothing radii of smoothed hankel
  !i         :rsmh must be specified for 1..ll(nlmh)+1
  !i   rsmg  :smoothing radius of gaussian
  !i   eh    :vector of l-dependent energies of smoothed Hankel
  !i         :eh must be specified for 1..ll(nlmh)+1
  !i   kmax  :polynomial cutoff
  !i   nlmh  :L-cutoff for smoothed Hankel functions being expanded
  !i   nlmg  :L-cutoff for P_kL expansion
  !i   k0    :leading dimension of coefficient array c
  !i   ndim  :second dimension of coefficient array c
  !i   cg    :Clebsch Gordon coefficients, stored in condensed form (scg.f)
  !i   indxcg:index for Clebsch Gordon coefficients
  !i   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
  !i   cy    :Normalization constants for spherical harmonics
  !i   slat  :struct containing information about the lattice
  !o Outputs
  !o   c     :cofficients c(k,M,L); see Remarks
  !o         :here k=0..kmax, M=1..nlmg, L=1..nlmh
  !r Remarks
  !r   Expansion is:  H_L(r-ph) = sum(kM) c(k,M,L) * P_kM(r-pg)
  !r   See J. Math. Phys. {\bf 39},3393 (1998), Eq. 12.17
  !r
  !r   As rsmg -> 0, expansion turns into a Taylor series of H_L.
  !r   As rsmg increases the error is spread out over a progressively
  !r   larger radius, thus reducing the error for larger r while
  !r   reducing accuracy at small r.
  !b Bugs
  !b   Doesn't properly handle case rsmh<0
  !u Updates
  !u   18 May 00 Made rsmh,eh l-dependent
  !u   24 Apr 00 Adapted from nfp hxp_bl.f
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: k0,kmax,ndim,nlmg,nlmh,jcg(1),indxcg(1)
  double precision :: eh(1),rsmg,rsmh(1),ph(3),pg(3),cg(1),cy(1),q(3)
  !     .slat(1)
  double complex c(0:k0,ndim,1)
  ! ... Local parameters
  integer :: ndim1,ndim2,ktop0,ilmg,ilmh,k,l,ll,lmaxg,m,nm
  double precision :: a,dfact,eg,fac,factk,fpi
  ! FCPP#if F90
  complex(8),allocatable:: s(:,:,:)
  ! FCPP#else
  ! FCPP      parameter (ndim1=49, ndim2=36, ktop0=20)
  ! FCPP      double complex s(ndim1,ndim2,0:ktop0)
  ! FCPP#endif

  if (nlmg == 0 .OR. nlmh == 0) return
  fpi = 16d0*datan(1d0)
  !     Memory allocation
  ! FCPP#if F90
  ndim1 = nlmg
  ndim2 = nlmh
  ktop0 = kmax
  allocate(s(ndim1,ndim2,0:ktop0))
  ! FCPP#endif

  if (kmax > ktop0) call rxi('hxpbl: increase ktop0, need',kmax)
  if (nlmg > ndim1) call rxi('hxpbl: increase ndim1, need',nlmg)
  if (nlmh > ndim2) call rxi('hxpbl: increase ndim2, need',nlmh)
  if (nlmg > ndim1) call rxi('hxpbl: increase ndim1, need',nlmg)

  ! ... Integrals of gaussians and smoothed Hankels
  eg = 0d0
  call ghibl(pg,ph,q,rsmg,rsmh,eg,eh,nlmg,nlmh,kmax,ndim1,ndim2,cg, &
       indxcg,jcg,cy,s)!,slat

  ! ... Scale to get coefficients of the P_kL
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
           fac = fpi / ((4*a*a)**k * a**l * factk * dfact)
           do  5  ilmh = 1, nlmh
              c(k,ilmg,ilmh) = s(ilmg,ilmh,k)*fac
5          enddo
           factk = factk*(k+1)
3       enddo
2    enddo
     dfact = dfact*(2*l+3)
1 enddo

  ! FCPP#if F90
  deallocate(s)
  ! FCPP#endif
end subroutine hxpbl

