module m_smhankel !Bloch sum of smooth Hankel, Gaussians.
  ! JMP39:
  ! Bott, E., M. Methfessel, W. Krabs, and P. C. Schmidt.
  ! “Nonsingular Hankel Functions as a New Basis for Electronic Structure Calculations.”
  ! Journal of Mathematical Physics 39, no. 6 (June 1, 1998): 3393–3425.
  ! https://doi.org/doi:10.1063/1.532437.
  public hxpbl, hxpgbl, hhigbl, hhibl,hhugbl, hgugbl, ggugbl  !*bl means blochsum 
  private
contains
subroutine hhugbl(mode,p1,p2,rsm1,rsm2,e1,e2,nlm1,nlm2,ndim1,ndim2,wk,dwk, s,ds)
  use m_lmfinit,only: rv_a_ocy,rv_a_ocg, iv_a_oidxcg, iv_a_ojcg
  use m_lattic,only:lat_vol
  !- Estatic energy integrals between Bloch Hankels, and gradients.
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   mode  :0 rsm1,rsm2,e1,e2 are scalars
  !i         :1 rsm1,rsm2,e1,e2 are l-dependent
  !i   p1    :first center
  !i   p2    :second center
  !i   rsm1  :smoothing radius of Hankels at p1 (l-dependent)
  !i   rsm2  :smoothing radius of Hankels at p2 (l-dependent)
  !i   e1    :energy  of Hankels at p1 (l-dependent)
  !i   e2    :energy  of Hankels at p2 (l-dependent)
  !i   nlm1  :L-max for  Hankels at p1
  !i   nlm2  :L-max for  Hankels at p2
  !i   ndim1 :leading dimensions of s,ds
  !i   ndim2 :second dimensions of s,ds
  !i   slat  :struct containing information about the lattice
  !i   wk    :work space of same size as s
  !i   dwk   :work space of same size as ds
  !r Remarks
  !r   Gradient is wrt p1; use -ds for grad wrt p2.
  !u Updates
  !u   23 May 00 Made rsm1,e1,rsm2,e2 l-dependent
  !u   22 Apr 00 Adapted from nfp hhug_bl.f
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: mode,nlm1,nlm2,ndim1,ndim2
  real(8):: rsm1(0:*) , rsm2(0:*) , e1(0:*) , e2(0:*) , p1(3) , p2(3)
  double complex s(ndim1,ndim2),ds(ndim1,ndim2,3), wk(ndim1,ndim2),dwk(ndim1,ndim2,3)
  integer:: kmax , kdim , i2 , i1 , lmxx , nlmx , l , ll
  parameter (lmxx=6,nlmx=(lmxx+1)**2)
  double precision :: add,pi,q(3),fpi,vol,gam1,gam2,xx1,xx2
  double precision :: zer(0:lmxx),bet1(nlmx),fac(nlmx)
  data q /0d0,0d0,0d0/
  data zer /lmxx*0d0,0d0/
  pi = 4d0*datan(1d0)
  vol=lat_vol
  kmax = 0
  kdim = 0
  if (nlm1 > nlmx) call rx('hhugbl: increase lmxx')
  call hhigbl ( mode , p1 , p2 , q , rsm1 , rsm2 , e1 , e2 , nlm1 &
       , nlm2 , kmax , ndim1 , ndim2 , kdim , rv_a_ocg , iv_a_oidxcg &
       , iv_a_ojcg , rv_a_ocy ,  s , ds ) !slat ,
  call hhigbl ( mode , p1 , p2 , q , rsm1 , rsm2 , e1 , zer , nlm1 &
       , nlm2 , kmax , ndim1 , ndim2 , kdim , rv_a_ocg , iv_a_oidxcg &
       , iv_a_ojcg , rv_a_ocy , wk , dwk ) !slat ,
  do  i2 = 1, nlm2
     if (mode == 0) then
        l = 0
     else
        l = ll(i2)
     endif
     bet1(i2) = dexp(e2(l)*rsm2(l)*rsm2(l)/4d0)
     fac(i2) = 8d0*pi/e2(l)
  enddo
  do  i2 = 1, nlm2
     do  i1 = 1, nlm1
        s(i1,i2)    = fac(i2)*(s(i1,i2)    - bet1(i2)*wk(i1,i2))
        ds(i1,i2,1) = fac(i2)*(ds(i1,i2,1) - bet1(i2)*dwk(i1,i2,1))
        ds(i1,i2,2) = fac(i2)*(ds(i1,i2,2) - bet1(i2)*dwk(i1,i2,2))
        ds(i1,i2,3) = fac(i2)*(ds(i1,i2,3) - bet1(i2)*dwk(i1,i2,3))
     enddo
  enddo
  ! ... Extra term for l1=l2=0
  fpi = 4d0*pi
  gam1 = 0.25d0*rsm1(0)*rsm1(0)
  gam2 = 0.25d0*rsm2(0)*rsm2(0)
  xx1 = fpi*dexp(gam1*e1(0))/(vol*e1(0))
  xx2 = fpi*dexp(gam2*e2(0))/(vol*e2(0))
  add = -2*xx1*xx2*vol/e2(0)
  s(1,1) = s(1,1) + add
end subroutine hhugbl

subroutine hgugbl(p1,p2,rsm1,rsm2,e1,nlm1,nlm2,ndim1,ndim2,  s,ds)
  use m_lmfinit,only: rv_a_ocy,rv_a_ocg, iv_a_oidxcg, iv_a_ojcg
  !- Estatic energy integrals between Bloch Hankels and gaussians, and grad
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   slat  :struct containing information about the lattice
  !i   p1    :first center
  !i   p2    :second center
  !i   rsm1  :smoothing radius of Hankels at p1
  !i   rsm2  :smoothing radius of gaussians at p2
  !i   e1    :energy  of Hankels at p1
  !i   nlm1  :L-max for  Hankels at p1
  !i   nlm2  :L-max for  gaussians at p2
  !i   ndim1 :leading dimensions of s,ds
  !i   ndim2 :second dimensions of s,ds
  !i   slat  :struct containing information about the lattice
  !o Outputs
  !o   s     :integrals between Bloch Hankels and gaussians
  !o   ds    :gradient of s; see Remarks
  !r Remarks
  !r   Gradient is wrt p1; use -ds for grad wrt p2.
  !u Updates
  !u   22 Apr 00 Adapted from nfp hgug_bl.f
  ! ----------------------------------------------------------------------
  implicit none
  integer :: nlm1,nlm2,ndim1,ndim2
  real(8):: rsm1 , rsm2 , p1(3) , p2(3) , e1
  complex(8):: s(ndim1,ndim2),ds(ndim1,ndim2,3)
  integer:: kmax , kdim , ilm2 , ilm1
  double precision :: q(3),e2
  q=[0d0,0d0,0d0]
  kmax = 0
  kdim = 0
  e2 = 0d0
  call hhigbl ( 0 , p1 , p2 , q , [rsm1] , [rsm2] , [e1] , [e2] , nlm1 &
       , nlm2 , kmax , ndim1 , ndim2 , kdim , rv_a_ocg , iv_a_oidxcg &
       , iv_a_ojcg , rv_a_ocy ,  s , ds )!slat ,
  do  ilm2 = 1, nlm2
     do  ilm1 = 1, nlm1
        s(ilm1,ilm2) = 2d0*s(ilm1,ilm2)
        ds(ilm1,ilm2,1) = 2d0*ds(ilm1,ilm2,1)
        ds(ilm1,ilm2,2) = 2d0*ds(ilm1,ilm2,2)
        ds(ilm1,ilm2,3) = 2d0*ds(ilm1,ilm2,3)
     enddo
  enddo
end subroutine hgugbl

subroutine hhigbl(mode,p1,p2,q,rsm1,rsm2,e1,e2,nlm1,nlm2,kmax, &
     ndim1,ndim2,k0,cg,indxcg,jcg,cy,s,ds) 
  !- Integrals between smooth hankels with k-th power of Laplace operator
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   mode  :1's digit
  !i         :0 rsm1,rsm2,e1,e2 are scalars
  !i         :1 rsm1,rsm2,e1,e2 are l-dependent vectors
  !i         :10's digit
  !i         :1: do not calculate strux for any (l1,l2) pairs for
  !i         :   if rsm(l1) or rsm(l2) is zero.
  !i   p1    :first center
  !i   p2    :second center
  !i   q     :wave number for Bloch sum
  !i   rsm1  :smoothing radii of Hankels at p1 (l-dependent)
  !i   rsm2  :smoothing radii of Hankels at p2 (l-dependent)
  !i   e1    :energies of smooth Hankels at p1 (l-dependent)
  !i   e2    :energies of smooth Hankels at p2 (l-dependent)
  !i   nlm1  :L-cutoff for functions at p1
  !i   nlm2  :L-cutoff for functions at p2
  !i   kmax  :cutoff in power of Laplace operator
  !i   ndim1 :leading dimensions of s,ds
  !i   ndim2 :second dimensions of s,ds
  !i   k0    :dimensions s,ds
  !i   cg    :Clebsch Gordon coefficients (scg.f)
  !i   indxcg:index for Clebsch Gordon coefficients
  !i   jcg   :L q.n. for the C.G. coefficients (scg.f)
  !i   cy    :Normalization constants for spherical harmonics
  !o Outputs
  !o   s     :integrals between smooth Bloch Hankels; see Remarks
  !o   ds    :gradient of s; see Remarks
  !r Remarks
  !r  s(L,M,k) contains integral of H_L^*(r-p1) (laplace)^k H_M(r-p2)
  !r  Gradient is wrt p1; use -ds for grad wrt p2.
  !u Updates
  !u   8 Jun 00  Added 10s digit mode.
  !u   23 May 00 Made rsm1,e1,rsm2,e2 l-dependent
  !u   22 Apr 00 Adapted from nfp hhig_bl.f
  !u   28 Apr 97 adapted to handle case q=0 and e1 or e2 zero.
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: mode,jcg(1),indxcg(1),nlm1,nlm2,kmax,ndim1,ndim2,k0
  real(8):: p1(3) , p2(3) , q(3) , cg(1) , cy(1) , rsm1(0:*) , &
       rsm2(0:*) , e1(0:*) , e2(0:*)

  double complex s(ndim1,ndim2,0:k0),ds(ndim1,ndim2,0:k0,3)
  ! ... Local parameters
  integer :: ilm,jlm,k,l1,l1t,l2,l2t,ll,lm11,lm12,lm21,lm22,lmx1,lmx2,m
  double precision :: dr(3)

  if (nlm1 == 0 .OR. nlm2 == 0) return

  do  1  m = 1, 3
     dr(m) = p1(m)-p2(m)
1 enddo
  lmx1 = ll(nlm1)
  lmx2 = ll(nlm2)

  do    k = 0, kmax
     do    jlm = 1, nlm2
        do    ilm = 1, nlm1
           s(ilm,jlm,k)    = dcmplx(0d0,0d0)
           ds(ilm,jlm,k,1) = dcmplx(0d0,0d0)
           ds(ilm,jlm,k,2) = dcmplx(0d0,0d0)
           ds(ilm,jlm,k,3) = dcmplx(0d0,0d0)
        enddo
     enddo
  enddo

  l1t = -1
  do  20  l1 = 0, lmx1
     if (l1 <= l1t) goto 20
     if (mod(mode,10) == 0) then
        l1t = lmx1
     else
        call gtbsl2(l1,lmx1,e1,rsm1,l1t)
     endif
     !       l1t = l1

     l2t = -1
     do  22  l2 = 0, lmx2
        if (l2 <= l2t) goto 22
        if (mod(mode,10) == 0) then
           l2t = lmx2
        else
           call gtbsl2(l2,lmx2,e2,rsm2,l2t)
        endif
        !         l2t = l2

        lm11 = l1**2+1
        lm12 = (l1t+1)**2
        lm21 = l2**2+1
        lm22 = (l2t+1)**2
        if (mode/10 == 1 .AND. rsm1(l1)*rsm2(l2) == 0) goto 22
        call phhigb(dr,q,rsm1(l1),rsm2(l2),e1(l1),e2(l2),lm11,lm12, &
             lm21,lm22,kmax,ndim1,ndim2,k0,cg,indxcg,jcg,cy,s,ds) !slat,
22   enddo
20 enddo
end subroutine hhigbl

subroutine phhigb(dr,q,rsm1,rsm2,e1,e2,mlm1,nlm1,mlm2,nlm2, &
     kmax,ndim1,ndim2,k0,cg,indxcg,jcg,cy,s,ds)!slat,
  use m_lattic,only:lat_vol
  !- Integrals between smooth hankels with k-th power of Laplace operator
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   dr    :p1-p2
  !i   q     :wave number for Bloch sum
  !i   rsm1  :smoothing radius of Hankels at p1
  !i   rsm2  :smoothing radius of Hankels at p2
  !i   e1    :energy of smooth Hankels at p1
  !i   e2    :energy of smooth Hankels at p2
  !i   mlm1  :lower L-cutoff for functions at p1
  !i   nlm1  :upper L-cutoff for functions at p1
  !i   mlm2  :lower L-cutoff for functions at p2
  !i   nlm2  :upper L-cutoff for functions at p2
  !i   kmax  :cutoff in power of Laplace operator
  !i   ndim1 :leading dimensions of s,ds
  !i   ndim2 :second dimensions of s,ds
  !i   k0    :dimensions s,ds
  !i   cg    :Clebsch Gordon coefficients (scg.f)
  !i   indxcg:index for Clebsch Gordon coefficients
  !i   jcg   :L q.n. for the C.G. coefficients (scg.f)
  !i   cy    :Normalization constants for spherical harmonics
  !i   slat  :struct containing information about the lattice
  !o Outputs
  !o   s     :integrals between smooth Bloch Hankels; see Remarks
  !o   ds    :gradient of s; see Remarks
  !r Remarks
  !r  s(L,M,k) contains integral of H_L^*(r-p1) (laplace)^k H_M(r-p2)
  !r  Gradient is wrt p1; use -ds for grad wrt p2.
  !u Updates
  !u   22 Apr 00 Adapted from nfp hhig_bl.f
  !u   28 Apr 97 adapted to handle case q=0 and e1 or e2 zero.
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: jcg(1),indxcg(1),nlm1,nlm2,kmax,ndim1,ndim2,k0,mlm1,mlm2
  real(8):: q(3) , cg(1) , cy(1) , dr(3) , rsm1 , rsm2 , e1 , e2

  double complex s(ndim1,ndim2,0:k0),ds(ndim1,ndim2,0:k0,3)
  ! ... Local parameters
  integer :: nlm0,ktop0,lmax1,ll,lmax2,lmaxx,nlmx,nlmxp1,ktop,ktopp1, &
       k,ilm,kz,kx1,kx2,ky1,ky2,ilm1,l1,ilm2,l2,ii,indx,icg1,icg2, &
       icg,lm,ip
  double precision :: gam1,fpi,gam2,gamx,rsmx,qq,fac1,fac2,e,cz,cx1, &
       cx2,cy1,cy2,fac,add,vol,dgets
  parameter (nlm0=100, ktop0=10)
  double complex hkl1(0:ktop0,nlm0),hkl2(0:ktop0,nlm0), &
       ghkl(0:ktop0,nlm0,3),hsm(nlm0),hsmp(nlm0)
  vol = lat_vol
  fpi = 16d0*datan(1.d0)
  gam1 = 0.25d0*rsm1*rsm1
  gam2 = 0.25d0*rsm2*rsm2
  gamx = gam1+gam2
  rsmx = 2d0*dsqrt(gamx)
  qq = dsqrt(q(1)*q(1)+q(2)*q(2)+q(3)*q(3))
  lmax1 = ll(nlm1)
  lmax2 = ll(nlm2)
  lmaxx = lmax1+lmax2
  nlmx = (lmaxx+1)**2
  nlmxp1 = (lmaxx+2)**2
  ktop = max0(lmax1,lmax2)+kmax
  ktopp1 = ktop+1
  if (nlmxp1 > nlm0)  call rxi('hhigbl: need nlm0 ge',nlmxp1)
  if (ktopp1 > ktop0) call rxi('hhigbl: need ktop0 ge',ktopp1)
  ! --- Set up functions for connecting vector p2-p1 ---
  if (dabs(e1-e2) > 1d-5) then
     fac1 = dexp(gam2*(e2-e1))/(e1-e2)
     fac2 = dexp(gam1*(e1-e2))/(e2-e1)
     if (dabs(e1) > 1d-6) then
        call hklbl(dr,rsmx,e1,q,ktopp1,nlmxp1,ktop0,cy,hkl1) !slat,
     else
        if (qq > 1d-6) call rx('hhigbl: e1=0 only allowed if q=0')
        call fklbl(dr,rsmx,ktopp1,nlmxp1,ktop0,cy,hkl1) !slat,
     endif
     if (dabs(e2) > 1d-6) then
        call hklbl(dr,rsmx,e2,q,ktopp1,nlmxp1,ktop0,cy,hkl2) !slat,
     else
        if (qq > 1d-6) call rx('hhigbl: e2=0 only allowed if q=0')
        call fklbl(dr,rsmx,ktopp1,nlmxp1,ktop0,cy,hkl2) !slat,
     endif
     do  4  ilm = 1, nlmxp1
        do  5  k = 0, ktopp1
           hkl1(k,ilm) = fac1*hkl1(k,ilm) + fac2*hkl2(k,ilm)
5       enddo
4    enddo
  else
     e = .5d0*(e1+e2)
     if (qq < 1d-6 .AND. dabs(e) < 1d-6) &
          call rx('hhigbl: case q=0 and e1=e2=0 not available')
     call hklbl(dr,rsmx,e,q,ktopp1,nlmxp1,ktop0,cy,hkl2) !slat,
     call hsmbl(dr,rsmx,e,q,lmaxx+1,cy,hsm,hsmp) !slat,
     do  6  ilm = 1, nlmxp1
        hkl1(0,ilm) = hsmp(ilm) - gamx*hsm(ilm)
        do  7  k = 1, ktopp1
           hkl1(k,ilm) = -e*hkl1(k-1,ilm) - hkl2(k-1,ilm)
7       enddo
6    enddo
  endif
  ! ... Make gradients using CGs for p functions
  call dpzero(ghkl, 2*(ktop0+1)*nlm0*3)
  do  ilm = 1, nlmx
     call scglp1(ilm,kz,cz,kx1,kx2,cx1,cx2,ky1,ky2,cy1,cy2)
     do  k = 0, ktop
        ghkl(k,ilm,1) = ghkl(k,ilm,1)-cx1*hkl1(k,kx1)-cx2*hkl1(k,kx2)
        ghkl(k,ilm,2) = ghkl(k,ilm,2)-cy1*hkl1(k,ky1)-cy2*hkl1(k,ky2)
        ghkl(k,ilm,3) = ghkl(k,ilm,3)-cz*hkl1(k,kz)
        if (ilm <= lmaxx*lmaxx) then
           ghkl(k,kx1,1) = ghkl(k,kx1,1) - cx1*hkl1(k+1,ilm)
           ghkl(k,kx2,1) = ghkl(k,kx2,1) - cx2*hkl1(k+1,ilm)
           ghkl(k,ky1,2) = ghkl(k,ky1,2) - cy1*hkl1(k+1,ilm)
           ghkl(k,ky2,2) = ghkl(k,ky2,2) - cy2*hkl1(k+1,ilm)
           ghkl(k,kz,3)  = ghkl(k,kz,3)  - cz *hkl1(k+1,ilm)
        endif
     enddo
  enddo
  ! ... Combine with Clebsch-Gordan coefficients
  do  1111  ilm1 = mlm1, nlm1
     l1 = ll(ilm1)
     do  111  ilm2 = mlm2, nlm2
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
              ds(ilm1,ilm2,ip,1) = ds(ilm1,ilm2,ip,1) + fac*ghkl(k+ip,ilm,1)
              ds(ilm1,ilm2,ip,2) = ds(ilm1,ilm2,ip,2) + fac*ghkl(k+ip,ilm,2)
              ds(ilm1,ilm2,ip,3) = ds(ilm1,ilm2,ip,3) + fac*ghkl(k+ip,ilm,3)
12         enddo
11      enddo
111  enddo
1111 enddo
  ! ... Extra term when e1 or e2 is zero (but not both)
  if (mlm1 == 1 .AND. mlm2 == 1) then
     add = 0d0
     if (dabs(e2) < 1d-6) add = fpi*dexp(gam1*e1)/(vol*e1*e1)
     if (dabs(e1) < 1d-6) add = fpi*dexp(gam2*e2)/(vol*e2*e2)
     s(1,1,0) = s(1,1,0) + add
  endif
end subroutine phhigb
  
subroutine hhibl(mode,p1,p2,q,rsm1,rsm2,e1,e2,nlm1,nlm2,kmax, &
     ndim1,ndim2,cg,indxcg,jcg,cy,s) !,slat
  !- Integrals between smooth Hankels with k-th power of Laplace operator.
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   mode  :1's digit (not implemented: always vectors)
  !i         :0 rsm1,rsm2,e1,e2 are scalars
  !i         :1 rsm1,rsm2,e1,e2 are l-dependent vectors
  !i         :10's digit
  !i         :1: do not calculate strux for any (l1,l2) pairs for
  !i         :   if rsm(l1) or rsm(l2) is zero.
  !i   p1    :first center
  !i   p2    :second center
  !i   q     :wave number for Bloch sum
  !i   rsm1  :smoothing radii of Hankels at p1 (l-dependent)
  !i   rsm2  :smoothing radii of Hankels at p2 (l-dependent)
  !i   e1    :energies of smooth Hankels at p1 (l-dependent)
  !i   e2    :energies of smooth Hankels at p2 (l-dependent)
  !i   nlm1  :L-cutoff for functions at p1
  !i   nlm2  :L-cutoff for functions at p2
  !i   kmax  :cutoff in power of Laplace operator
  !i   ndim1 :leading dimension of s
  !i   ndim2 :second dimension of s
  !i   cg    :Clebsch Gordon coefficients (scg.f)
  !i   indxcg:index for Clebsch Gordon coefficients
  !i   jcg   :L q.n. for the C.G. coefficients (scg.f)
  !i   cy    :Normalization constants for spherical harmonics
  !i   slat  :struct containing information about the lattice
  !o Outputs
  !o   s     :integrals; see Remarks
  !r Remarks
  !r   s(L,M) contains integral of H_L^*(r-p1) (laplace)^k H_M(r-p2)
  !r   Row L corresponds to p1 and col M corresponds to p2.
  !r   Strux s(L,M) is computed for dr=p1-p2
  !r   See JMP 39, 3393, Section 10
  !u Updates
  !u   8 Jun 00  Added 10s digit mode.  New arg list.
  !u   19 May 00 Made rsm1,e1,rsm2,e2 l-dependent.  Elements for which
  !u             rsm1 =0 or rsm2 = 0 are not computed.
  !u   18 May 00 Adapted from nfp hhi_bl.f
  ! ----------------------------------------------------------------------
  implicit none
  integer :: mode,jcg(1),indxcg(1),nlm1,nlm2,kmax,ndim1,ndim2
  double precision :: p1(3),p2(3),q(3),cg(1),cy(1),rsm1(0:*),rsm2(0:*),e1(0:*),e2(0:*)
  double complex s(ndim1,ndim2,0:kmax)
  integer :: m,lmx1,lmx2,ll,l1,l2,k,jlm,ilm,lm11,lm21,lm12,lm22,l1t,l2t
  double precision :: dr(3)
  if (nlm1 == 0 .OR. nlm2 == 0) return
  do  1  m = 1, 3
     dr(m) = p1(m)-p2(m)
1 enddo
  lmx1 = ll(nlm1)
  lmx2 = ll(nlm2)
  do    k = 0, kmax
     do    jlm = 1, nlm2
        do    ilm = 1, nlm1
           s(ilm,jlm,k) = dcmplx(0d0,0d0)
        enddo
     enddo
  enddo
  l1t = -1
  do  20  l1 = 0, lmx1
     if (l1 <= l1t) goto 20
     call gtbsl2(l1,lmx1,e1,rsm1,l1t)
     !       l1t = l1
     l2t = -1
     do  22  l2 = 0, lmx2
        if (l2 <= l2t) goto 22
        call gtbsl2(l2,lmx2,e2,rsm2,l2t)
        !         l2t = l2
        lm11 = l1**2+1
        lm12 = (l1t+1)**2
        lm21 = l2**2+1
        lm22 = (l2t+1)**2
        if (mode/10 == 1 .AND. rsm1(l1)*rsm2(l2) == 0) goto 22
        call phhibl(dr,q,rsm1(l1),rsm2(l2),e1(l1),e2(l2),lm11,lm12, &
             lm21,lm22,kmax,ndim1,ndim2,cg,indxcg,jcg,cy,s) !,slat
22   enddo
20 enddo
end subroutine hhibl

subroutine phhibl(dr,q,rsm1,rsm2,e1,e2,mlm1,nlm1,mlm2,nlm2,kmax, &
     ndim1,ndim2,cg,indxcg,jcg,cy,s) !,slat
  !- Integrals between smooth Hankels with k-th power of Laplace operator.
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   dr    :p1-p2
  !i   q     :wave number for Bloch sum
  !i   rsm1  :smoothing radius of Hankels at p1
  !i   rsm2  :smoothing radius of Hankels at p2
  !i   e1    :energy of smooth Hankels at p1
  !i   e2    :energy of smooth Hankels at p2
  !i   nlm1  :L-cutoff for functions at p1
  !i   nlm2  :L-cutoff for functions at p2
  !i   kmax  :cutoff in power of Laplace operator
  !i   ndim1 :leading dimension of s
  !i   ndim2 :second dimension of s
  !i   cg    :Clebsch Gordon coefficients (scg.f)
  !i   indxcg:index for Clebsch Gordon coefficients
  !i   jcg   :L q.n. for the C.G. coefficients (scg.f)
  !i   cy    :Normalization constants for spherical harmonics
  !i   slat  :struct containing information about the lattice
  !o Outputs
  !o   s     :integrals; see Remarks
  !r Remarks
  !r   s(L,M) contains integral of H_L^*(r-p1) (laplace)^k H_M(r-p2)
  !u Updates
  !u   18 May 00 Adapted from nfp hhi_bl.f
  ! ----------------------------------------------------------------------
  implicit none
  integer :: jcg(1),indxcg(1),mlm1,nlm1,mlm2,nlm2,kmax,ndim1,ndim2
  double precision :: dr(3),q(3),cg(1),cy(1),rsm1, rsm2,e1,e2
  double complex s(ndim1,ndim2,0:kmax)
  ! ... Local parameters
  integer :: nlm0,ktop0,icg,icg1,icg2,ii,ilm,ilm1,ilm2,indx,ip,k, &
       ktop,l1,l2,ll,lm,lmax1,lmax2,lmaxx,nlmx
  parameter( nlm0=100, ktop0=10 )
  double precision :: fpi,e,fac,fac1,fac2,gam1,gam2,gamx,rsmx
  double complex hkl1(0:ktop0,nlm0),hkl2(0:ktop0,nlm0), hsm(nlm0),hsmp(nlm0)
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
  if (nlmx > nlm0) call rxi('increase nlm0 in hhibl need',nlmx)
  if (ktop > ktop0) call rx('hhibl: increase ktop0')
  ! ... Set up functions for connecting vector p1-p2
  if (dabs(e1-e2) > 1d-5) then
     fac1 = dexp(gam2*(e2-e1))/(e1-e2)
     fac2 = dexp(gam1*(e1-e2))/(e2-e1)
     call hklbl(dr,rsmx,e1,q,ktop,nlmx,ktop0,cy, hkl1) !slat,
     call hklbl(dr,rsmx,e2,q,ktop,nlmx,ktop0,cy, hkl2) !slat,
     do  4  ilm = 1, nlmx
        do  5  k = 0, ktop
           hkl1(k,ilm) = fac1*hkl1(k,ilm) + fac2*hkl2(k,ilm)
5       enddo
4    enddo
  else
     e = .5d0*(e1+e2)
     call hklbl(dr,rsmx,e,q,ktop,nlmx,ktop0,cy, hkl2) !slat,
     call hsmbl(dr,rsmx,e,q,lmaxx,cy, hsm,hsmp) !slat,
     do  6  ilm = 1, nlmx
        hkl1(0,ilm) = hsmp(ilm) - gamx*hsm(ilm)
        do  7  k = 1, ktop
           hkl1(k,ilm) = -e*hkl1(k-1,ilm) - hkl2(k-1,ilm)
7       enddo
6    enddo
  endif
  ! ... Combine with Clebsch-Gordan coefficients
  do  1111  ilm1 = mlm1, nlm1
     l1 = ll(ilm1)
     do  111  ilm2 = mlm2, nlm2
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
12         enddo
11      enddo
111  enddo
1111 enddo
end subroutine phhibl

subroutine gklbl(p,rsm,e,q,kmax,nlm,k0,cy,gkl) !,slat
  use m_lattic,only: rv_a_oqlv,rv_a_odlv,plat=>lat_plat
  use m_lmfinit,only:alat=> lat_alat
  use m_lattic,only: qlat=>lat_qlat, vol=>lat_vol,awald=>lat_awald, lat_nkd, lat_nkq
  use m_shortn3_plat,only: shortn3_plat,nout,nlatout
  !- Bloch-sums of k,L-dependent gaussians
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   p     :Function is centered at p
  !i   rsm   :smoothing radius
  !i   e     :gkL scaled by exp(e*rsm**2/4)
  !i   q     :wave number for Bloch sum (units of 2*pi/alat)
  !i   kmax  :polynomial cutoff
  !i   nlm   :L-cutoff for gkl
  !i   k0    :leading dimension of gkl
  !i   cy    :Normalization constants for spherical harmonics
  !i   slat  :struct containing information about the lattice
  !o Outputs
  !o   gkl   :Bloch-summed Gaussians
  !u Updates
  !u   24 Apr 00 Adapted from nfp gkl_bl.f
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: k0,kmax,nlm
  real(8):: e , rsm , q(3) , p(3) , cy(1)
  double complex gkl(0:k0,nlm)
  ! ... Local parameters
  double complex cfac,phase
  integer:: ilm , k , ll , lmax , nkd , nkq , owk , l , m
  double precision :: pi,sp,p1(3),rwald,ppin(3)
  real(8),allocatable:: wk(:)
  if (nlm == 0) return
  pi = 4d0*datan(1d0)
  nkd=lat_nkd
  nkq=lat_nkq
  ! ... Shorten p, multiply by corresponding phase later
  !call shorbz(p,p1,plat,qlat)
  ppin=matmul(transpose(qlat),p) 
  call shortn3_plat(ppin) 
  p1= matmul(plat,ppin+nlatout(:,1))

  sp = 2*pi*(q(1)*(p(1)-p1(1))+q(2)*(p(2)-p1(2))+q(3)*(p(3)-p1(3)))
  phase = dcmplx(dcos(sp),dsin(sp))
  lmax = ll(nlm)
  rwald = 1d0/awald
  ! ... If the smoothing radius is larger than the Ewald parameter,
  !     make the summation in reciprocal space, else in real space
  if (rsm >= rwald) then
     call gklblq ( p1 , rsm , q , kmax , nlm , k0 , alat , rv_a_oqlv , nkq , vol , gkl )
  else
     allocate(wk((kmax+1)*(lmax+1)))
     call gklbld ( p1 , rsm , q , kmax , nlm , k0 , alat , rv_a_odlv , nkd , wk , gkl )
     deallocate(wk)
  endif
  cfac = (0d0,1d0)*dexp(0.25d0*e*rsm*rsm)*phase
  ilm = 0
  do    l = 0, lmax
     cfac = cfac*(0d0,-1d0)
     do    m = 1, 2*l+1
        ilm = ilm+1
        do    k = 0, kmax
           gkl(k,ilm) = cfac*cy(ilm)*gkl(k,ilm)
        enddo
     enddo
  enddo
end subroutine gklbl
! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
subroutine gklbld(p,rsm,q,kmax,nlm,k0,alat,dlv,nkd,wk,gkl)
  !- Evaluate gkl in real space
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   p     :Function is centered at p
  !i   rsm   :smoothing radius
  !i   q     :wave number for Bloch sum
  !i   kmax  :polynomial cutoff
  !i   nlm   :L-cutoff for gkl
  !i   k0    :leading dimension of gkl
  !i   alat  :length scale of lattice and basis vectors, a.u.
  !i   dlv   :direct-space lattice vectors
  !i   nkd   :number of dlv
  !i   wk    :work array of dimension (kmax+1)(lmax+1), lmax=ll(nlm)
  !o Outputs
  !o   gkl   :Bloch-summed Gaussians
  !r Remarks
  !u Updates
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: k0,kmax,nkd,nlm
  double precision :: p(3),q(3),alat,rsm,wk(0:kmax,0:*),dlv(3,nkd)
  double complex gkl(0:k0,nlm)
  ! ... Local parameters
  integer :: nlm0,ilm,ir,k,l,ll,lmax,m,nm
  parameter (nlm0=144)
  double precision :: yl(nlm0),r(3),qdotr,r1,r2,tpi
  double complex cfac
  if (nlm > nlm0) call rxi('increase nlm0 in gklbld; need',nlm)
  tpi = 8d0*datan(1d0)
  lmax = ll(nlm)
  do    ilm = 1, nlm
     do    k = 0, kmax
        gkl(k,ilm) = 0d0
     enddo
  enddo
  do  20  ir = 1, nkd
     r(1) = alat*(p(1)-dlv(1,ir))
     r(2) = alat*(p(2)-dlv(2,ir))
     r(3) = alat*(p(3)-dlv(3,ir))
     call sylm(r,yl,lmax,r2)
     r1 = dsqrt(r2)
     call radgkl(r1,rsm,kmax,lmax,kmax,wk)
     qdotr = tpi*(q(1)*dlv(1,ir)+q(2)*dlv(2,ir)+q(3)*dlv(3,ir))
     cfac = cdexp(dcmplx(0d0,qdotr))
     ilm = 0
     do  38  l = 0, lmax
        nm = 2*l+1
        do  39  m = 1, nm
           ilm = ilm+1
           do  40  k = 0, kmax
              gkl(k,ilm) = gkl(k,ilm) + yl(ilm)*cfac*wk(k,l)
40         enddo
39      enddo
        cfac = cfac*(0d0,1d0)
38   enddo
20 enddo
end subroutine gklbld
! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
subroutine gklblq(p,rsm,q,kmax,nlm,k0,alat,qlv,nkq,vol,gkl)
  !- Evaluate gkl in reciprocal space
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   p     :Function is centered at p
  !i   rsm   :smoothing radius
  !i   q     :wave number for Bloch sum
  !i   kmax  :polynomial cutoff
  !i   nlm   :L-cutoff for gkl
  !i   k0    :leading dimension of gkl
  !i   alat  :length scale of lattice and basis vectors, a.u.
  !i   dlv   :reciprocal-space lattice vectors
  !i   nkq   :number of qlv
  !i   vol   :cell volume
  !o Outputs
  !o   gkl   :Bloch-summed Gaussians
  !r Remarks
  !u Updates
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: k0,kmax,nkq,nlm
  double precision :: alat,rsm,vol,q(3),p(3),qlv(3,nkq)
  double complex gkl(0:k0,nlm)
  ! ... Local parameters
  integer :: ilm,ir,k,ll,lmax,nlm0
  parameter (nlm0=144)
  double precision :: a,gamma,r2,scalp,tpi,tpiba,vfac,r(3),yl(nlm0)
  double complex eiphi,add,add0
  if (nlm > nlm0) call rxi('increase nlm0 in gklblq; need',nlm)
  tpi = 8d0*datan(1d0)
  a = 1d0/rsm
  gamma = .25d0/(a*a)
  vfac = 1d0/vol
  tpiba = tpi/alat
  lmax = ll(nlm)
  do ilm = 1, nlm
     do k = 0, kmax
        gkl(k,ilm) = 0d0
     enddo
  enddo
  do ir = 1,nkq
     r(1) = tpiba*(q(1)+qlv(1,ir))
     r(2) = tpiba*(q(2)+qlv(2,ir))
     r(3) = tpiba*(q(3)+qlv(3,ir))
     scalp = alat*(r(1)*p(1)+r(2)*p(2)+r(3)*p(3))
     eiphi = dcmplx(dcos(scalp),dsin(scalp))
     call sylm(r,yl,lmax,r2)

     add0 = vfac*dexp(-gamma*r2)
     do ilm = 1, nlm
        add = add0
        do k = 0, kmax
           gkl(k,ilm) = gkl(k,ilm) + yl(ilm)*eiphi*add
           add = -r2*add
        enddo
     enddo

  enddo
end subroutine gklblq

subroutine hklbl(p,rsm,e,q,kmax,nlm,k0,cy,hkl)
  use m_lmfinit,only: lat_alat
  use m_lattic,only: rv_a_oqlv,rv_a_odlv,lat_plat
  use m_lattic,only: lat_qlat, lat_vol, lat_awald,lat_nkd,lat_nkq
  use m_shortn3_plat,only: shortn3_plat,nout,nlatout
  !- Bloch-sums of k,L-dependent smooth Hankel functions.
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   p     :Function is centered at p
  !i   rsm   :smoothing radius
  !i   e     :energy of smoothed Hankel
  !i   q     :wave number for Bloch sum
  !i   kmax  :polynomial cutoff
  !i   nlm   :L-cutoff for hkl
  !i   k0    :leading dimension of hkl
  !i   cy    :Normalization constants for spherical harmonics
  !o Outputs
  !o   hkl   :Bloch-summed smoothed Hankels
  !r Remarks
  !r   H_kL = laplace^k H_L
  !r   Uses the recursion relation H_k+1,L = -e*H_kL - 4*pi*G_kL
  !u Updates
  !u   24 Apr 00 Adapted from nfp hkl_bl.f
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: k0,kmax,nlm
  real(8):: e , rsm , q(3) , p(3) , cy(1)
  double complex hkl(0:k0,nlm)
  ! ... Local parameters
  integer:: nlm0 , ilm , job , k , ll , lmax , nkd , nkq , nrx &
       , owk , oyl
  parameter (nlm0=144)
  double precision :: alat,awald,fpi,pi,sp,vol,plat(3,3), &
       qlat(3,3),p1(3),ppin(3)
  double complex hsm(nlm0),hsmp(nlm0),phase,gklsav,gklnew
  real(8),allocatable:: wk(:),yl(:)
  double precision :: faca
  parameter (faca=1d0)

  if (nlm == 0) return
  if (nlm > nlm0) call rxi('increase nlm0 in hklbl need',nlm)
  pi = 4d0*datan(1d0)
  fpi = 4d0*pi
  lmax = ll(nlm)
  ! ... Shorten p, multiply by corresponding phase later
  alat=lat_alat
  plat=lat_plat
  qlat=lat_qlat
!  call shorbz(p,p1,plat,qlat)
  ppin=matmul(transpose(qlat),p) 
  call shortn3_plat(ppin) 
  p1= matmul(plat,ppin+nlatout(:,1))
  
  sp = 2*pi*(q(1)*(p(1)-p1(1))+q(2)*(p(2)-p1(2))+q(3)*(p(3)-p1(3)))
  phase = dcmplx(dcos(sp),dsin(sp))
  awald=lat_awald
!  tol=lat_tol
  vol=lat_vol
  nkd=lat_nkd
  nkq=lat_nkq
  nrx = max(nkd,nkq)
  allocate(wk(nrx*(2*lmax+10)),yl(nrx*(lmax+1)**2))
  call hsmq ( 1 , 0 , ll ( nlm ) , e , rsm , 0000 , q , p1 , nrx &
       , nlm , wk , yl , awald , alat , rv_a_oqlv , nkq , rv_a_odlv &
       , nkd , vol , hsm , hsmp )
  if (rsm > faca/awald) then
     call gklbl(p1,rsm,e,q,kmax-1,nlm,k0,cy, hkl) !slat,
  else
     job = 2
     call gklq ( lmax , rsm , q , p1 , e , kmax - 1 , k0 , alat , &
          rv_a_odlv , nkd , nrx , yl , wk , job , hkl )
  endif
  deallocate(wk,yl)
  ! --- H_kL by recursion ---
  do   ilm = 1, nlm
     gklsav = hkl(0,ilm)
     hkl(0,ilm) = hsm(ilm)
     do    k = 1, kmax
        gklnew = hkl(k,ilm)
        hkl(k,ilm) = -e*hkl(k-1,ilm) - fpi*gklsav
        gklsav = gklnew
     enddo
  enddo
  ! ... Put in phase to undo shortening
  if (sp /= 0) then
     do  20  ilm = 1,nlm
        do  22  k = 0,kmax
           hkl(k,ilm) = phase*hkl(k,ilm)
22      enddo
20   enddo
  endif
end subroutine hklbl

subroutine fklbl(p,rsm,kmax,nlm,k0,cy,fkl)
  use m_lattic,only: rv_a_odlv,rv_a_oqlv,plat=>lat_plat
  use m_lmfinit,only: alat=>lat_alat,tol=>lat_tol
  use m_lattic,only: qlat=>lat_qlat, vol=>lat_vol, awald=>lat_awald, nkd=>lat_nkd, nkq=>lat_nkq
  use m_shortn3_plat,only: shortn3_plat,nout,nlatout
  !- Bloch sum of smooth Hankels for e=0 and q=(0,0,0).
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   p     :Function is centered at p
  !i   rsm   :smoothing radius
  !i   kmax  :polynomial cutoff
  !i   nlm   :L-cutoff for gkl
  !i   k0    :leading dimension of gkl
  !i   cy    :Normalization constants for spherical harmonics
  !i   slat  :struct containing information about the lattice
  !o Outputs
  !o   fkl   :Bloch-summed Hankels for q=0 and e=0
  !r Remarks
  !r   For (k=0,l=0) f equals the limit of hklbl minus the avg value.
  !r   For all other cases f is the limit of hklbl as e goes to zero.
  !u Updates
  !u   24 Apr 00 Adapted from nfp fkl_bl.f
  ! ----------------------------------------------------------------------
  implicit none
  integer :: kmax,nlm,k0
  real(8):: p(3) , cy(*) , rsm,ppin(3)
  integer:: nlm0 , lmax , ll, nrx , owk , oyl , job, ilm , k
  parameter ( nlm0=196 )
  double precision :: q(3),p1(3),faca,fpi,y0,e
  complex(8):: fsm(nlm0),gklsav,gklnew,fkl(0:k0,nlm)
  parameter (faca=1d0)
  real(8),allocatable:: wk(:),yl(:)
  if (nlm == 0) return
  fpi = 16d0*datan(1d0)
  y0 = 1d0/dsqrt(fpi)
  if (nlm > nlm0) call rx('fklbl: increase nlm0')
  lmax = ll(nlm)
  e = 0d0
  q = 0d0
  !call shorbz(p,p1,plat,qlat)
  ppin=matmul(transpose(qlat),p) 
  call shortn3_plat(ppin) 
  p1= matmul(plat,ppin+nlatout(:,1))

  nrx = max(nkd,nkq)
  allocate(wk(nrx*(2*lmax+10)),yl(nrx*(lmax+1)**2))
  call hsmqe0 ( lmax , rsm , 0 , q , p1 , nrx , nlm , wk , yl , &
       awald , alat , rv_a_oqlv , nkq , rv_a_odlv , nkd , vol , fsm  )
  if (rsm > faca/awald) then
     call gklbl(p1,rsm,e,q,kmax-1,nlm,k0,cy, fkl) !slat,
  else
     job = 2
     call gklq ( lmax , rsm , q , p1 , e , kmax - 1 , k0 , alat , &
          rv_a_odlv , nkd , nrx , yl , wk , job , fkl )
  endif
  deallocate(wk,yl)
  ! ... Upward recursion in k: mainly sets fkl = -4*pi * g(k-1,l)
  do    ilm = 1, nlm
     gklsav = fkl(0,ilm)
     fkl(0,ilm) = fsm(ilm)
     do    k = 1, kmax
        gklnew = fkl(k,ilm)
        fkl(k,ilm) = -fpi*gklsav
        gklsav = gklnew
     enddo
  enddo
  fkl(1,1) = fkl(1,1) + fpi*y0/vol !! ... Add extra term to F(k=1,l=0)
end subroutine fklbl

subroutine gfigbl(pg,ph,rsmg,rsmh,nlmg,nlmh,kmax,ndim1,ndim2,kdim, &
     cg,indxcg,jcg,cy,s,ds)!,slat
  !- Integrals between smooth hankels(eh=0,q=0) and gaussians
  !  with some power of the laplace operator, and their gradients.
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   pg    :center of gaussians
  !i   ph    :center of smooth Hankels
  !i   rsmg  :smoothing radius of gaussians
  !i   rsmh  :smoothing radius of Hankels
  !i   nlmg  :L-max for gaussians
  !i   nlmh  :L-max for Hankels
  !i   kmax  :cutoff in power of Laplace operator
  !i   ndim1 :leading dimensions of s,ds
  !i   ndim2 :second dimensions of s,ds
  !i   kdim  :dimensions s,ds
  !i   cg    :Clebsch Gordon coefficients, stored in condensed form (scg.f)
  !i   indxcg:index for Clebsch Gordon coefficients
  !i   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
  !i   cy    :Normalization constants for spherical harmonics
  !i   vol   :cell volume
  !i   slat  :struct containing information about the lattice
  !o Outputs
  !o   s     :integrals between Hankels and gaussians; see Remarks
  !o   ds    :gradient of s; see Remarks
  !r Remarks
  !r   s(L,M,k) contains integral of G_L^*(r-pg) (laplace)^k F_M(r-ph)
  !r            s,ds are generated for L=1..nlmg and M=1..nlmh
  !r   ds is gradient of s wrt pg; use -ds for grad wrt ph.
  !u Updates
  !u   28 Apr 97 adapted to handle case q=0 and e1 or e2 zero.
  !u   22 Apr 00 Adapted from nfp hhig_bl.f
  ! ----------------------------------------------------------------------
  implicit none
  integer :: jcg(*),indxcg(*),nlmg,nlmh,kmax,ndim1,ndim2,kdim
  double precision :: ph(3),pg(3),cg(1),cy(1),rsmg,rsmh !,slat(1)
  double complex s(ndim1,ndim2,0:kdim),ds(ndim1,ndim2,0:kdim,3)
  integer :: nlm0,ktop0,m,lmaxh,ll,lmaxg,lmaxx,nlmx,nlmxp1,ktop,ktopp1, &
       ilm,kz,kx1,kx2,ky1,ky2,k,jlm,ilg,lg,ilh,lh,ii,indx,icg1,icg2, &
       icg,lm,ip,nlmxx
  double precision :: dr(3),gamh,gamg,rsmx,cz,cx1,cx2,cy1,cy2,fac
  complex(8),allocatable:: hkl(:,:),ghkl(:,:,:)
  if (nlmh == 0 .OR. nlmg == 0) return
  gamh = 0.25d0*rsmh*rsmh
  gamg = 0.25d0*rsmg*rsmg
  rsmx = 2d0*dsqrt(gamg+gamh)
  do  1  m = 1, 3
     dr(m) = pg(m)-ph(m)
1 enddo
  lmaxh = ll(nlmh)
  lmaxg = ll(nlmg)
  lmaxx = lmaxg+lmaxh
  nlmx = (lmaxx+1)**2
  nlmxp1 = (lmaxx+2)**2
  ktop = max0(lmaxg,lmaxh) + kmax
  ktopp1 = ktop+1
  nlm0 = nlmxp1
  ktop0 = ktopp1
  allocate(hkl(0:ktop0,nlm0),ghkl(0:ktop0,nlm0,3))
  if (nlmxp1 > nlm0)  call rxi('gfigbl: need nlm0 ge',nlmxp1)
  if (ktopp1 > ktop0) call rxi('gfigbl: need ktop0 ge',ktopp1)
  ! ... Make functions for connecting vector
  call fklbl(dr,rsmx,ktopp1,nlmxp1,ktop0,cy, hkl) 
  ! ... Make gradients using CGs for p functions
  call dpzero(ghkl, 2*(ktop0+1)*nlm0*3)
  nlmxx = lmaxx*lmaxx
  do  ilm = 1, nlmx
     call scglp1(ilm,kz,cz,kx1,kx2,cx1,cx2,ky1,ky2,cy1,cy2)
     do  k = 0, ktop
        ghkl(k,ilm,1) = ghkl(k,ilm,1)-cx1*hkl(k,kx1)-cx2*hkl(k,kx2)
        ghkl(k,ilm,2) = ghkl(k,ilm,2)-cy1*hkl(k,ky1)-cy2*hkl(k,ky2)
        ghkl(k,ilm,3) = ghkl(k,ilm,3)-cz*hkl(k,kz)
        if (ilm <= nlmxx) then
           ghkl(k,kx1,1) = ghkl(k,kx1,1) - cx1*hkl(k+1,ilm)
           ghkl(k,kx2,1) = ghkl(k,kx2,1) - cx2*hkl(k+1,ilm)
           ghkl(k,ky1,2) = ghkl(k,ky1,2) - cy1*hkl(k+1,ilm)
           ghkl(k,ky2,2) = ghkl(k,ky2,2) - cy2*hkl(k+1,ilm)
           ghkl(k,kz,3)  = ghkl(k,kz,3)  - cz *hkl(k+1,ilm)
        endif
     enddo
  enddo

  do    k = 0, kmax
     do    jlm = 1, nlmh
        do    ilm = 1, nlmg
           s(ilm,jlm,k)    = dcmplx(0d0,0d0)
           ds(ilm,jlm,k,1) = dcmplx(0d0,0d0)
           ds(ilm,jlm,k,2) = dcmplx(0d0,0d0)
           ds(ilm,jlm,k,3) = dcmplx(0d0,0d0)
        enddo
     enddo
  enddo
  ! --- Combine with Clebsch-Gordan coefficients ---
  do  1111  ilg = 1, nlmg
     lg = ll(ilg)
     do  111  ilh = 1, nlmh
        lh = ll(ilh)
        ii = max0(ilg,ilh)
        indx = (ii*(ii-1))/2 + min0(ilg,ilh)
        icg1 = indxcg(indx)
        icg2 = indxcg(indx+1)-1
        do  11  icg = icg1, icg2
           ilm = jcg(icg)
           lm = ll(ilm)
           k = (lg+lh-lm)/2
           fac = (-1d0)**lg*cg(icg)
           do  12  ip = 0, kmax
              s(ilg,ilh,ip) = s(ilg,ilh,ip) + fac*hkl(k+ip,ilm)
              ds(ilg,ilh,ip,1) = ds(ilg,ilh,ip,1)+fac*ghkl(k+ip,ilm,1)
              ds(ilg,ilh,ip,2) = ds(ilg,ilh,ip,2)+fac*ghkl(k+ip,ilm,2)
              ds(ilg,ilh,ip,3) = ds(ilg,ilh,ip,3)+fac*ghkl(k+ip,ilm,3)
12         enddo
11      enddo
111  enddo
1111 enddo
  deallocate(hkl,ghkl)
end subroutine gfigbl

subroutine ggugbl(p1,p2,rsm1,rsm2,nlm1,nlm2,ndim1,ndim2,s,ds) !slat,
  use m_lmfinit,only: rv_a_ocy,rv_a_ocg, iv_a_oidxcg, iv_a_ojcg
  !- Estatic energy integrals between Bloch gaussians, and gradients.
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   slat  :struct containing information about the lattice
  !i   p1    :first center
  !i   p2    :second center
  !i   rsm1  :smoothing radius of Gaussians at p1
  !i   rsm2  :smoothing radius of Gaussians at p2
  !i   e1    :energy  of Gaussians at p1
  !i   e2    :energy  of Gaussians at p2
  !i   nlm1  :L-max for  Gaussians at p1
  !i   nlm2  :L-max for  Gaussians at p2
  !i   ndim1 :leading dimensions of s,ds
  !i   ndim2 :second dimensions of s,ds
  !i   slat  :struct containing information about the lattice
  !o Outputs
  !o   s     :integrals between Bloch Gaussians
  !o   ds    :gradient of s; see Remarks
  !r Remarks
  !r   Gradient is wrt p1; use -ds for grad wrt p2.
  !u Updates
  !u   22 Apr 00 Adapted from nfp ggug_bl.f
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: nlm1,nlm2,ndim1,ndim2
  real(8):: rsm1 , rsm2 , p1(3) , p2(3)
  double complex s(ndim1,ndim2),ds(ndim1,ndim2,3)
  !! ... Local parameters
  integer:: kmax , kdim , ilm2 , ilm1
  kmax = 0
  kdim = 0
  call gfigbl ( p1 , p2 , rsm1 , rsm2 , nlm1 , nlm2 , kmax , ndim1 &
       , ndim2 , kdim , rv_a_ocg , iv_a_oidxcg , iv_a_ojcg , rv_a_ocy &
       ,  s , ds ) !slat ,
  do  ilm2 = 1, nlm2
     do  ilm1 = 1, nlm1
        s(ilm1,ilm2)    = 2d0*s(ilm1,ilm2)
        ds(ilm1,ilm2,1) = 2d0*ds(ilm1,ilm2,1)
        ds(ilm1,ilm2,2) = 2d0*ds(ilm1,ilm2,2)
        ds(ilm1,ilm2,3) = 2d0*ds(ilm1,ilm2,3)
     enddo
  enddo
end subroutine ggugbl

subroutine hklgbl(p,rsm,e,q,kmax,nlm,k0,nlm0,cy,hkl,ghkl)!slat,
  !- Bloch-sums of k,L-dependent smooth hankel functions and gradients
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   p     :Function is centered at p
  !i   rsm   :smoothing radius
  !i   e     :energy of smoothed Hankel
  !i   q     :wave number for Bloch sum
  !i   kmax  :polynomial cutoff
  !i   nlm   :L-cutoff for hkl
  !i   k0    :leading dimension of hkl
  !i   cy    :Normalization constants for spherical harmonics
  !i   slat  :struct containing information about the lattice
  !o Outputs
  !o   hkl   :Bloch-summed smoothed Hankels
  !o   ghkl  :gradients of hkl
  !r Remarks
  !r   H_kL = laplace^k H_L
  !r   Uses the recursion relation H_k+1,L = -e*H_kL - 4*pi*G_kL
  !r   H_kL are made to kmax+1, lmax+1 in order to assemble grads.
  !u Updates
  !u   25 May 00 Adapted from nfp hklg_bl.f
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: k0,kmax,nlm,nlm0
  double precision :: q(3),p(3),cy(1),e,rsm!slat(1),
  double complex hkl(0:k0,nlm0),ghkl(0:k0,nlm0,3)
  ! ... Local parameters
  integer :: ilm,k,kx1,kx2,ky1,ky2,kz,ll,lmax,m,nlm1
  double precision :: cx1,cx2,cy1,cy2,cz

  if (nlm == 0) return
  lmax = ll(nlm)
  nlm1 = (lmax+2)**2
  if (nlm1 > nlm0) call rxi('hklgbl: need nlm0 ge ',nlm1)
  if (kmax+1 > k0) call rxi('hklgbl: need k0 ge ',kmax+1)
  do  m = 1, 3
     do  ilm = 1, nlm
        do  k = 0, kmax
           ghkl(k,ilm,m) = (0d0,0d0)
        enddo
     enddo
  enddo

  ! ... Make Hkl's up to one higher in l and k
  call hklbl(p,rsm,e,q,kmax+1,nlm1,k0,cy,hkl)!,slat

  ! ... Assemble gradients using Clebsh-Gordans for p functions
  do  ilm = 1, nlm
     call scglp1(ilm,kz,cz,kx1,kx2,cx1,cx2,ky1,ky2,cy1,cy2)
     do  k = 0, kmax
        ghkl(k,ilm,1) = ghkl(k,ilm,1)-cx1*hkl(k,kx1)-cx2*hkl(k,kx2)
        ghkl(k,ilm,2) = ghkl(k,ilm,2)-cy1*hkl(k,ky1)-cy2*hkl(k,ky2)
        ghkl(k,ilm,3) = ghkl(k,ilm,3)-cz*hkl(k,kz)
        if (ilm <= lmax*lmax) then
           ghkl(k,kx1,1) = ghkl(k,kx1,1) - cx1*hkl(k+1,ilm)
           ghkl(k,kx2,1) = ghkl(k,kx2,1) - cx2*hkl(k+1,ilm)
           ghkl(k,ky1,2) = ghkl(k,ky1,2) - cy1*hkl(k+1,ilm)
           ghkl(k,ky2,2) = ghkl(k,ky2,2) - cy2*hkl(k+1,ilm)
           ghkl(k,kz,3) =  ghkl(k,kz,3)  - cz *hkl(k+1,ilm)
        endif
     enddo
  enddo
end subroutine hklgbl
  
subroutine ghibl(pg,ph,q,rsmg,rsmh,eg,eh,nlmg,nlmh,kmax, &
     ndim1,ndim2,cg,indxcg,jcg,cy,s)!,slat
  !- Block of integrals between smooth hankels and gaussians with some power
  !  of the laplace operator.
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   ph    :Function is centered at ph; see Remarks
  !i   pg    :Function it expansed is at pg; see Remarks
  !i   q     :wave number for Bloch sum
  !i   rsmg  :smoothing radius of gaussian
  !i   rsmh  :vector of l-dependent smoothing radii of smoothed hankel
  !i         :rsmh must be specified for 1..ll(nlmh)+1
  !i   eg    :gkL scaled by exp(e*rsm**2/4)
  !i   eh    :vector of l-dependent energies of smoothed Hankel
  !i         :eh must be specified for 1..ll(nlmh)+1
  !i   nlmg  :L-cutoff for P_kL expansion
  !i   nlmh  :L-cutoff for smoothed Hankel functions
  !i   kmax  :polynomial cutoff
  !i   ndim1 :leading dimension of s
  !i   ndim2 :second dimension of s
  !i   cg    :Clebsch Gordon coefficients, stored in condensed form (scg.f)
  !i   indxcg:index for Clebsch Gordon coefficients
  !i   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
  !i   cy    :Normalization constants for spherical harmonics
  !i   slat  :struct containing information about the lattice
  !o Outputs
  !o   s     :integrals of gaussian and Hankels; see Remarks
  !r Remarks
  !r   s(L,M,k) contains integral of G_L^*(r-pg) (lap)^k H_M(r-ph)
  !r   See J. Math. Phys. {\bf 39},3393 (1998), Eq. 8.4
  !u Updates
  !u   18 May 00 Made rsmh,eh l-dependent
  !u   24 Apr 00 Adapted from nfp ghi_bl.f
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: nlmg,nlmh,kmax,ndim1,ndim2,jcg(1),indxcg(1)
  double precision :: rsmg,rsmh(1),eg,eh(1)
  double precision :: ph(3),pg(3),cg(1),cy(1),q(3)!,slat(1)
  double complex s(ndim1,ndim2,0:*)
  ! ... Local parameters
  integer :: nlm0,ktop0,icg,icg1,icg2,ii,ilg,ilh,ilm,indx,ip,jlm,k, &
       ktop,lg,lh,ll,lm,lmaxg,lmaxh,lmaxx,m,nlmx,l1,l2,ilm1,ilm2
  ! FCPP#if F90
  complex(8),allocatable:: hkl(:,:)
  ! FCPP#else
  ! FCPP      parameter (nlm0=144, ktop0=21)
  ! FCPP      double complex hkl(0:ktop0,nlm0)
  ! FCPP#endif
  double precision :: ee,fac,gamg,gamh,rsmx,dr(3),e,rsm

  if (nlmh == 0 .OR. nlmg == 0) return

  ! ... rsmh- and eh- independent setup
  do  1  m = 1, 3
     dr(m) = pg(m)-ph(m)
1 enddo
  lmaxh = ll(nlmh)
  lmaxg = ll(nlmg)
  lmaxx = lmaxg+lmaxh
  nlmx = (lmaxx+1)**2
  ktop = max0(lmaxg,lmaxh)+kmax
  ktop0 = ktop
  nlm0  = nlmx
  allocate( hkl(0:ktop0,nlm0))
  do    k = 0, kmax
     do    jlm = 1, nlmh
        do    ilm = 1, nlmg
           s(ilm,jlm,k) = dcmplx(0d0,0d0)
        enddo
     enddo
  enddo

  ! --- Loop over sequences of l with a common rsm,e ---
  l2 = -1
  do  20  l1 = 0, lmaxh
     if (l1 <= l2) goto 20
     call gtbsl2(l1,lmaxh,eh,rsmh,l2)
     rsm  = rsmh(l1+1)
     e    = eh(l1+1)
     if (rsm <= 0 .OR. e > 0) goto 20
     ilm1 = l1**2+1
     ilm2 = (l2+1)**2
     lmaxx= lmaxg+l2
     nlmx = (lmaxx+1)**2
     gamh = 0.25d0*rsm*rsm
     gamg = 0.25d0*rsmg*rsmg
     rsmx = 2d0*dsqrt(gamg+gamh)
     ktop = max0(lmaxg,l2)+kmax
     call hklbl(dr,rsmx,e,q,ktop,nlmx,ktop0,cy,hkl)!,slat

     !   ... Combine with Clebsch-Gordan coefficients
     ee = dexp(gamg*(eg-e))
     do  1111  ilg = 1, nlmg
        lg = ll(ilg)
        do  111  ilh = ilm1, ilm2
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
12            enddo
11         enddo
111     enddo
1111 enddo
20 enddo
end subroutine ghibl
  
subroutine ghigbl(pg,ph,q,rsmg,rsmh,eg,eh,nlmg,nlmh,kmax, &
     ndim1,ndim2,k0,cg,indxcg,jcg,cy,s,ds)!,slat
  !- Block of integrals between smooth hankels and gaussians with some
  !  power of the laplace operator, and gradients.
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   pg    :Function it expansed is at pg; see Remarks
  !i   ph    :Function is centered at ph; see Remarks
  !i   q     :wave number for Bloch sum
  !i   rsmg  :smoothing radius of gaussian
  !i   rsmh  :vector of l-dependent smoothing radii of smoothed hankel
  !i         :rsmh must be specified for 1..ll(nlmh)+1
  !i   eg    :gkL scaled by exp(e*rsm**2/4)
  !i   eh    :vector of l-dependent energies of smoothed Hankel
  !i         :eh must be specified for 1..ll(nlmh)+1
  !i   nlmg  :L-cutoff for P_kL expansion
  !i   nlmh  :L-cutoff for smoothed Hankel functions
  !i   kmax  :polynomial cutoff
  !i   ndim1 :leading dimension of s,ds
  !i   ndim2 :second dimension of s,ds
  !i   k0    :third dimenson of s,ds
  !i   cg    :Clebsch Gordon coefficients, stored in condensed form (scg.f)
  !i   indxcg:index for Clebsch Gordon coefficients
  !i   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
  !i   cy    :Normalization constants for spherical harmonics
  !i   slat  :struct containing information about the lattice
  !o Outputs
  !o   s     :integrals of gaussian and Hankels; see Remarks
  !o   ds    :gradients of s; see Remarks
  !r Remarks
  !r   s(L,M,k) contains integral of G_L^*(r-pg) (lap)^k H_M(r-ph)
  !r  ds((L,M,k) contains grad s wrt ph; take negative for grad wrt pg.
  !u Updates
  !u   25 May 00 Made rsmh,eh l-dependent
  !u   25 May 00 Adapted from nfp ghig_bl.f
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: k0,kmax,ndim1,ndim2,nlmg,nlmh,jcg(1),indxcg(1)
  double precision :: rsmg,rsmh(1),eg,eh(1)
  double precision :: ph(3),pg(3),cg(1),cy(1),q(3)!,slat(1)
  double complex s(ndim1,ndim2,0:k0),ds(ndim1,ndim2,0:k0,3)
  ! ... Local parameters
  integer :: icg,icg1,icg2,ii,ilg,ilh,ilm,ilm1,ilm2,indx,ip,jlm,k,ktop, &
       ktop0,l1,l2,lg,lh,ll,lm,lmaxg,lmaxh,lmaxx,m,nlm0,nlmx
  double precision :: ee,fac,gamg,gamh,rsmx,dr(3),e,rsm
  complex(8),allocatable:: hkl(:,:),dhkl(:,:,:)
  if (nlmh == 0 .OR. nlmg == 0) return
  ! ... rsmh- and eh- independent setup
  do  1  m = 1, 3
     dr(m) = pg(m)-ph(m)
1 enddo
  lmaxh = ll(nlmh)
  lmaxg = ll(nlmg)
  lmaxx = lmaxg+lmaxh
  nlmx = (lmaxx+1)**2
  ktop = max0(lmaxg,lmaxh)+kmax
  ktop0 = ktop+1
  nlm0  = (lmaxx+2)**2
  allocate( hkl(0:ktop0,nlm0),dhkl(0:ktop0,nlm0,3))
  do    k = 0, kmax
     do    jlm = 1, nlmh
        do    ilm = 1, nlmg
           s(ilm,jlm,k) = dcmplx(0d0,0d0)
           ds(ilm,jlm,k,1) = dcmplx(0d0,0d0)
           ds(ilm,jlm,k,2) = dcmplx(0d0,0d0)
           ds(ilm,jlm,k,3) = dcmplx(0d0,0d0)
        enddo
     enddo
  enddo

  ! --- Loop over sequences of l with a common rsm,e ---
  l2 = -1
  do  20  l1 = 0, lmaxh
     if (l1 <= l2) goto 20
     call gtbsl2(l1,lmaxh,eh,rsmh,l2)
     rsm  = rsmh(l1+1)
     e    = eh(l1+1)
     if (rsm <= 0 .OR. e > 0) goto 20
     ilm1 = l1**2+1
     ilm2 = (l2+1)**2
     lmaxx= lmaxg+l2
     nlmx = (lmaxx+1)**2
     gamh = 0.25d0*rsm*rsm
     gamg = 0.25d0*rsmg*rsmg
     rsmx = 2d0*dsqrt(gamg+gamh)
     ktop = max0(lmaxg,l2)+kmax
     call hklgbl(dr,rsmx,e,q,ktop,nlmx,ktop0,nlm0,cy,hkl,dhkl)!,slat

     !   ... Combine with Clebsch-Gordan coefficients
     ee = dexp(gamg*(eg-e))
     do  1111  ilg = 1, nlmg
        lg = ll(ilg)
        do  111  ilh = ilm1, ilm2
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
                 s(ilg,ilh,ip)    = s(ilg,ilh,ip)    + fac*hkl(k+ip,ilm)
                 ds(ilg,ilh,ip,1) = ds(ilg,ilh,ip,1) + fac*dhkl(k+ip,ilm,1)
                 ds(ilg,ilh,ip,2) = ds(ilg,ilh,ip,2) + fac*dhkl(k+ip,ilm,2)
                 ds(ilg,ilh,ip,3) = ds(ilg,ilh,ip,3) + fac*dhkl(k+ip,ilm,3)
12            enddo
11         enddo
111     enddo
1111 enddo
20 enddo
  deallocate(hkl,dhkl)
end subroutine ghigbl

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
  deallocate(s)
end subroutine hxpbl
  
subroutine hxpgbl(ph,pg,q,rsmh,rsmg,eh,kmax,nlmh,nlmg,k0,ndimh, &
     ndimg,cg,indxcg,jcg,cy,c,dc)!,slat
  !- Coefficients to expand smooth bloch hankels and grads centered at ph
  !  into a sum of polynomials P_kL centered at pg.
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   ph    :Function is centered at ph; see Remarks
  !i   pg    :Function it expansed is at pg; see Remarks
  !i   q     :wave number for Bloch sum
  !i   rsmh  :smoothing radius of smoothed Hankel (l-dependent)
  !i   rsmg  :smoothing radius of gaussian
  !i   eh    :energies of smoothed Hankel (l-dependent)
  !i   kmax  :polynomial cutoff
  !i   nlmh  :L-cutoff for smoothed Hankel functions being expanded
  !i   nlmg  :L-cutoff for P_kL expansion
  !i   k0    :leading dimension of coefficient arrays c,dc
  !i   ndimg :second dimension of coefficient arrays c,dc
  !i   ndimh :htird dimension of coefficient arrays c,dc
  !i   cg    :Clebsch Gordon coefficients, stored in condensed form (scg.f)
  !i   indxcg:index for Clebsch Gordon coefficients
  !i   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
  !i   cy    :Normalization constants for spherical harmonics
  !i   slat  :struct containing information about the lattice
  !o Outputs
  !o   c     :cofficients c(k,M,L); see Remarks
  !o   dc    :cofficients dc(k,M,L,1..3) (gradients of c); see Remarks
  !o         :here k=0..kmax, M=1..nlmg, L=1..nlmh
  !r Remarks
  !r   Expansion is:  H_L(r-ph) = sum(kM)  c(k,M,L) * P_kM(r-pg)
  !r             grad H_L(r-ph) = sum(kM) dc(k,M,L) * P_kM(r-pg)
  !r   As rsmg -> 0, expansion turns into a Taylor series of H_L.
  !r   As rsmg increases the error is spread out over a progressively
  !r   larger radius, thus reducing the error for larger r while
  !r   reducing accuracy at small r.
  !r
  !r   Grads are wrt to ph; take negative for grad wrt to pg.
  !u Updates
  !u   25 May 00 Made rsmh,eh l-dependent
  !u   25 May 00 Adapted from nfp hxpg_bl.f
  ! ----------------------------------------------------------------------
  implicit none
  integer :: k0,kmax,ndimg,ndimh,nlmg,nlmh,jcg(1),indxcg(1)
  double precision :: eh(1),rsmg,rsmh(1),ph(3),pg(3),cg(1),cy(1),q(3)
  double complex c(0:k0,ndimg,ndimh),dc(0:k0,ndimg,ndimh,3)
  integer :: ndim1,ndim2,ktop0,ilmg,ilmh,k,l,ll,lmaxg,m,nm
  double precision :: a,dfact,eg,fac,factk,fpi
  complex(8),allocatable:: s(:,:,:),ds(:,:,:,:)
  if (nlmg == 0 .OR. nlmh == 0) return
  fpi = 16d0*datan(1d0)
  ndim1 = nlmg
  ndim2 = nlmh
  ktop0 = kmax
  allocate(s(ndim1,ndim2,0:ktop0),ds(ndim1,ndim2,0:ktop0,3))
  if (kmax > ktop0) call rxi('hxpgbl: increase ktop0, need',kmax)
  if (nlmg > ndim1) call rxi('hxpgbl: increase ndim1, need',nlmg)
  if (nlmh > ndim2) call rxi('hxpgbl: increase ndim2, need',nlmh)
  ! ... Integrals of Hankels with Gaussians
  eg = 0d0
  call ghigbl(pg,ph,q,rsmg,rsmh,eg,eh,nlmg,nlmh,kmax,ndim1,ndim2, &
       ktop0,cg,indxcg,jcg,cy,s,ds)!,slat
  ! ... Scale to get coefficients of the PkL
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
5          enddo
           factk = factk*(k+1)
3       enddo
2    enddo
     dfact = dfact*(2*l+3)
1 enddo
  deallocate(s,ds)
end subroutine hxpgbl

subroutine gklq(lmax,rsm,q,p,e,kmax,k0,alat,dlv,nkd,nrx,yl,wk,job, gkl)
  use m_ropyln,only: ropyln
  use m_lattic,only: qlat=>lat_qlat,plat=>lat_plat
  use m_shortn3_plat,only: shortn3_plat,nout,nlatout
  !- Bloch sum of k,L-dependent gaussians (vectorizes)
  ! ---------------------------------------------------------------
  !i Inputs:
  !i  lmax   :l-cutoff for gkl
  !i   rsm   :smoothing radius
  !i   q     :wave number for Bloch sum (units of 2*pi/alat)
  !i   p     :connecting vector (units of alat)
  !i   e     :G_kL scaled by exp(e*rsm**2/4)
  !i   alat  :length scale of lattice and basis vectors, a.u.
  !i  dlv,nkd:direct lattice vectors, and number
  !i   nrx   :leading dimension of wk,yl
  !i   yl    :work array, dimensioned nrx*(lmax+1)**2
  !i   wk    :work array, dimensioned at least nrx*(2*lmax+10)
  !i   nrx   :dimensions work arrays yl, and must be >= max(nkq,nkd)
  !i   job   :1s digit
  !i         :0, generate wk(1,2,3,4)
  !i          1  assume wk(1,3,4) and yl have been generated already
  !i          2  assume wk(1,2,3,4) and yl have been generated already
  !i         :10s digit
  !i         :0  use standard phase convention, do not shorten p
  !i         :1  use standard phase convention, but shorten p
  !i         :2  scale standard phase convention by exp(-i q . p)
  !i         :3  like 2, but also shorten p
  !i   k0    :leading dimension of gkl
  !o Outputs:
  !o   yl:  ylm(1..nkd,1..(lmax+1)**2) for points alat*(p-dlv(1..nkd))
  !o   wk:  (*,1) holds r**2
  !o        (*,2) holds Y0 exp(-(r/rsm)**2)
  !o        (*,3) holds cos(q.dlv)
  !o        (*,4) holds sin(q.dlv)
  !o   gkl: G_kL * exp(e*rsm**2/4) generated for (0:kmax,0:lmax)
  !e External routines required: ropyln
  !u Updates
  !u  15 Aug 00 extended to e>0; added 1000 digit job
  ! ---------------------------------------------------------------
  !     implicit none
  integer :: k0,kmax,lmax,nkd,nrx,job
  double precision :: alat,rsm,p(3),q(3),dlv(3,nkd),wk(nrx,*),yl(nrx,1)
  complex(8) :: gkl(0:k0,1),img=(0d0,1d0)
  real(8):: e
  ! Local variables
  integer :: ilm,ir,k,l,m,nlm,ik1,ik2,lc1,ls1,lc2,ls2,job0,job1
  double precision :: qdotr,pi,tpi,y0,ta2,x,y,a2,g0fac,xx1,xx2,x1,x2, &
       y2,p1(3),sp,cosp,sinp,pp(3)
  if (kmax < 0 .OR. lmax < 0 .OR. rsm == 0d0) return
  job0 = mod(job,10)
  job1 = mod(job/10,10)
  nlm = (lmax+1)**2
  pi  = 4*datan(1d0)
  tpi = 2*pi
  y0  = 1/dsqrt(4*pi)
  a2  = 1/rsm**2
  ta2 = 2*a2
  do    ilm = 1, nlm
     do    k = 0, kmax
        gkl(k,ilm) = 0d0
     enddo
  enddo
  ! ... Shorten connecting vector; need to adjust phase later
  if (job1 == 1 .OR. job1 == 3) then
     !call shortn(p,p1,dlv,nkd)
     pp=matmul(transpose(qlat),p)
     call shortn3_plat(pp)
     p1 = matmul(plat,pp+nlatout(:,1))
  else
     call dcopy(3,p,1,p1,1)
  endif
  ! --- Put ylm in yl and alat**2*(p-dlv)**2 in wk(1) ---
  if (job0 == 0) then
     do  20  ir = 1, nkd
        wk(ir,2) = alat*(p1(1)-dlv(1,ir))
        wk(ir,3) = alat*(p1(2)-dlv(2,ir))
        wk(ir,4) = alat*(p1(3)-dlv(3,ir))
20   enddo
     call ropyln(nkd,wk(1,2),wk(1,3),wk(1,4),lmax,nrx,yl,wk)
     ! ...   cos(q.dlv), sin(q.dlv) -> wk(3,4), Y0 exp(-(a dlv)**2) -> wk(2)
     do  22  ir = 1, nkd
        qdotr = 2*pi*(q(1)*dlv(1,ir)+ q(2)*dlv(2,ir)+ q(3)*dlv(3,ir))
        wk(ir,3) = dcos(qdotr)
        wk(ir,4) = dsin(qdotr)
        wk(ir,2) = y0*dexp(-wk(ir,1)*a2)
22   enddo
  elseif (job0 == 1) then
     do  24  ir = 1, nkd
        wk(ir,2) = y0*dexp(-wk(ir,1)*a2)
24   enddo
  endif
  lc1 = 5
  ls1 = 6
  lc2 = 7
  ls2 = 8
  ik1 = 9
  ik2 = 10+lmax
  ! --- Outer loop over k (in blocks of 2), and over l ---
  do  301  k = 0, kmax, 2
     do  30  l = 0, lmax
        g0fac = 1/rsm*ta2**(l+1)/pi * dexp(e*rsm*rsm/4)

        !   ... Make radial part of the G_kl(1..nkd) for k= 0, 1
        if (k == 0) then
           do  32  ir = 1, nkd
              xx1 = g0fac*wk(ir,2)
              xx2 = (ta2*wk(ir,1)-3-2*l)* ta2 * xx1
              wk(ir,ik1+l) = xx1
              wk(ir,ik2+l) = xx2
              wk(ir,lc1) = wk(ir,3)*xx1
              wk(ir,ls1) = wk(ir,4)*xx1
              wk(ir,lc2) = wk(ir,3)*xx2
              wk(ir,ls2) = wk(ir,4)*xx2
32         enddo
           !   ... Make radial part of the G_kl(1..nkd) for k, k+1 from k-1, k-2
           !       and cos(q.dlv) * G_kl and sin(q.dlv) * G_kl
        else
           x = 2*(k-1)*(2*k + 2*l-1)
           y = 4*k + 2*l-1
           x2 = 2*k*(2*(k+1) + 2*l-1)
           y2 = 4*(k+1) + 2*l-1
           do  34  ir = 1, nkd
              xx1 = ta2*((ta2*wk(ir,1)-y)*wk(ir,ik2+l) - x*ta2*wk(ir,ik1+l))
              xx2 = ta2*((ta2*wk(ir,1)-y2)*xx1         -x2*ta2*wk(ir,ik2+l))
              wk(ir,ik1+l) = xx1
              wk(ir,ik2+l) = xx2
              wk(ir,lc1) = wk(ir,3)*xx1
              wk(ir,ls1) = wk(ir,4)*xx1
              wk(ir,lc2) = wk(ir,3)*xx2
              wk(ir,ls2) = wk(ir,4)*xx2
34         enddo
        endif

        !   ... For each point, add G_kl Y_L exp(i q.dlv) into Bloch G_kL
        ilm = l*l
        if (k < kmax) then
           do  36  m = -l, l
              ilm = ilm+1
              do  38  ir = nkd, 1, -1
                 gkl(k,ilm)   = gkl(k,ilm)   + wk(ir,lc1)*yl(ir,ilm)+img*wk(ir,ls1)*yl(ir,ilm)
                 gkl(k+1,ilm) = gkl(k+1,ilm) + wk(ir,lc2)*yl(ir,ilm)+img*wk(ir,ls2)*yl(ir,ilm)
38            enddo
36         enddo
        else
           do  46  m = -l, l
              ilm = ilm+1
              do  48  ir = nkd, 1, -1
                 gkl(k,ilm) = gkl(k,ilm) + wk(ir,lc1)*yl(ir,ilm)+img*wk(ir,ls1)*yl(ir,ilm)
48            enddo
46         enddo
        endif
30   enddo
301 enddo
  ! ... Put in phase to undo shortening, or different phase convention
  sp = tpi*(q(1)*(p(1)-p1(1))+q(2)*(p(2)-p1(2))+q(3)*(p(3)-p1(3)))
  if (job1 >= 2) sp = sp-tpi*(q(1)*p1(1)+q(2)*p1(2)+q(3)*p1(3))
  if (sp /= 0d0) then
     cosp = dcos(sp)
     sinp = dsin(sp)
     do    ilm = 1, nlm
        do    k   = 0, kmax
           x1 = dreal(gkl(k,ilm))
           x2 = dimag(gkl(k,ilm))
           gkl(k,ilm) = (x1*cosp - x2*sinp)+img*( x2*cosp + x1*sinp)
        enddo
     enddo
  endif
end subroutine gklq

subroutine hsmbl(p,rsm,e,q,lmax,cy,hsm,hsmp) !,slat
  use m_lattic,only: rv_a_oqlv,rv_a_odlv,lat_plat
  use m_struc_def           !Cgetarg
  use m_lmfinit,only: lat_alat,lat_tol
  use m_lattic,only: lat_qlat
  use m_lattic,only: lat_vol
  use m_lattic,only: lat_awald
  use m_lattic,only: lat_nkd
  use m_lattic,only: lat_nkq
  use m_shortn3_plat,only: shortn3_plat,nout,nlatout
  !- Bloch-sum of smooth Hankel functions and energy derivatives
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   p     :Function is centered at p
  !i   rsm   :smoothing radius
  !i   e     :energy of smoothed Hankel
  !i   q     :wave number for Bloch sum
  !i   lmax  :l-cutoff for hsm
  !i   cy    :Normalization constants for spherical harmonics
  !i   slat  :struct containing information about the lattice
  !o Outputs
  !o   hsm   :Bloch-summed smoothed Hankels
  !o   hsmp  :Energy derivatives of hsm
  !r Remarks
  !r  Hankel functions evaluated by Ewald summation.
  !r  p in units of alat, qlv in units of 2 pi / alat
  !r  See also hsmq for a vectorized version.
  !u Updates
  !u   26 Jan 07 Works with positive energies e
  !u   1 May 00 Adapted from nfp hsmbl.f
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: lmax
  real(8):: rsm , e , q(3) , p(3) , cy(1), ppin(3)
  !      type(s_lat)::slat

  double complex hsm(1),hsmp(1)
  ! ... Local parameters
  integer:: nkd , nkq , ilm , l , m
  double precision :: alat,p1(3),plat(3,3),qlat(3,3),sp,pi,rwald,awald,asm, &
       tol,vol
  double complex cfac,phase
  alat=lat_alat
  plat=lat_plat
  qlat=lat_qlat

  !call shorbz(p,p1,plat,qlat)
  ppin=matmul(transpose(qlat),p) 
  call shortn3_plat(ppin) 
  p1= matmul(plat,ppin+nlatout(:,1))
  
  pi = 4d0*datan(1d0)
  sp = 2*pi*(q(1)*(p(1)-p1(1))+q(2)*(p(2)-p1(2))+q(3)*(p(3)-p1(3)))
  phase = dcmplx(dcos(sp),dsin(sp))
  awald=lat_awald
  tol=lat_tol
  vol=lat_vol
  nkd=lat_nkd
  nkq=lat_nkq
  rwald = 1d0/awald
  asm = 1d0/rsm
  if (rsm < rwald) then
     call hsmblq ( p1 , e , q , awald , lmax , alat , rv_a_oqlv , &
          nkq , vol , hsm , hsmp )
     call hsmbld ( p1 , rsm , e , q , awald , lmax , alat , rv_a_odlv &
          , nkd , hsm , hsmp )
  else
     call hsmblq ( p1 , e , q , asm , lmax , alat , rv_a_oqlv , nkq &
          , vol , hsm , hsmp )
  endif
  ! ... Multiply by phase to undo shortening
  cfac = dcmplx(0d0,1d0)*phase
  ilm = 0
  do    l = 0, lmax
     cfac = cfac*dcmplx(0d0,-1d0)
     do    m = 1, 2*l+1
        ilm = ilm+1
        hsm(ilm)  = cfac*cy(ilm)*hsm(ilm)
        hsmp(ilm) = cfac*cy(ilm)*hsmp(ilm)
     enddo
  enddo
end subroutine hsmbl

subroutine hsmblq(p,e,q,a,lmax,alat,qlv,nkq,vol,dl,dlp)
  !- k-space part of smooth hankel bloch sum
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   p     :Function is centered at p (units of alat)
  !i   e     :energy of smoothed Hankel
  !i   q     :wave number for Bloch sum
  !i   a     :Ewald smoothing parameter
  !i   lmax  :l-cutoff for dl,dlp
  !i   alat  :length scale of lattice and basis vectors, a.u.
  !i   qlv   :reciprocal lattice vectors, units of 2pi/alat
  !i   nkq   :number of r.l.v.
  !i   vol   :cell volume
  !o Outputs
  !i   dl    :k-summed smoothed Bloch hankels
  !i   dlp   :energy derivative of dl
  !u Updates
  ! ----------------------------------------------------------------------
  implicit none
  integer :: nkq,lmax
  double precision :: a,alat,e,vol
  double precision :: q(3),p(3),qlv(3,nkq)
  double complex dl(1),dlp(1)
  integer :: lmxx,nlm,ilm,ir
  parameter (lmxx=11)
  double precision :: r(3),tpi,gamma,fpibv,tpiba,scalp,r2,den0,den1
  double precision :: yl((lmxx+1)**2)
  double complex eiphi

  if (lmax > lmxx) call rxi('hsmblq: increase lmxx to',lmax)
  tpi = 8d0*datan(1d0)
  gamma = .25d0/(a*a)
  fpibv = 2d0*tpi/vol
  tpiba = tpi/alat
  nlm = (lmax+1)**2
  do  10  ilm = 1, nlm
     dl(ilm) = dcmplx(0d0,0d0)
     dlp(ilm) = dcmplx(0d0,0d0)
10 enddo
  do  26  ir = 1, nkq
     r(1) = tpiba*(q(1)+qlv(1,ir))
     r(2) = tpiba*(q(2)+qlv(2,ir))
     r(3) = tpiba*(q(3)+qlv(3,ir))
     scalp = alat*(r(1)*p(1)+r(2)*p(2)+r(3)*p(3))
     eiphi = dcmplx(dcos(scalp),dsin(scalp))
     call sylm(r,yl,lmax,r2)
     den0 = dexp(gamma*(e-r2))/(r2-e)
     den1 = den0/(r2-e)
     do  ilm = 1, nlm
        dl(ilm)  = dl(ilm) + yl(ilm)*eiphi*den0
        dlp(ilm) = dlp(ilm) + yl(ilm)*eiphi*den1
     enddo
26 enddo
  do  35  ilm = 1, nlm
     dl(ilm)  = fpibv*dl(ilm)
     dlp(ilm) = fpibv*dlp(ilm) + gamma*dl(ilm)
35 enddo
end subroutine hsmblq

subroutine hsmbld(p,rsm,e,q,a,lmax,alat,dlv,nkd,dl,dlp)
  !- Adds real space part of reduced structure constants (ewald).
  !u Updates
  !u   10 May 07 New protections to handle error functions underflow
  !     implicit none
  ! ... Passed parameters
  integer :: lmxx,ilm,ir,l,lmax,m,nkd,nm
  double precision :: q(3),p(3),dlv(3,nkd)
  double complex dl(1),dlp(1)
  ! ... Local parameters
  logical :: lpos,lzero
  parameter (lmxx=11)
  double precision :: a,a2,akap,alat,asm,asm2,cc,ccsm,derfc,e,emkr,gl, &
       qdotr,r1,r2,rsm,srpi,ta,ta2,tasm,tasm2,umins,uplus,tpi,kap
  double precision :: yl((lmxx+1)**2),chi1(-1:10),chi2(-1:10),r(3)
  !     double complex zchi1(-1:10),zchi2(-1:10)
  double complex cfac,zikap,expikr,zuplus,zerfc
  double precision :: srfmax,fmax
  parameter (srfmax=16d0, fmax=srfmax*srfmax)

  if (lmax > lmxx) call rxi('hsmbld: increase lmxx to',lmax)
  tpi = 8.d0*datan(1.d0)
  srpi = dsqrt(tpi/2.d0)
  lpos = e .gt. 0d0
  if (lpos) then
     !     Methfessl's 'akap' is sqrt(-e) = i kap by standard kap=sqrt(e)
     kap = dsqrt(e)
     zikap = dcmplx(0d0,1d0)*kap
  else
     akap = dsqrt(-e)
  endif
  ta = 2.d0*a
  a2 = a*a
  ta2 = 2.d0*a2
  cc = 4.d0*a2*a*dexp(e/(ta*ta))/srpi

  asm = 1d0/rsm
  tasm = 2.d0*asm
  asm2 = asm*asm
  tasm2 = 2.d0*asm2
  ccsm = 4d0*asm2*asm*dexp(e/(tasm*tasm))/srpi

  do  20  ir = 1, nkd
     r(1) = alat*(p(1)-dlv(1,ir))
     r(2) = alat*(p(2)-dlv(2,ir))
     r(3) = alat*(p(3)-dlv(3,ir))
     call sylm(r,yl,lmax,r2)
     r1 = dsqrt(r2)
     lzero = .false.

     ! --- Make the xi's from -1 to lmax ---
     if (r1 < 1d-6) then
        do  31  l = 1, lmax
           chi1(l) = 0d0
           chi2(l) = 0d0
31      enddo
        if (lpos) then
           !         do  l = 1, lmax
           !           zchi1(l) = 0d0
           !           zchi2(l) = 0d0
           !         enddo
           zuplus = zerfc(zikap/ta)
           chi1(0) = -zikap*zuplus &
                +2/dsqrt(4*datan(1d0))*a * cdexp(-(zikap/ta)**2)
           zuplus = zerfc(zikap/tasm)
           chi2(0) = -zikap*zuplus &
                +2/dsqrt(4*datan(1d0))*asm * cdexp(-(zikap/tasm)**2)
           chi1(-1) = -zerfc(-zikap/ta)/zikap
           chi2(-1) = -zerfc(-zikap/tasm)/zikap
        else
           chi1(0) = ta*dexp(e/(2d0*ta2))/srpi - akap*derfc(akap/ta)
           chi1(-1) = derfc(akap/ta)/akap
           chi2(0) = tasm*dexp(e/(2d0*tasm2))/srpi - akap*derfc(akap/tasm)
           chi2(-1) = derfc(akap/tasm)/akap
        endif
     else
        if (lpos) then
           expikr = exp(zikap*r1)
           zuplus = zerfc(zikap/ta+r1*a)*expikr
           chi1(0) = expikr/r1 - dble(zuplus)/r1
           chi1(-1) = expikr/zikap + dimag(zuplus)/kap
           !          zchi1(0) = expikr/r1 - dble(zuplus)/r1
           !          zchi1(-1) = expikr/zikap + dimag(zuplus)/kap
        else
           !         If large, these are -log(uplus),-log(umins); then chi->0
           !r        If both (akap*rsm/2 -/+ r/rsm) >> 1, we have
           !r        -log u(+/-) -> (akap*rsm/2 -/+ r/rsm)^2 -/+ akap*r
           !r                    =  (akap*rsm/2)^2 + (r/rsm)^2 >> 1
           !r         u(+/-)     -> exp[-(akap*rsm/2)^2 - (r/rsm)^2] -> 0
           !r        Also if akap*r >> 1,   chi < dexp(-akap*r1) -> 0
           emkr = dexp(-akap*r1)
           if (.5d0*akap/a+r1*a > srfmax .AND. &
                .5d0*akap/a-r1*a > srfmax .OR. akap*r1 > fmax) then
              lzero = .true.
              !            uplus = derfc(.5d0*akap/a+r1*a)/emkr
              !            umins = derfc(.5d0*akap/a-r1*a)*emkr
              !            print *, 'approx',(.5d0*akap/a+r1*a)**2 - akap*r1
              !            print *, '-log up',-dlog(uplus)
              !            print *, 'approx', (.5d0*akap/a-r1*a)**2 + akap*r1
              !            print *, '-log um', -dlog(umins)
              !            print *, r1*a
              !            stop
           else
              uplus = derfc(.5d0*akap/a+r1*a)/emkr
              umins = derfc(.5d0*akap/a-r1*a)*emkr
              chi1(0) = 0.5d0*(umins-uplus)/r1
              chi1(-1) = (umins+uplus)/(2.d0*akap)
           endif
        endif
        if (lzero) then
           do  30  l = -1, lmax
              chi1(l) = 0
30         enddo
           lzero = .false.
        else
           gl = cc*dexp(-a2*r2)/ta2
           do  32  l = 1, lmax
              chi1(l) = ((2*l-1)*chi1(l-1)-e*chi1(l-2)-gl)/r2
              gl = ta2*gl
32         enddo
        endif

        !        chi2 is complex; so is chi1, but the imaginary part
        !        is the same, so the difference is real
        if (lpos) then
           zuplus = zerfc(zikap/tasm+r1*asm)*expikr
           chi2(0) = expikr/r1 - dble(zuplus)/r1
           chi2(-1) = expikr/zikap + dimag(zuplus)/kap
           !          zchi2(0) = expikr/r1 - dble(zuplus)/r1
           !          zchi2(-1) = expikr/zikap + dimag(zuplus)/kap
        else
           if (.5d0*akap/asm+r1*asm > srfmax .AND. &
                .5d0*akap/asm-r1*asm > srfmax .OR. akap*r1 > fmax)then
              lzero = .true.
           else
              uplus = derfc(.5d0*akap/asm+r1*asm)/emkr
              umins = derfc(.5d0*akap/asm-r1*asm)*emkr
              chi2(0) = 0.5d0*(umins-uplus)/r1
              chi2(-1) = (umins+uplus)/(2d0*akap)
           endif
        endif
        if (lzero) then
           do  40  l = -1, lmax
              chi2(l) = 0
40         enddo
           lzero = .false.
        else
           gl = ccsm*dexp(-asm2*r2)/tasm2
           do  33  l = 1, lmax
              chi2(l) = ((2*l-1)*chi2(l-1)-e*chi2(l-2)-gl)/r2
              gl = tasm2*gl
33         enddo
        endif
     endif
     qdotr = tpi*(q(1)*dlv(1,ir)+q(2)*dlv(2,ir)+q(3)*dlv(3,ir))
     cfac = cdexp(dcmplx(0d0,qdotr))
     ilm = 0
     do  38  l = 0, lmax
        nm = 2*l+1
        do  39  m = 1, nm
           ilm = ilm+1
           dl(ilm) = dl(ilm) + yl(ilm)*(chi2(l)-chi1(l))*cfac
           dlp(ilm) = dlp(ilm) + yl(ilm)*0.5d0*(chi2(l-1)-chi1(l-1))*cfac
39      enddo
        cfac = cfac*dcmplx(0d0,1d0)
38   enddo
20 enddo

end subroutine hsmbld
end module m_smhankel



!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
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
  integer :: k0,kmax,nlmh
  double precision :: eh(1),rsmg,rsmh(1),c(0:k0,nlmh)
  integer :: ndim,ktop0,ilm,k,l,ll,lmax,m,nm
  double precision :: a,dfact,eg,fac,factk,fpi,sig
  real(8),allocatable:: s(:,:)
  if (nlmh == 0) return
  fpi = 16d0*datan(1d0)
  ktop0 = kmax
  ndim  = nlmh
  allocate(s(ndim,0:ktop0))
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
  deallocate(s)
end subroutine hxpos

subroutine ghios(rsmg,rsmh,eg,eh,nlmh,kmax,ndim,s)
  !- Integrals between 3D Gaussians and smooth Hankels at same site.
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   rsmg  :smoothing radius of gaussian
  !i   rsmh  :vector of l-dependent smoothing radii of smoothed hankel
  !i         :rsmh must be specified for 1..ll(nlmh)+1
  !i   eg    :gkL scaled by exp(eg*rsm**2/4)
  !i   eh    :vector of l-dependent energies of smoothed Hankel
  !i         :eh must be specified for 1..ll(nlmh)+1
  !i   nlmh  :L-cutoff for smoothed Hankel functions and P_kL
  !i   kmax  :polynomial cutoff
  !i   ndim  :nl*nl*nbas
  !o Outputs
  !o   s     :real-space structure constants c(k,M,L); see Remarks
  !o         :s(ilm,*) for rsmh(ll(ilm)) <= 0 are UNTOUCHED
  !r Remarks
  !r   s(L,M,k) contains integral of G_L^*(r) (lap)^k H_M(r)
  !r
  !r   Only diagonal elements ilm=jlm are nonzero and returned.
  !r   Equivalent to ghiml called for pg=ph.
  !u Updates
  !u   18 May 00 Made rsmh,eh l-dependent
  !u   24 Apr 00 Adapted from nfp ghi_os.f
  ! ----------------------------------------------------------------------
  integer :: kmax,ndim,nlmh
  double precision :: eg,eh(*),rsmg,rsmh(*),s(ndim,0:kmax)
  integer :: ktop0,ilm,k,ktop,l,ll,lmaxh,l1,l2,ilm1,ilm2
  parameter( ktop0=50 )
  double precision :: fac,fpi,gamg,gamh,rsmx,y0,h0k(0:ktop0),rsm,e
  if (nlmh == 0) return
  fpi = 16d0*datan(1d0)
  y0 = 1d0/dsqrt(fpi)
  lmaxh = ll(nlmh)
  l2 = -1
  do  20  l1  = 0, lmaxh ! --- Loop over sequences of l with a common rsm,e ---
     if (l1 <= l2) goto 20
     call gtbsl2(l1,lmaxh,eh,rsmh,l2)
     rsm  = rsmh(l1+1)
     e    = eh(l1+1)
     if (rsm <= 0 .OR. e > 0) goto 20
     ilm1 = l1**2+1
     ilm2 = (l2+1)**2
     gamh = 0.25d0*rsm*rsm
     gamg = 0.25d0*rsmg*rsmg
     rsmx = 2d0*dsqrt(gamg+gamh)
     !   ... Make hankels for l=0 and k=0..kmax
     ktop = l2+kmax
     if(ktop > ktop0) call rxi('ghios: increase ktop0, need',ktop)
     call hklos(rsmx,e,ktop,h0k)
     !   ... Evaluate what is left of Clebsch-Gordan sum
     fac = y0*dexp(gamg*(eg-e))
     do  ilm = ilm1, ilm2
        l = ll(ilm)
        do  k = 0, kmax
           s(ilm,k) = fac * (-1)**l * h0k(k+l)
        enddo
     enddo
20 enddo
end subroutine ghios

subroutine hklos(rsm,e,kmax,h0k)
  !- k,L-dependent smooth hankel functions at (0,0,0)
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   rsm   :smoothing radius
  !i   e     :energy of smoothed Hankel
  !i   kmax  :polynomial cutoff
  !o Outputs
  !i   h0k   :Smoothed Hankel at origin
  !r Remarks
  !r   Only the functions for l=0 are generated (remaining are zero)
  !r   Equivalent to hkl_ml for p=(0,0,0) and lmax=0.
  !u Updates
  !u   24 Apr 00 Adapted from nfp hkl_os.f
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: kmax
  double precision :: rsm,e,h0k(0:kmax)
  ! ... Local parameters
  integer :: k
  double precision :: akap,asm,cc,derfc,gam,gg,hh,pi,y0
  pi = 4d0*datan(1d0)
  y0 = 1d0/dsqrt(4*pi)
  gam = rsm*rsm/4d0
  asm = 0.5d0/dsqrt(gam)
  ! ... Make smooth Hankel at zero
  if (e > 0d0) call rx('hklos: e is positive')
  akap = dsqrt(dabs(e))
  cc = 4d0*asm**3 * dexp(gam*e) / dsqrt(pi)
  hh = cc/(2*asm*asm)-akap*derfc(akap/(2*asm))
  h0k(0) = hh*y0
  ! ... Upward recursion for laplace**k h0
  gg = dexp(gam*e) * (asm*asm/pi)**1.5d0
  do k = 1,kmax
     hh = -e*hh - 4*pi*gg
     h0k(k) = hh*y0
     gg = -2*asm*asm*(2*k+1)*gg
  enddo
end subroutine hklos

subroutine gtbsl2(l1,lmxh,eh,rsmh,l2)
  !- Returns the highest l with rsm,e common to those of a given l
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   l1    :current l
  !i   lmxh :basis l-cutoff
  !i   eh    :energy of smoothed Hankel
  !i   rsmh  :smoothing radius of smoothed hankel
  !o Outputs
  !o   l2    :large l for which eh and rsmh l1..l2 are in common
  !r Remarks
  !r   Routine used group functions, strux in blocks.
  !u Updates
  ! ----------------------------------------------------------------------
  !     implicit none
  integer :: l1,l2,lmxh
  double precision :: rsmh(0:lmxh),eh(0:lmxh)
  double precision :: e,rsm
  e = eh(l1)
  rsm = rsmh(l1)
  l2 = l1
10 if (l2 >= lmxh) return
  if (rsmh(l2+1) /= rsm .OR. eh(l2+1) /= e) return
  l2 = l2+1
  goto 10
end subroutine gtbsl2
