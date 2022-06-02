subroutine hhugbl(mode,p1,p2,rsm1,rsm2,e1,e2,nlm1,nlm2,ndim1, &
     ndim2,wk,dwk,s,ds) !slat,
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
  real(8):: rsm1(0:*) , rsm2(0:*) , e1(0:*) , e2(0:*) , p1(3) , &
       p2(3)
  double complex s(ndim1,ndim2),ds(ndim1,ndim2,3), &
       wk(ndim1,ndim2),dwk(ndim1,ndim2,3)
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


