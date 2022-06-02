subroutine hgugbl(p1,p2,rsm1,rsm2,e1,nlm1,nlm2,ndim1,ndim2,  s,ds)
  use m_lmfinit,only: rv_a_ocy,rv_a_ocg, iv_a_oidxcg, iv_a_ojcg
  use m_struc_def           !Cgetarg
  use m_lattic,only:lat_vol
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
  ! ... Passed parameters
  integer :: nlm1,nlm2,ndim1,ndim2
  real(8):: rsm1 , rsm2 , p1(3) , p2(3) , e1
  complex(8):: s(ndim1,ndim2),ds(ndim1,ndim2,3)
  integer:: kmax , kdim , ilm2 , ilm1
  double precision :: q(3),e2,vol
  q=[0d0,0d0,0d0]
  vol=lat_vol
  kmax = 0
  kdim = 0
  e2 = 0d0
  call hhigbl ( 0 , p1 , p2 , q , rsm1 , rsm2 , e1 , e2 , nlm1 &
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


