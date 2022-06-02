subroutine fklbl(p,rsm,kmax,nlm,k0,cy,fkl)
  use m_lattic,only: rv_a_odlv,rv_a_oqlv,lat_plat
  use m_lmfinit,only: lat_alat,lat_tol
  use m_lattic,only: lat_qlat
  use m_lattic,only: lat_vol
  use m_lattic,only: lat_awald
  use m_lattic,only: lat_nkd
  use m_lattic,only: lat_nkq
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
  ! ... Passed parameters
  integer :: kmax,nlm,k0
  real(8):: p(3) , cy(*) , rsm
  !      type(s_lat)::slat
  double complex fkl(0:k0,nlm)
  ! ... Local parameters
  integer:: nlm0 , lmax , ll , nkd , nkq , nrx , owk , oyl , job &
       , ilm , k
  parameter ( nlm0=196 )
  double precision :: q(3),alat,plat(3,3),qlat(3,3),p1(3)
  double precision :: faca,fpi,y0,e,vol,awald,tol
  double complex fsm(nlm0),gklsav,gklnew
  parameter (faca=1d0)
  real(8),allocatable:: wk(:),yl(:)

  if (nlm == 0) return
  fpi = 16d0*datan(1d0)
  y0 = 1d0/dsqrt(fpi)
  if (nlm > nlm0) call rx('fklbl: increase nlm0')
  lmax = ll(nlm)
  e = 0d0
  q(1) = 0d0
  q(2) = 0d0
  q(3) = 0d0

  alat=lat_alat
  !      i_copy_size=size(lat_plat)
  !      call dcopy(i_copy_size,lat_plat,1,plat,1)
  plat=lat_plat
  !      i_copy_size=size(lat_qlat)
  !      call dcopy(i_copy_size,lat_qlat,1,qlat,1)
  qlat=lat_qlat
  awald=lat_awald
  tol=lat_tol
  vol=lat_vol
  nkd=lat_nkd
  nkq=lat_nkq
  call shorbz(p,p1,plat,qlat)
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
  ! ... Add extra term to F(k=1,l=0)
  fkl(1,1) = fkl(1,1) + fpi*y0/vol
end subroutine fklbl


