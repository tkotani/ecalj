! FCPP#define F90 1
subroutine hklbl(p,rsm,e,q,kmax,nlm,k0,cy,hkl) !slat,
  use m_lattic,only: rv_a_oqlv,rv_a_odlv,lat_plat
  use m_struc_def  !Cgetarg
  use m_lmfinit,only: lat_tol, &
       lat_alat
  use m_lattic,only: lat_qlat
  use m_lattic,only: lat_vol
  use m_lattic,only: lat_awald
  use m_lattic,only: lat_nkd
  use m_lattic,only: lat_nkq
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
  !i   slat  :struct containing information about the lattice
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
  !      type(s_lat)::slat
  double complex hkl(0:k0,nlm)
  ! ... Local parameters
  integer:: nlm0 , ilm , job , k , ll , lmax , nkd , nkq , nrx &
       , owk , oyl
  parameter (nlm0=144)
  double precision :: alat,awald,fpi,pi,sp,tol,vol,plat(3,3), &
       qlat(3,3),p1(3)
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
  call shorbz(p,p1,plat,qlat)
  sp = 2*pi*(q(1)*(p(1)-p1(1))+q(2)*(p(2)-p1(2))+q(3)*(p(3)-p1(3)))
  phase = dcmplx(dcos(sp),dsin(sp))
  awald=lat_awald
  tol=lat_tol
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


