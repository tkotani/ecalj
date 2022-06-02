subroutine gklbl(p,rsm,e,q,kmax,nlm,k0,cy,gkl) !,slat
  use m_lattic,only: rv_a_oqlv,rv_a_odlv,lat_plat
  use m_struc_def           !Cgetarg
  use m_lmfinit,only: lat_alat
  use m_lattic,only: lat_qlat
  use m_lattic,only: lat_vol
  use m_lattic,only: lat_awald
  use m_lattic,only: lat_nkd
  use m_lattic,only: lat_nkq
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
  !      type(s_lat)::slat
  double complex gkl(0:k0,nlm)
  ! ... Local parameters
  double complex cfac,phase
  integer:: ilm , k , ll , lmax , nkd , nkq , owk , l , m
  double precision :: alat,awald,pi,sp,vol,plat(3,3),qlat(3,3),p1(3),rwald
  real(8),allocatable:: wk(:)
  if (nlm == 0) return
  pi = 4d0*datan(1d0)
  alat=lat_alat
  !      i_copy_size=size(lat_plat)
  !      call dcopy(i_copy_size,lat_plat,1,plat,1)
  !      i_copy_size=size(lat_qlat)
  !     call dcopy(i_copy_size,lat_qlat,1,qlat,1)
  plat=lat_plat
  qlat=lat_qlat
  awald=lat_awald
  vol=lat_vol
  nkd=lat_nkd
  nkq=lat_nkq
  ! ... Shorten p, multiply by corresponding phase later
  call shorbz(p,p1,plat,qlat)
  sp = 2*pi*(q(1)*(p(1)-p1(1))+q(2)*(p(2)-p1(2))+q(3)*(p(3)-p1(3)))
  phase = dcmplx(dcos(sp),dsin(sp))
  lmax = ll(nlm)
  rwald = 1d0/awald
  ! ... If the smoothing radius is larger than the Ewald parameter,
  !     make the summation in reciprocal space, else in real space
  if (rsm >= rwald) then
     call gklblq ( p1 , rsm , q , kmax , nlm , k0 , alat , rv_a_oqlv &
          , nkq , vol , gkl )
  else
     allocate(wk((kmax+1)*(lmax+1)))
     call gklbld ( p1 , rsm , q , kmax , nlm , k0 , alat , rv_a_odlv &
          , nkd , wk , gkl )

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


