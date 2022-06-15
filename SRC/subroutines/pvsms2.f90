subroutine pvsms2 ( ssite , sspec , rotm , nbas , nsp , sv_p_orhoat )
  use m_struc_def
  use m_lmfinit,only: stdo,slabl


  !- Rotate local densities by specified rotation
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   ssite :struct for site-specific information; see routine usite
  !i     Elts read: spec
  !i     Stored:    *
  !i     Passed to: *
  !i   sspec :struct for species-specific information; see routine uspec
  !i     Elts read: nr lmxl
  !i     Stored:    name
  !i     Passed to: spacks
  !i   rotm  :3x3 cartesian rotation matrix
  !i   nbas  :size of basis
  !i   nsp   :number of spin channels
  !i   orhoat:vector of offsets containing site density
  !o Outputs
  !o   orhoat:On output the different m-channels of rhoat(1) and rhoat(2)
  !o         :are mixed by the rotation
  !l Local variables
  !r Remarks
  !r   For a rotation matrix R, The density is stored in the 1-center form
  !r      rho_l(r) YL(rhat)
  !r   Given a rotation matrix R, this it transforms as
  !r      rho_l(r) YL(R rhat) = rho_l(r) rYL(rhat)
  !r   where rYL is made by ylmrtg
  !r
  !b Bugs
  !b   No ability is supplied when the Yl are true instead of real
  !b   spherical harmonics
  !u Updates
  !u   21 Dec 04 First created
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer:: nbas , nsp
  type(s_rv1) :: sv_p_orhoat(3,1)
  type(s_site)::ssite(*)
  type(s_spec)::sspec(*)
  double precision :: rotm(3,3)
  integer :: ib,i,j,is,lmxl,igetss,nr,nlml,ipr,nlx,nl2,nglob
  parameter (nlx=9, nl2=nlx*nlx)
  double precision :: rYL(nl2,nl2)
  character spid*8
  call getpr(ipr)
  ! ... Rotation matrix for real spherical harmonics
  call ylmrtg(nl2,rotm,rYL)
  ! --- For each site and l, rotate the m-components ---
  if (ipr >= 20) then
     call info0(20,0,0,' Rotate local densities using R=')
     write (stdo,350) ((rotm(i,j),j=1,3),i=1,3)
350  format(3f11.6)
  endif
  do  10  ib = 1, nbas
     is = ssite(ib)%spec
     spid=slabl(is) !sspec(is)%name
     nr  =sspec(is)%nr
     lmxl=sspec(is)%lmxl
     if (lmxl == -1) goto 10
     nlml = (lmxl+1)**2
     if (nlml > nl2) call rx('increase nl2 in pvsms2')
     call pvsms3 ( nr , nr , nlml , nsp , ryl , nl2 , sv_p_orhoat( 1 , ib )%v )
     call pvsms3 ( nr , nr , nlml , nsp , ryl , nl2 , sv_p_orhoat( 2 , ib )%v )
10 enddo
end subroutine pvsms2

subroutine pvsms3(nrx,nr,nlml,nsp,rYL,nl2,rho)
  !- Rotation of an l-dependent density
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nrx   :leading dimension of rho
  !i   nr    :number of radial mesh points
  !i   nlml  :L-cutoff for charge density on radial mesh
  !i   nsp   :2 for spin-polarized case, otherwise 1
  !l   rYL   :rotation matrix that rotates Y_lm
  !i   nl2   :leading dimension of rYL
  !o Outputs
  !o   rho   :On output the different m-channels of rho are
  !o         :mixed by rYL
  ! ----------------------------------------------------------------------
  implicit none
  integer :: nrx,nr,nlml,nsp,nl2
  double precision :: rho(nrx,nlml,nsp)
  double precision :: rYL(nl2,nl2)
  integer :: isp
  double precision :: rwk(nrx,nlml)
  if (nlml == 0) return
  do  isp = 1, nsp
     call dgemm('N','T',nr,nlml,nlml,1d0,rho(1,1,isp),nrx, &
          rYL,nl2,0d0,rwk,nrx)
     call dcopy(nrx*nlml,rwk,1,rho(1,1,isp),1)
  enddo
end subroutine pvsms3


