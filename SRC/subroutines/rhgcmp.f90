subroutine rhgcmp(mode,ib1,ib2,ssite,sspec,sv_p_orhoat,kmax,ng,cg)!slat,
  use m_supot,only: rv_a_ogv
  use m_struc_def  !Cgetarg
  !      use m_lmfinit,only: globalvariables


  !      use m_globalvariables
  use m_lmfinit,only: lat_alat,nspec,nsp
  use m_lattic,only: lat_qlat
  use m_lattic,only: lat_vol
  use m_supot,only: lat_nabc
  use m_supot,only: lat_ng


  use m_lattic,only:lat_plat
  !- Adds density of compensating gaussians to FT list
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   mode  : a compound of digits specifying what is to be included
  !i         : in the expansion coefficients
  !i         : 1s   digit = 1 add local density rho1-rho2
  !i         :              2 add local density rho1
  !i         :              3 add local density rho2
  !i         : 10s  digit = 1 add core density rhoc
  !i         :              2 add -1 * core density from sm-hankel
  !i         :                in the local density, restoring it
  !i         :                by adding the sm-hankel to the FT mesh
  !i         :              3 combination 1+2
  !i         : 100s digit = 1 add -1 * nuclear density Z delta(r)
  !i         :                In this mode, Z is smoothed into the G_kL
  !i         :              2 add -1 * nuclear density Z delta(r)
  !i         :                In this mode, Z is incporporated directly
  !i         :                in a PW expansion (Z is not smoothed).
  !i         :
  !i         :Examples:
  !i         :mode=130 include the core, the core tail and nuclear charges
  !i         :         This should make the system charge-neutral.
  !i         :mode=131 Like mode=130, but exclude nuclear charge.
  !i         :         The system should have net charge sum_z
  !i         :mode=2   Exclude all core charges, i.e. gaussian (qcorg-z)
  !i         :  and qcorh from the foca Hankel density.
  !i         :  The system should have the valence charge.
  !i         :3 Like 0, but include nuclear charge -Z delta(r)
  !i         :  directly in a PW expansion (Z is not smoothed).
  !i   ssite :struct for site-specific information; see routine usite
  !i     Elts read: spec pos
  !i     Stored:    *
  !i     Passed to: rhogkl
  !i   sspec :struct for species-specific information; see routine uspec
  !i     Elts read: lmxl rg
  !i     Stored:    *
  !i     Passed to: corprm rhogkl
  !i   slat  :struct for lattice information; see routine ulat
  !i     Elts read: alat plat qlat nabc ng ogv okv vol
  !i     Stored:    *
  !i     Passed to: *
  !i   w(orhat):vector of offsets to local site density arrays
  !i   ng    :number of G-vectors
  !o Outputs
  !o   cg    :FT of local densities is added to cg, depending on mode.
  !r Remarks
  !r   The local charges inside each augmentation sphere
  !r   (including -1 * the core tail) are smoothed by expanding
  !r   in a  G_kL expansion for k=0..kmax.  The latter is
  !r   subsequently converted into a PW expansion.
  !u Updates
  !u   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
  !u   23 Oct 01 rhgcmp now expands local densities in
  !u             GkL for k=0..kmax, l=1..nlml for each site
  !u             Recovers old rhgcmp for kmax=0.  New argument list.
  !u   09 Feb 01 Added mode
  !u   30 May 00 Adapted from nfp rho_gcomp.f
  ! ----------------------------------------------------------------------
  implicit none
  ! ... Passed parameters
  integer :: mode,ib1,ib2,ng,kmax,i_copy_size
  !      integer orhat
  type(s_rv1) :: sv_p_orhoat(3,1)
  type(s_site)::ssite(*)
  type(s_spec)::sspec(*)
  !      type(s_lat)::slat
  double complex cg(ng)
  ! ... Local parameters
  integer :: ib,is,iv0,igetss,lmxl,ltop,n1,n2,n3,ng1,nlm, &
       nlmtop,ngabc(3),lfoc,modgkl
  real(8) ,allocatable :: qkl_rv(:)
  real(8) ,allocatable :: cs_rv(:)
  real(8) ,allocatable :: g_rv(:)
  real(8) ,allocatable :: g2_rv(:)
  integer ,allocatable :: iv_iv(:)
  real(8) ,allocatable :: sn_rv(:)
  real(8) ,allocatable :: yl_rv(:)
  equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
  double precision :: alat,ceh,cofg,cofh,qcorg,qcorh,qsc,rfoc,rg, &
       vol,z,q0(3),df(0:20),plat(3,3),qlat(3,3),tau(3)
  external corprm,poppr,pshpr,rhgcm2,rhogkl,stdfac,suphas,suphs0,suylg,tcn,tcx
  data q0 /0d0,0d0,0d0/
  call tcn('rhgcmp')
  call stdfac(20,df)
  !      i_copy_size=size(lat_plat)
  !      call dcopy(i_copy_size,lat_plat,1,plat,1)
  !      i_copy_size=size(lat_qlat)
  !      call dcopy(i_copy_size,lat_qlat,1,qlat,1)
  !      i_copy_size=size(lat_nabc)
  !      call icopy(i_copy_size,lat_nabc,1,ngabc,1)
  alat=lat_alat
  plat=lat_plat
  qlat=lat_qlat
  ngabc=lat_nabc
  ng1=lat_ng
  vol=lat_vol
  !      nspec = globalvariables%nspec
  !      nsp   = globalvariables%nsp
  modgkl = mode
  if (mode >= 200) modgkl = mod(mode,100)
  ! --- Set up help arrays ---
  ltop = 0
  do  is = 1, nspec
     lmxl = int(sspec(is)%lmxl)
     ltop = max0(ltop,lmxl)
  enddo
  nlmtop = (ltop+1)**2
  allocate(yl_rv(ng*nlmtop))
  allocate(g2_rv(ng))
  allocate(g_rv(ng*3))
  call suylg ( ltop , alat , ng , rv_a_ogv , g_rv , g2_rv , yl_rv  )
  if (allocated(g_rv)) deallocate(g_rv)
  allocate(iv_iv(ng*3))
  call suphs0 ( plat , ng , rv_a_ogv , iv_iv )
  allocate(cs_rv(ng))
  allocate(sn_rv(ng))
  ! --- Loop over sites ---
  iv0 = 0
  do  ib = ib1, ib2
     is=ssite(ib)%spec
     !        i_copy_size=size(ssite(ib)%pos)
     !     call dcopy(i_copy_size,ssite(ib)%pos,1,tau,1)
     tau=ssite(ib)%pos
     lmxl=sspec(is)%lmxl
     rg=sspec(is)%rg
     if (lmxl == -1) goto 10
     call corprm(sspec,is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
     if (mode == 2) cofh = 0
     nlm = (lmxl+1)**2
     call suphas ( q0 , tau , ng , iv_iv , n1 , n2 , n3 , qlat  , cs_rv , sn_rv )
     allocate(qkl_rv(nlm*(kmax+1)))
     call pshpr(0)
     ! ino G_kL expansion of valence sphere densities
     call rhogkl ( ib , ib , nsp , modgkl , ssite , sspec , sv_p_orhoat &
          , kmax , qkl_rv )
     call poppr
     ! ino Convert G_kL expansion of function centered at a site to PW's
     ! ino for Gaussian(G)
     ! ino   rg: gam=1/4*rg**2
     ! ino      sqkl=sqkl + qkl(k,ilm)*fac
     ! ino   cg(i) = cg(i) + sqkl*cc*yl(i,ilm)
     ! ino for H(G)
     ! ino   rfoc: gamf=1/4*rfoc**2
     ! ino   ceh:  exp(gamf*(ceh-g2))/ceh-g2
     ! ino   cofh: cg(i)=cg(i)+cofh*aa*phase
     call rhgcm2 ( vol , rg , rfoc , ceh , cofh , kmax , mod ( mode &
          / 10 , 10 ) .ge.2 , qkl_rv , nlm , ng , g2_rv , yl_rv &
          , cs_rv , sn_rv , cg )
     ! ino PW expansion of Z * delta(r)
     if ( mode >= 200 ) call rhgcm3 ( - z , vol , ng , cs_rv &
          , sn_rv , cg )
     if (allocated(qkl_rv)) deallocate(qkl_rv)
     iv0 = iv0+nlm
10   continue
  enddo
  if (allocated(sn_rv)) deallocate(sn_rv)
  if (allocated(cs_rv)) deallocate(cs_rv)
  if (allocated(iv_iv)) deallocate(iv_iv)
  if (allocated(g2_rv)) deallocate(g2_rv)
  if (allocated(yl_rv)) deallocate(yl_rv)
  call tcx('rhgcmp')
end subroutine rhgcmp

subroutine rhgcm2(vol,rg,rfoc,ceh,cofh,kmax,lcor,qkl,nlm,ng,g2,yl, &
     cs,sn,cg)
  !- Convert G_kL expansion of function centered at a site to PW's
  !     implicit none
  ! ... Passed parameters
  integer :: ng,nlm,kmax
  logical :: lcor
  double precision :: ceh,cofh,rfoc,rg,vol,qkl(0:kmax,nlm)
  double precision :: g2(ng),yl(ng,1),cs(ng),sn(ng)
  double complex cg(ng)
  ! ... Local parameters
  integer :: i,ilm,l,ll,lmxl,m,k
  double precision :: aa,cfoc,cvol,gam,gamf,pi,y0,fac,sqkl
  double complex phase,cc

  if (nlm == 0) return
  lmxl = ll(nlm)
  pi = 4d0*datan(1d0)
  y0 = 1d0/dsqrt(4d0*pi)
  gam = 0.25d0*rg*rg
  gamf = 0.25d0*rfoc*rfoc
  cvol = 1d0/vol
  cfoc = -4d0*pi*y0/vol
  do  i = 1, ng
     phase = dcmplx(cs(i),sn(i))
     aa = dexp(-gam*g2(i))*cvol
     cc = aa*phase*(0d0,1d0)
     ilm = 0
     do  l = 0, lmxl
        cc = cc*(0d0,-1d0)
        do m = -l,l
           ilm = ilm+1
           fac = 1
           sqkl = 0
           do  k = 0, kmax
              sqkl = sqkl + qkl(k,ilm)*fac
              fac = -g2(i)*fac
           enddo
           cg(i) = cg(i) + sqkl*cc*yl(i,ilm)
        enddo
     enddo

     if (lcor) then
        aa = cfoc*dexp(gamf*(ceh-g2(i)))/(ceh-g2(i))
        cg(i) = cg(i) + cofh*aa*phase
     endif

  enddo

end subroutine rhgcm2

subroutine rhgcm3(z,vol,ng,cs,sn,cg)

  !- PW expansion of Z * delta(r)
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   z     :size of delta-function
  !i   vol   :cell volume
  !i   ng    :number of G-vectors
  !i   cs    :cos(-p*G)
  !i   sn    :cos(-p*G)
  !o Outputs
  !o   cg    :Clebsch Gordon coefficients, stored in condensed form (scg.f)
  !l Local variables
  !l         :
  !r Remarks
  !r
  !u Updates
  !u   26 Oct 01
  ! ----------------------------------------------------------------------

  !     implicit none
  ! ... Passed parameters
  integer :: ng
  double precision :: z,vol,cs(ng),sn(ng)
  double complex cg(ng)
  ! ... Local parameters
  integer :: i
  double complex phase

  do  i = 1, ng
     phase = dcmplx(cs(i),sn(i))
     cg(i) = cg(i) + z*phase/vol
  enddo

end subroutine rhgcm3

