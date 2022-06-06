subroutine rhgcmp(mode,ib1,ib2,ssite,sspec,sv_p_orhoat,kmax,ng,cg)
  use m_struc_def  
  use m_supot,only: rv_a_ogv
  use m_lmfinit,only: lat_alat,nspec,nsp
  use m_lattic,only: lat_qlat,lat_vol,lat_plat
  use m_supot,only: lat_nabc,lat_ng
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
  real(8) ,allocatable :: qkl_rv(:,:)
  real(8) ,allocatable :: cs_rv(:)
  real(8) ,allocatable :: g_rv(:)
  real(8) ,allocatable :: g2_rv(:)
  integer ,allocatable :: iv_iv(:)
  real(8) ,allocatable :: sn_rv(:)
  real(8) ,allocatable :: yl_rv(:,:)
  equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
  double precision :: alat,ceh,cofg,cofh,qcorg,qcorh,qsc,rfoc,rg, &
       vol,z,q0(3),df(0:20),plat(3,3),qlat(3,3),tau(3)
  external corprm,poppr,pshpr,rhogkl,stdfac,suphas,suphs0,suylg,tcn,tcx
  data q0 /0d0,0d0,0d0/
  call tcn('rhgcmp')
  call stdfac(20,df)
  alat=lat_alat
  plat=lat_plat
  qlat=lat_qlat
  ngabc=lat_nabc
  ng1=lat_ng
  vol=lat_vol
  modgkl = mode
  if (mode >= 200) modgkl = mod(mode,100)
  ! --- Set up help arrays ---
  ltop = 0
  do  is = 1, nspec
     lmxl = sspec(is)%lmxl
     ltop = max0(ltop,lmxl)
  enddo
  nlmtop = (ltop+1)**2
  allocate(yl_rv(ng,nlmtop))
  allocate(g2_rv(ng))
  allocate(g_rv(ng*3))
  call suylg ( ltop , alat , ng , rv_a_ogv , g_rv , g2_rv , yl_rv  )
  if (allocated(g_rv)) deallocate(g_rv)
  allocate(iv_iv(ng*3))
  call suphs0 ( plat , ng , rv_a_ogv , iv_iv )
  allocate(cs_rv(ng))
  allocate(sn_rv(ng))
  iv0 = 0
  do  ib = ib1, ib2
     is=ssite(ib)%spec
     tau=ssite(ib)%pos
     lmxl=sspec(is)%lmxl
     rg=sspec(is)%rg
     if (lmxl == -1) goto 10
     call corprm(sspec,is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)
     if (mode == 2) cofh = 0
     nlm = (lmxl+1)**2
     call suphas ( q0 , tau , ng , iv_iv , n1 , n2 , n3 , qlat  , cs_rv , sn_rv )
     allocate(qkl_rv(0:kmax,nlm))
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
!     call rhgcm2 ( vol , rg , rfoc , ceh , cofh , kmax , mod ( mode &
!          / 10 , 10 ) .ge.2 , qkl_rv , nlm , ng , g2_rv , yl_rv &
!          , cs_rv , sn_rv , cg )
     rhgcm2: block !!- Convert G_kL expansion of function centered at a site to PW's
       logical :: lcor
       complex(8):: phase(ng),cc
       integer :: i,ilm,l,m,k
       real(8):: aa,cfoc,cvol,gam,gamf,fac,sqkl
       real(8),parameter:: pi = 4d0*datan(1d0),y0 = 1d0/dsqrt(4d0*pi)
       if (nlm == 0) return
       !  lmxl = ll(nlm)
       lcor= mod(mode/ 10, 10) >=2
       gam = 0.25d0*rg*rg
       gamf = 0.25d0*rfoc*rfoc
       cvol = 1d0/vol
       cfoc = -4d0*pi*y0/vol
       phase = dcmplx(cs_rv(:),sn_rv(:))
       do  i = 1, ng
          cc = phase(i)*(0d0,1d0)*dexp(-gam*g2_rv(i))*cvol
          ilm = 0
          do  l = 0, lmxl
             cc = cc*(0d0,-1d0)
             do m = -l,l
                ilm = ilm+1
                fac = 1d0
                sqkl = 0d0
                do  k = 0, kmax
                   sqkl = sqkl + qkl_rv(k,ilm)*fac
                   fac = -g2_rv(i)*fac
                enddo
                cg(i) = cg(i) + sqkl*cc*yl_rv(i,ilm)
             enddo
          enddo
          if(lcor) cg(i) = cg(i) + cofh*cfoc*dexp(gamf*(ceh-g2_rv(i)))/(ceh-g2_rv(i))*phase(i)
       enddo
     end block rhgcm2
     if( mode >= 200 ) cg(:) = -z*dcmplx(cs_rv(:),sn_rv(:))/vol !- PW expansion of Z * delta(r)
     deallocate(qkl_rv)
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
