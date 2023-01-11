module m_makusq
  public makusq
  private
  contains
subroutine makusq(mode,nsites,isite, nev,isp,iq,q,evec,  aus)
  use m_lmfinit,only: ispec,sspec=>v_sspec,nbas,nlmax,nsp,nspc,nkapii,lhh,rsma
  use m_suham,only: ndham=>ham_ndham
  use m_igv2x,only: napw,ndimh,ndimhx,igvapw=>igv2x
  use m_mkpot,only: ppnl=>ppnl_rv
  use m_uspecb,only:uspecb
  use m_lattic,only: rv_a_opos
  !- Accumulate coefficients of (u,s) in all augmentation spheres at one k-pt
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   mode  :0 generate coefficients to values, slopes
  !i         :1 generate coefficients to phi,phidot
  !i   nbas  :number of basis atoms
  !i   nsites:If zero, coefficients are made all sites.
  !i         :If nonzero, coefficients are made just for a subset
  !i         :of sites (see isite); nsites is the number of sites
  !i   isite :sites at which to calculate coefficients; see nsites
  !i   nlmax :1st dimension of aus (maximum nlma over all sites)
  !i   ndham :dimensions aus
  !i   ndimh :dimensions evec
  !i   napw  :number of G vectors in PW basis (gvlst2.f)
  !i   igvapw:G vectors in PW basis, units of qlat (gvlst2.f)
  !i   nev   :number of eigenvectors for which to accumulate aus
  !i   nsp   :2 for spin-polarized case, otherwise 1
  !i   nspc  :2 for coupled spins; otherwise 1
  !i   isp   :spin channel, used only to address element in aus
  !i   iq    :qp index, used only to address element in aus
  !i   q     :Bloch vector
  !i   evec  :eigenvectors for this q
  !i   ppnl  :nmto-like pot pars
  !o Outputs
  !o   aus   :val,slo of w.f. at MT sphere surface added to aus; see Remarks
  !l Local variables
  !l   ispc  :the current spin index in the coupled spins case.
  !l         :Some quantities have no separate address space for each
  !l         :spin in the indepedent-spins case (evec,evl,ewgt) but do
  !l         :in the coupled-spins case.  A separate loop ispc=1..nspc
  !l         :must be added for the latter case
  !l         :ispc is the appropriate index for objects which distinguish
  !l         :spins in the spin-coupled case only
  !l   isp   :isp  is the appropriate index for objects which distinguish
  !l         :spins in the spin-uncoupled case only
  !l   ksp   :the current spin index in both independent and coupled
  !l         :spins cases.
  !l         :ksp is appropriate spin index for quantities that have
  !l         :separate address space for each spin in every case
  !l         :(potential- and density-like objects).
  !r Remarks
  !r   Makes coefficients for projection of wave function onto
  !r   augmented functions (u,s) which is valid inside the MT spheres.
  !r   u and s are linear combinations of and phi,phidot defined as:
  !r   u has val=1, slo=0 at rmax, s has val=0, slo=1
  !r
  !r   For example, for EELS matrix elements <nk|r|core> we will need
  !r    |nk> = \sum_L(au_nkL*u_l*Y_L + as_nkL*s_l*Y_L)
  !r
  !r   These are generated from the potential later (see vcdmel)
  !r   makusq returns the au_nkL and as_nkL at one spin and k-pt for
  !r   each of the sites in the unit cell, but ...
  !r   if nsites=nbas is passed then coeffs made for each site, otherwise
  !r   coeffs are made just for the nsites sites listed in isite (suclst)
  !u Updates
  !u   06 Jan 09 Adapt to include APW basis functions
  !u   08 Jul 08 Dimension aus separately from evec
  !u   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
  !u   23 Dec 04 Extended to the spin-coupled case
  !u   25 Aug 04 Adapted to extended local orbitals
  !u   21 Aug 03 Restrict to a list of sites (see Remarks and suclst)
  !u   12 Feb 02 Extended to local orbitals
  !u   28 Mar 01 (MvS) Added mode to generate coefficients to phi,phidot
  !u                   Some rearrangement of coefficients.
  !u   19 Feb 01 (MvS) shortened argument list
  !u   21 Nov 00 (ATP) Adapted from fp mkrout
  ! ----------------------------------------------------------------------
  implicit none
  real(8):: nlmax_,ndham_,ndimh_,napw_,igvapw_,nsp_,nspc_,ppnl_ !dummy input
  integer:: mode,isp,iq, nev,n0,nppn, nsites,isite(nsites)
  parameter (n0=10,nppn=12)
  real(8):: q(3)
  double complex evec(ndimh,nsp,nev), &
       aus(nlmax,ndham*nspc,3,nsp,nsites,iq)
  ! ... Local parameters
  integer :: nkap0
  parameter (nkap0=3)
  double precision :: eh(n0,nkap0),rsmh(n0,nkap0),a,rmt
  integer :: igetss,ib,nkapi,is,nr,kmax,lmxa,lmxl,lmxh,i
  real(8) ,allocatable :: rofi_rv(:)
  real(8) ,allocatable :: fh_rv(:)
  real(8) ,allocatable :: xh_rv(:)
  real(8) ,allocatable :: vh_rv(:)
  real(8) ,allocatable :: dh_rv(:)
  real(8) ,allocatable :: fp_rv(:)
  real(8) ,allocatable :: xp_rv(:)
  real(8) ,allocatable :: vp_rv(:)
  real(8) ,allocatable :: dp_rv(:)
  call tcn ('makusq')
  ! --- Start loop over atoms ---
  do  i = 1, nsites
     if (nsites == nbas) then
        ib = i
     else
        ib = isite(i)
     endif
     is = ispec(ib)
     lmxa=sspec(is)%lmxa
     lmxl=sspec(is)%lmxl
     kmax=sspec(is)%kmxt
     a=sspec(is)%a
     nr=sspec(is)%nr
     rmt=sspec(is)%rmt
     call uspecb(is,rsmh,eh)
     nkapi= nkapii(is)
     lmxh = sspec(is)%lmxb
     if (lmxa == -1) goto 10
     !   --- Set up all radial head and tail functions, and their BC's
     allocate(rofi_rv(nr))
     call radmsh ( rmt , a , nr , rofi_rv )
     allocate(fh_rv(nr*(lmxh+1)*nkapi))
     allocate(xh_rv(nr*(lmxh+1)*nkapi))
     allocate(vh_rv((lmxh+1)*nkapi))
     allocate(dh_rv((lmxh+1)*nkapi))
     call fradhd ( nkapi , eh , rsmh , lhh(:,is) , lmxh , nr , rofi_rv &
          , fh_rv , xh_rv , vh_rv , dh_rv )
     allocate(fp_rv(nr*(lmxa+1)*(kmax+1)))
     allocate(xp_rv(nr*(lmxa+1)*(kmax+1)))
     allocate(vp_rv((lmxa+1)*(kmax+1)))
     allocate(dp_rv((lmxa+1)*(kmax+1)))
     call fradpk(kmax,rsma(is),lmxa,nr,rofi_rv,fp_rv,xp_rv,vp_rv,dp_rv)
     !   --- Add to the coefficient for the projection onto (u,s) for this site
     call pusq1 ( mode , ib , isp , nspc , nlmax , lmxh &
          , nbas,  q , ndham , ndimh , napw , igvapw& ! & slat ,
          , nev , evec , vh_rv , dh_rv , vp_rv , dp_rv , ppnl ( 1 , 1 , &
          1 , ib ) , aus ( 1 , 1 , 1 , 1 , i , iq ) , aus ( 1 , 1 , 2 , &
          1 , i , iq ) , aus ( 1 , 1 , 3 , 1 , i , iq ) )
     if (allocated(dp_rv)) deallocate(dp_rv)
     if (allocated(vp_rv)) deallocate(vp_rv)
     if (allocated(xp_rv)) deallocate(xp_rv)
     if (allocated(fp_rv)) deallocate(fp_rv)
     if (allocated(dh_rv)) deallocate(dh_rv)
     if (allocated(vh_rv)) deallocate(vh_rv)
     if (allocated(xh_rv)) deallocate(xh_rv)
     if (allocated(fh_rv)) deallocate(fh_rv)
     if (allocated(rofi_rv)) deallocate(rofi_rv)
10   continue
  enddo
  call tcx('makusq')
end subroutine makusq
! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
subroutine pusq1(mode,ia,isp,nspc,nlmax,lmxh,nbas,q,ndham,ndimh,napw,igvapw,nev,evec,vh,dh,vp,dp,ppnl, &
     au,as,az)
  use m_lmfinit,only: nkapii,ispec,sspec=>v_sspec
  use m_uspecb,only: uspecb
  use m_bstrux,only: bstrux_set,bstr
  use m_lattic,only: rv_a_opos
  use m_rlocbl,only: rlocb1
  !- Add to the coefficient for the projection onto (u,s) for one site
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   mode  :0 generate coefficients to values, slopes
  !i         :1 generate coefficients to phi,phidot
  !i   ia    :augmentation sphere
  !i   isp   :current spin index for collinear case
  !i   nspc  :2 for coupled spins; otherwise 1
  !i   nlmax :dimensions au,as
  !i   lmxh  :basis l-cutoff
  !i   nbas  :size of basis
  !i   q     :bloch vector
  !i   ndham :dimensions au,as,az
  !i   ndimh :dimension of hamiltonian, evec
  !i   napw  :number of G vectors in PW basis (gvlst2.f)
  !i   igvapw:G vectors in PW basis, units of qlat (gvlst2.f)
  !i   nev   :number of eigenvectors to sum over
  !i   evec  :eigenvectors
  !i   vh    :value of head function in sphere ia
  !i   dh    :slope of head function in sphere ia
  !i   vp    :value of PkL expansion of tail function in sphere ia
  !i   dp    :slope of PkL expansion of tail function in sphere ia
  !i   ppnl  :NMTO pot pars (potpus.f)
  !l Local variables
  !l   ksp   :the current spin index in both independent and coupled
  !l         :spins cases.
  !o Outputs
  !o   au    :projection of this evec onto u function; see potpus.f
  !o         :If mode=1, au = projection of this evec onto phi function
  !o   as    :projection of this evec onto s function; see potpus.f
  !o         :If mode=1, au = projection of this evec onto phidot function
  !o   az    :projection of this evec onto local orbitals; see potpus.f
  !r Remarks
  !r   Adapted from augmbl
  !u Updates
  !u   23 Dec 04 Extended to the spin-coupled case
  !u    4 Jun 04 Relax condition nlmax>=nlma
  !u   10 Apr 02 Redimensionsed eh,rsmh to accommodate larger lmax
  !u   12 Feb 02 Extended to local orbitals
  ! ----------------------------------------------------------------------
  integer :: mode,ia,isp,nspc,lmxh,nlmax, &
       nbas,ndham,ndimh,napw,igvapw(3,napw),nev,nlmbx,n0,nppn
  parameter (nlmbx=25, n0=10, nppn=12)
  real(8):: ppnl(nppn,n0,2),q(3),vp(*),dp(*),vh(*),dh(*)
  complex(8):: evec(ndimh,nspc,ndimh), &
       au(nlmax,ndham*nspc,3,2), &
       as(nlmax,ndham*nspc,3,2), &
       az(nlmax,ndham*nspc,3,2)
  complex(8):: zdummy(1)
  integer :: nkap0,nlmxx
  parameter (nkap0=3,nlmxx=121)
  complex(8) ,allocatable :: a_zv(:) !b_zv(:),
  integer :: isa,lmxa,lmxha,kmax,nlma,ivec, ilm,k,ll,nkape,ksp,ispc,nlmto
  real(8) :: pa(3),rmt, phi,phip,dphi,dlphi,dphip,dlphip,det,rotpp(nlmxx,2,2),&
       eh(n0,nkap0),rsmh(n0,nkap0)
  call tcn ('pusq1')
  isa =ispec(ia)!%spec
  pa  =rv_a_opos(:,ia) !ssite(ia)%pos
  lmxa=sspec(isa)%lmxa
  if (lmxa == -1) return
  lmxha=sspec(isa)%lmxb
  kmax=sspec(isa)%kmxt
!  rsma=sspec(isa)%rsma
  rmt=sspec(isa)%rmt
  nlmto = ndimh-napw
  nlma  = (lmxa+1)**2
  call uspecb(isa,rsmh,eh)
  nkape = nkapii(isa)
  call bstrux_set(ia,q) !bstr
  allocate(a_zv((kmax+1)*nlma))
  !     In noncollinear case, isp=1 always => need internal ispc=1..2
  !     ksp is the current spin index in both cases:
  !     ksp = isp  in the collinear case
  !         = ispc in the noncollinear case
  !     whereas ispc=1 for independent spins, and spin index when nspc=2
  ispcloop: do  ispc = 1, nspc ! ... loop over noncollinear spins
     ksp = max(ispc,isp)
     ! if(mode == 1) then
     !    if (nlma > nlmxx) call rxi('makusq:  nlmxx < nlma=',nlma)
     !    do  ilm = 1, nlma
     !       k = ll(ilm)+1
     !       dlphi  = ppnl(3,k,ksp)/rmt
     !       dlphip = ppnl(4,k,ksp)/rmt
     !       phi    = ppnl(5,k,ksp)
     !       phip   = ppnl(6,k,ksp)
     !       dphi   = phi*dlphi/rmt
     !       dphip  = dlphip/rmt*phip
     !       det    = phi*dphip - dphi*phip
     !       rotpp(ilm,1,1) = dphip/det
     !       rotpp(ilm,1,2) = -dphi/det
     !       rotpp(ilm,2,1) = -phip/det
     !       rotpp(ilm,2,2) = phi/det
     !    enddo
     ! endif
     do  ivec = 1, nev ! --- Loop over eigenstates ---
        call rlocb1(ndimh, nlma, kmax, evec(1,ispc,ivec), bstr,a_zv) !b_zv ,a_zv ) !
        call pusq2 (mode, ia, nkape, kmax, lmxa , lmxh, nlmto , &
             min(nlma,nlmax),a_zv,rotpp,evec(1,ispc,ivec),vh,dh,vp,dp,ksp,&
             au(1,ivec,1,ksp), as(1,ivec,1,ksp), az(1,ivec,1,ksp) )
     enddo 
  enddo ispcloop
  deallocate(a_zv) 
  call tcx('pusq1')
end subroutine pusq1
subroutine pusq2(mode,ia,nkape,kmax,lmxa,lmxh,nlmto,nlma,cPkL,r,evec,vh,dh,vp,dp,ksp, au,as,az)
  use m_locpot,only: rotp
  use m_orbl,only: Orblib,ktab,ltab,offl,norb,ntab,blks
  !- Extract projection of eigenstate onto (u,s,z) for sphere at site ia
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   mode  :0 generate coefficients to values, slopes
  !i         :1 generate coefficients to phi,phidot
  !i   ia    :augmentation sphere
  !i   nkape :number of envelope function types which are joined to (u,s)
  !i         :Any ktab > nkape is a local orbital
  !i   kmax  :polynomial cutoff in P_kL expansion of envelope tails
  !i   lmxa  :augmentation l-cutoff
  !i   lmxh  :basis l-cutoff
  !i   nlmto :dimension of lmto component of basis
  !i   nlma  :number of L's in augmentation sphere = (lmxa+1)**2
  !o   cPkL  :coefficients to P_kL expansion of evec
  !i   r     :2x2 rotation matrices rotating (phi,phidot) to (u,s)
  !i   evec  :eigenvector
  !i   vh    :value of head function in sphere ia
  !i   dh    :slope of head function in sphere ia
  !i   vp    :value of PkL expansion of tail function in sphere ia
  !i   dp    :slope of PkL expansion of tail function in sphere ia
  !o Outputs
  !o   au    :projection of this evec onto u function; see potpus.f
  !o         :If mode=1, au = projection of this evec onto phi function
  !o   as    :projection of this evec onto s function; see potpus.f
  !o         :If mode=1, au = projection of this evec onto phidot function
  !o   az    :projection of this evec onto local orbitals; see potpus.f
  !l Local variables
  !l         :
  !r Remarks
  !r
  !u Updates
  !u   12 Feb 02 Extended to local orbitals
  ! ----------------------------------------------------------------------
  implicit none
  integer :: mode,ia,nkape,kmax,lmxa,lmxh,nlmto,nlma,ksp
  double precision :: vh(0:lmxh,1),dh(0:lmxh,1)
  double precision :: vp(0:lmxa,0:kmax),dp(0:lmxa,0:kmax)
  integer :: nlmxx
  parameter (nlmxx=121)
  double precision :: r(nlmxx,2,2)
  complex(8):: au(nlma),as(nlma),az(nlma), evec(nlmto),cPkL(0:kmax,nlma)
  integer :: n0,nkap0
  parameter (n0=10,nkap0=3)
!  integer :: blks(n0*nkap0),ntab(n0*nkap0)
  integer :: io1,l1,ik1,nlm11,nlm12,ilm1,i1,ilma,k
  integer :: l,ll
  double precision :: xx
  double complex wk(nlmxx,2)
  complex(8):: auas(2)
  !     call tcn('pusq2')
  if (nlmto == 0) return
  if (nlma > nlmxx) call rxi('makusq:  nlmxx < nlma=',nlma)
  ! --- Loop over all orbitals centered at this site ---
  call orblib(ia)!norb,ltab,ktab,xx,offl,xx)
  !     Block into groups of consecutive l
!  call gtbsl4(ia) !1(4,norb,ltab,ktab,xx,xx,ntab,blks)
  !     Contribution from head part
  do  io1 = 1, norb
     l1  = ltab(io1)
     ik1 = ktab(io1)
     nlm11 = l1**2+1
     nlm12 = nlm11 + blks(io1)-1
     i1 = offl(io1)-nlm11+1 !  i1 = hamiltonian offset for first orbital in block
     if (ik1 <= nkape) then
        do  ilm1 = nlm11, nlm12
           l = ll(ilm1)
           au(ilm1) = au(ilm1) + vh(l,ik1) * evec(ilm1+i1)
           as(ilm1) = as(ilm1) + dh(l,ik1) * evec(ilm1+i1)
        enddo
     else
        do  ilm1 = nlm11, nlm12
           az(ilm1) = az(ilm1) + evec(ilm1+i1)
        enddo
     endif
  enddo
  !     Contribution from tail part
  do  ilma = 1, nlma
     l = ll(ilma)
     do  k = 0, kmax
        au(ilma) = au(ilma) + vp(l,k) * cPkL(k,ilma)
        as(ilma) = as(ilma) + dp(l,k) * cPkL(k,ilma)
     enddo
  enddo
  !     Rotate to (phi,phidot)
   if(mode /= 0) then
      do  ilma = 1, nlma
         l=ll(ilma)
         auas = matmul([au(ilma),as(ilma)],rotp(l,ksp,:,:,ia))
         au(ilma) = auas(1)
         as(ilma) = auas(2)
      enddo
   endif
  !     call tcx('pusq2')
end subroutine pusq2
end module m_makusq
