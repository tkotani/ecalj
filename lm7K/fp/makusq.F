      subroutine makusq(mode,ssite,sspec,slat,sham,nbas,nsites,isite,
     .nlmax,ndham,ndimh,napw,igvapw,nev,nsp,nspc,isp,iq,q,evec,ppnl,
     .aus)

      use m_struc_def  !Cgetarg

C- Accumulate coefficients of (u,s) in all augmentation spheres at one k-pt
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 generate coefficients to values, slopes
Ci         :1 generate coefficients to phi,phidot
Ci   ssite :struct containing site-specific information
Ci     Elts read: spec pos
Ci     Passed to: pusq1
Ci   sspec :struct containing species-specific information
Ci     Elts read: rsma lmxa lmxl kmxt a nr rmt lmxb
Ci     Stored:
Ci     Passed to: uspecb pusq1
Ci   slat  :struct containing information about the lattice
Ci     Elts read: ocg ojcg oidxcg ocy
Ci     Stored:
Ci     Passed to: pusq1 hxpbl
Ci   sham  :struct containing information about the hamiltonian
Ci     Elts read: oindxo
Ci   nbas  :number of basis atoms
Ci   nsites:If zero, coefficients are made all sites.
Ci         :If nonzero, coefficients are made just for a subset
Ci         :of sites (see isite); nsites is the number of sites
Ci   isite :sites at which to calculate coefficients; see nsites
Ci   nlmax :1st dimension of aus (maximum nlma over all sites)
Ci   ndham :dimensions aus
Ci   ndimh :dimensions evec
Ci   napw  :number of G vectors in PW basis (gvlst2.f)
Ci   igvapw:G vectors in PW basis, units of qlat (gvlst2.f)
ci   nev   :number of eigenvectors for which to accumulate aus
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nspc  :2 for coupled spins; otherwise 1
Ci   isp   :spin channel, used only to address element in aus
Ci   iq    :qp index, used only to address element in aus
Ci   q     :Bloch vector
Ci   evec  :eigenvectors for this q
Ci   ppnl  :nmto-like pot pars
Co Outputs
Co   aus   :val,slo of w.f. at MT sphere surface added to aus; see Remarks
Cl Local variables
Cl   ispc  :the current spin index in the coupled spins case.
Cl         :Some quantities have no separate address space for each
Cl         :spin in the indepedent-spins case (evec,evl,ewgt) but do
Cl         :in the coupled-spins case.  A separate loop ispc=1..nspc
Cl         :must be added for the latter case
Cl         :ispc is the appropriate index for objects which distinguish
Cl         :spins in the spin-coupled case only
Cl   isp   :isp  is the appropriate index for objects which distinguish
Cl         :spins in the spin-uncoupled case only
Cl   ksp   :the current spin index in both independent and coupled
Cl         :spins cases.
Cl         :ksp is appropriate spin index for quantities that have
Cl         :separate address space for each spin in every case
Cl         :(potential- and density-like objects).
Cr Remarks
Cr   Makes coefficients for projection of wave function onto
Cr   augmented functions (u,s) which is valid inside the MT spheres.
Cr   u and s are linear combinations of and phi,phidot defined as:
Cr   u has val=1, slo=0 at rmax, s has val=0, slo=1
Cr
Cr   For example, for EELS matrix elements <nk|r|core> we will need
Cr    |nk> = \sum_L(au_nkL*u_l*Y_L + as_nkL*s_l*Y_L)
Cr
Cr   These are generated from the potential later (see vcdmel)
Cr   makusq returns the au_nkL and as_nkL at one spin and k-pt for
Cr   each of the sites in the unit cell, but ...
Cr   if nsites=nbas is passed then coeffs made for each site, otherwise
Cr   coeffs are made just for the nsites sites listed in isite (suclst)
Cu Updates
Cu   06 Jan 09 Adapt to include APW basis functions
Cu   08 Jul 08 Dimension aus separately from evec
Cu   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
Cu   23 Dec 04 Extended to the spin-coupled case
Cu   25 Aug 04 Adapted to extended local orbitals
Cu   21 Aug 03 Restrict to a list of sites (see Remarks and suclst)
Cu   12 Feb 02 Extended to local orbitals
Cu   28 Mar 01 (MvS) Added mode to generate coefficients to phi,phidot
Cu                   Some rearrangement of coefficients.
Cu   19 Feb 01 (MvS) shortened argument list
Cu   21 Nov 00 (ATP) Adapted from fp mkrout
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,nbas,nsp,nspc,isp,iq,
     .ndham,ndimh,napw,igvapw(3,napw),nev,nlmax,n0,nppn,
     .nsites,isite(nsites)
      parameter (n0=10,nppn=12)
      real(8):: q(3)
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)
      type(s_lat)::slat
      type(s_ham)::sham

      double precision ppnl(nppn,n0,nsp,nbas)
      double complex evec(ndimh,nsp,nev),
     .aus(nlmax,ndham*nspc,3,nsp,nsites,iq)
C ... Local parameters
      integer nkap0
      parameter (nkap0=3)
      integer lh(nkap0)
      double precision eh(n0,nkap0),rsmh(n0,nkap0),rsma,a,rmt
      integer igetss,ib,nkapi,is,nr,kmax,lmxa,lmxl,lmxh,i
ckino Dec.12.2011:         integer,pointer :: iv_p_oiprmb(:) =>NULL()

      real(8) ,allocatable :: rofi_rv(:)
      real(8) ,allocatable :: fh_rv(:)
      real(8) ,allocatable :: xh_rv(:)
      real(8) ,allocatable :: vh_rv(:)
      real(8) ,allocatable :: dh_rv(:)
      real(8) ,allocatable :: fp_rv(:)
      real(8) ,allocatable :: xp_rv(:)
      real(8) ,allocatable :: vp_rv(:)
      real(8) ,allocatable :: dp_rv(:)

C ... Heap

C ... Setup
C     stdo = lgunit(1)
C     ipr  = iprint()
      call tcn ('makusq')

C --- Start loop over atoms ---
      do  i = 1, nsites
        if (nsites .eq. nbas) then
          ib = i
        else
          ib = isite(i)
        endif
        is = int(ssite(ib)%spec)


ckino Dec.12.2011:           iv_p_oiprmb => sham%iv_p_oindxo 



        rsma=sspec(is)%rsma
        lmxa=sspec(is)%lmxa
        lmxl=sspec(is)%lmxl
        kmax=sspec(is)%kmxt


        a=sspec(is)%a
        nr=sspec(is)%nr
        rmt=sspec(is)%rmt

        call uspecb(0,1,sspec,is,is,lh,rsmh,eh,nkapi)

        lmxh=sspec(is)%lmxb

        if (lmxa .eq. -1) goto 10

C   --- Set up all radial head and tail functions, and their BC's
        allocate(rofi_rv(nr))

        call radmsh ( rmt , a , nr , rofi_rv )

        allocate(fh_rv(nr*(lmxh+1)*nkapi))

        allocate(xh_rv(nr*(lmxh+1)*nkapi))

        allocate(vh_rv((lmxh+1)*nkapi))

        allocate(dh_rv((lmxh+1)*nkapi))

        call fradhd ( nkapi , eh , rsmh , lh , lmxh , nr , rofi_rv 
     .  , fh_rv , xh_rv , vh_rv , dh_rv )

        allocate(fp_rv(nr*(lmxa+1)*(kmax+1)))

        allocate(xp_rv(nr*(lmxa+1)*(kmax+1)))

        allocate(vp_rv((lmxa+1)*(kmax+1)))

        allocate(dp_rv((lmxa+1)*(kmax+1)))

        call fradpk ( kmax , rsma , lmxa , nr , rofi_rv , fp_rv 
     .  , xp_rv , vp_rv , dp_rv )

C   --- Add to the coefficient for the projection onto (u,s) for this site
ckino Dec.12.2011:           call pusq1 ( mode , ib , isp , nspc , iv_p_oiprmb , nlmax , lmxh 
ckino Dec.27.2011:           call pusq1 ( mode , ib , isp , nspc , sham%iv_p_oindxo , nlmax , lmxh 
         call pusq1 ( mode , ib , isp , nspc , sham%iv_a_oindxo , nlmax , lmxh 
     .   , nbas , ssite , sspec , slat , q , ndham , ndimh , napw , igvapw 
     .   , nev , evec , vh_rv , dh_rv , vp_rv , dp_rv , ppnl ( 1 , 1 , 
     .   1 , ib ) , aus ( 1 , 1 , 1 , 1 , i , iq ) , aus ( 1 , 1 , 2 , 
     .   1 , i , iq ) , aus ( 1 , 1 , 3 , 1 , i , iq ) ) 


        if (allocated(dp_rv)) deallocate(dp_rv)
        if (allocated(vp_rv)) deallocate(vp_rv)
        if (allocated(xp_rv)) deallocate(xp_rv)
        if (allocated(fp_rv)) deallocate(fp_rv)
        if (allocated(dh_rv)) deallocate(dh_rv)
        if (allocated(vh_rv)) deallocate(vh_rv)
        if (allocated(xh_rv)) deallocate(xh_rv)
        if (allocated(fh_rv)) deallocate(fh_rv)
        if (allocated(rofi_rv)) deallocate(rofi_rv)

C   --- print
C        print *, 'isp,ib,nev,iq=',isp,ib,nev,iq
C        call zprm('val',2,aus(1,1,1,isp,ib,iq),nlmax,(lmxa+1)**2,nev)
C        call zprm('slo',2,aus(1,1,2,isp,ib,iq),nlmax,(lmxa+1)**2,nev)
C        call zprm('loc',2,aus(1,1,3,isp,ib,iq),nlmax,(lmxa+1)**2,nev)
   10   continue
      enddo
      call tcx('makusq')
      end subroutine makusq


      subroutine pusq1(mode,ia,isp,nspc,iprmb,nlmax,lmxh,nbas,ssite,
     .sspec,slat,q,ndham,ndimh,napw,igvapw,nev,evec,vh,dh,vp,dp,ppnl,
     .au,as,az)

      use m_struc_def  !Cgetarg

C- Add to the coefficient for the projection onto (u,s) for one site
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 generate coefficients to values, slopes
Ci         :1 generate coefficients to phi,phidot
Ci   ia    :augmentation sphere
Ci   isp   :current spin index for collinear case
Ci   nspc  :2 for coupled spins; otherwise 1
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   nlmax :dimensions au,as
Ci   lmxh  :basis l-cutoff
Ci   nbas  :size of basis
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec pos
Ci     Stored:    *
Ci     Passed to: *
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: lmxa lmxb kmxt rsma rmt
Ci     Stored:    *
Ci     Passed to: uspecb
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: ocg ojcg oidxcg ocy
Ci     Stored:    *
Ci     Passed to: hxpbl
Ci   q     :bloch vector
Ci   ndham :dimensions au,as,az
Ci   ndimh :dimension of hamiltonian, evec
Ci   napw  :number of G vectors in PW basis (gvlst2.f)
Ci   igvapw:G vectors in PW basis, units of qlat (gvlst2.f)
Ci   nev   :number of eigenvectors to sum over
Ci   evec  :eigenvectors
Ci   vh    :value of head function in sphere ia
Ci   dh    :slope of head function in sphere ia
Ci   vp    :value of PkL expansion of tail function in sphere ia
Ci   dp    :slope of PkL expansion of tail function in sphere ia
Ci   ppnl  :NMTO pot pars (potpus.f)
Cl Local variables
Cl   ksp   :the current spin index in both independent and coupled
Cl         :spins cases.
Co Outputs
Co   au    :projection of this evec onto u function; see potpus.f
Co         :If mode=1, au = projection of this evec onto phi function
Co   as    :projection of this evec onto s function; see potpus.f
Co         :If mode=1, au = projection of this evec onto phidot function
Co   az    :projection of this evec onto local orbitals; see potpus.f
Cr Remarks
Cr   Adapted from augmbl
Cu Updates
Cu   23 Dec 04 Extended to the spin-coupled case
Cu    4 Jun 04 Relax condition nlmax>=nlma
Cu   10 Apr 02 Redimensionsed eh,rsmh to accommodate larger lmax
Cu   12 Feb 02 Extended to local orbitals
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,ia,isp,nspc,lmxh,nlmax,
     .nbas,ndham,ndimh,napw,igvapw(3,napw),nev,nlmbx,n0,nppn
      parameter (nlmbx=25, n0=10, nppn=12)
      integer iprmb(ndimh)
      double precision ppnl(nppn,n0,2)
      double precision vp(*),dp(*),vh(*),dh(*)
      real(8):: q(3)
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)
      type(s_lat)::slat

      double complex evec(ndimh,nspc,ndimh),
     .au(nlmax,ndham*nspc,3,2),
     .as(nlmax,ndham*nspc,3,2),
     .az(nlmax,ndham*nspc,3,2)
C ... Local parameters
      integer nkap0,nlmxx
      parameter (nkap0=3,nlmxx=121)
      integer lh(nkap0)
ckino Dec.12.2011:            integer,pointer :: iv_p_ojcg(:) =>NULL()

ckino Dec.12.2011:           integer,pointer :: iv_p_oidxcg(:) =>NULL()

ckino Dec.12.2011:          real(8),pointer :: rv_p_ocy(:) =>NULL()

ckino Dec.12.2011:         real(8),pointer :: rv_p_ocg(:) =>NULL()

      complex(8) ,allocatable :: b_zv(:)
      complex(8) ,allocatable :: a_zv(:)

      integer isa,lmxa,lmxha,kmax,nlma,ivec,
     .ilm,k,ll,nkape,ksp,ispc,nlmto
      double precision eh(n0,nkap0),rsmh(n0,nkap0)
      double precision rsma,pa(3),rmt,
     .phi,phip,dphi,dlphi,dphip,dlphip,det,rotp(nlmxx,2,2)
C ... Heap
       integer ::iwdummy 


      call tcn ('pusq1')

ckino Dec.12.2011:         rv_p_ocg => slat%rv_p_ocg 

ckino Dec.12.2011:         iv_p_ojcg => slat%iv_p_ojcg 

ckino Dec.12.2011:         iv_p_oidxcg => slat%iv_p_oidxcg 

ckino Dec.12.2011:         rv_p_ocy => slat%rv_p_ocy 



      isa=ssite(ia)%spec
      i_copy_size=size(ssite(ia)%pos)
      call dcopy(i_copy_size,ssite(ia)%pos,1,pa,1)


      lmxa=sspec(isa)%lmxa
      lmxha=sspec(isa)%lmxb
      kmax=sspec(isa)%kmxt
      rsma=sspec(isa)%rsma


      rmt=sspec(isa)%rmt

      if (lmxa .eq. -1) return
      nlmto = ndimh-napw

C     nlmha = (lmxha+1)**2
      nlma  = (lmxa+1)**2
C     if (nlma .gt. nlmax) call rxi('makusq: need nlmax',nlma)
C     Count no. envelope functions connecting (phi,phidot) at site ia
      call uspecb(0,1,sspec,isa,isa,lh,rsmh,eh,nkape)

C --- Make strux to expand all orbitals at site ia ---
      allocate(b_zv((kmax+1)*nlma*ndimh))

ckino Dec.12.2011:             call bstrux ( 2 , slat , ssite , sspec , rv_p_ocg , iv_p_oidxcg 
ckino Dec.12.2011:             call bstrux ( 2 , slat , ssite , sspec , rv_p_ocg , slat%iv_p_oidxcg 
ckino Dec.27.2011:             call bstrux ( 2 , slat , ssite , sspec , slat%rv_p_ocg , slat%iv_p_oidxcg 
ckino Dec.27.2011:             call bstrux ( 2 , slat , ssite , sspec , slat%rv_p_ocg , slat%iv_a_oidxcg 
           call bstrux ( 2 , slat , ssite , sspec , slat%rv_a_ocg , slat%iv_a_oidxcg 
ckino Dec.12.2011:       .     , iv_p_ojcg , rv_p_ocy , iprmb , nbas , ia , pa , rsma , q , 
ckino Dec.12.2011:       .     , slat%iv_p_ojcg , rv_p_ocy , iprmb , nbas , ia , pa , rsma , q , 
ckino Dec.27.2011:       .     , slat%iv_p_ojcg , slat%rv_p_ocy , iprmb , nbas , ia , pa , rsma , q , 
ckino Dec.28.2011:       .     , slat%iv_p_ojcg , slat%rv_a_ocy , iprmb , nbas , ia , pa , rsma , q , 
     .     , slat%iv_a_ojcg , slat%rv_a_ocy , iprmb , nbas , ia , pa , rsma , q , 
     .     kmax , nlma , ndimh , napw , igvapw , b_zv , iwdummy ) 







C     In noncollinear case, isp=1 always => need internal ispc=1..2
C     ksp is the current spin index in both cases:
C     ksp = isp  in the collinear case
C         = ispc in the noncollinear case
C     whereas ispc=1 for independent spins, and spin index when nspc=2
      do  ispc = 1, nspc
        ksp = max(ispc,isp)

        if (mode .eq. 1) then
          if (nlma .gt. nlmxx) call rxi('makusq:  nlmxx.lt.nlma=',nlma)
          do  ilm = 1, nlma
            k = ll(ilm)+1
            dlphi  = ppnl(3,k,ksp)/rmt
            dlphip = ppnl(4,k,ksp)/rmt
            phi    = ppnl(5,k,ksp)
            phip   = ppnl(6,k,ksp)
            dphi   = phi*dlphi/rmt
            dphip  = dlphip/rmt*phip
            det    = phi*dphip - dphi*phip
            rotp(ilm,1,1) = dphip/det
            rotp(ilm,1,2) = -dphi/det
            rotp(ilm,2,1) = -phip/det
            rotp(ilm,2,2) = phi/det
          enddo
        endif

C --- Loop over eigenstates ---
        allocate(a_zv((kmax+1)*nlma))

        do  ivec = 1, nev
C       call pusq3(ndimh,nlma,kmax,evec(1,ispc,ivec),w(ob),w(oa))
          call rlocb1 ( ndimh , nlma , kmax , evec ( 1 , ispc , ivec ) 
     .    , b_zv , a_zv )

C       call zprm('cPkL',2,w(oa),kmax+1,kmax+1,nlma)
          call pusq2 ( mode , ia , nkape , kmax , lmxa , lmxh , nlmto , 
     .    min ( nlma , nlmax ) , iprmb , a_zv , rotp , evec ( 1 , ispc 
     .    , ivec ) , vh , dh , vp , dp , au ( 1 , ivec , 1 , ksp ) , as 
     .    ( 1 , ivec , 1 , ksp ) , az ( 1 , ivec , 1 , ksp ) )

        enddo
        if (allocated(a_zv)) deallocate(a_zv)


C  --- Debugging --
C      print *, 'ksp,ib,nev=',ksp,ib,nev
C      call zprm('val',2,au(1,1,1,ksp),nlmax,(lmxa+1)**2,nev)
C      call zprm('slo',2,as(1,1,1,ksp),nlmax,(lmxa+1)**2,nev)
C      call zprm('loc',2,az(1,1,3,ksp),nlmax,(lmxa+1)**2,nev)

C ... end loop over noncollinear spins
      enddo

      if (allocated(b_zv)) deallocate(b_zv)

      call tcx('pusq1')

      end subroutine pusq1


      subroutine pusq2(mode,ia,nkape,kmax,lmxa,lmxh,nlmto,nlma,iprmb,
     .cPkL,r,evec,vh,dh,vp,dp,au,as,az)

C- Extract projection of eigenstate onto (u,s,z) for sphere at site ia
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 generate coefficients to values, slopes
Ci         :1 generate coefficients to phi,phidot
Ci   ia    :augmentation sphere
Ci   nkape :number of envelope function types which are joined to (u,s)
Ci         :Any ktab > nkape is a local orbital
Ci   kmax  :polynomial cutoff in P_kL expansion of envelope tails
Ci   lmxa  :augmentation l-cutoff
Ci   lmxh  :basis l-cutoff
Ci   nlmto :dimension of lmto component of basis
Ci   nlma  :number of L's in augmentation sphere = (lmxa+1)**2
Ci   iprmb :permutations ordering orbitals in l+i+h blocks (makidx.f)
Co   cPkL  :coefficients to P_kL expansion of evec
Ci   r     :2x2 rotation matrices rotating (phi,phidot) to (u,s)
Ci   evec  :eigenvector
Ci   vh    :value of head function in sphere ia
Ci   dh    :slope of head function in sphere ia
Ci   vp    :value of PkL expansion of tail function in sphere ia
Ci   dp    :slope of PkL expansion of tail function in sphere ia
Co Outputs
Co   au    :projection of this evec onto u function; see potpus.f
Co         :If mode=1, au = projection of this evec onto phi function
Co   as    :projection of this evec onto s function; see potpus.f
Co         :If mode=1, au = projection of this evec onto phidot function
Co   az    :projection of this evec onto local orbitals; see potpus.f
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   12 Feb 02 Extended to local orbitals
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,ia,nkape,kmax,lmxa,lmxh,nlmto,nlma,iprmb(1)
      double precision vh(0:lmxh,1),dh(0:lmxh,1)
      double precision vp(0:lmxa,0:kmax),dp(0:lmxa,0:kmax)
      integer nlmxx
      parameter (nlmxx=121)
      double precision r(nlmxx,2,2)
      double complex au(nlma),as(nlma),az(nlma),
     .evec(nlmto),cPkL(0:kmax,nlma)
C Local
      integer n0,nkap0,norb
      parameter (n0=10,nkap0=3)
      integer ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0)
      integer blks(n0*nkap0),ntab(n0*nkap0)
      integer io1,l1,ik1,nlm11,nlm12,ilm1,i1,ilma,k
      integer l,ll
      double precision xx
      double complex wk(nlmxx,2)

C     call tcn('pusq2')
      if (nlmto .eq. 0) return
      if (nlma .gt. nlmxx) call rxi('makusq:  nlmxx.lt.nlma=',nlma)

C --- Loop over all orbitals centered at this site ---
      call orbl(ia,0,nlmto,iprmb,norb,ltab,ktab,xx,offl,xx)
C     Block into groups of consecutive l
      call gtbsl1(4,norb,ltab,ktab,xx,xx,ntab,blks)

C     Contribution from head part
      do  io1 = 1, norb
        l1  = ltab(io1)
        ik1 = ktab(io1)
        nlm11 = l1**2+1
        nlm12 = nlm11 + blks(io1)-1
C       i1 = hamiltonian offset for first orbital in block
        i1 = offl(io1)-nlm11+1
        if (ik1 .le. nkape) then
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

C     Contribution from tail part
      do  ilma = 1, nlma
        l = ll(ilma)
        do  k = 0, kmax
          au(ilma) = au(ilma) + vp(l,k) * cPkL(k,ilma)
          as(ilma) = as(ilma) + dp(l,k) * cPkL(k,ilma)
        enddo
      enddo

C     Rotate to (phi,phidot)
      if (mode .ne. 0) then
        call dcopy(2*nlma,au,1,wk(1,1),1)
        call dcopy(2*nlma,as,1,wk(1,2),1)
        do  ilma = 1, nlma
          au(ilma) = wk(ilma,1)*r(ilma,1,1) + wk(ilma,2)*r(ilma,2,1)
          as(ilma) = wk(ilma,1)*r(ilma,1,2) + wk(ilma,2)*r(ilma,2,2)
        enddo
      endif

C     call tcx('pusq2')
      end subroutine pusq2


C      subroutine pusq3(ndimh,nlma,kmax,evec,b,a)
CC- Add together coeffs to expand wavefct at this site
C see rlocb1
CCu  Adapted from pvrlakl in nfp, reversing index order in b
C      implicit none
CC ... Passed parameters
C      integer ndimh,nlma,kmax
CC      double complex b(ndimh,nlma,0:kmax),a(0:kmax,nlma),evec(ndimh)
C      double complex b(0:kmax,nlma,ndimh),a(0:kmax,nlma),evec(ndimh)
CC ... Local parameters
C      integer k,ilma,i
C      call dpzero(a, 2*(kmax+1)*nlma)
CC     call tcn('pusq3')
C      do k = 0, kmax
C        do ilma = 1 ,nlma
C          do i = 1, ndimh
CC            a(k,ilma) = a(k,ilma) + evec(i)*b(i,ilma,k)
C            a(k,ilma) = a(k,ilma) + evec(i)*b(k,ilma,i)
C          enddo
C        enddo
C      enddo
CC     call tcx('pusq3')
C      end

