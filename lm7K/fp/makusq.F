Cdelw1 var ['orofi', '1', 'orofi', 'nr', 'real(8)', 'rv_w_', 'makusq']
Cdelw1 var ['ofh', '1', 'ofh', 'nr*(lmxh+1)*nkapi', 'real(8)', 'rv_w_', 'makusq']
Cdelw1 var ['oxh', '1', 'oxh', 'nr*(lmxh+1)*nkapi', 'real(8)', 'rv_w_', 'makusq']
Cdelw1 var ['ovh', '1', 'ovh', '(lmxh+1)*nkapi', 'real(8)', 'rv_w_', 'makusq']
Cdelw1 var ['odh', '1', 'odh', '(lmxh+1)*nkapi', 'real(8)', 'rv_w_', 'makusq']
Cdelw1 var ['ofp', '1', 'ofp', 'nr*(lmxa+1)*(kmax+1)', 'real(8)', 'rv_w_', 'makusq']
Cdelw1 var ['oxp', '1', 'oxp', 'nr*(lmxa+1)*(kmax+1)', 'real(8)', 'rv_w_', 'makusq']
Cdelw1 var ['ovp', '1', 'ovp', '(lmxa+1)*(kmax+1)', 'real(8)', 'rv_w_', 'makusq']
Cdelw1 var ['odp', '1', 'odp', '(lmxa+1)*(kmax+1)', 'real(8)', 'rv_w_', 'makusq']
Cdelw1 var ['oiprmb', '-1', 'undef', '0', 'undef', 'x', 'makusq']
Cdelw1 var ['ob', '0', 'ob', '(kmax+1)*nlma*ndimh', 'complex(8)', 'zv_w_', 'pusq1']
Cdelw1 var ['ocg', '-1', 'undef', '0', 'undef', 'x', 'pusq1']
Cdelw1 var ['oidxcg', '-1', 'undef', '0', 'undef', 'x', 'pusq1']
Cdelw1 var ['ojcg', '-1', 'undef', '0', 'undef', 'x', 'pusq1']
Cdelw1 var ['ocy', '-1', 'undef', '0', 'undef', 'x', 'pusq1']
Cdelw1 var ['oa', '1', 'oa', '(kmax+1)*nlma', 'complex(8)', 'zv_w_', 'pusq1']
Cdelw1 undel [['undef', 'wref'], 'oiprmb', 'makusq']
Cdelw1 undel [['undef', 'wref'], 'ocg', 'pusq1']
Cdelw1 undel [['undef', 'wref'], 'oidxcg', 'pusq1']
Cdelw1 undel [['undef', 'wref'], 'ojcg', 'pusq1']
Cdelw1 undel [['undef', 'wref'], 'ocy', 'pusq1']
Cgetarg...info...           structure ['makusq', 'site', 'ssite']
Cgetarg...info...           structure ['makusq', 'ham', 'sham']
Cgetarg...info...           structure ['makusq', 'spec', 'sspec']
Cgetarg...info...           structure ['pusq1', 'lat', 'slat']
Cgetarg...info...           structure ['pusq1', 'site', 'ssite']
Cgetarg...info...           structure ['pusq1', 'spec', 'sspec']
Cgetarg...info...           use_to_add ['pusq1', 'm_struc_def']
Cgetarg...info...           use_to_add ['makusq', 'm_struc_def']
      subroutine makusq(mode,ssite,sspec,slat,sham,nbas,nsites,isite,
     .  nlmax,ndham,ndimh,napw,igvapw,nev,nsp,nspc,isp,iq,q,evec,ppnl,
     .  aus)
       
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
     .  ndham,ndimh,napw,igvapw(3,napw),nev,nlmax,n0,nppn,
     .  nsites,isite(nsites)
      parameter (n0=10,nppn=12)
Cgetarg       double precision ssite(1),sspec(1),slat(1),sham(1),q(3)
       real(8):: q(3) 
       type(s_site)::ssite(*)
       type(s_spec)::sspec(*)
       type(s_lat)::slat
       type(s_ham)::sham

      double precision ppnl(nppn,n0,nsp,nbas)
      double complex evec(ndimh,nsp,nev),
     .  aus(nlmax,ndham*nspc,3,nsp,nsites,iq)
C ... Local parameters
      integer nkap0
      parameter (nkap0=3)
      integer lh(nkap0)
      double precision eh(n0,nkap0),rsmh(n0,nkap0),rsma,a,rmt
      integer igetss,ib,nkapi,is,nr,kmax,lmxa,lmxl,lmxh,i
Cdelw1 do not change  [['undef', 'wref'], 'oiprmb', 'makusq']
Cdelw1       integer oiprmb,orofi,ofh,oxh,ovh,odh,ofp,oxp,ovp,odp
       integer:: oiprmb 
       real(8) ,allocatable :: rv_w_orofi(:)
       real(8) ,allocatable :: rv_w_ofh(:)
       real(8) ,allocatable :: rv_w_oxh(:)
       real(8) ,allocatable :: rv_w_ovh(:)
       real(8) ,allocatable :: rv_w_odh(:)
       real(8) ,allocatable :: rv_w_ofp(:)
       real(8) ,allocatable :: rv_w_oxp(:)
       real(8) ,allocatable :: rv_w_ovp(:)
       real(8) ,allocatable :: rv_w_odp(:)

C ... Heap
      integer w(1)
      common /w/ w

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
Cgetarg         is = igetss('site spec',ib,ssite)
         is = int(ssite(ib)%spec) 

Cgetarg...info...           integer :: oiprmb , integer(8) :: sham%oindxo 1
Cgetarg         call upack('ham oindxo',sham,oiprmb,0,0,0,0)
         
Cdelw1 warning(1) , probably  oiprmb  is not defined yet at linenumber= 140
Cdelw1 do not change oiprmb because of ['undef', 'wref']
         oiprmb=sham%oindxo

Cgetarg...info...           real(8) :: rsma , real(8) :: sspec%rsma 1
Cgetarg...info...           integer :: lmxa , integer(8) :: sspec%lmxa 1
Cgetarg...info...           integer :: lmxl , integer(8) :: sspec%lmxl 1
Cgetarg...info...           integer :: kmax , integer(8) :: sspec%kmxt 1
Cgetarg         call upack('spec rsma lmxa lmxl kmxt',sspec,is,rsma,lmxa,lmxl,
Cgetarg      .    kmax)
         
         rsma=sspec(is)%rsma
         lmxa=sspec(is)%lmxa
         lmxl=sspec(is)%lmxl
         kmax=sspec(is)%kmxt

Cgetarg...info...           real(8) :: a , real(8) :: sspec%a 1
Cgetarg...info...           integer :: nr , integer(8) :: sspec%nr 1
Cgetarg...info...           real(8) :: rmt , real(8) :: sspec%rmt 1
Cgetarg         call upack('spec a nr rmt',sspec,is,a,nr,rmt,0)
         
         a=sspec(is)%a
         nr=sspec(is)%nr
         rmt=sspec(is)%rmt

        call uspecb(0,1,sspec,is,is,lh,rsmh,eh,nkapi)
Cgetarg...info...           integer :: lmxh , integer(8) :: sspec%lmxb 1
Cgetarg         call upack('spec lmxb',sspec,is,lmxh,0,0,0)
         
         lmxh=sspec(is)%lmxb

        if (lmxa .eq. -1) goto 10

C   --- Set up all radial head and tail functions, and their BC's
Cdelw1         call defrr(orofi, nr)
         allocate(rv_w_orofi(nr))
         if (nr<0) rv_w_orofi(:)=0.0d0

Cdelw1         call radmsh(rmt,a,nr,w(orofi))
         call radmsh ( rmt , a , nr , rv_w_orofi ) 

Cdelw1         call defrr(ofh, nr*(lmxh+1)*nkapi)
         allocate(rv_w_ofh(nr*(lmxh+1)*nkapi))
         if (nr*(lmxh+1)*nkapi<0) rv_w_ofh(:)=0.0d0

Cdelw1         call defrr(oxh, nr*(lmxh+1)*nkapi)
         allocate(rv_w_oxh(nr*(lmxh+1)*nkapi))
         if (nr*(lmxh+1)*nkapi<0) rv_w_oxh(:)=0.0d0

Cdelw1         call defrr(ovh, (lmxh+1)*nkapi)
         allocate(rv_w_ovh((lmxh+1)*nkapi))
         if ((lmxh+1)*nkapi<0) rv_w_ovh(:)=0.0d0

Cdelw1         call defrr(odh, (lmxh+1)*nkapi)
         allocate(rv_w_odh((lmxh+1)*nkapi))
         if ((lmxh+1)*nkapi<0) rv_w_odh(:)=0.0d0

Cdelw1         call fradhd(nkapi,eh,rsmh,lh,lmxh,nr,w(orofi),w(ofh),
Cdelw1      .    w(oxh),w(ovh),w(odh))
         call fradhd ( nkapi , eh , rsmh , lh , lmxh , nr , rv_w_orofi 
     .   , rv_w_ofh , rv_w_oxh , rv_w_ovh , rv_w_odh ) 

Cdelw1         call defrr(ofp, nr*(lmxa+1)*(kmax+1))
         allocate(rv_w_ofp(nr*(lmxa+1)*(kmax+1)))
         if (nr*(lmxa+1)*(kmax+1)<0) rv_w_ofp(:)=0.0d0

Cdelw1         call defrr(oxp, nr*(lmxa+1)*(kmax+1))
         allocate(rv_w_oxp(nr*(lmxa+1)*(kmax+1)))
         if (nr*(lmxa+1)*(kmax+1)<0) rv_w_oxp(:)=0.0d0

Cdelw1         call defrr(ovp, (lmxa+1)*(kmax+1))
         allocate(rv_w_ovp((lmxa+1)*(kmax+1)))
         if ((lmxa+1)*(kmax+1)<0) rv_w_ovp(:)=0.0d0

Cdelw1         call defrr(odp, (lmxa+1)*(kmax+1))
         allocate(rv_w_odp((lmxa+1)*(kmax+1)))
         if ((lmxa+1)*(kmax+1)<0) rv_w_odp(:)=0.0d0

Cdelw1         call fradpk(kmax,rsma,lmxa,nr,w(orofi),w(ofp),w(oxp),w(ovp),
Cdelw1      .    w(odp))
         call fradpk ( kmax , rsma , lmxa , nr , rv_w_orofi , rv_w_ofp 
     .   , rv_w_oxp , rv_w_ovp , rv_w_odp ) 

C   --- Add to the coefficient for the projection onto (u,s) for this site
Cdelw1 warning(1) , probably  oiprmb  is not defined yet at linenumber= 187
Cdelw1 do not change oiprmb because of ['undef', 'wref']
Cdelw1 warning(1) , probably  oiprmb  is not defined yet at linenumber= 187
Cdelw1 do not change oiprmb because of ['undef', 'wref']
Cdelw1         call pusq1(mode,ib,isp,nspc,w(oiprmb),nlmax,lmxh,nbas,ssite,
Cdelw1      .    sspec,slat,q,ndham,ndimh,napw,igvapw,nev,evec,w(ovh),w(odh),
Cdelw1      .    w(ovp),w(odp),ppnl(1,1,1,ib),
Cdelw1      .    aus(1,1,1,1,i,iq),aus(1,1,2,1,i,iq),aus(1,1,3,1,i,iq))
         call pusq1 ( mode , ib , isp , nspc , w ( oiprmb ) , nlmax , 
     .   lmxh , nbas , ssite , sspec , slat , q , ndham , ndimh , napw 
     .   , igvapw , nev , evec , rv_w_ovh , rv_w_odh , rv_w_ovp , rv_w_odp 
     .   , ppnl ( 1 , 1 , 1 , ib ) , aus ( 1 , 1 , 1 , 1 , i , iq ) , 
     .   aus ( 1 , 1 , 2 , 1 , i , iq ) , aus ( 1 , 1 , 3 , 1 , i , iq 
     .   ) ) 

Cdelw1 rlse name= orofi old_list= orofi ofh oxh ovh odh ofp oxp ovp odp 
Cdelw1 rlse name= orofi new_list= (None)
Cdelw1         call rlse(orofi)
         if (allocated(rv_w_odp)) deallocate(rv_w_odp)
         if (allocated(rv_w_ovp)) deallocate(rv_w_ovp)
         if (allocated(rv_w_oxp)) deallocate(rv_w_oxp)
         if (allocated(rv_w_ofp)) deallocate(rv_w_ofp)
         if (allocated(rv_w_odh)) deallocate(rv_w_odh)
         if (allocated(rv_w_ovh)) deallocate(rv_w_ovh)
         if (allocated(rv_w_oxh)) deallocate(rv_w_oxh)
         if (allocated(rv_w_ofh)) deallocate(rv_w_ofh)
         if (allocated(rv_w_orofi)) deallocate(rv_w_orofi)

C   --- print
C        print *, 'isp,ib,nev,iq=',isp,ib,nev,iq
C        call zprm('val',2,aus(1,1,1,isp,ib,iq),nlmax,(lmxa+1)**2,nev)
C        call zprm('slo',2,aus(1,1,2,isp,ib,iq),nlmax,(lmxa+1)**2,nev)
C        call zprm('loc',2,aus(1,1,3,isp,ib,iq),nlmax,(lmxa+1)**2,nev)
   10   continue
      enddo
      call tcx('makusq')
Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [oiprmb]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine makusq 


      subroutine pusq1(mode,ia,isp,nspc,iprmb,nlmax,lmxh,nbas,ssite,
     .  sspec,slat,q,ndham,ndimh,napw,igvapw,nev,evec,vh,dh,vp,dp,ppnl,
     .  au,as,az)
       
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
     .  nbas,ndham,ndimh,napw,igvapw(3,napw),nev,nlmbx,n0,nppn
      parameter (nlmbx=25, n0=10, nppn=12)
      integer iprmb(ndimh)
      double precision ppnl(nppn,n0,2)
      double precision vp(*),dp(*),vh(*),dh(*)
Cgetarg       double precision q(3),ssite(1),sspec(1),slat(1)
       real(8):: q(3) 
       type(s_site)::ssite(*)
       type(s_spec)::sspec(*)
       type(s_lat)::slat

      double complex evec(ndimh,nspc,ndimh),
     .  au(nlmax,ndham*nspc,3,2),
     .  as(nlmax,ndham*nspc,3,2),
     .  az(nlmax,ndham*nspc,3,2)
C ... Local parameters
      integer nkap0,nlmxx
      parameter (nkap0=3,nlmxx=121)
      integer lh(nkap0)
Cdelw1 do not change  [['undef', 'wref'], 'ocg', 'pusq1']
Cdelw1 do not change  [['undef', 'wref'], 'ojcg', 'pusq1']
Cdelw1 do not change  [['undef', 'wref'], 'oidxcg', 'pusq1']
Cdelw1 do not change  [['undef', 'wref'], 'ocy', 'pusq1']
Cdelw1       integer ocg,ojcg,oidxcg,ocy,ob,oa
       integer:: ocg , ojcg , oidxcg , ocy 
       complex(8) ,allocatable :: zv_w_ob(:)
       complex(8) ,allocatable :: zv_w_oa(:)

      integer isa,lmxa,lmxha,kmax,nlma,ivec,
     .  ilm,k,ll,nkape,ksp,ispc,nlmto
      double precision eh(n0,nkap0),rsmh(n0,nkap0)
      double precision rsma,pa(3),rmt,
     .  phi,phip,dphi,dlphi,dphip,dlphip,det,rotp(nlmxx,2,2)
C ... Heap
      integer w(1)
      common /w/ w

      call tcn ('pusq1')
Cgetarg...info...           integer :: ocg , integer(8) :: slat%ocg 1
Cgetarg...info...           integer :: ojcg , integer(8) :: slat%ojcg 1
Cgetarg...info...           integer :: oidxcg , integer(8) :: slat%oidxcg 1
Cgetarg...info...           integer :: ocy , integer(8) :: slat%ocy 1
Cgetarg       call upack('lat ocg ojcg oidxcg ocy',slat,ocg,ojcg,oidxcg,ocy,0)
       
Cdelw1 warning(1) , probably  ocg  is not defined yet at linenumber= 302
Cdelw1 do not change ocg because of ['undef', 'wref']
       ocg=slat%ocg
Cdelw1 warning(1) , probably  ojcg  is not defined yet at linenumber= 303
Cdelw1 do not change ojcg because of ['undef', 'wref']
       ojcg=slat%ojcg
Cdelw1 warning(1) , probably  oidxcg  is not defined yet at linenumber= 304
Cdelw1 do not change oidxcg because of ['undef', 'wref']
       oidxcg=slat%oidxcg
Cdelw1 warning(1) , probably  ocy  is not defined yet at linenumber= 305
Cdelw1 do not change ocy because of ['undef', 'wref']
       ocy=slat%ocy

Cgetarg...info...           integer :: isa , integer(8) :: ssite%spec 1
Cgetarg...info...           real(8) :: pa(3) , real(8) :: ssite%pos 3
Cgetarg       call upack('site spec pos',ssite,ia,isa,pa,0,0)
       
       isa=ssite(ia)%spec
       i_copy_size=size(ssite(ia)%pos) 
       call dcopy(i_copy_size,ssite(ia)%pos,1,pa,1) 

Cgetarg...info...           integer :: lmxa , integer(8) :: sspec%lmxa 1
Cgetarg...info...           integer :: lmxha , integer(8) :: sspec%lmxb 1
Cgetarg...info...           integer :: kmax , integer(8) :: sspec%kmxt 1
Cgetarg...info...           real(8) :: rsma , real(8) :: sspec%rsma 1
Cgetarg       call upack('spec lmxa lmxb kmxt rsma',sspec,isa,lmxa,lmxha,kmax,
Cgetarg      .            rsma)
       
       lmxa=sspec(isa)%lmxa
       lmxha=sspec(isa)%lmxb
       kmax=sspec(isa)%kmxt
       rsma=sspec(isa)%rsma

Cgetarg...info...           real(8) :: rmt , real(8) :: sspec%rmt 1
Cgetarg       call upack2('spec rmt',sspec,isa,rmt)
       
       rmt=sspec(isa)%rmt

Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [ocg oidxcg ojcg ocy]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
      if (lmxa .eq. -1) return
      nlmto = ndimh-napw

C     nlmha = (lmxha+1)**2
      nlma  = (lmxa+1)**2
C     if (nlma .gt. nlmax) call rxi('makusq: need nlmax',nlma)
C     Count no. envelope functions connecting (phi,phidot) at site ia
      call uspecb(0,1,sspec,isa,isa,lh,rsmh,eh,nkape)

C --- Make strux to expand all orbitals at site ia ---
Cdelw1       call defcc (ob,   (kmax+1)*nlma*ndimh)
       allocate(zv_w_ob((kmax+1)*nlma*ndimh))
       if ((kmax+1)*nlma*ndimh<0) zv_w_ob(:)=0.0d0

Cdelw1 warning(1) , probably  ocg  is not defined yet at linenumber= 343
Cdelw1 do not change ocg because of ['undef', 'wref']
Cdelw1 warning(1) , probably  oidxcg  is not defined yet at linenumber= 343
Cdelw1 do not change oidxcg because of ['undef', 'wref']
Cdelw1 warning(1) , probably  ojcg  is not defined yet at linenumber= 343
Cdelw1 do not change ojcg because of ['undef', 'wref']
Cdelw1 warning(1) , probably  ocy  is not defined yet at linenumber= 343
Cdelw1 do not change ocy because of ['undef', 'wref']
Cdelw1 warning(1) , probably  ocg  is not defined yet at linenumber= 343
Cdelw1 do not change ocg because of ['undef', 'wref']
Cdelw1 warning(1) , probably  oidxcg  is not defined yet at linenumber= 343
Cdelw1 do not change oidxcg because of ['undef', 'wref']
Cdelw1 warning(1) , probably  ojcg  is not defined yet at linenumber= 343
Cdelw1 do not change ojcg because of ['undef', 'wref']
Cdelw1 warning(1) , probably  ocy  is not defined yet at linenumber= 343
Cdelw1 do not change ocy because of ['undef', 'wref']
Cdelw1       call bstrux(2,slat,ssite,sspec,w(ocg),w(oidxcg),w(ojcg),w(ocy),
Cdelw1      .  iprmb,nbas,ia,pa,rsma,q,kmax,nlma,ndimh,napw,igvapw,w(ob),w)
       call bstrux ( 2 , slat , ssite , sspec , w ( ocg ) , w ( oidxcg 
     . ) , w ( ojcg ) , w ( ocy ) , iprmb , nbas , ia , pa , rsma , 
     . q , kmax , nlma , ndimh , napw , igvapw , zv_w_ob , w ) 


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
Cdelw1       call defcc (oa,   (kmax+1)*nlma)
       allocate(zv_w_oa((kmax+1)*nlma))
       if ((kmax+1)*nlma<0) zv_w_oa(:)=0.0d0

      do  ivec = 1, nev
C       call pusq3(ndimh,nlma,kmax,evec(1,ispc,ivec),w(ob),w(oa))
Cdelw1         call rlocb1(ndimh,nlma,kmax,evec(1,ispc,ivec),w(ob),w(oa))
         call rlocb1 ( ndimh , nlma , kmax , evec ( 1 , ispc , ivec ) 
     .   , zv_w_ob , zv_w_oa ) 

C       call zprm('cPkL',2,w(oa),kmax+1,kmax+1,nlma)
Cdelw1         call pusq2(mode,ia,nkape,kmax,lmxa,lmxh,nlmto,min(nlma,nlmax),
Cdelw1      .    iprmb,w(oa),rotp,evec(1,ispc,ivec),vh,dh,vp,dp,
Cdelw1      .    au(1,ivec,1,ksp),as(1,ivec,1,ksp),az(1,ivec,1,ksp))
         call pusq2 ( mode , ia , nkape , kmax , lmxa , lmxh , nlmto , 
     .   min ( nlma , nlmax ) , iprmb , zv_w_oa , rotp , evec ( 1 , ispc 
     .   , ivec ) , vh , dh , vp , dp , au ( 1 , ivec , 1 , ksp ) , as 
     .   ( 1 , ivec , 1 , ksp ) , az ( 1 , ivec , 1 , ksp ) ) 

      enddo
Cdelw1 rlse name= oa old_list= ob oa 
Cdelw1 rlse name= oa new_list= ob 
Cdelw1       call rlse(oa)
       if (allocated(zv_w_oa)) deallocate(zv_w_oa)


C  --- Debugging --
C      print *, 'ksp,ib,nev=',ksp,ib,nev
C      call zprm('val',2,au(1,1,1,ksp),nlmax,(lmxa+1)**2,nev)
C      call zprm('slo',2,as(1,1,1,ksp),nlmax,(lmxa+1)**2,nev)
C      call zprm('loc',2,az(1,1,3,ksp),nlmax,(lmxa+1)**2,nev)

C ... end loop over noncollinear spins
      enddo

Cdelw1 rlse name= ob old_list= ob 
Cdelw1 rlse name= ob new_list= (None)
Cdelw1       call rlse(ob)
       if (allocated(zv_w_ob)) deallocate(zv_w_ob)

      call tcx('pusq1')

Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlistundel: [ocg oidxcg ojcg ocy]
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
       end subroutine pusq1 


      subroutine pusq2(mode,ia,nkape,kmax,lmxa,lmxh,nlmto,nlma,iprmb,
     .  cPkL,r,evec,vh,dh,vp,dp,au,as,az)

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
     .  evec(nlmto),cPkL(0:kmax,nlma)
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
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
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
Cgetarg       end
Cdelw1 w_varlist remains: (None)
Cdelw1 w_varlist (undel), remains: (None)
Cdelw1 w_varlist (del), remains: (None)
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

