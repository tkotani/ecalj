#define USEOPTIONAL
      module m_mkrou3
      contains
      subroutine mkrou3(mode,lmxa,nlml,nsp,pnz,dmatl,hab,sab,qsum,hsum)

C- l-decomposed charges and eigenvalue sum
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 do nothing
Ci         :1 make qsum only
Ci         :2 make qsum and hsum both
Ci   lmxa  :augmentation l-cutoff
Ci   dmatl :dmatl(l1,l2,mlm,i,j,isp) holds coefficients to a y_lm
Ci         :expansion of the function products f_i(l1) f_j(l2)
Ci         :where f_i,f_j form this table.
Ci         :    (uu  us  uz)
Ci         :    (su  ss  sz)
Ci         :    (zu  zs  zz)
Ci   hab   :<u,s | h | u,s> for each pair uu, us, su, ss; see remarks
Ci   sab   :<u,s | 1 | u,s>
Co Outputs
Co   qsum  :l-decomposed sphere charge
Co   hsum  :l-decomposed one-electron energy
Cr Remarks
Cr   qsum and hsum are used to find the band centers of gravity
Cr   u and s ae linear combinations radial wave functions defined as:
Cr   u has val=1, slo=1 at rmax, s has val=0, slo=1
Cu Updates
Cu   28 Aug 01 Extended to local orbitals.
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,lmxa,n0,nab,nlml,nsp
      parameter (n0=10,nab=9)
#ifdef USEOPTIONAL
      real(8),optional:: hab(nab,n0,nsp), hsum(n0,nsp)
      real(8):: qsum(n0,nsp) ,sab(nab,n0,nsp),
     .dmatl(0:lmxa,0:lmxa,nlml,3,3,nsp),pnz(n0,2)
#else
      double precision hab(nab,n0,nsp),sab(nab,n0,nsp),
     .qsum(n0,nsp),hsum(n0,nsp),
     .dmatl(0:lmxa,0:lmxa,nlml,3,3,nsp),pnz(n0,2)
#endif

C ... Local parameters
      integer l,isp,m
      double precision pi,srfpi,qz,hz

      if (mode .eq. 0) return
C     call tcn('mkrou3')

      pi = 4d0*datan(1d0)
      srfpi = dsqrt(4d0*pi)

      do  isp = 1, nsp
        do  l = 0, lmxa
          if (mode .ge. 1) qsum(l+1,isp) = 0d0
          if (mode .ge. 2) hsum(l+1,isp) = 0d0
        enddo
      enddo

      do  isp = 1, nsp
        do  l = 0, lmxa
          m = l+1
          qsum(m,isp) =
     .    + dmatl(l,l,1,1,1,isp)*sab(1,m,isp)*srfpi
     .    + dmatl(l,l,1,1,2,isp)*sab(2,m,isp)*srfpi
     .    + dmatl(l,l,1,2,1,isp)*sab(3,m,isp)*srfpi
     .    + dmatl(l,l,1,2,2,isp)*sab(4,m,isp)*srfpi
          if (mode .ge. 2) then
            hsum(m,isp) =
     .      + dmatl(l,l,1,1,1,isp)*hab(1,m,isp)*srfpi
     .      + dmatl(l,l,1,1,2,isp)*hab(2,m,isp)*srfpi
     .      + dmatl(l,l,1,2,1,isp)*hab(3,m,isp)*srfpi
     .      + dmatl(l,l,1,2,2,isp)*hab(4,m,isp)*srfpi
          endif
C         ... uz, sz, zu, zs, zz terms
          if (pnz(m,1) .ne. 0) then
            qz =
     .      + dmatl(l,l,1,1,3,isp)*sab(5,m,isp)*srfpi
     .      + dmatl(l,l,1,2,3,isp)*sab(6,m,isp)*srfpi
     .      + dmatl(l,l,1,3,1,isp)*sab(5,m,isp)*srfpi
     .      + dmatl(l,l,1,3,2,isp)*sab(6,m,isp)*srfpi
     .      + dmatl(l,l,1,3,3,isp)*sab(7,m,isp)*srfpi
            if (mode .ge. 2) then
              hz =
     .        + dmatl(l,l,1,1,3,isp)*hab(5,m,isp)*srfpi
     .        + dmatl(l,l,1,2,3,isp)*hab(6,m,isp)*srfpi
     .        + dmatl(l,l,1,3,1,isp)*hab(5,m,isp)*srfpi
     .        + dmatl(l,l,1,3,2,isp)*hab(6,m,isp)*srfpi
     .        + dmatl(l,l,1,3,3,isp)*hab(7,m,isp)*srfpi
            endif
          endif

        enddo
      enddo

C     call tcx('mkrou3')

      end subroutine mkrou3
      end module m_mkrou3


         subroutine mkrout ( ssite , sspec , slat , sham , nbas , nsp 
     .   , nlmto , lekkl , sv_p_oqkkl , sv_p_oeqkkl , sv_p_orhoat , hab 
     .   , sab , qbyl , hbyl , lrout ) 




      use m_struc_def  !Cgetarg
      use m_globalvariables
      use m_mkrou3
C- Assembles local output densities out of the qkkl, and core states
C ----------------------------------------------------------------------
Ci Inputs
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec pnu ov0 ov1 pz
Ci     Stored:    *
Ci     Passed to: *
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: rsma lmxa lmxl kmxt a nr rmt lmxb stc orhoc
Ci     Stored:    *
Ci     Passed to: corprm uspecb gtpcor
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: ocg ojcg oidxcg ocy
Ci     Stored:    *
Ci     Passed to: *
Ci   sham  :struct for parameters defining hamiltonian; see routine uham
Ci     Elts read: oindxo eterms
Ci     Stored:    eterms (see Outputs)
Ci     Passed to: *
Ci   nbas  :size of basis
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nlmto :dimension of lower block of LMTO part of hamiltonian
Cl   lekkl :1 Find band CG from eqkkl/qkkl
Ci   oqkkl :local density-matrix (rhocbl or comparable routine)
Co   oeqkkl:local part of energy-weighted density matrix
Ci   hab   :hamiltonian matrix elements of radial wave functions at each
Ci         :site.  See potpus for their definition.
Ci   sab   :overlap matrix elements of radial wave functions at each
Ci         :site.  See potpus for their definition.
Ci         :hab and sab are used here to find band cg
Ci   lrout :0 calculate core states part only
Co Outputs
Co   orhoat:vector of offsets containing site density
Co   sham->eterms various integrals for the total energy are stored:
Co         :(8)  sumec = sum of foca=0 core eigenvalues
Co         :(9)  sumtc = sum of all core kinetic energies
Co         :(12) sumt0 = sum of frozen core kinetic energies
Co   qbyl  :l-decomposed charge
Co   hbyl  :l-decomposed eigenvalue sum
Cr Remarks
Cr   u and s are linear combinations of and phi,phidot defined as:
Cr   u has val=1, slo=1 at rmax, s has val=0, slo=1
Cu Updates
Cu   30 Jul 08 (T. Kotani) Use ekkl to set hbyl
Cu   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
Cu   10 Apr 02 Redimensionsed eh,rsmh to accomodate larger lmax
Cu   28 Aug 01 Extended to local orbitals.
Cu   18 Jun 00 spin polarized
Cu   20 May 00 adapted from nfp mk_dens2
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
         integer:: n0 , nab , nbas , nsp , lekkl , nlmto , lrout 
         type(s_rv1) :: sv_p_orhoat(3,nbas)

        type(s_rv1) :: sv_p_oeqkkl(3,nbas)

       type(s_rv1) :: sv_p_oqkkl(3,nbas)

      parameter (n0=10,nab=9)
      real(8):: qbyl(n0,nsp,nbas) , hbyl(n0,nsp,nbas) , sab(nab,n0,nsp,nbas) 
     ., hab(nab,n0,nsp,nbas)
cki      type(s_site)::ssite(*)
      type(s_site)::ssite(nbas)
cki      type(s_spec)::sspec(*)
      type(s_spec)::sspec(nbas)
      type(s_ham)::sham
      type(s_lat)::slat

C ... Local parameters
                 integer:: ib , ipr , iprint , is , k , kcor , kmax , lcor , lfoc 
     .           , lgunit , lmxa , lmxh , lmxl , nlma , nlmh , nlml , nlml1 , 
     .           nr , ncore , stdo , igetss 
                 real(8),pointer :: rv_p_oqpp(:) =>NULL()

                real(8),pointer :: rv_p_oqhp(:) =>NULL()

               real(8),pointer :: rv_p_oqhh(:) =>NULL()

              real(8),pointer :: rv_p_orhoca(:) =>NULL()

             real(8),pointer :: rv_p_ov01(:) =>NULL()

            real(8),pointer :: rv_p_ov0(:) =>NULL()

           integer,pointer :: iv_p_ojcg(:) =>NULL()

          integer,pointer :: iv_p_oidxcg(:) =>NULL()

         real(8),pointer :: rv_p_ocy(:) =>NULL()

        real(8),pointer :: rv_p_ocg(:) =>NULL()

       integer,pointer :: iv_p_oiprmb(:) =>NULL()

      real(8) ,allocatable :: dmatl_rv(:)
      real(8) ,allocatable :: dh_rv(:)
      real(8) ,allocatable :: dp_rv(:)
      real(8) ,allocatable :: fh_rv(:)
      real(8) ,allocatable :: fp_rv(:)
      real(8) ,allocatable :: rofi_rv(:)
      real(8) ,allocatable :: rss_rv(:)
      real(8) ,allocatable :: rus_rv(:)
      real(8) ,allocatable :: ruu_rv(:)
      real(8) ,allocatable :: rwgt_rv(:)
      real(8) ,allocatable :: sl_rv(:)
      real(8) ,allocatable :: ul_rv(:)
      real(8) ,allocatable :: gz_rv(:)
      real(8) ,allocatable :: vh_rv(:)
      real(8) ,allocatable :: vp_rv(:)
      real(8) ,allocatable :: xh_rv(:)
      real(8) ,allocatable :: xp_rv(:)
      real(8) ,allocatable :: chh_rv(:)
      real(8) ,allocatable :: chp_rv(:)
      real(8) ,allocatable :: cpp_rv(:)

      double precision a,ceh,pi,qcor(2),rfoc,rmt,rsma,smec,smtc,stc0,
     .sum1,sum2,sums1,sums2,sumec,sumtc,sumt0,xx,y0,z,ddot,res,
     .rsml(n0),ehl(n0)
C ... For basis
      integer nkap0,nkapi,nkape,nkaph,nglob
      parameter (nkap0=3)
      integer lh(nkap0)
      double precision eh(n0,nkap0),rsmh(n0,nkap0)
      integer ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0)
      integer norb,ntab(n0*nkap0),blks(n0*nkap0)
      double precision pnu(n0,2),pnz(n0,2),qcorg,qcorh,qsc,cofg,cofh
C     this ordering must match sham->eterms; see uham
      double precision eterms(20)
      equivalence (eterms(8),sumec)
      equivalence (eterms(9),sumtc)
      equivalence (eterms(12),sumt0)
C ... Heap
       integer ::iwdummy 


C ... Setup
      stdo = lgunit(1)
      ipr  = iprint()
      pi   = 4d0*datan(1d0)
      y0   = 1d0/dsqrt(4d0*pi)
Changenglob      nkaph = nglob('nkaph')
      nkaph = globalvariables%nkaph
      sums1 = 0
      call tcn('mkrout')

       rv_p_ocg => slat%rv_p_ocg 

       iv_p_ojcg => slat%iv_p_ojcg 

       iv_p_oidxcg => slat%iv_p_oidxcg 

       rv_p_ocy => slat%rv_p_ocy 



       iv_p_oiprmb => sham%iv_p_oindxo 

      i_copy_size=size(sham%eterms)
      call dcopy(i_copy_size,sham%eterms,1,eterms,1)


C --- Start loop over atoms ---
      sumtc = 0d0
      sumec = 0d0
      sumt0 = 0d0
      if (ipr.ge.30 .and. lrout.gt.0) then
        if (nsp .eq. 1) write(stdo,201)
        if (nsp .eq. 2) write(stdo,202)
  201   format(/' mkrout:  Qtrue      sm,loc       local')
  202   format(/' mkrout:  Qtrue      sm,loc       local',8x,
     .  'true mm   smooth mm    local mm')
      endif
      do  ib = 1, nbas
        is = int(ssite(ib)%spec)


        rsma=sspec(is)%rsma
        lmxa=sspec(is)%lmxa
        lmxl=sspec(is)%lmxl
        kmax=sspec(is)%kmxt

        if (lmxa .eq. -1) goto 10


        a=sspec(is)%a
        nr=sspec(is)%nr
        rmt=sspec(is)%rmt


        lmxh=sspec(is)%lmxb
        stc0=sspec(is)%stc
         rv_p_orhoca => sspec(is)%rv_p_orhoc 


        call corprm(sspec,is,qcorg,qcorh,qsc,cofg,cofh,ceh,lfoc,rfoc,z)

        call uspecb(0,1,sspec,is,is,lh,rsmh,eh,nkapi)
        call uspecb(0,2,sspec,is,is,lh,rsmh,eh,nkape)
         call orbl ( ib , 0 , nlmto , iv_p_oiprmb , norb , ltab , ktab 
     .   , xx , offl , xx ) 

        call gtbsl1(4,norb,ltab,ktab,xx,xx,ntab,blks)


        is=ssite(ib)%spec
        i_copy_size=size(ssite(ib)%pnu)
        call dcopy(i_copy_size,ssite(ib)%pnu,1,pnu,1)
         rv_p_ov0 => ssite(ib)%rv_p_ov0 

         rv_p_ov01 => ssite(ib)%rv_p_ov1 



        i_copy_size=size(ssite(ib)%pz)
        call dcopy(i_copy_size,ssite(ib)%pz,1,pnz,1)

        call gtpcor(sspec,is,kcor,lcor,qcor)
        nlml = (lmxl+1)**2
        nlma = (lmxa+1)**2
        nlmh = (lmxh+1)**2
        allocate(rofi_rv(nr))

        allocate(rwgt_rv(nr))

        call radmsh ( rmt , a , nr , rofi_rv )

        call radwgt ( rmt , a , nr , rwgt_rv )

        if (lrout .ne. 0) then
          rv_p_oqhh => sv_p_oqkkl(3,ib)%v
          rv_p_oqhp => sv_p_oqkkl(2,ib)%v
          rv_p_oqpp => sv_p_oqkkl(1,ib)%v
           call dpzero ( sv_p_orhoat( 1 , ib )%v , nr * nlml * nsp ) 

           call dpzero ( sv_p_orhoat( 2 , ib )%v , nr * nlml * nsp ) 


C   --- Assemble rho1 = products of augmented functions ---
C   ... Set up all radial head and tail envelope functions, and their bc's
          allocate(fh_rv(nr*(lmxh+1)*nkaph))

          allocate(xh_rv(nr*(lmxh+1)*nkaph))

          allocate(vh_rv((lmxh+1)*nkaph))

          allocate(dh_rv((lmxh+1)*nkaph))

          call fradhd ( nkaph , eh , rsmh , lh , lmxh , nr , rofi_rv 
     .    , fh_rv , xh_rv , vh_rv , dh_rv )

          allocate(fp_rv(nr*(lmxa+1)*(kmax+1)))

          allocate(xp_rv(nr*(lmxa+1)*(kmax+1)))

          allocate(vp_rv((lmxa+1)*(kmax+1)))

          allocate(dp_rv((lmxa+1)*(kmax+1)))

          call fradpk ( kmax , rsma , lmxa , nr , rofi_rv , fp_rv 
     .    , xp_rv , vp_rv , dp_rv )


C   ... Augmented wave functions
          allocate(ul_rv(nr*(lmxa+1)*nsp))

          allocate(sl_rv(nr*(lmxa+1)*nsp))

          allocate(gz_rv(nr*(lmxa+1)*nsp))

          allocate(ruu_rv(nr*(lmxa+1)*2*nsp))

          allocate(rus_rv(nr*(lmxa+1)*2*nsp))

          allocate(rss_rv(nr*(lmxa+1)*2*nsp))

          call uspecb(0,4,sspec,is,is,lh,rsml,ehl,k)
           call makusp ( n0 , z , nsp , rmt , lmxa , rv_p_ov0 , a , nr , 
     .     xx , xx , pnu , pnz , rsml , ehl , ul_rv , sl_rv , gz_rv , ruu_rv 
     .     , rus_rv , rss_rv ) 



C   ... Contracted density matrix as coffs to products of (u,s,gz)
          k = max(nkaph,1+kmax)**2*(lmxa+1)**2*nlml*nsp
          allocate(chh_rv(k))
          chh_rv(:)=0d0
          allocate(chp_rv(k))
          chp_rv(:)=0d0
          allocate(cpp_rv(k))
          cpp_rv(:)=0d0
          allocate(dmatl_rv((lmxa+1)**2*nlml*nsp*9))
          dmatl_rv(:)=0d0
          call mkrou1 ( nsp , nlmh , nlma , nlml , kmax , rv_p_ocg , iv_p_ojcg 
     .     , iv_p_oidxcg , nkaph , nkapi , norb , ltab , ktab , blks , rv_p_oqhh 
     .     , rv_p_oqhp , rv_p_oqpp , vh_rv , dh_rv , vp_rv , dp_rv , chh_rv 
     .     , chp_rv , cpp_rv , dmatl_rv ) 



C        call prmx('dmatl',w(odmatl),(lmxa+1)**2,(lmxa+1)**2,nlml*nsp*9)

C   ... True local density for this sphere
           call mkrou2 ( nsp , lmxa , nlml , pnz , dmatl_rv , nr , ul_rv 
     .     , sl_rv , gz_rv , ruu_rv , rus_rv , rss_rv , sv_p_orhoat( 1 , ib )%v 
     .     ) 


C       call prrmsh('rout',w(orofi),w(orhoat(1,ib)),nr,nr,nlml*nsp)
C       call pshpr(40)
C       call wrhomt('rhoMT.','density',ib,w(orhoat(1,ib)),w(orofi),nr,
C    .    nlml,nsp)
C       call poppr

C   --- Assemble rho2 = unaugmented products Pkl*Pk'l' ---
C       H H product.  Include local orbitals w/ envelopes (nkapi->nkaph)
           call mkrou5 ( nsp , nr , nlml , nkaph , nkaph , fh_rv , lmxh 
     .     , nkaph , nkaph , fh_rv , lmxh , chh_rv , sv_p_orhoat( 2 , ib )%v 
     .     ) 


C       H Pkl product
           call mkrou5 ( nsp , nr , nlml , nkaph , nkaph , fh_rv , lmxh 
     .     , kmax + 1 , kmax + 1 , fp_rv , lmxa , chp_rv , sv_p_orhoat( 2 , ib )%v 
     .     ) 


C       Pkl Pkl product
           call mkrou5 ( nsp , nr , nlml , kmax + 1 , kmax + 1 , fp_rv , 
     .     lmxa , kmax + 1 , kmax + 1 , fp_rv , lmxa , cpp_rv , sv_p_orhoat( 2 , ib )%v 
     .     ) 


C       call prrmsh('sm rout',w(orofi),w(orhoat(2,ib)),nr,nr,nlml*nsp)

C   ... Site charge and eigenvalue sum decomposed by l
          call mkrou3 ( 2 , lmxa , nlml , nsp , pnz , dmatl_rv , hab 
     .    ( 1 , 1 , 1 , ib ) , sab ( 1 , 1 , 1 , ib ) , qbyl ( 1 , 1 , 
     .    ib ) , hbyl ( 1 , 1 , ib ) )

C       Remake hbyl from energy-weighted density matrix
          if (lekkl .ne. 0) then
            rv_p_oqhh => sv_p_oeqkkl(3,ib)%v
            rv_p_oqhp => sv_p_oeqkkl(2,ib)%v
            rv_p_oqpp => sv_p_oeqkkl(1,ib)%v
            nlml1 = 1
            call dpzero ( dmatl_rv , ( lmxa + 1 ) * * 2 * nlml1 * nsp 
     .      * 9 )

            k = max(nkaph,1+kmax)**2*(lmxa+1)**2*nlml1*nsp
            call dpzero ( chh_rv , k )

            call dpzero ( chp_rv , k )

            call dpzero ( cpp_rv , k )

                  call mkrou1 ( nsp , nlmh , nlma , nlml1 , kmax , rv_p_ocg , iv_p_ojcg 
     .            , iv_p_oidxcg , nkaph , nkapi , norb , ltab , ktab , blks , rv_p_oqhh 
     .            , rv_p_oqhp , rv_p_oqpp , vh_rv , dh_rv , vp_rv , dp_rv , chh_rv 
     .            , chp_rv , cpp_rv , dmatl_rv ) 






#ifdef USEOPTIONAL
             call mkrou3 ( mode=1 , lmxa=lmxa , nlml=nlml1 , nsp=nsp 
     .       ,pnz= pnz , dmatl=dmatl_rv 
     .       , sab=sab ( 1 , 1 , 1 , ib ) , qsum=hbyl ( 1 , 1 , ib )  ) 
#else
             call mkrou3 ( 1 , lmxa , nlml1 , nsp , pnz , dmatl_rv , iwdummy 
     .       , sab ( 1 , 1 , 1 , ib ) , hbyl ( 1 , 1 , ib ) , iwdummy ) 
#endif


          endif

C   --- Print charges for information ---
           call radsum ( nr , nr , nlml , nsp , rwgt_rv , sv_p_orhoat( 1 , ib )%v 
     .     , sum1 ) 
           call radsum ( nr , nr , nlml , nsp , rwgt_rv , sv_p_orhoat( 2 , ib )%v 
     .     , sum2 ) 

ccccccccccccccccccccc
c          print *,'ddd sum1 sum2=',sum1,sum(abs(sv_p_orhoat( 1 , ib )%v))
c     &    ,size(sv_p_orhoat( 1 , ib )%v),nr*nlml*nsp
ccccccccccccccccccccc
          sum1 = sum1/y0
          sum2 = sum2/y0
          if (nsp .eq. 2) then
             sums1 = sum1 - 2 * ddot ( nr , rwgt_rv , 1 , sv_p_orhoat( 1 , ib )%v 
     .       , 1 ) / y0 


             sums2 = sum2 - 2 * ddot ( nr , rwgt_rv , 1 , sv_p_orhoat( 2 , ib )%v 
     .       , 1 ) / y0 


          endif
          if (ipr.ge.30 .and. nsp .eq. 1)
     .    write(stdo,200) ib,sum1,sum2,sum1-sum2
          if (ipr.ge.30 .and. nsp .eq. 2)
     .    write(stdo,200) ib,sum1,sum2,sum1-sum2,
     .    -sums1,-sums2,-sums1+sums2
  200     format(i4,3f12.6,2x,3f12.6)
        endif

C       Contribution to mag outsite rmt: extrapolate tail to infinity
         call mkrou6 ( rofi_rv , sv_p_orhoat( 1 , ib )%v , nr , nlml , 
     .   nsp , xx , xx , res ) 


        if (ipr.ge.30 .and. res .ne. 0) then
          write(stdo,211) res,res-sums1
  211     format(7x,'contr. to mm extrapolated for r>rmt:',f11.6,
     .    ' est. true mm =',f9.6)
        endif

C   --- Make new core density and core eigenvalue sum ---
        if (lfoc .eq. 0) then
          call pshpr(ipr+11)
            call getcor ( 0 , z , a , pnu , pnz , nr , lmxa , rofi_rv , rv_p_ov01 
     .      , kcor , lcor , qcor , smec , smtc , sv_p_orhoat( 3 , ib )%v 
     .      , ncore , 0d0 , 0d0 ) 



          call poppr
C     ... here if v0 used to make core... then effectively frozen
C         call p1kden(nr,w(orwgt),w(orhoat(3,ib)),w(ov01),w(ov0),sum)
C         write(stdo,996) smec,sum,smec+sum
C 996     format('smec,sum,smec+sum=',3f14.6)
C         smec = smec+sum

          sumtc = sumtc + smtc
          sumec = sumec + smec
        else
          if (ipr .ge. 41) write(stdo,288) stc0
  288     format(' foca..  use smtc = ',f14.6)
          sumtc = sumtc + stc0
          sumt0 = sumt0 + stc0
            call dpcopy ( rv_p_orhoca , sv_p_orhoat( 3 , ib )%v , 1 , nr 
     .      * nsp , 1d0 ) 


        endif

        if (allocated(dmatl_rv)) deallocate(dmatl_rv)
        if (allocated(cpp_rv)) deallocate(cpp_rv)
        if (allocated(chp_rv)) deallocate(chp_rv)
        if (allocated(chh_rv)) deallocate(chh_rv)
        if (allocated(rss_rv)) deallocate(rss_rv)
        if (allocated(rus_rv)) deallocate(rus_rv)
        if (allocated(ruu_rv)) deallocate(ruu_rv)
        if (allocated(gz_rv)) deallocate(gz_rv)
        if (allocated(sl_rv)) deallocate(sl_rv)
        if (allocated(ul_rv)) deallocate(ul_rv)
        if (allocated(dp_rv)) deallocate(dp_rv)
        if (allocated(vp_rv)) deallocate(vp_rv)
        if (allocated(xp_rv)) deallocate(xp_rv)
        if (allocated(fp_rv)) deallocate(fp_rv)
        if (allocated(dh_rv)) deallocate(dh_rv)
        if (allocated(vh_rv)) deallocate(vh_rv)
        if (allocated(xh_rv)) deallocate(xh_rv)
        if (allocated(fh_rv)) deallocate(fh_rv)
        if (allocated(rwgt_rv)) deallocate(rwgt_rv)
        if (allocated(rofi_rv)) deallocate(rofi_rv)


   10   continue
      enddo

C --- Put sumec,sumtc into etot struct ---

      i_copy_size=size(sham%eterms)
      call dcopy(i_copy_size,eterms,1,sham%eterms,1)


      call tcx('mkrout')
      end subroutine mkrout


      subroutine mkrou1(nsp,nlmh,nlma,nlml,kmax,cg,jcg,indxcg,
     .nkaph,nkapi,norb,ltab,ktab,blks,qhh,qhp,qpp,vh,dh,vp,dp,
     .chh,chp,cpp,dmatl)

C- Contracted density matrix for one site
C ----------------------------------------------------------------------
Ci Inputs
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nlmh  :(lmxh+1)**2, where lmxh = L-cutoff in basis
Ci   nlma  :(lmxa+1)**2, where lmxa = L-cutoff in augmentation
Ci   nlml  :L-cutoff for charge density on radial mesh
Ci   kmax  :polynomial cutoff in augmentation
Ci   cg    :Clebsch Gordon coefficients, stored in condensed form (scg.f)
Ci   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
Ci   indxcg:index for Clebsch Gordon coefficients
Ci   nkaph :number of types of one l-quantum number in the MTO basis
Ci   nkapi :number of the nkaph functions that are envelope functions
Ci   norb  :number of orbital types for this site; see orbl.f
Ci   ltab  :table of l-quantum numbers for each type
Ci   ktab  :table of energy index for each type
Ci   blks  :blks(iorb) = size of contiguous block of orbitals for
Ci         :orbital type iorb and possibly iorb+1,iorb+2...
Ci   qhh   :head-head component of site part of density matrix (addrbl.f)
Ci   qhp   :head-tail component of site part of density matrix (addrbl.f)
Ci   qpp   :tail-tail component of site part of density matrix (addrbl.f)
Ci   vh    :values of head functions on MT sphere
Ci   dh    :slopes of tail functions on MT sphere
Ci   vp    :values of Pkl  functions on MT sphere
Ci   dp    :slopes of Pkl  functions on MT sphere
Co Outputs
Co   chh   :head-head product function coefficients
Co   chp   :head-tail product function coefficients
Co   cpp   :tail-tail product function coefficients
Co   dmatl :dmatl(l1,l2,mlm,i,j,isp) holds coefficients to a Y_lm
Co         :expansion of the function products f_i(l1) f_j(l2)
Co         :where f_i,f_j form this table.
Co         :    (uu  us  uz)
Co         :    (su  ss  sz)
Co         :    (zu  zs  zz)
Cr Remarks
Cr   Transforms density-matrix as generated by addrbl.f into
Cr   contracted density-matrix of wave function products.
Cr   Radial parts of chh,chp,cpp are
Cr     (psi1)_l (psi2)_l'
Cr   Radial parts of dmatl are:
Cr     (u, s, or gz)_l (u, s, or gz)_l'
Cu Updates
Cu   29 Jul 08 Adapted from mkrout.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nsp,nlmh,nlma,nlml,kmax,norb,nkaph,nkapi
      integer jcg(1),indxcg(1)
      integer ltab(norb),ktab(norb),blks(norb)
      double precision qhh(nkaph,nkaph,nlmh,nlmh,nsp)
      double precision qhp(nkaph,kmax+1,nlmh,nlma,nsp)
      double precision qpp(kmax+1,kmax+1,nlmh,nlma,nsp)
C     double precision dmatl(0:lmxa,0:lmxa,nlml,3,3,nsp)
C     double precision vh(0:lmxh,nkaph),dh(0:lmxh,nkaph),
C    .                 vp(0:lmxa,kmxa+1),dp(0:lmxa,kmxa+1)
C     chh = chh(nkaph,nkaph,0:lmxh,0:lmxh,nlml,nsp)
C     chp = chp(nkaph,kmax+1,0:lmxh,0:lmxa,nlml,nsp)
C     cpp = cpp(kmax+1,kmax+1,0:lmxa,0:lmxa,nlml,nsp)
      double precision dmatl(*),vh(*),dh(*),vp(*),dp(*),
     .chh(*),chp(*),cpp(*)

      double precision cg(1)
      integer ll,k,lmxa,lmxh
      integer n0,nkap0
      parameter (n0=10,nkap0=3)
      integer ltba(n0*nkap0),ktba(n0*nkap0),blka(n0*nkap0)

ccccccccccccccccccccc
c      print *,'sumcheck qhh,qhp,qpp=',sum(qhh),sum(qhp),sum(qpp)
c      stop 'test end'
ccccccccccccccccccccc

C     lmxl = ll(nlml)
      lmxa = ll(nlma)
      lmxh = ll(nlmh)
      do  k = 0, kmax
        ltba(k+1) = 0
        ktba(k+1) = k+1
        blka(k+1) = nlma
      enddo

C ... Contracted density-matrix chh,chp,cpp from qhh,qhp,qpp
C     H H product
      call mkrou4(nsp,nlml,cg,jcg,indxcg,
     .nkaph,norb,ltab,ktab,blks,lmxh,nlmh,
     .nkaph,norb,ltab,ktab,blks,lmxh,nlmh,
     .qhh,chh)
C     H Pkl product
      call mkrou4(nsp,nlml,cg,jcg,indxcg,
     .nkaph,norb,ltab,ktab,blks,lmxh,nlmh,
     .kmax+1,kmax+1,ltba,ktba,blka,lmxa,nlma,
     .qhp,chp)
C     Pkl Pkl product
      call mkrou4(nsp,nlml,cg,jcg,indxcg,
     .kmax+1,kmax+1,ltba,ktba,blka,lmxa,nlma,
     .kmax+1,kmax+1,ltba,ktba,blka,lmxa,nlma,
     .qpp,cpp)

C ... Contracted density matrix as coffs to products of (u,s,gz)
C     H H product
      call mkcfus(nsp,lmxa,nlml,
     .nkaph,nkapi,vh,dh,lmxh,
     .nkaph,nkapi,vh,dh,lmxh,
     .chh,dmatl)
C     H Pkl product
      call mkcfus(nsp,lmxa,nlml,
     .nkaph,nkapi,vh,dh,lmxh,
     .kmax+1,kmax+1,vp,dp,lmxa,
     .chp,dmatl)
C     Pkl Pkl product
      call mkcfus(nsp,lmxa,nlml,
     .kmax+1,kmax+1,vp,dp,lmxa,
     .kmax+1,kmax+1,vp,dp,lmxa,
     .cpp,dmatl)

      end subroutine mkrou1


      subroutine mkrou2(nsp,lmxa,nlml,pnz,dmatl,nr,ul,sl,gz,ruu,rus,rss,
     .rho)

C- Assemble true site density from product function coefficients
C ----------------------------------------------------------------------
Ci Inputs
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   lmxa  :augmentation l-cutoff
Ci   nlml  :l-cutoff for charge density
Ci   dmatl :dmatl(l1,l2,mlm,i,j,isp) holds coefficients to a y_lm
Ci         :expansion of the function products f_i(l1) f_j(l2)
Ci         :where f_i,f_j form this table.
Ci         :    (uu  us  uz)
Ci         :    (su  ss  sz)
Ci         :    (zu  zs  zz)
Ci   nr    :number of radial mesh points
Ci   ul    :r*radial wave functions; see remarks
Ci   sl    :r*radial wave functions; see remarks
Ci   gz    :r*semicore wave functions; see remarks
Ci   ruu   :l-diagonal (uu) product including small component (makusp)
Ci   rus   :l-diagonal (us) product including small component (makusp)
Ci   rss   :l-diagonal (ss) product including small component (makusp)
Co Outputs
Co   rho   :charge density assembled
Cr Remarks
Cr   This routine uses linear combinations (u,s) of phi,phidot
Cr   defined as : u has val=1, slo=1 at rmax, s has val=0, slo=1
Cr   ul and sl are returned as r * u and r * s, respectively.
Cr
Cr   Let phi_z be the w.f. corresponding to pnu_z.
Cr   A local orbital of the first type is defined as follows.
Cr      gz = r * ( phi_z - phi_z(rmax) u - (phi_z)'(rmax) s )
Cr   By construction, gz/r has both value = 0 and slope = 0 at rmax.
Cr   A local orbital of the second type is defined as gz=r*phi_z;
Cr   for r>rmax a smooth Hankel tail (spec'd by ehl,rsml) is attached.
Cr
Cr   Spherical density assumbled using ruu,rus,rss which include
Cr   small comonent of relativistic wave function.
Cu Updates
Cu   28 aug 01 extended to local orbitals.  altered argument list.
C ----------------------------------------------------------------------
C     implicit none
C ... passed parameters
      integer nsp,lmxa,nlml,nr,n0
      parameter (n0=10)
      double precision rho(nr,nlml,nsp),pnz(n0,nsp),
     .dmatl(0:lmxa,0:lmxa,nlml,3,3,nsp),
     .ul(nr,0:lmxa,nsp),ruu(nr,0:lmxa,2,nsp),
     .sl(nr,0:lmxa,nsp),rus(nr,0:lmxa,2,nsp),
     .gz(nr,0:lmxa,nsp),rss(nr,0:lmxa,2,nsp)
C ... local parameters
      logical lpz1,lpz2
      integer l,i,mlm,l1,l2,isp
      double precision xuu,xus,xsu,xss,xuz,xsz,xzu,xzs,xzz
C     double precision top

      call tcn('mkrou2')
ccccccccccc
c      print *,'ddd ', sum(abs(dmatl(:,:,:,1:2,1:2,:)))
c      print *,'ddd ', sum(abs(ul))
c      print *,'ddd ruu ', sum(abs(ruu))
c      print *,'ddd ', sum(abs(rus))
c      print *,'ddd ', sum(abs(rss))
c      print *,'ddd ul ', sum(abs(ul))
c      print *,'ddd ', sum(abs(sl))

C --- Full density as products (u,s) (u,s); no small component ---
      do  isp = 1, nsp
        do  mlm = 1, nlml
          do  l1 = 0, lmxa
            do  l2 = 0, lmxa
              lpz1 = pnz(l1+1,1) .ne. 0
              lpz2 = pnz(l2+1,1) .ne. 0

              xuu = dmatl(l1,l2,mlm,1,1,isp)
              xus = dmatl(l1,l2,mlm,1,2,isp)
              xsu = dmatl(l1,l2,mlm,2,1,isp)
              xss = dmatl(l1,l2,mlm,2,2,isp)

C           top = xuu*xuu + xus*xus + xsu*xsu + xss*xss
C           if (dsqrt(top).gt.1d-6)
C    .        write (6,700) mlm,l1,l2,xuu,xus,xsu,xss
C 700       format(3i5,4f14.8)

              do  i = 1, nr
                rho(i,mlm,isp) = rho(i,mlm,isp)
     .          + xuu * ul(i,l1,isp) * ul(i,l2,isp)
     .          + xsu * sl(i,l1,isp) * ul(i,l2,isp)
     .          + xus * ul(i,l1,isp) * sl(i,l2,isp)
     .          + xss * sl(i,l1,isp) * sl(i,l2,isp)
              enddo

              if (lpz1 .or. lpz2) then

                xuz = dmatl(l1,l2,mlm,1,3,isp)
                xsz = dmatl(l1,l2,mlm,2,3,isp)
                xzu = dmatl(l1,l2,mlm,3,1,isp)
                xzs = dmatl(l1,l2,mlm,3,2,isp)
                xzz = dmatl(l1,l2,mlm,3,3,isp)
                if (xuz.ne.0 .or. xsz.ne.0 .or. xzu.ne.0 .or.
     .          xzs.ne.0 .or. xzz.ne.0) then
                  do  i = 1, nr
                    rho(i,mlm,isp) = rho(i,mlm,isp)
     .              + xuz * ul(i,l1,isp) * gz(i,l2,isp)
     .              + xsz * sl(i,l1,isp) * gz(i,l2,isp)
     .              + xzu * gz(i,l1,isp) * ul(i,l2,isp)
     .              + xzs * gz(i,l1,isp) * sl(i,l2,isp)
     .              + xzz * gz(i,l1,isp) * gz(i,l2,isp)
                  enddo
                endif
              endif

            enddo
          enddo
        enddo
      enddo

C --- Remake spherical density including small component ---
c      print *, 'skip remaking spher. density'
c      return
      call dpzero(rho,nr)
      call dpzero(rho(1,1,nsp),nr)
      do  isp = 1, nsp
        do  l = 0, lmxa
          xuu = dmatl(l,l,1,1,1,isp)
          xus = dmatl(l,l,1,1,2,isp) + dmatl(l,l,1,2,1,isp)
          xss = dmatl(l,l,1,2,2,isp)
          xuz = dmatl(l,l,1,1,3,isp) + dmatl(l,l,1,3,1,isp)
          xsz = dmatl(l,l,1,2,3,isp) + dmatl(l,l,1,3,2,isp)
          xzz = dmatl(l,l,1,3,3,isp)

          do  i = 1, nr
            rho(i,1,isp) = rho(i,1,isp)
     .      +xuu*ruu(i,l,1,isp)+xus*rus(i,l,1,isp)+xss*rss(i,l,1,isp)
          enddo

          if (pnz(l+1,1) .ne. 0) then
            do  i = 1, nr
              rho(i,1,isp) = rho(i,1,isp)
     .        +xuz*ruu(i,l,2,isp)+xsz*rus(i,l,2,isp)+xzz*rss(i,l,2,isp)
            enddo
          endif

        enddo
      enddo

      call tcx('mkrou2')
      end subroutine mkrou2



      subroutine mkrou4(nsp,nlml,cg,jcg,indxcg,
     .nk1,norb1,ltab1,ktab1,blks1,lmx1,ndim1,
     .nk2,norb2,ltab2,ktab2,blks2,lmx2,ndim2,
     .qkk12,ckk)

C- Assemble contracted density-matrix (in Y_lm form)
C ----------------------------------------------------------------------
Ci Inputs
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nlml  :charge density L-cutoff
Ci   cg    :Clebsch Gordon coefficients, stored in condensed form (scg.f)
Ci   jcg   :L q.n. for the C.G. coefficients stored in condensed form (scg.f)
Ci   indxcg:index for Clebsch Gordon coefficients
Ci   nk1   :no. orbital types for a given L of first function
Ci         :nk1 merely dimensions qkk12,ckk
Ci   norb1 :total number of orbitals of the first type
Ci   ltab1 :table of l quantum numbers for the norb1 orbitals
Ci   ktab1 :table of k numbers (orbital type) for the norb1 orbitals
Ci   blks1 :block size for grouping orbitals into blocks (gtbls1)
Ci   lmx1  :dimensions ckk
Ci   ndim1 :dimensions qkk12; ndim1 should be (lmx1+1)**2
Ci   nk2   :no. orbital types for a given L of second function
Ci         :nk2 merely dimensions qkk12,ckk
Ci   norb2 :total number of orbitals of the first type
Ci   ltab2 :table of l quantum numbers for the norb1 orbitals
Ci   ktab2 :table of k numbers (orbital type) for the norb1 orbitals
Ci   blks2 :block size for grouping orbitals into blocks (gtbls1)
Ci   lmx2  :dimensions ckk
Ci   ndim2 :dimensions qkk12; ndim1 should be (lmx1+1)**2
Ci   qkk12 :density matrix between norb1 and norb2 orbitals
Co Outputs
Co   ckk   :product function coefficients
Cr Remarks
Cr   Wave function products have a radial part, written as products
Cr   of the radial functions (psi1)_l (psi2)_l'
Cr   The angular part Y_L Y_L' is contracted into a single index M
Cr   using Clebsch Gordan coefficients.  Thus ckk has indices
Cr      ckk(k1,k2,l1,l2,M),
Cr   where k1 ranges over the orbital types of the first function
Cr   and   k2 ranges over the orbital types of the second function.
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nsp,nlml,norb1,lmx1,ndim1,norb2,lmx2,ndim2,nk1,nk2
      integer jcg(1),indxcg(1)
      integer ltab1(norb1),ktab1(norb1),blks1(norb1),
     .ltab2(norb2),ktab2(norb2),blks2(norb2)
      double precision qkk12(nk1,nk2,ndim1,ndim2,nsp),cg(1),
     .ckk(nk1,nk2,0:lmx1,0:lmx2,nlml,nsp)
C ... Local parameters
      integer ilm1,io1,l1,nlm11,nlm12,k1,
     .ilm2,io2,l2,nlm21,nlm22,k2,
     .icg,ll,mlm,ix,isp
      double precision xx

C     call tcn('mkrou4')

      do  isp = 1, nsp
        do  io2 = 1, norb2
          if (blks2(io2) .ne. 0) then
C       k2,l2 = k and starting l index for this block
            l2 = ltab2(io2)
            k2 = ktab2(io2)
            nlm21 = l2**2+1
            nlm22 = nlm21 + blks2(io2)-1
            do  ilm2 = nlm21, nlm22
              l2 = ll(ilm2)
              do  io1 = 1, norb1
                if (blks1(io1) .ne. 0) then
C           k1,l1 = k and starting l index for this block
                  l1 = ltab1(io1)
                  k1 = ktab1(io1)
                  nlm11 = l1**2+1
                  nlm12 = nlm11 + blks1(io1)-1
                  do  ilm1 = nlm11, nlm12
                    l1 = ll(ilm1)

                    ix = max0(ilm1,ilm2)
                    ix = (ix*(ix-1))/2 + min0(ilm1,ilm2)
                    do icg = indxcg(ix),indxcg(ix+1)-1
                      mlm = jcg(icg)
                      if (mlm .le. nlml) then
                        xx = cg(icg)*qkk12(k1,k2,ilm1,ilm2,isp)
                        ckk(k1,k2,l1,l2,mlm,isp) = ckk(k1,k2,l1,l2,mlm,isp)+xx
                      endif
                    enddo
                  enddo
                endif
              enddo
            enddo
          endif
        enddo
      enddo

C     call tcx('mkrou4')
      end subroutine mkrou4


      subroutine mkcfus(nsp,lmxa,nlml,nf1,nf1s,val1,slo1,lmx1,
     .nf2,nf2s,val2,slo2,lmx2,ckk,dmatl)

C- Assemble contracted density matrix as coffs to products of (u,s,gz)
C ----------------------------------------------------------------------
Ci Inputs
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   lmxa  :leading dimension of dmatl
Ci   nlml  :charge density L-cutoff
Ci   nf1   :number of orbital types of first kind for each l
Ci   nf1s  :number of functions of first kind for each l, for which
Ci         :there is a smooth part to be subtracted (which also
Ci         :corresponds to the functions which connect to envelope
Ci         :functions)
Ci   val1  :function values at MT boundary, first function
Ci   slo1  :function slopes at MT boundary, first function
Ci   lmx1  :dimensions val1,slo1
Ci   nf2   :number of orbital types of second kind for each l
Ci   nf2s  :number of functions of second kind for each l, for which
Ci         :a smooth part is to be subtracted (which also
Ci         :corresponds to the functions which connect to envelope
Ci         :functions)
Ci   val2  :function values at MT boundary, second function
Ci   slo2  :function slopes at MT boundary, second function
Ci   lmx2  :dimensions val1,slo1
Ci   ckk   :density matrix between  and  orbitals
Co Outputs
Co   dmatl :dmatl(l1,l2,mlm,i,j,isp) holds coefficients to a Y_lm
Co         :expansion of the function products f_i(l1) f_j(l2)
Co         :where f_i,f_j form this table.
Co         :    (uu  us  uz)
Co         :    (su  ss  sz)
Co         :    (zu  zs  zz)
Cr Remarks
Cr   Transforms contracted density-function matrix into coefficients
Cr   of wave function products. Radial parts are
Cr     (u, s, or gz)_l (u, s, or gz)_l'
Cu Updates
Cu   28 Aug 01 Extended to local orbitals.  Altered argument list.
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nsp,nlml,lmxa,lmx1,lmx2,nf1,nf2,nf1s,nf2s
      double precision ckk(nf1,nf2,0:lmx1,0:lmx2,nlml,nsp),
     .val1(0:lmx1,nf1s),slo1(0:lmx1,nf1s),
     .val2(0:lmx2,nf2s),slo2(0:lmx2,nf2s)
      double precision dmatl(0:lmxa,0:lmxa,nlml,3,3,nsp)
C ... Local parameters
      integer l1,k1,l2,k2,mlm,isp
      double precision xx

C     call tcn('mkcfus')

      do  isp = 1, nsp
        do  mlm = 1, nlml
          do  k2 = 1, nf2s
            do  k1 = 1, nf1s
              do  l2 = 0, lmx2
                do  l1 = 0, lmx1
                  xx = ckk(k1,k2,l1,l2,mlm,isp)

                  dmatl(l1,l2,mlm,1,1,isp) = dmatl(l1,l2,mlm,1,1,isp)
     .            + xx * val1(l1,k1) * val2(l2,k2)
                  dmatl(l1,l2,mlm,1,2,isp) = dmatl(l1,l2,mlm,1,2,isp)
     .            + xx * val1(l1,k1) * slo2(l2,k2)
                  dmatl(l1,l2,mlm,2,1,isp) = dmatl(l1,l2,mlm,2,1,isp)
     .            + xx * slo1(l1,k1) * val2(l2,k2)
                  dmatl(l1,l2,mlm,2,2,isp) = dmatl(l1,l2,mlm,2,2,isp)
     .            + xx * slo1(l1,k1) * slo2(l2,k2)

                enddo
              enddo
            enddo
          enddo
        enddo
      enddo

C     call tcx('mkcfus')

C --- Products involving local orbitals ---
      if (nf1s .ge. nf1 .and. nf2s .ge. nf2) return
C     call tcn('mkcfus')

      do  isp = 1, nsp
        do  mlm = 1, nlml
          do  k2 = 1, nf2
            do  k1 = 1, nf1
              if (k1 .gt. nf1s .or. k2 .gt. nf2s) then
                do  l2 = 0, lmx2
                  do  l1 = 0, lmx1
                    xx = ckk(k1,k2,l1,l2,mlm,isp)

C               sc-sc product
                    if (k1 .gt. nf1s .and. k2 .gt. nf2s) then
                      dmatl(l1,l2,mlm,3,3,isp) = dmatl(l1,l2,mlm,3,3,isp)+xx

C               sc-valence product
                    elseif (k1 .gt. nf1s) then
                      dmatl(l1,l2,mlm,3,1,isp) = dmatl(l1,l2,mlm,3,1,isp)
     .                + xx * val2(l2,k2)
                      dmatl(l1,l2,mlm,3,2,isp) = dmatl(l1,l2,mlm,3,2,isp)
     .                + xx * slo2(l2,k2)

C               valence-sc product
                    elseif (k2 .gt. nf2s) then
                      dmatl(l1,l2,mlm,1,3,isp) = dmatl(l1,l2,mlm,1,3,isp)
     .                + val1(l1,k1) * xx
                      dmatl(l1,l2,mlm,2,3,isp) = dmatl(l1,l2,mlm,2,3,isp)
     .                + slo1(l1,k1) * xx

                    endif

C                if (l1.eq.2.and.l2.eq.2 .and. mlm.eq.1) then
C                  print 333, k1,k2,xx,
C     .              dmatl(l1,l2,mlm,1,3,isp),dmatl(l1,l2,mlm,3,1,isp)
C  333             format(2i4,3f12.6)
C                endif

                  enddo
                enddo
              endif
            enddo
          enddo
        enddo
      enddo

C     call tcx('mkcfus')

C      isp = 1
C      do  mlm = 1, nlml
C        print *, ' '
C        do  l2 = 0, lmx2
C          do  l1 = 0, lmx1
C            print 334, l1,l2,mlm,
C     .        dmatl(l1,l2,mlm,1,2,isp),
C     .        dmatl(l2,l1,mlm,2,1,isp),
C     .        dmatl(l1,l2,mlm,1,2,isp)-dmatl(l2,l1,mlm,2,1,isp)
C            print 334, l1,l2,mlm,
C     .        dmatl(l1,l2,mlm,1,3,isp),
C     .        dmatl(l2,l1,mlm,3,1,isp),
C     .        dmatl(l1,l2,mlm,1,3,isp)-dmatl(l2,l1,mlm,3,1,isp)
C            print 334, l1,l2,mlm,
C     .        dmatl(l1,l2,mlm,2,3,isp),
C     .        dmatl(l2,l1,mlm,3,2,isp),
C     .        dmatl(l1,l2,mlm,2,3,isp)-dmatl(l2,l1,mlm,3,2,isp)
C
C  334       format(3i4,2f12.6,f14.8)
C          enddo
C        enddo
C      enddo
C      pause

      end subroutine mkcfus


      subroutine mkrou5(nsp,nr,nlml,nf1,nf1s,f1,lmx1,nf2,nf2s,f2,lmx2,
     .ckk,rho)

C- Assemble smooth site density from contracted density matrix
C ----------------------------------------------------------------------
Ci Inputs
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nlml  :charge density L-cutoff
Ci   nf1   :number of orbital types of first kind for each l
Ci         :nf1 merely dimensions ckk
Ci   nf1s  :number of functions of first kind for each l, for which
Ci         :there are functions of f1 type defined
Ci   f1    :first function on a radial mesh
Ci   lmx1  :l-cutoff for f1 functions
Ci   nf2   :number of orbital types of second kind for each l
Ci         :nf2 merely dimensions ckk
Ci   nf2s  :number of functions of second kind for each l, for which
Ci         :there are functions of f2 type defined
Ci   f2    :second function on a radial mesh
Ci   lmx2  :l-cutoff for f2 functions
Ci   ckk   :density matrix between f1 and f2 orbitals
Co Outputs
Co   rho   :density assembled on mesh
Cr Remarks
Cu Updates
Cu   28 Aug 01 Extended to local orbitals.  Altered argument list.
Cu             Envelopes f1,f2 must be zero for all channels that
Cu             have no smooth counterparts to subtract.
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nsp,nr,nlml,lmx1,lmx2,nf1,nf2,nf1s,nf2s
      double precision ckk(nf1,nf2,0:lmx1,0:lmx2,nlml,nsp),
     .f1(nr,0:lmx1,nf1s),f2(nr,0:lmx2,nf2s),rho(nr,nlml,nsp)
C ... Local parameters
      integer l1,k1,l2,k2,mlm,isp,i
      double precision xx

      call tcn('mkrou5')

      do  isp = 1, nsp
        do  mlm = 1, nlml
          do  k2 = 1, nf2s
            do  k1 = 1, nf1s
              do  l2 = 0, lmx2
                do  l1 = 0, lmx1
                  xx = ckk(k1,k2,l1,l2,mlm,isp)
                  if (dabs(xx) .gt. 1d-16) then
                    do  i= 1, nr
                      rho(i,mlm,isp) = rho(i,mlm,isp)
     .                + xx*f1(i,l1,k1)*f2(i,l2,k2)
                    enddo
                  endif
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo

      call tcx('mkrou5')
      end subroutine mkrou5


      subroutine mkrou6(rofi,rho,nr,nlml,nsp,rho0,decay,res)

C- Fit tail of spin density; integrate charge beyond MT sphere
C ----------------------------------------------------------------------
Ci Inputs
Ci   rofi  :radial mesh points
Ci   rho   :spin-polarized charge density
Ci   nr    :number of radial mesh points
Ci   nlml  :L-cutoff for charge density on radial mesh
Co Outputs
Co   rho0  :fit density of form rho0*exp(-decay*r)
Co   decay :fit density of form rho0*exp(-decay*r)
Co   res   :integral of fit density from rofi(nr) to infinity
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   27 Sep 03  First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nr,nlml,nsp
      double precision rofi(nr),rho(nr,nlml,2),rho0,decay,res
C ... Local parameters
      integer ir
      double precision norm(2,2),tnorm(2,2),rhs(2),y,dy,fac,y0,r0,pi,a,b

      res = 0
      if (nr .lt. 10 .or. nsp .eq. 1) return
      call dpzero(norm,4)
      call dpzero(rhs,2)
      fac = 1
      if (rho(nr,1,1) .lt. rho(nr,1,2)) fac = -1
      do  ir = nr-5, nr
        y = fac*(rho(ir,1,1) - rho(ir,1,2))/rofi(ir)**2
C       If the spin density changes sign, nonsensical to try and fit
        if (y .le. 0) return
        dy = dlog(y)

        norm(1,1) = norm(1,1) + 1
        norm(1,2) = norm(1,2) + rofi(ir)
        norm(2,1) = norm(2,1) + rofi(ir)
        norm(2,2) = norm(2,2) + rofi(ir)**2

        rhs(1) = rhs(1) + dy
        rhs(2) = rhs(2) + rofi(ir)*dy
      enddo

      call dinv22(norm,tnorm)
      a = tnorm(1,1)*rhs(1) + tnorm(1,2)*rhs(2)
      b = tnorm(2,1)*rhs(1) + tnorm(2,2)*rhs(2)
      pi = 4d0*datan(1d0)
      y0 = 1d0/dsqrt(4d0*pi)
      a = fac*exp(a)/y0
      b = -b

C     Nonsensical if density not decaying fast enough
      if (b .lt. 1) return
      r0 = rofi(nr)

C      a = 1d6
C      b = 2.5d0
C      r0 = 3

C     Integral a*exp(-b*r)*r*r = (2+2*b*r0+b**2*r0*2)/b**3*exp(-b*r0)
      res = a*(2+2*b*r0+b**2*r0**2)/b**3*exp(-b*r0)
      decay = b
      rho0 = a

      end subroutine mkrou6


