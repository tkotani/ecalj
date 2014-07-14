      program hx0fp0
      use m_readqg
      use m_readeigen
      use m_read_bzdata,ngrp2=>ngrp
      use m_genallcf_v3
      use keyvalue
c-------------------------------------------------
ctm, Takashi Miyake, 0605;  maxloc mode (ixc=111)
c  Calculate x0, W-V, \epsilon
Cr
Cr Be careful for the indexing...
Cr      A routine idxlnmc(nindxv,nindxc,...  in index.f
Cr      specifies the order of the  (Core wave)+(Argumentation wave) in each MT.
Cr      The total number of the wave are mnl(ic)= mnlc(ic) + mnlv(ic). 
Cr      The indexing starts with core first and then valence on top of core 
Cr      So n-index in in for valence electron is different from inv. 
C
Cr the key routine is x0k_v2.
c
c-------------------------------------------------
      implicit none
      real(8):: q(3),  qgbin(3),qx(3)
      real(8):: ua=1d0 ! this is a dummy.

c worksize in megabytes (1 word = 4 bytes)
c      integer(4),parameter :: mbytes=50, mwords=mbytes/4,
c     &  iwksize=mwords * 1000 *1000
c      integer w
c      common /w/ w(iwksize)
c      integer(4),parameter:: ngnmax=10,mxclass=100
c      real(8):: gen(9,ngnmax)
c      character*6 clabl(mxclass)
c      character*120 symgrp

      integer(4):: ifrb(2),ifcb(2),ifrhb(2),ifchb(2) !,ifev(2)
      integer(4):: ndble=8
      integer(4):: nword

c> for genallc.f
c      real(8):: plat(3,3) !,qbas(3,3),ginv(3,3)
c     l          ifphiv(2),ifphic(2),
c     l          iphi(2),iphidot(2),ipgs(2),ipgsdot(2),
c     l          iphiv(2),iphivd(2),ipgv(2),ipgvd(2),
c     l          iphic(2),ipgc(2)
c
      real(8),allocatable:: vxcfp(:,:),
     &  wqt(:), wgt0(:,:),q0i(:,:)
      integer(4),allocatable :: ngvecpB(:,:,:),ngveccB(:,:) !,ngveccB(:,:,:)
     &           , ngvecp(:,:), ngvecc(:,:), !,ngpn(:),ngcni(:),iqib(:),
     &   nx(:,:),nblocha(:),lx(:),ifppb(:)   !ongveccBr(:,:,:)
      complex(8),allocatable:: geigB(:,:,:,:) ,geig(:,:),vcoul(:,:),
     &  zw(:,:),zw0(:,:),
     &  zxq(:,:,:),zxqi(:,:,:)
      real(8),allocatable :: ppbrd (:,:,:,:,:,:,:),cgr(:,:,:,:),eqt(:),
     &  ppbrdx(:,:,:,:,:,:,:),aaa(:,:),symope(:,:),
     &  ppb(:,:),pdb(:,:),dpb(:,:),ddb(:,:), 
     &  qbze(:,:),qibze(:,:) !,ecore(:,:)!,rcxq(:,:,:) !sf..13June
c     &  freqr(:),freqi(:) !rw(:,:),cw(:,:) --->zw
      complex(8),allocatable :: trwv(:),trwv2(:),rcxq(:,:,:,:),rcxqmean(:,:) !sf..13June

c  tetra 
      logical :: tetra,tmpwwk=.true.! If tmpwwk=.true., this use a temporary file tmp.wwk
                                    ! so as to reduce the memory usage.
      complex(8) :: fff,img=(0d0,1d0)
      complex(8),allocatable :: wwk(:,:,:)
      integer(4),allocatable :: 
     &         n1b(:,:,:),n2b(:,:,:),nbnb(:,:),nbnbtt(:,:),noccxvv(:)
c     &         idtetf(:,:),ib1bz(:) 
c 
      real(8),allocatable:: ekxx1(:,:),ekxx2(:,:)!,ekxx3(:,:) !,qbzw(:,: )
      real(8) ::qbzx(3),anfvec(3)
      logical :: anfexist 
      logical :: debug
      integer(4),allocatable:: ibasf(:)
      real(8),allocatable :: transaf(:,:)

      real(8),allocatable :: freq_r(:),freq_i(:),freqx(:),wx(:),expa(:)
     &                      ,frhis(:)
      logical :: realomega, imagomega
c      real(8),allocatable:: ppovl(:,:)
      complex(8),allocatable:: epsi(:,:),gbvec(:),zzr(:),x0mean(:)
      complex(8) :: epxxx,vcmean, vcmmmm
      complex(8),allocatable:: vcmmm(:)
      character*9 fileps
      character*11 fileps23
      character*15 filepsnolfc
      character*11  filele
      logical :: paralellx0=.true. !, hist
c      character(5) :: charnum5
c      character(8) :: xt
      character(20):: xxt

c faleev 
      real(8) :: Emin, Emax,emin2,emax2          
c      integer :: nw2
      real(8) :: omg2max,omg1max
      real(8), allocatable :: freqr2(:)  , ekxxx(:,:,:)
ctm, 080213
      real(8) :: emax_freqr,emin_freqr

c for tetwt5
      real(4),allocatable :: demin(:,:,:,:),demax(:,:,:,:)
      real(8),allocatable :: whw(:)

      integer(4),allocatable:: ihw(:,:,:),nhw(:,:,:),jhw(:,:,:),ibjb(:,:,:,:)
      
      logical::imagonly=.false.,realonly=.false. !,readgwinput
      integer(4)::iopen,maxocc2,iclose,
     & ixc,iqxini,iqxend,
c,nw,incwf,natom,nclass,ipos,igrp,
c     & iinvg,nspin,nl,nn,nnv,nnc,
c     o                   inindx,inindxv,inindxc,iiclass,             !l,n, dimensions 
c     d                   nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, !l,n,  dimensions
c     o                   izdummy, 
c     o   iil,iin,iim,   iilnm, i_mnl, ! l,n,m for Phi ! w(i_mnl)=> mnl(ic) for all electron
c     o   iilv,iinv,iimv,iilnmv,i_mnlv,! l,n,m for Phi
c     o   iilc,iinc,iimc,iilnmc,i_mnlc,! l,n,m for Phi
c     o   iecore,ikonf,iicore,incore,nctot,             !core
c     o   imagw_dummy,niw,idummy,
     &   nw_input,
     &   ifhbe,
     &   nprecb,mrecb,mrece,nlmtot,nqbzt,nband,
     &   nq0i,i,nq0ix,neps,ngrpmx,ngcmx,mxx,nqbze,nqibze,ini,ix,ngrpx
     &  ,mdimx,nbloch,nblochpmx,ifvcfpout,ndummy1,ndummy2,ifcphi,is,nwp,
     &   ifepscond,nxx,ifvxcpout,ifgb0vec
     &   ,nw0,iw,nwhis,ifinin,nw2,iw0,ifwwk,noccxv,noccx
     &   ,ifemesh,nprecx,mrecl,ifwd,ifrcwi,ifrcw,nspinmx,ifianf,ibas
     &   ,ibas1,irot,iq,ngb,iqixc2,ifepsdatnolfc,ifepsdat,ngbin,igc0
     &   ,kx,isf,kqxx,kp,job,nbnbx,nhwtot,noccxvx(2)=-9999,nwmax  !,ifev1,ifev2
     &   ,ihis,jhwtot,ik,ibib,ib1,ib2,ichkhis,ihww,j,imode
     &   ,ngpmx,  ifchipmlog

      real(8):: dum1,dum2,dum3,wqtsum,epsrng,dnorm,
     & dwry,dwh,omg_c,omg2

      integer(4)::nwin, incwfin,  verbose
      real(8)::efin

c for ixc==13
c      integer(4) ::iepsmode=0
c      logical :: ifgeigb

      integer(4):: bzcase, ngc,mrecg
      real(8):: qq(3),quu(3), deltaq(3),qqq(3)=0d0
      logical:: omitqbz=.false., noq0p

      logical,allocatable :: iwgt(:,:,:,:)
      complex(8),allocatable:: wgt(:,:,:)
      
      real(8),allocatable:: qbz2(:,:)
      logical :: qbzreg
      integer(4):: nbcut,nbcut2

      integer(4),allocatable:: nstibz(:) !Nov2004 Miyake's tote
      real(8),allocatable:: ecqw(:,:),wiw(:)
      real(8) :: erpaqw, trpvqw, trlogqw,rydberg,hartree
     &   ,pi,efz,qfermi,alpha,rs,voltot,ecelgas,efx,valn
      integer(4):: iqbz,iqindx,iflegas,nmx
     &   ,ifcor,nqitot,isx,ntot,ieclog,iww,iqq,ieceig,ecorr_on=-1
      real(8) :: eclda_bh,eclda_pz,wk4ec,faca
      logical ::legas

      real(8),allocatable::    evall(:)
      complex(8),allocatable:: ovlpc(:,:),evecc(:,:)
      integer(4):: nev !,  ifdpin


      real(8),allocatable:: ecut(:),ecuts(:) ,totexc(:), trpv(:),trlog(:)
      integer(4):: necut,iecut

c mode23
      integer(4):: ifv,lxx,ibasx,ilmx,ilm_r,nx_r,lb,nb,mb
      integer(4),allocatable:: nxx_r(:)
      real(8),allocatable:: svec(:),spinvec(:,:)
      character*3:: charnum3
      character*4:: charnum4
      complex(8),allocatable:: jcoup(:,:), mcm(:)
      complex(8):: aaax,aaax0,addx !bbbx,bbbx25,bbbx35
      real(8)::chg1,chg2,spinmom,schi
#ifdef COMMONLL
      integer(4)::ll(51**2)
      common/llblock/ll
#else
      integer(4) :: ll
      external ll
#endif
      complex(8),allocatable:: ovlp(:,:),evec(:,:),ovlpi(:,:)
      real(8),allocatable::eval(:)
      integer(4):: new,nmxx,ii,iy,ipl1,ixx

      complex(8),allocatable :: ppovl(:,:),oo(:,:)
      real(8)::qxx(3),mmnorm,ssm
c svd
      real(8),allocatable::SS(:),rwork(:),ss0(:)
      complex(8),allocatable:: UU(:,:),VT(:,:),work(:),zw0bk(:,:)
     & ,vtt(:,:),zzz(:,:),sqsvec(:),ooo(:,:),ppo(:,:),sqovlp(:,:),sqovlpi(:,:)
      integer(4)::lwork,info,imin,ifzxq
      complex(8)::ddd, x0meanx
      complex(8),allocatable:: UU0(:,:),VT0(:,:)

      logical :: ngczero=.false., sergeyv=.false., chipm=.false.,nolfco=.false.
     &  ,epsmode=.false.,normalm=.false., eiqr=.false.
      integer(4)::  ife, ifchipmn,ifchipm,isvd,idum4 !ingczero,
      real(8):: qs,qt,ww,muu, ddq(3)
      character*11 ::ttt

      logical ::smbasis
      integer(4):: nnmx,nomx

C... Feb2006 time-reversal=off case
      logical :: timereversal, testtimer,onceww
      integer(4):: npm ,jpm,ncc,nw_i
      real(8):: frr

ctm, 060517
      logical :: lmlw=.false.
      integer(4):: nwf,ifwd_full,ifrcw_full,ifrcwi_full,
     &             nprecx2,mrecl2,nblochpmx2,nwp2,niw2,nqibze2,
     &             nw_i2,nrec
      complex(8),allocatable:: zw_full(:,:),zw_full2(:,:)
ctm, 080311
      logical :: lwbnd=.false.
      integer(4),allocatable:: iwbnd(:,:,:)
ctm, 080721
      logical :: lent=.false.
      integer(4),allocatable:: ib_ent(:,:,:)

      hartree=2d0*rydberg()
      pi         = 4d0*datan(1d0)
!-------------------------------------------------------------------------
C --- computational mode select
      write(6,"(a)") 
     &      '--- Type numbers #1 #2 #3 [#2 and #3 are options] ---'  
      write(6,"(a)") ' #1:run mode'
      write(6,"(a)") '    1   : normal' 
      write(6,"(a)") '    2   : epsNolfc'
      write(6,"(a)") '    3   : eps'
      write(6,"(a)") '    22  : chi^+- NoLFC'
      write(6,"(a)") '    23  : chi^+- '
      write(6,"(a)") '    101 : total energy Miyake'
      write(6,"(a)") '    111 : maxloc-W     Miyake'
      write(6,"(a)") '    11  : normal    Sergey'
      write(6,"(a)") '    202 : eps       Sergey'
      write(6,"(a)") '    203 : eps NoLFC Sergey'
      write(6,"(a)") '    222 : chi^+- NoLFC Sergey'
      write(6,"(a)") '    223 : chi^+- Sergey'
      write(6,"(a)") '    12  : total energy Miyake Sergey'
      write(6,"(a)") '    -9999: just show version num'
      write(6,"(a)") ' #2=iqxini   #3=iqxend' !    #4=ngczero or not'
      write(6,"(a)") '  10222 : <e^{iqr}|chi^+-|e^{iqr}> NoLFC'
      write(6,"(a)") 
     &      '-------------------------------------------------------'
c Test mode   write(6,"(a)") '   tettest(4) Realtetwt4(5) Realtetwt5(6) '
      call readin6(ixc,iqxini,iqxend, idum4 ) !ingczero)
      call headver('hx0fp0',ixc)
c      if(ingczero==1) then
c        ngczero=.true.
c      endif
c      print *,' ngczero=',ngczero ! print *,' ixc,iqxini,iqxend, ingczero=', ixc,iqxini,iqxend, ingczero
      if(ixc==0) stop ' --- ixc=0 --- Choose computational mode!'
      if(iqxini ==0 ) then
        paralellx0=.false.
        write(6,"(' ixc = ',i4)") ixc
      else
        paralellx0=.true.
        write(6,"(' PARALELL mode: ixc iqxini iqxend =',3i5)")
     &  ixc, iqxini, iqxend
      endif

! original mode
      if    (ixc==1) then
        write(6,*) ' OK ixc=1. normal mode'
        normalm=.true.
!     -- eps modes
      elseif(ixc==2) then 
        write(6,*) ' OK ixc=2. eps mode. ONLY NoLFC' 
        omitqbz=.true.
      elseif(ixc==3) then
        write(6,*) ' OK ixc=3. eps mode'
        omitqbz=.true.
      elseif(ixc==22)then
        write(6,*) ' OK ixc=22. (chi^{+-})^{-1} NoLFC mode' 
        omitqbz=.true.
      elseif(ixc==23)then
        write(6,*) ' OK ixc=23. (chi^{+-})^{-1} mode' !May2005
        omitqbz=.true.
! test modes
      elseif(ixc==4) then
         write(6,*)            ' OK ixc=4. test. mode tetwt4.vs.tetwt5.'
      elseif(ixc==5) then
         write(6,*) ' OK ixc=5. Real-axis tetwt4 mode'

!total energy mode
      elseif(ixc==101) then
        write(6,*)' OK ixc=101. Miyake x0 mode'
        ecorr_on=901

! test mode
      elseif(ixc==6) then
         sergeyv=.true.
         write(6,*) ' OK ixc=6. Real-axis tetwt5 mode'
! Sergey (Hilbert-transformation) modes
      elseif(ixc==11) then; write(6,*) " OK ixc=11 s.f.'s normal mode " 
        sergeyv=.true. !ixc=11
        normalm=.true.
!     -- eps mode NoLFC
      elseif(ixc==202) then
        write(6,*) " OK ixc=202  sergey's eps mode Only NoLFC "
        sergeyv=.true. !ixc=11
        realonly=.true.
c        iepsmode=202
        omitqbz=.true.
!     -- eps mode with LFC 
      elseif(ixc==203) then
        write(6,*) " OK ixc=203 sergey's eps mode with LFC "
        sergeyv=.true. !ixc=11
        realonly=.true.
c        iepsmode=203
        omitqbz=.true.
! Total energy modes
      elseif(ixc==12) then 
        write(6,*) " OK ixc=12 Miyake's total energy Sergey"
        sergeyv=.true. !ixc=11
        imagonly=.true.
        ecorr_on=901
ctm, 060517
      elseif(ixc==111) then
        write(6,*) " OK ixc=111 Miyake's maxloc mode"
        sergeyv=.true. !ixc=11
        normalm=.true.
        lmlw = .true.
!     -- chipm mode NoLFC
      elseif(ixc==222) then
        write(6,*) " OK ixc=222    chipm sergey's "
        sergeyv=.true. 
        realonly=.true.
        omitqbz=.true.
!     -- chipm mode NoLFC
      elseif(ixc==10222) then
        write(6,*) " OK ixc=10222  <q|chipm|q> sergey"
        sergeyv=.true. 
        realonly=.true.
        omitqbz=.true.
        eiqr =.true.
!     -- eps mode with LFC 
      elseif(ixc==223) then
        write(6,*) " OK ixc=223  chipm with LFC sergey's  "
        sergeyv=.true. 
        realonly=.true.
        omitqbz=.true.
      else             
        stop' hx0fp0: mode ixc is not appropriate'
      endif
!
      if(    ixc==  2.or.ixc==  3.or.ixc== 22.or.ixc== 23
     &   .or.ixc==202.or.ixc==203.or.ixc==222.or.ixc==223
     &   .or.ixc==10222) then
        epsmode=.true. !        iepsmode = ixc
      endif

      if(epsmode) then
        if(mod(ixc,200)==22.or.mod(ixc,200)==23) chipm =.true.
        if(mod(ixc,10)==2)                       nolfco=.true.
        if(chipm) ifchipmlog = iopen ('ChiPM.log',1,3,0)
      endif
c
c      call getkeyvalue("GWinput","Ecorr",ecorr_on,default=0)
      if(ecorr_on > 0) then
         ieclog = 8155
         if(ecorr_on==901) then
           ieceig=8156
           open(ieceig,file='rpa_eigen.chk') 
           close(ieceig,status='delete')
         endif
         open(ieclog, file='ecorr.chk')
      endif
c        ifdpin = iopen('DPSIONin',0,-1,0) 

c---  readin BZDATA. See gwsrc/rwbzdata.f
c--------readin data set when you call read_BZDATA ---------------
c       integer(4)::ngrp,nqbz,nqibz,nqbzw,nteti,ntetf,
c     &   n_index_qbz  
c       integer(4):: n1,n2,n3
c       real(8):: qbas(3,3),ginv(3,3),qbasmc(3,3),dq_bzcase2
c       real(8),allocatable:: qbz(:,:),wbz(:),qibz(:,:)
c     &    ,wibz(:),qbzw(:,:)
c       integer(4),allocatable:: idtetf(:,:),ib1bz(:),idteti(:,:)
c     &    ,nstar(:),irk(:,:),index_qbz(:,:,:),nstbz(:)
c-----------------------------------------------------------------
      call read_BZDATA()

c--- Use regular mesh even for bzcase==2 and qbzreg()=T
      ddq = 0d0
      if(bzcase()==2) ddq= dq_bzcase2
      do iq = 1, nqbz
        qbz(1:3,iq) = qbz(1:3,iq) + ddq  
        ! This new qbz is regular mesh, which are identical in the both bzcase. 
      enddo
      if(.not.qbzreg()) then ! off-regular mesh case
        do i=1,nqbz
          qbz(:,i) = qbz(:,i) - dq_bzcase2
        enddo
      endif
      do i=1,nqbz
        write(6,"('i qbz=',i3,3f8.4)") i,qbz(:,i)
      enddo
      print *,' nqbz nqibz =',nqbz,nqibz

c      print *,  qbz 
c      print *,' irk=',irk
c      print *,' #### idtetf: ####'
c      print *, idtetf

C- oct2005 not implimented cases.
      if(smbasis()) then
        if(.not.sergeyv) 
     &  stop ' smbasis=T & sergeyv=F is not implimented yet'
        if(chipm) then
        print *, ' smbasis=T & chipm=T is not implimented yet.'//
     &       ' Supply consistent MixSpin for smbasis!'//
     &       ' MixSpin should be converted at the end of hvccfp0.'
        stop ' smbasis=T & chipm=T is not implimented yet.'
        endif
      endif

C--- readin GWIN and LMTO, then allocate and set datas.
c      call wkinit (iwksize)
      if(epsmode.and.sergeyv) then
         nwin = -999
      elseif(ixc==101.or.ixc==1.or.sergeyv) then
        nwin = 0      !Readin nw from NW file
      else
        nwin = -999   !Not readin NW file
      endif
c      if(iepsmode==202.or.iepsmode==203) nwin = -999

      incwfin= 0  !use ForX0 for core in GWIN
      efin =  0d0 !readin EFERMI
      call genallcf_v3(nwin,efin,incwfin) !in module m_genallcf_v3
      if(ngrp/= ngrp2) stop 'ngrp inconsistent: BZDATA and LMTO GWIN_V2'
c---  These are allocated and setted by genallcf_v3
c      integer(4)::  nclass,natom,nspin,nl,nn,nnv,nnc, ngrp,
c     o  nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, nctot,niw,nw
c      real(8) :: alat,ef, diw,dw,delta,deltaw,esmr   
c      character(120):: symgrp
c      character(6),allocatable :: clabl(:)
c      integer(4),allocatable:: iclass(:)
c     &  ,nindxv(:,:),nindxc(:,:),ncwf(:,:,:) ,
c     o    invg(:), il(:,:), in(:,:), im(:,:),   ilnm(:),  nlnm(:),
c     o    ilv(:),inv(:),imv(:),  ilnmv(:), nlnmv(:),
c     o    ilc(:),inc(:),imc(:),  ilnmc(:), nlnmc(:),
c     o    nindx(:,:),konf(:,:),icore(:,:),ncore(:),
c     &    occv(:,:,:),unoccv(:,:,:)   
c     &   ,occc(:,:,:),unoccc(:,:,:),
c     o    nocc(:,:,:),nunocc(:,:,:)
c      real(8), allocatable:: 
c     o  plat(:,:),pos(:,:),z(:),  ecore(:,:), freq(:), symgg(:,:,:) ! symgg=w(igrp)
c-----------------------------------------------------------------------
c      call genallcf_v2x(
c     o                   plat,alat,natom,nclass,ipos,   !structure
c     o                   symgrp,gen,ngnmax,ngrp,igrp,   !symmetry
c     o                   iinvg,ef,
c     o                   clabl,nspin,nl,nn,nnv,nnc,                  !l,n, dimensions
c     o                   inindx,inindxv,inindxc,iiclass,             !l,n, dimensions 
c     d                   nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, !l,n,  dimensions
c     o                   izdummy, 
c     o   iil,  iin,iim,   iilnm,i_mnl, ! l,n,m for Phi ! w(i_mnl)=> mnl(ic) for all electron
c     o   iilv,iinv,iimv,iilnmv, i_mnlv,! l,n,m for Phi
c     o   iilc,iinc,iimc,iilnmc, i_mnlc,! l,n,m for Phi
c     o   incwf,iecore,ikonf,iicore,incore,nctot,            !core
c     o   imagw_dummy,niw,dum1,nw,dw,delta,dum2,dum3,idummy) ! frequency

ctm, 060517
        if (niw==0) realonly = .true.
c ...
      if(chipm.and.nspin==1) stop 'chipm mode is for nspin=2'
      nw_input = nw ; print *, 'nw delta=',nw_input,delta
cc   if (nclass > mxclass) stop ' hx0fp0: increase mxclass'
      debug=.false.; if(verbose()>=100) debug=.true.
      if(debug) print *,' end of genallc'
!!!! WE ASSUME iclass(iatom)= iatom !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  We assume nclass = natom.     !!!!!!!!!!!!!!!!!!!!!!!!!
      if(nclass /= natom) stop ' nclass /= natom '

c ... tetra or not
      if(delta <= 0d0) then
        tetra =  .true.
        delta = -delta
        print *,' hx0fp0: tetrahedron mode delta=',delta
      else
        tetra = .false. ! switch for tetrahedron method for dielectric functions
      endif

c ... read dimensions of h,hb
      ifhbe      = iopen('hbe.d',1,0,0)
      read (ifhbe,*) nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg
      if(nlmto/=nlmtot) stop ' hx0fp0: nlmto/=nlmtot in hbe.d'
      if(nqbz /=nqbzt ) stop ' hx0fp0: nqbz /=nqbzt  in hbe.d'

c --- q near zero --------
      if(debug) print *, 'reading QOP'
      open (101,file='Q0P')
      read (101,"(i5)") nq0i
      write(6,*) ' *** nqibz nq0i=', nqibz,nq0i
      allocate( wqt(1:nq0i),q0i(1:3,1:nq0i) )
      do i=1,nq0i
        read (101, * ) wqt(i),q0i(1:3,i)
      enddo
      nq0ix = nq0i
      do i=1,nq0i
        if(wqt(i)==0d0 ) then
          nq0ix = i-1 
          exit
        endif
      enddo
      neps = nq0i - nq0ix  ! number of zero weight q0p which are used for ixc=2 or 3 mode.
      write( 6,*) ' num of zero weight q0p=',neps
      write(6,"(i3,f14.6,2x, 3f14.6)" )(i, wqt(i),q0i(1:3,i),i=1,nq0i)
      close(101)
      wqtsum = sum(abs(wqt(1:nq0i))) 
      call getkeyvalue("GWinput","TestNoQ0P",noq0p,default=.false.)
      if(normalm.and.abs(wqtsum-1d0) >1d-10.and.(.not.noq0p)) 
     &  stop ' wqtsum of Q0P /=1 '

c --- readin by rdpp ; Radial integrals ppbrd and plane wave part
      call getsrdpp2( nclass,nl,nxx)
      call readngmx('QGpsi',ngpmx)
      call readngmx('QGcou',ngcmx)
      print *,' ngcmx ngpmx=',ngcmx,ngpmx
! --- qibze(3,nqbze) qbze(3,nqibze)
      nqbze  = nqbz *(1 + nq0i)
      nqibze = nqibz + nq0i
      allocate( qbze(3, nqbze), qibze(3, nqibze))
      call dcopy(3*nqbz, qbz,  1, qbze,1)
      call dcopy(3*nqibz,qibz, 1, qibze,1)
      do i = 1,nq0i
        qibze(:,nqibz+i)  = q0i(:,i) 
        ini = nqbz*(1 + i -1)
        do ix=1,nqbz
          qbze (:,ini+ix)   = q0i(:,i) + qbze(:,ix)
        enddo
      enddo
      ngrpx = 1
      allocate( nx(0:2*(nl-1),nclass),  nblocha(nclass) ,lx(nclass),
     &   ppbrd ( 0:nl-1, nn, 0:nl-1,nn, 0:2*(nl-1),nxx, nspin*nclass),
     &   cgr(nl**2,nl**2,(2*nl-1)**2,ngrpx), symope(3,3))
      symope(1:3,1) = (/1d0,0d0,0d0/)
      symope(1:3,2) = (/0d0,1d0,0d0/)
      symope(1:3,3) = (/0d0,0d0,1d0/)
      call rdpp_v3(nxx, nl,ngrpx, nn, nclass, nspin,symope,qbas,
     o      nblocha, lx, nx,  ppbrd , mdimx, nbloch, cgr)
      if(.not.smbasis()) nblochpmx = nbloch + ngcmx 
      allocate(ngveccB(3,ngcmx)) ! work arry 
     
c... for legas test 
      legas = .false.
      INQUIRE (FILE = 'LEGAS', EXIST = legas)

c --- Check PARALELL.X0
      if(chipm.and.(.not.legas).and.(.not.eiqr) ) then
        iqxend = nqibz + nq0i
      elseif(paralellx0) then
        if(iqxini >99999) stop ' hx0fp0: iqxini >99999'
        if(iqxend >99999) stop ' hx0fp0: iqxend >99999'
        ifvcfpout = iopen('VCCFP.'
     &     //xxt(iqxini,iqxend)
     &  ,0,-1,0) 
        read(ifvcfpout) ndummy1, ndummy2
      else
        iqxend = nqibz + nq0i
        ifvcfpout = iopen('VCCFP',0,-1,0) 
        read(ifvcfpout) ndummy1, ndummy2
      endif
      print *, ' nqibz nqibze=',nqibz,nqibze

      call init_readeigen(ginv,nspin,nband,mrece) !initialization of readEigen
                                                  ! EVU EVD are readin in init_readeigen
      call init_readeigen2(mrecb,nlmto,mrecg)

c --- set realomega, imagomega tetra nw niw nwp ifgb0vec ------------
      if( ixc==101 ) then
        realomega =.false.
        imagomega =.true.
        nw = 0
        nwp= 0
      elseif( ixc==1 ) then !old imagw = 2 case
        realomega =.true.
        imagomega =.true.
      elseif( epsmode.and.(.not.sergeyv)) then !ixc==2.or.ixc==3.or.ixc==22.or.ixc==23 ) then
        print *, 'mode=',ixc
        realomega =.true.
        imagomega =.false.
        niw = 0
        call getkeyvalue("GWinput","EPSrange",epsrng )
        call getkeyvalue("GWinput","EPSdw",  dwry )
        epsrng=epsrng/2d0
        dw = dwry/2d0
        if(dw==0d0) then
          nw = 1
        else
          nw = (epsrng/2d0 - 1d-10)/(dw/2d0) + 2 !epsrng/2d0 corresponds to in a.u.
        endif
      elseif(ixc==4.or.ixc==5.or.ixc==6) then 
! ... These are test modes.
! ixc=4 tetrahedren weight test. tetwt5.vs.tetwt5. Write tethis.chk
! ixc=5 Spectrum function (Img part) along the Real axis with tetwt4
! ixc=6 Spectrum function (Img part) along the Real axis with tetwt5. Histgram method.
        realomega = .true.
        imagomega = .false.
        tetra     = .true.
        niw = 0 
! ---  For tetwt5 ---  the tetrahedron weight for spectrum function (imaginary part)
!   Histogram bins are specified by freq_r(1:nwp)
!     nwp=nw+1; frhis(1)=0
!     The 1st   bin  is     [frhis(1),  frhis(2)]   ...
!     The last  bin  is     [frhis(nw), frhis(nwp)].

! ... These parameters specifies a test histgram bins;Sergey's mesh just for test modes.
        nw0 = 200    !100    800
        dwh = 0.01d0 !0.02d0 0.0025d0 !in hartree 
        call findemaxmin(nband,qbze,nqbze,nspin,emax,emin)
c$$$ Even if we neglect this section, result is not so bad---anyway we may need to 
c$$$ fix this block if necessary. Oct2005
c$$$        if(bzcase()==2) then
c$$$          allocate(qbz2(3,nqbz))
c$$$          do iq=1,nqbz
c$$$           qbz2(:,iq)=qbz(:,iq)+dq_bzcase2
c$$$          enddo
c$$$          call findemaxmin(nband,qbz2,nqbz,nspin ,emax2,emin2)
c$$$          emax=max(emax,emax2)
c$$$          emin=min(emin,emin2)
c$$$          deallocate(qbz2)
c$$$        endif
        if (nctot .gt. 0) Emin = minval(ecore(:,1:nspin)) 
ctm, 080213 >
        call getkeyvalue("GWinput","emax_freqr",emax_freqr,default=-999d0)
        call getkeyvalue("GWinput","emin_freqr",emin_freqr,default=999d0)
        if (emax_freqr.gt.emin_freqr) then
           if (emax_freqr.gt.emax) emax = emax_freqr
           if (emin_freqr.lt.emin) emin = emin_freqr
           write(*,*)'hx0fp0: emax/emin from GWinput'
           write(*,*)'emax=',emax
           write(*,*)'emin=',emin
        endif   
ctm, 080213 <
        omg2max = (Emax-Emin)*.5+.2d0    !(in Hartree) covers  all relevant omega, +.2 for margin
        omg1max = dwh*(nw0-1)
        nwp = int(sqrt(omg2max*(2*nw0-1d0)/dwh-(nw0**2-3*nw0+1d0)))+1 ! + 1 for margin
        nw  = nwp-1
        write(*,*) Emax,Emin,nw0,nw  ! nwp is new max number in frequency array
        write(6,'(a32,2i7,2d15.3)')'hx0fp0: nw0,nw,omg1max,omg2max='
     &              , nw0,nw,  omg1max,omg2max
        if (nw <= nw0) stop 'hx0fp0:ixc==[456] nw2 <= nw' 
        allocate(freq_r(0:nw))
        do iw=0,nw  !This is a test mesh by Sergey.Faleev
          if(iw+1<=nw0) then;  freq_r(iw)=dwh*iw
          else;  freq_r(iw)=dwh*((iw+1)**2+nw0**2-3*nw0+1)/(2*nw0-1d0)
          endif  
        enddo !freq_r(iw) is linear for iw<=nw and quadratic for nw<iw<=nw2
               !freq_r(iw) chosen in such a way that it is continues with 
!!! nw nwp=nw+1 freq_r(1:nwp) are used after here.
        allocate(frhis(nwp))
        frhis=freq_r(0:nw)
        nwhis=nw

C ... Sergey Faleev's normal  mode
      elseif(sergeyv) then 
        realomega = .true.
        imagomega = .true.
        tetra     = .true.
        call getkeyvalue("GWinput","omg_c",omg_c )
        write(*,*)'dw, omg_c= ',dw, omg_c
C ... frhis        
        if(bzcase()==2) then
          write(6,'(" dq_bzcase2=",3f9.4)')dq_bzcase2
        endif
c        call findemaxmin(ifev,nband,nqbz,nspin
        call findemaxmin(nband,qbze,nqbze,nspin
     o  ,emax,emin)
c$$$ Even if we neglect this section, result is not so bad---anyway we may need to 
c$$$ fix this block if necessary. Oct2005
c$$$        if(bzcase()==2) then
c$$$          allocate(qbz2(3,nqbz))
c$$$          do iq=1,nqbz
c$$$           qbz2(:,iq)=qbz(:,iq)+dq_bzcase2
c$$$          enddo
c$$$          call findemaxmin(nband,qbz2,nqbz,nspin ,emax2,emin2)
c$$$          emax=max(emax,emax2)
c$$$          emin=min(emin,emin2)
c$$$          deallocate(qbz2)
c$$$        endif
        if (nctot > 0) Emin=minval(ecore(:,1:nspin)) 
        omg2max=(Emax-Emin)*.5d0+.2d0 !(in Hartree) covers  all relevant omega, +.2 for margin
        print *,' emin emax omg2max=',emin, emax,omg2max
        nw2=int(omg_c/dw*( sqrt(1.+2*omg2max/omg_c)-1. ) )+1+3 !+3 for margin
        allocate(freqr2(nw2))                                  !+1 b/c (iw-1)
        do iw=1,nw2
          freqr2(iw)=dw*(iw-1)+dw**2/2./omg_c*(iw-1)**2
        enddo !linear + quadratic term
        if (nw2 < 2 ) stop "hx0fp0: nw2 < 2"
!        if (dw*(nw-2) > freqr2(nw2-1)) stop
!     &                           "hx0fp0: dw*(nw-2) > freqr2(nw2-1)"
        if (dw*(nw-2) > freqr2(nw2-1)) then
          write(*,*)'dw=',dw 
          write(*,*)'nw=',nw 
          write(*,*)'nw2=',nw2
          write(*,*)'omg_c=',omg_c
          write(*,*)'omg2max,',omg2max
          write(*,*)'freqr2(nw2-1)',freqr2(nw2-1)
          write(*,*)'recommended nw<',int(freqr2(nw2-1)/dw) +1
!          do iw=1,nw2
!            write(*,*)iw,freqr2(iw)
!          enddo
          stop "hx0fp0: dw*(nw-2) > freqr2(nw2-1)"
        endif  
        nwhis = nw2-1
        allocate(frhis(1:nwhis+1)) 
        frhis = freqr2(1:nwhis+1)
        print *,' we set frhis nwhis=',nwhis
C ... for frhis_m
        nw=nw2-1   ! nw+1 is how many points of real omega we use 
                   ! for dressed coulomb line W(iw=0:nw) iw=0 corresponds omg=0
                   ! maximum nw=nw2-1 because nwhis=nw2-1 
        do iw=3,nw2-1  
        !nw is chosen from condition that frhis_m(nw-3)<dw*(nw-3) <frhis_m(nw-2).
        !Here frhis_m(iw)= (freqr2(iw)+freqr2(iw+1))/2d0   
        !nw was constructed such that omg=dw*(nw-2)> all relevant frequensies needed
        ! for correlation Coulomb Wc(omg),
        ! and one more point omg=dw*(nw-1) needed for extrapolation.    
        ! Now, frhis_m(nw-1)> all relevent frequensies for Wc(omg)
        ! and one more point omg=frhis_m(nw) needed for extropolation 
        ! used in subroutine alagr3z in  sxcf.f.
        omg2 = (freqr2(iw-2)+freqr2(iw-1))/2d0
         if (omg2 > dw*(nw_input-3)) then
              nw=iw  
              exit   ! 'nw_input' is only used to get maximum frequency for
                     !  dressed coulomb line
          endif                    
        enddo   
        if(epsmode.and.sergeyv) then !epsmode==202.or.iepsmode==203) then
          nw  = nwhis-1
          niw = 0
        endif  
        allocate(freq_r(0:nw))
        freq_r(0)=0d0
        do iw=1,nw
          freq_r(iw)=(frhis(iw)+frhis(iw+1))/2d0
        enddo
        nwp = nw+1
        delta=0d0 !only for writeemesh
      endif

C ... timereversal=F is implimented only for tetra=T and sergeyv=T
      npm=1
      nw_i=0
      if(.not.timereversal()) then
        write(6,"('TimeReversal off mode')")
        npm=2
        nw_i=-nw
        if(.not.tetra)      stop ' tetra=T for timereversal=off'
        if(.not.sergeyv)    stop ' sergeyv=T for timereversal=off'
        if(chipm)           stop ' chipm=F for timereversal=off'
      endif

! ...
      if(realomega.and.(.not.sergeyv)) then !     elseif (realomega.and.(ixc==22.or.ixc==23)) then
        nwp = nw+1 !Feb2006 I want to unify the meaning of nwp as nwp=nw+1
                   !With this change x0kf_v2z calculate one more freq_r along real axis.
                   ! at freq_r(nw). --->It WAS not calculated. It wat up to freq_r(nw-1).
                   ! But not tested well.
        allocate( freq_r(0:nw) ) 
        do iw = 0,nw
          freq_r(iw) = dw*iw ! freq_r
        enddo
      endif

!...get eigenvector corresponds to exp(iqr) (q is almost zero).
      if(epsmode) then !iepsmode/=0) then ;  print *, ' read in Mix0vec'
        allocate(epsi(nw_i:nw,neps)) !5July2005 nwp should be used after it is defined!
        if(paralellx0) then
          ifgb0vec = iopen ( "Mix0vec."//xxt(iqxini,iqxend),1,3,0)
        else
          ifgb0vec = iopen ( "Mix0vec",1,3,0) 
        endif
        print *, ' read end of Mix0vec'
      endif
! ...
      if(realomega) then
        open(UNIT=3111,file='freq_r') !write number of frequency
        !points nwp and frequensies in 'freq_r' file 
        write(3111,"(2i8,'  !(a.u.=2Ry)')") nwp, nw_i
        do iw= nw_i,-1
          write(3111,"(d23.15,2x,i6)") -freq_r(-iw),iw !This file is reffere by hsfp0 and so.
        enddo                     
        do iw= 0,nw
          write(3111,"(d23.15,2x,i6)") freq_r(iw),iw !This file is reffere by hsfp0 and so.
        enddo                     
        close(3111)
      endif  

C --- WVI only for imagonly for ixc==12
      if(imagonly.and.sergeyv) then
        realomega =.false.
        imagomega =.true.
      endif

C --- epsPP noLFC mode for ixc==13
      if(realonly.and.sergeyv) then
        realomega =.true.
        imagomega =.false.
      endif

C --- set freq_i 
      if (imagomega) then
        print *,' freqimg: niw =',niw
        allocate( freq_i(niw) ,freqx(niw),wx(niw),expa(niw) )                       
        call freq01 (niw,ua,  ! Generate gaussian frequencies x between (0,1) and w=(1-x)/x      
     o                freqx,freq_i,wx,expa)  
         allocate(wiw(niw))
         do iw=1,niw
           wiw(iw)=wx(iw)/(2d0*pi*freqx(iw)*freqx(iw))
         enddo
c        deallocate(freqx,wx,expa)! freqx,wx is used for Miyake's tote mode !Nov2004                       
      endif

c ... for Miyake tote mode Nov2004
      if(ecorr_on>0) then  !it was bzcase()==2 Was it bug?
        allocate(nstibz(nqibz))
        do iq=1,nqibz
          iqbz = iqindx(qibz(:,iq),ginv,qbz,nqbz)
          nstibz(iq) = nstbz(iqbz)
c          write(6,"(' iq qibz nstibz=',i5,3f9.4,i5)")iq,qibz(:,iq),nstibz(iq)
        enddo
      endif

ctm, 060517
      if (lmlw) then ! ixc=111
        print *,' nband before lmlw=',nband
        call init_readeigen_mlw_diag(nwf,nband,mrecb,mrecg)
        write(*,*)'readeigen_mlw_diag done'
        nctot = 0
ctm, 040919
c        open(3311,file='NOCC_MLW')
c        read(3311,*)
c        read(3311,*)
c        read(3311,*)
c        read(3311,*)ef
c        close(3311)
c        write(*,*)'ef=',ef
        print *,' nband after lmlw=',nband
      else ! lmlw
ctm, 080721
        nwf =  1 ! temporary
        call getkeyvalue("GWinput","crpa_ent",lent,default=.false.)
        if (lent) then 
          write(*,*)'cRPA_ent mode on'
          allocate(ib_ent(nband,nqbze,nspin))
          ib_ent = 0
          call init_readeigen_crpa_ent(ib_ent,nwf)
        else ! lent
ctm, 080311
!          call getkeyvalue("GWinput","wan_band",lwbnd,default=.false.)
        call getkeyvalue("GWinput","crpa_band",lwbnd,default=.false.)
          if (lwbnd) then
            write(*,*)'cRPA_band mode on'
            call get_nwf(nwf)
            allocate(iwbnd(nwf,nqbze,nspin))
            call choose_wanband(iwbnd,nwf,nqbze,nspin)
          else ! lwbnd
            write(*,*)'cRPA_window mode on'
          endif ! lwbnd
        endif ! lent
      endif ! lmlw
      if (.not.lwbnd) allocate(iwbnd(nwf,nqbze,nspin))
      if (.not.lent) allocate(ib_ent(nband,nqbze,nspin))

C --- tetra init
      call getkeyvalue("GWinput","TmpWWK",tmpwwk,default=.false.)
      if(tetra) then
        allocate( !wgt(nband+nctot,nband,nqbz), !noccxvv(nw+niw),
     &    nbnb(nqbz,npm),nbnbtt(nqbz,npm)
     &   ,ekxx1(nband,nqbz),ekxx2(nband,nqbz) !!! nband=nlmto 
     &   )
        if(tmpwwk)  ifwwk = iopen('tmp.wwk',0,-1,0)
      endif ;  if(debug) print *,' xxx1:'
c ... maximum no. occupied valence states
c      noccxv     = maxocc (ifev,nspin,ef,  nband,nqbze)
      noccxv     = maxocc2 (nspin,ef, nband, qbze,nqbze)
      write(*,*)'ef=',ef
ctm, 060517
      if(noccxv>nband) stop 'hx0fp0: all the bands filled! too large Ef'
c      if(noccxv>nband) then
c         write(*,*) 'hx0fp0: all the bands filled!'
c         noccxv = nband
c      endif   
      noccx      = noccxv + nctot

c --- write energy mesh ---
      if(paralellx0) then
        ifemesh = iopen('emesh.hx0fp0.'//xxt(iqxini,iqxend),1,-1,0)
      else
        ifemesh = iopen('emesh.hx0fp0',1,-1,0)
      endif
      call writeemesh(ifemesh,freq_i,niw,freq_r,nwp,delta)
      if(debug) print *,' xxx1aax:'

c ... ppb
      allocate( ppb(nlnmx*nlnmx*mdimx*nclass,nspin) )
c      nwblk  = nw /nwdiv
c      niwblk = niw/nwdiv
c      if ( nwdiv*nwblk  .ne. nw ) nwblk = nwblk + 1
c      if ( nwdiv*niwblk .ne. niw) niwblk = niwblk + 1
c      allocate(  freqr(nwblk), freqi(niwblk)   )
      nprecx = ndble  !We use double precision arrays only.

C... This is just to get nblochpmx
      if(smbasis()) then
      call getngbpomat(nqibz+nq0i,
     o nnmx,nomx)
      nblochpmx = nnmx 
      endif

      mrecl  = nprecx*2*nblochpmx*nblochpmx/nword()

ctm, 060517
      if (lmlw) then

C --- open WV* files
      if(paralellx0) then
        ifwd   = iopen('WV.d.maxloc.'//xxt(iqxini,iqxend),1,-1,0)
        ifwd_full = iopen('WV.d.'//xxt(iqxini,iqxend),1,0,0)
      else
        ifwd   = iopen('WV.d.maxloc',1,-1,0)
        ifwd_full   = iopen('WV.d',1,0,0)
      endif

c! July 2003. now WV.d contains true size of WVR.
      write (ifwd,"(1x,10i14)") !change format oct2005
     &   nprecx,mrecl,nblochpmx,nwp,niw,nqibz + nq0i-1,nw_i
ctm, 060517
      read (ifwd_full,"(1x,10i14)") 
     &   nprecx2,mrecl2,nblochpmx2,nwp2,niw2,nqibze2,nw_i2
      if (nprecx .ne. nprecx2) stop 'hx0fp0: nprecx error'
      if (mrecl .ne. mrecl2) stop 'hx0fp0: mrecl2 error'
      if (nblochpmx .ne. nblochpmx2) stop 'hx0fp0: nblochpmx error'
      if (nwp .ne. nwp2) stop 'hx0fp0: nwp error'
      if (niw .ne. niw2) stop 'hx0fp0: niwecx error'
      if (nqibz+nq0i-1 .ne. nqibze2) stop 'hx0fp0: nqibze error'
      if (nw_i .ne. nw_i2) stop 'hx0fp0: nw_i error'

      allocate(  zw(nblochpmx,nblochpmx) )
      if(debug) print *,' xxx1a2:'
      if    (ixc==101) then
         if(paralellx0) then
           ifrcwi = iopen('WVI.maxloc.'//xxt(iqxini,iqxend) ,0,-1,mrecl)
           ifrcwi_full = iopen('WVI.'
     &     //xxt(iqxini,iqxend) ,0,0,mrecl)
         else
           ifrcwi = iopen('WVI.maxloc.',0,-1,mrecl)
           ifrcwi_full = iopen('WVI',0,0,mrecl)
         endif
      elseif (normalm) then !ixc==1.or.sergeyv) then
         if(paralellx0) then
           ifrcw  = iopen('WVR.maxloc.'//xxt(iqxini,iqxend),0,-1,mrecl)
           ifrcwi = iopen('WVI.maxloc.'//xxt(iqxini,iqxend),0,-1,mrecl)
           ifrcw_full  = iopen('WVR.'//xxt(iqxini,iqxend),0,0,mrecl)
           ifrcwi_full = iopen('WVI.'//xxt(iqxini,iqxend),0,0,mrecl)
         else
           ifrcw  = iopen('WVR.maxloc',0,-1,mrecl)
           ifrcwi = iopen('WVI.maxloc',0,-1,mrecl)
           ifrcw_full  = iopen('WVR',0,0,mrecl)
           ifrcwi_full = iopen('WVI',0,0,mrecl)
         endif
c      elseif(chipm) then
c        ifrcw  = iopen('WVR.maxloc',0,-1,mrecl)
c        ifrcw_full  = iopen('WVR',0,0,mrecl)
      endif

      else ! .not.lmlw

C --- open WV* files
      if(paralellx0) then
        ifwd   = iopen('WV.d.'//xxt(iqxini,iqxend),1,-1,0)
      else
        ifwd   = iopen('WV.d',1,-1,0)
      endif

c! July 2003. now WV.d contains true size of WVR.
      write (ifwd,"(1x,10i14)") !change format oct2005
     &   nprecx,mrecl,nblochpmx,nwp,niw,nqibz + nq0i-1,nw_i
      allocate(  zw(nblochpmx,nblochpmx) )
      if(debug) print *,' xxx1a2:'
      if    (ixc==101) then
         if(paralellx0) then
           ifrcwi = iopen('WVI.'//xxt(iqxini,iqxend) ,0,-1,mrecl)
         else
           ifrcwi = iopen('WVI',0,-1,mrecl)
         endif
      elseif (normalm) then !ixc==1.or.sergeyv) then
         if(paralellx0) then
           ifrcw  = iopen('WVR.'//xxt(iqxini,iqxend),0,-1,mrecl)
           ifrcwi = iopen('WVI.'//xxt(iqxini,iqxend),0,-1,mrecl)
         else
           ifrcw  = iopen('WVR',0,-1,mrecl)
           ifrcwi = iopen('WVI',0,-1,mrecl)
         endif
c      elseif(chipm) then
c        ifrcw  = iopen('WVR',0,-1,mrecl)
      endif

      endif ! if (lmlw)

      nspinmx = nspin 

c --- READIN ANFcond------------------------------------
      inquire(file='ANFcond',EXIST=anfexist)
      print *,' ANFcond exist? =',anfexist
      if(anfexist) then
        allocate(ibasf(natom),transaf(1:3,natom))
        ifianf = 211
        open(ifianf,file='ANFcond')
        read(ifianf,*)
        read(ifianf,*)
        read(ifianf,*) anfvec(1:3)
        write(6,'(" ANFcond=",3f13.6)') anfvec
        do ibas=1,natom
          read(ifianf,*) ibas1, ibasf(ibas) 
          if(ibas1/=ibas) stop ' readin ANFcond error:ibas1/=ibas'
          write(6,"(' anf ',2i3,f13.6)") ibas1,ibasf(ibas)
c     &     ,transaf(1:3,ibas)
        enddo
        close(ifianf)
        nspinmx =1
      endif ;   if(debug) print *,' xxx1b:'
C... oct2005
      if(smbasis() .and. anfexist) then
       print *, ' smbais=T & anfexist=T is not yet.'//
     &       'At leaset anfx0k shold be changed'
       stop ' smbais=T & anfexist=T is not yet.'
      endif
      if((.not.timereversal()) .and. anfexist) then
       print *,' timereversal=F & anfexist=T is not yet examined.'//
     &       'At leaset you have to consider about anfx0k and so.'
       stop ' timereversal=F & anfexist=T is not yet.'
      endif

C... these are used x0k
      call getkeyvalue("GWinput","nbcutlow",nbcut, default=0 )
      call getkeyvalue("GWinput","nbcutlowto",nbcut2, default=0 )
      write(6,"(' nbcut nbcutlowto=',2i5)") nbcut,nbcut2

C --- ppb ---
!      allocate(ecore(nctot,nspin)) !core energies
      do  is = 1,nspin
! <Phi(R,L,n) Phi(R,L',n') B(R,i)> for each spin 
        if(debug) print *," goto ppbafp is=",is
        irot = 1
        call ppbafp_v2 (irot,ngrpx,is,nspin,
     i              il,in,im, nlnm,  
     d                 nl,nn,nclass,nlnmx,
     i                 mdimx,lx,nx,nxx, 
     i                 cgr, nl-1,       
     i                 ppbrd,           
     o              ppb(1,is)) 
      enddo
      if(debug) print *,' xxx1c:'

C... Set iqxini 
      if(omitqbz) then
        iqxini= nqibz + 1
      elseif(paralellx0) then
     &  !skip
      elseif(bzcase()==1) then   
        iqxini = 2
      else  
        iqxini = 1
      endif

C... check write 1st part for Ec mode to ecorr.chk Nov2004
      if(ecorr_on>0) then
          if(ecorr_on >0.and. (.not.imagomega) ) 
     &    stop ' hx0fp0: ecorr_on .and. (.not.imagomega)'
          write(ieclog,
     &          "('   iq                q                       wk')")
          do iqq = iqxini,iqxend
           if(iqq<=nqibz) then 
            qq = qibze(:,iqq) 
           else
            qq = 0d0 
           endif
           call getwk(iqq, wibz, wqt,nqbz,nqibz,nstibz,nq0i, wk4ec)  
           write(ieclog,"(i5,3x,3f12.8, f15.5)")
     &          iqq, qq, wk4ec 
          enddo
          write(ieclog,*)
          write(ieclog,"('   iw omega(Ry)       wiw')")
          do iww=1,niw
            write(ieclog,"(i5, f10.5, f10.5)") iww,2d0*freq_i(iww),wiw(iww)
          enddo
          write(ieclog,*)
          write(ieclog,"(' Note:IntWgt=wk*wiw.',
     &     '  Ec =\sum_{k,iw) IntWgt(k,iw)*ecqw(k,iw)')")
          close(ieclog)
          open(ieclog,file="ecorr.chk",access='append')

          call getkeyvalue("GWinput","necut_p",necut, default=1 )
          allocate(totexc(necut),trpv(necut),trlog(necut)) 
          totexc = 0d0
          trpv   = 0d0
          trlog  = 0d0
          if(.not.sergeyv) necut=1 
      else
        necut=1
      endif 
      allocate(ecut(necut),ecuts(necut))
      call getkeyvalue("GWinput","ecut_p" ,ecut, necut,default=(/1d10/) )
      call getkeyvalue("GWinput","ecuts_p",ecuts,necut,default=(/1d10/) )
      if(necut/=1) stop 'Sorry! this version is not for necut=1'

! I assume 1 is for majority for eiqr case.
      if(chipm.and.eiqr) then
        schi=1d0 !1d0 means Majority is isp=1. If Majority is isp=2, use schi=-1d0.
        allocate(jcoup(nw_i:nw,neps))
        mmnorm=1d0
      endif

      if(chipm.and.(.not.legas).and.(.not.eiqr)) then
       allocate(jcoup(nw_i:nw,neps))
       allocate(svec(nbloch))
       svec=0d0
c       do ibas=1,natom
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       ibas=1
       write(6,*)"NOTE !!! This version assume ibas=1 "//
     &           "is the only the magnetic atom in cell"
       if(smbasis()) then
         print *,'This version is not for spin sus for smbais=T.'//
     &         'You need to modify MixSpin.' 
         stop  'This version is not for spin sus for smbais=T'
       endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         open(ifv,file='MixSpin.'//charnum3(ibas))
         read(ifv,*) ibasx,lxx
         allocate(nxx_r(0:lxx))
         do i=0,lxx
           read(ifv,*) nxx_r(i)  !   write(6,"(2i5,d13.6)") nxx_r(i)
         enddo
         allocate(spinvec((lxx+1)**2,maxval(nxx_r)))
         spinvec=0d0
         do ilmx = 1, (lxx+1)**2
           lb = ll(ilmx )       !  print *,' lb=',lb,lxx,ilmx
         do ixx = 1, nxx_r(lb)  !  print *,' nn=',nn,nxx_r(lb)
           if(ilmx==1) then
             read(ifv,*) ilm_r, nx_r, spinvec(ilmx,ixx),chg1,chg2
           else
             read(ifv,*) ilm_r, nx_r, spinvec(ilmx,ixx)
           endif   
c           write(6,"(2i5,d13.6)") ilmx, ixx, spinvec(ilmx,ixx)
         enddo
         enddo

! Calculate ChiPM. So sign of omega should be correct.
         spinmom=(chg1-chg2)
         schi=1d0
         if(spinmom<0d0) then
           if(.not.sergeyv)  then 
              freq_r= -freq_r   
           else 
              schi  = -1d0   ! This affects to dpsion2. Obtained results 
                             ! should be the same in both mode.
           endif
         endif

C ... ReOrdering of spinvec in nbas ordering... 
         i=0
         if(ibas>1) i= sum(nblocha(1:ibas-1))
         do lb  = 0, lx (ibas)
         do nb  = 1, nx (lb,ibas)
         do mb  = -lb, lb
          i = i+1 
          ilmx = lb**2+ lb+ mb +1  
          svec(i) = spinvec(ilmx,nb)
          write(6,"(' i lb mb svec svec**2=',3i4,2d13.5)") i,lb,mb,svec(i),svec(i)**2
         enddo
         enddo
         enddo
         deallocate(nxx_r,spinvec)
         close(ifv)
         write(ifchipmlog,"( 'spinmom= ',f14.10)") chg1-chg2
c       enddo

       write(ifchipmlog,"(a)")  
     &' --- Tabel of mcm v.s. ss cutoff.'//
     & '  : mcm means <m|chi^{-1}|m>. '//
     & '  ss is the SVD values of chi^{-1}.'
       if(.not.nolfco) then
         write(ifchipmlog,"(13x,a)")'q              omega'//
     &   '   mcm(ss<ssmin+2Ry)    mcm(ss<ssmin+10Ry)'//
     &   '    mcm(ss<ssmin+5Ry)   mcm(ss<ssmin+100Ry)'//
     &   '    mcm(all)'
       endif
       mmnorm = sum(svec**2)       
       write (6,*)' svecsum=',sum(svec(:)**2)
       write (6,*)' mmnorm=<m|m>=' ,mmnorm
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      endif

      if(epsmode.and.sergeyv.and.nolfco) then
        allocate( rcxqmean(nwhis,npm))
      endif

C --- Calculate x0(q,iw) and W ------------------------------------------
c============ loop over q point ========================================
c=======================================================================
c      if(iepsmode==202) iqxini = nqibz+1
      do 1001 iq = iqxini,iqxend ! q=(0,0,0) is omitted!
        call cputid (0)
        q = qibze(:,iq)
        if(iq<=nqibz.or.ngczero) then ! iqx = iq
          qq=q 
        else  !    iqx = 1            ! corresponding q=0 
          qq=0d0 
        endif

        call readqg('QGcou', qq, ginv,  quu,ngc,ngveccB)



c ---readin coulomb matrix
        if(chipm.and.(.not.legas).and.(.not.eiqr)) then
         allocate(vcoul(1,1))!dummy
         ngb = nbloch + ngc !ngb is readin from vcoul 25jan2006
        else
         read(ifvcfpout) ngb       ! oct2005 !I added ngb to VCCFP
         allocate(vcoul(ngb,ngb))
         read(ifvcfpout) vcoul,qx  ! coulomb matrix
         if(sum(abs(q-qx))>1d-10) then
         write(6,"('iq q=',i5,3f9.4,' qx(VCCFP)=',3f9.4)")iq,q,qx
         stop 'hx0fp0: VCCFP is not compatible'
         endif
        endif
        write(6,"('do 1001 start----- iq q=',i5,3f9.4,' qq=',3f9.4)")iq,q,qq

c --- eps mode
        if(epsmode) then !iepsmode/=0) then      ! only calculate iq>nqibz
          if(iq>nqibz) then
            if(wqt(iq-nqibz)==0d0) then
              iqixc2 = iq- (nqibz+nq0ix)
              if(chipm) then
                ifchipmn=iopen('ChiPM'//charnum4(iqixc2)//'.nolfc.dat',1,3,0)
                if(eiqr) then
                  write(ifchipmn, "(10x,'q(1:3)',12x,'  w(Ry)',18x
     &            ,    '<e^iqr|chi^pm|e^iqr>^-1',23x
     &            ,8x ,'<e^iqr|chi^pm|e^iqr> ')")
                else
                  write(ifchipmn, "(10x,'q(1:3)',12x,'  w(Ry)',18x
     &            ,    '<m|m>**2/<m|chi|m> ',23x
     &            ,8x ,'<m|chi|m>/<m|m>**2 ')")
                endif
                if(.not.nolfco) then
                  ifchipm = iopen ('ChiPM'//charnum4(iqixc2)//'.dat',1,3,0)
                  write(ifchipm,"(10x,'q(1:3)',12x,'  w(Ry)',18x
     &            ,' <m|chi^{-1}|m> (ssmin only) '   ,15x
     &            ,' <m|chi^{-1}|m> (ss< ssmin+2Ry) ',8x
     &            ,' ratio=mcm(ssmin+2Ry)/nolfc')")
                endif
              else
               filepsnolfc ='EPS'//charnum4(iqixc2)//'.nolfc.dat'
               ifepsdatnolfc = iopen ( filepsnolfc,1,3,0)
               write(ifepsdatnolfc,"(a)")
     &          ' q(1:3)   w(Ry)   eps    epsi  --- NO LFC'
               if(.not.nolfco) then !ixc==3.or.iepsmode==203) then
                fileps = 'EPS'//charnum4(iqixc2)//'.dat'
                ifepsdat = iopen ( fileps,1,3,0)
                write(ifepsdat,"(a)") 
     &          ' q(1:3)   w(Ry)   eps  epsi --- LFC included. '
               endif
              endif
              goto 110
            endif
          endif
          deallocate(vcoul)
          cycle  !Cycle here
  110     continue
          if((.not.chipm).or.legas.or.eiqr) then
            read(ifgb0vec,*) qgbin(1:3),ngbin,igc0,dnorm
            write(6,"(' --- Readin Mix0vec: ',3d13.6,2i5,d18.8)")
     &      qgbin(1:3),ngbin,igc0,dnorm
            if(ngb/=ngbin) stop 'hx0fp0: ngb/=ngbgin'
            if(sum(abs(qgbin(1:3)-q)) >1d-8) stop 'hx0fp0: qgbin/=q'
            write(6,"(' ngb nwp niw=',3i8)")ngb,nwp,niw
            allocate(gbvec(ngb),zzr(ngb),x0mean(nw_i:nw))
            x0mean=0d0
            do i=1,ngb
             read(ifgb0vec,"(4d24.15)") gbvec(i), zzr(i)
            enddo
          elseif(chipm.and.nolfco) then
            allocate(zzr(ngb),x0mean(nw_i:nw))
            x0mean=0d0
            zzr   =0d0
            zzr(1:nbloch)=svec(1:nbloch)
          endif
        endif
   
        allocate( zw0(ngb,ngb),
     &    zxq (ngb,ngb,nw_i:nw), 
     &    zxqi(ngb,ngb,niw))
        write(6,*) "---  end of allocation zw0 zxq zxqi"
        if(ixc==6 ) then
           allocate( rcxq(ngb,ngb,nwp,npm) )
        elseif(epsmode.and.sergeyv.and.nolfco) then !iepsmode==202) then
          allocate( rcxq(1,1,1,npm) )
          rcxqmean=0d0
        elseif(sergeyv) then
          allocate( rcxq(ngb,ngb,nwhis,npm) )
        endif
        zxq=0d0; zxqi=0d0
        if(sergeyv) rcxq=(0d0,0d0)
c============ loop over spin============================================
c=======================================================================
        do 1003 is = 1,nspinmx  
          write(6,"(' *** ',2i4,' out of nqibz+n0qi nsp=',2i4,' *** ')") 
     &    iq, is, nqibz + nq0i,nspin
          if(debug) print *,' niw nw=',niw,nw

C.. spin chi_charge or chi_+-
          isf=is
          if(chipm) then
            write(6,*)" chi_+- mode ixc=",ixc
            if(is==1) isf=2
            if(is==2) isf=1
            if(sergeyv) then
               rcxq=0d0
               if(epsmode.and.nolfco) rcxqmean=0d0
            endif
          endif


c===========tetraini block tetra==.true.===============================1ini
c --- tetrahedron method --- preparatory part
c  nbnb    : total number of pairs (ib,jb) with non-zero weight.
c  ib=n1b  : band index for occ.   1\ge n1b \ge nband+nctot. 
c         Valence-core order(Core index follows valence index).
c  jb=n2b  : band index for unocc. 1\ge n2b \ge nband
c  wwk(ibib,...)  = (complex)weight for the pair for n1b(ibib...),n2b(ibib...).
          if(tetra) then
           write(6,"(' tetra mode niw nw=',2i5)") niw,nw  
C ... kqxx(kx) ekxx --- 
c    ekxx1  for  rk
c    ekxx2  for q+rk.  See tetwt4
            do kx = 1, nqbz
              call readeval(    qbz(:,kx)-qqq, is, ekxx1(1:nband, kx) )
              call readeval(q + qbz(:,kx)-qqq,isf, ekxx2(1:nband, kx) )
            enddo
C ... Get pairs (n1b n2b) with non-zero tetrahedron wieghts.
c     the pairs are not dependent on the energy omega 
c     in the denominator of the dielectric function.
            write(6,"(' -- First tetwt5 is to get size of array --')")
            job = 0
            if(npm==1) then
              ncc=0
            else
              ncc=nctot
            endif
            allocate( demin(nband+nctot,nband+ncc,nqbz,npm),
     &                demax(nband+nctot,nband+ncc,nqbz,npm) )        
            allocate( iwgt (nband+nctot,nband+ncc,nqbz,npm) )
            if(nctot==0) then
               deallocate(ecore)
               allocate(ecore(1,1)) !this is dummry
            endif   
            call tetwt5x_dtet4(  npm,ncc,
     i           q, ekxx1, ekxx2, qbas,ginv,ef, 
     d          ntetf,nqbzw,  nband,nqbz,
     i          nctot,ecore(1,is),idtetf,qbzw,ib1bz, 
     i          job,
     o          iwgt,nbnb,   !job=0
     o          demin,demax, !job=0
     i          frhis, nwhis,  ! job=1    not-used
     i          nbnbx,ibjb,nhwtot,  ! job=1    not-used !bug:kino-fix  nhwtot was missing though it was dummy
     i          ihw,nhw,jhw,   ! job=1    not-used
     o          whw,           ! job=1    not-used
     i        iq,is,isf,nqibz) 
            if(ixc/=4.and.ixc/=6.and.(.not.sergeyv)) deallocate(demin,demax)
            nbnbx = maxval(nbnb(1:nqbz,1:npm)) !nbnbx = nbnbxx
            if(debug) print *,' nbnbx=',nbnbx
            allocate(  n1b(nbnbx,nqbz,npm)
     &                ,n2b(nbnbx,nqbz,npm))
            do jpm=1,npm
              call rsvwwk00_4(jpm, iwgt(1,1,1,jpm),nqbz,nband,nctot,ncc, nbnbx, !2006
     o        n1b(1,1,jpm), n2b(1,1,jpm), noccxvx(jpm), nbnbtt(1,jpm))
            enddo
            if(sum(abs(nbnb-nbnbtt))/=0)stop'hx0fp0:sum(nbnb-nbnbtt)/=0'
            noccxv = maxval(noccxvx)
            noccx  = nctot + noccxv
            print *,' Tetra mode: nctot noccxv= ',nctot,noccxv
            deallocate(iwgt)
          endif
c=========end of tetraini block==========================================1end


c=======TetrahedronWeight_4 block. tettwt4  ixc<=5 ==============2ini
          if(tetra.and.(.not.sergeyv)) then !(ixc==101.or.ixc<=5.or.ixc==22.or.ixc==23)) then
            allocate( wgt(nband+nctot,nband,nqbz) ) 
            if(.not.tmpwwk) allocate( wwk(nbnbx,nqbz,niw+nw+1)) ! +1 is for nwp case 
            if(tmpwwk) rewind(ifwwk)
            nwmax = niw +nw +1
            if(ixc==4.or.ixc==5) nwmax = nw+1
            do iw = 1, nwmax
              if(ixc==4.or.ixc==5.or.ixc==6) then 
                fff = freq_r(iw-1)+img*delta 
                write(6,"(' *** iw on real= ',i3,' omega(Ry)=',2d15.6)")iw-1,fff*2d0
              elseif(iw<=niw) then
                fff = img * freq_i(iw)            ! along img axis
                write(6,"(' *** iw on imag= ',i3,' omega(Ry)=',2d15.6)")iw,fff*2d0
              endif
            enddo

            do 1010 iw = 1,nwmax
              if(ixc==4.or.ixc==5.or.ixc==6) then 
                fff = freq_r(iw-1)+img*delta 
c                write(6,"(' *** iw on real= ',i3,' omega(Ry)=',2d15.6)")
c     &          iw-1,fff*2d0
              else  
                if(iw<=niw) then
                  fff = img * freq_i(iw)            ! along img axis
c                 write(6,"(' *** iw on imag= ',i3,' omega(Ry)=',2d15.6)")
c     &           iw,fff*2d0
                else
cnw iw-niw ---> iw-niw-1
                  fff = freq_r(iw-niw-1) + img*delta ! delta is in a.u.  ! along real axis
c                 write(6,"(' *** iw on real= ',i3,' omega(Ry)=',2d15.6)") iw,fff*2d0
                endif  
              endif

              call tetwt4x_dev(fff,  
     i              q, ekxx1, ekxx2, qbas,ginv,ef,
     d          ntetf, nqbzw, nband,nqbz,
     i          nctot,ecore(1,is),idtetf,qbzw,ib1bz, chipm, !ixc,
     o          wgt, nbnbtt,
     i       iq,is,isf,nqibz) ! new input for devided-tetrahedron method.

!-------------------
              if(chipm.and.spinmom<0) wgt= dconjg(wgt)


              if(dreal(fff)==0d0) wgt = dreal(wgt) ! imag part zero clear for omg=0 !this is not necessary.

              if(tmpwwk) then
                print *,' keyx rewind ifwwk nwmax=',iw,nbnbx,nqbz
                call rsvwwk4_w( wgt,nqbz,nband,nctot,nbnbx,
     i          n1b, n2b, nbnb, ifwwk)
              else
                call rsvwwk4( wgt,nqbz,nband,nctot,nbnbx,
     i          n1b, n2b, nbnb,
     o          wwk(1,1,iw) )
              endif
ccccccccccccccccccccccccccccc
c             write(996,"(' --- iq iw=',2i3)" ),iq,iw
c             do ik  = 1,nqbz
c             do ibib= 1,nbnb(ik)
c               write(996,"(4i4,2d16.8)" )
c     &         ik,iw,n1b(ibib,ik),n2b(ibib,ik),wwk(ibib,ik,iw)
c             enddo
c             enddo
cccccccccccccccccccccccccccc
c--------------------------
c              if( maxval(nbnb(1:nqbz,iw)) > nbnbx) then ! security check
c               print *,' hx0fp0: maxval(nbnb(1:nqbz,iw)) > nbnbx ',
c     &          nbnbx,nbnb(1:nqbz,iw)           
c                stop ' hx0fp0: maxval(nbnb(1:nqbz,iw)) > nbnbx '
c              endif
c--------------------------
ccccccccccccccccccccccccccccccc
c            if(iq==2.and.iw==nwmax) write(996,*)' --- end of iq=2 ---'
c            if(iq==3.and.iw==1) then
c             write(996,"(' --- iq iw=',2i3)" ),iq,iw
c             do ik  = 1,nqbz
c             do ibib= 1,nbnb(ik)
c               write(996,"(' k ib jb=',3i4,' wwk=',2d16.8)" )
c     &         ik,n1b(ibib,ik),n2b(ibib,ik),wwk(ibib,ik,iw)
c             enddo
c             enddo
c             stop '----- test end ------------'
c           endif
ccccccccccccccccccccccccccccccc
 1010       continue
            print *,'end of 1010 loop'
            deallocate(wgt)
            if(tmpwwk) then
              allocate( wwk(nbnbx,nqbz,nwmax) )
              print *,' wwk allocated'
              ifwwk=iclose('tmp.wwk')
c              rewind(ifwwk)
              ifwwk = iopen('tmp.wwk',0,-1,0)
              print *,' keyx rewind ifwwk nwmax=',nwmax,nbnbx,nqbz,niw,nw
              do iw=1,nwmax
              read(ifwwk) wwk(1:nbnbx,1:nqbz,iw)
              enddo
              ifwwk=iclose('tmp.wwk')
              print *,' readend of wwk'
            endif
          endif
c=======TetrahedronWeight_4 block end.========================2end


c==============x0kf_v2 block original(.not.sergeyv) case =====================3ini
          if(tetra.and.((.not.sergeyv).and.ixc/=4) )then !ixc<=3.or.ixc==5.or.ixc==101.or.ixc==22.or.ixc==23)) then
           if(ixc==5) then
           wwk = dimag(wwk)
           endif
           print *,' --- goto x0kf_v2z'
           call x0kf_v2z(   ! lessl,
     i            wwk,n1b,n2b,nbnbx,nbnb, tetra, ! for tetra=T
     i              q, is,isf, nspin, !is,  !,ifcphi ifev1,ifev2 ifrb(is),ifcb(is),ifrhb(is),ifchb(is),
     i              qbas,ginv,ef,qbz,wbz,   !iindxk,
c     i            index_qbz,n_index_qbz, jan2004
     i              nblocha ,nlnm,nlnmv,nlnmc,iclass, !mdim=nblocha
     i              ppb(1,is), !,pdb(1,is),dpb(1,is),ddb(1,is),
     i              freq_r,freq_i,delta,realomega,imagomega, 
     i              icore,ncore,
     &              ecore(1,is),
     d            nlmto,nqbz,noccx,noccxv,
     d            nl,nclass,natom,nnc,
     d            nlnmx,mdimx,nbloch,nwp,  niw,
c     i    iq,ngpn,ngcni(iqx),ngpmx,ngcmx  , geigB(1,1,1,is), Jan2004
c     i    ngvecpB,ngveccB, ngveccB(1,1,iqx),
     i     iq,ngc,ngpmx,ngcmx  ,
c     i    ngveccB,  !ngveccB(1,1,iqx),
     i     nqbze,nband,nqibz, 
     o     zxq,zxqi,            ! zxq... are accumulated variable for spins 
     i     nolfco,zzr,    !for iepsmode
     o     x0mean, !   )  !for iepsmode
     i     vcoul,ngb ) ! for te
          endif
          if(epsmode.and.nolfco) 
     &      print *,' nolfco --- end of x0kf_v2 sum x0mean',sum(abs(x0mean))
c=============x0kf_v2 block end===============================3end


c=======TetrahedronWeight_5 block. tetwt5  ixc==,4,6,11 =======4ini
        if(ixc==4.or.sergeyv) then
C --- METHOD (tetwt5) for the tetrahedron weight
!   Histogram secstions are specified by frhis(1:nwp)
!     The 1st   bin  is     [frhis(1),  frhis(2)]   ...
!     The last  bin  is     [frhis(nw), frhis(nwp)].
!     nwp=nw+1; frhis(1)=0
! takao-feb/2002 
            if(frhis(1)/=0d0) stop ' hx0fp0: we assume frhis(1)=0d0'
            write(6,*)' --- Frequency  histogram  sections (a.u.)---- '
            if(onceww(1)) then
              do ihis= 1, nwhis
              write(6,"(' ihis Init  End=', i4,2f13.6)")
     &        ihis,frhis(ihis),frhis(ihis+1)  
              enddo
            endif  
            write(6,*)' ----------------nbnbx nqbz= ',nbnbx,nqbz
C ... make index sets
            allocate(ihw(nbnbx,nqbz,npm),nhw(nbnbx,nqbz,npm),jhw(nbnbx,nqbz,npm))
            jhwtot = 1
            do jpm =1,npm
            do ik   = 1,nqbz
            do ibib = 1,nbnb(ik,jpm)
c             print *,' ik ibib = ',ik,ibib
cc              ib1=  n1b(ibib,ik,jpm) 
cc              ib2 = n2b(ibib,ik,jpm) 
c              print *,' goto histrange=',ib1,ib2
c     &         ,demin(n1b(ibib,ik),n2b(ibib,ik),ik)
c     &         ,demax(n1b(ibib,ik),n2b(ibib,ik),ik)
              call hisrange( frhis, nwhis,  
     i         demin(n1b(ibib,ik,jpm),n2b(ibib,ik,jpm),ik,jpm),
     i         demax(n1b(ibib,ik,jpm),n2b(ibib,ik,jpm),ik,jpm),
     o         ihw(ibib,ik,jpm),nhw(ibib,ik,jpm))
              jhw(ibib,ik,jpm)= jhwtot
              jhwtot = jhwtot + nhw(ibib,ik,jpm)
            enddo
            enddo
            enddo

            nhwtot = jhwtot-1
            print *,' nhwtot=',nhwtot
            deallocate(demin,demax)
            allocate( whw(nhwtot),   ! histo-weight
     &                ibjb(nctot+nband,nband+ncc,nqbz,npm) )

            whw=0d0
            ibjb = 0
            do jpm=1,npm
            do ik   = 1,nqbz
            do ibib = 1,nbnb(ik,jpm)
              ib1  = n1b(ibib,ik,jpm) 
              ib2  = n2b(ibib,ik,jpm) 
              ibjb(ib1,ib2,ik,jpm) = ibib
            enddo
            enddo
            enddo
C ... Generate the histogram weights whw 
            job=1
            write(6,*) 'goto tetwt5x_dtet4 job=',job
            call tetwt5x_dtet4(  npm,ncc,
     i           q, ekxx1, ekxx2, qbas,ginv,ef, 
     d          ntetf,nqbzw, nband,nqbz,
     i          nctot,ecore(1,is),idtetf,qbzw,ib1bz, 
     i          job,
     o          iwgt,nbnb,      ! job=0
     o          demin,demax,   ! job=0
     i          frhis,nwhis,       ! job=1
     i          nbnbx,ibjb,nhwtot, ! job=1
     i          ihw,nhw,jhw,       ! job=1
     o          whw,               ! job=1
     i   iq,is,isf,nqibz) 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$            print *,' nqbz=',nqbz 
c$$$            do ik   = 1,nqbz
c$$$            print *,' nbnb=',ik,nbnb(ik) 
c$$$            do jpm  = 1,npm
c$$$            do ibib = 1,nbnb(ik,jpm)
c$$$             ib1  = n1b(ibib,ik,jpm) 
c$$$             ib2  = n2b(ibib,ik,jpm) 
c$$$             write(6,"('-sumchk whw-- jpm ik ib jb sum=',4i4,d13.6)") 
c$$$     &       jpm,ik,ib1,ib2,
c$$$     &       sum( whw(jhw(ibib,ik,jpm):jhw(ibib,ik,jpm)+nhw(ibib,ik,jpm)))
c$$$c              ihww = ihw(ibib,ik)
c$$$c              do iw = 1, nhw(ibib,ik)
c$$$c               write(3666,"('jpm ihis [a b] whw =',i3,2x,i4,2f8.5,f18.12,d14.6)")
c$$$c     &         jpm, ihww+iw-1, frhis(ihww+iw-1),frhis(ihww+iw),
c$$$c     &         whw(jhw(ibib,ik)+iw-1), whw(jhw(ibib,ik)+iw-1)
c$$$c              enddo     
c$$$            enddo
c$$$            enddo
c$$$            enddo
c$$$c            stop 'hx0fp0.f ---  testend xxx2'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c            write(1236,"('-----',5i8)") noccxvx
c            do ik  = 1,nqbz
c              write(1236,"('--- ',5i8)") ik, nbnb(ik),nbnbtt(ik)
c            do ibib= 1,nbnb(ik)
c              write(1236,"(2i6,2f12.6)") n1b(ibib,ik),n2b(ibib,ik)
c     &         ,demin(n1b(ibib,ik),n2b(ibib,ik),ik)
c     &         ,demax(n1b(ibib,ik),n2b(ibib,ik),ik)
c            enddo
c            enddo
c            stop 'hx0fp0.m.fal.f ---  testend'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        endif     
c=======TetrahedronWeight_5 block end ==================4end




c==============x0kf_v3h block. ixc==6,11 =========================5ini
c Instead of wwk, you have to use whw(nwhtot). 
c You will need to pass these arrays to x0kf
c   ihw(ibjb,kx): omega index, to specify the section of the histogram.
c   nhw(ibjb,kx): the number of histogram sections
c   jhw(ibjb,kx): pointer to whw 
c   whw( jhw(ibjb,kx) ) \to whw( jhw(ibjb,kx) + nhw(ibjb),kx)-1 ), where ibjb=ibjb(ib,jb,kx)
c     : histogram weights for given ib,jb,kx for histogram sections
c       from ihw(ibjb,kx) to ihw(ibjb,kx)+nhw(ibjb,kx)-1.
c------------------------------------------------------------------
        if(sergeyv) then ;  print *, ' --- goto x0kf_v3hz '
c- new version only along real axis for 
          call x0kf_v3hz(npm,ncc,   
     i          ihw,nhw,jhw,whw,nhwtot,  ! tetwt5
     i            n1b,n2b,nbnbx,nbnb,    ! use whw by tetwt5 ,
     i              q,  
     i              nspin,is,isf,
     i              qbas,ginv,ef,  qbz,wbz, 
     i              nblocha ,nlnm,nlnmv,nlnmc,iclass, 
     i              ppb(1,is), 
     i              icore,ncore,
     d            nlmto,nqbz,noccx,noccxv,
     d            nl,nclass,natom,nnc,
     d            nlnmx,mdimx,nbloch,  nwhis,
     i    iq,ngb,ngc,ngpmx,ngcmx,  !ngb/=ngc+nbloch for smbasis()=T oct2005
     i    nqbze,nband,nqibz, 
ctm, 080311
     i    qbze,iwbnd,nwf,lwbnd,
ctm, 080820
     i    ib_ent,lent,
     o    rcxq, ! See rcxq_zxq in x0kf_v3h.f; rcxq is equivalent with zxq; 
                !  rcxq is the accumulating variable for spins 
     i    nolfco,zzr,  !for nolfco
     o    rcxqmean   ) !for nolfco
        endif

cccccccccccccccccccc
c          do iw  = 1,nwhis 
c          do i=1,ngb
c          do j=1,ngb
c             if(i==j) rcxq(j,i,iw)= dreal(rcxq(j,i,iw))
c          enddo; enddo; enddo
cccccccccccccccccccccc

c=========== Check write mode4 block ixc==4 ============6ini
         if(ixc==4) then
C --- check write ---
           ichkhis = iopen ('tethis.chk',1,-1,0)
           do jpm   = 1,npm
           do ik   = 1,nqbz
           do ibib = 1,nbnb(ik,jpm)
              ib1  = n1b(ibib,ik,jpm) 
              ib2  = n2b(ibib,ik,jpm) 
              ihww = ihw(ibib,ik,jpm)
              write(ichkhis,"('--- ik ib jb jpm =',4i4)") ik,ib1,ib2,jpm
              do iw = 1, nhw(ibib,ik,jpm)
                write(ichkhis, "('his [min max]=',i6,2f7.4,
     &             ' whw wwk*h= ',2d12.5,'( =',d12.5,'+',d12.5,')')")
     &          ihww+iw-1, freq_r(ihww+iw-1), freq_r(ihww+iw),
!-------------------------------------------------------
     &          whw(jhw(ibib,ik,jpm)+iw-1),                     !weight for the histgram range. by tetwt5
!------------ v.s. -------------
     &          dimag(wwk(ibib,ik,ihww+iw-1)+wwk(ibib,ik,ihww+iw))/2d0
     &           * (freq_r(ihww+iw-1)- freq_r(ihww+iw)),    !weight=weight1+weight2 for each bin of the histgram range by tetwt4.fal
!-------------------------------------------------------
     &          dimag(wwk(ibib,ik,ihww+iw-1))/2d0
     &           * (freq_r(ihww+iw-1)- freq_r(ihww+iw)),    !weight1 for the histgram range by tetwt4
     &          dimag(wwk(ibib,ik,ihww+iw))/2d0
     &           * (freq_r(ihww+iw-1)- freq_r(ihww+iw))     !weight2 for the histgram range by tetwt4
              enddo     
           enddo
           enddo
           enddo
           ichkhis = iclose('tethis.chk')
           print *,' OK! hx0fp1_tetwt5_test end. Histgram method.'
           print *,' The final reulst wwk vs.whw is in tethis.chk!'
           print *,' This code is a routine to test tetwt5.'
           stop    ' ixc==4  Weight Histgram check. wwk vs.whw'
         endif  ! end of if(tetra)
c=========== Check write mode4 block end  ==============6end


         if(ixc==4.or.sergeyv)         deallocate(ihw,nhw,jhw, whw,ibjb)
         if(tetra.and.(.not.sergeyv) ) deallocate(wwk) ! !ixc<=5.or.ixc==101.or.ixc==22.or.ixc==23)) 
         if(tetra)                     deallocate( n1b,n2b)

         iecut=1

         if(sergeyv) then !This was in 1003 loop as if(ixc==11.and.is==nspinmx) then
          if(debug) write(6,"(a)") ' --- goto dpsion4 --- '
          if(is==nspinmx.or.chipm) then
            if(ixc==6) then ! ... convert from rcxq to zxq for ixc==6
              stop "--- !This path is not tested yet !!!!---"
              zxq(:,:,0)=0d0
              do jpm =1,npm 
              do iw  = 1,nw 
              do i=1,ngb
              do j=1,i  
                 ixx=1
                 if(jpm==2) ixx=-1
              zxq(j,i,iw*ixx )=rcxq(j,i,iw,jpm)
              zxq(i,j,iw*ixx )=dconjg(zxq(j,i,iw*ixx))
              enddo
              enddo
              enddo
              enddo
            else
cccccccccccccddddddccccccccccccccccccccccc
c          rcxq(:,:,:,1) =rcxq(:,:,:,2)
cccccccccccccccccccddddddccccccccccccccccc
             write(6,*)"goto dpsion4 nwhis nw_i nw_w=",nwhis,nw_i,nw
             call dpsion4(frhis,nwhis, freq_r, nw, freq_i,niw, realomega, imagomega, 
     i        ngb, rcxq, npm,nw_i,  ! rcxq is alterd---used as work
     o        zxq, zxqi, 
     i        nolfco,chipm, schi,is,  rcxqmean,ecut(iecut),ecuts(iecut), !iepsmode,
     o        x0mean)  !zxq x0mean is accumlating veriable for chipm mode
            endif
          endif  ;      print *,' --- end of dpsion4 ----'
          if(epsmode.and.nolfco) print *,' end of dpsion4 sum x0mean',sum(abs(x0mean))
         endif
 1003   continue  !end of spin loop =====
        if(allocated(rcxq) ) deallocate(rcxq)

 1033    continue
         print *," end of 1033 xxxxxxxxx"
c=======================================================================


c==============x0kf_v3h block end  =========================
c$$$- 13July2005 
c$$$    Takao needed to remove do 2005 loop so as to put dpsion2 into do 1003
c$$$    loop so as to make chipm mode works ok(do dpsion2 for eachs spin).
c$$$    If necessary we may need to recover do 2005 mode---this was for e_rpa test...
c$$$C... this loop is to calculate Ec for different ecut, ecuts.
c$$$        do 2005 iecut=1,necut !In usual mode, necut=1
c$$$c-----------------
c$$$        if(sergeyv) then !This was in 1003 loop as if(ixc==11.and.is==nspinmx) then
c$$$          if(debug) print *, ' --- goto dpsion2'
c$$$          call dpsion2(frhis,nwhis, freq_r, nw, freq_i,niw, realomega, imagomega, 
c$$$     i      ngb, rcxq,     
c$$$     o      zxq, zxqi, ! zxq... are Not the accumulating variable for spins 
c$$$     i      nolfco,chipm,is,  rcxqmean,ecut(iecut),ecuts(iecut), !iepsmode,
c$$$     o      x0mean)
c$$$          print *,' --- end of dpsion2'
c$$$          if(epsmode) print *,' end of dpsion sum x0mean',sum(abs(x0mean))
c$$$c          do iw=1,nwp
c$$$c            write(6,"(a,i5,2d13.5)")' sumcheck zxq=',iw,sum(abs(zxq(:,:,iw)))
c$$$c          enddo 
c$$$        endif
   
c       write(ifdpin)  nwhis,nw,niw, ngb, frhis,freq_r, freq_i,rcxq


C ---  RealOmega ====================================
        if (realomega) then
         if(chipm) then !ixc==22.or.ixc==23) then
           if (nspin==1) stop 'chipm modes are for nspin==2'
           if (anfexist) then 
! Not tested yet
!          if(anfexist) x0mean= 2d0*x0mean  !if paramagnetic, multiply x0 by 2
            write(6,*) 'I have not tested anfexist=t and ixc==22 23 '//
     &      '--- Maybe x0mean = 2d0* x0mean*2 will OK.'
            stop 'hx0fp0: tail output! Not yet for anfexist&ixc==22 23'
           endif   
         elseif(epsmode.and.nolfco.and.(.not.chipm)) then !ixc==2.or.iepsmode==202) then
           if (nspin==1) x0mean= 2d0*x0mean  !if paramagnetic, multiply x0 by 2
           if (anfexist) then 
! Not tested yet
!          if(anfexist) x0mean= 2d0*x0mean  !if paramagnetic, multiply x0 by 2
           write(6,*) 'I have not tested anfexist=t and ixc==2 '//
     &     '--- Maybe x0mean = 2d0* x0mean*2 will OK.'
           stop 'hx0fp0: tail output! Not yet for anfexist&ixc==2'
           endif   
         else 
          if (nspin == 1) zxq = 2d0*zxq   !if paramagnetic, multiply x0 by 2
          if (anfexist) then   !Antiferro accelaration.
c           nwmax = nw   
c           if(ixc==4.or.ixc==5.or.sergeyv) nwmax=nw
           do iw  = nw_i,nw !max 
           call anfx0k(natom,nclass,nblocha,iclass,pos
     i     ,nbloch,ngc, !ngcni(iqx),
c     i    q, ngveccB(1,1,iqx),qbas, ! for q+G
     i     q, ngveccB,qbas, ! for q+G
     i     anfvec,ibasf,             ! these are antiferro informations.
     i     zxq(1,1,iw))              ! i/o 
           enddo
          endif
         endif

c         write (ifxd,"(1x,3f10.4)") q(1),q(2),q(3)
c         write (ifrx) rxq,cxq
         if(epsmode) then
           if(nolfco) then
             ttt='without LFC'
           else
             ttt='with LFC'
           endif
           if(chipm) then
            write(6,*) '--- chi_{+-}}^{-1}      --- '//ttt
           else
            write(6,*) '--- dielectric constant --- '//ttt
           endif
         elseif(ixc==1.or.sergeyv) then
            write(6, *)" trace check for W-V"
         endif

C... Get O^{1/2}= sqrt(ovlp)
         if(chipm.and.(.not.nolfco)) then
           allocate(sqovlp(ngb,ngb),sqovlpi(ngb,ngb))
         write(6,*)"goto getsqovlp"
           call getsqovlp(q,ngc,ngb, sqovlp) 
           sqovlpi=sqovlp
         write(6,*)"goto matcinv"
           call matcinv(ngb,sqovlpi)
           allocate(zzz(ngb,ngb),UU(ngb,ngb),VT(ngb,ngb),SS(ngb))
           allocate(UU0(ngb,ngb),VT0(ngb,ngb),SS0(ngb))
         endif
         write(6,*)"zzzz nw_i nw=",nw_i,nw
         do iw  = nw_i,nw  !Feb2006. Before it was 1:nwp (nwp=nw+1).
                           !  So freq_r(iw-1) is shifted to freq_r(iw).
cnw all freq is from zero freq(0:nwp-1); nwp is the number of real flequencies.
            frr= dsign(freq_r(abs(iw)),dble(iw))
            if(.not.epsmode) then  !if(ixc==1.or.sergeyv) then
              imode = 1
              if(debug) write(6,"(' rrr sum zxq=',i5,4d23.14)") iw,sum(zxq(:,:,iw)),sum(abs(zxq(:,:,iw)))
              if (lmlw) then
ctm, 060517
                 allocate(zw_full(nblochpmx,nblochpmx)
     &                   ,zw_full2(ngb,ngb))
c                 read(ifrcw_full, rec=((iq-iqxini)*nwp +iw) ) zw_full
                 nrec = (iq-iqxini)*(nw-nw_i+1)+ iw-nw_i+1
                 read(ifrcw_full, rec=nrec ) zw_full   !  WP = vsc-v
                 zw_full2 = zw_full(1:ngb,1:ngb)
                 deallocate(zw_full)
                 call wcf_maxloc( ngb, vcoul, zxq(1,1,iw), zw_full2,
     &                            imode,
     &                            zw0)
                 deallocate(zw_full2)
              else
                 call wcf( ngb, vcoul, zxq(1,1,iw), imode,
     &                     zw0)
              endif
              if(debug) write(6,"(' rrr sum zw0=',i5,4d23.14)") iw,sum(zw0),sum(abs(zw0))
              zw(1:ngb,1:ngb) = zw0 ! zw(nblochpmx,nblochpmx)        
              write(ifrcw, rec=((iq-iqxini)*(nw-nw_i+1)+ iw-nw_i+1 ) ) zw   !  WP = vsc-v
             !----------------------------
              call tr_chkwrite("freq_r iq iw realomg trwv=", 
     &               zw, iw, frr,nblochpmx, nbloch,ngb,iq)
             !----------------------------
            elseif(epsmode.and.(.not.chipm)) then !.and.ixc/=23) then ! No LFC (local field correction). It's better to use echo 4| hbasfp0. 
              vcmean= sum( dconjg(gbvec) * matmul(vcoul,gbvec) ) 
              if(.not.nolfco) then
                x0mean(iw) = sum( dconjg(zzr)* matmul(zxq(:,:,iw),zzr)) 
              endif
              epsi(iw,iqixc2) = 1d0/(1- vcmean * x0mean(iw)) 
              write(6,'(" iq iw omega eps epsi noLFC=",2i6,f8.3,2d23.15,3x, 2d23.15,
     &        " vcmean x0mean =", 2d23.15,3x, 2d23.15)') iqixc2,iw,2*frr,
     &        1d0/epsi(iw,iqixc2),epsi(iw,iqixc2),vcmean,x0mean(iw)
              write(ifepsdatnolfc,'(3f12.8,2x,d12.4,2d23.15,2x,2d23.15)') 
     &        q, 2*frr, 1d0/epsi(iw,iqixc2),epsi(iw,iqixc2)
              if(.not.nolfco) then ! With LFC
                imode=2
                call wcf( ngb, vcoul, zxq(1,1,iw), imode,
     &                zw0) !  write(6,"('ssschk1=',3d13.5)") sum(abs(zw0)) sum(abs(gbvec))
                epsi(iw,iqixc2)= sum( dconjg(gbvec) * matmul(zw0,zzr) ) 
                write(6,'( " iq iw omega eps epsi  wLFC="
     &          ,2i6,f8.3,2d23.15,3x, 2d23.15)') 
     &          iqixc2,iw,2*frr,1d0/epsi(iw,iqixc2),epsi(iw,iqixc2)
                write(6,*)
                write(ifepsdat,'(3f12.8,2x,d12.4,2d23.15,2x,2d23.15)') 
     &          q, 2*frr,1d0/epsi(iw,iqixc2),epsi(iw,iqixc2)
              endif

            elseif(epsmode.and.chipm) then
              if(nolfco) then !ixc==22) then 
C... ChiPM mode without LFC
c              mmnorm = sum(svec**2);write(6,"('mmnorm=svec*svec=', d14.6)") mmnorm
                if(legas) then !test case
                 vcmean= sum( dconjg(gbvec) * matmul(vcoul,gbvec) ) 
                 write(ifchipmn,'(3f12.8,2x,f8.5,2x,2d23.15)') 
     &            q, 2*schi*frr, 1d0-vcmean*2*x0mean(iw)  !4*pi*alat**2/sum(q**2)/4d0/pi**2*x0mean(iw)
                else
c                   print *,'ttt: mmnorm=',mmnorm,q
c                   print *,'ttt: x0mean=',x0mean(iw)
c                   print *,'ttt: freq  =',frr
c                   print *,'ttt: schi  =',schi
                 x0mean(iw)       = x0mean(iw)/2d0/mmnorm**2 ! in Ry unit.
                 jcoup(iw,iqixc2) = 1d0/x0mean(iw) !x0mean= sum(svec*matmul(zxq(1:nbloch,1:nbloch,iw),svec))
                 write(ifchipmn,
     &            '(3f12.8,2x,f20.15,2x,2d23.15,2x,2d23.15)') 
     &            q, 2*schi*frr, jcoup(iw,iqixc2),x0mean(iw) 
c                write(ifchipmlog,'(3f12.8,x,f8.4," : ", (2f10.6,x))')
c     &         q, 2*freq_r(iw-1), jcoup(iw,iqixc2)
                endif
              else

C... ChiPM mode with LFC
C...... NoLFC part
c              mmnorm = sum(svec**2); write(6,"('mmnorm=svec*svec=', d14.6)") mmnorm
                zxq(1:ngb,1:ngb,iw) = zxq(1:ngb,1:ngb,iw)/2d0  ! in Ry.
                x0meanx= sum(svec*matmul(zxq(1:nbloch,1:nbloch,iw),svec))/mmnorm**2
!       x0meanx= <m|chi^+-(\omega)|m>/<m|m>**2
                jcoup(iw,iqixc2) = 1d0/x0meanx
                write(ifchipmn, '(3f12.8,2x,f20.15,2x,2d23.15
     &          ,3x,2d23.15)') q, 2*schi*frr, jcoup(iw,iqixc2), x0meanx 
!                <m|m>**2/<m|chipm|m> \approx  jcoup=<m|chipm^-1|m>

C...... With LFC
!  SVD of chi^-1: This procedure is necessary to remove inaccurarte high-value ones.
cx                zzz  = matmul(sqovlp, matmul(zw0, sqovlp) )
c                if(verbose()>50.or.iw<=2)
                 write(6,"(a,i5)")' ----SVD: chiinv --- iw=',iw
                call zgesvdnn(
     i            ngb, matmul(sqovlpi,matmul(zxq(1:ngb,1:ngb,iw),sqovlpi)),
     o            SS0,UU0,VT0)
                SS = 1d0/SS0
                VT = dconjg(transpose(UU0))
                UU = dconjg(transpose(VT0))

                ssm= ss(1)
                allocate(mcm(10));  mcm=0d0
                do ii= 1,ngb !ngb 
cx                ddd= ss(ii)*sum(svec(1:nbloch)*uu(1:nbloch,ii))*sum(vt(ii,1:nbloch)*svec(1:nbloch))
                  ddd= ss(ii)*sum(svec(1:nbloch)*uu(1:nbloch,ii)) 
     &                       *sum(vt(ii,1:nbloch)*svec(1:nbloch))
                  if(verbose()>50.or.iw<=2) then
                    write(6,"(' k eig_k <m|chiinv|m>=',
     &                      i4,d13.5,2d13.5 )" ) ii, ss(ii), ddd
                  endif
                  if(ii==1 ) mcm(1)=mcm(1) + ddd
                  if(ss(ii)< 2d0  +ssm) mcm(2)=mcm(2) + ddd
                  if(ss(ii)< 10d0 +ssm) mcm(3)=mcm(3) + ddd
                  if(ss(ii)< 50d0 +ssm) mcm(4)=mcm(4) + ddd
                  if(ss(ii)< 100d0+ssm) mcm(5)=mcm(5) + ddd
                                        mcm(6)=mcm(6) + ddd
                enddo
                if(iw==2) write(6,"(a)")
     &            '---SVD is not shown hereafter. Use Verbose>50---'
                write(ifchipm,
     &         '(3f12.8,2x,f20.15,2x,2d23.15,2x,2d23.15,2x,f9.5,x,f9.5)') 
     &           q, 2*schi*frr, mcm(1), mcm(2), jcoup(iw,iqixc2)/mcm(2)
                write(ifchipmlog,'(3f12.8,x,f8.4," : ", 7(2f10.6,x))')
     &           q, 2*schi*frr, mcm(2:6)
                deallocate(mcm)      
              endif
            endif
         enddo
         if(chipm.and.(.not.nolfco)) 
     &          deallocate(zzz,UU,VT,SS, sqovlp,sqovlpi,UU0,VT0,SS0)
         if( allocated(zzr)   ) deallocate(zzr)
         if( allocated(x0mean)) deallocate(x0mean)
         if( allocated(gbvec) ) deallocate(gbvec)
         if( ixc==5.or.ixc==6 ) then
            jpm=1
c            nwmax = nw
c            if(ixc==5) nwmax =nw
            allocate(trwv(nw_i:nw),trwv2(nw_i:nw))
            do iw = nw_i,nw !max ! trace check
              trwv(iw) = zxq(6,7,iw) 
              trwv2(iw) = 0d0
              do i = 1,ngb
              trwv2(iw) = trwv2(iw) + zxq(i,i,iw) 
              enddo
            enddo
            do iw= nw_i,nw-1
              if(ixc==5) 
     &         write(6,"('iq iw[min_max]=',2i5,2f7.4,' trwv by wwk*h= ',
     &         12d13.5)") iq, iw, freq_r(iw), freq_r(iw+1), 
     &         (trwv2(iw)+trwv2(iw+1))/2d0*(freq_r(iw)-freq_r(iw+1)),
     &         (trwv(iw)+trwv(iw+1))  /2d0*(freq_r(iw)-freq_r(iw+1))
                               !weight for the histgram range. by tetwt5
              if(ixc==6) 
     &         write(6,"('iq iw[min_max]=',2i5,2f7.4,' trwv by whw  = ',
     &         12d13.5)") iq,iw, freq_r(iw), freq_r(iw+1),
     &            trwv2(iw),trwv(iw) !weight for the histgram range. by tetwt5
            enddo  
            deallocate(trwv,trwv2)
         endif
        endif

C... Close files for epsmode
        if(epsmode) then !iepsmode/=0) then      ! only calculate iq>nqibz
          if(chipm) then
            ifchipmn = iclose ( 'ChiPM'//charnum4(iqixc2)//'.nolfc.dat')
            if(.not.nolfco) then
              ifchipm = iclose ('ChiPM'//charnum4(iqixc2)//'.dat')
            endif 
          else
            filepsnolfc ='EPS'//charnum4(iqixc2)//'.nolfc.dat'
            ifepsdatnolfc = iclose( filepsnolfc)
            if(.not.nolfco) then
              fileps = 'EPS'//charnum4(iqixc2)//'.dat'
                ifepsdat = iclose( fileps)
            endif
          endif
        endif
C--- RealOmega end ==============================


     
C--- ImagOmega ===================================
        if (imagomega) then
          print *,' goto imag omega'
          if (nspin == 1) zxqi = 2d0*zxqi   ! if paramagnetic, multiply x0 by 2
c ... Antiferro accelaration.
          if (anfexist) then
           do iw  = 1,niw 
           call anfx0k(natom,nclass,nblocha,iclass,pos
     i      ,nbloch,ngc, !ngcni(iqx),
c     i      q, ngveccB(1,1,iqx),qbas,    ! for q+G
     i      q, ngveccB,qbas,    ! for q+G
     i      anfvec,ibasf,  ! these are antiferro informations.
     i      zxqi(1,1,iw))  ! i/o 
           enddo
          endif
          if (ecorr_on>0)then !ixc==101.or.(sergeyv.and.imagonly)) then 
            imode=0
          else
            imode=1
          endif

          do iw  = 1,niw
            print *,'sss sumchk zxqi=',sum(abs(zxqi(:,:,iw)))
            if (lmlw) then
ctm, 060517
              allocate(zw_full(nblochpmx,nblochpmx)
     &                ,zw_full2(ngb,ngb))
              read(ifrcwi_full, rec=((iq-iqxini)*niw+iw) ) zw_full
              zw_full2 = zw_full(1:ngb,1:ngb)
              deallocate(zw_full)
              call wcf_maxloc( ngb, vcoul, zxqi(1,1,iw), zw_full2,
     &               imode,
     &               zw0)
              deallocate(zw_full2)
            else   
              call wcf( ngb, vcoul,zxqi(1,1,iw),imode, 
     &               zw0)
            endif  
            print *,'sss sumchk zw0=',sum(abs(zw0))
            zw(1:ngb,1:ngb) = zw0 ! zw(nblochpmx,nblochpmx)    
            write(ifrcwi, rec=((iq-iqxini)*niw+iw) ) zw    !  WP = vsc-v 

C--- for Miyake's total energy !Nov2004 ----------------
            if(ecorr_on>0) then
             if (debug) write(*,*)'ip,ix=',iq,iw,'  niw=',niw
             call getwk(iq, wibz, wqt,nqbz,nqibz,nstibz,nq0i, wk4ec)  
             call ecorq2 (vcoul, zw0, ngb, iq,iw,ieceig,
     o        erpaqw, trpvqw, trlogqw)
C--- Integration along imaginary axis.
c omit k and basis index for simplicity 
c wint = -(i/4pi) < [w'=-inf,inf] Q(w') >
c
c When w' ==> iw', w' is now real,
c   wint =  (1/2pi) < [w'=0,inf] Q(iw') >
c
c transform: x = 1/(1+w')
c this leads to a denser mesh in w' around 0 for equal mesh x
c which is desirable since Q is peaked around w'=0 
c    wint =  (1/2pi) < [x=0,1] Q(iw') / x^2 >
             faca   =  wk4ec* wiw(iw)
             trpv(iecut)    = trpv(iecut)  + faca* trpvqw
             trlog(iecut)   = trlog(iecut) + faca* trlogqw
             totexc(iecut)  = totexc(iecut)+ faca* erpaqw !  = trpv+ trlog
c            ecqw(iq,iw) = erpaqw
             if(iw==1) then
              write(ieclog,*)
             endif
             if(iw==1.and.iq==iqxini) then
               write(ieclog,
     &         "('   iq   iw   omega/i(Ry)        IntWgt',
     &         '    trpvqw(eV)      ecqw(eV)   ecqw*IntWgt',  
     &         ' :  ecut   ecuts')")
             endif
             write(ieclog,"( 2i5,3f14.6,3f14.6,2f8.3)") 
     &       iq,iw, 2d0*freq_i(iw), faca, trpvqw*hartree, erpaqw*hartree, 
     &       faca*erpaqw*hartree, ecut(iecut),ecuts(iecut)
             close(ieclog)
             open(ieclog,file="ecorr.chk",access='append')
ccccccccccccccccccccccccccccccccccccccccccccccccc
c           allocate( ovlpc(ngb,ngb),evall(ngb),
c     &          evecc(ngb,ngb))
c           evall=0d0
c           ovlpc=0d0
c           do i=1,ngb
c             ovlpc(i,i)=1d0
c           enddo
c           nmx=ngb
cc1           call diagcv(ovlpc,zw0/2d0+transpose(dconjg(zw0))/2d0,evecc,ngb, evall,nmx,1d99, nev)
c           call diagcv(ovlpc,zw0,evecc,ngb, evall,nmx,1d99, nev)
c           write(6,"('ngb nev=',2i5)") ngb,nev
c           write(6,"('chk eigen of zw0 Max Min=',2d13.6)")maxval(evall),minval(evall)
c           do i=1,3
c             write(6,*) i, evall(i)
c           enddo
c           do i=ngb-3,ngb
c            write(6,*) i, evall(i)
c           enddo
c           deallocate( ovlpc,evall,evecc)
ccccccccccccccccccccccccccccccccccccccccccccccc
            endif
            call tr_chkwrite("freq_i iq iw imgomg  trwv=",zw,iw,freq_i(iw),nblochpmx,nbloch,ngb,iq)
          enddo
        endif
C... ImagOmega end ===============================
c$$$ 2005 continue ! iecut loop


c 1002 continue  ! end of frequency block-loop
      deallocate(vcoul, zw0, zxq, zxqi)
 1001 continue  
C--- end of loop over q point =================================
c=======================================================================
 1100 continue


      if(epsmode.and.legas) then
         print *,' Find LEGAS. legas =',legas
         iflegas = 2101
         open (iflegas,file='LEGAS')
         read(iflegas,*)rs
         close(iflegas)
         alpha  = (9*pi/4d0)**(1d0/3d0)
         qfermi = alpha/rs
         efx  = qfermi**2
         valn = efx**1.5d0*voltot/3d0/pi**2
         write (6,*)'  #### egas test mode  legas=T #### given rs =',rs
         write (6,*)'     Exact Fermi momentum  qf  =', qfermi
         write (6,*)'     Exact Fermi energy    Ef  =', efx
         do iq = iqxini,iqxend ! q=(0,0,0) is omitted!
          if(iq<=nqibz) cycle
          print *,' iq=',iq
           iqixc2 = iq- (nqibz+nq0ix)
           filele ='EPSEG'//charnum4(iqixc2)//'.dat'
              ife = iopen ( filele,1,3,0)
              write(ife,"(a)")
     &          ' q(1:3)   w(Ry)   eps    epsi  --- NO LFC'
           q = qibze(:,iq)
           qt= sqrt(sum(qibze(1:,iq)**2))*2d0*pi/alat
           qs= qt/qfermi
           write(6,"(' qs qfermi=',2d13.5)"    ) qs,qfermi
           write(6,"(' q-q^2/2 q+q^2=',2d13.5)") qs-qs**2/2d0,qs+qs**2/2d0
         do iw  = nw_i,nw
            ww  = freq_r(iw)
            muu = ww/qfermi**2
            if(     qs<2d0 .and. muu < qs-qs**2/2d0) then
              x0meanx= -img*qfermi/(4*pi*qs)*2*muu
            elseif( qs<2d0 .and. muu < qs+qs**2/2d0) then
              x0meanx= -img*qfermi/(4*pi*qs)*( 1d0-(muu/qs-.5d0*qs)**2 )
            else
              x0meanx=0d0 
            endif
            vcmmmm= 4*pi/qt**2
            epsi(iw,iqixc2) = 1d0/(1- vcmmmm * x0meanx) 
c            epsi(iw,iqixc2) = 1d0/(1- vcmmm(iq) * x0meanx) 
           write(ife,'(3f12.8,2x,d12.4,2d23.15,2x,2d23.15)') 
     &        q, 2*ww,1d0/epsi(iw,iqixc2),epsi(iw,iqixc2)
         enddo
         enddo
         print *,' ----------legas end--------'
      endif

C--- Write TEECOR ecorr_on mode
      if(imagomega.and.ecorr_on>0) then
        hartree=2d0*rydberg()
        ifcor   = iopen('TEECORR2',1,-1,0) ! output files
        do iecut=1,necut
         write(6,"( ' RPA Ec =' 3f23.15,'   ecut ecuts (Ry)=',2d12.4)") 
     &   totexc(iecut)*hartree,trpv(iecut)*hartree, trlog(iecut)*hartree 
     &    ,ecut(iecut),ecuts(iecut)
         write(ifcor,*) '============================'
         write(ifcor,*) 'Correlation energy Erpa (eV)'
         write(ifcor,*) '============================'
         write(ifcor,*)' *** '
         write(ifcor,"(5d23.15)") 
     &     totexc(iecut)*hartree,trpv(iecut)*hartree,trlog(iecut)*hartree
     &    ,ecut(iecut),ecuts(iecut)
        enddo
C... output ecqw !    write(ifcor,*)'*** ecqw(q,w) ***'
        write(ifcor,*)' nqibz =',nqibz
        write(ifcor,*)' nq0i  =',nq0i
        write(ifcor,*)' niw   =',niw
        write(ifcor,*)' --- See details of Ec in ecor.chk ---'
c        nqitot = nqibz + nq0i
c        call wecqw(ifcor,
c     d           nqibz,nqbz,nq0i,nqitot,niw,
c     o           wibz,wqt,wx,freqx,ecqw)
C... Write electron gas correlation energy 
        legas = .false.
        INQUIRE (FILE = 'LEGAS', EXIST = legas)
        if(legas) then !!! test for electron gas case. 
         print *,' find LEGAS. legas =',legas
         iflegas = 2101
         open (iflegas,file='LEGAS')
         read(iflegas,*)rs
         close(iflegas)
         alpha = (9*pi/4d0)**(1d0/3d0)
         qfermi = alpha/rs
         efx  = qfermi**2
         valn = efx**1.5d0*voltot/3d0/pi**2
         write (6,*)'  #### egas test mode  legas=T #### given rs =',rs
         write (6,*)' egas  Exact Fermi momentum  qf  =', qfermi
         write (6,*)' egas  Exact Fermi energy    Ef  =', efx
         if(tetra) stop 'legas You have to give ef of  tetrahedron'
         efz=(ntot*3*pi**2/voltot)**(2d0/3d0) ! ef is calculated from ntot.
         qfermi= dsqrt(efz)
         alpha = (9*pi/4d0)**(1d0/3d0)
         rs    = alpha/qfermi
         write (ifcor,*)' --- electron gas ---'
         write (ifcor,*)' density parameter rs= ', rs
         write (ifcor,*)' kf= ',qfermi
         write (ifcor,*)' *** Barth-Hedin formula'
         ecelgas = eclda_bh(rs) * hartree * ntot
         write (ifcor,*)ecelgas
         write (ifcor,*)' *** Perdew-Zunger formula'
         ecelgas = eclda_pz(rs) * hartree * ntot
         write (ifcor,*)ecelgas
         write (ifcor,*)' *** Gell-Mann and Brueckner formula'
         ecelgas = (-0.0311d0 * dlog(rs) -0.048d0) * hartree * ntot
         write (ifcor,*)ecelgas
        endif
      endif
C...close files
c        close(ifdpin)
      is = iclose('hbe.d')
      call cputid(0)
      if(ixc==101) stop ' OK! hx0fp0 mode=101 Ecor '
      if(ixc==111) stop ' OK! hx0fp0 mode=111 maxloc-W '
      if(ixc==1)   stop ' OK! hx0fp0 mode=1 normal'
      if(ixc==2)   stop ' OK! hx0fp0 mode=2  eps NoLFC'
      if(ixc==22)  stop ' OK! hx0fp0 mode=22 chi+- NoLFC'
      if(ixc==23)  stop ' OK! hx0fp0 mode=23 chi+- LFC  '
      if(ixc==3)   stop ' OK! hx0fp0 mode=3 eps with LFC'
      if(ixc==202) stop ' OK! hx0fp0 mode=202 sergeyv epsPP NoLFC'
      if(ixc==203) stop ' OK! hx0fp0 mode=203 sergeyv eps LFC '
      if(ixc==222) stop ' OK! hx0fp0 mode=222 chi+- NoLFC sergeyv'
      if(ixc==10222)stop' OK! hx0fp0 mode=10222 <q|chi+-|q> NoLFC sergeyv'
      if(ixc==223) stop ' OK! hx0fp0 mode=223 chi+- LFC   sergeyv'
      if(ixc==11)  stop ' OK! hx0fp0 mode=11 normal sergeyv'
      if(ixc==12)  stop ' OK! hx0fp0 mode=12 Ecor sergeyv mode'
      if(ixc==4)   stop ' OK! hx0fp0 ixc=4 test'
      if(ixc==5)   stop ' OK! hx0fp0 ixc=5 test'
      if(ixc==6)   stop ' OK! hx0fp0 ixc=6 test'
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c--------------------------------------------------------------------
      real*8 function eclda_bh(rs)
      real(8) :: rs,cp,rp,z

      cp       = 0.0504d0*0.5d0 ! 0.5 changes unit from Ry to Hartree
      rp       = 30.d0
      z        = rs / rp
      eclda_bh = -cp * ( (1.d0+z**3)*dlog(1.d0+1.d0/z) 
     .                 + 0.5d0*z - z**2 - 0.33333333d0 )
      end
c--------------------------------------------------------------------
      real*8 function eclda_pz(rs)
      real(8) :: rs
      if (rs.ge.1.d0) then
      eclda_pz = -0.1423d0 / (1.d0 + 1.0529d0*dsqrt(rs) + 0.334d0*rs)
      else
      eclda_pz = -0.0480d0 + 0.0311d0*dlog(rs) - 0.0116d0 * rs 
     .         + 0.0020d0*rs*dlog(rs)
      endif
      end
c--------------------------------------------------------------------
      subroutine wecqw(ifcor,
     d                 nqibz,nqbz,nq0i,nqitot,niw,
     o                 wibz,wqt,wx,freqx,ecqw)
      
      implicit double precision (a-h,o-z)

      dimension   wibz(nqibz),wqt(nq0i),wx(niw),
     .            freqx(niw),ecqw(nqitot,niw)
      real(8):: rydberg
      write(ifcor,*)'*** ecqw(q,w) ***'
      write(ifcor,*)'nqibz =',nqibz
      write(ifcor,*)'nq0i  =',nq0i
      write(ifcor,*)'niw   =',niw

      do ip = 2,nqitot
        if (ip <= nqibz) then 
        wk = wibz(ip)*0.5d0 ! 0.5 for the normalization of wibz
        else
c        wk = wqt(ip-nqibz)*wibz(1)*0.5d0 ! 0.5 for the normalization of wibz
        wk = wqt(ip-nqibz)* 1d0/dble(nqbz)
        endif
        write(ifcor,*)'*** iq,wq = ',ip,wk

      sume=0d0
      do ix = 1,niw
      write(ifcor,*)freqx(ix),ecqw(ip,ix),wx(ix)
      sume=sume+  wx(ix)/(freqx(ix)*freqx(ix)) * ecqw(ip,ix)
      end do
      write(ifcor,*) '  sum ecqw*wx=', wk*sume*2d0*rydberg()

c end of ip-loop
      end do
      return
      end

c---------------------------------------------------------------------
       subroutine getsqovlp(q,ngc,ngb,sqovlp)
C  ... Get ppovl
       implicit none
       real(8)::q(3)
       integer(4):: ngc,ngb,nbloch,i,nmxx,ix,iy,nev
       complex(8):: sqovlp(ngb,ngb)
       complex(8),allocatable:: ooo(:,:),ppo(:,:),sqovlpi(:,:),ppovl(:,:)
       complex(8),allocatable:: ovlp(:,:),evec(:,:)
       real(8),allocatable:: eval(:)
       nbloch = ngb-ngc
       allocate(ppovl(1:ngc,1:ngc))
       call readppovl0(q,ngc,ppovl)
       allocate(ooo(ngc,ngc),ppo(ngc,ngc),evec(ngc,ngc),eval(ngc))
       ooo= 0d0
       do ix=1,ngc
         ooo(ix,ix)=1d0
       enddo   
       ppo = ppovl
       deallocate(ppovl)
       nmxx = ngc
       evec = 0d0
       eval = 0d0
       call diagcv(ooo, ppo,
     &     evec, ngc, eval, nmxx, 1d99, nev)
c         print *,' diagcv overlap ngc nev=',ngc,nev
       sqovlp=0d0
       do i=1,nbloch
         sqovlp(i,i)=1d0
       enddo   
       do i=1,ngc
         if(eval(i)<0d0) then
            stop 'getsqovlp:  eval(i) <0d0'
         endif   
         do ix=1,ngc;  do iy=1,ngc
           sqovlp(ix+nbloch,iy+nbloch)= 
     &     sqovlp(ix+nbloch,iy+nbloch)
     &      + evec(ix,i)* sqrt(eval(i))* dconjg(evec(iy,i))
         enddo ;      enddo
       enddo  
       deallocate(ooo,ppo,evec,eval)
c         sqovlpi = sqovlp
c         call matcinv(ngb,sqovlp)     !  inverse
c         ovlpi=ovlp
c         deallocate(ppovl,ovlp)
       end
c---------------------------------------------------------------------
       subroutine readppovl0(q,ngc,ppovl)
       implicit none
       integer(4):: ngc_r,ippovl0,ixx,ngc
       complex(8):: ppovl(ngc,ngc)
       real(8):: q(3),qx(3)
       ippovl0=2301
       open(ippovl0,file='PPOVL0',form='unformatted')
       do 
        read(ippovl0) qx,ngc_r
        if(sum(abs(qx-q))<1d-6) then
c           print *,' q qx ngc ngc_r=',q, qx, ngc, ngc_r
           if(ngc_r/=ngc) stop 'readin ppovl: ngc_r/=ngc'
           read(ippovl0) ppovl
c           print *,' readin PPOVL0 OK q=',q
           exit
        endif  
      enddo
      close(ippovl0)
      end

      subroutine zgesvdnn(ngb,zzz, SS,UU,VT)
C--- SVD zzz= matmul(UU,matmul(SS,VT)) ------------
c$$$C--- SVD of chipm test !SVD procedure is not required to calculate <m|chi|m>
c$$$      lwork=4*ngb
c$$$      allocate(work(LWORK),rwork(5*ngb),zzz(ngb,ngb),UU(ngb,ngb),VT(ngb,ngb),VTT(ngb,ngb),ss0(ngb))
c$$$      zzz= matmul(transpose(conjg(ovlpi)), matmul(zxq(1:ngb,1:ngb,iw), ovlpi) )
c$$$      call zgesvd('A','A',ngb,ngb,zzz,ngb,SS0,UU,ngb,VT,ngb,work,lwork,rwork,info) 
c$$$      write(6,*)' -------SVD: Oinv*chipm*Oinv ----------------'
c$$$      aaax = 0d0
c$$$      do i=1,ngb
c$$$        addx= sum(svec(1:nbloch)*uu(1:nbloch,i)) *ss0(i)* sum(VT(i,1:nbloch)*svec(1:nbloch)) 
c$$$        write(6,"(' SVD OcO: eig_k <m|chi|m>_k=',i4,2x, d13.5,2x,2d14.6)")i,SS0(i),addx
c$$$        if(i<25) aaax= aaax+ addx
c$$$      enddo 
c$$$      aaax= mmnorm**2/aaax
c$$$      deallocate(work,rwork,zzz,uu,vt,vtt)
c$$$      deallocate(ovlpi)
      implicit none
      integer(4)::lwork,info,ngb,i
      complex(8):: zzz(ngb,ngb),UU(ngb,ngb),VT(ngb,ngb)
      real(8):: ss(ngb)
      real(8),allocatable:: rwork(:)
      complex(8),allocatable:: work(:),zw0bk(:,:),vtt(:,:)
      lwork=4*ngb
      allocate(zw0bk(ngb,ngb))
      allocate(work(LWORK),rwork(5*ngb)) !,VTT(ngb,ngb))
      zw0bk = zzz
!      print *,' zgesvdnn: singular value decomp '
      call zgesvd('A','A',ngb,ngb,zzz,ngb,SS,UU,ngb,VT,ngb,work,lwork,rwork,info) 
!      do i=1,ngb
!         write(6,"(' i ss=',i4,' ', d13.5 )")i,SS(i) !    write(6,"(' i ss=',i4,'  ', d13.5,' ss0*ss=',d13.5 )")i,SS(i),ss(i)*ss0(ngb-i+1)
!         vtt(i,:)=ss(i)*vt(i,:)
!      enddo
!      write(6,"('sumcheck zzz  zzz-uu*s*vt=',d13.5,d13.5)")
!     &  sum(abs(zw0bk)), sum(abs(zw0bk - matmul(uu,vtt)))
!      if(abs(sum(abs(zw0bk - matmul(uu,vtt))))>1d-8*sum(abs(zw0bk))) 
!     &  stop 'sumcheck zzz  zzz-uu*s*vt= error'
!      deallocate(vtt)
      end

C===================================================================
              subroutine tr_chkwrite(tagname,zw,iw,freqq,nblochpmx,nbloch,ngb,iq)
              implicit none
              integer(4):: nblochpmx,nbloch,ngb,iw,i,iq
              complex(8):: zw(nblochpmx,nblochpmx),trwv,trwv2
              real(8):: freqq
              logical :: smbasis
              character*(*)::tagname
              trwv=0d0
              if(.not.smbasis()) then
                do i = 1,nbloch
                trwv = trwv + zw(i,i) 
                enddo
              endif
              trwv2 = 0d0
              do i = 1,ngb
              trwv2 = trwv2 + zw(i,i) 
              enddo !  write(6,'(" realomg trwv=",2i6,4d22.14)') iq,iw,trwv(iw),trwv2(iw)
              write(6,'(a,f10.4,2i5,4d22.14)')tagname,freqq,iq,iw,trwv,trwv2
              end
C===================================================================
      subroutine get_nwf(nwf)
      implicit none
      integer(4):: nwf,nqbz,iko_ix,iko_fx,ifi
      integer(4):: iopen,iclose

         ifi  = iopen('MLWU',0,0,0)
         read(ifi)nqbz,nwf,iko_ix,iko_fx
         ifi = iclose('MLWU')
         end
C===================================================================
      subroutine choose_wanband(iwbnd,nwf,nqbze,nspin)
      implicit none
      integer(4):: nwf,nqbze,nspin,iwbnd(nwf,nqbze,nspin)
      integer(4):: iopen,iclose

c local
      integer(4):: is,ifmlw,ifuu,nq0i,nwf2,nqbz,nqbz2,iko_ix,iko_fx,
     &             iko_ix2,iko_fx2,iqbz,iqbz2,iq0i,iq0i2,j1,j2,
     &             ikp,ib,iwf,nbnd
      integer(4),allocatable:: isort(:)
      real(8):: q(3)
      real(8),allocatable:: wbnd(:)
      complex(8),allocatable:: dnk(:,:,:,:),uum(:,:,:,:,:),cbwf(:,:)


c --- Readin MLWU/D and UUq0U/D
      do is = 1,nspin

c fileopen
         if (is.eq.1) then
           ifmlw  = iopen('MLWU',0,0,0)
           ifuu   = iopen('UUq0U',0,0,0)
         else ! is
           ifmlw  = iopen('MLWD',0,0,0)
           ifuu   = iopen('UUq0D',0,0,0)
         endif ! is

c nqbz mesh-points
         read(ifmlw)nqbz,nwf2,iko_ix,iko_fx
         if (nwf2.ne.nwf) stop 'choose_wanband: nwf error'
!         if (nqbz.ne.nqbze) stop 'choose_wanband: nqbz error'
         if (is.eq.1) allocate(dnk(iko_ix:iko_fx,nwf,nqbz,nspin))
         do iqbz = 1,nqbz
            read(ifmlw)iqbz2,q(1:3)
            if (iqbz2 .ne. iqbz) stop 'choose_wanband: iqbz error'
            read(ifmlw)dnk(iko_ix:iko_fx,1:nwf,iqbz,is)
         enddo ! iqbz

c shifted mesh points
         read(ifuu)
         read(ifuu)nqbz2,nq0i,iko_ix2,iko_fx2
         if (is.eq.1) 
     >     allocate(uum(iko_ix:iko_fx,iko_ix:iko_fx,nqbz,nq0i,nspin))
         if (nqbz2.ne.nqbz) stop "choose_wanband: nqbz2 error"
         if (iko_ix2.ne.iko_ix)
     >      stop "choose_wanband: iko_ix2 error"
         if (iko_fx2.ne.iko_fx)
     >      stop "choose_wanband: iko_fx2 error"
         do iqbz = 1,nqbz
         do iq0i =1,nq0i  
            read(ifuu)
            read(ifuu)iqbz2,iq0i2
            if (iqbz2 .ne. iqbz) stop 'choose_wanband: iqbz error'
            if (iq0i2 .ne. iq0i) stop 'choose_wanband: iq0i error'
            read(ifuu)
     &      ((uum(j1,j2,iqbz,iq0i,is),
     &         j1=iko_ix,iko_fx),j2=iko_ix,iko_fx)
         enddo ! iq0i
         enddo ! iqbz

c fileclose
         if (is.eq.1) then
            ifmlw  = iclose('MLWU')
            ifuu   = iclose('UUq0U')
         else ! is
            ifmlw  = iclose('MLWD')
            ifuu   = iclose('UUq0D')
         endif ! is

      enddo ! is

      allocate(cbwf(iko_ix:iko_fx,nwf),wbnd(iko_ix:iko_fx),
     &         isort(iko_fx-iko_ix+1))
      do ikp = 1,nqbze
         iqbz = mod(ikp,nqbz)
         if (iqbz .eq. 0) iqbz = nqbz
         iq0i = (ikp - iqbz)/nqbz
      do is = 1,nspin
         cbwf = 0d0
         if (iq0i.eq.0) then
             cbwf(:,:) = dnk(:,:,iqbz,is)
         else ! iq0i
c   <psi(k+q0,n) | psi(k+q0,m)^W> 
c = S[l] <psi(k+q0,n) |e^(iq0.r)| psi(k,l)>
c      * <psi(k,l) |e^(-iq0.r)| psi(k+q0,m)^W>  
c ~ S[l] <psi(k+q0,n) |e^(iq0.r)| psi(k,l)> <psi(k,l) |psi(k,m)^W>  
c
c psi^W : bloch fn. in the Wannier gauge
           do ib = iko_ix,iko_fx
           do iwf= 1,nwf
             cbwf(ib,iwf) =
     >            sum( conjg(uum(iko_ix:iko_fx,ib,iqbz,iq0i,is))
     >                      *dnk(iko_ix:iko_fx,iwf,iqbz,is) )
           enddo ! iwf
           enddo ! ib
         endif ! iq0i

c choose bands
         wbnd = 0.0d0
         nbnd = iko_fx - iko_ix + 1
         do ib = iko_ix,iko_fx
         do iwf = 1,nwf
           wbnd(ib) = wbnd(ib) + abs(cbwf(ib,iwf))**2
         enddo ! iwf
         enddo ! ib   
         call sortbnd(wbnd(iko_ix:iko_fx),iko_fx-iko_ix+1,
     &                isort)
         do iwf = 1,nwf
            iwbnd(iwf,ikp,is) = isort(iwf) + iko_ix - 1
         enddo ! iwf
      enddo ! is
      enddo ! ikp
      deallocate(dnk,uum,cbwf,wbnd,isort)

      end
C===================================================================
      subroutine sortbnd(rin,n,isort)
      implicit none
      integer(4):: n,isort(n)
      integer(4):: i,j,itmp
      real(8):: rin(n),r(n),rtmp

      do j = 1,n
         isort(j) = j
      enddo ! j   

      r(:) = rin(:)
      do j = 2,n
        rtmp = r(j)
        do i = j-1,1,-1
           if (r(i).ge.rtmp) goto 999
           r(i+1) = r(i)
           isort(i+1) = isort(i)
        enddo ! i
        i = 0
 999    continue
        r(i+1) = rtmp
        isort(i+1) = j
      enddo ! j

check
!      do j = 1,n
!         write(*,*)j,isort(j),rin(isort(j))
!      enddo
      do j = 1,n-1
         if (rin(isort(j)).lt.rin(isort(j+1)))
     &   stop "hx0fp0: sortbnd error"
      enddo   

      end
C===================================================================
