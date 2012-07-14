      program hx0fp0_sc
c-------------------------------------------------
c  Calculate x0, W-V, \epsilon
Cr
Cr Be careful for the indexing...
Cr      A routine idxlnmc(nindxv,nindxc,...  in index.f
Cr      specifies the order of the  (Core wave)+(Argumentation wave) in each MT.
Cr      The total number of the wave are mnl(ic)= mnlc(ic) + mnlv(ic).
Cr      The indexing starts with core first and then valence on top of core
Cr      So n-index in "in" for valence electron is different from "inv".
Cr      Search in and inv.
Cr the key routine is x0k_v2.
c
c-------------------------------------------------
      use m_readqg,only: readngmx,readqg
      use m_readeigen,only: init_readeigen,init_readeigen2,readeval
      use m_read_bzdata,only: read_bzdata,
     &  ngrp2=>ngrp,nqbz,nqibz,nqbzw,nteti,ntetf,n1,n2,n3,qbas,ginv,qbasmc,
     &  dq_bzcase2,qbz,wbz,qibz,wibz,qbzw,
     &  idtetf,ib1bz,idteti,
     &  nstar,irk,nstbz
      use m_genallcf_v3,only: genallcf_v3,
     &  nclass,natom,nspin,nl,nn,nnv,nnc, ngrp,
     &  nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, nctot,niw,nw,
     &  alat,ef, diw,dw,delta,deltaw,esmr,symgrp,clabl,iclass,nindxv,nindxc,ncwf,
     &  invg, il, in, im, ilnm, nlnm, ilv,inv,imv,  ilnmv, nlnmv,
     &  ilc,inc,imc,  ilnmc, nlnmc,nindx,konf,icore,ncore,
     &  occv,unoccv ,occc,unoccc, nocc, nunocc, plat, pos,z,ecore, freq, symgg ! symgg=w(igrp)
      use keyvalue,only: getkeyvalue

      implicit none
      integer(4)::iopen,maxocc2,iclose,
     & ixc,iqxini,iqxend,
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
     &  ,    nw_w,nwmin

      real(8):: dum1,dum2,dum3,wqtsum,epsrng,dnorm,
     & dwry,dwh,omg_c,omg2

      real(8):: q(3),  qgbin(3),qx(3)
      real(8):: ua=1d0 ! this is a dummy.

c-------------------------------------------------
c worksize in megabytes (1 word = 4 bytes)
c      parameter (mbytes=50)
c      parameter (mwords=mbytes/4)
c      parameter (iwksize=mwords * 1000 *1000)
c work array w is an integer array
c      integer w
c      common /w/ w(iwksize)
c class parameters
c      parameter (mxclass=100)
c      character*6 clabl(mxclass)
c symmetry group
c      parameter (ngnmax=10)
c      dimension gen(9,ngnmax)
c      character*120 symgrp

      integer(4):: ifrb(2),ifcb(2),ifrhb(2),ifchb(2),ifev(2)
      integer(4):: ndble=8
      integer(4):: nword

c> for genallc.f
c      dimension plat(3,3) !,qbas(3,3),ginv(3,3)
c     l          ifphiv(2),ifphic(2),
c     l          iphi(2),iphidot(2),ipgs(2),ipgsdot(2),
c     l          iphiv(2),iphivd(2),ipgv(2),ipgvd(2),
c     l          iphic(2),ipgc(2)
c
      real(8),allocatable:: vxcfp(:,:),
     &  wqt(:), wgt0(:,:),q0i(:,:)
      integer(4),allocatable :: ngveccB(:,:)
     &           , ngvecc(:,:),iqib(:),
     &   nx(:,:),nblocha(:),lx(:),ifppb(:)
      complex(8),allocatable:: geigB(:,:,:,:) ,geig(:,:),vcoul(:,:),
     &  zw(:,:),zw0(:,:),
     &  zxq(:,:,:),zxqi(:,:,:)
c
c      complex(4), allocatable :: zw4(:,:)
c      complex(8), allocatable :: zw4(:,:) !takao
c
      real(8),allocatable :: ppbrd (:,:,:,:,:,:,:),cgr(:,:,:,:),eqt(:),
     &  ppbrdx(:,:,:,:,:,:,:),aaa(:,:),symope(:,:), !qibz(:,:),
     &  ppb(:,:),pdb(:,:),dpb(:,:),ddb(:,:), 
     &  qbze(:,:),qibze(:,:)  !,ecore(:,:)
c     &  freqr(:),freqi(:) !rw(:,:),cw(:,:) --->zw
      complex(8),allocatable :: rcxq(:,:,:,:)

c  tetra
      logical :: tetra,tmpwwk=.false.! If tmpwwk=.true., this use a temporary file tmp.wwk
      ! so as to reduce the memory usage.
      complex(8) :: fff,img=(0d0,1d0)
      complex(8),allocatable :: wwk(:,:,:)
      integer(4),allocatable :: 
     &         n1b(:,:,:),n2b(:,:,:),nbnb(:,:),nbnbtt(:,:),noccxvv(:)
c     &         idtetf(:,:),ib1bz(:)
c
      real(8),allocatable:: ekxx1(:,:),ekxx2(:,:) !qbzw(:,: )
      real(8) ::qbzx(3),anfvec(3)
      logical :: anfexist
      logical :: debug
      integer(4),allocatable:: ibasf(:)
      real(8),allocatable :: transaf(:,:)

      real(8),allocatable :: freq_r(:),freq_i(:),freqx(:),wx(:),expa(:)
     &                      ,frhis(:)
      logical :: realomega, imagomega
      real(8),allocatable:: ppovl(:,:)
      complex(8),allocatable:: zzr(:)
      complex(8) :: epxxx,vcmean
      character*9 fileps
      character*15 filepsnolfc
      logical :: paralellx0=.true. !, hist
c      character(5) :: charnum5
c      character(8) :: xt
      character(20):: xxt

c faleev
      real(8) :: Emin, Emax      ,emax2,emin2
      integer(4) :: iSigma_en  !sf..21May02  !iSigma_en is integer
      !parameter stored in GWIN_V2
      !which determines approximation for  self-energy.
      !Self-energy should be made hermitian for energies to be real
      !iSigma_en==0 SE_nn'(ef)+img integral:delta_nn'([SE_nn(e_n)+c.c.]/2-SE_nn(ef))
      !iSigma_en==1 SE_nn'(ef)+delta_nn'([SE_nn(e_n)+c.c.]/2-SE_nn(ef))
      !iSigma_en==2 [SE_nn'((e_n+e_n')/2)+h.c.]/2
      !iSigma_en==3 [(SE_nn'(e_n)+SE_nn'(e_n'))/2+h.c.]/2
      real(8) :: omg2max,omg1max
      real(8), allocatable :: freqr2(:)  , ekxxx(:,:,:)

c for tetwt5
c      real(8),allocatable :: demin(:,:,:),demax(:,:,:),whw(:)
      real(4),allocatable :: demin(:,:,:,:),demax(:,:,:,:)
      real(8),allocatable :: whw(:)

      integer(4),allocatable:: ihw(:,:,:),nhw(:,:,:),jhw(:,:,:),ibjb(:,:,:,:)

      logical::imagonly=.false. , noq0p !,readgwinput

      integer(4)::nwin, incwfin, verbose
      real(8)::efin
c      integer(4) ::iepsmode=0
      logical :: nolfco=.false.
      complex(8) :: x0mean,rcxqmean !dummy

      integer(4):: bzcase, isp1,isp2, ngc,mrecg
      real(8)::  qq(3),quu(3),deltaq(3),qqq(3)=0d0

      logical,allocatable :: iwgt(:,:,:,:)
      complex(8),allocatable:: wgt(:,:,:)

      real(8),allocatable:: qbz2(:,:)
      logical :: qbzreg
      character(2):: soflag
      integer(4):: nbcut,nbcut2

      logical ::smbasis
      integer(4):: ifpomat,nnmx,ikpo,nn_,noo,iqxxx,nomx
      real(8):: q_r(3)
      complex(8),allocatable:: pomat(:,:)


C... Feb2006 time-reversal=off case
      logical   :: timereversal,onceww
      integer(4):: npm ,jpm,ncc,nw_i
      real(8)   :: frr

      real(8):: sciss

c      logical :: negative_testtr,testtr
!-------------------------------------------------------------------------
C --- computational mode select
      write(6,*) ' --- hx0fp0_sc Choose omodes below ----------------'
      write(6,*) ' Only for ixc==11 (ixc==1 may work---not checked)'
      write(6,*) ' --- Put number above ! -----------------'
      call readin5(ixc,iqxini,iqxend)
      if(ixc==0) stop ' --- ixc=0 --- Choose computational mode!'
      if(iqxini ==0 ) then
        paralellx0=.false.
        write(6,"(' ixc = ',i4)") ixc
      else
        paralellx0=.true.
        write(6,"(' PARALELL mode: ixc iqxini iqxend =',3i5)")
     &  ixc, iqxini, iqxend
      endif
      call headver('hx0fp0_sc',ixc)

      call getkeyvalue("GWinput","ScaledGapX0",sciss,default=1d0)
      write(6,"(' ScaledGapX0=',f5.3)") sciss

!      print *, xt(iqxini), 'to',xt(iqxend); print *, xxt(iqxini,iqxend), 'xxt '
c      if    (ixc==1) then; write(6,*) ' OK ixc=1. normal mode'
ccccccccccccccccccccccccc sf 21May02 ccccccccccccccccccccccccccccc
c      elseif(ixc==2) then; write(6,*) ' OK ixc=2. eps mode. ONLY NoLFC'
c      elseif(ixc==3) then; write(6,*) ' OK ixc=3. eps mode'
c      elseif(ixc==4) then
c         write(6,*)            ' OK ixc=4. test. mode tetwt4.vs.tetwt5.'
c      elseif(ixc==5) then; write(6,*) ' OK ixc=5. Real-axis tetwt4 mode'
c      elseif(ixc==6) then; write(6,*) ' OK ixc=6. Real-axis tetwt5 mode'
c      elseif(ixc==11) then; write(6,*) " OK ixc=11 s.f.'s normal mode "
c      elseif(ixc==12) then
c        write(6,*) " OK ixc=12 s.f.'s normal mode only---WVI "
c        ixc= 11
c        imagonly=.true.
      if    (ixc==1)  then; write(6,*) ' OK ixc=1. normal mode'
      elseif(ixc==11) then; write(6,*) " OK ixc=11 s.f.'s normal mode "
      endif

c---  readin BZDATA. See gwsrc/rwbzdata.f
c--------readin data set when you call read_BZDATA ---------------
c       integer(4)::ngrp,nqbz,nqibz,nqbzw,nteti,ntetf,
c     &   n_index_qbz
c       integer(4):: n1,n2,n3
c       real(8):: qbas(3,3),ginv(3,3),qbasmc(3,3),dq_bzcase2
c       real(8),allocatable:: qbz(:,:),wbz(:),qibz(:,:)
c     &    ,wibz(:),qbzw(:,:)
c       integer(4),allocatable:: idtetf(:,:),ib1bz(:),idteti(:,:)
c     &    ,nstar(:),irk(:,:),index_qbz(:,:,:)
c-----------------------------------------------------------------
      call read_BZDATA()

c--- Use regular mesh even for bzcase==2 and qbzreg()=T
c         off-regular mesh for bzcase==1 and qbzreg()=F
      if( ( bzcase()==2.and.qbzreg() )       .or.
     &    ( bzcase()==1.and.(.not.qbzreg()))      ) then
        deltaq= qbas(:,1)/n1 + qbas(:,2)/n2 +qbas(:,3)/n3
        do i=1,nqbz
          qbz(:,i) = qbz(:,i) - deltaq/2d0
          write(6,"('i qbz=',i3,3f8.4)") i,qbz(:,i)
        enddo
      endif
      print *,' nqbz qbz =',nqbz
      print *,' nqibz ngrp=',nqibz,ngrp
c      print *,  qbz
c      print *,' irk=',irk
c      print *,' #### idtetf: ####'
c      print *, idtetf

C --- readin compuational conditions
c      call wkinit (iwksize)
      nwin   = 0  !Readin nw from NW file
      if(ixc/=1.and.ixc/=11) nwin = -999 !Not readin NW file
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
c------------------------------------------------------------------
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


      call readd_iSigma_en(ifinin,iSigma_en) !sf..21May02
c ...
      nw_input = nw ; print *, 'nw delta=',nw_input,delta
c      if (nclass > mxclass) stop ' hx0fp0: increase mxclass'

      debug=.false.
      if(verbose()>=100) debug=.true.

      if(debug) print *,' end of genallc'


!!!! WE ASSUME iclass(iatom)= iatom !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  We assume nclass = natom.     !!!!!!!!!!!!!!!!!!!!!!!!!
      if(nclass /= natom) stop ' nclass /= natom '

c --- tetra or not
      if(delta <= 0d0) then
        tetra =  .true.
        delta = -delta
        print *,' hx0fp0: tetrahedron mode delta=',delta
      else
        tetra = .false. ! switch for tetrahedron method for dielectric functions
      endif

c --- read dimensions of h,hb
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
c      if(ixc==1.and.abs(wqtsum-1d0) >1d-10) stop ' wqtsum of Q0P /=1 '
      call getkeyvalue("GWinput","TestNoQ0P",noq0p,default=.false.)
      if(ixc==1.and.abs(wqtsum-1d0) >1d-10.and.(.not.noq0p)) 
     &  stop ' wqtsum of Q0P /=1 '

c --- readin by rdpp ; Radial integrals ppbrd and plane wave part
c      call getsrdpp( nclass,nl,
c     o               ngpmx,ngcmx,nxx )
      call getsrdpp2( nclass,nl,nxx)
      call readngmx('QGpsi',ngpmx)
      call readngmx('QGcou',ngcmx)
      print *,' ngcmx ngpmx=',ngcmx,ngpmx

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
      allocate( nx(0:2*(nl-1),nclass), nblocha(nclass) ,lx(nclass),
     &   ppbrd ( 0:nl-1, nn, 0:nl-1,nn, 0:2*(nl-1),nxx, nspin*nclass),
     &   cgr(nl**2,nl**2,(2*nl-1)**2,ngrpx), symope(3,3))

      symope(1:3,1) = (/1d0,0d0,0d0/)
      symope(1:3,2) = (/0d0,1d0,0d0/)
      symope(1:3,3) = (/0d0,0d0,1d0/)
c      call rdpp_v2( ngpmx,ngcmx,nxx,  qibze,nqibze, qbze,nqbze,
c     i      nband, nl,ngrpx, nn,  nclass, nspin, symope,   qbas,
c     o      nblocha, lx, nx, ppbrd ,
c     o      mdimx, nbloch, cgr,
c     o      nblochpmx, ngpn,geigB,ngvecpB,  ngcni,ngveccB )
      call rdpp_v3(nxx, nl,ngrpx, nn, nclass, nspin,symope,qbas,
     o      nblocha, lx, nx,  ppbrd , mdimx, nbloch, cgr)
      if(.not.smbasis()) nblochpmx = nbloch + ngcmx

c      call rdpp_pln(ngpmx,ngcmx, qibze,nqibze, qbze,nqbze,nband,nspin,
c     o      ngpn,geigB,ngvecpB,ngcni,ngveccB)

c     &   geigB  (ngpmx,nband, nqbze,nspin), ngpn(nqbze),
c     &   ngvecpB(3,ngpmx, nqbze) , ngcni(nqibze), ! IBZ !
c     &   ngveccB(3,ngcmx,nqibze), !ngveccBr(3,ngcmx,nqibze),

c work arry --------------------------
      allocate(ngveccB(3,ngcmx))

c --- Check PARALELL.X0
      if(paralellx0) then
        if(iqxini >99999) stop ' hx0fp0: iqxini >99999'
        if(iqxend >99999) stop ' hx0fp0: iqxend >99999'
        ifvcfpout = iopen('VCCFP.'
     &     //xxt(iqxini,iqxend)
     &  ,0,-1,0)
      else
        iqxini= 2
        iqxend = nqibz + nq0i
        ifvcfpout = iopen('VCCFP',0,-1,0)
      endif
      read(ifvcfpout) ndummy1, ndummy2

c ... check write  ! mdim =nblocha
c      do i=1,nqbze
c       write(6,"(' i ngpn  gsum=',2i5,2d15.6)") i,ngpn(i),
c     &           sum(geigB(1:ngpn(i),1:nband, i))
c      enddo
c      do i=1,nqibze
c        call readqg('QGcou',qibze(1:3,i),ginv,  quu,ngc, ngveccB)
c        write(6,"(' i ngcni nsum=',2i5,i8)") i,ngc, sum(ngveccB(1:3,1:ngc))
c      enddo

C --- Open files EV* and CPHI ; eigenvlaues and eigenfunctions
c      ifev(1)    = iopen('EVU', 0,0,mrece)
c      if(nspin == 2) then
c        ifev(2)    = iopen('EVD', 0,0,mrece)
c      endif

      call init_readeigen(ginv,nspin,nband,mrece) !initialization of readEigen
      ! EVU EVD are readin in init_readeigen
      call init_readeigen2(mrecb,nlmto,mrecg)
c      ifcphi  = iopen('CPHI',0,0,mrecb) ! Coeeficients of eigenfunction

c --- ecore ---
c      allocate(ecore(nctot,nspin)) !core energies
c      do  is = 1,nspin
c        if (nctot .gt. 0) then
c           call catch1 (w(iecore),is,nctot,2,ecore(:,is)) !core energies
c           print *,' ecore is=',is,ecore(:,is)
c        endif
c      enddo

c --- set realomega, imagomega tetra nw niw nwp ifgb0vec ---------------------
!         nwp, freq_r,  frhis(1:nwhis+1)
      if   ( ixc==1 ) then !old imagw = 2 case
        realomega =.true.
        imagomega =.true.
        stop 'hsfp0sc: ixc==1 is not implimented'
cccccccccccccccccccc faleev 21May02,  use only ixc=1,11 modes cccccccccc
c      elseif( ixc==2.or.ixc==3 ) then
c        realomega =.true.
c        imagomega =.false.
c        niw = 0
c        ifepscond = 2102
c        open (ifepscond,file='EPScond')
c        read (ifepscond,*) epsrng, dwry !epsrng dw in Ry
c        dw = dwry/2d0
c        close(ifepscond)
c        if(dw==0d0) then
c          nw = 1
c        else
c          nw = (epsrng/2d0 - 1d-10)/(dw/2d0) + 2 !epsrng/2d0 corresponds to in a.u.
c        endif
c        allocate(epsi(nw,neps))
c        if(paralellx0) then
c          ifgb0vec = iopen ( "Mix0vec."
c     &     //xxt(iqxini,iqxend),1,3,0)
cc          "//charnum5(iqxini)//'to'//charnum5(iqxend),1,3,0)
c        else
c          ifgb0vec = iopen ( "Mix0vec",1,3,0)
c        endif
c      elseif(ixc==4.or.ixc==5.or.ixc==6) then
c! ... These are test modes.
c! ixc=4 tetrahedren weight test. tetwt5.vs.tetwt5. Write tethis.chk
c! ixc=5 Spectrum function (Img part) along the Real axis with tetwt4
c! ixc=6 Spectrum function (Img part) along the Real axis with tetwt5. Histgram method.
c        realomega = .true.
c        imagomega = .false.
c        tetra     = .true.
c        niw = 0
c! ---  For tetwt5 ---  the tetrahedron weight for spectrum function (imaginary part)
c!   Histogram bins are specified by freq_r(1:nwp)
c!     nwp=nw+1; frhis(1)=0
c!     The 1st   bin  is     [frhis(1),  frhis(2)]   ...
c!     The last  bin  is     [frhis(nw), frhis(nwp)].
c
c! ... These parameters specifies a test histgram bins;Sergey's mesh just for test modes.
c        nw0 = 200    !100    800
c        dwh = 0.01d0 !0.02d0 0.0025d0 !in hartree
c! ...
c        call findemaxmin(ifev,nband,nqbz,nspin,emax,emin)
c        if (nctot .gt. 0) Emin = minval(ecore)
c        omg2max = (Emax-Emin)*.5+.2d0    !(in Hartree) covers  all relevant omega, +.2 for margin
c        omg1max = dwh*(nw0-1)
c        nwp = int(sqrt(omg2max*(2*nw0-1d0)/dwh-(nw0**2-3*nw0+1d0)))+1 ! + 1 for margin
c        nw  = nwp-1
c        write(*,*) Emax,Emin,nw0,nw  ! nwp is new max number in frequency array
c        write(6,'(a32,2i7,2d15.3)')'hx0fp1: nw0,nw,omg1max,omg2max='
c     &              , nw0,nw,  omg1max,omg2max
c        if (nw <= nw0) stop 'hx0fp0:ixc==[456] nw2 <= nw'
c        allocate(freq_r(nwp))
c        do iw=1,nwp  !This is a test mesh by Sergey.Faleev
c          if(iw<=nw0) then;  freq_r(iw)=dwh*(iw-1)
c          else;  freq_r(iw)=dwh*(iw**2+nw0**2-3*nw0+1)/(2*nw0-1d0)
c          endif
c        enddo !freq_r(iw) is linear for iw<=nw and quadratic for nw<iw<=nw2
c               !freq_r(iw) chosen in such a way that it is continues with
c!!! nw nwp=nw+1 freq_r(1:nwp) are used after here.
c        allocate(frhis(nwp))
c        frhis=freq_r(1:nwp)
c        nwhis=nw
cC ... ixc=11 Sergey Faleev's normal  mode
      elseif(ixc==11) then
        realomega = .true.
        imagomega = .true.
        tetra     = .true.
c        if(readgwinput() ) then
        call getkeyvalue("GWinput","omg_c",omg_c )
c        else
c          call readomgc(ifinin,omg_c)
c        endif
        write(*,*)'dw, omg_c= ',dw, omg_c
C ... frhis
c        call findemaxmin(ifev,nband,nqbz,nspin
        call findemaxmin(nband,qbze,nqbze,nspin
     o  ,emax,emin)
        if(bzcase()==2) then
          allocate(qbz2(3,nqbz))
          do iq=1,nqbz
            qbz2(:,iq)=qbz(:,iq)+dq_bzcase2
          enddo
          call findemaxmin(nband,qbz2,nqbz,nspin ,emax2,emin2)
          emax=max(emax,emax2)
          emin=min(emin,emin2)
          deallocate(qbz2)
        endif

        print *,' emin emax=',emin, emax
        if (nctot .gt. 0) Emin=minval(ecore(:,1:nspin))
        omg2max=(Emax-Emin)*.5d0+.2d0 !(in Hartree) covers  all relevant omega, +.2 for margin
        print *,' omega2max=',omg2max
!             if (.false.) then

ccccccccccccccccccccccccccccccccc  histogram divisions cccccccccccccccccc

        nw2=int(omg_c/dw*( sqrt(1.+2*omg2max/omg_c)-1. ) )+1+3 !+3 for margin
        allocate(freqr2(nw2))                                  !+1 b/c (iw-1)
        do iw=1,nw2
          freqr2(iw)=dw*(iw-1)+dw**2/2./omg_c*(iw-1)**2
        enddo !linear + quadratic term
!             endif
        if (nw2 < 2 ) stop "hx0fp0: nw2 < 2"
        if (dw*(nw-2) > freqr2(nw2-1)) stop
     &                           "hx0fp0: dw*(nw-2) > freqr2(nw2-1)"
        nwhis = nw2-1
        allocate(frhis(1:nwhis+1))
        frhis = freqr2(1:nwhis+1)
        print *,' we set frhis nwhis=',nwhis
C ... for frhis_m
        nw=nw2-1   ! nw+1 is how many points of real omega we use
        ! for dressed coulomb line W(iw=0:nw) iw=0 corresponds omg=0
        ! maximum nw=nw2-1 because nwhis=nw2-1
        do iw=3,nw2-1
          !nw is chosen from condition that frhis_m(nw-3)<dw*(nw_input-3) <frhis_m(nw-2).
          !Here frhis_m(iw)= (freqr2(iw)+freqr2(iw+1))/2d0
          !nw was constructed such that omg=dw*(nw-2)> all relevant frequensies needed
          ! for correlation Coulomb Wc(omg),
          ! and one more point omg=dw*(nw-1) needed for extrapolation.
          ! Now, frhis_m(nw-1)> all relevent frequensies for Wc(omg)
          ! and one more point omg=frhis_m(nw) needed for extropolation
          ! used in subroutine alagr3z in  sxcf.f.
          omg2 = (freqr2(iw-2)+freqr2(iw-1))/2d0
cccccccccccccccccc
          print *,'iw omg2=',iw,omg2,dw,nw_input
cccccccccccccccccc
          if (omg2 > dw*(nw_input-3)) then
            nw=iw
            exit       ! 'nw_input' is only used to get maximum frequency for
            !  dressed coulomb line
          endif
        enddo
        allocate(freq_r(0:nw))
        freq_r(0)=0d0
        do iw=1,nw
          freq_r(iw)=(frhis(iw)+frhis(iw+1))/2d0
        enddo
        nwp = nw+1

        delta=0d0 !only for writeemesh

C--------- Apr2007
c        write(6,"(' ixc=11: niw nw nwhis=',3i5)") niw,nw,nwhis
c        stop 'xxxxxxxxxxxxxxx'


      endif

C ... timereversal=F is implimented only for tetra=T and sergeyv=T
      npm=1
      nw_i=0
      if(.not.timereversal()) then
        write(6,"('TimeReversal off mode')")
        npm=2
        nw_i=-nw
        if(.not.tetra)      stop ' tetra=T for timereversal=off'
      endif

! ...
c      if (realomega.and.ixc==1) then
c        nwp=nw
ccnw        allocate( freq_r(nwp) )
c        allocate( freq_r(0:nwp-1) )
c        do iw = 0,nwp-1
c          freq_r(iw) = dw*iw ! freq_r
c        enddo
c      endif
      if(realomega.and.ixc==1) then !
        nwp = nw+1 !Feb2006 I want to unify the meaning of nwp as nwp=nw+1
        !With this change x0kf_v2z calculate one more freq_r along real axis.
        ! at freq_r(nw). --->It WAS not calculated. It wat up to freq_r(nw-1).
        ! But not tested well.
        allocate( freq_r(0:nw) )
        do iw = 0,nw
          freq_r(iw) = dw*iw ! freq_r
        enddo
      endif
c-------------------------------------

      if(iSigma_en==0) nwp=1

c        open(UNIT=3111,file='freq_r') !write number of frequency
c        !points nwp and frequensies in 'freq_r' file
c        write(3111,*)nwp
ccnw          iw0=0            !for ixc=11
cccnw          if(ixc==1) iw0=1 !for ixc==1
c        do iw=0,nwp-1
c        write(3111,*)freq_r(iw)
c        enddo
c        close(3111)
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
ccccccccccccccccccccccccccccccccccccccccc faleev 21May02
c      if(imagonly.and.ixc==11) then
c        realomega =.false.
c        imagomega =.true.
c      endif

C --- set freq_i
      if (imagomega) then
        print *,' freqimg: niw =',niw
        allocate( freq_i(niw) ,freqx(niw),wx(niw),expa(niw) )
        call freq01 (niw,ua,  ! Generate gaussian frequencies x between (0,1) and w=(1-x)/x      
     o                freqx,freq_i,wx,expa)
        deallocate(freqx,wx,expa)
      endif

C --- tetra init
      if(tetra) then
        allocate( !wgt(nband+nctot,nband,nqbz), !noccxvv(nw+niw),
     &    nbnb(nqbz,npm),nbnbtt(nqbz,npm)
     &   ,ekxx1(nband,nqbz),ekxx2(nband,nqbz)) !, !!! nband=nlmto
c     &      idtetf(0:3,6*n1*n2*n3),ib1bz((n1+1)*(n2+1)*(n3+1)),
c     &      qbzw(3,(n1+1)*(n2+1)*(n3+1)) )
c        call TETFBZF(qbas, N1,N2,N3,qbz,nqbz,
c     o    idtetf,qbzw,ib1bz) ! make tetrahedron.
        if(tmpwwk)  ifwwk = iopen('tmp.wwk',0,-1,0)
      endif
      if(debug) print *,' xxx1:'
c ... maximum no. occupied valence states
      noccxv     = maxocc2 (nspin,ef, nband, qbze,nqbze)
      if(noccxv>nband) stop 'hx0fp0: all the bands filled! too large Ef'
      noccx      = noccxv + nctot

c --- write energy mesh ---
      if(paralellx0) then
        ifemesh = iopen('emesh.hx0fp0.'//xxt(iqxini,iqxend),1,-1,0)
      else
        ifemesh = iopen('emesh.hx0fp0',1,-1,0)
      endif
c      call writeemesh(ifemesh,freq_i,niw,freq_r,nwp,delta)

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
     o  nnmx,nomx)
        nblochpmx = nnmx
      endif
      mrecl  = nprecx*2*nblochpmx*nblochpmx/nword()

C --- open WV* files
      if(paralellx0) then
        ifwd   = iopen('WV.d.'//xxt(iqxini,iqxend),1,-1,0)
      else
        ifwd   = iopen('WV.d',1,-1,0)
      endif
c july 2003!
c      write (ifwd,*) !"(1x,i3,i8,i5,5i4)")
c     &   nprecx,mrecl,nblochpmx,nw_input,niw,nqibz + nq0i-1
      write (ifwd,"(1x,10i14)") !"(1x,i3,i8,i5,5i4)") 
     &   nprecx,mrecl,nblochpmx,nwp,niw,nqibz + nq0i-1,nw_i

cccccccccccccccccccccccccccc
c      stop 'test end 11111111111111111'
ccccccccccccccccccccccccccccccc
      allocate(  zw(nblochpmx,nblochpmx) )
      if (ixc==1.or.ixc==11) then
        if(paralellx0) then
c           ifrcw  = iopen('WVR.'
c     &     //xxt(iqxini,iqxend),0,-1,mrecl)
c
c           ifrcwi = iopen('WVI.'
c     &     //xxt(iqxini,iqxend) ,0,-1,mrecl)
        else
ctakao
          ifrcw  = iopen('WVR',0,-1,mrecl)
          ifrcwi = iopen('WVI',0,-1,mrecl)
c           ifrcw  = iopen('WVR',0,-1,mrecl/2)
c           ifrcwi = iopen('WVI',0,-1,mrecl/2)
          !mrecl/2 because  we use complex(4) format in WVR WVI!
          ! faleev 15 Okt2002
        endif
      endif

C --- nspinmx
      nspinmx = nspin

c --- READIN ANFcond------------------------------------
      inquire(file='ANFcond',EXIST=anfexist)
      print *,' ANFcond exist? =',anfexist
      if(anfexist) then
        allocate(ibasf(natom),transaf(1:3,natom))
        ifianf = 211
        open(ifianf,file='ANFcond')
        read(ifianf,*) soflag
        if(soflag=='SO') then !SO case is not implimented
          print *,' ANFcond found. but SO case is not implimented. '
     &           ,'Not using ANF symm'
          anfexist=.false.
          goto 889
        endif
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
      endif
 889  continue
C... oct2005
      if(smbasis() .and. anfexist) then
        print *,' smbais=T & anfexist=T is not yet. '//
     &        ' At leaset anfx0k shold be changed'
        stop ' smbais=T & anfexist=T is not yet. '
      endif

C --- ppb ---
!      allocate(ecore(nctot,nspin)) !core energies
      do  is = 1,nspin
! <Phi(R,L,n) Phi(R,L',n') B(R,i)> for each spin
        if(debug) print *," goto ppbafp is=",is
        irot = 1
        call ppbafp_v2 (irot,ngrpx,is,nspin,
     i              il,in,im, nlnm,
     d                 nl,nn,nclass,nlnmx,
     i                 mdimx,lx,nx,nxx,  !Bloch wave    
     i                 cgr, nl-1,        !rotated CG
     i                 ppbrd,            !radial integrals
     o              ppb(1,is)) !,pdb(1,is),dpb(1,is),ddb(1,is) )
c       if (nctot .gt. 0) call catch1  (w(iecore),is,nctot,2,ekc(:,is)) !core energies
c        if (nctot .gt. 0) then
c           call catch1 (w(iecore),is,nctot,2,ecore(:,is)) !core energies
c           print *,' ecore is=',is,ecore(:,is)
c        endif
      enddo

      print *, ' end of ppbafp_v2'

C ---- these are used in x0k
      call getkeyvalue("GWinput","nbcutlow",nbcut, default=0 )
      call getkeyvalue("GWinput","nbcutlowto",nbcut2, default=0 )
      write(6,"(' nbcut nbcutlowto=',2i5)") nbcut,nbcut2

C --- Calculate x0(q,iw) and W ------------------------------------------
c============ loop over q point ========================================
c=======================================================================
      if(paralellx0) then
     &  !skip
      elseif(bzcase()==1) then
        iqxini = 2
      else
        iqxini = 1
      endif

cccccccccccccccccccccccccccccddddcccccccccccccccccccccccccccccccccc
      do iq = iqxini,iqxend ! q=(0,0,0) is omitted!
        q = qibze(:,iq)
        write(6,"(' q list =',i4,3f14.6)")iq, q
      enddo
c      if (testtr()) then
c      iqxini = nqibz + 1
c      iqxend = nqibz + nq0i
c      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



cccccccccccccccccccccccccccccccccccccccccc
      do 1001 iq = iqxini,iqxend ! q=(0,0,0) is omitted!
        call cputid (0)
        q = qibze(:,iq)
ccccccccccccccccccccccccccccccccccccccccccccccc
c        if (negative_testtr()) then
c           q=-q
c        endif
ccccccccccccccccccccccccccccccccccccccccccccccc
        if(iq<=nqibz) then
c           iqx = iq
          qq=q
        else
c           iqx = 1              ! corresponding q=0
          qq=0d0
        endif

        write(6,"(' qq vec =',3f14.6)") qq
        call readqg('QGcou', qq, ginv,  quu,ngc,ngveccB)
c       ngb = nbloch + ngc !ngb is readin from vcoul

c ---readin coulomb matrix
c        allocate(vcoul(ngb,ngb))
c        read(ifvcfpout) vcoul          ! coulomb matrix


ccccccccccccccccccccccccccccccccccccccccccccccc
c        if (negative_testtr()) then
c           q=-q
c        endif
ccccccccccccccccccccccccccccccccccccccccccccccc
        do
          read(ifvcfpout,end=1005) ngb       ! oct2005 !I added ngb to VCCFP
          allocate(vcoul(ngb,ngb))
          read(ifvcfpout,end=1005) vcoul,qx  ! coulomb matrix
          if(sum(abs(q-qx))<1d-10) goto 1006
        enddo
 1005   continue
        write(6,"('iq q=',i5,3f9.4,' qx(VCCFP)=',3f9.4)")iq,q,qx
        stop 'hx0fp0_sc: can not find q in VCCFP'
 1006   continue
ccccccccccccccccccccccccccccccccccccccccccccccc
c        if (negative_testtr()) then
c           q=-q
c        endif
ccccccccccccccccccccccccccccccccccccccccccccccc

        write(6,"(' q vec =',3f14.6)") q
c        stop 'testxxxxxxxxxxxxxx'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ctest
c        if(iq<5) then
c           print *, ' --- skip iq=',iq
c          deallocate(vcoul)
c          cycle
c        endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccc faleev 21May2002 no use of ixc=2,3
c        if(ixc==2.or.ixc==3) then      ! only calculate iq>nqibz for ixc=2 mode
c          if(iq>nqibz) then; if(wqt(iq-nqibz)==0d0) then
c            iqixc2 = iq- (nqibz+nq0ix)
c            filepsnolfc ='EPS'
c     &      //char(48+iqixc2/10)//char(48+mod(iqixc2,10))//'.nolfc.dat'
c            ifepsdatnolfc = iopen ( filepsnolfc,1,3,0)
c            write(ifepsdatnolfc,*)
c     &                 ' q(1:3)   w(Ry)   eps    epsi  --- NO LFC'
c            if(ixc==3) then
c              fileps = 'EPS'
c     &        //char(48+iqixc2/10)//char(48+mod(iqixc2,10))//'.dat'
c              ifepsdat = iopen ( fileps,1,3,0)
c              write(ifepsdat,*) ' q(1:3)   w(Ry)   eps  epsi'//
c     &                  ' --- LFC included. '
c            endif
c            goto 110
c          endif;  endif
c          deallocate(vcoul)
c          cycle
c  110     continue
c        endif

c      print *," xxx1:"

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      print *," trace check for vcou"
c      trwv = 0d0
c      do i = 1,nbloch
c        trwv = trwv + vcoul(i,i)
c      enddo
c      write(6,'(" trvcou1=",i6,2d22.14)') iq,trwv
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        allocate( zw0(ngb,ngb) )

c     if(ixc==6 ) allocate( rcxq(ngb,ngb,nwp) )    !sf 21May02

        if(ixc==11) then
          print *," xxx2: memsize 8*ngb*ngb*nwhis=", 8*ngb*ngb*nwhis,
     &  ' ngb nwhis=',ngb,nwhis
          allocate( rcxq(ngb,ngb,nwhis,npm) )
          rcxq=(0d0,0d0)
        else
          if(onceww(2))
     &  print *," xxx2:allocate zxq zxqi memsize 16*ngb*ngb*(nwp+niw)=",
     &  16*ngb*ngb*(1+nwp+niw),' ngb nwp niw=',ngb,nwp,niw
          allocate(
     &    zxq (ngb,ngb,nw_i:nw),   !,nwp) feb2006
     &    zxqi(ngb,ngb,niw))
          zxq=0d0; zxqi=0d0
        endif

c     &    zxq (ngb,ngb,nwblk),
c     &    zxqi(ngb,ngb,niwblk))
c      do 1002 iwb = 1,nwdiv   !loop over frequency block
c      print *," xxx3"

c============ loop over spin============================================
ccccccccccccccccccccccccccccccccccccccccccc
c      nspinmx=1 !test
ccccccccccccccccccccccccccccccccccccccccccc
c=======================================================================
        do 1003 is = 1,nspinmx
          write(6,"(' *** ',2i4,' out of nqibz+n0qi nsp=',2i4,' *** ')") 
     &    iq, is, nqibz + nq0i,nspin
          if(debug) print *,' niw nw=',niw,nw

c===========tetraini block tetra==.true.===============================1ini
c --- tetrahedron method --- preparatory part
c  nbnb = total number of weight.
c  n1b  = band index for occ.   1\ge n1b \ge nband+nctot.
c         Valence-core order(Core index follows valence index).
c  n2b  = band index for unocc. 1\ge n2b \ge nband
c  wwk(ibib,...)  = (complex)weight for the pair for n1b(ibib...),n2b(ibib...).
          if(tetra) then
            write(6,"(' tetra mode niw nw nwhis=',3i5)") niw,nw,nwhis
cccccccccccccccccccccccc
c           qqq=q
cccccccccccccccccccccccc
C ... kqxx(kx) ekxx ---
c    ekxx1  for rk
c    ekxx2 for q+rk See tetwt4
            do kx = 1, nqbz
c             call rwdd1 (ifev(is), kx, nband, ekxx1(1:nband, kx) )
              call readeval(qbz(:,kx)-qqq,  is, ekxx1(1:nband, kx) )
            enddo
            do kx = 1, nqbz
              call readeval(q+qbz(:,kx)-qqq,is, ekxx2(1:nband, kx) )
            enddo

c            if(iq<=nqibz) then
c              do kx = 1, nqbz
cc                call fbz2 ( q(1:3) + qbze(1:3,kx),
cc     i           ginv ,index_qbz,n_index_qbz,qbz,nqbz,
cc     o           qbzx,kqxx)
cc               ekxx2(1:nband, kx) = ekxx1(1:nband, kqxx) !c    write(1211,*) kx, kqxx
c              call readeval(q + qbz(:,kx),is, ekxx2(1:nband, kx) )
c              enddo
c            else
c              do kx = 1, nqbz;  kp  = nqbz *(iq - nqibz) + kx
cc              call rwdd1 (ifev(is), kp, nband, ekxx2(1:nband, kx) )
c               call readeval(q+qbz(:,kx),is, ekxx2(1:nband, kx) )
c              enddo
c            endif

ctakao-feb/2002 I replaced tetwt4 (1d30) with tetwt5(job=0) -----
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
c            allocate( demin(nband+nctot,nband,nqbz),
c     &                demax(nband+nctot,nband,nqbz) )
c            allocate(iwgt(nband+nctot,nband,nqbz) )
! wgt, demin, demax may require too much memory in epsilon mode.
! We will have to remove these memory allocations in future.
! tetwt5x_dtet2 can be very slow because of these poor memory allocation.
            if(nctot==0) then
              deallocate(ecore)
              allocate(ecore(1,2)) !this is dummry
            endif
            call tetwt5x_dtet4(npm,ncc,  
     i           q, ekxx1, ekxx2, qbas,ginv,ef, 
     d          ntetf,nqbzw, nband,nqbz,
     i          nctot,ecore(1,is),idtetf,qbzw,ib1bz, 
     i          job,
     o          iwgt,nbnb,      !job=0
     o          demin,demax,   !job=0
     i          frhis, nwhis,  ! job=1    not-used
     i          nbnbx,ibjb,nhwtot,  ! job=1    not-used
     i          ihw,nhw,jhw, ! job=1    not-used
     o          whw,         ! job=1    not-used
     i   iq,is,is,nqibz) ! new input for devided-tetrahedron method.
            if(ixc==1) deallocate(demin,demax) !sf 21May02
            nbnbx = maxval(nbnb(1:nqbz,1:npm)) !nbnbx = nbnbxx
            if(debug) print *,' nbnbx=',nbnbx
            allocate(  n1b(nbnbx,nqbz,npm)
     &                ,n2b(nbnbx,nqbz,npm))
            do jpm=1,npm
              call rsvwwk00_4(jpm, iwgt(1,1,1,jpm),nqbz,nband,nctot,ncc, nbnbx,
     o        n1b(1,1,jpm), n2b(1,1,jpm), noccxvx(jpm), nbnbtt(1,jpm))
            enddo

c ... check
            if(debug) then
              do kx  = 1, nqbz
                do jpm = 1, npm
c        do ix  = 1, nbnb(kx,jpm)
c        write(6,"('jpm kx ixx n1b n2b=',5i7)")jpm,kx,ix,n1b(ix,kx,jpm),n2b(ix,kx,jpm)
c        enddo
                  write(6,"('jpm kx  minval n1b n2b=',4i5)")jpm,kx,
     &        minval(n1b(1:nbnb(kx,jpm),kx,jpm)),
     &        minval(n2b(1:nbnb(kx,jpm),kx,jpm))
                enddo
              enddo
            endif
c
            if(sum(abs(nbnb-nbnbtt))/=0)then
              do ik=1,nqbz
                write(6,*)
                write(6,*)"nbnb  =",nbnb(ik,:)
                write(6,*)"nbnbtt=",nbnbtt(ik,:)
              enddo
              stop'hx0fp0:sum(nbnb-nbnbtt)/=0'
            endif
            noccxv = maxval(noccxvx)
            noccx  = nctot + noccxv
            print *,' Tetra mode: nctot noccxv= ',nctot,noccxv
            deallocate(iwgt)
          endif
c=========end of tetraini block==========================================1end



c=======TetrahedronWeight_4 block. tettwt4  ixc==1 ==============2ini
          if(tetra.and.ixc==1) then
            stop 'hx0fp0sc: ixc==1 tetra is not impliemented'
          endif !ixc==1



c=======TetrahedronWeight_5 block. tetwt5  ixc==,4,6,11 =======4ini
          if(ixc==11) then !sf 21May02
C --- METHOD (tetwt5) for the tetrahedron weight
!   Histogram secstions are specified by frhis(1:nwp)
!     The 1st   bin  is     [frhis(1),  frhis(2)]   ...
!     The last  bin  is     [frhis(nw), frhis(nwp)].
!     nwp=nw+1; frhis(1)=0
! takao-feb/2002
            if(frhis(1)/=0d0) stop ' hx0fp0: we assume frhis(1)=0d0'
            write(6,*)' --- Frequency  histogram  sections ---- '
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
c              ib1=  n1b(ibib,ik)
c              ib2 = n2b(ibib,ik)
cccccccccccc
c              print *,' goto histrange=',ib1,ib2,jpm
c     &         ,demin(n1b(ibib,ik,jpm),n2b(ibib,ik,jpm),ik,jpm)
c     &         ,demax(n1b(ibib,ik,jpm),n2b(ibib,ik,jpm),ik,jpm)
cccccccccccc
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
     i   iq,is,is,nqibz) ! new input for devided-tetrahedron method.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c            do ik   = 1,nqbz
c            do ibib = 1,nbnb(ik)
c              ib1  = n1b(ibib,ik)
c              ib2  = n2b(ibib,ik)
c              write(3666,"('--- ik ib jb =',3i4)") ik,ib1,ib2
c              ihww = ihw(ibib,ik)
c              do iw = 1, nhw(ibib,ik)
c               write(3666,"('ihis [a b] whw =',i4,2f8.5,f18.12,d14.6)")
c     &         ihww+iw-1, frhis(ihww+iw-1),frhis(ihww+iw),
c     &         whw(jhw(ibib,ik)+iw-1), whw(jhw(ibib,ik)+iw-1)
c              enddo
c            enddo
c            enddo
c            stop 'hx0fp0.f ---  testend xxx2'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            ! endif                        !ixc==11   !sf 21May02
c=======TetrahedronWeight_5 block end ==================4end




c==============x0kf_v2h block. ixc==6,11 =========================5ini
c Instead of wwk, you have to use whw(nwhtot).
c You will need to pass these arrays to x0kf
c   ihw(ibjb,kx): omega index, to specify the section of the histogram.
c   nhw(ibjb,kx): the number of histogram sections
c   jhw(ibjb,kx): pointer to whw
c   whw( jhw(ibjb,kx) ) \to whw( jhw(ibjb,kx) + nhw(ibjb),kx)-1 ), where ibjb=ibjb(ib,jb,kx)
c     : histogram weights for given ib,jb,kx for histogram sections
c       from ihw(ibjb,kx) to ihw(ibjb,kx)+nhw(ibjb,kx)-1.
c------------------------------------------------------------------
            !if(ixc==6.or.ixc==11) then    !ixc==11   !sf 21May02
            print *, ' goto x0kf_v3hz'


ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c          goto 8888
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        do k=1,nqbz
c        n1bmax  = 0; n2bmin = nband + nctot +1000000  ! for sanity check
c        do ibib = 1, nbnb(k)
c          if( n1b(ibib,k)<=nband) n1bmax = max(n1b(ibib,k),n1bmax)
c          n2bmin = min(n2b(ibib,k),n2bmin)
c        enddo
c          print *,' n2bmin=',n2bmin
c        enddo
c        stop 'test end'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c- new version only along real axis for
            call x0kf_v3hz( npm,ncc,
     i          ihw,nhw,jhw,whw,nhwtot,  ! tetwt5
     i            n1b,n2b,nbnbx,nbnb,    ! use whw by tetwt5 ,
     i              q,  
     i              nspin,is,is, 
     i              qbas,ginv,  qbz,wbz, 
     i              nblocha ,nlnm,nlnmv,nlnmc,iclass, 
     i              ppb(1,is), 
     i              icore,ncore,
     d            nlmto,nqbz,noccx,noccxv,
     d            nl,nclass,natom,nnc,
     d            nlnmx,mdimx,nbloch,  nwhis,
     i    iq,ngb,ngc,ngpmx,ngcmx,  !ngb/=ngc+nbloch for smbasis()=T oct2005
     i    nqbze,nband,nqibz, 
     o    rcxq, ! See rcxq_zxq in x0kf_v3h.f; rcxq is equivalent with zxq;
            !  rcxq is the accumulating variable for spins 
     i    nolfco,zzr,  !for nolfco
     o    rcxqmean   ) !for nolfco

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 8888     continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

          endif !ixc==11
c==============x0kf_v2h block end  ===========================5end

c=========== Check write mode4 block ixc==4 ============6ini
cccccccccccccccccc faleev 21May02 ixc==4 not used
c        if(ixc==4) then
cC --- check write ---
c          ichkhis = iopen ('tethis.chk',1,-1,0)
c          do ik   = 1,nqbz
c          do ibib = 1,nbnb(ik)
c              ib1  = n1b(ibib,ik)
c              ib2  = n2b(ibib,ik)
c              ihww = ihw(ibib,ik)
c              write(ichkhis,"('--- ik ib jb =',3i4)") ik,ib1,ib2
c              do iw = 1, nhw(ibib,ik)
c                write(ichkhis, "('his [min max]=',i6,2f7.4,
c     &             ' whw wwk*h= ',2d12.5,'( =',d12.5,'+',d12.5,')')")
c     &          ihww+iw-1, freq_r(ihww+iw-1), freq_r(ihww+iw),
c!-------------------------------------------------------
c     &          whw(jhw(ibib,ik)+iw-1),                     !weight for the histgram range. by tetwt5
c!------------ v.s. -------------
c     &          dimag(wwk(ibib,ik,ihww+iw-1)+wwk(ibib,ik,ihww+iw))/2d0
c     &           * (freq_r(ihww+iw-1)- freq_r(ihww+iw)),    !weight=weight1+weight2 for each bin of the histgram range by tetwt4.fal
c!-------------------------------------------------------
cc
c     &          dimag(wwk(ibib,ik,ihww+iw-1))/2d0
c     &           * (freq_r(ihww+iw-1)- freq_r(ihww+iw)),    !weight1 for the histgram range by tetwt4
c     &          dimag(wwk(ibib,ik,ihww+iw))/2d0
c     &           * (freq_r(ihww+iw-1)- freq_r(ihww+iw))     !weight2 for the histgram range by tetwt4
c              enddo
c          enddo
c          enddo
c          ichkhis = iclose('tethis.chk')
c          print *,' OK! hx0fp1_tetwt5_test end. Histgram method.'
c          print *,' The final reulst wwk vs.whw is in tethis.chk!'
c          print *,' This code is a routine to test tetwt5.'
c          stop ' ixc==4  Weight Histgram check. wwk vs.whw'
c         endif  ! end of if(tetra)
c=========== Check write mode4 block end  ==============6end

          if(ixc==11) 
     &       deallocate(ihw,nhw,jhw, whw,ibjb )
          if(tetra.and.ixc==1) deallocate(wwk)
          if(tetra)            deallocate( n1b,n2b)

 1003   continue;print *, 'end of spin-loop nwp=',nwp  !end of spin-loop
c===========end of spin loop============================================
c=======================================================================


c--- Hilbert transform.  Genrerate Real part from Imaginary part. ======
        if(ixc==11) then
          write(6,'(" xxx4: allocate zxq zxqi"//
     &    " =16*ngb*ngb*(nwp+niw)=",i8," ngb nwp niw=",3i7)')
     &    16*ngb*ngb*(1+nwp+niw), ngb, nwp, niw
          allocate(
     &      zxq (ngb,ngb,nw_i:nw),  !,nwp), 
     &      zxqi(ngb,ngb,niw))
          zxq =0d0
          zxqi=0d0
          nw_w=nw
          if(iSigma_en==0) nw_w=0 !i think this is ok , but not tested... feb2006
          write(6,*)"goto dpsion4 nwhis nw_i nw_w=",nwhis,nw_i,nw
          call dpsion4(frhis,nwhis, freq_r, nw_w, freq_i,niw, realomega, imagomega, 
     i        ngb, rcxq, npm,nw_i,  ! rcxq is alterd---used as work
     o        zxq, zxqi, 
     i        nolfco,.false., 0d0,is, rcxqmean,1d50,1d50, 
     o        x0mean)  !zxq x0mean is accumlating veriable for chipm mode
c          call dpsion3(frhis,nwhis, freq_r, nw_w, freq_i,niw,
c     i                             realomega,    imagomega,
c     i    ngb, rcxq,
c     o    zxq, zxqi, ! zxq... are Not the accumulating variable for spins
c     i    nolfco,.false.,0d0,is, rcxqmean,   1d50,1d50,
c     o    x0mean)
cc         do iw=1,nwp
cc         write(6,"(a,i5,2d13.5)")' zzz=',iw,sum(zxq(:,5,iw))
cc         enddo
        endif

c ... convert from rcxq to zxq for ixc==6
c        if(ixc==6) then
c          do iw  = 1,nw
c          call rcxq_zxq(rcxq(1,1,iw),zxq(1,1,iw),ngb,1)
c          enddo
c        endif

c ---  RealOmega ====================================
        if (realomega) then
          if (nspin == 1) zxq = 2d0*zxq   !if paramagnetic, multiply x0 by 2
          nwmax = nw
          nwmin = nw_i
c          if(ixc==11) nwmax=nwp
          if(iSigma_en==0) then
            nwmin=0
            nwmax=0
          endif
          if (anfexist) then   !Antiferro accelaration.
            do iw  = nwmin,nwmax
              call anfx0k(natom,nclass,nblocha,iclass,pos
     i     ,nbloch,ngc,
     i      q, ngveccB,qbas, ! for q+G
     i      anfvec,ibasf,             ! these are antiferro informations.
     i      zxq(1,1,iw))              ! i/o
            enddo
          endif
c          write (6,"(' end of anfx0k... q='3f10.4)") q(1),q(2),q(3)
          write(6, *)" trace check for W-V nwmin nwmax=",nwmin,nwmax
          do iw  = nwmin,nwmax
            frr= dsign(freq_r(abs(iw)),dble(iw))
            imode = 1
            if(debug) write(6,"(' rrr sum zxq=',i5,4d23.14)") iw,sum(zxq(:,:,iw)),sum(abs(zxq(:,:,iw)))
            call wcf( ngb, vcoul, zxq(1,1,iw), imode,
     &                zw0)
            if(debug) write(6,"(' rrr sum zw0=',i5,4d23.14)") iw,sum(zw0),sum(abs(zw0))
            zw(1:ngb,1:ngb) = zw0 ! zw(nblochpmx,nblochpmx)
            write(ifrcw, rec=((iq-iqxini)*(nwmax-nwmin+1)+ iw-nwmin+1)) zw   !  WP = vsc-v
            !----------------------------
            call tr_chkwrite("freq_r iq iw realomg trwv=", 
     &               zw, iw, frr,nblochpmx, nbloch,ngb,iq)
          enddo !iw
        endif !realomega
c --- RealOmega end ==============================




C --- ImagOmega ===================================
        if (imagomega) then
          print *,' goto imag omega'
          if (nspin == 1) zxqi = 2d0*zxqi   ! if paramagnetic, multiply x0 by 2
c ... Antiferro accelaration.
          if (anfexist) then
            do iw  = 1,niw
              call anfx0k(natom,nclass,nblocha,iclass,pos
     i      ,nbloch,ngc,
c     i      q, ngveccB(1,1,iqx),qbas,    ! for q+G 12jan2004
     i      q, ngveccB,qbas,    ! for q+G
     i      anfvec,ibasf,  ! these are antiferro informations.
     i      zxqi(1,1,iw))  ! i/o
            enddo
          endif
          imode=1
          do iw  = 1,niw
            print *,'sss sumchk zxqi=',sum(abs(zxqi(:,:,iw)))
            call wcf( ngb, vcoul,zxqi(1,1,iw),imode, 
     &              zw0)
            print *,'sss sumchk zw0=',sum(abs(zw0))
            zw(1:ngb,1:ngb) = zw0 ! zw(nblochpmx,nblochpmx)
            write(ifrcwi, rec=((iq-iqxini)*niw+iw) ) zw    !  WP = vsc-v
            call tr_chkwrite("freq_i iq iw imgomg  trwv=",zw,iw,freq_i(iw),nblochpmx,nbloch,ngb,iq)
          enddo
        endif !imagomega
C --- ImagOmega end ===============================

c 1002 continue  ! end of frequency block-loop
        deallocate(vcoul, zw0, zxq, zxqi)
        if(ixc==11) deallocate(rcxq)
 1001 continue
c============ end of loop over q point =================================
c=======================================================================

c close files
      is = iclose('hbe.d')
c      is = iclose('EVU')
c      if(nspin == 2) is = iclose('EVD')
      call cputid(0)
      if(ixc==1) stop ' OK! hx0fp0_sc mode=1'
c      if(ixc==2) stop ' OK! hx0fp0 eps without LFC mode=2'
c      if(ixc==3) stop ' OK! hx0fp0 eps with LFC mode=3'
c      if(ixc==4) stop ' OK! hx0fp0 ixc=4'
c      if(ixc==5) stop ' OK! hx0fp0 ixc=5'
c      if(ixc==6) stop ' OK! hx0fp0 ixc=6'
      if(ixc==11) stop ' OK! hx0fp0_sc ixc=11 Sergey F. mode'
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
c      do i = 1,ngb
c        write(6,'("iii i=",i4,a,f10.4,2i5,4d22.14)')i,tagname,freqq,iq,iw,zw(i,i)
c      enddo
      end
