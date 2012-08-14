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
      use m_pbindex,only: PBindex,norbt,l_tbl,k_tbl,ibas_tbl,offset_tbl,offset_rev_tbl
c      use m_spacegrouprot,only: Spacegrouprot
      use m_readqgcou,only: readqgcou
c      use m_rotMPB,only: rotMPB

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
     & ,rcxqmean(:,:,:,:) !takao sep2006

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
c      real(8),allocatable:: ppovl(:,:)
      complex(8),allocatable:: zzr(:,:),ppovl(:,:),ppovlz(:,:),ppovlzinv(:,:)
      complex(8) :: epxxx,vcmean
      character*9 fileps
      character*15 filepsnolfc
      logical :: paralellx0=.true. !, hist
c      character(5) :: charnum5
c      character(8) :: xt
      character(20):: xxt

c faleev
      real(8) :: Emin, Emax      ,emax2,emin2
c      integer(4) :: iSigma_en  !sf..21May02  !iSigma_en is integer
      !parameter stored in GWIN_V2
      !which determines approximation for  self-energy.
      !Self-energy should be made hermitian for energies to be real
cxxx      !iSigma_en==0 SE_nn'(ef)+img integral:delta_nn'([SE_nn(e_n)+c.c.]/2-SE_nn(ef))
cxxx      !iSigma_en==1 SE_nn'(ef)+delta_nn'([SE_nn(e_n)+c.c.]/2-SE_nn(ef))
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
      complex(8) :: x0meann !dummy

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

!!
      logical :: newaniso,newaniso2 !apr2012
      integer :: ngb0,ifvcoud,idummy,ifepstinv,igb1,igb2,ngb_in,nmbas1,nmbas2,iq0,ifisk,iqx,ig
      complex(8),allocatable:: zcousq(:,:),epstinv(:,:),epstilde(:,:),zcousqrsum(:,:,:),zcousqr(:,:)
      real(8),allocatable:: vcousq(:)
      real(8):: fourpi,sqfourpi,tpioa,absq,vcou1,vcou1sq

!! Eq.(40) in PRB81 125102
c      complex(8),allocatable::sk(:,:,:),sks(:,:,:),skI(:,:,:),sksI(:,:,:),
c     &  w_k(:,:,:),w_ks(:,:,:),w_kI(:,:,:),w_ksI(:,:,:), llw(:,:), llwI(:,:),
      complex(8),allocatable::sk(:),sks(:),skI(:),sksI(:), 
     &  w_k(:),w_ks(:),w_kI(:), w_ksI(:), s_vc(:),vw_k(:),vw_ks(:)
      complex(8),allocatable:: llw(:,:), llwI(:,:),w0(:),w0i(:),aaamat(:,:)
      real(8),allocatable:: dmlx(:,:),epinvq0i(:,:),epinv(:,:,:),epinvq0i_m1(:,:),wklm(:),qeibz(:,:,:)
      integer:: lxklm,nlxklm,ifidmlx,ifrcwx,iq0xx,ircw,nini,nend,iwxx,nw_ixxx,nwxxx,niwxxx,iwx,icc1,icc2
      complex(8):: vc1vc2
      integer,allocatable:: neibz(:),nwgt(:,:),ngrpt(:),igx(:,:,:),igxt(:,:,:),eibzsym(:,:,:)

      real(8),allocatable:: aik(:,:,:,:)
      integer,allocatable:: aiktimer(:,:)
      integer:: l2nl, nmbas_in,schi=-9999!dummy
      logical:: eibz4x0,tiii,iprintx,chipm=.false.
      real(8)::qvv(3),ecut,ecuts,hartree,rydberg,pi
      
c      integer(4),allocatable :: iclasst(:)
!-------------------------------------------------------------------------
      hartree=2d0*rydberg()
      pi         = 4d0*datan(1d0)
      fourpi=4d0*pi
      sqfourpi=sqrt(fourpi)
      call cputid (0)
C --- computational mode select
      write(6,*) ' --- hx0fp0_sc Choose omodes below ----------------'
      write(6,*) ' Only for ixc==11 '
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
c      if    (ixc==1)  then; write(6,*) ' OK ixc=1. normal mode'
      if(ixc==11) then; write(6,*) " OK ixc=11 s.f.'s normal mode "
      else
         print *,'we only allow ixc==11. given ixc=',ixc
         stop 'error:we only allow ixc==11.'
      endif

!! ====newaniso2====
      newaniso2=.false.
      if(newaniso()) then !.and.ixc==11) then
        newaniso2=.true.
      endif
!!
      call getkeyvalue("GWinput","ecut_p" ,ecut, default=1d10 )
      call getkeyvalue("GWinput","ecuts_p",ecuts,default=1d10 )

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


c      call readd_iSigma_en(ifinin,iSigma_en) !sf..21May02
c ...
      nw_input = nw ; print *, 'nw delta=',nw_input,delta
c      if (nclass > mxclass) stop ' hx0fp0: increase mxclass'

      debug=.false.
      if(verbose()>=100) debug=.true.
      if(debug) print *,' end of genallc'

      tpioa=2d0*pi/alat

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
      if(.not.newaniso2) then
        wqtsum = sum(abs(wqt(1:nq0i)))
c       if(ixc==1.and.abs(wqtsum-1d0) >1d-10) stop ' wqtsum of Q0P /=1 '
        call getkeyvalue("GWinput","TestNoQ0P",noq0p,default=.false.)
c        if(ixc==1.and.abs(wqtsum-1d0) >1d-10.and.(.not.noq0p)) 
c     &  stop ' wqtsum of Q0P /=1 '
      endif 

c --- readin by rdpp ; Radial integrals ppbrd and plane wave part
c      call getsrdpp( nclass,nl,
c     o               ngpmx,ngcmx,nxx )
      call getsrdpp2( nclass,nl,nxx)
      call readngmx('QGpsi',ngpmx)
      call readngmx('QGcou',ngcmx)
      print *,' ngcmx ngpmx=',ngcmx,ngpmx
! qibze(3,nqbze) qbze(3,nqibze)
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
      l2nl=2*(nl-1)
      allocate( nx(0:l2nl,nclass),  nblocha(nclass) ,lx(nclass),
     &   ppbrd ( 0:nl-1, nn, 0:nl-1,nn, 0:l2nl,nxx, nspin*nclass),
     &   cgr(nl**2,nl**2,(2*nl-1)**2,ngrpx), symope(3,3))
      symope(1:3,1) = (/1d0,0d0,0d0/)
      symope(1:3,2) = (/0d0,1d0,0d0/)
      symope(1:3,3) = (/0d0,0d0,1d0/)
      call rdpp_v3(nxx, nl,ngrpx, nn, nclass, nspin,symope,qbas,
     o      nblocha, lx, nx,  ppbrd , mdimx, nbloch, cgr)
      if(.not.smbasis()) nblochpmx = nbloch + ngcmx
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
c      if   ( ixc==1 ) then !old imagw = 2 case
c        realomega =.true.
c        imagomega =.true.
c        stop 'hsfp0sc: ixc==1 is not implimented'
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
c      elseif(ixc==11) then
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

c      endif

C ... timereversal=F is implimented only for tetra=T and sergeyv=T
      npm=1
      nw_i=0
      if(.not.timereversal()) then
        write(6,"('TimeReversal off mode')")
        npm=2
        nw_i=-nw
        if(.not.tetra)      stop ' tetra=T for timereversal=off'
      endif
      print *,'Timereversal=',Timereversal()


! ...
c      if (realomega.and.ixc==1) then
c        nwp=nw
ccnw        allocate( freq_r(nwp) )
c        allocate( freq_r(0:nwp-1) )
c        do iw = 0,nwp-1
c          freq_r(iw) = dw*iw ! freq_r
c        enddo
c      endif

c      if(realomega.and.ixc==1) then !
c        nwp = nw+1 !Feb2006 I want to unify the meaning of nwp as nwp=nw+1
c        !With this change x0kf_v2z calculate one more freq_r along real axis.
c        ! at freq_r(nw). --->It WAS not calculated. It wat up to freq_r(nw-1).
c        ! But not tested well.
c        allocate( freq_r(0:nw) )
c        do iw = 0,nw
c          freq_r(iw) = dw*iw ! freq_r
c        enddo
c      endif
c-------------------------------------

c      if(iSigma_en==0) nwp=1

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
      if(newaniso2) then
         iqxini=1
      elseif(paralellx0) then
     &  !skip
      elseif(bzcase()==1) then
        iqxini = 2
      else
        iqxini = 1
      endif

cccccccccccccccccccccccccccccddddcccccccccccccccccccccccccccccccccc
c      do iq = iqxini,iqxend ! q=(0,0,0) is omitted!
c        q = qibze(:,iq)
c        write(6,"(' q list =',i4,3f14.6)")iq, q
c      enddo
c      if (testtr()) then
c      iqxini = nqibz + 1
c      iqxend = nqibz + nq0i
c      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!! === Use of symmetry. EIBZ procedure PRB81,125102 ===
!!  For rotation of zcousq.  See readeigen.F rotwv.F ppbafp.fal.F(for index of product basis).
      if(eibz4x0()) then
        allocate( nwgt(nqbz,iqxini:nqibz), !qeibz(3,nqbz,iqxini:nqibz),neibz(iqxini:nqibz),
     &    igx(ngrp*2,nqbz,iqxini:nqibz),igxt(ngrp*2,nqbz,iqxini:nqibz),
     &    eibzsym(ngrp,-1:1,iqxini:nqibz))
         iprintx=.false.
         write(6,"('TimeRevesal switch = ',l1)") timereversal()
         call eibzgen(nqibz,symgg,ngrp,qibze,iqxini,nqibz,qbz,nqbz,timereversal(),ginv,iprintx,
     o   nwgt,igx,igxt,eibzsym)
!! Check timereversal is required for symmetrization operation or not. If not tiii=timereversal=F is used.
!! this is because the symmetrization is a little time-consuming.
        tiii=timereversal() 
        if(minval(igxt)==1) tiii=.false.
        iprintx=.true.
ccccccccccccccccccccccccc
c        tiii=.true.
cccccccccccccccccccc
        write(6,"('=== goto eibzgen === used timereversal=',l1)")tiii
        call eibzgen(nqibz,symgg,ngrp,qibze,iqxini,nqibz,qbz,nqbz,tiii,ginv,iprintx,
     o    nwgt,igx,igxt,eibzsym)
        call PBindex(natom,lx,l2nl,nx) !all input. this returns requied index stored in arrays in m_pbindex.
            ! PBindex: index for product basis.  We will unify this system; still similar is used in ppbafp_v2.
        call readqgcou() !no input. Read QGcou and store date into variables.
c       call Spacegrouprot(symgg,ngrp,plat,natom,pos) ! all inputs.
      endif
c      stop 'xxxxxxxxxxxxxxxxxx'

      allocate( llw(nw_i:nw,iqxini:iqxend), llwI(niw,iqxini:iqxend) )
      ifisk=iopen('WKcomponent',0,-1,0)


cccccccccccccccccccccccccccccccccccccccccc
      do 1001 iq = iqxini,iqxend ! q=(0,0,0) is omitted!
        call cputid (0)
        q = qibze(:,iq)
        if(iq<=nqibz) then
          qq=q
        else  ! iqx = 1    corresponding q=0
          qq=0d0  
        endif
        write(6,"(' qq vec =',3f14.6)") qq
        call readqg('QGcou', qq, ginv,  quu,ngc,ngveccB)
        if( newaniso2.and.iq==1 ) then ! *sanity check
          if(sum(q**2)>1d-10) then
            stop ' hx0fp0.sc: sanity check. |q(iqx)| /= 0'
          endif
        endif

!! ==== readin Coulomb matrix ====
        if(newaniso2) then
           allocate(vcoul(1,1))
           ngb=nbloch+ngc
        else
          do
          read(ifvcfpout,end=1005) ngb       ! oct2005 !I added ngb to VCCFP
          allocate(vcoul(ngb,ngb))
          read(ifvcfpout,end=1005) vcoul,qx  ! coulomb matrix
          if(sum(abs(q-qx))<1d-10) goto 1006
          enddo
 1005     continue
          write(6,"('iq q=',i5,3f9.4,' qx(VCCFP)=',3f9.4)")iq,q,qx
          stop 'hx0fp0_sc: can not find q in VCCFP'
 1006     continue
        endif  

        write(6,"(' q vec =',3f14.6)") q

!! === readin diagonalized Coulomb interaction ===
!! zcousq: E(\nu,I), given in PRB81,125102; vcousq: sqrt(v), as well.

        if(newaniso2) then
!! qq is used.
          if(sum(abs(qq))<1d-6) then !reopen Vcoud
            idummy  = iclose('Vcoud') !close and open again. This is because first data in Voud is for q=0
            ifvcoud = iopen('Vcoud',0,0,0)
          endif
          read(ifvcoud) ngb0
          read(ifvcoud) qvv
c         print *,'qvv=',qvv
          if(sum(abs(qvv-qq))>1d-10) then
            print *,'qvv =',qvv
            stop 'hx0fp0: qvv/=0 hvcc is not consistent'
          endif
          if(allocated(zcousq)) deallocate( zcousq,vcousq )
          allocate( zcousq(ngb0,ngb0),vcousq(ngb0))
          read(ifvcoud) vcousq
c         print *,'qvv =',qvv
c         do igb1=1,ngb0
c           print *,'qvv sum vcousq ',igb1,vcousq(igb1)
c         enddo
          read(ifvcoud) zcousq
          vcousq=sqrt(vcousq)
        endif
!!
        if(newaniso2.and. iq>nqibz  ) then !for mmmw
cc        if(.false.) then !comment out for mmmz
          if( ngb0/=ngb ) then
            stop 'hx0fp0.m.f:ngb0/=ngb'
          endif
          nolfco =.true.
          nmbas_in=ngb
          nmbas1 = 2  !nmbas1=2 is for newaniso modes. For P1mu and Pmu1 components.
          nmbas2 = ngb
          if(allocated(rcxqmean)) deallocate(rcxqmean)
c          if(allocated(x0mean)  ) deallocate(x0mean)
          if(allocated(zzr))      deallocate(zzr)
          allocate( rcxqmean(nwhis,npm,nmbas1,nmbas2) )
c          allocate(   x0mean(nw_i:nw,  nmbas1,nmbas2) )
          allocate( zzr(ngb,ngb) )
          rcxqmean= 0d0
        elseif(newaniso2) then !.and.iq==1) then
          nmbas1=ngb !meaningful in dpsion5, not used in x0kf_v4h.
          nmbas2=ngb !
        endif

        if(newaniso2) then
          if(allocated(ppovlz)) deallocate(ppovlz)
          if(allocated(ppovlzinv)) deallocate(ppovlzinv)
          if(allocated(ppovl)) deallocate(ppovl)
          allocate(ppovl(ngc,ngc),ppovlz(ngb,ngb),   ppovlzinv(ngb,ngb))
          call readppovl0(qq,ngc,ppovl)
          ppovlz(1:nbloch,:) = zcousq(1:nbloch,:)
          ppovlz(nbloch+1:nbloch+ngc,:) 
     &      = matmul(ppovl,zcousq(nbloch+1:nbloch+ngc,:))
          if(allocated(zzr)) zzr = ppovlz
c          deallocate(ppovl)
          write(*,*)'nnnnn',nbloch+ngc,ngb
        endif
  

c$$$        if(ixc==11) then
c$$$          print *," xxx2: memsize 8*ngb*ngb*nwhis=", 8*ngb*ngb*nwhis,' ngb nwhis=',ngb,nwhis
c$$$          allocate( rcxq(ngb,ngb,nwhis,npm) )
c$$$          rcxq=(0d0,0d0)
c$$$        else
c$$$          if(onceww(2)) print *," xxx2:allocate zxq zxqi memsize 16*ngb*ngb*(nwp+niw)=",
c$$$     &    16*ngb*ngb*(1+nwp+niw),' ngb nwp niw=',ngb,nwp,niw
c$$$          allocate(
c$$$     &    zxq (ngb,ngb,nw_i:nw),   !,nwp) feb2006
c$$$     &    zxqi(ngb,ngb,niw))
c$$$          zxq=0d0; zxqi=0d0
c$$$        endif

c     &    zxq (ngb,ngb,nwblk),
c     &    zxqi(ngb,ngb,niwblk))
c      do 1002 iwb = 1,nwdiv   !loop over frequency block
c      print *," xxx3"

!! takao apr2012
c        if(nolfco) then
c          if(allocated(rcxq)) deallocate(rcxq)
c          if(allocated(zxq) ) deallocate(zxq)
c          if(allocated(zxqi) ) deallocate(zxqi)
          allocate( rcxq(nmbas1,nmbas2,nwhis,npm) )
c          allocate( zxq (nmbas1,nmbas2,nw_i:nw), zxqi (nmbas1,nmbas2,niw))
c        else
          allocate( zw0(ngb,ngb)) !, zxq (ngb,ngb,nw_i:nw), zxqi(ngb,ngb,niw) )
c        endif
c        zxq=0d0;  zxqi=0d0;  
        rcxq = 0d0


c============ loop over spin============================================
ccccccccccccccccccccccccccccccccccccccccccc
c      nspinmx=1 !test
ccccccccccccccccccccccccccccccccccccccccccc
c=======================================================================
        do 1003 is = 1,nspinmx
          write(6,"(' *** ',2i4,' out of nqibz+n0qi nsp=',2i4,' *** ')") 
     &    iq, is, nqibz + nq0i,nspin
          if(debug) print *,' niw nw=',niw,nw
          isf=is
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



cc=======TetrahedronWeight_4 block. tettwt4  ixc==1 ==============2ini
c          if(tetra.and.ixc==1) then
c            stop 'hx0fp0sc: ixc==1 tetra is not impliemented'
c          endif !ixc==1



c=======TetrahedronWeight_5 block. tetwt5  ixc==,4,6,11 =======4ini
c          if(ixc==11) then !sf 21May02
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

ctakao july2012 
c$$$c- new version only along real axis for
c$$$          print *, ' --- goto x0kf_v3hz ---- newaniso= ',newaniso2
c$$$            call x0kf_v3hz( npm,ncc,
c$$$     i          ihw,nhw,jhw,whw,nhwtot,  ! tetwt5
c$$$     i            n1b,n2b,nbnbx,nbnb,    ! use whw by tetwt5 ,
c$$$     i              q,  
c$$$     i              nspin,is,is, 
c$$$     i              qbas,ginv,  qbz,wbz, 
c$$$     i              nblocha ,nlnm,nlnmv,nlnmc,iclass, 
c$$$     i              ppb(1,is), 
c$$$     i              icore,ncore,
c$$$     d            nlmto,nqbz,noccx,noccxv,
c$$$     d            nl,nclass,natom,nnc,
c$$$     d            nlnmx,mdimx,nbloch,  nwhis,
c$$$     i    iq,ngb,ngc,ngpmx,ngcmx,  !ngb/=ngc+nbloch for smbasis()=T oct2005
c$$$     i    nqbze,nband,nqibz, 
c$$$     o    rcxq, ! See rcxq_zxq in x0kf_v3h.f; rcxq is equivalent with zxq;
c$$$            !  rcxq is the accumulating variable for spins 
c$$$     i    nolfco,zzr,  !for nolfco
c$$$     o    rcxqmean   ) !for nolfco
c$$$
          print *, ' --- goto x0kf_v4hz ---- newaniso= ',newaniso2
          call x0kf_v4hz(npm,ncc,   
     i        ihw,nhw,jhw,whw,nhwtot,  ! tetwt5
     i        n1b,n2b,nbnbx,nbnb,    ! use whw by tetwt5 ,
     i              q,  
     i              nspin,is,isf,
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
     &    !  rcxq is the accumulating variable for spins 
     i    nolfco,zzr,nmbas_in,nmbas1,nmbas2, ppovlz, ppovl,zcousq,
cooo
c     i    .false.,!for nolfco Add nmbas Sep2006
     i    newaniso2,newaniso2,!for nolfco Add nmbas Sep2006
     o    rcxqmean,
     i    nwgt(:,iq),igx(:,:,iq),igxt(:,:,iq),ngrp, eibzsym(:,:,iq))
          print *,' end of x0kf_v4h sum rcxq=',sum(abs(rcxq))
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 8888     continue
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c          endif !ixc==11
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

c          if(ixc==11) 
            deallocate(ihw,nhw,jhw, whw,ibjb )
c          if(tetra.and.ixc==1) deallocate(wwk)
          if(tetra)            deallocate( n1b,n2b)
 1003   continue;print *, 'end of spin-loop nwp=',nwp  !end of spin-loop
c===========end of spin loop============================================
c=======================================================================

!! --- Hilbert transform.  Genrerate Real part from Imaginary part. ======
c$$$        if(ixc==11) then
c$$$          write(6,'(" xxx4: allocate zxq zxqi"//
c$$$     &    " =16*ngb*ngb*(nwp+niw)=",i8," ngb nwp niw=",3i7)')
c$$$     &    16*ngb*ngb*(1+nwp+niw), ngb, nwp, niw
c$$$          zxq =0d0
c$$$          zxqi=0d0
c$$$          nw_w=nw
c$$$          if(iSigma_en==0) nw_w=0 !i think this is ok , but not tested... feb2006
c$$$          write(6,*)"goto dpsion4 nwhis nw_i nw_w=",nwhis,nw_i,nw
c$$$          call dpsion4(frhis,nwhis, freq_r, nw_w, freq_i,niw, realomega, imagomega, 
c$$$     i        ngb, rcxq, npm,nw_i,  ! rcxq is alterd---used as work
c$$$     o        zxq, zxqi, 
c$$$     i        nolfco,.false., 0d0,is, rcxqmean,1d50,1d50, 
c$$$     o        x0mean)  !zxq x0mean is accumlating veriable for chipm mode
c$$$c          call dpsion3(frhis,nwhis, freq_r, nw_w, freq_i,niw,
c$$$c     i                             realomega,    imagomega,
c$$$c     i    ngb, rcxq,
c$$$c     o    zxq, zxqi, ! zxq... are Not the accumulating variable for spins
c$$$c     i    nolfco,.false.,0d0,is, rcxqmean,   1d50,1d50,
c$$$c     o    x0mean)
c$$$cc         do iw=1,nwp
c$$$cc         write(6,"(a,i5,2d13.5)")' zzz=',iw,sum(zxq(:,5,iw))
c$$$cc         enddo
c$$$        endif

          if(allocated(zxq) ) deallocate(zxq,zxqi)
          allocate(zxq (nmbas1,nmbas2,nw_i:nw), zxqi(nmbas1,nmbas2,niw))

          write(6,'("goto dpsion5: "
     &        //"nwhis nw_i niw nw_w nmbas1 nmbas2=",6i5)') nwhis,nw_i,nw,niw,nmbas1,nmbas2
!!              write(6,*)'sumchk rcxq goto dpsion5=', sum(abs(rcxq))
!! rcxqmean--->rcxq index conversion
        if(nolfco) then
            do iw=1,nwhis
                do jpm=1,npm
                  rcxq(:,:,iw,jpm)=rcxqmean(iw,jpm,:,:)
                enddo
            enddo
        endif
!! nolfco mode in dpsion5 is removed.
        print *,' -------- nmbas1,nmbas2=', nmbas1,nmbas2
        call dpsion5(frhis,nwhis, freq_r, nw, freq_i,niw, realomega, imagomega, 
     i        rcxq, npm,nw_i, nmbas1,nmbas2, ! rcxq is alterd---used as work
     o        zxq, zxqi,
     i        chipm, schi,is,  ecut,ecuts)
        print *,' --- end of dpsion5 ----',sum(abs(zxq)),sum(abs(zxqi))
        if(allocated(rcxq) ) deallocate(rcxq)

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
c          if(iSigma_en==0) then
c            nwmin=0
c            nwmax=0
c          endif
          if (anfexist) then   !Antiferro accelaration.
            do iw  = nwmin,nwmax
              call anfx0k(natom,nclass,nblocha,iclass,pos
     i     ,nbloch,ngc,
     i      q, ngveccB,qbas, ! for q+G
     i      anfvec,ibasf,             ! these are antiferro informations.
     i      zxq(1,1,iw))              ! i/o
            enddo
          endif

!! prepare for iq0.
          iq0 = iq - nqibz
          if(newaniso2) then
            if( iq==1 ) then
              ifepstinv = iopen('EPS0inv',0,-1,0)
              write(ifepstinv) ngb
            elseif(iq0>0) then
              write(ifisk) ngb,nw_i,nw,niw
              allocate( sk( ngb),  sks(ngb), w_k(ngb), w_ks(ngb) )
              allocate( skI( ngb),  sksI(ngb), w_kI(ngb), w_ksI(ngb) )
c            allocate( s_vc(ngb-1))
            endif
c            print *,' dddddddddd allocate epstilde dddddddddd'
cccc           if(iq<=nqibz) then
            allocate(epstilde(ngb,ngb))
ccccc           endif
            if(iq0 >= 1) then
              ifepstinv = iopen('EPS0inv',0,0,0)
              read(ifepstinv) ngb
            endif
            allocate(epstinv(ngb,ngb))
          endif

c          write (6,"(' end of anfx0k... q='3f10.4)") q(1),q(2),q(3)
          write(6, *)" === trace check for W-V === nwmin nwmax=",nwmin,nwmax
          do 1015 iw  = nwmin,nwmax
            frr= dsign(freq_r(abs(iw)),dble(iw))
            imode = 1
              if(newaniso2.and.iq<=nqibz) then !for mmmw
cc            if(.true.) then !for mmmz
cmmm

                if(iq==1) then
cc              if(iq==1.or.iq>nqibz) then
                  ix=1
                  zw0(:,1)=0d0
                  zw0(1,:)=0d0
                else
                  ix=0
                endif

!!  Eqs.(37),(38) in PRB81 125102
                do igb1=ix+1,ngb
                  do igb2=ix+1,ngb
                    epstilde(igb1,igb2)= -vcousq(igb1)*zxq(igb1,igb2,iw)*vcousq(igb2)
                    if(igb1==igb2) epstilde(igb1,igb2)=1+epstilde(igb1,igb2)
                  enddo
                enddo
                epstinv(ix+1:ngb,ix+1:ngb)=epstilde(ix+1:ngb,ix+1:ngb)
                call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
c$$$cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$cmmm direct inversion vs. block inversion
c$$$                if(iq>nqibz) then
c$$$c direct inversion
c$$$                ix=0
c$$$                do igb1=ix+1,ngb
c$$$                do igb2=ix+1,ngb
c$$$                  epstilde(igb1,igb2)= -vcousq(igb1)*zxq(igb1,igb2,iw)*vcousq(igb2)
c$$$                  if(igb1==igb2) epstilde(igb1,igb2)=1+epstilde(igb1,igb2)
c$$$                enddo
c$$$                enddo
c$$$                epstinv(ix+1:ngb,ix+1:ngb)=epstilde(ix+1:ngb,ix+1:ngb)
c$$$                call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
c$$$                do igb1=1+ix,ngb
c$$$                do igb2=1+ix,ngb
c$$$                  zw0(igb1,igb2)= vcousq(igb1)*epstinv(igb1,igb2)*vcousq(igb2)
c$$$                  if(igb1==igb2) zw0(igb1,igb2)= zw0(igb1,igb2)-vcousq(igb1)*vcousq(igb2)
c$$$                enddo
c$$$                enddo
c$$$c              write(*,"('mmmmzp99x  ',i3,10(2d13.5,2x))") iw,zw0(1,1),zw0(2:10:3,1),zw0(63:70:3,1)
c$$$              write(*,"('mmmmzp99x  ',i3,10(2d13.5,2x))") iw,1d0/epstinv(1,1),zw0(2:10:3,1),zw0(63:70:3,1)
c$$$c             write(*,"('mmmmzp99x  ',i3,10(2d13.5,2x))") iw,zw0(1,1),zw0(1,2:10:3),zw0(1,63:70:3)
c$$$c block inversion
c$$$                ix=1
c$$$                do igb1=ix+1,ngb
c$$$                do igb2=ix+1,ngb
c$$$                  epstilde(igb1,igb2)= -vcousq(igb1)*zxq(igb1,igb2,iw)*vcousq(igb2)
c$$$                  if(igb1==igb2) epstilde(igb1,igb2)=1+epstilde(igb1,igb2)
c$$$                enddo
c$$$                enddo
c$$$                epstinv(ix+1:ngb,ix+1:ngb)=epstilde(ix+1:ngb,ix+1:ngb)
c$$$                call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
c$$$                absq=sqrt(sum(q**2*tpioa**2))
c$$$                sk(  1:ngb)= zxq(1,1:ngb,iw)
c$$$                sks( 1:ngb)= zxq(1:ngb,1,iw)
c$$$                w_k(1) =0d0
c$$$                w_ks(1)=0d0
c$$$                w_k( 2:ngb)= vcousq(2:ngb)*vcousq(1)*matmul(vcousq(1)*sk(2:ngb)*vcousq(2:ngb),epstinv(2:ngb,2:ngb))
c$$$                w_ks(2:ngb)= vcousq(2:ngb)*vcousq(1)*matmul(epstinv(2:ngb,2:ngb),vcousq(1)*sks(2:ngb)*vcousq(2:ngb))
c$$$                llw(iw,iq0)=
c$$$     &             1d0
c$$$     &            -vcousq(1)*sk(1)*vcousq(1) ! sk(1,1,iw)=sks(1,1,iw)=H of Eq.(40).
c$$$     &            -vcousq(1)*vcousq(1)* sum( vcousq(2:ngb)*sk(2:ngb) * matmul(epstinv(2:ngb,2:ngb),sks(2:ngb)*vcousq(2:ngb)))
c$$$                write(*,"('mmmmzwp99x ',i3,10(2d13.5,2x))") iw,llw(iw,iq0), !(1d0/llw(iw,iq0)-1d0)*vcousq(1)**2,
c$$$c     &                  w_k(2:10:3)/llw(iw,iq0), w_k(63:70:3)/llw(iw,iq0)
c$$$     &                  w_ks(2:10:3)/llw(iw,iq0), w_ks(63:70:3)/llw(iw,iq0)
c$$$                write(*,"('mmmmzwp99x ')")
c$$$                endif
c$$$ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

                do igb1=1+ix,ngb
                  do igb2=1+ix,ngb
                    zw0(igb1,igb2)= vcousq(igb1)*epstinv(igb1,igb2)*vcousq(igb2)
                    if(igb1==igb2) zw0(igb1,igb2)= zw0(igb1,igb2)-vcousq(igb1)*vcousq(igb2)
                  enddo
                enddo
c                print *,'qqqqq sumcheck epstinv',sum(abs(epstinv(ix+1:ngb,ix+1:ngb)))
                if(iq==1) write(ifepstinv) epstinv(ix+1:ngb,ix+1:ngb),iq,iw
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cmmmm
c       if(iq>nqibz) then
c         write(*,"('mmmmz99x ',i3,10(2d13.5,2x))") iw,zw0(1,1)+vcousq(1)**2,zw0(2:10:3,1),zw0(63:70:3,1)
c       endif
c       if(iq==1.or.iq>nqibz) then
c                write(*,"('mmmz0  ',i3,10(2d13.5,2x))") iw,zxq(1,1,iw),zxq(1,2:10:3,iw),zxq(1,63:70:3,iw)
c                write(*,"('mmmz0* ',i3,10(2d13.5,2x))") iw,zxq(1,1,iw),zxq(2:10:3,1,iw),zxq(63:70:3,1,iw)
c         write(*,"('mmmmz99x ',i3,10(2d13.5,2x))") iw,zw0(1,1)+vcousq(1)**2,zw0(1,2:10:3),zw0(1,63:70:3)
c       write(*,"('mmmzx  ',2i3,10(2d13.5,2x))") iq,iw,zxq(2,1,iw),zxq(2,2:10:3,iw),zxq(2,63:70:3,iw)
c       write(*,"('mmmzx  ',2i3,10(2d13.5,2x))") iq,iw,zxq(3,1,iw),zxq(3,2:10:3,iw),zxq(3,63:70:3,iw)
c       write(*,"('mmmzxs ',2i3,10(2d13.5,2x))") iq,iw,zxq(1,1,iw),zxq(2:10:3,1,iw),zxq(63:70:3,1,iw)
c       write(*,"('mmmzxs ',2i3,10(2d13.5,2x))") iq,iw,zxq(1,2,iw),zxq(2:10:3,2,iw),zxq(63:70:3,2,iw)
c       write(*,"('mmmmzee',2i3,10(2d13.5,2x))")iq,iw,epstilde(2,2),epstilde(2,2:10:3),epstilde(2,63:70:3)
c       write(*,"('mmmmzee',2i3,10(2d13.5,2x))")iq,iw,epstilde(3,2),epstilde(3,2:10:3),epstilde(3,63:70:3)
c       endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
              endif

              if(newaniso2.and.iq>nqibz) then
cmmm3
                read(ifepstinv) epstinv(2:ngb,2:ngb),iqx,iwx
c                print *,'sumcheck epstinv=',sum(abs(epstinv(2:ngb,2:ngb)))
                if(iw/=iwx) stop 'hx0fp0: iw/=iwx'  !sanity check
cmmm3
                sk  (1:ngb)= zxq(1,1:ngb,iw)
                sks (1:ngb)= zxq(2,1:ngb,iw)
cc                sk(  1:ngb)= zxq(1,1:ngb,iw)
cc                sks( 1:ngb)= zxq(1:ngb,1,iw)

                w_k(1) =0d0
                w_ks(1)=0d0
                vcou1= fourpi/sum(q**2*tpioa**2) !test --> vcousq(1)**2!  !fourpi/sum(q**2*tpioa**2-eee)
                vcou1sq= sqrt(vcou1)   ! only vcousq(1) should be replaced.
                w_k( 2:ngb)= vcou1sq*matmul( sk(2:ngb)*vcousq(2:ngb), epstinv(2:ngb,2:ngb) )
                w_ks(2:ngb)= vcou1sq*matmul( epstinv(2:ngb,2:ngb), sks(2:ngb)*vcousq(2:ngb))
cmmm epsPP mode - vcou1sq*sum( sk(2:ngb) * w_ks(2:ngb)*vcousq(2:ngb) )
                llw(iw,iq0)=  1d0 -vcou1*sk(1) !- vcou1sq*sum( sk(2:ngb) * w_ks(2:ngb)*vcousq(2:ngb) )
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cmmmm
c                write(*,"('mmmw0  ',i3,10(2d13.5,2x))") iw,zxq(1,1,iw),zxq(1,2:10:3,iw),zxq(1,63:70:3,iw)
c                write(*,"('mmmw0* ',i3,10(2d13.5,2x))") iw,zxq(1,1,iw),zxq(2,2:10:3,iw),zxq(2,63:70:3,iw)
c                write(*,"('mmmmw99x ',i3,10(2d13.5,2x))") iw,fourpi/sum(q**2*tpioa**2)/llw(iw,iq0),
c     &                  w_k(2:10:3)/llw(iw,iq0),w_k(63:70:3)/llw(iw,iq0)
c                write(*,"('mmmmw99x ',i3,10(2d13.5,2x))") iw,llw(iw,iq0),
c     &                  w_ks(2:10:3)/llw(iw,iq0),w_ks(63:70:3)/llw(iw,iq0)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                write(*,*) 'epsPP iq iw',iq,iw, 1d0 - fourpi* sk(1)/sum(q**2*tpioa**2)
c                write(*,*) 'uuuu iq iw llw vc x0mean',iq,iw,fourpi/sum(q**2*tpioa**2),sk(1)
                write(ifisk) iw,iq0,q
                write(ifisk) vcousq(2:ngb)*w_k( 2:ngb),vcousq(2:ngb)*w_ks( 2:ngb)
              endif

              if(.not.newaniso2) then               ! Original mode
                call wcf( ngb, vcoul, zxq(1,1,iw), imode, zw0)
              endif

c$$$!!... a debug mode
c$$$                write(6,"('hhh --- EigenValues for Im( W) --------')")
c$$$                allocate(ebb(ngb))
c$$$                call diagcvh2( (zw0-transpose(dconjg(zw0)))/2d0/img, ngb, ebb)
c$$$                do ii=1,ngb
c$$$                  if( abs(ebb(ii))>1d-8 .and. ebb(ii)>0) then
c$$$                    write(6, "('hhhIWq : iw ii eb=',2i4,d13.5)") iw, ii, ebb(ii)
c$$$                  else
c$$$                    write(6, "('hhhIWqxxx : iw ii eb=',2i4,d13.5)") iw, ii, ebb(ii)
c$$$                  endif
c$$$                enddo
c$$$                deallocate(ebb)

cccccccccccccccccccccccccccccccccc
cooo
c            print *,'ddddddddddddd original like dddddddddddd iq ngb=',iq,ngb
c            if(iq==1) then
c               zw0=0d0
c            endif
cccccccccccccccccccccccccccccccccc

ccccccccccccccc
cooo
c              print *,'ddddd2222 original like mode dddddddddddddddd'
c              if(.false.) then
            if(newaniso2.and.iq>nqibz) then
c                zw(1:ngb,1:ngb) = 0d0
c                write(ifrcw, rec=((iq-iqxini)*(nw-nw_i+1)+ iw-nw_i+1 ) ) zw   !  WP = vsc-v
            else
                zw(1:ngb,1:ngb) = zw0
                write(ifrcw, rec=((iq-iqxini)*(nw-nw_i+1)+ iw-nw_i+1 ) ) zw   !  WP = vsc-v
                call tr_chkwrite("freq_r iq iw realomg trwv=", zw, iw, frr,nblochpmx, nbloch,ngb,iq)
            endif
 1015   continue  !iw
          if(newaniso2) then
            if(allocated(sk)) deallocate(sk,sks,w_k,w_ks)
          endif
          if( allocated(zzr)   ) deallocate(zzr)
c          if( allocated(x0mean)) deallocate(x0mean)
c          if( allocated(gbvec) ) deallocate(gbvec)
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
          do 1016 iw  = 1,niw
            if(newaniso2 .and. iq<=nqibz ) then
cc             if(.true.) then
!!  Eqs.(37),(38) in PRB81 125102
cmmm3
              if(iq==1) then
cc              if(iq==1.or.iq>nqibz) then
                ix=1
                zw0(:,1)=0d0
                zw0(1,:)=0d0
              else
                ix=0
              endif
!!  Eqs.(37),(38) in PRB81 125102
              do igb1=ix+1,ngb
                do igb2=ix+1,ngb
                  epstilde(igb1,igb2)= -vcousq(igb1)*zxqi(igb1,igb2,iw)*vcousq(igb2)
                  if(igb1==igb2) epstilde(igb1,igb2)=1+epstilde(igb1,igb2)
                enddo
              enddo
              epstinv=epstilde
              call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
              do igb1=ix+1,ngb
                do igb2=ix+1,ngb
                  zw0(igb1,igb2)= vcousq(igb1)*epstinv(igb1,igb2)*vcousq(igb2)
                  if(igb1==igb2) zw0(igb1,igb2)= zw0(igb1,igb2)-vcousq(igb1)*vcousq(igb2)
                enddo
              enddo
              if(iq==1) write(ifepstinv) epstinv(ix+1:ngb,ix+1:ngb),iq,iw
            endif
cmmm3
            if(newaniso2.and.iq>nqibz) then
              read(ifepstinv) epstinv(2:ngb,2:ngb),iqx,iwx
              if(iw/=iwx) then
                print *,'iw iwx=',iw,iwx
                stop 'hx0fp0: iw/=iwx' !sanity check
              endif
cmmm3
cc                skI  (1:ngb)= zxqi(1,1:ngb,iw)
cc                sksI (1:ngb)= zxqi(1:ngb,1,iw)
              skI  (1:ngb)= zxqi(1,1:ngb,iw)
              sksI (1:ngb)= zxqi(2,1:ngb,iw)
              w_kI(1) = 0d0
              w_ksI(1)= 0d0
              vcou1  = fourpi/sum(q**2*tpioa**2) ! test-->vcousq(1)**2 !fourpi/sum(q**2*tpioa**2-eee)
              vcou1sq= sqrt(vcou1)
c                w_kI( 2:ngb)= vcou1*vcousq(2:ngb)*matmul( skI(2:ngb)*vcousq(2:ngb), epstinv(2:ngb,2:ngb) )
c                w_ksI(2:ngb)= vcou1*vcousq(2:ngb)*matmul( epstinv(2:ngb,2:ngb), sksI(2:ngb)*vcousq(2:ngb))
              w_kI( 2:ngb)= vcou1sq*matmul( skI(2:ngb)*vcousq(2:ngb), epstinv(2:ngb,2:ngb) )
              w_ksI(2:ngb)= vcou1sq*matmul( epstinv(2:ngb,2:ngb), sksI(2:ngb)*vcousq(2:ngb))
cmmm epsPP mode ---> no - vcou1sq*sum( skI(2:ngb) * w_ksI(2:ngb)*vcousq(2:ngb) )
              llwI(iw,iq0)=  1d0 -vcou1*skI(1) !- vcou1sq*sum( skI(2:ngb) * w_ksI(2:ngb)*vcousq(2:ngb) )
ccccccc cmmm2 xxxxxx
              write(*,*) 'iq iw llwI',iq,iw,llwI(iw,iq0)
              write(ifisk) iw,iq0,q
              write(ifisk) vcousq(2:ngb)*w_kI(2:ngb),vcousq(2:ngb)*w_ksI( 2:ngb)
            endif

            if(.not.newaniso2) then               ! original mode
              call wcf( ngb, vcoul,zxqi(1,1,iw),imode,  zw0)
            endif

            if(newaniso2.and.iq>nqibz) then
c              zw(1:ngb,1:ngb) = 0d0 ! zw(nblochpmx,nblochpmx)
c              write(ifrcwi, rec=(iq-iqxini)*niw + iw)  zw    !  WP = vsc-v
            else
              zw(1:ngb,1:ngb) = zw0 ! zw(nblochpmx,nblochpmx)
              write(ifrcwi, rec=(iq-iqxini)*niw + iw)  zw    !  WP = vsc-v
              call tr_chkwrite("freq_i iq iw imgomg trwv=",zw,iw,freq_i(iw),nblochpmx,nbloch,ngb,iq)
            endif
 1016    continue
          if(newaniso2) then
            if(iq==1)    ifepstinv = iclose('EPS0inv') !iq==1 close write mode.
            if(iq>nqibz) ifepstinv = iclose('EPS0inv') !iq>nqibz --> readin mode
            deallocate(epstinv)
            if(allocated(epstilde)) deallocate(epstilde)
            if(allocated(skI)) deallocate(skI,sksI,w_kI,w_ksI)
            if(allocated(epstilde)) deallocate(epstilde)
          endif
c$$$          enddo
        endif !imagomega
C --- ImagOmega end ===============================

c 1002 continue  ! end of frequency block-loop
        deallocate(vcoul)
        if(allocated(zw0)) deallocate(zw0)
        if(allocated(zxq )) deallocate(zxq)
        if(allocated(zxqi)) deallocate(zxqi)
 1001 continue
c============ end of loop over q point =================================
c=======================================================================

!! == W(0) divergent part and W(0) non-analytic constant part.==
      if(newaniso2) then
        print *
        print *,' ==== newaniso2 mode W(0) divergent part ==== '
!! == W(0) divergent part ==
!! getw0 routine...
        ifidmlx = iopen('EPSwklm',0,0,0)
        read(ifidmlx) nq0i,lxklm
        allocate( dmlx(nq0i,9))
        allocate( epinvq0i(nq0i,nq0i),epinv(3,3,nq0i))
        nlxklm=(lxklm+1)**2
        allocate( wklm(nlxklm))
        read(ifidmlx) dmlx, epinv,epinvq0i
        read(ifidmlx) wklm
        ifidmlx = iclose(ifidmlx)
!! starting from llw(iw,iq0),llwI(iw,iq0)
!! === <e|L|e> (eq.36 in Friedrich paper) is expanded in YL -->stored in llwyl. ===
        allocate(w0(nw_i:nw),w0i(niw))
c       print *,' goto getw0 nq0i epinvq0i=',nq0i,epinvq0i
!! wbz(1) is the weight for q=0 = 1/(n1*n2*n3)
        call getw0(llw, nw_i,nw,nq0i,dmlx,epinvq0i,wklm,wbz(1), lxklm,  q0i,epinv,w0)
        call getw0(llwI,1,niw  ,nq0i,dmlx,epinvq0i,wklm,wbz(1), lxklm,  q0i,epinv,w0i)
c       print *,'sumcheck w0,w0i=',sum(abs(w0)),sum(abs(w0i))
!! === w0,w0i are stored to zw for q=0 ===
!! === w_ks*wk are stored to zw for iq >nqibz ===
        ifisk = iclose('WKcomponent')
        ifisk = iopen('WKcomponent',0,-1,0)
        do iq = 1,1 !iq=1 only 4pi/k**2 /eps part only ! iq = iqxini,iqxend
          q = qibze(:,iq)
          if(iq>nqibz) then
            iq0 = iq - nqibz
            read(ifisk) ngb,nw_ixxx,nwxxx,niwxxx
            allocate(vw_k(ngb),vw_ks(ngb))
          endif
          do ircw=1,2
            if    (ircw==1) then;  nini=nw_i;   nend=nw ;   ifrcwx=ifrcw
            elseif(ircw==2) then;  nini=1;      nend=niw;   ifrcwx=ifrcwi
            endif
            do iw=nini,nend
              if(iq<=nqibz) read(ifrcwx, rec=((iq-iqxini)*(nend-nini+1)+ iw-nini+1 ) ) zw !(1:ngb,1:ngb)
              if( iq==1 ) then
                if(ircw==1) zw(1,1) = w0(iw)
                if(ircw==2) zw(1,1) = w0i(iw)
              endif
              write(ifrcwx,rec=((iq-iqxini)*(nend-nini+1)+ iw-nini+1 ) ) zw !(1:ngb,1:ngb)
            enddo
          enddo
          if(iq>nqibz) deallocate(vw_k,vw_ks)
        enddo
      endif

c close files
      is = iclose('hbe.d')
c      is = iclose('EVU')
c      if(nspin == 2) is = iclose('EVD')
      call cputid(0)
c      if(ixc==1) stop ' OK! hx0fp0_sc mode=1'
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
