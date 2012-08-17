      program hsfp0_v2
c-------------------------------------------------------------
c     Calculates the  self-energy \Sigma within the
c     GW approximation at the LDA energy and two neighbouring points:
c
c     SEx(q,t,t) = <psi(q,t) |SEx| psi(q,t)>
c     SEc(q,t,t) = <psi(q,t) |SEc| psi(q,t)>
c
c     SEc(r,r';w) = (i/2pi) < [w'=-inf,inf] G(r,r';w+w') Wc(r,r';w') >
c     See details in sxcf_v2.f
c
c     This routine requirs an input from standard IO.
c     In scripts, you can do it like, prompt>echo mode|../exec/hsfp0 >lsc,
c     where mode is 1 or 2.
c
c     mode= 1: exchange    mode SEx, the exchange part of the self-energy
c     mode= 2: correlation mode SEc, the correlated part of the self-energy
c     mode= 3: core exchange mode SEXcore
c     mode= 4: plot spectrum function ---See manual
c     mode= 5: exchange    mode Exx, the exchange part of the total energy
c     mode= 6: core exchange mode Exxc, Exx(core-valence)
c
c     May 2002 Takashi Miyake. Total energy calc.
c     Apr 2002 takao kotani. multiple argumentation wave per l.
c     This hsfp0 is build from hsec10.f by F.Aryasetiawan.
c------------------------------------------------------------
      use m_readqg,only: readqg,readngmx
      use m_readeigen,only: init_readeigen,readeval,lowesteval,init_readeigen2
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
c     real(8),parameter ::
c     &  ua    = 1d0    ! constant in w(0)exp(-ua^2*w'^2) to take care of peak around w'=0
c------------------------------------
c     test switches to calculate the self-energy based on an another separation of \Sigma.
!     \Sigma = \Sigma_{sx} + \Sigma_{coh} + \Sigma_{img axis} + \Sigma_{pole} by Hedin PR(1965)A785
!     I found COH term has inevitably poor accuracy.
      logical ::tetra, tetra_hsfp0,         
     &     screen  = .false.,   ! \Sigma_{sx} for mode 1 and
      ! \Sigma_{img axis} + \Sigma_{pole} for mode 2
     &     cohtest = .false.,   ! \Sigma_{coh}. mode swich is not required.
     &     tetraex = .false.    ! This switch is only meaningful for mode=1,5,6
      ! If you want to calculate exchange, use tetraex=T .
      ! Note that you have to supply EFERMI by the tetrahedon method.
c------------------------------------
      real(8)    :: esmr2,shtw
      integer(4)::
     &     ixc,iopen,ifhbed, nprecb,mrecb,mrece,nlmtot,nqbzt, nband,
     &     ibas,ibasx,ngpmx,nxx,ngcmx,nbloch,ifqpnt,ifwd,
     &     nprecx,mrecl,nblochpmx2,nwp=0,niwt, nqnum,mdimx,nblochpmx,
!     nwp=0 for initialization March26-2007
     &     ifrcw,ifrcwi,  noccxv,maxocc2,noccx,ifvcfpout,iqall,iaf,ntq,
     &     i,k,nspinmx, nq,is,ip,iq,idxk,ifoutsex,iclose,nq0i,ig,
     &     mxkp,nqibzxx,ntet,nene,iqi, ix,iw,
     &     nlnx4,niwx,irot,invr,invrot,ivsum, ifoutsec,ntqx,
     &     ifsec(2)             !,ifcphi
     &     ,ifxc(2),ifsex(2), ifphiv(2),ifphic(2),ifec,ifexsp(2),
     &     ifsecomg(2),ifexx,ndble=8
      real(8) :: pi,tpia,vol,voltot,rs,alpha,
     &     qfermi,efx,valn,efnew,edummy,efz,qm,xsex,egex,
     &     zfac1,zfac2,dscdw1,dscdw2,dscdw,zfac,ef2=1d99,exx,exxq,exxelgas
      logical lqall,laf
      integer(4),allocatable :: itq(:)
      real(8),allocatable    :: q(:,:)

c     takao
      integer(4),allocatable :: ngvecpB(:,:,:), 
     &     ngvecp(:,:), ngvecc(:,:),iqib(:),    
     &     kount(:,:), nx(:,:),nblocha(:),lx(:)
      real(8),allocatable:: vxcfp(:,:,:),
     &     wqt(:), wgt0(:,:),q0i(:,:),
     &     ppbrd (:,:,:,:,:,:,:),cgr(:,:,:,:),eqt(:),
     &     ppbrdx(:,:,:,:,:,:,:),aaa(:,:), 
     &     ppb(:), eq(:),       
     &     eqx(:,:,:),eqx0(:,:,:),ekc(:),coh(:,:)
      complex(8),allocatable:: geigB(:,:,:,:) ,zsec(:,:,:)
c
      logical :: exchange, legas, tote
      real(8) ::  rydberg,hartree
      real(8):: qreal(3), ntot,nocctotg2,tripl,xxx(3,3)
      logical ::nocore

c     space group infermation
      integer(4),allocatable :: iclasst(:), invgx(:), miat(:,:)
      real(8),allocatable    :: tiat(:,:,:),shtvg(:,:)

c     tetra
      real(8),allocatable :: qz(:,:),qbzxx(:),wbzxx(:),wtet(:,:,:,:),
     &     eband(:,:,:), ene(:) !,ecore(:,:)
      integer(4),allocatable ::idtetx(:,:),idtet(:,:),ipq(:)
     &     ,iene(:,:,:),ibzx(:) ! ,nstar(:)

      integer(4) ::ib,iqx,igp,iii,ivsumxxx,isx,iflegas, iqpntnum
c
      real(8),allocatable   :: eex1(:,:,:),exsp1(:,:,:),qqex1(:,:,:,:)
      integer(4),allocatable:: nspex(:,:),ieord(:),itex1(:,:,:)
      real(8)    :: qqex(1:3), eex,exsp,eee, exwgt,deltax0
      integer(4) :: itmx,ipex,itpex,itex,nspexmx,nnex,isig,iex,ifexspx
     &     ,ifexspxx ,ifefsm, nq0ix,ifemesh,nz
      character(3)  :: charnum3,sss
      character(12) :: filenameex
      logical :: exspwrite=.false.
      character*8 xt

      integer(4) :: iwini,iwend
      real(8),allocatable:: omega(:,:)
      real(8) ::  omegamax,dwplot,omegamaxin
c     logical :: sergeys

      integer(4)::nqbze,ini,nq0it,idummy
      real(8),allocatable:: qbze(:,:)

      real(8)   :: ebmx(2)
      integer(4):: nbmx(2)

      real(8):: volwgt

      integer(4)::nwin, incwfin
      real(8)::efin,ddw
      integer(4),allocatable::imdim(:)
      real(8),allocatable::freqx(:),freqw(:),wwx(:) !,expa(:)

      logical:: GaussSmear      !readgwinput,
      integer(4)::ret
      character*(150):: ddd


      integer(4):: bzcase,  ngpn1,mrecg,verbose,ngcn1,nwxx
      real(8)   :: wgtq0p,quu(3)

      real(8),allocatable:: freq_r(:)

      logical ::smbasis
      integer(4):: ifpomat,nkpo,nnmx,nomx,ikpo,nn_,no,nss(2)
      real(8):: q_r(3)
      real(8),allocatable:: qrr(:,:)
      integer(4),allocatable:: nnr(:),nor(:)
      complex(8),allocatable:: pomatr(:,:,:),pomat(:,:)

      real(8)::sumimg
      logical allq0i            !S.F.Jan06

      integer(4):: nw_i

c     logical :: imgonly=.false.
      complex(8):: zseciip(-1:1)
      logical :: debug=.false.

!! newaniso mode
      real(8),allocatable:: vcousq(:),dmlx(:,:),epinvq0i(:,:),wklm(:),vcoud(:)
      complex(8),allocatable:: zcousq(:,:)
      logical:: newaniso
      integer:: ifvcoud,lxklm,ifidmlx

C...  mode switch. --------------
      if(cohtest) then
        screen = .true.
        ixc = 2; nz=0
        open(671,file='COH')
      else
        write(6,*) ' --- Choose omodes below ----------------'
        write(6,*) '  Sx(1) Sc(2) ScoreX(3) Spectrum(4) '
        write(6,*) '  EXX_val-val (5)  Exx_core-val(6) '
        write(6,*) '  Sx_sf(11) Sc_sf(12) '
        write(6,*) '  [option --- (+ QPNT.{number} ?)] '
        write(6,*) ' --- Put number above ! -----------------'
        call readin5(ixc,nz,idummy)
      endif
      call headver('hsfp0',ixc)
      call cputid (0) !initialization

      write(6,*) ' ixc nz=',ixc, nz
      if(ixc==0) stop ' --- ixc=0 --- Choose computational mode!'

      hartree=2d0*rydberg()
      tetraex=tetra_hsfp0()
      iii=verbose()
      print *,' verbose=',iii

C...  tetraex is only for Ex.
      if(ixc==1.or.ixc==5.or.ixc==6) then;
      else; tetraex=.false.
      endif


c     Sergey Fallev's sxcf_fal2
c     sergeys =.false.
      if(ixc==11.or.ixc==12) then
c     sergeys=.true.
        ixc=ixc-10
c     elseif(ixc==22) then
c     imgonly=.true.
c     ixc=2
      endif
c---  readin BZDATA. See gwsrc/rwbzdata.f
c--------readin data set when you call read_BZDATA ---------------
c     integer(4)::ngrp,nqbz,nqibz,nqbzw,nteti,ntetf
ccccc ! &   ,n_index_qbz
c     integer(4):: n1,n2,n3
c     real(8):: qbas(3,3),ginv(3,3),qbasmc(3,3)
c     real(8),allocatable:: qbz(:,:),wbz(:),qibz(:,:)
c     &    ,wibz(:),qbzw(:,:)
c     integer(4),allocatable:: idtetf(:,:),ib1bz(:),idteti(:,:)
c     &    ,nstar(:),irk(:,:),nstbz(:)          !,index_qbz(:,:,:)
c-----------------------------------------------------------------
      call read_BZDATA()
      print *,' nqbz  =',nqbz
c     print *,  qbz
      print *,' nqibz ngrp=',nqibz,ngrp
c     print *,' irk=',irk
c     print *,' #### idtetf: ####'
c     print *, idtetf

c     set up work array
c     call wkinit (iwksize)
      call pshprt(60)


C---  readin GWIN and LMTO, then allocate and set datas.
      nwin =-999                !not readin NW file
      efin =-999d0              !not readin EFERMI
      if    (ixc==3) then;  incwfin= -2 !core exchange mode
      elseif(ixc==5) then;  incwfin= -4 ! valence-valence Ex energy mode. ! See rgwinf called from genallcf_*
      elseif(ixc==6) then;  incwfin= -3 ! core-valence Ex energy mode.
      else               ;  incwfin= -1 !use 7th colmn for core at the end section of GWIN
      endif
      call genallcf_v3(nwin,efin,incwfin) !in module m_genallcf_v3
      if(ngrp/= ngrp2) stop 'ngrp inconsistent: BZDATA and LMTO GWIN_V2'
c---  These are allocated and setted.
c     integer(4)::  nclass,natom,nspin,nl,nn,nnv,nnc, ngrp,
c     o  nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, nctot,niw, !not readin nw
c     real(8) :: alat,ef, diw,dw,delta,deltaw,esmr
c     character(120):: symgrp
c     character(6),allocatable :: clabl(:)
c     integer(4),allocatable:: iclass(:)
c     &  ,nindxv(:,:),nindxc(:,:),ncwf(:,:,:) ,
c     o    invg(:), il(:,:), in(:,:), im(:,:),   ilnm(:),  nlnm(:),
c     o    ilv(:),inv(:),imv(:),  ilnmv(:), nlnmv(:),
c     o    ilc(:),inc(:),imc(:),  ilnmc(:), nlnmc(:),
c     o    nindx(:,:),konf(:,:),icore(:,:),ncore(:),
c     &    occv(:,:,:),unoccv(:,:,:)
c     &   ,occc(:,:,:),unoccc(:,:,:),
c     o    nocc(:,:,:),nunocc(:,:,:)
c     real(8), allocatable::
c     o  plat(:,:),pos(:,:),z(:),  ecore(:,:),  symgg(:,:,:) ! symgg=w(igrp),freq(:)
c-----------------------------------------------------------------------

c---  Get maximums takao 18June03
      call getnemx8(nbmx,ebmx)
      write(6,"('  nbmx ebmx from GWinput=',2i8,2d13.5)") nbmx,ebmx

c-------------------------------------------------------------------
c     if (nclass > mxclass) stop ' hsfp0: increase mxclass'
c!!!! WE ASSUME iclass(iatom)= iatom !!!!!!!!!!!!!!!!!!!!!!!!!
      if (nclass /= natom ) stop ' hsfp0: nclass /= natom ' ! We assume nclass = natom.
      print *,' hsfp0: end of genallcf2'
c
      call pshprt(30)
      pi   = 4d0*datan(1d0)
      tpia = 2d0*pi/alat

      shtw = 0d0
c     Jan2008 I don't think this is nesssary now... takao.
c     if(esmr<1d-5) shtw=0.01d0 ! Ferdi's shift to avoid resonance effect(maybe)
c
      call dinv33(plat,1,xxx,vol)
      voltot = dabs(vol)*(alat**3)

c---  ef is taken as rs for the empty-sphere test case of legas=T case -------------
      legas = .false.
      INQUIRE (FILE = 'LEGAS', EXIST = legas)
      if(legas) then            !!! test for electron gas case.
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
        if(tetraex) stop 'legas You have to give ef of  tetrahedron'
      endif

c
      ifexsp=0
      if(ixc==1) then
        exchange=.true.
        tote=.false.
        write(6,*) ' --- Exchange mode --- '
        ifxc(1)  = iopen('XCU'//xt(nz),1,-1,0)
        ifsex(1) = iopen('SEXU'//xt(nz),1,-1,0)
        if (nspin == 2) then
          ifxc(2)  = iopen('XCD'//xt(nz),1,-1,0)
          ifsex(2) = iopen('SEXD'//xt(nz),1,-1,0)
        endif

        INQUIRE (FILE = 'EXspTEST', EXIST = exspwrite)
        if(exspwrite) then
          print *,'--- Find EXspTEST ExspectrumWrite=',exspwrite
          print *,'--- esmr is chosen to be 2d0 Ry'
          esmr= 2d0
          do is=1,nspin
            ifexsp(is)  = iopen('EXSP.'//char(48+is),1,-1,0)
          enddo
        endif
      elseif(ixc==2) then
        exchange=.false.
        tote=.false.
        write(6,*) ' --- Correlation mode --- '
        if(cohtest) write(6,*) ' COH calculation mode. Results in COH'
        ifsec(1)   = iopen('SECU'//xt(nz),1,-1,0) ! output files
        if (nspin == 2) 
     .        ifsec(2)   = iopen('SECD'//xt(nz),1,-1,0)
      elseif(ixc==3) then
        exchange=.true.
        tote=.false.
        esmr=0d0
        write(6,*) ' --- CORE Exchange mode --- '
        ifsex(1)   = iopen('SEXcoreU'//xt(nz),1,-1,0)
        if (nspin == 2) then
          ifsex(2)   = iopen('SEXcoreD'//xt(nz),1,-1,0)
        endif
      elseif(ixc==4) then
        write(6,*) ' --- Spectrum function Sigma(omega) mode --- '
        exchange=.false.
        tote=.false.
c     ifsecomg(1) = iopen('SEComgU'//xt(nz),1,-1,0) ! output files
c     if (nspin == 2)
c     .  ifsecomg(2) = iopen('SEComgD'//xt(nz),1,-1,0)
      elseif(ixc==5) then
        exchange=.true.
        tote=.true.
        write(6,*) ' --- Exx mode valence-valence --- '
        ifexx = iopen('TEEXXvv',1,-1,0)
        esmr2 = esmr           !takao June2002
      elseif(ixc==6) then
        exchange=.true.
        tote=.true.
        esmr2 = esmr           !takao June2002
        esmr  = 0d0            !takao June2002
        write(6,*) ' --- CORE Exx mode core-valence --- '
        ifexx = iopen('TEEXXcv',1,-1,0)
      else
        stop ' hsfp0: Need input (std input) 1-6!'
      endif

c---  Neglect core is NoCore exists ------------
c     inquire(file='NoCore',exist=nocore)
c     if(nocore) nctot=0

c---
      write(6, *) ' --- computational conditions --- '
      write(6,'("    deltaw  =",f13.6)') deltaw
c     write(6,'("    ua      =",f13.6)') ua
      write(6,'("    esmr    =",f13.6)') esmr
      if(ixc==6) write(6,'("    esmr2    =",f13.6)') esmr2
      write(6,'("    alat voltot =",2f13.6)') alat, voltot
      if(nw/=-999) stop "nw here should be -999==nwin"
      write(6,'("    niw dw   =",i7,f13.6)') niw,dw

c>>   read dimensions of wc,b,hb
      ifhbed     = iopen('hbe.d',1,0,0)
      read (ifhbed,*) nprecb,mrecb,mrece,nlmtot,nqbzt, nband,mrecg
      if (nprecb == 4) stop 'hsfp0: b,hb in single precision'

      call init_readeigen(ginv,nspin,nband,mrece) !initialization of readEigen

c     --- get space group information ---------------------------------
c     true class information in order to determine the space group -----------
c     because the class in the generated GW file is dummy.(iclass(ibas)=ibas should be kept).
      open (102,file='CLASS')
      allocate(iclasst(natom),invgx(ngrp)
     &     ,miat(natom,ngrp),tiat(3,natom,ngrp),shtvg(3,ngrp))
      print *,'  --- Readingin CLASS info ---'
      do ibas = 1,natom
        read(102,*) ibasx, iclasst(ibas)
        write(6, "(2i10)") ibasx, iclasst(ibas)
      enddo

c     Get space-group transformation information. See header of mptaouof.
      call mptauof(symgg,ngrp,plat,natom,pos,iclasst
     o     ,miat,tiat,invgx,shtvg )
c     write (*,*)  'tiat=', tiat(1:3,1:natom,invr),invr

c     Get array size to call rdpp
c     call getsrdpp( nclass,nl,
c     o               ngpmx,ngcmx,nxx )
      call getsrdpp2( nclass,nl,nxx)
      call readngmx('QGpsi',ngpmx)
      call readngmx('QGcou',ngcmx)
      print *,' ngcmx ngpmx=',ngcmx,ngpmx
      allocate( nx(0:2*(nl-1),nclass), nblocha(nclass) ,lx(nclass),
     &     ppbrd ( 0:nl-1, nn, 0:nl-1,nn, 0:2*(nl-1),nxx, nspin*nclass),
     &     cgr(nl**2,nl**2,(2*nl-1)**2,ngrp))

c     - readin plane wave parts, and Radial integrals ppbrd.
!     ppbrd = radial integrals
!     cgr   = rotated cg coeffecients.
!     geigB = eigenfunction's coefficiens for planewave.
!     ngvecpB (in 1stBZ) contains G vector for eigen function.
!     ngveccB (in IBZ)   contains G vector for Coulomb matrix.
c     call rdpp_v2( ngpmx,ngcmx,nxx,  qibz,nqibz, qbz,nqbz,
c     i      nband, nl,ngrp, nn,  nclass, nspin, symgg,    qbas,
c     o      nblocha, lx, nx, ppbrd ,
c     o      mdimx, nbloch, cgr,
c     o      nblochpmx, ngpn,geigB,ngvecpB,  ngcni,ngveccB )
      call rdpp_v3(nxx, nl, ngrp, nn, nclass, nspin, symgg,qbas,
     o     nblocha, lx, nx,  ppbrd , mdimx, nbloch, cgr)
c     nblochpmx = nbloch + ngcmx !oct2005

c     allocate(ngcni(nqibz)) !, ngveccB(3,ngcmx,nqibz)) !, ngveccBr(3,ngcmx,nqibz))
c     geigB(ngpmx,nband,nqbz,nspin),ngpn(nqbz),ngvecpB(3,ngpmx,nqbz),
c     &   )  ! in IBZ
c
c     call rdpp_pln(ngpmx,ngcmx, qibz,nqibz, qbz,nqbz,nband,nspin,
c     o      ngpn,geigB,ngvecpB,ngcni,ngveccB)
c
c     do iq = 1,nqibz
c     call readqg('QGcou',qibz(1:3,iq),ginv,  quu,ngcni(iq), ngveccB(1,1,iq))
c     write(6,"('--From QGcou  qibz quu ngc=',3f9.4,'  ',3f9.4,i5)")
c     &       qibz(1:3,iq),quu,ngcni(iq)
c     enddo
c     for info
      allocate(ngvecp(3,ngpmx),ngvecc(3,ngcmx))
      call readqg('QGpsi',qibz(1:3,1),ginv, quu,ngpn1, ngvecp)
      call readqg('QGcou',qibz(1:3,1),ginv, quu,ngcn1, ngvecc)
      deallocate(ngvecp,ngvecc)

      print *, ' end of read QGcou'

cccccccccccccccccccccccccccccccccccccccccccccccccccc
c     iqx=1
c     do ib=1,nband
c     write(6,'("  iband iqx sumgeigB=",2i3,12d12.3)')
c     &   ib,iqx, sum(geigB(1:ngpn(1),ib,iqx))
c     enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     ib =1
c     iqx=1
c     do igp =1,nbloch + ngpn(iqx)
c     write(6,'("  igb ib iqx geigB=",3i3,1x,3i2,12d12.3)')
c     &   igp,ib,iqx, ngvecpB(1:3,igp,iqx), geigB(igp,ib,iqx)
c     enddo
c     stop "xxxxxx zzz"
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c----------------------------------------------
      call pshprt(60)


C---  Readin WV.d
      if(.not.exchange.or.(exchange.and.screen)) then !screened exchange
        ifwd      = iopen('WV.d',1,-1,0) !direct access files WVR and WVI which include W-V.
        read (ifwd,*) nprecx,mrecl,nblochpmx,nwp,niwt, nqnum, nw_i
        !nw_i added at Feb2006
c     nblochpmx from WV.d oct2005
        write(6,"(' Readin WV.d =', 10i5)") 
     &        nprecx, mrecl, nblochpmx, nwp, niwt, nqnum, nw_i
        call checkeq(nprecx,ndble)
c     call checkeq(nblochpmx,nblochpmx2)
c     if (nwp /= nw)   stop 'hsfp0: wrong nw' !july 2003
        nw=nwp-1
        if (niwt /= niw) stop 'hsfp0: wrong niw'
        ifrcw     = iopen('WVR',0,-1,mrecl)
        if(.not.exchange) ifrcwi = iopen('WVI',0,-1,mrecl)
C...  reading general energy mesh from file 'freq_r'
        open(UNIT=3111,file='freq_r') !this is in a.u.
        read(3111,*)nwxx       !number of energy points
        if(nwxx/= nw+1) stop ' freq_r nw /=nw'
        allocate(freq_r(nw_i:nw)) !freq_r(1)=0d0
        do iw = nw_i,nw
          read(3111,*)freq_r(iw)
        enddo
        close(3111)
      endif
c     write(6,*) freq_r

c...  Readin eigen functions
c     ifev(1)   = iopen('EVU', 0,0,mrece)
c     if (nspin==2) ifev(2) = iopen('EVD', 0,0,mrece)


c     -tetraex section ------------------------
      if(tetraex) then
        print *, ' !!!!  tetraex=T is not tested well !!!!!!!!!'
        ifefsm  = iopen('EFERMI',1,-1,0)
        read(ifefsm,*) ef
        write(6,"(' ef is readin from EFERMI',d23.15)")ef
c     --- get tetrahedron
c     mxkp = n1*n2*n3
c     allocate( qbzxx(3*mxkp),wbzxx(mxkp),ipq(mxkp) )
c     call bzmesh (plat,qbasmc,n1,n2,n3,w(igrp),ngrp,ipq,
c     .               qbzxx,wbzxx,nqibzxx,mxkp)
c     allocate(idtetx(0:4,mxkp*6))
c     call tetirr(qbasmc,n1,n2,n3,ipq,nqibz,ntet,
c     .              idtetx)
c     allocate(idtet(0:4,ntet))
c     idtet(0:4,1:ntet) = idtetx(0:4,1:ntet)
c     deallocate(idtetx,qbzxx,wbzxx,ipq)

c
c     nene=0
c     nene = ntq*nq*nspin ! for energy points.
c     if(exchange) nene=0
c     allocate(wtet(nband,nspin,nqibz,0:3*nene),
c     &    eband(nband,nspin,nqibz), qz(3,nqibz),nstar(nqibz),
c     &    iene(3*ntq,nq,nspin), ene(0:3*nene) ) ! pointer for
        allocate(wtet(nband,nspin,nqibz,0:0),
     &        eband(nband,nspin,nqibz), qz(3,nqibz)) ! ,nstar(nqibz))
        call dcopy (3*nqibz,qibz,1,qz,1)
        do  is  = 1,nspin      !Readin eband
          do  iqi = 1,nqibz
c     iq    = idxk (qz(1:3,iqi),qbz,nqbz)
c     call rwdd1 (ifev(is), iq, nband, eband(:,is,iqi))
            call readeval(qz(1:3,iqi),is, eband(:,is,iqi))
          enddo
        enddo
c     wtet(nband,nsp,nqibz,iene) where
c     the energy pointer as iene(itp,ip,ispin) corresponding its energy value.
c     ene(0) = ef
c     if(.not.exchange) then
c     ix =0
c     do is = 1,nspin
c     do ip = 1,nq
c     do i  = 1,ntq
c     do iw = -1,1
c     ix  = ix+1
c     iene(3*i+iw-1,ip,is) = ix
c     ene(ix) = eqx0(i,ip,is) + 2.d0*(dble(iw)-shtw)*deltaw
c     enddo
c     enddo
c     enddo
c     enddo
c     endif
c     do ix = 0,3*nene
c     ene(ix) = ene(ix)-1d-15  ! to avoid coincidence
c     call bzints2(n1,n2,n3,eband,wtet(:,:,:,ix),nqibz,nband,nband,
c     .                nspin,edummy,edummy,edummy,1,ene(ix),2,ntet,idtet)
c     enddo
c     ene(ix) = ene(ix)-1d-15  ! to avoid coincidence
ccccccccccccccc
c     call bzints2(n1,n2,n3,eband,wtet(:,:,:,0),nqibz,nband,nband,
c     .               nspin,edummy,edummy,edummy,1,ef,2,ntet,idtet)
        volwgt = (3d0 - nspin) / ntetf ! ntetf was =6*n1*n2*n3
        call bzints2x(volwgt,eband,wtet(:,:,:,0),nqibz,nband,nband,
c     ,dum,nkp,ldim,nbmax,
     .        nspin,edummy,edummy,edummy,1,ef,2,nteti,idteti)

        ntot= sum(wtet)
        if(legas) then
          write(6,"(' tetra=T ef ntot nexact ratio=',15f12.6)") ef,ntot
     &           , ef**1.5d0/3d0/pi**2*voltot, ef**1.5d0 /3d0/pi**2*voltot/ntot
        else
          write(6,"(' tetra=T ef nvalence)=',15f12.6)") ef,ntot
        endif
c
        if(nspin==1) wtet = wtet/2d0
c     call icopy(nqibz,w(instar),nstar)
        do iqi = 1,nqibz
          wtet(:,:,iqi,:) = wtet(:,:,iqi,:)/nstar(iqi)
        enddo
c     do ix = 1,3*nene
c     if( ene(ix)> ef) then
c     wtet(:,:,:,ix) = wtet(:,:,:,ix) - wtet(:,:,:,0)
c     else
c     wtet(:,:,:,ix) = wtet(:,:,:,0)  - wtet(:,:,:,ix)
c     endif
c     enddo
        deallocate( eband,qz)  !, ene ) ! pointer for
c     -- ibzx denote the index of k{FBZ for given k{1BZ.
        allocate(ibzx(nqbz))
c     call invkibzx(w(iirk),nqibz,ngrp,nqbz,
        call invkibzx(irk,nqibz,ngrp,nqbz,
     o        ibzx)

c     print *,' *** wtet sum 1=', sum( wtet)
c     stop ' -------------test end------------- '
c     call cputid(0)
        if (tote) then
          do i=1,nband
            if( sum(abs(wtet(i,:,:,0))) == 0d0 ) exit
          enddo
          if(  sum(abs(wtet(i+1:nband,:,:,0)))/=0d0) 
     &           stop ' hsfp0: wtetef sum err1'
          noccxv = i-1
        endif
c
c=======================
        goto 201
c=======================
      else
        allocate(wtet(1,1,1,1),ibzx(1)) !dummy
      endif
c     -end of tetra section --------------------------------------------

      if((ixc==3).or.(ixc==6)) then
        ef = lowesteval() -1d-3 !lowesteigen(nspin,nband,qbz,nqbz) - 1d-3 !lowesteb was
        if(maxval(ecore(:,1:nspin))>ef) then
          stop 'hsfp0 ixc=3/6:  ecore>evalence '
        endif

        if(ixc==6) then
          call efsimplef2a(nspin,wibz,qibz,ginv, 
     i           nband,nqibz
     i           ,konf,z,nl,natom,iclass,nclass
     i           ,valn, legas, esmr2, !!! valn is input for legas=T, output otherwise.
c     
     i           qbz,nqbz       !index_qbz, n_index_qbz,
     o           ,efnew)
          ef2 = efnew
          print *,' end of efsimple ef2 esmr2=',ef2,esmr2
        endif
c
      else                      ! if(esmr/=0d0) then
c     --- determine Fermi energy ef for given valn (legas case) or corresponding charge given by z and konf.
!     When esmr is negative, esmr is geven automatically by efsimplef.
        call efsimplef2a(nspin,wibz,qibz,ginv, 
     i        nband,nqibz
     i        ,konf,z,nl,natom,iclass,nclass
     i        ,valn, legas, esmr, !!! valn is input for legas=T, output otherwise.
c     
     i        qbz,nqbz          !index_qbz, n_index_qbz,
     o        ,efnew)
c
c     print *,' end of efsimple'
        ef = efnew
        if (ixc==5) ef2 = ef
c     - check total ele number -------
        ntot  = nocctotg2(nspin, ef,esmr, qbz,wbz, nband,nqbz) !wbz
        print *,' ef    =',ef
        print *,' esmr  =',esmr
        print *,' valn  =',valn
        print *,' ntot  =',ntot
c
c     if(ixc==1) then
c     ifefsm  = iopen('EFERMI.hsfp0.mode1',1,-1,0)
c     write(ifefsm,*) ef,'!ef by smearing. by hsfp0'
c     ifefsm  = iclose('EFERMI.hsfp0.mode1')
c     elseif(ixc==2) then
c     ifefsm  = iopen('EFERMI.hsfp0.mode2',1,-1,0)
c     write(ifefsm,*) ef,'!ef by smearing. by hsfp0'
c     ifefsm  = iclose('EFERMI.hsfp0.mode2')
c     endif
c     c
c     c        if(abs(valn-ntot)>1d-6) stop ' abs(valn-ntot)>1d-6'  !20001 May
c     c
cxx   elseif(legas) then
cxx   ef = efx
      endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(tote) then
        ddw= .5d0
        if(GaussSmear()) ddw= 10d0
        if(ixc==6) then
c     noccxv = maxocc (ifev,nspin, ef2+ ddw*esmr2, nband,nqbz)
          noccxv = maxocc2 (nspin, ef2+ ddw*esmr2, nband,qbz,nqbz)
        else
          noccxv = maxocc2 (nspin, ef + ddw*esmr,  nband, qbz,nqbz)
        endif
      endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     noccxv = maxocc (ifev,nspin, ef+0.5d0*esmr, nband,nqbz)  ! maximum no. of occupied valence states
!     maxocc seems to give (the maxmum number of occ + 1).


c-------------------------
 201  continue
c-------------------------

!
c     noccx  = noccxv + nctot
c     write(6,"(' --- noccx = noccxv + nctot    =',3i4)")
c     &          noccx,noccxv, nctot

c     write(6,"(' --- Total number of electrons =',d24.16)") ntot
c
c     if(legas) then
c     efx   = (ntot*3*pi**2/voltot)**(2d0/3d0)
c     qfermi= dsqrt(efx)
c     alpha = (9*pi/4d0)**(1d0/3d0)
c     write (6,*)' -g-  density parameter rs = ', alpha/qfermi
c     endif

c     ifcphi  = iopen('CPHI',0,0,mrecb)

      call init_readeigen2(mrecb,nlmto,mrecg) !initialize m_readeigen

c     the coulomb matrix for exchange.
      if(exchange) ifvcfpout = iopen('VCCFP',0,-1,0)

c     QPNT data
      if(.not.tote) then
        if(nz==0) then
c     if(readgwinput()) then
          call getkeyvalue("GWinput","<QPNT>",unit=ifqpnt,status=ret)
c     else
c     ifqpnt = iopen('QPNT',1,0,0)
c     endif
        else
          ifqpnt  = iopen('QPNT'//xt(nz),1,0,0)
        endif
      endif
      print *,' ifqpnt ret=',ifqpnt,ret
c     do i=1,ret
c     read(ifqpnt,'(a)') ddd
c     write(6,'(a)') ddd
c     enddo



c     read q-points and states
      if(ixc==4) then
        call readxx(ifqpnt)    !skip to ***** for q point for spectrum function.
        omegamaxin = 1d70
        read (ifqpnt,*,err=2038,end=2038) dwplot,omegamaxin
 2038   continue
        omegamax = 2*freq_r(nw) !This is in Ry.
c     print *,' nwp dw(a.u.) omegamax (Ry)=',nwp,dw,omegamax
        if(omegamaxin < omegamax) then
          write(6,*)' --Use readin dwplot and omegamaxin from <QPNT>'
          omegamax = omegamaxin
        endif
        if( omegamax <0) stop 'hsfp0 :strange omegamax <0 '
        iwini =  -int( omegamax / dwplot )
        iwend =   int( omegamax/  dwplot )
        print *,' iwini:iwend omegamax(Ry)=',iwini,iwend,omegamax
c     elseif(imgonly) then
c     iwini=0
c     iwend=0
      else
        iwini = -1
        iwend = 1
      endif
c     ---
      if (tote) then
        lqall      = .true.
      else
        lqall      = .false.
        laf        = .false.
        call readx   (ifqpnt,10)
        read (ifqpnt,*) iqall,iaf
        if (iqall == 1) lqall = .true.
        if (iaf   == 1)   laf = .true.
        call readx   (ifqpnt,100)
      endif
!     states
      if (tote) then
        ntq    = noccxv
        allocate( itq(ntq) )
        do i   = 1,ntq
          itq(i) = i
        enddo
      else
        read (ifqpnt,*) ntq
        allocate( itq(ntq) )
        read (ifqpnt,*) (itq(i),i=1,ntq)
      endif
!     q-points
      if (lqall) then           !all q-points case
        nq         = nqibz
        allocate(q(3,nq))
        call dcopy   (3*nqibz,qibz,1,q,1)
      else
        call readx   (ifqpnt,100)
        read (ifqpnt,*) nq
        allocate(q(3,nq))
        do       k = 1,nq
          read (ifqpnt,*) i,q(1,k),q(2,k),q(3,k)
          write(6,'(i3,3f13.6)') i,q(1,k),q(2,k),q(3,k)
        enddo
      endif
c
      nspinmx = nspin
      if (laf) nspinmx =1
      close(ifqpnt)

C...
      call getkeyvalue("GWinput","QPNT_nbandrange",nss,2 !nov2005
     &     ,default=(/-99997,-99997/) )
      if(nss(2)/=-99997) then
        if(allocated(itq)) deallocate(itq)
        ntq=nss(2)-nss(1)+1
        allocate( itq(ntq) )
        do i=max(1,nss(1)),min(nss(2),nband)
          itq(i-nss(1)+1) = i
        enddo
      endif

C...
      allocate(omega(ntq,iwini:iwend))
      if(ixc==4) then
        do iw = iwini,iwend
          omega(1:ntq,iw) =  dwplot* iw + ef
        enddo
      endif

C     ---  q near zero
      print *, 'reading QOP'
      open (101,file='Q0P')
      read (101,"(i5)") nq0i
      if(.not.exchange) call checkeq(nqibz+nq0i-1, nqnum)
      write(6,*) ' *** nqibz nq0i_total=', nqibz,nq0i
      nq0it = nq0i
      allocate( wqt(1:nq0i),q0i(1:3,1:nq0i) )
c     read (101,"(d24.16,3x, 3d24.16)" )( wqt(i),q0i(1:3,i),i=1,nq0i)
      nq0ix = nq0i
      do i=1,nq0i
        read (101,* ) wqt(i),q0i(1:3,i)
        if(wqt(i)==0d0 ) nq0ix = i-1
      enddo
      nq0i = nq0ix              ! New nq0i July 2001
      write(6,*) ' Used k number in Q0P =', nq0i
      write(6,"(i3,f14.6,2x, 3f14.6)" )(i, wqt(i),q0i(1:3,i),i=1,nq0i)
      close(101)
      allocate( wgt0(nq0i,ngrp) )
C     Sergey's 1stFeb2005
c     call q0iwgt2(symgg,ngrp,wqt,q0i,nq0i,
c     o            wgt0)
      call getkeyvalue("GWinput","allq0i",allq0i,default=.false.) !S.F.Jan06
      call q0iwgt3(allq0i,symgg,ngrp,wqt,q0i,nq0i, !S.F.Jan06
     o     wgt0)                ! added allq0i argument
C--------------------------

      if(nq0i/=0) write(6,*) ' *** tot num of q near 0   =', 1/wgt0(1,1)
      write(6,"('  sum(wgt0) from Q0P=',d14.6)")sum(wgt0)
      if(bzcase()==2) then
        wgt0= wgt0*wgtq0p()/dble(nqbz)
        write(6,"('bzcase=2:  sum(wgt0_modified )=',d14.6)")sum(wgt0)
      endif

c     --- qbze(3,nqibze)
      nqbze  = nqbz *(1 + nq0it)
      allocate( qbze(3, nqbze) )
      call dcopy(3*nqbz, qbz, 1, qbze,1)
      do i = 1,nq0it
        ini = nqbz*(1 + i -1)
        do ix=1,nqbz
          qbze (:,ini+ix)   = q0i(:,i) + qbze(:,ix)
        enddo
      enddo

c     --- read LDA eigenvalues
c     ntp0=ntq
      allocate(eqx(ntq,nq,nspin),eqx0(ntq,nq,nspin),eqt(nband))
      do      is = 1,nspin
        do      ip = 1,nq
c     iq       = idxk (q(1,ip),qbze,nqbze)
c     call rwdd1   (ifev(is), iq, nband, eqt) !direct access read b,hb and e(q,t)
          call readeval(q(1,ip),is,eqt)
c     print *,' eqt=',eqt
          eqx0(1:ntq,ip,is) = eqt(itq(1:ntq))
          eqx (1:ntq,ip,is) = rydberg()*(eqt(itq(1:ntq))- ef)
        enddo
      enddo
      deallocate(eqt)

c     --- info
      if(tote) then
        call winfo2(6,nspin,nq,ntq,is,nbloch
     &        ,ngpn1,ngcn1,nqbz,nqibz,ef,ef2,deltaw,alat,esmr,esmr2)
      else
        call winfo(6,nspin,nq,ntq,is,nbloch
     &        ,ngpn1,ngcn1,nqbz,nqibz,ef,deltaw,alat,esmr)
      endif

!! === readin Vcoud and EPSwklm for newaniso()=T ===
      if(newaniso()) then
        ifidmlx = iopen('EPSwklm',0,0,0)
        read(ifidmlx) nq0ix,lxklm
        if(nq0i/=nq0ix) then
          print *,'nq0i from EPSwklm /= nq0i',nq0i,nq0ix
          stop 'nq0i from EPSwklm /= nq0i'
        endif
        allocate( dmlx(nq0i,9))
        allocate( epinvq0i(nq0i,nq0i) )
        allocate( wklm((lxklm+1)**2))
        read(ifidmlx) dmlx, epinvq0i
        read(ifidmlx) wklm
        ifidmlx = iclose(ifidmlx)
      endif


c-------------------------
c     LDA exchange-correlation
c-------------------------
c     - takao Vxc is readin from VXCFP
      if(ixc==1) then
        allocate(  vxcfp(ntq,nq,nspin) )
        call rsexx2(nspin,itq,q,ntq,nq, ginv, symgg,ngrp, vxcfp) !add ginv july2011
c     loop over spins
        do is = 1,nspinmx
          write (ifxc(is),*) '==================================='
          write (ifxc(is),"(' LDA exchange-correlation : is=',i3)")is
          write (ifxc(is),*) '==================================='
          call winfo(ifxc(is),nspin,nq,ntq,is,nbloch
     &           ,ngpn1,ngcn1,nqbz,nqibz,ef,deltaw,alat,esmr)
          write (ifxc(is),*)' ***'
          write (ifxc(is),"(a)") ' jband   iq ispin                  
     &           qvec 
     &           eigen-Ef (in eV)
     &           LDA XC (in eV)'
          ifoutsex = ifxc(is)
          write(6,*)
          do ip = 1,nq
            do i  = 1,ntq
              write(ifoutsex,"(3i5,3d24.16,3x,d24.16,3x,d24.16)")
     &                 itq(i),ip,is, q(1:3,ip), eqx(i,ip,is), 
     &                 vxcfp(i,ip,is)
              write(6,"(' j iq isp=' i3,i4,i2,'  q=',3f8.4,
     &                 '  eig=',f10.4,'  Sxc(LDA)=',f10.4)") 
     &                 itq(i),ip,is, q(1:3,ip), eqx(i,ip,is),
     &                 vxcfp(i,ip,is)
            end do
          end do
c     end of spin-loop
          if(is==1) isx = iclose('XCU'//xt(nz))
          if(is==2) isx = iclose('XCD'//xt(nz))
        end do
        deallocate(vxcfp)
      endif


c     pointer to optimal product basis
      allocate(imdim(natom))
      call indxmdm (nblocha,nclass,
     i     iclass,natom,
     o     imdim )

      if(niw/=0) then           ! generate gaussian frequencies x between (0,1) and w=(1-x)/x
        allocate(freqx(niw),freqw(niw),wwx(niw)) !,expa(niw))
        call freq01x  (niw,    !ua,
     o        freqx,freqw,wwx)  !, expa)
      endif

c     ------ write energy mesh ----------
c     if(.not.sergeys) then
c     write(6,*)"niw=, nwp",niw,nwp
      ifemesh = iopen('emesh.hsfp0'//xt(nz),1,-1,0)
      deltax0 = 0d0
      call writeemesh(ifemesh,freqw,niw,freq_r,nwp,deltax0)
c
      iii=ivsumxxx(irk,nqibz*ngrp)
      write(6,*) " sum of nonzero iirk=",iii, nqbz
C...  Read pomatr
      if(smbasis()) then
        print *,' smooth mixed basis : augmented zmel'
        call getngbpomat(nqibz+nq0i,
     o        nnmx,nomx)
        nkpo = nqibz+nq0i
        ifpomat = iopen('POmat',0,-1,0) !oct2005
        allocate( pomatr(nnmx,nomx,nkpo),qrr(3,nkpo),nor(nkpo),nnr(nkpo) )
        do ikpo=1,nkpo
          read(ifpomat) qrr(:,ikpo),nn_,no,iqx !readin reduction matrix pomat
c     write(6,"('smbasis: ikp q no nn=',i5,3f8.4,4i5)") ikp,qrr(:,ikpo),no,nn_
          nnr(ikpo)=nn_
          nor(ikpo)=no
          read(ifpomat) pomatr(1:nn_,1:no,ikpo)
        enddo
        isx = iclose("POmat")
        write(6,*)"Read end of POmat ---"
      else                      !dummy
        nkpo = 1
        nnmx =1
        nomx =1
        allocate( pomatr(nnmx,nomx,nkpo), qrr(3,nkpo),nor(nkpo),nnr(nkpo) )
      endif

c-----------------------------------------------------------
c     calculate the correlated part of the self-energy SEc(qt,w)
c-----------------------------------------------------------
      nlnx4    = nlnx**4
      niwx     = max0 (nw+1,niw)
      allocate( ppb(nlnmx*nlnmx*mdimx*nclass),  eq(nband), 
     &     kount(nqibz,nq), zsec(iwini:iwend,ntq,nq),  
     &     coh(ntq,nq) )


c     loop over spin ----------------------------------------------------
      if (tote) exx = 0d0
      do 2000 is = 1,nspinmx

        if((exchange).and.(.not.tote)) then
          write(ifsex(is),*) '======================================='
          write(ifsex(is),"('Self-energy exchange SEx(q,t): is=',i3)") is
          write(ifsex(is),*) '======================================='

          call winfo(ifsex(is),nspin,nq,ntq,is,nbloch
     &           ,ngpn1,ngcn1,nqbz,nqibz,ef,deltaw,alat,esmr)

          write (ifsex(is),*)' *** '
          write (ifsex(is),"(a)") ' jband   iq ispin                  
     &           qvec 
     &           eigen-Ef (in eV)
     &           exchange (in eV)'

        elseif(ixc==2) then
          write(ifsec(is),*) '=========================================='
          write(ifsec(is),"('Self-energy correlated SEc(qt,w): is=',i3)")
     &           is
          write(ifsec(is),*) '=========================================='

          call winfo(ifsec(is),nspin,nq,ntq,is,nbloch
     &           ,ngpn1,ngcn1,nqbz,nqibz,ef,deltaw,alat,esmr)

          write (ifsec(is),*)' *** '
          write (ifsec(is),"(a)") ' jband   iq ispin                  
     &           qvec 
     &           eigen-Ef (in eV)
     &           Re(Sc) 3-points (in eV)                        
     &           In(Sc) 3-points (in eV)     
     &           Zfactor'
        endif

c     initialise secq and kount
        zsec  = 0d0
        coh   = 0d0
        kount = 0
c
!make dummy SEXcore
        if(ixc==3.and.nctot==0) goto 2001

c     loop over rotations -------------------------------
        do 1000 irot = 1,ngrp
c          print *,'uuuuuuuu takao do 1000 ',sum(abs( irk(:,irot) )),sum(abs( wgt0(:,irot)))
          if( sum(abs( irk(:,irot) )) ==0 .and.
     &           sum(abs( wgt0(:,irot))) == 0d0 ) then
            cycle
          endif
          write (6,"(i3,'  out of ',i3,'  rotations ',$)") irot,ngrp
          call cputid (0)
c     if (irot == 1 .or. irot == ngrp) then
c     call cputid(0); write(*,*)' ppba '
c     endif

c     rotate atomic positions invrot*R = R' + T
          invr       = invrot (irot,invg,ngrp)

c     -- ppb= <Phi(SLn,r) Phi(SL'n',r) B(S,i,Rr)>
          if(debug) print *,'### goto ppbafp_v2 ###'
          call ppbafp_v2 (irot,ngrp,is,nspin,
     i           il,in,im, nlnm,
     d           nl,nn,nclass,nlnmx,
     i           mdimx,lx,nx,nxx, !Bloch wave    
     i           cgr, nl-1,     !rotated CG
     i           ppbrd,         !radial integrals
     o           ppb)

c     c -- Rotated gvecc
c     call rotgvec(symgg(:,:,irot), nqibz,
c     i    ngcmx,ngcni,qbas,ngveccB,
c     o    ngveccBr)

c------------------------------------------------------
c     calculate the correlated part of the self-energy within GW
c     ntqx = 0
c     if(tetra.and.(.not.exchange)) then
c     ntqx =3*ntq
c     endif
c     loop over q
          do 1001 ip = 1,nq   !
            if(debug) write (*,*) ip,'  out of ',nq,'  k-points ' ! call cputid  (0)
            if(ixc==2) then
c     iq  = idxk (q(1,ip),qbz,nqbz)
c     call rwdd1 (ifev(is),iq, nband,eq)
              call readeval(q(1,ip),is,eq)
              do iw = iwini,iwend
                do i  = 1,ntq
                  omega(i,iw) = eq(itq(i)) + 2d0*(dble(iw)-shtw)*deltaw
                enddo
              enddo
            endif

            if(nctot==0) then
              deallocate(ecore)
              allocate(ecore(1,2)) !this is dummry
            endif

c     if(sergeys) then
            if(debug) print *,'### goto sxcf_fal3z ### nlmto=',nlmto
            call sxcf_fal3z (kount, ip,q(1,ip),irot,itq,ntq,ef,ef2,esmr,esmr2,
     i              tiat(1:3,1:natom,invr),miat( 1:natom,invr),
c     v2
     i              nspin,is,   !ifcphi,ifrb(is),ifcb(is),ifrhb(is),ifchb(is),
     i              ifrcw,ifrcwi,
     i              qbas,ginv,qibz,qbz,wbz,nstbz, wibz, !iindxk,
     i              nstar,irk,  !kount,
c
c     i        iiclass,nblocha,i_mnlv,i_mnlc,iicore,incore,iimdim,
     i              iclass,nblocha,nlnmv, nlnmc, !w(i_mnlv),w(i_mnlc)
     i              icore,ncore, imdim,
     i              ppb,        ! pdb,dpb,ddb,
     i              freq_r,freqx, wwx, !expa,ua,
     i              dw,         !deltaw,
     i              ecore(:,is),
c     
     d              nlmto,nqibz,nqbz,nctot,
c     i        index_qbz, n_index_qbz,
     d              nl,nnc,nclass,natom,
     d              nlnmx,mdimx,nbloch,ngrp, nw_i,nw,  niw,niwx,nq,
!     nw_i,nw instead of nwp feb2006
c
c     i     nblochpmx, ngpn,ngcni,ngpmx,ngcmx,
     i              nblochpmx,ngpmx,ngcmx,
c     i     geigB(1,1,1,is), ngvecpB,ngveccBr,
c     i     ngveccBr,
     i              wgt0,nq0i,q0i, symgg(:,:,irot),alat,
     i              matmul(symgg(:,:,irot),shtvg(:,invr)),nband,
     i              ifvcfpout,
c     i     shtw, 
     i              exchange, tote, screen, cohtest, ifexsp(is),
     i              omega, iwini,iwend,
     i              nbmx(2),ebmx(2), !takao 18June2003 
     i              pomatr, qrr,nnr,nor,nnmx,nomx,nkpo, !imgonly,
c     i    pos,wklm,lxklm,zcousq,vcoud,ngb0, !mar2012
     i    wklm,lxklm, !mar2012
            !oct2005 for pomat 
     o              zsec,coh,exx) ! acuumulation variable
            print *,' end of sxcf_fal3z xxxxxxxxxxxxxxx'

c     else                                                 !june2002 add esmr2 takao
c     call sxcf_v2z (kount, ip,q(1,ip),irot,itq,ntq,ef,ef2,esmr,esmr2,
c     i              tiat(1:3,1:natom,invr),miat( 1:natom,invr),
c     i              nspin,is, !ifcphi,ifrb(is),ifcb(is),ifrhb(is),ifchb(is),
c     i              ifrcw,ifrcwi,
c     i              qbas,ginv,qibz,qbz,wbz,nstbz, wibz,
c     i              nstar,  irk,     !iindxk,kount,
c     i        iclass,nblocha,nlnmv, nlnmc, !w(i_mnlv),w(i_mnlc)
c     i        icore,ncore, imdim,
c     i              ppb, ! pdb,dpb,ddb,
c     i              freq, freqx, wwx, expa,
c     i              ua,dw,  !deltaw,
c     i              ecore(:,is),
c     d        nlmto,nqibz,nqbz,nctot,
c     d              nl,nnc,nclass,natom,
c     d              nlnmx,mdimx,nbloch,ngrp,nwp,niw,niwx,nq,
c     i     nblochpmx,ngpmx,ngcmx,
c     i     wgt0,nq0i,q0i, symgg(:,:,irot),alat,
c     i     matmul(symgg(:,:,irot),shtvg(:,invr)),nband,
c     i     ifvcfpout,
c     i     exchange, tote, screen, cohtest, ifexsp(is),
c     &     wtet(:,is,:,0),
c     &     ibzx,tetraex,
c     &     ixc,
c     c
c     i        omega, iwini,iwend,
c     i     nbmx,ebmx, !takao 18June2003
c     o        zsec,coh,exx) ! acuumulation variable
c     endif
cccccccccccccccccccccccccccccccccccccccccc
c     iii = ivsum(kount,nqibz*nq)
c     print *," sumkount 2=",nqibz,nq,iii
c     stop "--- kcount test end --- "
cccccccccccccccccccccccccccccccccccccccccc
 1001     continue            !end of q-loop
 1000   continue               !end of rotation-loop
cccccccccccccccccccccccc

c     check that all k { FBZ have been included
c     if (ivsum(kount,nqibz*nq) /= nqbz*nq) then
c     iii = ivsum(kount,nqibz*nq)
c     print *," ivsum=",iii, nqbz*nq
c     stop 'hsfp0: missing k-pts'
c     endif

c---------------------------------
c     electron gas bare exchange (exact)
        if (legas.and.exchange.and.(.not.tote)) then
          efz=(ntot*3*pi**2/voltot)**(2d0/3d0) ! ef is calculated from ntot.
          pi         = 4.d0*datan(1.d0)
          tpia       = 2.d0*pi/alat
          qfermi= dsqrt(efz)
          alpha = (9*pi/4d0)**(1d0/3d0)
          write (6,*)' --- exact electron gas bare exchange --- '
          write (6,*)' density parameter rs= ', alpha/qfermi
          write (6,*)' kf= ',qfermi
          do      ip = 1,nq
            qreal =  tpia*q(1:3,ip)
            qm    = dsqrt ( sum(qreal**2) )
            xsex  = hartree * egex (qm,efz)
            write (6,*)
            write (6,"(' True qm-ef Sx=',2f14.6,' q/qf=',f14.6)")
     &              rydberg()*(qm**2-efz), xsex, qm/qfermi
            write (6,"(' Num  qm-ef Sx=',2f14.6)") 
     &              eqx(1,ip,is),        hartree*dreal(zsec(iwini,1,ip))
            write (6,"(' === diff     =',2f14.6)") 
     &              rydberg()*(qm**2-efz)-eqx(1,ip,is)
     &              , xsex - hartree*dreal(zsec(iwini,1,ip))

            write (661,"(' qm True qm-ef Sx=',3f14.6)") 
     &              qm,rydberg()*(qm**2-efz), xsex
            write (662,"(' qm Num  qm-ef Sx=',3f14.6)") 
     &              qm,eqx(1,ip,is),     hartree*dreal(zsec(iwini,1,ip))
ccc   write (ifsex(is),6600) qreal(1),qreal(2),qreal(3),xsex
ccc   write (6,6600) qreal(1),qreal(2),qreal(3),xsex
ccc   6600   format (' qreal =',3f8.4,'   SEx(q) =',d13.5)
            write (663,"(2f14.6)") qm/qfermi, qfermi
          end do
        endif

 2001   continue
c---------------------------------

c     write SEx(q,t)
c     call cv      (hartree,sexq,ntq*nq,sexq) !Because vcoul is in a.u.

        if(ixc==5.or.ixc==6) then
!     nothing
        elseif(exchange) then
c     do       ip = 1,nq
c     call wsex2 (ifsex(is),q(1:3,ip),itq(1:ntq),eqx(1:ntq,ip,is),
c     &              hartree*dreal(zsec(-1,1:ntq,ip)), !sexq(1:ntq,ip),
c     d              ntq)
c     end do
c
c-------------------------
          ifoutsex=ifsex(is)
          write(6,*)
          do ip = 1,nq
            do i  = 1,ntq
              write(ifoutsex,"(3i5,3d24.16,3x,d24.16,3x,d24.16)")
     &                 itq(i),ip,is, q(1:3,ip), eqx(i,ip,is), 
     &                 hartree*dreal(zsec(iwini,i,ip))
              write(6,"(' j iq isp=' i3,i4,i2,'  q=',3f8.4,
     &                 '  eig=',f10.4,'  Sx=',f10.4)")
c     &     '  eig=',f10.4,'  Sx=',d15.7)") 
     &                 itq(i),ip,is, q(1:3,ip), eqx(i,ip,is),
     &                 hartree*dreal(zsec(iwini,i,ip))
            end do
          end do
c-------------------------
        elseif(ixc==2) then
c     write SEc(qt,w)
c     write (ifsec(is),*)' '
c     write (ifsec(is),6500) is,nbloch,deltaw,ef
c     do      ip = 1,nq
c     call wsec2 (ifsec(is),is,q(1,ip),itq(1:ntq),
c     i              eqx(1:ntq,ip,is),deltaw,nbloch,
c     i              dreal(zsec(:,:,ip)), dimag(zsec(:,:,ip)),
c     d              3,ntq,nband)
c     end do
c
c-------------------------
          ifoutsec=ifsec(is)
          do ip = 1,nq
            do i  = 1,ntq
              zseciip = 0d0
c     if(imgonly) zseciip(iwini:iwend) = zsec(0, i, ip)
              zseciip(iwini:iwend) = zsec(iwini:iwend, i, ip)

              dscdw1  = dreal( zseciip(1) - zseciip(0)  ) /deltaw
              zfac1   = 1d0/(1d0-dscdw1)
              dscdw2  = dreal( zseciip(0) - zseciip(-1) ) /deltaw
              zfac2   = 1d0/(1d0-dscdw2)
              dscdw   = dreal( zseciip(1) - zseciip(-1) ) /(2d0*deltaw)
              zfac    = 1d0/(1d0-dscdw)
              write(6,"(' j iq isp=' i3,i4,i2,'  q=',3f8.4,
     &                 '  eig=',f8.4,'  Re(Sc) =',3f8.4,'  Img(Sc) =',3f8.4,
     &                 '  zfac=',f8.4,'  zfac1&2 =',2f8.4)") 
     &                 itq(i),ip,is, q(1:3,ip), eqx(i,ip,is), 
     &                 hartree*dreal(zseciip(-1:1)), 
     &                 hartree*dimag(zseciip(-1:1)),zfac,zfac1,zfac2
              write(ifoutsec,"(3i5,3d24.16,3x,d24.16,3x,3d24.16,
     &                 3x,3d24.16,3x,3d24.16)")
     &                 itq(i),ip,is, q(1:3,ip), eqx(i,ip,is), 
     &                 hartree*dreal(zseciip(-1:1)), 
     &                 hartree*dimag(zseciip(-1:1)),zfac,zfac1,zfac2
ccccccccccccccccccccc
              if(cohtest) then
                write(671,"( i3,i4,i2, 3f8.4,
     &                    '  eig=',f8.4,'  coh =',f10.4)") 
     &                    itq(i),ip,is, q(1:3,ip), eqx(i,ip,is), 
     &                    hartree*coh(i,ip)
              endif
ccccccccccccccccccccc
            end do
          end do
c-------------------------
        elseif(ixc==4) then    !spectrum mode
c     ifoutsec=ifsecomg(is)
          ifoutsec=9300
          if(is==1) sss='.UP'
          if(is==2) sss='.DN'
          open(ifoutsec,file='SEComg'//sss)
          do ip = 1,nq
            do i  = 1,ntq
              write(6,"(' --- j iq isp=' i3,i4,i2,' q=',3f8.4,' eig=',f8.4)") 
     &                 itq(i),ip,is, q(1:3,ip), eqx(i,ip,is)

c...  for previous version
c     open(ifoutsec,file='SEComg'//'_ib'//charnum3(itq(i))//'_ip'//charnum3(ip)//sss)
c     sumimg=0d0
              do iw = iwini,iwend
                if(mod(iw,50)==1) 
     &                    write(6,"(' omega-ef=',f9.4,'  Sc=',2f9.4)") 
     &                    (omega(i,iw)-ef)*rydberg(), hartree*zsec(iw,i,ip)
c     sumimg = sumimg + abs(dimag(zsec(iw,i,ip))*dwplot)

c...  for previous version
c     write(ifoutsec,"(4i5,3f10.6,3x,f10.6,2x,f16.8,x,3f16.8)")
c     &      iw,itq(i),ip,is, q(1:3,ip), eqx(i,ip,is),
c     &      (omega(i,iw)-ef)*rydberg(),  hartree*zsec(iw,i,ip) !,sumimg
                write(ifoutsec,"(4i5,3f10.6,3x,f10.6,2x,f16.8,x,3f16.8)")
     &                    iw,itq(i),ip,is, q(1:3,ip),  eqx(i,ip,is),
     &                    (omega(i,iw)-ef)*rydberg(),  hartree*zsec(iw,i,ip) !,sumimg
              end do
c     write(ifoutsec,*)
              write(6,*)
            end do
          end do
          close(ifoutsec)
        endif
 2000 continue                  !end of spin-loop

      if (tote) then
        write(ifexx,*) '======================================='
        if(ixc==5) write(ifexx,*) 
     &        '  Exchange energy valence-valence    Exx (eV)   '
        if(ixc==6) write(ifexx,*) 
     &        '  Exchange energy core-valence       Exx (eV)   '
c     write(ifexx,*) '    Exchange energy     Exx (eV)       '
        write(ifexx,*) '======================================='

        write (ifexx,*)' *** '
        if (nspinmx == 1) exx = exx * 2.d0
        write(ifexx,*)exx*hartree

        call winfo2(ifexx,nspin,nq,ntq,is,nbloch
     &        ,ngpn1,ngcn1,nqbz,nqibz,ef,ef2,deltaw,alat,esmr,esmr2)

        if (legas) then
          pi         = 4.d0*datan(1.d0)
          efz=(ntot*3*pi**2/voltot)**(2d0/3d0) ! ef is calculated from ntot.
          qfermi= dsqrt(efz)
          alpha = (9*pi/4d0)**(1d0/3d0)
          rs    = alpha/qfermi
          write (ifexx,*)' --- electron gas ---'
          write (ifexx,*)' density parameter rs= ', rs
          write (ifexx,*)' kf= ',qfermi
          write (ifexx,*)' *** exact exchange'
          exxelgas = -0.4582 * hartree /rs * ntot
          write (ifexx,*)exxelgas
        endif
      endif

C...  close files
      isx = iclose ('wc.d')
      isx = iclose ('wci.d')
      isx = iclose ('hbe.d')


C     --- EXspectrum ----------------------------------------------------
c     This section is similar with efsimplef.f
      if(sum(ifexsp(1:nspin))/=0) then
        do is = 1,nspin
          print *,' --- Goto ExSpectrum section --- is=',is
          rewind (ifexsp(is))
          itmx = 0
          do
            read(ifexsp(is),*,end=1215)ipex,itpex,itex,qqex(1:3), eex,exsp
            if(itex>itmx) itmx=itex
          enddo
 1215     continue
          nspexmx = itmx*(nqbz+nq0i*ngrp) !Get marimum value of the number of the ex spectrum
c
          allocate( eex1(nspexmx,ntq,nq), exsp1(nspexmx,ntq,nq),
     &           nspex(ntq,nq) ,
     &           itex1(nspexmx,ntq,nq),
     &           qqex1(3,nspexmx,ntq,nq) )
          print *,' nspexmx =',nspexmx
c
          rewind (ifexsp(is))
          nspex = 0
          do
            read(ifexsp(is),*,end=1216) ipex,itpex,itex,qqex(1:3),eex,exsp
            nspex(itpex,ipex) = nspex(itpex,ipex)+1
            iex  = nspex(itpex,ipex)
            eex1  (iex,itpex,ipex) = eex
            exsp1 (iex,itpex,ipex) = exsp
            itex1 (iex,itpex,ipex) = itex
            qqex1(:,iex,itpex,ipex)= qqex
          enddo
 1216     continue            !Get eex1(1:nspex) exsp1(1:nspex) for itp ip.
          print *,' nspex(1 1)=',nspex(1,1)
c
          do ipex = 1,nq
            do itpex=1,ntq
              write(6,*)' is itq ip =',is,itq,ip
              nnex = nspex(itpex,ipex)
              allocate( ieord(1:nnex) )
              call sortea( eex1(1:nnex,itpex,ipex),ieord, nnex,isig)
              eex1 (1:nnex,itpex,ipex)  = eex1  (ieord(1:nnex),itpex,ipex)
              exsp1 (1:nnex,itpex,ipex) = exsp1 (ieord(1:nnex),itpex,ipex)
              itex1 (1:nnex,itpex,ipex) = itex1 (ieord(1:nnex),itpex,ipex)
              qqex1(:,1:nnex,itpex,ipex)= qqex1 (:,ieord(1:nnex),itpex,ipex)

              filenameex = 'EXSP'//charnum3(ipex)//charnum3(itpex)
     &                 //'.'//char(48+is)
              ifexspx=4111
              open(ifexspx,file=filenameex)

              filenameex = 'EXSS'//charnum3(ipex)//charnum3(itpex)
     &                 //'.'//char(48+is)
              ifexspxx=4112
              open(ifexspxx,file=filenameex)

              do i=1,nnex
                write(ifexspx, "(2d14.6, i4, 3f14.6)")
     &                    eex1  (i,itpex,ipex), exsp1 (i,itpex,ipex),
     &                    itex1 (i,itpex,ipex), qqex1 (1:3,i,itpex,ipex)
              enddo
c
              eee  =-1d99
              exwgt= 0d0
              do i=1,nnex
                if(eex1(i,itpex,ipex) > eee+1d-4 .or. i==nnex) then
                  if(i/=1) write(ifexspxx, "(2d23.15)")
     &                       eee, exwgt*hartree
                  eee  = eex1(i,itpex,ipex)
                  exwgt= exsp1 (i,itpex,ipex)
                else
                  exwgt= exwgt + exsp1 (i,itpex,ipex)
                endif
              enddo
c
              deallocate( ieord )
              close(ifexspx)
              close(ifexspxx)
            enddo
          enddo
          deallocate( eex1, exsp1, nspex, itex1, qqex1 )
        enddo
        print *,' End of ExSpectrum section ---'
      endif
c-------------------------------------------------------------------------------
      call cputid(0)
c     if(ixc==1.and.sergeys)stop " OK! hsfp0: Sergey's Exchange mode"
c     if(ixc==2.and.sergeys)stop " OK! hsfp0: Sergey's Correlation mode"
      if(ixc==1 ) stop ' OK! hsfp0: Exchange mode'
      if(ixc==2 ) stop ' OK! hsfp0: Correlation mode'
      if(ixc==3 ) stop ' OK! hsfp0: Core-exchange mode'
      if(ixc==4 ) stop ' OK! hsfp0: spectrum mode'
      if(ixc==5 ) stop ' OK! hsfp0: Exx mode val-val  See TEEXXvv'
      if(ixc==6 ) stop ' OK! hsfp0: Exx mode core-val See TEEXXcv'
      end

