      program hsfp0_v2_sc
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
c     !mode= 4: plot spectrum function ---See manual
c
c     iSigma_en parameter which determines approximation for  self-energy,
c     stored in GWIN_V2.
c     iSigma_en==0 SE_nn'(ef)+image integr:delta_nn'(SE_nn(e_n)-SE_nn(ef))
c     iSigma_en==1 SE_nn'(ef)+delta_nn'(SE_nn(e_n)-SE_nn(ef))
c     iSigma_en==2 SE_nn'((e_n+e_n')/2)
c     iSigma_en==3 (SE_nn'(e_n)+SE_nn'(e_n'))/2
c     iSigma_en==5 delta_nn' SE_nn(e_n)
c     output file should contain hermitean part of SE for energies to be real
c     (for example, hermitean conjunction of SE_nn'(e_n) means SE_n'n(e_n')^* )
c
c     Apr 2002 takao kotani. multiple argumentation wave per l.
c     This hsfp0 is build from hsec10.f by F.Aryasetiawan.
c------------------------------------------------------------
      use keyvalue
      use m_readqg
      use m_readeigen,only: init_readeigen,init_readeigen2,readeval,lowesteval
      use m_read_bzdata,ngrp2=>ngrp
      use m_genallcf_v3
      implicit none
c     real(8),parameter ::
c     &  ua    = 1d0    ! constant in w(0)exp(-ua^2*w'^2) to take care of peak around w'=0
c------------------------------------
c     test switches to calculate the self-energy based on an another separation of \Sigma.
!     \Sigma = \Sigma_{sx} + \Sigma_{coh} + \Sigma_{img axis} + \Sigma_{pole} by Hedin PR(1965)A785
!     I found COH term has inevitably poor accuracy.
      logical ::tetra, tetra_hsfp0,         
     &     screen = .false.,    ! \Sigma_{sx} for mode 1 and
      ! \Sigma_{img axis} + \Sigma_{pole} for mode 2
     &     cohtest= .false.     ! \Sigma_{coh}. mode swich is not required.
c
c     &  , tetra  = .false.  ! test switch for tetrahedron method test.
c     ! tetra=T is only effective for exchange=T case.
c     ! Tetrahedron mehod for correlation is a bit
      ! difficult and I gave up.
c     ! If you want to calculate with tetra=T for exchange, you
      ! have to uncomment tetra related part in
      ! sxcf.f, and a part calling sxcf in this routine. Note wtet wtetef!
      ! They sometimes cause array destruction if you run tetra=T without comment them.
c------------------------------------
      real(8) :: shtw
c     integer(4) :: mxclass,ngnmax,mbytes,mwords,iwksize,
c     &   natom,nclass,ipos,ngrp,igrp,
c     &   iinvg,
c     o   nspin,nl,nn,nnv,nnc,
c     o   inindx,inindxv,inindxc,iiclass,
c     d   nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc,
c     o   iz,
c     o   iil,iin,iim,iilnm,i_mnl,
c     o   iilv,iinv,iimv,iilnmv,i_mnlv,
c     o   iilc,iinc,iimc,iilnmc,i_mnlc,
c     o   incwf,iecore,ikonf,iicore,incore,nctot,
c     o   imagw,niw,nw,ifreq,
      integer(4)::
     &     ixc,iopen,ifhbed, nprecb,mrecb,mrece,nlmtot,nqbzt, nband,
     &     ibas,ibasx,ngpmx,nxx,ngcmx,nbloch,ifqpnt,ifwd,
     &     nprecx,mrecl,nblochpmx2,nwt,niwt, nqnum,mdimx,nblochpmx,
     &     ifrcw,ifrcwi,  noccxv,maxocc,noccx,ifvcfpout,iqall,iaf,ntq,
     &     i,k,nspinmx, nq,is,ip,iq,idxk,ifoutsex,iclose,nq0i,ig,
     &     mxkp,nqibzxx,ntet,nene,iqi, ix,iw,
     &     nlnx4,niwx,irot,invr,invrot,ivsum, ifoutsec,ntqx,
c
c     &   ifrb(2),ifcb(2),ifrhb(2),ifchb(2)
c     &    ifev(2),
     &     ifsec(2)             !ifcphi
     &     ,ifxc(2),ifsex(2), ifphiv(2),ifphic(2),ifec,ifexsp(2),
     &     ifsex2(2),ifsec2(2), !out S_nn'
     &     ifsecomg(2),ndble=8
c
c     real(8) :: alat,ef,diw,dw,delta,pi,tpia,vol,voltot,rs,alpha,
      real(8) :: pi,tpia,vol,voltot,rs,alpha,
     &     qfermi,efx,valn,efnew,edummy,efz,qm,xsex,egex,
     &     zfac1,zfac2,dscdw1,dscdw2,dscdw,zfac
c
c     &   lowesteval !defined in readeigen

ccc   qbas(3,3),ginv(3,3) plat(3,3),

      logical lqall,laf
c     character*120 symgrp

c     class parameters
c     parameter (mxclass=100)
c     character*6 clabl(mxclass)
c     symmetry group
c     parameter (ngnmax=10)
c     real(8) :: gen(9,ngnmax)
c
      integer(4),allocatable :: itq(:)
      real(8),allocatable    :: q(:,:)

c     takao
      integer(4),allocatable :: ngvecpB(:,:,:), !ngveccB(:,:,:),
     &     ngvecp(:,:), ngvecc(:,:),iqib(:), !ngcni(:),
     &     kount(:,:), nx(:,:),nblocha(:),lx(:) !,ngveccBr(:,:,:)
      real(8),allocatable:: vxcfp(:,:,:),
     &     wqt(:), wgt0(:,:),q0i(:,:),
     &     ppbrd (:,:,:,:,:,:,:),cgr(:,:,:,:),eqt(:),
     &     ppbrdx(:,:,:,:,:,:,:), !aaa(:,:), symope(:,:,:), !qibz(:,:),
     &     ppb(:), eq(:),       !,pdb(:),dpb(:),ddb(:)
     &     eqx(:,:,:),eqx0(:,:,:),ekc(:),coh(:,:)
      complex(8),allocatable:: geigB(:,:,:,:) ,zsec(:,:,:)
c
      logical :: exchange, legas
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
     &     ,iene(:,:,:),ibzx(:) !,nstar(:)
c     real(8) :: qbasmc(3,3)

c     worksize in megabytes (1 word = 4 bytes)
c     parameter (mbytes=60)
c     parameter (mwords=mbytes/4)
c     parameter (iwksize=mwords * 1000 *1000)
c     integer w
c     common /w/ w(iwksize)

      integer(4) ::ib,iqx,igp,iii,ivsumxxx,isx,iflegas, iqpntnum
c
      real(8),allocatable   :: eex1(:,:,:),exsp1(:,:,:),qqex1(:,:,:,:)
      integer(4),allocatable:: nspex(:,:),ieord(:),itex1(:,:,:)
      real(8)    :: qqex(1:3), eex,exsp,eee, exwgt,deltax0
      integer(4) :: itmx,ipex,itpex,itex,nspexmx,nnex,isig,iex,ifexspx
     &     ,ifexspxx ,ifefsm, nq0ix,ifemesh,nz
      character(3)  :: charnum3
      character(12) :: filenameex
      logical :: exspwrite=.false.
      character*8 xt

      integer(4) :: iSigma_en,ifinin ,idummy !iwini,iwend       !sf 21May02

      real(8),allocatable:: omega(:) !sf 21May02
c     real(8) ::  omegamax,dwplot,omegamaxin
c     logical :: sergeys     !sf 21may02

c     real(8)   :: ebmx
c     integer(4):: nbmx
      real(8)   :: ebmx(2)
      integer(4):: nbmx(2)

      real(8):: volwgt

      integer(4)::nwin, incwfin
      real(8)::efin
      integer(4),allocatable::imdim(:)
      real(8),allocatable::freqx(:),freqw(:),wwx(:) !,expa(:)

      integer(4):: bzcase,  ngpn1,mrecg,ngcn1
      real(8)   :: wgtq0p,quu(3)

      character(2):: soflag
      integer(4):: ifianf


      logical ::smbasis
      integer(4):: ifpomat,nkpo,nnmx,nomx,ikpo,nn_,no
      real(8):: q_r(3)
      real(8),allocatable:: qrr(:,:)
      integer(4),allocatable:: nnr(:),nor(:)
      complex(8),allocatable:: pomatr(:,:,:),pomat(:,:)

      logical allq0i            !S.F.Jan06

      integer(4):: nw_i
      logical:: exonly
      real(8):: wex

!! newaniso mode
      real(8),allocatable:: vcousq(:),dmlx(:,:),epinvq0i(:,:),wklm(:),vcoud(:)
      complex(8),allocatable:: zcousq(:,:)
      logical:: newaniso
      integer:: ifvcoud,lxklm,ifidmlx

c---------------------------------------
      hartree=2d0*rydberg()

c     mode switch. --------------
      if(cohtest) then
        screen = .true.
        ixc = 2; nz=0
        open(671,file='COH')
      else
        write(6,*) ' --- Choose modes below ------------'
        write(6,*) '  Sx(1) Sc(2) ScoreX(3) '
        write(6,*) '  [option --- (+ QPNT.{number} ?)] '
        write(6,*) ' --- Put number above ! ------------'
        call readin5(ixc,nz,idummy)
      endif
      call headver('hsfp0_sc',ixc)
      write(6,*) ' ixc nz=',ixc, nz
      if(ixc==0) stop ' --- ixc=0 --- Choose computational mode!'

c---  readin BZDATA. See gwsrc/rwbzdata.f
c--------readin data set when you call read_BZDATA ---------------
c     integer(4)::ngrp,nqbz,nqibz,nqbzw,nteti,ntetf
c     !  &   ,n_index_qbz
c     integer(4):: n1,n2,n3
c     real(8):: qbas(3,3),ginv(3,3),qbasmc(3,3)
c     real(8),allocatable:: qbz(:,:),wbz(:),qibz(:,:)
c     &    ,wibz(:),qbzw(:,:)
c     integer(4),allocatable:: idtetf(:,:),ib1bz(:),idteti(:,:)
c     &    ,nstar(:),irk(:,:),nstbz(:)           !,index_qbz(:,:,:)
c-----------------------------------------------------------------
      call read_BZDATA()
      print *,' nqbz =',nqbz
c     print *,  qbz
      print *,' nqibz ngrp=',nqibz,ngrp2
c     print *,' irk=',irk
c     print *,' #### idtetf: ####'
c     print *, idtetf

c
c     set up work array
c     call wkinit (iwksize)
      call pshprt(60)

C---  readin GWIN and LMTO, then allocate and set datas.
      nwin   = 0                !Readin nw from NW file
      efin=-999d0               !not readin EFERMI
      if(ixc==3) then; incwfin= -2 !core exchange mode
      else           ; incwfin= -1 !use 7th colmn for core at the end section of GWIN
      endif
      call genallcf_v3(nwin,efin,incwfin) !in module m_genallcf_v3
      if(ngrp/= ngrp2) stop 'ngrp inconsistent: BZDATA and LMTO GWIN_V2'
c---  These are allocated and setted by genallcf_v3
c     integer(4)::  nclass,natom,nspin,nl,nn,nnv,nnc, ngrp,
c     o  nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, nctot,niw,nw
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
c     o  plat(:,:),pos(:,:),z(:),  ecore(:,:), freq(:), symgg(:,:,:) ! symgg=w(igrp)
c-----------------------------------------------------------------------
c     call genallcf_v2x(
c     c> structure
c     o                   plat,alat,natom,nclass,ipos,
c     c> symmetry
c     o                   symgrp,gen,ngnmax,ngrp,igrp,
c
c     o                   iinvg,ef,
c
c>>   l,n and dimensions
c     o                   clabl,nspin,nl,nn,nnv,nnc,
c     o                   inindx,inindxv,inindxc,iiclass,
c     d                   nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc,
c     c z value
c     o                   iz,
c>>   l,n,m indices for Phi
c     o                   iil,iin,iim,   iilnm  ,i_mnl,
c     o                   iilv,iinv,iimv,iilnmv ,i_mnlv,
c     o                   iilc,iinc,iimc,iilnmc ,i_mnlc,
c>>   core
c     o                   incwf,iecore,ikonf,iicore,incore,nctot,
c     > frequency
c     o                   imagw,niw,diw,nw,dw,delta,deltaw,esmr,ifreq )


c-----------reading self-energy mode parameter from file 'GWIN_V2'
      call readd_iSigma_en(ifinin,iSigma_en) !sf..21May02

c     allocate(ecore(nctot,nspin)) !core energies
c     do is = 1,nspin
c     if (nctot > 0) call catch1 (w(iecore),is,nctot,2,ecore(:,is)) !core energies
c     enddo

c---  Get maximums takao 18June03
      call getnemx8(nbmx,ebmx)
      write(6,"('  nbmx ebmx from GWinput=',2i8,2d13.5)") nbmx,ebmx
!     See sxcf_fal2z--------
!     nbmx1 ebmx1: to set how many bands of <i|sigma|j>  do you calculate.
!     nbmx2 ebmx2: to restrict num of bands of G to calculate G \times W

c-------------------------------------------------------------------
c     if (nclass > mxclass) stop ' hsfp0: increase mxclass'
c!!!! WE ASSUME iclass(iatom)= iatom !!!!!!!!!!!!!!!!!!!!!!!!!
      if (nclass /= natom ) stop ' hsfp0: nclass /= natom ' ! We assume nclass = natom.
      print *,' hsfp0: end of genallcf_v3'
c
      call pshprt(30)
      pi   = 4d0*datan(1d0)
      tpia = 2d0*pi/alat

      shtw = 0d0
      if(esmr<1d-5) shtw=0.01d0 ! Ferdi's shift to avoid resonance effect(maybe)
c
      call dinv33(plat,1,xxx,vol)
      voltot = dabs(vol)*(alat**3)

      tetra= tetra_hsfp0()

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
        if(tetra) stop 'legas You have to give ef of  tetrahedron'
      endif

c
      ifexsp=0
      if(ixc==1) then
        exchange=.true.
        write(6,*) ' --- Exchange mode --- '
        ifxc(1)  = iopen('XCU'//xt(nz),1,-1,0)
        ifsex(1) = iopen('SEXU'//xt(nz),1,-1,0)
        ifsex2(1)= iopen('SEX2U',0,-1,0) !out SEX_nn'
        if (nspin == 2) then
          ifxc(2)  = iopen('XCD'//xt(nz),1,-1,0)
          ifsex(2) = iopen('SEXD'//xt(nz),1,-1,0)
          ifsex2(2)= iopen('SEX2D',0,-1,0) !out SEX_nn'
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
        write(6,*) ' --- Correlation mode --- '
        if(cohtest) write(6,*) ' COH calculation mode. Results in COH'
        ifsec(1)   = iopen('SECU'//xt(nz),1,-1,0) ! output files
        ifsec2(1)= iopen('SEC2U',0,-1,0) !out SEC_nn'
        if (nspin == 2) 
     .        ifsec(2)   = iopen('SECD'//xt(nz),1,-1,0)
        ifsec2(2)= iopen('SEC2D',0,-1,0) !out SEC_nn'
      elseif(ixc==3) then
        exchange=.true.
        esmr=0d0
        write(6,*) ' --- CORE Exchange mode --- '
        ifsex(1)   = iopen('SEXcoreU'//xt(nz),1,-1,0)
        ifsex2(1)= iopen('SEXcore2U',0,-1,0) !out SEXcore_nn'
        if (nspin == 2) then
          ifsex(2)   = iopen('SEXcoreD'//xt(nz),1,-1,0)
          ifsex2(2)= iopen('SEXcore2D',0,-1,0) !out SEXcore_nn'
        endif
ccccccccccccccccccccccccccccccccccccccccccfaleev 21May02 do not use ixc==4
c     elseif(ixc==4) then
c     write(6,*) ' --- Spectrum function Sigma(\omega) mode --- '
c     exchange=.false.
c     ifsecomg(1) = iopen('SEComgU'//xt(nz),1,-1,0) ! output files
c     if (nspin == 2)
c     .  ifsecomg(2) = iopen('SEComgD'//xt(nz),1,-1,0)
      else
        stop ' hsfp0: Need input (std input) 1(Sx) 2(Sc) or 3(ScoreX)!'
      endif

c---  Neglect core is NoCore exists ------------
c     inquire(file='NoCore',exist=nocore)
c     if(nocore) nctot=0

c---
      write(6, *) ' --- computational conditions --- '
      write(6,'("    deltaw  =",f13.6)') deltaw
c     write(6,'("    ua      =",f13.6)') ua
      write(6,'("    esmr    =",f13.6)') esmr
      write(6,'("    alat voltot =",2f13.6)') alat, voltot
      write(6,'("    niw nw dw   =",2i3,f13.6)') niw,nw,dw

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
c     call mptauof(w(igrp),ngrp,plat,natom,pos,iclasst
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
c     nblochpmx = nbloch + ngcmx !not use this oct2005

c     allocate(ngcni(nqibz)) !, ngveccB(3,ngcmx,nqibz), ngveccBr(3,ngcmx,nqibz))
c     geigB(ngpmx,nband,nqbz,nspin),ngpn(nqbz),ngvecpB(3,ngpmx,nqbz),
c     &   )  ! in IBZ
c
c     call rdpp_pln(ngpmx,ngcmx, qibz,nqibz, qbz,nqbz,nband,nspin,
c     o      ngpn,geigB,ngvecpB,ngcni,ngveccB)
c
c     do iq = 1,nqibz
c     call readqg('QGcou',qibz(1:3,iq),ginv,  quu,ngcni(iq)), ngveccB(1,1,iq))
c     write(6,"('--From QGcou  qibz quu ngc=',3f9.4,'  ',3f9.4,i5)")
c     &       qibz(1:3,iq),quu,ngcni(iq)
c     enddo
c     for info
      allocate(ngvecp(3,ngpmx),ngvecc(3,ngcmx))
      call readqg('QGpsi',qibz(1:3,1),ginv, quu,ngpn1, ngvecp)
      call readqg('QGcou',qibz(1:3,1),ginv, quu,ngcn1, ngvecc)
      deallocate(ngvecp,ngvecc)
c     do i=1,nqibz
c     write (6,*) i, ngcni(i)
c     enddo
c     stop 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx '

      print *, ' end of read QGcou'

cccccccccccccccccccccccccccccccccccccccccccccccccccc
c     iqx=1
c     do ib=1,nband
c     write(6,'("  iband iqx sumgeigB=",2i3,12d12.3)')
c     &   ib,iqx, sum(geigB(1:ngpn1,ib,iqx,1))
c     enddo
c     stop 'test end'
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

c------------
c     open files
c------------
c     > input files
c     ifqpnt    = iopen('QPNT'//xt(nz),1,0,0)

c     screened exchange
      if(.not.exchange.or.(exchange.and.screen)) then
        !screened exchange
c     direct access files WVR and WVI which include W-V.
        ifwd      = iopen('WV.d',1,-1,0)
        read (ifwd,*) nprecx,mrecl,nblochpmx,nwt,niwt, nqnum,nw_i
        write(6,"(' Readin WV.d =', 10i5)") 
     &        nprecx,mrecl,nblochpmx,nwt,niwt, nqnum, nw_i
c     nblochpmx from WV.d oct2005
        call checkeq(nprecx,ndble)
        nw=nwt
        if (niwt /= niw) stop 'hsfp0: wrong niw'
        ifrcw     = iopen('WVR',0,-1,mrecl)
        if(.not.exchange) ifrcwi = iopen('WVI',0,-1,mrecl)
      endif

c     eigen functions
c     ifev(1)   = iopen('EVU', 0,0,mrece)
c     if (nspin == 2) ifev(2) = iopen('EVD', 0,0,mrece)

c-------------------------
      if(tetra) goto 201
c-------------------------

      if(ixc==3) then
c     ef = lowesteb(ifev,nspin, nband,nqbz)-1d-3
c     ef = lowesteigen(nspin, nband,qbz,nqbz)-1d-3
        ef = lowesteval() -1d-3 !lowesteigen(nspin,nband,qbz,nqbz) - 1d-3 !lowesteb was

c-------------------
!     exonly 03June2006
        call getkeyvalue("GWinput","EXonly",wex,default=0d0)
        if(wex==0d0) then
          exonly=.false.
        else
          exonly=.true.
          print *,' exonly=T ecore shift: ecore---> ecore-100'
          ecore = ecore-100.0
        endif

        if(maxval(ecore(:,1:nspin))>ef) then
          print *,' ef nspin=',ef,nspin,nctot
          do is=1,nspin
            print *,' maxval( ecore) nctot=', is,nctot
            do ix=1,nctot
              write(6,"(i4,d13.5)") ix,ecore(ix, is)
            enddo
          enddo
          stop 'hsfp0 ixc=3:  ecore>evalence '
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
     i        qbz,nqbz          ! index_qbz, n_index_qbz,
     o        ,efnew)

c     call efsimplef(ifev,nspin,w(iwibz),w(iindxk),w(iqibz)
c     i       ,n1,n2,n3,qbas,ginv, nband,nqibz
c     i       ,w(ikonf),w(iz),nl,natom,w(iiclass),nclass
c     i       ,valn, legas, esmr  !!! valn is input for legas=T, output otherwise.
c     o       ,efnew)
c
c     print *,' end of efsimple'
        ef = efnew
c     - check total ele number -------
        ntot  = nocctotg2(nspin, ef,esmr, qbz,wbz, nband,nqbz)
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

cccccccccccccccccccccccccccccccccccccc
c     ef=-0.22
c     print *, " !!! test fermi energy=",ef
cccccccccccccccccccccccccccccccccccccc

c--
c     noccxv = maxocc (ifev,nspin, ef+0.5d0*esmr, nband,nqbz)  ! maximum no. of occupied valence states
!     maxocc seems to give (the maxmum number of occ + 1).
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

c-------------------------
 201  continue
c-------------------------

c     ifcphi  = iopen('CPHI',0,0,mrecb)

      call init_readeigen2(mrecb,nlmto,mrecg) !initialize m_readeigen

c     the coulomb matrix for exchange.
      if(exchange) ifvcfpout = iopen('VCCFP',0,-1,0)


c     read q-points and states
cccccccccccccccccccccccccccccccccfaleev 21May02 do not use ixc==4
c     if(ixc==4) then
c     call readxx(ifqpnt) !skip to ***** for q point for spectrum function.
c     omegamaxin = 1d70
c     read (ifqpnt,*,err=2038,end=2038) dwplot,omegamaxin
c     2038   continue
c     omegamax = (nw-2)*2d0*dw - 0.5d0*esmr !This is in Ry.
c     if(omegamaxin < omegamax) then
c     write(6,*) ' --- Use readin dwplot and omegamaxin from QPNT.'
c     omegamax = omegamaxin
c     else
c     2          write(6,*) ' --- omegamax by hx0fp >omegamaxin'
c     write(6,*) ' --- Use dwplot from QPNT and omegamax by hx0fp.'
c     endif
c     if( omegamax <0) stop 'hsfp0 :strange omegamax <0 '
c     iwini =  -int( omegamax / dwplot )
c     iwend =   int( omegamax/  dwplot )
c     print *,' Readin dwplot from QPNT: iwini iwend',iwini,iwend
c     2        dwidth = nw*dw
c     else
ccccccccccccccccccccccccccccccfaleev 21May02 do not use iwini, iwend
c     iwini= -1
c     iwend= 1
c     endif
c     ---
c     lqall      = .false.
c     laf        = .false.
c     call readx   (ifqpnt,10)
c     read (ifqpnt,*) iqall,iaf
c     if (iqall == 1) lqall = .true.
c     if (iaf   == 1)   laf = .true.
c     call readx   (ifqpnt,100)
!     states
c     read (ifqpnt,*) ntq
c     if (ntq<0) then ! all states !takao May 2003
      ntq =  min(nbmx(1),nband) !takao 19June2003 nband
      allocate (itq(ntq))
      do i = 1, ntq
        itq(i) = i
      enddo
c     else              ! number states from GWIN0
c     allocate( itq(ntq) )
c     read (ifqpnt,*) (itq(i),i=1,ntq)
c     endif

!     q-points
c     if (lqall) then !all q-points case
      if(bzcase()==1)then
        nq = nqibz
        allocate(q(3,nq))
        call dcopy   (3*nqibz,qibz,1,q,1)
      else
c     nq = nqibz+nqibz_r
c     allocate(q(3,nq))
c     do iq= 1,nqibz_r
c     q(:,iq) = qibz_r(:,iq)
c     enddo
c     do iq= 1,nqibz
c     q(:,iq+nqibz_r)=qibz(:,iq)
c     enddo
        nq = nqibz_r
        allocate(q(3,nq))
        do iq= 1,nqibz_r
          q(:,iq)=qibz_r(:,iq)
        enddo
      endif
      do iq=1,nq
        write(6,'(" Target iq q=",i5,3f9.4)')iq,q(:,iq)
      enddo
c     else
c     call readx   (ifqpnt,100)
c     read (ifqpnt,*) nq
c     allocate(q(3,nq))
c     do       k = 1,nq
c     read (ifqpnt,*) i,q(1,k),q(2,k),q(3,k)
c     write(6,'(i3,3f13.6)') i,q(1,k),q(2,k),q(3,k)
c     enddo
c     endif
cccc  if (ntq /= nband) stop 'hsfp0: SC GW: ntq /= nband'
cccc  if (nq /= nqibz)  stop 'hsfp0: SC GW: nq /= nqibz'
cccc  last statement should be commented if  less then nqibz
cccc  momentums in output SEU(D) wanted
c
c
      nspinmx = nspin
      inquire(file='ANFcond',EXIST=laf)
      if(laf) then
        ifianf = 211
        open(ifianf,file='ANFcond')
        read(ifianf,*) soflag
        if(soflag=='SO') then
          print *,' ANFcond found. but SO case is not implimented. Not using ANF symm'
          goto 889            !SO case is not implimented
        endif
        print *, 'Find ANFcond--- up only'
        nspinmx =1
      endif
 889  continue

c------------'
      allocate(omega(ntq))      !sf 21May02
cccccccccccccccccccccccccccccccccfaleev 21May02 do not use ixc==4
c     if(ixc==4) then
c     do iw = iwini,iwend
c     omega(1:ntq,iw) =  dwplot* iw + ef
c     enddo
c     endif

c     --- read LDA eigenvalues
c     ntp0=ntq
      allocate(eqx(ntq,nq,nspin),eqx0(ntq,nq,nspin),eqt(nband))
      do      is = 1,nspin
        do      ip = 1,nq
c     iq       = idxk (q(1,ip),qbz,nqbz)
c     call rwdd1   (ifev(is), iq, nband, eqt) !direct access read b,hb and e(q,t)
          call readeval(q(1,ip),is,eqt)
          eqx0(1:ntq,ip,is) = eqt(itq(1:ntq))
          eqx (1:ntq,ip,is) = rydberg()*(eqt(itq(1:ntq))- ef)
        enddo
      enddo
      deallocate(eqt)
c     --- info
      call winfo(6,nspin,nq,ntq,is,nbloch
     &     ,ngpn1,ngcn1,nqbz,nqibz,ef,deltaw,alat,esmr)

c-------------------------
c     LDA exchange-correlation
c-------------------------
c     - takao Vxc is readin from VXCFP
      if(ixc==1) then
        allocate(  vxcfp(ntq,nq,nspin) )
        call rsexx(nspin,itq,q,ntq,nq, ginv, vxcfp) !add ginv july2011
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
              if(eqx(i,ip,is) <1d20.and.vxcfp(i,ip,is)/=0d0) then  !takao june2009. See lmf2gw (evl_d=1d20; in Ry.. but eqx is in eV. no problem for inequality).
                write(6,"(' j iq isp=' i3,i4,i2,'  q=',3f8.4,
     &                 '  eig=',f10.4,'  Sxc(LDA)=',f10.4)") 
     &                 itq(i),ip,is, q(1:3,ip), eqx(i,ip,is),
     &                 vxcfp(i,ip,is)
              endif
            end do
          end do
c     end of spin-loop
          if(is==1) isx = iclose('XCU'//xt(nz))
          if(is==2) isx = iclose('XCD'//xt(nz))
        end do
        deallocate(vxcfp)
      endif
c


c---  takao not modifed weights -----------------
ccccccccccccccccccccccccccccccccccccc
c     call defi    (ikfbz,ngrp)
c     call modwt02 (w(iwbz),w(iqbz),w(iqibz),qbas,w(iirk),w(instar),
c     w              w(ikfbz),
c     d              nqibz,nqbz,ngrp)
ccccccccccccccccccccccccccccccccccccc


c     q near zero
      print *, 'reading QOP'
      open (101,file='Q0P')
      read (101,"(i5)") nq0i
      if(.not.exchange) call checkeq(nqibz+nq0i-1, nqnum)
      write(6,*) ' *** nqibz nq0i_total=', nqibz,nq0i
      allocate( wqt(1:nq0i),q0i(1:3,1:nq0i) )
c     read (101,"(d24.16,3x, 3d24.16)" )( wqt(i),q0i(1:3,i),i=1,nq0i)
      nq0ix = nq0i

      do i=1,nq0i
        read (101,* ) wqt(i),q0i(1:3,i)
        if(wqt(i)==0d0 ) nq0ix = i-1
      enddo
      nq0i = nq0ix              ! New nq0i July 2001
      write(6,*) ' Used k number in Q0P =', nq0i
c
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

      if (nq0i/=0 ) write(6,*) ' *** tot num of q near 0   =', 1/wgt0(1,1)
      write(6,"('  sum(wgt0) from Q0P=',d14.6)")sum(wgt0)
      if(bzcase()==2) then
        wgt0= wgt0*wgtq0p()/dble(nqbz)
        write(6,"('bzcase=2:  sum(wgt0_modified )=',d14.6)")sum(wgt0)
      endif
c
c     do i =1,nq0i
c     do ig=1,ngrp
c     write(66,'(2i3,f12.5)')i,ig,wgt0(i,ig)
c     enddo
c     enddo


c     pointer to optimal product basis
c     call defi    (iimdim,natom)
      allocate(imdim(natom))
      call indxmdm (nblocha,nclass,
     i     iclass,natom,
     o     imdim )
      if(niw/=0) then
c     generate gaussian frequencies x between (0,1) and w=(1-x)/x
c     call defdr   (ifreqx,niw)
c     call defdr   (ifreqw,niw)
c     call defdr   (iwx,niw)
c     call defdr   (iexpa,niw)
        allocate(freqx(niw),freqw(niw),wwx(niw)) !,expa(niw))
        call freq01x  (niw,    !ua,
     o        freqx,freqw,wwx)  !,expa)
      endif
c     ------ write energy mesh ----------
c     ifemesh = iopen('emesh.hsfp0'//xt(nz),1,-1,0)
c     deltax0 = 0d0
c     call writeemesh(ifemesh,freqw,niw,freq,nw,deltax0)


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

c---------------------------------------------
      if(tetra) then
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
c     nene = ntq*nq*nspin ! for energy points.
c     if(exchange) nene=0
c     allocate(wtet(nband,nspin,nqibz,0:3*nene),
c     &    eband(nband,nspin,nqibz), qz(3,nqibz),nstar(nqibz),
c     &    iene(3*ntq,nq,nspin), ene(0:3*nene) ) ! pointer for
        allocate(wtet(nband,nspin,nqibz,0:0),
     &        eband(nband,nspin,nqibz), qz(3,nqibz) ) ! pointer for
        call dcopy (3*nqibz,qibz,1,qz,1)
        do  is    = 1,nspin    !Readin eband
          do  iqi = 1,nqibz
c     iq       = idxk (qz(1:3,iqi),qbz,nqbz)
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
        volwgt = (3d0 - nspin) / ntetf ! ntetf was =6*n1*n2*n3
        call bzints2x(volwgt,eband,wtet(:,:,:,0),nqibz,nband,nband,
c     ,dum,nkp,ldim,nbmax,
     .        nspin,edummy,edummy,edummy,1,ef,2,nteti,idteti)

c
        ntot= sum(wtet)
        if(legas) then
          write(6,"(' tetra=T ef ntot nexact ratio=',15f12.6)") ef,ntot
     &           , ef**1.5d0/3d0/pi**2*voltot, ef**1.5d0 /3d0/pi**2*voltot/ntot
        else
          write(6,"(' tetra=T ef nvalence)=',15f12.6)") ef,ntot
        endif
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
        deallocate( eband, qz, ene ) ! pointer for
c     -- ibzx denote the index of k{FBZ for given k{1BZ.
        allocate(ibzx(nqbz))
c     call invkibzx(w(iirk),nqibz,ngrp,nqbz,
        call invkibzx(irk,nqibz,ngrp,nqbz,
     o        ibzx)
c     print *,' *** wtet sum 1=', sum( wtet)
c     stop ' -------------test end------------- '
c     call cputid(0)
      else
        allocate(wtet(1,1,1,1), iene(1,1,1)) !dummy
      endif
c     -end of tetra section --------------------------------------------




C...  Read pomatr
      if(smbasis()) then
        print *,' smooth mixed basis : augmented zmel'
        ifpomat = iopen('POmat',0,-1,0) !oct2005
        nkpo = nqibz+nq0i
        nnmx=0
        nomx=0
        do ikpo=1,nkpo
          read(ifpomat) q_r,nn_,no,iqx !readin reduction matrix pomat
          if(nn_>nnmx) nnmx=nn_
          if(no>nomx) nomx=no
          allocate( pomat(nn_,no) )
          read(ifpomat) pomat
          deallocate(pomat)
        enddo
        isx = iclose("POmat")
        ifpomat = iopen('POmat',0,-1,0) !oct2005
        allocate( pomatr(nnmx,nomx,nkpo),qrr(3,nkpo),nor(nkpo),nnr(nkpo) )
        do ikpo=1,nkpo
          read(ifpomat) qrr(:,ikpo),nn_,no,iqx !readin reduction matrix pomat
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


c-------------------------------------
      iii=ivsumxxx(irk,nqibz*ngrp)
      write(6,*) " sum of nonzero iirk=",iii, nqbz

c-----------------------------------------------------------
c     calculate the correlated part of the self-energy SEc(qt,w)
c-----------------------------------------------------------
c     arrays for sxcf.f
      nlnx4    = nlnx**4
      niwx     = max0 (nw+1,niw) !nw --->nw+1 feb2006
      allocate( ppb(nlnmx*nlnmx*mdimx*nclass),
     &     eq(nband), 
     &     kount(nqibz,nq),zsec(ntq,ntq,nq), 
     &     coh(ntq,nq) )


c     loop over spin ----------------------------------------------------
      do 2000 is = 1,nspinmx

        if(exchange) then
          write(ifsex2(is)) nspin, nq, ntq,nqbz,nqibz, n1,n2,n3
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
          write(ifsec2(is)) nspin, nq, ntq ,nqbz,nqibz  ,n1,n2,n3
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
     &           Zfactor(=1)'   !sf 21May02
        endif

c     initialise secq and kount
        zsec  = 0d0
        coh   = 0d0
        kount = 0
c
        if(ixc==3.and.nctot==0) goto 2001 !make dummy SEXcore

c     loop over rotations -------------------------------
        do 1000 irot = 1,ngrp
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
c     call rotra   (invr,w(igrp),w(ipos),plat,w(iiclass),
c     d              natom,ngrp,
c     o              w(itr),w(iiatomp) )


c     -- ppb= <Phi(SLn,r) Phi(SL'n',r) B(S,i,Rr)>
          call ppbafp_v2 (irot,ngrp,is,nspin,
     i           il,in,im,nlnm, !w(i_mnl),
     d           nl,nn,nclass,nlnmx,
     i           mdimx,lx,nx,nxx, !Bloch wave    
     i           cgr, nl-1,     !rotated CG
     i           ppbrd,         !radial integrals
     o           ppb)

c     c -- Rotated gvecc
c     c        call rotgvec(symgg(:,:,irot), nqibz,
c     c     i    ngcmx,ngcni,qbas,ngveccB,
c     c     o    ngveccBr)

c------------------------------------------------------
c     calculate the correlated part of the self-energy within GW
          ntqx = 0
          if(tetra.and.(.not.exchange)) then
            ntqx =3*ntq
          endif
c     loop over q
          do 1001 ip = 1,nq   !;write (*,*) ip,'  out of ',nq,'  k-points ' ! call cputid  (0)
ccc   if(ixc==2) then
c     iq  = idxk (q(1,ip),qbz,nqbz)
c     call rwdd1 (ifev(is),iq, nband,eq)
            call readeval(q(1,ip),is,eq)
            ! do iw = iwini,iwend     !sf 21May02
            do i  = 1,ntq
              omega(i) = eq(itq(i)) !sf 21May02     + 2d0*(dble(iw)-shtw)*deltaw
            enddo
            ! enddo
ccc   endif
c     if(sergeys) then   !sf 21May02
c     call sxcf_fal2_scy (kount, ip,q(1,ip),irot,itq,ntq,ef,esmr,
c     print *, 'nblochpmx=',nblochpmx

            if(nctot==0) then !this is to override -Mbounds check
              deallocate(ecore)
              allocate(ecore(1,2)) !this is dummry
            endif
cccccccccccccccccccccccccccccccc
        write(6,*) "goto sxcf_fal3_scz"
ccccccccccccccccccccccccccccccc
            call sxcf_fal3_scz (kount, ip,q(1,ip),irot,itq,ntq,ef,esmr,
     i              tiat(1:3,1:natom,invr),miat( 1:natom,invr),
c     v2
     i              nspin,is,   ! ifcphi,ifrb(is),ifcb(is),ifrhb(is),ifchb(is), 12jan2004
     i              ifrcw,ifrcwi,
     i              qbas,ginv,qibz,qbz ,wbz, nstbz, !iindxk,
     i              nstar,irk,  !kount,
c
c     i       iiclass,nblocha,i_mnlv,i_mnlc,iicore,incore,iimdim,
     i              iclass,nblocha,nlnmv,nlnmc, !w(i_mnlv),w(i_mnlc)
     i              icore,ncore,imdim,
     i              ppb,        ! pdb,dpb,ddb,
     i              freqx,wwx,  !expa,ua, !sf 21May02 (no w(ifreq))
     i              dw,         !deltaw,
     i              ecore(:,is),
c     
     d              nlmto,nqibz,nqbz,nctot,
c     i        index_qbz, n_index_qbz, !12jan2004
     d              nl,nnc,nclass,natom,
     d              nlnmx,mdimx,nbloch,ngrp,niw,nq, !niwx removed in sxcf_fal2_scz2 feb2006
c     i     nblochpmx, ngpn,ngcni,ngpmx,ngcmx, !12jan2004
     i              nblochpmx, ngpmx,ngcmx,
c     i     geigB(1,1,1,is), ngvecpB,ngveccBr, !12jan2004
c     i     ngveccBr,
     i              wgt0,nq0i,q0i, symgg(:,:,irot),alat,
     i              matmul(symgg(:,:,irot),shtvg(:,invr)),nband,
     i              ifvcfpout,
c     i     shtw, 
     i              exchange, screen, cohtest, ifexsp(is),
     i              omega,iSigma_en, ! iwini,iwend,            !sf 21May02
     i              nbmx,ebmx,  !takao 18June2003 
     i              pomatr, qrr,nnr,nor,nnmx,nomx,nkpo, !oct2005 for pomat 
     i    wklm,lxklm, !mar2012
     o              zsec,coh)   ! acuumulation variable
cccccccccccccccccccccccccccccccc
        write(6,*) "end of sxcf_fal3_scz"
ccccccccccccccccccccccccccccccc

ccccccccccccccccccccccfaleev 21May02 we use general mesh in sxcf_fal1, so
ccccccccccccccccccccccwe do not need sxce_v2
c     else
c     call sxcf_v2 (kount, ip,q(1,ip),irot,itq,ntq,ef,esmr,
c     i              tiat(1:3,1:natom,invr),miat( 1:natom,invr),
c     v2
c     i              ifev(is),ifcphi,nspin,is, !ifrb(is),ifcb(is),ifrhb(is),ifchb(is),
c     i              ifrcw,ifrcwi,
c     i              qbas,ginv,w(iqibz),w(iqbz),w(iwbz),iindxk,
c     i              w(instar),w(iirk),     !kount,
c     1
c     i        iiclass,nblocha,i_mnlv,i_mnlc,iicore,incore,iimdim,
c     i              ppb, ! pdb,dpb,ddb,
c     i              w(ifreq),w(ifreqx),w(iwx),w(iexpa),
c     i              ua,dw,  !deltaw,
c     i              ecore(:,is),
c     1
c     d              n1,n2,n3,nlmto,nqibz,nqbz,nctot,
c     d              nl,nnc,nclass,natom,
c     d              nlnmx,mdimx,nbloch,ngrp,nw,niw,niwx,nq,
c     1
c     i     nblochpmx, ngpn,ngcni,ngpmx,ngcmx,
c     i     geigB(1,1,1,is), ngvecpB,ngveccBr,
c     i     wgt0,nq0i,q0i, symope(:,:,irot),alat,
c     i     matmul(symope(:,:,irot),shtvg(:,invr)),nband,
c     i     ifvcfpout,
c     1     i     shtw,
c     i     exchange, screen, cohtest, ifexsp(is),
c     i        omega, iwini,iwend,
c     o        zsec,coh) ! acuumulation variable
c     endif
cccccccccccccccccccccccccccccccccccccccccc
c     iii = ivsum(kount,nqibz*nq)
c     print *," sumkount 2=",nqibz,nq,iii
c     stop "--- kcount test end --- "
cccccccccccccccccccccccccccccccccccccccccc
 1001     continue            !end of q-loop
 1000   continue               !end of rotation-loop


c     check that all k { FBZ have been included
c     if (ivsum(kount,nqibz*nq) /= nqbz*nq) then
c     iii = ivsum(kount,nqibz*nq)
c     print *," ivsum=",iii, nqbz*nq
c     stop 'hsfp0: missing k-pts'
c     endif

c---------------------------------
c     electron gas bare exchange (exact)
        if (legas.and.exchange) then
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
     &              eqx(1,ip,is),        hartree*dreal(zsec(1,1,ip)) !sf 21May02
            write (6,"(' === diff     =',2f14.6)") 
     &              rydberg()*(qm**2-efz)-eqx(1,ip,is)
     &              , xsex - hartree*dreal(zsec(1,1,ip)) !sf 21May02

            write (661,"(' qm True qm-ef Sx=',3f14.6)") 
     &              qm,rydberg()*(qm**2-efz), xsex
            write (662,"(' qm Num  qm-ef Sx=',3f14.6)") 
     &              qm,eqx(1,ip,is),     hartree*dreal(zsec(1,1,ip)) !sf 21May02
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
        if(exchange) then
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
     &                 hartree*dreal(zsec(i,i,ip)) !sf 21May02
              if( eqx(i,ip,is)<1d20.and.abs(zsec(i,i,ip))/=0d0 ) then  !takao june2009
                write(6,"(' j iq isp=' i3,i4,i2,'  q=',3f8.4,
     &                 '  eig=',f10.4,'  Sx=',f10.4)")
     &                 itq(i),ip,is, q(1:3,ip), eqx(i,ip,is),
     &                 hartree*dreal(zsec(i,i,ip)) !sf 21May02
              endif
            end do
            write(ifsex2(is)) is, q(1:3,ip), zsec(1:ntq,1:ntq,ip) !SEC_nn' out
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
ccccccccccccccccccccccccccccccfaleev 21May02
c     dscdw1  = (dreal(zsec(1,i,ip)) - dreal(zsec(0,i,ip)))
c     &             /deltaw
c     zfac1   = 1d0/(1d0-dscdw1)
c     dscdw2  = (dreal(zsec(0,i,ip)) - dreal(zsec(-1,i,ip)))
c     &             /deltaw
c     zfac2   = 1d0/(1d0-dscdw2)
c     dscdw  = (dreal(zsec(1,i,ip)) - dreal(zsec(-1,i,ip)))
c     &             /(2d0*deltaw)
c     zfac   = 1d0/(1d0-dscdw)
              if( eqx(i,ip,is)<1d20.and.abs(zsec(i,i,ip))/=0d0 ) then  !takao june2009
                write(6,"(' j iq isp=' i3,i4,i2,'  q=',3f8.4,
     &          '  eig=',f8.4,'  Re(Sc) =',f8.4,'  Img(Sc) =',f8.4 )") !sf 21May02
     &          itq(i),ip,is, q(1:3,ip), eqx(i,ip,is), 
     &          hartree*dreal(zsec(i,i,ip)), !sf 21May02
     &          hartree*dimag(zsec(i,i,ip)) !,zfac,zfac1,zfac2  !sf 21May02
              endif
              write(ifoutsec,"(3i5,3d24.16,3x,d24.16,3x,d24.16, 
     &                 3x,d24.16)")
     &                 itq(i),ip,is, q(1:3,ip), eqx(i,ip,is), 
     &                 hartree*dreal(zsec(i,i,ip)), !sf 21May02
     &                 hartree*dimag(zsec(i,i,ip)) !,zfac,zfac1,zfac2  !sf 21May02
cccccccccccccccccccccfaleev 21May02
c     if(cohtest) then
c     write(671,"( i3,i4,i2, 3f8.4,
c     &     '  eig=',f8.4,'  coh =',f10.4)")
c     &     itq(i),ip,is, q(1:3,ip), eqx(i,ip,is),
c     &     hartree*coh(i,ip)
c     endif
ccccccccccccccccccccc
            end do
            write(ifsec2(is)) is, q(1:3,ip), zsec(1:ntq,1:ntq,ip) !SEC_nn' out
          end do
c-------------------------
cccccccccccccccccccccccccfaleev 21May02 do not use ixc==4
c     elseif(ixc==4) then !spectrum mode
c     ifoutsec=ifsecomg(is)
c     do ip = 1,nq
c     do i  = 1,ntq
c     write(6,"(' --- j iq isp=' i3,i4,i2,' q=',3f8.4,
c     &     ' eig=',f8.4)")
c     &      itq(i),ip,is, q(1:3,ip), eqx(i,ip,is)
c     do iw = iwini,iwend
c     write(6,"(' omega-ef=',f8.4,'  Sc=',2f8.4)")
c     &     (omega(i,iw)-ef)*rydberg(), hartree*zsec(iw,i,ip)
c
c     write(ifoutsec,"(4i5,3d24.16,3x,d24.16,
c     &     2x,d24.16,x,2d24.16)")
c     &     iw,itq(i),ip,is, q(1:3,ip), eqx(i,ip,is)
c     &     ,(omega(i,iw)-ef)*rydberg(),  hartree*zsec(iw,i,ip)
c     end do
c     write(ifoutsec,*)
c     write(6,*)
c     end do
c     end do
        endif                  !ixc
 2000 continue                  !end of spin-loop

c------------
c     close files
c------------
      isx = iclose ('wc.d')
      isx = iclose ('wci.d')
      isx = iclose ('hbe.d')
c     isx = iclose ('RBU')
c     isx = iclose ('CBU')
c     isx = iclose ('RHBU')
c     isx = iclose ('CHBU')
c     isx = iclose ('EVU')
c     isx = iclose ('RBD')
c     isx = iclose ('CBD')
c     isx = iclose ('RHBD')
c     isx = iclose ('CHBD')
c     isx = iclose ('EVD')
c     isx = iclose ('CPHI')

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


c     formats
c     6500 format (' spin =',i2,'   nbloch =',i4,'  deltaw =',f7.4,' a.u.',
c     .'   ef=',f10.5,' Rydberg' )
c
      call cputid(0)
c     if(ixc==1.and.sergeys)stop " OK! hsfp0: Sergey's Exchange mode"
c     if(ixc==2.and.sergeys)stop " OK! hsfp0: Sergey's Correlation mode"
      if(ixc==1 ) stop ' OK! hsfp0_sc: Exchange mode'
      if(ixc==2 ) stop ' OK! hsfp0_sc: Correlation mode'
      if(ixc==3 ) stop ' OK! hsfp0_sc: Core-exchange mode'
c     if(ixc==4 ) stop ' OK! hsfp0: spectrum mode'
      end

