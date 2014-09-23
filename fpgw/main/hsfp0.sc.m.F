      program hsfp0_sc
      use m_readqg,only: readqg,readngmx
      use m_readeigen,only: init_readeigen,init_readeigen2,readeval,lowesteval
      use m_read_bzdata,only: nqbz,nqibz,nqbzw,nteti,ntetf
     &     ,n1,n2,n3,qbas,ginv,qbasmc,qbz,wbz,qibz,wibz,qbzw,idtetf,ib1bz,idteti
     &     ,nstar,irk,nstbz,ngrp2=>ngrp,qibz_r,nqibz_r,  read_bzdata
      use m_genallcf_v3,only: genallcf_v3,
     &     nclass,natom,nspin,nl,nn,nnv,nnc, ngrp,
     &     nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, nctot,niw,nw,
     &     alat,ef, diw,dw,delta,deltaw,esmr,symgrp,clabl,iclass,nindxv,nindxc,ncwf,
     &     invg, il, in, im, ilnm, nlnm, ilv,inv,imv,  ilnmv, nlnmv,
     &     ilc,inc,imc,  ilnmc, nlnmc,nindx,konf,icore,ncore,
     &     occv,unoccv ,occc,unoccc, nocc, nunocc, plat, pos,z,ecore, freq, symgg ! symgg=w(igrp)
      use keyvalue,only: getkeyvalue
      use m_mpi  ! MIZUHO-IR
      implicit none
!     ! =  Calculates the  self-energy \Sigma in GW approximation, including Off-diagonal components. =
!     !
c     SEx(q,itp,itpp) = <psi(q,itp) |SEx| psi(q,itpp)>
c     SEc(q,itp,itpp) = <psi(q,itp) |SEc| psi(q,itpp)>
c     
c     where, SEc(r,r';w) = (i/2pi) < [w'=-inf,inf] G(r,r';w+w') Wc(r,r';w') >
c     
c     This routine requirs a number from standard Input.
c     Look into script calling this routine ecalj/fpgw/exec/gwsc. E.g, this gwsc calls echo 2|../exec/hsfp0_sc >lsc when mode=2.
c     
c     mode= 1: exchange    mode SEx, the exchange part of the self-energy
c     mode= 2: correlation mode SEc, the correlated part of the self-energy
c     mode= 3: core exchange mode SEXcore
c     !mode= 4: plot spectrum function ---See manual ---> this is in hsfp0.
c     
c     iSigMode parameter which determines approximation for self-energy is given by GWinput file as iSigMode.
c     iSigMode=iSigMode
c     iSigMode==0 SE_nn'(ef)+image integr:delta_nn'(SE_nn(e_n)-SE_nn(ef))
c     iSigMode==1 SE_nn'(ef)+delta_nn'(SE_nn(e_n)-SE_nn(ef))
c     iSigMode==2 SE_nn'((e_n+e_n')/2)
c     iSigMode==3 (SE_nn'(e_n)+SE_nn'(e_n'))/2  <--- this is mainly used
c     iSigMode==5 delta_nn' SE_nn(e_n)
c     output file should contain hermitean part of SE for energies to be real
c     (for example, hermitean conjunction of SE_nn'(e_n) means SE_n'n(e_n')^* )
c     
c     Apr 2002 takao kotani. multiple argumentation wave per l.
c     This hsfp0 is build from hsec10.f by F.Aryasetiawan.
!     ! -----------------------------------------------------------------------------------

c------------------------------------------------------------
c     real(8),parameter :: ua  = 1d0    ! constant in w(0)exp(-ua^2*w'^2) to take care of peak around w'=0
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
                                ! difficult and I gave up for a while.
                                ! If you want to calculate with tetra=T for exchange, you
                                ! have to uncomment tetra related part in
                                ! sxcf.f, and a part calling sxcf in this routine. Note wtet wtetef!
                                ! They sometimes cause array destruction if you run tetra=T without comment them.
c------------------------------------
      real(8) :: shtw
c     integer :: mxclass,ngnmax,mbytes,mwords,iwksize,
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
      integer::
     &     ixc,iopen,ifhbed, nprecb,mrecb,mrece,nlmtot,nqbzt, nband,
     &     ibas,ibasx,ngpmx,nxx,ngcmx,nbloch,ifqpnt,ifwd,
     &     nprecx,mrecl,nblochpmx2,nwt,niwt, nqnum,mdimx,nblochpmx,
     &     noccxv,maxocc,noccx,ifvcfpout,iqall,iaf,ntq, !ifrcw,ifrcwi,  
     &     i,k,nspinmx, nq,is,ip,iq,idxk,ifoutsex,iclose,nq0i,ig,
     &     mxkp,nqibzxx,ntet,nene,iqi, ix,iw,
     &     nlnx4,irot,invr,invrot,ivsum, ifoutsec, !niwx,
c     
c     &   ifrb(2),ifcb(2),ifrhb(2),ifchb(2), ifev(2),
     &     ifsec(2)             
     &     ,ifxc(2),ifsex(2), ifphiv(2),ifphic(2),ifec,ifexsp(2),
     &     ifsex2(2),ifsec2(2), !out S_nn'
     &     ifsecomg(2),ndble=8
c     
c     real(8) :: alat,ef,diw,dw,delta,pi,tpia,vol,voltot,rs,alpha,
      real(8) :: pi,tpia,vol,voltot,rs,alpha,
     &     qfermi,efx,valn,efnew,edummy,efz,qm,xsex,egex,edummyd(1),
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
      integer,allocatable :: itq(:)
      real(8),allocatable    :: q(:,:)

c     takao
      integer,allocatable :: ngvecpB(:,:,:), !ngveccB(:,:,:),
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
      integer,allocatable :: iclasst(:), invgx(:), miat(:,:)
      real(8),allocatable    :: tiat(:,:,:),shtvg(:,:)

c     tetra
      real(8),allocatable :: qz(:,:),qbzxx(:),wbzxx(:),wtet(:,:,:,:),
     &     eband(:,:,:), ene(:) !,ecore(:,:)
      integer,allocatable ::idtetx(:,:),idtet(:,:),ipq(:)
     &     ,iene(:,:,:),ibzx(:) !,nstar(:)
c     real(8) :: qbasmc(3,3)
      integer ::ib,iqx,igp,iii,ivsumxxx,isx,iflegas, iqpntnum
c     
      real(8),allocatable   :: eex1(:,:,:),exsp1(:,:,:),qqex1(:,:,:,:)
      integer,allocatable:: nspex(:,:),ieord(:),itex1(:,:,:)
      real(8)    :: qqex(1:3), eex,exsp,eee, exwgt,deltax0
      integer :: itmx,ipex,itpex,itex,nspexmx,nnex,isig,iex,ifexspx
     &     ,ifexspxx ,ifefsm, nq0ix,ifemesh,nz
      character(3)  :: charnum3
      character(12) :: filenameex
      logical :: exspwrite=.false.
      character*8 xt

      integer :: iSigMode,ifinin ,idummy !iwini,iwend       !sf 21May02

      real(8),allocatable:: omega(:) !sf 21May02
c     real(8) ::  omegamax,dwplot,omegamaxin
c     logical :: sergeys     !sf 21may02

c     real(8)   :: ebmx
c     integer:: nbmx
      real(8)   :: ebmx(2)
      integer:: nbmx(2)

      real(8):: volwgt

      integer::nwin, incwfin
      real(8)::efin
      integer,allocatable::imdim(:)
      real(8),allocatable::freqx(:),freqw(:),wwx(:) !,expa(:)

      integer:: bzcase,  ngpn1,mrecg,ngcn1
      real(8)   :: wgtq0p,quu(3)

      character(2):: soflag
      integer:: ifianf

      logical ::smbasis
      integer:: ifpomat,nkpo,nnmx,nomx,ikpo,nn_,no
      real(8):: q_r(3)
      real(8),allocatable:: qrr(:,:)
      integer,allocatable:: nnr(:),nor(:)
      complex(8),allocatable:: pomatr(:,:,:),pomat(:,:)
      logical allq0i       
      integer:: nw_i
      logical:: exonly
      real(8):: wex
!     ! newaniso mode
      real(8),allocatable:: vcousq(:),dmlx(:,:),epinvq0i(:,:),wklm(:),vcoud(:)
      complex(8),allocatable:: zcousq(:,:)
      logical:: newaniso
      integer:: ifvcoud,lxklm,ifidmlx,iqq

      integer,allocatable:: irkip_all(:,:,:,:),irkip(:,:,:,:)

!     !
      integer,allocatable:: nrkip_all(:,:,:,:),nrkip(:,:,:,:)
      integer,allocatable:: neibz(:),nwgt(:,:),ngrpt(:),igx(:,:,:),igxt(:,:,:),eibzsym(:,:,:)
      integer:: iqxend,iqxini
      integer:: l2nl,igrp,kx,kr
      logical :: iprintx,tiii,timereversal, eibz4sig,tiiiout

      integer:: ifevec,ifvxc,nsym,nhdim,it,nblk,iband,napw,ldim,ierr,ispx,nbsize,nbsizemx
     &     ,iblk1,iblk2,ii1,ii2,ie1,ie2,ne1,ne2,iqxx, ndimhx, nspx,nnnx
      integer,allocatable::iblki(:),iblke(:),nbandmx(:,:)
      complex(8),allocatable:: evec(:,:),evec_inv(:,:),evecrot(:,:),rmatjj(:,:,:)
      real(8),allocatable::evaliq(:)
      real(8)::tolry=1d-4,qqqx(3),qtarget(3)
      complex(8),allocatable::zsect(:,:)

      logical :: selectqp=.false.,diagonly=.false.
      integer:: ret,dest,nnn !,ifiwqfac
c      real(8),allocatable :: qbzwww(:,:),wqfac(:),wbzz(:)
      character(128) :: ixcc
      real(8):: eftrue,esmref !jan2013
      real(4):: time_red1,time_red2
      integer:: timevalues(8) ,ibz

      integer:: procid,nrankv,ifigwb_,ifigwx1_,ifigwx2_,ifvxc_,ifevec_,ifiproc,iqqxx,
     &  isp,ixx,ixxx,nqixx,nspxx,ispxx,iqbz
      character*256:: extn,ext
      character*256,allocatable:: extp(:)
      integer,allocatable:: ifevec__(:),ifvxc__(:),iprocq(:,:)
      integer:: verbose
c-----------------------------------------------------------------------

c---------------------------------------
      call MPI__Initialize()      ! MIZUHO-IR
      call date_and_time(values=timevalues)
      write(*,'(a,9i5)')'dateandtime1=',mpi__rank,timevalues(1:8)

!$      call omp_set_nested(1)
!TIME0
!TIME0 
      hartree=2d0*rydberg()
!     !     mode switch. --------------
      if(cohtest) then
         screen = .true.
         ixc = 2; nz=0
         open(671,file='COH')
      elseif(mpi__root) then
         write(6,*) ' --- Choose modes below ------------'
         write(6,*) '  Sx(1) Sc(2) ScoreX(3) '
         write(6,*) '  [option --- (+ QPNT.{number} ?)] '
         write(6,*) ' Add 1000, eg, 1001 is diagonal only mode for one-shot Z=1'
         write(6,*) ' --- Put number above ! ------------'
         call readin5(ixc,nz,idummy)
         if(ixc>1000) then      !selected QP
            ixc=mod(ixc,1000)
            selectqp=.true.
            diagonly=.true.
         endif
      endif
      call MPI__Broadcast(ixc)
      call MPI__Broadcast(nz)
c     call MPI__Broadcast(idummy)

! END   MIZUHO-IR
      if(mpi__root) call headver('hsfp0_sc',ixc)
      write(ixcc,"('.mode=',i4.4)")ixc
      call MPI__consoleout('hsfp0_sc'//trim(ixcc))
      write(6,*) ' ixc nz=',ixc, nz
Cstop2rx 2013.08.09 kino      if(ixc==0) stop ' --- ixc=0 --- Choose computational mode!'
      if(ixc==0) call rx( ' --- ixc=0 --- Choose computational mode!')

!     ! ===  readin BZDATA. See gwsrc/rwbzdata.f ===
!     ! See use m_read_bzdata,only: at the top of this routine
c--------readin data set when you call read_BZDATA ---------------
c     integer::ngrp,nqbz,nqibz,nqbzw,nteti,ntetf
c     integer:: n1,n2,n3
c     real(8):: qbas(3,3),ginv(3,3),qbasmc(3,3)
c     real(8),allocatable:: qbz(:,:),wbz(:),qibz(:,:)
c     &    ,wibz(:),qbzw(:,:)
c     integer,allocatable:: idtetf(:,:),ib1bz(:),idteti(:,:)
c     &    ,nstar(:),irk(:,:),nstbz(:)    
c-----------------------------------------------------------------
      call read_BZDATA()
      print *,' nqbz =',nqbz
      print *,' nqibz ngrp=',nqibz,ngrp2
      call pshprt(60)
!     ! -----------------------------------------------------------

!     ! === readin GWIN and LMTO, then allocate and set datas. ===
!     ! See use m_genallcf_v3,only: at the top of this routine
      nwin   = 0                !Readin nw from NW file
      efin=-999d0               !not readin EFERMI
      if(ixc==3) then; incwfin= -2 !core exchange mode
      else           ; incwfin= -1 !use 7th colmn for core at the end section of GWIN
      endif
      call genallcf_v3(nwin,efin,incwfin) !in module m_genallcf_v3
Cstop2rx 2013.08.09 kino      if(ngrp/= ngrp2) stop 'ngrp inconsistent: BZDATA and LMTO GWIN_V2'
      if(ngrp/= ngrp2) call rx( 'ngrp inconsistent: BZDATA and LMTO GWIN_V2')
c---  These are allocated and setted by genallcf_v3
c     integer::  nclass,natom,nspin,nl,nn,nnv,nnc, ngrp,
c     o  nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, nctot,niw,nw
c     real(8) :: alat,ef, diw,dw,delta,deltaw,esmr
c     character(120):: symgrp
c     character(6),allocatable :: clabl(:)
c     integer,allocatable:: iclass(:)
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
!     ! -----------------------------------------------------------
      esmref=esmr

!     ! reading self-energy mode parameter from file 'GWinput'
      call readd_iSigma_en(ifinin,iSigMode) 
c     allocate(ecore(nctot,nspin)) !core energies
c     do is = 1,nspin
c     if (nctot > 0) call catch1 (w(iecore),is,nctot,2,ecore(:,is)) !core energies
c     enddo

c---  Get maximums takao 18June03
      call getnemx8(nbmx,ebmx)

ccccccccccccccccccccccccccccc
!! forget ebmx2 nbmx2 
      nbmx(2)=9999999
      ebmx(2)=1d10
ccccccccccccccccccccccccccccc
      write(6,"('  nbmx ebmx from GWinput=',2i8,2d13.5)") nbmx,ebmx
!     !     See sxcf_fal2z--------
!     !     nbmx1 ebmx1: to set how many bands of <i|sigma|j>  do you calculate.
!     !     nbmx2 ebmx2: to restrict num of bands of G to calculate G \times W

c-------------------------------------------------------------------
c     if (nclass > mxclass) stop ' hsfp0: increase mxclass'
!!!!  WE ASSUME iclass(iatom)= iatom !!!!!!!!!!!!!!!!!!!!!!!!!
Cstop2rx 2013.08.09 kino      if (nclass /= natom ) stop ' hsfp0: nclass /= natom ' ! We assume nclass = natom.
      if (nclass /= natom ) call rx( ' hsfp0: nclass /= natom ')
      print *,' hsfp0: end of genallcf_v3'

      call pshprt(30)
      pi   = 4d0*datan(1d0)
      tpia = 2d0*pi/alat
      call dinv33(plat,1,xxx,vol)
      voltot = dabs(vol)*(alat**3)
      shtw = 0d0
c     if(esmr<1d-5) shtw=0.01d0 ! Ferdi's shift to avoid resonance effect(maybe) !I had this until sep2012
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
Cstop2rx 2013.08.09 kino         if(tetra) stop 'legas You have to give ef of  tetrahedron'
         if(tetra) call rx( 'legas You have to give ef of  tetrahedron')
      endif

c      ifexsp=0
      if(ixc==1) then
         exchange=.true.
         write(6,*) ' --- Exchange mode --- '
         if(mpi__root) then
            ifxc(1)  = iopen('XCU'//xt(nz),1,-1,0)
            ifsex(1) = iopen('SEXU'//xt(nz),1,-1,0)
            ifsex2(1)= iopen('SEX2U',0,-1,0) !out SEX_nn'
            if (nspin == 2) then
               ifxc(2)  = iopen('XCD'//xt(nz),1,-1,0)
               ifsex(2) = iopen('SEXD'//xt(nz),1,-1,0)
               ifsex2(2)= iopen('SEX2D',0,-1,0) !out SEX_nn'
            endif
         endif
c         INQUIRE (FILE = 'EXspTEST', EXIST = exspwrite)
c         if(exspwrite) then
c            print *,'--- Find EXspTEST ExspectrumWrite=',exspwrite
c            print *,'--- esmr is chosen to be 2d0 Ry'
c            esmr= 2d0
c            do is=1,nspin
c               ifexsp(is)  = iopen('EXSP.'//char(48+is),1,-1,0)
c            enddo
c         endif
      elseif(ixc==2) then
         exchange=.false.
         write(6,*) ' --- Correlation mode --- '
         if(cohtest) write(6,*) ' COH calculation mode. Results in COH'
         if(mpi__root) then
            ifsec(1)   = iopen('SECU'//xt(nz),1,-1,0) ! output files
            ifsec2(1)= iopen('SEC2U',0,-1,0) !out SEC_nn'
            if (nspin == 2) 
     .           ifsec(2)   = iopen('SECD'//xt(nz),1,-1,0)
            ifsec2(2)= iopen('SEC2D',0,-1,0) !out SEC_nn'
         endif
      elseif(ixc==3) then
         exchange=.true.
         esmr=0d0
         write(6,*) ' --- CORE Exchange mode --- '
         if(mpi__root) then
            ifsex(1)   = iopen('SEXcoreU'//xt(nz),1,-1,0)
            ifsex2(1)= iopen('SEXcore2U',0,-1,0) !out SEXcore_nn'
            if (nspin == 2) then
               ifsex(2)   = iopen('SEXcoreD'//xt(nz),1,-1,0)
               ifsex2(2)= iopen('SEXcore2D',0,-1,0) !out SEXcore_nn'
            endif
         endif   
ccccccccccccccccccccccccccccccccccccccccccfaleev 21May02 do not use ixc==4
c     elseif(ixc==4) then
c     write(6,*) ' --- Spectrum function Sigma(\omega) mode --- '
c     exchange=.false.
c     ifsecomg(1) = iopen('SEComgU'//xt(nz),1,-1,0) ! output files
c     if (nspin == 2)
c     .  ifsecomg(2) = iopen('SEComgD'//xt(nz),1,-1,0)
      else
Cstop2rx 2013.08.09 kino         stop ' hsfp0: Need input (std input) 1(Sx) 2(Sc) or 3(ScoreX)!'
         call rx( ' hsfp0: Need input (std input) 1(Sx) 2(Sc) or 3(ScoreX)!')
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
      write(6,'("    niw nw dw   =",2i6,f13.6)') niw,nw,dw

c>>   read dimensions of wc,b,hb
      ifhbed     = iopen('hbe.d',1,0,0)
      read (ifhbed,*) nprecb,mrecb,mrece,nlmtot,nqbzt, nband,mrecg
Cstop2rx 2013.08.09 kino      if (nprecb == 4) stop 'hsfp0: b,hb in single precision'
      if (nprecb == 4) call rx( 'hsfp0: b,hb in single precision')

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
      if(verbose()>=40) write (*,*)' hsfp0.sc.m: end of mptauof'

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

!     ! WVR, WVI
      if(.not.exchange.or.(exchange.and.screen)) then !screen means screened exchange case
         ifwd = iopen('WV.d',1,-1,0)
         read (ifwd,*) nprecx,mrecl,nblochpmx,nwt,niwt, nqnum,nw_i
         ifwd = iclose('WV.d')
         write(6,"(' Readin WV.d =', 10i5)") nprecx,mrecl,nblochpmx,nwt,niwt, nqnum, nw_i
         ifwd =iclose('WV.d')
         call checkeq(nprecx,ndble)
         nw = nwt
Cstop2rx 2013.08.09 kino         if(niwt /= niw) stop 'hsfp0: wrong niw'
         if(niwt /= niw) call rx( 'hsfp0: wrong niw')
c!!   direct access files WVR and WVI which include W-V.
c     if(.not.exchange) then
c     ifrcw  = iopen('WVR',0,-1,mrecl)
c     ifrcwi = iopen('WVI',0,-1,mrecl)
c     endif 
      endif

      if(tetra) goto 201

!! == Determine Fermi energy ef for given valn (legas case) or corresponding charge given by z and konf.==
!!    When esmr is negative, esmr is geven automatically by efsimplef.
      call efsimplef2a(nspin,wibz,qibz,ginv, 
     i     nband,nqibz
     i     ,konf,z,nl,natom,iclass,nclass
     i     ,valn, legas, esmref,  !!! valn is input for legas=T, output otherwise.
     i     qbz,nqbz             ! index_qbz, n_index_qbz,
     o     ,efnew)
      if(ixc/=3) ef = efnew
      eftrue = efnew
!! ==== check total ele number =====
      ntot  = nocctotg2(nspin, ef,esmr, qbz,wbz, nband,nqbz)
      print *,' ef    =',ef
      print *,' esmr  =',esmr
      print *,' valn  =',valn
      print *,' ntot  =',ntot
!! == Core-exchange case. ef means just below the valence eigenvalue (to take only core in sxcf).==
      if(ixc==3) then
         ef = lowesteval() -1d-3 !lowesteigen(nspin,nband,qbz,nqbz) - 1d-3 !lowesteb was
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
Cstop2rx 2013.08.09 kino            stop 'hsfp0 ixc=3: ecore>evalence. '
            call rx( 'hsfp0 ixc=3: ecore>evalence. ')
         endif
      endif

c-------------------------
 201  continue
c-------------------------
      call init_readeigen2(mrecb,nlmto,mrecg) !initialize m_readeigen
c     the coulomb matrix for exchange.
      if((.not.newaniso()).and.exchange) ifvcfpout = iopen('VCCFP',0,-1,0)

!!    read q-points and states
      if(selectqp .and. mpi__root) then
         call getkeyvalue("GWinput","<QPNT>",unit=ifqpnt,status=ret)
         lqall      = .false.
         laf        = .false.
         call readx   (ifqpnt,10)
         read (ifqpnt,*) iqall,iaf
         if (iqall == 1) lqall = .true.
         if (iaf   == 1)   laf = .true.
         call readx   (ifqpnt,100)
c         read (ifqpnt,*) ntq 
c         allocate( itq(ntq) )
c         read (ifqpnt,*) (itq(i),i=1,ntq)
         if (lqall) then        !all q-points case
            nq         = nqibz
            allocate(q(3,nq))
            call dcopy   (3*nqibz,qibz,1,q,1)
         else
            call readx   (ifqpnt,100)
            read (ifqpnt,*) nq
            allocate(q(3,nq))
            do       k = 1,nq
               read (ifqpnt,*) i,q(1,k),q(2,k),q(3,k)
            enddo
         endif
         nspinmx = nspin
         if (laf) nspinmx =1
         close(ifqpnt)
      else
!     q-points
         if(bzcase()==1)then
            nq = nqibz
            allocate(q(3,nq))
            call dcopy   (3*nqibz,qibz,1,q,1)
         else
            nq = nqibz_r
            allocate(q(3,nq))
            do iq= 1,nqibz_r
               q(:,iq)=qibz_r(:,iq)
            enddo
         endif
      endif
!
      call MPI__Broadcast(nq)
      if(mpi__root) then
        do dest=1,mpi__size-1
          call MPI__REAL8send(q,3*nq,dest)
        enddo
      else  
        call MPI__REAL8recv(q,3*nq,0)
      endif  

!! == Determine ntq.  See also ntqxx in sxcf_fal.sc.F ntq should be common for all ixc modes.==
      if(mpi__root) then
         allocate(eqt(nband))
         ntq=0
         do is = 1,nspin
            do ip = 1,nq
               call readeval(qibz(1,ip),is, eqt)
               do iband=1,nband
                  ntq = max(iband,ntq)
                  if(eqt(iband)-eftrue>ebmx(1)) exit
               enddo  
            enddo
         enddo
         ntq = min(ntq, nbmx(1))
         deallocate(eqt)
         write(*,"(' nband ntq ef eftrue=', 2i5,2f9.4)")nband,ntq,ef,eftrue
      endif
      call MPI__Broadcast(ntq)

      allocate (itq(ntq))
      do i = 1, ntq
         itq(i) = i !itq is used also in hsfp0.m.F
      enddo
      do iq=1,nq
         write(6,'(" Target iq q=",i5,3f9.4)')iq,q(:,iq)
      enddo

!! ANF mode ---> not so much used now... need to revise it in future.
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

!!
      allocate(omega(ntq))
!     ! read LDA eigenvalues
      allocate(eqx(ntq,nq,nspin),eqx0(ntq,nq,nspin),eqt(nband))
      do is = 1,nspin
         do ip = 1,nq
            call readeval(q(1,ip),is,eqt)
            eqx0(1:ntq,ip,is) = eqt(itq(1:ntq))
            eqx (1:ntq,ip,is) = rydberg()*(eqt(itq(1:ntq))- eftrue)
         enddo
      enddo
      deallocate(eqt)
      call winfo(6,nspin,nq,ntq,is,nbloch,ngpn1,ngcn1,nqbz,nqibz,ef,deltaw,alat,esmr)

c-------------------------
c     LDA exchange-correlation
c-------------------------
c     - takao Vxc is readin from VXCFP
      if(ixc==1) then
         allocate(  vxcfp(ntq,nq,nspin) )
         call rsexx(nspin,itq,q,ntq,nq, ginv, vxcfp) !add ginv july2011
         if(mpi__root) then
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
                  if(eqx(i,ip,is) <1d20.and.vxcfp(i,ip,is)/=0d0) then !takao june2009. See lmf2gw (evl_d=1d20; in Ry.. but eqx is in eV. no problem for inequality).
                     write(6,"(' j iq isp=' i3,i4,i2,'  q=',3f8.4,
     &                    '  eig=',f10.4,'  Sxc(LDA)=',f10.4)") 
     &                    itq(i),ip,is, q(1:3,ip), eqx(i,ip,is),
     &                    vxcfp(i,ip,is)
                  endif
               end do
            end do
            if(is==1) isx = iclose('XCU'//xt(nz))
            if(is==2) isx = iclose('XCD'//xt(nz))
         enddo                  !     end of spin-loop
         endif !mpi__root
         deallocate(vxcfp)
      endif

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
         allocate(freqx(niw),freqw(niw),wwx(niw)) !,expa(niw))
         call freq01x  (niw,    !ua,
     o        freqx,freqw,wwx)  !,expa)
      endif
c     ------ write energy mesh ----------
c     ifemesh = iopen('emesh.hsfp0'//xt(nz),1,-1,0)
c     deltax0 = 0d0
c     call writeemesh(ifemesh,freqw,niw,freq,nw,deltax0)


!     ! === readin Vcoud and EPSwklm for newaniso()=T ===
      if(newaniso()) then
         ifidmlx = iopen('EPSwklm',0,0,0)
         read(ifidmlx) nq0ix,lxklm
         if(nq0i/=nq0ix) then
            print *,'nq0i from EPSwklm /= nq0i',nq0i,nq0ix
Cstop2rx 2013.08.09 kino            stop 'nq0i from EPSwklm /= nq0i'
            call rx( 'nq0i from EPSwklm /= nq0i')
         endif
         allocate( dmlx(nq0i,9))
         allocate( epinvq0i(nq0i,nq0i) )
         allocate( wklm((lxklm+1)**2))
         read(ifidmlx) dmlx, epinvq0i
         read(ifidmlx) wklm
         ifidmlx = iclose('EPSwklm')
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
     .        nspin,edummy,edummy,edummyd,1,ef,2,nteti,idteti)

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
c     calculate the the self-energy SEx(ip) or SEc(ip)
c-----------------------------------------------------------
!     ! EIBZ symmetrization
      if(eibz4sig()) then
         allocate(nwgt(nqbz,1:nq),igx(ngrp*2,nqbz,nq))
         allocate(igxt(ngrp*2,nqbz,nq), eibzsym(ngrp,-1:1,nq))
         iqxini=1
         iqxend=nq
c         write(6,"('TimeRevesal switch = ',l1)") timereversal()
c         call eibzgen(nq,symgg,ngrp,q(:,iqxini:iqxend),
c     &        iqxini,iqxend,qbz,nqbz,timereversal(),ginv,iprintx,
c     o        nwgt,igx,igxt,eibzsym,tiii)
!     ! Check timereversal is required for symmetrization operation or not. If tiii=timereversal=F is enforced,
!     ! the symmetrization procedure in x0kf_v4h becomes a little time-consuming.
cccccccccccccccccccc
         tiii=.false.    !Enforce no time reversal. time reversal not yet...
         print *,'NOTE:TimeReversal not yet implemented in hsfp0.sc.m.F'
         write(6,"('=== goto eibzgen === used timereversal=',l1)")tiii
         iprintx=.false.
         if(mpi__root) iprintx=.true.
         call eibzgen(nq,symgg,ngrp,q(:,iqxini:iqxend),
     &        iqxini,iqxend,qbz,nqbz,tiii,ginv,iprintx,
     o        nwgt,igx,igxt,eibzsym,tiiiout)
c     call PBindex(natom,lx,l2nl,nx) !all input. this returns requied index stored in arrays in m_pbindex.
                                ! PBindex: index for product basis.  We will unify this system; still similar is used in ppbafp_v2.
c     call readqgcou() !no input. Read QGcou and store date into variables.
c     call Spacegrouprot(symgg,ngrp,plat,natom,pos) ! all inputs.

C         do iq=iqxini,iqxini
C         do ibz=1,200
C           if(nwgt(ibz,iq)/=0) then
C             write(6,"('yyy1: ',i8,2x,25(i3,i2))") ibz,(igx(i,ibz,iq),igxt(i,ibz,iq),i=1,nwgt(ibz,iq))
C           endif
C         enddo
C         enddo
      endif

!     ! == irkip control paralellization  ==
!     ! We have to distribute non-zero irkip into processes (nrank).
!     ! When irkip(nqibz,ngrp,nq,nspinmx)/=0, we expect grain-size 
!     ! for each job of (iqibz,igrp,iq,isp) is almost the same.
!     ! Our pupose is to calculate zsec(itp,itpp,iq).
!     ! Thus we need to set up communicator (grouping) MPI__COMM_iqisp(iq,isp) to do all_reduce.
!     ! (for given zsec(iq,isp), we take sum on zsec for (iqibz,igrp) by all_reduce.)
!     ! ---
!     ! NOTE: in future, we will further extend irkip for itp and itpp
      if(eibz4sig()) then
         allocate(irkip_all(nspinmx,nqibz,ngrp,nq)) !this is global 
         allocate(nrkip_all(nspinmx,nqibz,ngrp,nq)) !this is global 
         allocate(nrkip(nspinmx,nqibz,ngrp,nq)) !this is global 
         nrkip_all=0
         irkip_all=0
         is=1                   ! not spin dependent
         do iqq=1,nq
c     irkip_all(is,:,:,iqq)=irk
            do kx=1,nqibz
               do igrp=1,ngrp
                  kr = irk(kx,igrp) !ip_all(is,kx,igrp,iqq) !kr is index for qbz (for example, nonzero # of kr is 64 for 4x4x4)
                  if(kr==0) cycle
                  if(nwgt(kr,iqq)/=0) then
                     irkip_all(is,kx,igrp,iqq)= irk(kx,igrp)
                     nrkip_all(is,kx,igrp,iqq)= nwgt(kr,iqq)
                  endif   
c     print *,' iqq kr irk =',iqq,kr,irkip_all(is,kx,igrp,iqq),nrkip_all(is,kx,igrp,iqq)
               enddo
            enddo
         enddo  
C          do iqq=1,nq
C             write(6,"('iq=',i4,' # of EIBZ: Used(TimeR 1 or -1)=',i3,'=',i3,'+',i3)")
C      &           iqq,sum(eibzsym(:,:,iqq)),sum(eibzsym(:,1,iqq)),sum(eibzsym(:,-1,iqq))
C             write(6,"('eibz: iqq sum(nrkip_all)=nqbz  ',i3,3f11.5,3i8)")
C      &           iqq,q(:,iqq),sum(nrkip_all(is,:,:,iqq)),nqbz
C             do kx=1,nqibz
C                do igrp=1,ngrp
C                   kr = irkip_all(is,kx,igrp,iqq) !kr is index for qbz
C                   if(kr/=0) write(6,"('      ',i8,3f11.5,i8,2x,25(i4,i2))")
C      &                 kr,qbz(:,kr),nrkip_all(is,kx,igrp,iqq)
C      &                 ,(igx(i,kr,iqq),igxt(i,kr,iqq),i=1,nwgt(kr,iqq))
C                enddo
C             enddo
C!     !   Probably partial group symmetrization is enough. But it may not reduce computational time so much.
C         enddo  
         if(nspinmx==2) then
            irkip_all(2,:,:,:)=irkip_all(1,:,:,:)
            nrkip_all(2,:,:,:)=nrkip_all(1,:,:,:)
         endif
!     ! not eibz4sig
      else 
         do is = 1,nspinmx
            do iqq=1,nq
               irkip_all(is,:,:,iqq)=irk
            enddo
         enddo
      endif



cccccccccccccccccccccccccccccccccccccc
c      call MPI__finalize
c      stop 'test end xxxxxxxxxxxxxx'
ccccccccccccccccccccccccccccccccccccccc




cccccccccccccccccccccccccccccccccccccccc
c     do is = 1,nspinmx
c     do iqq=1,nq
c     irkip_all(is,:,:,iqq)=irk
c     nrkip_all(is,:,:,iqq)=1
c     enddo
c     enddo
c     nrkip=nrkip_all
cccccccccccccccccccccccccccccccccccccccc

      allocate(irkip(nspinmx,nqibz,ngrp,nq)) !local
      call MPI__sxcf_rankdivider(irkip,irkip_all,nspinmx,nqibz,ngrp,nq) ! MIZUHO-IR

!     ! nrkip is weight correspoinding to irkip for a node.
      nrkip = nrkip_all         !we don't need to change this for MPI case. It just need to distribute non-zero irkip.
!     ! ----------------------------------------
      nlnx4    = nlnx**4        ! niwx     = max0 (nw+1,niw) !nw --->nw+1 feb2006
      allocate( ppb(nlnmx*nlnmx*mdimx*nclass),
     &     kount(nqibz,nq),zsec(ntq,ntq,nq), coh(ntq,nq),nbandmx(nq,nspinmx) ) !eq(nband),
!TIME1 "before 2000 loop"
c     loop over spin ----------------------------------------------------
      do 2000 is = 1,nspinmx
!TIME0 
         if(mpi__root) then
         if(exchange) then
            write(ifsex2(is)) nspin, nq, ntq,nqbz,nqibz, n1,n2,n3
            write(ifsex(is),*) '======================================='
            write(ifsex(is),"('Self-energy exchange SEx(q,t): is=',i3)") is
            write(ifsex(is),*) '======================================='
            call winfo(ifsex(is),nspin,nq,ntq,is,nbloch,ngpn1,
     &           ngcn1,nqbz,nqibz,ef,deltaw,alat,esmr)
            write (ifsex(is),*)' *** '
            write (ifsex(is),"(a)") ' jband   iq ispin                  '//
     &           '           qvec            eigen-Ef (in eV)           exchange (in eV)'
         elseif(ixc==2) then
            write(ifsec2(is)) nspin, nq, ntq ,nqbz,nqibz  ,n1,n2,n3
            write(ifsec(is),*) '=========================================='
            write(ifsec(is),"('Self-energy correlated SEc(qt,w): is=',i3)") is
            write(ifsec(is),*) '=========================================='
            call winfo(ifsec(is),nspin,nq,ntq,is,nbloch,ngpn1,
     &           ngcn1,nqbz,nqibz,ef,deltaw,alat,esmr)
            write (ifsec(is),*)' *** '
            write (ifsec(is),"(a)") ' jband   iq ispin                  '//
     &           '           qvec            eigen-Ef (in eV)           '//
     &           'Re(Sc) 3-points (in eV)                        '//
     &           '           In(Sc) 3-points (in eV)                Zfactor(=1)'
         endif
         endif
         zsec  = 0d0
         coh   = 0d0
         kount = 0
         nbandmx(:,is) = 0 ! MIZUHO-IR  moved here for nctot=0 3rdNov2012
         if(ixc==3.and.nctot==0) goto 2001 !make dummy SEXcore

!TIME1 "before sxcf_fal3_scz"
!TIME0

!! dummy to overlaid -check bounds sep2014
         if(size(ecore)==0) then
            deallocate(ecore)
            allocate(ecore(1,2))
         endif

         call sxcf_fal3_scz (kount,q,itq,ntq,ef,esmr,
     i        tiat,miat,
c     v2
     i        nspin,is,   
c     i              ifrcw,ifrcwi,
     i        qbas,ginv,qibz,qbz ,wbz, nstbz, 
     i        nstar,irkip(is,:,:,:),nrkip(is,:,:,:),
c     
c     i       iiclass,nblocha,i_mnlv,i_mnlc,iicore,incore,iimdim,
     i        iclass,nblocha,nlnmv,nlnmc, !w(i_mnlv),w(i_mnlc)
     i        icore,ncore,imdim,
     i        ppb,              ! pdb,dpb,ddb,
     i        freqx,wwx,        !expa,ua, !sf 21May02 (no w(ifreq))
     i        dw,               !deltaw,
     i        ecore(:,is),
c     
     d        nlmto,nqibz,nqbz,nctot,
c     i        index_qbz, n_index_qbz, !12jan2004
     d        nl,nnc,nclass,natom,
     d        nlnmx,mdimx,nbloch,ngrp,niw,nq, !niwx removed in sxcf_fal2_scz2 feb2006
c     i     nblochpmx, ngpn,ngcni,ngpmx,ngcmx, !12jan2004
     i        nblochpmx, ngpmx,ngcmx,
c     i     geigB(1,1,1,is), ngvecpB,ngveccBr, !12jan2004
c     i     ngveccBr,
     i        wgt0,nq0i,q0i,symgg,alat,
     i        shtvg,nband,
     i        ifvcfpout,
c     i     shtw, 
     i        exchange, screen, cohtest, ifexsp(is),
c     i              omega,iSigMode, ! iwini,iwend,            
     i        iSigMode,         ! iwini,iwend,            sep2012takao
     i        nbmx,ebmx,        !takao 18June2003 
     i        pomatr, qrr,nnr,nor,nnmx,nomx,nkpo, !oct2005 for pomat 
     i        wklm,lxklm,       !mar2012
     i        invg,il,in,im,nlnm,nn, lx,nx,nxx,cgr,ppbrd, !sep 2012
     i      diagonly, eftrue, !diagonly nov2012, eftrue jan2013
     o        zsec,coh,nbandmx(1,is)) ! acuumulation variable
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c         print *,'---end of sxcf_fal3--- is mpi__rank',is,mpi__rank
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c         deallocate(wbzz)
!TIME1 "after sxcf_fal3_scz"
      call date_and_time(values=timevalues)
      write(*,'(a,9i5)')'dateandtime2=',mpi__rank,timevalues(1:8) 
      call cpu_time(time_red1)
!TIME0 

!! CAUITION! Allreduce wait all cpu jobs done here. 
!! Before nov2013, MPI__sxcf_rankdivider was stpid---> half of cores assigned for isp=2
!! was just waiting here!
      call MPI__AllreduceMax( nbandmx(:,is), nq ) ! MIZUHO-IR
      call cpu_time(time_red2)
!TIME1 "reduction nbandmx end"
      write(*,*) mpi__rank,'time(MPI__AllreduceMax)=',time_red2-time_red1

c---------------------------------
!TIME0
!!     electron gas bare exchange (exact)
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
!TIME1 "end legas exchange"
 2001    continue
!TIME0 

!     ! --- symmetrize zsec for eibz4sig mode. Sep2012 -----------------
!     !     Rotate zsec(itp,itp)
         allocate( zsect(ntq,ntq))


!! === readin lmfgw_kdivider, and get extensions === apr2013
         ifiproc=4001
         open(unit=ifiproc,file='lmfgw_kdivider',status='old')
         read(ifiproc,*) ext
         read(ifiproc,*) nqixx, nspxx, nrankv
         if(allocated(iprocq)) deallocate(iprocq)
         allocate(iprocq(nqixx,nspxx))
         do isp=1,nspxx
            do iqq=1,nqixx
               read(ifiproc,*) iqqxx, ispxx, ixxx
Cstop2rx 2013.08.09 kino               if(iqqxx/=iqq) stop 'iqqxx/=iqq'
               if(iqqxx/=iqq) call rx( 'iqqxx/=iqq')
Cstop2rx 2013.08.09 kino               if(ispxx/=isp) stop 'ispxx/=isp'
               if(ispxx/=isp) call rx( 'ispxx/=isp')
               iprocq(iqq,isp) = ixxx
               write(6,"('iqq isp irank=',i8,i2,i6)") iqq,isp, iprocq(iqq,isp)
            enddo
         enddo
         close(ifiproc)
!! for multiple files.
         if(allocated(extp)) deallocate(extp,ifvxc__,ifevec__)
         allocate(extp(0:nrankv-1),ifvxc__(0:nrankv-1),ifevec__(0:nrankv-1))
         extp(0) = trim(ext)
         write(6,"('  0 ext= ',a,a)") trim(extp(0)),' ----------'
         do procid=1,nrankv-1
            write(extn,"(i10)") procid
            extp(procid)=trim(adjustl(ext))//'_'//trim(adjustl(extn))
            write(6,"(i3,' ext= ',a,a)") procid,trim(extp(procid)),' ----------'
         enddo
         ixx=5001 
         do procid=0,nrankv-1
            ixx=ixx+1
            open(ixx,file='vxc'//extp(procid),form='unformatted')
            ifvxc__(procid)=ixx
            ixx=ixx+1
            open(ixx,file='evec'//extp(procid),form='unformatted')
            ifevec__(procid)=ixx
         enddo
!!-----------------------------------------------------------


         ifvxc_ = ifvxc__(0)  !0 is root
         ifevec_= ifevec__(0)
         read(ifevec_)  ndimhx, nspx,nnnx 
         read(ifvxc_)            !skip ndimh, nsp,nnn
         allocate(evaliq(nband),iblki(nband),iblke(nband))

!TIME1 "end of allocate zsect"
!TIME0
         iqq=0  !iqq is to read multiple vxc.* evec.*
         do 3020 iq=1,nq !nq means iq for which we will calculate sigma
            iqq=iqq+1
         do 3030 ispx=1,nspinmx !ispx loop is to find isx=is
            ifvxc_  = ifvxc__ (iprocq(iqq,ispx)) 
            ifevec_ = ifevec__(iprocq(iqq,ispx)) 
               if(ispx==is) then 
                 !this if-block is due to evec and v_xc file-->they shall be divieded into spin files.
                  read(ifvxc_)  nhdim,ldim
                  read(ifvxc_)
                  allocate( evec(nhdim,nhdim), evec_inv(nhdim,nhdim),evecrot(nhdim,nhdim))
                  read(ifevec_) qqqx(1:3), evec(1:nhdim,1:nhdim)
                  zsect = 0d0
               else             !skip isx/=is. Need to get access sequential files evec and v_xc.
                  read(ifvxc_)  
                  read(ifvxc_)  
                  read(ifevec_)
                  cycle
               endif  
               do i=1,nnnx  !nq     !qqqx from evec v_xc.
                  if(sum(abs(qqqx-q(:,i)))<1d-6) then
                     iqxx=i
                     goto 3011
                  endif  
               enddo  
               deallocate(evec,evec_inv,evecrot)
Cstop2rx 2013.08.09 kino               stop 'hsfp0_sc: bug:qqqx can not find ...'
               call rx( 'hsfp0_sc: bug:qqqx can not find ...')
 3011          continue

Cstop2rx 2013.08.09 kino               if(tiii) stop 'timereversal is not yet implemented'
               if(tiii) call rx( 'timereversal is not yet implemented')
               call readeval(q(:,iqxx), is, evaliq)
               evec_inv = evec
               call matcinv(nhdim,evec_inv)
               nsym = sum(eibzsym(:,:,iqxx))

               do it=1,1        !no-time reversal yet !it=1,-1,-2 !c.f. x0kf_v4h
                  do igrp=1,ngrp !A-rotator
                     if( eibzsym(igrp,it,iqxx)==0) cycle
                     nblk=0
                     iblki=0
                     iblke=0
                     iblki(1)=1
!!  degeneracy divider for evaliq. See How to apply EIBZ to 
!! Is this procedure really make speed up so much?
                     tolry= 0.2d0 !Degeneracy tol. if tolry is large, 
!! larger tolry is safer, although a little inefficient.
!! If tolry is too small to divide degenerated values to different blocks --> then we have wrong results.
                                !(NOTE that Hamiltonian can be not so symmetric in some reasons) 
                     nbsizemx=0
                     do iband=2,nbandmx(iqxx,is) 
                                ! nbandmx is the number of bands for which we calculate self-energy.
                                ! We assume nbandmx(iqxx,is) is well separated for degeneracy.
                        if(evaliq(iband) > evaliq(iband-1)+tolry 
     &                       .or.iband==nbandmx(iqxx,is)) then
                           nblk=nblk+1
                           if(nblk>=2) iblki(nblk)=iblke(nblk-1)+1
                           if(iband==nbandmx(iqxx,is)) then
                              iblke(nblk)=iband
                           else
                              iblke(nblk)=iband-1
                           endif  
                           nbsize = iblke(nblk)- iblki(nblk)+1
                           if( nbsize>nbsizemx ) nbsizemx = nbsize
                        endif
                     enddo      ! iband
!     ! rotation of evec. Generate evecrot. (Within degenerated block, evec are mapped).e
                     allocate(rmatjj(nbsizemx,nbsizemx,nblk))
                     napw=nhdim-ldim
                     do iblk1=1,nblk  
                        ii1=iblki(iblk1)
                        ie1=iblke(iblk1) 
                        ne1=ie1-ii1+1
                        call rotwvigg(igrp,q(:,iqxx),q(:,iqxx),nhdim,
     &                       napw,ne1,evec(:,ii1:ie1),evecrot(:,ii1:ie1),ierr )
                        rmatjj(1:ne1,1:ne1,iblk1) = 
     &                       matmul(evec_inv(ii1:ie1,:),evecrot(:,ii1:ie1))
                     enddo      ! iblk1
                     do iblk1=1,nblk  
                        do iblk2=1,nblk  
                           ii1=iblki(iblk1)
                           ie1=iblke(iblk1) 
                           ne1=ie1-ii1+1
                           ii2=iblki(iblk2)
                           ie2=iblke(iblk2) 
                           ne2=ie2-ii2+1
                           zsect(ii1:ie1,ii2:ie2)= zsect(ii1:ie1,ii2:ie2) 
     &                          + matmul( dconjg(transpose(rmatjj(1:ne1,1:ne1,iblk1))), 
     &                          matmul(zsec(ii1:ie1,ii2:ie2,iqxx),
     &                          rmatjj(1:ne2,1:ne2,iblk2)) )
                        enddo   ! iblk2
                     enddo      ! iblk1
                     deallocate(rmatjj)
                  enddo         ! igrp
               enddo            ! it
               deallocate(evec, evec_inv, evecrot)
               zsec(:,:,iqxx) = zsect(:,:)/dble(nsym)
               call MPI__AllreduceSum( zsec(:,:,iqxx),ntq*ntq ) ! MIZUHO-IR
 3030       continue            ! ispx
 3020    continue               ! iq
!TIME1 "end 3020 loop"

      ixx=5001
      do procid=0,nrankv-1
        ixx=ixx+1
        close(unit=ixx)
        ixx=ixx+1
        close(unit=ixx)
      enddo
      deallocate(iblki,iblke)
      deallocate(evaliq)
      deallocate(zsect)

!TIME0 
         if(mpi__root) then
         if(exchange) then
            ifoutsex=ifsex(is)
            write(6,*)
            do ip = 1,nq
               do i  = 1,ntq
                  write(ifoutsex,"(3i5,3d24.16,3x,d24.16,3x,d24.16)")
     &                 itq(i),ip,is, q(1:3,ip), eqx(i,ip,is), 
     &                 hartree*dreal(zsec(i,i,ip)) !sf 21May02
                  if( eqx(i,ip,is)<1d20.and.abs(zsec(i,i,ip))/=0d0 ) then !takao june2009
                     write(6,"(' j iq isp=' i3,i4,i2,'  q=',3f8.4,
     &                    '  eig=',f10.4,'  Sx=',f10.4)")
     &                    itq(i),ip,is, q(1:3,ip), eqx(i,ip,is),
     &                    hartree*dreal(zsec(i,i,ip)) !sf 21May02
                  endif
               end do
               write(ifsex2(is)) is, q(1:3,ip), zsec(1:ntq,1:ntq,ip) !SEC_nn' out
            end do
         elseif(ixc==2) then
            ifoutsec=ifsec(is)
            do ip = 1,nq
               do i  = 1,ntq
                  if( eqx(i,ip,is)<1d20.and.abs(zsec(i,i,ip))/=0d0 ) then !takao june2009
                     write(6,"(' j iq isp=' i3,i4,i2,'  q=',3f8.4,
     &                    '  eig=',f8.4,'  Re(Sc) =',f8.4,'  Img(Sc) =',f8.4 )") !sf 21May02
     &                    itq(i),ip,is, q(1:3,ip), eqx(i,ip,is), 
     &                    hartree*dreal(zsec(i,i,ip)), !sf 21May02
     &                    hartree*dimag(zsec(i,i,ip)) !,zfac,zfac1,zfac2  !sf 21May02
                  endif
                  write(ifoutsec,"(3i5,3d24.16,3x,d24.16,3x,d24.16, 3x,d24.16)")
     &                 itq(i),ip,is, q(1:3,ip), eqx(i,ip,is), 
     &                 hartree*dreal(zsec(i,i,ip)), !sf 21May02
     &                 hartree*dimag(zsec(i,i,ip)) !,zfac,zfac1,zfac2  !sf 21May02
               end do
               write(ifsec2(is)) is, q(1:3,ip), zsec(1:ntq,1:ntq,ip) !SEC_nn' out
            end do
         endif                  !ixc
         endif  !mpi__root
!TIME1 "before 2000 continue"
 2000 continue                  !end of spin-loop

c$$$!TIME0 
c$$$!!   --- EXspectrum ----------------------------------------------------
c$$$c     This section is similar with efsimplef.f
c$$$      if(sum(ifexsp(1:nspin))/=0) then
c$$$         do is = 1,nspin
c$$$            print *,' --- Goto ExSpectrum section --- is=',is
c$$$            rewind (ifexsp(is))
c$$$            itmx = 0
c$$$            do
c$$$               read(ifexsp(is),*,end=1215)ipex,itpex,itex,qqex(1:3), eex,exsp
c$$$               if(itex>itmx) itmx=itex
c$$$            enddo
c$$$ 1215       continue
c$$$            nspexmx = itmx*(nqbz+nq0i*ngrp) !Get marimum value of the number of the ex spectrum
c$$$c     
c$$$            allocate( eex1(nspexmx,ntq,nq), exsp1(nspexmx,ntq,nq),
c$$$     &           nspex(ntq,nq) ,
c$$$     &           itex1(nspexmx,ntq,nq),
c$$$     &           qqex1(3,nspexmx,ntq,nq) )
c$$$            print *,' nspexmx =',nspexmx
c$$$c     
c$$$            rewind (ifexsp(is))
c$$$            nspex = 0
c$$$            do
c$$$               read(ifexsp(is),*,end=1216) ipex,itpex,itex,qqex(1:3),eex,exsp
c$$$               nspex(itpex,ipex) = nspex(itpex,ipex)+1
c$$$               iex  = nspex(itpex,ipex)
c$$$               eex1  (iex,itpex,ipex) = eex
c$$$               exsp1 (iex,itpex,ipex) = exsp
c$$$               itex1 (iex,itpex,ipex) = itex
c$$$               qqex1(:,iex,itpex,ipex)= qqex
c$$$            enddo
c$$$ 1216       continue            !Get eex1(1:nspex) exsp1(1:nspex) for itp ip.
c$$$            print *,' nspex(1 1)=',nspex(1,1)
c$$$c     
c$$$            do ipex = 1,nq
c$$$               do itpex=1,ntq
c$$$                  write(6,*)' is itq ip =',is,itq,ip
c$$$                  nnex = nspex(itpex,ipex)
c$$$                  allocate( ieord(1:nnex) )
c$$$                  call sortea( eex1(1:nnex,itpex,ipex),ieord, nnex,isig)
c$$$                  eex1 (1:nnex,itpex,ipex)  = eex1  (ieord(1:nnex),itpex,ipex)
c$$$                  exsp1 (1:nnex,itpex,ipex) = exsp1 (ieord(1:nnex),itpex,ipex)
c$$$                  itex1 (1:nnex,itpex,ipex) = itex1 (ieord(1:nnex),itpex,ipex)
c$$$                  qqex1(:,1:nnex,itpex,ipex)= qqex1 (:,ieord(1:nnex),itpex,ipex)
c$$$
c$$$                  filenameex = 'EXSP'//charnum3(ipex)//charnum3(itpex)
c$$$     &                 //'.'//char(48+is)
c$$$                  ifexspx=4111
c$$$                  open(ifexspx,file=filenameex)
c$$$
c$$$                  filenameex = 'EXSS'//charnum3(ipex)//charnum3(itpex)
c$$$     &                 //'.'//char(48+is)
c$$$                  ifexspxx=4112
c$$$                  open(ifexspxx,file=filenameex)
c$$$
c$$$                  do i=1,nnex
c$$$                     write(ifexspx, "(2d14.6, i4, 3f14.6)")
c$$$     &                    eex1  (i,itpex,ipex), exsp1 (i,itpex,ipex),
c$$$     &                    itex1 (i,itpex,ipex), qqex1 (1:3,i,itpex,ipex)
c$$$                  enddo
c$$$c     
c$$$                  eee  =-1d99
c$$$                  exwgt= 0d0
c$$$                  do i=1,nnex
c$$$                     if(eex1(i,itpex,ipex) > eee+1d-4 .or. i==nnex) then
c$$$                        if(i/=1) write(ifexspxx, "(2d23.15)")
c$$$     &                       eee, exwgt*hartree
c$$$                        eee  = eex1(i,itpex,ipex)
c$$$                        exwgt= exsp1 (i,itpex,ipex)
c$$$                     else
c$$$                        exwgt= exwgt + exsp1 (i,itpex,ipex)
c$$$                     endif
c$$$                  enddo
c$$$c     
c$$$                  deallocate( ieord )
c$$$                  close(ifexspx)
c$$$                  close(ifexspxx)
c$$$               enddo
c$$$            enddo
c$$$            deallocate( eex1, exsp1, nspex, itex1, qqex1 )
c$$$         enddo
c$$$         print *,' End of ExSpectrum section ---'
c$$$      endif
c$$$!TIME1 "end of exspectrum"
!TIME1 "end program"
!TIMESHOW 
      isx = iclose ('wc.d')
      isx = iclose ('wci.d')
      isx = iclose ('hbe.d')
      call cputid(0)
      call MPI__Finalize 
      if(ixc==1 ) call rx0( ' OK! hsfp0_sc: Exchange mode')
      if(ixc==2 ) call rx0( ' OK! hsfp0_sc: Correlation mode')
      if(ixc==3 ) call rx0( ' OK! hsfp0_sc: Core-exchange mode')
      end program hsfp0_sc 

