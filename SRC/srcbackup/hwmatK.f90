program hwmat
  !!-------------------------------------------------------------
  !! Calculates the bare/screened interaction W
  !! See note in hmatK_MPI.F
  !!------------------------------------------------------------
  use m_readqg,only: Readngmx2,ngcmx,ngpmx,Readqg0,Readqg
  use m_hamindex,only:   Readhamindex,ngrp,symgg=>symops,invg=>invgx
  use m_read_bzdata,only: Read_bzdata,qibz,irk,ginv,n1,n2,n3,nqbz,nqibz,nstar,nstbz,qbas=>qlat,qbz,wibz,wbz
  &  ,  nq0i=>nq0ix,q0i,wqt=>wt
  use m_readeigen,only: onoff_write_pkm4crpa,init_readeigen &
       ,init_readeigen2, init_readeigen_mlw_noeval, init_readeigen_phi_noeval,  nwf
  use m_genallcf_v3,nctot2=>nctot,niwg=>niw
  use m_keyvalue,only: getkeyvalue
  use m_readhbe,only: Readhbe, nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg
  use m_zmel,only: ppbafp_v2

  implicit none
  real(8),parameter :: &
       ua    = 1d0    ! constant in w(0)exp(-ua^2*w'^2) to take care of peak around w'=0
  !------------------------------------
  real(8)    :: esmr2,shtw
  !      integer(4) :: mxclass,ngnmax,mbytes,mwords,iwksize,
  !     &   natom,nclass,ipos,igrp,
  !     &   iqibz,
  !     &   iqbz,
  !     &   iinvg,
  !     o   nspin,nl,nn,nnv,nnc,
  !     o   inindx,inindxv,inindxc,iiclass,
  !     d   nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc,
  !     o   iz,
  !     o   iil,iin,iim,iilnm,i_mnl,
  !     o   iilv,iinv,iimv,iilnmv,i_mnlv,
  !     o   iilc,iinc,iimc,iilnmc,i_mnlc,
  !     o   incwf,iecore,ikonf,iicore,incore,nctot,
  !     o   imagw,niw,nw,ifreq,
  integer(4):: &
       ixc,iopen, ! & ifhbed, nprecb,mrecb,mrece,nlmtot,nqbzt, nband,
  ibas,ibasx,nxx,nbloch,ifqpnt,ifwd,ifmloc, &
       nprecx,mrecl,nblochpmx2,nwp,niwt, nqnum,mdimx,nblochpmx, &
       ifrcw,ifrcwi,  noccxv,maxocc2,noccx,ifvcfpout,iqall,iaf,ntq, &
       i,k,nspinmx, nq,is,ip,iq,idxk,ifoutsex,iclose,ig, &
       mxkp,nqibzxx,ntet,nene,iqi, ix,iw, &
       nlnx4,niwx,irot,invr,invrot,ivsum, ifoutsec,ntqx, &
       
       !     &   ifrb(2),ifcb(2),ifrhb(2),ifchb(2)
       !     &    ifev(2),
       ifwmat(2) ! & ,ifcphi
  ,ifxc(2),ifsex(2), ifphiv(2),ifphic(2),ifec,ifexsp(2), &
       ifsecomg(2),ifexx,ndble=8

  !      real(8) :: alat,ef,diw,dw,delta,pi,tpia,vol,voltot,rs,alpha,
  real(8) :: pi,tpia,vol,voltot,rs,alpha, &
       qfermi,efx,valn,efnew,edummy,efz,qm,xsex,egex, &
       zfac1,zfac2,dscdw1,dscdw2,dscdw,zfac,ef2=1d99,exx,exxq,exxelgas

  !     &   lowesteval !defined in readeigen
  ! c   qbas(3,3),ginv(3,3) plat(3,3),
  logical :: lqall,laf
  !      character*120 symgrp

  ! class parameters
  !      parameter (mxclass=100)
  !      character*6 clabl(mxclass)
  ! symmetry group
  !      parameter (ngnmax=10)
  !      real(8) :: gen(9,ngnmax)

  integer(4),allocatable :: itq(:)
  real(8),allocatable    :: q(:,:)

  ! takao
  integer(4),allocatable :: ngvecpB(:,:,:),! & ngveccB(:,:,:),
  ngvecp(:,:), ngvecc(:,:),iqib(:), ! & ,ngpn(:)ngcni(:)
  kount(:), nx(:,:),nblocha(:),lx(:) !ngveccBr(:,:,:)
  real(8),allocatable:: vxcfp(:,:,:), &
       wgt0(:,:), &
       ppbrd (:,:,:,:,:,:,:),cgr(:,:,:,:),eqt(:), &
       ppbrdx(:,:,:,:,:,:,:),aaa(:,:), ! & symope(:,:,:)=symgg, ! qibz(:,:),
  ppb(:), eq(:), ! & ,pdb(:),dpb(:),ddb(:)
  eqx(:,:,:),eqx0(:,:,:),ekc(:),coh(:,:) &
       , rw_w(:,:,:,:,:,:),cw_w(:,:,:,:,:,:), &
       rw_iw(:,:,:,:,:,:),cw_iw(:,:,:,:,:,:)
  complex(8),allocatable:: geigB(:,:,:,:)

  logical :: screen, exchange, cohtest, legas, tote, lueff
  real(8) ::  rydberg,hartree
  real(8):: qreal(3), ntot,nocctotg2,tripl!,xxx(3,3)
  logical ::nocore

  ! space group infermation
  integer(4),allocatable :: iclasst(:), invgx(:), miat(:,:)
  real(8),allocatable    :: tiat(:,:,:),shtvg(:,:)

  ! tetra
  real(8),allocatable :: qz(:,:),qbzxx(:),wbzxx(:),wtet(:,:,:,:), &
       eband(:,:,:), ene(:) !,ecore(:,:)
  integer(4),allocatable ::idtetx(:,:),idtet(:,:),ipq(:) &
       ,iene(:,:,:),ibzx(:) ! ,nstar(:)
  !      real(8) :: qbasmc(3,3)

  ! worksize in megabytes (1 word = 4 bytes)
  !      integer(4) :: mbytes,mwords,iwksize
  !      parameter (mbytes=60)
  !      parameter (mwords=mbytes/4)
  !      parameter (iwksize=mwords * 1000 *1000)
  !      integer w
  !      common /w/ w(iwksize)

  integer(4) ::ib,iqx,igp,iii,ivsumxxx,isx,iflegas, iqpntnum

  real(8),allocatable   :: eex1(:,:,:),exsp1(:,:,:),qqex1(:,:,:,:)
  integer(4),allocatable:: nspex(:,:),ieord(:),itex1(:,:,:)
  real(8)    :: qqex(1:3), eex,exsp,eee, exwgt,deltax0
  integer(4) :: itmx,ipex,itpex,itex,nspexmx,nnex,isig,iex,ifexspx &
       ,ifexspxx ,ifefsm, ifemesh,nz
  character(3)  :: charnum3,sss
  character(12) :: filenameex
  logical :: exspwrite=.false.
  character(8) :: xt

  integer(4) :: iwini,iwend
  real(8),allocatable:: omega(:,:)
  real(8) ::  omegamax,dwplot,omegamaxin
  !      logical :: sergeys

  integer(4)::nqbze,ini,nq0it,idummy
  real(8),allocatable:: qbze(:,:)

  real(8)   :: ebmx(2)
  integer(4):: nbmx(2)

  real(8):: volwgt

  integer(4)::nwin, incwfin
  real(8)::efin,ddw,dwdummy
  integer(4),allocatable::imdim(:)
  real(8),allocatable::freqx(:),freqw(:),wwx(:),expa(:)

  logical:: GaussSmear !readgwinput,
  integer(4)::ret
  character*(150):: ddd


  integer(4)::   ngpn1,verbose,ngcn1,nwxx !bzcase,mrecg,
  real(8)   :: wgtq0p,quu(3)

  real(8),allocatable:: freq_r(:)

  logical ::smbasis
  integer(4):: ifpomat,nkpo,nnmx,nomx,ikpo,nn_,no,nss(2)
  real(8):: q_r(3)
  real(8),allocatable:: qrr(:,:)
  integer(4),allocatable:: nnr(:),nor(:)
  complex(8),allocatable:: pomatr(:,:,:),pomat(:,:)

  real(8)::sumimg
  logical :: allq0i                                             !S.F.Jan06

  integer(4):: nw_i,ifile_handle,if101,if3111

  logical:: latomic,lfull,lstatic,lwssc
  logical:: l1d
  real(8):: rsite(3),rcut1,rcut2,ef
  real(8),allocatable :: rws(:,:),drws(:),rws1(:,:),rws2(:,:)
  integer(4):: nrws,nrws1,nrws2,ir1,ir2,ir3,ir,nrw,ia
  integer(4),allocatable:: irws(:),irws1(:),irws2(:)
  integer:: nw, nctot0,niw,if102
  !---------------------------------------
  hartree=2d0*rydberg()
  iii=verbose()
  write(6,*)' verbose=',iii

  ! mode switch. --------------
  write(6,*) ' --- Choose omodes below ----------------'
  write(6,*) '  V (1) or W (2) or U(3)'
  write(6,*) '  [option --- (+ QPNT.{number} ?)] '
  write(6,*) ' --- Put number above ! -----------------'
  call readin5(ixc,nz,idummy)
  write(6,*) ' ixc nz=',ixc, nz
  if(ixc==0) stop ' --- ixc=0 --- Choose computational mode!'

  !---  readin BZDATA. See gwsrc/rwbzdata.f
  !--------readin data set when you call read_BZDATA ---------------
  !       integer(4)::ngrp,nqbz,nqibz,nqbzw,nteti,ntetf
  ! ccc    ! &   ,n_index_qbz
  !       integer(4):: n1,n2,n3
  !       real(8):: qbas(3,3),ginv(3,3),qbasmc(3,3)
  !       real(8),allocatable:: qbz(:,:),wbz(:),qibz(:,:)
  !     &    ,wibz(:),qbzw(:,:)
  !       integer(4),allocatable:: idtetf(:,:),ib1bz(:),idteti(:,:)
  !     &    ,nstar(:),irk(:,:),nstbz(:)          !,index_qbz(:,:,:)
  !-----------------------------------------------------------------
  call read_BZDATA()
  write(6,*)' nqbz  =',nqbz
  !      write(6,*)  qbz
  write(6,*)' nqibz ngrp=',nqibz,ngrp
  !      write(6,*)' irk=',irk
  !      write(6,*)' #### idtetf: ####'
  !      write(6,*) idtetf

  ! set up work array
  !      call wkinit (iwksize)
  call pshprt(60)


  !--- readin GWIN and LMTO, then allocate and set datas.
  !      nwin =-999    !not readin NW file
  !      efin =-999d0  !not readin EFERMI
  incwfin= -1  !use 7th colmn for core at the end section of GWIN
  call genallcf_v3(incwfin) !in module m_genallcf_v3
  niw=niwg
  !      if(ngrp/= ngrp2) stop 'ngrp inconsistent: BZDATA and LMTO GWIN_V2'
  ef=1d99
  !---  These are allocated and setted.
  !      integer(4)::  nclass,natom,nspin,nl,nn,nnv,nnc, ngrp,
  !     o  nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, nctot,niw, !not readin nw
  !      real(8) :: alat,ef, diw,dw,delta,deltaw,esmr
  !      character(120):: symgrp
  !      character(6),allocatable :: clabl(:)
  !      integer(4),allocatable:: iclass(:)
  !     &  ,nindxv(:,:),nindxc(:,:),ncwf(:,:,:) ,
  !     o    invg(:), il(:,:), in(:,:), im(:,:),   ilnm(:),  nlnm(:),
  !     o    ilv(:),inv(:),imv(:),  ilnmv(:), nlnmv(:),
  !     o    ilc(:),inc(:),imc(:),  ilnmc(:), nlnmc(:),
  !     o    nindx(:,:),konf(:,:),icore(:,:),ncore(:),
  !     &    occv(:,:,:),unoccv(:,:,:)
  !     &   ,occc(:,:,:),unoccc(:,:,:),
  !     o    nocc(:,:,:),nunocc(:,:,:)
  !      real(8), allocatable::
  !     o  plat(:,:),pos(:,:),z(:),  ecore(:,:),  symgg(:,:,:) ! symgg=w(igrp),freq(:)
  !-----------------------------------------------------------------------

  !--- Get maximums takao 18June03
  !      call getnemx(nbmx,ebmx,8,.true.) !8+1 th line of GWIN0
  call getnemx8(nbmx,ebmx)
  write(6,"('  nbmx ebmx from GWinput=',2i8,2d13.5)") nbmx,ebmx

  !-------------------------------------------------------------------
  !      if (nclass > mxclass) stop ' hsfp0: increase mxclass'
!!!!! WE ASSUME iclass(iatom)= iatom !!!!!!!!!!!!!!!!!!!!!!!!!
  if (nclass /= natom ) stop ' hsfp0: nclass /= natom ' ! We assume nclass = natom.
  write(6,*)' hsfp0: end of genallcf2'

  call pshprt(30)
  pi   = 4d0*datan(1d0)
  tpia = 2d0*pi/alat

  shtw = 0d0
  if(esmr<1d-5) shtw=0.01d0 ! Ferdi's shift to avoid resonance effect(maybe)

  !      call dinv33(plat,1,xxx,vol)
  !      voltot = dabs(vol)*(alat**3)
  voltot = abs(alat**3*tripl(plat,plat(1,2),plat(1,3)))

  legas = .false.
  latomic = .false.
  l1d = .false.
  inquire(file='UU1dU',exist=l1d)

  if (ixc==1) then
     write(6,*)' --- bare Coulomb mode --- '
     exchange =.true.
     lueff = .false.
     ifwmat(1) = iopen('VMATU',1,-1,0)
     if (nspin == 2) ifwmat(2) = iopen('VMATD',1,-1,0)
  elseif (ixc ==2) then
     write(6,*)' --- screening (Wc) mode --- '
     exchange =.false.
     lueff = .false.
     ifwmat(1) = iopen('WcMATU',1,-1,0)
     if (nspin == 2) ifwmat(2) = iopen('WcMATD',1,-1,0)
  elseif (ixc ==3) then
     write(6,*)' --- Ueff mode --- '
     exchange =.false.
     lueff = .true.
     ifwmat(1) = iopen('UMATU',1,-1,0)
     if (nspin == 2) ifwmat(2) = iopen('UMATD',1,-1,0)
  elseif (ixc==4) then
     write(6,*)' --- bare Coulomb mode --- '
     exchange =.true.
     lueff = .false.
     latomic = .true.
     ixc=1
     ifwmat(1) = iopen('VMATaU',1,-1,0)
     if (nspin == 2) ifwmat(2) = iopen('VMATaD',1,-1,0)
  elseif (ixc ==5) then
     write(6,*)' --- screening (Wc) mode --- '
     exchange =.false.
     lueff = .false.
     latomic = .true.
     ixc=2
     ifwmat(1) = iopen('WcMATaU',1,-1,0)
     if (nspin == 2) ifwmat(2) = iopen('WcMATaD',1,-1,0)
  elseif (ixc ==6) then
     write(6,*)' --- Ueff mode --- '
     exchange =.false.
     lueff = .true.
     latomic = .true.
     ixc=3
     ifwmat(1) = iopen('UMATaU',1,-1,0)
     if (nspin == 2) ifwmat(2) = iopen('UMATaD',1,-1,0)
  else
     stop 'ixc error'
  endif

  !---
  write(6, *) ' --- computational conditions --- '
  write(6,'("    deltaw  =",f13.6)') deltaw
  write(6,'("    ua      =",f13.6)') ua
  write(6,'("    esmr    =",f13.6)') esmr
  write(6,'("    alat voltot =",2f13.6)') alat, voltot
  !      write(6,'("    niw nw dw   =",2i3,f13.6)') niw,nw,dw

  !>> read dimensions of wc,b,hb
  !      ifhbed     = iopen('hbe.d',1,0,0)
  !      read (ifhbed,*) nprecb,mrecb,mrece,nlmtot,nqbzt, nband,mrecg
  !      if (nprecb == 4) stop 'hsfp0: b,hb in single precision'
  call Readhbe()
  call Readhamindex()
  call init_readeigen()     !nband,mrece) !initialization of readEigen

  ! --- get space group information ---------------------------------
  ! true class information in order to determine the space group -----------
  ! because the class in the generated GW file is dummy.(iclass(ibas)=ibas should be kept).
  if102=ifile_handle()
  open (if102,file='CLASS')
  allocate(iclasst(natom),invgx(ngrp) &
       ,miat(natom,ngrp),tiat(3,natom,ngrp),shtvg(3,ngrp))
  write(6,*)'  --- Readingin CLASS info ---'
  do ibas = 1,natom
     read(if102,*) ibasx, iclasst(ibas)
     write(6, "(2i10)") ibasx, iclasst(ibas)
  enddo
  close(if102)

  ! Get space-group transformation information. See header of mptaouof.
  call mptauof(symgg,ngrp,plat,natom,pos,iclasst &
       ,miat,tiat,invgx,shtvg )
  !        write (*,*)  'tiat=', tiat(1:3,1:natom,invr),invr

  ! Get array size to call rdpp
  !      call getsrdpp( nclass,nl,
  !     o               ngpmx,ngcmx,nxx )
  call getsrdpp2( nclass,nl,nxx)
  !      call readngmx('QGpsi',ngpmx)
  !      call readngmx('QGcou',ngcmx)
  call readngmx2()
  write(6,*)' ngcmx ngpmx=',ngcmx,ngpmx
  allocate( nx(0:2*(nl-1),nclass), nblocha(nclass) ,lx(nclass), &
       ppbrd ( 0:nl-1, nn, 0:nl-1,nn, 0:2*(nl-1),nxx, nspin*nclass), &
       cgr(nl**2,nl**2,(2*nl-1)**2,ngrp))

  !- readin plane wave parts, and Radial integrals ppbrd.
  ! ppbrd = radial integrals
  ! cgr   = rotated cg coeffecients.
  ! geigB = eigenfunction's coefficiens for planewave.
  ! ngvecpB (in 1stBZ) contains G vector for eigen function.
  ! ngveccB (in IBZ)   contains G vector for Coulomb matrix.
  !      call rdpp_v2( ngpmx,ngcmx,nxx,  qibz,nqibz, qbz,nqbz,
  !     i      nband, nl,ngrp, nn,  nclass, nspin, symgg,    qbas,
  !     o      nblocha, lx, nx, ppbrd ,
  !     o      mdimx, nbloch, cgr,
  !     o      nblochpmx, ngpn,geigB,ngvecpB,  ngcni,ngveccB )
  call rdpp_v3(nxx, nl, ngrp, nn, nclass, nspin, symgg,qbas, &
       nblocha, lx, nx,  ppbrd , mdimx, nbloch, cgr)
  !      nblochpmx = nbloch + ngcmx

  !      allocate(ngcni(nqibz)) !, ngveccB(3,ngcmx,nqibz)) !, ngveccBr(3,ngcmx,nqibz))
  !   geigB(ngpmx,nband,nqbz,nspin),ngpn(nqbz),ngvecpB(3,ngpmx,nqbz),
  !     &   )  ! in IBZ

  !      call rdpp_pln(ngpmx,ngcmx, qibz,nqibz, qbz,nqbz,nband,nspin,
  !     o      ngpn,geigB,ngvecpB,ngcni,ngveccB)

  !      do iq = 1,nqibz
  !        call readqg('QGcou',qibz(1:3,iq),ginv,  quu,ngcni(iq), ngveccB(1,1,iq))
  !        write(6,"('--From QGcou  qibz quu ngc=',3f9.4,'  ',3f9.4,i5)")
  !     &       qibz(1:3,iq),quu,ngcni(iq)
  !      enddo
  ! for info
  allocate(ngvecp(3,ngpmx),ngvecc(3,ngcmx))
  call readqg('QGpsi',qibz(1:3,1), quu,ngpn1, ngvecp)
  call readqg('QGcou',qibz(1:3,1), quu,ngcn1, ngvecc)
  deallocate(ngvecp,ngvecc)

  write(6,*) ' end of read QGcou'

  ! cccccccccccccccccccccccccccccccccccccccccccccccccc
  !      iqx=1
  !      do ib=1,nband
  !        write(6,'("  iband iqx sumgeigB=",2i3,12d12.3)')
  !     &   ib,iqx, sum(geigB(1:ngpn(1),ib,iqx))
  !      enddo
  ! ccccccccccccccccccccccccccccccccccccccccccccccccccc
  !      ib =1
  !      iqx=1
  !      do igp =1,nbloch + ngpn(iqx)
  !         write(6,'("  igb ib iqx geigB=",3i3,1x,3i2,12d12.3)')
  !     &   igp,ib,iqx, ngvecpB(1:3,igp,iqx), geigB(igp,ib,iqx)
  !      enddo
  !      stop "xxxxxx zzz"
  ! cccccccccccccccccccccccccccccccccccccccccccccccccc
  !----------------------------------------------
  call pshprt(60)


  !--- Readin WV.d
  if ( .NOT. exchange) then
     if (lueff) then
        ifwd      = iopen('WV.d.maxloc',1,-1,0)
     else
        ifwd      = iopen('WV.d',1,-1,0)
     endif
     read (ifwd,*) nprecx,mrecl,nblochpmx,nwp,niwt,nqnum,nw_i
     ! nblochpmx from WV.d oct2005
     write(6,"(' Readin WV.d =', 10i5)") &
          nprecx, mrecl, nblochpmx, nwp, niwt, nqnum, nw_i
     call checkeq(nprecx,ndble)
     !       call checkeq(nblochpmx,nblochpmx2)
     !       if (nwt /= nw)   stop 'hwmat: wrong nw'
     !       nw = nwt
     nw=nwp-1
     if (niwt /= niw) stop 'hwmat: wrong niw'
     ! m 050518
     niw = 0
     niwt = 0
     if (lueff) then
        ifrcw     = iopen('WVR.maxloc',0,-1,mrecl)
        !          ifrcwi = iopen('WVI.maxloc',0,-1,mrecl)
     else
        !          ifrcw     = iopen('WVR',0,-1,mrecl)
        !          ifrcwi = iopen('WVI',0,-1,mrecl)
     endif
     !... reading general energy mesh from file 'freq_r'
     !        open(UNIT=3111,file='freq_r') !this is in a.u.
     !        read(3111,*)nwxx             !number of energy points
     !        if(nwxx/=nw) stop ' freq_r nw /=nw'
     !        allocate(freq_r(nw))       !freq_r(1)=0d0
     !        read(3111,*)freq_r
     !     close(3111)
     if3111=ifile_handle()
     open(UNIT=if3111,file='freq_r') !this is in a.u.
     read(if3111,*)nwxx             !number of energy points
     if(nwxx/= nw+1) stop ' freq_r nw /=nw'
     allocate(freq_r(nw_i:nw))       !freq_r(1)=0d0
     do iw = nw_i,nw
        read(if3111,*)freq_r(iw)
     enddo
     close(if3111)

  else
     ifvcfpout = iopen('VCCFP',0,-1,0)
     allocate(freq_r(1)); freq_r=1d99
     !       nw = 1
     nw = 0
  endif

  nrw = nw
  call getkeyvalue("GWinput","wmat_static",lstatic,default=.false.)
  if (lstatic) nrw = 0

  !... Readin eigen functions
  !      ifev(1)   = iopen('EVU', 0,0,mrece)
  !      if (nspin==2) ifev(2) = iopen('EVD', 0,0,mrece)

  !$$$c --- determine Fermi energy ef for given valn (legas case) or corresponding charge given by z and konf.
  !$$$! When esmr is negative, esmr is geven automatically by efsimplef.
  !$$$        write(6,*)'nnn1 nband=',nband
  !$$$        call efsimplef2a(nspin,wibz,qibz,ginv,
  !$$$     i        nband,nqibz
  !$$$     i       ,konf,z,nl,natom,iclass,nclass
  !$$$     i       ,valn, legas, esmr,  !!! valn is input for legas=T, output otherwise.
  !$$$c
  !$$$     i        qbz,nqbz !index_qbz, n_index_qbz,
  !$$$     o       ,efnew)
  !$$$        write(6,*)'nnn2 nband=',nband
  !$$$c
  !$$$c        write(6,*)' end of efsimple'
  !$$$c        ef = efnew
  !$$$c- check total ele number -------
  !$$$        ntot  = nocctotg2(nspin, ef,esmr, qbz,wbz, nband,nqbz) !wbz
  !$$$        write(6,*)' ef    =',ef
  !$$$        write(6,*)' esmr  =',esmr
  !$$$        write(6,*)' valn  =',valn
  !$$$        write(6,*)' ntot  =',ntot

  !      ifcphi  = iopen('CPHI',0,0,mrecb)

  call init_readeigen2() !mrecb,nlmto,mrecg) !initialize m_readeigen
  write(*,*)'nnn3 nband =',nband,latomic

  !! We get nwf after one of following init_... are called
  if (latomic) then
     call init_readeigen_phi_noeval() !nwf,nband,mrecb,mrecg)
  else
     !     if (l1d) then
     !     call init_readeigen_mlw_noeval1D(nwf,nband,mrecb,mrecg)
     !     else
     call init_readeigen_mlw_noeval() !nwf,nband,mrecb,mrecg)
     !     endif
  endif

  write(*,*)'Caution! evals are zero hereafter.'
  write(*,*)'nwf =',nwf
  !      write(*,*)'nwf =',nwf
  !      write(*,*)'mrecb =',mrecb
  !      write(*,*)'mrecg =',mrecg
  write(*,*)'init_readeigen_mlw: done'

  ! QPNT data
  if(nz==0) then
     !        if(readgwinput()) then
     call getkeyvalue("GWinput","<QPNT>",unit=ifqpnt,status=ret)
     !        else
     !         ifqpnt = iopen('QPNT',1,0,0)
     !        endif
  else
     ifqpnt  = iopen('QPNT'//xt(nz),1,0,0)
  endif
  write(6,*)' ifqpnt ret=',ifqpnt,ret

  ! read q-points and states
  ! ---
  ! read QPNT
  lqall      = .true.
  laf        = .false.
  call readx   (ifqpnt,10)
  read (ifqpnt,*) iqall,iaf
  if (iaf   == 1)   laf = .TRUE. 
  call readx   (ifqpnt,100)
  read (ifqpnt,*) ntq
  allocate (itq(ntq))
  read (ifqpnt,*) (itq(i),i=1,ntq)
  deallocate(itq)
  nq         = nqibz
  allocate(q(3,nq))
  call dcopy   (3*nqibz,qibz,1,q,1)

  nspinmx = nspin
  if (laf) nspinmx =1

  close(ifqpnt)

  ! m, 070521
  call getkeyvalue("GWinput","wmat_all",lfull,default=.false.)
  if (lfull) then
     call getkeyvalue("GWinput","wmat_rcut1",rcut1, default=0.01d0 )
     call getkeyvalue("GWinput","wmat_rcut2",rcut2, default=0.01d0 )
     ! m, 070814
     call getkeyvalue("GWinput","wmat_WSsuper",lwssc,default=.true.)
     if (lwssc) then
        allocate(irws(n1*n2*n3*8),rws(3,n1*n2*n3*8),drws(n1*n2*n3*8))
        call wigner_seitz(alat,plat,n1,n2,n3,nrws,rws,irws,drws)
        write(*,*)'*** Wigner-Seitz Super cell'
        do i=1,nrws
           write(*,"(i5,4f12.6,i5)")i,rws(1,i),rws(2,i),rws(3,i), &
                drws(i),irws(i)
        enddo
     else
        allocate(irws(n1*n2*n3),rws(3,n1*n2*n3),drws(n1*n2*n3))
        call super_cell(alat,plat,n1,n2,n3,nrws,rws,irws,drws)
        write(*,*)'*** Super cell (Not Wigner-Seitz super cell)'
        do i=1,nrws
           write(*,"(i5,4f12.6,i5)")i,rws(1,i),rws(2,i),rws(3,i), &
                drws(i),irws(i)
        enddo
     endif
     nrws1 = 1
     do i=1,nrws
        if (drws(i) <= rcut1) nrws1=i
     enddo
     nrws2 = 1
     do i=1,nrws
        if (drws(i) <= rcut2) nrws2=i
     enddo
     allocate(rws1(3,nrws1),rws2(3,nrws2),irws1(nrws1),irws2(nrws2))
     rws1(:,1:nrws1) = rws(:,1:nrws1)
     rws2(:,1:nrws2) = rws(:,1:nrws2)
     irws1(1:nrws1) = irws(1:nrws1)
     irws2(1:nrws2) = irws(1:nrws2)
     nrws = nrws1*nrws2*nrws2
     deallocate(irws,rws,drws)
  else
     rcut1 = 0.0d0
     rcut2 = 0.0d0
     call getkeyvalue("GWinput","wmat_rsite", rsite,3, &
          default=(/0.0d0,0.0d0,0.0d0/),status=ret)
     nrws1 = 1
     nrws2 = 1
     allocate(rws1(3,nrws1),rws2(3,nrws2),irws1(nrws1),irws2(nrws2))
     rws1(:,1) = rsite(:)
     rws2 = 0d0
     irws1(1) = 1
     irws2(1) = 1
     nrws = nrws1*nrws2*nrws2
  endif
  write(*,'(a14,i5,f12.6)')'nrws1, rcut1 =',nrws1,rcut1
  write(*,'(a14,i5,f12.6)')'nrws2, rcut2 =',nrws2,rcut2
  write(*,'(a7,i7)')'nrws  =',nrws

  !$$$c ---  q near zero
  !$$$      write(6,*) 'reading QOP'
  !$$$      if101=ifile_handle()
  !$$$      open (if101,file='Q0P')
  !$$$      read (if101,"(i5)") nq0i
  !$$$      if(.not.exchange) call checkeq(nqibz+nq0i-1, nqnum)
  !$$$      write(6,*) ' *** nqibz nq0i_total=', nqibz,nq0i
  !$$$      nq0it = nq0i
  !$$$      allocate( wqt(1:nq0i),q0i(1:3,1:nq0i) )
  !$$$c      read (if101,"(d24.16,3x, 3d24.16)" )( wqt(i),q0i(1:3,i),i=1,nq0i)
  !$$$      nq0ix = nq0i
  !$$$      do i=1,nq0i
  !$$$      read (if101,* ) wqt(i),q0i(1:3,i)
  !$$$      if(wqt(i)==0d0 ) nq0ix = i-1
  !$$$      enddo
  !$$$      nq0i = nq0ix ! New nq0i July 2001
  write(6,*) ' Used k number in Q0P =', nq0i
  write(6,"(i3,f14.6,2x, 3f14.6)" )(i, wqt(i),q0i(1:3,i),i=1,nq0i)
  close(if101)
  !      allocate( wgt0(nq0i,ngrp) )
  ! Sergey's 1stFeb2005
  !      call q0iwgt2(symgg,ngrp,wqt,q0i,nq0i,
  !     o            wgt0)
  call getkeyvalue("GWinput","allq0i",allq0i,default=.false.)!S.F.Jan06
  call q0iwgt3(allq0i,symgg,ngrp,wqt,q0i,nq0i,               ! & S.F.Jan06
  wgt0)                   ! added allq0i argument
  !--------------------------
  if(nq0i/=0) write(6,*) ' *** tot num of q near 0   =', 1/wgt0(1,1)
  write(6,"('  sum(wgt0) from Q0P=',d14.6)")sum(wgt0)
  !$$$      if(bzcase()==2) then
  !$$$        wgt0= wgt0*wgtq0p()/dble(nqbz)
  !$$$        write(6,"('bzcase=2:  sum(wgt0_modified )=',d14.6)")sum(wgt0)
  !$$$      endif

  ! --- qbze(3,nqibze)
  nqbze  = nqbz *(1 + nq0it)
  allocate( qbze(3, nqbze) )
  call dcopy(3*nqbz, qbz, 1, qbze,1)
  do i = 1,nq0it
     ini = nqbz*(1 + i -1)
     do ix=1,nqbz
        qbze (:,ini+ix)   = q0i(:,i) + qbze(:,ix)
     enddo
  enddo

  ! --- read LDA eigenvalues
  !     ntp0=ntq
  !      allocate(eqx(ntq,nq,nspin),eqx0(ntq,nq,nspin),eqt(nband))
  !      do      is = 1,nspin
  !      do      ip = 1,nq
  !c        iq       = idxk (q(1,ip),qbze,nqbze)
  !c        call rwdd1   (ifev(is), iq, nband, eqt) !direct access read b,hb and e(q,t)
  !        call readeval(q(1,ip),is,eqt)
  !c        write(6,*)' eqt=',eqt
  !        eqx0(1:ntq,ip,is) = eqt(itq(1:ntq))
  !        eqx (1:ntq,ip,is) = rydberg()*(eqt(itq(1:ntq))- ef)
  !      enddo
  !      enddo
  !      deallocate(eqt)


  ! --- info
  call winfo(6,nspin,nq,ntq,is,nbloch &
       ,ngpn1,ngcn1,nqbz,nqibz,ef,deltaw,alat,esmr)

  ! pointer to optimal product basis
  allocate(imdim(natom)) !12may2015 bugfix
  !      call indxmdm (nblocha,nclass,
  !     i              iclass,natom,
  !     o              imdim )
  do ia = 1,natom
     imdim(ia)  = sum(nblocha(iclass(1:ia-1)))+1
  enddo

  if(niw/=0) then ! generate gaussian frequencies x between (0,1) and w=(1-x)/x
     allocate(freqx(niw),freqw(niw),wwx(niw))!,expa(niw))
     call freq01x  (niw, ! & ua,
     freqx,freqw,wwx) !, expa)
  endif
  allocate(expa(1));expa=1d99
  ! c ------ write energy mesh ----------
  ! c      if(.not.sergeys) then
  !c      ifemesh = iopen('emesh.hwmat'//xt(nz),1,-1,0)
  !c      deltax0 = 0d0
  !c      call writeemesh(ifemesh,freqw,niw,freq_r,nwp,deltax0)
  ! c
  iii=ivsumxxx(irk,nqibz*ngrp)
  write(6,*) " sum of nonzero iirk=",iii, nqbz

  !... Read pomatr
  if(smbasis()) then
     write(6,*)' smooth mixed basis : augmented zmel'
     call getngbpomat(nqibz+nq0i, &
          nnmx,nomx)
     nkpo = nqibz+nq0i
     ifpomat = iopen('POmat',0,-1,0) !oct2005
     allocate( pomatr(nnmx,nomx,nkpo),qrr(3,nkpo),nor(nkpo),nnr(nkpo) )
     do ikpo=1,nkpo
        read(ifpomat) qrr(:,ikpo),nn_,no,iqx !readin reduction matrix pomat
        !         write(6,"('smbasis: ikp q no nn=',i5,3f8.4,4i5)") ikp,qrr(:,ikpo),no,nn_
        nnr(ikpo)=nn_
        nor(ikpo)=no
        read(ifpomat) pomatr(1:nn_,1:no,ikpo)
     enddo
     isx = iclose("POmat")
     write(6,*)"Read end of POmat ---"
  else !dummy
     nkpo = 1
     nnmx =1
     nomx =1
     allocate( pomatr(nnmx,nomx,nkpo), qrr(3,nkpo),nor(nkpo),nnr(nkpo) )
  endif

  !-----------------------------------------------------------
  ! calculate the correlated part of the self-energy SEc(qt,w)
  !-----------------------------------------------------------
  ! arrays for sxcf.f
  nlnx4    = nlnx**4
  niwx     = max0 (nw,niw)
  allocate( ppb(nlnmx*nlnmx*mdimx*nclass),  eq(nwf), &
       kount(nqibz), &
       rw_w(nwf,nwf,nwf,nwf,nrws,0:nrw), &
       cw_w(nwf,nwf,nwf,nwf,nrws,0:nrw), &
       rw_iw(nwf,nwf,nwf,nwf,nrws,niw), &
       cw_iw(nwf,nwf,nwf,nwf,nrws,niw))



  ! loop over spin ----------------------------------------------------
  do 2000 is = 1,nspinmx

     ! initialise secq and kount
     kount = 0
     rw_w = 0d0
     cw_w = 0d0
     rw_iw = 0d0
     cw_iw = 0d0


     ! loop over rotations -------------------------------
     ! m
     call chkrot() !ngrp)
     do 1000 irot = 1,ngrp
        !      do 1000 irot = 1,1
        if( sum(abs( irk(:,irot) )) ==0 .AND. &
             sum(abs( wgt0(:,irot))) == 0d0 ) then
           cycle
        endif
        write (6,"(i3,'  out of ',i3,'  rotations ',$)") irot,ngrp
        call cputid (0)
        !        if (irot == 1 .or. irot == ngrp) then
        !        call cputid(0); write(*,*)' ppba '
        !        endif

        ! rotate atomic positions invrot*R = R' + T
        invr       = invrot (irot,invg,ngrp)

        ! -- ppb= <Phi(SLn,r) Phi(SL'n',r) B(S,i,Rr)>
        call ppbafp_v2 (irot,ngrp,is,nspin, &
             il,in,im, nlnm, &
             nl,nn,nclass,nlnmx, &
             mdimx,lx,nx,nxx,  ! & Bloch wave
        cgr, nl-1,        ! & rotated CG
        ppbrd,            ! & radial integrals
        ppb)

        !c -- Rotated gvecc
        !        call rotgvec(symgg(:,:,irot), nqibz,
        !     i    ngcmx,ngcni,qbas,ngveccB,
        !     o    ngveccBr)

        !------------------------------------------------------
        ! calculate the correlated part of the self-energy within GW
        !        ntqx = 0
        !        if(tetra.and.(.not.exchange)) then
        !        ntqx =3*ntq
        !        endif
        ! loop over q
        !        do ip = 1,nq   !;write (*,*) ip,'  out of ',nq,'  k-points ' ! call cputid  (0)
        !          iq  = idxk (q(1,ip),qbz,nqbz)
        !          call rwdd1 (ifev(is),iq, nwf,eq)
        !          call readeval(q(1,ip),is,eq)

        ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! nctot=0 in this version
        nctot0=0
        ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        !         write(*,*) 'wmatq in'

        call wmatqK (kount, irot,ef,ef2,esmr,esmr2, &
                                ! m, 070501
                                !     i              tiat(1:3,1:natom,invr),miat( 1:natom,invr), rsite,
             tiat(1:3,1:natom,invr),miat( 1:natom,invr), &
                                !     i              rws,irws,nrws,
             rws1,rws2,nrws1,nrws2,nrws, &
                                ! 2
             nspin,is, ! & ifcphi,ifrb(is),ifcb(is),ifrhb(is),ifchb(is),
        ifrcw,ifrcwi, &
             qbas,ginv,qibz,qbz,wbz,nstbz, wibz, ! & iindxk,
        nstar,irk,     ! & kount,

        !     i        iiclass,nblocha,i_mnlv,i_mnlc,iicore,incore,iimdim,
        iclass,nblocha,nlnmv, nlnmc,   ! & w(i_mnlv),w(i_mnlc)
        icore,ncore, imdim, &
             ppb, ! &  pdb,dpb,ddb,
        freq_r,freqx, wwx, expa, &
             ua,dwdummy,  ! & deltaw,
        ecore(:,is), &
             
             nlmto,nqibz,nqbz,nctot0, &
                                !     i        index_qbz, n_index_qbz,
             nl,nnc,nclass,natom, &
             nlnmx,mdimx,nbloch,ngrp,nw_i,nw,nrw,niw,niwx,nq, &
             
                                !     i     nblochpmx, ngpn,ngcni,ngpmx,ngcmx,
             nblochpmx,ngpmx,ngcmx, &
                                !     i     geigB(1,1,1,is), ngvecpB,ngveccBr,
                                !     i     ngveccBr,
             wgt0,wqt,nq0i,q0i, symgg(:,:,irot),alat, &
             matmul(symgg(:,:,irot),shtvg(:,invr)),nwf, &
             ifvcfpout, &
                                !     i     shtw,
             exchange, ! & tote, screen, cohtest, ifexsp(is),
        !     i        omega, iwini,iwend,
        nbmx(2),ebmx(2), ! & takao 18June2003
        pomatr, qrr,nnr,nor,nnmx,nomx,nkpo,          ! & oct2005 for pomat
        nwf, &
             rw_w,cw_w,rw_iw,cw_iw) ! acuumulation variable
        !         write(*,*) 'wmatq out'
        ! cccccccccccccccccccccccccccccccccccccccc
        !        iii = ivsum(kount,nqibz*nq)
        !        write(6,*)" sumkount 2=",nqibz,nq,iii
        !        stop "--- kcount test end --- "
        ! cccccccccccccccccccccccccccccccccccccccc
        !< end of q-loop
        !        enddo

        !< end of rotation-loop
        !        enddo
        continue !end of rotation-loop
1000 enddo

     !      if (debug) write(*,*) 'do-rot out'

     ! check that all k { FBZ have been included
     !      if (ivsum(kount,nqibz*nq) /= nqbz*nq) then
     !        iii = ivsum(kount,nqibz*nq)
     !        write(6,*)" ivsum=",iii, nqbz*nq
     !        stop 'hsfp0: missing k-pts'
     !      endif

     !---------------------------------

2001 continue
     !---------------------------------

     ! write <p p | W | p p>
     if (exchange) then
        call       wvmat (is,ifwmat(is),nwf, &
             rws1,rws2,irws1,irws2,nrws1,nrws2,nrws, &
             alat,rcut1,rcut2,rw_w(:,:,:,:,:,0),cw_w(:,:,:,:,:,0))
     else
        call       wwmat (is,ifwmat(is),nw_i,nrw+1,nwf, &
             rws1,rws2,irws1,irws2,nrws1,nrws2,nrws, &
             alat,rcut1,rcut2, &
             freq_r(0:nrw), &
             rw_w,cw_w)
     endif
     continue !end of spin-loop
2000 enddo

  !------------
  ! close files
  !------------
  isx = iclose ('wc.d')
  isx = iclose ('wci.d')
  isx = iclose ('hbe.d')
  isx = iclose ('RBU')
  isx = iclose ('CBU')
  isx = iclose ('RHBU')
  isx = iclose ('CHBU')
  isx = iclose ('EVU')
  isx = iclose ('RBD')
  isx = iclose ('CBD')
  isx = iclose ('RHBD')
  isx = iclose ('CHBD')
  isx = iclose ('EVD')
  !-----------------------------------------------------------------------
  call cputid(0)
  call rx0s(' OK! hwmat')
END PROGRAM hwmat

!-----------------------------------------------------------------------
subroutine wwmat (is,ifwmat,nw_i,nw,nwf, &
     rws1,rws2,irws1,irws2,nrws1,nrws2,nrws, &
     alat,rcut1,rcut2, &
     freq, &
     rw_w,cw_w)

  implicit real*8(a-h,o-z)
  implicit integer (i-n)
  integer(4) :: irws1(nrws1),irws2(nrws2)
  real(8) :: rws1(3,nrws1),rws2(3,nrws2)
  real(8) :: rydberg,hartree
  real(8) :: freq(1:nw),freq2(1:nw)      ! iw is shifted by one
  real(8) :: rw_w(nwf,nwf,nwf,nwf,nrws,1:nw), ! &  iw is shifted by one
  cw_w(nwf,nwf,nwf,nwf,nrws,1:nw)  ! iw is shifted by one
  !     &          ,rw_iw(nwf,nwf,nwf,nwf,nrws,niw),
  !     &           cw_iw(nwf,nwf,nwf,nwf,nrws,niw)

  hartree=2d0*rydberg()

  call cv      (hartree,freq,nw,freq2)
  !      call cv      (hartree,freqw,niw,freqw2)
  call cv      (hartree,rw_w,nwf*nwf*nwf*nwf*nrws*nw,rw_w)
  call cv      (hartree,cw_w,nwf*nwf*nwf*nwf*nrws*nw,cw_w)
  !      call cv      (hartree,rw_iw,nwf*nwf*nwf*nwf*nrws*niw,rw_iw)
  !      call cv      (hartree,cw_iw,nwf*nwf*nwf*nwf*nrws*niw,cw_iw)


  write(ifwmat,*)'*** nwf,nw,alat'
  write(ifwmat,*)nwf,nw,alat
  write(ifwmat,*)'*** rcut1,rcut2'
  write(ifwmat,*)rcut1,rcut2
  write(ifwmat,*)'*** nrws1,nrws2,nrws'
  write(ifwmat,*)nrws1,nrws2,nrws
  write(ifwmat,*)'*** w along the real-axis'
  write(ifwmat,*)freq2
  !      write(ifwmat,*)'*** w along the imaginary-axis'
  !      write(ifwmat,*)freqw
  write(ifwmat,*)'*** rws1,irws1'
  write(ifwmat,*)rws1,irws1
  write(ifwmat,*)'*** rws2,irws2'
  write(ifwmat,*)rws2,irws2
  write(ifwmat,*)'*** rw_w'
  write(ifwmat,*)rw_w
  write(ifwmat,*)'*** cw_w'
  write(ifwmat,*)cw_w
  !      write(ifwmat,*)'*** rw_iw'
  !      write(ifwmat,*)rw_iw
  !      write(ifwmat,*)'*** cw_iw'
  !      write(ifwmat,*)cw_iw

  write(*,*)'Screening part: diagonal,real'
  do i = 1,nwf
     write(*,"('  Wannier',2i5,f12.6,' eV ',2f12.6,' eV')") &
          is,i,freq(1),rw_w(i,i,i,i,1,1),cw_w(i,i,i,i,1,1)
  enddo

  return
end subroutine wwmat
!-----------------------------------------------------------------------
subroutine wvmat (is,ifwmat,nwf, &
     rws1,rws2,irws1,irws2,nrws1,nrws2,nrws, &
     alat,rcut1,rcut2,rw_w,cw_w)

  implicit real*8(a-h,o-z)
  implicit integer (i-n)
  integer(4) :: irws1(nrws1),irws2(nrws2)
  real(8) :: rydberg,hartree
  real(8) :: rws1(3,nrws1),rws2(3,nrws2)
  real(8) :: rw_w(nwf,nwf,nwf,nwf,nrws),cw_w(nwf,nwf,nwf,nwf,nrws)

  hartree=2d0*rydberg()

  call cv      (hartree,rw_w,nwf*nwf*nwf*nwf*nrws,rw_w)
  call cv      (hartree,cw_w,nwf*nwf*nwf*nwf*nrws,cw_w)

  write(ifwmat,*)'*** nwf,alat'
  write(ifwmat,*)nwf,alat
  write(ifwmat,*)'*** rcut1,rcut2'
  write(ifwmat,*)rcut1,rcut2
  write(ifwmat,*)'*** nrws1,nrws2,nrws'
  write(ifwmat,*)nrws1,nrws2,nrws
  write(ifwmat,*)'*** rws1,irws1'
  write(ifwmat,*)rws1,irws1
  write(ifwmat,*)'*** rws2,irws2'
  write(ifwmat,*)rws2,irws2
  write(ifwmat,*)'*** vcoul: Re'
  write(ifwmat,*)rw_w
  write(ifwmat,*)'*** vcoul: Im'
  write(ifwmat,*)cw_w

  write(*,*)'Coulomb interaction: diagonal'
  do i = 1,nwf
     write(*,"('  Wannier',2i5,2f12.6,' eV')")is,i, &
          rw_w(i,i,i,i,1),cw_w(i,i,i,i,1)
  enddo

  return
end subroutine wvmat
!-----------------------------------------------------------------------
subroutine chkrot () !ngrp)
  use m_hamindex,only:   Readhamindex, symops, ngrp
  implicit real*8(a-h,o-z)
  implicit integer (i-n)
  parameter (eps = 1d-6)
  real(8) :: symope(3,3,ngrp)
  !$$$      call readhamindex()
  !$$$      ifi = ifile_handle()
  !$$$      open(ifi,file='SYMOPS')
  !$$$      read(ifi,*)ngrp2
  !$$$      if (ngrp .ne. ngrp2) stop 'chkrot: ngrp error'
  !$$$
  !$$$      do ig = 1,ngrp
  !$$$         read(ifi,*)ig2
  !$$$         do i = 1,3
  !$$$            read(ifi,*)symope(i,1:3,ig)
  !$$$         enddo
  !$$$      enddo
  !$$$
  !$$$      close(ifi)
  symope=symops
  do i = 1,3
     symope(i,i,1) = symope(i,i,1) - 1d0
  enddo
  do i = 1,3
     do j = 1,3
        as = dabs(symope(i,j,1))
        if (as > eps) stop 'chkrot: irot=1 error'
     enddo
  enddo

  return
end subroutine chkrot
!-----------------------------------------------------------------------
