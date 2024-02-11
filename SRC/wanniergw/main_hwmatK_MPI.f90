subroutine hwmatK_MPI()
  !!== Calculates the bare/screened interaction W ===
  !!
  !! W(w) = <phi(n1,dR) phi(n2,0) |W(w)| phi(n3,R') phi(n4,R'+dR')>
  !! phi(n,R): maximally localized Wannier orbital
  !!    n : band index
  !!    R : site
  !!
  !! We included a cRPA method by juelich's group.
  !!
  !! This routine requirs an input from standard IO.
  !! In scripts, you can do it like, prompt>echo mode|../exec/hwmat >lwmat,
  !! where mode is 1 or 2, or 10011.
  !!
  !!  mode= 1: bare Coulomb mode, V, <phi phi | V  | phi phi>
  !!  mode= 2: screening mode, Wc,   <phi phi | Wc | phi phi>
  !!  mode=11: <phi phi | V(omega=0)  | phi phi>
  !!  mode=12: <phi phi | Wc(omega=0) | phi phi>
  !!  mode=10011: cRPA mode
  !!
  !!  nwf: total number of wannier within the primitive cell.
  !!  nrws1: index for phi1(n1,dR) phi(n2,0), where range of dR is by dR < rcut1
  !!  nrws2: index for phi(n3,R'), where range of R' is            by R' < rcut2
  !!         This is also for R'+dR'.
  !!  nrws = nrws1*nrws2*nrws2 : total spatial index for W.
  !!  Thus irws=1,nrws where specifies a set (dR,R',dR').
  !!
  ! cccccccccccccc old history ----
  !x  mode= 3: effective U mode,   <phi phi | U | phi phi>
  !x  mode= 4: similar to mode=1, but phi is atomic orbital+gaussian tail
  !x  mode= 5: similar to mode=2, but phi is atomic orbital+gaussian tail
  !x  mode= 6: similar to mode=3, but phi is atomic orbital+gaussian tail

  ! Mar 2008 Takashi Miyake, MPI version (parallelized along the outer q only)
  ! May 2007 Takashi Miyake, full matrix elements
  ! May 2006 Takashi Miyake, updated for new fpgw
  ! Sep 2004 Takashi Miyake, off-site W
  ! Jul 2004 Takashi Miyake. from hsfp0.m.f
  ! May 2002 Takashi Miyake. Total energy calc.
  ! Apr 2002 takao kotani. multiple argumentation wave per l.
  ! This hsfp0 is build from hsec10.f by F.Aryasetiawan.
  !------------------------------------------------------------
  use m_readqg,only: readngmx2,ngcmx,ngpmx,readqg0,readqg
  use m_hamindex,only:   Readhamindex,symgg=>symops,ngrp,invg=>invgx
  use m_read_bzdata,only: Read_bzdata,qibz,irk,ginv,n1,n2,n3,nqbz,nqibz,nstar,nstbz,qbas=>qlat,qbz,wibz,wbz &
       ,nq0i=>nq0ix,wqt=>wt,q0i
  use m_readeigen,only: onoff_write_pkm4crpa,init_readeigen,init_readeigen2, &
       init_readeigen_mlw_noeval,  nwf !,init_readeigen_phi_noeval
  use m_genallcf_v3,niwg=>niw
  use m_keyvalue,only: getkeyvalue
  use m_readhbe,only: Readhbe, nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg
  use m_zmel_old,only: ppbafp_v2
  use m_hamindex0,only: readhamindex0,iclasst
  ! RS: MPI module
  use rsmpi
  use rsmpi_rotkindex
  use m_mksym_util,only:mptauof
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
  integer:: ixc,iopen, ibas,ibasx,nxx,nbloch,ifqpnt,ifwd,ifmloc, &
       nprecx,mrecl,nblochpmx2,nwp,niwt, nqnum,mdimx,nblochpmx, &
       ifrcw,ifrcwi,  noccxv,maxocc2,noccx,ifvcfpout,iqall,iaf,ntq, &
       i,k,nspinmx, nq,is,ip,iq,idxk,ifoutsex,iclose,ig, &
       mxkp,nqibzxx,ntet,nene,iqi, ix,iw, &
       nlnx4,niwx,irot,invr,invrot,ivsum, ifoutsec,ntqx,&
       ifwmat(2) ,ifxc(2),ifsex(2), ifphiv(2),ifphic(2),ifec,ifexsp(2), &
       ifsecomg(2),ifexx,ndble=8
  real(8) :: pi,tpia,vol,voltot,rs,alpha, &
       qfermi,efx,valn,efnew,edummy,efz,qm,xsex,egex, &
       zfac1,zfac2,dscdw1,dscdw2,dscdw,zfac,ef2=1d99,exx,exxq,exxelgas
  logical :: lqall,laf
  integer(4),allocatable :: itq(:)
  real(8),allocatable    :: q(:,:)
  ! takao
  integer(4),allocatable :: ngvecpB(:,:,:),ngvecp(:,:), ngvecc(:,:),iqib(:),&
       kount(:), nx(:,:),nblocha(:),lx(:) !ngveccBr(:,:,:)
  real(8),allocatable:: vxcfp(:,:,:), &
       wgt0(:,:), &
       ppbrd (:,:,:,:,:,:,:),cgr(:,:,:,:),eqt(:), &
       ppbrdx(:,:,:,:,:,:,:),aaa(:,:),& ! & symope(:,:,:)=symgg, ! qibz(:,:),
       ppb(:), eq(:), &! & ,pdb(:),dpb(:),ddb(:)
       eqx(:,:,:),eqx0(:,:,:),ekc(:),coh(:,:) &
       , rw_w(:,:,:,:,:,:),cw_w(:,:,:,:,:,:), &
       rw_iw(:,:,:,:,:,:),cw_iw(:,:,:,:,:,:)
  complex(8),allocatable:: geigB(:,:,:,:)

  logical :: screen, exchange, cohtest, legas, tote, lueff
  logical :: lcrpa
  real(8) ::  rydberg,hartree
  real(8):: qreal(3), ntot,nocctotg2,tripl,xxx(3,3)
  logical ::nocore

  ! space group infermation
  integer(4),allocatable :: invgx(:), miat(:,:)
  real(8),allocatable    :: tiat(:,:,:),shtvg(:,:)

  ! tetra
  real(8),allocatable :: qz(:,:),qbzxx(:),wbzxx(:),wtet(:,:,:,:), &
       eband(:,:,:), ene(:) !,ecore(:,:)
  integer(4),allocatable ::idtetx(:,:),idtet(:,:),ipq(:) &
       ,iene(:,:,:),ibzx(:) ! ,nstar(:)
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


  integer(4)::  ngpn1,verbose,ngcn1,nwxx !bzcase, mrecg,
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

  integer(4):: nw_i,ifile_handle,if102,if3111,if101

  logical:: latomic,lfull,lstatic,lwssc
  logical:: l1d,lll
  real(8):: rsite(3),rcut1,rcut2
  real(8),allocatable :: rws(:,:),drws(:),rws1(:,:),rws2(:,:)
  integer(4):: nrws,nrws1,nrws2,ir1,ir2,ir3,ir,nrw
  integer(4),allocatable:: irws(:),irws1(:),irws2(:)

  ! RS: variables for MPI
  integer(4) :: input3(3),irot_local,ip_local,iq_local, nq0ixxx
  integer,allocatable :: nq_local(:),iqx_index(:,:)
  real(8),allocatable:: &
       rw_w_sum(:,:,:,:,:,:),rw_iw_sum(:,:,:,:,:,:), &
       cw_w_sum(:,:,:,:,:,:),cw_iw_sum(:,:,:,:,:,:)

  integer :: iwf1, iwf2, iwf3, iwf4, ia
  integer :: ifcou, ifscr
  real(8),allocatable:: &
       rv_w(:,:,:,:,:),cv_w(:,:,:,:,:), &
       rv_iw(:,:,:,:,:),cv_iw(:,:,:,:,:)
  real(8)::ef
  ! For hmagnon (only omega=0 is used)
  integer::nw,nctot0,niw
  logical:: lomega0
  call RSMPI_Init()

  hartree=2d0*rydberg()

  iii=verbose()
  if (Is_IO_Root_RSMPI())write(6,*)' verbose=',iii

  ! mode switch. --------------
  if (Is_IO_Root_RSMPI()) then
     write(6,*) ' --- Choose omodes below ----------------'
     write(6,*) '  V (1) or W (2) or U(3)'
     write(6,*) '  V_omega=0 (11) or W_omega=0 (12)'
     write(6,*) '  [option --- (+ QPNT.{number} ?)] '
     write(6,*) ' --- Put number above ! -----------------'
     call readin5(ixc,nz,idummy)

     input3(1) = ixc
     input3(2) = nz
     input3(3) = idummy
     write(6,*) ' ixc nz=',ixc, nz
     if(ixc==0) stop ' --- ixc=0 --- Choose computational mode!'
  endif

  call MPI_Bcast(input3,3,MPI_INTEGER,io_root_rsmpi, &
       MPI_COMM_WORLD,ierror_rsmpi)
  call RSMPI_Check("MPI_Bcast(input3)",ierror_rsmpi)
  ixc=input3(1)
  nz=input3(2)
  idummy=input3(3)

  lomega0=.false.
  if (ixc==11) then
     ixc=1
     lomega0=.true.
  elseif(ixc==12) then
     ixc=2
     lomega0=.true.
  endif

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
  !      write(6,"(a,9f9.4)")'hwmatK_MPI:ginv=',ginv

  if (Is_IO_Root_RSMPI()) then
     write(6,*)' nqbz  =',nqbz
     !      write(6,*)  qbz
     write(6,*)' nqibz ngrp=',nqibz,ngrp
     !      write(6,*)' irk=',irk
     !      write(6,*)' #### idtetf: ####'
     !      write(6,*) idtetf
  endif

  ! set up work array
  !      call wkinit (iwksize)
  call pshpr(60)

  !--- readin GWIN and LMTO, then allocate and set datas.
  !      nwin =-999    !not readin NW file
  !      efin =-999d0  !not readin EFERMI
  incwfin= -1  !use 7th colmn for core at the end section of GWIN
  call genallcf_v3(incwfin) !in module m_genallcf_v3
  niw=niwg
  ef=1d99
  !      if(ngrp/= ngrp2)
  !     &  call RSMPI_Stop( 'ngrp inconsistent: BZDATA and LMTO GWIN_V2')
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
  !      call getnemx8(nbmx,ebmx)
  !      if (Is_IO_Root_RSMPI())
  !     &  write(6,"('  nbmx ebmx from GWinput=',2i8,2d13.5)") nbmx,ebmx

  !-------------------------------------------------------------------
  !      if (nclass > mxclass) stop ' hsfp0: increase mxclass'
!!!!! WE ASSUME iclass(iatom)= iatom !!!!!!!!!!!!!!!!!!!!!!!!!
  if (nclass /= natom ) call RSMPI_Stop( ' hsfp0: nclass /= natom ') ! We assume nclass = natom.
  if (Is_IO_Root_RSMPI()) write(6,*)' hsfp0: end of genallcf2'

  call pshpr(30)
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
     if (Is_IO_Root_RSMPI()) &
          write(6,*)' --- bare Coulomb mode --- '
     exchange =.true.
     lueff = .false.
     if (Is_IO_Root_RSMPI()) then
        ifwmat(1) = iopen('VMATU',1,-1,0)
        if (nspin == 2) ifwmat(2) = iopen('VMATD',1,-1,0)
     endif
  elseif (ixc ==2) then
     if (Is_IO_Root_RSMPI()) &
          write(6,*)' --- screening (Wc) mode --- '
     exchange =.false.
     lueff = .false.
     lcrpa = .false.
     print *, "lcrpa =", lcrpa
     print *, "lueff=", lueff
     if (Is_IO_Root_RSMPI()) then
        ifwmat(1) = iopen('WcMATU',1,-1,0)
        if (nspin == 2) ifwmat(2) = iopen('WcMATD',1,-1,0)
     endif
  elseif (ixc == 100) then
     if (Is_IO_Root_RSMPI()) &
          write(6,*)' --- screening_crpa (Wc) mode --- '
     exchange =.false.
     lueff = .false.
     lueff = .false.
     lcrpa = .true.
     ixc = 2
     if (Is_IO_Root_RSMPI()) then
        ifwmat(1) = iopen('WcMATU',1,-1,0)
        if (nspin == 2) ifwmat(2) = iopen('WcMATD',1,-1,0)
     endif
     !$$$      elseif (ixc ==3) then
     !$$$        if (Is_IO_Root_RSMPI())
     !$$$     &   write(6,*)' --- Ueff mode --- '
     !$$$        exchange =.false.
     !$$$        lueff = .true.
     !$$$        if (Is_IO_Root_RSMPI()) then
     !$$$          ifwmat(1) = iopen('UMATU',1,-1,0)
     !$$$          if (nspin == 2) ifwmat(2) = iopen('UMATD',1,-1,0)
     !$$$        endif
     !$$$      elseif (ixc==4) then
     !$$$        if (Is_IO_Root_RSMPI())
     !$$$     &   write(6,*)' --- bare Coulomb mode --- '
     !$$$        exchange =.true.
     !$$$        lueff = .false.
     !$$$        latomic = .true.
     !$$$        ixc=1
     !$$$        if (Is_IO_Root_RSMPI()) then
     !$$$          ifwmat(1) = iopen('VMATaU',1,-1,0)
     !$$$          if (nspin == 2) ifwmat(2) = iopen('VMATaD',1,-1,0)
     !$$$        endif
     !$$$      elseif (ixc ==5) then
     !$$$        if (Is_IO_Root_RSMPI())
     !$$$     &   write(6,*)' --- screening (Wc) mode --- '
     !$$$        exchange =.false.
     !$$$        lueff = .false.
     !$$$        latomic = .true.
     !$$$        ixc=2
     !$$$        if (Is_IO_Root_RSMPI()) then
     !$$$          ifwmat(1) = iopen('WcMATaU',1,-1,0)
     !$$$          if (nspin == 2) ifwmat(2) = iopen('WcMATaD',1,-1,0)
     !$$$        endif
     !$$$      elseif (ixc ==6) then
     !$$$        if (Is_IO_Root_RSMPI())
     !$$$     &   write(6,*)' --- Ueff mode --- '
     !$$$        exchange =.false.
     !$$$        lueff = .true.
     !$$$        latomic = .true.
     !$$$        ixc=3
     !$$$        if (Is_IO_Root_RSMPI()) then
     !$$$          ifwmat(1) = iopen('UMATaU',1,-1,0)
     !$$$          if (nspin == 2) ifwmat(2) = iopen('UMATaD',1,-1,0)
     !$$$        endif
  elseif(ixc==10011) then
     write(6,*) 'ixc=10011 pkm4crpa mode'
     exchange =.false.
     lueff = .false.
  else
     call RSMPI_Stop( 'ixc error')
  endif

  !---
  if (Is_IO_Root_RSMPI()) then
     write(6, *) ' --- computational conditions --- '
     write(6,'("    deltaw  =",f13.6)') deltaw
     write(6,'("    ua      =",f13.6)') ua
     write(6,'("    esmr    =",f13.6)') esmr
     write(6,'("    alat voltot =",2f13.6)') alat, voltot
     !      write(6,'("    niw nw dw   =",2i5,f13.6)') niw,nw,dw
  endif

  !>> read dimensions of wc,b,hb
  !      ifhbed     = iopen('hbe.d',1,0,0)
  !      read (ifhbed,*) nprecb,mrecb,mrece,nlmtot,nqbzt, nband,mrecg
  !      if (nprecb == 4)
  !     & call RSMPI_Stop( 'hsfp0: b,hb in single precision')
  call Readhbe()
  call Readhamindex()
  call init_readeigen()!nband,mrece) !initialization of readEigen

  ! --- get space group information ---------------------------------
  ! true class information in order to determine the space group -----------
  !     because the class in the generated GW file is dummy.(iclass(ibas)=ibas should be kept).
  !      if102=ifile_handle()
  !      open (if102,file='CLASS')
  allocate(invgx(ngrp) &
       ,miat(natom,ngrp),tiat(3,natom,ngrp),shtvg(3,ngrp))
  !      if (Is_IO_Root_RSMPI())
  !     & write(6,*)'  --- Readingin CLASS info ---'
  !      do ibas = 1,natom
  !        read(if102,*) ibasx, iclasst(ibas)
  !        if (Is_IO_Root_RSMPI())
  !     &   write(6, "(2i10)") ibasx, iclasst(ibas)
  !      enddo
  !      close(if102)
  call readhamindex0()

  ! Get space-group transformation information. See header of mptaouof.
  call mptauof(symgg,ngrp,plat,natom,pos,iclasst &
       ,miat,tiat,invgx,shtvg )
  !        write (*,*)  'tiat=', tiat(1:3,1:natom,invr),invr

  ! Get array size to call rdpp
  !      call getsrdpp( nclass,nl,
  !     o               ngpmx,ngcmx,nxx )
  call getsrdpp2( nclass,nl,nxx)
  call readngmx2()
  !      call readngmx('QGpsi',ngpmx)
  !      call readngmx('QGcou',ngcmx)
  if (Is_IO_Root_RSMPI()) &
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
  call rdpp_v3(nxx, nl, ngrp, nn, nclass, nspin, symgg, nblocha, lx, nx, ppbrd, mdimx, nbloch, cgr)
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

  if (Is_IO_Root_RSMPI()) write(6,*) ' end of read QGcou'

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
  call pshpr(60)

  if(ixc==10011) goto 1018

  !--- Readin WV.d
  if ( .NOT. exchange) then
     if (lueff) then
        call RSMPI_Stop('Ueff mode not implimented in the MPI version')
        ifwd      = iopen('WV.d.maxloc',1,-1,0)
     else
        ifwd      = iopen('WV.d',1,-1,0)
     endif
     read (ifwd,*) nprecx,mrecl,nblochpmx,nwp,niwt,nqnum,nw_i
     ! nblochpmx from WV.d oct2005
     if (Is_IO_Root_RSMPI()) write(6,"(' Readin WV.d =', 10i5)") &
          nprecx, mrecl, nblochpmx, nwp, niwt, nqnum, nw_i
     !call checkeq(nprecx,ndble)
     if(nprecx/=ndble)call rx("hwmatK_MPI: WVR and WVI not compatible")!call checkeq(nprecx,ndble)
     !       call checkeq(nblochpmx,nblochpmx2)
     !       if (nwt /= nw)   stop 'hwmatK: wrong nw'
     !       nw = nwt
     nw=nwp-1
     if (niwt /= niw) call RSMPI_Stop( 'hwmatK: wrong niw')
     ! m 050518
     niw = 0
     niwt = 0
     if (lueff) then
        ifrcw     = iopen('WVR.maxloc',0,-1,mrecl)
        !          ifrcwi = iopen('WVI.maxloc',0,-1,mrecl)
     else
        ifrcw     = iopen('WVR',0,-1,mrecl)
        !          ifrcwi = iopen('WVI',0,-1,mrecl)
     endif
     !... reading general energy mesh from file 'freq_r'
     open(newunit=if3111,file='freq_r') !this is in a.u.
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
  if (Is_IO_Root_RSMPI()) &
       call getkeyvalue("GWinput","wmat_static",lstatic,default= .FALSE. )
  call MPI_Bcast(lstatic,1,MPI_LOGICAL,io_root_rsmpi, &
       MPI_COMM_WORLD,ierror_rsmpi)
  call RSMPI_Check("MPI_Bcast(lstatic)",ierror_rsmpi)
  if (lstatic) nrw = 0

  !... Readin eigen functions
  !      ifev(1)   = iopen('EVU', 0,0,mrece)
  !     if (nspin==2) ifev(2) = iopen('EVD', 0,0,mrece)

  !$$$c     --- determine Fermi energy ef for given valn (legas case) or corresponding charge given by z and konf.
  !$$$!     When esmr is negative, esmr is geven automatically by efsimplef.
  !$$$      call efsimplef2a_RSMPI(nspin,wibz,qibz,ginv,
  !$$$     i     nband,nqibz
  !$$$     i     ,konf,z,nl,natom,iclass,nclass
  !$$$     i     ,valn, legas, esmr,  !!! valn is input for legas=T, output otherwise.
  !$$$c
  !$$$     i     qbz,nqbz             !index_qbz, n_index_qbz,
  !$$$     o     ,efnew)
  !$$$c
  !$$$c     write(6,*)' end of efsimple'
  !$$$c     ef = efnew
  !$$$c     - check total ele number -------
  !$$$      ntot  = nocctotg2(nspin, ef,esmr, qbz,wbz, nband,nqbz) !wbz
  !$$$      write(6,*)' ef    =',ef
  !$$$      write(6,*)' esmr  =',esmr
  !$$$      write(6,*)' valn  =',valn
  !$$$      write(6,*)' ntot  =',ntot

  !      ifcphi  = iopen('CPHI',0,0,mrecb)
1018 continue
  call init_readeigen2()!mrecb,nlmto,mrecg) !initialize m_readeigen
  write(*,*)'nband =',nband
  !! jan2015
  lll=.false.
  if(ixc==10011 .AND. Is_IO_Root_RSMPI()) lll= .TRUE. 
  call onoff_write_pkm4crpa(lll)
  ! his is for writing pkm4crpa in init_readeigen_mlw_noeval.

!  if (latomic) then
!     call init_readeigen_phi_noeval()!nwf,nband,mrecb,mrecg)
!  else
     !         if (l1d) then
     !           call init_readeigen_mlw_noeval1D(nwf,nband,mrecb,mrecg)
     !         else
     call init_readeigen_mlw_noeval()!nwf,nband,mrecb,mrecg)
     !         endif
!  endif
  if (Is_IO_Root_RSMPI()) then
     write(*,*)'Caution! evals are zero hereafter.'
     write(*,*)'nwf =',nwf
     !      write(*,*)'nband =',nband
     !      write(*,*)'mrecb =',mrecb
     !      write(*,*)'mrecg =',mrecg
     write(*,*)'init_readeigen_mlw: done'
  endif

  !! pkm4crpa mode. generated by init_readeigen_mlw_noeval
  if(ixc==10011) then
     call RSMPI_Finalize()
     if (Is_IO_Root_RSMPI()) call rx0s(' OK! hwmatK_MPI ixc=10011')
     stop
  endif

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
  if (Is_IO_Root_RSMPI()) write(6,*)' ifqpnt ret=',ifqpnt,ret

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
  if (Is_IO_Root_RSMPI()) &
       call getkeyvalue("GWinput","wmat_all",lfull,default= .FALSE. )
  call MPI_Bcast(lfull,1,MPI_LOGICAL,io_root_rsmpi, &
       MPI_COMM_WORLD,ierror_rsmpi)
  call RSMPI_Check("MPI_Bcast(lfull)",ierror_rsmpi)

  print *, "Here!!!!!!!!!!!!!", lfull, lwssc, nrws

  if (lfull) then
     if (Is_IO_Root_RSMPI()) then
        call getkeyvalue("GWinput","wmat_rcut1",rcut1, default=0.01d0 )
        call getkeyvalue("GWinput","wmat_rcut2",rcut2, default=0.01d0 )
        ! m, 070814
        call getkeyvalue("GWinput","wmat_WSsuper",lwssc,default=.true.)
     endif ! Is_IO_Root_RSMPI
     call MPI_Bcast(rcut1,1,MPI_DOUBLE_PRECISION,io_root_rsmpi, &
          MPI_COMM_WORLD,ierror_rsmpi)
     call RSMPI_Check("MPI_Bcast(rcut1)",ierror_rsmpi)
     call MPI_Bcast(rcut2,1,MPI_DOUBLE_PRECISION,io_root_rsmpi, &
          MPI_COMM_WORLD,ierror_rsmpi)
     call RSMPI_Check("MPI_Bcast(rcut2)",ierror_rsmpi)
     call MPI_Bcast(lwssc,1,MPI_LOGICAL,io_root_rsmpi, &
          MPI_COMM_WORLD,ierror_rsmpi)
     call RSMPI_Check("MPI_Bcast(lwssc)",ierror_rsmpi)
     if (lwssc) then
        allocate(irws(n1*n2*n3*8),rws(3,n1*n2*n3*8),drws(n1*n2*n3*8))
        call wigner_seitz(alat,plat,n1,n2,n3,nrws,rws,irws,drws)
        if (Is_IO_Root_RSMPI()) then
           write(*,*)'*** Wigner-Seitz Super cell'
           do i=1,nrws
              write(*,"(i5,4f12.6,i5)")i,rws(1,i),rws(2,i),rws(3,i), &
                   drws(i),irws(i)
           enddo
        endif
     else
        allocate(irws(n1*n2*n3),rws(3,n1*n2*n3),drws(n1*n2*n3))
        call super_cell(alat,plat,n1,n2,n3,nrws,rws,irws,drws)
        if (Is_IO_Root_RSMPI()) then
           write(*,*)'*** Super cell (Not Wigner-Seitz super cell)'
           do i=1,nrws
              write(*,"(i5,4f12.6,i5)")i,rws(1,i),rws(2,i),rws(3,i), &
                   drws(i),irws(i)
           enddo
        endif
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
     if (Is_IO_Root_RSMPI()) &
          call getkeyvalue("GWinput","wmat_rsite", rsite,3, &
          default=(/0.0d0,0.0d0,0.0d0/),status=ret)
     call MPI_Bcast(rsite,3,MPI_DOUBLE_PRECISION,io_root_rsmpi, &
          MPI_COMM_WORLD,ierror_rsmpi)
     call RSMPI_Check("MPI_Bcast(rsite)",ierror_rsmpi)
     rcut1 = 0.0d0
     rcut2 = 0.0d0
     nrws1 = 1
     nrws2 = 1
     allocate(rws1(3,nrws1),rws2(3,nrws2),irws1(nrws1),irws2(nrws2))
     rws1(:,1) = rsite(:)
     rws2 = 0d0
     irws1(1) = 1
     irws2(1) = 1
     nrws = nrws1*nrws2*nrws2
  endif
  print *, "Here!!!!!!!!!!!!!", lfull, lwssc, nrws
  if (Is_IO_Root_RSMPI()) then
     write(*,'(a14,i5,f12.6)')'nrws1, rcut1 =',nrws1,rcut1
     write(*,'(a14,i5,f12.6)')'nrws2, rcut2 =',nrws2,rcut2
     write(*,'(a7,i7)')'nrws  =',nrws
  endif

  !$$$c ---  q near zero
  !$$$      write(6,*) 'reading QOP'
  !$$$      if101=ifile_handle()
  !$$$      open (if101,file='Q0P')
  !$$$      read (if101,"(i5)") nq0i
  !$$$      if(.not.exchange) call checkeq(nqibz+nq0i-1, nqnum)
  !$$$      if (Is_IO_Root_RSMPI())
  !$$$     & write(6,*) ' *** nqibz nq0i_total=', nqibz,nq0i
  !$$$      nq0it = nq0i
  !$$$      allocate( wqt(1:nq0i),q0i(1:3,1:nq0i) )
  !$$$c      read (101,"(d24.16,3x, 3d24.16)" )( wqt(i),q0i(1:3,i),i=1,nq0i)
  !$$$      nq0ix = nq0i
  !$$$      do i=1,nq0i
  !$$$      read (if101,* ) wqt(i),q0i(1:3,i)
  !$$$      if(wqt(i)==0d0 ) nq0ix = i-1
  !$$$      enddo
  !$$$      nq0i = nq0ix ! New nq0i July 2001
  !$$$      if (Is_IO_Root_RSMPI()) then
  write(6,*) ' Used k number in Q0P =', nq0i
  write(6,"(i3,f14.6,2x, 3f14.6)" )(i, wqt(i),q0i(1:3,i),i=1,nq0i)
  !$$$      endif
  !$$$      close(if101)
  allocate( wgt0(nq0i,ngrp) )
  ! Sergey's 1stFeb2005
  !      call q0iwgt2(symgg,ngrp,wqt,q0i,nq0i,
  !     o            wgt0)
  if (Is_IO_Root_RSMPI()) &
       call getkeyvalue("GWinput","allq0i",allq0i,default= .FALSE. )!S.F.Jan06
  call MPI_Bcast(allq0i,1,MPI_LOGICAL,io_root_rsmpi, &
       MPI_COMM_WORLD,ierror_rsmpi)
  call RSMPI_Check("MPI_Bcast(allq0i)",ierror_rsmpi)
  call q0iwgt3(allq0i,symgg,ngrp,wqt,q0i,nq0i,     wgt0)                   ! added allq0i argument
  !--------------------------
  if (Is_IO_Root_RSMPI()) then
     if(nq0i/=0) write(6,*) ' *** tot num of q near 0   =', 1/wgt0(1,1)
     write(6,"('  sum(wgt0) from Q0P=',d14.6)")sum(wgt0)
  endif
  !$$$      if(bzcase()==2) then
  !$$$        wgt0= wgt0*wgtq0p()/dble(nqbz)
  !$$$        if (Is_IO_Root_RSMPI())
  !$$$     &   write(6,"('bzcase=2:  sum(wgt0_modified )=',d14.6)")sum(wgt0)
  !$$$      endif

  ! --- qbze(3,nqibze)
  nqbze  = nqbz *(1 + nq0i)
  allocate( qbze(3, nqbze) )
  call dcopy(3*nqbz, qbz, 1, qbze,1)
  do i = 1,nq0i
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
  if (Is_IO_Root_RSMPI()) &
       call winfo(6,nspin,nq,ntq,nspin,nbloch &
       ,ngpn1,ngcn1,nqbz,nqibz,ef,deltaw,alat,esmr)

  ! pointer to optimal product basis
  allocate(imdim(natom)) !bugfix 12may2015
  !      call indxmdm (nblocha,nclass,
  !     i              iclass,natom,
  !     o              imdim )
  do ia = 1,natom
     imdim(ia)  = sum(nblocha(iclass(1:ia-1)))+1
  enddo
  if(niw/=0) then ! generate gaussian frequencies x between (0,1) and w=(1-x)/x
     allocate(freqx(niw),freqw(niw),wwx(niw))!,expa(niw))
     call freq01x  (niw, freqx,freqw,wwx) 
  endif
  ! c ------ write energy mesh ----------
  ! c      if(.not.sergeys) then
  !c      ifemesh = iopen('emesh.hwmat'//xt(nz),1,-1,0)
  !c      deltax0 = 0d0
  !c      call writeemesh(ifemesh,freqw,niw,freq_r,nwp,deltax0)
  ! c
  iii=count(irk/=0) !ivsumxxx(irk,nqibz*ngrp)
  if (Is_IO_Root_RSMPI()) write(6,*) " sum of nonzero iirk=",iii, nqbz
  !... Read pomatr
  if(smbasis()) then
     call RSMPI_Stop('hwmatK_MPI: smbasis notimplemented!')
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
  allocate( ppb(nlnmx*nlnmx*mdimx*nclass),  eq(nband), &
       kount(nqibz), &
       rw_w(nwf,nwf,nwf,nwf,nrws,0:nrw), &
       cw_w(nwf,nwf,nwf,nwf,nrws,0:nrw), &
       rw_iw(nwf,nwf,nwf,nwf,nrws,niw), &
       cw_iw(nwf,nwf,nwf,nwf,nrws,niw), &
       rv_w(nwf,nwf,nwf,nwf,nrws), &
       cv_w(nwf,nwf,nwf,nwf,nrws))

  ! RS: set MPI parameters
  ! RS: see gwsrc/RSMPI_rotkindex_mod.F
  call MPI_Barrier(MPI_COMM_WORLD,ierror_rsmpi)
  call RSMPI_Check("MPI_Barrier",ierror_rsmpi)
  !      call setup_rotkindex(ngrp,irk,wgt0,bzcase(),nqibz,nq0i,nq)
  nq0ixxx=0
  !      call setup_rotkindex(ngrp,irk,wgt0,bzcase(),nqibz,nq0ixxx,1) ! nq=1
  call setup_rotkindex(ngrp,irk,wgt0,1,nqibz,nq0ixxx,1) ! nq=1

  if (Is_IO_Root_RSMPI()) then ! debug
     if(nq0i/=0) then
        write(6,*) 'RS: total number of k-points should be', &
             nqbz +  1/wgt0(1,1)   - 2 + 1 !bzcase()
     else
        write(6,*) 'RS: total number of k-points should be', &
             nqbz  - 2 + 1 !bzcase()
     endif
  endif
  call MPI_Barrier(MPI_COMM_WORLD,ierror_rsmpi)
  call RSMPI_Check("MPI_Barrier",ierror_rsmpi)
  ! RS: openlogfile for each process
  if (ixc == 1) then
     ifile_rsmpi = iopen ('lwt_v.MPI'//myrank_id_rsmpi,1,3,0)
  else if (ixc == 2) then
     ifile_rsmpi = iopen ('lwt_wc.MPI'//myrank_id_rsmpi,1,3,0)
  else if (ixc == 3) then
     ifile_rsmpi = iopen ('lwt_u.MPI'//myrank_id_rsmpi,1,3,0)
  else if (ixc == 4) then
     ifile_rsmpi = iopen ('lwt_v_phi.MPI'//myrank_id_rsmpi,1,3,0)
  else if (ixc == 5) then
     ifile_rsmpi = iopen ('lwt_wc_phi.MPI'//myrank_id_rsmpi,1,3,0)
  else if (ixc == 6) then
     ifile_rsmpi = iopen ('lwt_u_phi.MPI'//myrank_id_rsmpi,1,3,0)
  elseif( ixc==10011) then
  else
     call RSMPI_Stop("unknown ixc")
  endif
  ! RS: print how symmetry operations and k-points are devided ..
  write(ifile_rsmpi,*) "rank : ", myrank_id_rsmpi
  write(ifile_rsmpi,*) "nrotk_local:",nk_local_qkgroup
  write(ifile_rsmpi,*) "nrot_local :",nrot_local_rotk
  if (nrot_local_rotk > 0) then
     write(ifile_rsmpi,*) &
          "irot_index :",irot_index_rotk(1:nrot_local_rotk)
  endif
  write(ifile_rsmpi,*) "nk_local(1:ngrp) :"
  write(ifile_rsmpi,*) nk_local_rotk(:)
  do irot=1,ngrp
     if (nk_local_rotk(irot) > 0) then
        write(ifile_rsmpi,*) "> irot,nk_local(irot) = ", &
             irot, nk_local_rotk(irot)
        write(ifile_rsmpi,*) "   ik_index : ", &
             ik_index_rotk(irot,1:nk_local_rotk(irot))
     endif
  enddo
  ! ccccccccccccc
  call MPI_Barrier(MPI_COMM_WORLD,ierror_rsmpi)
  call RSMPI_Check("MPI_Barrier",ierror_rsmpi)

  if (Is_IO_Root_RSMPI()) then
     write(6,*) "RS: loop over spin --"
  endif


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
     do 1000 irot_local = 1,nrot_local_rotk
        irot = irot_index_rotk(irot_local)
        if( sum(abs( irk(:,irot) )) ==0 .AND. &
             sum(abs( wgt0(:,irot))) == 0d0 ) then
           ! RS: this should not happen... see gwsrc/RSMPI_rotindex_mod.F
           call RSMPI_Stop("hwmatK_RSMPI, cylce occurs in do 1000 -loop!")
           cycle
        endif
        write (6,"(i3,'  out of ',i3,'  rotations ',$)") irot,ngrp
        call cputid (0)

        ! rotate atomic positions invrot*R = R' + T
        !        invr       = invrot (irot,invg,ngrp)
        invr     = invg(irot)
        ! -- ppb= <Phi(SLn,r) Phi(SL'n',r) B(S,i,Rr)>
        call ppbafp_v2 (irot,ngrp,is, &
             mdimx,lx,nx,nxx,   & !Bloch wave
             cgr, nl-1,         & !rotated CG
             ppbrd,             & !radial integrals
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
        !          call rwdd1 (ifev(is),iq, nband,eq)
        !          call readeval(q(1,ip),is,eq)

        ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! nctot=0 in this version
        nctot0=0
        ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        write(*,*) 'wmatq in',irot_local,nrot_local_rotk

        call wmatqk_MPI (kount, irot,ef,ef2,esmr,esmr2, &
                                ! m, 070501
                                !     i              tiat(1:3,1:natom,invr),miat( 1:natom,invr), rsite,
             tiat(1:3,1:natom,invr),miat( 1:natom,invr), &
                                !     i              rws,irws,nrws,
             rws1,rws2,nrws1,nrws2,nrws, &
                                ! 2
             nspin,is,  & !ifcphi,ifrb(is),ifcb(is),ifrhb(is),ifchb(is),
             ifrcw,ifrcwi, &
             qbas,ginv,qibz,qbz,wbz,nstbz, wibz, &! & iindxk,
             nstar,irk,   &  ! & kount,

             !     i        iiclass,nblocha,i_mnlv,i_mnlc,iicore,incore,iimdim,
             iclass,nblocha,nlnmv, nlnmc, &  ! & w(i_mnlv),w(i_mnlc)
             icore,ncore, imdim, &
             ppb, &! &  pdb,dpb,ddb,
             freq_r,freqx, wwx, expa, &
             ua,dwdummy,  &! & deltaw,
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
             matmul(symgg(:,:,irot),shtvg(:,invr)),nband, &
             ifvcfpout, &
                                !     i     shtw,
             exchange, &! & tote, screen, cohtest, ifexsp(is),
        !     i        omega, iwini,iwend,
        !     i     nbmx(2),ebmx(2), !takao 18June2003
             pomatr, qrr,nnr,nor,nnmx,nomx,nkpo, &         ! & oct2005 for pomat
             nwf, &
             rw_w,cw_w,rw_iw,cw_iw) ! acuumulation variable
        write(*,*) 'wmatq out',irot_local,nrot_local_rotk
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


     print *,'xxxxxxxxbbbbbbbbbbbbbbbbbbbb'
     ! RS: accumulate zw
     allocate( rw_w_sum(nwf,nwf,nwf,nwf,nrws,0:nrw), &
          cw_w_sum(nwf,nwf,nwf,nwf,nrws,0:nrw), &
          rw_iw_sum(nwf,nwf,nwf,nwf,nrws,niw), &
          cw_iw_sum(nwf,nwf,nwf,nwf,nrws,niw))
     call MPI_AllReduce(rw_w,rw_w_sum,(nrw+1)*nwf**4*nrws, &
          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror_rsmpi)
     call RSMPI_Check("MPI_AllReduce,rw_w",ierror_rsmpi)
     rw_w = rw_w_sum
     call MPI_AllReduce(cw_w,cw_w_sum,(nrw+1)*nwf**4*nrws, &
          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror_rsmpi)
     call RSMPI_Check("MPI_AllReduce,cw_w",ierror_rsmpi)
     cw_w = cw_w_sum

     if (niw > 0) then
        call MPI_AllReduce(rw_iw,rw_iw_sum,niw*nwf**4*nrws, &
             MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror_rsmpi)
        call RSMPI_Check("MPI_AllReduce,rw_iw",ierror_rsmpi)
        rw_iw = rw_iw_sum
        call MPI_AllReduce(cw_iw,cw_iw_sum,niw*nwf**4*nrws, &
             MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror_rsmpi)
        call RSMPI_Check("MPI_AllReduce,cw_iw",ierror_rsmpi)
        cw_iw = cw_iw_sum
     endif                   ! niw
     deallocate(rw_w_sum,cw_w_sum,rw_iw_sum,cw_iw_sum)

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
     if (Is_IO_Root_RSMPI()) then
        if (exchange) then
           call       wvmat (is,ifwmat(is),nwf, &
                rws1,rws2,irws1,irws2,nrws1,nrws2,nrws, &
                alat,rcut1,rcut2,rw_w(:,:,:,:,:,0),cw_w(:,:,:,:,:,0), lcrpa, lomega0)
           rv_w = rw_w(:,:,:,:,:,0)
           cv_w = cw_w(:,:,:,:,:,0)
           !        write(*,*) "This one!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!"
           !        write(*,*) rw_w
           !        write(*,*) cw_w
           !        write(*,*) rw_w(:,:,:,:,:,0)
           !        write(*,*) cw_w(:,:,:,:,:,0)

        else
           call       wwmat (is,ifwmat(is),nw_i,nrw+1,nwf, &
                rws1,rws2,irws1,irws2,nrws1,nrws2,nrws, &
                alat,rcut1,rcut2, &
                freq_r(0:nrw), &
                rw_w,cw_w,rv_w,cv_w, &
                lcrpa, lomega0)
        endif
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
  call cputid(ifile_rsmpi)
  call cputid(0)
  call RSMPI_Finalize()
  if (Is_IO_Root_RSMPI()) call rx0s(' OK! hwmatK_MPI')
end subroutine hwmatK_MPI

!-----------------------------------------------------------------------
subroutine wwmat (is,ifwmat,nw_i,nw,nwf, &
     rws1,rws2,irws1,irws2,nrws1,nrws2,nrws, &
     alat,rcut1,rcut2, &
     freq, &
     rw_w,cw_w,rv_w,cv_w, &
     lcrpa, lomega0)
  implicit real*8(a-h,o-z)
  implicit integer (i-n)
  integer(4) :: irws1(nrws1),irws2(nrws2)
  real(8) :: rws1(3,nrws1),rws2(3,nrws2)
  real(8) :: rydberg,hartree
  real(8) :: freq(1:nw),freq2(1:nw)   
  real(8) :: rw_w(nwf,nwf,nwf,nwf,nrws,1:nw), cw_w(nwf,nwf,nwf,nwf,nrws,1:nw)  
  integer:: iwf1, iwf2, iwf3, iwf4, ifreq2
  real(8) :: rv_w(nwf,nwf,nwf,nwf,nrws),cv_w(nwf,nwf,nwf,nwf,nrws)
  logical:: lcrpa, lomega0

  hartree=2d0*rydberg()

  !      call cv      (hartree,freq,nw,freq2)
  !      call cv      (hartree,rw_w,nwf*nwf*nwf*nwf*nrws*nw,rw_w)
  !      call cv      (hartree,cw_w,nwf*nwf*nwf*nwf*nrws*nw,cw_w)
  freq2 = hartree*freq
  rw_w  = hartree*rw_w
  cw_w  = hartree*cw_w

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


  !        write(*,*) rw_w
  !        write(*,*) cw_w
  !        write(*,*) rw_w(:,:,:,:,:,0)
  !        write(*,*) cw_w(:,:,:,:,:,0)

  write(*,*)'Writing Screened Couloumb interaction (W-v) : Real'
  ifscr = ifile_handle()

  if ((is==1) .AND. (lcrpa .eqv. .FALSE. )) then
     open(ifscr,file="Screening_W-v.UP")
  else if ((is==2) .AND. (lcrpa .eqv. .FALSE. )) then
     open(ifscr,file="Screening_W-v.DN")
  else if ((is==1) .AND. (lcrpa .eqv. .TRUE. )) then
     open(ifscr,file="Screening_W-v_crpa.UP")
  else if ((is==2) .AND. (lcrpa .eqv. .TRUE. )) then
     open(ifscr,file="Screening_W-v_crpa.DN")
  end if
  if (lomega0) then !only omega=0
     nrws1=1
     nw=1
  endif
  do ir1 = 1,nrws1
     do ifreq2 = 1,nw
        do iwf1 = 1,nwf
           do iwf2 = 1,nwf
              do iwf3 = 1, nwf
                 do iwf4 = 1, nwf
                    write(ifscr,"(' Wannier ',2i5, 3f12.6, 5i5,4f12.6)") &
                         ir1, irws1(ir1), rws1(:,ir1) &
                         ,is,iwf1,iwf2,iwf3,iwf4,freq(ifreq2), freq2(ifreq2) &
                         ,rw_w(iwf1,iwf2,iwf3,iwf4,ir1,ifreq2) &
                         ,cw_w(iwf1,iwf2,iwf3,iwf4,ir1,ifreq2)
                 enddo
              enddo
           enddo
        enddo
        !          write(ifscr,*)''
     enddo
  enddo
  close(ifscr)
  print *, "lcrpa, lomega0 =", lcrpa, lomega0
  print *, "lueff=", lueff
  if (lcrpa .eqv. .FALSE. ) then
     write(*,*) 'See "Screening_W-v.UP" and "Screening_W-v.DN" files'
  else if (lcrpa .eqv. .TRUE. ) then
     write(*,*) 'See "Screening_W-v_crpa.UP" and "Screening_W-v_crpa.DN" files'
  end if

end subroutine wwmat
!-----------------------------------------------------------------------
subroutine wvmat (is,ifwmat,nwf, &
     rws1,rws2,irws1,irws2,nrws1,nrws2,nrws, &
     alat,rcut1,rcut2,rw_w,cw_w, &
     lcrpa, lomega0)

  implicit real*8(a-h,o-z)
  implicit integer (i-n)
  integer(4) :: irws1(nrws1),irws2(nrws2)
  real(8) :: rydberg,hartree
  real(8) :: rws1(3,nrws1),rws2(3,nrws2)
  real(8) :: rw_w(nwf,nwf,nwf,nwf,nrws),cw_w(nwf,nwf,nwf,nwf,nrws)
  integer:: iwf1, iwf2, iwf3, iwf4,ir1
  logical:: lcrpa, lomega0

  hartree=2d0*rydberg()

  !      call cv      (hartree,rw_w,nwf*nwf*nwf*nwf*nrws,rw_w)
  !      call cv      (hartree,cw_w,nwf*nwf*nwf*nwf*nrws,cw_w)
  rw_w= hartree*rw_w
  cw_w= hartree*cw_w

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


  write(*,*)'Coulomb interaction (v) : '
  ifcou = ifile_handle()
  if (is==1) then
     open(ifcou,file="Coulomb_v.UP")
  else if (is==2) then
     open(ifcou,file="Coulomb_v.DN")
  end if
  if (lomega0) then !only omega=0
     nrws1=1
     nw=1
  endif
  do ir1 = 1,nrws1
     do iwf1 = 1,nwf
        do iwf2 = 1,nwf
           do iwf3 = 1, nwf
              do iwf4 = 1, nwf
                 write(ifcou,"(' Wannier ',2i5, 3f12.6, 5i5,2f12.6)") &
                      ir1, irws1(ir1), rws1(:,ir1), is, iwf1, iwf2, iwf3, iwf4 &
                      ,rw_w(iwf1,iwf2,iwf3,iwf4,ir1),cw_w(iwf1,iwf2,iwf3,iwf4,ir1)
              enddo
           enddo
        enddo
     enddo
     write(*,*)
  enddo
  close(ifcou)

  return
end subroutine wvmat
!-----------------------------------------------------------------------
subroutine chkrot ()
  use m_hamindex,only:   Readhamindex, symops, ngrp
  implicit integer (i-n)
  implicit real*8(a-h,o-z)
  parameter (eps = 1d-6)
  real(8) :: symope(3,3,ngrp)

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
