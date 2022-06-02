program hef
  !-------------------------------------------------------------
  ! hef is from hsfp0 but it is only to calcuate ef and Ex test.
  ! so this routine contains unused things so much.
  use m_hamindex,only:   Readhamindex
  use m_readeigen,only:init_readeigen,readeval
  use m_read_bzdata,only: Read_bzdata,ginv,nqibz,nqbz,qibz,wibz,qbz,wbz
  use m_genallcf_v3,only: genallcf_v3, &
       nclass,natom,nspin,nl,nn,nnv,nnc, &
       nlmto,nlnmx, nctot,niw,! & nw_input=>nw,ef
  alat, delta,deltaw,esmr,clabl,iclass, ! & diw,dw,
  il,in,im,nlnm, &
       plat, pos,z,ecore,konf,nlnx
  use m_keyvalue,only: getkeyvalue
  use m_readhbe,only: Readhbe, nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg
  use m_lgunit,only:m_lgunit_init

  implicit none

  !      real(8) :: shtw
  integer(4) :: mxclass,ngnmax &
       ! ,mbytes,mwords,iwksize,
       !     &   natom,nclass,ipos,ngrp,igrp,
       !     &   iinvg,instar,iirk,
       !     o   nspin,nl,nn,nnv,nnc,
       !     o   inindx,inindxv,inindxc,iiclass,
       !     d   nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc,
       !     o   iz,
       !     o   iil,iin,iim,iilnm,inlnm,
       !     o   iilv,iinv,iimv,iilnmv,inlnmv,
       !     o   iilc,iinc,iimc,iilnmc,inlnmc,
       !     o   incwf,iecore,ikonf,iicore,incore,nctot,
       !     o   imagw,niw,nw,ifreq,
       ixc,! & ifhbed, nprecb,mrecb,mrece,nlmtot,nqbzt, nband, !iopen,
  ibas,ibasx,ngpmx,nxx,ngcmx,nbloch,ifqpnt,ifwd, &
       nprecx,nblochpmx2,nwt,niwt, nqnum,mdimx,nblochpmx, &
       ifrcw,ifrcwi,  noccxv,maxocc,noccx,ifvcfpout,iqall,iaf,ntq, &
       i,k,nspinmx, nq,is,ip,iq,idxk,ifoutsex,iclose,nq0i,ig,iimdim, &
       ifreqx,ifreqw,iwx,iexpa,mxkp,nqibzxx,ntet,nene,iqi, ix,iw, &
       nlnx4,niwx,irot,invr,invrot,ivsum, ifoutsec,ntqx, &
       
       !     &   ifrb(2),ifcb(2),ifrhb(2),ifchb(2),ifev(2),ifsec(2)
       ifev(2), &
       !     &             ,ifxc(2),ifsex(2), ifphiv(2),ifphic(2),ifec,
       ndble=8

  !      real(8) :: alat,ef,diw,dw,delta,pi,tpia,vol,voltot,rs,alpha,
  real(8) :: pi,tpia,vol,voltot,rs,alpha, &
       qfermi,efx,valn,efnew,edummy,efz,qm,xsex, &
       zfac1,zfac2,dscdw1,dscdw2,dscdw,zfac

  !     &   plat(3,3)

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
  integer(4),allocatable :: ngvecpB(:,:,:),ngveccB(:,:,:), &
       ngvecp(:,:), ngvecc(:,:),ngpn(:),ngcni(:),iqib(:), &
       kount(:,:), nx(:,:),nblocha(:),lx(:),ngveccBr(:,:,:)
  real(8),allocatable:: vxcfp(:,:,:), &
       wqt(:), wgt0(:,:),q0i(:,:), &
       ppbrd (:,:,:,:,:,:,:),cgr(:,:,:,:),eqt(:), &
       ppbrdx(:,:,:,:,:,:,:),aaa(:,:),symope(:,:,:), &
       ppb(:),pdb(:),dpb(:),ddb(:), eq(:,:), &
       eqx(:,:,:),eqx0(:,:,:),ekc(:),coh(:,:)
  complex(8),allocatable:: geigB(:,:,:,:) ,zsec(:,:,:)

  logical :: exchange, legas=.false.
  real(8):: qreal(3), tripl !ntot,nocctotg2,
  logical ::nocore

  ! space group infermation
  integer(4),allocatable :: iclasst(:), invgx(:), miat(:,:)
  real(8),allocatable    :: tiat(:,:,:),shtvg(:,:)

  ! tetra
  real(8),allocatable :: qz(:,:),qbzxx(:),wbzxx(:),wtet(:,:,:,:), &
       eband(:,:,:), ene(:) !,ecore(:,:)
  integer(4),allocatable ::idtetx(:,:),idtet(:,:),ipq(:) &
       ,iene(:,:,:),ibzx(:)
  !      real(8) :: qbasmc(3,3)

  ! worksize in megabytes (1 word = 4 bytes)
  !      parameter (mbytes=60)
  !      parameter (mwords=mbytes/4)
  !      parameter (iwksize=mwords * 1000 *1000)
  !      integer w
  !      common /w/ w(iwksize)

  integer(4) ::ib,iqx,igp,iii,isx
  integer(4) :: ipex,itpex,ifexspxx,ifexspxxw
  character(12) :: filenameex
  character(3) :: charnum3
  real(8) :: ex,eee_dummy,exwgt_dummy,eee,exwgt,wfac,wfacx
  logical :: exsptest

  integer(4):: incwfin,ret,ifile_handle
  real(8)::ef
  !      logical:: readgwinput

  !      call headver('hef',0)

  !-----------------------------------------------------------------------
  !---  readin BZDATA. See gwsrc/rwbzdata.f
  !--------readin data set when you call read_BZDATA ---------------
  !       integer(4)::ngrp,nqbz,nqibz,nqbzw,nteti,ntetf,
  !     &   n_index_qbz
  !       real(8):: qbas(3,3),ginv(3,3)
  !       real(8),allocatable:: qbz(:,:),wbz(:),qibz(:,:)
  !     &    ,wibz(:),qbzw(:,:)
  !       integer(4),allocatable:: idtetf(:,:),ib1bz(:),idteti(:,:)
  !     &    ,nstar(:),irk(:,:),index_qbz(:,:,:)
  !-----------------------------------------------------------------
  call m_lgunit_init()
  call read_BZDATA()
  write(6,*)' nqbz qbz =',nqbz
  !      write(6,*)' nqibz ngrp=',nqibz,ngrp2
  call pshprt(60)
  !-------------------------------
  ! generate all quantities needed
  !-------------------------------
  incwfin=  0     ! dummy
  call genallcf_v3(incwfin) !in module m_genallcf_v3
  !      if(ngrp/= ngrp2) call rx( 'ngrp inconsistent: BZDATA and LMTO GWIN_V2')
  !-------------------------------------------------------------------
  if (nclass /= natom ) call rx( ' hsfp0: nclass /= natom ')
  write(6,*)' hsfp0: end of genallcf2'
  call pshprt(30)
  pi   = 4d0*datan(1d0)
  tpia = 2d0*pi/alat
  voltot = abs(alat**3*tripl(plat,plat(1,2),plat(1,3)))
  call Readhbe()
  call Readhamindex()
  call init_readeigen()!nband,mrece) !initialization of readEigen

  ! --- determine Fermi energy ef for given valn (legas case) or corresponding charge given by z and konf.
  ! When esmr is negative, esmr is geven automatically by efsimplef.
  call efsimplef2a(nspin,wibz,qibz,ginv, &
       nband,nqibz &
       ,konf,z,nl,natom,iclass,nclass &
       ,valn, legas, esmr,  ! & !! valn is input for legas=T, output otherwise.

  qbz,nqbz ! & index_qbz, n_index_qbz,
  ,efnew)
  ef = efnew
  !- check total ele number -------
  !      ntot  = nocctotg2(nspin, ef,esmr, qbz,wbz, nband,nqbz)
  write(6,*)' ef    =',ef
  write(6,*)' esmr  =',esmr
  write(6,*)' valn  =',valn
  !      write(6,*)' ntot  =',ntot
  ! 2001 May
  !        if(abs(valn-ntot)>1d-6) stop ' abs(valn-ntot)>1d-6'


  ! --- Ex test --- Calculate weights.
  INQUIRE (FILE = 'EXspTEST', EXIST = exsptest)
  if( .NOT. exsptest) goto 999
  ! -------------------
  !> input files
  !      if(readgwinput()) then
  call getkeyvalue("GWinput","<QPNT>",unit=ifqpnt,status=ret)
  !      else
  !      ifqpnt    = iopen('QPNT',1,0,0)
  !      endif

  ! read q-points and states
  lqall      = .false.
  laf        = .false.
  call readx   (ifqpnt,10)
  read (ifqpnt,*) iqall,iaf
  if (iqall == 1) lqall = .TRUE. 
  if (iaf   == 1)   laf = .TRUE. 
  call readx   (ifqpnt,100)
  ! states
  read (ifqpnt,*) ntq
  allocate( itq(ntq) )
  read (ifqpnt,*) (itq(i),i=1,ntq)
  ! q-points
  if (lqall) then !all q-points case
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
  ! -------------------
  do is=1,nspin
     do ipex = 1,nq
        do itpex=1,ntq
           filenameex = 'EXSS'//charnum3(ipex)//charnum3(itpex) &
                //'.'//char(48+is)
           ifexspxx =ifile_handle()
           open(ifexspxx, file=filenameex)
           filenameex='EXWT'//charnum3(ipex)//charnum3(itpex) &
                //'.'//char(48+is)
           ifexspxxw=ifile_handle()
           open(ifexspxxw,file=filenameex)
           ex=0d0
           do
              read(ifexspxx,*,end=1013) &
                   eee, exwgt
              wfac = wfacx(-1d99, ef, eee, esmr)
              ex=ex + wfac*exwgt
              write(ifexspxxw,"(4f15.8)") eee, wfac, &
                   exwgt,wfac*exwgt
              if( wfac==0d0 ) exit
           enddo
1013       write(6,"(' j iq isp=' i3,i4,i2,'  q=',3f8.4, &
                '  Sx(eV)=',f10.4)") &
                itq(itpex),ipex,is, q(1:3,ipex), &
                ex
           close(ifexspxx)
           close(ifexspxxw)
        enddo
     enddo
  enddo
  ! top2rx 2013.08.09 kino 999  stop '--- Efermi and EX test --- '
999 call rx( '--- Efermi and EX test --- ')
END PROGRAM hef

