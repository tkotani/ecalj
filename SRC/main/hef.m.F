      program hef
c-------------------------------------------------------------
c hef is from hsfp0 but it is only to calcuate ef and Ex test.
c so this routine contains unused things so much.
      use m_hamindex,only:   Readhamindex
      use m_readeigen,only:init_readeigen,readeval
      use m_read_bzdata,only: Read_bzdata,ginv,nqibz,nqbz,qibz,wibz,qbz,wbz
      use m_genallcf_v3,only: genallcf_v3,
     & nclass,natom,nspin,nl,nn,nnv,nnc, 
     & nlmto,nlnmx, nctot,niw,!nw_input=>nw,ef
     & alat, delta,deltaw,esmr,clabl,iclass, !diw,dw,
     & il,in,im,nlnm, 
     & plat, pos,z,ecore,konf,nlnx
      use m_keyvalue,only: getkeyvalue
      use m_readhbe,only: Readhbe, nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg

      implicit none
c
c      real(8) :: shtw
      integer(4) :: mxclass,ngnmax
c ,mbytes,mwords,iwksize,
c     &   natom,nclass,ipos,ngrp,igrp,
c     &   iinvg,instar,iirk,
c     o   nspin,nl,nn,nnv,nnc,
c     o   inindx,inindxv,inindxc,iiclass,
c     d   nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc,
c     o   iz,
c     o   iil,iin,iim,iilnm,inlnm,
c     o   iilv,iinv,iimv,iilnmv,inlnmv,
c     o   iilc,iinc,iimc,iilnmc,inlnmc,
c     o   incwf,iecore,ikonf,iicore,incore,nctot,
c     o   imagw,niw,nw,ifreq,  
     &   ixc,!ifhbed, nprecb,mrecb,mrece,nlmtot,nqbzt, nband, !iopen,
     &   ibas,ibasx,ngpmx,nxx,ngcmx,nbloch,ifqpnt,ifwd,
     &   nprecx,nblochpmx2,nwt,niwt, nqnum,mdimx,nblochpmx,
     &   ifrcw,ifrcwi,  noccxv,maxocc,noccx,ifvcfpout,iqall,iaf,ntq,
     &   i,k,nspinmx, nq,is,ip,iq,idxk,ifoutsex,iclose,nq0i,ig,iimdim,
     &   ifreqx,ifreqw,iwx,iexpa,mxkp,nqibzxx,ntet,nene,iqi, ix,iw,
     &   nlnx4,niwx,irot,invr,invrot,ivsum, ifoutsec,ntqx,
c
c     &   ifrb(2),ifcb(2),ifrhb(2),ifchb(2),ifev(2),ifsec(2)
     &   ifev(2),
c     &             ,ifxc(2),ifsex(2), ifphiv(2),ifphic(2),ifec,
     &   ndble=8
c
c      real(8) :: alat,ef,diw,dw,delta,pi,tpia,vol,voltot,rs,alpha,
      real(8) :: pi,tpia,vol,voltot,rs,alpha,
     &   qfermi,efx,valn,efnew,edummy,efz,qm,xsex,
     &   zfac1,zfac2,dscdw1,dscdw2,dscdw,zfac
c
c     &   plat(3,3)

      logical lqall,laf
c      character*120 symgrp

c class parameters
c      parameter (mxclass=100)
c      character*6 clabl(mxclass)
c symmetry group
c      parameter (ngnmax=10)
c      real(8) :: gen(9,ngnmax)
c
      integer(4),allocatable :: itq(:)
      real(8),allocatable    :: q(:,:)

c takao
      integer(4),allocatable :: ngvecpB(:,:,:),ngveccB(:,:,:),
     &  ngvecp(:,:), ngvecc(:,:),ngpn(:),ngcni(:),iqib(:),
     &  kount(:,:), nx(:,:),nblocha(:),lx(:),ngveccBr(:,:,:)
      real(8),allocatable:: vxcfp(:,:,:),
     &  wqt(:), wgt0(:,:),q0i(:,:),
     &  ppbrd (:,:,:,:,:,:,:),cgr(:,:,:,:),eqt(:),
     &  ppbrdx(:,:,:,:,:,:,:),aaa(:,:),symope(:,:,:), 
     &  ppb(:),pdb(:),dpb(:),ddb(:), eq(:,:),
     &  eqx(:,:,:),eqx0(:,:,:),ekc(:),coh(:,:)
      complex(8),allocatable:: geigB(:,:,:,:) ,zsec(:,:,:)
c
      logical :: exchange, legas=.false.
      real(8):: qreal(3), ntot,nocctotg2,tripl
      logical ::nocore

c space group infermation
      integer(4),allocatable :: iclasst(:), invgx(:), miat(:,:)
      real(8),allocatable    :: tiat(:,:,:),shtvg(:,:)

c tetra
      real(8),allocatable :: qz(:,:),qbzxx(:),wbzxx(:),wtet(:,:,:,:),
     &     eband(:,:,:), ene(:) !,ecore(:,:)
      integer(4),allocatable ::idtetx(:,:),idtet(:,:),ipq(:)
     &      ,iene(:,:,:),ibzx(:)
c      real(8) :: qbasmc(3,3)

c worksize in megabytes (1 word = 4 bytes)
c      parameter (mbytes=60)
c      parameter (mwords=mbytes/4)
c      parameter (iwksize=mwords * 1000 *1000)
c      integer w
c      common /w/ w(iwksize)

      integer(4) ::ib,iqx,igp,iii,isx
      integer(4) :: ipex,itpex,ifexspxx,ifexspxxw
      character(12) :: filenameex
      character(3) :: charnum3
      real(8) :: ex,eee_dummy,exwgt_dummy,eee,exwgt,wfac,wfacx
      logical :: exsptest

      integer(4):: incwfin,ret,ifile_handle
      real(8)::ef
c      logical:: readgwinput

c      call headver('hef',0)

c-----------------------------------------------------------------------
c---  readin BZDATA. See gwsrc/rwbzdata.f
c--------readin data set when you call read_BZDATA ---------------
c       integer(4)::ngrp,nqbz,nqibz,nqbzw,nteti,ntetf,
c     &   n_index_qbz
c       real(8):: qbas(3,3),ginv(3,3)
c       real(8),allocatable:: qbz(:,:),wbz(:),qibz(:,:)
c     &    ,wibz(:),qbzw(:,:)
c       integer(4),allocatable:: idtetf(:,:),ib1bz(:),idteti(:,:)
c     &    ,nstar(:),irk(:,:),index_qbz(:,:,:)
c-----------------------------------------------------------------
      call read_BZDATA()
      write(6,*)' nqbz qbz =',nqbz
c      write(6,*)' nqibz ngrp=',nqibz,ngrp2
      call pshprt(60)
c-------------------------------
c generate all quantities needed
c-------------------------------
      incwfin=  0     ! dummy
      call genallcf_v3(incwfin) !in module m_genallcf_v3
c      if(ngrp/= ngrp2) call rx( 'ngrp inconsistent: BZDATA and LMTO GWIN_V2')
c-------------------------------------------------------------------
      if (nclass /= natom ) call rx( ' hsfp0: nclass /= natom ')
      write(6,*)' hsfp0: end of genallcf2'
      call pshprt(30)
      pi   = 4d0*datan(1d0)
      tpia = 2d0*pi/alat
      voltot = abs(alat**3*tripl(plat,plat(1,2),plat(1,3)))
      call Readhbe()
      call Readhamindex()
      call init_readeigen()!nband,mrece) !initialization of readEigen

c --- determine Fermi energy ef for given valn (legas case) or corresponding charge given by z and konf.
! When esmr is negative, esmr is geven automatically by efsimplef.
      call efsimplef2a(nspin,wibz,qibz,ginv, 
     i        nband,nqibz
     i       ,konf,z,nl,natom,iclass,nclass
     i       ,valn, legas, esmr,  !!! valn is input for legas=T, output otherwise.
c
     i        qbz,nqbz !index_qbz, n_index_qbz,
     o       ,efnew)
      ef = efnew
c- check total ele number -------
      ntot  = nocctotg2(nspin, ef,esmr, qbz,wbz, nband,nqbz)
      write(6,*)' ef    =',ef
      write(6,*)' esmr  =',esmr
      write(6,*)' valn  =',valn
      write(6,*)' ntot  =',ntot
c 2001 May
c        if(abs(valn-ntot)>1d-6) stop ' abs(valn-ntot)>1d-6'


c --- Ex test --- Calculate weights.
      INQUIRE (FILE = 'EXspTEST', EXIST = exsptest)
      if(.not.exsptest) goto 999
c -------------------
c> input files
c      if(readgwinput()) then
      call getkeyvalue("GWinput","<QPNT>",unit=ifqpnt,status=ret)
c      else
c      ifqpnt    = iopen('QPNT',1,0,0)
c      endif

c read q-points and states
      lqall      = .false.
      laf        = .false.
      call readx   (ifqpnt,10)
      read (ifqpnt,*) iqall,iaf
      if (iqall == 1) lqall = .true.
      if (iaf   == 1)   laf = .true.
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
c -------------------
      do is=1,nspin
        do ipex = 1,nq
          do itpex=1,ntq
            filenameex = 'EXSS'//charnum3(ipex)//charnum3(itpex)
     &    //'.'//char(48+is)
            ifexspxx =ifile_handle()
            open(ifexspxx, file=filenameex)
            filenameex='EXWT'//charnum3(ipex)//charnum3(itpex)
     &    //'.'//char(48+is)
            ifexspxxw=ifile_handle()
            open(ifexspxxw,file=filenameex)
            ex=0d0
            do
              read(ifexspxx,*,end=1013)
     &    eee, exwgt
              wfac = wfacx(-1d99, ef, eee, esmr)
              ex=ex + wfac*exwgt
              write(ifexspxxw,"(4f15.8)") eee, wfac,
     &    exwgt,wfac*exwgt
              if( wfac==0d0 ) exit
            enddo
 1013       write(6,"(' j iq isp=' i3,i4,i2,'  q=',3f8.4,
     &  '  Sx(eV)=',f10.4)") 
     &  itq(itpex),ipex,is, q(1:3,ipex),
     &  ex
            close(ifexspxx)
            close(ifexspxxw)
          enddo
        enddo
      enddo
Cstop2rx 2013.08.09 kino 999  stop '--- Efermi and EX test --- '
 999  call rx( '--- Efermi and EX test --- ')
      end

