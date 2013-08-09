      program hef
c-------------------------------------------------------------
c hef is from hsfp0 but it is only to calcuate ef and Ex test.
c so this routine contains unused things so much.
      use m_readeigen,only:init_readeigen,readeval
      use m_read_bzdata,ngrp2=>ngrp
      use m_genallcf_v3
      use keyvalue
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
     &   ixc,iopen,ifhbed, nprecb,mrecb,mrece,nlmtot,nqbzt, nband,
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
      logical :: exchange, legas
      real(8):: qreal(3), ntot,nocctotg2,tripl,xxx(3,3)
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

      integer(4):: nwin,incwfin,ret
      real(8)::efin
c      logical:: readgwinput

      call headver('hef',0)

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
      print *,' nqbz qbz =',nqbz
c      print *,  qbz
      print *,' nqibz ngrp=',nqibz,ngrp2
c      print *,' irk=',irk
c      print *,' #### idtetf: ####'
c      print *, idtetf

c set up work array
c      call wkinit (iwksize)
      call pshprt(60)
c-------------------------------
c generate all quantities needed
c-------------------------------
      nwin   = -999   ! dummy
      incwfin=  0     ! dummy
      efin=  -999d0   ! dummy
C--- readin GWIN and LMTO, then allocate and set datas.
      call genallcf_v3(nwin,efin,incwfin) !in module m_genallcf_v3
Cstop2rx 2013.08.09 kino      if(ngrp/= ngrp2) stop 'ngrp inconsistent: BZDATA and LMTO GWIN_V2'
      if(ngrp/= ngrp2) call rx( 'ngrp inconsistent: BZDATA and LMTO GWIN_V2')
c------------------------------
c      integer(4)::
c     o  nclass,natom,nspin,nl,nn,nnv,nnc,
c     o  ngrp,
c     o  nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc,
c     o  nctot,
c     o  niw,nw
c      real(8) :: alat,ef, diw,dw,delta,deltaw,esmr
c      character(120):: symgrp
c      character(6),pointer :: clabl(:)
c      integer(4),pointer::
c     &   iclass(:)
c     &  ,nindxv(:,:),nindxc(:,:)
c     &  ,ncwf(:,:,:) ,
c     o    invg(:),
c     o    il(:,:), in(:,:), im(:,:),   ilnm(:),  nlnm(:),
c     o    ilv(:),inv(:),imv(:),  ilnmv(:), nlnmv(:),
c     o    ilc(:),inc(:),imc(:),  ilnmc(:), nlnmc(:),
c     o    nindx(:,:),konf(:,:),icore(:,:),ncore(:),
c     &    occv(:,:,:),unoccv(:,:,:)
c     &   ,occc(:,:,:),unoccc(:,:,:),
c     o    nocc(:,:,:),nunocc(:,:,:)
c      real(8), pointer::
c     o  plat(:,:),pos(:,:),z(:), symgg(:,:,:), !w(igrp)
c     o  ecore(:,:), freq(:)
c--------------------------------
c
c$$$      call genallcf_v2x(
c$$$c> structure
c$$$     o                   plat,alat,natom,nclass,ipos,
c$$$c> symmetry
c$$$     o                   symgrp,gen,ngnmax,ngrp,igrp,
c$$$c> BZ
c$$$c     o                   n1,n2,n3,qbas,ginv,iqibz,iwibz,nqibz,
c$$$c     o                   iqbz,iwbz,nqbz,   iindxk,
c$$$c     o                   iinvg,  instar,iirk,   ef,
c$$$     o                   iinvg,ef,
c$$$c>> file units
c$$$c     o                   ifphiv,ifphic,ifec,
c$$$c>> l,n and dimensions
c$$$     o                   clabl,nspin,nl,nn,nnv,nnc,
c$$$     o                   inindx,inindxv,inindxc,iiclass,
c$$$     d                   nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc,
c$$$c z value
c$$$     o                   iz,
c$$$c>> l,n,m indices for Phi
c$$$     o                   iil,iin,iim,iilnm,inlnm,
c$$$     o                   iilv,iinv,iimv,iilnmv,inlnmv,
c$$$     o                   iilc,iinc,iimc,iilnmc,inlnmc,
c$$$c>> core
c$$$     o                   incwf,iecore,ikonf,iicore,incore,nctot,
c$$$c> frequency
c$$$     o                   imagw,niw,diw,nw,dw,delta,deltaw,esmr,ifreq )
c$$$c

c      allocate(ecore(nctot,nspin)) !core energies
c      do is = 1,nspin
c        if (nctot > 0) call catch1 (w(iecore),is,nctot,2,ecore(:,is)) !core energies
c      enddo

c-------------------------------------------------------------------
c      if (nclass > mxclass) stop ' hsfp0: increase mxclass'
c!!!! WE ASSUME iclass(iatom)= iatom !!!!!!!!!!!!!!!!!!!!!!!!!
Cstop2rx 2013.08.09 kino      if (nclass /= natom ) stop ' hsfp0: nclass /= natom ' ! We assume nclass = natom.
      if (nclass /= natom ) call rx( ' hsfp0: nclass /= natom ')
      print *,' hsfp0: end of genallcf2'
c
      call pshprt(30)
      pi   = 4d0*datan(1d0)
      tpia = 2d0*pi/alat

c      shtw = 0d0
c      if(esmr<1d-5) shtw=0.01d0 ! Ferdi's shift to avoid resonance effect(maybe)
c
      call dinv33(plat,1,xxx,vol)
      voltot = dabs(vol)*(alat**3)

c>> read dimensions of wc,b,hb
      ifhbed     = iopen('hbe.d',1,0,0)
      read (ifhbed,*) nprecb,mrecb,mrece,nlmtot,nqbzt, nband
Cstop2rx 2013.08.09 kino      if (nprecb == 4) stop 'hsfp0: b,hb in single precision'
      if (nprecb == 4) call rx( 'hsfp0: b,hb in single precision')



      call init_readeigen(ginv,nspin,nband,mrece) !initialization of readEigen

c------------
c open files
c------------
c eigen functions
c      ifev(1)   = iopen('EVU', 0,0,mrece)
c      if (nspin == 2) ifev(2) = iopen('EVD', 0,0,mrece)

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
      print *,' ef    =',ef
      print *,' esmr  =',esmr
      print *,' valn  =',valn
      print *,' ntot  =',ntot
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
            ifexspxx =4112
            open(ifexspxx, file=filenameex)
            filenameex='EXWT'//charnum3(ipex)//charnum3(itpex)
     &    //'.'//char(48+is)
            ifexspxxw=4113
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

