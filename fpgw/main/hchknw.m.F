      program hchknw
c---     check the number of nw
      use m_readeigen,only: init_readeigen,readeval
      use m_read_bzdata,ngrp2=>ngrp
      use m_genallcf_v3
      use keyvalue
      implicit none
      logical :: tetra=.false.
      real(8) :: shtw=0d0
      integer(4) :: mxclass,ngnmax,mbytes,mwords,iwksize,
c     &   natom,nclass,ipos,igrp,
c     &   iindxk,
c     &   iinvg,iirk,
c     o   nspin,nl,nn,nnv,nnc,
c     o   inindx,inindxv,inindxc,iiclass,
c     d   nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc,
c     o   iz,
c     o   iil,iin,iim,iilnm,inlnm,
c     o   iilv,iinv,iimv,iilnmv,inlnmv,
c     o   iilc,iinc,iimc,iilnmc,inlnmc,
c     o   incwf,iecore,ikonf,iicore,incore,nctot,
c     o   imagw,niw,nw,ifreq,  
     &  incwf,imagw,
c 
     &   iii,iopen,ifhbed, nprecb,mrecb,mrece,nlmtot,nqbzt, nband,
     &   ibas,ibasx,ngpmx,nxx,ngcmx,nbloch,ifqpnt,ifwd,
     &   nprecx,nblochpmx2,nwt,niwt, nqnum,mdimx,nblochpmx,
     &   ifrcw,ifrcwi,  noccxv,maxocc,noccx,ifvcfpout,iqall,iaf,ntq,
     &   i,k,nspinmx, nq,is,ip,iq,idxk,ifoutsex,iclose,nq0i,ig,iimdim,
     &   ifreqx,ifreqw,iwx,iexpa,mxkp,nqibzxx,ntet,nene,iqi, ix,iw,
     &   nlnx4,niwx,irot,invr,invrot,ivsum, ifoutsec,ntqx,
c
     &   ifrb(2),ifcb(2),ifrhb(2),ifchb(2),ifsec(2)           !,ifev(2)
     &             ,ifxc(2),ifsex(2),ifphi_v2, ifphic(2),ifec,
     &   ndble=8
c
      real(8) :: pi,tpia,vol,voltot,rs,alpha,
     &   qfermi,efx,valn,efnew,edummy,efz,qm,xsex,egex,
     &   zfac1,zfac2,dscdw1,dscdw2,dscdw,zfac
c
c     &   plat(3,3)

      logical lqall,laf
c      character*120 symgrp

c class parameters
      parameter (mxclass=100)
c      character*6 clabl(mxclass)
c symmetry group
      parameter (ngnmax=10)
      real(8) :: gen(9,ngnmax)
c
      integer(4),allocatable :: itq(:)
      real(8),allocatable    :: q(:,:)

c takao
      integer(4),allocatable :: ngvecpB(:,:,:),ngveccB(:,:,:),
     &  ngvecp(:,:), ngvecc(:,:),ngpn(:),ngcni(:),iqib(:),
     &  kount(:), nx(:,:),nblocha(:),lx(:),ngveccBr(:,:,:)
      real(8),allocatable:: vxcfp(:,:),
     &  wqt(:), wgt0(:,:),q0i(:,:),
     &  ppbrd (:,:,:,:,:,:,:),cgr(:,:,:,:),eqt(:),
     &  ppbrdx(:,:,:,:,:,:,:),symope(:,:,:) ,   !qibz(:,:),
     &  ppb(:),pdb(:),dpb(:),ddb(:), eq(:,:),
     &  eqx(:,:,:),eqx0(:,:,:),ekc(:)
      complex(8),allocatable:: geigB(:,:,:) ,geig(:,:),zsec(:,:,:)
c
      logical :: exchange, legas = .false. ! test forelectron gas exchange check write for electron gas test only for exchange.
      real(8):: qreal(3), ntot,nocctotg2,tripl,xxx(3,3)

c space group infermation
      integer(4),allocatable :: iclasst(:), invgx(:), miat(:,:)
      real(8),allocatable    :: tiat(:,:,:),shtvg(:,:)

c tetra
      real(8),allocatable :: qz(:,:),qbzxx(:),wbzxx(:),wtet(:,:,:,:),
     &     eband(:,:,:), ene(:) !,ecore(:,:)
      integer(4),allocatable ::idtetx(:,:),idtet(:,:),ipq(:)
     &      ,iene(:,:,:),ibzx(:)  !,nstar(:)
c      real(8) :: qbasmc(3,3)

c worksize in megabytes (1 word = 4 bytes)
c      parameter (mbytes=0)
c      parameter (mwords=mbytes/4)
c      parameter (iwksize=mwords * 1000 *1000)
c      integer w
c      common /w/ w(iwksize)

      real(8) :: eee,omegav,omegac,eminc,emaxv,we
      integer(4) ::it,ib,ipx1,itx1,ipx2,itx2,iflegas, ifphiv_2

      real(8),allocatable:: qbze(:,:)
      integer(4):: nq0it,nqbze,ini
      real(8):: wdummy
      logical:: GaussSmear !external

      real(8):: ffac
      logical ::qpntexist
      integer(4)::ixc,nz,idummy

      integer(4)::nwin, incwfin,ret
      real(8)::efin

      integer(4):: nss(2)
      logical:: debug=.true.

c---------------------------------------------------------------
      INQUIRE (FILE = 'LEGAS', EXIST = legas)
      write(6,*) ' --- Choose modes below ------------'
      write(6,*) '  From QPNT (0 or return). Fullband(1) '
      write(6,*) ' --- Put number above ! ------------'
      call readin5(ixc,nz,idummy)
      call headver('hchknw',ixc)
      if(ixc/=0.and.ixc/=1) stop "hchknw:Wrong mode!"

c---  readin BZDATA. See gwsrc/rwbzdata.f
      call read_BZDATA()
c--------readin data set when you call read_BZDATA ---------------
c       integer(4)::ngrp,nqbz,nqibz,nqbzw,nteti,ntetf,
c     &   n_index_qbz
c       integer(4):: n1,n2,n3
c       real(8):: qbas(3,3),ginv(3,3),qbasmc(3,3)
c       real(8),allocatable:: qbz(:,:),wbz(:),qibz(:,:)
c     &    ,wibz(:),qbzw(:,:)
c       integer(4),allocatable:: idtetf(:,:),ib1bz(:),idteti(:,:)
c     &    ,nstar(:),irk(:,:),index_qbz(:,:,:)
c-----------------------------------------------------------------
      print *,' nqbz qbz =',nqbz
c      print *,  qbz
      print *,' nqibz ngrp=',nqibz,ngrp
c      print *,' qout=',qout
c      stop 'test end qout'
c      print *,' irk=',irk
c      print *,' #### idtetf: ####'
c      print *, idtetf
c------------------------------------------------------------------

c set up work array
c      call wkinit (iwksize)
      call pshprt(60)


C--- readin GWIN and LMTO, then allocate and set datas.
      nwin    = -999 !Not readin NW file
      efin    = 0d0 !readin EFERMI
      if(legas) efin=-999d0! not readin EFERMI
      incwfin = 0 !dummy
      call genallcf_v3(nwin,efin,incwfin) !in module m_genallcf_v3
      if(ngrp/= ngrp2) stop 'ngrp inconsistent: BZDATA and LMTO GWIN_V2'
c--- These are allocated and setted.
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
c-----------------------------------------------------------------------
c$$$      call genallcf_v2x(
c$$$c> structure
c$$$     o                   plat,alat,natom,nclass,ipos,
c$$$c> symmetry
c$$$     o                   symgrp,gen,ngnmax,ngrp,igrp,
c$$$     o                   iinvg,ef,
c$$$c>> l,n and dimensions
c$$$     o                   clabl,nspin,nl,nn,nnv,nnc,
c$$$     o                   inindx,inindxv,inindxc,iiclass,
c$$$     d                   nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc,
c$$$     o                   iz,
c$$$c>> l,n,m indices for Phi
c$$$     o                   iil,iin,iim,iilnm,inlnm,
c$$$     o                   iilv,iinv,iimv,iilnmv,inlnmv,
c$$$     o                   iilc,iinc,iimc,iilnmc,inlnmc,
c$$$c>> core
c$$$     o                   incwf,iecore,ikonf,iicore,incore,nctot,
c$$$c> frequency
c$$$     o                   imagw,niw,diw,nw,dw,delta,deltaw,esmr,ifreq )


c=== dump check ===============================
      gen  =0d0 !dummy
      imagw=0   !dummy
      call genallcf_v3_dump(
     o          plat,alat,natom,nclass, pos,             !structure
     o          symgrp,gen, ngnmax,ngrp, symgg, !igrp,   !symmetry
     o          invg, ef,
c>> l,n and dimensions
     o           clabl,nspin,nl,nn,nnv,nnc,
     o           nindx, nindxv, nindxc, iclass,
     d           nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc,
     o           z,
c>> l,n,m indices for Phi
     o                   il, in, im, ilnm, nlnm,
     o                   ilv,inv,imv,ilnmv,nlnmv,
     o                   ilc,inc,imc,ilnmc,nlnmc,
     o      ncwf, ecore, konf,icore,ncore,nctot,         !core
     o      imagw,niw,diw,nw,dw,delta,deltaw,esmr,freq ) ! frequency
c=== dump test end ==================================================


      if (nclass > mxclass) stop ' hchknw: increase mxclass'
      if (nclass /= natom ) stop ' hchknw: nclass /= natom ' ! We assume nclass = natom.
      print *,' hchknw: end of genallc'
c
      pi   = 4d0*datan(1d0)
      tpia = 2d0*pi/alat
      call dinv33(plat,1,xxx,vol)
      voltot = dabs(vol)*(alat**3)

c--- !!! ef is taken as rs for the empty-sphere test case of legas=T case -------------
      if(legas) then
        iflegas=2101
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
      endif

c------------
c open files
c------------
c eigen values ------------------------------
c>> read dimensions of wc,b,hb
      ifhbed     = iopen('hbe.d',1,0,0)
      read (ifhbed,*) nprecb,mrecb,mrece,nlmtot,nqbzt, nband
      if (nprecb == 4) stop 'hsfp0: b,hb in single precision'
c
c      ifev(1)   = iopen('EVU', 0,0,mrece)
c      if (nspin == 2) ifev(2) = iopen('EVD', 0,0,mrece)


c read q-points and states ----------------------------
c      INQUIRE (FILE ='QPNT', EXIST = qpntexist)
c      if(.not.qpntexist) then
      if(ixc==1) then
        write(6,*) ' Take all bands ank qpnt'
        lqall  = .true.
        INQUIRE (FILE ='ANFcond', EXIST = laf)
        if(laf) print *,' Find ANFcond--- laf=T only UP band'
        ntq=-1
      else
c        ifqpnt = iopen('QPNT',1,0,0)
c        if(readgwinput()) then
        call getkeyvalue("GWinput","<QPNT>",unit=ifqpnt,status=ret)
        print *,' end of getkeyvalue------ ifqpnt=',ifqpnt
c        else
c        ifqpnt    = iopen('QPNT',1,0,0)
c        endif

        lqall      = .false.
        laf        = .false.
        call readx   (ifqpnt,10)
        read (ifqpnt,*) iqall,iaf
        if (iqall == 1) lqall = .true.
        if (iaf   == 1)   laf = .true.
        call readx   (ifqpnt,100)
        read (ifqpnt,*) ntq
      endif

! states
      if(ntq<0) then  !takao May 2003
        ntq = nband
        allocate (itq(ntq))
        do i = 1, ntq
          itq(i) = i
        enddo
      else
        allocate( itq(ntq) )
        read (ifqpnt,*) (itq(i),i=1,ntq)
      endif

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
      if(ixc/=1) close(ifqpnt) !Walter fixed this Jul2008
      nspinmx = nspin
      if (laf) nspinmx =1

C...
      call getkeyvalue("GWinput","QPNT_nbandrange",nss,2 !nov2005
     &               ,default=(/-99997,-99997/) )
      if(nss(2)/=-99997) then
        if(allocated(itq)) deallocate(itq)
        ntq=nss(2)-nss(1)+1
        allocate( itq(ntq) )
        do i=max(1,nss(1)),min(nss(2),nband)
          itq(i-nss(1)+1) = i
        enddo
      endif

      if(debug)print *,' goto init_readeigen'
      call init_readeigen(ginv,nspin,nband,mrece) !initialization of readEigen
      if(debug)print *,' end of init_readeigen'
c----------------------------------------------
      if(.not.tetra) then
        if(esmr/=0d0) then
          call efsimplef2a(nspin,wibz,qibz,ginv, 
     i        nband,nqibz
c     i       ,w(ikonf),w(iz),nl,natom,w(iiclass),nclass
     i       ,konf,z,nl,natom,iclass,nclass
     i       ,valn, legas, esmr,  !!! valn is input for legas=T, output otherwise.
c
     i       qbz,nqbz ! index_qbz, n_index_qbz,
     o       ,efnew)
          ef = efnew
c- check total ele number -------
          ntot  = nocctotg2(nspin, ef,esmr, qbz,wbz, nband,nqbz)
          print *,' ef    =',ef
          print *,' esmr=',esmr
          print *,' valn  =',valn
          print *,' ntot  =',ntot
ccc
ccc        if(abs(valn-ntot)>1d-6) stop ' abs(valn-ntot)>1d-6'  !20001 May
c      else
c        if(legas) ef = efx
        endif
      endif



c ---  q near zero
      print *, 'reading QOP'
      open (101,file='Q0P')
      read (101,"(i5)") nq0it
      write(6,*) ' *** nqibz nq0i_total=', nqibz,nq0it
      allocate( q0i(1:3,1:nq0it) )
      do i=1,nq0it
        read (101,* ) wdummy,q0i(1:3,i)
      enddo
      write(6,*) ' k in Q0P =', nq0it
      write(6,"(i3, 3f14.6)" )(i, q0i(1:3,i),i=1,nq0it)
      close(101)
c --- qibze(3,nqbze) qbze(3,nqibze)
      nqbze  = nqbz *(1 + nq0it)
      allocate(  qbze(3,nqbze) )
      call dcopy(3*nqbz, qbz, 1, qbze,1)
      do i = 1,nq0it
        ini = nqbz*(1 + i -1)
        do ix=1,nqbz
          qbze (:,ini+ix) = q0i(:,i) + qbze(:,ix)
        enddo
      enddo
c

c --- read LDA eigenvalues for occ
c     ntp0=ntq
      allocate(eqt(nband))
      omegav =  0d0 !1d99 fixed oct.2003 takao
      omegac =  0d0 !-1d99
      do is = 1,nspinmx
        do ip = 1,nq
c        iq    = idxk (q(1,ip),qbze,nqbze)
c        call rwdd1   (ifev(is), iq, nband, eqt) !direct access read b,hb and e(q,t)
          call readeval(q(1,ip),is,eqt)
          do it=1,ntq
            eee = eqt(itq(it)) - 2d0*deltaw  !2.d0*(-1d0-shtw)*deltaw
            print *,' is ip it eee=',eee,eqt(itq(it))
            if(eee>=1d20-1d10) cycle !takao jun2009
            if( eee <ef .and. eee< omegav )  then
              omegav = eee
              ipx1 = ip
              itx1 = it
            endif
            eee = eqt(itq(it)) + 2d0*deltaw  !+ 2.d0*(1d0-shtw)*deltaw
c
            print *,' is ip it eee=',eee,eqt(itq(it))
            if( eee >ef .and. eee> omegac )  then
              omegac = eee
              ipx2 = ip
              itx2 = it
            endif
          enddo
        enddo
      enddo

      if(GaussSmear()) then
        ffac=10d0 !This is OK?
      else
        ffac=0.5d0
      endif
      emaxv =  0d0 !-1d99 fixed oct.2003 takao
      eminc =  0d0 !1d99
      do is = 1, nspinmx
        do iq = 1, nqbz
c       call rwdd1 (ifev(is), iq, nband, eqt)
          call readeval(qbz(1,iq),is,eqt)
          print *
          do ib=1,nband
            eee = eqt(ib)
            if(ib<6) print *,' iq ib eee=',iq,ib,eee
            if( eee <ef+ffac*esmr .and. eee>emaxv) emaxv = eee
            if( eee >ef-ffac*esmr .and. eee<eminc) eminc = eee
          enddo
        enddo
      enddo
      deallocate(eqt)
      we  = max(abs(emaxv - ef), abs(omegav-ef)
     &         ,abs(omegac- ef) , abs(eminc-ef) )
      we= we+ffac*esmr
      nw  = idint (we/2d0/dw) + 3
      print *
      print *,' --------------------------------'
      print *,' emaxv= ',emaxv
      print *,' eminc= ',eminc
      print *,' omegav ip it=',omegav,ipx1,itx1
      print *,' omegac ip it=',omegac,ipx2,itx2
      print *,' we max for <ef', emaxv - omegav
      print *,' we max for >ef', omegac- eminc
      print *,' *** required we nw=  ', we,nw
      print *,"write nw to NW"
      open(1101,file='NW')
      write(1101,*) nw
      close(1101)

c
      write(6, *) ' --- computational conditions --- '
      write(6,'("    deltaw  =",f13.6)') deltaw
      write(6,'("    esmr    =",f13.6)') esmr
      write(6,'("    alat voltot =",2f13.6)') alat, voltot
      write(6,'("    niw nw dw   =",2i6,f13.6)') niw,nw,dw
      stop " OK! hchknw: write nw to NW"
c
      end
