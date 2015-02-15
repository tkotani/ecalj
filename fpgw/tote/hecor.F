      program hecor
c-------------------------------------------------------------
cm revised for fpgw025, May 07, 2002, Takashi Miyake
c
c Calculates the RPA correlation part in the total energy
c from gwsrc/hsfp0.f, July 05, 2001, Takashi Miyake
c------------------------------------------------------------
c      implicit real*8(a-h,o-z)
      use m_readeigen,only: init_readeigen,init_readeigen2,readeval,lowesteval
      use m_readqg
      use m_read_bzdata,ngrp2=>ngrp
      use m_genallcf_v3
      use keyvalue
      implicit none
c parameter
      real(8),parameter :: 
     &  ua    = 1d0    ! constant in w(0)exp(-ua^2*w'^2) to take care of peak around w'=0
c------------------------------------
      logical ::         
     &    tetra  = .false.  ! test switch for tetrahedron method test.
c                           ! tetra=T is only effective for exchange=T case.
c                           ! Tetrahedron mehod for correlation is a bit
      ! difficult and I gave up.
c                           ! If you want to calculate with tetra=T for exchange, you
      ! have to uncomment tetra related part in
      ! sxcf.f, and a part calling sxcf in this routine. Note wtet wtetef!
      ! They sometimes cause array destruction if you run tetra=T without comment them.
c------------------------------------
      integer(4) ::
     o   nqitot,
     &   ixc,iopen,ifhbed, nprecb,mrecb,mrece,nlmtot,nqbzt, nband,
     &   ibas,ibasx,ngpmx,nxx,ngcmx,nbloch,ifqpnt,ifwd,
     &   nprecx,mrecl,nblochpmx2,nwt,niwt, nqnum,mdimx,nblochpmx,
     &   ifrcw,ifrcwi,  noccxv,maxocc,noccx,ifvcfpout,iqall,iaf,ntq,
     &   i,k,nspinmx, nq,is,ip,iq,idxk,ifoutsex,iclose,nq0i,ig,iimdim,
     &   ifreqx,ifreqw,iwx,iexpa,mxkp,nqibzxx,ntet,nene,iqi, ix,iw,
     &   nlnx4,niwx,irot,invr,invrot,ivsum, ifoutsec,ntqx,
c
c     &   ifrb(2),ifcb(2),ifrhb(2),ifchb(2),
c     &   ifev(2),
c     &   ifsec(2),ifxc(2),ifsex(2), 
     &   ifphiv(2),ifphic(2),ifec,ifexsp(2),ifcor,iflegas,
     &   ndble=8,ndummy1,ndummy2,ngb
c
      real(8) :: pi,tpia,vol,voltot,rs,alpha,
     &   qfermi,efx,valn,efnew,edummy,efz,qm,xsex,egex,
     &   zfac1,zfac2,dscdw1,dscdw2,dscdw,zfac,
c
c     &   plat(3,3),qbas(3,3),ginv(3,3),lowesteb,
     &   totexc,trpv,trlog,ecelgas

c      logical lqall,laf
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
      integer(4),allocatable :: ngvecpB(:,:,:),
     &  ngvecp(:,:), iqib(:),
     &  kount(:,:), nx(:,:),nblocha(:),lx(:)
      real(8),allocatable:: vxcfp(:,:,:),
     &  wqt(:), wgt0(:,:),q0i(:,:),
     &  ppbrd (:,:,:,:,:,:,:),cgr(:,:,:,:),eqt(:),
     &  ppbrdx(:,:,:,:,:,:,:),aaa(:,:),symope(:,:), qibze(:,:),
     &  ppb(:),pdb(:),dpb(:),ddb(:), eq(:,:),
     &  eqx(:,:,:),eqx0(:,:,:),ekc(:),coh(:,:)
      complex(8),allocatable:: geigB(:,:,:,:) ,zsec(:,:,:)
c
      logical :: legas
      real(8) :: rydberg,hartree
      real(8):: qreal(3), ntot,nocctotg2,tripl,xxx(3,3)
      logical ::nocore

c space group infermation
      integer(4),allocatable :: iclasst(:), invgx(:), miat(:,:)
      real(8),allocatable    :: tiat(:,:,:),shtvg(:,:)

c tetra
      real(8),allocatable :: qz(:,:),qbzxx(:),wbzxx(:),wtet(:,:,:,:),
     &     eband(:,:,:), ene(:) !,ecore(:,:)
      integer(4),allocatable ::idtetx(:,:),idtet(:,:),ipq(:)
     &      ,iene(:,:,:),ibzx(:) !,nstar(:)
c      real(8) :: qbasmc(3,3)

ctm
      complex(8),allocatable :: zv(:,:),zw(:,:)
      real(8),allocatable :: rv(:,:),cv(:,:),rpv(:,:),cpv(:,:),
     &                       rw(:,:),cw(:,:),
     &                       wpv(:,:,:),ovlp(:,:,:),wdiag(:),
     &                       evec(:,:,:),eval(:),eval2(:),ecqw(:,:)
      integer(4),allocatable :: iwdiag(:)



c worksize in megabytes (1 word = 4 bytes)
c      parameter (mbytes=60)
c      parameter (mwords=mbytes/4)
c      parameter (iwksize=mwords * 1000 *1000)
c      integer w
c      common /w/ w(iwksize)


      integer(4) ::ib,iqx,igp,isx
c
      real(8),allocatable   :: eex1(:,:,:),exsp1(:,:,:),qqex1(:,:,:,:)
      integer(4),allocatable:: nspex(:,:),ieord(:),itex1(:,:,:)
      real(8)    :: qqex(1:3), eex,exsp,eee, exwgt
      integer(4) :: itmx,ipex,itpex,itex,nspexmx,nnex,isig,iex,ifexspx
     & ,ifexspxx ,ifefsm
      character(12) :: filenameex
      logical :: exspwrite=.false.

ctm 010705
      real(8) :: eclda_bh,eclda_pz

      real(8),allocatable :: freq_r(:),freq_i(:),freqx(:),wx(:),expa(:)
     &                      ,frhis(:)
      integer(4)::nwin, incwfin,  verbose,nqibze,ngc,iqxini !,bzcase
      real(8)::efin, qq(3),quu(3)
      integer(4),allocatable:: ngvecc(:,:)


      integer(4):: iqindx,iqbz, irecwq
      integer(4),allocatable:: nstibz(:)
      real(8):: erpaq, trpvq, trlogq !accumlating variable
      real(8):: erpaqw, trpvqw, trlogqw
c---------------------------------------

c set up work array
c      call wkinit (iwksize)
      call pshprt(60)


c---  readin BZDATA. See gwsrc/rwbzdata.f
c--------readin data set when you call read_BZDATA ---------------
c       integer(4)::ngrp,nqbz,nqibz,nqbzw,nteti,ntetf,
c     &   n_index_qbz
c       integer(4):: n1,n2,n3
c       real(8):: qbas(3,3),ginv(3,3),qbasmc(3,3)
c       real(8),allocatable:: qbz(:,:),wbz(:),qibz(:,:)
c     &    ,wibz(:),qbzw(:,:)
c       integer(4),allocatable:: idtetf(:,:),ib1bz(:),idteti(:,:)
c     &    ,nstar(:),irk(:,:),nstbz(:),  !index_qbz(:,:,:)
c-----------------------------------------------------------------
      call read_BZDATA()
      print *,' nqbz qbz =',nqbz
c      print *,  qbz
      print *,' nqibz ngrp=',nqibz,ngrp
c      print *,' BZmesh=',bzcase()
c      print *,' irk=',irk
c      print *,' #### idtetf: ####'
c      print *, idtetf
c$$$      if(bzcase()==2) then
c$$$        allocate(nstibz(nqibz))
c$$$c        do iq=1,nqbz
c$$$c          write(6,"(' iq qbz nstibz=',i5,3f9.4,i5)")iq,qbz(:,iq),nstbz(iq)
c$$$c        enddo
c$$$        do iq=1,nqibz
c$$$          iqbz = iqindx(qibz(:,iq),ginv,qbz,nqbz)
c$$$          nstibz(iq) = nstbz(iqbz)
c$$$          write(6,"(' iq qibz nstibz=',i5,3f9.4,i5)")iq,qibz(:,iq),nstibz(iq)
c$$$        enddo
c$$$c        stop
c$$$      endif
c-------------------------------
c generate all quantities needed
c-------------------------------
      nwin = -999d0   !Not readin nw
      efin = -999d0  !not readin EFERMI
      incwfin= 0  !use 7th colmn for core at the end section of GWIN
c
      call genallcf_v3(nwin,efin,incwfin) !in module m_genallcf_v3
Cstop2rx 2013.08.09 kino      if(ngrp/= ngrp2) stop 'ngrp inconsistent: BZDATA and LMTO GWIN_V2'
      if(ngrp/= ngrp2) call rx( 'ngrp inconsistent: BZDATA and LMTO GWIN_V2')
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
c-----------------------------------------------------------------------
c$$$      call genallcf_v2(
c$$$c> structure
c$$$     o                   plat,alat,natom,nclass,ipos,
c$$$c> symmetry
c$$$     o                   symgrp,gen,ngnmax,ngrp,igrp,
c$$$c> BZ
c$$$     o                   n1,n2,n3,qbas,ginv,iqibz,iwibz,nqibz,
c$$$     o                   iqbz,iwbz,nqbz,iindxk,
c$$$     o                   iinvg,instar,iirk,ef,
c$$$c>> file units
c$$$c     o                   ifphiv,ifphic,ifec,
c$$$c>> l,n and dimensions
c$$$     o                   clabl,nspin,nl,nn,nnv,nnc,
c$$$     o                   inindx,inindxv,inindxc,iiclass,
c$$$     d                   nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc,
c$$$c z value
c$$$     o                   iz,
c$$$c>> l,n,m indices for Phi
c$$$     o                   iil,iin,iim,   iilnm  ,i_mnl,
c$$$     o                   iilv,iinv,iimv,iilnmv ,i_mnlv,
c$$$     o                   iilc,iinc,iimc,iilnmc ,i_mnlc,
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
c      if (nclass > mxclass) stop ' hecor: increase mxclass'
c!!!! WE ASSUME iclass(iatom)= iatom !!!!!!!!!!!!!!!!!!!!!!!!!
c      if (nclass /= natom ) stop ' hecor: nclass /= natom ' ! We assume nclass = natom.
      print *,' hecor: end of genallc_v3'
c
      call pshprt(30)
      pi   = 4d0*datan(1d0)
      tpia = 2d0*pi/alat

c      shtw = 0d0
c      if(esmr<1d-5) shtw=0.01d0 ! Ferdi's shift to avoid resonance effect(maybe)
c
      call dinv33(plat,1,xxx,vol)
      voltot = dabs(vol)*(alat**3)

c--- ef is taken as rs for the empty-sphere test case of legas=T case ----
      legas = .false.
      INQUIRE (FILE = 'LEGAS', EXIST = legas)
      if(legas) then !!! test for electron gas case.
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
Cstop2rx 2013.08.09 kino        if(tetra) stop 'legas You have to give ef of  tetrahedron'
        if(tetra) call rx( 'legas You have to give ef of  tetrahedron')
      endif

c---
      write(6, *) ' --- computational conditions --- '
      write(6,'("    ua      =",f13.6)') ua
      write(6,'("    esmr    =",f13.6)') esmr
      write(6,'("    alat voltot =",2f13.6)') alat, voltot
c      write(6,'("    niw nw dw   =",2i3,f13.6)') niw,nw,dw

c>> read dimensions of wc,b,hb
      ifhbed     = iopen('hbe.d',1,0,0)
      read (ifhbed,*) nprecb,mrecb,mrece,nlmtot,nqbzt, nband
Cstop2rx 2013.08.09 kino      if (nprecb == 4) stop 'hecor: b,hb in single precision'
      if (nprecb == 4) call rx( 'hecor: b,hb in single precision')

      call init_readeigen(ginv,nspin,nband,mrece) !initialization of readEigen

c --- get space group information ---------------------------------
c true class information in order to determine the space group -----------
c because the class in the generated GW file is dummy.(iclass(ibas)=ibas should be kept).
      open (102,file='CLASS')
      allocate(iclasst(natom),invgx(ngrp)
     &         ,miat(natom,ngrp),tiat(3,natom,ngrp),shtvg(3,ngrp))
      print *,'  --- Readingin CLASS info ---'
      do ibas = 1,natom
        read(102,*) ibasx, iclasst(ibas)
        write(6, "(2i10)") ibasx, iclasst(ibas)
      enddo

c Get space-group transformation information. See header of mptaouof.
      call mptauof(symgg,ngrp,plat,natom,pos,iclasst
     o              ,miat,tiat,invgx,shtvg )
c        write (*,*)  'tiat=', tiat(1:3,1:natom,invr),invr

c Get array size to call rdpp
      call getsrdpp2( nclass,nl,nxx)
      call readngmx('QGpsi',ngpmx)
      call readngmx('QGcou',ngcmx)
      print *,' ngcmx ngpmx=',ngcmx,ngpmx
c      call getsrdpp( nclass,nl,
c     o               ngpmx,ngcmx,nxx )


      allocate( nx(0:2*(nl-1),nclass),
     &   nblocha(nclass) ,lx(nclass),
cv2
     &   ppbrd ( 0:nl-1, nn, 0:nl-1,nn, 0:2*(nl-1),nxx, nspin*nclass),
c     &   ppbrd ( 0:nl-1, nn, 0:nl-1,nn, 0:2*(nl-1),nxx, 4*nclass*nspin),
     &   cgr(nl**2,nl**2,(2*nl-1)**2,1),
     &          symope(3,3))

c     &   geigB(ngpmx,nband,nqbz,nspin),ngpn(nqbz),ngvecpB(3,ngpmx,nqbz),
c     &   ngcni(nqibz),ngveccB(3,ngcmx,nqibz), ngveccBr(3,ngcmx,nqibz), ! in IBZ !
c      call dcopy(ngrp*9,symgg,1,symope,1)
c- readin plane wave parts, and Radial integrals ppbrd.
! ppbrd = radial integrals
! cgr   = rotated cg coeffecients.
! geigB = eigenfunction's coefficiens for planewave.
! ngvecpB (in 1stBZ) contains G vector for eigen function.
! ngveccB (in IBZ)   contains G vector for Coulomb matrix.

c      call rdpp_v2( ngpmx,ngcmx,nxx,  qibz,nqibz, qbz,nqbz,
c     i      nband, nl,ngrp, nn,  nclass, nspin, symgg,    qbas,
c     o      nblocha, lx, nx, ppbrd ,
c     o      mdimx, nbloch, cgr,
c     o      nblochpmx, ngpn,geigB,ngvecpB,  ngcni,ngveccB )
ccc these are dummy.
      symope(1:3,1) = (/1d0,0d0,0d0/)
      symope(1:3,2) = (/0d0,1d0,0d0/)
      symope(1:3,3) = (/0d0,0d0,1d0/)
      ngrp=1
ccc
      call rdpp_v3(nxx, nl, ngrp, nn, nclass, nspin,symope,qbas,
     o      nblocha, lx, nx,  ppbrd , mdimx, nbloch, cgr)
      nblochpmx = nbloch + ngcmx

c dummy
      allocate(ngvecc(3,ngcmx))

      print *, ' end of rdpp_v3'

c----------------------------------------------
      call pshprt(60)

c------------
c open files
c------------
c> input files
c      ifqpnt    = iopen('QPNT',1,0,0)

c direct access files WVR and WVI which include W-V.
      ifwd      = iopen('WV.d',1,-1,0)
      read (ifwd,*) nprecx,mrecl,nblochpmx2,nwt,niwt, nqnum
      write(6,*) ' WV.d =', nprecx,mrecl,nblochpmx2,nwt,niwt, nqnum
      call checkeq(nprecx,ndble)
      call checkeq(nblochpmx,nblochpmx2)
c     if (nwt /= nw)   stop 'hecor: wrong nw'
Cstop2rx 2013.08.09 kino      if (niwt /= niw) stop 'hecor: wrong niw'
      if (niwt /= niw) call rx( 'hecor: wrong niw')
c     ifrcw     = iopen('WVR',0,-1,mrecl)
      ifrcwi    = iopen('WVI',0,-1,mrecl)

c eigen functions
c      ifev(1)   = iopen('EVU', 0,0,mrece)
c      if (nspin == 2) ifev(2) = iopen('EVD', 0,0,mrece)

c --- determine Fermi energy ef for given valn (legas case) or corresponding charge given by z and konf.
! When esmr is negative, esmr is geven automatically by efsimplef.
c        call efsimplef(ifev,nspin,wibz,w(iindxk),qibz
c     i       ,n1,n2,n3,qbas,ginv, nband,nqibz
c     i       ,w(ikonf),z,nl,natom,w(iiclass),nclass
c     i       ,valn, legas, esmr  !!! valn is input for legas=T, output otherwise.
c     o       ,efnew)
      call efsimplef2a(nspin,wibz,qibz,ginv, 
     i        nband,nqibz
     i       ,konf,z,nl,natom,iclass,nclass
     i       ,valn, legas, esmr,  !!! valn is input for legas=T, output otherwise.
c
     i       qbz,nqbz ! index_qbz, n_index_qbz,
     o       ,efnew)
c
      ef = efnew
c- check total ele number -------
      ntot  = nocctotg2(nspin, ef,esmr, qbz,wbz, nband,nqbz)
      print *,' ef    =',ef
      print *,' esmr  =',esmr
      print *,' valn  =',valn
      print *,' ntot  =',ntot

c        ifefsm  = iopen('EFERMI.hecor',1,-1,0)
c        write(ifefsm,*) ef,'!ef by smearing. by hecor'
c        ifefsm  = iclose('EFERMI.hecor')
c
c        if(abs(valn-ntot)>1d-6) stop ' abs(valn-ntot)>1d-6'  !20001 May

c--

      ifvcfpout = iopen('VCCFP',0,-1,0)
      read(ifvcfpout) ndummy1, ndummy2

c read q-points and states
c      lqall      = .true.
c      laf        = .false.
c      call readx   (ifqpnt,10)
c      read (ifqpnt,*) iqall,iaf
c      if (iaf   == 1)   laf = .true.
c      call readx   (ifqpnt,100)
! q-points
      nq         = nqibz
      allocate(q(3,nq))
      call dcopy   (3*nqibz,qibz,1,q,1)
c
c      nspinmx = nspin
c      if (laf) nspinmx =1

c      call winfo(6,nspin,nq,ntq,is,nbloch
c     &    ,ngpn(1),ngcni(1),nqbz,nqibz,ef,deltaw,alat,esmr)

c q near zero
      print *, 'reading QOP'
      open (101,file='Q0P')
      read (101,"(i5)") nq0i
      call checkeq(nqibz+nq0i-1, nqnum)
      write(6,*) ' *** nqibz nq0i=', nqibz,nq0i
      allocate( wqt(1:nq0i),q0i(1:3,1:nq0i) )
      do i=1,nq0i
        read (101,* ) wqt(i),q0i(1:3,i)
      enddo
      write(6,"(i3,f14.6,2x, 3f14.6)" )(i, wqt(i),q0i(1:3,i),i=1,nq0i)
      close(101)

c --- qibze(3,nqbze)
      nqibze = nqibz + nq0i
      allocate( qibze(3, nqibze))
      call dcopy(3*nqibz,qibz, 1, qibze,1)
      do i = 1,nq0i
        qibze(:,nqibz+i)  = q0i(:,i)
      enddo

c generate gaussian frequencies x between (0,1) and w=(1-x)/x
c      call defdr   (ifreqx,niw)
c      call defdr   (ifreqw,niw)
c      call defdr   (iwx,niw)
c      call defdr   (iexpa,niw)
      allocate( freq_i(niw) ,freqx(niw),wx(niw),expa(niw) )
      call freq01  (niw,ua,
     o              freqx,freq_i,wx,expa)

c-----------------------------------------------------------
c do loop over q { IBZ
      nqitot  = nqibz + nq0i
      totexc  = 0.d0
      trpv    = 0.d0
      trlog   = 0.d0
      allocate ( ecqw(nqitot,niw) )


c      if(bzcase()==1) then
c        iqxini = 2
c      else
c        iqxini = 1
c      endif
      iqxini=1

      allocate( zw(nblochpmx2, nblochpmx2) )
      do iq = iqxini,nqitot ! q=(0,0,0) is omitted!
c        if(iq<=nqibz) then
c          iqx = iq
c        else
c          iqx = 1                !corresponding q=0 !BUG fix at July 20th
c        endif
c        q = qibz(:,iq)
        if(iq<=nqibz) then
          qq= qibze(:,iq)
        else
          qq=0d0
        endif
        write(6,"('iq qq =',3f9.4,i5,' out of',i5)") qq, iq,nqibz+nq0i

        call readqg('QGcou', qq, ginv, quu,ngc,ngvecc)
        ngb = nbloch + ngc

c allocate
c        allocate( zw(nblochpmx2,nblochpmx2), zv(ngb,ngb),
c     &            rv(ngb,ngb), cv(ngb,ngb),
c     &            rpv(ngb,ngb), cpv(ngb,ngb),
c     &            rw(ngb,ngb), cw(ngb,ngb),
c     &            wpv(ngb,ngb,2), ovlp(ngb,ngb,2),
c     &            wdiag(11*ngb), iwdiag(ngb),
c     &            evec(ngb,ngb,2), eval(ngb), eval2(ngb) )

        allocate(zv(ngb,ngb)) !, rv(ngb,ngb), cv(ngb,ngb))

c Coulomb matrix
        read(ifvcfpout) zv(1:ngb,1:ngb)

c integration over the imaginary axis
c        erpaq      = 0d0
c        trpvq      = 0d0
c        trlogq     = 0d0
        do ix=1,niw
c          if(bzcase()==2) then
c            irecwq = (iq-1)*niw + ix
c          else
            irecwq = (iq-2)*niw + ix
c          endif
          read(ifrcwi, rec=irecwq) zw
          call ecorq2 (iq,ix,
     i     zv, zw(1:ngb,1:ngb),
     i     wibz,wqt,wx,freqx,
     d     nqibz,nq0i,ngb,niw, nstibz,nqbz,!takao
     o     erpaqw, trpvqw, trlogqw,
     o     ecqw(iq,ix))

c contribution for q
c         erpaq  = erpaq  + erpaqw
c         trpvq  = trpvq  + trpvqw
c         trlogq = trlogq + trlogqw
c total
          totexc  = totexc+ erpaqw
          trpv    = trpv  + trpvqw
          trlog   = trlog + trlogqw
        enddo

c
c        call ecorq (iq,
c     i              ifrcwi,zw,rv,cv,
c     i              wibz,wqt,wx,freqx,
c     w              rpv,cpv,rw,cw,
c     w              wpv,ovlp,wdiag,iwdiag,evec,eval,eval2,
c     d              nqibz,nq0i,nqitot,ngb,niw,nblochpmx2,   nstibz,nqbz,!takao
c     o totexc,trpv,trlog,ecqw)

c        deallocate(zw,rv,cv,rpv,cpv,rw,cw,
c     &             wpv,ovlp,wdiag,iwdiag,evec,eval,eval2)
        deallocate(zv)
c< end of q-loop
      end do

c---------------------------------

c------------
c output
c------------
      ifcor   = iopen('TEECORR',1,-1,0) ! output files
      write(ifcor,*) '============================'
      write(ifcor,*) 'Correlation energy Erpa (eV)'
      write(ifcor,*) '============================'

c        call winfo(ifsec(is),nspin,nq,ntq,is,nbloch
c     &    ,ngpn(1),ngcni(1),nqbz,nqibz,ef,deltaw,alat,esmr)

      write (ifcor,*)' *** '
      hartree=2d0*rydberg()
      write(ifcor,*)totexc*hartree,trpv*hartree,trlog*hartree

c---------------------------------
c electron gas correlation energy
      if (legas) then
        pi         = 4.d0*datan(1.d0)
        efz=(ntot*3*pi**2/voltot)**(2d0/3d0) ! ef is calculated from ntot.
        qfermi= dsqrt(efz)
        alpha = (9*pi/4d0)**(1d0/3d0)
        rs    = alpha/qfermi
        write (ifcor,*)' --- electron gas ---'
        write (ifcor,*)' density parameter rs= ', rs
        write (ifcor,*)' kf= ',qfermi
        write (ifcor,*)' *** Barth-Hedin formula'
        ecelgas = eclda_bh(rs) * hartree * ntot
        write (ifcor,*)ecelgas
        write (ifcor,*)' *** Perdew-Zunger formula'
        ecelgas = eclda_pz(rs) * hartree * ntot
        write (ifcor,*)ecelgas
        write (ifcor,*)' *** Gell-Mann and Brueckner formula'
        ecelgas = (-0.0311d0 * dlog(rs) -0.048d0) * hartree * ntot
        write (ifcor,*)ecelgas
      endif

c---------------------------------
c output ecqw
      call wecqw(ifcor,
     d           nqibz,nqbz,nq0i,nqitot,niw,
     o           wibz,wqt,wx,freqx,ecqw)

c------------
c close files
c------------
      isx = iclose ('wc.d')
      isx = iclose ('wci.d')
      isx = iclose ('hbe.d')

      call cputid(0)
Cstop2rx 2013.08.09 kino      stop  ' OK! hecor'
      call rx0( ' OK! hecor')
      end
c--------------------------------------------------------------------
      subroutine winfo(ifi,nspin,nq,ntq,is,nbloch,ngp,ngc,nqbz,nqibz,ef
     &    ,deltaw,alat,esmr)
      implicit none
      integer(4) :: ifi,nspin,nq,ntq,is,nbloch,ngp,ngc,nqbz,nqibz
      real(8) :: ef,deltaw,alat,esmr
      write (ifi,*)' ***'
      write (ifi,6700) nspin,nq,ntq
      write (ifi,6501) is,nbloch,ngp,ngc,nqbz,nqibz,ef
     &  ,deltaw,alat,ef,esmr
 6501 format (' spin =',i2,'   nbloch ngp ngc=',3i4
     &        ,'  nqbz =',i6,'  nqibz =',
     .        i6,'   ef=', f10.4,' Rydberg'
     &        ,/,d23.16,' <= deltaw(Hartree)'
     &        ,/,d23.16,' <= alat'
     &        ,/,d23.16,' <= ef '
     &        ,/,d23.16,' <= esmr')
 6700 format (1x,3i4,'  nspin  nq  ntq')
      end
c------------------------------------------------------------------
      subroutine checkeq(i,j)
Cstop2rx 2013.08.09 kino      if(i/=j) stop 
      if(i/=j) call rx0( ''//
     & " checkeq in hsfp0: dim of WVR and WVI not compatible")
      end
c--------------------------------------------------------------------
      real*8 function eclda_bh(rs)
      real(8) :: rs,cp,rp,z

      cp       = 0.0504d0*0.5d0 ! 0.5 changes unit from Ry to Hartree
      rp       = 30.d0
      z        = rs / rp
      eclda_bh = -cp * ( (1.d0+z**3)*dlog(1.d0+1.d0/z) 
     .                 + 0.5d0*z - z**2 - 0.33333333d0 )
      end
c--------------------------------------------------------------------
      real*8 function eclda_pz(rs)
      real(8) :: rs
      if (rs.ge.1.d0) then
        eclda_pz = -0.1423d0 / (1.d0 + 1.0529d0*dsqrt(rs) + 0.334d0*rs)
      else
        eclda_pz = -0.0480d0 + 0.0311d0*dlog(rs) - 0.0116d0 * rs 
     .         + 0.0020d0*rs*dlog(rs)
      endif
      end
c--------------------------------------------------------------------
      subroutine wecqw(ifcor,
     d                 nqibz,nqbz,nq0i,nqitot,niw,
     o                 wibz,wqt,wx,freqx,ecqw)

      implicit double precision (a-h,o-z)

      dimension   wibz(nqibz),wqt(nq0i),wx(niw),
     .            freqx(niw),ecqw(nqitot,niw)
      real(8):: rydberg
      write(ifcor,*)'*** ecqw(q,w) ***'
      write(ifcor,*)'nqibz =',nqibz
      write(ifcor,*)'nq0i  =',nq0i
      write(ifcor,*)'niw   =',niw

      do ip = 2,nqitot
        if (ip <= nqibz) then
          wk = wibz(ip)*0.5d0 ! 0.5 for the normalization of wibz
        else
c        wk = wqt(ip-nqibz)*wibz(1)*0.5d0 ! 0.5 for the normalization of wibz
          wk = wqt(ip-nqibz)* 1d0/dble(nqbz)
        endif
        write(ifcor,*)'*** iq,wq = ',ip,wk

        sume=0d0
        do ix = 1,niw
          write(ifcor,*)freqx(ix),ecqw(ip,ix),wx(ix)
          sume=sume+  wx(ix)/(freqx(ix)*freqx(ix)) * ecqw(ip,ix)
        end do
        write(ifcor,*) '  sum ecqw*wx=', wk*sume*2d0*rydberg()

c end of ip-loop
      end do

      return
      end
c---------------------------------------------------------------------
