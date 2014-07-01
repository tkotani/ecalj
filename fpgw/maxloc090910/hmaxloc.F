      program hmaxloc
c-------------------------------------------------------------
c construct maximally localized Wannier functions
c
c References
c [1] N. Marzari and D.Vanderbilt, PRB56,12847(1997)
c [2] I. Souza, N. Marzari and D.Vanderbilt, PRB65,035109(2002)
c
c mode 1:  determine parameters for <u(m,k)|u(n,k+b> (uu-matrix)
c mode 2:  main part
c    Step 1: choose Hilbert space  (Ref.[2])
c    Step 2: maximally localize Wannier functions (Ref.[1])
c    Step 3: construct effective Hamiltonian and interpolate bands (Ref.[2])
c
cm Oct 2008 Takashi Miyake, updated
cm Aug 2007 Takashi Miyake, berry connection in the Wannier gauge
c  May 2004 Takashi Miyake, from hwmat.f
c------------------------------------------------------------
      use m_readqg
      use m_readeigen
      use m_read_bzdata,ngrp2=>ngrp
      use m_genallcf_v3
      use keyvalue
      implicit none
c------------------------------------
      real(8)    :: esmr2,shtw
      integer(4)::
     &   ixc,iopen,ifhbed, nprecb,mrecb,mrece,nlmtot,nqbzt, nband,
     &   ibas,ibasx,ngpmx,nxx,ngcmx,nbloch,ifqpnt,ifwd,ifbb,
     &   nprecx,mrecl,nblochpmx2,nwt,niwt, nqnum,mdimx,nblochpmx,
     &   ifrcw,ifrcwi,  noccxv,maxocc2,noccx,ifvcfpout,iqall,iaf,ntq,
     &   i,j,k,nspinmx, nq,is,ip,iq,idxk,ifoutsex,iclose,nq0i,ig,
     &   mxkp,nqibzxx,ntet,nene,iqi, ix,iw,
     &   nlnx4,niwx,irot,invr,invrot,ivsum, ifoutsec,ntqx,
     &    ifmlw(2),ifmlwe(2) !,ifcphi
     &   ,ifxc(2),ifsex(2), ifphiv(2),ifphic(2),ifec,ifexsp(2),
     &   ifsecomg(2),ifexx,ifwand,ndble=8
      real(8) :: pi,tpia,vol,voltot,rs,alpha,
     & qfermi,efx,valn,efnew,edummy,efz,qm,xsex,egex,
     & zfac1,zfac2,dscdw1,dscdw2,dscdw,zfac,ef2=1d99,exx,exxq,exxelgas
      logical lqall,laf

      integer(4),allocatable :: itq(:)
      real(8),allocatable    :: q(:,:)

c takao
      integer(4),allocatable :: ngvecpB(:,:,:),!ngveccB(:,:,:),
     &  ngvecp(:,:), ngvecc(:,:),iqib(:), !,ngpn(:)ngcni(:)
     &  kount(:,:), nx(:,:),nblocha(:),lx(:) !ngveccBr(:,:,:)
      real(8),allocatable:: vxcfp(:,:,:),
     &  wqt(:), wgt0(:,:),q0i(:,:),
     &  ppbrd (:,:,:,:,:,:,:),cgr(:,:,:,:),eqt(:),
     &  ppbrdx(:,:,:,:,:,:,:),aaa(:,:), !symope(:,:,:)=symgg, ! qibz(:,:),
     &  ppb(:), eq(:), !,pdb(:),dpb(:),ddb(:)
     &  eqx(:,:,:),eqx0(:,:,:),ekc(:),coh(:,:)
     &         , rw_w(:,:,:,:,:),cw_w(:,:,:,:,:),
     &           rw_iw(:,:,:,:,:),cw_iw(:,:,:,:,:)
      complex(8),allocatable:: geigB(:,:,:,:)
c 
      logical :: screen, exchange, cohtest, legas, tote 
      real(8) ::  rydberg,hartree 
      real(8):: qreal(3), ntot,nocctotg2,tripl,xxx(3,3)
      logical ::nocore

c space group infermation 
      integer(4),allocatable :: iclasst(:), invgx(:), miat(:,:)
      real(8),allocatable    :: tiat(:,:,:),shtvg(:,:)

c
      real(8),allocatable   :: eex1(:,:,:),exsp1(:,:,:),qqex1(:,:,:,:)
      integer(4),allocatable:: nspex(:,:),ieord(:),itex1(:,:,:)
      real(8)    :: qqex(1:3), eex,exsp,eee, exwgt,deltax0
      integer(4) :: itmx,ipex,itpex,itex,nspexmx,nnex,isig,iex,ifexspx
     & ,ifexspxx ,ifefsm, nq0ix,ifemesh,nz
      character(3)  :: charnum3,sss
      character(12) :: filenameex
      logical :: exspwrite=.false.
      character*8 xt


      integer(4)::nqbze,ini,nq0it,idummy
      real(8),allocatable:: qbze(:,:)

      real(8)   :: ebmx
      integer(4):: nbmx
      
      real(8):: volwgt

      integer(4)::nwin, incwfin
      real(8)::efin,ddw
      integer(4),allocatable::imdim(:)
      real(8),allocatable::freqx(:),freqw(:),wwx(:),expa(:)

      logical:: GaussSmear !readgwinput, 
      integer(4)::ret
      character*(150):: ddd


      integer(4):: bzcase,  ngpn1,mrecg,verbose,ngcn1,nwxx
      real(8)   :: wgtq0p,quu(3)

      integer(4):: iii,isx,ivsumxxx

c for maxloc
      real(8)   :: wbb(12),wbbsum,bb(3,12),
     c             eomin,eomax,eimin,eimax,
     c             qwf0(3),dqwf0(3),qks(3),q0(3)
      complex(8),allocatable:: uumat(:,:,:,:),evecc(:,:),eveccs(:,:),
     c                         amnk(:,:,:),cnk(:,:,:),umnk(:,:,:)
      real(8),allocatable:: ku(:,:),kbu(:,:,:),eunk(:,:),eval(:),evals(:),
     c                      eks(:),rt(:,:),rt8(:,:,:),qbz0(:,:),r0g(:,:),
     c                      wphi(:,:)
      integer(4):: nbb,isc,nwf,ifmloc,ifq0p,
     c             nox,iko_ix,iko_fx,
     c             noxs(2),iko_ixs(2),iko_fxs(2),
     c             ieo_swt,iei_swt,itin_i,itin_f,itout_i,itout_f,nphix,
     c             nbbelow,nbabove
      integer(4),allocatable:: ikbidx(:,:)
      integer(4),allocatable:: iki_i(:),iki_f(:),
     c                         ikbi_i(:,:),ikbi_f(:,:),
     c                         iko_i(:),iko_f(:),
     c                         ikbo_i(:,:),ikbo_f(:,:),
     c                         iphi(:,:),iphidot(:,:),
     c                         nphi(:)
      logical :: leout,lein,lbin,lq0p,lsyml,lbnds
      logical :: debug=.false.
!
      integer(4):: nlinex,ntmp
      parameter (nlinex=100)
      integer(4)::nline,np(nlinex)
      real(8):: qi(3,nlinex),qf(3,nlinex)
c step 1
      complex(8),allocatable:: cnq0(:,:),
     c                         upu(:,:,:,:),cnk2(:,:,:),
     c                         zmn(:,:)
      complex(8):: ctmp
      real(8),allocatable:: omgik(:)
      real(8)   :: omgi,omgiold,conv1,alpha1,domgi,qtmp(3)
      integer(4):: nsc1,ndz,nin,ifhoev,ifuu0,ifpsig
c step 2
      complex(8),allocatable:: mmn(:,:,:,:),mmn0(:,:,:,:),
     c                         rmn(:,:),smn(:,:),amn(:,:),
     c                         tmn(:,:),dwmn(:,:)
      real(8),allocatable:: rn(:,:),qn(:)
      real(8)   :: omgd,omgod,omgdod,omgidod,omgdodold,domgdod,
     c             conv2,alpha2
      integer(4):: nsc2,ibb,ii,ij,ik
      logical :: lrmn,lmmn
c step 3
      complex(8),allocatable:: hrotk(:,:,:),hrotr(:,:,:),hrotkp(:,:)
     c                       , hrotkps(:,:)
      real(8):: e1,e2,rcut
      integer(4):: iband,ifbnd,iftb,ifsh,nsh,nsh1,nsh2
      logical :: lsh
      real(8),allocatable :: rws(:,:),drws(:)
      integer(4),allocatable:: irws(:)
      integer(4):: nrws,ifham

c ixc=3
      character(20)::filename
      complex(8),allocatable::  hrotrcut(:,:,:)
      integer:: ifh
      real(8):: heps ,r_v
c---------------------------------------
      hartree=2d0*rydberg()

      iii=verbose()
      print *,' verbose=',iii

c mode switch. --------------
      write(6,*) ' --- Choose omodes below ----------------'
      write(6,*) '  bb vectors (1) or Wannier fn. (2) or TB Hamiltonian (3)'
      write(6,*) ' --- Put number above ! -----------------'
      call readin5(ixc,nz,idummy)
      write(6,*) ' ixc=',ixc
      if(ixc==0) stop ' --- ixc=0 --- Choose computational mode!'

c---  readin BZDATA. See gwsrc/rwbzdata.f
c--------readin data set when you call read_BZDATA ---------------
c       integer(4)::ngrp,nqbz,nqibz,nqbzw,nteti,ntetf
ccccc    ! &   ,n_index_qbz  
c       integer(4):: n1,n2,n3
c       real(8):: qbas(3,3),ginv(3,3),qbasmc(3,3)
c       real(8),allocatable:: qbz(:,:),wbz(:),qibz(:,:)
c     &    ,wibz(:),qbzw(:,:)
c       integer(4),allocatable:: idtetf(:,:),ib1bz(:),idteti(:,:)
c     &    ,nstar(:),irk(:,:),nstbz(:)          !,index_qbz(:,:,:)
c-----------------------------------------------------------------
      call read_BZDATA()
      print *,' nqbz  =',nqbz
c      print *,  qbz 
      print *,' nqibz ngrp=',nqibz,ngrp
c      print *,' irk=',irk
c      print *,' #### idtetf: ####'
c      print *, idtetf

c set up work array
c      call wkinit (iwksize)
      call pshprt(60)


C--- readin GWIN and LMTO, then allocate and set datas.
      nwin =-999    !not readin NW file
      efin =-999d0  !not readin EFERMI
c      efin = 0d0    !readin EFERMI
      incwfin= -1  !use 7th colmn for core at the end section of GWIN
      call genallcf_v3(nwin,efin,incwfin) !in module m_genallcf_v3
      if(ngrp/= ngrp2) stop 'ngrp inconsistent: BZDATA and LMTO GWIN_V2'
c---  These are allocated and setted.
c      integer(4)::  nclass,natom,nspin,nl,nn,nnv,nnc, ngrp,
c     o  nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, nctot,niw, !not readin nw
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
c     o  plat(:,:),pos(:,:),z(:),  ecore(:,:),  symgg(:,:,:) ! symgg=w(igrp),freq(:)
c-----------------------------------------------------------------------

c--- Get maximums takao 18June03
      call getnemx(nbmx,ebmx,8,.true.) !8+1 th line of GWIN0

c-------------------------------------------------------------------
c      if (nclass > mxclass) stop ' hsfp0: increase mxclass'
c!!!! WE ASSUME iclass(iatom)= iatom !!!!!!!!!!!!!!!!!!!!!!!!!
      if (nclass /= natom ) stop ' hsfp0: nclass /= natom ' ! We assume nclass = natom.
      print *,' hsfp0: end of genallcf2'
c
      call pshprt(30)
      pi   = 4d0*datan(1d0)
      tpia = 2d0*pi/alat

      call dinv33(plat,1,xxx,vol)
      voltot = dabs(vol)*(alat**3)

      ifmlw(1) = iopen('MLWU',0,-1,0)
      ifmlwe(1)= iopen('MLWEU',0,-1,0)
      if (nspin == 2) then 
         ifmlw(2) = iopen('MLWD',0,-1,0)
         ifmlwe(2)= iopen('MLWED',0,-1,0)
      endif

c>> read dimensions of wc,b,hb
      ifhbed     = iopen('hbe.d',1,0,0)
      read (ifhbed,*) nprecb,mrecb,mrece,nlmtot,nqbzt, nband,mrecg
      if (nprecb == 4) stop 'hsfp0: b,hb in single precision'

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
      
c----------------------------------------------
      call pshprt(60)

c... Readin eigen functions
c      ifev(1)   = iopen('EVU', 0,0,mrece)
c      if (nspin==2) ifev(2) = iopen('EVD', 0,0,mrece)

! read EF from 'BNDS' if exists
      lbnds=.false.
      inquire(file='BNDS',exist=lbnds)      
      if (lbnds) then
        write(*,*)'Read EF from BNDS'
        open(99,file='BNDS',status='old')
        read(99,*)ntmp,ef
        close(99)
      else ! lbnds   
c --- determine Fermi energy ef for given valn (legas case) or corresponding charge given by z and konf.
! When esmr is negative, esmr is geven automatically by efsimplef.
        write(*,*)'Calculate EF in efsimplef2a'
        legas = .false.
        call efsimplef2a(nspin,wibz,qibz,ginv, 
     i        nband,nqibz
     i       ,konf,z,nl,natom,iclass,nclass
     i       ,valn, legas, esmr,  !!! valn is input for legas=T, output otherwise. 
c
     i        qbz,nqbz !index_qbz, n_index_qbz,
     o       ,efnew)
c
c        print *,' end of efsimple'
        ef = efnew
      endif ! lbnds
c- check total ele number -------
      ntot  = nocctotg2(nspin, ef,esmr, qbz,wbz, nband,nqbz) !wbz
      print *,' ef    =',ef
      print *,' esmr  =',esmr
      print *,' valn  =',valn
      print *,' ntot  =',ntot

c      ifcphi  = iopen('CPHI',0,0,mrecb)

      call init_readeigen2(mrecb,nlmto,mrecg) !initialize m_readeigen

!c QPNT data
ctm, 080222
! read QPNT from 'SYML' if exists
      lsyml=.false.
      inquire(file='SYML',exist=lsyml)      
      if (lsyml) then
        write(*,*)'Read k points for bands from SYML'
        lqall      = .false.
        laf        = .false.
        open(99,file='SYML',status='old')
        do i = 1,nlinex
          read(99,*)np(i),qi(1:3,i),qf(1:3,i)
          if (np(i).eq.0) exit
          nline = i
        enddo
        if (nline.eq.nlinex) stop 'hmaxloc: too many lines in SYML'
        close(99)
        nq = 0
        do i = 1,nline
           nq = nq + np(i)
        enddo ! i
        allocate(q(3,nq))
        iq = 0
        do i = 1,nline
        do j = 0,np(i)-1
          iq = iq + 1
          q(:,iq) = qi(:,i) + (qf(:,i)-qi(:,i))*dble(j)/dble(np(i)-1)
        enddo ! j
        enddo ! i
      else ! lsyml
        write(*,*)'Read k points for bands from GWinput'
        call getkeyvalue("GWinput","<QPNT>",unit=ifqpnt,status=ret) 
        print *,' ifqpnt ret=',ifqpnt,ret
c
        lqall      = .false.
        laf        = .false.
        call readx   (ifqpnt,10)
        read (ifqpnt,*) iqall,iaf
        if (iqall == 1) lqall = .true.
        if (iaf   == 1)   laf = .true.
        call readx   (ifqpnt,100)
ctm 040622
        read (ifqpnt,*)
        read (ifqpnt,*)

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
        endif ! lqall
        close(ifqpnt)
      endif ! syml
c
      nspinmx = nspin
      if (laf) nspinmx =1
c------------
c input parameters specific to MAXLOC
      call getkeyvalue("GWinput","<MLWF>",unit=ifmloc,status=ret) 
      print *,' ifmloc ret=',ifmloc,ret
      read (ifmloc,*) nwf
      allocate (nphi(nwf))
      read (ifmloc,*) (nphi(i),i=1,nwf)
      nphix = 0
      do i = 1,nwf
         if(nphi(i).gt.nphix)nphix = nphi(i)
      enddo
      allocate (r0g(nphix,nwf),iphi(nphix,nwf),iphidot(nphix,nwf),
     &          wphi(nphix,nwf))
      do i=1,nwf
         do j=1,nphi(i)
            read(ifmloc,*)iphi(j,i),iphidot(j,i),r0g(j,i),wphi(j,i)
         enddo
      enddo
      close(ifmloc)

      call wan_input(leout,lein,lbin,ieo_swt,iei_swt,
     &    eomin,eomax,itout_i,itout_f,nbbelow,nbabove,
     &    eimin,eimax,itin_i,itin_f,
     &    nsc1,nsc2,conv1,conv2,alpha1,alpha2,rcut)
c

cskino
      r_v=rcut
      call getkeyvalue("GWinput",'wan_tbcut_rcut',heps,default=r_v) 
      call getkeyvalue("GWinput",'wan_tbcut_heps',heps,default=0.0d0) 
      write(*,*) 'mloc.heps ', heps 
cekino


cc --- read LDA eigenvalues
      ntq = nwf
cc     ntp0=ntq
c      allocate(eqx(ntq,nq,nspin),eqx0(ntq,nq,nspin),eqt(nband))
c      do      is = 1,nspin
c      do      ip = 1,nq
cc        iq       = idxk (q(1,ip),qbze,nqbze)
cc        call rwdd1   (ifev(is), iq, nband, eqt) !direct access read b,hb and e(q,t)
c        call readeval(q(1,ip),is,eqt)
cc        print *,' eqt=',eqt
c        eqx0(1:ntq,ip,is) = eqt(itq(1:ntq))
c        eqx (1:ntq,ip,is) = rydberg()*(eqt(itq(1:ntq))- ef)
c      enddo
c      enddo
c      deallocate(eqt)

c --- info
      call winfo(6,nspin,nq,ntq,is,nbloch
     &    ,0,0,nqbz,nqibz,ef,deltaw,alat,esmr)

c
      iii=ivsumxxx(irk,nqibz*ngrp)
      write(6,*) " sum of nonzero iirk=",iii, nqbz

c-------------------------------------------------------------
c debug:
c      allocate(eqt(nband))
c      do ip = 1,nqbz
c        call readeval(qbz(1,ip),1,eqt)         
c        write(80,"('***',3f10.5)")qbz(:,ip)
c        do is=1,nband
c           write(80,"(i5,f12.6)")is,eqt(is)
c        enddo
c      enddo
   
c Rt vectors
      allocate (rt(3,nqbz),rt8(3,8,nqbz),qbz0(3,nqbz))
      call getrt(qbz,qbas,plat,n1,n2,n3,nqbz,
     o           rt,rt8,qbz0)

c b vectors
      call getbb(plat,alat,n1,n2,n3,
     o           nbb,wbb,wbbsum,bb)

c index for k and k+bb
      allocate (ku(3,nqbz),kbu(3,nbb,nqbz),ikbidx(nbb,nqbz))

      call kbbindx(qbz,ginv,bb,
     d             nqbz,nbb,
     o             ikbidx,ku,kbu)


      allocate (iko_i(nqbz),iko_f(nqbz),
     &          iki_i(nqbz),iki_f(nqbz),
     &          ikbo_i(nbb,nqbz),ikbo_f(nbb,nqbz),
     &          ikbi_i(nbb,nqbz),ikbi_f(nbb,nqbz))


c ixc = 1
      if (ixc.eq.1) then
         do is = 1,nspin
         call ewindow(is,ieo_swt,iei_swt,itout_i,itout_f,itin_i,itin_f,
     i                eomin,eomax,eimin,eimax,ef,qbz,ikbidx,
     i                nbbelow,nbabove,
     d                nqbz,nbb,nband,nwf,nspin,
     o                iko_i,iko_f,iki_i,iki_f,
     o                ikbo_i,ikbo_f,ikbi_i,ikbi_f,
     o                iko_ixs(is),iko_fxs(is),noxs(is),
     o                leout,lein)
         enddo

c write bb vectors to 'BBVEC'
         call writebb(ifbb,wbb(1:nbb),bb(1:3,1:nbb),
     i                ikbidx,ku,kbu,
     i                iko_ixs,iko_fxs,noxs,
     d                nspin,nqbz,nbb)

ctm, 060923 !!!
         ifwand = iopen('wan.d',1,-1,0)
         iko_ix = iko_ixs(1)
         iko_fx = iko_fxs(1)
         if (nspin.eq.2) then
            if (iko_ixs(2).lt.iko_ix) iko_ix = iko_ixs(2)
            if (iko_fxs(2).gt.iko_fx) iko_fx = iko_fxs(2)
         endif
         write(ifwand,*)nqbz,nwf,iko_ix,iko_fx
         write(ifwand,*)nspin
         do is = 1,nspin
           write(ifwand,*)nqbz,nwf,iko_ixs(is),iko_fxs(is)
         enddo
         isx = iclose('wan.d')
         stop 'hmaxloc: ixc=1 ok'
      endif

c loop over spin
      do is = 1,nspin
         write(*,*)'is =',is,'  out of',nspin

c energy window
      call ewindow(is,ieo_swt,iei_swt,itout_i,itout_f,itin_i,itin_f,
     i             eomin,eomax,eimin,eimax,ef,qbz,ikbidx,
     i             nbbelow,nbabove,
     d             nqbz,nbb,nband,nwf,nspin,
     o             iko_i,iko_f,iki_i,iki_f,
     o             ikbo_i,ikbo_f,ikbi_i,ikbi_f,
     o             iko_ix,iko_fx,nox,
     o             leout,lein)
!      call chk_ewindow(ifbb,is,nspin,nqbz,nbb,iko_ix,iko_fx)

c read uu-matrix
      allocate (uumat(iko_ix:iko_fx,iko_ix:iko_fx,nbb,nqbz))
      call readuu(is,iko_ix,iko_fx,ikbidx,
     d            nqbz,nbb,
     o            uumat)
      call chkuu(is,iko_ix,iko_fx,ikbidx,uumat,
     d           nqbz,nbb)

c step 1  -- choose Hilbert space --
      write(*,*)'Step 1: Hilbert space branch'

      allocate (amnk(iko_ix:iko_fx,nwf,nqbz),
     &          upu(iko_ix:iko_fx,iko_ix:iko_fx,nbb,nqbz),
     &          cnk(iko_ix:iko_fx,nwf,nqbz),
     &          cnk2(iko_ix:iko_fx,nwf,nqbz),
     &          omgik(nqbz))

      call init_unkg(is,qbz,ginv,ef,lein,
     i              iko_ix,iko_fx,iko_i,iko_f,
     i              iki_i,iki_f,
     d              nwf,nband,nqbz,
     o              amnk,cnk)

!      call chk_amnkweight(qbz,iko_ix,iko_fx,amnk,
!     &     nqbz,nwf,nband,nlmto)
!      call chk_cnkweight(qbz,iko_ix,iko_fx,cnk,
!     &     nqbz,nwf,nband,nlmto)

      do isc = 1,nsc1

         do iq = 1,nqbz

            call dimz(lein,iko_i(iq),iko_f(iq),iki_i(iq),iki_f(iq),
     o                ndz,nin)

            if (nwf.gt.nin) then
              if (ndz.lt.1) stop 'ndz < 1'

c (1-2) <u_mk | P_k+b | u_nk>
               call getupu(isc,
     i                     uumat(:,:,:,iq),cnk,
     i                     lein,alpha1,iq,ikbidx(:,iq),
     i                     iko_ix,iko_fx,
     i                     iko_i(iq),iko_f(iq),
     i                     iki_i(iq),iki_f(iq),
     i                     ikbo_i(:,iq),ikbo_f(:,iq),
     i                     ikbi_i(:,iq),ikbi_f(:,iq),
     d                     nwf,nbb,nqbz,
     u                     upu(:,:,:,iq))

c (1-3) Zmn(k) > phi,eval
               allocate (zmn(ndz,ndz),evecc(ndz,ndz),eval(ndz))

               call getzmn(upu(:,:,:,iq),wbb,lein,
     i                     iko_ix,iko_fx,
     i                     iko_i(iq),iko_f(iq),
     i                     iki_i(iq),iki_f(iq),
     d                     nwf,nbb,nqbz,ndz,
     o                     zmn)

               call chk_hm(zmn,ndz)
               call diag_hm(zmn,ndz,eval,evecc)

               call new_cnk(cnk(:,:,iq),evecc,iq,
     i                      iko_ix,iko_fx,
     i                      iko_i(iq),iko_f(iq),
     i                      iki_i(iq),iki_f(iq),
     d                      nwf,ndz,
     o                      cnk2(:,:,iq))

c (1-3) w_I(k)  eq.(18)
               call chk_eval(wbb,eval,nbb,ndz)
               call get_omgik(wbb,eval,
     i                        iko_i(iq),iko_f(iq),
     i                        iki_i(iq),iki_f(iq),
     d                        nbb,nwf,ndz,
     o                        omgik(iq))

               deallocate (zmn,evecc,eval)

            else
               omgik(iq) = 0d0
               cnk2(:,:,iq) = cnk(:,:,iq)

c end if (ndz>1)
            endif

c end of iq-loop
         enddo

c (1-5) w_I(k) > Omaga_I  eq.(11)
         omgi = sum(omgik(:)*wbz(:))

c (1-6) check self-consistency
         write(*,*)'#SC-loop, conv.',isc,omgi 
         if (isc.ge.2) then
            domgi = dabs((omgiold - omgi) / omgiold)
            if (domgi .lt. conv1) then
               write(*,*) 'step1: converged!'
               goto 810
            endif
         endif

c update
         omgiold = omgi
         cnk     = cnk2

c end of self-consistent loop
      enddo
      write(*,*)'step1: not converged'

 810  continue

      deallocate(upu,cnk2)

!      call chk_cnkweight(qbz,iko_ix,iko_fx,cnk,
!     &     nqbz,nwf,nband,nlmto)

c step 2 -- localize Wannier fn. --
      write(*,*)'Step 2: Wannier fn. branch'

      allocate (mmn(nwf,nwf,nbb,nqbz),mmn0(nwf,nwf,nbb,nqbz),
     &          umnk(nwf,nwf,nqbz),
     &          rmn(nwf,nwf),amn(nwf,nwf),smn(nwf,nwf),
     &          rn(3,nwf),qn(nwf),tmn(nwf,nwf),dwmn(nwf,nwf),
     &          eunk(nwf,nqbz))

c (2-0) construct initlal u~ from u
      call diag_unk(is,qbz,
     i              iko_ix,iko_fx,iko_i,iko_f,
     d              nband,nwf,nqbz,
     u              cnk,
     o              eunk)

!      call chk_cnkweight(qbz,iko_ix,iko_fx,cnk,
!     &     nqbz,nwf,nband,nlmto)
!
! check ortho-normality of u~'s
!      call chk_cnk(cnk,
!     i             iko_ix,iko_fx,iko_i,iko_f,
!     d             nband,nwf,nqbz)
!
! check: eunk vs. KS energy
!      call chk_eunk(is,qbz,eunk,ef,
!     d              nqbz,nband,nwf)


c (2-1) initial: uumat -> Mmn
      call init_mmn(cnk,uumat,ikbidx,
     i              iko_ix,iko_fx,iko_i,iko_f,ikbo_i,ikbo_f,
     d              nwf,nqbz,nbb,
     o              mmn0)


c (2-2) initial U
      call init_Umnk(amnk,cnk,
     i               iko_ix,iko_fx,iko_i,iko_f,
     d               nwf,nqbz,
     o               umnk)

!      call chk_umn(cnk,umnk,qbz,
!     i             iko_ix,iko_fx,iko_i,iko_f,
!     d             nwf,nqbz,nband,nlmto)

      call updt_mmn(umnk,mmn0,ikbidx,
     d              nwf,nqbz,nbb,
     u              mmn)

      do isc = 1, nsc2

c <r_n> ([1] eq.31)
         call get_rn(mmn,bb,wbb,wbz,
     d               nwf,nqbz,nbb,
     o               rn)

         do iq = 1,nqbz

            dwmn = (0d0,0d0)

            do ibb = 1,nbb
c (2-3) A[R] matrix
               call getrmn(mmn(:,:,ibb,iq),
     d                     nwf,
     o                     rmn)
               call getamn(rmn,
     d                     nwf,
     o                     amn)

c (2-4) S[T] matrix
               call gettmn(rn,mmn(:,:,ibb,iq),bb(:,ibb),
     d                     nwf,
     o                     qn,tmn)
               call getsmn(tmn,
     d                     nwf,
     o                     smn)

c DW(k) ([1] eq.57)
               dwmn(:,:) =  dwmn(:,:)
     &       + wbb(ibb) * (amn(:,:) - smn(:,:)) * alpha2 / wbbsum

c end of ibb-loop
            enddo

c (2-5) DW(k) -> U(k) ([1] eq.60)
            call updt_uk(dwmn,
     d                   nwf,
     u                   umnk(:,:,iq))
c            call chk_um(umnk(:,:,iq),nwf)

c end of iq-loop
         enddo   

c update Mmn ([1] eq.61)
         call updt_mmn(umnk,mmn0,ikbidx,
     d                 nwf,nqbz,nbb,
     u                 mmn)

c (2-6) Omeg_I, Omega_D and Omega_OD ([1] eq.34,35,36)
         call getOmg(mmn,rn,bb,wbb,wbz,
     d               nwf,nqbz,nbb,
     o               omgi,omgd,omgod,omgdod,omgidod)

c check self-consistency
         write(*,*)'#SC-loop, conv.',isc,omgdod
         write(*,950)'Omg: I, OD, D ',omgi,omgod,omgd
         if (isc.ge.2) then
            domgdod = dabs((omgdodold - omgdod) / omgdodold)
            if (domgdod .lt. conv2) then
               write(*,*) 'step2: converged!'
               goto 820
            endif
         endif
         omgdodold = omgdod

c end of self-consistent loop
      enddo
      write(*,*)'step2: not converged'
 820  continue

!      call chk_dnk(is,eunk,qbz,
!     i             umnk,cnk,
!     i             iko_ix,iko_fx,iko_i,iko_f,
!     d             nband,nwf,nqbz)
!
!      call chk_umn(cnk,umnk,qbz,
!     i             iko_ix,iko_fx,iko_i,iko_f,
!     d             nwf,nqbz,nband,nlmto)

c output
      call wmaxloc(ifmlw(is),ifmlwe(is),
     i             qbz,umnk,cnk,eunk,
     i             iko_ix,iko_fx,iko_i,iko_f,
     d             nwf,nqbz,nband,nlmto)
      call writeOmg(is,mmn,rn,bb,wbb,wbz,tpia,
     d             nwf,nqbz,nbb)
c 070824
      call getkeyvalue("GWinput","wan_write_rmn",lrmn,default=.false.)
      if (lrmn)
     & call writermn(is,mmn,bb,wbb,qbz,qbz0,wbz,rt,
     d              nwf,nqbz,nbb,n1,n2,n3)
c 070830
      call getkeyvalue("GWinput","wan_write_mmn",lmmn,default=.false.)
      if (lmmn)
     & call writemmn(is,mmn,bb,wbb,qbz,wbz,rt,
     d              nwf,nqbz,nbb,n1,n2,n3)

      deallocate(uumat,amnk,omgik,mmn,mmn0,
     &           rmn,amn,smn,rn,qn,tmn,dwmn)

c step 3 -- reduced Hamiltonian --
      write(*,*)'Step 3: reduced Hamiltonian branch'

c open file
      if (is .eq. 1) then
         ifbnd = iopen('bnds.maxloc.up',1,-1,0)         
         iftb  = iopen('bnds.tb.up',1,-1,0)         
      else
         ifbnd = iopen('bnds.maxloc.dn',1,-1,0)         
         iftb  = iopen('bnds.tb.dn',1,-1,0)         
      endif
      write(ifbnd,*)nq
      write(ifbnd,*)nwf
      write(iftb,*)nq
      write(iftb,*)nwf

c allocate
      allocate (hrotk(nwf,nwf,nqbz), ! hrotr(nwf,nwf,nqbz),
     o          hrotkp(nwf,nwf),evecc(nwf,nwf),eval(nwf))

c for small Hamiltonian
      call getkeyvalue("GWinput","wan_small_ham",lsh,default=.false.)
      if (lsh) then
         call getkeyvalue("GWinput","wan_nsh1",nsh1, default=1 )
         call getkeyvalue("GWinput","wan_nsh2",nsh2, default=2 )
         write(*,*)'SmallHam on',nsh1,nsh2
         nsh = nsh2 - nsh1 + 1
         if (is .eq. 1) then 
            ifsh = iopen('bnds.sh.up',1,-1,0)
         else                                 
            ifsh = iopen('bnds.sh.dn',1,-1,0)
         endif
         write(ifsh,*)nq
         write(ifsh,*)nsh
         allocate (hrotkps(nsh,nsh),eveccs(nsh,nsh),evals(nsh))
      endif

c (3-1) ~H(k) -> Hrot(k)
      call rot_hmnk(umnk,eunk,
     d              nwf,nqbz,
     o              hrotk)
c      stop 'rot_hmnk done'


c (3-2) Hrot_mn(R)
      allocate(irws(n1*n2*n3*8),rws(3,n1*n2*n3*8),drws(n1*n2*n3*8))
      call wigner_seitz(alat,plat,n1,n2,n3,nrws,rws,irws,drws)
      allocate(hrotr(nwf,nwf,nrws))

        if (ixc.eq.2) then

      call get_hrotr_ws(hrotk,qbz,wbz,
     i                  rws,irws,drws,
     d                  nwf,nqbz,nrws,
     o                  hrotr)
cskino
c write hrotr and *rws 
          if (is .eq. 1) then
            ifh = iopen('hrotr.up',1,-1,0)
          else
            ifh = iopen('hrotr.dn',1,-1,0)
          endif

          call write_hrotr(ifh, hrotr,
     i                  rws,irws,drws,
     d                  nwf,nrws )

          close (ifh) 

cekino

cskino
         else if (ixc.eq.3) then

          if (is .eq. 1) then
            filename='hrotr.up'
          else
            filename = 'hrotr.dn'
          endif
       
            call read_hrotr(filename,nwf,nrws,
     o          hrotr) 

         endif


          if (is .eq. 1) then
            ifh = iopen('hrotr.cut.up',1,-1,0)
          else
            ifh = iopen('hrotr.cut.dn',1,-1,0)
          endif

          allocate(hrotrcut(nwf,nwf,nrws))

          call make_hrotrcut( hrotr,
     i                  rws,irws,drws,
     i                  rcut,heps,
     d                  nwf,nrws,
     o                  hrotrcut )

          call write_hrotr(ifh, hrotrcut,
     i                  rws,irws,drws,
     d                  nwf,nrws )


          close (ifh)


cekino


c k-point mesh
      call get_nqbze(nqbz,nqbze)
      allocate(qbze(3,nqbze))
      call get_qbze(qbz,nqbz,
     o              qbze,nqbze)
      write(ifmlw(is))nqbze,nwf
      write(ifmlwe(is))nqbze,nwf
      do iq = 1,nqbze
          call get_hrotkp_ws(hrotr,rws,drws,irws,qbz(:,iq),
     d                     nwf,nqbz,nrws,
     o                     hrotkp)
c          call q2q0(qbze(:,iq),plat,q0(:))

          call diag_hm(hrotkp,nwf,eval,evecc)
!          call chk_diag(qbze(:,iq),hrotkp,nwf,eval,evecc)
!          if (iq.le.nqbz)
!     &    call chk_umnk(qbz(:,iq),hrotkp,nwf,eval,umnk(:,:,iq))
!          if (iq.le.nqbz)
!     &    call cmp_umn_evec(qbz(:,iq),umnk(:,:,iq),evecc,eval,nwf)

          call wmaxloc_diag(ifmlw(is),ifmlwe(is),
     i                 iq,qbze(1:3,iq),umnk,cnk,eunk,evecc,eval,
     i                 iko_ix,iko_fx,iko_i,iko_f,
     d                 nwf,nqbz)

c          write(6900,*)'***',iq
c          write(6900,*)evecc
cc          write(*,990)'iq =',iq,qbz(1:3,iq)
c          if (iq.le.nqbz) then
c          do iband = 1,nwf
c             e1 = (eval(iband)   -ef)*rydberg()
c             e2 = (eunk(iband,iq)-ef)*rydberg()
c             write(*,*)iband,e1,e2
c          enddo
c          endif
      enddo
      deallocate(qbze)

c other k-points
      write(ifbnd,*)ef,' ef'
      write(iftb,*)ef,' ef'
      if (lsh) write(ifsh,*)ef,' ef'
      do iq = 1,nq
c (3-3) Hrot_mn(k')
          call get_hrotkp_ws(hrotr,rws,drws,irws,q(:,iq),
     d                     nwf,nqbz,nrws,
     o                     hrotkp)
c          call q2q0(q(:,iq),plat,q0(:))

c (3-4) diagonalize
          call diag_hm(hrotkp,nwf,eval,evecc)

c          write(ifbnd,990)'iq =',iq,q(1:3,iq)
          write(ifbnd,*)'iq =',iq
          write(ifbnd,990)q(1:3,iq)
          do iband = 1,nwf
             write(ifbnd,*)iband,eval(iband)
             e1 = (eval(iband)   -ef)*rydberg()
c             e1 = (dreal(hrotkp(iband,iband))-ef)*rydberg()
c             write(ifbnd,*)iband,e1
          enddo

c (3-4) diagonalize  -- Small Hamiltonian --
          if (lsh) then
          hrotkps(1:nsh,1:nsh) = hrotkp(nsh1:nsh2,nsh1:nsh2)
          call diag_hm(hrotkps,nsh,evals,eveccs)

c          write(ifbnd,990)'iq =',iq,q(1:3,iq)
          write(ifsh,*)'iq =',iq
          write(ifsh,990)q(1:3,iq)
          do iband = 1,nsh
             write(ifsh,*)iband,evals(iband)
             e1 = (evals(iband)   -ef)*rydberg()
c             e1 = (dreal(hrotkps(iband,iband))-ef)*rydberg()
c             write(ifsh,*)iband,e1
          enddo
          endif ! lsh

c (3-3) Hrot_mn(k')  -- Tight-binding ---
         call get_hrotkp_tb_ws(rcut,plat,alat,
     i                     hrotr,rws,drws,irws,q(:,iq),
     d                     nwf,nqbz,nrws,
     o                     hrotkp)

c (3-4) diagonalize -- Tight-binding --
          call diag_hm(hrotkp,nwf,eval,evecc)

c          write(iftb,990)'iq =',iq,q(1:3,iq)
          write(iftb,*)'iq =',iq
          write(iftb,990)q(1:3,iq)
          do iband = 1,nwf
             write(iftb,*)iband,eval(iband)
             e1 = (eval(iband)   -ef)*rydberg()
c             write(iftb,*)iband,e1
          enddo
c end-of iq-loop
      enddo

c 080707
      call writeham(ifham,is,ef,alat,plat,pos,qbz,wbz,rws,irws,hrotk,
     &              nspin,natom,nwf,nqbz,nrws)
      deallocate(cnk,umnk,eunk,hrotk,hrotr,hrotkp,evecc,eval,
     &           irws,rws,drws)
      if (lsh) deallocate(hrotkps,eveccs,evals)

      close(ifbnd)

c end of loop over spin
      enddo

 950  format(a14,3f23.16)
c 990  format(a4,i5,3f12.6)
 990  format(3f12.6)

c-----------------------------------------------------------------------
      call cputid(0)
      stop 'hmaxloc: ixc=2 ok'
      end

c-----------------------------------------------------------------------
      subroutine chk_amnkweight(qbz,iko_ix,iko_fx,amnk,
     &     nqbz,nwf,nband,nlmto)
      use m_readqg
      use m_readeigen
      use keyvalue
      implicit real*8(a-h,o-z)

      complex(8) :: amnk(iko_ix:iko_fx,nwf,nqbz)
      complex(8),allocatable:: cphi1(:,:),cphi2(:,:)
      real(8) :: qbz(3,nqbz),q(3),quu(3)
      real(8),allocatable:: wbas(:,:)
      integer(4),allocatable:: 
     &  m_indx(:),n_indx(:),l_indx(:),ibas_indx(:)

c --- Readin nlam index
      ifoc = iopen('@MNLA_CPHI',1,0,0)
      ldim2 = nlmto
      read(ifoc,*)
      allocate(m_indx(ldim2),n_indx(ldim2),l_indx(ldim2),ibas_indx(ldim2))
      do ix =1,ldim2
        read(ifoc,*)m_indx(ix),n_indx(ix),l_indx(ix),ibas_indx(ix),ixx
        if(ixx/=ix) stop  'failed to readin @MNLA_CPHI'
      enddo   

      nbas = ibas_indx(nlmto)
      allocate(cphi1(nlmto,nband),cphi2(nlmto,nwf),wbas(nbas,nwf))
      wbas = 0d0
      cphi2=0d0
      do iq = 1,nqbz
      q = qbz(:,iq)
      call readcphi(q,nlmto,1,quu,cphi1)

       do iwf=1,nwf
         do ib=iko_ix,iko_fx
            cphi2(:,iwf) = cphi2(:,iwf) + cphi1(:,ib)*amnk(ib,iwf,iq)
         enddo
      enddo

      enddo ! iq

      do iwf=1,nwf
         do ia=1,nlmto
            ibas = ibas_indx(ia)
            wbas(ibas,iwf) = wbas(ibas,iwf) + 
     &                 conjg(cphi2(ia,iwf))*cphi2(ia,iwf)
         enddo ! ia
      enddo ! iwf
      wbas = wbas / dble(nqbz**2)

      write(*,*)'*** ibas,iwf,wbas'
      do iwf=1,nwf
      do ibas=1,nbas
         write(*,*)ibas,iwf,wbas(ibas,iwf)
      enddo
      write(*,*)
      enddo   
      write(*,*)'*** ibas,wbas'
      do ibas=1,nbas
         w = 0d0
         do iwf=1,nwf
            w = w + wbas(ibas,iwf)
         enddo
         write(*,*)ibas,w
      enddo
   
      deallocate(cphi1,cphi2,wbas,m_indx,l_indx,n_indx,ibas_indx)
      ix = iclose('@MNLA_CPHI')

      end
c-----------------------------------------------------------------------
      subroutine chk_cnkweight(qbz,iko_ix,iko_fx,cnk,
     &     nqbz,nwf,nband,nlmto)
      use m_readqg
      use m_readeigen
      use keyvalue
      implicit real*8(a-h,o-z)

      complex(8) :: cnk(iko_ix:iko_fx,nwf,nqbz)
      complex(8),allocatable:: cphi1(:,:),cphi2(:,:)
      real(8) :: qbz(3,nqbz),q(3),quu(3)
      real(8),allocatable:: wbas(:,:)
      integer(4),allocatable:: 
     &  m_indx(:),n_indx(:),l_indx(:),ibas_indx(:)

c --- Readin nlam index
      ifoc = iopen('@MNLA_CPHI',1,0,0)
      ldim2 = nlmto
      read(ifoc,*)
      allocate(m_indx(ldim2),n_indx(ldim2),l_indx(ldim2),ibas_indx(ldim2))
      do ix =1,ldim2
        read(ifoc,*)m_indx(ix),n_indx(ix),l_indx(ix),ibas_indx(ix),ixx
        if(ixx/=ix) stop  'failed to readin @MNLA_CPHI'
      enddo   

      nbas = ibas_indx(nlmto)
      allocate(cphi1(nlmto,nband),cphi2(nlmto,nwf),wbas(nbas,nwf))
      wbas = 0d0
      cphi2=0d0

      do iq = 1,nqbz
      q = qbz(:,iq)
      call readcphi(q,nlmto,1,quu,cphi1)

      do iwf=1,nwf
         do ib=iko_ix,iko_fx
            cphi2(:,iwf) = cphi2(:,iwf) + cphi1(:,ib)*cnk(ib,iwf,iq)
         enddo
      enddo

      enddo ! iq

      do iwf=1,nwf
         do ia=1,nlmto
            ibas = ibas_indx(ia)
            wbas(ibas,iwf) = wbas(ibas,iwf) + 
     &                 conjg(cphi2(ia,iwf))*cphi2(ia,iwf)
         enddo ! ia
      enddo ! iwf
      wbas = wbas / dble(nqbz*nqbz)

      write(*,*)'*** ibas,iwf,wbas'
      do iwf=1,nwf
      do ibas=1,nbas
         write(*,*)ibas,iwf,wbas(ibas,iwf)
      enddo
      write(*,*)
      enddo   
      write(*,*)'*** ibas,wbas'
      do ibas=1,nbas
         w = 0d0
         do iwf=1,nwf
            w = w + wbas(ibas,iwf)
         enddo
         write(*,*)ibas,w
      enddo
   
      deallocate(cphi1,cphi2,wbas,m_indx,l_indx,n_indx,ibas_indx)
      ix = iclose('@MNLA_CPHI')

      end
c-----------------------------------------------------------------------
      subroutine chk_umn(cnk,umnk,qbz,
     i                  iko_ix,iko_fx,iko_i,iko_f,
     d                  nwf,nqbz,nband,nlmto)
      use m_readqg
      use m_readeigen
      use keyvalue
      implicit real*8(a-h,o-z)

      complex(8) :: cnk(iko_ix:iko_fx,nwf,nqbz),
     &              dnk(iko_ix:iko_fx,nwf,nqbz),
     &              umnk(nwf,nwf,nqbz)
      real(8) :: qbz(3,nqbz)
      integer(4) :: iko_i(nqbz),iko_f(nqbz)

      dnk = (0d0,0d0)
      do iq = 1,nqbz
         do imp = iko_i(iq),iko_f(iq)
         do in = 1,nwf
            do im = 1,nwf
               dnk(imp,in,iq) = dnk(imp,in,iq) 
     &                 + umnk(im,in,iq) * cnk(imp,im,iq)
            enddo ! im
         enddo ! in
         enddo ! imp      
      enddo ! iq  

      call chk_cnkweight(qbz,iko_ix,iko_fx,dnk,
     &     nqbz,nwf,nband,nlmto)

      end
c-----------------------------------------------------------------------
