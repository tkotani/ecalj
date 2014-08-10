      program hibz2fbz0

c-------------------------------------------------
c q(IBZ) to k(FBZ)
c SE_nm(k) = <kn| SE |km> = O(kn,qi) SE_ij(q) conj(O(km,qj))
c Vxc_nm(k) = <kn| Vxc |km> = O(kn,qi) Vxc_ij(q) conj(O(km,qj))
c O(kn,qn') = < psi_kn(r) | psi_qi(R^-1 r) >
c k = Rq (R: symmetry natrix)
c
c Takahi Miyake, May, 2008
c from hibz2fbz.F
c-------------------------------------------------
      use m_readqg
      use m_readeigen
      use m_read_bzdata,ngrp2=>ngrp
      use m_genallcf_v3
      use keyvalue

      implicit none
      real(8):: q(3),  qgbin(3),qx(3)
      integer(4),allocatable :: ngvecpB(:,:,:),ngveccB(:,:) !,ngveccB(:,:,:)
     & , ngvecpf1(:,:), ngvecpf2(:,:), 
     &   nx(:,:),nblocha(:),ifppb(:) !ongveccBr(:,:,:)
      real(8),allocatable :: ppbrd (:,:,:,:,:,:,:),cg(:,:,:),symope(:,:),
     &phij(:),psij(:),rprodx(:,:),rphiphi(:),q0i(:,:),
     &qbzs(:,:),qbz2(:,:)
      complex(8),parameter:: img=(0d0,1d0)
c,nw,incwf,natom,nclass,ipos,igrp,
c     & iinvg,nspin,nl,nn,nnv,nnc,
c     o                   inindx,inindxv,inindxc,iiclass,             !l,n, dimensions 
c     d                   nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, !l,n,  dimensions
c     o                   izdummy, 
c     o   iil,iin,iim,   iilnm, i_mnl, ! l,n,m for Phi ! w(i_mnl)=> mnl(ic) for all electron
c     o   iilv,iinv,iimv,iilnmv,i_mnlv,! l,n,m for Phi
c     o   iilc,iinc,iimc,iilnmc,i_mnlc,! l,n,m for Phi
c     o   iecore,ikonf,iicore,incore,nctot,             !core
c     o   imagw_dummy,niw,idummy,
      integer(4)
     &   nw_input,
     &   ifhbe,
     &   nprecb,mrecb,mrece,nlmtot,nqbzt,nband,
     &   nq0i,i,nq0ix,neps,ngrpmx,ngcmx,mxx,nqbze,nqibze,ini,ix,ngrpx
     &  ,mdimx,nbloch,nblochpmx,ifvcfpout,ndummy1,ndummy2,ifcphi,is,nwp,
     &   ifepscond,nxx,ifvxcpout,ifgb0vec
     &   ,nw0,iw,nwhis,ifinin,nw2,iw0,ifwwk,noccxv,noccx
     &   ,ifemesh,nprecx,mrecl,ifwd,ifrcwi,ifrcw,nspinmx,ifianf,ibas
     &   ,ibas1,irot,iq,ngb,iqixc2,ifepsdatnolfc,ifepsdat,ngbin,igc0
     &   ,kx,isf,kqxx,kp,job,nbnbx,nhwtot,noccxvx,nwmax  !,ifev1,ifev2
     &   ,ihis,jhwtot,ik,ibib,ib1,ib2,ichkhis,ihww,j,imode
     &   ,ngpmx

      real(8):: dum1,dum2,dum3,wqtsum,epsrng,dnorm,dini,
     & dwry,dwh,omg_c,omg2,xxx
      integer(4)::nwin, incwfin,  verbose
      real(8)::efin
      integer(4):: bzcase, mrecg,ifphi,
     & nbas,nradmx,ncoremx,nrx,ic,icx,isp,l,n,irad,ifoc,
     & ldim2,ixx,ngp1,ngp2,nq0it
      real(8):: qq(3),quu(3), deltaq(3),q1x(3),q2x(3)
      real(8),parameter::  pi =     3.14159265358979323846d0
      real(8),parameter::  fpi =    4d0*pi

c      logical:: test_qbzreg
      logical:: qbzreg
!-------------------------------------------------------------------------
      integer(4),allocatable:: ncindx(:,:),
     &           lcindx(:,:),
     &           nrad(:),
     &           nindx_r(:,:),
     &           lindx_r(:,:),
     &           nc_max(:,:),
     &  m_indx(:),n_indx(:),l_indx(:),ibas_indx(:), nrofi(:)
      real(8),allocatable:: phitoto(:,:,:,:,:), aa(:),rr(:,:)
     &                     ,phitotr(:,:,:,:,:),
     &        bb(:),zz(:),rmax(:),cy(:),yl(:)


      complex(8),allocatable:: geig1(:,:),geig2(:,:),cphi1(:,:),cphi2(:,:)
     & ,uum(:,:)
      complex(8):: ppj,phaseatom
      real(8)   :: q1(3),q2(3),dq(3),absqg2,absdq,r2s,absqg
      integer(4):: j1,j2,j1max,j2max,j1min,j2min,ispin
     & ,l1,l2,lm1,lm2,ibas2,lm3,ig1,ig2,ir,ia,ia1,ma,ia2,m2,l3,m1,lxx
     &, iopen,ico,lxd,lx !,ll
      real(8):: ylk
#ifdef COMMONLL
      integer(4):: ll(51**2)
      common/llblock/ll
#else
      integer(4) ll
#endif  
ctm
      integer(4),allocatable:: kpk(:),kprot(:),itq(:),itq2(:)
      integer(4) :: ixc,idummy,idummy2,i1,i2,i3,ifpp,ifkkr
     &              ,iqibz,iqbz,iqbz2,itmp,itmp2,iti,itf
     &              ,nqibz2,nqbz2,ntq,ntq2,nq2,nspin2,iqtmp,nadd(3),invr
     &              ,lmx,iii,ikp,ikr,kount,mrecp
     &              ,ifo,ifmlw,ifmlwe,iko_ix,iko_fx,nwf
     &              ,iko_ixs,iko_fxs,ifwand,nwf2
     &              ,iko_ix_wan,iko_fx_wan
      integer(4) :: nword,invrot
      real(8) :: uunorm,qdiff(3),add(3),qbasinv(3,3),det
     &          ,qibz_k(3),qbz_kr(3)
      complex(8),allocatable:: dnk(:,:,:)
      logical:: lallb,lwan,lsc

c space group infermation 
      integer(4) :: ibasx
      integer(4),allocatable :: iclasst(:), invgx(:), miat(:,:)
      real(8),allocatable    :: tiat(:,:,:),shtvg(:,:)

c
      integer(4) :: mrect,ntau,ntau2,itau,iflda,ifvxc,ifsxc,ifsx,ifsc,
     &              mrecs,nrec,isx,iclose,ngrp3,ikp2,ii,ij,ii2,ij2,ntmp
     &             ,ifv_xc,ifevec,ndimh,nqqq,nqqqx,ndummy,itp,itpp
     &             ,ip,it1,it2,iv,nfreq,ifvxc2,ifsxc2,ifsx2,ifsc2
     &             ,iwini,iwend
      real(8) :: beta,eflda,qqqx(3),wn,vn,rcos,rsin,ccos,csin,
     &           rsum,csum,rydberg,ef2,dwplot,rtmp
      complex(8) :: ztmp
      real(8),allocatable:: tau(:),qqq(:,:),elda(:,:,:),qqqxx(:,:),
     &               qibz2(:,:),eq2(:,:,:)
      complex(8),allocatable:: zsec2(:,:,:,:),zw(:,:),zw2(:,:)
     &               ,zw3(:,:),zw4(:,:),w1(:,:,:),w2(:),zwexp(:),
     &                vxc(:,:), sx(:,:),sxc(:,:),dsx(:,:),sc(:,:)
     &               ,vxcf(:,:,:,:),v_xc(:,:,:,:),evec(:,:,:,:)
c-------------------------

c mode switch. --------------
c      write(6,*) ' --- Choose omodes below -------------------'
c      write(6,*) '  (1) q to k,  (2) project d out, (3) ???' 
c      call readin5(ixc,idummy,idummy2)
c      write(6,*) ' ixc=',ixc
c      if(ixc==0) stop ' --- ixc=0 --- Choose computational mode!'


c---  readin BZDATA. See gwsrc/rwbzdata.f
c--------readin data set when you call read_BZDATA ---------------
c       integer(4)::ngrp,nqbz,nqibz,nqbzw,nteti,ntetf,
c     &   n_index_qbz  
c       integer(4):: n1,n2,n3
c       real(8):: qbas(3,3),ginv(3,3),qbasmc(3,3),dq_bzcase2
c       real(8),allocatable:: qbz(:,:),wbz(:),qibz(:,:)
c     &    ,wibz(:),qbzw(:,:)
c       integer(4),allocatable:: idtetf(:,:),ib1bz(:),idteti(:,:)
c     &    ,nstar(:),irk(:,:),index_qbz(:,:,:)
c-----------------------------------------------------------------
      call read_BZDATA()
      call dinv33(qbas,0,qbasinv,det) 

c--- Use regular mesh even for bzcase==2
      allocate(qbzs(3,nqbz))
      if(bzcase()==2.and.qbzreg()) then
      deltaq= qbas(:,1)/n1 + qbas(:,2)/n2 +qbas(:,3)/n3
      do i=1,nqbz
        qbzs(:,i) = qbz(:,i) -deltaq/2d0
        write(6,"('i qbz=',i3,3f8.3)") i,qbz(:,i)
      enddo
ctm 060916
c      stop "hibz2fbz: bzcase=2 is not implemented."
      endif
      print *,' ======== nqbz qbz  =',nqbz
      print *,' ======== nqibz ngrp=',nqibz,ngrp
      print *,  qbz 
      print *,'============================'
      print *

C--- readin GWIN and LMTO, then allocate and set datas.
      nwin = -999    !readin condition. Not readin NW file
      incwfin= 0     !readin condition. use ForX0 for core in GWIN
      efin =  -999d0 !readin condition. Not readin EFERMI
      call genallcf_v3(nwin,efin,incwfin) !in module m_genallcf_v3
      if(ngrp/= ngrp2) stop 'ngrp inconsistent: BZDATA and LMTO GWIN_V2'
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

c --- read dimensions of h,hb
      ifhbe      = iopen('hbe.d',1,0,0)
      read (ifhbe,*) nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg
      if(nlmto/=nlmtot) stop ' hibz2fbz: nlmto/=nlmtot in hbe.d'
      if(nqbz /=nqbzt ) stop ' hibz2fbz: nqbz /=nqbzt  in hbe.d'

c --- input file
c      ifmlw  = iopen('MLWU',0,0,0)
c      read(ifmlw)nqbz2,nwf,iko_ix,iko_fx
c      isx = iclose('MLWU')
!      call getkeyvalue("GWinput","Wannier",lwan,default=.false.)
      call getkeyvalue("GWinput","SelfConsistent",lsc,default=.false.)
!      if (lwan) then
         write(*,*)'Read wan.d'
         ifwand  = iopen('wan.d',-1,0,0)
         read(ifwand,*)nqbz2,nwf,iko_ix_wan,iko_fx_wan
         isx = iclose('wan.d')
         if (nqbz2.ne.nqbz) stop 'hibz2fbz: wrong nqbz in wan.d'
!      else
!         stop 'hibzfbz: not .wannier.'
!      endif
      if (lsc) then
         write(*,*)'Self-consistent calc.'
         iko_ix = 1
         iko_fx = nband
      else
         iko_ix = iko_ix_wan
         iko_fx = iko_fx_wan
      endif   
      ntq = iko_fx - iko_ix + 1
      allocate( itq(ntq) )
      do i   = 1,ntq
        itq(i) = i + (iko_ix - 1)
      enddo

c      mrecp = 16*nband*nband/nword()
      mrecp = 16*ntq*ntq/nword()
      ifpp     = iopen('PSIPSI',0,0,mrecp)

      ifkkr     = iopen('kkrot',1,0,0)
      read(ifkkr,*) nqbz2,nqibz2,ngrp3
      if (nqbz2.ne.nqbz) stop 'hibz2fbz: wrong nqbz'
      if (nqibz2.ne.nqibz) stop 'hibz2fbz: wrong nqibz'
      if (ngrp3.ne.ngrp) stop 'hibz2fbz: wrong ngrp'
      allocate (kpk(nqbz),kprot(nqbz))
      do ikp = 1,nqbz
        read(ifkkr,*) ikp2,kpk(ikp),kprot(ikp)
      enddo  

!c --- read tau mesh from file 'TAU'
!      open(UNIT=3111,file='TAU')
!      read(3111,*)ntau
!      allocate(tau(ntau))
!      do itau=1,ntau
!        read(3111,*)tau(itau)
!      enddo
!      read(3111,*)beta
!      close(3111)

!c --- read HLDA
!      iflda  = iopen('HLDA',0,0,0)
!      read(iflda)ntq2,nqibz2,nqbz2,nspin2,eflda
!      if (ntq2.ne.ntq) stop 'hibz2fbz: wrong ntq in HLDA'
!      if (nqibz2.ne.nqibz) stop 'hibz2fbz: wrong nqibz in HLDA'
!      if (nqbz2.ne.nqbz) stop 'hibz2fbz: wrong nqbz in HLDA'
!      if (nspin2.ne.nspin) stop 'hibz2fbz: wrong nspin in HLDA'
!      allocate(qqq(3,nqibz),itq2(ntq),elda(ntq,nqibz,nspin))
!      read(iflda)qqq
!      read(iflda)itq2
!      read(iflda)elda
!      do ii = 1,ntq
!        if (itq(ii).ne.itq2(ii)) stop 'hibz2fbz: wrong itq in HLDA'
!      enddo
!      do iq = 1,nqibz
!        call checkagree(qqq(:,iq),qibz(:,iq),' q error in HLDA')
!      enddo
!      deallocate(qqq)

c --- input files
      mrecs = 16*ntq*ntq/nword()
      ifvxc  = iopen('XC2U',1,0,0)
      ifsx   = iopen('SEX2U',1,0,0)
      ifsc   = iopen('SECW2U',1,0,0)
      ifsxc  = iopen('SEXcore2U',1,0,0)

      ifvxc2  = iopen('XC2Uda',0,-1,mrecs)
      ifsx2   = iopen('SEX2Uda',0,-1,mrecs)
      ifsc2   = iopen('SECW2Uda',0,-1,mrecs)
      ifsxc2  = iopen('SEXcore2Uda',0,-1,mrecs)
      allocate(vxc(ntq,ntq),sx(ntq,ntq),
     &         sxc(ntq,ntq),sc(ntq,ntq),dsx(ntq,ntq))
      allocate( uum(ntq,ntq) )

c --- output file
      ifo     = iopen('GW0',1,-1,0)

c-----------------
c XC2U
      read(ifvxc,*)ntq2,nqibz2,nqbz2,nspin2,ef2
      if (ntq2.ne.ntq) stop 'hibz2fbz0: ntq error'
      if (nqibz2.ne.nqibz) stop 'hibz2fbz0: nqibz error'
      if (nqbz2.ne.nqbz) stop 'hibz2fbz0: nqbz error'
      if (nspin2.ne.nspin) stop 'hibz2fbz0: nspin error'
      allocate(qibz2(3,nqibz2),itq2(ntq2),eq2(ntq2,nqibz2,nspin2))
      read(ifvxc,*)qibz2
      read(ifvxc,*)itq2
      read(ifvxc,*)eq2
      do ik = 1,nqibz2
        read(ifvxc,*)
        do ij = 1,ntq
        do ii = 1,ntq
          read(ifvxc,*)ii2,ij2,vxc(ii,ij)
        enddo  
        enddo  
!        nrec = nqibz*(ispin-1)+ik
        nrec = ik
        write(ifvxc2,rec=nrec) vxc
      end do ! ik

c SEXcore2U
      read(ifsxc,*)ntq2,nqibz2,nqbz2,nspin2,ef2
      if (ntq2.ne.ntq) stop 'hibz2fbz0: ntq error'
      if (nqibz2.ne.nqibz) stop 'hibz2fbz0: nqibz error'
      if (nqbz2.ne.nqbz) stop 'hibz2fbz0: nqbz error'
      if (nspin2.ne.nspin) stop 'hibz2fbz0: nspin error'
      read(ifsxc,*)qibz2
      read(ifsxc,*)itq2
      read(ifsxc,*)eq2
      do ii=1,ntq
         if (itq2(ii).ne.itq(ii)) stop 'hibz2fbz0: itq error'
      enddo
      do ik = 1,nqibz2
        read(ifsxc,*)
        do ij = 1,ntq2
        do ii = 1,ntq2
           read(ifsxc,*)ii2,ij2,sxc(ii,ij)
        enddo   
        enddo   
!        nrec = nqibz*(ispin-1)+ik
        nrec = ik
        write(ifsxc2,rec=nrec) sxc
      end do ! ik

c SEX2U
      read(ifsx,*)ntq2,nqibz2,nqbz2,nspin2,ef2
      if (ntq2.ne.ntq) stop 'hibz2fbz0: ntq error'
      if (nqibz2.ne.nqibz) stop 'hibz2fbz0: nqibz error'
      if (nqbz2.ne.nqbz) stop 'hibz2fbz0: nqbz error'
      if (nspin2.ne.nspin) stop 'hibz2fbz0: nspin error'
      read(ifsx,*)qibz2
      read(ifsx,*)itq2
      read(ifsx,*)eq2
      do ii=1,ntq
         if (itq2(ii).ne.itq(ii)) stop 'hibz2fbz0: itq error'
      enddo
      do ik = 1,nqibz2
        read(ifsx,*)
        do ij = 1,ntq2
        do ii = 1,ntq2
           read(ifsx,*)ii2,ij2,sx(ii,ij)
        enddo   
        enddo   
!        nrec = nqibz*(ispin-1)+ik
        nrec = ik
        write(ifsx2,rec=nrec) sx
      end do ! ik

c SECW2U
      read(ifsc,*)ntq2,nqibz2,nqbz2,nspin2,ef2,dwplot,iwini,iwend
      if (ntq2.ne.ntq) stop 'hibz2fbz0: ntq error'
      if (nqibz2.ne.nqibz) stop 'hibz2fbz0: nqibz error'
      if (nqbz2.ne.nqbz) stop 'hibz2fbz0: nqbz error'
      if (nspin2.ne.nspin) stop 'hibz2fbz0: nspin error'
      read(ifsc,*)qibz2
      read(ifsc,*)itq2
      read(ifsc,*)eq2
      do ii=1,ntq
         if (itq2(ii).ne.itq(ii)) stop 'hibz2fbz0: itq error'
      enddo
      allocate(zsec2(ntq,ntq,nqibz,iwini:iwend))
      read(ifsc,*)zsec2
      nfreq = iwend-iwini+1
      do iw = 1,nfreq
      do ik = 1,nqibz2
!        nrec = nfreq*(nqibz*(ispin-1)+ik-1)+iw
        nrec = nfreq*(ik-1)+iw
        write(ifsc2,rec=nrec) zsec2(:,:,ik,iwini+iw-1)
      end do ! ik
      end do ! iw

c reopen
      ii = iclose('XC2Uda')
      ii = iclose('SEX2Uda')
      ii = iclose('SECW2Uda')
      ii = iclose('SEXcore2Uda')
      ifvxc2  = iopen('XC2Uda',0,0,mrecs)
      ifsx2   = iopen('SEX2Uda',0,0,mrecs)
      ifsc2   = iopen('SECW2Uda',0,0,mrecs)
      ifsxc2  = iopen('SEXcore2Uda',0,0,mrecs)

c --- initiallization to get eigenfunctions
      call init_readeigen(ginv,nspin,nband,mrece) !initialization of readEigen
      call init_readeigen2(mrecb,nlmto,mrecg)
      call readngmx('QGpsi',ngpmx)
      allocate( geig1(ngpmx,nband),geig2(ngpmx,nband))
      print *, 'end of initialization'

c output: header
      write(ifo,*)nwf,nqbz,n1,n2,n3,iwini,iwend,dwplot,nspin,ef2,alat
      write(ifo,*)plat
      write(ifo,*)qbas
      write(ifo,*)qbz
!      write(ifo,*)tau

      do 1010 ispin = 1,nspin
        if (ispin.eq.2) stop'spin polarized case not implemented'

c open maxloc file
        if (ispin.eq.1) then
          ifmlw  = iopen('MLWU',0,0,0)
        else
          ifmlw  = iopen('MLWD',0,0,0)
        endif
        read(ifmlw)nqbz2,nwf2,iko_ixs,iko_fxs
        if (nqbz2.ne.nqbz) stop 'hibz2fbz: wrong nqbz'
        if (nwf2.ne.nwf) stop 'hibz2fbz: wrong nwf'
c        if (iko_ixs.lt.itq(1)) stop 'hibz2fbz: wrong iko_ixs'
c        if (iko_fxs.gt.itq(ntq)) stop 'hibz2fbz: wrong iko_fxs'
        if (iko_ixs.ne.iko_ix_wan) stop 'hibz2fbz: wrong iko_ixs'
        if (iko_fxs.ne.iko_fx_wan) stop 'hibz2fbz: wrong iko_fxs'
        allocate(dnk(ntq,nwf,nqbz),zw(nwf,nwf),zw2(iko_ixs:iko_fxs,nwf))
        dnk = 0d0

      do 1070 ikp = 1,nqbz ! k'
c        write(*,*)'ikp =',ikp, 'out of',nqbz

        irot = kprot(ikp)
        ik = kpk(ikp)
        ikr = irk(ik,irot)
        if (ikp.ne.ikr) stop 'hibz2fbz: kp .ne. rot*k'
        qbz_kr= qbz (:,ikr)
        qibz_k= qibz(:,ik)

c read maxloc wannier
        read(ifmlw)iqbz2,q(1:3)
c        read(ifmlw)dnk(1:ntq,1:nwf,ikp)
        read(ifmlw)zw2(iko_ixs:iko_fxs,1:nwf)
        ntmp = iko_ix - 1
        do ii = iko_ixs,iko_fxs
          dnk(ii-ntmp,1:nwf,ikp) = zw2(ii,1:nwf)
        enddo
        call checkagree(q,qbz(:,ikp),' q error in MLW')

c read psipsi
        nrec = nqbz*(ispin-1)+ikp
        read(ifpp,rec=nrec) uum(1:ntq,1:ntq)
c debug: >
!        write(*,*)'iko_ixs,iko_fxs,ntq'
!        write(*,*)iko_ixs,iko_fxs,ntq
!        write(81,"('***',i5,6f12.6)")ikp,qbz_kr,qibz_k
!        write(82,"('***',i5,6f12.6)")ikp,qbz_kr,qibz_k
!        do ii = 1,10
!           write(81,"(10f10.4)")real(uum(ii,1:10))
!           write(82,"(10f10.4)")imag(uum(ii,1:10))
!        enddo   
        write(83,"('***',i5,6f12.6)")ikp,qbz_kr,qibz_k
        do ii = 1,3
           rtmp = sum(abs(uum(ii,:))**2)
           write(83,*)ii,rtmp
        enddo
c debug: <

c Hlda
        call ldanm(eq2(:,ik,ispin),dnk(:,:,ikp),ntq,nwf,zw)
        write(ifo,*)zw
c debug: 
c        write(999,"(i5,3f12.6)")ikp,eq2(1:3,ik,ispin)

c Vxc, SEXcore, SEX
        nrec = nqibz*(ispin-1)+ik
        read(ifvxc2,rec=nrec) vxc
c debug: >
        vxc = 0d0
        do ii = 1,ntq
        do ij = 1,ntq
           vxc(ii,ii) = 1d0
        enddo   
        enddo   
c debug: <
        call senm(vxc,uum,dnk(:,:,ikp),ntq,nwf,zw)

c debug: >
        if (abs(qbz_kr(3)).lt.1.d-4) then
        write(81,"('***',i5,6f12.6)")ikp,qbz_kr,qibz_k
        write(82,"('***',i5,6f12.6)")ikp,qbz_kr,qibz_k
        do ii = 1,3
           write(81,"(3f10.4)")real(zw(ii,1:3))
           write(82,"(3f10.4)")imag(zw(ii,1:3))
        enddo   
        endif
c debug: <


        zw = zw*rydberg()*2d0
        write(ifo,*)zw
        read(ifsxc2,rec=nrec) sxc
        read(ifsx2,rec=nrec) sx
c        dsx = sx + sxc - vxc
        dsx = sx + sxc
        call senm(dsx,uum,dnk(:,:,ikp),ntq,nwf,zw)
        zw = zw*rydberg()*2d0
        write(ifo,*)zw

c SEC
        nfreq = iwend-iwini+1
        do iw = 1,nfreq
          nrec = nfreq*(nqibz*(ispin-1)+ik-1)+iw
          read(ifsc2,rec=nrec) sc
          call senm(sc,uum,dnk(:,:,ikp),ntq,nwf,zw)
          zw = zw*rydberg()*2d0
          write(ifo,*)zw
        enddo

 1070 continue ! ikp

      if (ispin.eq.1) then
         isx  = iclose('MLWU')
      else
         isx  = iclose('MLWD')
      endif
      deallocate(dnk,zw,zw2)

 1010 continue ! ispin


      print *, ' ====== end ========================================'
      end
c-----------------------------------------------------------------------
      subroutine checkagree(a,b,char)
      real(8):: a(3),b(3) 
      character*(*) :: char
      if(sum(abs(a-b))>1d-6) then
         print *,' Error in checkagree:',char
         stop ' Error in checkagree:'
      endif
      end
c-----------------------------------------------------------------------
      subroutine ldanm(eig,evec,ni,no,zw)

      implicit none
      integer ni,no,in,im,ii
      complex(8) :: evec(ni,no), zw(no,no)
      real(8) :: eig(ni),eps
      parameter (eps=1d-4)

      zw = 0d0
      do in = 1,no
      do im = 1,no
         zw(in,im) = sum(dconjg(evec(:,in))*evec(:,im)*eig(:))
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine senm(se,uum,evec,ni,no,zw)

      implicit none
      complex(8) :: se(ni,ni),uum(ni,ni),evec(ni,no), zw(no,no)
     &            ,sek(ni,ni),zw2(ni,ni)
      real(8) :: eps
      integer ni,no,in,im,ii
      parameter (eps=1d-4)

c IBZ to FBZ
c      sek = 0d0
c      do in = 1,ni
c      do im = 1,ni
c         sek(in,im) = sum(uum(1:ni,in)
c    &      *matmul(se(1:ni,1:ni),dconjg(uum(1:ni,im))) )
c      enddo
c      enddo

      zw2(:,:) = matmul(uum(:,:),se(:,:))
      sek = 0d0
      do in = 1,ni
      do im = 1,ni
        do ii = 1,ni
          sek(in,im) = sek(in,im) + zw2(in,ii)*dconjg(uum(im,ii))
        enddo
      enddo
      enddo

c KS to wannier
      zw = 0d0
      do in = 1,no
      do im = 1,no
         zw(in,im) = sum(dconjg(evec(1:ni,in))
     &      *matmul(sek(1:ni,1:ni),evec(1:ni,im)) )
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine senm2(se,uum,evec,ni,no,zw)

      implicit none
      complex(8) :: se(ni,ni),uum(ni,ni),evec(ni,no), zw(no,no)
     &            ,sek(ni,ni),zw2(ni,ni),zsum
      real(8) :: eps
      integer ni,no,in,im,ii
      parameter (eps=1d-4)

c IBZ to FBZ
c      sek = 0d0
c      do in = 1,ni
c      do im = 1,ni
c         sek(in,im) = sum(uum(1:ni,in)
c    &      *matmul(se(1:ni,1:ni),dconjg(uum(1:ni,im))) )
c      enddo
c      enddo

      zw2(:,:) = matmul(uum(:,:),se(:,:))
      sek = 0d0
      do in = 1,ni
      do im = 1,ni
        do ii = 1,ni
          sek(in,im) = sek(in,im) + zw2(in,ii)*dconjg(uum(im,ii))
        enddo
      enddo
      enddo

c KS to wannier
      zw = 0d0
      do in = 1,no
      do im = 1,no
         zw(in,im) = sum(dconjg(evec(1:ni,in))
     &      *matmul(sek(1:ni,1:ni),evec(1:ni,im)) )
      enddo
      enddo

c debug: >
      do in = 1,no
         write(61,"(i5,6f12.6)")in,zw(in,in),sek(in,in)
      enddo
      zsum = 0d0
      do in = 1,no
         zsum = zsum + zw(in,in)
      enddo   
      write(61,"(2f12.6)")zsum

c      write(61,*)'==='
c
c      do in = 1,no
c         write(61,"(i5,2f12.6)")in,sum(evec(:,in)*dconjg(evec(:,in)))
c      enddo
c      write(61,*)'==='
c debug: <

      return
      end
c-----------------------------------------------------------------------
