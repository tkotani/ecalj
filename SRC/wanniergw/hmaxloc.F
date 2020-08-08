
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
c      use m_readqg,only: readngmx,readqg
      use m_readeigen,only: init_readeigen,init_readeigen2,readeval
      use m_hamindex,only: qtt,nqtt
      use m_read_bzdata,only: read_bzdata,
     &     ngrp2=>ngrp,nqbz,nqibz,nqbzw,nteti,ntetf,n1,n2,n3,qbas,ginv, !qbasmc,
     &     dq_,qbz,wbz,qibz,wibz,qbzw,
     &     idtetf,ib1bz,idteti,
     &     nstar,irk,nstbz
      use m_genallcf_v3,only: genallcf_v3,
     &     nclass,natom,nspin,nl,nn,ngrp,
     &     nlmto,nlnmx, nctot,niw, !nw_input=>nw,ef
     &     alat,delta,deltaw,esmr,symgrp,clabl,iclass,!,diw,dw
     &     invg, il, in, im, nlnm, 
     &     plat, pos, ecore, symgg , konf,z,
     &     spid
      use m_read_Worb,only: s_read_Worb, s_cal_Worb,
     &  nwf, nclass_mlwf, cbas_mlwf, nbasclass_mlwf, 
     &   classname_mlwf, iclassin,
     &  iphi, iphidot, nphi, nphix
      use m_keyvalue,only: getkeyvalue
      use m_readhbe,only:Readhbe,nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg

      implicit none
c------------------------------------
      real(8),allocatable:: r0g(:,:), wphi(:,:)
      real(8)    :: esmr2,shtw
      integer :: iclass2
      integer(4)::
     &   ixc,iopen, !ifhbed, nprecb,mrecb,mrece,nlmtot,nqbzt, nband,
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
      real(8):: qreal(3), ntot,nocctotg2,tripl!,xxx(3,3)
      real(8):: qlat(3,3)
      logical ::nocore

c     okumura Sep,2017
      integer(4):: ifdorb,iiwf
      integer(4),allocatable::idorb(:),idorb_in(:) !check
      
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


      integer(4):: bzcase,  ngpn1,verbose,ngcn1,nwxx !mrecg,
      real(8)   :: wgtq0p,quu(3)

      integer(4):: iii,isx,ivsumxxx

c for maxloc
      real(8)   :: wbb(12),wbbsum,bb(3,12),
     c             eomin,eomax,eimin,eimax,
     c             qwf0(3),dqwf0(3),qks(3),q0(3)
      complex(8),allocatable:: uumat(:,:,:,:),evecc(:,:),eveccs(:,:),
     c     amnk(:,:,:),cnk(:,:,:),umnk(:,:,:),evecc1(:,:,:),evecc2(:,:,:)
      real(8),allocatable:: ku(:,:),kbu(:,:,:),eunk(:,:),eval(:),evals(:),
     c                      eks(:),rt(:,:),rt8(:,:,:),qbz0(:,:)
      integer(4):: nbb,isc,ifq0p,
     c             nox,iko_ix,iko_fx,
     c             noxs(2),iko_ixs(2),iko_fxs(2),
     c             ieo_swt,iei_swt,itin_i,itin_f,itout_i,itout_f,
     c             nbbelow,nbabove
      integer(4),allocatable:: ikbidx(:,:)
      integer(4),allocatable:: iki_i(:),iki_f(:),
     c                         ikbi_i(:,:),ikbi_f(:,:),
     c                         iko_i(:),iko_f(:),
     c                         ikbo_i(:,:),ikbo_f(:,:)
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
      integer(4):: iband,ifbnd,iftb,ifsh,nsh,nsh1,nsh2,iffb
      logical :: lsh
      real(8),allocatable :: rws(:,:),drws(:)
      integer(4),allocatable:: irws(:)
      integer(4):: nrws,ifham

c ixc=3
      character(20)::filename
      complex(8),allocatable::  hrotrcut(:,:,:),evecc_w(:,:,:,:)
      integer:: ifh
      real(8):: heps ,r_v

      real(8)::qold(3)
      real(8),allocatable:: xq(:),eval1(:,:),eval2(:,:),eval3(:,:),eval_w(:,:,:)

      integer::npin 
      real(8):: qiin(3),qfin(3)

      integer(4),allocatable:: 
     &  m_indx(:),n_indx(:),l_indx(:),ibas_indx(:),ibasiwf(:)
      integer:: ifoc,iwf,ldim2,ixx,ifile_handle

      real(8):: enwfmax,qxx(3),eeee,enwfmaxi, ef
      integer:: inii,if102
      logical:: leauto,leinauto

cccccccccccccccccccccccccccccccccccccccccxxxxxxx
c      open(1107,file='xxx1')
c      open(1108,file='xxx2')
cccccccccccccccccccccccccccccccccccccccccccc

c---------------------------------------
      hartree=2d0*rydberg()

      iii=verbose()
      write(6,*)' verbose=',iii

c mode switch. --------------
      write(6,*) ' --- Choose omodes below ----------------'
      write(6,*) '  bb vectors (1) or Wannier fn. (2) or TB Hamiltonian (3)'
      write(6,*) ' --- Put number above ! -----------------'
      call readin5(ixc,nz,idummy)
      write(6,*) ' ixc=',ixc
      if(ixc<1.or.ixc>3) call rx(' --- ixc=0 --- Choose computational mode!')

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
      call read_Bzdata()
      write(6,*)' nqibz ngrp=',nqibz,ngrp
      write(6,*)' nqbz  =',nqbz
c      write(6,*)  qbz 
c      write(6,*)' irk=',irk
c      write(6,*)' #### idtetf: ####'
c      write(6,*) idtetf

c set up work array
c      call wkinit (iwksize)
      call pshprt(60)

C--- readin GWIN and LMTO, then allocate and set datas.
c      nwin =-999    !not readin NW file
c      efin =-999d0  !not readin EFERMI
c      efin = 0d0    !readin EFERMI
c      call readefermi()
      incwfin= -1  !use 7th colmn for core at the end section of GWIN
      call genallcf_v3(incwfin) !in module m_genallcf_v3
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

ccccccccccccccccccccccccccccccccccccccccccccc
c$$$      do i=1,natom
c$$$       print *,'  iatom, spid= ',i,spid(i)
c$$$      enddo
ccccccccccccccccccccccccccccccccccccccccccccc




c--- Get maximums takao 18June03
      call getnemx(nbmx,ebmx,8,.true.) !8+1 th line of GWIN0

c-------------------------------------------------------------------
c      if (nclass > mxclass) stop ' hsfp0: increase mxclass'
c!!!! WE ASSUME iclass(iatom)= iatom !!!!!!!!!!!!!!!!!!!!!!!!!
      if (nclass /= natom ) stop ' hsfp0: nclass /= natom ' ! We assume nclass = natom.
      write(6,*)' hsfp0: end of genallcf2'
c
      call pshprt(30)
      pi   = 4d0*datan(1d0)
      tpia = 2d0*pi/alat

c      call dinv33(plat,1,xxx,vol)
c      call dinv33(plat,1,qlat,vol)
c      voltot = dabs(vol)*(alat**3)
      call minv33tp(plat,qlat)
      voltot = abs(alat**3*tripl(plat,plat(1,2),plat(1,3)))

      ifmlw(1) = iopen('MLWU',0,-1,0)
      ifmlwe(1)= iopen('MLWEU',0,-1,0)
      if (nspin == 2) then 
         ifmlw(2) = iopen('MLWD',0,-1,0)
         ifmlwe(2)= iopen('MLWED',0,-1,0)
      endif

c>> read dimensions of wc,b,hb
c      ifhbed     = iopen('hbe.d',1,0,0)
c      read (ifhbed,*) nprecb,mrecb,mrece,nlmtot,nqbzt, nband,mrecg
c      if (nprecb == 4) stop 'hsfp0: b,hb in single precision'
      call Readhbe()
      call init_readeigen()!nband,mrece) !initialization of readEigen
      
c --- get space group information ---------------------------------
c true class information in order to determine the space group -----------
c because the class in the generated GW file is dummy.(iclass(ibas)=ibas should be kept).
      if102=ifile_handle()
      open (if102,file='CLASS')
      allocate(iclasst(natom),invgx(ngrp)
     &         ,miat(natom,ngrp),tiat(3,natom,ngrp),shtvg(3,ngrp))
      write(6,*)'  --- Readingin CLASS info ---'
      do ibas = 1,natom
        read(if102,*) ibasx, iclasst(ibas)
        write(6, "(2i10)") ibasx, iclasst(ibas)
      enddo
      close(if102)

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
        ifh=ifile_handle()
        open(ifh,file='BNDS',status='old')
        read(ifh,*)ntmp,ef
        close(ifh)
      else ! lbnds   
        call rx('you have to perform job_band in advance')
      endif
c$$$c --- determine Fermi energy ef for given valn (legas case) or corresponding charge given by z and konf.
c$$$! When esmr is negative, esmr is geven automatically by efsimplef.
c$$$        write(*,*)'Calculate EF in efsimplef2a'
c$$$        legas = .false.
c$$$        call efsimplef2a(nspin,wibz,qibz,ginv, 
c$$$     i        nband,nqibz
c$$$     i       ,konf,z,nl,natom,iclass,nclass
c$$$     i       ,valn, legas, esmr,  !!! valn is input for legas=T, output otherwise. 
c$$$c
c$$$     i        qbz,nqbz !index_qbz, n_index_qbz,
c$$$     o       ,efnew)
c$$$c
c$$$c        write(6,*)' end of efsimple'
c$$$        ef = efnew
c$$$      endif ! lbnds
c- check total ele number -------
      ntot  = nocctotg2(nspin, ef,esmr, qbz,wbz, nband,nqbz) !wbz
      write(6,*)' ef    =',ef
      write(6,*)' esmr  =',esmr
      write(6,*)' valn  =',valn
      write(6,*)' ntot  =',ntot

c      ifcphi  = iopen('CPHI',0,0,mrecb)

      call init_readeigen2()!mrecb,nlmto,mrecg) !initialize m_readeigen

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
        nline=0
        do i = 1,nlinex
          read(99,*,err=551,end=552)npin,qiin,qfin
          if (npin==0) exit
          nline = nline+1
          np(nline)=npin
          qi(1:3,nline)=qiin
          qf(1:3,nline)=qfin
 551      continue
        enddo
 552    continue
        if (nline.eq.nlinex) call rx('hmaxloc: too many lines in SYML')
        close(99)
        nq = 0
        do i = 1,nline
           nq = nq + np(i)
        enddo ! i
        allocate(q(3,nq),xq(nq))
        iq = 0
        xq=0d0
        qold=q(:,1)
        do i = 1,nline
        do j = 0,np(i)-1
          iq = iq + 1
          q(:,iq) = qi(:,i) + (qf(:,i)-qi(:,i))*dble(j)/dble(np(i)-1)
          if(iq>1) then
             !!! xq may be too long if q is outside of 1stBZ
             if(0.2 < dsqrt(sum((q(:,iq)-qold)**2))) then !!maybe same point (must be checked)
                xq(iq)= xq(iq-1)
             else
                xq(iq)= xq(iq-1) + dsqrt( sum((q(:,iq)-qold)**2) )
             endif
          endif
          qold=q(:,iq)
        enddo ! j
        enddo ! i
      else ! lsyml
        write(*,*)'Read k points for bands from GWinput'
        call getkeyvalue("GWinput","<QPNT>",unit=ifqpnt,status=ret) 
        write(6,*)' ifqpnt ret=',ifqpnt,ret
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
        allocate(xq(nq))
        xq=0d0
      endif ! syml
c
      nspinmx = nspin
      if (laf) nspinmx =1
c------------
c input parameters specific to MAXLOC
      call s_read_Worb()

      allocate(idorb(nwf)) !!okumura
      allocate(idorb_in(nwf))
      idorb=-1
      iiwf=0
      
      do iclass2=1,nclass_mlwf
        write(*,*)'output:',iclassin(iclass2), nwf
     & ,trim(classname_mlwf(iclass2)),cbas_mlwf(1:nbasclass_mlwf(iclass2),iclass2)
ccc okumura
        do iwf=1,nbasclass_mlwf(iclass2)
           iiwf=iiwf+1
           if (cbas_mlwf(iwf,iclass2) <= 1) then
              idorb(iiwf)=0     !s
           elseif (cbas_mlwf(iwf,iclass2) <= 4) then
              idorb(iiwf)=1     !p
           elseif (cbas_mlwf(iwf,iclass2) <= 9) then
              idorb(iiwf)=2     !d
           elseif (cbas_mlwf(iwf,iclass2) <= 16) then
              idorb(iiwf)=3     !f
           elseif (cbas_mlwf(iwf,iclass2) <= 25) then
              idorb(iiwf)=4     !g
           endif
           print *,"checkorb: iiwf,idorb(iiwf)",iiwf,idorb(iiwf)
        enddo

      enddo
      ifdorb=ifile_handle()
c      open(ifdorb,file="worb_list.d",form="unformatted")
      open(ifdorb,file="Worb2lorb.d",form="unformatted")
      write(ifdorb) idorb(1:nwf)
      write(ifdorb) nclass_mlwf
      write(ifdorb) nbasclass_mlwf(1:nclass_mlwf)
      close(ifdorb)

c$$$      ifdorb=ifile_handle()
c$$$      open(ifdorb,file="Worb2lorb.d",form="formatted")
c$$$      read(ifdorb,"(I4)",ADVANCE='NO') idorb_in(1:nwf)
c$$$      do iwf=1,nwf
c$$$         print *,"idorb check, iwf,idorb",iwf,idorb_in(iwf)
c$$$      enddo
c$$$      close(ifdorb)

ccc end okumura
      
      call s_cal_Worb()

      allocate (r0g(nphix,nwf), wphi(nphix,nwf))

      r0g = 2d0
      wphi = 1d0


!! Read parameters in GWinput. expand wan_input
c      call wan_input(leout,lein,lbin,ieo_swt,iei_swt,
c     &    eomin,eomax,itout_i,itout_f,nbbelow,nbabove,
c     &    eimin,eimax,itin_i,itin_f,
c     &    nsc1,nsc2,conv1,conv2,alpha1,alpha2,rcut)
c
      ieo_swt = 0
      eomin   = 0d0
      eomax   = 0d0
      itout_i = 0
      itout_f = 0
      iei_swt = 0
      eimin   = 0d0
      eimax   = 0d0
      itin_i  = 0
      itin_f  = 0
      call getkeyvalue("GWinput","wan_out_ewin",leout,default=.true.)
      call getkeyvalue("GWinput","wan_in_ewin",lein,default=.false.)
      call getkeyvalue("GWinput","wan_in_bwin",lbin,default=.false.)
      if (leout) then
         call getkeyvalue("GWinput","wan_out_emin",eomin,default=999d0 )  
         call getkeyvalue("GWinput","wan_out_emax",eomax,default=-999d0 )  
         if (eomin.gt.eomax) call rx('hmaxloc: eomin > eomax')
         ieo_swt = 1
      else
         call getkeyvalue("GWinput","wan_out_bmin",itout_i,default=999 )  
         call getkeyvalue("GWinput","wan_out_bmax",itout_f,default=-999 )  
         if (itout_i.gt.itout_f) call rx('hmaxloc: itout_i > itout_f')
      endif
      if (lein) then
         call getkeyvalue("GWinput","wan_in_emin",eimin,default=999d0 )  
         call getkeyvalue("GWinput","wan_in_emax",eimax,default=-999d0 )  
         if (eimin.gt.eimax) call rx('hmaxloc: eimin > eimax')
         iei_swt = 1
      endif
      if (lbin) then
         call getkeyvalue("GWinput","wan_in_bmin",itin_i,default=999 )  
         call getkeyvalue("GWinput","wan_in_bmax",itin_f,default=-999 )  
         if (itin_i.gt.itin_f) call rx('hmaxloc: itin_i > itin_f')
         iei_swt = 2
      endif
      call getkeyvalue("GWinput","wan_maxit_1st",nsc1,default=100)  
      call getkeyvalue("GWinput","wan_conv_1st",conv1,default=1d-5)  
      call getkeyvalue("GWinput","wan_mix_1st",alpha1,default=0.1d0)  
      call getkeyvalue("GWinput","wan_maxit_2nd",nsc2,default=100)  
      call getkeyvalue("GWinput","wan_conv_2nd",conv2,default=1d-5)  
      call getkeyvalue("GWinput","wan_mix_2nd",alpha2,default=0.1d0)  
      call getkeyvalue("GWinput","wan_tb_cut",rcut,default=1.01d0)  
      call getkeyvalue("GWinput","wan_nb_below",nbbelow,default=0)  
      call getkeyvalue("GWinput","wan_nb_above",nbabove,default=0)  
!     
      
c     skino
      r_v=rcut
      call getkeyvalue("GWinput",'wan_tbcut_rcut',heps,default=r_v) 
      call getkeyvalue("GWinput",'wan_tbcut_heps',heps,default=0.0d0) 
      write(*,*) 'mloc.heps ', heps 
c     ekino
      
      
cc --- read LDA eigenvalues
      ntq = nwf
cc     ntp0=ntq
c      allocate(eqx(ntq,nq,nspin),eqx0(ntq,nq,nspin),eqt(nband))
c      do      is = 1,nspin
c      do      ip = 1,nq
cc        iq       = idxk (q(1,ip),qbze,nqbze)
cc        call rwdd1   (ifev(is), iq, nband, eqt) !direct access read b,hb and e(q,t)
c        call readeval(q(1,ip),is,eqt)
cc        write(6,*)' eqt=',eqt
c        eqx0(1:ntq,ip,is) = eqt(itq(1:ntq))
c        eqx (1:ntq,ip,is) = rydberg()*(eqt(itq(1:ntq))- ef)
c      enddo
c      enddo
c      deallocate(eqt)

c --- info
      call winfo(6,nspin,nq,ntq,is,nbloch
     &     ,0,0,nqbz,nqibz,ef,deltaw,alat,esmr)

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
c      write(6,"(a,9f9.4)")'qbas=',qbas
c      write(6,"(a,9f9.4)")'plat=',plat
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

!! takao list eigen -----
      enwfmax =-1d9
      enwfmaxi=1d9
      allocate(eqt(1:nband))
      do is = 1,nspin
      do iq = 1,nqbz
         qxx = qbz(:,iq)
         call readeval(qxx,is,eqt)
         ini=1
         do i=1,nband
c           write(6,*)'eqeq',eqt(i),eomin,eqt(nwf)
c            if (eqt(i)>eomin) then
            if( (eqt(i)-ef)*rydberg() >eomin) then
             inii=i
             exit
           endif
         enddo
         eeee= (eqt(nwf+inii-1)-ef)*rydberg()
         write(6,"('elist: q iq is nwfi nwfe e(nwf)= ',3f9.4,i5,i2,2i5,f10.3)") qxx,iq,is,inii,nwf+inii-1, eeee
         if (enwfmax < eeee) enwfmax  = eeee
         if (enwfmaxi >eeee) enwfmaxi = eeee
      enddo   
      enddo   
      deallocate(eqt)
      write(6,"('elist max enwf enwfmaxi=',2f13.5)") enwfmax,enwfmaxi
      call getkeyvalue("GWinput","wan_out_emax_auto",leauto,default=.false.)
      if(leauto) then
         eomax= enwfmax + 1d-4
         write(6,*)
         write(6,"(' WE USE wan_out_emax_auto on ==> +1d-3 ==> eomax=',3f13.5)") eomax
      endif
      call getkeyvalue("GWinput","wan_in_emax_auto",leinauto,default=.false.)
      if(leinauto) then
         eimax= enwfmaxi + 1d-4
         write(6,*)
         write(6,"(' WE USE  wan_in_emax_auto on ==> +1d-3 ==> eimax=',3f13.5)") eimax
      endif
      
c      stop 'qqqqqqqqqqqqqqqqqq'

!! ixc = 1 ----------------
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
         call rx0('hmaxloc: ixc=1 ok')
      endif

!! loop over spin -----------------------
      do 1000 is = 1,nspin
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

ccccccccccccccccccccccc
      do iq=1,nqbz
      write(6,"('iiii: iq and internal window region=',3i5)")iq,iki_i(iq),iki_f(iq)
      enddo
cccccccccccccccccccccc

c read uu-matrix
      allocate (uumat(iko_ix:iko_fx,iko_ix:iko_fx,nbb,nqbz))
      call readuu(is,iko_ix,iko_fx,ikbidx,
     d            nqbz,nbb,
     o            uumat)
      call chkuu(is,iko_ix,iko_fx,ikbidx,uumat,
     d           nqbz,nbb)

!! step 1  -- choose Hilbert space -- determine cnk
      write(*,*)'Step 1: Hilbert space branch'
      write(6,*)' iko_ix iko_fx=',iko_ix,iko_fx
      allocate (amnk(iko_ix:iko_fx,nwf,nqbz),
     &          upu(iko_ix:iko_fx,iko_ix:iko_fx,nbb,nqbz),
     &          cnk(iko_ix:iko_fx,nwf,nqbz),
     &          cnk2(iko_ix:iko_fx,nwf,nqbz),
     &          omgik(nqbz))
! amnk appered in Eq.22 in Ref.II. <psi|Gaussian>
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
              if (ndz.lt.1) call rx('ndz < 1')
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
         write(*,"('#SC-loop, conv.',i5,d13.5)")isc,omgi 
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

c      call chk_cnkweight(qbz,iko_ix,iko_fx,cnk,
c     &     nqbz,nwf,nband,nlmto)

!! NOTE: cnk is the final results of step 1
!!   cnk(iko_ix:iko_fx,nwf,nqbz)
!!   cnk(iko_i(iq):iko_f(iq),nwf,iq) gives nwf-dimentional space. 
!!   step 1 (minimization of Omega_I)


!! === step 2 -- localize Wannier fn. ================
      write(*,*)'Step 2: Wannier fn. branch'

      allocate (mmn(nwf,nwf,nbb,nqbz),mmn0(nwf,nwf,nbb,nqbz),
     &          umnk(nwf,nwf,nqbz),
     &          rmn(nwf,nwf),amn(nwf,nwf),smn(nwf,nwf),
     &          rn(3,nwf),qn(nwf),tmn(nwf,nwf),dwmn(nwf,nwf),
     &          eunk(nwf,nqbz))

!! (2-0) construct initlal u~ from u
!! eunk(= e~) of {H~}_mn): eigenvalue within the nwf-dimentional Hilbert space
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

!! (2-1) initial: uumat -> M_mn(0) Eq.58 in Ref.[1]
      call init_mmn(cnk,uumat,ikbidx,
     i              iko_ix,iko_fx,iko_i,iko_f,ikbo_i,ikbo_f,
     d              nwf,nqbz,nbb,
     o              mmn0)

!! (2-2) initial U
!!    umnk= U(m,n) = ( A S^{-1/2} )_mn. See Eq.23 in Ref.II.
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

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c         do i=1,nbb
c            write(1106+is,"(a,i4,13f13.5)")'bbbb',i,bb(1:3,i),wbb(i)
c         enddo
c         do i=1,nqbz
c            write(1106+is,"(a,i4,13f13.5)")'www',i,wbz(i)
c         enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do isc = 1, nsc2
cccccccccccccccccccccccc
c      mmn=mmn+(0d0,1d-8)
ccccccccccccccccccccccc

c <r_n> ([1] eq.31)
         call get_rn(mmn,bb,wbb,wbz,
     d               nwf,nqbz,nbb,
     o               rn)
ccccccccccccccccccccccccccccccccccccccccccccc
c         do i=1,nwf
c            write(1106+is,"(a,3f13.5)")'rrrrrn',rn(1:3,i)
c         enddo
ccccccccccccccccccccccccccccccccccccccccccccccccc

         do iq = 1,nqbz
            dwmn = (0d0,0d0)
            do ibb = 1,nbb
cccccccccccccccccccccccccccccccccccccccccccccc
c      do i = 1,nwf
c      do j = 1,nwf
c         write(1106+is,"(a,4i5,2f13.3)")' mmmmm ',i,j,ibb,iq,mmn(i,j,ibb,iq)+(0d0,0.0001)
c      enddo
c      enddo
c      do i = 1,nwf
c      do j = 1,nwf
c         write(1106+is,"(a,4i5,2f13.3)")' nnnnnn ',i,j,ibb,iq,mmn0(i,j,ibb,iq)+(0d0,0.0001)
c      enddo
c      enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccc

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

cccccccccccccccccccccccccccccccccc
c               smn=0d0
c               amn=0d0
cccccccccccccccccccccccccccccccccc

c DW(k) ([1] eq.57)
               dwmn(:,:) =  dwmn(:,:)
     &       + wbb(ibb) * (amn(:,:) - smn(:,:)) * alpha2 / wbbsum

c end of ibb-loop
            enddo

ccccccccccccccccccccccccc
c            dwmn=0d0
ccccccccccccccccccccccccc
c (2-5) DW(k) -> U(k) ([1] eq.60)
            call updt_uk(dwmn,
     d                   nwf,
     u                   umnk(:,:,iq))
c            call chk_um(umnk(:,:,iq),nwf)

ccccccccccccccccccccccccccccccccccccccccccccccccccc
c      do i = 1,nwf
c      do j = 1,nwf
c         write(1106+is,"(a,3i5,2f13.3)")' zzzzz',iq,i,j,umnk(i,j,iq)
c      enddo
c      enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccc

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
c         write(*,*)'#SC-loop, conv.',isc,omgdod
c         write(*,950)'Omg: I, OD, D ',omgi,omgod,omgd
         write(*,"('#SC-loop, conv.',i6,e13.5,' Omg:_I,_OD,_D= ',3f17.10)")
     &   isc,omgdod,omgi,omgod,omgd
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
      write(*,*)"---------- wlaxloc isp =",is
      call wmaxloc(ifmlw(is),ifmlwe(is),
     i             qbz,umnk,cnk,eunk,
     i             iko_ix,iko_fx,iko_i,iko_f,
     d             nwf,nqbz,nband,nlmto, is)
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


!! step 3 -- reduced Hamiltonian ------------------------------
      write(*,*)'Step 3: reduced Hamiltonian branch'
c open file
      if (is .eq. 1) then
         ifbnd = iopen('bnds.maxloc.up',1,-1,0)         
         iftb  = iopen('bnds.tb.up',1,-1,0)         
         iffb  = iopen('bnds.fb.up',1,-1,0)
      else
         ifbnd = iopen('bnds.maxloc.dn',1,-1,0)         
         iftb  = iopen('bnds.tb.dn',1,-1,0)
         iffb  = iopen('bnds.fb.dn',1,-1,0)         
      endif
      write(ifbnd,*)nq
      write(ifbnd,*)nwf
      write(iftb,*)nq
      write(iftb,*)nwf
      write(iffb,*)'#',nq
      write(iffb,*)'#',nwf
c allocate
      if(allocated(hrotk)) deallocate(hrotk,hrotkp,evecc,eval)
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
c (3-1) ~H(k) -> Hrot(k): note eunk is eigenvalues in the basis of cnk
      call rot_hmnk(umnk,eunk,
     d              nwf,nqbz,
     o              hrotk) !rotated Hamiltonian in MLW basis. 
c (3-2) Hrot_mn(R)
      if(allocated(irws)) deallocate(irws,rws,drws)
      allocate(irws(n1*n2*n3*8),rws(3,n1*n2*n3*8),drws(n1*n2*n3*8))
      call wigner_seitz(alat,plat,n1,n2,n3,nrws,rws,irws,drws)
      if(allocated(hrotr)) deallocate(hrotr)
      allocate(hrotr(nwf,nwf,nrws)) !real space Hamiltonian in Wannier funciton basis
c      write(*,*) 'xxxxxxxxxx1'
      if (ixc.eq.2) then
         call get_hrotr_ws(hrotk,qbz,wbz,
     i        rws,irws,drws,
     d        nwf,nqbz,nrws,
     o        hrotr)
c     skino
c     write hrotr and *rws 
         if (is .eq. 1) then
            ifh = iopen('hrotr.up',1,-1,0)
         else
            ifh = iopen('hrotr.dn',1,-1,0)
         endif

         call write_hrotr(ifh, hrotr,
     i        rws,irws,drws,
     d        nwf,nrws )

         close (ifh) 
c     ekino
c     skino
      else if (ixc.eq.3) then
         if (is .eq. 1) then
            filename='hrotr.up'
         else
            filename = 'hrotr.dn'
         endif
         call read_hrotr(filename,nwf,nrws,
     o        hrotr) 
         if (is .eq. 1) then
            ifh = iopen('hrotr.cut.up',1,-1,0)
         else
            ifh = iopen('hrotr.cut.dn',1,-1,0)
         endif
         allocate(hrotrcut(nwf,nwf,nrws))
         call make_hrotrcut( hrotr,
     i        rws,irws,drws,
     i        rcut,heps,
     d        nwf,nrws,
     o        hrotrcut )
         call write_hrotr(ifh, hrotrcut,
     i        rws,irws,drws,
     d        nwf,nrws )
         close (ifh)
         deallocate(hrotrcut)
c     ekino
      endif
c      write(*,*) 'xxxxxxxxxx2'

!! ----------------------------------------------------------
!! k-point mesh
      call get_nqbze(nqbz,nqbze)
      allocate(qbze(3,nqbze))
      call get_qbze(qbz,nqbz,
     o              qbze,nqbze)
      write(ifmlw(is))nqbze,nwf
      write(ifmlwe(is))nqbze,nwf
      do iq = 1,nqbze
c         write(*,*)'goto get_hrotkp_ws iq=',iq,nqbze
          call get_hrotkp_ws(hrotr,rws,drws,irws,qbze(:,iq), !july2014 qbz->qbze
     d                     nwf,nqbz,nrws,
     o                     hrotkp)
          call diag_hm(hrotkp,nwf,eval,evecc)
          call wmaxloc_diag(ifmlw(is),ifmlwe(is),
     i                 iq,qbze(1:3,iq),umnk,cnk,eunk,evecc,eval,
     i                 iko_ix,iko_fx,iko_i,iko_f,
     d                 nwf,nqbz)
          write(6,*)"iq is",iq,is
       enddo
c      write(6,*)'eeeeeeeee'
      deallocate(qbze)
ccc          write(*,990)'iq =',iq,qbz(1:3,iq)
cc          if (iq.le.nqbz) then
cc          do iband = 1,nwf
cc             e1 = (eval(iband)   -ef)*rydberg()
cc             e2 = (eunk(iband,iq)-ef)*rydberg()


!! --------------------------------------------------------------
c --- Readin nlam index
      ifoc = iopen('@MNLA_CPHI',1,0,0)
      ldim2 = nlmto
      read(ifoc,*)
      if(allocated(m_indx)) deallocate(m_indx,n_indx,l_indx,ibas_indx,ibasiwf)
      allocate(m_indx(ldim2),n_indx(ldim2),l_indx(ldim2),ibas_indx(ldim2))
      do ix =1,ldim2
        read(ifoc,*)m_indx(ix),n_indx(ix),l_indx(ix),ibas_indx(ix),ixx
        if(ixx/=ix) call rx('failed to readin @MNLA_CPHI')
      enddo   
      ix = iclose('@MNLA_CPHI')
      allocate(ibasiwf(nwf))
      do iwf=1,nwf
        ibasiwf(iwf) = ibas_indx(iphi(1,iwf))
      enddo


!! write HrotRS
      ifh=ifile_handle()
      if(is==1) open(ifh,file='HrotRS.up',form='unformatted')
      if(is==2) open(ifh,file='HrotRS.dn',form='unformatted')
      write(ifh)alat,plat,natom
      write(ifh)pos
      write(ifh)ef 
      write(ifh)nwf,nrws,n1,n2,n3
      write(ifh) irws,rws,drws,hrotr, ibasiwf !drws added by okumura Aug28,2017
      close(ifh)

      ifh = ifile_handle()

      call write_hopping_output(is, ifh, hrotr,
     &                  rws,irws,alat,plat,qlat,pos,natom,
     &                  ibasiwf, nwf,nrws,spid , m_indx, l_indx,
     &                  nphix, iphi, ldim2)

      close(ifh)

!! TEST okumura: iq=1,nq, q->qbz?   (2017/06/10)
      !!! qtt -> q, nqtt -> nqbz     
      ! data list for wannier
      ifh=ifile_handle()
      open(ifh,file="wan4chi.d",form="unformatted")
      write(ifh) nwf,nspin,nqbz
      close(ifh)

      ! generate eigenvalue and eigenvector of Wannier Hamiltonian
      ! Index:: evecc_w (orbital,band,q-point,spin) 
      write(6,*) 
      if (is==1) allocate(eval_w(nwf,nqbz,nspin),evecc_w(nwf,nwf,nqbz,nspin))
      do iq = 1,nqbz
         write(6,*)' got get_hrotkp_ws iq =',iq
         call get_hrotkp_ws(hrotr,rws,drws,irws,qbz(:,iq),
     d                     nwf,nqbz,nrws,
     o                     hrotkp)
         call diag_hm(hrotkp,nwf,eval,evecc)
         eval_w(1:nwf,iq,is)=eval
         evecc_w(1:nwf,1:nwf,iq,is)=evecc
      enddo

      if(is==2) then
         ifh=ifile_handle()
         open(ifh,file='EValue_w',form='unformatted')
         write(ifh) nwf,nqbz,nspin
         write(ifh) eval_w(1:nwf,1:nqbz,1:nspin)
         close(ifh)
         ifh=ifile_handle()
         open(ifh,file='EVec_w',form='unformatted')
         write(ifh) nwf,nqbz,nspin
         write(ifh) qbz(1:3,1:nqbz)
         write(ifh) evecc_w(1:nwf,1:nwf,1:nqbz,1:nspin)
         close(ifh)
      endif
     !! end okumura
      
!! other k-points
      write(ifbnd,*)ef,' ef'
      write(iftb,*)ef,' ef'
      write(iffb,*)'#',ef,' ef'
      if (lsh) write(ifsh,*)ef,' ef'
      allocate(eval1(nwf,nq),eval3(nwf,nq),evecc1(nwf,nwf,nq))
      if(lsh) allocate(eval2(nwf,nq),evecc2(nwf,nwf,nq))

      do iq = 1,nq
c     write(6,*)' got get_hrotkp_ws iq =',iq
c (3-3) Hrot_mn(k')
          call get_hrotkp_ws(hrotr,rws,drws,irws,q(:,iq),
     d                     nwf,nqbz,nrws,
     o                     hrotkp)
c (3-4) diagonalize
          call diag_hm(hrotkp,nwf,eval,evecc)
          eval1(1:nwf,iq)=eval
          evecc1(1:nwf,1:nwf,iq)=evecc
c     (3-4) diagonalize  -- Small Hamiltonian --
          if (lsh) then
             hrotkps(1:nsh,1:nsh) = hrotkp(nsh1:nsh2,nsh1:nsh2)
             call diag_hm(hrotkps,nsh,evals,eveccs)
             write(ifsh,*)'iq =',iq
             write(ifsh,990)q(1:3,iq)
             eval2(1:nsh,iq)= evals(1:nsh)
             evecc2(1:nwf,1:nwf,iq)=eveccs
          endif                 ! lsh
c (3-3) Hrot_mn(k')  -- Tight-binding ---
          call get_hrotkp_tb_ws(rcut,plat,alat,
     i         hrotr,rws,drws,irws,q(:,iq),  ibasiwf,pos,natom,
     d         nwf,nqbz,nrws,
     o         hrotkp)

c     (3-4) diagonalize -- Tight-binding --
          call diag_hm(hrotkp,nwf,eval,evecc)
          eval3(1:nwf,iq)=eval

       enddo
       do iband = 1,nwf
         do iq = 1,nq
            write(ifbnd,"(i5,3f13.5,'  ',f13.6,f13.6,i5,' !eee! x eval-ef(ev) iband' )")
     &           iq,q(1:3,iq),  xq(iq),(eval1(iband,iq)-ef)*rydberg(),iband
            write(iftb,"(i5,3f13.5,'  ',f13.6,f13.6,i5,' !eee! x eval-ef(ev) iband' )")
     &           iq,q(1:3,iq),  xq(iq),(eval3(iband,iq)-ef)*rydberg(),iband

            write(iffb,"(i5,3f13.5,'  ',f13.6,f13.6,i5,' ')",ADVANCE='NO') 
     &           iq,q(1:3,iq),  xq(iq),(eval1(iband,iq)-ef)*rydberg(),iband
            do iwf=1,nwf
               write(iffb,"(f13.6)",ADVANCE='NO')
     &              (abs(evecc1(iwf,iband,iq)))**2 
            enddo
            write(iffb,*)            
         enddo
         write(ifbnd,*)
         write(iftb,*)
         write(iffb,*)
         write(iffb,*)
      enddo

      deallocate(eval1,eval3,evecc1)

      if(lsh) then
         do iband = 1,nsh
            do iq = 1,nq
               write(ifsh,"(i5,3f13.5,'  ',f13.6,f13.6,i5,' !eee! x eval-ef(ev) iband' )") 
     &              iq,q(1:3,iq),  xq(iq),(eval2(iband,iq)-ef)*rydberg(),iband
            enddo
         enddo
      endif
      call writeham(ifham,is,ef,alat,plat,pos,qbz,wbz,rws,irws,hrotk,nspin,natom,nwf,nqbz,nrws)
      deallocate(cnk,umnk,eunk,hrotk,hrotr,hrotkp,evecc,eval,irws,rws,drws, 
     &           ibasiwf,m_indx,n_indx,l_indx,ibas_indx)
      if (lsh) deallocate(hrotkps,eveccs,evals,evecc2)
      close(ifbnd)
      close(iftb)
      close(iffb)
c     end of loop over spin
 1000 continue
 950  format(a14,3f23.16)
 990  format(3f12.6)
      call cputid(0)
      call rx0('hmaxloc: ixc=2 ok')
      end

c$$$c-----------------------------------------------------------------------
c$$$      subroutine chk_amnkweight(qbz,iko_ix,iko_fx,amnk,
c$$$     &     nqbz,nwf,nband,nlmto)
c$$$      use m_readqg
c$$$      use m_readeigen
c$$$      implicit real*8(a-h,o-z)
c$$$
c$$$      complex(8) :: amnk(iko_ix:iko_fx,nwf,nqbz)
c$$$      complex(8),allocatable:: cphi1(:,:),cphi2(:,:)
c$$$      real(8) :: qbz(3,nqbz),q(3),quu(3)
c$$$      real(8),allocatable:: wbas(:,:)
c$$$      integer(4),allocatable:: 
c$$$     &  m_indx(:),n_indx(:),l_indx(:),ibas_indx(:)
c$$$
c$$$c --- Readin nlam index
c$$$      ifoc = iopen('@MNLA_CPHI',1,0,0)
c$$$      ldim2 = nlmto
c$$$      read(ifoc,*)
c$$$      allocate(m_indx(ldim2),n_indx(ldim2),l_indx(ldim2),ibas_indx(ldim2))
c$$$      do ix =1,ldim2
c$$$        read(ifoc,*)m_indx(ix),n_indx(ix),l_indx(ix),ibas_indx(ix),ixx
c$$$        if(ixx/=ix) call rx('failed to readin @MNLA_CPHI')
c$$$      enddo   
c$$$
c$$$      nbas = ibas_indx(nlmto)
c$$$      allocate(cphi1(nlmto,nband),cphi2(nlmto,nwf),wbas(nbas,nwf))
c$$$      wbas = 0d0
c$$$      cphi2=0d0
c$$$      do iq = 1,nqbz
c$$$      q = qbz(:,iq)
c$$$      call readcphi(q,nlmto,1,quu,cphi1)
c$$$
c$$$       do iwf=1,nwf
c$$$         do ib=iko_ix,iko_fx
c$$$            cphi2(:,iwf) = cphi2(:,iwf) + cphi1(:,ib)*amnk(ib,iwf,iq)
c$$$         enddo
c$$$      enddo
c$$$
c$$$      enddo ! iq
c$$$
c$$$      do iwf=1,nwf
c$$$         do ia=1,nlmto
c$$$            ibas = ibas_indx(ia)
c$$$            wbas(ibas,iwf) = wbas(ibas,iwf) + 
c$$$     &                 conjg(cphi2(ia,iwf))*cphi2(ia,iwf)
c$$$         enddo ! ia
c$$$      enddo ! iwf
c$$$      wbas = wbas / dble(nqbz**2)
c$$$
c$$$      write(*,*)'*** ibas,iwf,wbas'
c$$$      do iwf=1,nwf
c$$$      do ibas=1,nbas
c$$$         write(*,*)ibas,iwf,wbas(ibas,iwf)
c$$$      enddo
c$$$      write(*,*)
c$$$      enddo   
c$$$      write(*,*)'*** ibas,wbas'
c$$$      do ibas=1,nbas
c$$$         w = 0d0
c$$$         do iwf=1,nwf
c$$$            w = w + wbas(ibas,iwf)
c$$$         enddo
c$$$         write(*,*)ibas,w
c$$$      enddo
c$$$   
c$$$      deallocate(cphi1,cphi2,wbas,m_indx,l_indx,n_indx,ibas_indx)
c$$$      ix = iclose('@MNLA_CPHI')
c$$$
c$$$      end
c$$$c-----------------------------------------------------------------------
c$$$      subroutine chk_cnkweight(qbz,iko_ix,iko_fx,cnk,
c$$$     &     nqbz,nwf,nband,nlmto)
c$$$      use m_readqg
c$$$      use m_readeigen
c$$$      implicit real*8(a-h,o-z)
c$$$
c$$$      complex(8) :: cnk(iko_ix:iko_fx,nwf,nqbz)
c$$$      complex(8),allocatable:: cphi1(:,:),cphi2(:,:)
c$$$      real(8) :: qbz(3,nqbz),q(3),quu(3)
c$$$      real(8),allocatable:: wbas(:,:)
c$$$      integer(4),allocatable:: 
c$$$     &  m_indx(:),n_indx(:),l_indx(:),ibas_indx(:)
c$$$
c$$$c --- Readin nlam index
c$$$      ifoc = iopen('@MNLA_CPHI',1,0,0)
c$$$      ldim2 = nlmto
c$$$      read(ifoc,*)
c$$$      allocate(m_indx(ldim2),n_indx(ldim2),l_indx(ldim2),ibas_indx(ldim2))
c$$$      do ix =1,ldim2
c$$$        read(ifoc,*)m_indx(ix),n_indx(ix),l_indx(ix),ibas_indx(ix),ixx
c$$$        if(ixx/=ix) call rx('failed to readin @MNLA_CPHI')
c$$$      enddo   
c$$$
c$$$      nbas = ibas_indx(nlmto)
c$$$      allocate(cphi1(nlmto,nband),cphi2(nlmto,nwf),wbas(nbas,nwf))
c$$$      wbas = 0d0
c$$$      cphi2=0d0
c$$$
c$$$      do iq = 1,nqbz
c$$$      q = qbz(:,iq)
c$$$      call readcphi(q,nlmto,1,quu,cphi1)
c$$$
c$$$      do iwf=1,nwf
c$$$         do ib=iko_ix,iko_fx
c$$$            cphi2(:,iwf) = cphi2(:,iwf) + cphi1(:,ib)*cnk(ib,iwf,iq)
c$$$         enddo
c$$$      enddo
c$$$
c$$$      enddo ! iq
c$$$
c$$$      do iwf=1,nwf
c$$$         do ia=1,nlmto
c$$$            ibas = ibas_indx(ia)
c$$$            wbas(ibas,iwf) = wbas(ibas,iwf) + 
c$$$     &                 conjg(cphi2(ia,iwf))*cphi2(ia,iwf)
c$$$         enddo ! ia
c$$$      enddo ! iwf
c$$$      wbas = wbas / dble(nqbz*nqbz)
c$$$
c$$$      write(*,*)'*** ibas,iwf,wbas'
c$$$      do iwf=1,nwf
c$$$      do ibas=1,nbas
c$$$         write(*,*)ibas,iwf,wbas(ibas,iwf)
c$$$      enddo
c$$$      write(*,*) 
c$$$      enddo   
c$$$
c$$$      write(*,*)'*** ibas,wbas'
c$$$      do ibas=1,nbas
c$$$         w = 0d0
c$$$         do iwf=1,nwf
c$$$            w = w + wbas(ibas,iwf)
c$$$         enddo
c$$$         write(*,*)ibas,w
c$$$      enddo
c$$$      deallocate(cphi1,cphi2,wbas,m_indx,l_indx,n_indx,ibas_indx)
c$$$      ix = iclose('@MNLA_CPHI')
c$$$      end
c-----------------------------------------------------------------------
c$$$      subroutine chk_umn(cnk,umnk,qbz,
c$$$     i                  iko_ix,iko_fx,iko_i,iko_f,
c$$$     d                  nwf,nqbz,nband,nlmto)
c$$$      use m_readqg
c$$$      use m_readeigen
c$$$
c$$$      implicit real*8(a-h,o-z)
c$$$
c$$$      complex(8) :: cnk(iko_ix:iko_fx,nwf,nqbz),
c$$$     &              dnk(iko_ix:iko_fx,nwf,nqbz),
c$$$     &              umnk(nwf,nwf,nqbz)
c$$$      real(8) :: qbz(3,nqbz)
c$$$      integer(4) :: iko_i(nqbz),iko_f(nqbz)
c$$$
c$$$      dnk = (0d0,0d0)
c$$$      do iq = 1,nqbz
c$$$         do imp = iko_i(iq),iko_f(iq)
c$$$         do in = 1,nwf
c$$$            do im = 1,nwf
c$$$               dnk(imp,in,iq) = dnk(imp,in,iq) 
c$$$     &                 + umnk(im,in,iq) * cnk(imp,im,iq)
c$$$            enddo ! im
c$$$         enddo ! in
c$$$         enddo ! imp      
c$$$      enddo ! iq  
c$$$
c$$$      call chk_cnkweight(qbz,iko_ix,iko_fx,dnk,
c$$$     &     nqbz,nwf,nband,nlmto)
c$$$
c$$$      end
