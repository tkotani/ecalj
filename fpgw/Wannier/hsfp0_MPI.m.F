      program hsfp0_v2_MPI
c-------------------------------------------------------------
c Calculates the  self-energy \Sigma within the
c GW approximation at the LDA energy and two neighbouring points:
c
c SEx(q,t,t) = <psi(q,t) |SEx| psi(q,t)>  
c SEc(q,t,t) = <psi(q,t) |SEc| psi(q,t)>
c 
c SEc(r,r';w) = (i/2pi) < [w'=-inf,inf] G(r,r';w+w') Wc(r,r';w') >
c See details in sxcf_v2.f
c
c This routine requirs an input from standard IO.
c In scripts, you can do it like, prompt>echo mode|../exec/hsfp0 >lsc,
c where mode is 1 or 2.
c
c  mode= 1: exchange    mode SEx, the exchange part of the self-energy
c  mode= 2: correlation mode SEc, the correlated part of the self-energy
c  mode= 3: core exchange mode SEXcore
c  mode= 4: plot spectrum function ---See manual
c  mode= 5: exchange    mode Exx, the exchange part of the total energy
c  mode= 6: core exchange mode Exxc, Exx(core-valence)
c
c Oct 2008 Takashi Miyake, from hsfp0_RSMPI.m.F
c May 2002 Takashi Miyake. Total energy calc.
c Apr 2002 takao kotani. multiple argumentation wave per l. 
c This hsfp0 is build from hsec10.f by F.Aryasetiawan.
c------------------------------------------------------------
c
c    May 2007 MPI version by Rei Sakuma (RS)
c      Parallelized with respect to q-points.
c      Summation over rotations and k-points (in irreducible k-points)
c      is also parallelized
c      if number of processes are greater than that of q-points.
c      See gwsrc/RSMPI_rotkindex_mod.F, gwsrc/RSMPI_qkgroup_mod.F
c------------------------------------------------------------
      use m_readqg
      use m_readeigen
      use m_read_bzdata,ngrp2=>ngrp
      use m_genallcf_v3
      use keyvalue
c RS: MPI module
      use rsmpi
      use rsmpi_rotkindex
      implicit none
c      real(8),parameter :: 
c     &  ua    = 1d0    ! constant in w(0)exp(-ua^2*w'^2) to take care of peak around w'=0 
c------------------------------------
c test switches to calculate the self-energy based on an another separation of \Sigma.
! \Sigma = \Sigma_{sx} + \Sigma_{coh} + \Sigma_{img axis} + \Sigma_{pole} by Hedin PR(1965)A785
! I found COH term has inevitably poor accuracy.
      logical ::tetra, tetra_hsfp0,         
     &  screen = .false., ! \Sigma_{sx} for mode 1 and 
                            ! \Sigma_{img axis} + \Sigma_{pole} for mode 2
     &   cohtest= .false.   ! \Sigma_{coh}. mode swich is not required.
     & , tetraex  = .false. ! This switch is only meaningful for mode=1,5,6
                            ! If you want to calculate exchange, use tetraex=T .
                            ! Note that you have to supply EFERMI by the tetrahedon method.
c------------------------------------
      real(8)    :: esmr2,shtw
c      integer(4) :: mxclass,ngnmax,mbytes,mwords,iwksize,
c     &   natom,nclass,ipos,igrp,
c     &   iqibz,
c     &   iqbz,
c     &   iinvg,
c     o   nspin,nl,nn,nnv,nnc,
c     o   inindx,inindxv,inindxc,iiclass,
c     d   nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc,
c     o   iz, 
c     o   iil,iin,iim,iilnm,i_mnl,
c     o   iilv,iinv,iimv,iilnmv,i_mnlv,
c     o   iilc,iinc,iimc,iilnmc,i_mnlc,
c     o   incwf,iecore,ikonf,iicore,incore,nctot,
c     o   imagw,niw,nw,ifreq,  
      integer(4)::
     &   ixc,iopen,ifhbed, nprecb,mrecb,mrece,nlmtot,nqbzt, nband,
     &   ibas,ibasx,ngpmx,nxx,ngcmx,nbloch,ifqpnt,ifwd,
     &   nprecx,mrecl,nblochpmx2,nwp,niwt, nqnum,mdimx,nblochpmx,
     &   ifrcw,ifrcwi,  noccxv,maxocc2,noccx,ifvcfpout,iqall,iaf,ntq,
     &   i,k,nspinmx, nq,is,ip,iq,idxk,ifoutsex,iclose,nq0i,ig,
     &   mxkp,nqibzxx,ntet,nene,iqi, ix,iw,
     &   nlnx4,niwx,irot,invr,invrot,ivsum, ifoutsec,ntqx,
c
c     &   ifrb(2),ifcb(2),ifrhb(2),ifchb(2)
c     &    ifev(2),
     &    ifsec(2) !,ifcphi
     &   ,ifxc(2),ifsex(2), ifphiv(2),ifphic(2),ifec,ifexsp(2),
     &   ifsecomg(2),ifexx,ndble=8
c
c      real(8) :: alat,ef,diw,dw,delta,pi,tpia,vol,voltot,rs,alpha,
      real(8) :: pi,tpia,vol,voltot,rs,alpha,
     & qfermi,efx,valn,efnew,edummy,efz,qm,xsex,egex,
     & zfac1,zfac2,dscdw1,dscdw2,dscdw,zfac,ef2=1d99,exx,exxq,exxelgas
c
c     &   lowesteval !defined in readeigen
ccc   qbas(3,3),ginv(3,3) plat(3,3),
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
      integer(4),allocatable :: ngvecpB(:,:,:),!ngveccB(:,:,:),
     &  ngvecp(:,:), ngvecc(:,:),iqib(:), !,ngpn(:)ngcni(:)
     &  kount(:,:), nx(:,:),nblocha(:),lx(:) !ngveccBr(:,:,:)
      real(8),allocatable:: vxcfp(:,:,:),
     &  wqt(:), wgt0(:,:),q0i(:,:),
     &  ppbrd (:,:,:,:,:,:,:),cgr(:,:,:,:),eqt(:),
     &  ppbrdx(:,:,:,:,:,:,:),aaa(:,:), !symope(:,:,:)=symgg, ! qibz(:,:),
     &  ppb(:), eq(:), !,pdb(:),dpb(:),ddb(:)
     &  eqx(:,:,:),eqx0(:,:,:),ekc(:),coh(:,:)
      complex(8),allocatable:: geigB(:,:,:,:) ,zsec(:,:,:)
c 
      logical :: exchange, legas, tote 
      real(8) ::  rydberg,hartree 
      real(8):: qreal(3), ntot,nocctotg2,tripl,xxx(3,3)
      logical ::nocore

c space group infermation 
      integer(4),allocatable :: iclasst(:), invgx(:), miat(:,:)
      real(8),allocatable    :: tiat(:,:,:),shtvg(:,:)

c tetra 
      real(8),allocatable :: qz(:,:),qbzxx(:),wbzxx(:),wtet(:,:,:,:),
     &     eband(:,:,:), ene(:) !,ecore(:,:)
      integer(4),allocatable ::idtetx(:,:),idtet(:,:),ipq(:)
     &      ,iene(:,:,:),ibzx(:) ! ,nstar(:)
c      real(8) :: qbasmc(3,3)

c worksize in megabytes (1 word = 4 bytes)
c      integer(4) :: mbytes,mwords,iwksize
c      parameter (mbytes=60)
c      parameter (mwords=mbytes/4)
c      parameter (iwksize=mwords * 1000 *1000)
c      integer w
c      common /w/ w(iwksize)

      integer(4) ::ib,iqx,igp,iii,ivsumxxx,isx,iflegas, iqpntnum
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

      integer(4) :: iwini,iwend,itini,itfin
      real(8),allocatable:: omega(:,:)
      real(8) ::  omegamax,dwplot,omegamaxin,eqtmin,eqtmax,etmp
c      logical :: sergeys

      integer(4)::nqbze,ini,nq0it,idummy
      real(8),allocatable:: qbze(:,:)

      real(8)   :: ebmx(2)
      integer(4):: nbmx(2)
      
      real(8):: volwgt

      integer(4)::nwin, incwfin
      real(8)::efin,ddw
      integer(4),allocatable::imdim(:)
      real(8),allocatable::freqx(:),freqw(:),wwx(:) !,expa(:)

      logical:: GaussSmear !readgwinput, 
      integer(4)::ret
      character*(150):: ddd


      integer(4):: bzcase,  ngpn1,mrecg,verbose,ngcn1,nwxx
      real(8)   :: wgtq0p,quu(3)

      real(8),allocatable:: freq_r(:)

      logical ::smbasis
      integer(4):: ifpomat,nkpo,nnmx,nomx,ikpo,nn_,no,nss(2)
      real(8):: q_r(3)
      real(8),allocatable:: qrr(:,:)
      integer(4),allocatable:: nnr(:),nor(:)
      complex(8),allocatable:: pomatr(:,:,:),pomat(:,:)

      real(8)::sumimg
      logical allq0i                                             !S.F.Jan06

      integer(4):: nw_i

ctm, 070130
      integer(4):: ifdk
      real(8):: dk(3)
      logical:: ldk
      real(8),allocatable :: qtmp(:,:)
ctm, 070621
      integer(4):: ifv_xc,ifevec,ifvxc2,ifse2,mrecs,ndim,ndimh,ndimh2
     &            ,nqqq,nqqqx,nspin2,ndummy,itp,itpp,itp2,it1,it2
     &            ,ii,ij
      real(8):: qqqx(3)
      real(8),allocatable :: qqq(:,:)
      complex(8), allocatable:: zsec2(:,:,:,:)
     &           ,v_xc(:,:,:,:),evec(:,:,:,:),vxc(:,:,:,:),vxc2(:,:)
      logical:: loffd

ctm, 081020
      logical:: lwang ! Wannier gauge
      integer(4) :: nwf,iko_ix,iko_fx
      complex(8),allocatable::umn(:,:)
      real(8):: rtmp

c RS: variables for MPI
      integer(4) :: input3(3),irot_local,ip_local
      complex(8),allocatable:: zsec_sum(:,:,:),zsec2_sum(:,:,:,:)

c RS: program begins here
      call RSMPI_Init()

c mode switch. --------------
      if(cohtest) then
c RS: disabled
c        screen = .true.
c        ixc = 2; nz=0
c        open(671,file='COH')
        call RSMPI_Stop("cohtest==T is not implemented!")
      else
        if (Is_IO_Root_RSMPI()) then
          write(6,*) ' --- Choose omodes below ----------------'
          write(6,*) '  Sx(1) Sc(2) ScoreX(3) Spectrum(4) ' 
          write(6,*) '  EXX_val-val (5)  Exx_core-val(6) '
          write(6,*) '  Sx_sf(11) Sc_sf(12) '
          write(6,*) '  [option --- (+ QPNT.{number} ?)] '
          write(6,*) ' --- Put number above ! -----------------'
          call readin5(ixc,nz,idummy)
          input3(1)=ixc
          input3(2)=nz
          input3(3)=idummy
        endif
      endif

      call MPI_Bcast(input3,3,MPI_INTEGER,io_root_rsmpi,
     &     MPI_COMM_WORLD,ierror_rsmpi)
      call RSMPI_Check("MPI_Bcast(input3)",ierror_rsmpi)
      ixc=input3(1)
      nz=input3(2)
      idummy=input3(3)

      if (Is_IO_Root_RSMPI()) then
        call headver('hsfp0_RSMPI',ixc)
        write(6,*) ' ixc nz=',ixc, nz
      endif

      if(ixc==0) then
        call RSMPI_Stop(' --- ixc=0 --- Choose computational mode!')
      endif

      hartree=2d0*rydberg()
      tetraex=tetra_hsfp0() ! always false
      iii=verbose()

      if (Is_IO_Root_RSMPI()) then
        write(6,*) ' verbose=',iii
      endif

c---tetraex is only for Ex.
      if(ixc==1.or.ixc==5.or.ixc==6) then;
      else; tetraex=.false. 
      endif


c Sergey Fallev's sxcf_fal2
c      sergeys =.false.
      if(ixc==11.or.ixc==12) then
c        sergeys=.true.
        ixc=ixc-10
      endif  
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
      if (Is_IO_Root_RSMPI()) then
        write(6,*) ' nqbz  =',nqbz
c      print *,  qbz 
        write(6,*) ' nqibz ngrp=',nqibz,ngrp
c      print *,' irk=',irk
c      print *,' #### idtetf: ####'
c      print *, idtetf
      endif
c set up work array
c      call wkinit (iwksize)
      call pshprt(60)


C--- readin GWIN and LMTO, then allocate and set datas.
      nwin =-999    !not readin NW file
      efin =-999d0  !not readin EFERMI
      if    (ixc==3) then;  incwfin= -2  !core exchange mode
      elseif(ixc==5) then;  incwfin= -4  ! valence-valence Ex energy mode. ! See rgwinf called from genallcf_*
      elseif(ixc==6) then;  incwfin= -3  ! core-valence Ex energy mode.
      else               ;  incwfin= -1  !use 7th colmn for core at the end section of GWIN
      endif
      call genallcf_v3(nwin,efin,incwfin) !in module m_genallcf_v3
      if(ngrp/= ngrp2) then
        call RSMPI_Stop('ngrp inconsistent: BZDATA and LMTO GWIN_V2')
      endif
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
      call getnemx8(nbmx,ebmx)
      if (Is_IO_Root_RSMPI()) then
        write(6,"('  nbmx ebmx from GWinput=',2i8,2d13.5)") nbmx,ebmx
      endif


c-------------------------------------------------------------------
c      if (nclass > mxclass) stop ' hsfp0: increase mxclass'
c!!!! WE ASSUME iclass(iatom)= iatom !!!!!!!!!!!!!!!!!!!!!!!!!
      if (nclass /= natom ) then
        call RSMPI_Stop("hsfp0: nclass /= natom")
      endif
      if (Is_IO_Root_RSMPI()) write(6,*) ' hsfp0: end of genallcf2'
c
      call pshprt(30)
      pi   = 4d0*datan(1d0)
      tpia = 2d0*pi/alat

      shtw = 0d0
      if(esmr<1d-5) shtw=0.01d0 ! Ferdi's shift to avoid resonance effect(maybe)
c
      call dinv33(plat,1,xxx,vol)
      voltot = dabs(vol)*(alat**3)

c--- ef is taken as rs for the empty-sphere test case of legas=T case -------------
      legas = .false.
      INQUIRE (FILE = 'LEGAS', EXIST = legas)
      if(legas) then !!! test for electron gas case. 
c RS: legas=T is disabled..
        call RSMPI_Stop("legas=T")
      endif     

ctm, 081020
      if (Is_IO_Root_RSMPI()) then
        call getkeyvalue("GWinput","sig_wan_gauge",lwang,default=.false.)
        write(*,*)'Wannier gauge mode =',lwang
      endif ! Is_IO_Root_RSMPI
      call MPI_Bcast(lwang,1,MPI_LOGICAL,io_root_rsmpi,
     &     MPI_COMM_WORLD,ierror_rsmpi)
      call RSMPI_Check("MPI_Bcast(lwang)",ierror_rsmpi)
c
c RS: only io_root opens output file
      ifexsp=0
      if(ixc==1) then
        exchange=.true.
        tote=.false.
        if (Is_IO_Root_RSMPI()) then
          write(6,*) ' --- Exchange mode --- '
          if (lwang) then
          ifxc(1)  = iopen('XCU.wan'//xt(nz),1,-1,0)
          ifsex(1) = iopen('SEXU.wan'//xt(nz),1,-1,0)
          else ! lwang
          ifxc(1)  = iopen('XCU'//xt(nz),1,-1,0)
          ifsex(1) = iopen('SEXU'//xt(nz),1,-1,0)
          endif ! lwang
        endif
        if (nspin == 2) then
          if (Is_IO_Root_RSMPI()) then
           if (lwang) then
            ifxc(1)  = iopen('XCU.wan'//xt(nz),1,-1,0)
            ifsex(1) = iopen('SEXU.wan'//xt(nz),1,-1,0)
           else ! lwang
            ifxc(2)  = iopen('XCD'//xt(nz),1,-1,0)
            ifsex(2) = iopen('SEXD'//xt(nz),1,-1,0)
           endif ! lwang 
          endif
        endif
        
        INQUIRE (FILE = 'EXspTEST', EXIST = exspwrite)
        if(exspwrite) then
c RS: exspwrite=T is disabled..
          call RSMPI_Stop("exspwrite=T")
        endif
      elseif(ixc==2) then
        exchange=.false.
        tote=.false.
        if (Is_IO_Root_RSMPI()) then
          write(6,*) ' --- Correlation mode --- '
          if(cohtest) write(6,*) ' COH calculation mode. Results in COH'
          if (lwang) then
          ifsec(1) = iopen('SECU.wan'//xt(nz),1,-1,0) ! output files
          if (nspin == 2) 
     .         ifsec(2) = iopen('SECD.wan'//xt(nz),1,-1,0)

          else ! lwang
          ifsec(1) = iopen('SECU'//xt(nz),1,-1,0) ! output files
          if (nspin == 2) 
     .         ifsec(2) = iopen('SECD'//xt(nz),1,-1,0)
          endif ! lwang
        endif
      elseif(ixc==3) then
        exchange=.true.
        tote=.false.
        esmr=0d0
        if (Is_IO_Root_RSMPI()) then
          write(6,*) ' --- CORE Exchange mode --- '
          if (lwang) then
          ifsex(1)   = iopen('SEXcoreU.wan'//xt(nz),1,-1,0)
          if (nspin == 2) then
            ifsex(2)   = iopen('SEXcoreD.wan'//xt(nz),1,-1,0)
          endif
          else ! lwang   
          ifsex(1)   = iopen('SEXcoreU.wan'//xt(nz),1,-1,0)
          if (nspin == 2) then
            ifsex(2)   = iopen('SEXcoreD.wan'//xt(nz),1,-1,0)
          endif
          endif ! lwang
        endif
c RS: ixc==4 for calculating sigma_c
      elseif(ixc==4) then
        exchange=.false.
        tote=.false.
        if (Is_IO_Root_RSMPI()) then
          write(6,*) ' --- Spectrum function Sigma(omega) mode --- '
        endif
      else
c RS: ixc==5,6 is disabled..
        call RSMPI_Stop(' hsfp0_RSMPI: Need input (std input) 1-4!')
      endif

c---
      if (Is_IO_Root_RSMPI()) then
        write(6, *) ' --- computational conditions --- '
        write(6,'("    deltaw  =",f13.6)') deltaw 
c      write(6,'("    ua      =",f13.6)') ua
        write(6,'("    esmr    =",f13.6)') esmr
        if(ixc==6) write(6,'("    esmr2    =",f13.6)') esmr2
        write(6,'("    alat voltot =",2f13.6)') alat, voltot
c     RS: this line is moved out of this block
c     if(nw/=-999) stop "nw here should be -999==nwin"
        write(6,'("    niw dw   =",i7,f13.6)') niw,dw
      endif
      if(nw/=-999) call RSMPI_Stop("nw here should be -999==nwin")


c>> read dimensions of wc,b,hb
      ifhbed     = iopen('hbe.d',1,0,0)
      read (ifhbed,*) nprecb,mrecb,mrece,nlmtot,nqbzt, nband,mrecg
      if (nprecb == 4) then
        call RSMPI_Stop('hsfp0_RSMPI: b,hb in single precision')
      endif

      call init_readeigen(ginv,nspin,nband,mrece) !initialization of readEigen

c --- get space group information ---------------------------------
c true class information in order to determine the space group -----------
c because the class in the generated GW file is dummy.(iclass(ibas)=ibas should be kept).
      open (102,file='CLASS')
      allocate(iclasst(natom),invgx(ngrp)
     &         ,miat(natom,ngrp),tiat(3,natom,ngrp),shtvg(3,ngrp))
      if (Is_IO_Root_RSMPI()) write(6,*) '  --- Readingin CLASS info ---'
      do ibas = 1,natom
        read(102,*) ibasx, iclasst(ibas)
        if (Is_IO_Root_RSMPI()) write(6, "(2i10)") ibasx, iclasst(ibas)
      enddo

c Get space-group transformation information. See header of mptaouof.
      call mptauof(symgg,ngrp,plat,natom,pos,iclasst
     o              ,miat,tiat,invgx,shtvg )
c        write (*,*)  'tiat=', tiat(1:3,1:natom,invr),invr
      
c Get array size to call rdpp
c      call getsrdpp( nclass,nl,  
c     o               ngpmx,ngcmx,nxx )
      call getsrdpp2( nclass,nl,nxx)
      call readngmx('QGpsi',ngpmx)
      call readngmx('QGcou',ngcmx)
      if (Is_IO_Root_RSMPI()) write(6,*) ' ngcmx ngpmx=',ngcmx,ngpmx

      allocate( nx(0:2*(nl-1),nclass), nblocha(nclass) ,lx(nclass),
     &   ppbrd ( 0:nl-1, nn, 0:nl-1,nn, 0:2*(nl-1),nxx, nspin*nclass),
     &   cgr(nl**2,nl**2,(2*nl-1)**2,ngrp))

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
      call rdpp_v3(nxx, nl, ngrp, nn, nclass, nspin, symgg,qbas,
     o      nblocha, lx, nx,  ppbrd , mdimx, nbloch, cgr)
c      nblochpmx = nbloch + ngcmx !oct2005

c      allocate(ngcni(nqibz)) !, ngveccB(3,ngcmx,nqibz)) !, ngveccBr(3,ngcmx,nqibz))
c   geigB(ngpmx,nband,nqbz,nspin),ngpn(nqbz),ngvecpB(3,ngpmx,nqbz),  
c     &   )  ! in IBZ 
c
c      call rdpp_pln(ngpmx,ngcmx, qibz,nqibz, qbz,nqbz,nband,nspin,
c     o      ngpn,geigB,ngvecpB,ngcni,ngveccB)
c
c      do iq = 1,nqibz
c        call readqg('QGcou',qibz(1:3,iq),ginv,  quu,ngcni(iq), ngveccB(1,1,iq))
c        write(6,"('--From QGcou  qibz quu ngc=',3f9.4,'  ',3f9.4,i5)")
c     &       qibz(1:3,iq),quu,ngcni(iq)
c      enddo
c for info
      allocate(ngvecp(3,ngpmx),ngvecc(3,ngcmx))
      call readqg('QGpsi',qibz(1:3,1),ginv, quu,ngpn1, ngvecp)
      call readqg('QGcou',qibz(1:3,1),ginv, quu,ngcn1, ngvecc)
      deallocate(ngvecp,ngvecc)

      if(Is_IO_Root_RSMPI()) then
        write(6,*) ' end of read QGcou' 
      endif

cccccccccccccccccccccccccccccccccccccccccccccccccccc
c      iqx=1
c      do ib=1,nband
c        write(6,'("  iband iqx sumgeigB=",2i3,12d12.3)') 
c     &   ib,iqx, sum(geigB(1:ngpn(1),ib,iqx))
c      enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      ib =1
c      iqx=1
c      do igp =1,nbloch + ngpn(iqx)
c         write(6,'("  igb ib iqx geigB=",3i3,1x,3i2,12d12.3)') 
c     &   igp,ib,iqx, ngvecpB(1:3,igp,iqx), geigB(igp,ib,iqx)
c      enddo
c      stop "xxxxxx zzz"      
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c----------------------------------------------
      call pshprt(60)


C--- Readin WV.d
c RS: for exchange..
      nwp = 0
      if(.not.exchange.or.(exchange.and.screen)) then   !screened exchange
        ifwd      = iopen('WV.d',1,-1,0) !direct access files WVR and WVI which include W-V.
        read (ifwd,*) nprecx,mrecl,nblochpmx,nwp,niwt, nqnum, nw_i 
                                                            !nw_i added at Feb2006
c nblochpmx from WV.d oct2005
        if (Is_IO_Root_RSMPI()) write(6,"(' Readin WV.d =', 10i5)") 
     &        nprecx, mrecl, nblochpmx, nwp, niwt, nqnum, nw_i 

c RS: checkeq() just checks whether nprecx == ndble  or not..
c       call checkeq(nprecx,ndble)
        if (nprecx .ne. ndble) then
          call RSMPI_Stop("dim of WVR and WVI not compatible")
        endif

c       call checkeq(nblochpmx,nblochpmx2) 
c       if (nwp /= nw)   stop 'hsfp0: wrong nw' !july 2003
        nw=nwp-1
        if (niwt /= niw) call RSMPI_Stop('hsfp0_RSMPI: wrong niw' )
        ifrcw     = iopen('WVR',0,-1,mrecl)
        if(.not.exchange) ifrcwi = iopen('WVI',0,-1,mrecl)
C... reading general energy mesh from file 'freq_r'
        open(UNIT=3111,file='freq_r') !this is in a.u.
        read(3111,*)nwxx        !number of energy points
        if(nwxx/= nw+1) call RSMPI_Stop(' freq_r nw /=nw')
        allocate(freq_r(nw_i:nw)) !freq_r(1)=0d0
        do iw = nw_i,nw
          read(3111,*)freq_r(iw)
        enddo
        close(3111)
      endif                     ! if (.not.exchange .or. ...)
c      write(6,*) freq_r
c      stop 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'

c... Readin eigen functions
c      ifev(1)   = iopen('EVU', 0,0,mrece)
c      if (nspin==2) ifev(2) = iopen('EVD', 0,0,mrece)


c-tetraex section ------------------------
c RS: tetraex = T is not implemented, so this section is removed.
      if(tetraex) then
        call RSMPI_Stop("tetraex=T is not implemented")
      else
        allocate(wtet(1,1,1,1),ibzx(1)) !dummy
      endif
c-end of tetra section --------------------------------------------

      if((ixc==3).or.(ixc==6)) then 
        ef = lowesteval() -1d-3 !lowesteigen(nspin,nband,qbz,nqbz) - 1d-3 !lowesteb was
        if(maxval(ecore(:,1:nspin))>ef) then
          call RSMPI_Stop('hsfp0 ixc=3 or 6:  ecore>evalence ')
        endif

c RS: ixc==6 (Exxc) is disabled
        if(ixc==6) then
          call RSMPI_Stop("ixc==6")
        endif

        if (Is_IO_Root_RSMPI())
     &   call getkeyvalue("GWinput","sig_shiftk",ldk,default=.false.)
        call MPI_Bcast(ldk,1,MPI_LOGICAL,io_root_rsmpi,
     &     MPI_COMM_WORLD,ierror_rsmpi)
        call RSMPI_Check("MPI_Bcast(ldk)",ierror_rsmpi)
c
      else                      ! if(esmr/=0d0) then
c --- determine Fermi energy ef for given valn (legas case) or corresponding charge given by z and konf.
! When esmr is negative, esmr is geven automatically by efsimplef.
        call efsimplef2a_RSMPI(nspin,wibz,qibz,ginv, 
     i        nband,nqibz
     i        ,konf,z,nl,natom,iclass,nclass
     i        ,valn, legas, esmr, !!! valn is input for legas=T, output otherwise. 
c     
     i        qbz,nqbz          !index_qbz, n_index_qbz,
     o        ,efnew)
c
c     print *,' end of efsimple'
ctm, 070130
!      ldk = .false.
!      INQUIRE (FILE = 'DK', EXIST = ldk)
      if (Is_IO_Root_RSMPI()) then
      call getkeyvalue("GWinput","sig_shiftk",ldk,default=.false.)
      if(ldk) then 
        ifdk = 2101
        open (ifdk,file='EFERMI')
        read(ifdk,*)efnew
        close(ifdk)
      endif     
      endif ! Is_IO_Root_RSMPI
      call MPI_Bcast(ldk,1,MPI_LOGICAL,io_root_rsmpi,
     &     MPI_COMM_WORLD,ierror_rsmpi)
      call RSMPI_Check("MPI_Bcast(ldk)",ierror_rsmpi)
      call MPI_Bcast(efnew,1,MPI_DOUBLE_PRECISION,io_root_rsmpi,
     &     MPI_COMM_WORLD,ierror_rsmpi)
      call RSMPI_Check("MPI_Bcast(efnew)",ierror_rsmpi)
      
         ef = efnew
         if (ixc==5) ef2 = ef
c RS: to avoid messy output...
         call MPI_Barrier(MPI_COMM_WORLD,ierror_rsmpi)
         call RSMPI_Check("MPI_Barrier",ierror_rsmpi)
c- check total ele number -------
         ntot  = nocctotg2(nspin, ef,esmr, qbz,wbz, nband,nqbz) !wbz
         call MPI_Barrier(MPI_COMM_WORLD,ierror_rsmpi)
         call RSMPI_Check("MPI_Barrier",ierror_rsmpi)

         if (Is_IO_Root_RSMPI()) then
           write(6,*) ' ef    =',ef
           write(6,*) ' esmr  =',esmr
           write(6,*) ' valn  =',valn
           write(6,*) ' ntot  =',ntot
         endif
c
c        if(ixc==1) then
c          ifefsm  = iopen('EFERMI.hsfp0.mode1',1,-1,0)
c          write(ifefsm,*) ef,'!ef by smearing. by hsfp0'
c          ifefsm  = iclose('EFERMI.hsfp0.mode1')
c        elseif(ixc==2) then
c          ifefsm  = iopen('EFERMI.hsfp0.mode2',1,-1,0)
c          write(ifefsm,*) ef,'!ef by smearing. by hsfp0'
c          ifefsm  = iclose('EFERMI.hsfp0.mode2')
c        endif
cc
cc        if(abs(valn-ntot)>1d-6) stop ' abs(valn-ntot)>1d-6'  !20001 May
cc
cxx      elseif(legas) then
cxx         ef = efx
      endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(tote) then
        call RSMPI_Stop("tote==T")
      endif 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c      noccxv = maxocc (ifev,nspin, ef+0.5d0*esmr, nband,nqbz)  ! maximum no. of occupied valence states
! maxocc seems to give (the maxmum number of occ + 1).


c-------------------------
  201 continue
c-------------------------
      call init_readeigen2(mrecb,nlmto,mrecg) !initialize m_readeigen

ctm, 081020
      if (lwang) then
        call init_readeigen_mlw(nwf,iko_ix,iko_fx)
        if (Is_IO_Root_RSMPI()) then
          write(*,*)'nwf =',nwf
          write(*,*)'iko_ix =',iko_ix
          write(*,*)'iko_fx =',iko_fx
        endif ! IO  
        allocate(umn(iko_ix:iko_fx,nwf))
      else
        iko_ix=1
        iko_fx=1
        nwf = 1
        allocate(umn(iko_ix:iko_fx,nwf)) ! dummy
      endif ! lwang  

c the coulomb matrix for exchange.
      if(exchange) ifvcfpout = iopen('VCCFP',0,-1,0) 

c QPNT data
      if(.not.tote) then
        if(nz==0) then
c          if(readgwinput()) then
          call getkeyvalue("GWinput","<QPNT>",unit=ifqpnt,status=ret) 
c          else
c           ifqpnt = iopen('QPNT',1,0,0)
c          endif
        else
          ifqpnt  = iopen('QPNT'//xt(nz),1,0,0)
        endif
      endif
      if (Is_IO_Root_RSMPI()) write(6,*) ' ifqpnt ret=',ifqpnt,ret
c     do i=1,ret
c         read(ifqpnt,'(a)') ddd
c         write(6,'(a)') ddd
c      enddo



c read q-points and states
      if(ixc==4) then
        call readxx(ifqpnt) !skip to ***** for q point for spectrum function.
        omegamaxin = 1d70
        read (ifqpnt,*,err=2038,end=2038) dwplot,omegamaxin
 2038   continue
c RS: for safety..
c        omegamax = 2*freq_r(nw) !This is in Ry.
        omegamax = 2*freq_r(nw-1) !This is in Ry.
c        print *,' nwp dw(a.u.) omegamax (Ry)=',nwp,dw,omegamax
        if(omegamaxin < omegamax) then
          if (Is_IO_Root_RSMPI()) then
            write(6,*)' --Use readin dwplot and omegamaxin from <QPNT>'
          endif
          omegamax = omegamaxin
        endif
        if( omegamax <0) call RSMPI_Stop('hsfp0 :strange omegamax <0 ')
        iwini =  -int( omegamax / dwplot )
        iwend =   int( omegamax/  dwplot )
        if (Is_IO_Root_RSMPI()) then
          write(6,*) ' iwini:iwend omegamax(Ry)=',iwini,iwend,omegamax
        endif
      else
        iwini = -1
        iwend = 1 
      endif
c ---
      if (tote) then 
        lqall      = .true.
      else
        lqall      = .false.
        laf        = .false.
        call readx   (ifqpnt,10)
        read (ifqpnt,*) iqall,iaf
        if (iqall == 1) lqall = .true.
        if (iaf   == 1)   laf = .true.
        call readx   (ifqpnt,100)
      endif
! states
      if (tote) then
        ntq    = noccxv
        allocate( itq(ntq) )
        do i   = 1,ntq
        itq(i) = i
        enddo
      else  
        read (ifqpnt,*) ntq
        allocate( itq(ntq) )
        read (ifqpnt,*) (itq(i),i=1,ntq)
      endif
! q-points
      if (lqall) then !all q-points case
        nq         = nqibz
        allocate(q(3,nq),qtmp(3,nq))
        call dcopy   (3*nqibz,qibz,1,q,1)
      else
        call readx   (ifqpnt,100)
        read (ifqpnt,*) nq
        allocate(q(3,nq),qtmp(3,nq))
        do       k = 1,nq
          read (ifqpnt,*) i,q(1,k),q(2,k),q(3,k)
          if (Is_IO_Root_RSMPI()) then
            write(6,'(i3,3f13.6)') i,q(1,k),q(2,k),q(3,k)
          endif
        enddo
      endif
c
      nspinmx = nspin
      if (laf) nspinmx =1
      close(ifqpnt)

ctm, 081026
      if (lwang) then
        deallocate(itq)
        call getkeyvalue("GWinput","wan_ntq",ntq,default=-999)
        if (ntq.le.0.or.ntq.gt.nwf) then
          if (Is_IO_Root_RSMPI()) write(*,*)'ntq,nwf',ntq,nwf
          call RSMPI_Stop('hsfp0: wrong ntq')
        endif  
        allocate(itq(ntq))
        call getkeyvalue("GWinput","wan_itq",itq,ntq)
        if (Is_IO_Root_RSMPI()) then
           write(*,*)'ntq =',ntq
           write(*,*)'itq =',(itq(i),i=1,ntq)
        endif ! IO   
      endif

ctm, 081002 sig_emin and sig_emax
      if (Is_IO_Root_RSMPI()) then
      call getkeyvalue("GWinput","sig_emin",eqtmin,default=-999d0)
      call getkeyvalue("GWinput","sig_emax",eqtmax,default=999d0)
      if (eqtmin.gt.-998d0.and.Is_IO_Root_RSMPI())
     &    write(*,*)'eqtmin=',eqtmin
      if (eqtmax.lt.998d0.and.Is_IO_Root_RSMPI())
     &    write(*,*)'eqtmax=',eqtmax
      eqtmin = eqtmin / rydberg()
      eqtmax = eqtmax / rydberg()

ctm, 070131
      dk = 0d0
      if(ldk) then 
        ifdk = 2101
        open (ifdk,file='DK')
        read(ifdk,*)dk(1),dk(2),dk(3)
        close(ifdk)
      endif     
      endif ! Is_IO_Root_RSMPI
      call MPI_Bcast(eqtmin,1,MPI_DOUBLE_PRECISION,io_root_rsmpi,
     &     MPI_COMM_WORLD,ierror_rsmpi)
      call RSMPI_Check("MPI_Bcast(eqtmin)",ierror_rsmpi)
      call MPI_Bcast(eqtmax,1,MPI_DOUBLE_PRECISION,io_root_rsmpi,
     &     MPI_COMM_WORLD,ierror_rsmpi)
      call RSMPI_Check("MPI_Bcast(eqtmax)",ierror_rsmpi)
      call MPI_Bcast(dk,3,MPI_DOUBLE_PRECISION,io_root_rsmpi,
     &     MPI_COMM_WORLD,ierror_rsmpi)
      call RSMPI_Check("MPI_Bcast(dk)",ierror_rsmpi)

      do iq = 1,nq
        q(:,iq) = q(:,iq) + dk(:)
      enddo
C...
      if (Is_IO_Root_RSMPI())
     & call getkeyvalue("GWinput","QPNT_nbandrange",nss,2 !nov2005
     &               ,default=(/-99997,-99997/) )
      call MPI_Bcast(nss,2,MPI_INTEGER,io_root_rsmpi,
     &     MPI_COMM_WORLD,ierror_rsmpi)
      call RSMPI_Check("MPI_Bcast(nse)",ierror_rsmpi)

      if(nss(2)/=-99997) then
        if(allocated(itq)) deallocate(itq)
        ntq=nss(2)-nss(1)+1
        allocate( itq(ntq) )
        do i=max(1,nss(1)),min(nss(2),nband)
          itq(i-nss(1)+1) = i
        enddo
      endif

C...
      allocate(omega(ntq,iwini:iwend))
      if(ixc==4) then
        do iw = iwini,iwend      
          omega(1:ntq,iw) =  dwplot* iw + ef
        enddo
      endif
      
C ---  q near zero
      if (Is_IO_Root_RSMPI()) then
        write(6,*) 'reading QOP'
      endif

      open (101,file='Q0P')
      read (101,"(i5)") nq0i
      if(.not.exchange) then        !call checkeq(nqibz+nq0i-1, nqnum)
        if (nqibz+nq0i-1 .ne. nqnum) then
          call RSMPI_Stop("dim of WVR and WVI not compatible")
        endif
      endif

      if(Is_IO_Root_RSMPI()) then
        write(6,*) ' *** nqibz nq0i_total=', nqibz,nq0i
      endif
      nq0it = nq0i
      allocate( wqt(1:nq0i),q0i(1:3,1:nq0i) )
c      read (101,"(d24.16,3x, 3d24.16)" )( wqt(i),q0i(1:3,i),i=1,nq0i)
      nq0ix = nq0i
      do i=1,nq0i
        read (101,* ) wqt(i),q0i(1:3,i)
        if(wqt(i)==0d0 ) nq0ix = i-1
      enddo
      nq0i = nq0ix ! New nq0i July 2001
      if (Is_IO_Root_RSMPI()) then
        write(6,*) ' Used k number in Q0P =', nq0i
        write(6,"(i3,f14.6,2x, 3f14.6)" )(i, wqt(i),q0i(1:3,i),i=1,nq0i)
      endif
      close(101)
      allocate( wgt0(nq0i,ngrp) )
C Sergey's 1stFeb2005
c      call q0iwgt2(symgg,ngrp,wqt,q0i,nq0i, 
c     o            wgt0)
      if (Is_IO_Root_RSMPI())
     & call getkeyvalue("GWinput","allq0i",allq0i,default=.false.)!S.F.Jan06
      call MPI_Bcast(allq0i,1,MPI_LOGICAL,io_root_rsmpi,
     &     MPI_COMM_WORLD,ierror_rsmpi)
      call RSMPI_Check("MPI_Bcast(allq0i)",ierror_rsmpi)
      call q0iwgt3(allq0i,symgg,ngrp,wqt,q0i,nq0i,               !S.F.Jan06
     o            wgt0)                   ! added allq0i argument
C--------------------------

      if(nq0i/=0) then
        if (Is_IO_Root_RSMPI()) then
          write(6,*) ' *** tot num of q near 0   =', 1/wgt0(1,1)      
        endif
      endif
      if (Is_IO_Root_RSMPI()) then
        write(6,"('  sum(wgt0) from Q0P=',d14.6)")sum(wgt0)
      endif
      if(bzcase()==2) then
        wgt0= wgt0*wgtq0p()/dble(nqbz)
        if (Is_IO_Root_RSMPI()) then
           write(6,"('bzcase=2:  sum(wgt0_modified )=',d14.6)")sum(wgt0)
        endif
      endif

c --- qbze(3,nqibze)
      nqbze  = nqbz *(1 + nq0it)
      allocate( qbze(3, nqbze) )
      call dcopy(3*nqbz, qbz, 1, qbze,1)
      do i = 1,nq0it
        ini = nqbz*(1 + i -1)
        do ix=1,nqbz
          qbze (:,ini+ix)   = q0i(:,i) + qbze(:,ix)
        enddo
      enddo

c --- read LDA eigenvalues
c     ntp0=ntq
      allocate(eqx(ntq,nq,nspin),eqx0(ntq,nq,nspin),eqt(nband))
      eqx0 = 0d0
      do      is = 1,nspin
        do      ip = 1,nq
c        iq       = idxk (q(1,ip),qbze,nqbze)
c        call rwdd1   (ifev(is), iq, nband, eqt) !direct access read b,hb and e(q,t)
          call readeval(q(1,ip),is,eqt)
c        print *,' eqt=',eqt
         if (lwang) then
          call readumn_mlw(q(1,ip),is, umn)
          do i = 1,ntq
          do itp = iko_ix,iko_fx
            eqx0(i,ip,is) = eqx0(i,ip,is)
     &                    + eqt(itp)*abs(umn(itp,itq(i)))**2
            eqx (1:ntq,ip,is) = rydberg()*(eqx0(1:ntq,ip,is)- ef)
          enddo ! itp  
          enddo ! i
         else ! lwang
          eqx0(1:ntq,ip,is) = eqt(itq(1:ntq))
          eqx (1:ntq,ip,is) = rydberg()*(eqt(itq(1:ntq))- ef)
         endif
        enddo
      enddo
      deallocate(eqt)

c --- info
      if(tote) then
        call winfo2(6,nspin,nq,ntq,is,nbloch
     &    ,ngpn1,ngcn1,nqbz,nqibz,ef,ef2,deltaw,alat,esmr,esmr2)
      else
        if (Is_IO_Root_RSMPI()) then
          call winfo(6,nspin,nq,ntq,is,nbloch
     &           ,ngpn1,ngcn1,nqbz,nqibz,ef,deltaw,alat,esmr)
        endif
      endif

ctm, 070621
      if (Is_IO_Root_RSMPI()) then
      call getkeyvalue("GWinput","sig_fullmat",loffd,default=.false.)
      if (loffd) write(6,*)"Calculates full matrix elements of Sigma"
      endif ! Is_IO_Root_RSMPI
      call MPI_Bcast(loffd,1,MPI_LOGICAL,io_root_rsmpi,
     &     MPI_COMM_WORLD,ierror_rsmpi)
      call RSMPI_Check("MPI_Bcast(loffd)",ierror_rsmpi)

c-------------------------
c LDA exchange-correlation
c-------------------------
c- takao Vxc is readin from VXCFP
      if(ixc==1) then
ctm, 070621
       if ((lwang.or.loffd).and.Is_IO_Root_RSMPI()) then
c read Vxc in lmto basis
        open(UNIT=3007,file='v_xc',form='UNFORMATTED',status='OLD')
        ifv_xc=3007
        open(UNIT=3008,file='evec',form='UNFORMATTED',status='OLD')
        ifevec=3008
        read(ifv_xc) ndimh, nspin2,nqqq
        if (nspin2 .ne. nspinmx)
     &   call RSMPI_Stop('hsfp0: wrong nspin in v_xc')
c        if (nqqq .ne. nq) stop 'hsfp0: wrong nq in v_xc'
        if (ndimh .ne. nband)
     &   call RSMPI_Stop('hsfp0: wrong ndimh in v_xc')
        read(ifevec) ndimh2, nspin2,nqqqx
        if (nspin2 .ne. nspinmx)
     &   call RSMPI_Stop('hsfp0: wrong nspin in evec')
        if (nqqqx .ne. nqqq)
     &   call RSMPI_Stop('hsfp0: wrong nq in evec')
        if (ndimh2 .ne. ndimh)
     &   call RSMPI_Stop('hsfp0: wrong ndimh in evec')
        allocate(v_xc(ndimh,ndimh,nqqq,nspin),
     &           evec(ndimh,ndimh,nqqq,nspin),qqq(3,nqqq))
        do iq = 1,nqqq
        do is = 1,nspinmx
          read(ifv_xc) qqq(1:3,iq),v_xc(1:ndimh,1:ndimh,iq,is)
          read(ifevec) qqqx(1:3),evec(1:ndimh,1:ndimh,iq,is)
          if (sum(abs( qqq(:,iq)-qqqx(:) )) .gt. 1d-10)
     &       call RSMPI_Stop('hsfp0: qqq/=qqqx')
        end do ! is
        end do ! iq
        close(ifv_xc)
        close(ifevec)

c calc Vxc in KS basis
        allocate(vxc(ntq,ntq,nq,nspinmx))
        if (lwang)allocate(vxc2(iko_ix:iko_fx,iko_ix:iko_fx))
        do ip = 1,nq
        ndummy = 0
        do iq = 1,nqqq
         if (sum(abs( qqq(:,iq)-q(:,ip) )) .gt. 1d-10) cycle
         ndummy = ndummy + 1
         do is  = 1,nspinmx
         if (lwang) then
           vxc2 = 0d0
           do itp = iko_ix,iko_fx
           do itpp= iko_ix,iko_fx
             vxc2(itp,itpp) = sum(dconjg(evec(:,itp,iq,is))*
     &          matmul(v_xc(:,:,iq,is),evec(:,itpp,iq,is)))
           enddo ! itpp
           enddo ! itp
           call readumn_mlw(q(1,ip),is, umn)
           do itp  = 1,ntq
           do itpp = 1,ntq
             it1 = itq(itp)
             it2 = itq(itpp)
             vxc(itp,itpp,ip,is) = sum(dconjg(umn(:,it1))*
     &          matmul(vxc2(:,:),umn(:,it2)))
          enddo ! itpp
          enddo ! itp
          vxc(:,:,ip,is) = rydberg()*vxc(:,:,ip,is)
         else ! lwang   
         do itp = 1,ntq
         do itpp= 1,ntq
          it1 = itq(itp)
          it2 = itq(itpp)
          vxc(itp,itpp,ip,is) = rydberg()*sum(dconjg(evec(:,it1,iq,is))*
     &          matmul(v_xc(:,:,iq,is),evec(:,it2,iq,is)))
         enddo
         enddo
         endif ! lwang
         enddo ! is
        enddo ! iq
        if (ndummy .ne. 1) call RSMPI_Stop('hsfp0: qqq .ne. q')
        enddo ! ip
       if (loffd) then ! output
c output:
        do is = 1,nspinmx
         if (lwang) then
         if (is.eq.1) then
            ifvxc2  = iopen('XC2U.wan',1,-1,0)
         else
            ifvxc2  = iopen('XC2D.wan',1,-1,0)
         endif   
         else ! lwang
         if (is.eq.1) then
            ifvxc2  = iopen('XC2U',1,-1,0)
         else
            ifvxc2  = iopen('XC2D',1,-1,0)
         endif   
         endif ! lwang
         write(ifvxc2,"(4i5,f12.6)")ntq,nq,nqbz,nspin,ef
         write(ifvxc2,*)q
         write(ifvxc2,*)itq
         write(ifvxc2,*)eqx
         do ip = 1,nq
c          write(ifvxc2,*)vxc(:,:,ip,is)/rydberg()/2d0
          write(ifvxc2,*)'ip=',ip
          do ij = 1,ntq
          do ii = 1,ntq
            write(ifvxc2,*)ii,ij,vxc(ii,ij,ip,is)/rydberg()/2d0
          enddo  
          enddo  
         end do ! ip
         if (lwang) then
         if (is.eq.1) then
            isx = iclose('XC2U.wan')
         else   
            isx = iclose('XC2U.wan')
         endif   
         else ! lwang
         if (is.eq.1) then
            isx = iclose('XC2U')
         else   
            isx = iclose('XC2U')
         endif   
         endif ! lwang
        end do ! is
       endif ! loffd  ! output
!        deallocate(v_xc,evec,qqq,vxc)
        deallocate(v_xc,evec,qqq)
        if (lwang)deallocate(vxc2)
       endif ! (lwang.or.loffd).and.Is_IO_ROOT_RSMPI()

        allocate(  vxcfp(ntq,nq,nspin) )
       if (lwang) then
          if (Is_IO_Root_RSMPI()) then
          do i = 1,ntq
             vxcfp(i,:,:) = real(vxc(i,i,:,:))
          enddo ! i   
          endif ! IO_Root
       else ! lwang
c        call rsexx(nspin,itq,q,ntq,nq, vxcfp)
        call rsexx2(ginv,nspin,itq,q,ntq,nq, vxcfp)
       endif ! lwang
c loop over spins
c RS: only io_root prints infomation
        if (Is_IO_Root_RSMPI()) then
          do is = 1,nspinmx
            write (ifxc(is),*) '==================================='
            write (ifxc(is),"(' LDA exchange-correlation : is=',i3)")is
            write (ifxc(is),*) '==================================='
            call winfo(ifxc(is),nspin,nq,ntq,is,nbloch
     &           ,ngpn1,ngcn1,nqbz,nqibz,ef,deltaw,alat,esmr)
            write (ifxc(is),*)' ***'
            write (ifxc(is),"(a)") ' jband   iq ispin                  
     &                 qvec 
     &                                       eigen-Ef (in eV)
     &          LDA XC (in eV)'  
            ifoutsex = ifxc(is)
            write(6,*)
            do ip = 1,nq
              do i  = 1,ntq
                write(ifoutsex,"(3i5,3d24.16,3x,d24.16,3x,d24.16)")
     &               itq(i),ip,is, q(1:3,ip), eqx(i,ip,is), 
     &               vxcfp(i,ip,is) 
                write(6,"(' j iq isp=' i3,i4,i2,'  q=',3f8.4,
     &               '  eig=',f10.4,'  Sxc(LDA)=',f10.4)") 
     &               itq(i),ip,is, q(1:3,ip), eqx(i,ip,is),
     &               vxcfp(i,ip,is) 
              end do
            end do
c end of spin-loop
            if(is==1) isx = iclose('XCU'//xt(nz))
            if(is==2) isx = iclose('XCD'//xt(nz))
          end do
        endif                   ! if (Is_IO_Root_RSMPI())
        deallocate(vxcfp)

      endif ! ixc==1

c pointer to optimal product basis
      allocate(imdim(natom))
      call indxmdm (nblocha,nclass,
     i              iclass,natom,
     o              imdim )

      if(niw/=0) then           ! generate gaussian frequencies x between (0,1) and w=(1-x)/x
        allocate(freqx(niw),freqw(niw),wwx(niw)) !,expa(niw))
        call freq01x  (niw,     !ua,
     o       freqx,freqw,wwx)   !, expa)
      endif

c ------ write energy mesh ----------
c      if(.not.sergeys) then
      if(.not.exchange.or.(exchange.and.screen)) then   !screened exchange
      if (Is_IO_Root_RSMPI()) ifemesh = iopen('emesh.hsfp0'//xt(nz),1,-1,0)
      deltax0 = 0d0
      if (Is_IO_Root_RSMPI()) then
        call writeemesh(ifemesh,freqw,niw,freq_r,nwp,deltax0)
c subroutine writeemesh(ifi,freqi,niw,freqr,nnw,delta)
      endif
      endif
c
      iii=ivsumxxx(irk,nqibz*ngrp)
      if (Is_IO_Root_RSMPI()) write(6,*) " sum of nonzero iirk=",iii, nqbz

C... Read pomatr
c RS: smbasis==T is not implemented
      if(smbasis()) then
         call RSMPI_Stop("smbasis==T is not implemented!")
      else                      !dummy
         nkpo = 1
         nnmx =1
         nomx =1 
         allocate( pomatr(nnmx,nomx,nkpo), qrr(3,nkpo),nor(nkpo),nnr(nkpo) )
      endif

c-----------------------------------------------------------
c calculate the correlated part of the self-energy SEc(qt,w)
c-----------------------------------------------------------
c arrays for sxcf.f
      nlnx4    = nlnx**4
      niwx     = max0 (nw+1,niw)
      allocate( ppb(nlnmx*nlnmx*mdimx*nclass),  eq(nband), 
     &     kount(nqibz,nq), zsec(iwini:iwend,ntq,nq),  
     &     coh(ntq,nq) ) 
ctm, 070621
c      if (loffd) allocate(zsec2(ntq,ntq,nq))
      allocate(zsec2(ntq,ntq,nq,iwini:iwend))      

c RS: set MPI parameters
c RS: see gwsrc/RSMPI_rotkindex_mod.F
      call MPI_Barrier(MPI_COMM_WORLD,ierror_rsmpi)
      call RSMPI_Check("MPI_Barrier",ierror_rsmpi)
      call setup_rotkindex(ngrp,irk,wgt0,bzcase(),nqibz,nq0i,nq)

      if (Is_IO_Root_RSMPI()) then ! debug
        if(nq0i/=0) then
          write(6,*) 'RS: total number of k-points should be',
     &         nqbz +  1/wgt0(1,1)   - 2 + bzcase()   
        else
          write(6,*) 'RS: total number of k-points should be',
     &         nqbz  - 2 + bzcase()   
        endif
      endif
      call MPI_Barrier(MPI_COMM_WORLD,ierror_rsmpi)
      call RSMPI_Check("MPI_Barrier",ierror_rsmpi)
c RS: openlogfile for each process
      if (ixc .eq. 3) then
        ifile_rsmpi = iopen ('lsxC.RSMPI'//myrank_id_rsmpi,1,3,0) 
      else if (ixc .eq. 1) then
        ifile_rsmpi = iopen ('lsx_sf.RSMPI'//myrank_id_rsmpi,1,3,0) 
      else if (ixc .eq. 2) then         
        ifile_rsmpi = iopen ('lsc_sf.RSMPI'//myrank_id_rsmpi,1,3,0) 
      else if (ixc .eq. 4) then         
        ifile_rsmpi = iopen ('lsigm.RSMPI'//myrank_id_rsmpi,1,3,0) 
      else
        call RSMPI_Stop("unknown ixc")
      endif
c RS: print how symmetry operations and k-points are devided ..
      write(ifile_rsmpi,*) "rank : ", myrank_id_rsmpi
      write(ifile_rsmpi,*) "nrotk_local:",nk_local_qkgroup 
      write(ifile_rsmpi,*) "nrot_local :",nrot_local_rotk
      if (nrot_local_rotk .gt. 0) then
        write(ifile_rsmpi,*)
     &        "irot_index :",irot_index_rotk(1:nrot_local_rotk)
      endif
      write(ifile_rsmpi,*) "nk_local(1:ngrp) :"
      write(ifile_rsmpi,*) nk_local_rotk(:)
      do irot=1,ngrp
        if (nk_local_rotk(irot).gt.0) then
          write(ifile_rsmpi,*) "> irot,nk_local(irot) = ",
     &         irot, nk_local_rotk(irot)
          write(ifile_rsmpi,*) "   ik_index : ",
     &         ik_index_rotk(irot,1:nk_local_rotk(irot))
        endif
      enddo
ccccccccccccccc

      call MPI_Barrier(MPI_COMM_WORLD,ierror_rsmpi)
      call RSMPI_Check("MPI_Barrier",ierror_rsmpi)

      if (Is_IO_Root_RSMPI()) then
        write(6,*) "RS: loop over spin --"
      endif
c loop over spin ----------------------------------------------------
      if (tote) exx = 0d0

      do 2000 is = 1,nspinmx
        if (Is_IO_Root_RSMPI()) then
          if((exchange).and.(.not.tote)) then
            write(ifsex(is),*) '======================================='
            write(ifsex(is),"('Self-energy exchange SEx(q,t): is=',i3)") is
            write(ifsex(is),*) '======================================='

            call winfo(ifsex(is),nspin,nq,ntq,is,nbloch
     &           ,ngpn1,ngcn1,nqbz,nqibz,ef,deltaw,alat,esmr)

            write (ifsex(is),*)' *** '
            write (ifsex(is),"(a)") ' jband   iq ispin                  
     &           qvec 
     &                                       eigen-Ef (in eV)
     &           exchange (in eV)'  

          elseif(ixc==2) then
            write(ifsec(is),*) '=========================================='
            write(ifsec(is),"('Self-energy correlated SEc(qt,w): is=',i3)")
     &           is
            write(ifsec(is),*) '=========================================='

            call winfo(ifsec(is),nspin,nq,ntq,is,nbloch
     &           ,ngpn1,ngcn1,nqbz,nqibz,ef,deltaw,alat,esmr)

            write (ifsec(is),*)' *** '
            write (ifsec(is),"(a)") ' jband   iq ispin                  
     &           qvec 
     &           eigen-Ef (in eV)
     &           Re(Sc) 3-points (in eV)                        
     &           In(Sc) 3-points (in eV)     
     &           Zfactor'   
          endif
        endif                   !if (ioroot)
c initialise secq and kount
        zsec  = 0d0
        coh   = 0d0
        kount = 0
ctm, 070621
        zsec2 = 0d0
c 
        if(ixc==3.and.nctot==0) goto 2001 !make dummy SEXcore

        if(Is_IO_Root_RSMPI()) write(6,*) "loop over rotations"
        call RSMPI_Print_WTime()
c loop over rotations -------------------------------
c         do 1000 irot = 1,ngrp
c RS: local index
        do 1000 irot_local = 1,nrot_local_rotk
c RS: correct index
          irot = irot_index_rotk(irot_local)

          if( sum(abs( irk(:,irot) )) ==0 .and.
     &         sum(abs( wgt0(:,irot))) == 0d0 ) then
c RS: this should not happen... see gwsrc/RSMPI_rotindex_mod.F
            call RSMPI_Stop("hsfp0_RSMPI,
     & cylce occurs in do 1000 -loop!")
            cycle
          endif

          write (ifile_rsmpi,"(i3,'  out of ',i3,'  rotations ',$)")
     &           irot,ngrp
          call cputid (ifile_rsmpi)

c        if (irot == 1 .or. irot == ngrp) then
c        call cputid(0); write(*,*)' ppba '
c        endif

c rotate atomic positions invrot*R = R' + T
          invr       = invrot (irot,invg,ngrp)

c -- ppb= <Phi(SLn,r) Phi(SL'n',r) B(S,i,Rr)> 
          call ppbafp_v2 (irot,ngrp,is,nspin,
     i           il,in,im, nlnm,
     d           nl,nn,nclass,nlnmx,
     i           mdimx,lx,nx,nxx, !Bloch wave    
     i           cgr, nl-1,     !rotated CG
     i           ppbrd,         !radial integrals
     o           ppb)

cc -- Rotated gvecc
c        call rotgvec(symgg(:,:,irot), nqibz, 
c     i    ngcmx,ngcni,qbas,ngveccB, 
c     o    ngveccBr)

c------------------------------------------------------
c calculate the correlated part of the self-energy within GW
c        ntqx = 0
c        if(tetra.and.(.not.exchange)) then
c        ntqx =3*ntq
c        endif

c
c loop over q
c            do 1001 ip = 1,nq
c RS: local index
          do 1001 ip_local = 1,nq_local_qkgroup
c RS: correct index
            ip = iq_index_qkgroup(ip_local)

            if(ixc==2) then
c          iq  = idxk (q(1,ip),qbz,nqbz)
c          call rwdd1 (ifev(is),iq, nband,eq)
              call readeval(q(1,ip),is,eq)
ctm 020621
             itini = 1
             itfin = ntq
             if (.not.lwang) then
              do i  = 1,ntq
                if (eq(itq(i))-ef < eqtmin) itini = i+1
              enddo
              do i  = ntq,1,-1
                if (eq(itq(i))-ef > eqtmax) itfin = i-1
              enddo
             endif ! lwang
              if (Is_IO_Root_RSMPI()) then
                write(*,*)'ip,itini,itfin'
                write(*,*)ip,itini,itfin
              endif
c
              if (lwang) then
               call readumn_mlw(q(1,ip),is, umn)               
               do iw = iwini,iwend
               do i  = 1,ntq
                rtmp = 0d0
                do itp = iko_ix,iko_fx
                  rtmp = rtmp + eq(itp)*abs(umn(itp,itq(i)))**2
                enddo   
                omega(i,iw) = rtmp + 2d0*(dble(iw)-shtw)*deltaw
               enddo
               enddo
              else ! lwang   
              do iw = iwini,iwend
              do i  = 1,ntq
                omega(i,iw) = eq(itq(i)) + 2d0*(dble(iw)-shtw)*deltaw
              enddo
              enddo
              endif ! lwang
            elseif(ixc==4) then
              itini = 1
              itfin = ntq
            endif ! ixc

            if(nctot==0) then
              deallocate(ecore)
              allocate(ecore(1,2)) !this is dummry
            endif   

c       if(sergeys) then
c RS:
            call sxcf_fal3z_RSMPI(kount, ip,q(1,ip),irot,itq,ntq,
     i              itini,itfin,ef,ef2,esmr,esmr2,
     i              tiat(1:3,1:natom,invr),miat( 1:natom,invr), 
cv2
     i              nspin,is, !ifcphi,ifrb(is),ifcb(is),ifrhb(is),ifchb(is),
     i              ifrcw,ifrcwi,
     i              qbas,ginv,qibz,qbz,wbz,nstbz, wibz, !iindxk,
     i              nstar,irk,     !kount,
c
c     i        iiclass,nblocha,i_mnlv,i_mnlc,iicore,incore,iimdim,
     i        iclass,nblocha,nlnmv, nlnmc,   !w(i_mnlv),w(i_mnlc)
     i        icore,ncore, imdim,
     i              ppb, ! pdb,dpb,ddb,
     i              freq_r,freqx, wwx,  !expa,ua,
     i              dw,  !deltaw,
     i              ecore(:,is), 
c
     d        nlmto,nqibz,nqbz,nctot,
c     i        index_qbz, n_index_qbz,
     d              nl,nnc,nclass,natom,
     d              nlnmx,mdimx,nbloch,ngrp, nw_i,nw,  niw,niwx,nq,
!                                            nw_i,nw instead of nwp feb2006
c
c     i     nblochpmx, ngpn,ngcni,ngpmx,ngcmx,
     i     nblochpmx,ngpmx,ngcmx,
c     i     geigB(1,1,1,is), ngvecpB,ngveccBr,
c     i     ngveccBr,
     i     wgt0,nq0i,q0i, symgg(:,:,irot),alat,
     i     matmul(symgg(:,:,irot),shtvg(:,invr)),nband,
     i     ifvcfpout,
c     i     shtw, 
     i     exchange, tote, screen, cohtest, ifexsp(is),
     i        omega, iwini,iwend,
     i     nbmx(2),ebmx(2), !takao 18June2003 
     i   pomatr, qrr,nnr,nor,nnmx,nomx,nkpo,          !oct2005 for pomat 
ctm 070621
     i        loffd,
ctm 081021
     i        lwang,iko_ix,iko_fx,nwf,umn,
     o        zsec,zsec2,coh,exx) ! acuumulation variable
c        else                                                 !june2002 add esmr2 takao
c         call sxcf_v2z (kount, ip,q(1,ip),irot,itq,ntq,ef,ef2,esmr,esmr2,
c     i              tiat(1:3,1:natom,invr),miat( 1:natom,invr), 
c     i              nspin,is, !ifcphi,ifrb(is),ifcb(is),ifrhb(is),ifchb(is),
c     i              ifrcw,ifrcwi,
c     i              qbas,ginv,qibz,qbz,wbz,nstbz, wibz,
c     i              nstar,  irk,     !iindxk,kount,
c     i        iclass,nblocha,nlnmv, nlnmc, !w(i_mnlv),w(i_mnlc)
c     i        icore,ncore, imdim,
c     i              ppb, ! pdb,dpb,ddb,
c     i              freq, freqx, wwx, expa,
c     i              ua,dw,  !deltaw,
c     i              ecore(:,is), 
c     d        nlmto,nqibz,nqbz,nctot,
c     d              nl,nnc,nclass,natom,
c     d              nlnmx,mdimx,nbloch,ngrp,nwp,niw,niwx,nq,
c     i     nblochpmx,ngpmx,ngcmx,
c     i     wgt0,nq0i,q0i, symgg(:,:,irot),alat,
c     i     matmul(symgg(:,:,irot),shtvg(:,invr)),nband,
c     i     ifvcfpout,
c     i     exchange, tote, screen, cohtest, ifexsp(is),
c     &     wtet(:,is,:,0),
c     &     ibzx,tetraex,
c     &     ixc,    
cc
c     i        omega, iwini,iwend,
c     i     nbmx,ebmx, !takao 18June2003 
c     o        zsec,coh,exx) ! acuumulation variable
c        endif
cccccccccccccccccccccccccccccccccccccccccc
c        iii = ivsum(kount,nqibz*nq)
c        print *," sumkount 2=",nqibz,nq,iii
c        stop "--- kcount test end --- "
cccccccccccccccccccccccccccccccccccccccccc
 1001     continue              !end of q-loop
 1000   continue                !end of rotation-loop
cccccccccccccccccccccccc


C RS: accumulate zsec
        allocate(zsec_sum(iwini:iwend,ntq,nq))
        call MPI_AllReduce(zsec,zsec_sum,(iwend-iwini+1)*ntq*nq,
     &         MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierror_rsmpi)
        call RSMPI_Check("MPI_AllReduce,zsec",ierror_rsmpi)
        zsec(:,:,:) = zsec_sum(:,:,:)
        deallocate(zsec_sum)
c TM: accumulate zsec2
        allocate(zsec2_sum(ntq,ntq,nq,iwini:iwend))
        call MPI_AllReduce(zsec2,zsec2_sum,(iwend-iwini+1)*ntq*ntq*nq,
     &         MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierror_rsmpi)
        call RSMPI_Check("MPI_AllReduce,zsec",ierror_rsmpi)
        zsec2(:,:,:,:) = zsec2_sum(:,:,:,:)
        deallocate(zsec2_sum)

c check that all k { FBZ have been included
c      if (ivsum(kount,nqibz*nq) /= nqbz*nq) then
c        iii = ivsum(kount,nqibz*nq)
c        print *," ivsum=",iii, nqbz*nq 
c        stop 'hsfp0: missing k-pts'
c      endif

c---------------------------------
c electron gas bare exchange (exact)
        if (legas.and.exchange.and.(.not.tote)) then
          call RSMPI_Stop("legas==T")
        endif

 2001   continue
c---------------------------------
      
c write SEx(q,t)
c      call cv      (hartree,sexq,ntq*nq,sexq) !Because vcoul is in a.u.

        if(ixc==5.or.ixc==6) then
! nothing
        elseif(exchange) then
c        do       ip = 1,nq
c          call wsex2 (ifsex(is),q(1:3,ip),itq(1:ntq),eqx(1:ntq,ip,is),
c     &              hartree*dreal(zsec(-1,1:ntq,ip)), !sexq(1:ntq,ip),
c     d              ntq)
c        end do
c
c-------------------------
          if (Is_IO_Root_RSMPI()) then
            ifoutsex=ifsex(is)
            write(6,*)
            do ip = 1,nq
              do i  = 1,ntq
                write(ifoutsex,"(3i5,3d24.16,3x,d24.16,3x,d24.16)")
     &               itq(i),ip,is, q(1:3,ip), eqx(i,ip,is), 
     &               hartree*dreal(zsec(iwini,i,ip)) 
                write(6,"(' j iq isp=' i3,i4,i2,'  q=',3f8.4,
     &               '  eig=',f10.4,'  Sx=',f10.4)") 
c     &     '  eig=',f10.4,'  Sx=',d15.7)") 
     &                    itq(i),ip,is, q(1:3,ip), eqx(i,ip,is),
     &                    hartree*dreal(zsec(iwini,i,ip)) 
              end do
            end do
ctm, 070621
            if (loffd) then
              if (ixc.eq.1) then
               if (lwang) then
                if (is.eq.1) then
                  ifse2= iopen('SEX2U.wan',1,-1,0)
                else
                  ifse2= iopen('SEX2D.wan',1,-1,0)
                endif
               else ! lwang   
               if (is.eq.1) then
                 ifse2= iopen('SEX2U',1,-1,0)
               else
                 ifse2= iopen('SEX2D',1,-1,0)
               endif
               endif ! lwang
              else ! ixc
               if (lwang) then
                if (is.eq.1) then
                  ifse2= iopen('SEXcore2U.wan',1,-1,0)
                else
                  ifse2= iopen('SEX2coreD.wan',1,-1,0)
                endif
               else ! lwang 
               if (is.eq.1) then
                 ifse2= iopen('SEXcore2U',1,-1,0)
               else
                 ifse2= iopen('SEX2coreD',1,-1,0)
               endif
               endif ! lwang
              endif ! ixc
              write(ifse2,"(4i5,f12.6)")ntq,nq,nqbz,nspin,ef
              write(ifse2,*)q
              write(ifse2,*)itq
              write(ifse2,*)eqx
              do ip = 1,nq
c               write(ifse2,*)zsec2(1:ntq,1:ntq,ip)
               write(ifse2,*)'ip=',ip
               do ij = 1,ntq
               do ii = 1,ntq
                 write(ifse2,*)ii,ij,zsec2(ii,ij,ip,iwini)
               enddo   
               enddo   
              end do ! ip                         
              if (ixc.eq.1) then
               if (is.eq.1) then
                 isx= iclose('SEX2U')
               else
                 isx= iclose('SEX2D')
               endif
              else
               if (is.eq.1) then
                 isx= iclose('SEXcore2U')
               else
                 isx= iclose('SEXcore2D')
               endif
              endif 
            endif ! loffd
         endif ! Is_IO_Root_RSMPI()
c-------------------------
        elseif(ixc==2) then
c write SEc(qt,w)
c        write (ifsec(is),*)' '
c        write (ifsec(is),6500) is,nbloch,deltaw,ef
c        do      ip = 1,nq
c          call wsec2 (ifsec(is),is,q(1,ip),itq(1:ntq),
c     i              eqx(1:ntq,ip,is),deltaw,nbloch,
c     i              dreal(zsec(:,:,ip)), dimag(zsec(:,:,ip)), 
c     d              3,ntq,nband) 
c        end do
c
c-------------------------
          if (Is_IO_Root_RSMPI()) then
            ifoutsec=ifsec(is)
            do ip = 1,nq
              do i  = 1,ntq
                dscdw1  = (dreal(zsec(1,i,ip)) - dreal(zsec(0,i,ip)))
     &               /deltaw
                zfac1   = 1d0/(1d0-dscdw1)
                dscdw2  = (dreal(zsec(0,i,ip)) - dreal(zsec(-1,i,ip)))
     &               /deltaw
                zfac2   = 1d0/(1d0-dscdw2)
                dscdw  = (dreal(zsec(1,i,ip)) - dreal(zsec(-1,i,ip)))
     &               /(2d0*deltaw)
                zfac   = 1d0/(1d0-dscdw)
                write(6,"(' j iq isp=' i3,i4,i2,'  q=',3f8.4,
     &               '  eig=',f8.4,'  Re(Sc) =',3f8.4,'  Img(Sc) =',3f8.4,
     &               '  zfac=',f8.4,'  zfac1&2 =',2f8.4)") 
     &               itq(i),ip,is, q(1:3,ip), eqx(i,ip,is), 
     &               hartree*dreal(zsec(-1:1,i,ip)), 
     &               hartree*dimag(zsec(-1:1,i,ip)),zfac,zfac1,zfac2
                write(ifoutsec,"(3i5,3d24.16,3x,d24.16,3x,3d24.16,
     &               3x,3d24.16,3x,3d24.16)")
     &               itq(i),ip,is, q(1:3,ip), eqx(i,ip,is), 
     &               hartree*dreal(zsec(-1:1,i,ip)), 
     &               hartree*dimag(zsec(-1:1,i,ip)),zfac,zfac1,zfac2
ccccccccccccccccccccc
                if(cohtest) then
                  write(671,"( i3,i4,i2, 3f8.4,
     &                 '  eig=',f8.4,'  coh =',f10.4)") 
     &                 itq(i),ip,is, q(1:3,ip), eqx(i,ip,is), 
     &                 hartree*coh(i,ip)
                endif
ccccccccccccccccccccc
              end do            ! do i
            end do              !do ip
ctm, 070621
            if (loffd) then
              if (lwang) then
               if (is.eq.1) then
                 ifse2= iopen('SEC2U.wan',1,-1,0)
               else
                 ifse2= iopen('SEC2D.wan',1,-1,0)
               endif
              else ! lwang   
              if (is.eq.1) then
                ifse2= iopen('SEC2U',1,-1,0)
              else
                ifse2= iopen('SEC2D',1,-1,0)
              endif
              endif ! lwang
              write(ifse2,"(4i5,f12.6)")ntq,nq,nqbz,nspin,ef
              write(ifse2,*)q
              write(ifse2,*)itq
              write(ifse2,*)eqx
              do ip = 1,nq
c               write(ifse2,*)zsec2(1:ntq,1:ntq,ip)
               write(ifse2,*)'ip=',ip
               do ij=1,ntq
               do ii=1,ntq
                write(ifse2,*)ii,ij,zsec2(ii,ij,ip,0)
               enddo
               enddo
              end do ! ip                         
              if (lwang) then
               if (is.eq.1) then
                 isx= iclose('SEC2U.wan')
               else
                 isx= iclose('SEC2D.wan')
               endif
              else ! lwang   
              if (is.eq.1) then
                isx= iclose('SEC2U')
              else
                isx= iclose('SEC2D')
              endif
              endif ! lwang
            endif ! loffd
          endif ! Is_IO_Root_RSMPI()
c-------------------------
        elseif(ixc==4) then     !spectrum mode
          ifoutsec=9300
          if(is==1) sss='.UP'
          if(is==2) sss='.DN'
          if (Is_IO_Root_RSMPI()) then
            if (lwang) then
              open(ifoutsec,file='SEComg'//sss//'.wan') 
            else ! lwang   
            open(ifoutsec,file='SEComg'//sss) 
            endif ! lwang
            write(ifoutsec,*)'*** iwini,iwend,ntq,nq,dwplot,ef'
            write(ifoutsec,*)iwini,iwend,ntq,nq,dwplot,ef
            do ip = 1,nq
              do i  = 1,ntq
                write(6,"(' --- j iq isp=' i3,i4,i2,' q=',3f8.4,' eig=',f8.4)") 
     &               itq(i),ip,is, q(1:3,ip), eqx(i,ip,is) 

                do iw = iwini,iwend
                  if(mod(iw,50)==1) 
     &                 write(6,"(' omega-ef=',f9.4,'  Sc=',2f9.4)") 
     &                 (omega(i,iw)-ef)*rydberg(), hartree*zsec(iw,i,ip)

                  write(ifoutsec,"(4i5,3f10.6,3x,f10.6,2x,f16.8,x,3f16.8)")
     &                 iw,itq(i),ip,is, q(1:3,ip),  eqx(i,ip,is),
     &                 (omega(i,iw)-ef)*rydberg(),  hartree*zsec(iw,i,ip) !,sumimg
                end do
                write(6,*) 
              end do
            end do
            close(ifoutsec)  
ctm, 080521                                                                            
            if (loffd) then
              if (lwang) then
               if (is.eq.1) then
                 ifse2= iopen('SECW2U.wan',1,-1,0)
               else
                 ifse2= iopen('SECW2D.wan',1,-1,0)
               endif
              else ! lwang
              if (is.eq.1) then
                 ifse2= iopen('SECW2U',1,-1,0)
              else
                 ifse2= iopen('SECW2D',1,-1,0)
              endif
              endif ! lwang
              write(ifse2,"(4i5,2f12.6,2i7)")ntq,nq,nqbz,nspin,ef,
     &         dwplot,iwini,iwend
              write(ifse2,*)q
              write(ifse2,*)itq
              write(ifse2,*)eqx
              write(ifse2,*)zsec2
              if (lwang) then
               if (is.eq.1) then
                 isx= iclose('SECW2U.wan')
               else
                 isx= iclose('SECW2D.wan')
               endif
              else ! lwang   
              if (is.eq.1) then
                 isx= iclose('SECW2U')
              else
                 isx= iclose('SECW2D')
              endif
              endif ! lwang
            endif ! loffd
          endif ! Is_IO_Root_RSMPI
        endif ! ixc
 2000 continue                  !end of spin-loop


C... close files
c      isx = iclose ('wc.d')
c      isx = iclose ('wci.d')
      isx = iclose ('hbe.d')


C --- EXspectrum ----------------------------------------------------
c This section is similar with efsimplef.f
c RS: this mode is disabled
      if(sum(ifexsp(1:nspin))/=0) then
        call RSMPI_Stop("sum(ifexsp)/=0. This cannot happen")
      endif
c-------------------------------------------------------------------------------
 3000 continue
      call cputid(ifile_rsmpi)
      call cputid(0)

c      if(ixc==1.and.sergeys)stop " OK! hsfp0: Sergey's Exchange mode"
c      if(ixc==2.and.sergeys)stop " OK! hsfp0: Sergey's Correlation mode"

      call RSMPI_Finalize()
      write(ifile_rsmpi,*) "RSMPI_Finalize is called."

      if (Is_IO_Root_RSMPI()) then
         if(ixc==1) then
            write(0,*) ' OK! hsfp0: Exchange mode'
            write(6,*) ' OK! hsfp0: Exchange mode'
         elseif(ixc==2) then
            write(0,*) ' OK! hsfp0: Correlation mode'
            write(6,*) ' OK! hsfp0: Correlation mode'
         elseif(ixc==3 ) then
            write(0,*) ' OK! hsfp0: Core-exchange mode'
            write(6,*) ' OK! hsfp0: Core-exchange mode'
         elseif(ixc==4 ) then
            write(0,*) ' OK! hsfp0: spectrum mode'
            write(6,*) ' OK! hsfp0: spectrum mode'
         else
            write(0,*) "This cannot happen!"
         endif
      endif

      end

