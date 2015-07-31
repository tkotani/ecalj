      program hsfp0
      use m_readqg,only: readqg,readngmx
      use m_readeigen,only: init_readeigen,readeval,lowesteval,init_readeigen2
      use m_read_bzdata,only: read_bzdata,
     &  ngrp2=>ngrp,nqbz,nqibz,nqbzw,nteti,ntetf,n1,n2,n3,qbas,ginv,qbasmc,
     &  dq_bzcase2,qbz,wbz,qibz,wibz,qbzw,
     &  idtetf,ib1bz,idteti,
     &  nstar,irk,nstbz
      use m_genallcf_v3,only: genallcf_v3,
     &  nclass,natom,nspin,nl,nn,nnv,nnc, ngrp,
     &  nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, nctot,niw,nw,
     &  alat,ef, diw,dw,delta,deltaw,esmr,symgrp,clabl,iclass,nindxv,nindxc,ncwf,
     &  invg, il, in, im, ilnm, nlnm, ilv,inv,imv,  ilnmv, nlnmv,
     &  ilc,inc,imc,  ilnmc, nlnmc,nindx,konf,icore,ncore,
     &  occv,unoccv ,occc,unoccc, nocc, nunocc, plat, pos,z,ecore, freq, symgg ! symgg=w(igrp)
      use keyvalue,only: getkeyvalue

      use m_rdpp,only: rdpp, !"call rdpp" generate following data set.
     &   nblocha,lx,nx,ppbrd,mdimx,nbloch,cgr

      use m_zmel,only: !these data are stored in this program.
     &   nband,itq,imdim,ngcmx,ngpmx,
     &   miat,tiat,shtvg, ntq, ppbir

      use m_mpi
      implicit none
!! = Calculate the  self-energy \Sigma within the GW approximation. Some other functions =
!!
c     SEx(q,t,t) = <psi(q,t) |SEx| psi(q,t)>
c     SEc(q,t,t) = <psi(q,t) |SEc| psi(q,t)>
c
c     where SEc(r,r';w) = (i/2pi) < [w'=-inf,inf] G(r,r';w+w') Wc(r,r';w') >
c
c     This routine requirs an input from standard IO.
c     In scripts, you can do it like, prompt>echo mode|../exec/hsfp0 >lsc,
c     where mode is 1 or 2.
c
c     This routine requirs a number from standard Input. ixc below.
c     Look into script calling this routine ecalj/fpgw/exec/gw_lmfh. E.g, this gwsc calls echo 12|../exec/hsfp0 >lsc when mode=12
c
c     mode= 11: exchange    mode SEx, the exchange part of the self-energy
c     mode= 12: correlation mode SEc, the correlated part of the self-energy
c     mode= 3: core exchange mode SEXcore
c     mode= 4: plot spectrum function ---See manual
c     mode= 5: exchange    mode Exx, the exchange part of the total energy
c     mode= 6: core exchange mode Exxc, Exx(core-valence)
c
c     May 2002 Takashi Miyake. Total energy calc.
c     This hsfp0 is build from hsec10.f by F.Aryasetiawan.
!! ---------------------------------------------------------------

c------------------------------------------------------------
c     real(8),parameter ::  ua = 1d0    ! constant in w(0)exp(-ua^2*w'^2) to take care of peak around w'=0
c------------------------------------
c     test switches to calculate the self-energy based on an another separation of \Sigma.
!     \Sigma = \Sigma_{sx} + \Sigma_{coh} + \Sigma_{img axis} + \Sigma_{pole} by Hedin PR(1965)A785
!     I found COH term has inevitably poor accuracy.
      logical ::tetra, tetra_hsfp0,         
     &     screen  = .false.,   ! \Sigma_{sx} for mode 1 and
      ! \Sigma_{img axis} + \Sigma_{pole} for mode 2
     &     cohtest = .false.,   ! \Sigma_{coh}. mode swich is not required.
     &     tetraex = .false.    ! This switch is only meaningful for mode=1,5,6
      ! If you want to calculate exchange, use tetraex=T .
      ! Note that you have to supply EFERMI by the tetrahedon method.
c------------------------------------
      real(8):: esmr2,shtw
      integer::
     &     ixc,iopen,ifhbed, nprecb,mrecb,mrece,nlmtot,nqbzt, !nband,
     &     ibas,ibasx,nxx,ifqpnt,ifwd, !ngpmx,ngcmx,nbloch,
     &     nprecx,mrecl,nblochpmx2,nwp=0,niwt, nqnum,nblochpmx, !mdimx,
!     nwp=0 for initialization March26-2007
     &     noccxv,maxocc2,noccx,ifvcfpout,iqall,iaf, !ntq, !ifrcw,ifrcwi,  
     &     i,k,nspinmx, nq,is,ip,iq,idxk,ifoutsex,iclose,nq0i,ig,
     &     mxkp,nqibzxx,ntet,nene,iqi, ix,iw,
     &     nlnx4,niwx,irot,invr,invrot,ivsum, ifoutsec,
     &     ifsec(2),             !,ifcphi
     &     ifxc(2),ifsex(2), ifphiv(2),ifphic(2),ifec,ifexsp(2),
     &     ifsecomg(2),ifexx,ndble=8
      real(8) :: pi,tpia,vol,voltot,rs,alpha,
     &     qfermi,efx,valn,efnew,edummy,efz,qm,xsex,egex,edummyd(1),
     &     zfac1,zfac2,dscdw1,dscdw2,dscdw,zfac,ef2=1d99,exx,exxq,exxelgas
      logical lqall,laf
c      integer,allocatable :: itq(:)
      real(8),allocatable    :: q(:,:)
      integer,allocatable :: ngvecp(:,:), ngvecc(:,:),iqib(:), kount(:,:) 
        !lx(:) !nx(:,:),nblocha(:), !ngvecpB(:,:,:), 
      real(8),allocatable:: vxcfp(:,:,:),
     &     wqt(:), wgt0(:,:),q0i(:,:),
     &     eqt(:), !ppbrd (:,:,:,:,:,:,:),cgr(:,:,:,:),
     &     ppbrdx(:,:,:,:,:,:,:),aaa(:,:), 
     &     ppb(:), eq(:),       
     &     eqx(:,:,:),eqx0(:,:,:),ekc(:),coh(:,:)
      complex(8),allocatable:: geigB(:,:,:,:) ,zsec(:,:,:)

      logical :: exchange, legas, tote
      real(8) ::  rydberg,hartree
      real(8):: qreal(3), ntot,nocctotg2,tripl,xxx(3,3)
      logical ::nocore

c     space group infermation
      integer,allocatable :: iclasst(:), invgx(:) !, miat(:,:)
c      real(8),allocatable    :: shtvg(:,:) !tiat(:,:,:),

c     tetra
      real(8),allocatable :: qz(:,:),qbzxx(:),wbzxx(:),wtet(:,:,:,:),
     &     eband(:,:,:), ene(:) !,ecore(:,:)
      integer(4),allocatable ::idtetx(:,:),idtet(:,:),ipq(:)
     &     ,iene(:,:,:),ibzx(:) ! ,nstar(:)

      integer ::ib,iqx,igp,iii,ivsumxxx,isx,iflegas, iqpntnum
c
      real(8),allocatable   :: eex1(:,:,:),exsp1(:,:,:),qqex1(:,:,:,:)
      integer,allocatable:: nspex(:,:),ieord(:),itex1(:,:,:)
      real(8)    :: qqex(1:3), eex,exsp,eee, exwgt,deltax0
      integer :: itmx,ipex,itpex,itex,nspexmx,nnex,isig,iex,ifexspx
     &     ,ifexspxx ,ifefsm, nq0ix,ifemesh,nz
      character(3)  :: charnum3,sss
      character(12) :: filenameex
      logical :: exspwrite=.false.
      character(8):: xt

      integer :: iwini,iwend
      real(8),allocatable:: omega(:,:)
      real(8) ::  omegamax,dwplot,omegamaxin
c     logical :: sergeys

      integer::nqbze,ini,nq0it,idummy
c      real(8),allocatable:: qbze(:,:)

      real(8)   :: ebmx(2)
      integer:: nbmx(2)

      real(8):: volwgt

      integer::nwin, incwfin
      real(8)::efin,ddw
c      integer,allocatable::imdim(:)
      real(8),allocatable::freqx(:),freqw(:),wwx(:) !,expa(:)

c      logical:: GaussSmear=.true.      !readgwinput,
      integer::ret
      character*(150):: ddd


      integer:: bzcase=1,  ngpn1,mrecg,verbose,ngcn1,nwxx
      real(8)   :: wgtq0p,quu(3)

      real(8),allocatable:: freq_r(:)

      logical ::smbasis
      integer:: ifpomat,nkpo,nnmx,nomx,ikpo,nn_,no,nss(2)
      real(8):: q_r(3)
      real(8),allocatable:: qrr(:,:)
      integer,allocatable:: nnr(:),nor(:)
      complex(8),allocatable:: pomatr(:,:,:),pomat(:,:)

      real(8)::sumimg
      logical allq0i            !S.F.Jan06

      integer:: nw_i

c     logical :: imgonly=.false.
      complex(8):: zseciip(-1:1)
      logical :: debug=.false.

!! newaniso mode
      real(8),allocatable:: vcousq(:),dmlx(:,:),epinvq0i(:,:),wklm(:),vcoud(:)
      complex(8),allocatable:: zcousq(:,:)
c      logical:: newaniso=.true.
      integer:: ifvcoud,lxklm,ifidmlx,iqq

      integer,allocatable:: irkip_all(:,:,:,:),irkip(:,:,:,:)
      integer:: irank,nrank

      integer:: ififr,ifile_handle
      integer:: timevalues(8)  ,isp,dest
      character(128) :: ixcc

C...  mode switch. --------------
      call MPI__Initialize()      ! MIZUHO-IR
      call date_and_time(values=timevalues)
      write(6,'(a,9i5)')'dateandtime1=',mpi__rank,timevalues(1:8)

      hartree=2d0*rydberg()
      if(cohtest) then
        screen = .true.
        ixc = 2; nz=0
        open(671,file='COH')
      elseif(mpi__root) then
        write(6,*) ' --- Choose omodes below ----------------'
        write(6,*) '  Sx(1) Sc(2) ScoreX(3) Spectrum(4) '
        write(6,*) '  EXX_val-val (5)  Exx_core-val(6) '
        write(6,*) '  Sx_sf(11) Sc_sf(12) '
        write(6,*) '  [option --- (+ QPNT.{number} ?)] '
        write(6,*) ' --- Put number above ! -----------------'
        call readin5(ixc,nz,idummy)
      endif
      call MPI__Broadcast(ixc)
      call MPI__Broadcast(nz)
      if(mpi__root) call headver('hsfp0',ixc)
      write(ixcc,"('.mode=',i4.4)")ixc
      call MPI__consoleout('hsfp0'//trim(ixcc))
      write(6,*) ' ixc nz=',ixc, nz
      if(ixc==0) call rx( ' --- ixc=0 --- Choose computational mode!')

!!  tetraex is only for Ex.
      tetraex=tetra_hsfp0()
      iii=verbose()
      write(6,*)' verbose=',iii
      if(ixc==1.or.ixc==5.or.ixc==6) then;
      else; tetraex=.false.
      endif
!! These are mainly used now.
      if(ixc==11.or.ixc==12) then
        ixc=ixc-10
      endif

!! ===  readin BZDATA. See gwsrc/rwbzdata.f ===
!! See use m_read_bzdata,only: at the top of this routine
      call read_BZDATA()
      write(6,*)' nqbz  =',nqbz
      write(6,*)' nqibz ngrp=',nqibz,ngrp
      call pshprt(60)

!! === readin GWIN and LMTO, then allocate and set datas. ===
!! See use m_genallcf_v3,only: at the top of this routine
      nwin =-999                !not readin NW file
      efin =-999d0              !not readin EFERMI
      if    (ixc==3) then;  incwfin= -2 !core exchange mode
      elseif(ixc==5) then;  incwfin= -4 ! valence-valence Ex energy mode. ! See rgwinf called from genallcf_*
      elseif(ixc==6) then;  incwfin= -3 ! core-valence Ex energy mode.
      else               ;  incwfin= -1 !use 7th colmn for core at the end section of GWIN
      endif
      call genallcf_v3(nwin,efin,incwfin) ! module m_genallcf_v3. See use m_genallcf in this rouitine
      if(ngrp/= ngrp2) call rx( 'ngrp inconsistent: BZDATA and LMTO GWIN_V2')

!! Get maximums 
      call getnemx8(nbmx,ebmx)
!!     nbmx1 ebmx1: to set how many bands of <i|sigma|j>  do you calculate.
!!     nbmx2 ebmx2: to restrict num of bands of G to calculate G \times W
!! ebmx2 nbmx2 are dummy now.
      nbmx(2)=9999999
      ebmx(2)=1d10
      write(6,"('  nbmx ebmx from GWinput=',2i8,2d13.5)") nbmx,ebmx

!!!!  WE ASSUME iclass(iatom)= iatom (because of historical reason)
      if (nclass /= natom ) call rx( ' hsfp0: nclass /= natom ')
      call pshprt(30)
      pi   = 4d0*datan(1d0)
      tpia = 2d0*pi/alat
      call dinv33(plat,1,xxx,vol)
      voltot = dabs(vol)*(alat**3)
      shtw = 0d0
!! if(esmr<1d-5) shtw=0.01d0 ! Ferdi's shift to avoid resonance effect(maybe), I used this until sep2012

!! ef is taken as rs for the empty-sphere test case of legas=T case 
!! HOMOGENIOUS GAS code. Usually not used.Need fixing if necessary.
      legas = .false.
      INQUIRE (FILE = 'LEGAS', EXIST = legas)
      if(legas) then            !!! test for electron gas case.
        write(6,*)' find LEGAS. legas =',legas
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
        if(tetraex) call rx( 'legas You have to give ef of  tetrahedron')
      endif

!!
      ifexsp=0
      if(ixc==1) then
        exchange=.true.
        tote=.false.
        write(6,*) ' --- Exchange mode --- '
        if(mpi__root) then
          ifxc(1)  = iopen('XCU'//xt(nz),1,-1,0)
          ifsex(1) = iopen('SEXU'//xt(nz),1,-1,0)
c          ifsex2(1)= iopen('SEX2U',0,-1,0) !out SEX_nn'
          if (nspin == 2) then
            ifxc(2)  = iopen('XCD'//xt(nz),1,-1,0)
            ifsex(2) = iopen('SEXD'//xt(nz),1,-1,0)
c            ifsex2(2)= iopen('SEX2D',0,-1,0) !out SEX_nn'
          endif
          INQUIRE (FILE = 'EXspTEST', EXIST = exspwrite)
          if(exspwrite) then
            write(6,*)'--- Find EXspTEST ExspectrumWrite=',exspwrite
            write(6,*)'--- esmr is chosen to be 2d0 Ry'
            esmr= 2d0
            do is=1,nspin
              ifexsp(is)  = iopen('EXSP.'//char(48+is),1,-1,0)
            enddo
          endif
        endif
      elseif(ixc==2) then
        exchange=.false.
        tote=.false.
        write(6,*) ' --- Correlation mode --- '
        if(cohtest) write(6,*) ' COH calculation mode. Results in COH'
         if(mpi__root) then
           ifsec(1)   = iopen('SECU'//xt(nz),1,-1,0) ! output files
           if (nspin == 2) 
     .      ifsec(2)   = iopen('SECD'//xt(nz),1,-1,0)
         endif  
      elseif(ixc==3) then
        exchange=.true.
        tote=.false.
        esmr=0d0
        write(6,*) ' --- CORE Exchange mode --- '
        if(mpi__root) then
          ifsex(1)   = iopen('SEXcoreU'//xt(nz),1,-1,0)
          if (nspin == 2) then
            ifsex(2)   = iopen('SEXcoreD'//xt(nz),1,-1,0)
          endif
        endif  
!! ixc=4,5,6 not checked now... NEED MPI and so on...
      elseif(ixc==4) then
        write(6,*) ' --- Spectrum function Sigma(omega) mode --- '
        exchange=.false.
        tote=.false.
c     ifsecomg(1) = iopen('SEComgU'//xt(nz),1,-1,0) ! output files
c     if (nspin == 2)
c     .  ifsecomg(2) = iopen('SEComgD'//xt(nz),1,-1,0)
        call rx('need to test this mode again')
      elseif(ixc==5) then
        exchange=.true.
        tote=.true.
        write(6,*) ' --- Exx mode valence-valence --- '
        ifexx = iopen('TEEXXvv',1,-1,0)
        esmr2 = esmr           !takao June2002
        call rx('need to test this mode again')
      elseif(ixc==6) then
        exchange=.true.
        tote=.true.
        esmr2 = esmr           !takao June2002
        esmr  = 0d0            !takao June2002
        write(6,*) ' --- CORE Exx mode core-valence --- '
        ifexx = iopen('TEEXXcv',1,-1,0)
        call rx('need to test this mode again')
      else
        call rx( ' hsfp0: Need input (std input) 1-6!')
      endif

c---  Neglect core is NoCore exists ------------
c     inquire(file='NoCore',exist=nocore)
c     if(nocore) nctot=0

c---
      write(6, *) ' --- computational conditions --- '
      write(6,'("    deltaw  =",f13.6)') deltaw
c     write(6,'("    ua      =",f13.6)') ua
      write(6,'("    esmr    =",f13.6)') esmr
      if(ixc==6) write(6,'("    esmr2    =",f13.6)') esmr2
      write(6,'("    alat voltot =",2f13.6)') alat, voltot
      write(6,'("    niw dw   =",i7,f13.6)') niw,dw

!!   read dimensions of wc,b,hb from hbe.d
      ifhbed     = iopen('hbe.d',1,0,0)
      read (ifhbed,*) nprecb,mrecb,mrece,nlmtot,nqbzt, nband,mrecg
      isx = iclose ('hbe.d')
      if (nprecb == 4) call rx( 'hsfp0: b,hb in single precision')
      call init_readeigen(ginv,nspin,nband,mrece) !initialization of readEigen

!!  === get space group information ===
!!  True class information in order to determine the space group -----------
!!  because the class in the generated GW file is dummy.(iclass(ibas)=ibas should be kept).
      open (102,file='CLASS')
      allocate(iclasst(natom),invgx(ngrp)
     &     ,miat(natom,ngrp),tiat(3,natom,ngrp),shtvg(3,ngrp))
      write(6,*)'  --- Readingin CLASS info ---'
      do ibas = 1,natom
        read(102,*) ibasx, iclasst(ibas)
        write(6, "(2i10)") ibasx, iclasst(ibas)
      enddo
!! ===  Get space-group transformation information. See header of mptaouof. ===
      call mptauof(symgg,ngrp,plat,natom,pos,iclasst
     o     ,miat,tiat,invgx,shtvg )

!! ====  Get array size to call rdpp ====
      call getsrdpp2( nclass,nl,nxx)
      call readngmx('QGpsi',ngpmx)
      call readngmx('QGcou',ngcmx)
      write(6,*)' max number of G for QGpsi and QGcou: ngcmx ngpmx=',ngcmx,ngpmx
      allocate(ngvecp(3,ngpmx),ngvecc(3,ngcmx))
      call readqg('QGpsi',qibz(1:3,1),ginv, quu,ngpn1, ngvecp)
      call readqg('QGcou',qibz(1:3,1),ginv, quu,ngcn1, ngvecc)
      deallocate(ngvecp,ngvecc)
      if(debug) write(6,*) ' end of read QGcou'

c      allocate( nx(0:2*(nl-1),nclass), nblocha(nclass) ,lx(nclass),
c     &     ppbrd ( 0:nl-1, nn, 0:nl-1,nn, 0:2*(nl-1),nxx, nspin*nclass),
c     &     cgr(nl**2,nl**2,(2*nl-1)**2,ngrp))

!! === readin plane wave parts, and Radial integrals ppbrd. ===
!     ppbrd = radial integrals
!     cgr   = rotated cg coeffecients.
!     geigB = eigenfunction's coefficiens for planewave.
!     ngvecpB (in 1stBZ) contains G vector for eigen function.
!     ngveccB (in IBZ)   contains G vector for Coulomb matrix.
c      call rdpp_v3(nxx, nl, ngrp, nn, nclass, nspin, symgg,qbas,
c     o     nblocha, lx, nx,  ppbrd , mdimx, nbloch, cgr)
      call rdpp(nxx, nl, ngrp, nn, nclass, nspin, symgg,qbas)
      !     output ==>     nblocha, lx, nx,  ppbrd , mdimx, nbloch, cgr
      call pshprt(60)

!! Readin WV.d
      if(.not.exchange.or.(exchange.and.screen)) then !screen means screened exchange case
        ifwd = iopen('WV.d',1,-1,0) !direct access files WVR and WVI which include W-V.
        read (ifwd,*) nprecx,mrecl,nblochpmx,nwp,niwt, nqnum, nw_i
        write(6,"(' Readin WV.d =', 10i5)") nprecx,mrecl,nblochpmx, nwp, niwt, nqnum, nw_i
        ifwd =iclose('WV.d')
        call checkeq(nprecx,ndble)
        nw=nwp-1
        if (niwt /= niw) call rx( 'hsfp0: wrong niw')
!! Energy mesh; along real axis. Read 'freq_r' 
!! NOTE nw_i=nw for non-timereversal case.
!!      nw_i=0 for time-reversal case.
!!  NOTE: We assume freq_r(i) == -freq_r(-i) in this code. feb2006
!!  NOTE: this program assumes freq_r(iw)=freq_r(-iw). freq_r(iw <0) is redundant.
        ififr=ifile_handle()
        open(UNIT=ififr,file='freq_r')
        read(ififr,*)nwxx       !number of energy points
        if(nwxx/= nw+1) call rx( ' freq_r nw /=nw')
        allocate(freq_r(nw_i:nw)) !freq_r(1)=0d0
        do iw = nw_i,nw
          read(ififr,*)freq_r(iw)
        enddo
        close(ififr)
        if(nw_i/=0) then
          if(nw/= -nw_i)        call rx( "sxcf_fal3: nw/=-nw_i")
          if(freq_r(0)/=0d0)    call rx( "sxcf_fal3: freq_r(0)/=0")
          if( sum(abs( freq_r(1:nw)+freq_r(-1:-nw:-1)))/=0) 
     &     call rx( "sxcf_fal3_scz: freq_r /= -freq_r")
        endif
      endif

      if(tetraex) goto 201

!!== Determine Fermi energy ef for given valn (legas case), or corresponding charge given by z and konf.==
!!    When esmr is negative, esmr is geven automatically by efsimplef.
      if((ixc==3).or.(ixc==6)) then
        ef = lowesteval() -1d-3 !lowesteigen(nspin,nband,qbz,nqbz) - 1d-3 !lowesteb was
        write(*,*)'efffffffffff=',ef
        if(maxval(ecore(:,1:nspin))>ef) then
          call rx( 'hsfp0 ixc=3/6:  ecore>evalence ')
        endif
        if(ixc==6) then
          call efsimplef2a(nspin,wibz,qibz,ginv, 
     i           nband,nqibz
     i           ,konf,z,nl,natom,iclass,nclass
     i           ,valn, legas, esmr2, !!! valn is input for legas=T, output otherwise.
c     
     i           qbz,nqbz       !index_qbz, n_index_qbz,
     o           ,efnew)
          ef2 = efnew
          write(6,*)' end of efsimple ef2 esmr2=',ef2,esmr2
        endif
      else                      ! if(esmr/=0d0) then
c     --- determine Fermi energy ef for given valn (legas case) or corresponding charge given by z and konf.
!     When esmr is negative, esmr is geven automatically by efsimplef.
        call efsimplef2a(nspin,wibz,qibz,ginv, 
     i        nband,nqibz
     i        ,konf,z,nl,natom,iclass,nclass
     i        ,valn, legas, esmr, !!! valn is input for legas=T, output otherwise.
     i        qbz,nqbz          !index_qbz, n_index_qbz,
     o        ,efnew)
        ef = efnew
        if (ixc==5) ef2 = ef
c     - check total ele number -------
        ntot  = nocctotg2(nspin, ef,esmr, qbz,wbz, nband,nqbz) !wbz
        write(6,*)' ef    =',ef
        write(6,*)' esmr  =',esmr
        write(6,*)' valn  =',valn
        write(6,*)' ntot  =',ntot
      endif
!!
 201  continue

      if(tote) then
        ddw= .5d0
c        if(GaussSmear) ddw= 10d0
        ddw= 10d0
        if(ixc==6) then
c     noccxv = maxocc (ifev,nspin, ef2+ ddw*esmr2, nband,nqbz)
          noccxv = maxocc2 (nspin, ef2+ ddw*esmr2, nband,qbz,nqbz)
        else
          noccxv = maxocc2 (nspin, ef + ddw*esmr,  nband, qbz,nqbz)
        endif
      endif
c     noccxv = maxocc (ifev,nspin, ef+0.5d0*esmr, nband,nqbz)  ! maximum no. of occupied valence states
!     maxocc seems to give (the maxmum number of occ + 1).

      call init_readeigen2(mrecb,nlmto,mrecg) !initialize m_readeigen

      if(.not.tote) then
        if(nz==0) then
          call getkeyvalue("GWinput","<QPNT>",unit=ifqpnt,status=ret)
        else
          ifqpnt  = iopen('QPNT'//xt(nz),1,0,0)
        endif
      endif
c      write(6,*)' ifqpnt ret=',ifqpnt,ret

!!  read q-points and states
      if(ixc==4) then
        call readxx(ifqpnt)    !skip to ***** for q point for spectrum function.
        omegamaxin = 1d70
        read (ifqpnt,*,err=2038,end=2038) dwplot,omegamaxin
 2038   continue
        omegamax = 2*freq_r(nw) !This is in Ry.
        if(omegamaxin < omegamax) then
          write(6,*)' --Use readin dwplot and omegamaxin from <QPNT>'
          omegamax = omegamaxin
        endif
        if( omegamax <0) call rx( 'hsfp0 :strange omegamax <0 ')
        iwini =  -int( omegamax / dwplot )
        iwend =   int( omegamax/  dwplot )
        write(6,*)' iwini:iwend omegamax(Ry)=',iwini,iwend,omegamax
      else
        iwini = -1
        iwend = 1
      endif

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

!     states
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

!     q-points
      if (lqall) then           !all q-points case
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
c
      nspinmx = nspin
      if (laf) nspinmx =1
      close(ifqpnt)
!
      call MPI__Broadcast(nq)
      if(mpi__root) then
        do dest=1,mpi__size-1
          call MPI__REAL8send(q,3*nq,dest)
        enddo
      else  
        call MPI__REAL8recv(q,3*nq,0)
      endif  

C...
      call getkeyvalue("GWinput","QPNT_nbandrange",nss,2 !nov2005
     &     ,default=(/-99997,-99997/) )
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

!! read LDA eigenvalues
      allocate(eqx(ntq,nq,nspin),eqx0(ntq,nq,nspin),eqt(nband))
      do      is = 1,nspin
        do      ip = 1,nq
          call readeval(q(1,ip),is,eqt)
          eqx0(1:ntq,ip,is) = eqt(itq(1:ntq))
          eqx (1:ntq,ip,is) = rydberg()*(eqt(itq(1:ntq))- ef)
        enddo
      enddo
      deallocate(eqt)
      if(tote) then
        call winfo2(6,nspin,nq,ntq,is,nbloch,ngpn1,ngcn1,nqbz,nqibz,ef,ef2,deltaw,alat,esmr,esmr2)
      else
        call winfo(6,nspin,nq,ntq,is,nbloch,ngpn1,ngcn1,nqbz,nqibz,ef,deltaw,alat,esmr)
      endif
!!-------------------------
!!     LDA exchange-correlation
!!-------------------------
      if(ixc==1) then
        allocate(  vxcfp(ntq,nq,nspin) )
        call rsexx2(nspin,itq,q,ntq,nq, ginv, symgg,ngrp, vxcfp) !add ginv july2011
c     loop over spins
        if(mpi__root) then
        do is = 1,nspinmx
          write (ifxc(is),*) '==================================='
          write (ifxc(is),"(' LDA exchange-correlation : is=',i3)")is
          write (ifxc(is),*) '==================================='
          call winfo(ifxc(is),nspin,nq,ntq,is,nbloch
     &           ,ngpn1,ngcn1,nqbz,nqibz,ef,deltaw,alat,esmr)
          write (ifxc(is),*)' ***'
          write (ifxc(is),"(a)") ' jband   iq ispin                  
     &           qvec 
     &           eigen-Ef (in eV)
     &           LDA XC (in eV)'
          ifoutsex = ifxc(is)
          write(6,*)
          do ip = 1,nq
            do i  = 1,ntq
              write(ifoutsex,"(3i5,3d24.16,3x,d24.16,3x,d24.16)")
     &                 itq(i),ip,is, q(1:3,ip), eqx(i,ip,is), 
     &                 vxcfp(i,ip,is)
              write(6,"(' j iq isp=' i3,i4,i2,'  q=',3f8.4,
     &                 '  eig=',f10.4,'  Sxc(LDA)=',f10.4)") 
     &                 itq(i),ip,is, q(1:3,ip), eqx(i,ip,is),
     &                 vxcfp(i,ip,is)
            end do
          end do
c     end of spin-loop
          if(is==1) isx = iclose('XCU'//xt(nz))
          if(is==2) isx = iclose('XCD'//xt(nz))
        enddo
        endif !mpi__root
        deallocate(vxcfp)
      endif

!!  Offset Gamma point Q0P
      write(6,*) 'reading QOP'
      open (101,file='Q0P')
      read (101,"(i5)") nq0i
      if(.not.exchange) call checkeq(nqibz+nq0i-1, nqnum)
      write(6,*) ' *** nqibz nq0i_total=', nqibz,nq0i
      nq0it = nq0i
      allocate( wqt(1:nq0i),q0i(1:3,1:nq0i) )
c     read (101,"(d24.16,3x, 3d24.16)" )( wqt(i),q0i(1:3,i),i=1,nq0i)
      nq0ix = nq0i
      do i=1,nq0i
        read (101,* ) wqt(i),q0i(1:3,i)
        if(wqt(i)==0d0 ) nq0ix = i-1
      enddo
      nq0i = nq0ix              ! New nq0i July 2001
      write(6,*) ' Used k number in Q0P =', nq0i
      write(6,"(i3,f14.6,2x, 3f14.6)" )(i, wqt(i),q0i(1:3,i),i=1,nq0i)
      close(101)
      allocate( wgt0(nq0i,ngrp) )
C     Sergey's 1stFeb2005
c     call q0iwgt2(symgg,ngrp,wqt,q0i,nq0i,
c     o            wgt0)
      call getkeyvalue("GWinput","allq0i",allq0i,default=.false.) !S.F.Jan06
      call q0iwgt3(allq0i,symgg,ngrp,wqt,q0i,nq0i, !S.F.Jan06
     o     wgt0)                ! added allq0i argument
C--------------------------
      if(nq0i/=0) write(6,*) ' *** tot num of q near 0   =', 1/wgt0(1,1)
      write(6,"('  sum(wgt0) from Q0P=',d14.6)")sum(wgt0)
c$$$      if(bzcase==2) then
c$$$        wgt0= wgt0*wgtq0p()/dble(nqbz)
c$$$        write(6,"('bzcase=2:  sum(wgt0_modified )=',d14.6)")sum(wgt0)
c$$$      endif

c$$$c     --- qbze(3,nqibze)
c$$$      nqbze  = nqbz *(1 + nq0it)
c$$$      allocate( qbze(3, nqbze) )
c$$$      call dcopy(3*nqbz, qbz, 1, qbze,1)
c$$$      do i = 1,nq0it
c$$$        ini = nqbz*(1 + i -1)
c$$$        do ix=1,nqbz
c$$$          qbze (:,ini+ix)   = q0i(:,i) + qbze(:,ix)
c$$$        enddo
c$$$      enddo

!! Pointer to optimal product basis
      allocate(imdim(natom))
      call indxmdm (nblocha,nclass,
     i     iclass,natom,
     o     imdim )
      if(niw/=0) then
!! Generate gaussian frequencies x between (0,1) and w=(1-x)/x
         allocate(freqx(niw),freqw(niw),wwx(niw)) !,expa(niw))
         call freq01x  (niw,    !ua,
     o        freqx,freqw,wwx)  !,expa)
      endif

c     ------ write energy mesh ----------
      ifemesh = iopen('emesh.hsfp0'//xt(nz),1,-1,0)
      deltax0 = 0d0
      if(mpi__root) then
      call writeemesh(ifemesh,freqw,niw,freq_r,nwp,deltax0)
      endif
      iii=ivsumxxx(irk,nqibz*ngrp)
      write(6,*) " sum of nonzero iirk=",iii, nqbz


!! === readin Vcoud and EPSwklm for newaniso()=T ===
c      if(newaniso) then
        ifidmlx = iopen('EPSwklm',0,0,0)
        read(ifidmlx) nq0ix,lxklm
        if(nq0i/=nq0ix) then
          write(6,*)'nq0i from EPSwklm /= nq0i',nq0i,nq0ix
          call rx( 'nq0i from EPSwklm /= nq0i')
        endif
        allocate( dmlx(nq0i,9))
        allocate( epinvq0i(nq0i,nq0i) )
        allocate( wklm((lxklm+1)**2))
        read(ifidmlx) dmlx, epinvq0i
        read(ifidmlx) wklm
        ifidmlx = iclose('EPSwklm') !jan2013 bugfix ifidmlx)
c      endif


c------------------------------------------------------------------
!! tetraex section. not checked well recentrly(sep2012) -------
!! this can be recovered in future.
      if(tetraex) then
        write(6,*) ' !!!!  tetraex=T is not tested well !!!!!!!!!'
        ifefsm  = iopen('EFERMI',1,-1,0)
        read(ifefsm,*) ef
        write(6,"(' ef is readin from EFERMI',d23.15)")ef
c     --- get tetrahedron
c     mxkp = n1*n2*n3
c     allocate( qbzxx(3*mxkp),wbzxx(mxkp),ipq(mxkp) )
c     call bzmesh (plat,qbasmc,n1,n2,n3,w(igrp),ngrp,ipq,
c     .               qbzxx,wbzxx,nqibzxx,mxkp)
c     allocate(idtetx(0:4,mxkp*6))
c     call tetirr(qbasmc,n1,n2,n3,ipq,nqibz,ntet,
c     .              idtetx)
c     allocate(idtet(0:4,ntet))
c     idtet(0:4,1:ntet) = idtetx(0:4,1:ntet)
c     deallocate(idtetx,qbzxx,wbzxx,ipq)

c     nene=0
c     nene = ntq*nq*nspin ! for energy points.
c     if(exchange) nene=0
c     allocate(wtet(nband,nspin,nqibz,0:3*nene),
c     &    eband(nband,nspin,nqibz), qz(3,nqibz),nstar(nqibz),
c     &    iene(3*ntq,nq,nspin), ene(0:3*nene) ) ! pointer for
        allocate(wtet(nband,nspin,nqibz,0:0),
     &        eband(nband,nspin,nqibz), qz(3,nqibz)) ! ,nstar(nqibz))
        call dcopy (3*nqibz,qibz,1,qz,1)
        do  is  = 1,nspin      !Readin eband
          do  iqi = 1,nqibz
c     iq    = idxk (qz(1:3,iqi),qbz,nqbz)
c     call rwdd1 (ifev(is), iq, nband, eband(:,is,iqi))
            call readeval(qz(1:3,iqi),is, eband(:,is,iqi))
          enddo
        enddo
c     wtet(nband,nsp,nqibz,iene) where
c     the energy pointer as iene(itp,ip,ispin) corresponding its energy value.
c     ene(0) = ef
c     if(.not.exchange) then
c     ix =0
c     do is = 1,nspin
c     do ip = 1,nq
c     do i  = 1,ntq
c     do iw = -1,1
c     ix  = ix+1
c     iene(3*i+iw-1,ip,is) = ix
c     ene(ix) = eqx0(i,ip,is) + 2.d0*(dble(iw)-shtw)*deltaw
c     enddo
c     enddo
c     enddo
c     enddo
c     endif
c     do ix = 0,3*nene
c     ene(ix) = ene(ix)-1d-15  ! to avoid coincidence
c     call bzints2(n1,n2,n3,eband,wtet(:,:,:,ix),nqibz,nband,nband,
c     .                nspin,edummy,edummy,edummy,1,ene(ix),2,ntet,idtet)
c     enddo
c     ene(ix) = ene(ix)-1d-15  ! to avoid coincidence
ccccccccccccccc
c     call bzints2(n1,n2,n3,eband,wtet(:,:,:,0),nqibz,nband,nband,
c     .               nspin,edummy,edummy,edummy,1,ef,2,ntet,idtet)
        volwgt = (3d0 - nspin) / ntetf ! ntetf was =6*n1*n2*n3
        call bzints2x(volwgt,eband,wtet(:,:,:,0),nqibz,nband,nband,
     .        nspin,edummy,edummy,edummyd,1,ef,2,nteti,idteti)
        ntot= sum(wtet)
        if(legas) then
          write(6,"(' tetra=T ef ntot nexact ratio=',15f12.6)") ef,ntot
     &           , ef**1.5d0/3d0/pi**2*voltot, ef**1.5d0 /3d0/pi**2*voltot/ntot
        else
          write(6,"(' tetra=T ef nvalence)=',15f12.6)") ef,ntot
        endif
        if(nspin==1) wtet = wtet/2d0
        do iqi = 1,nqibz
          wtet(:,:,iqi,:) = wtet(:,:,iqi,:)/nstar(iqi)
        enddo
        deallocate( eband,qz)  !, ene ) ! pointer for
!!     -- ibzx denote the index of k{FBZ for given k{1BZ.
        allocate(ibzx(nqbz))
        call invkibzx(irk,nqibz,ngrp,nqbz,
     o        ibzx)
        if (tote) then
          do i=1,nband
            if( sum(abs(wtet(i,:,:,0))) == 0d0 ) exit
          enddo
          if(  sum(abs(wtet(i+1:nband,:,:,0)))/=0d0) 
     &           call rx( ' hsfp0: wtetef sum err1')
          noccxv = i-1
        endif
      else
        allocate(wtet(1,1,1,1),ibzx(1)) !dummy
      endif
!!  end of tetra section --------------------------------------------
      iii=ivsumxxx(irk,nqibz*ngrp)
      write(6,*) " sum of nonzero iirk=",iii, nqbz


c-----------------------------------------------------------
c     calculate the the self-energy SEx(ip) or SEc(ip)
c-----------------------------------------------------------
!! == irkip control paralellization  ==
!! We have to distribute non-zero irkip_all into processes (nrank).
!! When irkip_all(nqibz,ngrp,nq,nspinmx)/=0, we expect grain-size 
!! for each job of (iqibz,igrp,iq,isp) is almost the same.
!! Our pupose is to calculate zsec(itp,itpp,iq).
!! Thus we need to set up communicator (grouping) MPI_COMM_iqisp(iq,isp) to do all_reduce.
!! (for given zsec(iq,isp), we take sum on zsec for (iqibz,igrp) by all_reduce.)
!! ---
!! NOTE: in future, we will further extend irkip for itp and itpp
      allocate(irkip_all(nspinmx,nqibz,ngrp,nq)) !this is global
      do is = 1,nspinmx
      do iqq=1,nq
        irkip_all(is,:,:,iqq)=irk
      enddo
      enddo
 

!! -- ppb= <Phi(SLn,r) Phi(SL'n',r) B(S,i,Rr)>
      allocate( ppbir(nlnmx*nlnmx*mdimx*nclass,ngrp,nspin))
      do irot = 1,ngrp
      do isp = 1,nspin
        call ppbafp_v2 (irot,ngrp,isp,nspin,
     i              il,in,im,nlnm, !w(i_mnl),
     d              nl,nn,nclass,nlnmx,
     i              mdimx,lx,nx,nxx, !Bloch wave    
     i              cgr, nl-1,  !rotated CG
     i              ppbrd,      !radial integrals
     o              ppbir(:,irot,isp))
      enddo  
      enddo  

!! == job divider for MPI ==
c      nrank=1 !total number of rank
c      irank=0 !rank for local process
c      irkip=irkip_all
      allocate(irkip(nspinmx,nqibz,ngrp,nq)) !local
      call MPI__hsfp0_rankdivider(irkip,irkip_all,nspinmx,nqibz,ngrp,nq) ! MIZUHO-IR
c      nrkip = nrkip_all         !we don't need to change this for MPI case. It just need to distribute non-zero irkip.
!! ----------------------------------------
      nlnx4    = nlnx**4
      niwx     = max0 (nw+1,niw)
      allocate( ppb(nlnmx*nlnmx*mdimx*nclass),  eq(nband), 
     &     kount(nqibz,nq), zsec(iwini:iwend,ntq,nq),  
     &     coh(ntq,nq) )

c     loop over spin ----------------------------------------------------
      if (tote) exx = 0d0
      do 2000 is = 1,nspinmx
        if(mpi__root) then
        if((exchange).and.(.not.tote)) then
          write(ifsex(is),*) '======================================='
          write(ifsex(is),"('Self-energy exchange SEx(q,t): is=',i3)") is
          write(ifsex(is),*) '======================================='

          call winfo(ifsex(is),nspin,nq,ntq,is,nbloch
     &           ,ngpn1,ngcn1,nqbz,nqibz,ef,deltaw,alat,esmr)

          write (ifsex(is),*)' *** '
          write (ifsex(is),"(a)") ' jband   iq ispin                  
     &           qvec 
     &           eigen-Ef (in eV)
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
        endif
c     initialise secq and kount
        zsec  = 0d0
        coh   = 0d0
        kount = 0
!make dummy SEXcore
        if(ixc==3.and.nctot==0) goto 2001

        call sxcf_fal3z (kount, ixc,deltaw,shtw,q,itq,ntq,ef,ef2,esmr,esmr2,
c     i        tiat,miat,
     i   nspin,is,
     i   qbas,ginv,qibz,qbz,wbz,nstbz, wibz, 
     i   nstar,irkip(is,:,:,:),  
c     i              iclass,nblocha,nlnmv, nlnmc, !w(i_mnlv),w(i_mnlc)
c     i              icore,ncore, imdim,
c     i              ppb,        
     i   freq_r,freqx, wwx,     !expa,ua,
     i   dw, ecore(:,is),
     d   nlmto,nqibz,nqbz,nctot,
     d   nl,nnc,nclass,natom,
     d   mdimx,nbloch,ngrp, nw_i,nw,  niw,niwx,nq, !nlnmx,
     i   nblochpmx,ngpmx,ngcmx,
     i   wgt0,nq0i,q0i, symgg,alat, nband, ifvcfpout, !shtvg,
     i   exchange, tote, screen, cohtest, ifexsp(is),
     i   iwini,iwend,
     i   nbmx(2),ebmx(2),
c     i   pomatr, qrr,nnr,nor,nnmx,nomx,nkpo, 
     i   wklm,lxklm,            !mar2012
     i   invg,il,in,im,nn, lx,nx,nxx, !ppbrd, !sep 2012 cgr,nlnm,
     o   zsec,coh,exx)          ! acuumulation variable

c        write(6,*)'### end of sxcf_fal3z ### '

c 1001   continue         !end of q-loop
c---------------------------------
c     electron gas bare exchange (exact)
         if (legas.and.exchange.and.(.not.tote)) then
          efz=(ntot*3*pi**2/voltot)**(2d0/3d0) ! ef is calculated from ntot.
          pi         = 4.d0*datan(1.d0)
          tpia       = 2.d0*pi/alat
          qfermi= dsqrt(efz)
          alpha = (9*pi/4d0)**(1d0/3d0)
          write (6,*)' --- exact electron gas bare exchange --- '
          write (6,*)' density parameter rs= ', alpha/qfermi
          write (6,*)' kf= ',qfermi
          do      ip = 1,nq
            qreal =  tpia*q(1:3,ip)
            qm    = dsqrt ( sum(qreal**2) )
            xsex  = hartree * egex (qm,efz)
            write (6,*)
            write (6,"(' True qm-ef Sx=',2f14.6,' q/qf=',f14.6)")
     &              rydberg()*(qm**2-efz), xsex, qm/qfermi
            write (6,"(' Num  qm-ef Sx=',2f14.6)") 
     &              eqx(1,ip,is),        hartree*dreal(zsec(iwini,1,ip))
            write (6,"(' === diff     =',2f14.6)") 
     &              rydberg()*(qm**2-efz)-eqx(1,ip,is)
     &              , xsex - hartree*dreal(zsec(iwini,1,ip))

            write (661,"(' qm True qm-ef Sx=',3f14.6)") 
     &              qm,rydberg()*(qm**2-efz), xsex
            write (662,"(' qm Num  qm-ef Sx=',3f14.6)") 
     &              qm,eqx(1,ip,is),     hartree*dreal(zsec(iwini,1,ip))
ccc   write (ifsex(is),6600) qreal(1),qreal(2),qreal(3),xsex
ccc   write (6,6600) qreal(1),qreal(2),qreal(3),xsex
ccc   6600   format (' qreal =',3f8.4,'   SEx(q) =',d13.5)
            write (663,"(2f14.6)") qm/qfermi, qfermi
          end do
         endif
 2001   continue
!!
        call MPI__AllreduceSum( zsec(:,:,:),(iwend-iwini+1)*ntq*nq ) ! MIZUHO-IR

c     write SEx(q,t)
c     call cv      (hartree,sexq,ntq*nq,sexq) !Because vcoul is in a.u.


!! write mpi case
        if(mpi__root) then
        if(ixc==5.or.ixc==6) then
!     nothing
        elseif(exchange) then
          ifoutsex=ifsex(is)
          write(6,*)
          do ip = 1,nq
            do i  = 1,ntq
              write(ifoutsex,"(3i5,3d24.16,3x,d24.16,3x,d24.16)")
     &                 itq(i),ip,is, q(1:3,ip), eqx(i,ip,is), 
     &                 hartree*dreal(zsec(iwini,i,ip))
              write(6,"(' j iq isp=' i3,i4,i2,'  q=',3f8.4,
     &                 '  eig=',f10.4,'  Sx=',f10.4)")
c     &     '  eig=',f10.4,'  Sx=',d15.7)") 
     &                 itq(i),ip,is, q(1:3,ip), eqx(i,ip,is),
     &                 hartree*dreal(zsec(iwini,i,ip))
            end do
          end do
c-------------------------
        elseif(ixc==2) then
          ifoutsec=ifsec(is)
          do ip = 1,nq
            do i  = 1,ntq
              zseciip = 0d0
              zseciip(iwini:iwend) = zsec(iwini:iwend, i, ip)
              dscdw1  = dreal( zseciip(1) - zseciip(0)  ) /deltaw
              zfac1   = 1d0/(1d0-dscdw1)
              dscdw2  = dreal( zseciip(0) - zseciip(-1) ) /deltaw
              zfac2   = 1d0/(1d0-dscdw2)
              dscdw   = dreal( zseciip(1) - zseciip(-1) ) /(2d0*deltaw)
              zfac    = 1d0/(1d0-dscdw)
              write(6,"(' j iq isp=' i3,i4,i2,'  q=',3f8.4,
     &                 '  eig=',f8.4,'  Re(Sc) =',3f8.4,'  Img(Sc) =',3f8.4,
     &                 '  zfac=',f8.4,'  zfac1&2 =',2f8.4)") 
     &                 itq(i),ip,is, q(1:3,ip), eqx(i,ip,is), 
     &                 hartree*dreal(zseciip(-1:1)), 
     &                 hartree*dimag(zseciip(-1:1)),zfac,zfac1,zfac2
              write(ifoutsec,"(3i5,3d24.16,3x,d24.16,3x,3d24.16,
     &                 3x,3d24.16,3x,3d24.16)")
     &                 itq(i),ip,is, q(1:3,ip), eqx(i,ip,is), 
     &                 hartree*dreal(zseciip(-1:1)), 
     &                 hartree*dimag(zseciip(-1:1)),zfac,zfac1,zfac2
              if(cohtest) then
                write(671,"( i3,i4,i2, 3f8.4,
     &                    '  eig=',f8.4,'  coh =',f10.4)") 
     &                    itq(i),ip,is, q(1:3,ip), eqx(i,ip,is), 
     &                    hartree*coh(i,ip)
              endif
            end do
          end do
c-------------------------
        elseif(ixc==4) then    !spectrum mode
          ifoutsec=9300
          if(is==1) sss='.UP'
          if(is==2) sss='.DN'
          open(ifoutsec,file='SEComg'//sss)
          do ip = 1,nq
            do i  = 1,ntq
              write(6,"(' --- j iq isp=' i3,i4,i2,' q=',3f8.4,' eig=',f8.4)") 
     &                 itq(i),ip,is, q(1:3,ip), eqx(i,ip,is)
              do iw = iwini,iwend
                if(mod(iw,50)==1) 
     &                    write(6,"(' omega-ef=',f9.4,'  Sc=',2f9.4)") 
     &                    (omega(i,iw)-ef)*rydberg(), hartree*zsec(iw,i,ip)
                write(ifoutsec,"(4i5,3f10.6,3x,f10.6,2x,f16.8,x,3f16.8)")
     &                    iw,itq(i),ip,is, q(1:3,ip),  eqx(i,ip,is),
     &                    (omega(i,iw)-ef)*rydberg(),  hartree*zsec(iw,i,ip) !,sumimg
              end do
              write(6,*)
            end do
          end do
          close(ifoutsec)
        endif
        endif
 2000 continue                  !end of spin-loop
!!
      if(mpi__root) then
      if (tote) then
        write(ifexx,*) '======================================='
        if(ixc==5) write(ifexx,*) 
     &        '  Exchange energy valence-valence    Exx (eV)   '
        if(ixc==6) write(ifexx,*) 
     &        '  Exchange energy core-valence       Exx (eV)   '
c     write(ifexx,*) '    Exchange energy     Exx (eV)       '
        write(ifexx,*) '======================================='

        write (ifexx,*)' *** '
        if (nspinmx == 1) exx = exx * 2.d0
        write(ifexx,*)exx*hartree

        call winfo2(ifexx,nspin,nq,ntq,is,nbloch
     &        ,ngpn1,ngcn1,nqbz,nqibz,ef,ef2,deltaw,alat,esmr,esmr2)

        if (legas) then
          pi         = 4.d0*datan(1.d0)
          efz=(ntot*3*pi**2/voltot)**(2d0/3d0) ! ef is calculated from ntot.
          qfermi= dsqrt(efz)
          alpha = (9*pi/4d0)**(1d0/3d0)
          rs    = alpha/qfermi
          write (ifexx,*)' --- electron gas ---'
          write (ifexx,*)' density parameter rs= ', rs
          write (ifexx,*)' kf= ',qfermi
          write (ifexx,*)' *** exact exchange'
          exxelgas = -0.4582 * hartree /rs * ntot
          write (ifexx,*)exxelgas
        endif
      endif  
      endif  

C     --- EXspectrum ----------------------------------------------------
c     This section is similar with efsimplef.f
      if(mpi__root) then
      if(sum(ifexsp(1:nspin))/=0) then
        do is = 1,nspin
          write(6,*)' --- Goto ExSpectrum section --- is=',is
          rewind (ifexsp(is))
          itmx = 0
          do
            read(ifexsp(is),*,end=1215)ipex,itpex,itex,qqex(1:3), eex,exsp
            if(itex>itmx) itmx=itex
          enddo
 1215     continue
          nspexmx = itmx*(nqbz+nq0i*ngrp) !Get marimum value of the number of the ex spectrum
c
          allocate( eex1(nspexmx,ntq,nq), exsp1(nspexmx,ntq,nq),
     &           nspex(ntq,nq) ,
     &           itex1(nspexmx,ntq,nq),
     &           qqex1(3,nspexmx,ntq,nq) )
          write(6,*)' nspexmx =',nspexmx
c
          rewind (ifexsp(is))
          nspex = 0
          do
            read(ifexsp(is),*,end=1216) ipex,itpex,itex,qqex(1:3),eex,exsp
            nspex(itpex,ipex) = nspex(itpex,ipex)+1
            iex  = nspex(itpex,ipex)
            eex1  (iex,itpex,ipex) = eex
            exsp1 (iex,itpex,ipex) = exsp
            itex1 (iex,itpex,ipex) = itex
            qqex1(:,iex,itpex,ipex)= qqex
          enddo
 1216     continue            !Get eex1(1:nspex) exsp1(1:nspex) for itp ip.
          write(6,*)' nspex(1 1)=',nspex(1,1)
c
          do ipex = 1,nq
            do itpex=1,ntq
              write(6,*)' is itq ip =',is,itq,ip
              nnex = nspex(itpex,ipex)
              allocate( ieord(1:nnex) )
              call sortea( eex1(1:nnex,itpex,ipex),ieord, nnex,isig)
              eex1 (1:nnex,itpex,ipex)  = eex1  (ieord(1:nnex),itpex,ipex)
              exsp1 (1:nnex,itpex,ipex) = exsp1 (ieord(1:nnex),itpex,ipex)
              itex1 (1:nnex,itpex,ipex) = itex1 (ieord(1:nnex),itpex,ipex)
              qqex1(:,1:nnex,itpex,ipex)= qqex1 (:,ieord(1:nnex),itpex,ipex)

              filenameex = 'EXSP'//charnum3(ipex)//charnum3(itpex)//'.'//char(48+is)
              ifexspx=4111
              open(ifexspx,file=filenameex)
              filenameex = 'EXSS'//charnum3(ipex)//charnum3(itpex)//'.'//char(48+is)
              ifexspxx=4112
              open(ifexspxx,file=filenameex)
              do i=1,nnex
                write(ifexspx, "(2d14.6, i4, 3f14.6)")
     &                    eex1  (i,itpex,ipex), exsp1 (i,itpex,ipex),
     &                    itex1 (i,itpex,ipex), qqex1 (1:3,i,itpex,ipex)
              enddo

              eee  =-1d99
              exwgt= 0d0
              do i=1,nnex
                if(eex1(i,itpex,ipex) > eee+1d-4 .or. i==nnex) then
                  if(i/=1) write(ifexspxx, "(2d23.15)")
     &                       eee, exwgt*hartree
                  eee  = eex1(i,itpex,ipex)
                  exwgt= exsp1 (i,itpex,ipex)
                else
                  exwgt= exwgt + exsp1 (i,itpex,ipex)
                endif
              enddo
c
              deallocate( ieord )
              close(ifexspx)
              close(ifexspxx)
            enddo
          enddo
          deallocate( eex1, exsp1, nspex, itex1, qqex1 )
        enddo
        write(6,*)' End of ExSpectrum section ---'
      endif
      endif
      call cputid(0)
      write(6,*) '--- end of hsfp0_sc --- irank=',mpi__rank
      call flush(6)
      call MPI__Finalize 
      if(ixc==1 ) call rx0( ' OK! hsfp0: Exchange mode')
      if(ixc==2 ) call rx0( ' OK! hsfp0: Correlation mode')
      if(ixc==3 ) call rx0( ' OK! hsfp0: Core-exchange mode')
      if(ixc==4 ) call rx0( ' OK! hsfp0: spectrum mode')
      if(ixc==5 ) call rx0( ' OK! hsfp0: Exx mode val-val  See TEEXXvv')
      if(ixc==6 ) call rx0( ' OK! hsfp0: Exx mode core-val See TEEXXcv')
      end program hsfp0
