subroutine hsfp0()
  use m_ReadEfermi,only: Readefermi,ef_read=>ef
  use m_readqg,only: Readqg0,Readngmx2,ngpmx,ngcmx
  use m_hamindex,only:   Readhamindex
  use m_readeigen,only: Init_readeigen,Readeval,Lowesteval,Init_readeigen2
  use m_read_bzdata,only: Read_bzdata, &
       nqbz,nqibz,nqbzw,nteti,ntetf,n1,n2,n3,ginv,qlat, &
       dq_,qbz,wbz,qibz,wibz,qbzw, idtetf,ib1bz,idteti, &
       nstar,irk,nstbz, lxklm,dmlx,epinvq0i,wklm, wqt=>wt,q0i,nq0i
  use m_hamindex,only: ngrp, symgg=>symops
  use m_genallcf_v3,only: Genallcf_v3, &
       nclass,natom,nspin,nl,nn, nlmto,nlnmx, nctot,niw, &
       alat,delta,deltaw,esmr_in=>esmr,clabl,iclass, il,in,im,nlnm, &
       plat, pos,z,ecore,  konf,nlnx
  use m_keyvalue,only: Getkeyvalue
  use m_rdpp,only: Rdpp, nblocha,lx,nx,ppbrd,mdimx,nbloch,cgr,nxx
  use m_zmel,only: Mptauof_zmel
  use m_itq,only: Setitq_hsfp0, itq,ntq
  use m_mpi,only: &
       MPI__Initialize,MPI__real8send,MPI__real8recv,MPI__send_iv,MPI__recv_iv,MPI__sxcf_rankdivider, &
       MPI__Finalize,MPI__root,MPI__Broadcast,MPI__rank,MPI__size,MPI__allreducesum, &
       MPI__consoleout, &
       MPI__barrier
  use m_readhbe,only: Readhbe, nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg
  use m_lgunit,only: m_lgunit_init
  use m_freq,only: freq01
  use m_anf,only:  Anfcond, laf
  use m_bzints,only: bzints2x
  implicit none
  !! = Calculate the diagonal part of self-energy \Sigma within the GW approximation. And some other functions =
  !  See document at the top of hsfp0.sc.m.F

  !     SEx(q,t,t) = <psi(q,t) |SEx| psi(q,t)>
  !     SEc(q,t,t) = <psi(q,t) |SEc| psi(q,t)>

  !     where SEc(r,r';w) = (i/2pi) < [w'=-inf,inf] G(r,r';w+w') Wc(r,r';w') >

  !    For given ixc (read(5,*), hsfp0 yields variety of results.
  !    (but some are already obsolate).
  !     mode= 11: exchange    mode SEx, the exchange part of the self-energy
  !     mode= 12: correlation mode SEc, the correlated part of the self-energy
  !     mode= 3: core exchange mode SEXcore
  !     mode= 4: plot spectrum function ---See manual
  !     mode= 5: exchange    mode Exx, the exchange part of the total energy
  !     mode= 6: core exchange mode Exxc, Exx(core-valence)

  !     May 2002 Takashi Miyake. Total energy calc.
  !     This hsfp0 is build from hsec10.f by F.Aryasetiawan.
  !-----------------------------------------
  !     real(8),parameter ::  ua = 1d0    ! constant in w(0)exp(-ua^2*w'^2) to take care of peak around w'=0
  !------------------------------------
  !  \Sigma = \Sigma_{sx} + \Sigma_{coh} + \Sigma_{img axis} + \Sigma_{pole} by Hedin PR(1965)A785
  !  I found COH term has inevitably poor accuracy.
  logical ::tetra, tetra_hsfp0, screen  = .false., &
       cohtest = .false., &  ! &  \Sigma_{coh}. mode swich is not required.
       tetraex = .false.    ! This switch is only meaningful for mode=1,5,6
  ! If you want to calculate exchange, use tetraex=T .
  ! Note that you have to supply EFERMI by the tetrahedon method.
  !------------------------------------
  real(8):: esmr2,shtw,esmr
  integer::  ixc,  ibas,ibasx,ifqpnt,ifwd, &! & ngpmx,ngcmx,!nbloch,
       nprecx,mrecl,nblochpmx2,nwp=0,niwt, nqnum,nblochpmx, &! & mdimx,
       noccxv,maxocc2,noccx,ifvcfpout,iqall,iaf, &! & ntq, !ifrcw,ifrcwi,
       i,k,nspinmx, nq,is,ip,iq,idxk,ifoutsex,ig, &
       mxkp,nqibzxx,ntet,nene,iqi, ix,iw, &
       nlnx4,niwx,irot,invr,invrot,ivsum, ifoutsec, &
       ifsec(2),  ifxc(2),ifsex(2), ifphiv(2),ifphic(2),ifec,ifexsp(2), &
       ifsecomg(2),ifexx,ndble=8
  real(8) :: pi,tpia,vol,voltot,rs,alpha, &
       qfermi,efx,valn,efnew,edummy,efz,qm,xsex,egex,edummyd(1), &
       zfac1,zfac2,dscdw1,dscdw2,dscdw,zfac,ef2=1d99,exx,exxq,exxelgas
  logical :: lqall!,laf
  real(8),allocatable    :: q(:,:)
  integer,allocatable :: ngvecp(:,:), ngvecc(:,:),iqib(:), kount(:,:)
  real(8),allocatable:: vxcfp(:,:,:), &
       wgt0(:,:), &
       eqt(:), &
       ppbrdx(:,:,:,:,:,:,:),aaa(:,:), &
       ppb(:), eq(:), &
       eqx(:,:,:),eqx0(:,:,:),ekc(:),coh(:,:),eqxx(:,:)
  complex(8),allocatable:: geigB(:,:,:,:) ,zsec(:,:,:)

  logical :: exchange, legas, tote
  real(8) ::  rydberg,hartree
  real(8):: qreal(3), ntot,tripl,xxx(3,3) !nocctotg2,
  logical ::nocore
  real(8),allocatable :: qz(:,:),qbzxx(:),wbzxx(:),wtet(:,:,:,:), &
       eband(:,:,:), ene(:) !,ecore(:,:)
  integer(4),allocatable ::idtetx(:,:),idtet(:,:),ipq(:) &
       ,iene(:,:,:),ibzx(:) ! ,nstar(:)

  integer ::ib,iqx,igp,iii,ivsumxxx,isx,iflegas, iqpntnum

  real(8),allocatable   :: eex1(:,:,:),exsp1(:,:,:),qqex1(:,:,:,:)
  integer,allocatable:: nspex(:,:),ieord(:),itex1(:,:,:)
  real(8)    :: qqex(1:3), eex,exsp,eee, exwgt,deltax0
  integer :: itmx,ipex,itpex,itex,nspexmx,nnex,isig,iex,ifexspx &
       ,ifexspxx ,ifefsm, nq0ix,ifemesh,nz
  character(3)  :: charnum3,sss
  character(12) :: filenameex
  logical :: exspwrite=.false.
  character(8):: xt

  integer :: iwini,iwend
  real(8),allocatable:: omega(:,:)
  real(8) ::  omegamax,dwplot,omegamaxin

  integer::nqbze,ini,nq0it,idummy

  real(8)   :: ebmx(2)
  integer:: nbmx(2)

  real(8):: volwgt

  integer:: incwfin
  real(8):: ddw
  real(8),allocatable::freqx(:),freqw(:),wwx(:),expa(:)

  !      logical:: GaussSmear=.true.      !readgwinput,
  integer::ret
  character*(150):: ddd
  integer:: bzcase=1,  ngpn1,verbose,ngcn1,nwxx !mrecg,
  real(8)   :: wgtq0p,quu(3)
  real(8),allocatable:: freq_r(:)
  logical ::smbasis
  integer:: ifpomat,nkpo,nnmx,nomx,ikpo,nn_,no,nss(2)
  real(8):: q_r(3)
  real(8),allocatable:: qrr(:,:)
  integer,allocatable:: nnr(:),nor(:)
  complex(8),allocatable:: pomatr(:,:,:),pomat(:,:)
  real(8)::sumimg
  logical :: allq0i            !S.F.Jan06
  integer:: nw_i
  complex(8):: zseciip(-1:1)
  logical :: debug=.false.
  real(8),allocatable:: vcousq(:),vcoud(:)
  complex(8),allocatable:: zcousq(:,:)
  !      logical:: newaniso=.true.
  integer:: ifvcoud,ifidmlx,iqq

  integer,allocatable:: irkip_all(:,:,:,:),irkip(:,:,:,:)
  integer:: irank,nrank
  integer:: ififr
  integer:: timevalues(8)  ,isp,dest,ificlass,ifiq0p
  character(128) :: ixcc
  integer:: nw,ifcoh, ixx(2)
  real(8)::dwdummy, ef
  logical:: cmdopt2
  character(20):: outs=''
  !...  mode switch. --------------
  call MPI__Initialize()
  call M_lgunit_init()
  call date_and_time(values=timevalues)
  write(6,'(a,9i5)')'dateandtime1=',MPI__rank,timevalues(1:8)
  hartree=2d0*rydberg()
  if(cohtest) then
     screen = .true.
     ixc = 2; nz=0
     open(newunit=ifcoh,file='COH')
  elseif(MPI__root) then
     if(cmdopt2('--job=',outs)) then
        read(outs,*) ixc
     else
        write(6,*) ' --- Choose omodes below ----------------'
        write(6,*) '  Sx(1) Sc(2) ScoreX(3) Spectrum(4) '
        write(6,*) '  EXX_val-val (5)  Exx_core-val(6) '
        write(6,*) '  Sx_sf(11) Sc_sf(12) '
        write(6,*) '  [option --- (+ QPNT.{number} ?)] '
        write(6,*) ' --- Put number above ! -----------------'
        read(5,*) ixc           !call readin5(ixc,nz,idummy)
        write(*,*) ' ixc=', ixc !computational mode index
     endif
     nz=0
  endif
  call MPI__Broadcast(ixc)
  call MPI__Broadcast(nz)
  write(ixcc,"('.mode=',i4.4)")ixc
  call MPI__consoleout('hsfp0'//trim(ixcc))
  write(6,*) ' ixc nz=',ixc, nz
  if(ixc==0) call rx( ' --- ixc=0 --- Choose computational mode!')
  !!  tetraex is only for Ex.
  tetraex=tetra_hsfp0()
  iii=verbose()
  write(6,*)' verbose=',iii
  if(ixc==1 .OR. ixc==5 .OR. ixc==6) then;
  else; tetraex=.false.
  endif
  !! These are mainly used now.
  if(ixc==11 .OR. ixc==12) then
     ixc=ixc-10
  endif

  !! ===  readin BZDATA. See gwsrc/rwbzdata.f ===
  !! See use m_read_bzdata,only: at the top of this routine
  call read_BZDATA()
  write(6,*)' nqbz nqibz =',nqbz,nqibz
  call pshpr(60)

  !! === readin GWIN and LMTO, then allocate and set datas. ===
  !! See use m_genallcf_v3,only: at the top of this routine
  !      nwin =-999                !not readin NW file
  !      efin =-999d0              !not readin EFERMI
  if    (ixc==3) then;  incwfin= -2 !core exchange mode
  elseif(ixc==5) then;  incwfin= -4 ! valence-valence Ex energy mode. ! See rgwinf called from genallcf_*
  elseif(ixc==6) then;  incwfin= -3 ! core-valence Ex energy mode.
  else               ;  incwfin= -1 !use 7th colmn for core at the end section of GWIN
  endif
  call genallcf_v3(incwfin) ! module m_genallcf_v3. See use m_genallcf in this rouitine

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
  call pshpr(30)
  pi   = 4d0*datan(1d0)
  tpia = 2d0*pi/alat
  !      call dinv33(plat,1,xxx,vol)
  !      voltot = dabs(vol)*(alat**3)
  voltot = abs(alat**3*tripl(plat,plat(1,2),plat(1,3)))
  shtw = 0d0
  !! if(esmr<1d-5) shtw=0.01d0 ! Ferdi's shift to avoid resonance effect(maybe), I used this until sep2012
  esmr=esmr_in
  !$$$!! ef is taken as rs for the empty-sphere test case of legas=T case
  !$$$!! HOMOGENIOUS GAS code. Usually not used.Need fixing if necessary.
  !$$$      INQUIRE (FILE = 'LEGAS', EXIST = legas)
  !$$$      if(legas) then            !!! test for electron gas case.
  !$$$        write(6,*)' find LEGAS. legas =',legas
  !$$$        iflegas = 2101
  !$$$        open (iflegas,file='LEGAS')
  !$$$        read(iflegas,*)rs
  !$$$        close(iflegas)
  !$$$        alpha = (9*pi/4d0)**(1d0/3d0)
  !$$$        qfermi = alpha/rs
  !$$$        efx  = qfermi**2
  !$$$        valn = efx**1.5d0*voltot/3d0/pi**2
  !$$$        write (6,*)'  #### egas test mode  legas=T #### given rs =',rs
  !$$$        write (6,*)' egas  Exact Fermi momentum  qf  =', qfermi
  !$$$        write (6,*)' egas  Exact Fermi energy    Ef  =', efx
  !$$$        if(tetraex) call rx( 'legas You have to give ef of  tetrahedron')
  !$$$      endif

  !!
  ifexsp=0
  if(ixc==1) then
     exchange=.true.
     tote=.false.
     write(6,*) ' --- Exchange mode --- '
     if(MPI__root) then
        open(newunit=ifxc(1),file='XCU'//xt(nz))
        open(newunit=ifsex(1),file='SEXU'//xt(nz))
        if (nspin == 2) then
           open(newunit=ifxc(2),file='XCD'//xt(nz))
           open(newunit=ifsex(2),file='SEXD'//xt(nz))
        endif
        INQUIRE (FILE = 'EXspTEST', EXIST = exspwrite)
        if(exspwrite) then
           write(6,*)'--- Find EXspTEST ExspectrumWrite=',exspwrite
           write(6,*)'--- esmr is chosen to be 2d0 Ry'
           esmr= 2d0
           do is=1,nspin
              open(newunit=ifexsp(is),file='EXSP.'//char(48+is))
           enddo
        endif
     endif
  elseif(ixc==2) then
     exchange=.false.
     tote=.false.
     write(6,*) ' --- Correlation mode --- '
     if(cohtest) write(6,*) ' COH calculation mode. Results in COH'
     if(MPI__root) then
        open(newunit=ifsec(1),file='SECU'//xt(nz))
        if (nspin == 2) &
             open(newunit=ifsec(2),file='SECD'//xt(nz))
     endif
  elseif(ixc==3) then
     exchange=.true.
     tote=.false.
     esmr=0d0
     write(6,*) ' --- CORE Exchange mode --- '
     if(MPI__root) then
        open(newunit=ifsex(1),file='SEXcoreU'//xt(nz))
        if(nspin == 2) open(newunit=ifsex(2),file='SEXcoreD'//xt(nz))
     endif
     !! ixc=4,5,6 not checked now... NEED MPI and so on...
  elseif(ixc==4) then
     write(6,*) ' --- Spectrum function Sigma(omega) mode --- '
     exchange=.false.
     tote=.false.
     !        open(newunit=ifsecomg(1),file='SEComgU'//xt(nz))
     !        if (nspin == 2) open(newunit=ifsecomg(2),file='SEComgD'//xt(nz))
     !        call rx('need to test this mode again')
  elseif(ixc==5) then
     exchange=.true.
     tote=.true.
     write(6,*) ' --- Exx mode valence-valence --- '
     open(newunit=ifexx,file='TEEXXvv')
     esmr2 = esmr           !takao June2002
     call rx('need to test this mode again')
  elseif(ixc==6) then
     exchange=.true.
     tote=.true.
     esmr2 = esmr           !takao June2002
     esmr  = 0d0            !takao June2002
     write(6,*) ' --- CORE Exx mode core-valence --- '
     open(newunit=ifexx,file='TEEXXcv')
     call rx('need to test this mode again')
     write(6,'("    esmr2    =",f13.6)') esmr2
  else
     call rx( ' hsfp0: Need input (std input) 1-6!')
  endif
  write(6, *) ' --- computational conditions --- '
  write(6,'("    deltaw  =",f13.6)') deltaw
  write(6,'("    esmr    =",f13.6)') esmr
  write(6,'("    alat voltot =",2f13.6)') alat, voltot
  call Readhbe()    !Read dimensions of h,hb
  call Readhamindex()
  call init_readeigen()!nband,mrece) !initialization of readEigen
  call Mptauof_zmel(symgg,ngrp) !Get space-group transformation information. in m_zmel
  call Readngmx2() !return ngpmx and ngcmx in m_readqg
  call readqg0('QGpsi',qibz(1:3,1), quu,ngpn1)
  call readqg0('QGcou',qibz(1:3,1), quu,ngcn1)
  write(6,*) ' max number of G for QGpsi and QGcou: ngcmx ngpmx=',ngcmx,ngpmx
  if(debug) write(6,*) ' end of read QGcou'
  !! === readin plane wave parts, and Radial integrals ppbrd. ===
  !     ppbrd = radial integrals
  !     cgr   = rotated cg coeffecients.
  !     geigB = eigenfunction's coefficiens for planewave.
  !     ngvecpB (in 1stBZ) contains G vector for eigen function.
  !     ngveccB (in IBZ)   contains G vector for Coulomb matrix.
  !      call Rdpp(ngrp,symgg) =>moved to m_zmel
  call pshpr(60)
  !! Readin WV.d
  if( .NOT. exchange .OR. (exchange .AND. screen)) then !screen means screened exchange case
     open(newunit=ifwd,file='WV.d') !direct access files WVR and WVI which include W-V.
     read(ifwd,*) nprecx,mrecl,nblochpmx,nwp,niwt, nqnum, nw_i
     write(6,"(' Readin WV.d =', 10i8)") nprecx,mrecl,nblochpmx, nwp, niwt, nqnum, nw_i
     close(ifwd)
     if(nprecx/=ndble)call rx("hsfp0: dim of WVR and WVI not compatible")!call checkeq(nprecx,ndble)
     nw=nwp-1
     if (niwt /= niw) call rx( 'hsfp0: wrong niw')
     !! Energy mesh; along real axis. Read 'freq_r'
     !! NOTE nw_i=nw for non-timereversal case.
     !!      nw_i=0 for time-reversal case.
     !!  NOTE: We assume freq_r(i) == -freq_r(-i) in this code. feb2006
     !!  NOTE: this program assumes freq_r(iw)=freq_r(-iw). freq_r(iw <0) is redundant.
     open(newUNIT=ififr,file='freq_r')
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
        if( sum(abs( freq_r(1:nw)+freq_r(-1:-nw:-1)))/=0) &
             call rx( "sxcf_fal3_scz: freq_r /= -freq_r")
     endif
  endif
  !! efermi by tetrahedron. this can be overwritten
  call readefermi()
  ef=ef_read
  if(tetraex) goto 201

  !!== Determine Fermi energy ef for given valn (legas case), or corresponding charge given by z and konf.==
  !!    When esmr is negative, esmr is geven automatically by efsimplef.
  legas=.false.
  if((ixc==3) .OR. (ixc==6)) then
     ef = lowesteval() -1d-3 !lowesteigen(nspin,nband,qbz,nqbz) - 1d-3 !lowesteb was
     write(*,*)'efffffffffff=',ef
     if(maxval(ecore(:,1:nspin))>ef) then
        call rx( 'hsfp0 ixc=3/6:  ecore>evalence ')
     endif
     if(ixc==6) then
        call efsimplef2ax(legas,esmr2, valn,efnew)!Get num of val electron valn and Fermi energy ef. legas=T give ef for given valn.
!        call efsimplef2a(nspin,wibz,qibz,ginv, &
!             nband,nqibz &
!             ,konf,z,nl,natom,iclass,nclass &
!             ,valn, legas, esmr2, qbz,nqbz,efnew)
        ef2 = efnew
        write(6,*)' end of efsimple ef2 esmr2=',ef2,esmr2
     endif
  else                      ! if(esmr/=0d0) then
     !     --- determine Fermi energy ef for given valn (legas case) or corresponding charge given by z and konf.
     !     When esmr is negative, esmr is geven automatically by efsimplef.
        call efsimplef2ax(legas,esmr, valn,efnew)!Get num of val electron valn and Fermi energy ef. legas=T give ef for given valn.
!     call efsimplef2a(nspin,wibz,qibz,ginv, nband,nqibz &
!          ,konf,z,nl,natom,iclass,nclass,valn, legas, esmr, qbz,nqbz,efnew)
     ef = efnew
     if (ixc==5) ef2 = ef
     !     - check total ele number -------
     !        ntot  = nocctotg2(nspin, ef,esmr, qbz,wbz, nband,nqbz) !wbz
     write(6,*)' ef    =',ef
     write(6,*)' esmr  =',esmr
     write(6,*)' valn  =',valn
     !        write(6,*)' ntot  =',ntot
  endif
  !!
201 continue
  if(tote) then
     ddw= .5d0
     ddw= 10d0
     block
       real(8):: ekt(nband,nqbz,nspin),eff
       if(ixc==6) then
          eff= ef2+ ddw*esmr2
       else
          eff= ef + ddw*esmr
       endif
       do is = 1,nspin
          do iq = 1,nqbz
             ekt(:,iq,is)= readeval(qbz(:,iq),is)
          enddo
       enddo
       noccxv = maxval(count(ekt(1:nband,1:nqbz,1:nspin)<eff,1))
     endblock
  endif
  !     noccxv = maxocc (ifev,nspin, ef+0.5d0*esmr, nband,nqbz)  ! maximum no. of occupied valence states
  !     maxocc seems to give (the maxmum number of occ + 1).

  call init_readeigen2()!mrecb,nlmto,mrecg) !initialize m_readeigen
  call Getkeyvalue("GWinput","QPNT_nbandrange",nss,2,default=(/-99997,-99997/) )

  if( .NOT. tote) then
     if(nz==0) then
        call getkeyvalue("GWinput","<QPNT>",unit=ifqpnt,status=ret)
     else
        open(newunit=ifqpnt,file='QPNT'//xt(nz))
     endif
  endif
  !   write(6,*)' ifqpnt ret=',ifqpnt,ret
  call anfcond()
!  lqall=.true.
  if (tote) then
     lqall      = .true.
  else
     lqall      = .false.
     call readx   (ifqpnt,10)
     read (ifqpnt,*) iqall!,iaf
     if (iqall == 1) lqall = .TRUE. 
     call readx   (ifqpnt,100)
  endif
  call setitq_hsfp0(ngcmx,ngpmx,tote,ifqpnt,noccxv,nss)
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

  nspinmx = nspin
  if (laf) nspinmx =1

  ! necessary ?
  call MPI__Broadcast(nq)
  if(MPI__root) then
     do dest=1,MPI__size-1
        call MPI__REAL8send(q,3*nq,dest)
     enddo
  else
     call MPI__REAL8recv(q,3*nq,0)
  endif

  !!  read q-points and states
  if(ixc==4) then
     call readxx(ifqpnt)    !skip to ***** for q point for spectrum function.
     omegamaxin = 1d70
     read (ifqpnt,*,err=2038,end=2038) dwplot,omegamaxin
2038 continue
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
  close(ifqpnt)
  !! omega
  !! this is also in sxcf.
  if(ixc==4) then
     allocate(omega(ntq,iwini:iwend))
     do iw = iwini,iwend
        omega(1:ntq,iw) =  dwplot* iw + ef
     enddo
  endif

  !! read LDA eigenvalues
  allocate(eqx(ntq,nq,nspin),eqx0(ntq,nq,nspin),eqt(nband))
  do      is = 1,nspin
     do      ip = 1,nq
        eqt = readeval(q(1,ip),is)
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
     !     loop over spins
     if(MPI__root) then
        do is = 1,nspinmx
           write (ifxc(is),*) '==================================='
           write (ifxc(is),"(' LDA exchange-correlation : is=',i3)")is
           write (ifxc(is),*) '==================================='
           call winfo(ifxc(is),nspin,nq,ntq,is,nbloch &
                ,ngpn1,ngcn1,nqbz,nqibz,ef,deltaw,alat,esmr)
           write (ifxc(is),*)' ***'
           write (ifxc(is),"(a)") ' jband   iq ispin &
                qvec &
                eigen-Ef (in eV) &
                LDA XC (in eV)'
           ifoutsex = ifxc(is)
           write(6,*)
           do ip = 1,nq
              do i  = 1,ntq
                 write(ifoutsex,"(3i5,3d24.16,3x,d24.16,3x,d24.16)") &
                      itq(i),ip,is, q(1:3,ip), eqx(i,ip,is), &
                      vxcfp(i,ip,is)
                 write(6,"(' j iq isp=' i3,i4,i2,'  q=',3f8.4, &
                      '  eig=',f10.4,'  Sxc(LDA)=',f10.4)") &
                      itq(i),ip,is, q(1:3,ip), eqx(i,ip,is), &
                      vxcfp(i,ip,is)
              end do
           end do
           !     end of spin-loop
           close(ifxc(is)) !isx = iclose('XCU'//xt(nz))
           !          if(is==2) isx = iclose('XCD'//xt(nz))
        enddo
     endif !MPI__root
     deallocate(vxcfp)
  endif
  !if( .NOT. exchange) call checkeq(nqibz+nq0i-1, nqnum)
  if(.not.exchange) then
     if(nqibz+nq0i-1/=nqnum) call rx("hsfp0: nqibz+nq0i-1/=nqnum")
  endif
  write(6,*) ' *** nqibz nq0i_total=', nqibz,nq0i
  write(6,*) ' Used k number in Q0P =', nq0i
  write(6,"(i3,f14.6,2x, 3f14.6)" )(i, wqt(i),q0i(1:3,i),i=1,nq0i)
  allocate( wgt0(nq0i,ngrp) )
  call getkeyvalue("GWinput","allq0i",allq0i,default=.false.) 
  call q0iwgt3(allq0i,symgg,ngrp,wqt,q0i,nq0i, wgt0)   
  !--------------------------
  if(nq0i/=0) write(6,*) ' *** tot num of q near 0   =', 1/wgt0(1,1)
  write(6,"('  sum(wgt0) from Q0P=',d14.6)")sum(wgt0)
  if(niw/=0) then !! Generate gaussian frequencies x between (0,1) and w=(1-x)/x
     allocate(freqx(niw),freqw(niw),wwx(niw),expa(niw))
     call freq01(niw, 1d0, freqx,freqw,wwx,expa) !ua=1d0 is dummy 
  endif
  iii= count(irk/=0) !iii=ivsumxxx(irk,nqibz*ngrp)
  write(6,*) " sum of nonzero iirk=",iii, nqbz
  !! === readin Vcoud and EPSwklm for newaniso()=T ===
  !        open(newunit=ifidmlx,file='EPSwklm',form='unformatted')
  !        read(ifidmlx) nq0ix,lxklm
  !        if(nq0i/=nq0ix) then
  !          write(6,*)'nq0i from EPSwklm /= nq0i',nq0i,nq0ix
  !          call rx( 'nq0i from EPSwklm /= nq0i')
  !        endif
  !        allocate( dmlx(nq0i,9))
  !        allocate( epinvq0i(nq0i,nq0i) )
  !        allocate( wklm((lxklm+1)**2))
  !        read(ifidmlx) dmlx, epinvq0i
  !        read(ifidmlx) wklm
  !        close(ifidmlx) !jan2013 bugfix ifidmlx)

  !------------------------------------------------------------------
  !! tetraex section. not checked well recentrly(sep2012) -------
  !! this can be recovered in future.
  if(tetraex) then
     write(6,*) ' !!!!  tetraex=T is not tested well !!!!!!!!!'
     open(newunit=ifefsm,file='EFERMI')
     read(ifefsm,*) ef
     write(6,"(' ef is readin from EFERMI',d23.15)")ef
     allocate(wtet(nband,nspin,nqibz,0:0), &
          eband(nband,nspin,nqibz), qz(3,nqibz)) ! ,nstar(nqibz))
     qz=qibz !call dcopy (3*nqibz,qibz,1,qz,1)
     do  is  = 1,nspin      !Readin eband
        do  iqi = 1,nqibz
           eband(:,is,iqi) = readeval(qz(1:3,iqi),is)
        enddo
     enddo
     volwgt = (3d0 - nspin) / ntetf ! ntetf was =6*n1*n2*n3
     call bzints2x(volwgt,eband,wtet(:,:,:,0),nqibz,nband,nband, &
          nspin,edummy,edummy,edummyd,1,ef,2,nteti,idteti)
     ntot= sum(wtet)
     !$$$        if(legas) then
     !$$$          write(6,"(' tetra=T ef ntot nexact ratio=',15f12.6)") ef,ntot
     !$$$     &           , ef**1.5d0/3d0/pi**2*voltot, ef**1.5d0 /3d0/pi**2*voltot/ntot
     !$$$        else
     !$$$          write(6,"(' tetra=T ef nvalence)=',15f12.6)") ef,ntot
     !$$$        endif
     write(6,"(' tetra=T ef nvalence)=',15f12.6)") ef,ntot
     if(nspin==1) wtet = wtet/2d0
     do iqi = 1,nqibz
        wtet(:,:,iqi,:) = wtet(:,:,iqi,:)/nstar(iqi)
     enddo
     deallocate( eband,qz)  !, ene ) ! pointer for
     !!     -- ibzx denote the index of k{FBZ for given k{1BZ.
     allocate(ibzx(nqbz))
     do iqx  = 1,nqbz
        ixx = findloc(irk(:,:)-iqx,value=0)
        ibzx(iqx)= ixx(1)
     enddo
     if (tote) then
        do i=1,nband
           if( sum(abs(wtet(i,:,:,0))) == 0d0 ) exit
        enddo
        if(  sum(abs(wtet(i+1:nband,:,:,0)))/=0d0) &
             call rx( ' hsfp0: wtetef sum err1')
        noccxv = i-1
     endif
  else
     allocate(wtet(1,1,1,1),ibzx(1)) !dummy
  endif
  !!  end of tetra section --------------------------------------------
  iii=count(irk/=0) !ivsumxxx(irk,nqibz*ngrp)
  write(6,*) " sum of nonzero iirk=",iii, nqbz
  !-----------------------------------------------------------
  !     calculate the the self-energy SEx(ip) or SEc(ip)
  !-----------------------------------------------------------
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
  !! == job divider for MPI ==
  !      nrank=1 !total number of rank
  !      irank=0 !rank for local process
  !      irkip=irkip_all
  allocate(irkip(nspinmx,nqibz,ngrp,nq)) !local
  call MPI__sxcf_rankdivider(irkip_all,nspinmx,nqibz,ngrp,nq, irkip)
  !nrkip = nrkip_all !we don't need to change this for MPI case. It just need to distribute non-zero irkip.
  !! ----------------------------------------
  niwx     = max0 (nw+1,niw)
  allocate( ppb(nlnmx*nlnmx*mdimx*nclass),  eq(nband), &
       kount(nqibz,nq), zsec(iwini:iwend,ntq,nq), coh(ntq,nq) )
  if (tote) exx = 0d0
  isploop: do 2000 is = 1,nspinmx
     if(MPI__root) call hswrite1() !internal routine. Write initial part.
     zsec  = 0d0
     coh   = 0d0
     kount = 0
     if(ixc==3 .AND. nctot==0) goto 2001 !!make dummy SEXcore
     call sxcf_fal3z (kount, ixc,deltaw,shtw,q,itq,ntq,ef,ef2,esmr,esmr2, &
          nspin,is, &
          qlat,ginv,qibz,qbz,wbz,nstbz, wibz, &
          nstar,irkip(is,:,:,:), &
          freq_r,freqx, wwx, &
          dwdummy, ecore(:,is), &
          nlmto,nqibz,nqbz,nctot, &
          nbloch,ngrp, nw_i,nw,  niw,niwx,nq, &
          nblochpmx,ngpmx,ngcmx, &
          wgt0,nq0i,q0i, symgg,alat, nband, ifvcfpout, &
          exchange, tote, screen, cohtest, ifexsp(is), &
          iwini,iwend, &
          nbmx(2),ebmx(2), &
          wklm,lxklm,  dwplot, &
          zsec,coh,exx)          ! acuumulation variable
     !---------------------------------
     !$$$!! Electron gas bare exchange (exact)
     !$$$         if (legas.and.exchange.and.(.not.tote)) then
     !$$$          efz=(ntot*3*pi**2/voltot)**(2d0/3d0) ! ef is calculated from ntot.
     !$$$          pi         = 4.d0*datan(1.d0)
     !$$$          tpia       = 2.d0*pi/alat
     !$$$          qfermi= dsqrt(efz)
     !$$$          alpha = (9*pi/4d0)**(1d0/3d0)
     !$$$          write (6,*)' --- exact electron gas bare exchange --- '
     !$$$          write (6,*)' density parameter rs= ', alpha/qfermi
     !$$$          write (6,*)' kf= ',qfermi
     !$$$          do      ip = 1,nq
     !$$$            qreal =  tpia*q(1:3,ip)
     !$$$            qm    = dsqrt ( sum(qreal**2) )
     !$$$            xsex  = hartree * egex (qm,efz)
     !$$$            write (6,*)
     !$$$            write (6,"(' True qm-ef Sx=',2f14.6,' q/qf=',f14.6)")
     !$$$     &              rydberg()*(qm**2-efz), xsex, qm/qfermi
     !$$$            write (6,"(' Num  qm-ef Sx=',2f14.6)")
     !$$$     &              eqx(1,ip,is),        hartree*dreal(zsec(iwini,1,ip))
     !$$$            write (6,"(' === diff     =',2f14.6)")
     !$$$     &              rydberg()*(qm**2-efz)-eqx(1,ip,is)
     !$$$     &              , xsex - hartree*dreal(zsec(iwini,1,ip))
     !$$$
     !$$$            write (661,"(' qm True qm-ef Sx=',3f14.6)")
     !$$$     &              qm,rydberg()*(qm**2-efz), xsex
     !$$$            write (662,"(' qm Num  qm-ef Sx=',3f14.6)")
     !$$$     &              qm,eqx(1,ip,is),     hartree*dreal(zsec(iwini,1,ip))
     !$$$ccc   write (ifsex(is),6600) qreal(1),qreal(2),qreal(3),xsex
     !$$$ccc   write (6,6600) qreal(1),qreal(2),qreal(3),xsex
     !$$$ccc   6600   format (' qreal =',3f8.4,'   SEx(q) =',d13.5)
     !$$$            write (663,"(2f14.6)") qm/qfermi, qfermi
     !$$$          end do
     !$$$         endif
2001 continue
     call MPI__AllreduceSum( zsec(:,:,:),(iwend-iwini+1)*ntq*nq )
     if(MPI__root) call Hswrite2() ! main output. internal subroutine
2000 enddo isploop
  if(MPI__root .AND. tote) call Hswrite3() ! total energy mode, obsolate now...
  if(MPI__root .AND. sum(ifexsp(1:nspin))/=0) call Hswrite4() !EXspectrum
  write(6,*) '--- end of hsfp0_sc --- irank=',MPI__rank
  call cputid(0)
  call flush(6)
  if(ixc==1 ) call rx0( ' OK! hsfp0: Exchange mode')
  if(ixc==2 ) call rx0( ' OK! hsfp0: Correlation mode')
  if(ixc==3 ) call rx0( ' OK! hsfp0: Core-exchange mode')
  if(ixc==4 ) call rx0( ' OK! hsfp0: spectrum mode')
  if(ixc==5 ) call rx0( ' OK! hsfp0: Exx mode val-val  See TEEXXvv')
  if(ixc==6 ) call rx0( ' OK! hsfp0: Exx mode core-val See TEEXXcv')
  stop
contains !followings are only for writing files.
  ! SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
  subroutine hswrite1()
    if((exchange) .AND. ( .NOT. tote)) then
       write(ifsex(is),*) '======================================='
       write(ifsex(is),"('Self-energy exchange SEx(q,t): is=',i3)") is
       write(ifsex(is),*) '======================================='
       call winfo(ifsex(is),nspin,nq,ntq,is,nbloch,ngpn1,ngcn1,nqbz,nqibz,ef,deltaw,alat,esmr)
       write (ifsex(is),*)' *** '
       write (ifsex(is),"(a)") ' jband   iq ispin &
            qvec &
            eigen-Ef (in eV) &
            exchange (in eV)'
    elseif(ixc==2) then
       write(ifsec(is),*) '=========================================='
       write(ifsec(is),"('Self-energy correlated SEc(qt,w): is=',i3)") is
       write(ifsec(is),*) '=========================================='
       call winfo(ifsec(is),nspin,nq,ntq,is,nbloch,ngpn1,ngcn1,nqbz,nqibz,ef,deltaw,alat,esmr)
       write (ifsec(is),*)' *** '
       write (ifsec(is),"(a)") ' jband   iq ispin &
            qvec &
            eigen-Ef (in eV) &
            Re(Sc) 3-points (in eV) &
            In(Sc) 3-points (in eV) &
            Zfactor'
    endif
  end subroutine hswrite1
  ! SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
  subroutine Hswrite2()
    if(ixc==5 .OR. ixc==6) then
       continue 
    elseif(exchange) then
       ifoutsex=ifsex(is)
       write(6,*)
       do ip = 1,nq
          do i  = 1,ntq
             write(ifoutsex,"(3i5,3d24.16,3x,d24.16,3x,d24.16)") &
                  itq(i),ip,is, q(1:3,ip), eqx(i,ip,is), hartree*dreal(zsec(iwini,i,ip))
             write(6,"(' j iq isp=' i3,i4,i2,'  q=',3f8.4,'  eig=',f10.4,'  Sx=',f10.4)")&
                  itq(i),ip,is, q(1:3,ip), eqx(i,ip,is), hartree*dreal(zsec(iwini,i,ip))
          enddo
       enddo
       !-------------------------
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
             write(6,"(' j iq isp=' i3,i4,i2,'  q=',3f8.4, &
                  '  eig=',f8.4,'  Re(Sc) =',3f8.4,'  Img(Sc) =',3f8.4, &
                  '  zfac=',f8.4,'  zfac1&2 =',2f8.4)") &
                  itq(i),ip,is, q(1:3,ip), eqx(i,ip,is), &
                  hartree*dreal(zseciip(-1:1)), &
                  hartree*dimag(zseciip(-1:1)),zfac,zfac1,zfac2
             write(ifoutsec,"(3i5,3d24.16,3x,d24.16,3x,3d24.16,3x,3d24.16,3x,3d24.16)") &
                  itq(i),ip,is, q(1:3,ip), eqx(i,ip,is), &
                  hartree*dreal(zseciip(-1:1)), &
                  hartree*dimag(zseciip(-1:1)),zfac,zfac1,zfac2
             if(cohtest) then
                write(671,"( i3,i4,i2, 3f8.4,'  eig=',f8.4,'  coh =',f10.4)") &
                     itq(i),ip,is, q(1:3,ip), eqx(i,ip,is),hartree*coh(i,ip)
             endif
          end do
       end do
       !-------------------------
    elseif(ixc==4) then       !spectrum mode
       ifoutsec=9300
       if(is==1) sss='.UP'
       if(is==2) sss='.DN'
       open(ifoutsec,file='SEComg'//sss)
       do ip = 1,nq
          do i  = 1,ntq
             write(6,"(' --- j iq isp=' i3,i4,i2,' q=',3f8.4,' eig=',f8.4)") &
                  itq(i),ip,is, q(1:3,ip), eqx(i,ip,is)
             do iw = iwini,iwend
                if(mod(iw-iwini+1,50)==1) &
                     write(6,"(' omega-ef=',f9.4,'  Sc=',2f9.4)") &
                     (omega(i,iw)-ef)*rydberg(), hartree*zsec(iw,i,ip)
                write(ifoutsec,"(4i5,3f10.6,3x,f10.6,2x,2f16.8,x,3f16.8)") &
                     iw,itq(i),ip,is, q(1:3,ip), wibz(ip), eqx(i,ip,is), &
                     (omega(i,iw)-ef)*rydberg(),  hartree*zsec(iw,i,ip) !,sumimg
             end do
             write(ifoutsec,*)
             write(ifoutsec,*)
          end do
       end do
       close(ifoutsec)
    endif
  end subroutine Hswrite2
  ! SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
  subroutine Hswrite3()
    write(ifexx,*) '======================================='
    if(ixc==5) write(ifexx,*)'  Exchange energy valence-valence    Exx (eV)   '
    if(ixc==6) write(ifexx,*)'  Exchange energy core-valence       Exx (eV)   '
    write(ifexx,*) '======================================='
    write (ifexx,*)' *** '
    if (nspinmx == 1) exx = exx * 2.d0
    write(ifexx,*)exx*hartree
    call winfo2(ifexx,nspin,nq,ntq,is,nbloch,ngpn1,ngcn1,nqbz,nqibz,ef,ef2,deltaw,alat,esmr,esmr2)
    !$$$  if (legas) then
    !$$$  pi         = 4.d0*datan(1.d0)
    !$$$  efz=(ntot*3*pi**2/voltot)**(2d0/3d0) ! ef is calculated from ntot.
    !$$$  qfermi= dsqrt(efz)
    !$$$  alpha = (9*pi/4d0)**(1d0/3d0)
    !$$$  rs    = alpha/qfermi
    !$$$  write (ifexx,*)' --- electron gas ---'
    !$$$  write (ifexx,*)' density parameter rs= ', rs
    !$$$  write (ifexx,*)' kf= ',qfermi
    !$$$  write (ifexx,*)' *** exact exchange'
    !$$$  exxelgas = -0.4582 * hartree /rs * ntot
    !$$$  write (ifexx,*)exxelgas
    !$$$  endif
  end subroutine Hswrite3
  ! SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
  subroutine Hswrite4()
    character(3)  :: charnum3
    do is = 1,nspin
       write(6,*)' --- Goto ExSpectrum section --- is=',is
       rewind (ifexsp(is))
       itmx = 0
       do
          read(ifexsp(is),*,end=1215)ipex,itpex,itex,qqex(1:3), eex,exsp
          if(itex>itmx) itmx=itex
       enddo
1215   continue
       nspexmx = itmx*(nqbz+nq0i*ngrp) !Get marimum value of the number of the ex spectrum
       allocate( eex1(nspexmx,ntq,nq), exsp1(nspexmx,ntq,nq),nspex(ntq,nq),itex1(nspexmx,ntq,nq), &
            qqex1(3,nspexmx,ntq,nq) )
       write(6,*)' nspexmx =',nspexmx
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
1216   continue               !Get eex1(1:nspex) exsp1(1:nspex) for itp ip.
       write(6,*)' nspex(1 1)=',nspex(1,1)
       do ipex = 1,nq
          do itpex=1,ntq
             write(6,*)' is itq ip =',is,itq,ip
             nnex = nspex(itpex,ipex)
             allocate( ieord(1:nnex) )
             call sortea( eex1(1:nnex,itpex,ipex),ieord, nnex,isig)
             eex1 (1:nnex,itpex,ipex) = eex1 (ieord(1:nnex),itpex,ipex)
             exsp1(1:nnex,itpex,ipex) = exsp1(ieord(1:nnex),itpex,ipex)
             itex1(1:nnex,itpex,ipex) = itex1(ieord(1:nnex),itpex,ipex)
             qqex1(:,1:nnex,itpex,ipex)= qqex1 (:,ieord(1:nnex),itpex,ipex)
             open(newunit=ifexspx, file='EXSP'//charnum3(ipex)//charnum3(itpex)//'.'//char(48+is))
             open(newunit=ifexspxx,file='EXSS'//charnum3(ipex)//charnum3(itpex)//'.'//char(48+is))
             do i=1,nnex
                write(ifexspx,"(2d14.6, i4, 3f14.6)")eex1(i,itpex,ipex),exsp1(i,itpex,ipex),&
                     itex1 (i,itpex,ipex),qqex1 (1:3,i,itpex,ipex)
             enddo
             eee  =-1d99
             exwgt= 0d0
             do i=1,nnex
                if(eex1(i,itpex,ipex) > eee+1d-4 .OR. i==nnex) then
                   if(i/=1) write(ifexspxx, "(2d23.15)") eee, exwgt*hartree
                   eee  = eex1(i,itpex,ipex)
                   exwgt= exsp1 (i,itpex,ipex)
                else
                   exwgt= exwgt + exsp1 (i,itpex,ipex)
                endif
             enddo
             deallocate( ieord )
             close(ifexspx)
             close(ifexspxx)
          enddo
       enddo
       deallocate( eex1, exsp1, nspex, itex1, qqex1 )
    enddo
    write(6,*)' End of ExSpectrum section ---'
  end subroutine Hswrite4
end subroutine hsfp0

subroutine rsexx2 (nspin, itq, q, ntq,nq,ginv, symgg,ng, vxco)
  implicit real*8 (a-h,o-z)
  implicit integer (i-n)
  dimension vxco(ntq,nq,nspin),q(3,nq),itq(ntq) !itq is not dependent on q, right?
  real(8),allocatable :: qqq(:,:),vxcfpx(:,:,:)
  logical ::nocore,lfind
  real(8)::  rydberg,tolq=1d-5,qx(3),ginv(3,3),qr(3),symgg(3,3,ng),sym(3,3)
  integer:: ikpx=999999
  write(6,*)' OPEN VXCFP '
  open(newunit=ifvxcfp,file='VXCFP',form='unformatted')
  read(ifvxcfp) ldim,nqbz
  write(6,*)' rsexx ldim,nqbz',ldim,nqbz
  allocate(qqq(3,nqbz),vxcfpx(ldim,nqbz,nspin))
  do ikp = 1,nqbz
     read(ifvxcfp) qqq(1:3,ikp),vxcfpx(1:ldim,ikp,1:nspin)
     write(6,"(i5,100d13.5)") ikp,qqq(1:3,ikp)
  enddo
  close(ifvxcfp)
  do iq=1,nq
     do ikp=1,nqbz
        do ig = 1,ng
           sym = symgg(:,:,ig)
           qr=matmul(sym,qqq(1:3,ikp))
           lfind=.false.
           if(sum( (qr-q(1:3,iq))**2) <tolq) then
              lfind=.true.
           else
              call rangedq( matmul(ginv,q(1:3,iq)-qr), qx)
              if(sum(abs(qx))< tolq) lfind= .TRUE. 
           endif
           if(lfind) then
              ikpx=ikp
              goto 100
           endif
        enddo
     enddo
     call rx( ' rsexx: not find ikp')
100  continue
     vxco(1:ntq,iq,1:nspin)=rydberg()*vxcfpx(itq(1:ntq),ikpx,1:nspin)
  enddo
end subroutine rsexx2