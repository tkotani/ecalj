program lmfham2
  ! MLO version of hmaxloc.
  ! 
  !   
  !-------------------------------------------------------------
  ! construct maximally localized Wannier functions
  
  ! References
  ! [1] N. Marzari and D.Vanderbilt, PRB56,12847(1997)
  ! [2] I. Souza, N. Marzari and D.Vanderbilt, PRB65,035109(2002)
  
  ! mode 1:  determine parameters for <u(m,k)|u(n,k+b> (uu-matrix)
  ! mode 2:  main part
  !    Step 1: choose Hilbert space  (Ref.[2])
  !    Step 2: maximally localize Wannier functions (Ref.[1])
  !    Step 3: construct effective Hamiltonian and interpolate bands (Ref.[2])
  
  !m Oct 2008 Takashi Miyake, updated
  !m Aug 2007 Takashi Miyake, berry connection in the Wannier gauge
  !  May 2004 Takashi Miyake, from hwmat.f
  !------------------------------------------------------------
  !     use m_readqg,only: readngmx,readqg
  use m_ftox
  use m_lgunit,only: stdo
!  use m_hamindex,only:   Readhamindex,symgg=>symops,ngrp
!  use m_hamindex,only: qtt,nqtt
!  use m_hamindex0,only: readhamindex0,iclasst
!  use m_readeigen,only: init_readeigen,init_readeigen2,readeval
  
!  use m_read_bzdata,only: read_bzdata, nqbz,n1,n2,n3,qbz !,wbz
  
!  use m_qbze,only: Setqbze, nqbze,qbze
  use m_genallcf_v3,only: genallcf_v3, nclass,natom,nspin,nl,nn, &
       nlmto,nlnmx, nctot,niw, alat,delta,deltaw,esmr,clabl,iclass, il, in, im, nlnm, &
       plat, pos, ecore, konf,z, spid
  use m_read_Worb,only: s_read_Worb, s_cal_Worb, &
       nwf, nclass_mlwf, cbas_mlwf, nbasclass_mlwf, &
       classname_mlwf, iclassin, &
       iphi, iphidot, nphi, nphix
  use m_keyvalue,only: getkeyvalue
!  use m_readhbe,only:Readhbe,nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg
  use m_HamPMT,  only: npair,nlat,nqwgt, ldim, nkp, qbz=>qplist,ib_table, ReadHamPMTInfo,n1=>nkk1,n2=>nkk2,n3=>nkk3
  use m_HamRsMTO,only: hammr,ovlmr,ndimMTO,  ReadHamRsMTO,ixa=>ix,npairmx,nspx
  use m_zhev,only:zhev_tk4
  use m_MPItk,only:    m_MPItk_init, m_MPItk_finalize, nsize, master_mpi
  implicit none
  real(8),allocatable:: r0g(:,:), wphi(:,:)
  real(8)    :: esmr2,shtw
  integer :: iclass2
  integer(4):: &
       ixc,iopen,ibas,ibasx,ngpmx,nxx,ngcmx,nbloch,ifqpnt,ifwd,ifbb, &
       nprecx,mrecl,nblochpmx2,nwt,niwt, nqnum,mdimx,nblochpmx, &
       ifrcw,ifrcwi,  noccxv,maxocc2,noccx,ifvcfpout,iqall,iaf,ntq, &
       i,j,k,nspinmx, nq,is,ip,iq,idxk,ifoutsex,iclose,nq0i,ig, &
       mxkp,ntet,nene,iqi, ix,iw, &
       nlnx4,niwx,irot,invr,invrot,ivsum, ifoutsec,ntqx, &
       ifmlw(2),ifmlwe(2),ifxc(2),ifsex(2), ifphiv(2),ifphic(2),ifec,ifexsp(2), &
       ifsecomg(2),ifexx,ifwand,ndble=8
  real(8) :: tpia,vol,voltot,rs,alpha, &
       qfermi,efx,valn,efnew,edummy,efz,qm,xsex,egex, &
       zfac1,zfac2,dscdw1,dscdw2,dscdw,zfac,ef2=1d99,exx,exxq,exxelgas
  logical :: lqall,laf
  real(8),parameter:: pi = 4d0*datan(1d0)
  integer,allocatable :: itq(:)
  real(8),allocatable    :: q(:,:)
  integer(4),allocatable :: ngvecpB(:,:,:),ngvecp(:,:), ngvecc(:,:),iqib(:),&
       kount(:,:), nx(:,:),nblocha(:),lx(:) !ngveccBr(:,:,:)
  real(8),allocatable:: vxcfp(:,:,:), &
       wqt(:), wgt0(:,:),q0i(:,:), &
       ppbrd (:,:,:,:,:,:,:),cgr(:,:,:,:),eqt(:), &
       ppbrdx(:,:,:,:,:,:,:),aaa(:,:), ppb(:), eq(:), eqx(:,:,:),eqx0(:,:,:),ekc(:),coh(:,:) &
       , rw_w(:,:,:,:,:),cw_w(:,:,:,:,:), rw_iw(:,:,:,:,:),cw_iw(:,:,:,:,:)
  complex(8),allocatable:: geigB(:,:,:,:)
  logical :: screen, exchange, cohtest, legas, tote
  real(8) ::  rydberg,hartree
  real(8):: qreal(3), ntot,tripl!,xxx(3,3) !,nocctotg2
  real(8):: qlat(3,3)
  logical ::nocore

  !     okumura Sep,2017
  integer(4):: ifdorb,iiwf
  integer(4),allocatable::idorb(:),idorb_in(:) !check

  ! space group infermation
  integer(4),allocatable ::  invgx(:), miat(:,:)
  real(8),allocatable    :: tiat(:,:,:),shtvg(:,:)


  real(8),allocatable   :: eex1(:,:,:),exsp1(:,:,:),qqex1(:,:,:,:)
  integer(4),allocatable:: nspex(:,:),ieord(:),itex1(:,:,:)
  real(8)    :: qqex(1:3), eex,exsp,eee, exwgt,deltax0
  integer(4) :: itmx,ipex,itpex,itex,nspexmx,nnex,isig,iex,ifexspx &
       ,ifexspxx ,ifefsm, nq0ix,ifemesh,nz
  character(3)  :: sss
  character(12) :: filenameex
  logical :: exspwrite=.false.
  character(8) :: xt


  integer(4)::ini,nq0it,idummy
  !      real(8),allocatable:: qbze(:,:)

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

  ! for maxloc
  real(8)   :: wbb(12),wbbsum,bb(3,12), &
       eomin,eomax,eimin,eimax, &
       qwf0(3),dqwf0(3),qks(3),q0(3)
  complex(8),allocatable:: uumat(:,:,:,:),evecc(:,:),eveccs(:,:), &
       amnk(:,:,:),cnk(:,:,:),umnk(:,:,:),evecc1(:,:,:),evecc2(:,:,:)
  real(8),allocatable:: ku(:,:),kbu(:,:,:),eunk(:,:),eval(:),evals(:), &
       eks(:)!,rt(:,:),rt8(:,:,:),qbz0(:,:)
  integer(4):: nbb,isc,ifq0p, &
       nox,iko_ix,iko_fx, &
       noxs(2),iko_ixs(2),iko_fxs(2), &
       ieo_swt,iei_swt,itin_i,itin_f,itout_i,itout_f, &
       nbbelow,nbabove
  integer(4),allocatable:: ikbidx(:,:)
  integer(4),allocatable:: iki_i(:),iki_f(:), &
       ikbi_i(:,:),ikbi_f(:,:), &
       iko_i(:),iko_f(:), &
       ikbo_i(:,:),ikbo_f(:,:)
  logical :: leout,lein,lbin,lq0p,lsyml,lbnds
  logical :: debug=.false.

  integer(4):: nlinex!,ntmp
  parameter (nlinex=100)
  integer(4)::nline,np(nlinex)
  real(8):: qi(3,nlinex),qf(3,nlinex)
  ! step 1
  complex(8),allocatable:: cnq0(:,:), &
       upu(:,:,:,:),cnk2(:,:,:), &
       zmn(:,:)
  complex(8):: ctmp
  real(8),allocatable:: omgik(:)
  real(8)   :: omgi,omgiold,conv1,alpha1,domgi,qtmp(3)
  integer(4):: nsc1,ndz,nin,ifhoev,ifuu0,ifpsig
  ! step 2
  complex(8),allocatable:: mmn(:,:,:,:),mmn0(:,:,:,:), &
       rmn(:,:),smn(:,:),amn(:,:), &
       tmn(:,:),dwmn(:,:)
  real(8),allocatable:: rn(:,:),qn(:)
  complex(8),allocatable::r_nm(:,:,:)
  real(8)   :: omgd,omgod,omgdod,omgidod,omgdodold,domgdod, &
       conv2,alpha2
  integer(4):: nsc2,ibb,ii,ij,ik
  logical :: lrmn,lmmn
  ! step 3
  complex(8),allocatable:: hrotk(:,:,:),hrotr(:,:,:),hrotkp(:,:) &
       , hrotkps(:,:)
  real(8):: e1,e2,rcut
  integer(4):: iband,ifbnd,iftb,ifsh,nsh,nsh1,nsh2,iffb
  logical :: lsh
  real(8),allocatable :: rws(:,:),drws(:)
  integer(4),allocatable:: irws(:)
  integer(4):: nrws,ifham

  ! ixc=3
  character(20)::filename
  complex(8),allocatable::  hrotrcut(:,:,:),evecc_w(:,:,:,:)
  integer:: ifh
  real(8):: heps ,r_v

  real(8)::qold(3)
  real(8),allocatable:: xq(:),eval1(:,:),eval2(:,:),eval3(:,:),eval_w(:,:,:)

  integer::npin
  real(8):: qiin(3),qfin(3)

  integer(4),allocatable:: &
       m_indx(:),n_indx(:),l_indx(:),ibas_indx(:),ibasiwf(:)
  integer:: ifoc,iwf,ldim2,ixx,ifile_handle

  real(8):: enwfmax,qxx(3),eeee,enwfmaxi, ef,epsovl=1d-8,ginv(3,3)
  integer:: inii,if102,iwf2,ib,itmp,itmp2,nqbz2,nspin2,ib1,ib2,iqb,iqbz,it,jsp,nmx,nev,isyml,nband,nqbz!,n1,n2,n3
  logical:: leauto,leinauto,iprint,cmdopt2
  complex(8):: cccx(7)
  character(20):: outs=''
  
  complex(8),allocatable::ovlm(:,:),ovlmx(:,:),hamm(:,:),evec(:,:),ovec(:,:),emat(:,:)
  real(8),allocatable:: evl(:,:),ovl(:)
  real(8),allocatable:: bbv(:,:),wbz(:)
  complex(8),parameter:: img=(0d0,1d0)
  !  integer,allocatable:: ikbidx(:,:)
  complex(8),allocatable:: hmlor(:,:,:,:),omlor(:,:,:,:)
  call m_MPItk_init('lmfham2') ! mpi initialization
  hartree=2d0*rydberg()
  write(6,*) ' verbose=',verbose()
  write(6,*) '  Oneof --job=? 1:preparation, 2:main '
  if(cmdopt2('--job=',outs)) then
     read(outs,*) ixc
  else
     write(6,*)'set --job=1 or 2'
     call exit(-1)
  endif
  write(6,*) ' ixc=',ixc
  if(ixc<1 .OR. ixc>3) call rx(' --- ixc=0 --- Choose computational mode!')
  call ReadHamPMTinfo()
  call ReadHamRsMTO()
  nqbz=nkp
!  call read_Bzdata() !  write(6,*)' nqibz ngrp=',nqibz,ngrp,' nqbz  =',nqbz
!  call Setqbze()    ! extented BZ points list
  call pshpr(60)
  incwfin= -1  !use 7th colmn for core at the end section of GWIN
  call genallcf_v3(incwfin) !in module m_genallcf_v3
  call getnemx(nbmx,ebmx,8,.true.) !8+1 th line of GWIN0
  if (nclass /= natom ) stop ' hsfp0: nclass /= natom ' ! We assume nclass = natom.
  write(6,*)' hsfp0: end of genallcf2'
  call pshpr(30)
  tpia = 2d0*pi/alat
  call minv33tp(plat,qlat)
  ginv = transpose(plat)
  voltot = abs(alat**3*tripl(plat,plat(1,2),plat(1,3)))
  efermiread: block
    integer:: ifi
    open(newunit=ifi,file='efermi.lmf')
    read(ifi,*) ef
    close(ifi)
  endblock efermiread
  readsyml: block
    write(*,*)'Read k points for bands from SYML'     !      lqall      = .false.    !      laf        = .false.
    open(newunit=isyml,file='syml',status='old')
    nline=0
    do i = 1,nlinex
       read(isyml,*,err=551,end=552)npin,qiin,qfin
       if (npin==0) exit
       nline = nline+1
       np(nline)=npin
       qi(1:3,nline)=qiin
       qf(1:3,nline)=qfin
551    continue
    enddo
552 continue
    if (nline == nlinex) call rx('hmaxloc: too many lines in SYML')
    close(isyml)
    nq = 0
    do i = 1,nline
       nq = nq + np(i)
    enddo
    allocate(q(3,nq),xq(nq))
    iq = 0
    xq=0d0
    qold=q(:,1)
    nlineloop: do i = 1,nline
       nploop: do j = 0,np(i)-1
          iq = iq + 1
          q(:,iq) = qi(:,i) + (qf(:,i)-qi(:,i))*dble(j)/dble(np(i)-1)
          if(iq>1) then  !!! xq may be too long if q is outside of 1stBZ
             if(0.2 < dsqrt(sum((q(:,iq)-qold)**2))) then !!maybe same point (must be checked)
                xq(iq)= xq(iq-1)
             else
                xq(iq)= xq(iq-1) + dsqrt( sum((q(:,iq)-qold)**2) )
             endif
          endif
          qold=q(:,iq)
       enddo nploop
    enddo nlineloop
  endblock readsyml
  call s_read_Worb() ! input parameters specific to MAXLOC
  call s_cal_Worb()
  allocate (r0g(nphix,nwf), wphi(nphix,nwf))
  r0g = 2d0
  wphi = 1d0
  !! Read parameters in GWinput. expand wan_input
  !      call wan_input(leout,lein,lbin,ieo_swt,iei_swt,
  !     &    eomin,eomax,itout_i,itout_f,nbbelow,nbabove,
  !     &    eimin,eimax,itin_i,itin_f,
  !     &    nsc1,nsc2,conv1,conv2,alpha1,alpha2,rcut)
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
     if (eomin > eomax) call rx('hmaxloc: eomin > eomax')
     ieo_swt = 1
  else
     call getkeyvalue("GWinput","wan_out_bmin",itout_i,default=999 )
     call getkeyvalue("GWinput","wan_out_bmax",itout_f,default=-999 )
     if (itout_i > itout_f) call rx('hmaxloc: itout_i > itout_f')
  endif
  if (lein) then
     call getkeyvalue("GWinput","wan_in_emin",eimin,default=999d0 )
     call getkeyvalue("GWinput","wan_in_emax",eimax,default=-999d0 )
     if (eimin > eimax) call rx('hmaxloc: eimin > eimax')
     iei_swt = 1
  endif
  if (lbin) then
     call getkeyvalue("GWinput","wan_in_bmin",itin_i,default=999 )
     call getkeyvalue("GWinput","wan_in_bmax",itin_f,default=-999 )
     if (itin_i > itin_f) call rx('hmaxloc: itin_i > itin_f')
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
  r_v=rcut
  call getkeyvalue("GWinput",'wan_tbcut_rcut',heps,default=r_v)
  call getkeyvalue("GWinput",'wan_tbcut_heps',heps,default=0.0d0)
  write(*,*) 'mloc.heps ', heps   !c --- read LDA eigenvalues
!  ntq = nwf
!  call winfo(6,nspin,nq,ntq,is,nbloch ,0,0,nqbz,nqibz,ef,deltaw,alat,esmr)
! Rt vectors
!  allocate (rt(3,nqbz),rt8(3,8,nqbz),qbz0(3,nqbz))
!  call getrt(qbz,qbas,plat,n1,n2,n3,nqbz, rt,rt8,qbz0)
  
  call getbb(plat,alat,n1,n2,n3, nbb,wbb,wbbsum,bb) ! b vectors
  allocate (ku(3,nqbz),kbu(3,nbb,nqbz),ikbidx(nbb,nqbz))
  call kbbindx(qbz,ginv,bb, nqbz,nbb, ikbidx,ku,kbu) ! index for k and k+bb
  allocate(iko_i(nqbz),iko_f(nqbz), iki_i(nqbz),iki_f(nqbz), ikbo_i(nbb,nqbz),ikbo_f(nbb,nqbz), ikbi_i(nbb,nqbz),ikbi_f(nbb,nqbz))
  nband= ndimMTO
  if (ixc == 1) then     ! write bb vectors to 'BBVEC'
     call writebb(ifbb,wbb(1:nbb),bb(1:3,1:nbb), ikbidx,ku,kbu, &
          iko_ixs,iko_fxs,noxs, &
          nspin,nqbz,nbb)
     write(6,"(/,a)") 'OK! lmfham2: --job=1'
     call exit(0)
  endif
  allocate(bbv,source=bb)
  ! readbbvec: block
  !   open(newunit=ifbb,file='BBVEC')
  !   read(ifbb,*)
  !   read(ifbb,*)nbb,nqbz2
  !   if (nqbz /= nqbz2) stop 'readbb: nqbz is wrong!'
  !   allocate (bbv(3,nbb))!,ikbidx(nbb,nqbz))
  !   do i = 1,nbb
  !      read(ifbb,*)bbv(1,i),bbv(2,i),bbv(3,i)
  !   enddo
  !   close(ifbb)
  ! endblock readbbvec

  allocate(hmlor(nwf,nwf,npairmx,nspin),omlor(nwf,nwf,npairmx,nspin),source=(0d0,0d0))
  ispinloop: do 1000 is = 1,nspin
     write(*,*)'is =',is,'  out of',nspin
     ! energy window
     iko_i=1; iko_f=36! outer 
     iki_i=1; iki_f=6 ! inner
     ikbo_i=1; ikbo_f=36 !k+b for projection
     lein=.true.
     
     iko_ix=minval(iko_i)
     iko_fx=maxval(iko_f)
     ikbi_i=9999 !dummy
     ikbi_f=9999 !
     nox = iko_fx - iko_ix + 1
     
     allocate(ovlm(1:ndimMTO,1:ndimMTO),ovlmx(1:ndimMTO,1:ndimMTO),hamm(1:ndimMTO,1:ndimMTO))
     allocate(evec(ndimMTO,ndimMTO),evl(ndimMTO,nqbz), ovec(ndimMTO,ndimMTO),ovl(ndimMTO))
     allocate(amnk(iko_ix:iko_fx,nwf,nqbz))
     allocate(wbz(nqbz),source=1d0/nqbz)
     write(6,*)'goto uumatambk block ndimMTO=',ndimMTO,iko_ix,iko_fx
     uumatambk: block !Get uumat and amnk
       real(8):: qp(3)
       complex(8):: emat(1:ndimMTO,1:ndimMTO),osq(1:ndimMTO,1:ndimMTO),o2al(1:ndimMTO,1:ndimMTO,nqbz),phase,ovlmm(ndimMTO,nwf),&
            evec(ndimMTO,ndimMTO,nqbz),evecx(1:ndimMTO,1:ndimMTO)
       ! integer,parameter:: nred=18
       ! complex(8):: ovlmi(1:ndimMTO,1:ndimMTO), ovv(1:ndimMTO,1:nred),hammm(nred,nred),ovvv(nred,nred)
       ! real(8):: evlx(nred)
       ! integer:: ired(nred)
       ! open(newunit=iband,file='band_lmfham2init.dat')
       ! ired=[(i,i=1,9),(i,i=19,27)]
       ! do iq=1,nq
       !    qp=q(:,iq)
       !    ovlm = 0d0
       !    hamm = 0d0
       !    do i=1,ndimMTO
       !       do j=1,ndimMTO
       !          ib1 = ib_table(ixa(i)) !atomic-site index in the primitive cell
       !          ib2 = ib_table(ixa(j))
       !          do it =1,npair(ib1,ib2)
       !             phase= 1d0/nqwgt(it,ib1,ib2)*exp(-img*2d0*pi* sum(qp*matmul(plat,nlat(:,it,ib1,ib2))))
       !             hamm(i,j)= hamm(i,j)+ hammr(i,j,it,is)*phase
       !             ovlm(i,j)= ovlm(i,j)+ ovlmr(i,j,it,is)*phase
       !          enddo
       !       enddo
       !    enddo
       !    ovlmi=ovlm
       !    call matcinv(ndimMTO,ovlmi)
       !    ovv(:,:)= matmul(ovlmi,ovlm(:,ired(:)))
       !    hammm = matmul(dconjg(transpose(ovv)),matmul(hamm,ovv))
       !    ovvv  = matmul(dconjg(transpose(ovv)),matmul(ovlm,ovv))
       !    nmx=nred
       !    call zhev_tk4(nmx,hammm,ovvv,nmx,nev, evlx, evecx, epsovl)!Diangonale (hamm- evl ovlm) z=0
       !    do i=1,nev
       !       write(iband,ftox) ftof(xq(iq)),ftof(evlx(i)), is,i
       !    enddo
       ! enddo
       ! close(iband)
       ! stop 'vvvvvvvvvvv'
       
       emat=0d0
       forall(i=1:ndimMTO) emat(i,i)=1d0
       do iqbz=1,nqbz
          qp=qbz(:,iqbz)
          ovlm = 0d0
          hamm = 0d0
          do i=1,ndimMTO
             do j=1,ndimMTO
                ib1 = ib_table(ixa(i)) !atomic-site index in the primitive cell
                ib2 = ib_table(ixa(j))
                do it =1,npair(ib1,ib2)
                   phase= 1d0/nqwgt(it,ib1,ib2)*exp(-img*2d0*pi* sum(qp*matmul(plat,nlat(:,it,ib1,ib2))))
                   hamm(i,j)= hamm(i,j)+ hammr(i,j,it,is)*phase
                   ovlm(i,j)= ovlm(i,j)+ ovlmr(i,j,it,is)*phase
                enddo
             enddo
          enddo
          nmx=ndimMTO
          ovlmx=ovlm
          call zhev_tk4(ndimMTO,ovlm,emat,nmx,nev, ovl, ovec, epsovl)!Diangonale (ovlm - e ) alp=0
          ovlm=ovlmx
          call zhev_tk4(ndimMTO,hamm,ovlm,nmx,nev, evl(:,iqbz), evec(:,:,iqbz), epsovl)!Diangonale (hamm- evl ovlm) z=0
!          do i=1,nwf
!             write(6,*)'   e0=',i,evl(i,iqbz)
!          enddo
          do concurrent (i=1:ndimMTO,j=1:ndimMTO)
             osq(i,j)=sum(ovec(i,:)*ovl(:)**0.5d0*dconjg(ovec(j,:)))
          enddo   
          o2al(:,:,iqbz) = matmul(osq, evec(:,:,iqbz)) !o2al(basis index, band index, iqbz index)
!          o2al(:,:,iqbz) = matmul(ovlmx,evec(:,:,iqbz)) !matmul(osq, evec(:,:,iqbz)) !o2al(basis index, band index, iqbz index)
          forall(i=1:ndimMTO) ovlmm(i,:) = [ovlmx(i,1:9),ovlmx(i,19:27)]
          amnk(:,:,iqbz)= matmul(transpose(dconjg(evec(1:ndimMTO,iko_ix:iko_fx,iqbz))),ovlmm) !amnk= <psi|MTO> minimum basis MTO =9+9
       enddo
       allocate (uumat(iko_ix:iko_fx,iko_ix:iko_fx,nbb,nqbz))
       do iqbz=1,nqbz
          do ibb=1,nbb
             iqb = ikbidx(ibb,iqbz)             !q1(:) = qbz(:,iqbz)             !q2(:) = q1(:) + bbv(:,ibb)
             do concurrent(ib1=iko_ix:iko_fx, ib2=iko_ix:iko_fx) !ib1,ib2 band index of outer-inner window
                uumat(ib1,ib2,ibb,iqbz)= sum(dconjg(o2al(1:ndimMTO,ib1,iqbz))*o2al(1:ndimMTO,ib2,iqb)) ! define connection <q ib1| q+b ib2>
!uumat(ib1,ib2,ibb,iqbz)= sum(dconjg(o2al(1:ndimMTO,ib1,iqbz))*evec(1:ndimMTO,ib2,iqb)) !o2al(1:ndimMTO,ib2,iqb)) ! define connection <q ib1| q+b ib2>evec
             enddo
          enddo
       enddo
     endblock uumatambk
     
     !! step 1  -- choose Hilbert space -- determine cnk
     write(*,*)'Step 1: Hilbert space branch'
     write(6,*)' iko_ix iko_fx=',iko_ix,iko_fx
     allocate (& !amnk(iko_ix:iko_fx,nwf,nqbz), &
          upu(iko_ix:iko_fx,iko_ix:iko_fx,nbb,nqbz), &
          cnk(iko_ix:iko_fx,nwf,nqbz), &
          cnk2(iko_ix:iko_fx,nwf,nqbz), &
          omgik(nqbz))
     ! amnk appered in Eq.22 in Ref.II. <psi|Gaussian>
     cnk  = 0d0
     ! amnk = 0d0
     ! read psig(it,iwf,iqbz) = < psi(it,iqbz) | g(iwf) >
     !print *,'goto readpsig'
     !call readpsig(is,  iko_ix,iko_fx, nqbz,nwf,   amnk) !amnk given already.
     call amnk2unk(amnk,iko_ix,iko_fx,iko_i,iko_f, nwf,nqbz,  cnk)
     print *,'goto init_iew=',lein,iko_ix,iko_fx !,iko_i,iko_f,iki_i,iki_f
     if(lein)  call init_iew(iko_ix,iko_fx,iko_i,iko_f, iki_i,iki_f, nwf,nband,nqbz, cnk) ! inner energy window
     write(6,*)'end of amnk2'
     
!     call init_unkg(is,qbz,ginv,ef,lein, &
!          iko_ix,iko_fx,iko_i,iko_f, &
!          iki_i,iki_f, &
!          nwf,nband,nqbz, &
!          amnk,cnk)
  
     !      call chk_amnkweight(qbz,iko_ix,iko_fx,amnk,
     !     &     nqbz,nwf,nband,nlmto)
     !      call chk_cnkweight(qbz,iko_ix,iko_fx,cnk,
     !     &     nqbz,nwf,nband,nlmto)
     write(*,*) 'gotto step1looop: iko_ix:iko_fx,nwf,nqbz=',iko_ix,iko_fx,nwf,nqbz
     Step1loop: do isc = 1,nsc1
        iqloop: do iq = 1,nqbz
           call dimz(lein,iko_i(iq),iko_f(iq),iki_i(iq),iki_f(iq), ndz,nin)
           !           write(6,*)'iq =',iq,'iko=',iko_i(iq),iko_f(iq),iki_i(iq),iki_f(iq),'ndz=',ndz
           if (nwf > nin) then
              if (ndz < 1) call rx('ndz < 1')
              ! (1-2) <u_mk | P_k+b | u_nk>
              call getupu(isc, &
                   uumat(:,:,:,iq),cnk, &
                   lein,alpha1,iq,ikbidx(:,iq), &
                   iko_ix,iko_fx, &
                   iko_i(iq),iko_f(iq), &
                   iki_i(iq),iki_f(iq), &
                   ikbo_i(:,iq),ikbo_f(:,iq), &
                   ikbi_i(:,iq),ikbi_f(:,iq), &
                   nwf,nbb,nqbz, &
                   upu(:,:,:,iq))
              ! (1-3) Zmn(k) > phi,eval
              allocate (zmn(ndz,ndz),evecc(ndz,ndz),eval(ndz))
              call getzmn(upu(:,:,:,iq),wbb,lein, &
                   iko_ix,iko_fx, &
                   iko_i(iq),iko_f(iq), &
                   iki_i(iq),iki_f(iq), &
                   nwf,nbb,nqbz,ndz, &
                   zmn)
              call chk_hm(zmn,ndz)
              call diag_hm(zmn,ndz,eval,evecc)
              call new_cnk(cnk(:,:,iq),evecc,iq, &
                   iko_ix,iko_fx, &
                   iko_i(iq),iko_f(iq), &
                   iki_i(iq),iki_f(iq), &
                   nwf,ndz, &
                   cnk2(:,:,iq))
              ! (1-3) w_I(k)  eq.(18)
              !call chk_eval(wbb,eval,nbb,ndz)
              call get_omgik(wbb,eval, &
                   iko_i(iq),iko_f(iq), &
                   iki_i(iq),iki_f(iq), &
                   nbb,nwf,ndz, &
                   omgik(iq))
              deallocate (zmn,evecc,eval)
           else
              omgik(iq) = 0d0
              cnk2(:,:,iq) = cnk(:,:,iq)
              ! end if (ndz>1)
           endif
        enddo iqloop
        ! (1-5) w_I(k) > Omaga_I  eq.(11)
        omgi = sum(omgik(:)*wbz(:))
        ! (1-6) check self-consistency
        write(*,"('#SC-loop, conv.',i5,d13.5)")isc,omgi
        if (isc >= 2) then
           domgi = dabs((omgiold - omgi) / omgiold)
           if (domgi < conv1) then
              write(*,*) 'step1: converged!'
              goto 810
           endif
        endif
        ! update
        omgiold = omgi
        cnk     = cnk2
     enddo Step1loop! end of self-consistent loop
     write(*,*)'step1: not converged'
810  continue
     deallocate(upu,cnk2)
     !! NOTE: cnk is the final results of step 1
     !!   cnk(iko_ix:iko_fx,nwf,nqbz)
     !!   cnk(iko_i(iq):iko_f(iq),nwf,iq) gives nwf-dimentional space.
     !!   step 1 (minimization of Omega_I)

     write(6,*)'nnnn nkp nqbz=',nkp,nqbz
     eigenw: block ! 
       real(8):: qq(3)
       integer:: il,im,in,ib1,ib2,jsp
       complex(8):: ham(nwf,nwf),ovlx(nwf,nwf),phase,proj(iko_ix:iko_fx,iko_ix:iko_fx),pa(iko_ix:iko_fx,nwf)
       jsp=is
       do iqbz = 1,nqbz
          qq(:) = qbz(:,iqbz)
!          write(6,*)' xxxx iq q=',iq,qq
          do concurrent(i=iko_ix:iko_fx, j=iko_ix:iko_fx)     !outer bandindex
             proj(i,j) = sum(cnk(i,:,iqbz)*dconjg(cnk(j,:,iqbz))) !sum for MLOindex
          enddo
          pa(iko_ix:iko_fx,1:nwf) = matmul(proj,amnk(iko_ix:iko_fx,1:nwf,iqbz))
          do concurrent(i=1:nwf,j=1:nwf)
             ham(i,j)  = sum(dconjg(pa(:,i))*evl(:,iqbz)*pa(:,j)) !sum(dconjg(pa(:,i))*pa(:,j))!
             ovlx(i,j) = sum(dconjg(pa(:,i))*pa(:,j))
          enddo

          ! ham=0d0
          ! ovlx=0d0
          ! forall(i=1:nwf)
          !    ham(i,i)=1d0
          !    ovlx(i,i)=1d0
          ! end forall
          !          write(6,*)'sum ham=',iqbz,sum(abs(ham))
          !          write(6,*)'sum ovl=',iqbz,sum(abs(ovl))
          
          do concurrent(i=1:nwf,j=1:nwf) !Real space Hamiltonian. H(k)->H(T) FT to real space
             ib1 = merge(1,2,i<=9) !test ib_table(i) 
             ib2 = merge(1,2,j<=9) !     ib_table(j) 
             do it =1,npair(ib1,ib2)! hammr_ij (T)= \sum_k hamm(k) exp(ikT). it is the index for T
                phase = 1d0/dble(nkp)* exp(img*2d0*pi* sum(qq*matmul(plat,nlat(:,it,ib1,ib2))))
                hmlor(i,j,it,jsp)= hmlor(i,j,it,jsp)+ ham(i,j)*phase
                omlor(i,j,it,jsp)= omlor(i,j,it,jsp)+ ovlx(i,j)*phase
             enddo
          enddo
       enddo
     endblock eigenw
     write(6,*)' get hamlor. Goto band_lmfham2.dat ---------'
     bandplotMLO: block
       real(8):: qp(3),evlm(nwf,nq)
       complex(8):: phase,hamm(nwf,nwf),ovlm(nwf,nwf)
       integer:: iband
       open(newunit=iband,file='band_lmfham2.dat')
       do iq= 1,nq! nq !syml
          qp = q(:,iq)
          ovlm = 0d0
          hamm = 0d0
          do concurrent(i=1:nwf,j=1:nwf)
             ib1 = merge(1,2,i<=9) !test ib_table(i) 
             ib2 = merge(1,2,j<=9) !     ib_table(j) 
             !                ib1 = ib_table(ixa(i)) !atomic-site index in the primitive cell
             !                ib2 = ib_table(ixa(j))
             do it =1,npair(ib1,ib2)
                phase= 1d0/nqwgt(it,ib1,ib2)*exp(-img*2d0*pi* sum(qp*matmul(plat,nlat(:,it,ib1,ib2))))
                hamm(i,j)= hamm(i,j)+ hmlor(i,j,it,is)*phase
                ovlm(i,j)= ovlm(i,j)+ omlor(i,j,it,is)*phase
             enddo
          enddo
!          write(6,*)"sumham",sum(abs(hamm)),sum(abs(ovlm))
          nmx =nwf
          call zhev_tk4(nwf,hamm,ovlm,nmx,nev, evlm(:,iq), evec, epsovl)! Diangonale (hamm - evl ovlm ) evec=0
          write(stdo,"(' iq q=',i3,*(a))") iq,' ',ftof(qp,3),' e=',ftof(evlm(1:12,iq),6)
          do i=1,nev
             write(iband,ftox)  ftof(xq(iq)),ftof(evlm(i,iq)), is,i
          enddo
       enddo
       close(iband)
     endblock bandplotMLO
     call rx0('mmmm end of bandplotMTO mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm')



     
     
     !! === step 2 -- localize Wannier fn. ================
     write(*,*)'Step 2: Wannier fn. branch'
     allocate (mmn(nwf,nwf,nbb,nqbz),mmn0(nwf,nwf,nbb,nqbz), &
          umnk(nwf,nwf,nqbz), &
          rmn(nwf,nwf),amn(nwf,nwf),smn(nwf,nwf), &
          rn(3,nwf),qn(nwf),tmn(nwf,nwf),dwmn(nwf,nwf),r_nm(3,nwf,nwf), &
          eunk(nwf,nqbz))
     !! (2-0) construct initlal u~ from u
     !! eunk(= e~) of {H~}_mn): eigenvalue within the nwf-dimentional Hilbert space
     ! call diag_unk(is,qbz, &
     !      iko_ix,iko_fx,iko_i,iko_f, &
     !      nband,nwf,nqbz, &
     !      cnk, &
     !      eunk)
     
     !! (2-1) initial: uumat -> M_mn(0) Eq.58 in Ref.[1]
     call init_mmn(cnk,uumat,ikbidx, &
          iko_ix,iko_fx,iko_i,iko_f,ikbo_i,ikbo_f, &
          nwf,nqbz,nbb, &
          mmn0)
     !! (2-2) initial U
     !!    umnk= U(m,n) = ( A S^{-1/2} )_mn. See Eq.23 in Ref.II.
     call init_Umnk(amnk,cnk, &
          iko_ix,iko_fx,iko_i,iko_f, &
          nwf,nqbz, &
          umnk)
     call updt_mmn(umnk,mmn0,ikbidx, &
          nwf,nqbz,nbb, &
          mmn)
     
     ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     !         do i=1,nbb
     !            write(1106+is,"(a,i4,13f13.5)")'bbbb',i,bb(1:3,i),wbb(i)
     !         enddo
     !         do i=1,nqbz
     !            write(1106+is,"(a,i4,13f13.5)")'www',i,wbz(i)
     !         enddo
     ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     ! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      step2loop: do isc = 1, nsc2
!         ! <r_n> ([1] eq.31)
!         call get_rn(mmn,bb,wbb,wbz, &
!              nwf,nqbz,nbb, &
!              rn)
!         !         do i=1,nwf
!         !            write(6,ftox)'wwww: i,i, <wannier i|r|wannier i> ',i,ftof(rn(1:3,i))
!         !         enddo
!         do iq = 1,nqbz
!            dwmn = (0d0,0d0)
!            do ibb = 1,nbb
!               ! cccccccccccccccccccccccccccccccccccccccccccc
!               !      do i = 1,nwf
!               !      do j = 1,nwf
!               !         write(1106+is,"(a,4i5,2f13.3)")' mmmmm ',i,j,ibb,iq,mmn(i,j,ibb,iq)+(0d0,0.0001)
!               !      enddo
!               !      enddo
!               !      do i = 1,nwf
!               !      do j = 1,nwf
!               !         write(1106+is,"(a,4i5,2f13.3)")' nnnnnn ',i,j,ibb,iq,mmn0(i,j,ibb,iq)+(0d0,0.0001)
!               !      enddo
!               !      enddo
!               ! ccccccccccccccccccccccccccccccccccccccccccccccccc

!               ! (2-3) A[R] matrix
!               call getrmn(mmn(:,:,ibb,iq), &
!                    nwf, &
!                    rmn)
!               call getamn(rmn, &
!                    nwf, &
!                    amn)

!               ! (2-4) S[T] matrix
!               call gettmn(rn,mmn(:,:,ibb,iq),bb(:,ibb), &
!                    nwf, &
!                    qn,tmn)
!               call getsmn(tmn, &
!                    nwf, &
!                    smn)

!               ! cccccccccccccccccccccccccccccccc
!               !               smn=0d0
!               !               amn=0d0
!               ! cccccccccccccccccccccccccccccccc

!               ! DW(k) ([1] eq.57)
!               dwmn(:,:) =  dwmn(:,:) &
!                    + wbb(ibb) * (amn(:,:) - smn(:,:)) * alpha2 / wbbsum

!               ! end of ibb-loop
!            enddo

!            ! ccccccccccccccccccccccc
!            !            dwmn=0d0
!            ! ccccccccccccccccccccccc
!            ! (2-5) DW(k) -> U(k) ([1] eq.60)
!            call updt_uk(dwmn, &
!                 nwf, &
!                 umnk(:,:,iq))
!            !            call chk_um(umnk(:,:,iq),nwf)

!            ! ccccccccccccccccccccccccccccccccccccccccccccccccc
!            !      do i = 1,nwf
!            !      do j = 1,nwf
!            !         write(1106+is,"(a,3i5,2f13.3)")' zzzzz',iq,i,j,umnk(i,j,iq)
!            !      enddo
!            !      enddo
!            ! ccccccccccccccccccccccccccccccccccccccccccccccccc

!            ! end of iq-loop
!         enddo



!         ! update Mmn ([1] eq.61)
!         call updt_mmn(umnk,mmn0,ikbidx, &
!              nwf,nqbz,nbb, &
!              mmn)

!         !     (2-6) Omeg_I, Omega_D and Omega_OD ([1] eq.34,35,36)
!         iprint=.false.
!         call getOmg(mmn,rn,bb,wbb,wbz, &
!              nwf,nqbz,nbb, &
!              omgi,omgd,omgod,omgdod,omgidod,alat,iprint)

!         ! check self-consistency
!         !         write(*,*)'#SC-loop, conv.',isc,omgdod
!         !         write(*,950)'Omg: I, OD, D ',omgi,omgod,omgd
!         write(*,"('#SC-loop, spread(bohr**2)',i6,' Omg_I Omg_OD Omg_D=',3f17.10)") &
!                                 !     &   isc,omgdod,omgi,omgod,omgd
!              isc, omgi/tpia**2,omgod/tpia**2,omgd/tpia**2
!         if (isc >= 2) then
!            domgdod = dabs((omgdodold - omgdod) / omgdodold)
!            if (domgdod < conv2) then
!               write(*,*) 'step2: converged!'
!               goto 820
!            endif
!         endif
!         omgdodold = omgdod
!         ! end of self-consistent loop
!      enddo step2loop
!      write(*,*)'step2: not converged'
! 820  continue

!      iprint=.true.
!      call getOmg(mmn,rn,bb,wbb,wbz, nwf,nqbz,nbb, omgi,omgd,omgod,omgdod,omgidod,alat,iprint)
!      !      call chk_dnk(is,eunk,qbz,
!      !     i             umnk,cnk,
!      !     i             iko_ix,iko_fx,iko_i,iko_f,
!      !     d             nband,nwf,nqbz)

!      !      call chk_umn(cnk,umnk,qbz,
!      !     i             iko_ix,iko_fx,iko_i,iko_f,
!      !     d             nwf,nqbz,nband,nlmto)

!      write(*,*)"---------- wlaxloc isp =",is
!      block
!        complex(8),allocatable :: dnk(:,:,:)
!        allocate(dnk(iko_ix:iko_fx,nwf,nqbz))
!        call wmaxloc(ifmlw(is),ifmlwe(is), &
!             qbz,umnk,cnk,eunk, &
!             iko_ix,iko_fx,iko_i,iko_f, &
!             nwf,nqbz,nband,nlmto, is,dnk)
!        !$$$cccccccccccccccccccccccccccccccccccccccccccccccccc
!        !$$$      write(6,ftox)' uuuu'
!        !$$$      do iq=1,nqbz
!        !$$$      do iwf=1,nwf
!        !$$$c         write(6,ftox)' uuuu umnk  ',is,iq,ftof(umnk(iwf,1:nwf,iq),2)
!        !$$$c         write(6,ftox)' uuuu   amnk',is,iq,ftof(amnk(iko_ix:iko_fx,iwf,iq),2)
!        !$$$         do iwf2=1,nwf
!        !$$$           cccx(iwf2)= sum( dconjg(amnk(iko_ix:iko_fx,iwf,iq))*cnk(iko_ix:iko_fx,iwf2,iq) )
!        !$$$         enddo
!        !$$$         write(6,ftox)' uuuu   <G|psi*c  >',is,iq,ftof(cccx,3)
!        !$$$         write(6,ftox)' uuuu   <G|psi*c*u>',is,iq,ftof(matmul(cccx,umnk(:,:,iq)),3)
!        !$$$         do iwf2=1,nwf
!        !$$$            cccx(iwf2)= sum( dconjg(amnk(iko_ix:iko_fx,iwf,iq))*dnk(iko_ix:iko_fx,iwf2,iq))
!        !$$$         enddo
!        !$$$         write(6,ftox)' uuuu   <G|ib>dnk',is,iq,ftof(cccx,3)
!        !$$$      enddo
!        !$$$      enddo
!        !$$$cccccccccccccccccccccccccccccccccccccccccccccccccc
!        deallocate(dnk)
!      endblock
!      call writeOmg(is,mmn,rn,bb,wbb,wbz,tpia,    nwf,nqbz,nbb)     ! 070824
!      call getkeyvalue("GWinput","wan_write_rmn",lrmn,default=.false.)
!      if (lrmn) call writermn(is,mmn,bb,wbb,qbz,qbz0,wbz,rt, nwf,nqbz,nbb,n1,n2,n3)  ! 070830
!      call getkeyvalue("GWinput","wan_write_mmn",lmmn,default=.false.)
!      if (lmmn) call writemmn(is,mmn,bb,wbb,qbz,wbz,rt, nwf,nqbz,nbb,n1,n2,n3)
!      ! ccccccccccccccccccccccccccccccccccccccccccc
!      call get_rnm(mmn,bb,wbb,wbz, &
!           nwf,nqbz,nbb, &
!           r_nm)
!      do i=1,nwf
!         do j=1,nwf
!            if(i==j) cycle
!            write(6,ftox)'dddd: i,j, <wannier i|r|wannier j> ',i,j,ftof(alat/(2d0*pi)*r_nm(1:3,i,j)) &
!                 ,'  ',ftod(alat/(2d0*pi)*sqrt(sum(abs(r_nm(1:3,i,j)**2))))
!         enddo
!      enddo
!      ! ccccccccccccccccccccccccccccccccccccccccccccccc
!      deallocate(uumat,amnk,omgik,mmn,mmn0, rmn,amn,smn,rn,qn,tmn,dwmn,r_nm)

     !! step 3 -- reduced Hamiltonian ------------------------------
     write(*,*)'Step 3: reduced Hamiltonian branch'
     if (is == 1) then
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
     if(allocated(hrotk)) deallocate(hrotk,hrotkp,evecc,eval)
     allocate (hrotk(nwf,nwf,nqbz),hrotkp(nwf,nwf),evecc(nwf,nwf),eval(nwf))! for small Hamiltonian
     call getkeyvalue("GWinput","wan_small_ham",lsh,default=.false.)
     if (lsh) then
        call getkeyvalue("GWinput","wan_nsh1",nsh1, default=1 )
        call getkeyvalue("GWinput","wan_nsh2",nsh2, default=2 )
        write(*,*)'SmallHam on',nsh1,nsh2
        nsh = nsh2 - nsh1 + 1
        if (is == 1) then
           ifsh = iopen('bnds.sh.up',1,-1,0)
        else
           ifsh = iopen('bnds.sh.dn',1,-1,0)
        endif
        write(ifsh,*)nq
        write(ifsh,*)nsh
        allocate (hrotkps(nsh,nsh),eveccs(nsh,nsh),evals(nsh))
     endif
     ! (3-1) ~H(k) -> Hrot(k): note eunk is eigenvalues in the basis of cnk
     call rot_hmnk(umnk,eunk, &
          nwf,nqbz, &
          hrotk) !rotated Hamiltonian in MLW basis.
     ! (3-2) Hrot_mn(R)
     if(allocated(irws)) deallocate(irws,rws,drws)
     allocate(irws(n1*n2*n3*8),rws(3,n1*n2*n3*8),drws(n1*n2*n3*8))
     call wigner_seitz(alat,plat,n1,n2,n3,nrws,rws,irws,drws)
     if(allocated(hrotr)) deallocate(hrotr)
     allocate(hrotr(nwf,nwf,nrws)) !real space Hamiltonian in Wannier funciton basis
     if (ixc == 2) then
        call get_hrotr_ws(hrotk,qbz,wbz, &
             rws,irws,drws, &
             nwf,nqbz,nrws, &
             hrotr)        !     write hrotr and *rws
        if (is == 1) then
           ifh = iopen('hrotr.up',1,-1,0)
        else
           ifh = iopen('hrotr.dn',1,-1,0)
        endif
        call write_hrotr(ifh, hrotr, &
             rws,irws,drws, &
             nwf,nrws )
        close (ifh)
     else if (ixc == 3) then
        if (is == 1) then
           filename='hrotr.up'
        else
           filename = 'hrotr.dn'
        endif
        call read_hrotr(filename,nwf,nrws, &
             hrotr)
        if (is == 1) then
           ifh = iopen('hrotr.cut.up',1,-1,0)
        else
           ifh = iopen('hrotr.cut.dn',1,-1,0)
        endif
        allocate(hrotrcut(nwf,nwf,nrws))
        call make_hrotrcut( hrotr, &
             rws,irws,drws, &
             rcut,heps, &
             nwf,nrws, &
             hrotrcut )
        call write_hrotr(ifh, hrotrcut, &
             rws,irws,drws, &
             nwf,nrws )
        close (ifh)
        deallocate(hrotrcut)
     endif

     
     ! write(ifmlw(is))nqbze,nwf
     ! write(ifmlwe(is))nqbze,nwf
     ! do iq = 1,nqbze
     !    !         write(*,*)'goto get_hrotkp_ws iq=',iq,nqbze
     !    call get_hrotkp_ws(hrotr,rws,drws,irws,qbze(:,iq), nwf,nqbz,nrws, hrotkp)
     !    call diag_hm(hrotkp,nwf,eval,evecc)
     !    call wmaxloc_diag(ifmlw(is),ifmlwe(is), &
     !         iq,qbze(1:3,iq),umnk,cnk,eunk,evecc,eval, &
     !         iko_ix,iko_fx,iko_i,iko_f, &
     !         nwf,nqbz)        !          write(6,*)"iq is",iq,is
     ! enddo

     ! --- Readin nlam index
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

     !      ifh = ifile_handle()
     call write_hopping_output(is, hrotr, &
          rws,irws,alat,plat,qlat,pos,natom, &
          ibasiwf, nwf,nrws,spid , m_indx, l_indx, &
          nphix, iphi, ldim2)
     !      close(ifh)

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
        if(iq<5 .OR. iq>nqbz-3)write(6,*)' got get_hrotkp_ws iq =',iq
        if(iq==5)write(6,*)' ...'
        call get_hrotkp_ws(hrotr,rws,drws,irws,qbz(:,iq), &
             nwf,nqbz,nrws, &
             hrotkp)
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
     ! end okumura

     !! other k-points
     write(ifbnd,*)ef,' ef'
     write(iftb,*)ef,' ef'
     write(iffb,*)'#',ef,' ef'
     if (lsh) write(ifsh,*)ef,' ef'
     allocate(eval1(nwf,nq),eval3(nwf,nq),evecc1(nwf,nwf,nq))
     if(lsh) allocate(eval2(nwf,nq),evecc2(nwf,nwf,nq))
     do iq = 1,nq
        !     write(6,*)' got get_hrotkp_ws iq =',iq
        ! (3-3) Hrot_mn(k')
        call get_hrotkp_ws(hrotr,rws,drws,irws,q(:,iq), &
             nwf,nqbz,nrws, &
             hrotkp)
        ! (3-4) diagonalize
        call diag_hm(hrotkp,nwf,eval,evecc)
        eval1(1:nwf,iq)=eval
        evecc1(1:nwf,1:nwf,iq)=evecc
        !     (3-4) diagonalize  -- Small Hamiltonian --
        if (lsh) then
           hrotkps(1:nsh,1:nsh) = hrotkp(nsh1:nsh2,nsh1:nsh2)
           call diag_hm(hrotkps,nsh,evals,eveccs)
           write(ifsh,*)'iq =',iq
           write(ifsh,990)q(1:3,iq)
           eval2(1:nsh,iq)= evals(1:nsh)
           evecc2(1:nwf,1:nwf,iq)=eveccs
        endif                 ! lsh
        ! (3-3) Hrot_mn(k')  -- Tight-binding ---
        call get_hrotkp_tb_ws(rcut,plat,alat, &
             hrotr,rws,drws,irws,q(:,iq),  ibasiwf,pos,natom, &
             nwf,nqbz,nrws, &
             hrotkp)
        !     (3-4) diagonalize -- Tight-binding --
        call diag_hm(hrotkp,nwf,eval,evecc)
        eval3(1:nwf,iq)=eval
     enddo
     do iband = 1,nwf
        do iq = 1,nq
           write(ifbnd,"(i5,3f13.5,'  ',f13.6,f13.6,i5,' !eee! x eval-ef(ev) iband' )") &
                iq,q(1:3,iq),  xq(iq),(eval1(iband,iq)-ef)*rydberg(),iband
           write(iftb,"(i5,3f13.5,'  ',f13.6,f13.6,i5,' !eee! x eval-ef(ev) iband' )") &
                iq,q(1:3,iq),  xq(iq),(eval3(iband,iq)-ef)*rydberg(),iband
           write(iffb,"(i5,3f13.5,'  ',f13.6,f13.6,i5,' ')",ADVANCE='NO') &
                iq,q(1:3,iq),  xq(iq),(eval1(iband,iq)-ef)*rydberg(),iband
           do iwf=1,nwf
              write(iffb,"(f13.6)",ADVANCE='NO') (abs(evecc1(iwf,iband,iq)))**2
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
              write(ifsh,"(i5,3f13.5,'  ',f13.6,f13.6,i5,' !eee! x eval-ef(ev) iband' )") &
                   iq,q(1:3,iq),  xq(iq),(eval2(iband,iq)-ef)*rydberg(),iband
           enddo
        enddo
     endif
     call writeham(ifham,is,ef,alat,plat,pos,qbz,wbz,rws,irws,hrotk,nspin,natom,nwf,nqbz,nrws)
     deallocate(cnk,umnk,eunk,hrotk,hrotr,hrotkp,evecc,eval,irws,rws,drws, &
          ibasiwf,m_indx,n_indx,l_indx,ibas_indx)
     if (lsh) deallocate(hrotkps,eveccs,evals,evecc2)
     close(ifbnd)
     close(iftb)
     close(iffb)
1000 enddo ispinloop
950 format(a14,3f23.16)
990 format(3f12.6)
  call cputid(0)
  call rx0s('lmfham2: ixc=2 ok')
END PROGRAM lmfham2

! real(8) function nocctotg2(ispin, ef,esmr,qbz,wbz, nband,nqbz)
!   use m_readeigen, only: readeval
!   ! Count the total number of electrons under Ef.
!   ! use readeval
!   ! ispin   = 1, paramagnetic
!   !           2, ferromagnetic
!   ! ef      = fermi level
!   ! nband   = no. states
!   ! nqbz    = no. k-points
!   implicit none
!   integer ::it,is,k,ispin,nqbz,nband
!   real(8):: wbz(nqbz),qbz(3,nqbz),ekt(nband),esmr,ef,wgt,wiocc
!   nocctotg2 = 0d0
!   wgt       = 0d0
!   do is = 1,ispin
!      do k   = 1,nqbz
!         ekt= readeval(qbz(:,k),is)
!         wiocc = 0d0
!         do it = 1,nband
!            if(    ekt(it)  + 0.5d0*esmr < ef  ) then
!               wiocc = wiocc  + 1d0
!            elseif(ekt(it) - 0.5d0*esmr < ef  ) then
!               wiocc  = wiocc + (ef- (ekt(it)-0.5d0*esmr))/esmr
!            endif
!         enddo
!         nocctotg2 = nocctotg2 + wbz(k)* wiocc
!         wgt       = wgt       + wbz(k)
!      enddo
!   enddo
!   if(ispin==1) nocctotg2 = nocctotg2*2
!   write(6,*)' Ef=',ef
!   write(6,*)' wgt nocc=',wgt,nocctotg2
! END function nocctotg2
!------------------------------------------------------------------

!$$$c-----------------------------------------------------------------------
!$$$      subroutine chk_amnkweight(qbz,iko_ix,iko_fx,amnk,
!$$$     &     nqbz,nwf,nband,nlmto)
!$$$      use m_readqg
!$$$      use m_readeigen
!$$$      implicit real*8(a-h,o-z)
!$$$
!$$$      complex(8) :: amnk(iko_ix:iko_fx,nwf,nqbz)
!$$$      complex(8),allocatable:: cphi1(:,:),cphi2(:,:)
!$$$      real(8) :: qbz(3,nqbz),q(3),quu(3)
!$$$      real(8),allocatable:: wbas(:,:)
!$$$      integer(4),allocatable::
!$$$     &  m_indx(:),n_indx(:),l_indx(:),ibas_indx(:)
!$$$
!$$$c --- Readin nlam index
!$$$      ifoc = iopen('@MNLA_CPHI',1,0,0)
!$$$      ldim2 = nlmto
!$$$      read(ifoc,*)
!$$$      allocate(m_indx(ldim2),n_indx(ldim2),l_indx(ldim2),ibas_indx(ldim2))
!$$$      do ix =1,ldim2
!$$$        read(ifoc,*)m_indx(ix),n_indx(ix),l_indx(ix),ibas_indx(ix),ixx
!$$$        if(ixx/=ix) call rx('failed to readin @MNLA_CPHI')
!$$$      enddo
!$$$
!$$$      nbas = ibas_indx(nlmto)
!$$$      allocate(cphi1(nlmto,nband),cphi2(nlmto,nwf),wbas(nbas,nwf))
!$$$      wbas = 0d0
!$$$      cphi2=0d0
!$$$      do iq = 1,nqbz
!$$$      q = qbz(:,iq)
!$$$      call readcphi(q,nlmto,1,quu,cphi1)
!$$$
!$$$       do iwf=1,nwf
!$$$         do ib=iko_ix,iko_fx
!$$$            cphi2(:,iwf) = cphi2(:,iwf) + cphi1(:,ib)*amnk(ib,iwf,iq)
!$$$         enddo
!$$$      enddo
!$$$
!$$$      enddo ! iq
!$$$
!$$$      do iwf=1,nwf
!$$$         do ia=1,nlmto
!$$$            ibas = ibas_indx(ia)
!$$$            wbas(ibas,iwf) = wbas(ibas,iwf) +
!$$$     &                 conjg(cphi2(ia,iwf))*cphi2(ia,iwf)
!$$$         enddo ! ia
!$$$      enddo ! iwf
!$$$      wbas = wbas / dble(nqbz**2)
!$$$
!$$$      write(*,*)'*** ibas,iwf,wbas'
!$$$      do iwf=1,nwf
!$$$      do ibas=1,nbas
!$$$         write(*,*)ibas,iwf,wbas(ibas,iwf)
!$$$      enddo
!$$$      write(*,*)
!$$$      enddo
!$$$      write(*,*)'*** ibas,wbas'
!$$$      do ibas=1,nbas
!$$$         w = 0d0
!$$$         do iwf=1,nwf
!$$$            w = w + wbas(ibas,iwf)
!$$$         enddo
!$$$         write(*,*)ibas,w
!$$$      enddo
!$$$
!$$$      deallocate(cphi1,cphi2,wbas,m_indx,l_indx,n_indx,ibas_indx)
!$$$      ix = iclose('@MNLA_CPHI')
!$$$
!$$$      end
!$$$c-----------------------------------------------------------------------
!$$$      subroutine chk_cnkweight(qbz,iko_ix,iko_fx,cnk,
!$$$     &     nqbz,nwf,nband,nlmto)
!$$$      use m_readqg
!$$$      use m_readeigen
!$$$      implicit real*8(a-h,o-z)
!$$$
!$$$      complex(8) :: cnk(iko_ix:iko_fx,nwf,nqbz)
!$$$      complex(8),allocatable:: cphi1(:,:),cphi2(:,:)
!$$$      real(8) :: qbz(3,nqbz),q(3),quu(3)
!$$$      real(8),allocatable:: wbas(:,:)
!$$$      integer(4),allocatable::
!$$$     &  m_indx(:),n_indx(:),l_indx(:),ibas_indx(:)
!$$$
!$$$c --- Readin nlam index
!$$$      ifoc = iopen('@MNLA_CPHI',1,0,0)
!$$$      ldim2 = nlmto
!$$$      read(ifoc,*)
!$$$      allocate(m_indx(ldim2),n_indx(ldim2),l_indx(ldim2),ibas_indx(ldim2))
!$$$      do ix =1,ldim2
!$$$        read(ifoc,*)m_indx(ix),n_indx(ix),l_indx(ix),ibas_indx(ix),ixx
!$$$        if(ixx/=ix) call rx('failed to readin @MNLA_CPHI')
!$$$      enddo
!$$$
!$$$      nbas = ibas_indx(nlmto)
!$$$      allocate(cphi1(nlmto,nband),cphi2(nlmto,nwf),wbas(nbas,nwf))
!$$$      wbas = 0d0
!$$$      cphi2=0d0
!$$$
!$$$      do iq = 1,nqbz
!$$$      q = qbz(:,iq)
!$$$      call readcphi(q,nlmto,1,quu,cphi1)
!$$$
!$$$      do iwf=1,nwf
!$$$         do ib=iko_ix,iko_fx
!$$$            cphi2(:,iwf) = cphi2(:,iwf) + cphi1(:,ib)*cnk(ib,iwf,iq)
!$$$         enddo
!$$$      enddo
!$$$
!$$$      enddo ! iq
!$$$
!$$$      do iwf=1,nwf
!$$$         do ia=1,nlmto
!$$$            ibas = ibas_indx(ia)
!$$$            wbas(ibas,iwf) = wbas(ibas,iwf) +
!$$$     &                 conjg(cphi2(ia,iwf))*cphi2(ia,iwf)
!$$$         enddo ! ia
!$$$      enddo ! iwf
!$$$      wbas = wbas / dble(nqbz*nqbz)
!$$$
!$$$      write(*,*)'*** ibas,iwf,wbas'
!$$$      do iwf=1,nwf
!$$$      do ibas=1,nbas
!$$$         write(*,*)ibas,iwf,wbas(ibas,iwf)
!$$$      enddo
!$$$      write(*,*)
!$$$      enddo
!$$$
!$$$      write(*,*)'*** ibas,wbas'
!$$$      do ibas=1,nbas
!$$$         w = 0d0
!$$$         do iwf=1,nwf
!$$$            w = w + wbas(ibas,iwf)
!$$$         enddo
!$$$         write(*,*)ibas,w
!$$$      enddo
!$$$      deallocate(cphi1,cphi2,wbas,m_indx,l_indx,n_indx,ibas_indx)
!$$$      ix = iclose('@MNLA_CPHI')
!$$$      end
!-----------------------------------------------------------------------
!$$$      subroutine chk_umn(cnk,umnk,qbz,
!$$$     i                  iko_ix,iko_fx,iko_i,iko_f,
!$$$     d                  nwf,nqbz,nband,nlmto)
!$$$      use m_readqg
!$$$      use m_readeigen
!$$$
!$$$      implicit real*8(a-h,o-z)
!$$$
!$$$      complex(8) :: cnk(iko_ix:iko_fx,nwf,nqbz),
!$$$     &              dnk(iko_ix:iko_fx,nwf,nqbz),
!$$$     &              umnk(nwf,nwf,nqbz)
!$$$      real(8) :: qbz(3,nqbz)
!$$$      integer(4) :: iko_i(nqbz),iko_f(nqbz)
!$$$
!$$$      dnk = (0d0,0d0)
!$$$      do iq = 1,nqbz
!$$$         do imp = iko_i(iq),iko_f(iq)
!$$$         do in = 1,nwf
!$$$            do im = 1,nwf
!$$$               dnk(imp,in,iq) = dnk(imp,in,iq)
!$$$     &                 + umnk(im,in,iq) * cnk(imp,im,iq)
!$$$            enddo ! im
!$$$         enddo ! in
!$$$         enddo ! imp
!$$$      enddo ! iq
!$$$
!$$$      call chk_cnkweight(qbz,iko_ix,iko_fx,dnk,
!$$$     &     nqbz,nwf,nband,nlmto)
!$$$
!$$$      end
