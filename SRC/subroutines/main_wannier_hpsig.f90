subroutine hpsig_MPI()
  !-------------------------------------------------
  ! from huumat.f
  ! Calculate <psi|g>
  ! psi : Kohn Sham orbital
  ! g   : c1*phi + c2*phidot (MT)
  !       gaussian           (otherwise)

  ! Approximation
  ! in the MTs except for the central MT, phi and phidot are
  ! replaced with plane-waves

  ! Takashi Miyake, Mar 2008, parallelized
  ! Takashi Miyake, Jun 2004
  !-------------------------------------------------
  use mpi
  use m_readqg,only: readngmx,ngcmx,readqg0,readqg
  use m_hamindex,only:   readhamindex,ngrp
  use m_readeigen,only:init_readeigen,init_readeigen2,readcphif,readgeigf
  use m_read_bzdata,only: read_bzdata,qbz,nqbz,nqibz,nqbz,qbas=>qlat,nq0i,q0i
  use m_genallcf_v3, ncore2=>ncore,nrxx=>nrx !, nprecb,mrecb,mrece,ndimat,nqbzt,nband,mrecg
  use m_keyvalue,only: getkeyvalue
  use m_read_Worb,only: s_read_Worb, s_cal_Worb, &
       nwf, nclass_mlwf, cbas_mlwf, nbasclass_mlwf, &
       classname_mlwf, iclassin, &
       iphi, iphidot, nphi, nphix
!  use m_readhbe,only: readhbe, nprecb,mrecb,mrece,ndimat,nqbzt,nband,mrecg
  ! RS: MPI module
  !      use rsmpi
  !      use rsmpi_rotkindex
!  use m_scg,only: rotcg

  implicit none
  character(8):: xt
  integer :: iclass2
  real(8):: q(3),  qgbin(3),qx(3)
  integer(4),allocatable :: ngvecpB(:,:,:),ngveccB(:,:), ngvecpf1(:,:), ngvecpf2(:,:), &
       nx(:,:),nblocha(:),ifppb(:) !ongveccBr(:,:,:)
  real(8),allocatable :: ppbrd (:,:,:,:,:,:,:),cg(:,:,:),symope(:,:), &
       phij(:),psij(:),rprodx(:,:),rphiphi(:)!,q0i(:,:)
  complex(8),parameter:: img=(0d0,1d0)
  ! nw,incwf,natom,natom,ipos,igrp,
  !     & iinvg,nspin,nl,nn,nnv,nnc,
  !     o                   inindx,inindxv,inindxc,iiclass,             !l,n, dimensions
  !     d                   ndima,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, !l,n,  dimensions
  !     o                   izdummy,
  !     o   iil,iin,iim,   iilnm, i_mnl, ! l,n,m for Phi ! w(i_mnl)=> mnl(ic) for all electron
  !     o   iilv,iinv,iimv,iilnmv,i_mnlv,! l,n,m for Phi
  !     o   iilc,iinc,iimc,iilnmc,i_mnlc,! l,n,m for Phi
  !     o   iecore,ikonf,iicore,incore,nctot,             !core
  !     o   imagw_dummy,niw,idummy,
  integer(4) &
       nw_input, &
       !     &   ifhbe,
       !     &   nprecb,mrecb,mrece,ndimat,nqbzt,nband,
       i,nq0ix,ngrpmx,mxx,nqbze,nqibze,ini,ix,ngrpx &
       ,mdimx,nbloch,nblochpmx,ifvcfpout,ndummy1,ndummy2,ifcphi,is,nwp, &
       ifepscond,nxx,ifvxcpout,ifgb0vec &
       ,nw0,iw,nwhis,ifinin,nw2,iw0,ifwwk,noccxv,noccx &
       ,ifemesh,nprecx,mrecl,ifwd,ifrcwi,ifrcw,nspinmx,ifianf,ibas &
       ,ibas1,irot,iq,ngb,iqixc2,ifepsdatnolfc,ifepsdat,ngbin,igc0 &
       ,kx,isf,kqxx,kp,job,nbnbx,nhwtot,noccxvx,nwmax &
       ,ihis,jhwtot,ik,ibib,ib1,ib2,ichkhis,ihww,j,imode &
       ,ngpmx

  real(8):: dum1,dum2,dum3,wqtsum,epsrng,dnorm,dini, &
       dwry,dwh,omg_c,omg2,xxx
  integer(4)::nwin, incwfin,  verbose
  real(8)::efin
  integer(4):: ifphi, &
       nbas,nradmx,ncoremx,nrx,ic,icx,isp,l,n,irad,ifoc, &
       ldim2,ixx,ngp1,ngp2,nq0it
  real(8):: qq(3),quu(3), deltaq(3),q1x(3),q2x(3)
  real(8),parameter::  pi =     3.14159265358979323846d0
  real(8),parameter::  fpi =    4d0*pi

  !      logical:: test_qbzreg
  logical:: qbzreg
  !-------------------------------------------------------------------------
  integer(4),allocatable:: ncindx(:,:), &
       lcindx(:,:), &
       nrad(:), &
       nindx_r(:,:), &
       lindx_r(:,:), &
       nc_max(:,:), &
       m_indx(:),n_indx(:),l_indx(:),ibas_indx(:), nrofi(:)
  real(8),allocatable:: phitoto(:,:,:,:,:), aa(:),rr(:,:) &
       ,phitotr(:,:,:,:,:), &
       bb(:),zz(:),rmax(:),cy(:),yl(:)
  complex(8),allocatable:: geig1(:,:),cphi1(:,:),uum(:,:,:), ppovl(:,:)
  complex(8):: ppj,phaseatom
  real(8)   :: q1(3),q2(3),dq(3),absqg2,absdq,r2s,absqg
  integer(4):: j1,j2,j1max,j2max,j1min,j2min,ispin &
       ,l1,l2,lm1,lm2,ibas2,lm3,ig1,ig2,ir,ia1,ma,ia2,m2,l3,m1,lxx &
       ,ico,lxd,lx,ll
  real(8):: ylk

  ! m
  integer(4) :: ixc,idummy,idummy2,i1,i2,i3,nbbloop, &
       ifpsig,ifmloc,ret, &
       ifbb,nbb,iko_ixs(2),iko_fxs(2),noxs(2), &
       iqibz,iqbz,ibb,itmp,itmp2,iti,itf, &
       nqibz2,nqbz2,iqb,ibb2,iqtmp,ibbtmp, &
       ia,iwf!,ifile_handle!,if101
  integer(4),allocatable:: ikidx(:),ikbidx(:,:)
  real(8),allocatable :: bbv(:,:),r0g(:,:),c1(:,:,:),c2(:,:,:), &
       phig(:,:,:,:),wphi(:,:)
  real(8) :: pgnorm,wgt,ndg(3),sij,wphis
  complex(8),allocatable :: psig(:,:,:),qgg(:,:,:)
  logical :: ghead,tailt,debug=.true.

  ! m, MPI
!  include 'mpif.h'
  integer(4):: istatus(MPI_STATUS_SIZE)
  integer(4):: ierr,iftmp
  integer(4):: myproc,nproc,nproc1,nproc2,nq_proc,ii,jj,kk
  integer(4),allocatable:: iq_proc(:)
  character(4) charnum4
  integer,allocatable:: ncore(:)
  !-------------------------

  ! m, MPI
  ! initialize MPI
  call mpi_init(ierr)
  call mpi_COMM_RANK( MPI_COMM_WORLD, myproc, ierr )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, nproc, ierr)
  !      call RSMPI_Init()
  !      myproc = myrank_rsmpi
  !      nproc = nproc_rsmpi

  ! m
  ! mode switch. --------------
  !      write(6,*) ' --- Choose modes below ----------------'
  !      write(6,*) '  ????????????????????? '
  !      write(6,*) ' --- Put number above ! -----------------'
  !      call readin5(ixc,idummy,idummy2)
  !      write(6,*) ' ixc=',ixc
  !      if(ixc==0) stop ' --- ixc=0 --- Choose computational mode!'


  !---  readin BZDATA. See gwsrc/rwbzdata.f
  !--------readin data set when you call read_BZDATA ---------------
  !       integer(4)::ngrp,nqbz,nqibz,nqbzw,nteti,ntetf,
  !     &   n_index_qbz
  !       integer(4):: n1,n2,n3
  !       real(8):: qbas(3,3),ginv(3,3),qbasmc(3,3),dq_bzcase2
  !       real(8),allocatable:: qbz(:,:),wbz(:),qibz(:,:)
  !     &    ,wibz(:),qbzw(:,:)
  !       integer(4),allocatable:: idtetf(:,:),ib1bz(:),idteti(:,:)
  !     &    ,nstar(:),irk(:,:),index_qbz(:,:,:)
  !-----------------------------------------------------------------
  call read_BZDATA()

  !$$$c--- Use regular mesh even for bzcase==2
  !$$$      if(bzcase()==2.and.qbzreg()) then
  !$$$         deltaq= qbas(:,1)/n1 + qbas(:,2)/n2 +qbas(:,3)/n3
  !$$$         do i=1,nqbz
  !$$$            qbz(:,i) = qbz(:,i) -deltaq/2d0
  !$$$            if (myproc.eq.0) write(6,"('i qbz=',i3,3f8.3)") i,qbz(:,i)
  !$$$         enddo
  !$$$      endif

  if (myproc == 0) then
     write(6,*)' ======== nqbz qbz  =',nqbz
     !      write(6,*)' ======== nqibz ngrp ngrp2=',nqibz,ngrp,ngrp2
     do i=1,nqbz
        write(6,"(3f10.5)") qbz(1:3,i)
     enddo
     write(6,*)'============================'
     write(6,*)
  endif

  !--- readin GWIN and LMTO, then allocate and set datas.
  !      nwin = -999    !readin condition. Not readin NW file
  incwfin= 0     !readin condition. use ForX0 for core in GWIN
  !      efin =  -999d0 !readin condition. Not readin EFERMI
  call genallcf_v3(incwfin) !in module m_genallcf_v3
  !      if(ngrp/= ngrp2) stop 'ngrp inconsistent: BZDATA and LMTO GWIN_V2'
  !---  These are allocated and setted by genallcf_v3
  !      integer(4)::  natom,natom,nspin,nl,nn,nnv,nnc, ngrp,
  !     o  ndima,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc, nctot,niw,nw
  !      real(8) :: alat,ef, diw,dw,delta,deltaw,esmr
  !      character(120):: symgrp
  !      character(6),allocatable :: clabl(:)
  !      integer(4),allocatable:: iclass(:)
  !     &  ,nindxv(:,:),nindxc(:,:),ncwf(:,:,:) ,
  !     o    invg(:), il(:,:), in(:,:), im(:,:),   ilnm(:),  nlnm(:),
  !     o    ilv(:),inv(:),imv(:),  ilnmv(:), nlnmv(:),
  !     o    ilc(:),inc(:),imc(:),  ilnmc(:), nlnmc(:),
  !     o    nindx(:,:),konf(:,:),icore(:,:),ncore(:),
  !     &    occv(:,:,:),unoccv(:,:,:)
  !     &   ,occc(:,:,:),unoccc(:,:,:),
  !     o    nocc(:,:,:),nunocc(:,:,:)
  !      real(8), allocatable::
  !     o  plat(:,:),pos(:,:),z(:),  ecore(:,:), freq(:), symgg(:,:,:) ! symgg=w(igrp)

!!!! WE ASSUME iclass(iatom)= iatom !!!!!!!!!!!!!!!!!!!!!!!!!
!  if(natom /= natom) stop ' natom /= natom '

  ! --- read dimensions of h,hb
  !      ifhbe      = iopen('hbe.d',1,0,0)
  !      read (ifhbe,*) nprecb,mrecb,mrece,ndimat,nqbzt,nband,mrecg
  !      if(ndima/=ndimat) stop ' hx0fp0: ndima/=ndimat in hbe.d'
  !      if(nqbz /=nqbzt ) stop ' hx0fp0: nqbz /=nqbzt  in hbe.d'
!  call Readhbe()

  ! --- read by rdpp ; Radial integrals ppbrd and plane wave part
  call getsrdpp2(natom,nl,nxx)
  call readngmx('QGpsi',ngpmx)
  if (myproc == 0) write(6,*)' ngpmx=',ngpmx

  ! --- read radial functions PHIVC   (taken from hasfp0)
  if (myproc == 0) write(6,*)' Go to readining phivc'
  open(newunit=ifphi,file='__PHIVC',form='unformatted')
  read(ifphi) nbas, nradmx, ncoremx,nrx
  if( nbas/=natom ) stop ' nbas(PHIVC) /= natom '
  !      deallocate(ncore)
  allocate(  ncindx(ncoremx,nbas), &
       lcindx(ncoremx,nbas), &
       nrad(nbas), &
       nindx_r(1:nradmx,1:nbas), &
       lindx_r(1:nradmx,1:nbas), &
       aa(nbas),bb(nbas),zz(nbas), rr(nrx,nbas), nrofi(nbas) , &
       phitoto(nrx,0:nl-1,nn,nbas,nspin), &
       phitotr(nrx,0:nl-1,nn,nbas,nspin), &
       nc_max(0:nl-1,nbas),ncore(nbas),rmax(nbas) )
  if (myproc == 0) write(6,*)' end of allocation'
  read(ifphi) nrad(1:nbas)
  read(ifphi) nindx_r(1:nradmx,1:nbas),lindx_r(1:nradmx,1:nbas)
  nc_max=0
  do ibas=1,nbas
     ic = ibas
     if (myproc == 0) &
          write(6,*)' --- read PHIVC of ibas nrad=',ibas,nrad(ic)
     read(ifphi) ncore(ic), ncoremx                            !core
     read(ifphi) ncindx(1:ncoremx,ibas),lcindx(1:ncoremx,ibas) !core
     if (myproc == 0) write(6,*)' xxx0'
     read(ifphi) icx,zz(ic),nrofi(ic),aa(ic),bb(ic)
     if (myproc == 0) &
          write(6,*) 'ic icx=',ic,icx,zz(ic),nrofi(ic),aa(ic),bb(ic)
     if(ic/=icx) then
        stop ' h_uu: ic/=icx'
     endif
     if (myproc == 0) write(6,*)' xxx1 ncoremx ncore(ic)=',ncoremx,ncore(ic)
     read(ifphi) rr(1:nrofi(ic),ic)
     if (myproc == 0) write(6,*)' xxx2 ncoremx ncore(ic)=',ncoremx,ncore(ic)
     if (myproc == 0) write(6,*)' xxx2 nspin=',nspin
     rmax(ic) = rr(nrofi(ic),ic)
     do isp = 1, nspin
        if (myproc == 0) write(6,*)'          ---  isp nrad ncore(ic)=',isp, nrad(ic),ncore(ic)
        do ico = 1, ncore(ic) !core
           l =  lcindx(ico,ic)
           n =  ncindx(ico,ic)
           read(ifphi) phitoto(1:nrofi(ic),l,n, ic,isp)   !core orthogonal
           phitotr(1:nrofi(ic),l,n, ic,isp)=     &         ! & core raw= core orthgonal
           phitoto(1:nrofi(ic),l,n, ic,isp)               !
           if(n>nc_max(l,ic)) nc_max(l,ic)=n
           if (myproc == 0) &
                write(6,*)' sss1c=',sum(abs(phitoto(1:nrofi(ic),l,n, ic,isp)))
        enddo
        do irad = 1, nrad(ic)   !valence
           l = lindx_r (irad,ic)
           n = nindx_r (irad,ic) + nc_max(l,ic)
           read(ifphi) phitoto(1:nrofi(ic),l,n, ic,isp) !valence orthogonal
           read(ifphi) phitotr(1:nrofi(ic),l,n, ic,isp) !valence raw
           if (myproc == 0) then
              write(6,*)' sss1=',sum(abs(phitoto(1:nrofi(ic),l,n, ic,isp)))
              write(6,*)' sss2=',sum(abs(phitotr(1:nrofi(ic),l,n, ic,isp)))
           endif
        enddo
     enddo
  enddo
  !--- cg coefficient.  y = cg y y ; y is the real spherical harmonics
  ngrpx=1
  allocate( cg(nl**2,nl**2,(2*nl-1)**2), symope(3,3) )
  symope(1:3,1) = (/1d0,0d0,0d0/)
  symope(1:3,2) = (/0d0,1d0,0d0/)
  symope(1:3,3) = (/0d0,0d0,1d0/)
  cg = 0d0 !for sanity check
  call rotcg(nl-1,symope,ngrpx,cg)
  !! initiallization to get eigenfunctions
  call Readhamindex()
  call init_readeigen()!nband,mrece) !initialization of readEigen
  call init_readeigen2()!mrecb,ndima,mrecg)
  call readngmx('QGpsi',ngpmx)
  allocate( geig1(ngpmx*nspc,nband) )
  if (myproc == 0) write(6,*) 'end of initialization'
  ! --- Readin nlam index
  !  ifoc = iopen('@MNLA_CPHI',1,0,0)
  open(newunit=ifoc,file='@MNLA_CPHI')
  ldim2 = ndima
  read(ifoc,*)
  allocate(m_indx(ldim2),n_indx(ldim2),l_indx(ldim2),ibas_indx(ldim2))
  do ix =1,ldim2
     read(ifoc,*)m_indx(ix),n_indx(ix),l_indx(ix),ibas_indx(ix),ixx
     if(ixx/=ix) stop  'failed to readin @MNLA_CPHI'
  enddo
  ! ---  q near zero
  if (myproc == 0) write(6,*) 'reading QOP'
  !$$$      if101=ifile_handle()
  !$$$      open (if101,file='Q0P')
  !$$$      read (if101,"(i5)") nq0i
  !$$$!      if(.not.exchange) call checkeq(nqibz+nq0i-1, nqnum)
  !$$$      if (myproc.eq.0) write(6,*) ' *** nqibz nq0i_total=', nqibz,nq0i
  !$$$      nq0it = nq0i
  !$$$      allocate( q0i(1:3,1:nq0i) ) !wqt(1:nq0i),
  !$$$!      read (if101,"(d24.16,3x, 3d24.16)" )( wqt(i),q0i(1:3,i),i=1,nq0i)
  !$$$      nq0ix = nq0i
  !$$$      do i=1,nq0i
  !$$$      read (if101,* ) xxx,q0i(1:3,i)
  !$$$      if(xxx==0d0 ) nq0ix = i-1
  !$$$      enddo
  !$$$      nq0i = nq0ix ! New nq0i July 2001
  if (myproc == 0) then
     write(6,*) ' Used k number in Q0P =', nq0i
     write(6,"(i3,2x, 3f14.6)" )(i,q0i(1:3,i),i=1,nq0i)
  endif
!  close(if101)

  ! m
!  ifbb = iopen('BBVEC',1,0,0)
  open(newunit=ifbb,file='BBVEC')
  read(ifbb,*)
  read(ifbb,*)nbb,nqbz2
  if (nqbz /= nqbz2) stop 'readbb: nqbz is wrong!'
  allocate (bbv(3,nbb),ikbidx(nbb,nqbz))
  call  readbb(ifbb,nqbz,nspin,nbb, &
       bbv, &
       ikbidx, &
       iko_ixs,iko_fxs,noxs)
  ! GWinput data
  call s_read_Worb()
  do iclass2=1,nclass_mlwf
     write(*,*)'output:',iclassin(iclass2), nwf &
          ,trim(classname_mlwf(iclass2)),cbas_mlwf(1:nbasclass_mlwf(iclass2) &
          ,iclass2)
  enddo

  call s_cal_Worb()

  allocate (r0g(nphix,nwf), wphi(nphix,nwf))
  r0g = 2d0
  wphi = 1d0

  ! 061004
  if (myproc == 0) then
     call getkeyvalue("GWinput","wan_gauss_head",ghead,default=.false.)
     call getkeyvalue("GWinput","wan_truncate",tailt,default=.false.)
  endif ! myproc
  call MPI_Bcast(ghead,1,MPI_LOGICAL,0, &
       MPI_COMM_WORLD,ierr)
  !      call RSMPI_Check("MPI_Bcast(ghead)",ierr)
  call MPI_Bcast(tailt,1,MPI_LOGICAL,0, &
       MPI_COMM_WORLD,ierr)
  !      call RSMPI_Check("MPI_Bcast(tailt)",ierr)

  ! normalize wphi
  do i = 1,nwf
     wphis = dsqrt(sum(wphi(1:nphi(i),i)**2))
     wphi(1:nphi(i),i) = wphi(1:nphi(i),i)/wphis
  enddo
  ! matching condition at MT radius
  allocate(c1(nphix,nwf,nspin),c2(nphix,nwf,nspin))
  if ( .NOT. ghead) &
       call getc1c2(phitotr,r0g, &
       iphi,iphidot,n_indx,l_indx,m_indx,ibas_indx, &
       rr,nrofi,rmax,nphi,nc_max, &
       nwf,ndima,nrx,nl,nn,nbas,nspin,nphix, &  !ndima =>ndima
       c1,c2)

  if (myproc == 0) then
     write(*,*)'Gaussian Head =',ghead
     write(*,*)'# gaussian =',nwf
     write(*,*)'        No.   n    l    m  ibas'
     do i = 1,nwf
        write(*,*)'iwf,',i,nphi(i)
        do j = 1,nphi(i)
           ix = iphi(j,i)
           write(*,"('phi   ',5i5)") ix,n_indx(ix),l_indx(ix),m_indx(ix),ibas_indx(ix)
           ix = iphidot(j,i)
           write(*,"('phidot',5i5)") ix,n_indx(ix),l_indx(ix),m_indx(ix),ibas_indx(ix)
        enddo
        write(*,*)
     enddo
  endif ! myproc

  if (myproc == 0) then
     open(newunit=ifpsig,file='PSIGU.'//charnum4(0), form='unformatted')
     write(ifpsig)nqbz,iko_ixs(1),iko_fxs(1),nwf
     close(ifpsig)
     open(newunit=ifpsig,file='PSIGD.'//charnum4(0), form='unformatted')
     write(ifpsig)nqbz,iko_ixs(2),iko_fxs(2),nwf
     close(ifpsig)
  endif                     ! myproc

  !======================================================================
  ! --- Set q1(j1range)
  !======================================================================
  ! Note that the true q when we generate eigenfunctions are q1x and q2x.
  ! q1-q1x should be a G vector.

  ! --- I inserted checkagree to make sure that q1=q1x and q2=q2x ...

  ! m
  j1min = iko_ixs(1)
  j1max = iko_fxs(1)
  if (nspin == 2) then
     if (iko_ixs(2) < j1min) j1min = iko_ixs(2)
     if (iko_fxs(2) > j1max) j1max = iko_fxs(2)
  endif
  j2min = 1
  j2max = nwf
  !======================================================================

  allocate( psig(j1min:j1max,nwf,nspin) )
  ! m, MPI
  nq_proc = nqbz / nproc
  nproc1 = nqbz - nq_proc*nproc
  nproc2 = nproc - nproc1
  allocate(iq_proc(nq_proc+1))
  iq_proc = 0
  !iftmp = 100 + myproc
  open(newunit=iftmp,file='myproc'//xt(myproc))
  write(iftmp,*)'*** my proc,i,q(i)'
  kk = 0
  do jj = 1,nproc
     do ii = 1,nq_proc+1
        if (jj > nproc1 .AND. ii > nq_proc) cycle
        kk = kk + 1
        if (myproc == jj-1) then
           iq_proc(ii) = kk
           write(iftmp,*)myproc,ii,iq_proc(ii)
        endif
     enddo ! ii
  enddo ! jj
  do 1070 ii = 1,nq_proc+1
     if(debug) write(6,*)'do1070start ii=',ii
     iqbz = iq_proc(ii)
     if (iqbz == 0) cycle
     write(*,*)    'iqbz =',iqbz, 'out of',nqbz
     write(iftmp,*)'iqbz =',iqbz, 'out of',nqbz
     q1(:) = qbz(:,iqbz)
     ! --- q1x
     call readqg0('QGpsi',q1,  q1x, ngp1)     !      call checkagree(q1,q1x,' q1 ne q1x')
     if (myproc == 0) write(6,"(' q1 q1x=',3f9.4,3x,3f9.4)") q1,q1x
     ! ... lxx and allocations
     lxx=2*(nl-1)
     ! --- < q+G | g >
     allocate( ngvecpf1(3,ngp1),qgg(ngp1,nphix,nwf) )
     if (tailt) then
        qgg = 0d0
     else
        call readqg('QGpsi',q1, q1x, ngp1, ngvecpf1)        !      call checkagree(q1,q1x,' q1 ne q1x xxx2')
        call qggmat(q1,alat,plat,qbas, &
             iphi,r0g,pos, &
             rr,nrofi,rmax,aa,bb, &
             n_indx,l_indx,m_indx,ibas_indx, &
             ngvecpf1,nphi, &
             lxx,ngp1,nwf,ndima,nbas,nrx,nphix, &
             qgg)
        if(debug) write(*,*)'qggmat done'
     endif
     ! --- < phi | g >
     allocate ( phig(ndima,nphix,nwf,nspin) )
     call phigmat(iphi,iphidot,c1,c2,phitoto,phitotr, &
          rr,nrofi,rmax,aa,bb,r0g,ghead, &
          n_indx,l_indx,m_indx,ibas_indx,nphi,nc_max, &
          nwf,ndima,nn,nl,nbas,nrx,nspin,nphix, &
          phig)
     if(debug) write(*,*)'phigmat done'
     ! --- Calcuate <psi{q1x j1} | g_{iwf}>
     psig = 0d0
     do 1050 ispin=1,nspin
        ! ... MT part
        !r   ldim2 = ndima
        !r   n_indx   (1;ldim2) : n index (phi=1 phidto=2 localorbital=3)
        !r   l_indx   (1:ldim2) : l index
        !r   ibas_indx(1:ldim2) : ibas index.
        allocate( cphi1 (ndima*nspc,nband) )
        cphi1= readcphif(q1,ispin) !    call readcphi(q1, ndima, ispin, quu, cphi1)
        !   call checkagree(q1,q1x,' q1 ne quu')
        do 1020 ia = 1,ndima
           do j1= iko_ixs(ispin),iko_fxs(ispin)
              do iwf = 1,nwf
                 do j   = 1,nphi(iwf)
                    psig(j1,iwf,ispin) =  psig(j1,iwf,ispin) &
                         + wphi(j,iwf)*dconjg(cphi1(ia,j1))*phig(ia,j,iwf,ispin)
                 enddo
              enddo
           enddo
1020    enddo
        if(debug)  write(*,*)'MT part done'
        ! ... Interstitial Plane Wave part
        geig1=readgeigf(q1,ispin) !call readgeig(q1, ngpmx, ispin, quu, geig1)
        !       call checkagree(q1,quu,' q1 ne quu eig')
        !       do j1=j1min,j1max
        do j1  = iko_ixs(ispin),iko_fxs(ispin)
           do iwf = 1,nwf
              do j   = 1,nphi(iwf)
                 psig(j1,iwf,ispin) =  psig(j1,iwf,ispin) &
                      +  wphi(j,iwf)*sum( dconjg(geig1(1:ngp1,j1)) * qgg(1:ngp1,j,iwf))
              enddo
           enddo
        enddo
        deallocate(cphi1)
1050 enddo
     if(debug) write(6,*)'end of do1050'
     !! write result
     do ispin = 1,nspin
        if(ispin==1) open(newunit=ifpsig,file='PSIGU.'//charnum4(iqbz),form='unformatted')
        if(ispin==2) open(newunit=ifpsig,file='PSIGD.'//charnum4(iqbz),form='unformatted')
        iti = iko_ixs(ispin)
        itf = iko_fxs(ispin)
        if(debug) write(6,*)'wwwwwww', iqbz,ispin,iti,itf
        write(ifpsig) iqbz
        write(ifpsig) ((psig(j1,j2,ispin),j1=iti,itf),j2=1,nwf)
        close(ifpsig)
     enddo
     ! eck write
     !      if (mod(iqbz,10).eq.1) then
     !         do ispin = 1,nspin
     !            iti = iko_ixs(ispin)
     !            itf = iko_fxs(ispin)
     !            do iwf = 1,nwf
     !               pgnorm = 0d0
     !               do j1 = iti,itf
     !                  pgnorm = pgnorm + abs(psig(j1,iwf,ispin))**2
     !c                  write(98,"(3i5,2f12.6)")iqbz,iwf,j1,
     !c     &                  abs(psig(j1,iwf,ispin))**2
     !               enddo
     !               write(*,*)ispin,iwf,pgnorm
     !            enddo
     !         enddo
     !      endif
     !      do ispin = 1,nspin
     !         iti = iko_ixs(ispin)
     !         itf = iko_fxs(ispin)
     !         do j1 = iti,itf
     !         do j2 = iti,itf
     !            sij = 0d0
     !            do iwf = 1,nwf
     !               sij = sij + dconjg(psig(j1,iwf,ispin))*psig(j2,iwf,ispin)
     !            enddo
     !            if (j1.eq.j2) then
     !               write(98,"(2i5,f12.6)")iqbz,j1,sij
     !            else
     !               write(97,"(3i5,f12.6)")iqbz,j1,j2,sij
     !            endif
     !         enddo
     !         enddo
     !      enddo
     deallocate(ngvecpf1, qgg, phig)
     ! end of qbz-loop
     if(debug) write(6,*)'do1070end ii=',ii
1070 enddo
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  ! finalize MPI
  call mpi_finalize(ierr)
  if (myproc == 0) write(6,*) ' ====== end ========================================'
  if (myproc == 0) write(*,*) 'hpsig: ok'
end subroutine hpsig_MPI

! !====================================================================
! subroutine checkagree(a,b,char)
!   real(8):: a(3),b(3)
!   character*(*) :: char
!   if(sum(abs(a-b))>1d-6) then
!      write(6,*)' Error in checkagree:',char
!      stop ' Error in checkagree:'
!   endif
! end subroutine checkagree
! !-----------------------------------------------------------------------
! subroutine  readbb(ifbb,nqbz,nspin,nbb, &
!      bbv, &
!      ikbidx, &
!      iko_ixs,iko_fxs,noxs)
!   implicit integer (i-n)
!   implicit real*8(a-h,o-z)
!   parameter (eps = 1d-4)

!   real(8) :: u(3),bbv(3,nbb)
!   integer :: iopen, &
!        iko_ixs(2),iko_fxs(2),noxs(2)
!   integer:: ikbidx(nbb,nqbz)
!   !      integer(4),allocatable:: ikidx(:),ikbidx(:,:)

!   !      ifbb = iopen('BBVEC',1,0,0)
!   !      read(ifbb,*)
!   !      read(ifbb,*)nbb,nqibz2,nqbz2
!   !      if (nqibz.ne.nqibz2) stop 'readbb: nqibz is wrong!'
!   !      if (nqbz.ne.nqbz2) stop 'readbb: nqbz is wrong!'

!   !      allocate (ikidx(nqibz),ikbidx(nbb,nqibz))

!   do i = 1,nbb
!      read(ifbb,*)bbv(1,i),bbv(2,i),bbv(3,i),dummy4
!   enddo
!   do iq = 1,nqbz
!      read(ifbb,*)itmp,u(1:3)
!      do ib = 1,nbb
!         read(ifbb,*)itmp,itmp2,ikbidx(ib,iq),u(1:3)
!      enddo
!   enddo
!   read(ifbb,*)
!   read(ifbb,*)nspin2
!   if (nspin /= nspin2) stop 'nspin is wrong!'
!   do is = 1,nspin
!      read(ifbb,*)iko_ixs(is),iko_fxs(is),noxs(is)
!   enddo


!   return
! end subroutine readbb
! !-----------------------------------------------------------------------
subroutine  getc1c2(phitotr,r0g, &
     iphi,iphidot,n_indx,l_indx,m_indx,ibas_indx, &
     rr,nrofi,rmax,nphi,nc_max, &
     nwf,ndima,nrx,nl,nn,nbas,nspin,nphix, &
     c1,c2)
  implicit integer (i-n)
  implicit real*8(a-h,o-z)
  real(8) :: phitotr(nrx,0:nl-1,nn,nbas,nspin),r0g(nphix,nwf), &
       rr(nrx,nbas),rmax(nbas), &
       c1(nphix,nwf,nspin),c2(nphix,nwf,nspin)
  integer(4) :: iphi(nphix,nwf),iphidot(nphix,nwf), &
       n_indx(ndima),l_indx(ndima),m_indx(ndima), &
       ibas_indx(ndima),nrofi(nbas),nphi(nwf), &
       nc_max(0:nl-1,nbas)


  pi = 4d0*datan(1d0)

  c1 = 0d0
  c2 = 0d0

  do is  = 1,nspin
     do iwf = 1,nwf
        do j   = 1,nphi(iwf)

           ip = iphi(j,iwf)
           ibp = ibas_indx(ip)
           ilp = l_indx(ip)
           imp = m_indx(ip)
           inp = n_indx(ip) + nc_max(ilp,ibp)

           id = iphidot(j,iwf)
           ibd = ibas_indx(id)
           ild = l_indx(id)
           imd = m_indx(id)
           ind = n_indx(id) + nc_max(ild,ibd)

           if (ibp /= ibd) stop 'getc1c2: ibp /= ibd'
           if (ilp /= ild) stop 'getc1c2: ilp /= ild'
           if (imp /= imd) stop 'getc1c2: imp /= imd'

           ! g(r) = A exp(-(r/r0)**2)
           ! A = (128 / (pi*r0))**(1/4)
           Ag = (128d0 / (pi*r0g(j,iwf)**6))**0.25d0
           rmt = rr(nrofi(ibp),ibp)
           g  = Ag * dexp(-(rmt/r0g(j,iwf))**2)
           dg = g * (-2d0*rmt/r0g(j,iwf)**2)

           ! phi
           rmt  = rr(nrofi(ibp),ibp)
           rmt1 = rr(nrofi(ibp)-1,ibp)
           p  =  phitotr(nrofi(ibp),ilp,inp,ibp,is)/rmt
           dp = (phitotr(nrofi(ibp),  ilp,inp,ibp,is)/rmt &
                - phitotr(nrofi(ibp)-1,ilp,inp,ibp,is)/rmt1) &
                / (rmt - rmt1)

           ! phidot
           d  =  phitotr(nrofi(ibd),ild,ind,ibd,is)/rmt
           dd = (phitotr(nrofi(ibd),  ild,ind,ibd,is)/rmt &
                - phitotr(nrofi(ibd)-1,ild,ind,ibd,is)/rmt1) &
                / (rmt - rmt1)

           ! c1 *  p + c2 *  d =  g
           ! c1 * dp + c2 * dd = dg
           detinv = 1d0 / (p * dd - dp * d)
           c1(j,iwf,is) = detinv * (dd  * g - d * dg)
           c2(j,iwf,is) = detinv * (-dp * g + p * dg)

           ! eck
           gtmp = c1(j,iwf,is) *  p + c2(j,iwf,is) *  d
           dgtmp = c1(j,iwf,is) * dp + c2(j,iwf,is) * dd
           !         write(*,*)'iwf,j =',iwf,j
           write(*,*)'c1,c2             ',is,c1(j,iwf,is),c2(j,iwf,is)
           !         write(*,*)'phi,phidot        ',p,d
           !         write(*,*)'dphi/dr,dphidot/dr',dp,dd
           !         write(*,*)'g,gtmp            ',g,gtmp
           !         write(*,*)'dg,dgtmp          ',dg,dgtmp

        enddo
     enddo
  enddo

  return
end subroutine getc1c2
!--------------------------------------------------------------
subroutine qggmat(q,alat,plat,qlat, &
     iphi,r0g,bas, &
     rr,nrofi,rmax,aa,bb, &
     n_indx,l_indx,m_indx,ibas_indx, &
     ngvec,nphi, &
     lxx,ng,nwf,ndima,nbas,nrx,nphix, &
     qgg)
  ! < q+G | g> where q+G denotes IPW, zero within MT sphere   ! g(r) = Ag exp(-(r/r0)**2)
  implicit real*8 (a-h,o-z)
  implicit integer (i-n)
  parameter (nmt = 4)
  parameter (eps = 1d-5)
  complex(8) :: qgg(ng,nphix,nwf),eiqgr
  complex(8) :: img =(0d0,1d0)
  real(8) :: q(3),plat(3,3),qlat(3,3),r0g(nphix,nwf),bas(3,nbas), &
       rr(nrx,nbas),rmax(nbas),aa(nbas),bb(nbas), &
       g(3),qg(3),ag(nphix,nwf)
  integer(4) :: iphi(nphix,nwf),nrofi(nbas),ntail(nphix,nwf), &
       n_indx(ndima),l_indx(ndima),m_indx(ndima),ibas_indx(ndima), &
       ngvec(3,ng),nphi(nwf)
  integer(4) :: ng,nwf,ndima,nbas,nrx
  real(8) :: absqg,tripl,pi,rbas(3,nphix,nwf), alat,tpiba,tpibaqlat(3,3),voltot
  real(8),allocatable :: cy(:),yl(:),rtail(:,:,:),rg(:,:,:), phij(:),psij(:),rjl(:)
  !-----------------------------------------------------
  !r True q is given by
  !r    True_q(1:3)     = 2*pi/alat * q(1:3)
  !r True G is given by
  !r    True_G(1:3,igp) = 2*pi/alat * matmul(qlat * ngvec(1:3,ig)) ,ig=1,ng
  !-----------------------------------------------------
  write(6,*)' qggmat:'
  pi        = 4d0*datan(1d0)
  voltot    = abs(alat**3*tripl(plat,plat(1,2),plat(1,3)))
  tpibaqlat =  2*pi/alat *qlat
  tpiba     =  2*pi/alat

  qgg = 0d0
  ntailx = 0
  do iwf = 1,nwf
     do j   = 1,nphi(iwf)
        ilmto = iphi(j,iwf)
        ibas  = ibas_indx(ilmto)
        rbas(1:3,j,iwf) =  bas(1:3,ibas)*alat
        ag(j,iwf) = (128d0 / (pi*r0g(j,iwf)**6))**0.25d0
        ntail(j,iwf) = int(log(dble(nmt)) / aa(ibas)) + nrofi(ibas) + 1
        if (ntailx < ntail(j,iwf)) ntailx = ntail(j,iwf)
     enddo
  enddo
  allocate(cy((lxx+1)**2),yl((lxx+1)**2))
  allocate(phij(0:lxx),psij(0:lxx),rjl(ntailx))
  write(6,*)'rrrrrrrr',ntailx,nphix,nwf
  allocate(rtail(ntailx,nphix,nwf),rg(ntailx,nphix,nwf))
  do iwf = 1,nwf
     do j   = 1,nphi(iwf)
        ilmto = iphi(j,iwf)
        ibas =  ibas_indx(ilmto)
        do ir = 1,ntail(j,iwf)
           rtmp = bb(ibas)*exp(aa(ibas)*dble(ir-1))
           rtail(ir,j,iwf) = rtmp
           ! rg = r * gaussian
           rg(ir,j,iwf) = rtmp*ag(j,iwf)*dexp(-(rtmp/r0g(j,iwf))**2)
        enddo
     enddo
  enddo


  do ig  = 1, ng
     qg(1:3) = tpiba * (q(1:3)+ matmul(qlat, ngvec(1:3,ig)))
     absqg2  = sum(qg**2)
     absqg   = sqrt(absqg2)

     if (absqg > eps) then

        do iwf = 1,nwf
           do j   = 1,nphi(iwf)

              ilmto = iphi(j,iwf)
              ibas  = ibas_indx(ilmto)
              il    = l_indx(ilmto)
              im    = m_indx(ilmto)
              eiqgr = exp( -img* sum(qg(1:3)*rbas(1:3,j,iwf)) )
              a     = aa(ibas)
              b     = bb(ibas)

              ! Ylm(qg)
              call sylmnc(cy,lxx)
              call sylm(qg/absqg,yl,lxx,r2s) !spherical factor Y(dq)
              ilm   = il**2 + il+1 + im
              ylmqg = cy(ilm)*yl(ilm)

              ! r \times j_l(|qg|r)  !bessel function
              do ir =1,ntail(j,iwf)
                 rtmp = rtail(ir,j,iwf)
                 call bessl(absqg2*rtmp**2,lxx,phij,psij)
                 !  phij(lx) \approx 1/(2l+1)!! for small absqg*rr(ir,ibas).
                 if(rtmp==0d0) then
                    rjl(ir)=0d0
                 else
                    rjl(ir) = rtmp* phij(il)* (absqg*rtmp)**il
                 endif
              enddo

              ! I[rmt:\infty] dr r^2 jl(|qg|r) gaussian(r)
              nr = nrofi(ibas)
              call gintxx(rg(1:nr,j,iwf),rjl(1:nr),a,b,nr,sum1)
              nr = ntail(j,iwf)
              call gintxx(rg(1:nr,j,iwf),rjl(1:nr),a,b,nr,sum2)
              r2jlg = sum2 - sum1

              ! < qg | gaussian >
              qgg(ig,j,iwf) = 4d0*pi*(-img)**il * ylmqg * r2jlg*eiqgr

           enddo
           ! end of iwf-loop
        enddo

        ! if (absqg=0)
     else

        do iwf = 1,nwf
           do j   = 1,nphi(iwf)
              ilmto = iphi(j,iwf)
              ibas  = ibas_indx(ilmto)
              il    = l_indx(ilmto)
              im    = m_indx(ilmto)
              a     = aa(ibas)
              b     = bb(ibas)

              if (il /= 0) then
                 qgg(ig,j,iwf) = 0d0
              else
                 nr = nrofi(ibas)
                 call gintxx(rg(1:nr,j,iwf),rtail(1:nr,j,iwf),a,b,nr,sum1)
                 nr = ntail(j,iwf)
                 call gintxx(rg(1:nr,j,iwf),rtail(1:nr,j,iwf),a,b,nr,sum2)
                 qgg(ig,j,iwf) = (sum2-sum1) * dsqrt(4d0*pi)
              endif
           enddo
        enddo

        ! end of if (absqg=0)
     endif

     ! end of ig-loop
  enddo
  deallocate(cy,yl,rtail,rg,phij,psij,rjl)
  !      write(6,*)' qggmat: done '
end subroutine qggmat

!--------------------------------------------------------------
subroutine phigmat(iphi,iphidot,c1,c2,phitoto,phitotr, &
     rr,nrofi,rmax,aa,bb,r0g,ghead, &
     n_indx,l_indx,m_indx,ibas_indx,nphi,nc_max, &
     nwf,ndima,nn,nl,nbas,nrx,nspin,nphix, &
     phig)
  ! < phi | g>
  ! g = c1*phi + c2*phidot, which is conneted to the gaussian tail
  !     at r=rmt
  implicit integer (i-n)
  implicit real*8 (a-h,o-z)

  real(8) :: phig(ndima,nphix,nwf,nspin), &
       phitoto(nrx,0:nl-1,nn,nbas,nspin), &
       phitotr(nrx,0:nl-1,nn,nbas,nspin)
  real(8) :: rr(nrx,nbas),rmax(nbas),aa(nbas),bb(nbas), &
       c1(nphix,nwf,nspin),c2(nphix,nwf,nspin), &
       rphi(nrx),rhead(nrx),r0g(nphix,nwf)
  integer(4) :: iphi(nphix,nwf),iphidot(nphix,nwf),nrofi(nbas), &
       n_indx(ndima),l_indx(ndima),m_indx(ndima),ibas_indx(ndima), &
       nphi(nwf),nc_max(0:nl-1,nbas)
  logical :: ghead

  pi   = 4d0*datan(1d0)
  phig = 0d0

  do is  = 1,nspin
     do iwf = 1,nwf
        do j   = 1,nphi(iwf)

           if (ghead) then
              inlmp = iphi(j,iwf)
              ibasp = ibas_indx(inlmp)
              ilp   = l_indx(inlmp)
              imp   = m_indx(inlmp)
              inp   = n_indx(inlmp) + nc_max(ilp,ibasp)
              ! g(r) = A exp(-(r/r0)**2)
              ! A = (128 / (pi*r0))**(1/4)
              Ag = (128d0 / (pi*r0g(j,iwf)**6))**0.25d0
              rhead(1) = 0d0
              do ir = 2,nrofi(ibasp)
                 rmt = rr(ir,ibasp)
                 g  = Ag * dexp(-(rmt/r0g(j,iwf))**2)
                 rhead(ir) = rmt * g
              enddo

           else ! ghead
              ! r * phi
              inlmp = iphi(j,iwf)
              ibasp = ibas_indx(inlmp)
              ilp   = l_indx(inlmp)
              imp   = m_indx(inlmp)
              inp   = n_indx(inlmp) + nc_max(ilp,ibasp)

              ! r * phidot
              inlmd = iphidot(j,iwf)
              ibasd = ibas_indx(inlmd)
              ild    = l_indx(inlmd)
              imd    = m_indx(inlmd)
              ind    = n_indx(inlmd) + nc_max(ild,ibasd)

              if (ibasp /= ibasd) stop 'phigmat: bas(phi) /= bas(phidot)'
              if (ilp /= ild) stop 'phigmat: l(phi) /= l(phidot)'
              if (imp /= imd) stop 'phigmat: m(phi) /= m(phidot)'

              rhead(1) = 0d0
              do ir = 2,nrofi(ibasp)
                 rhead(ir) = c1(j,iwf,is) * phitotr(ir,ilp,inp,ibasp,is) &
                      + c2(j,iwf,is) * phitotr(ir,ild,ind,ibasd,is)
              enddo

           endif ! ghead

           do inlm = 1,ndima
              ibas = ibas_indx(inlm)
              il   = l_indx(inlm)
              im   = m_indx(inlm)
              in   = n_indx(inlm) + nc_max(il,ibas)
              if (ibas /= ibasp) cycle
              if (il /= ilp) cycle
              if (im /= imp) cycle
              call gintxx(phitoto(1:nrofi(ibas),il,in,ibas,is), &
                   rhead,aa(ibas),bb(ibas),nrofi(ibas),sum)
              phig(inlm,j,iwf,is) = sum
              !            write(*,*)inlm,iwf,sum
           enddo
        enddo
        ! end of iwf-loop
     enddo
     ! end of is-loop
  enddo
end subroutine phigmat
subroutine  readbb(ifbb,nqbz,nspin,nbb, bbv, ikbidx, iko_ixs,iko_fxs,noxs)
   implicit integer (i-n)
   implicit real*8(a-h,o-z)
   parameter (eps = 1d-4)
   real(8) :: u(3),bbv(3,nbb)
   integer :: iko_ixs(2),iko_fxs(2),noxs(2)
   integer:: ikbidx(nbb,nqbz)
   do i = 1,nbb
      read(ifbb,*)bbv(1,i),bbv(2,i),bbv(3,i),dummy4
   enddo
   do iq = 1,nqbz
      read(ifbb,*)itmp,u(1:3)
      do ib = 1,nbb
         read(ifbb,*)itmp,itmp2,ikbidx(ib,iq),u(1:3)
      enddo
   enddo
   read(ifbb,*)
   read(ifbb,*)nspin2
   if (nspin /= nspin2) call rx('nspin is wrong!')
   do is = 1,nspin
      read(ifbb,*)iko_ixs(is),iko_fxs(is),noxs(is)
   enddo
end subroutine readbb
