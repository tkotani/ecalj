module util_hwmatK
  use mpi
  use m_maxloc0,only: sortvec2
  contains
subroutine wwmat (is,nw_i,nw,nwf, &
     rws1,rws2,irws1,irws2,nrws1,nrws2,nrws, &
     alat,rcut1,rcut2, &
     freq, &
     rw_w,cw_w,rv_w,cv_w, &
     lcrpa, lomega0)
  implicit real*8(a-h,o-z)
  implicit integer (i-n)
  integer :: irws1(nrws1),irws2(nrws2)
  real(8) :: rws1(3,nrws1),rws2(3,nrws2)
  real(8) :: rydberg,hartree
  real(8) :: freq(1:nw),freq2(1:nw)   
  real(8) :: rw_w(nwf,nwf,nwf,nwf,nrws,1:nw), cw_w(nwf,nwf,nwf,nwf,nrws,1:nw)  
  integer:: iwf1, iwf2, iwf3, iwf4, ifreq2
  real(8) :: rv_w(nwf,nwf,nwf,nwf,nrws),cv_w(nwf,nwf,nwf,nwf,nrws)
  logical:: lcrpa, lomega0
  hartree=2d0*rydberg()
  freq2 = hartree*freq
  rw_w  = hartree*rw_w
  cw_w  = hartree*cw_w

  write(*,*)'Writing Screened Couloumb interaction (W-v) : Real'
  if ((is==1) .AND. (lcrpa .eqv. .FALSE. )) then
     open(newunit=ifscr,file="Screening_W-v.UP")
  else if ((is==2) .AND. (lcrpa .eqv. .FALSE. )) then
     open(newunit=ifscr,file="Screening_W-v.DN")
  else if ((is==1) .AND. (lcrpa .eqv. .TRUE. )) then
     open(newunit=ifscr,file="Screening_W-v_crpa.UP")
  else if ((is==2) .AND. (lcrpa .eqv. .TRUE. )) then
     open(newunit=ifscr,file="Screening_W-v_crpa.DN")
  end if
  if (lomega0) then !only omega=0
     nrws1=1
     nw=1
  endif
  do ir1 = 1,nrws1
     do ifreq2 = 1,nw
        do iwf1 = 1,nwf
           do iwf2 = 1,nwf
              do iwf3 = 1, nwf
                 do iwf4 = 1, nwf
                    write(ifscr,"(' Wannier ',2i5, 3f12.6, 5i5,4f12.6)") &
                         ir1, irws1(ir1), rws1(:,ir1) &
                         ,is,iwf1,iwf2,iwf3,iwf4,freq(ifreq2), freq2(ifreq2) &
                         ,rw_w(iwf1,iwf2,iwf3,iwf4,ir1,ifreq2) &
                         ,cw_w(iwf1,iwf2,iwf3,iwf4,ir1,ifreq2)
                 enddo
              enddo
           enddo
        enddo
        !          write(ifscr,*)''
     enddo
  enddo
  close(ifscr)
  print *, "lcrpa, lomega0 =", lcrpa, lomega0
  print *, "lueff=", lueff
  if (lcrpa .eqv. .FALSE. ) then
     write(*,*) 'See "Screening_W-v.UP" and "Screening_W-v.DN" files'
  else if (lcrpa .eqv. .TRUE. ) then
     write(*,*) 'See "Screening_W-v_crpa.UP" and "Screening_W-v_crpa.DN" files'
  end if
end subroutine wwmat
subroutine wvmat (is,nwf, &
     rws1,rws2,irws1,irws2,nrws1,nrws2,nrws, &
     alat,rcut1,rcut2,rw_w,cw_w, &
     lcrpa, lomega0)
  implicit real*8(a-h,o-z)
  implicit integer (i-n)
  integer :: irws1(nrws1),irws2(nrws2)
  real(8) :: rydberg,hartree
  real(8) :: rws1(3,nrws1),rws2(3,nrws2)
  real(8) :: rw_w(nwf,nwf,nwf,nwf,nrws),cw_w(nwf,nwf,nwf,nwf,nrws)
  integer:: iwf1, iwf2, iwf3, iwf4,ir1
  logical:: lcrpa, lomega0
  write(6,*) 'start wvmat'
  hartree=2d0*rydberg()
  rw_w= hartree*rw_w
  cw_w= hartree*cw_w
  write(*,*)'Coulomb interaction (v) : '
  if (is==1) then
     open(newunit=ifcou,file="Coulomb_v.UP")
  else if (is==2) then
     open(newunit=ifcou,file="Coulomb_v.DN")
  end if
  if (lomega0) then !only omega=0
     nrws1=1
     nw=1
  endif
  do ir1 = 1,nrws1
     do iwf1 = 1,nwf
        do iwf2 = 1,nwf
           do iwf3 = 1, nwf
              do iwf4 = 1, nwf
                 write(ifcou,"(' Wannier ',2i5, 3f12.6, 5i5,2f12.6)") &
                      ir1, irws1(ir1), rws1(:,ir1), is, iwf1, iwf2, iwf3, iwf4 &
                      ,rw_w(iwf1,iwf2,iwf3,iwf4,ir1),cw_w(iwf1,iwf2,iwf3,iwf4,ir1)
              enddo
           enddo
        enddo
     enddo
     write(*,*)
  enddo
  close(ifcou)
end subroutine wvmat
subroutine chkrot ()
  use m_hamindex,only:   Readhamindex, symops, ngrp
  implicit integer (i-n)
  implicit real*8(a-h,o-z)
  parameter (eps = 1d-6)
  real(8) :: symope(3,3,ngrp)
  symope=symops
  do i = 1,3
     symope(i,i,1) = symope(i,i,1) - 1d0
  enddo

  do i = 1,3
     do j = 1,3
        as = dabs(symope(i,j,1))
        if (as > eps) stop 'chkrot: irot=1 error'
     enddo
  enddo
end subroutine chkrot
!-----------------------------------------------------------------------
subroutine wigner_seitz(alat,plat,n1,n2,n3,nrws,rws,irws,drws)
  implicit real*8(a-h,o-z)
  implicit integer (i-n)
  integer :: n1,n2,n3,nrws
  real(8) :: alat,plat(3,3)
  integer :: irws(n1*n2*n3*8)
  real(8) :: rws(3,n1*n2*n3*8),drws(n1*n2*n3*8)
  integer(4):: ii0(3,8),isort(8), &
       iwork1(n1*n2*n3*8),iwork2(n1*n2*n3*8)
  real(8):: rr(3,8),dd(8)
  parameter (tol=1.d-6)
  nrws = 0
  do i1=0,n1-1
     do i2=0,n2-1
        do i3=0,n3-1
           n = 0
           do j1=0,1
              do j2=0,1
                 do j3=0,1
                    n = n+1
                    ii0(1,n) = i1 - j1*n1
                    ii0(2,n) = i2 - j2*n2
                    ii0(3,n) = i3 - j3*n3
                 enddo ! j3
              enddo ! j2
           enddo ! j1
           do n=1,8
              rr(1:3,n) =  ( plat(1:3,1)*dble(ii0(1,n)) &
                   +   plat(1:3,2)*dble(ii0(2,n)) &
                   +   plat(1:3,3)*dble(ii0(3,n)) )
           enddo
           call sortvec2(8,rr,dd,isort)
           ndegen = 1
           do n=2,8
              if ((dd(n)-dd(1)) <= tol) ndegen = n
           enddo
           do n=1,ndegen
              nrws = nrws + 1
              rws(1:3,nrws) = rr(1:3,n)
              drws(nrws) = dd(n)
              irws(nrws) = ndegen
           enddo
        enddo ! i3
     enddo ! i2
  enddo ! i1
  call sortvec2(nrws,rws,drws,iwork1)
  iwork2(1:nrws) = irws(1:nrws)
  do n=1,nrws
     irws(n) = iwork2(iwork1(n))
  enddo
end subroutine wigner_seitz
subroutine super_cell(alat,plat,n1,n2,n3,nrws,rws,irws,drws)
  implicit real*8(a-h,o-z)
  implicit integer (i-n)
  integer(4):: n1,n2,n3,nrws
  real(8):: alat,plat(3,3)
  integer(4) :: irws(n1*n2*n3)
  real(8) :: rws(3,n1*n2*n3),drws(n1*n2*n3)
  real(8):: rr(3),dd
  parameter (tol=1.d-6)
  nrws = 0
  do i1=0,n1-1
     do i2=0,n2-1
        do i3=0,n3-1
           rr(1:3) =  ( plat(1:3,1)*dble(i1) &
                + plat(1:3,2)*dble(i2) &
                + plat(1:3,3)*dble(i3) )
           nrws = nrws + 1
           rws(1:3,nrws) = rr(1:3)
           drws(nrws) = dsqrt(sum(rr(:)**2))
           irws(nrws) = 1
        enddo ! i3
     enddo ! i2
  enddo ! i1
  if (nrws /= n1*n2*n3) stop 'super_cell: nrws error!'
  return
end subroutine super_cell
end module util_hwmatK
subroutine hwmatK_MPI() !== Calculates the bare/screened interaction W ===
  !! W(w) = <phi(n1,dR) phi(n2,0) |W(w)| phi(n3,R') phi(n4,R'+dR')>
  !! phi(n,R): maximally localized Wannier orbital
  !!    n : band index
  !!    R : site
  !!
  !! We included a cRPA method by juelich's group.
  !!
  !! This routine requirs an input from standard IO.
  !! In scripts, you can do it like, prompt>echo mode|../exec/hwmat >lwmat,
  !! where mode is 1 or 2, or 10011.
  !!
  !!  mode= 1: bare Coulomb mode, V, <phi phi | V  | phi phi>
  !!  mode= 2: screening mode, Wc,   <phi phi | Wc | phi phi>
  !!  mode=11: <phi phi | V(omega=0)  | phi phi>
  !!  mode=12: <phi phi | Wc(omega=0) | phi phi>
  !!  mode=10011: cRPA mode
  !!
  !!  nwf: total number of wannier within the primitive cell.
  !!  nrws1: index for phi1(n1,dR) phi(n2,0), where range of dR is by dR < rcut1
  !!  nrws2: index for phi(n3,R'), where range of R' is            by R' < rcut2
  !!         This is also for R'+dR'.
  !!  nrws = nrws1*nrws2*nrws2 : total spatial index for W.
  !!  Thus irws=1,nrws where specifies a set (dR,R',dR').
  !!
  ! cccccccccccccc old history ----
  !x  mode= 3: effective U mode,   <phi phi | U | phi phi>
  !x  mode= 4: similar to mode=1, but phi is atomic orbital+gaussian tail
  !x  mode= 5: similar to mode=2, but phi is atomic orbital+gaussian tail
  !x  mode= 6: similar to mode=3, but phi is atomic orbital+gaussian tail

  ! Mar 2008 Takashi Miyake, MPI version (parallelized along the outer q only)
  ! May 2007 Takashi Miyake, full matrix elements
  ! May 2006 Takashi Miyake, updated for new fpgw
  ! Sep 2004 Takashi Miyake, off-site W
  ! Jul 2004 Takashi Miyake. from hsfp0.m.f
  ! May 2002 Takashi Miyake. Total energy calc.
  ! Apr 2002 takao kotani. multiple argumentation wave per l.
  ! This hsfp0 is build from hsec10.f by F.Aryasetiawan.
  !------------------------------------------------------------
  use mpi
  use util_hwmatK,only: super_cell,chkrot,wwmat,wvmat
  use m_wmatqk,only: wmatqk_mpi
  use m_maxloc0,only:wigner_seitz
  use m_readqg,only: readngmx2,ngcmx,ngpmx,readqg0,readqg
  use m_hamindex,only:   Readhamindex,symgg=>symops,ngrp,invg=>invgx
  use m_read_bzdata,only: Read_bzdata,qibz,irkin=>irk,ginv,n1,n2,n3,nqbz,nqibz,nstar,nstbz,qbas=>qlat,qbz,wibz,wbz &
       ,nq0i=>nq0ix,wqt=>wt,q0i
  use m_readeigen,only: onoff_write_pkm4crpa,init_readeigen,init_readeigen2, &
       init_readeigen_mlw_noeval,  nwf !,init_readeigen_phi_noeval
  use m_genallcf_v3,only:niwg=>niw,alat,deltaw,esmr,icore,natom,nl,nlnmc,nlnmv,nlnmc,nlnmx,nlnx,laf
  use m_genallcf_v3,only: genallcf_v3,ncore,nn,nnc,nspin,pos,plat, nprecb,mrecb,mrece,nqbzt,nband,mrecg,ndima
  use m_keyvalue,only: getkeyvalue
  use m_zmel_old,only: ppbafp_v2
  use m_hamindex0,only: readhamindex0,iclasst

  use m_mksym_util,only:mptauof

  ! RS: MPI module
!  use rsmpi,only: rsmpi_init,mpi_comm_world,mpi_double_precision,mpi_integer,mpi_sum
!  use rsmpi_rotkindex,only: setup_rotkindex, nrot_local_rotk,irot_index_rotk
  
  implicit none
  real(8),parameter :: &
       ua    = 1d0    ! constant in w(0)exp(-ua^2*w'^2) to take care of peak around w'=0
  !------------------------------------
  real(8)    :: esmr2,shtw
  !      integer :: mxclass,ngnmax,mbytes,mwords,iwksize,
  !     &   natom,natom,ipos,igrp,
  !     &   iqibz,
  !     &   iqbz,
  !     &   iinvg,
  !     o   nspin,nl,nn,nnv,nnc,
  !     o   inindx,inindxv,inindxc,iiclass,
  !     d   nlmto,nlnx,nlnxv,nlnxc,nlnmx,nlnmxv,nlnmxc,
  !     o   iz,
  !     o   iil,iin,iim,iilnm,i_mnl,
  !     o   iilv,iinv,iimv,iilnmv,i_mnlv,
  !     o   iilc,iinc,iimc,iilnmc,i_mnlc,
  !     o   incwf,iecore,ikonf,iicore,incore,nctot,
  !     o   imagw,niw,nw,ifreq,
  integer:: ixc, ibas,ibasx,nxx,nbloch,ifqpnt,ifwd,ifmloc, &
       nprecx,mrecl,nblochpmx2,nwp,niwt, nqnum,mdimx,nblochpmx, &
       ifrcw,ifrcwi,  noccxv,maxocc2,noccx,ifvcfpout,iqall,iaf,ntq, &
       i,k,nspinmx, nq,is,ip,iq,idxk,ifoutsex,iclose,ig, &
       mxkp,nqibzxx,ntet,nene,iqi, ix,iw, &
       nlnx4,niwx,irot,invr,invrot,ivsum, ifoutsec,ntqx,&
       ifxc(2),ifsex(2), ifphiv(2),ifphic(2),ifec,ifexsp(2), &
       ifsecomg(2),ifexx,ndble=8
  real(8) :: pi,tpia,vol,voltot,rs,alpha, &
       qfermi,efx,valn,efnew,edummy,efz,qm,xsex,egex, &
       zfac1,zfac2,dscdw1,dscdw2,dscdw,zfac,ef2=1d99,exx,exxq,exxelgas
!  logical :: lqall
!  integer,allocatable :: itq(:)
  real(8),allocatable    :: q(:,:)
  ! takao
  integer,allocatable :: ngvecpB(:,:,:),ngvecp(:,:), ngvecc(:,:),iqib(:),&
       kount(:), nx(:,:),nblocha(:),lx(:) !ngveccBr(:,:,:)
  real(8),allocatable:: vxcfp(:,:,:), &
       wgt0(:,:), &
       ppbrd (:,:,:,:,:,:,:),cgr(:,:,:,:),eqt(:), &
       ppbrdx(:,:,:,:,:,:,:),aaa(:,:),& ! & symope(:,:,:)=symgg, ! qibz(:,:),
       ppb(:), eq(:), &! & ,pdb(:),dpb(:),ddb(:)
       eqx(:,:,:),eqx0(:,:,:),ekc(:),coh(:,:) &
       , rw_w(:,:,:,:,:,:),cw_w(:,:,:,:,:,:), &
       rw_iw(:,:,:,:,:,:),cw_iw(:,:,:,:,:,:)
  complex(8),allocatable:: geigB(:,:,:,:)

  logical :: screen, exchange, cohtest, legas, tote, lueff
  logical :: lcrpa
  real(8) ::  rydberg,hartree
  real(8):: qreal(3), ntot,nocctotg2,tripl,xxx(3,3)
  logical ::nocore

  ! space group infermation
  integer,allocatable :: invgx(:), miat(:,:)
  real(8),allocatable    :: tiat(:,:,:),shtvg(:,:)

  ! tetra
  real(8),allocatable :: qz(:,:),qbzxx(:),wbzxx(:),wtet(:,:,:,:), &
       eband(:,:,:), ene(:) !,ecore(:,:)
  integer,allocatable ::idtetx(:,:),idtet(:,:),ipq(:) &
       ,iene(:,:,:),ibzx(:) ! ,nstar(:)
  integer ::ib,iqx,igp,iii,ivsumxxx,isx,iflegas, iqpntnum,ndiv

  real(8),allocatable   :: eex1(:,:,:),exsp1(:,:,:),qqex1(:,:,:,:)
  integer,allocatable:: nspex(:,:),ieord(:),itex1(:,:,:)
  real(8)    :: qqex(1:3), eex,exsp,eee, exwgt,deltax0
  integer :: itmx,ipex,itpex,itex,nspexmx,nnex,isig,iex,ifexspx &
       ,ifexspxx ,ifefsm, ifemesh,nz
  character(3)  :: charnum3,sss
  character(12) :: filenameex
  logical :: exspwrite=.false.
  character(8) :: xt

  integer :: iwini,iwend
  real(8),allocatable:: omega(:,:)
  real(8) ::  omegamax,dwplot,omegamaxin
  !      logical :: sergeys

  integer::nqbze,ini,nq0it,idummy
  real(8),allocatable:: qbze(:,:)

  real(8)   :: ebmx(2)
  integer:: nbmx(2)

  real(8):: volwgt

  integer::nwin, incwfin
  real(8)::efin,ddw,dwdummy
  integer,allocatable::imdim(:)
  real(8),allocatable::freqx(:),freqw(:),wwx(:),expa(:)

  logical:: GaussSmear !readgwinput,
  integer::ret
  character*(150):: ddd


  integer::  ngpn1,verbose,ngcn1,nwxx !bzcase, mrecg,
  real(8)   :: wgtq0p,quu(3)

  real(8),allocatable:: freq_r(:)

  integer:: ifpomat,nkpo,nnmx,nomx,ikpo,nn_,no,nss(2)
  real(8):: q_r(3)
  real(8),allocatable:: qrr(:,:)
  integer,allocatable:: nnr(:),nor(:)
  complex(8),allocatable:: pomatr(:,:,:),pomat(:,:)

  real(8)::sumimg
  logical :: allq0i                                             !S.F.Jan06

  integer:: nw_i,if102,if3111,if101

  logical:: latomic,lfull,lstatic,lwssc
  logical:: l1d,lll
  real(8):: rsite(3),rcut1,rcut2
  real(8),allocatable :: rws(:,:),drws(:),rws1(:,:),rws2(:,:)
  integer:: nrws,nrws1,nrws2,ir1,ir2,ir3,ir,nrw
  integer,allocatable:: irws(:),irws1(:),irws2(:)

  ! RS: variables for MPI
  integer :: input3(3),irot_local,ip_local,iq_local, nq0ixxx
  integer,allocatable :: nq_local(:),iqx_index(:,:)
  real(8),allocatable:: &
       rw_w_sum(:,:,:,:,:,:),rw_iw_sum(:,:,:,:,:,:), &
       cw_w_sum(:,:,:,:,:,:),cw_iw_sum(:,:,:,:,:,:)

  integer :: iwf1, iwf2, iwf3, iwf4, ia
  integer :: ifcou, ifscr
  real(8),allocatable:: &
       rv_w(:,:,:,:,:),cv_w(:,:,:,:,:), &
       rv_iw(:,:,:,:,:),cv_iw(:,:,:,:,:)
  real(8)::ef,shtv(3)   ! For hmagnon (only omega=0 is used)
  integer::nw,nctot0,niw
  logical:: lomega0
  integer:: ierr,procid,master=0,comm,nrank,irr,iqibz
  integer,allocatable::irkall(:,:),irk(:,:)
  logical:: master_mpi
!  include "mpif.h"
  comm= mpi_comm_world
  call mpi_init(ierr)
  call mpi_comm_size(comm, nrank, ierr)
  call MPI_COMM_RANK(comm, procid, ierr )
  master_mpi = procid == master
!  call RSMPI_Init()
  hartree=2d0*rydberg()
  iii=verbose()
  if(master_mpi)write(6,*)' verbose=',iii
  if(master_mpi) then
     write(6,*) ' --- Choose omodes below ----------------'
     write(6,*) '  V (1) or W (2) or U(3)'
     write(6,*) '  V_omega=0 (11) or W_omega=0 (12)'
     write(6,*) '  [option --- (+ QPNT.{number} ?)] '
     write(6,*) ' --- Put number above ! -----------------'
     call readin5(ixc,nz,idummy)
     input3(1) = ixc
     input3(2) = nz
     input3(3) = idummy
     write(6,*) ' ixc nz=',ixc, nz
     if(ixc==0) stop ' --- ixc=0 --- Choose computational mode!'
  endif
  call MPI_Bcast(input3,3,MPI_INTEGER,master, MPI_COMM_WORLD,ierr)
  ixc=input3(1)
  nz=input3(2)
  idummy=input3(3)
  lomega0=.false.
  if (ixc==11) then
     ixc=1
     lomega0=.true.
  elseif(ixc==12) then
     ixc=2
     lomega0=.true.
  endif
  call read_BZDATA()
  if (master_mpi) then
     write(6,*)' nqbz  =',nqbz
     write(6,*)' nqibz ngrp=',nqibz,ngrp
  endif
  call pshpr(60)
  incwfin= -1  !use 7th colmn for core at the end section of GWIN
  call genallcf_v3(incwfin) !in module m_genallcf_v3
  niw=niwg
  ef=1d99
  if (master_mpi) write(6,*)' hsfp0: end of genallcf2'
  call pshpr(30)
  pi   = 4d0*datan(1d0)
  tpia = 2d0*pi/alat
  shtw = 0d0
  if(esmr<1d-5) shtw=0.01d0 ! Ferdi's shift to avoid resonance effect(maybe)
  voltot = abs(alat**3*tripl(plat,plat(1,2),plat(1,3)))
  legas = .false.
  latomic = .false.
  l1d = .false.
  inquire(file='UU1dU',exist=l1d)
  if (ixc==1) then
     if (master_mpi) &
          write(6,*)' --- bare Coulomb mode --- '
     exchange =.true.
     lueff = .false.
  elseif (ixc ==2) then
     if (master_mpi) &
          write(6,*)' --- screening (Wc) mode --- '
     exchange =.false.
     lueff = .false.
     lcrpa = .false.
     print *, "lcrpa =", lcrpa
     print *, "lueff=", lueff
  elseif (ixc == 100) then
     if (master_mpi) &
          write(6,*)' --- screening_crpa (Wc) mode --- '
     exchange =.false.
     lueff = .false.
     lueff = .false.
     lcrpa = .true.
     ixc = 2
  elseif(ixc==10011) then
     write(6,*) 'ixc=10011 pkm4crpa mode'
     exchange =.false.
     lueff = .false.
  else
     call rx( 'ixc error')
  endif
  if (master_mpi) then
     write(6, *) ' --- computational conditions --- '
     write(6,'("    deltaw  =",f13.6)') deltaw
     write(6,'("    ua      =",f13.6)') ua
     write(6,'("    esmr    =",f13.6)') esmr
     write(6,'("    alat voltot =",2f13.6)') alat, voltot
  endif
  call Readhamindex()
  call init_readeigen()!nband,mrece) !initialization of readEigen
  allocate(invgx(ngrp),miat(natom,ngrp),tiat(3,natom,ngrp),shtvg(3,ngrp))
  call readhamindex0()
  call mptauof(symgg,ngrp,plat,natom,pos,iclasst,miat,tiat,invgx,shtvg )
  call getsrdpp2( natom,nl,nxx)
  call readngmx2()
  if (master_mpi) write(6,*)' ngcmx ngpmx=',ngcmx,ngpmx
  allocate( nx(0:2*(nl-1),natom), nblocha(natom) ,lx(natom), &
       ppbrd ( 0:nl-1, nn, 0:nl-1,nn, 0:2*(nl-1),nxx, nspin*natom), &
       cgr(nl**2,nl**2,(2*nl-1)**2,ngrp))
  call rdpp_v3(nxx, nl, ngrp, nn, natom, nspin, symgg, nblocha, lx, nx, ppbrd, mdimx, nbloch, cgr)
  allocate(ngvecp(3,ngpmx),ngvecc(3,ngcmx))
  call readqg('QGpsi',qibz(1:3,1), quu,ngpn1, ngvecp)
  call readqg('QGcou',qibz(1:3,1), quu,ngcn1, ngvecc)
  deallocate(ngvecp,ngvecc)
  if (master_mpi) write(6,*) ' end of read QGcou'
  call pshpr(60)
  if(ixc==10011) goto 1018
  !--- Readin WV.d
  if ( .NOT. exchange) then
     if (lueff) then
        call rx('Ueff mode not implimented in the MPI version')
        open(newunit=ifwd,file='WV.d.maxloc')
     else
        open(newunit=ifwd,file='WV.d')
     endif
     read (ifwd,*) nprecx,mrecl,nblochpmx,nwp,niwt,nqnum,nw_i
     if (master_mpi) write(6,"(' Readin WV.d =', 10i5)") &
          nprecx, mrecl, nblochpmx, nwp, niwt, nqnum, nw_i
     if(nprecx/=ndble)call rx("hwmatK_MPI: WVR and WVI not compatible")!call checkeq(nprecx,ndble)
     nw=nwp-1
     if (niwt /= niw) call rx( 'hwmatK: wrong niw')
     niw = 0
     niwt = 0
     if (lueff) then
        open(newunit=ifrcw, file='__WVR.maxloc',form='unformatted',access='direct',recl=mrecl)
     else
        open(newunit=ifrcw, file='__WVR',form='unformatted',access='direct',recl=mrecl)
     endif
     !... reading general energy mesh from file 'freq_r'
     open(newunit=if3111,file='freq_r') !this is in a.u.
     read(if3111,*)nwxx             !number of energy points
     if(nwxx/= nw+1) stop ' freq_r nw /=nw'
     allocate(freq_r(nw_i:nw))       !freq_r(1)=0d0
     do iw = nw_i,nw
        read(if3111,*)freq_r(iw)
     enddo
     close(if3111)
  else
     allocate(freq_r(1)); freq_r=1d99
     nw = 0
  endif
  nrw = nw
  call getkeyvalue("GWinput","wmat_static",lstatic,default= .FALSE. )
  if (lstatic) nrw = 0
1018 continue
  call init_readeigen2()!mrecb,nlmto,mrecg) !initialize m_readeigen
  write(*,*)'nband =',nband
  lll=.false.
  if(ixc==10011 .AND. master_mpi) lll= .TRUE. 
  call onoff_write_pkm4crpa(lll)
  call init_readeigen_mlw_noeval()!nwf,nband,mrecb,mrecg)
  if (master_mpi) then
     write(*,*)'Caution! evals are zero hereafter.'
     write(*,*)'nwf =',nwf
     write(*,*)'init_readeigen_mlw: done'
  endif
  if(ixc==10011) then
     call mpi_finalize(ierr)
     if (master_mpi) call rx0s(' OK! hwmatK_MPI ixc=10011')
  endif
  nq         = nqibz
  allocate(q(3,nq))
  call dcopy   (3*nqibz,qibz,1,q,1)
  nspinmx = nspin
  if (laf) nspinmx =1
  call getkeyvalue("GWinput","wmat_all",lfull,default= .FALSE. )
  if(lfull)then
    write(6,*) 'NEEDtoExamin code main_hwmatK_MPI again ! '
    write(6,*) ' Because of the bug in nvfortran24.1 we can not pass nrws2 to wmatqk_MPI (kount,irot,nrws1,nrws2,nrws.'
    write(6,*) ' Thus we call   wmatqk_MPI (kount,irot,1,1,1 , which means fixed value is passoed to.'
    write(6,*) ' If you need wmat_all, need to fix this part. or use fixed code with ifort/gfortran'
    call rx('wmat_all is not implemented because of the bug in nvfortran24.1')
  endif
  print *, "Here!!!!!!!!!!!!!", lfull, lwssc!, nrws
  if (lfull) then
     call getkeyvalue("GWinput","wmat_rcut1",rcut1, default=0.01d0 )
     call getkeyvalue("GWinput","wmat_rcut2",rcut2, default=0.01d0 )
     call getkeyvalue("GWinput","wmat_WSsuper",lwssc,default=.true.)
     if (lwssc) then
        allocate(irws(n1*n2*n3*8),rws(3,n1*n2*n3*8),drws(n1*n2*n3*8))
        call wigner_seitz(alat,plat,n1,n2,n3,nrws,rws,irws,drws)
        if (master_mpi) then
           write(*,*)'*** Wigner-Seitz Super cell'
           do i=1,nrws
              write(*,"(i5,4f12.6,i5)")i,rws(1,i),rws(2,i),rws(3,i), &
                   drws(i),irws(i)
           enddo
        endif
     else
        allocate(irws(n1*n2*n3),rws(3,n1*n2*n3),drws(n1*n2*n3))
        call super_cell(alat,plat,n1,n2,n3,nrws,rws,irws,drws)
        if (master_mpi) then
           write(*,*)'*** Super cell (Not Wigner-Seitz super cell)'
           do i=1,nrws
              write(*,"(i5,4f12.6,i5)")i,rws(1,i),rws(2,i),rws(3,i), drws(i),irws(i)
           enddo
        endif
     endif
     nrws1 = 1
     do i=1,nrws
        if (drws(i) <= rcut1) nrws1=i
     enddo
     nrws2 = 1
     do i=1,nrws
        if (drws(i) <= rcut2) nrws2=i
     enddo
     allocate(rws1(3,nrws1),rws2(3,nrws2),irws1(nrws1),irws2(nrws2))
     rws1(:,1:nrws1) = rws(:,1:nrws1)
     rws2(:,1:nrws2) = rws(:,1:nrws2)
     irws1(1:nrws1) = irws(1:nrws1)
     irws2(1:nrws2) = irws(1:nrws2)
     nrws = nrws1*nrws2*nrws2
     deallocate(irws,rws,drws)
  else
     call getkeyvalue("GWinput","wmat_rsite", rsite,3, default=(/0.0d0,0.0d0,0.0d0/),status=ret)
     rcut1 = 0.0d0
     rcut2 = 0.0d0
     nrws1 = 1
     nrws2 = 1
     allocate(rws1(3,nrws1),rws2(3,nrws2),irws1(nrws1),irws2(nrws2))
     rws1(:,1) = rsite(:)
     rws2 = 0d0
     irws1(1) = 1
     irws2(1) = 1
     nrws = nrws1*nrws2*nrws2
  endif
  print *, "Here!!!!!!!!!!!!!", lfull, lwssc, nrws!,' Is_IO',master_mpi
  write(*,'(a14,i5,f12.6)')'nrws1, rcut1 =',nrws1,rcut1
  write(*,'(a14,i5,f12.6)')'nrws2, rcut2 =',nrws2,rcut2
  write(*,*)'rrrrrrr rws1 rws1 =',rws1,' rrrrrr2',rws2
  write(6,*) ' Used k number in Q0P =', nq0i
  write(6,"(i3,f14.6,2x, 3f14.6)" )(i, wqt(i),q0i(1:3,i),i=1,nq0i)
  allocate( wgt0(nq0i,ngrp) )
  call getkeyvalue("GWinput","allq0i",allq0i,default= .FALSE. )!S.F.Jan06
  call q0iwgt3(allq0i,symgg,ngrp,wqt,q0i,nq0i,     wgt0)                   ! added allq0i argument
  if (master_mpi) then
     if(nq0i/=0) write(6,*) ' *** tot num of q near 0   =', 1/wgt0(1,1)
     write(6,"('  sum(wgt0) from Q0P=',d14.6)")sum(wgt0)
  endif
  nqbze  = nqbz *(1 + nq0i)
  allocate( qbze(3, nqbze) )
  call dcopy(3*nqbz, qbz, 1, qbze,1)
  do i = 1,nq0i
     ini = nqbz*(1 + i -1)
     do ix=1,nqbz
        qbze (:,ini+ix)   = q0i(:,i) + qbze(:,ix)
     enddo
  enddo
  if (master_mpi) call winfo(6,nspin,nq,ntq,nspin,nbloch,ngpn1,ngcn1,nqbz,nqibz,ef,deltaw,alat,esmr)
  allocate(imdim(natom)) !bugfix 12may2015
  do ia = 1,natom
     imdim(ia)  = sum(nblocha(1:ia-1))+1
  enddo
  if(niw/=0) then ! generate gaussian frequencies x between (0,1) and w=(1-x)/x
     allocate(freqx(niw),freqw(niw),wwx(niw))!,expa(niw))
     call freq01x  (niw, freqx,freqw,wwx) 
  endif
  iii=count(irkin/=0) !ivsumxxx(irk,nqibz*ngrp)
  !  write(6,*) 'sssh shape irk=',shape(irk)
  ! write(6,*) "sssh sum of nonzero iirk=",procid,iii
  ! irk divider
  allocate(irk,mold=irkin)
  ndiv= iii/nrank
  if(mod(iii,nrank)/=0) ndiv=ndiv+1
  irk=0
  ix=0
  do iqibz=1,nqibz
     do irot=1,ngrp
        irr=irkin(iqibz,irot)
        if(irr==0) cycle
        ix=ix+1
        if(ndiv*procid<ix.and.ix<=ndiv*(procid+1)) irk(iqibz,irot)=irr
     enddo
  enddo
  ! irk divider !  call setup_rotkindex(ngrp,irk,wgt0,1,nqibz,nq0ixxx,1) ! nq=1
  
  nkpo = 1
  nnmx = 1
  nomx =1
  allocate( pomatr(nnmx,nomx,nkpo), qrr(3,nkpo),nor(nkpo),nnr(nkpo) )
  nlnx4    = nlnx**4
  niwx     = max0 (nw,niw)
  allocate( ppb(nlnmx*nlnmx*mdimx*natom),  eq(nband), &
       kount(nqibz), &
       rw_w(nwf,nwf,nwf,nwf,nrws,0:nrw), &
       cw_w(nwf,nwf,nwf,nwf,nrws,0:nrw), &
       rw_iw(nwf,nwf,nwf,nwf,nrws,niw), &
       cw_iw(nwf,nwf,nwf,nwf,nrws,niw), &
       rv_w(nwf,nwf,nwf,nwf,nrws), &
       cv_w(nwf,nwf,nwf,nwf,nrws))
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  nq0ixxx=0
  

  
  if (master_mpi) then ! debug
     if(nq0i/=0) then
        write(6,*) ' total number of k-points should be',  nqbz +  1/wgt0(1,1)   - 2 + 1 !bzcase()
     else
        write(6,*) ' total number of k-points should be',  nqbz  - 2 + 1 !bzcase()
     endif
  endif
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  
  ! ! RS: openlogfile for each process
  ! if (ixc == 1) then
  !    ifile_rsmpi = iopen ('lwt_v.MPI'//myrank_id_rsmpi,1,3,0)
  ! else if (ixc == 2) then
  !    ifile_rsmpi = iopen ('lwt_wc.MPI'//myrank_id_rsmpi,1,3,0)
  ! else if (ixc == 3) then
  !    ifile_rsmpi = iopen ('lwt_u.MPI'//myrank_id_rsmpi,1,3,0)
  ! else if (ixc == 4) then
  !    ifile_rsmpi = iopen ('lwt_v_phi.MPI'//myrank_id_rsmpi,1,3,0)
  ! else if (ixc == 5) then
  !    ifile_rsmpi = iopen ('lwt_wc_phi.MPI'//myrank_id_rsmpi,1,3,0)
  ! else if (ixc == 6) then
  !    ifile_rsmpi = iopen ('lwt_u_phi.MPI'//myrank_id_rsmpi,1,3,0)
  ! elseif( ixc==10011) then
  ! else
  !    call rx("unknown ixc")
  ! endif
  
! RS: print how symmetry operations and k-points are devided ..
!  write(ifile_rsmpi,*) "rank : ", myrank_id_rsmpi
!  write(ifile_rsmpi,*) "nrotk_local:",nk_local_qkgroup
!  write(ifile_rsmpi,*) "nrot_local :",nrot_local_rotk
!  if (nrot_local_rotk > 0) write(ifile_rsmpi,*) "iiiiii irot_index :",irot_index_rotk(1:nrot_local_rotk)
!  write(ifile_rsmpi,*) "nk_local(1:ngrp) :"
!  write(ifile_rsmpi,*) nk_local_rotk(:)
!  do irot=1,ngrp
!     if (nk_local_rotk(irot) > 0) then
!        write(ifile_rsmpi,*) "> irot,nk_local(irot) = ", irot, nk_local_rotk(irot)
!        write(ifile_rsmpi,*) "   ik_index : ",           ik_index_rotk(irot,1:nk_local_rotk(irot))
!     endif
!  enddo
  
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  if (master_mpi) write(6,*) "RS: loop over spin --"
  ! loop over spin ----------------------------------------------------
  spinloop: do 2000 is = 1,nspinmx
     write(6,*)' ssssss spinloop',is,nspinmx
     ! initialise secq and kount
     kount = 0
     rw_w = 0d0
     cw_w = 0d0
     rw_iw = 0d0
     cw_iw = 0d0
     call chkrot()
!     do ix=1,nrank
!        write(6,"('xxxirk=',2i5,' xxx ',255i3)") procid,ix,irk(ix,:),nrot_local_rotk
!     enddo
!     write(6,"('mmmxxx=',i5,' xxx ',244i3)") procid, [(irot_index_rotk(irot_local),irot_local = 1,nrot_local_rotk)]
     rotloop: do 1000 irot=1,ngrp  !irot_local = 1,nrot_local_rotk
!        irot = irot_index_rotk(irot_local)
!        if( sum(abs( irk(:,irot) )) ==0 .AND. sum(abs( wgt0(:,irot))) == 0d0 ) then
!           call rx("hwmatK_RSMPI, cylce occurs in do 1000 -loop!")
!           cycle
!        endif
        write (6,"(i3,'  out of ',i3,'  rotations ',$)") irot,ngrp
        call cputid (0)
        ! rotate atomic positions invrot*R = R' + T         !        invr       = invrot (irot,invg,ngrp)
        invr     = invg(irot)
        ! -- ppb= <Phi(SLn,r) Phi(SL'n',r) B(S,i,Rr)>
        call ppbafp_v2 (irot,ngrp,is,mdimx,lx,nx,nxx,cgr,nl-1,ppbrd, ppb)
        nctot0=0
!        write(*,*) 'wmatq in',irot_local,nrot_local_rotk
        shtv = matmul(symgg(:,:,irot),shtvg(:,invr))
        call wmatqk_MPI (kount,irot,     1,   1,   1,  tiat(1:3,1:natom,invr),miat(1:natom,invr), &
             rws1,rws2, nspin,is,  & !ifcphi,ifrb(is),ifcb(is),ifrhb(is),ifchb(is),
             ifrcw,ifrcwi, qbas,ginv,qibz,qbz,wbz,nstbz, wibz,nstar,irk,  &! & iindxk,
             nblocha,nlnmv, nlnmc,  icore,ncore, imdim, &
             ppb,    freq_r,freqx, wwx, expa, ua, dwdummy,  &! & deltaw,
             ndima,nqibz,nqbz,nctot0, &
             nl,nnc,natom,natom, &
             nlnmx,mdimx,nbloch,ngrp,nw_i,nw,nrw,niw,niwx,nq, &
             nblochpmx,ngpmx,ngcmx, &
             wgt0,wqt,nq0i,q0i, symgg(:,:,irot),alat, &
             shtv,nband, & !ifvcfpout, &
             exchange, pomatr, qrr,nnr,nor,nnmx,nomx,nkpo, nwf,  rw_w,cw_w,rw_iw,cw_iw) ! acuumulation variable
        
!        write(*,*) 'wmatq out',irot_local,nrot_local_rotk
1000 enddo rotloop
     allocate( rw_w_sum(nwf,nwf,nwf,nwf,nrws,0:nrw), &
          cw_w_sum(nwf,nwf,nwf,nwf,nrws,0:nrw), &
          rw_iw_sum(nwf,nwf,nwf,nwf,nrws,niw), &
          cw_iw_sum(nwf,nwf,nwf,nwf,nrws,niw))
!     write(6,*)'sssssss sumcheck rw_w... ',procid,sum(rw_w),sum(cw_w)
     call MPI_AllReduce(rw_w,rw_w_sum,(nrw+1)*nwf**4*nrws, MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
     rw_w = rw_w_sum
     call MPI_AllReduce(cw_w,cw_w_sum,(nrw+1)*nwf**4*nrws, MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
     cw_w = cw_w_sum
     if (niw > 0) then
        call MPI_AllReduce(rw_iw,rw_iw_sum,niw*nwf**4*nrws, MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        rw_iw = rw_iw_sum
        call MPI_AllReduce(cw_iw,cw_iw_sum,niw*nwf**4*nrws, MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        cw_iw = cw_iw_sum
     endif                   ! niw
     deallocate(rw_w_sum,cw_w_sum,rw_iw_sum,cw_iw_sum)
2001 continue      ! write <p p | W | p p>
     if (master_mpi) then
        if (exchange) then
           call wvmat (is,nwf, &
                rws1,rws2,irws1,irws2,nrws1,nrws2,nrws, &
                alat,rcut1,rcut2,rw_w(:,:,:,:,:,0),cw_w(:,:,:,:,:,0), lcrpa, lomega0)
           rv_w = rw_w(:,:,:,:,:,0)
           cv_w = cw_w(:,:,:,:,:,0)
        else
           call       wwmat (is,nw_i,nrw+1,nwf, &
                rws1,rws2,irws1,irws2,nrws1,nrws2,nrws, &
                alat,rcut1,rcut2, &
                freq_r(0:nrw), &
                rw_w,cw_w,rv_w,cv_w, &
                lcrpa, lomega0)
        endif
     endif
2000 enddo spinloop
  ! isx = iclose ('wc.d')
  ! isx = iclose ('wci.d')
  ! isx = iclose ('hbe.d')
  ! isx = iclose ('RBU')
  ! isx = iclose ('CBU')
  ! isx = iclose ('RHBU')
  ! isx = iclose ('CHBU')
  ! isx = iclose ('EVU')
  ! isx = iclose ('RBD')
  ! isx = iclose ('CBD')
  ! isx = iclose ('RHBD')
  ! isx = iclose ('CHBD')
  ! isx = iclose ('EVD')
!  call cputid(ifile_rsmpi)
  call cputid(0)
  call mpi_finalize(ierr)
  if (master_mpi) call rx0s(' OK! hwmatK_MPI')
end subroutine hwmatK_MPI
!-----------------------------------------------------------------------
