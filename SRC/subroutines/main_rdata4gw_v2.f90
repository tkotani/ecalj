!> Generate input files used in GW calculations
!!  Read gwb.* and module variables in m_lmf2gw (set by call lmf2gw() for inputs.
!module m_rdata4gw
!  contains
subroutine rdata4gw()
  use m_nvfortran,only:findloc
  use m_hamindex0,only: nclass,iclass=>iclasst,nindx,lindx,ibasindx,nphimx
  use m_lmfinit,only: nsp,nbas,zz=>z,iantiferro,alat, slabl,rmt,spec_a,lmxa,ispec,nspec,nr,iantiferro,pos,lmxax
  use m_lattic,only:  plat=>lat_plat,qlat=>lat_qlat
  use m_sugw,only: ndima,ncoremx,nqirr,konfig,ncores,ecore,gcore,gval,konf0
  use m_suham,only: nbandmx=>ham_ndham
  use m_read_bzdata,only: Read_bzdata, nqibz,qibz, nq0i,nq0iadd,wt,q0i,iq0pin
  use m_lgunit,only: m_lgunit_init
  use m_pwmat,only: mkppovl2
  use m_lgunit,only:stdo
  use m_qplist,only: qplist,ngplist
  use m_mkpot,only: qval
  use m_igv2x,only: ndimhall
  use m_rdata1,only: nradmx,nrmx,nnc,nrad,nindx_r,lindx_r,iord,nvmax,nrc,mindx,gval_n,gcore_n,aac,bbc,gval_orth,zzpi
  use m_ftox
  ! gvec: open(newunit=ippovlgg,file= "PPOVLGG",form='unformatted')
  ! overlap of IPW:   open(newunit=ippovl0,form='unformatted',file='PPOVL0')
  ! overlap of IPW:   open(newunit=ippovlg,file= "PPOVLG."//charnum3(iqi),form='unformatted')
  ! overlap of IPW:   open(newunit=ippovli,file= "PPOVLI."//charnum3(iqi),form='unformatted')
  ! open(newunit=ifigwin,file='LMTO',form='unformatted')
  !
  !  nbas            : the number atom in the primitive cell
  !  bas(1:3,1:nbas) : atomic posion in Cartesian coordinate (alat unit),.
  !  lmxa(1:nbas)    : maximum l number for each atom for augmentation.
  !  alat        : lattice constant in a.u.
  !  plat(1:3,1) : 1st primitive translation vector in the unit of alat
  !  plat(1:3,2) : 2nd primitive translation vector
  !  plat(1:3,3) : 3rd primitive translation vector
  !  evl (1:nbandmx)          :  eigenvalue
  !  cphi(1:ndima,1:nbandmx) : the coefficienets of eigenfunction for phi(augmentation wave)
  !  geig(1:ngp,  1:nbandmx) : the coefficienets of eigenfunction for IPW
  !  nindx   (1;ndima) : n index (phi=1 phidot=2 localorbital=3)
  !  lindx   (1:ndima) : l index
  !  ibasindx(1:ndima) : ibas index . These are used to re-ordering cphi.
  !                     :  mindx(1:ndima)  is generated under the asuumption that
  !                        m=-l:l is successively ordered.
  !
  !  nphi,     ! number of augmentation nRl channels, as distinct from:
  !  ndima,       ! number of nRLm channels = 2*ldim for phi+phidot case.
  !  nphimx       ! Maxmum number of phi for all l ---  2 for phi+phidot case.
  !  vxclda (1:nbandmx) : lda XC potential <psi|Vxc(n)|psi> for each eigenfunctions
  !
  !-  Atomic data for all the atom in the cell, ibas=1,nbas
  !   Z:
  !   nr a,b; mesh is r(i)=b*(exp(a(i-1))-1)
  !   ncore: number of core.
  !   konf(0:lmxa): principle quantum number of valence electron.
  !   ec(1:ncore)
  !   gx : radial wave function phi.
  use m_keyvalue,only: getkeyvalue
!  use m_mpi,only: MPI__Initialize
  !use m_mathlib,only: rs
  implicit none
  integer :: &
       ifigw0,ifiqc,i,ifi &
       ,ifnlax,ibas,ic,k,icore,kkkmx,kkk,n,mmm,ix,irr,irad &
       ,ifhbed,ndble,mrecb,mrece,ifcphi,ifec,ife0,if0c,lx,nx & !,ifcphi0
       ,mx,irad1,irad2,l1,n1,l2,n2,nn,ierr,i1,i2,nnv &
       ,iband,ib,m,nm,iqqisp,ngp,ngc,ifhvccfp,ngrp,ig,ifiqibz &
       ,iqibz,ifnqibz,iqi,iq,ifigwin ,l,ifoc,ir,isp,  iqi2
  real(8):: ovvxx,ovvx,qpgcut_psi2,dummy
  integer,parameter:: nspx=2
  real(8) :: qm(3,3) !,ginv(3,3)
  real(8):: qqq(3),qqqa(3),qqqb(3),qx(3),ecx(2)
  integer::    ifphi,ifev,ifv,ifplane
  real(8),allocatable:: rofi(:), evl(:,:,:),vvv1(:),vvv2(:),vvv3(:)
  integer,allocatable:: nrGG(:), indxk(:),ipq(:) ,lcutmx(:) 
  real(8),allocatable ::gen(:,:),wgt(:) ,symgr(:,:,:),qibze(:,:)
  integer :: idxk
  complex(8),allocatable :: ppovl(:,:),ppx(:,:),ppovlinv(:,:)
  integer:: icore1,nnnn, mrecg!,nqi
  integer,allocatable:: invgx(:)
  real(8),allocatable:: shtvg(:,:)
  integer,allocatable:: ncindx(:),lcindx(:)
  complex(8),allocatable::  cphi(:,:,:,:), geig(:,:,:,:), cphix(:,:),zzpx(:,:)
  complex(8),allocatable::  cphi0(:,:,:,:), geig0(:,:,:,:), cphix0(:,:)
!  integer,allocatable::  iord(:,:,:,:)!,nvmax(:,:)
  logical:: core_orth,cccxxx !readgwinput,
  real(8)::ovcv
  character(20) :: fcore
  character(3) :: charnum3
  integer,allocatable:: ifcore(:)
  real(8)::xx(5)
  real(8),allocatable:: sumc(:,:,:)
  integer:: inorm ,ib1,ib2,itest,inormo
  complex(8),allocatable:: cphir(:,:),geigr(:,:,:)
  complex(8),allocatable:: cphir0(:,:),geigr0(:,:,:)
  integer::ificg
  integer::iff,iffx
  integer,allocatable:: iffg(:)
  integer:: ifgeig,normcheck,verbose !,bzcase,ifgeig0
  integer::version
  character(8):: xt
  integer:: l_given,m_given,ic_given,isp_given,iffr,ixx
  complex(8),allocatable :: gg(:)
  real(8),allocatable:: www(:,:,:)
  logical(8):: ppovl0l=.true.
  integer:: ippovl0
  logical::  ngczero=.false.,qreduce
  integer:: nqnumt,imx     ,ix1     ,idum1,idum2,idum3, idum4 !ingczero
  real(8),allocatable::qsave(:,:)
  logical :: debug=.false.
  integer:: nqini
  integer:: ippovlc=501, zvztest
  logical:: rmeshrefine
  real(8):: aa_n,bb_n,rmaxx,ab_n  !,delr,delrset,qu(3)
  integer:: igp1,igp2,   icold,nev,iread !ifnband,
  integer:: igg,iorb, igc
  real(8):: platt(3,3),qtarget(3)
  logical:: newaniso, noo
  integer:: ngcB,iqx,iqindx,iq0i
  real(8):: quu(3),dQpG,abq0i,dQQ
  integer,allocatable:: ngveccB(:,:)
  integer:: idummy11(1,1)
  integer:: ngggmx,ngcgpmx,ippovlg,ippovli,ippovlgg
  real(8):: QpGcutggg,QpGcutgcgp,tolq=1d-8
  complex(8),allocatable:: ggg(:)
  integer,allocatable:: nvggg(:,:),nvgcgp(:,:)
  integer:: nggg,ngcgp,ixyz !q independent
  logical:: cmdopt
  character(20):: outs=''
  integer:: nl,iqtt,iqq
  real(8),parameter:: pi = 4d0*datan(1d0)
  integer,pointer:: ngvecp(:,:),ngvecc(:,:)
  integer:: ifigwb_,ifigwb0_,ndimh !ifigwx1_,ifigwx2_,
  character*256:: extn,aaa,fname
  real(8),allocatable :: qirr(:,:),vxclda(:,:,:)
!  integer::  is,ispecc(nclass),ibasf(nbas)
!  logical::laf
  integer,allocatable,target:: ngvecptt(:,:,:),ngvecctt(:,:,:),ngptt(:),ngctt(:),iqindex(:)
  real(8),allocatable:: qtt(:,:),rmax(:)
  integer::ibasx,ifiqg,ifiqgc,irrq, nqtt, nqnum,ngpmx,nqnumc,ngcmx,nqbz
  real(8):: QpGcut_psi,QpGcut_cou,qxx(3) !  real(8),allocatable :: gval(:,:,:,:,:)!,gcore(:,:,:,:) !ec (:,:,:),
  if(verbose()>50 ) debug= .TRUE.
! Reading q+G and bzdata
  open(newunit=ifiqg ,file='QGpsi',form='unformatted')
  open(newunit=ifiqgc,file='QGcou',form='unformatted')
  read(ifiqg ) nqnum, ngpmx,QpGcut_psi,nqbz,nqirr
  read(ifiqgc) nqnumc,ngcmx,QpGcut_cou
  nqtt = nqnum
  allocate(qtt(3,nqtt),ngvecptt(3,ngpmx,nqtt),ngvecctt(3,ngpmx,nqtt),ngptt(nqtt),ngctt(nqtt),iqindex(nqtt))
  irrq=0
  do iq=1,nqtt
     read(ifiqg)  qtt(1:3,iq), ngptt(iq) , irr
     read(ifiqg)  ngvecptt(1:3,1:ngptt(iq),iq)
     read(ifiqgc) qxx, ngctt(iq)
     read(ifiqgc) ngvecctt(1:3,1:ngctt(iq),iq)
     if(irr==1) then
        irrq=irrq+1
        iqindex(irrq)=iq
     endif
     if(sum(abs(qtt(:,iq)-qxx))>1d-8) call rx('QGpsi QGcou q/=q')
     write(*,"(' qtt =',i4,3f9.5)")iq, qtt(1:3,iq)
  enddo
  close(ifiqg)
  close(ifiqgc)
  call read_bzdata()
  write(6,*)'QpGcut_psi QpGcutCou =',QpGcut_psi,QpGcut_Cou
  qirr = qplist
  allocate( rmax(nbas))
  rmax=rmt(ispec(1:nbas))
  write(stdo,ftox)' iclass=',iclass,' nqirr nbas',nqirr,nbas
  do iq=1,nqirr;  write(stdo,ftox)'iq qirr=',iq,ftof(qirr(1:3,iq)); enddo
! Write HVCCIN  
  open(newunit=ifhvccfp,file='HVCCIN',form='unformatted')
  write(ifhvccfp) alat, plat,qlat,nqirr, nbas,nbandmx
  write(ifhvccfp) qirr(:,1:nqirr), pos, rmax
  write(ifhvccfp) nqibz
  write(ifhvccfp) qibz(1:3,1:nqibz)
  close(ifhvccfp)
! Generate ppovl  
  write(stdo,*)
  write(stdo,"('  ngrp   = ',i3)")ngrp
  write(stdo,'("  qibz=",i3,3f12.5)')(i,qibz(1:3,i),i=1,nqibz)
  dQpG=maxval(sum(q0i(1:3,1:nq0i+nq0iadd)**2,dim=1))**.5
  allocate(qibze(3,nqibz + nq0i+nq0iadd)) !+iadd)) !nq0i+1 is for gamma point for bzcase==2
  qibze(1:3,1:nqibz)=qibz(1:3,1:nqibz)
  qibze(1:3,nqibz+1:nqibz+nq0i+nq0iadd)=q0i(1:3,1:nq0i+nq0iadd)
  do iqi =1,nqibz + nq0i+nq0iadd; write(stdo,"(' iqi qibze=',i4,3f9.5)")iqi,qibze(:,iqi); enddo
  nqnumt= nqibz+nq0i+nq0iadd !+iadd
  write(stdo,"(' nqnumt nqibz nq0i+nq0iadd =',4i6)")nqnumt,nqibz,nq0i+nq0iadd !,iadd
  QpGcutggg = (2d0+1d-2)*QpGcut_psi + QpGcut_cou+ 2d0*pi/alat*dQpG
  if(QpGcutggg<1d-9) QpGcutggg=0d0
  dQQ= 2d0*maxval(sum(qibze(1:3,:)**2,dim=1)**.5)
  QpGcutgcgp = (1d0+1d-2)*QpGcut_psi + QpGcut_cou+ 2d0* 2d0*pi/alat*dQQ
  if(QpGcutgcgp<1d-9) QpGcutgcgp=0d0
  qx = 0d0
  nggg=1 !dummy
  ngcgp=1 !dummy
  call getgv2( alat,plat,qlat,qx, QpGcutggg, 1,nggg,  idummy11)
  call getgv2( alat,plat,qlat,qx, QpGcutgcgp,1,ngcgp, idummy11)
  if(debug) write(stdo,"(' iqi qx nggg ngcgp=',i5,3f8.3,2i8)") iqi,qx,nggg,ngcgp
  allocate( nvggg(3,nggg) )
  allocate(nvgcgp(3,ngcgp))
  !! Here the range of G for nvggg is |Gc+Gp+Gp|< |Gcou|+ |Gphi|+ |Gphi| !triangle inequality.
  qx = 0d0
  call getgv2( alat,plat,qlat,qx, QpGcutggg,  2, nggg,  nvggg  )
  call getgv2( alat,plat,qlat,qx, QpGcutgcgp, 2, ngcgp, nvgcgp )
  write(stdo,"(' getgv2--> iqi qx nggg ngcgp=',i5,3f9.3,2i7)") iqi,qx,nggg,ngcgp
  allocate(ggg(nggg))
  call mkppovl2(alat,plat,qlat, 1,  (/0,0,0/),  nggg,  nvggg, nbas, rmax, pos, ggg)
  !! ... IPW(interstitial-plane-wave) overlap matrix
  nqini =1
  if(iq0pin==2) nqini=nqibz+1
  write(stdo,*)
  do ix=1,nqirr; write(stdo,"(' qirr  = ',3f10.4)") qirr(1:3,ix); enddo
  if(debug) print *,' --- goto PPOVLG section ---',nggg,ngcgp,nqnumt-nqini+1
!Write PPOVLGG === make <Gc-Gp1+Gp2> matrix ===
  open(newunit=ippovlgg,file= "PPOVLGG",form='unformatted')
  write(ippovlgg) nggg, ngcgp, nqnumt-nqini+1,nqini,nqnumt
  write(ippovlgg) nvgcgp(1:3,1:ngcgp)
  write(ippovlgg) nvggg(1:3,1:nggg)
  write(ippovlgg) ggg(1:nggg)
  deallocate(ggg,nvggg,nvgcgp)
  close(ippovlgg)
!Write PPOVL0,PPOVLG, PPOVLI 
  print *,' --- Write PPOVLG ----'
  if(ppovl0l) open(newunit=ippovl0,form='unformatted',file='PPOVL0')
  iqiloop: do 2010 iqi = nqini, nqnumt    !nqibz + nq0i !+ iadd
     open(newunit=ippovlg,file= "PPOVLG."//charnum3(iqi),form='unformatted')
     open(newunit=ippovli,file= "PPOVLI."//charnum3(iqi),form='unformatted')
     qx  = qibze(1:3,iqi)
     iqx= findloc([(sum(abs(qx(:)-qtt(:,iqtt)))<tolq,iqtt=1,nqtt)],dim=1,value=.true.)
     ngvecp =>ngvecptt(1:3,1:ngptt(iqx),iqx)
     ngvecc =>ngvecctt(1:3,1:ngctt(iqx),iqx)
     ngp=ngptt(iqx)
     ngc=ngctt(iqx)
     write(ippovli) qx,ngc
     write(ippovlg) qx,ngc
     write(stdo,"(' iqi qx ngc iqx=',i5,3f8.4,14i5)" ) iqi,qx,ngc,iqx,sum(abs(ngvecp)),sum(abs(ngvecc)),ngc,ngp
     if(ppovl0l) write(ippovl0)   qx,ngc
     if(ngc==0) cycle
     allocate(ppovl(ngc,ngc),ppovlinv(ngc,ngc)) !This is necessary for matcinv
     call mkppovl2(alat,plat,qlat, ngc,ngvecc, ngc,ngvecc, nbas,rmax,pos, ppovl)
     if(ppovl0l)  write(ippovl0) ppovl(1:ngc,1:ngc)
     ppovlinv = ppovl
     call matcinv(ngc,ppovlinv)
     deallocate(ppovl)
     !! ggg= < exp(i G r) > integral in the interstitial region.
     if(ngc/=0) write(ippovlg) ngvecc(1:3,1:ngc)
     if(ngc/=0) write(ippovli) ppovlinv(1:ngc,1:ngc)
     deallocate(ppovlinv)
     close(ippovlg)
     close(ippovli)
2010 enddo iqiloop
  if(ppovl0l) close(ippovl0)
  write(stdo,*)" OK! end of rdata4gw "
end subroutine rdata4gw
