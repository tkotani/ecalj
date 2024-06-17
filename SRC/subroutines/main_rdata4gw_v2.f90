!> Generate input files used in GW calculations
!!  Read gwb.* and module variables in m_lmf2gw (set by call lmf2gw() for inputs.
!module m_rdata4gw
!  contains
subroutine rdata4gw() !bind(C)
  use m_nvfortran,only:findloc
  use m_hamindex0,only: nclass,iclass=>iclasst,nindx,lindx,ibasindx,nphimx
  use m_lmfinit,only: nsp,nbas,zz=>z,iantiferro,alat, slabl,rmt,spec_a,lmxa,ispec,nspec,nr,iantiferro,pos !bas=>pos,
  use m_lattic,only:  plat=>lat_plat,qlat=>lat_qlat
  use m_sugw,only: ndima,lmxamx,ncoremx,nqirr,konfig,ncores,ecore,gcore,gval
  use m_suham,only: nbandmx=>ham_ndham
  use m_read_bzdata,only: Read_bzdata, nqibz,qibz, nq0i,nq0iadd,wt,q0i,iq0pin
  use m_lgunit,only: m_lgunit_init
  use m_pwmat,only: mkppovl2
  use m_lgunit,only:stdo
  use m_ftox
  use m_qplist,only: qplist,ngplist
  use m_mkpot,only: qval
  use m_igv2x,only: ndimhall
!  use m_lmf2gw,only: lmf2gw, ec=>ec_d, gval=>gx_d,gcore=>gcore_d   !zz,aa,efermi,
  !check core index               : open(newunit=ifnlax,file='@MNLA_core.chk')
  !datasize                       : open(newunit=ifhbed,file='hbe.d') 
  !eigencore                      : open(newunit=ifec, file='ECORE')
  !valence phi phidot, local core: open(newunit=ifphi,file='PHIVC',form='unformatted')
  !Ording of index for CPHI.(eigen within MT): open(newunit=ifoc,file='@MNLA_CPHI')
  !eigenf within MT               : open(newunit=ifcphi,file='CPHI',form='unformatted',access='direct',recl=mrecb)
  !eigenf PW                      : open(newunit=ifgeig,file='GEIG',form='unformatted',access='direct',recl=mrecg)
  !LDA Vxc.                       : open(newunit=ifv,file='VXCFP',form='unformatted')
  !valence eigen value            : open(newunit=ifev,file='EValue',form='unformatted')
  !input for hvccfp0              : open(newunit=ifhvccfp,file='HVCCIN',form='unformatted')
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
       ,nxx,lxx,ibxx &
       ,ifnlax,inum,ibas,ic,k,icore,kkkmx,kkk,n,mmm,ix,irr,nradmx,irad &
       ,ifhbed,ndble,mrecb,mrece,ifcphi,ifec,ife0,if0c,nmax,lx,nx & !,ifcphi0
       ,mx,nnc,lmxax,ic1,irad1,irad2,l1,n1,l2,n2,nn,ierr,i1,i2,nnv &
       ,iband,ib,m,nm,iqqisp,ngp,ngc,ifhvccfp,ngrp,ig,ifiqibz &
       ,iqibz,ifnqibz,iqi,iq,ifigwin ,l,ifoc,ir,isp,  iqi2
  real(8):: ovvxx,ovvx,qpgcut_psi2,dummy
  integer,parameter:: nspx=2
  real(8) :: qm(3,3) !,ginv(3,3)
  real(8):: qqq(3),qqqa(3),qqqb(3),qx(3),ecx(2)
  integer::    ifphi,ifev,ifv,ifplane
  real(8),allocatable:: rofi(:), evl(:,:,:),vvv1(:),vvv2(:),vvv3(:)
  integer,allocatable:: nncx(:,:),nrGG(:), indxk(:),ipq(:) ,lcutmx(:) 
  real(8),allocatable :: rmax(:),gen(:,:),wgt(:) ,symgr(:,:,:),qibze(:,:)
  integer :: idxk
  complex(8),allocatable :: ppovl(:,:),ppx(:,:),ppovlinv(:,:)
  integer:: icore1,nnnn, mrecg!,nqi
  integer,allocatable:: invgx(:)
  real(8),allocatable:: shtvg(:,:)
  integer,allocatable::  konf0(:,:),     mindx(:)
  real(8),allocatable :: gxcopy (:,:,:,:,:), gx_orth (:,:,:,:,:),  gx_in(:,:,:,:,:), &
       ovv(:,:,:,:,:),zzp(:,:,:,:,:),zzpi(:,:,:,:,:),eb(:)
  complex(8),allocatable::  cphi(:,:,:,:), geig(:,:,:,:), cphix(:,:),zzpx(:,:)
  complex(8),allocatable::  cphi0(:,:,:,:), geig0(:,:,:,:), cphix0(:,:)
  integer,allocatable:: nrad(:),nindx_r(:,:),lindx_r(:,:),ncindx(:),lcindx(:)
  integer,allocatable::  iord(:,:,:,:),nvmax(:,:)
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
  integer:: nr_n
  real(8):: aa_n,bb_n,rmaxx,ab_n,delr,delrset,qu(3)
  real(8),allocatable:: rofi_n(:), gval_n (:,:,:,:,:), gcore_n  (:,:,:,:) , aann(:),bbnn(:)
  integer,allocatable:: nrnn(:)
  integer:: igp1,igp2,   icold,nev,iread !ifnband,
  integer:: igg,iorb, igc
  real(8):: platt(3,3),qtarget(3)
  logical:: newaniso, noo
  integer:: ngcB,iqx,iqindx,iq0i
  real(8):: quu(3),dQpG,abq0i,dQQ
  integer,allocatable:: ngveccB(:,:)
  integer:: idummy11(1,1)
  integer:: nxxx(1:3)=0,ngggmx,ngcgpmx,ippovlg,ippovli,ippovlgg
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
  real(8):: aac(nclass), bbc(nclass), bb(nspec)
  integer:: nrc(nclass), is,ispecc(nclass),ibasf(nbas)
  logical::laf
  integer,allocatable,target:: ngvecptt(:,:,:),ngvecctt(:,:,:),ngptt(:),ngctt(:),iqindex(:)
  real(8),allocatable:: qtt(:,:)
  integer::ibasx,ifiqg,ifiqgc,irrq, nqtt, nqnum,ngpmx,nqnumc,ngcmx,nqbz
  real(8):: QpGcut_psi,QpGcut_cou,qxx(3) !  real(8),allocatable :: gval(:,:,:,:,:)!,gcore(:,:,:,:) !ec (:,:,:),
  integer::icors(nsp),icor1,ifigwa,ispx,kkkdummy,ldummy,nr_A,icorex,nrmx
  if(verbose()>50 ) debug= .TRUE.
  
! @MNLA_core.chk
  allocate(konf0(0:lmxamx,nclass), ovv(nphimx,nphimx,0:lmxamx, nclass,nsp) )
  nxx= -999; lxx= -999; ibxx=-999
  allocate(nrad(nbas), nindx_r(ndima,nbas),lindx_r(ndima,nbas))
  nrad = 0
  open(newunit=ifnlax,file='@MNLA_core.chk')
  write(ifnlax,"(a)") '    m    n    l  icore ibas   ' ! Index for core
  write(ifnlax,'(" ------- core ------------")')
  inum = 0
  do ibas = 1,nbas
     is   = ispec(ibas)
     ic   = iclass(ibas)
     konf0(:,ic) = [(mod(konfig(l,ibas),10),l=0,lmxa(is))] 
     icore = 0             
     do l  = 0, lmxa(is)
        do kkk = l+1,konf0(l,ic)-1 !kkk is the quantum principle number
           icore = icore+1
           n  = kkk - l   ! n is starting from 1 for any l.  
           do mmm=-l,l
              inum = inum+1
              write(ifnlax,'(10i5)') mmm, n, l, icore, ibas !MNLA index. magnetic radial l numcore, ibas
           enddo
        enddo
     enddo
     if(ncores(is)/=icore) call rx('rdata xxx1: ncores(is)/=icore '//xt(ncores(is))//xt(icore) )
  enddo
  close(ifnlax)
!  
  inum=0
  allocate(mindx(ndima))
  do ix = 1,ndima
     if(ix>1 .AND. nxx == nindx(ix) .AND. lxx == lindx(ix) .AND. ibxx== ibasindx(ix)) then
        mmm=mmm+1
     else
        lxx = lindx(ix)
        ibxx= ibasindx(ix)
        nxx = nindx(ix)
        nrad(ibxx)  = nrad(ibxx)+1
        irr = nrad(ibxx)
        nindx_r(irr,ibxx) = nxx
        lindx_r(irr,ibxx) = lxx
        mmm = -lxx
     endif
     mindx(ix) = mmm
     inum = inum+1
  enddo
  nradmx= maxval(nrad)
  write(stdo,*)
  write(stdo,*) " --- Radial function indexing --- "
  write(stdo,"(' nradmx=',i5)") nradmx
  do ibas=1,nbas
     write(stdo,*)' ---- ibas nrad(ibas) =', ibas, nrad(ibas)
     do irad = 1,nrad(ibas)
        write(stdo,'("      irad=",i3," nindx_r lindx_r=",2i3)')irad, nindx_r(irad,ibas), lindx_r(irad,ibas)
     enddo
  enddo
! ECORE PHIVC 
  open(newunit=ifec, file='ECORE')
  open(newunit=ifphi,file='PHIVC',form='unformatted')
  !! --- augmentation wave numbering--- idxnlmc in index.f order
  !     m,n,l,ic order. On the other hand, fplmto gives m l n ic  for valence!
  !     This ordering is important for psi2b_v2 and psibc_v2 in x0kf_v2 and sxcf_v2, and so on.
  open(newunit=ifoc,file='@MNLA_CPHI')
  write(ifoc,"('    m    n    l ibas')")
  allocate(nvmax(0:lmxamx,nclass))
  nvmax=0
  do ix =1,ndima
     nxx = nindx(ix)
     l   = lindx(ix)
     ibas =ibasindx(ix)
     ic = iclass(ibas)
     if( nxx> nvmax(l,ic) ) nvmax(l,ic) = nxx
  enddo
  nmax = maxval(nvmax)
  allocate( iord(-lmxamx:lmxamx,nmax,0:lmxamx,nbas))
  ix=0
  iorb=0
  do ibas = 1,nbas
     ic = iclass(ibas)
     is = ispec(ibas)
     do lx = 0,lmxa(is)
        do nx = 1,nvmax(lx,ic)
           iorb=iorb+1
           do mx = -lx,lx
              ix = ix+1
              iord(mx,nx,lx,ibas)=ix
              write(ifoc,"(10i6)")mx,nx,lx,ibas,ix,iorb
           enddo
        enddo
     enddo
  enddo
  if(ix/=ndima) call rx( 'rdata4gw:ix/=ndima')
  do ibas=1,nbas !iclass is crytalographically equivalent atoms. 
     ispecc(iclass(ibas))=ispec(ibas)!note iclass(ibas) ispec(ibas). The same iclass should have the same ispec.
  enddo
  !     Get new !     nrmx, gval, gcore, aa(ic),bb(ic),nr(ic)
  !     These are replaced during this if-block.
  ! 
  !     For given two conditions;
  !       a. dr/dI (delrset() in switches.F) and the
  !       b. keeping dr/dI at r=0 (= a*b),
  !     we can deternie required nr(ic), a(ic), b(ic).
  !  if(rmeshrefine()) then
  print *,'rmeshrefine:'
  delr = delrset()       !delr is dr/di at rmax in a.u.
  write(stdo,*)' meshrefine : delr nclass=',delr,nclass
  allocate(nrnn(nclass),aann(nclass),bbnn(nclass))
  do ic = 1,nclass
     if(minval(abs(iclass-ic))/=0) cycle !jan2008
     is= ispecc(ic)
     bb(is) = rmt(is)/(exp(spec_a(is)*(nr(is)-1))-1d0)
     rmaxx = bb(is)*(exp((nr(is)-1)*spec_a(is))-1d0)
     aa_n= (delr-spec_a(is)*bb(is))/rmaxx
     bb_n= spec_a(is)*bb(is)/aa_n
     write(stdo,"(' ic aa bb=',i5,2d13.6,i5)") ic, spec_a(is),bb(is),nr(is)
     write(stdo,"(' rmaxx aa_n bb_n=',3d13.6)")rmaxx,aa_n,bb_n
     ir =0
     do
        ir = ir + 1
        if( bb_n*( exp(aa_n*(ir-1))-1d0 ) >rmaxx ) exit
     enddo
     nrnn(ic) = (ir/2)*2 + 1 !odd for simplson integral later on.
     aann(ic) = aa_n
  enddo
  nrmx = maxval(nrnn(1:nclass))
  write(stdo,*) ' New nrmx =',nrmx
  allocate( gval_n(nrmx,0:lmxamx, nphimx, nsp,nclass), gcore_n(nrmx, ncoremx, nsp,nclass) )
  do ic = 1,nclass
     if(minval(abs(iclass-ic))/=0) cycle !jan2008
     is=ispecc(ic)
     write(stdo,"('  input  nr a b =',i5,3d13.6)") nr(is),spec_a(is),bb(is)
     allocate(rofi,source = [(bb(is)*(exp((ir-1)*spec_a(is))-1d0),ir=1,nr(is))])
     nr_n = nrnn(ic)
     aa_n = aann(ic)
     bb_n = rmt(is)/(exp(aa_n*(nr_n-1))-1d0)
     allocate(rofi_n,source=[(bb_n*(exp((ir-1)*aa_n)-1d0),ir=1,nr_n)])
     do isp = 1, nsp
        do lx = 0,lmxa(is)
           do nx = 1,nvmax(lx,ic)
              call rrefine(rofi,nr(is),rofi_n,nrnn(ic), gval(1, lx,nx,isp,ic), gval_n(1, lx, nx, isp,ic) )
           enddo
        enddo
        do icore = 1,ncores(is)
           !              write(stdo,ftox)' ggg2 gcore: icore2 isp ic=',icore, isp,ic, gcore(100,icore,isp,ic)
           call rrefine(rofi,nr(is),rofi_n,nrnn(ic), gcore(1,icore,isp,ic), gcore_n(1,icore, isp,ic) )
        enddo
     enddo
     aac(ic) = aa_n
     bbc(ic) = bb_n
     nrc(ic) = nr_n
     deallocate(rofi,rofi_n)
     write(stdo,"(' output  nr a b =',i5,3d13.6)") nrc(ic),aac(ic),bbc(ic)
  enddo   ! scaled gval to avoid degeneracy of overalp OrthoNormalized
  allocate(gx_in(nrmx,0:lmxamx,nphimx,nsp,nclass), gx_orth(nrmx,0:lmxamx,nphimx,nsp,nclass) )
  nnc = 0
  lmxax = lmxamx
  allocate( zzp(nmax,nmax,0:lmxax,nclass,nsp),zzpi(nmax,nmax,0:lmxax,nclass,nsp),eb(nmax),nncx(0:lmxax,nbas),rmax(nbas))
  do ibas=1,nbas
     ic   = iclass(ibas)
     do l  = 0,lmxa(is)
        nncx(l,ibas) = konf0(l,ic) -1 -(l+1) +1
        nnc          = max(nnc,nncx(l,ibas))
     enddo
  enddo
  allocate(ncindx(ncoremx),lcindx(ncoremx))
  ncindx=-9999;lcindx=-9999
  write(ifphi) nbas, nradmx,ncoremx,nrmx
  write(ifphi) nrad(1:nbas)
  write(ifphi) nindx_r(1:nradmx,1:nbas),lindx_r(1:nradmx,1:nbas)
  ibasloop: do 1010 ibas = 1,nbas
     ic    = iclass(ibas)
     ic1   = ibas
     is    = ispec(ibas)
     allocate(rofi(nr(ic)))
     rofi = [(bbc(ic)*(exp((ir-1)*aac(ic))-1d0), ir=1,nrc(ic))]
     rmax(ibas) = rofi(nrc(ic))
     write(stdo,*)
     write(stdo,*)' ### ibas ic =',ibas,ic
     write(stdo,"(4i4,2d14.6)")  nrc(ic),lmxa(is), nsp , ncores(is), aac(ic), bbc(ic)
     write(ifec,*)            !ECORE
     write(ifec,*) slabl(is) !spid(ibas) !ECORE
     write(ifec,*) ' z,atom=class,nr,a,b,nsp ' !ECORE
     write(ifec,"(1x,f5.1,2i10,f13.5,d14.6,i4)") zz(is),ic1,nrc(ic),aac(ic),bbc(ic),nsp !ECORE
     write(ifec,*)' configuration'!   !!! LocalOrbital 2=upper 1=lower' !ECORE
     write(ifec,ftox)(konf0(l,ic),l=0,lmxa(is)) !principl quantum  number of valence minimum
     ! related to LocalOrbital part lower(=1) upper(=2).
     write(ifec,*)' l,n, ecore(up), ecore(down) ' !ECORE
     icore = 0
     do l  = 0,lmxa(is)
        do kkk = l+1 ,konf0(l,ic)-1
           icore = icore+1
           n    = kkk - l
           ncindx(icore)= n
           lcindx(icore)= l
           write(ifec,ftox) l,n,ftod(ecore(icore,1:nsp,ic),16) !,ftod(ec(icore,ic,1:nsp),16) !ECORE 
        enddo
     enddo
     write(ifphi) ncores(is), ncoremx !core
     write(ifphi) ncindx,lcindx !core
     write(ifphi) ic1,zz(is),nrc(ic),aac(ic),bbc(ic)
     write(ifphi) rofi(1:nrc(ic))
     ! --- This section it to keep the numerical stability when we have degeneracy.
     !     (mainly in the case of orthnormalized input)
     do isp = 1, nsp
        do l1  = 0,lmxa(is)
           do nn = 1, nvmax(l1,ic)
              gx_in (1:nrc(ic),l1,nn,isp,ic)=gval_n(1:nrc(ic),l1,nn,isp,ic)*sqrt(1d0+0.1d0*nn)
           enddo
        enddo
     enddo
     ovv=0d0
     do isp = 1, nsp ! ... Get overlap matrix ovv of radial functions.
        do irad1 = 1,nrad(ibas)
           do irad2 = 1,nrad(ibas)
              l1 = lindx_r (irad1,ibas)
              n1 = nindx_r (irad1,ibas)
              l2 = lindx_r (irad2,ibas)
              n2 = nindx_r (irad2,ibas)
              if(l1/=l2) cycle
              call gintxx( gx_in(1,l1,n1,isp,ic),gx_in(1,l1,n2,isp,ic),aac(ic),bbc(ic), nrc(ic), ovv(n1,n2,l1,ic,isp) )
           enddo
        enddo
     enddo
     !     ------------------------------------------------------------------
     !     ovv(1:nm,1:nm,l1,ic,isp) ---> zzp (1:nm,1:nm,l1,ic,isp),
     !     where nm = nvmax(l1,ic)
     !     \sum_j ovv(i,j)*zzp0(j,k)= e_k zzp0(i,k)
     !     zzp(i,j) = zzp0(i,k)  /sqrt(e_k)
     !     zzpi = inverse of zzp
     !     ------------------------------------------------------------------
     do isp = 1, nsp
        do l1  = 0,lmxa(is)
           n1 = nvmax(l1,ic) 
           !Get zzp : eigenfunctions of ovv
           call rss(n1, ovv(1:n1,1:n1,l1,ic,isp), eb, zzp(1:n1,1:n1,l1,ic,isp), ierr) !rs=>rss 2022-6-13
           write(stdo,"(' eb=',10f12.6)") eb(1:n1)
           if(ierr/=0) call rx( ' rdata4gw: error in rs ')
           do i1=1,n1
              do i2=1,n1
                 zzp(i1,i2,l1,ic,isp) = zzp(i1,i2,l1,ic,isp) /sqrt(eb(i2))
              enddo
           enddo
           ! ... Get zzpi : inverse of zzp
           allocate(zzpx(1:n1,1:n1))
           zzpx = zzp(1:n1,1:n1,l1,ic,isp)
           call matcinv(n1,zzpx)
           zzpi(1:n1,1:n1,l1,ic,isp) = dreal(zzpx)
           deallocate(zzpx)
           do ir=1,nrc(ic)
              gx_orth(ir,l1,1:n1,isp,ic)= matmul(gx_in(ir,l1,1:n1,isp,ic),zzp(1:n1,1:n1,l1,ic,isp))
           enddo
        enddo
     enddo
     do isp = 1, nsp
        do icore = 1, ncores(is)
           write(ifphi) gcore_n(1:nrc(ic),icore, isp,ic) ! core
        enddo
        do irad = 1,nrad(ibas)
           l = lindx_r (irad,ibas)
           n = nindx_r(irad,ibas)
           write(ifphi) gx_orth(1:nrc(ic),l, n, isp,ic) ! valence orthogonalized
           write(ifphi) gval_n (1:nrc(ic),l, n, isp,ic) ! valence raw
        enddo
        nnv = maxval(nindx(1:ndima))
     enddo
     deallocate(rofi)
     continue                  ! end of ibas loop
1010 enddo ibasloop
  close(ifphi)
  close(ifec)
! LMTO file. basic part of crystal structure.
  write(stdo,*) " --- Write LMTO file(crystal structure and so on) ---"
  ibasf=-999
  do ibas=1,nbas
     do ibasx=ibas+1,nbas !is this fine?
        if(abs(iantiferro(ibas))/=0 .AND. iantiferro(ibas)+iantiferro(ibasx)==0) then
           ibasf(ibas)=ibasx
           exit
        endif
     enddo
     if(ibasf(ibas)/=-999) write(6,"(a,2i5)")' AF pair: ibas ibasf(ibas)=',ibas,ibasf(ibas)
  enddo
  laf=.false.
  if(sum(abs(iantiferro))/=0) laf= .TRUE.
  write(stdo,ftox)' iantiferro=',iantiferro(1:nbas)
  nl=lmxax+1
  open(newunit=ifigwin,file='LMTO',form='unformatted')
  write(ifigwin) nbas,alat,plat,nsp,nl,nnv,nnc,nrmx,qval
  write(ifigwin) pos,zz(ispec(1:nbas)),slabl(ispec(1:nbas)) 
  write(ifigwin) laf,ibasf
  close(ifigwin)
  
! reading q+G and bzdata
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
! hbe.d size file  
  open(newunit=ifhbed,file='hbe.d')
  write(stdo,'( " ndima nbandmx=",3i5)') ndima, nbandmx
  ndble = 8
  mrecb = 2*ndima*nbandmx *ndble !byte size !Use -assume byterecl for ifort recognize the recored in the unit of bytes.
  mrece = nbandmx         *ndble 
  mrecg = 2*ngpmx*nbandmx *ndble 
  write(ifhbed,"(*(g0,x))") ndble,mrecb,mrece,ndima,nqbz,nbandmx,mrecg
  write(ifhbed,*)' precision, mrecl of b, mrecl of eval, ndima(p+d+l)  nqbz  nbandmx mrecg'
  close(ifhbed)
! CPHI GEIG VXCFP  
  write(stdo,*) '########## goto Eigfun part ############## '
  allocate(cphix(ndima,nbandmx), cphir(ndima,nbandmx),geigr(1:ngpmx,1:nbandmx,1:nsp))
  open(newunit=ifcphi,file='CPHI',form='unformatted',access='direct',recl=mrecb)
  open(newunit=ifgeig,file='GEIG',form='unformatted',access='direct',recl=mrecg)
  open(newunit=ifv,file='VXCFP',form='unformatted')
  write(ifv) nbandmx,nqirr
  allocate(qirr(3,nqirr),evl(nbandmx, nqirr, nsp),vxclda(nbandmx, nqirr, nsp),vvv1(nbandmx),vvv2(nbandmx),vvv3(nbandmx))
  irreducibleqloop: do 1200 iqq = 1,nqirr !irreducible points. At qirr, we calculated eigenfunctions.
     if(mod(iqq,10)==1 .OR. iqq>nqirr-5) write(stdo,*) ' iqq=',iqq
     ndimh      = ndimhall(iqq)
     ngp        = ngplist(iqq)
     qirr(:,iqq)= qplist(:,iqq)
     isploop:do 1201 isp =1,nsp
        open(newunit=ifigwb_, file='gwb'//trim(xt(iqq))//trim(xt(isp)),form='unformatted')
        geigr(1:ngpmx,1:ndimh,isp)=0d0
        evl(:,iqq,isp)=1d20
        read(ifigwb_) evl(1:ndimh,iqq,isp),cphir(1:ndima,1:ndimh),geigr(1:ngp,1:ndimh,isp),vxclda(1:ndimh,iqq,isp),nev !nev is # of eigenfunctions.
        do ibas=1,nbas
           do ix = 1,ndima
              if(ibasindx(ix)==ibas) cphir (ix,1:nev) = cphir(ix, 1:nev)/sqrt(1d0+0.1d0*nindx(ix))
           enddo
        enddo
        cphix=0d0 !     Augmentation wave part
        do iband = 1,nev
           do ix= 1,ndima
              l  = lindx(ix)
              ib = ibasindx(ix)
              n  = nindx(ix)
              m  = mindx(ix)
              ic = iclass(ib)
              nm = nvmax(l,ic)
              cphix (iord(m,1:nm,l,ib),iband) = cphix (iord(m,1:nm,l,ib),iband) + zzpi(1:nm,n,l,ic,isp)*cphir(ix,iband)
           enddo
        enddo
        iqqisp= isp + nsp*(iqq-1)
        write(ifcphi,  rec=iqqisp)  cphix(1:ndima,1:nbandmx)
        close(ifigwb_)
        iqqisp= isp + nsp*(iqq-1)
        if(ngpmx/=0) write(ifgeig,  rec=iqqisp)  geigr(1:ngpmx,1:nbandmx,isp)
1201 enddo isploop
     write(ifv) qirr(1:3,iqq), vxclda(1:nbandmx,iqq,1:nsp) ! VXCFP
1200 enddo irreducibleqloop
  close(ifgeig); print *,' end of eigensection-----'

! Evalue
  open(newunit=ifev,file='EValue',form='unformatted')
  write(ifev) nbandmx, nqirr, nsp
  write(ifev) qirr(1:3,1:nqirr) !qirr
  write(ifev) evl(1:nbandmx, 1:nqirr, 1:nsp )
  close(ifev)
  write(stdo,ftox)' iclass=',iclass,' nqirr nbas',nqirr,nbas
  do iq=1,nqirr;  write(stdo,ftox)'iq qirr=',iq,ftof(qirr(1:3,iq)); enddo
  open(newunit=ifhvccfp,file='HVCCIN',form='unformatted')
  write(ifhvccfp) alat, plat,qlat,nqirr, nbas,nbandmx
  write(ifhvccfp) qirr(:,1:nqirr), pos, rmax
  write(ifhvccfp) nqibz
  write(ifhvccfp) qibz(1:3,1:nqibz)
  write(stdo,*)
  write(stdo,"('  ngrp   = ',i3)")ngrp
  write(stdo,'("  qibz=",i3,3f12.5)')(i,qibz(1:3,i),i=1,nqibz)
  
! PPOVLG
  !! === make <Gc-Gp1+Gp2> matrix ===
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
  open(newunit=ippovlgg,file= "PPOVLGG",form='unformatted')
  write(ippovlgg) nggg, ngcgp, nqnumt-nqini+1,nqini,nqnumt
  write(ippovlgg) nvgcgp(1:3,1:ngcgp)
  write(ippovlgg) nvggg(1:3,1:nggg)
  write(ippovlgg) ggg(1:nggg)
  deallocate(ggg,nvggg,nvgcgp)
  close(ippovlgg)
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
     write(stdo,ftox)'mmmmmmmmmatcinv',sum(abs(ppovlinv))
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
  !  call rx0( ' OK! rdata4gw')
end subroutine rdata4gw
subroutine rrefine(rofio,nro,rofin,nrn,go, gn )
  implicit none 
  intent(in):: rofio,nro,rofin,nrn,go
  intent(out):: gn
  integer:: nro,nrn,ir
  real(8):: polinta,rofio(nro),rofin(nrn),go(nro),gn(nrn)
  gn = [(polinta(rofin(ir), rofio,go,nro),ir=1,nrn)]
end subroutine rrefine

!endmodule m_rdata4gw
