!! Generate input files used in GW calculations
!!  Read gwb.* and module variables in m_lmf2gw (set by call lmf2gw() for inputs.
!! Output files are (I think too much divided...)
!      open(newunit=ifnlax,file='@MNLA_core')
!      open(newunit=ifhbed,file='hbe.d')
!      open(newunit=ifec, file='ECORE')
!      open(newunit=ifphi,file='PHIVC',form='unformatted')
!      open(newunit=ifoc,file='@MNLA_CPHI')
!            fcore='Core_ibas'//charnum3(ibas)//'_l'//char(48+l2)//'.chk'
!            open(newunit=ifcore(l2),file=fcore)
!            open(newunit=inorm, file = "normchk.dia.chk")
!            open(newunit=inormo,file = "normchk.off.chk")
!      open(newunit=ifcphi,file='CPHI',form='unformatted',access='direct',recl=mrecb)
!      open(newunit=ifgeig,file='GEIG',form='unformatted',access='direct',recl=mrecg)
!      open(newunit=ifv,file='VXCFP',form='unformatted')
!      open(newunit=ifev,file='EValue',form='unformatted')
!      open(newunit=ifhvccfp,file='HVCCIN',form='unformatted')
!      open(newunit=ifnqibz,file='NQIBZ')
!      open(newunit=ippovlgg,file= "PPOVLGG",form='unformatted')
!      open(newunit=ippovl0,form='unformatted',file='PPOVL0')
!      open(newunit=ippovlg,file= "PPOVLG."//charnum3(iqi),form='unformatted')
!      open(newunit=ippovli,file= "PPOVLI."//charnum3(iqi),form='unformatted')
!      open(newunit=ifigwin,file='LMTO',form='unformatted')
program rdata4gw_v2
  use m_lmf2gw,only:Lmf2gw, nindx,lindx,ibasindx,nbandmx,nphimx, &
       nsp,nbas,nclass,ncoremx,lmxamx,ldim2, &
       iantiferro,spid,iclass,lmxa=>lmxa_d,nr,konf=>konf_d,ncore=>ncore_d, &
       zz,aa,bb,ec=>ec_d, gx_raw=>gx_d,gcore=>gcore_d,plat, alat, efermi, &
       bas,laf,ibasf,nqbz,  nqnum,ngpmx,ngcmx,nqnumc,nqtt, &
       QpGcut_psi,QpGcut_cou,qtt, &
       ngvecptt,ngvecctt,ngptt,ngctt, ngplist,ndimhall,iqindex,qplist,nqirr
  use m_read_bzdata,only: Read_bzdata, &
       nqibz,qibz, nq0i,nq0iadd,wt,q0i,iq0pin
  use m_lgunit,only: m_lgunit_init

  ! xx nphi,        ! number of augmentation nRl channels, as distinct from:
  ! ldim2,       ! number of nRLm channels = 2*ldim for phi+phidot case.
  ! nphimx       ! Maxmum number of phi for all l ---  2 for phi+phidot case.
  !o  hbe.d   :  datasize
  !o  CPHI    :  Coefficients of eigenfun for augmentation wave in each MT
  !o  EV[UD]  : valence eigen value
  !o  PHIVC   : valence phi phidot, local core
  !o  ECORE   : core eigenvalue
  !o  VXCFP   : LDA Vxc. VXCFPV which contains Vxc(nval) is generated for NoCore case.
  !o  LMTO    : basic date for the crystal.
  !o  HVCCIN  : required inputs for hvccfp0. Informations in this files
  !o            are also in GWIN (duplicated--it shuold be improved in future.)
  !o
  !o   NQIBZ             : number of q point ---used only for hparainfo.
  !o @MNLA_CPHI  : Checkwrite. check write. Ording of index for CPHI.
  !o @MNLA_core  : Checkwrite. index for core


  !r--- From DATA4GW_V2 (old document)-----
  !r   nbas            : the number atom in the primitive cell
  !r   bas(1:3,1:nbas) : atomic posion in Cartesian coordinate (alat unit),.
  !r   lmxa(1:nbas)    : maximum l number for each atom for augmentation.

  !r  alat        : lattice constant in a.u.
  !r  plat(1:3,1) : 1st primitive translation vector in the unit of alat
  !r  plat(1:3,2) : 2nd primitive translation vector
  !r  plat(1:3,3) : 3rd primitive translation vector

  !--   eigenfunctions for all q points.
  !r  evl (1:nbandmx)          :  eigenvalue
  !r  cphi(1:ldim2,1:nbandmx) : the coefficienets of eigenfunction for phi(augmentation wave)
  !r  geig(1:ngp,  1:nbandmx) : the coefficienets of eigenfunction for IPW
  !r
  !r   nindx   (1;ldim2) : n index (phi=1 phidot=2 localorbital=3)
  !r   lindx   (1:ldim2) : l index
  !r   ibasindx(1:ldim2) : ibas index . These are used to re-ordering cphi.
  !r                     :  mindx(1:ldim2)  is generated under the asuumption that
  !r                        m=-l:l is successively ordered.
  !r
  !r  vxclda (1:nbandmx) : lda XC potential <psi|Vxc(n)|psi> for each eigenfunctions

  !--   Atomic data for all the atom in the cell, ibas=1,nbas
  !r   Z:
  !r   nr a,b; mesh is r(i)=b*(exp(a(i-1))-1)
  !r   ncore: number of core.
  !r   konf(0:lmxa): principle quantum number of valence electron.
  !r   ec(1:ncore)
  !r   gx : radial wave function phi.
  !r
  use m_keyvalue,only: getkeyvalue
  use m_mpi,only: MPI__Initialize
  !use m_mathlib,only: rs
  implicit none
  integer :: &
       ifigw0,ifiqc,i,ifi &
       ,nxx,lxx,ibxx &
       ,ifnlax,inum,ibas,ic,k,icore,kkkmx,kkk,n,mmm,ix,irr,nradmx,irad &
       ,ifhbed,ndble,mrecb,mrece,ifcphi,ifec,ife0,if0c,nmax,lx,nx &
       ,mx,nnc,lmxax,ic1,irad1,irad2,l1,n1,l2,n2,nn,ierr,i1,i2,nnv &
       ,iband,ib,m,nm,iqqisp,ngp,ngc,ifhvccfp,ngrp,ig,ifiqibz &
       ,iqibz,ifnqibz,iqi,iq,ifigwin ,l,ifoc,ir,isp,  iqi2
  real(8):: ovvxx,ovvx,qpgcut_psi2,dummy
  integer,parameter:: nspx=2
  real(8) :: qlat(3,3) ,qm(3,3) !,ginv(3,3)
  real(8):: qqq(3),qqqa(3),qqqb(3),qx(3),ecx(2)
  integer::    ifphi,ifev,ifv,ifplane
  real(8),allocatable:: rofi(:), evl(:,:,:),vvv1(:),vvv2(:),vvv3(:)
  integer,allocatable:: lmxaa(:),nncx(:,:),nrGG(:), indxk(:),ipq(:) ,lcutmx(:) 
  real(8),allocatable :: rmax(:),gen(:,:),wgt(:) ,symgr(:,:,:),qibze(:,:)
  integer :: idxk,nwordr
  complex(8),allocatable :: ppovl(:,:),ppx(:,:),ppovlinv(:,:)
  integer:: icore1,nnnn, mrecg!,nqi
  integer,allocatable:: invgx(:)
  real(8),allocatable:: shtvg(:,:)
  integer,allocatable::  vkonf(:,:),     mindx(:)
  real(8),allocatable :: gxcopy (:,:,:,:,:), gx_orth (:,:,:,:,:),  gx_in(:,:,:,:,:), &
       ovv(:,:,:,:,:),zzp(:,:,:,:,:),zzpi(:,:,:,:,:),eb(:)
  complex(8),allocatable::  cphi(:,:,:,:), &
       geig(:,:,:,:), cphix(:,:),zzpx(:,:)
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
  integer(4)::ificg
  integer(4)::iff,iffx
  integer(4),allocatable:: iffg(:)
  integer(4):: ifgeig,iadd,normcheck,verbose !,bzcase
  integer(4)::version
  character(8):: xt
  integer(4):: l_given,m_given,ic_given,isp_given,iffr,ixx
  complex(8),allocatable :: gg(:)
  real(8),allocatable:: www(:,:,:)
  logical(8):: ppovl0l=.true.
  integer:: ippovl0,ifile_handle
  logical::  ngczero=.false.,qreduce
  integer:: nqnumt,imx     ,ix1     ,idum1,idum2,idum3, idum4 !ingczero
  real(8),allocatable::qsave(:,:)
  logical :: debug=.false.
  integer:: nqini
  integer:: ippovlc=501, zvztest
  logical:: rmeshrefine
  integer:: nr_n
  real(8):: aa_n,bb_n,rmaxx,ab_n,delr,delrset,qu(3)
  real(8),allocatable:: rofi_n(:), gx_raw_n (:,:,:,:,:), gcore_n  (:,:,:,:) , aann(:)
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
  logical:: cmdopt2
  character(20):: outs=''
  integer:: nrmx,nl,iqtt,iqq
  real(8),parameter:: pi = 4d0*datan(1d0)
  integer,pointer:: ngvecp(:,:),ngvecc(:,:)
  integer:: ifigwb_,ndimh !ifigwx1_,ifigwx2_,
  character*256:: extn,aaa,fname
  real(8),allocatable :: qirr(:,:),vxclda(:,:,:)

  !!
  call MPI__Initialize()
  call m_lgunit_init()
  call lmf2gw() !set variables, and CphiGeig
  call read_bzdata()
  if(verbose()>50 ) debug= .TRUE. 
  allocate(vkonf(0:lmxamx,nclass), ovv(nphimx,nphimx,0:lmxamx, nclass,nsp) )
  allocate(lmxaa(1:nbas))
  lmxaa(1:nbas) = lmxa(iclass(1:nbas))
  call minv33tp(plat,qlat)
  !! Make index-file @MLNA  ------------------------------------
  nxx= -999; lxx= -999; ibxx=-999
  allocate(nrad(nbas), nindx_r(ldim2,nbas),lindx_r(ldim2,nbas))
  nrad = 0
  open(newunit=ifnlax,file='@MNLA_core.chk')
  write(ifnlax,"(a)") '    m    n    l  icore ibas   ' ! Index for core
  write(ifnlax,'(" ------- core ------------")')
  inum = 0
  do ibas = 1,nbas;  ic  = iclass(ibas)
     do l = 0, lmxa(ic)
        vkonf(l,ic) = mod(konf(l,ic),10)
        ! nly 1th digit of konf has original meaning of atomic configulation
        ! ark use 10th digit for local orbital
     enddo
  enddo
  do ibas = 1,nbas
     ic   = iclass(ibas)
     icore = 0                ! --- core
     do l   = 0, lmxa(ic)
        kkkmx = vkonf(l,ic) - 1
        do kkk = l+1,kkkmx  ! kkkmx is the principle quantum number for core.
           icore = icore+1
           n    = kkk - l   ! n is starting from 1 for any l.
           ! kkk is the usual quantum principle number in text book
           do mmm=-l,l
              inum = inum+1
              write(ifnlax,'(10i5)') mmm, n, l, icore, ibas
           enddo
        enddo
     enddo
     if(ncore(ic)/=icore) then
        write(6,*)ncore(ic),icore
        call rx( ' rdata xxx1: ncore(ic)/=icore')
     endif
  enddo
  close(ifnlax)
  !! -------------------------------------------------------------
  inum=0
  allocate(mindx(ldim2))
  do ix = 1,ldim2
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
  !! radial function indexing  ----
  nradmx= maxval(nrad)
  write(6,*)
  write(6,*) " --- Radial function indexing --- "
  write(6,"(' nradmx=',i5)") nradmx
  do ibas=1,nbas
     write(6,*)' ---- ibas nrad(ibas) =', ibas, nrad(ibas)
     do irad = 1,nrad(ibas)
        write(6,'("      irad=",i3," nindx_r lindx_r=",2i3)') &
             irad, nindx_r(irad,ibas), lindx_r(irad,ibas)
     enddo
  enddo
  !!
  open(newunit=ifhbed,file='hbe.d')
  ndble = 8
  write(6,'( " ldim2 nbandmx=",3i5)') ldim2, nbandmx
  mrecb = 2*ldim2*nbandmx *ndble/nwordr()
  mrece = nbandmx         *ndble/nwordr()
  mrecg = 2*ngpmx*nbandmx *ndble/nwordr()
  write(ifhbed,"(1x,i3,6i12)") ndble,mrecb,mrece,ldim2,nqbz,nbandmx,mrecg
  write(ifhbed,*)' precision, mrecl of b, mrecl of eval, ldim2(p+d+l)  nqbz  nbandmx mrecg'
  close(ifhbed)
  open(newunit=ifec, file='ECORE')
  open(newunit=ifphi,file='PHIVC',form='unformatted')
  !! --- augmentation wave numbering--- idxnlmc in index.f order
  !     m,n,l,ic order. On the other hand, fplmto gives m l n ic  for valence!
  !     This ordering is important for psi2b_v2 and psibc_v2 in x0kf_v2 and sxcf_v2, and so on.
  open(newunit=ifoc,file='@MNLA_CPHI')
  write(ifoc,"('    m    n    l ibas')")
  allocate(nvmax(0:lmxamx,nclass))
  nvmax=0
  do ix =1,ldim2
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
     do lx = 0,lmxa(ic)
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
  if(ix/=ldim2) call rx( 'rdata4gw_v2:ix/=ldim2')
  !! --- Mesh refining  ! nov2005
  !     Get new
  !     nrmx, gx_raw, gcore, aa(ic),bb(ic),nr(ic)
  !     These are replaced during this if-block.

  !     For given two conditions;
  !     a. dr/dI (delrset() in switches.F) and the
  !     b. keeping dr/dI at r=0 (= a*b),
  !     we can deternie required nr(ic), a(ic), b(ic).
  if(rmeshrefine()) then
     print *,'rmeshrefine:'
     delr = delrset()       !delr is dr/di at rmax in a.u.
     print *,' delr nclass=',delr,nclass
     allocate(nrnn(nclass),aann(nclass))
     do ic = 1,nclass
        if(minval(abs(iclass-ic))/=0) cycle !jan2008
        rmaxx = bb(ic)*(exp((nr(ic)-1)*aa(ic))-1d0)
        aa_n= (delr-aa(ic)*bb(ic))/rmaxx
        bb_n= aa(ic)*bb(ic)/aa_n
        write(6,"(' ic aa bb=',i5,2d13.6,i5)") ic, aa(ic),bb(ic),nr(ic)
        write(6,"(' rmaxx aa_n bb_n=',3d13.6)")rmaxx,aa_n,bb_n
        ir =0
        do
           ir = ir + 1
           if( bb_n*( exp(aa_n*(ir-1))-1d0 ) >rmaxx ) exit
        enddo
        nrnn(ic) = (ir/2)*2 + 1 !odd for simplson integral later on.
        aann(ic) = aa_n
     enddo
     nrmx = maxval(nrnn(1:nclass))
     write(6,*) ' New nrmx =',nrmx
     allocate( gx_raw_n (nrmx, 0:lmxamx, nphimx, nclass,nsp), &
          gcore_n  (nrmx, ncoremx, nclass,nsp) )
     do ic = 1,nclass
        if(minval(abs(iclass-ic))/=0) cycle !jan2008
        write(6,"('  input  nr a b =',i5,3d13.6)") nr(ic),aa(ic),bb(ic)
        allocate(rofi(nr(ic)))
        do ir=1,nr(ic)
           rofi(ir) = bb(ic)*(exp((ir-1)*aa(ic))-1d0)
        enddo
        nr_n = nrnn(ic)
        aa_n = aann(ic)
        bb_n = rofi(nr(ic))/(exp(aa_n*(nr_n-1))-1d0)
        !     aa_n = aa(ic) !initial guess
        !     bb_n = bb(ic) !initial guess
        !     call getab(nr_n,nr(ic),  aa_n,bb_n) !return new a b
        allocate(rofi_n(1:nr_n))
        do ir=1,nr_n
           rofi_n(ir) = bb_n*(exp((ir-1)*aa_n)-1d0)
        enddo
        do isp = 1, nsp
           do lx = 0,lmxa(ic)
              do nx = 1,nvmax(lx,ic)
                 call rrefine(rofi,nr(ic),rofi_n,nr_n, &
                      gx_raw(1, lx, nx, ic,isp), &
                      gx_raw_n(1, lx, nx, ic,isp) )
              enddo
           enddo
           do icore = 1,ncore(ic)
              call rrefine(rofi,nr(ic),rofi_n,nr_n, &
                   gcore(1,icore, ic,isp), &
                   gcore_n(1,icore, ic,isp) )
           enddo
        enddo
        aa(ic) = aa_n
        bb(ic) = bb_n
        nr(ic) = nr_n
        deallocate(rofi,rofi_n)
        write(6,"(' output  nr a b =',i5,3d13.6)") nr(ic),aa(ic),bb(ic)
     enddo
     deallocate( gcore )
     deallocate( gx_raw)
     allocate  ( gx_raw (nrmx, 0:lmxamx, nphimx, nclass,nsp), &
          gcore  (nrmx, ncoremx, nclass,nsp) )
     gx_raw = gx_raw_n
     gcore  = gcore_n
     print *,'rmeshrefine: end'
  endif
  allocate(gx_in  (nrmx, 0:lmxamx, nphimx, nclass,nsp), &! & scaled gx_raw to avoid degeneracy of overalp
       gx_orth(nrmx, 0:lmxamx, nphimx, nclass,nsp)  ) ! OrthoNormalized
  ! ATOMIC PART ic = ibas scheme -------------------------------
  write(6,*) '########## goto atomic part ##############'
  nnc = 0
  lmxax = lmxamx
  allocate( zzp(nmax,nmax,0:lmxax,nclass,nsp),zzpi(nmax,nmax,0:lmxax,nclass,nsp), &
       eb(nmax),nncx(0:lmxax,nbas),rmax(nbas))
  do ibas=1,nbas
     ic   = iclass(ibas)
     do l  = 0,lmxa(ic)
        nncx(l,ibas) = vkonf(l,ic) -1 -(l+1) +1
        nnc          = max(nnc,nncx(l,ibas))
     enddo
  enddo
  allocate(ncindx(ncoremx),lcindx(ncoremx))
  ncindx=-9999;lcindx=-9999
  ! ... PHIVC initial
  write(ifphi) nbas, nradmx,ncoremx,nrmx
  write(ifphi) nrad(1:nbas)
  write(ifphi) nindx_r(1:nradmx,1:nbas),lindx_r(1:nradmx,1:nbas)
  do 1010 ibas = 1,nbas
     ic    = iclass(ibas)
     ic1   = ibas
     allocate(rofi(nr(ic)))
     rofi = [(bb(ic)*(exp((ir-1)*aa(ic))-1d0), ir=1,nr(ic))]
     rmax(ibas) = rofi(nr(ic))
     write(6,*)
     write(6,*)' ### ibas ic =',ibas,ic
     write(6,"(4i4,2d14.6)")  nr(ic),lmxa(ic), nsp , ncore(ic), aa(ic), bb(ic)
     write(ifec,*)          !ECORE
     write(ifec,*) spid(ibas) !ECORE
     write(ifec,*) ' z,atom=class,nr,a,b,nsp ' !ECORE
     write(ifec,6200) zz(ic),ic1,nr(ic),aa(ic),bb(ic),nsp !ECORE
6200 format(1x,f5.1,2i10,f13.5,d14.6,i4) !ECORE 31jan2006 a and b can be negative
     write(ifec,*)' configuration   !!! LocalOrbital 2=upper 1=lower' !ECORE
     write(ifec,"($,1x,10i3)") (vkonf(l,ic),  l=0,lmxa(ic)) !ECORE
     ! principl quantum  number of valence minimum
     write(ifec,"('      ',10i3)") (konf(l,ic)/10,l=0,lmxa(ic)) !ECORE
     ! related to LocalOrbital part lower(=1) upper(=2).
     write(ifec,*)' l,n, ecore(up), ecore(down) ' !ECORE
     icore = 0
     do l  = 0,lmxa(ic)
        do kkk = l+1 ,vkonf(l,ic)-1
           icore = icore+1
           n    = kkk - l
           ncindx(icore)= n
           lcindx(icore)= l
           write(ifec,"(1x,2i3,2d17.6)") l,n,ec(icore,ic,1:nsp) !ECORE
        enddo
     enddo
     write(ifphi) ncore(ic), ncoremx !core
     write(ifphi) ncindx,lcindx !core
     write(ifphi) ic1,zz(ic),nr(ic),aa(ic),bb(ic)
     write(ifphi) rofi(1:nr(ic))
     ! --- This section it to keep the numerical stability when we have degeneracy.
     !     (mainly in the case of orthnormalized input)
     do isp = 1, nsp
        do l1  = 0,lmxa(ic)
           do nn = 1, nvmax(l1,ic)
              gx_in (1:nr(ic),l1,nn,ic,isp)=gx_raw(1:nr(ic),l1,nn,ic,isp)*sqrt(1d0+0.1d0*nn)
           enddo
        enddo
     enddo
     ! ... Get overlap matrix ovv of radial functions.
     ovv=0d0
     do isp = 1, nsp
        do irad1 = 1,nrad(ibas)
           do irad2 = 1,nrad(ibas)
              l1 = lindx_r (irad1,ibas)
              n1 = nindx_r (irad1,ibas)
              l2 = lindx_r (irad2,ibas)
              n2 = nindx_r (irad2,ibas)
              if(l1/=l2) cycle
              call gintxx( gx_in(1, l1, n1, ic,isp), &
                   gx_in(1, l1, n2, ic,isp),aa(ic),bb(ic), nr(ic), ovv(n1,n2,l1,ic,isp) )
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
        do l1  = 0,lmxa(ic)
           n1 = nvmax(l1,ic)
           ! ... Get zzp : eigenfunctions of ovv
           call rss(n1, ovv(1:n1,1:n1,l1,ic,isp), eb, &
                zzp(1:n1,1:n1,l1,ic,isp), ierr) !!rs==>rss 2022-6-13
           write(6,"(' eb=',10f12.6)") eb(1:n1)
           if(ierr/=0) call rx( ' rdata4gw_v2: error in rs ')
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
           do ir=1,nr(ic)
              gx_orth(ir,l1,1:n1,ic,isp) = &
                   matmul(gx_in(ir,l1,1:n1,ic,isp),zzp(1:n1,1:n1,l1,ic,isp))
           enddo
        enddo
     enddo
     do isp = 1, nsp
        do icore = 1, ncore(ic)
           write(ifphi) gcore(1:nr(ic),icore, ic,isp) ! core
        enddo
        do irad = 1,nrad(ibas)
           l = lindx_r (irad,ibas)
           n = nindx_r(irad,ibas)
           write(ifphi) gx_orth(1:nr(ic),l, n, ic,isp) ! valence orthogonalized
           write(ifphi) gx_raw (1:nr(ic),l, n, ic,isp) ! valence raw
        enddo
        nnv = maxval(nindx(1:ldim2))
     enddo
     deallocate(rofi)
     continue                  ! end of ibas loop
1010 enddo
  close(ifphi)
  close(ifec)
  ! --- end of atomic part --------------------------------------------

  !! Band part
  write(6,*) '########## goto Eigfun part ############## '
  !$$$      if(checknorm) then
  !$$$         open(newunit=inorm, file = "normchk.dia.chk")
  !$$$         open(newunit=inormo,file = "normchk.off.chk")
  !$$$         write(inorm,849)
  !$$$         write(inormo,849)
  !$$$ 849     format('#       IPW         IPW(diag)   Onsite(tot)',
  !$$$     .      '   Onsite(phi)      Total  #check in rdata4gw ')
  !$$$         allocate(sumc(nbandmx,nbandmx,nsp))
  !$$$      endif
  allocate(cphix(ldim2,nbandmx), cphir(ldim2,nbandmx),geigr(1:ngpmx,1:nbandmx,1:nsp))
  open(newunit=ifcphi,file='CPHI',form='unformatted',access='direct',recl=mrecb)
  open(newunit=ifgeig,file='GEIG',form='unformatted',access='direct',recl=mrecg)
  open(newunit=ifv,file='VXCFP',form='unformatted')
  write(ifv) nbandmx,nqirr
  allocate( qirr(3,nqirr),evl(nbandmx, nqirr, nsp),vxclda(nbandmx, nqirr, nsp), &
       vvv1(nbandmx),vvv2(nbandmx),vvv3(nbandmx)) !,cphi(ldim2, nbandmx))!, evl(nbandmx) )
  do 1200 iqq = 1,nqirr !irreducible points. At qirr, we calculated eigenfunctions.
     if(mod(iqq,10)==1 .OR. iqq>nqirr-5) write(6,*) ' iqq=',iqq
     do 1201 isp =1,nsp
        open(newunit=ifigwb_,file='gwb'//trim(xt(iqq))//trim(xt(isp)),form='unformatted')
        ndimh  =  ndimhall(iqq)
        ngp       =  ngplist(iqq)
        qirr(:,iqq)= qplist(:,iqq)
        geigr(1:ngpmx,1:ndimh,isp)=0d0
        evl(:,iqq,isp)=1d20
        ! Read band data.
        read(ifigwb_) evl(1:ndimh,iqq,isp),cphir(1:ldim2,1:ndimh),geigr(1:ngp,1:ndimh,isp), &
             vxclda(1:ndimh, iqq,isp),nev !nev is the obtained size of eigenfunctions.
        do ibas=1,nbas
           do ix = 1,ldim2
              if(ibasindx(ix)==ibas) &
                   cphir (ix,1:nev)= cphir(ix, 1:nev)/sqrt(1d0+0.1d0*nindx(ix))
           enddo
        enddo
        !     Augmentation wave part
        cphix=0d0
        do iband = 1,nev
           do ix= 1,ldim2
              l  = lindx(ix)
              ib = ibasindx(ix)
              n  = nindx(ix)
              m  = mindx(ix)
              ic = iclass(ib)
              nm = nvmax(l,ic)
              cphix(iord(m,1:nm,l,ib),iband) = &
                   cphix(iord(m,1:nm,l,ib),iband) +  zzpi(1:nm,n,l,ic,isp) * cphir(ix,iband)
           enddo
        enddo
        iqqisp= isp + nsp*(iqq-1)
        write(ifcphi, rec=iqqisp) cphix(1:ldim2,1:nbandmx) !cphi(1:ldim2,1:nband,iqq,isp)
        !$$$        if(checknorm) then
        !$$$           do ib1=1,nev
        !$$$             do ib2=1,nev
        !$$$               if(ib1/=ib2.and.checkdia) cycle
        !$$$               sumc(ib1,ib2,isp) = sum(dconjg(cphix(1:ldim2,ib1))*cphix(1:ldim2,ib2))
        !$$$             enddo
        !$$$           enddo
        !$$$  endif
        close(ifigwb_)
1201 enddo
     !!--- Plane wave part ---
     iqx=-99999
     do iqtt=1,nqtt
        if(sum(abs(qirr(:,iqq)-qtt(:,iqtt)))<tolq) then
           iqx=iqtt
           exit
        endif
     enddo
     ngvecp =>ngvecptt(1:3,1:ngptt(iqx),iqx)
     ngvecc =>ngvecctt(1:3,1:ngctt(iqx),iqx)
     ngp=ngptt(iqx)
     ngc=ngctt(iqx)
     if(debug) print *,' iqq ngp ngc qirr =',iqq,ngp,qirr(:,iqq)
     do isp=1,nsp
        iqqisp= isp + nsp*(iqq-1)
        if(ngpmx/=0) write(ifgeig, rec=iqqisp) geigr(1:ngpmx,1:nbandmx,isp)
     enddo
     write(ifv) qirr(1:3,iqq), vxclda(1:nbandmx,iqq,1:nsp) ! VXCFP
     !$$$        if(checknorm) then
     !$$$           xx=0d0
     !$$$           qx   = qirr(1:3,iqq)
     !$$$           allocate(ppovl(ngp,ngp))
     !$$$!     ppovl = <exp (i (q+G')) | exp (i (q+G))>
     !$$$!     I think P_{k+G} for eieigenfunction is completele zero (no high-l contribution)
     !$$$!     Thus we need to use mkppovl2 for normalization check.
     !$$$           call mkppovl2(alat,plat,qlat,
     !$$$     &        ngp, ngvecp, ngp, ngvecp,  nbas, rmax, bas,     ppovl)
     !$$$           do isp=1,nsp
     !$$$           do ib1=1,nev
     !$$$           do ib2=1,nev
     !$$$             if(ib1/=ib2.and.checkdia) cycle
     !$$$             xx(3)=sum( dconjg(geigr(1:ngp,ib1,isp))*matmul(ppovl,geigr(1:ngp,ib2,isp)))
     !$$$             xx(1)=sumc(ib1,ib2,isp)
     !$$$             if(ib1==ib2) write(inorm, '(5f14.6)') xx(3),xx(4),xx(1),xx(2),xx(1)+xx(3)
     !$$$             if(ib1/=ib2) write(inormo,'(2i4,5f14.6)') ib1,ib2,xx(3),xx(1),xx(1)+xx(3)
     !$$$           enddo
     !$$$           enddo
     !$$$           enddo
     !$$$           write(inorm,*)
     !$$$           write(inormo,*)
     !$$$           deallocate(ppovl)
     !$$$        endif
1200 enddo
  print *,' end of eigensection-----'
  print *

  !!===================================================
  open(newunit=ifev,file='EValue',form='unformatted')
  write(ifev) nbandmx, nqirr, nsp
  write(ifev) qirr(1:3,1:nqirr) !qirr
  write(ifev) evl(1:nbandmx, 1:nqirr, 1:nsp )
  close(ifev)

  print *,' iclass=',iclass,' nqirr nbas',nqirr,nbas
  write(6,"(a,100i4)") ' iantiferro=',iantiferro(1:nbas)
  do iq=1,nqirr
     write(*,"(' qirr=',i4,3f9.5)")iq, qirr(1:3,iq)
  enddo
  !! --- required Input for hvccfp0 ---
  open(newunit=ifhvccfp,file='HVCCIN',form='unformatted')
  write(ifhvccfp) alat, plat,qlat,nqirr, nbas,nbandmx
  write(ifhvccfp) qirr(:,1:nqirr), bas, rmax
  write(ifhvccfp) nqibz
  write(ifhvccfp) qibz(1:3,1:nqibz)
  write(6,*)
  write(6,"('  ngrp   = ',i3)")ngrp
  write(6,'("  qibz=",i3,3f12.5)')(i,qibz(1:3,i),i=1,nqibz)

  !! === make <Gc-Gp1+Gp2> matrix ===
  iadd=0
  dQpG=0d0
  do iq0i =1,nq0i+nq0iadd
     abq0i= sum(q0i(1:3,iq0i)**2)**.5
     if(abq0i>dQpG ) dQpG=abq0i
  enddo
  allocate(qibze(3,nqibz + nq0i+nq0iadd+iadd)) !nq0i+1 is for gamma point for bzcase==2
  do iqi =1,nqibz + nq0i+nq0iadd
     if( iqi <=nqibz) then
        qibze(1:3,iqi)  = qibz(1:3,iqi)
     elseif( iqi > nqibz) then
        qibze(1:3,iqi)  = q0i (1:3,iqi-nqibz)
     endif
     write(6,"(' iqi qibze=',i4,3f9.5)")iqi,qibze(:,iqi)
  enddo
  nqnumt= nqibz+nq0i+nq0iadd+iadd
  write(6,"(' nqnumt nqibz nq0i+nq0iadd iadd =',4i6)")nqnumt,nqibz,nq0i+nq0iadd,iadd
  ! cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  QpGcutggg = (2d0+1d-2)*QpGcut_psi + QpGcut_cou+ 2d0*pi/alat*dQpG
  if(QpGcutggg<1d-9) QpGcutggg=0d0
  dQQ=0d0
  do iqi =1,nqnumt
     do iqi2=1,nqnumt
        dQQ= max( sum(qibze(1:3,iqi)**2)**.5 + sum(qibze(1:3,iqi2)**2)**.5, dQQ)
     enddo
  enddo
  QpGcutgcgp = (1d0+1d-2)*QpGcut_psi + QpGcut_cou+ 2d0* 2d0*pi/alat*dQQ
  if(QpGcutgcgp<1d-9) QpGcutgcgp=0d0
  qx = 0d0
  nggg=1 !dummy
  ngcgp=1 !dummy
  call getgv2( alat,plat,qlat,qx, QpGcutggg, 1,nggg,  idummy11)
  call getgv2( alat,plat,qlat,qx, QpGcutgcgp,1,ngcgp, idummy11)
  if(debug) write(6,"(' iqi qx nggg ngcgp=',i5,3f8.3,2i8)") iqi,qx,nggg,ngcgp
  allocate( nvggg(3,nggg) )
  allocate(nvgcgp(3,ngcgp))
  !! Here the range of G for nvggg is |Gc+Gp+Gp|< |Gcou|+ |Gphi|+ |Gphi| !triangle inequality.
  qx = 0d0
  call getgv2( alat,plat,qlat,qx, QpGcutggg,  2, nggg,  nvggg  )
  call getgv2( alat,plat,qlat,qx, QpGcutgcgp, 2, ngcgp, nvgcgp )
  write(6,"(' getgv2--> iqi qx nggg ngcgp=',i5,3f9.3,2i7)") iqi,qx,nggg,ngcgp
  allocate(ggg(nggg))
  call mkppovl2(alat,plat,qlat, 1,  (/0,0,0/),  nggg,  nvggg, nbas, rmax, bas, ggg)
  !! ... IPW(interstitial-plane-wave) overlap matrix
  nqini =1
  if(iq0pin==2) nqini=nqibz+1
  write(6,*)
  do ix=1,nqirr
     write(6,"(' qirr  = ',3f10.4)") qirr(1:3,ix)
  enddo
  !----------------------
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
  do 2010 iqi = nqini, nqnumt    !nqibz + nq0i + iadd
     open(newunit=ippovlg,file= "PPOVLG."//charnum3(iqi),form='unformatted')
     open(newunit=ippovli,file= "PPOVLI."//charnum3(iqi),form='unformatted')
     qx  = qibze(1:3,iqi)
     iqx=-99999
     do iqtt=1,nqtt
        if(sum(abs(qx(:)-qtt(:,iqtt)))<tolq) then
           iqx=iqtt
           exit
        endif
     enddo
     ngvecp =>ngvecptt(1:3,1:ngptt(iqx),iqx)
     ngvecc =>ngvecctt(1:3,1:ngctt(iqx),iqx)
     ngp=ngptt(iqx)
     ngc=ngctt(iqx)
     write(ippovli) qx,ngc
     write(ippovlg) qx,ngc
     write(6,"(' iqi qx ngc=',i5,3f8.4,4i5)" ) iqi,qx,ngc
     if(ppovl0l) write(ippovl0)   qx,ngc
     if(ngc==0) then
        cycle
     endif
     allocate(ppovl(ngc,ngc),ppovlinv(ngc,ngc)) !This is necessary for matcinv
     call mkppovl2(alat,plat,qlat, &
          ngc,  ngvecc, &
          ngc,  ngvecc, &
          nbas, rmax, bas, &
          ppovl)
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
2010 enddo
  if(ppovl0l) close(ippovl0)
  write(6,*) " --- Write LMTO file(crystal structure and so on) ---"
  nl=lmxax+1
  open(newunit=ifigwin,file='LMTO',form='unformatted')
  write(ifigwin) nbas,alat,plat,nsp,nl,nnv,nnc,nrmx!,n1q,n2q,n3q
  write(ifigwin) bas,zz(iclass(1:nbas)),spid(1:nbas)
  write(ifigwin) laf,ibasf
  close(ifigwin)
  print *, " OK! end of rdata4gw_v2 "
  call rx0( ' OK! rdata4gw_v2')
END PROGRAM rdata4gw_v2

! ssssssssssssssssssssssssssssssssssssssssssssssssssss
subroutine rrefine(rofio,nro,rofin,nrn,go, gn )
  implicit none
  integer:: nro,nrn,ir
  real(8):: polinta,rofio(nro),rofin(nrn),go(nro),gn(nrn)
  gn = [(polinta(rofin(ir), rofio,go,nro),ir=1,nrn)]
end subroutine rrefine
