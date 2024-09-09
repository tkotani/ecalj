module m_lmf2gw !Note this is now only for wanplot for backward comatibility. Not in gwsc. So may unused varirables.
!> lmf2gw() set variables to module variables by reading files from following input files.
  integer,allocatable,protected:: nindx(:),lindx(:),ibasindx(:),mnla(:,:) !,iantiferro(:)
  integer,protected :: nbandmx,nphimx,&
  nsp,       &  !=1 or 2, corresponding to para or ferro.
  nbas,      &  !Number of atom in the primitive cell
  nclass,    &  !Number of atomic class (or type) for in the primitive cell
  nrmx,      &  != maxval(nr(1:nclass))  Maximum number of nr
  ncoremx,   &  != maxval(ncore(1:nclass))
  lmxamx,    &  != maxval(lmxa(1:nclass))
  ngpmx,     &  !Maximum number of G vector.
  ldim2    ! = total number of augmentation functions nRlm
  character(8),allocatable,protected:: spid(:)
  integer,allocatable:: &
       iclass(:),lmxa_d  (:),nr(:),konf_d(:,:),ncore_d(:),ibasf(:)
  real(8),allocatable :: &
       zz(:),aa(:),bb(:),ec_d (:,:,:),evl_d(:,:,:), &
       gx_d(:,:,:,:,:),gcore_d(:,:,:,:), bas(:,:)
  real(8):: plat(3,3), alat, efermi,qval
  complex(8),allocatable:: cphi_d(:,:,:,:)
  complex(8),allocatable:: geig_d(:,:,:,:)
  !      real(8),allocatable :: qirr(:,:)
  logical,protected:: laf   !! - laf: antiferro switch
  integer,protected:: ngcmx,nqnum,nqnumc,nqtt,nq0i,iq0pin,nq0iadd,nqbz,nqibz,nqbzx
  real(8):: QpGcut_psi,QpGcut_cou
  real(8),allocatable :: wt(:),q0i(:,:)
  integer,allocatable,target:: ngvecptt(:,:,:),ngvecctt(:,:,:),ngptt(:),ngctt(:)
  real(8),allocatable:: qtt(:,:)
!  integer,allocatable:: ngplist(:),ndimhall(:),iqindex(:)
!  real(8),allocatable:: qplist(:,:)
  integer:: nqirr     ! = Number of q points for irr=1 (see m_qplist, output of qg4gw).
  real(8),allocatable,protected :: qibz(:,:)
contains
  subroutine lmf2gw() !read atomic part wanplotatom.dat written in sugw.f90. This is not clean historically.
    !wanplot is expected to be unsuppported 2024-6-18
    use m_keyvalue,only: getkeyvalue
    use m_hamindex0,only: readhamindex0,nclass_in=>nclass,iclass_in=>iclasst
    use m_hamindex0,only: nindx_in=>nindx,lindx_in=>lindx,ibasindx_in=>ibasindx,nphimx_in=>nphimx
    implicit none
    integer:: iq0p
    integer:: ldim,      & ! = sum ( (lmxa(1:nbas)+1)**2 )
         nband    ! Number of bands given by GWIN0
    integer :: icor1,icorex,i,i1,i2,ibas,ibasx,ibx,ic,icore, &
         ifichkv,ifigw0,ifigwa,ifigwb,ifigwx1,ifigwx2,ifigwx3,isp,ispx, &
         ispxx,ix,kkk,kkkdummy,l,ldummy,lmxa,lxx,nclassx,m,n, &
         ncore,ndimh,ngp,nnc,nspdummy,IKP,NR_A !takao feb2012 ngc,ngcmx,
    real(8) ovv(20),ef0,z,a_a,b_a,rofi_anr
    character(120) ::  ext0(256), ext(256)
    real(8):: qqq(3),qxx(3), vvvv(18)
    integer:: ifi,ifefclass,icors(2)
    complex(8),allocatable:: zegf(:,:) ,geig(:,:)
    complex(8),allocatable:: cphi(:,:)
    real(8),allocatable:: evl(:), vvv1(:),vvv2(:),vvv3(:),rofi_A(:),gcore_A(:,:), ec_A(:)
    integer,allocatable:: konf(:,:),nncx(:,:),ngvecp(:,:),lmxaa(:)
    real(8),parameter ::  rydberg=13.6058d0
    ! nocore is obtained by inquire(file='NoCore',exist=nocore) in the upper routine.
    ! If nocore exist. you have to supply
    !  <psi|Vxc(n_valence)|psi>  to  vxclda (nband, nqirr).
    ! If not, <psi|Vxc(n_total)|psi>  to  vxclda.
    !----------------------------------------------
    integer:: ificg
    integer:: procid,nrank,ifigwb_,ifigwx1_,ifigwx2_
    integer:: iq,iqq,iqqx,nxxx,ifibz
    character*256:: extn,aaa,fname
    integer,parameter :: nsize= 1000000
    integer:: ifiproc,nqixx,nspxx,numprocxx,ixxx,ifiqibz
    integer::  id,nsizex,iqqxx,ib,ii,ipqn,nn,nnn(3),ifiqg,ifiqgc,irr,irrq,iqibz
    !! =================================================================
    open(newunit=ifigwa,file='wanplotatom.dat',form='unformatted')
    read (ifigwa) nbas,nsp,ldim2,nbandmx,lmxamx,ncoremx,nrmx,plat,alat!,nqirr
    allocate(bas(3,nbas),lmxaa(nbas))!,qplist(3,nqirr),ngplist(nqirr),ndimhall(nqirr))
    read(ifigwa) bas,lmxaa!,qplist,ngplist,ndimhall,qval
    call readhamindex0()
    nclass=nclass_in
    allocate(nindx(ldim2),lindx(ldim2),ibasindx(ldim2))
    nindx=nindx_in
    lindx=lindx_in
    ibasindx=ibasindx_in
    nphimx=nphimx_in
    allocate( iclass(nbas) )
    iclass=iclass_in
    allocate(lmxa_d(nclass), nr(nclass), ncore_d(nclass), konf_d(0:lmxamx,nclass), zz(nclass),aa(nclass),bb(nclass) )
    lmxa_d(iclass(1:nbas)) = lmxaa(1:nbas)
    !! ATOMIC PART ic = ibas scheme ==,  GET nrxx and ncoremx ----------------------
!    open(newunit=ifigwa,file='gwa',form='unformatted')
    allocate(nncx(0:lmxamx,nbas),konf(lmxamx+1,nbas),spid(nbas),ec_d(ncoremx, nclass, nsp),&
         gx_d(nrmx,0:lmxamx,nphimx,nclass,nsp), gcore_d(nrmx,ncoremx,nclass,nsp)  )
    do 3001 ibas = 1, nbas
       read(ifigwa) z, nr_A, a_A, b_A, rofi_Anr,lmxa,nspdummy,ncore,spid(ibas)
       allocate(rofi_A(nr_A), gcore_A(nr_A,ncore),ec_A(ncore))
       read(ifigwa) konf(1:lmxa+1,ibas)
       read(ifigwa) rofi_A(1:nr_A)
       write(6,"(' site',i3,'  z=',f5.1,'  rmax=',f8.5,'  lmax=',i1,'  konf=',10i1)")ibas,z,rofi_A(nr_A),lmxa,konf(1:lmxa+1,ibas)
       ic = iclass(ibas)
       zz(ic)= z
       aa(ic)= a_A
       bb(ic)= b_A
       nr(ic)= nr_A
       ncore_d(ic) = ncore/nsp
       konf_d(0:lmxa,ic) = konf(1:lmxa+1,ibas)
       write(6,"('  l    g(rmax)    gp(rmax)',4x,'<g g>',9x,'<gp gp>',9x,'<g gp>')")
       do  l = 0, lmxa
          do  isp = 1, nsp
             read(ifigwa) lxx,ispxx
             if(lxx /= l .OR. isp /=ispxx) call rx('lmf2gw:lxx or isp wrong')
             read(ifigwa) gx_d(1:nr_A,l,1,ic,isp) !phi
             read(ifigwa) gx_d(1:nr_A,l,2,ic,isp) !phidot
             if (konf_d(l,ic) >= 10) read(ifigwa) gx_d(1:nr_A,l,3,ic,isp) !phiz
          enddo
       enddo
       if(ncore/=0) write(6,'(''  l  k isp       ecore      gc(rmax)     <gc gc>'')')! core part
       icore = 0
       icors = 0
       do isp = 1, nsp
          do l = 0, lmxa
             nncx(l,ibas) = mod(konf(l+1,ibas),10)-1 -(l+1) +1
             nnc          = max(nnc,nncx(l,ibas))
             do kkk = l+1, mod(konf(l+1,ibas),10)-1
                icore = icore+1
                icors(isp) = icors(isp) +1
                icor1=icors(isp)
                read(ifigwa) icorex,ldummy,ispx,kkkdummy,ec_A(icore)
                if(icore/=icorex)  call rx('lmf2gw:icore/=icorex')
                read(ifigwa) gcore_A(1:nr_A,icore) ! gcore
                ec_d(icor1, ic, isp) = ec_A(icore)
                gcore_d(1:nr_A,icor1,ic,isp)  = gcore_A(1:nr_A,icore)
             enddo
          enddo
       enddo
       deallocate(rofi_A,gcore_A,ec_A)
3001 enddo
!    allocate(iantiferro(nbas))
!    read(ifigwa)iantiferro(1:nbas) !iantiferro may2015
    close(ifigwa)
  end subroutine lmf2gw
end module m_lmf2gw

subroutine wanplot()
  !! == Wannier function plot. Wannier function is expanded in the PW (spacial mesh). ==
  !! NOTE: Because os lazyness, not yet MPI. In cases, it may be useful...
  !! === Usage: ===
  !!  We first need to run maxloc generation as
  !!   echo 0|hbasfp0
  !!   echo 1|$nfpgw/hmaxloc   >lmaxloc1
  !!   $nfpgw/hpsig            >lpsig
  !!   echo 2|$nfpgw/huumat    >luumat2
  !!   echo 2|$nfpgw/hmaxloc   >lmaxloc2
  !!
  !! === Remarks ===
  !!   iclass means equivalent sites.
  !!   Use Xcrysden to see xsf file. Need to figure out automatic controll.
  !!   Currently a little inconvenient.
  !!   plot Wannier functions. Do this after Wannier matrix dnk is generated.
  !!
  !! ==== History ====
  !! takao modified jul2014 from calc_wannier6 by H.Kino.
  !! 080603, For plotting arbitary region, from calc_wannier
  !! 071114, Takashi Miyake, from calc_wfn.F
  !! -------------------------------------------------------------------------
  !      program calc_wannier6
  !      use m_LMTO
  !      use m_MLWF
  !      use m_wfrho_abc
  use m_keyvalue,only: getkeyvalue
  use m_readqg,only: readqg,readngmx
  use m_hamindex,only:   Readhamindex

  use m_readeigen,only: init_readeigen,init_readeigen2,readeval,lowesteval,readcphif,readgeigf,readcphifq 
  use m_read_bzdata,only: read_bzdata, nqbz,nqibz,nqbzw,nteti,ntetf &
       ,n1,n2,n3,qbas=>qlat,ginv,qbz,wbz,qibz,wibz,qbzw,idtetf,ib1bz,idteti &
       ,nstar,irk,nstbz,ngrp2=>ngrp !,qibz_r,nqibz_r
  use m_lmf2gw,only: lmf2gw,iclass,nclass,zz,alat,nbas,nsp,plat,ldim2,bas !set_mnla,
  use m_qg,only: read_qg,ngp
  use m_cubeformat
  use m_xsfformat
  use m_expand_mesh
!  use m_readhbe,only:Readhbe,nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg
  use m_genallcf_v3,only: genallcf_v3, nprecb,mrecb,mrece,nlmtot,nqbzt,nband,mrecg

  implicit none
  integer :: nwf,iko_ix,iko_fx,nband_wfn
  real(8) :: a,qdum(3)
  integer :: nq_wfn,tlat(3),nfac
  integer,allocatable :: bindx_wfn(:),bindx_wfn0(:)
  real(8),allocatable:: qreg(:,:),q_wfn(:,:)
  complex(8),allocatable :: &
       phipw0(:,:,:,:,:,:), phiaug0(:,:,:,:,:,:),phitot0(:,:,:,:,:,:), &
       phipw(:,:,:,:,:,:),  phiaug(:,:,:,:,:,:), phitot(:,:,:,:,:,:), &
       phipw1(:,:,:,:,:,:), phiaug1(:,:,:,:,:,:),phitot1(:,:,:,:,:,:), &
       wanpw(:,:,:,:,:),  wanaug(:,:,:,:,:),wantot(:,:,:,:,:), &
       wanpw1(:,:,:,:,:), wanaug1(:,:,:,:,:),wantot1(:,:,:,:,:), &
       dnk(:,:,:,:)
  integer:: i_rini(3),i_rfin(3)
  integer :: npw(3),mesh(3),mesh0(3),mesh1(3),meshrange(2,3)
  integer :: iq,ib,  ngpmx
  integer:: i,id,j,ifi,iqbz2,ifile_handle
  integer :: isp,iqbz,ikp,iwf,iopen!,nprecb,mrecb,mrece,nlmtot,nqbzt,ifhbed, mrecg,nband
  real(8):: r_rini0(3),r_rfin0(3),ang,alat_ang
  real(8):: r_rini(3), r_rfin(3),r
  real(8):: r_rini1(3), r_rfin1(3)
  logical :: lrho,flag
  character(10):: vis_unit='none'
  character(20):: outputformat ='none'
  character(20)::inputfile='GWinput'
  real(8)::quu(3),det,qlat(3,3),q(3)
  complex(8),allocatable :: cphi(:,:,:,:),geig(:,:,:,:),geig2(:,:),cphi2(:,:)
  character*4::fname
  logical:: debug=.false.,vis_skip
  integer::nqbzx,incwfin

  !! MPI dummy
  include 'mpif.h'
  integer:: ierr
  call mpi_init(ierr)
  !-----------------------------------------------------
  call getkeyvalue(inputfile,'vis_skip',  vis_skip,  default=.false. )
  if(vis_skip) call rx0s('wanplot: we found vis_skip on. Do nothing and Quit!')
  !! Readin all data
  call read_qg()
  call read_BZDATA()
  call lmf2gw()
  !      call set_mnla()
  call minv33tp (plat,qlat)  !inverse and transpose
  do isp=1,nsp
     ifi = ifile_handle()
     if (isp == 1) fname='MLWU'
     if (isp == 2) fname='MLWD'
     open(ifi,file=fname,form='unformatted',status='old', action='read')
     read(ifi)nqbzx,nwf,iko_ix,iko_fx
     if(nqbz/=nqbzx)  call rx('wanplot:nqbz/=nqbzx')
     if (isp == 1) allocate(dnk(iko_ix:iko_fx,nwf,nqbz,nsp),qreg(3,nqbz))
     do iqbz = 1,nqbz
        read(ifi)iqbz2,q(1:3)
        if(debug)write(6,"(i5,3f13.5)") iqbz,q(:)
        qreg(:,iqbz) = q
        read(ifi)dnk(iko_ix:iko_fx,1:nwf,iqbz,isp)
     enddo
     close(ifi)
  enddo                     ! isp
  write(6,*)'read end of MLWU/D ...'
  !!
  incwfin= -1  !use 7th colmn for core at the end section of GWIN
  call genallcf_v3(incwfin) !in module m_genallcf_v3
!  call Readhbe()
  call Readhamindex()
  call init_readeigen()!nband,mrece) !initialization of readEigen
  call init_readeigen2()!mrecb,ldim2,mrecg) !initialize m_readeigen

  !! replace geig and cphi with those for WF. Converted by the dnk matrix.
  !! geig2 = cphi2 = 0 for ik > nqbz
  call readngmx('QGpsi',ngpmx)
  write(6,*)'ngpmx nband ldim2=',ngpmx,nband,ldim2
  write(6,*)'nwf   nqbz  nsp  =',nwf,nqbz,nsp
  allocate(geig2(ngpmx,nband))
  allocate(cphi2(ldim2,nband))
  allocate(geig(ngpmx,nwf,nqbz,nsp))
  allocate(cphi(ldim2,nwf,nqbz,nsp))
  geig = 0d0
  cphi = 0d0
  do ikp = 1,nqbz
     do isp = 1,nsp
        geig2 = readgeigf(qreg(:,ikp),isp) !    call readgeig(qreg(:,ikp),ngpmx,isp, quu, geig2)
        cphi2 = readcphif(qreg(:,ikp),isp)
        quu   = readcphifq()
        if(sum(abs(qreg(:,ikp)-quu))>1d-6) call rx('wanplot: mmlf222eeeee')
        !! may2015 use zaxpy. This can avoid bug? when wanplot by qsub.
        do iwf = 1,nwf
           do ib  = iko_ix,iko_fx
              call zaxpy(ngp(ikp),dnk(ib,iwf,ikp,isp),geig2(1,ib),1,geig(1,iwf,ikp,isp),1)
              call zaxpy(ldim2,   dnk(ib,iwf,ikp,isp),cphi2(1,ib),1,cphi(1,iwf,ikp,isp),1)
              !           geig(:,iwf,ikp,isp) = geig(:,iwf,ikp,isp) + geig2(:,ib)*dnk(ib,iwf,ikp,isp)
              !           cphi(:,iwf,ikp,isp) = cphi(:,iwf,ikp,isp) + cphi2(:,ib)*dnk(ib,iwf,ikp,isp)
           enddo ! ib
        enddo ! iwf
     enddo ! isp
  enddo ! ikp
  deallocate(geig2,cphi2,dnk)
  write(6,*) '### ib,bas(1:3,ib) ############'
  do ib=1,nbas
     write(*,"(i5,3f12.6)")ib,bas(1:3,ib)
  enddo

  !! NOTE: nq_wfn = nqbz
  nq_wfn = nqbz
  allocate(q_wfn(3,nq_wfn))
  q_wfn(1:3,1:nqbz) = qbz(1:3,1:nqbz)

  !! == Readin vis_* settings in GWinput
  call getkeyvalue(inputfile,'vis_wan_band_n',nband_wfn,default=nwf)
  if (nband_wfn > nband) call rx('wanplot: nband_wfn > nband !')
  write(6,"(a,2i5)") '### nq_wfn, nband_wfn =',nq_wfn,nband_wfn
  allocate(bindx_wfn0(nband_wfn),bindx_wfn(nband_wfn))
  do ib=1,nband_wfn
     bindx_wfn0(ib) = ib
  enddo
  call getkeyvalue(inputfile,'vis_wan_band_id',bindx_wfn,default=bindx_wfn0,size=nband_wfn)
  do ib=1,nband_wfn
     write(*,"(a,2i5)") 'ib bndinx=',ib,bindx_wfn(ib)
  enddo
  !!
  call getkeyvalue(inputfile, 'vis_wan_tvec',   tlat,size=3,default=(/0,0,0/))
  write(*,"(a,3i5)")'### tlat',tlat
  call getkeyvalue(inputfile, 'vis_wan_interpolation', nfac,default=1)  !FFT
  vis_unit='abc'
  write(6,*)' CAUTION: range of ubond and lbound are in abc(cell) unit'
  !      call getkeyvalue(inputfile, 'vis_wan_unit', vis_unit)
  !      if ( trim(vis_unit).ne.'abc') then
  !         write(*,*) 'support only vis.wan.unit=abc'
  !         stop 'support only vis.wan.unit=abc'
  !      endif
  call calc_npw(nfac,npw)
  call getkeyvalue(inputfile,'vis_wan_mesh',  mesh0,  size=3, default=(/10,10,10/) )
  ! esh size 0:mesh0(1),...
  call getkeyvalue(inputfile,'vis_wan_lbound',r_rini0,size=3,default=(/-1d0,-1d0,-1d0/))!lower bound
  call getkeyvalue(inputfile,'vis_wan_ubound',r_rfin0,size=3,default=(/1d0,1d0,1d0/))   !upper bound
  write(*,*) ' mesh=',mesh0
  write(*,*) ' lbound=',r_rini0
  write(*,*) ' ubound=',r_rfin0
  !      call getkeyvalue(inputfile,'vis_wan_outputformat',outputformat,default='xsf')

  do i=1,3
     i_rini(i)= floor(r_rini0(i))
     i_rfin(i)= ceiling(r_rfin0(i))
  enddo
  !! for plot mesh
  r_rini= i_rini
  r_rfin = i_rfin
  do i=1,3
     mesh(i)= (i_rfin(i)-i_rini(i))*mesh0(i)
  enddo
  write(6,*)'mmm: i_rfin=',i_rini
  write(6,*)'mmm: i_rini=',i_rfin
  write(6,*)'mmm:   mesh=',mesh

  allocate(phipw0 (mesh0(1)+1,mesh0(2)+1,mesh0(3)+1,   nband_wfn,nq_wfn,nsp))
  allocate(phiaug0(mesh0(1)+1,mesh0(2)+1,mesh0(3)+1,  nband_wfn,nq_wfn,nsp))
  allocate(phitot0(mesh0(1)+1,mesh0(2)+1,mesh0(3)+1,  nband_wfn,nq_wfn,nsp))

  !! == Generate phi, which is the real-space rep. of the Bloch functions (on mesh0, real mesh points).
  !! Time consuming part.
  call calc_phiall_abc2(nq_wfn,nband_wfn,q_wfn,bindx_wfn, &
       npw,mesh0,nsp,nband,ldim2,ngpmx, &
       geig,cphi,nwf, &
       phipw0,phiaug0,phitot0)
  write(6,*)'sumchk 000 =',sum(abs(phipw0)),sum(abs(phiaug0)),sum(abs(phitot0))

  allocate(phipw(mesh(1)+1,mesh(2)+1,mesh(3)+1,   nband_wfn,nq_wfn,nsp))
  allocate(phiaug(mesh(1)+1,mesh(2)+1,mesh(3)+1,  nband_wfn,nq_wfn,nsp))
  allocate(phitot(mesh(1)+1,mesh(2)+1,mesh(3)+1,  nband_wfn,nq_wfn,nsp))

  !! phipw,phiaug,phitot are the real-space rep. on extended mesh points (mesh).
  call expand_mesh(plat, &
       nq_wfn,nband_wfn,q_wfn,nsp, &
       i_rini,i_rfin, &
       mesh0, phipw0,phiaug0,phitot0, &
       mesh, phipw,phiaug,phitot )
  !      write(6,*)'sumchk 22222=',sum(abs(phipw)),sum(abs(phiaug)),sum(abs(phitot))

  !! from Bloch to Wannier
  allocate(wanpw(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nsp) )
  allocate(wanaug(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nsp))
  allocate(wantot(mesh(1)+1,mesh(2)+1,mesh(3)+1,nband_wfn,nsp))

  call b2w(nq_wfn,nband_wfn,q_wfn,bindx_wfn,tlat,mesh,phipw,       wanpw )
  call b2w(nq_wfn,nband_wfn,q_wfn,bindx_wfn,tlat,mesh,phiaug,      wanaug)
  call b2w(nq_wfn,nband_wfn,q_wfn,bindx_wfn,tlat,mesh,phitot,      wantot)
  write(*,*) 'mesh in b2w'
  write(*,*) 'mesh=',mesh
  write(*,*) 'r_rini=',r_rini
  write(*,*) 'r_rfin=',r_rfin
  write(*,*)' '

  !!..... calculate range  from rini0 and rfin0 !abc (fractional basis only)
  !!    r = rini + (rfin-rini)*(i-1)/mesh
  meshrange=0
  r_rini1=0
  r_rfin1=0
  do id=1,3
     do i=1,mesh(id)+1
        r =  r_rini(id) + (r_rfin(id)-r_rini(id))*(i-1)/mesh(id)
        if ( r> r_rini0(id) ) then
           j=i-1
           meshrange(1,id)=j
           r_rini1(id) = r_rini(id) + (r_rfin(id)-r_rini(id))*(j-1)/mesh(id)
           exit
        endif
     enddo
     do i=1,mesh(id)+1
        r =  r_rini(id) + (r_rfin(id)-r_rini(id))*(i-1)/mesh(id)
        if ( r> r_rfin0(id) .OR. i==mesh(id)+1) then
           j=i
           meshrange(2,id)=j
           r_rfin1(id) = r_rini(id) + (r_rfin(id)-r_rini(id))*(j-1)/mesh(id)
           exit
        endif
     enddo
  enddo
  print *,'meshrange2=',meshrange(2,:) !upper limits
  print *,'meshrange1=',meshrange(1,:) !lower limits

  mesh1(:)=meshrange(2,:)-meshrange(1,:)
  allocate(wanpw1 (mesh1(1)+1,mesh1(2)+1,mesh1(3)+1,nband_wfn,nsp))
  allocate(wanaug1(mesh1(1)+1,mesh1(2)+1,mesh1(3)+1,nband_wfn,nsp))
  allocate(wantot1(mesh1(1)+1,mesh1(2)+1,mesh1(3)+1,nband_wfn,nsp))
  allocate(phipw1 (mesh1(1)+1,mesh1(2)+1,mesh1(3)+1,nband_wfn,nq_wfn,nsp))
  allocate(phiaug1(mesh1(1)+1,mesh1(2)+1,mesh1(3)+1,nband_wfn,nq_wfn,nsp))
  allocate(phitot1(mesh1(1)+1,mesh1(2)+1,mesh1(3)+1,nband_wfn,nq_wfn,nsp))

  write(*,*)'range in inputfile'
  write(*,*) 'rini0=',r_rini0
  write(*,*) 'rfin0=',r_rfin0
  write(*,*)' '
  write(*,*)'cutted mesh'
  write(*,*) 'mesh=',mesh1
  write(*,*) 'rini=',r_rini1
  write(*,*) 'rfin=',r_rfin1
  write(*,*)' '

  wanpw1(1:mesh1(1)+1,1:mesh1(2)+1,1:mesh1(3)+1, :,:) = &
       wanpw(meshrange(1,1):meshrange(1,1)+mesh1(1), &
       meshrange(1,2):meshrange(1,2)+mesh1(2), &
       meshrange(1,3):meshrange(1,3)+mesh1(3),:,: )
  wanaug1(1:mesh1(1)+1,1:mesh1(2)+1,1:mesh1(3)+1, :,:) = &
       wanaug(meshrange(1,1):meshrange(1,1)+mesh1(1), &
       meshrange(1,2):meshrange(1,2)+mesh1(2), &
       meshrange(1,3):meshrange(1,3)+mesh1(3),:,: )
  wantot1(1:mesh1(1)+1,1:mesh1(2)+1,1:mesh1(3)+1, :,:) = &
       wantot(meshrange(1,1):meshrange(1,1)+mesh1(1), &
       meshrange(1,2):meshrange(1,2)+mesh1(2), &
       meshrange(1,3):meshrange(1,3)+mesh1(3),:,: )

  phipw1(1:mesh1(1)+1,1:mesh1(2)+1,1:mesh1(3)+1, :,:,:) = &
       phipw(meshrange(1,1):meshrange(1,1)+mesh1(1), &
       meshrange(1,2):meshrange(1,2)+mesh1(2), &
       meshrange(1,3):meshrange(1,3)+mesh1(3),:,:,: )
  phiaug1(1:mesh1(1)+1,1:mesh1(2)+1,1:mesh1(3)+1, :,:,:) = &
       phiaug(meshrange(1,1):meshrange(1,1)+mesh1(1), &
       meshrange(1,2):meshrange(1,2)+mesh1(2), &
       meshrange(1,3):meshrange(1,3)+mesh1(3),:,:,: )
  phitot1(1:mesh1(1)+1,1:mesh1(2)+1,1:mesh1(3)+1, :,:,:) = &
       phitot(meshrange(1,1):meshrange(1,1)+mesh1(1), &
       meshrange(1,2):meshrange(1,2)+mesh1(2), &
       meshrange(1,3):meshrange(1,3)+mesh1(3),:,:,: )
  ! rini -> r_rini1
  ! rfin -> r_rfin1
  ! mesh -> mesh1

  !! Dump phi and wannier functions -------------------
  write(*,*) '-- dump phi(bloch funciton) and wan(Wannier funciton) --'
  qdum = 0.0d0
  ang = 0.529177d0
  alat_ang=alat*ang
  !! opendx and cube need to be fixed...
  !$$$      if (outputformat.eq.'opendx') then
  !$$$         call wfn2dx_2(alat_ang,plat,nsp,1,nband_wfn,qdum,bindx_wfn,
  !$$$     &        mesh1,r_rini1,r_rfin1,wanpw1,wanaug1,wantot1)
  !$$$         call crystal2dx_2(alat_ang,plat,r_rini1,r_rfin1,
  !$$$     &        nbas,bas,nclass,iclass,zz)
  !$$$      else if (outputformat.eq.'cube') then
  !$$$         call wrt_cube(
  !$$$     i        'wan',
  !$$$     i        alat,plat,nsp,1,nband_wfn,q_wfn,bindx_wfn,
  !$$$c     i     mesh,rini,rfin,phipw,phiaug,phitot  ! for bloch orbital
  !$$$     i        mesh1,r_rini1,r_rfin1,wanpw1,wanaug1,wantot1, ! for wannier function
  !$$$     i        nbas,bas,nclass,iclass,zz )
  !$$$         call wrt_cube(
  !$$$     i        'phi',
  !$$$     i        alat,plat,nsp,1,nband_wfn,q_wfn,bindx_wfn,
  !$$$     i        mesh1,r_rini1,r_rfin1,phipw1,phiaug1,phitot1, ! for bloch orbital
  !$$$c     i     mesh,rini,rfin,wanpw,wanaug,wantot,  ! for wannier function
  !$$$     i        nbas,bas,nclass,iclass,zz )
  !$$$      else                      !--- if(outputformat.eq.'xsf') then, default
  write(6,*) 'Writing xsf (Xcrysden) file...'
  call wrt_xsf( &
       'wan',vis_unit, &
       alat_ang,plat,nsp,1,nband_wfn,q_wfn,bindx_wfn, &
       mesh1,r_rini1,r_rfin1,wanpw1,wanaug1,wantot1,  nbas,bas,nclass,iclass,zz )
  call wrt_xsf( &
       'phi',vis_unit, &
       alat_ang,plat,nsp,1,nband_wfn,q_wfn,bindx_wfn, &
       mesh1,r_rini1,r_rfin1,phipw1,phiaug1,phitot1,  nbas,bas,nclass,iclass,zz )
  !$$$      endif

  ! dump rho or not ----------
  lrho=.false.
  write(6,"(a,l)") 'dump rho? [T/F] (if needed set lrho=T in wanplot.F)=',lrho
  if (lrho) call calc_rho_2(alat_ang,nq_wfn,nband_wfn,mesh,r_rini,r_rfin,wanpw,wanaug,wantot)
  !      if (lrho) call calc_rho(nq_wfn,nband_wfn,npw,phipw,phiaug,phitot)
  !      if (lrho) call calc_rho(1,nband_wfn,npw,wanpw,wanaug,wantot)
  !      call cputid(0)
  call rx0s('wanplot: ok')
end subroutine wanplot

