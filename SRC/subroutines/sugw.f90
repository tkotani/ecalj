!> Generate all the inputs for GW calculation. Need q+G info from QGpsi and QGcou which are generated a qg4gw.
module m_sugw
  use m_mpiio,only: openm,writem,closem
  real(8),allocatable,public::ecore(:,:,:),gcore(:,:,:,:),gval(:,:,:,:,:)
  integer,public::   ndham, nqirr,nqibz    !ndima, ncoremx,
!  integer,allocatable,public::  konf0(:,:) !konfig(:,:),ncores(:),
  private
  public:: m_sugw_init
contains
  subroutine m_sugw_init (socmatrix,eferm,vmag,qval) !Driver for GW calculation
    use m_lgunit,only:stdo
    use m_lmfinit,only: zz=>z,nris=>nr,lmxa,rmt,spec_a, konfig,ncores,ndimaa,ncoremx,ndima,konf0
    use m_ext,only:   sname
    use m_lmfinit,only: ham_pwmode,pwemax,ham_oveps,lrsig=>ham_lsig,nlmto,lso,nspx !nspx=nsp/nspc
    use m_lmfinit,only: ham_scaledsigma,alat=>lat_alat,nsp,nspc,ispec,nspec,pos
    use m_lmfinit,only: nbas,n0,nppn,nkap0,slabl,nmcorex=>nmcore,iantiferro,lmxax,nrmx, lmxax
    use m_lattic,only: plat=>lat_plat, qlat=>lat_qlat,bas=>rv_a_opos
    use m_supot,only: n1,n2,n3, gmax=>lat_gmax
    use m_rdsigm2,only: getsenex, senex,dsene
    use m_mkpot,only: smpot=>osmpot, vconst !, vesrmt
    use m_mkpot,only: osig, otau, oppi, ohsozz,ohsopm, oppix,spotx
    use m_MPItk,only: numproc=>nsize,procid,master,master_mpi,comm
    use m_igv2x,only: napw,ndimh,ndimhx,igv2x,m_Igv2x_setiq,ndimhall,nbandmx
    use m_elocp,only: rsmlss=>rsml, ehlss=>ehl
    use m_qplist,only: qplist,ngplist,ngvecp,iqibzmax,niqisp,iqproc,isproc
    use m_hamindex0,only: Readhamindex0, nlindx,nclass,iclass=>iclasst,lindx,nindx,nphimx
    use m_density,only: v0pot,pnuall,pnzall
    use m_augmbl,only: aughsoc
    use m_makusq,only: makusq
    use m_pwmat,only: pwmat
    use m_ftox
    use m_zhev,only: zhev_tk4
    use m_hambl,only: hambl
    use m_rdata1,only:rdata1init,nradmx,nnc,nrad,nindx_r,lindx_r,iord,nvmax,nrc,mindx,&
         gcore_n,aac,bbc,gval_orth,zzpi,nrmxe=>nrmx ,gval_n
    use m_blas,only: zmm => zmm_h, m_op_T, m_op_C
    use m_ppj,only: m_ppj_init,ppj
    implicit none
    intent(in)::          socmatrix,eferm,vmag,qval
    !  qval: valence charge
    !  osig,otau,oppi  augmentation matrices, s_rv1
    !  senex: real space Sigma_vxc
    !  hamm: Hamiltonian
    !  ovlm: Overlapmatrix
    !  nbas  :size of basis
    !  smpot :smooth potential on uniform mesh
    !  vconst:constant to be added to potential
    !
    !  evl(1:ndimh,isp): eigenvalue
    !  cphi(:,:,isp):MT part of　eigenfunciton
    !  pwz: IPW part of eigenfunciton
    !  vxclda(1:ndimh): LDA/GGA exchange correlation
    !  nev: number of calculated eigenfuncitons
    !  gval(:,:,1:3) radial w.f. gval: phi, phidot phiz for 1:3
    !  evec: eigenvectors.
    !  vxc:  matrix elements of LDA XC potential
    !---------------
    !For SO=1,  nsp=2, nspc=2, nspx=1,   ndimhx=ndimh*2
    !For SO/=1,        nspc=1, nspx=nsp, ndimhx=ndimh (nsp=1 or 2)
    integer :: lchk=1,i,i1,i2,iat,ib,ibr,icore,ierr,idat,ifievec,ifvxcevec,ifiqg,iflband(2),ifqeigen,ifsyml,iix,iline, &
         im1,im2,ipb(nbas),ipqn,ipr,iprint,iq,is,isp,ispc,j,job,k1, &
         k2,k3,konf,l,ldim,loldpw, lsig,mx,mxint, ncore,nevl,nev,nglob,ngp,ngp_p, &
         ngpmx,nline,nlinemax,nlmax,nmx,nn1,nn2,nnn, kkk,mmm,n, &
         npqn,nqbz,nqnum,nqnumx,nqtot,nr,iqibz,imx,ifigwb,ifinormchk,ifigw1,ifildima,ifigwn,ifigwbhead,&
         ispSS,ispEE,ispx,iqbk=-999,icore2,icore2o,ic,nmcore,ifnlax,ifigwa,ifgeigm,ifcphim,ixx,ix,nn
    real(8):: rsml(n0), ehl(n0) ,eferm,qval, vmag,vnow, QpGcut_psi,QpGcut_cou,dum,xx(5),a,z,vshft, qp(3),qpos,q_p(3), epsovl,dq(3),qpx(3)
    real(8),allocatable:: rofi(:),rwgt(:) !, cphiw(:,:) 
    real(8),pointer:: pnu(:,:),pnz(:,:)
    integer,allocatable :: konft(:,:,:),iiyf(:),ibidx(:,:),nqq(:), m_indx(:),n_indx(:),l_indx(:),ibas_indx(:)
    complex(8),allocatable :: aus_zv(:,:,:,:,:), hamm(:,:,:,:),ovlm(:,:,:,:),ovlmtoi(:,:),ovliovl(:,:) ,hammhso(:,:,:)
    complex(8),allocatable:: evec(:,:),evec0(:,:),vxc(:,:,:,:),ppovl(:,:),phovl(:,:),pwh(:,:),pwz(:,:,:),pzovl(:,:,:), pwz0(:,:),&
         testcc(:,:),testc(:,:,:),testcd(:,:),ppovld(:),cphi(:,:,:),cphi0(:,:,:),cphi_p(:,:,:),geig(:,:,:),geig_p(:,:,:),sene(:,:),ppovli(:,:)
    logical :: lwvxc,cmdopt0, emptyrun, magexist, debug=.false.,sigmamode,wanatom=.false.,once=.true.
    logical,optional:: socmatrix 
    character(8) :: xt
    character(256):: ext,sprocid,extn
    complex(8),allocatable::  geigr(:,:,:), cphix(:,:,:)
    integer:: mrecb,mrece,mrecg,ndble,ifv,iqq,ifev,konfigk,konfz
    real(8),allocatable::evl(:,:,:),vxclda(:,:,:) !,evl(:,:),vxclda(:)!qirr(:,:), !    complex(8),allocatable:: ppj(:,:,:)
    integer :: istat,ifoc,ldim2,lxx,nl,nrx
    logical :: blas_mode = .true. !    real(8):: rmax(nclass)
    include "mpif.h"
    call tcn ('m_sugw_init')
    debug=cmdopt0('--debugsugw')
    call getpr(ipr)
    if(lmxax<=0) call rx('sugw: lmxax>0 for gw mode')
    emptyrun  = cmdopt0('--emptyrun')
    sigmamode = mod(lrsig,10) .ne. 0
    magexist  = abs(vmag)>1d-6
    if(master_mpi) write(stdo,"(a)") 'm_sugw_init: start'
    if(master_mpi) write(stdo,"(' MagField added to Hailtonian -vmag/2 for isp=1, +vmag/2 for isp=2: vmag(Ry)=',d13.6)") vmag
    call Readhamindex0() ! ==== Read file NLAindx ====
    if(cmdopt0('--wanatom').and.master_mpi) wanatom=.true. 
    if(wanatom) then ! 'wanplotatom.dat' is originally a part of gwa and gwb.head. only for wanplot which will be unsupported.
      open(newunit=ifigwa,file='wanplotatom.dat',form='unformatted') 
      write(ifigwa)nbas,nsp,ndima,nbandmx,maxval(lmxa),ncoremx,nrmx,plat,alat!,nqirr,nqibz
      write(ifigwa)bas,lmxa(ispec(1:nbas)) !,ndimaa(1:nbas) !,qplist,ngplist,ndimhall,qval
    endif
    if(master_mpi) write(stdo,ftox)' Generate radial wave functions ncoremx,nphimx=',ncoremx,nphimx
    allocate(gval(nrmx,0:lmxax,nphimx,nsp,nclass), ecore(ncoremx,nsp,nclass),gcore(nrmx,ncoremx,nsp,nclass),source=0d0)
    ibmain: do 1150 ib = 1, nbas
      if(master_mpi) write(stdo,ftox)' ibas=',ib
      is=ispec(ib)  !spec index
      ic=iclass(ib) !class index
      pnu=>pnuall(:,1:nsp,ib)
      pnz=>pnzall(:,1:nsp,ib)
      a = spec_a(is)
      nr= nris(is)
      z = zz(is)
      nmcore= nmcorex(is)
      CreateAugmentedWaveFunctions: block
        real(8):: gvall(nr*2,0:lmxa(is),nphimx,nsp),gcorel(nr*2,ncoremx*nsp),ecorel(ncoremx*nsp)
        allocate(rofi(nr),rwgt(nr))
        call radmsh(rmt(is),a,nr,rofi)
        call radwgt(rmt(is),a,nr,rwgt)
        rsml= rsmlss(:,is)
        ehl = ehlss(:,is)
        RadialWaveFunctions: block 
          use m_rhocor,only: getcor
          use m_atwf,only: makrwf,wf2lo
          real(8) :: rhoc(nr,2),gp(2*nr*4),e, phi,dphi,phip,dphip,p,phz,dphz,phzp,dphzp, sumtc,sumec,ez
          do  l = 0, lmxa(is)
            do  isp = 1, nsp
              call makrwf(z,rofi(nr),l,v0pot(ib)%v(1,isp),a,nr,rofi,pnuall(1,isp,ib),4,gvall(1,l,1,isp),  gp,e,phi,dphi,phip,dphip,p)
              gvall(:,l,2,isp)=gp(1:2*nr) 
              if (floor(pnuall(l+1,1,ib)) /= konfig(l,ib)) then
                call makrwf(z,rofi(nr),l,v0pot(ib)%v(1,isp),a,nr,rofi,pnzall(1,isp,ib),2,gvall(1,l,3,isp),gp,ez,phz,dphz,phzp,dphzp,p)
                call wf2lo(l,a,nr,rofi,rwgt,phi,dphi,phip,dphip,phz,dphz, &
                     phzp,dphzp,pnzall(1,isp,ib),rsml,ehl, gvall(1,l,1,isp),gvall(1,l,2,isp),gvall(1,l,3,isp))
              endif
            enddo
          enddo
          call getcor(1,z,a,pnu,pnz,nr,lmxa(is),rofi,v0pot(ib)%v,0,0,[0d0,0d0],sumec,sumtc,rhoc,ncore,ecorel,gcorel,nmcore) !nmcore is non-magnetic core
        endblock RadialWaveFunctions
        if(wanatom) then
          write(ifigwa) z, nr, a, rmt(is)/(dexp(a*nr-a)-1d0), rmt(is), lmxa(is), nsp, ncore,slabl(ib) !b=rmt(is)/(dexp(a*nr-a)-1d0)
          write(ifigwa) konfig(0:lmxa(is),ib)
          write(ifigwa) rofi
          do  l = 0, lmxa(is)
            do  i = 1, nsp
              write(ifigwa) l,i
              write(ifigwa) gvall(1:nr,l,1,i)
              write(ifigwa) gvall(1:nr,l,2,i)
              if (konfig(l,ib) >= 10 .AND. master_mpi) write(ifigwa) gvall(1:nr,l,3,i)
            enddo
          enddo
          do  l = 0, lmxa(is)
            do  isp = 1, nsp
              do  konf = l+1, mod(konfig(l,ib),10)-1
                icore = icore+1
                write(ifigwa) icore, l, isp, konf, ecorel(icore)
                write(ifigwa) gcorel(1:nr,icore)
              enddo
            enddo
          enddo
        endif
        do  l = 0, lmxa(is) !!  Write orthonormalized valence wave functions for this atom
          do  i = 1, nsp
            gval(1:nr,l,1,i,ic)=  gvall(1:nr,l,1,i)
            gval(1:nr,l,2,i,ic)=  gvall(1:nr,l,2,i)
            if (konfig(l,ib) >= 10) gval(1:nr,l,3,i,ic)=gvall(1:nr,l,3,i)
          enddo
        enddo
        icore = 0 !  Core wave functions for this atom          !vshft = vesrmt(ib) !         vshft = 0
        do  isp = 1, nsp
          icore2=0
          do  l = 0, lmxa(is)
            konfigk = floor(pnuall(l+1,1,ib))           !take isp=1
            konfz   = floor(mod(pnzall(l+1,1,ib),10d0))
            if(konfz == 0) konfz = konfigk
            do  konf = l+1, min(konfz,konfigk)-1 !mod(konfig(l,ib),10)-1
              icore2= icore2+1
              icore = icore+1
              ecore     (icore2,isp,ic) = ecorel(icore)
              gcore(1:nr,icore2,isp,ic) = gcorel(1:nr,icore)
            enddo
          enddo
        enddo
      endblock CreateAugmentedWaveFunctions
      deallocate(rofi,rwgt)
!      if(ib==nbas) stop 'xxxxxxxxxxxxxx222'
1150 enddo ibmain
    
    if(wanatom) close(ifigwa)
    call rdata1init(ncores,ndima,ncoremx,konf0,gval,gcore) ! Write refined mesh and indexes to m_rdata1
    ! IPW part  Main loop for eigenfunction generation ==
    open(newunit=ifiqg,file='QGpsi',form='unformatted')
    read(ifiqg ) nqnum, ngpmx ,QpGcut_psi,nqbz,nqirr,imx,nqibz
    close(ifiqg)
!    if (lchk>=1 ) then !normcheck file
!      open(newunit=ifinormchk,file='norm.'//'procid'//trim(xt(procid))//'.chk')
!      write(ifinormchk,"(a)") '#         eval          IPW        IPW(diag)    Onsite(tot)      Total ! lmfgw'
!    endif
    if(ham_scaledsigma/=1d0 .AND. sigmamode) write(stdo,*)' Scaled Sigma method: ScaledSigma=',ham_scaledsigma
    ndble = 8
    mrecb = 2*ndima* nspc* nbandmx *ndble !byte size  !Use -assume byterecl for ifort, so that ifort recognizes the recored in the unit of bytes.
    mrece = nbandmx          *ndble 
    mrecg = 2*ngpmx*nbandmx  *ndble
    WriteGWfiles: if(master_mpi) then
      WriteGWfilesB: block
        integer,allocatable:: ncindx(:),lcindx(:)
        integer:: iorb,lx,nx,ifoc,ibas,ifec,irad,ifphi,ir,ibasf(nbas),ibasx,ifigwin,nnv,ifhbed
        logical:: laf
        ! @MNLA_core.chk index for core
        open(newunit=ifnlax,file='@MNLA_core.chk')
        write(ifnlax,"(a)") '    m    n    l  icore ibas   ' ! Index for core
        write(ifnlax,'(" ------- core ------------")')
        do ibas = 1,nbas
          is   = ispec(ibas)
          ic   = iclass(ibas)
          icore = 0             
          do l  = 0, lmxa(is)
            do kkk = l+1,konf0(l,ibas)-1 !kkk is the quantum principle number
              icore = icore+1
              n  = kkk - l   ! n is starting from 1 for any l.  
              do mmm=-l,l
                write(ifnlax,'(10i5)') mmm, n, l, icore, ibas !MNLA index. magnetic radial l numcore, ibas
              enddo
            enddo
          enddo
        enddo
        close(ifnlax)
        ! @MNLA_CPHI index for cphi 
        open(newunit=ifoc,file='@MNLA_CPHI')
        write(ifoc,"('    m    n    l ibas')")
        iorb=0
        do ibas = 1,nbas
          ic = iclass(ibas)
          is = ispec(ibas)
          do lx = 0,lmxa(is)
            do nx = 1,nvmax(lx,ic)
              iorb=iorb+1
              do mx = -lx,lx
                write(ifoc,"(10i6)")mx,nx,lx,ibas,iord(mx,nx,lx,ibas),iorb
              enddo
            enddo
          enddo
        enddo
        close(ifoc)
        write(stdo,ftox)" === Radial function indexing === nradmx=", nradmx
        do ibas=1,nbas
          write(stdo,ftox)' ---- ibas nrad(ibas) =', ibas, nrad(ibas)
          do irad = 1,nrad(ibas)
            write(stdo,'("      irad=",i3," nindx_r lindx_r=",2i3)')irad, nindx_r(irad,ibas), lindx_r(irad,ibas)
          enddo
        enddo
        !ECORE
        write(stdo,ftox)" === Write ECORE ==="
        open(newunit=ifec, file='ECORE')
        ibasloopc: do ibas = 1,nbas
          ic    = iclass(ibas)
          is    = ispec(ibas)
          write(ifec,*)            !ECORE
          write(ifec,*) slabl(is) !spid(ibas) !ECORE
          write(ifec,*) ' z,atom=class,nr,a,b,nsp ' !ECORE
          write(ifec,"(1x,f5.1,2i10,f13.5,d14.6,i4)") zz(is),ibas,nrc(ic),aac(ic),bbc(ic),nsp !ECORE
          write(ifec,*)' configuration'!   !!! LocalOrbital 2=upper 1=lower' !ECORE
          write(ifec,ftox)(konf0(l,ibas),l=0,lmxa(is)) !principl quantum  number of valence minimum
          write(ifec,*)' l,n, ecore(up), ecore(down) ' !ECORE ! related to LocalOrbital part lower(=1) upper(=2).
          icore = 0
          do l  = 0,lmxa(is)
            do kkk = l+1 ,konf0(l,ibas)-1
              icore = icore+1
              n    = kkk - l
              write(ifec,ftox) l,n,ftod(ecore(icore,1:nsp,ic),16) !,ftod(ec(icore,ic,1:nsp),16) !ECORE 
            enddo
          enddo
        enddo ibasloopc
        close(ifec)
        !PHIVC  
        write(stdo,ftox)" === Write PHIVC ==="
        open(newunit=ifphi,file='PHIVC',form='unformatted')
        write(ifphi) nbas, nradmx,ncoremx,nrmxe !extented for nrc
        write(ifphi) nrad(1:nbas)
        write(ifphi) nindx_r(1:nradmx,1:nbas),lindx_r(1:nradmx,1:nbas)
        allocate(ncindx(ncoremx),lcindx(ncoremx),source=-9999) !        allocate(nc_max(0:lmxax,nbas))
        ibasloopw: do ibas = 1,nbas
          ic    = iclass(ibas)
          is    = ispec(ibas)
          allocate(rofi(nrc(is)))
          rofi = [(bbc(ic)*(exp((ir-1)*aac(ic))-1d0), ir=1,nrc(ic))]
          icore = 0
          do l  = 0,lmxa(is) !lmxax !nl-1 !lmxa(is)
            do kkk = l+1 ,konf0(l,ibas)-1
              icore = icore+1
              n    = kkk - l
              ncindx(icore)= n
              lcindx(icore)= l
            enddo
          enddo
          write(ifphi) ncores(is), ncoremx !core
          write(ifphi) ncindx,lcindx !core !          if(ncindx>nc_max(l,ibas)) nc_max(l,ibas)=ncindx !nc_max
          write(ifphi) ibas,zz(is),nrc(ic),aac(ic),bbc(ic)
          write(ifphi) rofi(1:nrc(ic))
          do isp = 1, nsp 
            do icore = 1, ncores(is)
              write(ifphi) gcore_n(1:nrc(ic),icore, isp,ic) ! core
            enddo
            do irad = 1,nrad(ibas)
              l = lindx_r(irad,ibas)
              n = nindx_r(irad,ibas)
              write(ifphi) gval_orth(1:nrc(ic),l, n, isp,ic)  ! valence orthogonalized
              write(ifphi) gval_n (1:nrc(ic),l, n, isp,ic)  ! valence raw
            enddo
          enddo
          deallocate(rofi)
        enddo ibasloopw
        close(ifphi)
        ! MTOindex file. basic part of crystal structure.
        write(stdo,ftox)" === Write MTOindex file(crystal structure info and so on) ==="
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
        laf= sum(abs(iantiferro))/=0
        nnv = maxval(nindx(1:ndima))
        write(stdo,ftox)' iantiferro=',iantiferro(1:nbas)
        open(newunit=ifigwin,file='MTOindex',form='unformatted')    
        write(ifigwin) nbas,alat,plat,nsp,lmxax+1,nnv,nnc,nrmxe,qval,nspc,nlmto
        write(ifigwin) pos,zz(ispec(1:nbas)),slabl(ispec(1:nbas)),lmxa(ispec(1:nbas))
        write(ifigwin) ndble,mrecb,mrece,ndima,nqbz,nbandmx,mrecg
        write(ifigwin) laf,ibasf 
        close(ifigwin)
        open(newunit=ifhbed,file='hbe.d.chk') !human check only 
        write(stdo,'( " ndima nbandmx=",3i5)') ndima, nbandmx
        write(ifhbed,"('hbe output=',*(g0,x))") ndble,mrecb,mrece,ndima,nqbz,nbandmx,mrecg,nspc
        write(ifhbed,*)' precision, mrecl of b, mrecl of eval, ndima(p+d+l)  nqbz  nbandmx mrecg nspc'
        close(ifhbed)
      endblock WriteGWfilesB
    endif WriteGWfiles
    if(cmdopt0('--skipCPHI')) goto 1011
    ! ppj: ovalap matrix within MT. Gramschmidt2
    call mpi_barrier(comm,ierr)
    call m_ppj_init() 
    ! CPHI GEIG. We use mpi-io from 2024-9-26
    allocate(cphix(ndima,nspc,nbandmx),geigr(ngpmx,nspc,nbandmx))
    i=openm(newunit=ifcphim,file='CPHI',recl=mrecb)
    i=openm(newunit=ifgeigm,file='GEIG',recl=mrecg)
    allocate(evl(nbandmx, nqirr, nspx),vxclda(nbandmx, nqirr, nspx),source=0d0)!nqirr: # ofirreducible q points
    iqisploop: do 1001 idat=1,niqisp !iq = iqini,iqend ! iqini:iqend for this procid
      write(stdo,ftox) 'do 1001 idat=',idat
      iq  = iqproc(idat) ! iq index
      isp = isproc(idat) ! spin index: Note isp=1:nspx, where nspx=nsp/nspc.  isp=1 nspc=2 only for lso=1 if(debug)write(stdo,ftox)' iqisploop',iq,isp  
      qp  = qplist(:,iq) ! q vector containing nqirr
      ngp = ngplist(iq)  ! number of planewaves for PMT basis
      call m_igv2x_setiq(iq) ! Set napw ndimh ndimhx, and igv2x !note ndimhx is given here.
      allocate(hamm(ndimh,nspc,ndimh,nspc),ovlm(ndimh,nspc,ndimh,nspc)) !Spin-offdiagonal block included since nspc=2 for lso=1.
      allocate(evec(ndimhx,ndimhx),vxc(ndimh,nspc,ndimh,nspc),cphi(ndima,ndimhx,nspc))!,cphiw(ndimhx,nspc))
      if(iqbk==iq) then
        continue
      elseif( lso/=0 .OR. socmatrix) then
        if(allocated(hammhso)) deallocate(hammhso)
        allocate(hammhso(ndimh,ndimh,3))
        call aughsoc(qp, ohsozz,ohsopm, ndimh, hammhso)
        iqbk=iq !q index for hammhso
      endif
      if(emptyrun) then !set dummy to avoid error exit
        evec=1d0 
        evl(:,iq,isp)=[(i*0.1,i=1,ndimhx)]
        nev=ndimh
        goto 1212
      endif;       if(debug)write(stdo,ftox)' iqisploop333'
      GetHamiltonianAndDiagonalize: block
        integer:: iprint
        if(lso==1) then !L.S case nspc=2
          vxc=0d0 !zeroclear offdiagonal parts
          ovlm=0d0
          do ispc=1,nspc  ! nspc==2 ! LDA part of Hamiltonian and overlap matrices for this qp ---
            call hambl(ispc,qp,spotx,vconst,osig,otau,oppix, vxc(:,ispc,:,ispc),ovlm(:,ispc,:,ispc))
            call hambl(ispc,qp,smpot,vconst,osig,otau,oppi, hamm(:,ispc,:,ispc),ovlm(:,ispc,:,ispc))
            vxc(:,ispc,:,ispc) = hamm(:,ispc,:,ispc) - vxc(:,ispc,:,ispc) ! vxc(LDA) part
            hamm(:,ispc,:,ispc)= hamm(:,ispc,:,ispc) + hammhso(:,:,ispc) !spin-diag SOC elements (1,1), (2,2) added
          enddo
          hamm(:,1,:,2)= hammhso(:,:,3)                    !spin-offdiagonal SOC elements (1,2) block added
          hamm(:,2,:,1)= transpose(dconjg(hammhso(:,:,3))) !                              (2,1) block
          if(sigmamode) then !Add  Vxc(QSGW)-Vxc 
            do ispc=1,nspc
              call getsenex(qp, ispc, ndimh, ovlm(:,ispc,:,ispc)) !bugfix at 2024-4-24 obata: ispc was 1 when 2023-9-20
              hamm(:,ispc,:,ispc) = hamm(:,ispc,:,ispc) + ham_scaledsigma*senex !sene= Vxc(QSGW)-Vxc(LDA)
              call dsene()
            enddo
          endif
        else ! lso=0 (No SO) or lso=2(Lz.Sz)  Spin Diagonal case.spin diagonal) nspc=1 only
          call hambl(isp,qp,spotx,vconst,osig,otau,oppix, vxc(:,1,:,1),  ovlm(:,1,:,1)) !vxc=<F_i|H(LDA)-vxc(LDA)|F_j>
          call hambl(isp,qp,smpot,vconst,osig,otau,oppi,  hamm(:,1,:,1), ovlm(:,1,:,1)) !ham=<F_i|H(LDA)|F_j> and ovl=<F_i|F_j>
          if(lso==2) hamm(:,1,:,1) = hamm(:,1,:,1) + hammhso(:,:,isp) !diagonal part of SOC matrix added for Lz.Sz mode.
          vxc(:,1,:,1) = hamm(:,1,:,1) - vxc(:,1,:,1) ! vxc(LDA) part
          if(sigmamode) then !Add  Vxc(QSGW)-Vxc 
            call getsenex(qp,isp,ndimh,ovlm(:,1,:,1))
            hamm(:,1,:, 1) = hamm(:,1,:,1) + ham_scaledsigma*senex !senex= Vxc(QSGW)-Vxc(LDA)
            call dsene()
          endif
        endif;    if(debug)write(stdo,ftox)'sumcheck hamm=',sum(abs(hamm)),sum(abs(ovlm))
        if (mod(iq,10) /= 1) call pshpr(iprint()-6)  
        epsovl = ham_oveps
        evec=-1d99
        evl(:,iq,isp)=1d99
        AddExternelMagneticField: if(magexist) then !
          if(nspc==2) then
            hamm(:,1,:,1)= hamm(:,1,:,1) - vmag/2d0*ovlm(:,1,:,1)
            hamm(:,2,:,2)= hamm(:,2,:,2) + vmag/2d0*ovlm(:,2,:,2)
          else
            if(isp==1) hamm(:,1,:,1)= hamm(:,1,:,1) - vmag/2d0*ovlm(:,1,:,1)
            if(isp==2) hamm(:,1,:,1)= hamm(:,1,:,1) + vmag/2d0*ovlm(:,1,:,1)
          endif
        endif AddExternelMagneticField;            if(debug)write(stdo,ftox)' iqisploop666'
        call zhev_tk4(ndimhx,hamm,ovlm,ndimhx,nev,evl(1,iq,isp),evec,epsovl) ! Diagonalization. nev:Calculated number of eigenvec
      endblock GetHamiltonianAndDiagonalize;       if(debug)write(stdo,ftox)' iqisploop777 1212'
1212  continue
      lwvxc = (socmatrix .or. iq<=iqibzmax).and.(.not.cmdopt0('--novxc'))
      if(lwvxc) then
        open(newunit=ifvxcevec, file= 'vxcevec'//trim(xt(iq))//trim(xt(isp)),form='unformatted')
        write(ifvxcevec) qp,ndimhx,nev
        write(ifvxcevec) vxc(:,1:nspc,:,1:nspc)
        write(ifvxcevec) evec(1:ndimhx,1:ndimhx),evl(1:ndimhx,iq,isp)
        if(lso/=0.or.socmatrix) write(ifvxcevec) hammhso 
        close(ifvxcevec)
      endif
      if(emptyrun) then
        allocate(pwz(ngp,nspc,ndimh)) !dummy
        goto 1214
      endif
      write(stdo,ftox)' sugw: kpt isp=',iq,isp,'of',nqnum,'k= ',ftof(qp,5),'ndimh=',ndimh,'irank=',procid,'lwvxc=',lwvxc,'nev=',nev
      write(stdo,"(9f8.4)",advance='no') (evl(i,iq,isp),i=1,min(18,nev))
      write(stdo,ftox)' ...'
      evl(1+nev:nbandmx,iq,isp)=1d20 !padding
      if(mod(iq,10) /= 1) call poppr;  if(debug) write(stdo,"(' sugw:procid iq isp lwvxc= ',3i3,' ',l)")procid, iq,isp,lwvxc
      nlmax = (lmxax+1)**2
      CPHIpart: block
        use m_locpot,only: rotp !        use m_mkpot,only : sab_rv
        integer:: ilm,im,iv
        complex(8):: auasaz(3),aus_zv(nlmax,nbandmx,3,nsp,nbas),usz(3)
        !i   nlmax :leading dimension of aus
        !i   nev   :number of eigenvectors to accumulate cphi
        !i   nlindx: offset that set index in cphi for element (ipqn,l,ib)
        !i   aus   : values of (phi,phidot,pz) at MT sphere boundary; see makusq
        !o Outputs
        !o   cphi(ichan,iv) : coefficients to phi,phidot,phiz 
        !o      ichan = 1:ndima, orbital channel, i.e. one of phi,phidot,phiz (phiz is local orbital for a given site, l, and m; see nlindx.
        !o      iv = eigenvector
        !o   cphiw: diagonal matrix elements, one for each eigenvector. only for check  : cphiw(1,iv) = <cphi(iv) | overlap | cphi(iv)>
        call makusq(nbas,[-999], nev,  isp, 1,qp,reshape(evec(1:ndimhx,1:nev),[ndimh,nspc,nev]), aus_zv ) !    cphiw=0d0
        ispcc: do ispc=1,nspc
          if(lso==1) ispx=ispc
          if(lso/=1) ispx=isp
          ibb: do ib = 1, nbas
            do iv = 1, nev
              ilm  = 0
              do l = 0, lmxa(ispec(ib))
                do im = 1, 2*l+1
                  ilm = ilm+1
                  auasaz=aus_zv(ilm,iv,1:3,ispx,ib) !coefficient for (u,s,gz)
                  ! cphi are coefficients for augmented functions {phi,phidot,gz(val=slo=0)}, which are not orthnormal.
                  ! TK checked that sab_rv is almost the same as zzpi*zzpi as for nio_gwsc
                  cphi(nlindx(1:2,l,ib)+im,iv,ispc)= matmul(auasaz(1:2),rotp(l,ispx,:,:,ib))
                  if (nlindx(3,l,ib) >= 0) cphi(nlindx(3,l,ib) + im,iv,ispc) = auasaz(3)
                  !NOTE: cphiw(iv,ispc) =cphiw(iv,ispc)+sum(dconjg(auasaz)*matmul(sab_rv(:,:,l+1,ispx,ib),auasaz)) ! \approx 1d0 c.f. ovv below
                enddo
              enddo
            enddo
          enddo ibb
        enddo ispcc
      endblock CPHIpart
      GEIGpart: if(ngp > 0) then !IPW expansion of eigenfunctions pwz 
        !  ppovl: = O_{G1,G2} = <IPW_G1 | IPW_G2>
        !  phovl: <IPW_G1 | basis function = smooth Hankel or APW >   
        !    pwz:  IPW expansion of eigen function    
        allocate(ppovl(ngp,ngp),pwz(ngp,nspc,ndimhx),phovl(ngp,ndimh),ppovli(ngp,ngp))
        call pwmat(nbas,ndimh,napw,igv2x,qp,ngp,nlmax,ngvecp(1,1,iq),gmax, ppovl, phovl )
        ! MO added blas_mode to replaced matmul by a BLAS call 2024-11-09. This is because matmul in mic(intel) was very slow.
        if(blas_mode) then
          do ispc=1, nspc
            istat = zmm(phovl, evec(ndimh*(ispc-1)+1,1), pwz(1,ispc,1), m=ngp, n=ndimhx, k=ndimh, ldb=ndimh*nspc, ldc=ngp*nspc)
          enddo
        else
          pwz(1:ngp,1,1:ndimhx) = matmul(phovl(1:ngp,1:ndimh),evec(1:ndimh,1:ndimhx))
          if(lso==1) pwz(1:ngp,2,1:ndimhx) = matmul(phovl(1:ngp,1:ndimh),evec(ndimh+1:2*ndimh,1:ndimhx))
        endif
        deallocate(phovl)
        !        if (lchk >= 1) then   
        !          allocate(pzovl,source=pwz)
        !          allocate(ppovld(ngp)) ! extract diagonal before ppovl overwritten
        !          forall(i = 1:ngp) ppovld(i) = ppovl(i,i)
        !        endif
        ppovli=ppovl
        call matcinv(ngp,ppovli)! inversion of hermitian ppovl
        ! MO added blas_mode to replaced matmul by a BLAS call 2024-11-09.
        if(blas_mode) then
          block
          complex(8) :: ppovl_pwz(ngp,ndimhx)
          do ispc=1, nspc
            istat = zmm(ppovli, pwz(1,ispc,1), ppovl_pwz, m=ngp, n=ndimhx, k=ngp, ldb=ngp*nspc)
            pwz(1:ngp,ispc,1:ndimhx)=ppovl_pwz(1:ngp,1:ndimhx)
          enddo
          endblock
        else
          pwz(1:ngp,1,1:ndimhx) = matmul(ppovli,pwz(1:ngp,1,1:ndimhx)) !pwz= O^-1 *phovl * evec  ! IPW expansion of eigenfunction
          if(lso==1) pwz(1:ngp,2,1:ndimhx) = matmul(ppovli,pwz(1:ngp,2,1:ndimhx))
        endif
        deallocate(ppovli)
      !   if (lchk >= 1) then
      !     allocate(testc(ndimhx,ndimhx,nspc),testcd(ndimhx,nspc))
      !     do ispc=1,nspc
      !       if(lso/=1) ispx=isp
      !       if(lso==1) ispx=ispc
      !       ! MO added blas_mode to replaced matmul by a BLAS call 2024-11-10.
      !       if(blas_mode) then
      !         istat = zmm(pzovl(1,ispc,1), pwz(1,ispc,1), testc(1,1,ispc), m=ndimhx, n=ndimhx, k=ngp, &
      !                     opA=m_op_C, ldA=ngp*nspc, ldB=ngp*nspc)
      !         ! testc(:,:,ispc)= matmul(transpose(dconjg(pzovl(:,ispc,:))),pwz(1:ngp,ispc,1:ndimhx))
      !         testcd(:,ispc) = [(sum(dconjg(pwz(:,ispc,i))*ppovld(:)*pwz(:,ispc,i)),i=1,ndimhx)]
      !       else
      !       associate(pwz1=>pwz(1:ngp,ispc,1:ndimhx),pzo=>dconjg(pzovl(:,ispc,:)))
      !         testc(:,:,ispc)= matmul(transpose(pzo),pwz1)
      !         ! testc(:,:,ispc)= matmul(transpose(dconjg(pzo)),pwz1)
      !         testcd(:,ispc) = [(sum(dconjg(pwz1(:,i))*ppovld*pwz1(:,i)),i=1,ndimhx)] !dimhx)]
      !       end associate
      !       endif
      !     enddo
      !     write(ifinormchk,"('# iq',i5,'   q',3f12.6,' ndimhx nev',2i7)") iq,qp,ndimhx,nev
      !     ! xx(1) = sum over all augmentation w.f.  cphi+ ovl cphi
      !     ! xx(3) = IPW contribution to phi+ phi.   
      !     ! xx(4) = IPW contribution to phi+ phi, using diagonal part only]
      !     do i1 = 1, nev !nev: number of calculated eigenvals(eigenfunctions).
      !       i2=i1
      !       xx(1) = sum(cphiw(i1,1:nspc))
      !       xx(3) = sum(testc(i1,i2,:))
      !       xx(4) = sum(testcd(i1,:))    !if(i1==i2)write(ifinormchk,'(f12.5,5f14.6)')evl(i1,isp),xx(3),xx(4),xx(1),xx(1)+xx(3)
      !       write(ifinormchk,'(i4,f12.5,5f14.6)')i1,evl(i1,iq,isp),xx(3),xx(4),xx(1),xx(1)+xx(3) !xx(1)+xx(3) should be close to unity.
      !     enddo
      !     deallocate(testc,testcd)
      !     write(ifinormchk,*)
      !     deallocate(pzovl)
      !     deallocate(ppovld)
      !   endif
      endif GEIGpart;    if(debug) write (stdo,"('q ndimh=',3f10.5,i10)") qp, ndimh
      VXCmat: block
        allocate(testcc(1:nev,ndimhx),source=matmul(transpose(dconjg(evec(:,1:nev))),reshape(vxc,[ndimhx,ndimhx])))
        forall(i1 = 1:nev) vxclda(i1,iq,isp) = sum(testcc(i1,1:ndimhx) * evec(1:ndimhx,i1))  !<i|Vxc^lda|i>
        vxclda(nev+1:,iq,isp)=0d0 !padding
        deallocate(testcc)
      endblock VXCmat
1214  continue;  if(debug)write(stdo,ftox)'1214 goto writechpigeig'   
      WriteCphiGeig: block
        use m_readqg,only: readngmx,ngcmx,readqg0,readqg
        use m_hamindex0,only: nindx,ibasindx !    use m_mkpot,only: sab_rv
        use m_locpot,only: rotp
        use m_pwmat,only: mkppovl2
        real(8)::add,zzz(2,2),epscheck=1d-8
        complex(8)::ccc(3),nnn,rrr(3,3),mmm(3,3),ovv
        integer::iband,ibas,iqqisp,ix,m,nm,i,ilm,im,iv,ngvecpf1(3,ngp),ndg1(3)
        do ispc=1,nspc
          geigr(1:ngp,      ispc,1:ndimhx)=pwz(1:ngp,ispc,1:ndimhx)
          geigr(ngp+1:ngpmx,ispc,1:ndimhx)=0d0
        enddo ! skip cphi(ix,1:nev,1:nspc) = cphi(ix, 1:nev,1:nspc) /sqrt(1d0+0.1d0*nindx(ix)) here because zzpi includes this factor 2025-5-7
        cphix=0d0 !   Augmentation wave part. cphix is coefficients for the orthogonalized functions gval_ortrh
        if(debug)write(stdo,ftox)' writechpigeig 1111'
        do ispc=1,nspc
          ispx = merge(ispc,isp,lso==1)
          do iband = 1,nev
            do ix= 1,ndima
              ib= ibasindx(ix); l= lindx(ix); n= nindx(ix); m= mindx(ix)
              ic= iclass(ib)
              nm= nvmax(l,ic)
              cphix      (iord(m,1:nm,l,ib),ispc,iband) = & !chipx is for coeffieients of orthogonalized atomic funcitons gval_orth*Ylm
                   cphix (iord(m,1:nm,l,ib),ispc,iband) + zzpi(1:nm,n,l,ic,ispx)*cphi(ix,iband,ispc) 
            enddo
          enddo
        enddo;        if(debug)write(stdo,ftox)' writechpigeig 2222'
        GramSchmidtCphiGeig :block
          if(.not.cmdopt0('--skipGS')) &
               call GramSchmidt2(nspc,nev,ndima,ngp,ngpmx, ppj(1:ndima,1:ndima,isp),ppovl, cphix,geigr) !Improve Orthogonalization
          ncheckw: do ispc=1,nspc !Normalization check of MT+IPW division of eigenfunctions 
            do   i=1,nev
              do j=1,nev
                ovv= sum( dconjg(cphix(1:ndima,ispc,i))*matmul(ppj(:,:,isp), cphix(1:ndima,ispc,j))) + & !MT parts
                     sum( dconjg(geigr(1:ngp,  ispc,i))*matmul(ppovl,      geigr(1:ngp,  ispc,j)))     !IPW parts
                if(i/=j.and.abs(ovv)    >epscheck) write(stdo,ftox)'oooovlap=',i,j,ispc,ftod(abs(ovv))
                if(i==j.and.abs(ovv-1d0)>epscheck) write(stdo,ftox)'oooovlap=',i,j,ispc,ftod(abs(ovv))
              enddo
            enddo
          enddo ncheckw
          deallocate(ppovl)
        endblock GramSchmidtCphiGeig
        cphix(1:ndima,1:nspc,nev+1:nbandmx)=1d20 !padding 
        iqqisp= isp + nsp*(iq-1)
        i=writem(ifcphim,rec=iqqisp,data=cphix(1:ndima,1:nspc,1:nbandmx)) 
        if(ngpmx/=0) geigr(1:ngpmx,1:nspc,nev+1:nbandmx)=1d20   ! padding
        if(ngpmx/=0) i=writem(ifgeigm,rec=iqqisp,data=geigr(1:ngpmx,1:nspc,1:nbandmx));  if(debug)write(stdo,ftox)'end of writechpigeig'
      endblock WriteCphiGeig; if(debug)write(stdo,ftox)' writechpigeig 1001'  
      deallocate(pwz,hamm,ovlm,evec,vxc,cphi)!,cphiw)
1001 enddo iqisploop
    i=closem(ifcphim) !mpi-io
    i=closem(ifgeigm)
    call mpi_barrier(comm,ierr)
    call mpibc2_real(evl,   nbandmx*nqirr*nspx,'evl')
    call mpibc2_real(vxclda,nbandmx*nqirr*nspx,'vxclda')
    WriteGWfiles2: if(master_mpi) then
      open(newunit=ifv,file='VXCFP',form='unformatted')
      write(ifv) nbandmx,nqirr
      do iqq = 1,nqirr 
         write(ifv) qplist(1:3,iqq), vxclda(1:nbandmx,iqq,1:nspx) ! VXCFP
      enddo
      open(newunit=ifev,file='EValue',form='unformatted')
      write(ifev) nbandmx, nqirr, nsp,nspc
      write(ifev) qplist(1:3,1:nqirr) !qirr
      write(ifev) evl(1:nbandmx, 1:nqirr, 1:nspx )
      close(ifev)
    endif WriteGWfiles2
    if(master_mpi) then
      rdata4gwblock: block
        use m_nvfortran,only:findloc
        use m_read_bzdata,only: Read_bzdata, nqibz,qibz, nq0i,nq0iadd,q0i,iq0pin
        use m_pwmat,only: mkppovl2
        use m_qplist,only: qirr=>qplist
        real(8),parameter:: pi = 4d0*datan(1d0)
        integer:: ifhvccfp,i,ngp,ngc,iq,iqi,irr,ix, idummy11(1,1),ippovlg,ippovli,ippovlgg
        integer:: nggg,ngcgp, ifiqg,ifiqgc,irrq, nqtt, nqnum,ngpmx,nqnumc,ngcmx!,nqbz
        integer:: ippovl0, nqnumt,iqx,nqini,iqtt
        real(8):: qpgcut_psi2,qx(3),dQpG,dQQ,QpGcut_psi,QpGcut_cou,qxx(3), QpGcutggg,QpGcutgcgp, tolq=1d-8
        logical:: ppovl0l=.true.
        character(3) :: charnum3
        real(8),allocatable ::qibze(:,:),qsave(:,:),qtt(:,:),rmax(:) 
        integer,allocatable:: nvggg(:,:),nvgcgp(:,:), ngveccB(:,:)
        integer,allocatable,target:: ngvecptt(:,:,:),ngvecctt(:,:,:),ngptt(:),ngctt(:),iqindex(:)
        complex(8),allocatable :: ppovl(:,:),ppx(:,:),ppovlinv(:,:),ggg(:)
        integer,pointer:: ngvecp(:,:),ngvecc(:,:)
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
          if(iq<4.or.iq>nqtt-3) write(stdo,ftox)' qtt =',iq,ftof(qtt(1:3,iq),5)
          if(iq==4) write(stdo,ftox)' ...'
        enddo
        close(ifiqg)
        close(ifiqgc)
        call read_bzdata()
        write(stdo,ftox)' QpGcut_psi QpGcutCou =',ftof(QpGcut_psi),ftof(QpGcut_Cou),'nqirr nbas',nqirr,nbas
        allocate(rmax,source=rmt(ispec(1:nbas))) ; if(debug) write(stdo,'("  qibz=",i3,3f12.5)')(i,qibz(1:3,i),i=1,nqibz)
        dQpG=maxval(sum(q0i(1:3,1:nq0i+nq0iadd)**2,dim=1))**.5
        allocate(qibze(3,nqibz + nq0i+nq0iadd)) !+iadd)) !nq0i+1 is for gamma point for bzcase==2
        qibze(1:3,1:nqibz)=qibz(1:3,1:nqibz)
        qibze(1:3,nqibz+1:nqibz+nq0i+nq0iadd)=q0i(1:3,1:nq0i+nq0iadd)
        if(debug) then
          do iqi =1,nqibz + nq0i+nq0iadd
            write(stdo,"(' iqi qibze=',i4,3f9.5)")iqi,qibze(:,iqi)
          enddo
        endif
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
        call getgv2( alat,plat,qlat,qx, QpGcutgcgp,1,ngcgp, idummy11) ;if(debug) write(stdo,"(' iqi qx nggg ngcgp=',i5,3f8.3,2i8)") iqi,qx,nggg,ngcgp
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
        if(debug) then
          do ix=1,nqirr
            write(stdo,"(' qirr= ',3f10.4)") qirr(1:3,ix)
          enddo
        endif
        !Write PPOVLGG === make <Gc-Gp1+Gp2> matrix ===
        write(stdo,ftox)' === Write PPOVLGG PPOVLG PPOVLI ==='
        write(stdo,ftox)' nggg ngcgp=',nggg,ngcgp,'nqnumt nqini nqtt',nqnumt,nqini,nqtt
        open(newunit=ippovlgg,file= "PPOVLGG",form='unformatted')
        write(ippovlgg) nggg, ngcgp, nqnumt-nqini+1,nqini,nqnumt
        write(ippovlgg) nvgcgp(1:3,1:ngcgp)
        write(ippovlgg) nvggg(1:3,1:nggg)
        write(ippovlgg) ggg(1:nggg)
        deallocate(ggg,nvggg,nvgcgp)
        close(ippovlgg)
        !Write PPOVL0,PPOVLG, PPOVLI 
        if(ppovl0l) open(newunit=ippovl0,form='unformatted',file='PPOVL0')
        iqiloop: do iqi = nqini, nqnumt    !nqibz + nq0i !+ iadd
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
          write(stdo,"(' iqi qx iqx=',i5,3f8.4,i5,' ngc ngp=',4i5)")iqi,qx,iqx,ngc,ngp !sum(abs(ngvecp)),sum(abs(ngvecc))
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
        enddo iqiloop
        if(ppovl0l) close(ippovl0)
        write(stdo,*)" end of rdata4gw "
      endblock rdata4gwblock
    endif
1011 continue !skipcphi
    call tcx('m_sugw_init')
  end subroutine m_sugw_init
end module m_sugw
