!> Generate all the inputs for GW calculation
module m_sugw
  use m_lgunit,only:stdo
  use m_lmfinit,only: zz=>z,nris=>nr,lmxa,rmt,spec_a
  real(8),allocatable,public::ecore(:,:,:),gcore(:,:,:,:),gval(:,:,:,:,:)
  integer,public::  ndima, ndham, ncoremx,nqirr,nqibz
  integer,allocatable,public:: konfig(:,:),ncores(:), konf0(:,:)
  private
  public:: m_sugw_init
contains
  subroutine m_sugw_init (socmatrix,eferm,vmag,qval) !Driver for GW calculation
    use m_ext,only:   sname
    use m_suham,only: ndham=>ham_ndham !max dimension of hamiltonian +napwad (for so=0,2)
    use m_lmfinit,only: ham_pwmode,pwemax,ham_oveps,lrsig=>ham_lsig,nlmto,lso,nspx !nspx=nsp/nspc
    use m_lmfinit,only: ham_scaledsigma,alat=>lat_alat,nsp,nspc,ispec,nspec
    use m_lmfinit,only: nbas,n0,nppn,nkap0,slabl,nmcorex=>nmcore,iantiferro,lmxax,nrmx
    use m_lattic,only: plat=>lat_plat, qlat=>lat_qlat,bas=>rv_a_opos
    use m_supot,only: n1,n2,n3, gmax=>lat_gmax
    use m_rdsigm2,only: getsenex, senex,dsene
    use m_mkpot,only: smpot=>osmpot, vconst !, vesrmt
    use m_mkpot,only: osig, otau, oppi, ohsozz,ohsopm, oppix,spotx
    use m_MPItk,only: numproc=>nsize,procid,master,master_mpi,comm
    use m_igv2x,only: napw,ndimh,ndimhx,igv2x,m_Igv2x_setiq,ndimhall
    use m_elocp,only: rsmlss=>rsml, ehlss=>ehl
    use m_qplist,only: qplist,ngplist,ngvecp,iqibzmax,niqisp,iqproc,isproc
    use m_hamindex0,only: Readhamindex0, nlindx,nclass,iclass=>iclasst,lindx
    use m_density,only: v0pot,pnuall,pnzall
    use m_augmbl,only: aughsoc
    use m_makusq,only: makusq
    use m_pwmat,only: pwmat
    use m_ftox
    use m_zhev,only: zhev_tk4
    use m_hambl,only: hambl
    use m_rdata1,only:rdata1init,nradmx,nnc,nrad,nindx_r,lindx_r,iord,nvmax,nrc,mindx,gval_n,gcore_n,aac,bbc,gx_orth,zzpi,rmax
    use m_suham,only: nbandmx=>ham_ndham
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
    !For SO=1,  nsp=2, nspc=2, nspx=1,   ndimhx=ndimh*2
    !For SO/=1,        nspc=1, nspx=nsp, ndimhx=ndimh (nsp=1 or 2)
    integer :: lchk=1,i,i1,i2,iat,ib,ibr,icore,ierr,idat,ifievec,ifiv,ifiqg,iflband(2),ifqeigen,ifsyml,igets,igetss,iix,iline, &
         im1,im2,ipb(nbas),ipqn,ipr,iprint,iq,is,isp,ispc,j,job,k1, &
         k2,k3,konf,l,ldim,loldpw, lsig,mx,mxint, ncore,nevl,nev,nglob,ngp,ngp_p, &
         ngpmx,nline,nlinemax,nlmax,nmx,nn1,nn2,nnn, kkk,mmm,n, &
         nphimx,npqn,nqbz,nqnum,nqnumx,nqtot,nr,iqibz,imx,ifigwb,ifinormchk,ifigw1,ifildima,ifigwn,ifigwbhead,&
         ispSS,ispEE,ispx,iqbk,konfigk,konfz,icore2,icore2o,ic,nmcore,ifnlax,ifigwa,ifcphi,ifgeig
    real(8):: rsml(n0), ehl(n0) ,eferm,qval, vmag,vnow, QpGcut_psi,QpGcut_cou,dum,xx(5),a,z,vshft, qp(3),qpos,q_p(3), epsovl
    real(8),allocatable:: rofi(:),rwgt(:), cphiw(:,:) 
    real(8),pointer:: pnu(:,:),pnz(:,:)
    integer,allocatable :: konft(:,:,:),iiyf(:),ibidx(:,:),nqq(:)
    complex(8),allocatable :: aus_zv(:,:,:,:,:), hamm(:,:,:,:),ovlm(:,:,:,:),ovlmtoi(:,:),ovliovl(:,:) ,hammhso(:,:,:)
    complex(8),allocatable:: evec(:,:),evec0(:,:),vxc(:,:),ppovl(:,:),phovl(:,:),pwh(:,:),pwz(:,:),pzovl(:,:), pwz0(:,:),&
         testc(:,:),testcd(:),ppovld(:),cphi(:,:,:),cphi0(:,:,:),cphi_p(:,:,:),geig(:,:,:),geig_p(:,:,:),sene(:,:)
    logical :: lwvxc,cmdopt0, emptyrun, magexist, debug=.false.,sigmamode,wanatom=.false.
    logical,optional:: socmatrix 
    character(8) :: xt
    character(256):: ext,sprocid,extn
    complex(8),allocatable::  geigr(:,:,:), cphix(:,:)
    integer:: mrecb,mrecg,ndble,ifv,iqq,ifev
    real(8),allocatable::qirr(:,:),evl(:,:,:),vxclda(:,:,:)!,evl(:,:),vxclda(:)
    include "mpif.h"
    call tcn ('m_sugw_init')
    call getpr(ipr)
    if(lmxax<=0) call rx('sugw: lmxax>0 for gw mode')
    emptyrun  = cmdopt0('--emptyrun')
    sigmamode = mod(lrsig,10) .ne. 0
    magexist  = abs(vmag)>1d-6
    lwvxc     = .not. cmdopt0('--novxc')
    if(master_mpi) write(stdo,"(a)") 'm_sugw_init: start'
    if(master_mpi) write(stdo,"('MagField added to Hailtonian -vmag/2 for isp=1, +vmag/2 for isp=2: vmag(Ry)=',d13.6)") vmag
    call Readhamindex0() ! ==== Read file NLAindx ====
! MT part -----------------------
    allocate(konfig(0:lmxax,nbas))
    nphimx=2
    allocate(ncores(nspec))
    do ib=1,nbas
       is=ispec(ib)
       ncore=0
       do l = 0, lmxa(is)
          konfigk = floor(pnuall(l+1,1,ib))           !take isp=1 since spin-independent
          konfz   = floor(mod(pnzall(l+1,1,ib),10d0))
          if(konfz == 0) konfz = konfigk
          ncore = ncore+min(konfz,konfigk)-1 -l
          if(konfz /= konfigk) nphimx=3 !pz mode local orbital is for nphi=3
          if (konfz < konfigk) then
             konfig(l,ib) = konfz + 10
          elseif (konfz > konfigk) then
             konfig(l,ib) = konfigk + 20
          else
             konfig(l,ib)= konfigk
          endif
       enddo
       ncores(is)=ncore
    enddo
    ncoremx=maxval(ncores)
    allocate(konf0(0:lmxax,nclass))
    do ib = 1,nbas
       is   = ispec(ib)
       ic   = iclass(ib)
       konf0(:,ic) = [(mod(konfig(l,ib),10),l=0,lmxa(is))]
    enddo
    ndima = 0
    do  ib = 1, nbas
       pnz=>pnzall(:,1:nsp,ib)
       if(sum(abs(pnz(:,1:nsp)-pnzall(:,1:nsp,ib)))>1d-9) call rx('sugw xxx1aaa')
       ndima = ndima + sum([((2*l+1)*merge(3,2,pnz(l+1,1)/=0), l=0,lmxa(ispec(ib)))])
    enddo
    ! 'wanplotatom.dat' is originally a part of gwa and gwb.head. only for wanplot which will be unsupported.
    if(cmdopt0('--wanatom').and.master_mpi) wanatom=.true.
    if(wanatom) then 
       open(newunit=ifigwa,file='wanplotatom.dat',form='unformatted') 
       write(ifigwa)nbas,nsp,ndima,ndham,maxval(lmxa),ncoremx,nrmx,plat,alat!,nqirr,nqibz
       write(ifigwa)bas,lmxa(ispec(1:nbas))!,qplist,ngplist,ndimhall,qval
    endif
    write(stdo,ftox)' ... Generate radial wave functions ncoremx,nphimx=',ncoremx,nphimx
    allocate(gval(nrmx,0:lmxax,nphimx,nsp,nclass), ecore(ncoremx,nsp,nclass),gcore(nrmx,ncoremx,nsp,nclass),source=0d0)
    do ib = 1, nbas
       is=ispec(ib)
       ic=iclass(ib)
       pnu=>pnuall(:,1:nsp,ib)
       pnz=>pnzall(:,1:nsp,ib)
       a= spec_a(is)
       nr=nris(is)
       z= zz(is)
       nmcore= nmcorex(is)
       CreateAugmentedWaveFunctions: block
         real(8):: gvall(nr*2,0:lmxa(is),nphimx,nsp),gcorel(nr*2,ncoremx*nsp),ecorel(ncoremx*nsp)
         allocate(rofi(nr),rwgt(nr))
         call radmsh(rmt(is),a,nr,rofi)
         call radwgt(rmt(is),a,nr,rwgt)
         rsml= rsmlss(:,is)
         ehl = ehlss(:,is)
!      call atwf(a,lmxa(is),nr,nsp,pnu,pnz,rsml,ehl,rmt(ib),z,v0pot(ib)%v,nphimx,ncore,konfig(0,ib),ecorel, gcorel,gvall,nmcore)
         RadialWaveFunctions: block 
           use m_rhocor,only: getcor
           use m_atwf,only: makrwf,wf2lo
           real(8) :: rhoc(nr,2),gp(2*nr*4),e, phi,dphi,phip,dphip,p,phz,dphz,phzp,dphzp, sumtc,sumec,ez
           do  l = 0, lmxa(is)
              do  isp = 1, nsp
                 call makrwf(z,rofi(nr),l,v0pot(ib)%v(1,isp),a,nr,rofi,pnuall(1,isp,ib),4,gvall(1,l,1,isp),&
                      gp,e,phi,dphi,phip,dphip,p)
                 gvall(:,l,2,isp)=gp(1:2*nr) 
                 if (floor(pnuall(l+1,1,ib)) /= konfig(l,ib)) then
                    call makrwf(z,rofi(nr),l,v0pot(ib)%v(1,isp),a,nr,rofi,pnzall(1,isp,ib),2,gvall(1,l,3,isp),&
                         gp,ez,phz,dphz,phzp,dphzp,p)
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
    enddo
    if(wanatom) close(ifigwa)
! Write refined mesh and indexes to m_rdata1
    call rdata1init(ncores,ndima,ncoremx,konf0,gval,gcore)
! IPW part ! Main loop for eigenfunction generation ==
    open(newunit=ifiqg,file='QGpsi',form='unformatted')
    read(ifiqg ) nqnum, ngpmx ,QpGcut_psi,nqbz,nqirr,imx,nqibz
    close(ifiqg)
    if (lchk>=1 ) then !normcheck file
       open(newunit=ifinormchk,file='norm.'//'procid'//trim(xt(procid))//'.chk')
       write(ifinormchk,"(a)") '#     eval          IPW        IPW(diag)    Onsite(tot)      Total ! lmfgw'
    endif
    if(ham_scaledsigma/=1d0 .AND. sigmamode) write(stdo,*)' Scaled Sigma method: ScaledSigma=',ham_scaledsigma
    ndble = 8
    mrecb = 2*ndima*nbandmx *ndble !byte size !Use -assume byterecl for ifort recognize the recored in the unit of bytes.
    !mrece = nbandmx         *ndble 
    mrecg = 2*ngpmx*nbandmx *ndble 
    allocate(cphix(ndima,nbandmx),geigr(1:ngpmx,1:nbandmx,1:nsp))
    open(newunit=ifcphi,file='CPHI',form='unformatted',access='direct',recl=mrecb)
    open(newunit=ifgeig,file='GEIG',form='unformatted',access='direct',recl=mrecg)
    allocate(qirr(3,nqirr),evl(nbandmx, nqirr, nsp),vxclda(nbandmx, nqirr, nsp))
    iqisploop: do 1001 idat=1,niqisp !iq = iqini,iqend ! iqini:iqend for this procid
       iq  = iqproc(idat) ! iq index
       isp = isproc(idat) ! spin index: isp=1:nspx=nsp/nspc
       qp  = qplist(:,iq) ! q vector
       ngp = ngplist(iq)  ! number of planewaves for PMT basis
       lwvxc = iq<=iqibzmax
       if (cmdopt0('--novxc')) lwvxc = .FALSE. 
       if (socmatrix) lwvxc = .TRUE. 
       call m_Igv2x_setiq(iq) ! Set napw ndimh ndimhx, and igv2x
       allocate(hamm(ndimh,nspc,ndimh,nspc),ovlm(ndimh,nspc,ndimh,nspc)) !Spin-offdiagonal block included since nspc=2 for lso=1.
       allocate(evec(ndimhx,ndimhx),vxc(ndimh,ndimh))
       allocate(cphi(ndima,ndimhx,nsp),cphiw(ndimhx,nsp))
!       allocate(evl(ndimhx,nspx),vxclda(ndimhx))
       if(iqbk==iq) then
          continue
       elseif( lso/=0 .OR. socmatrix) then
          deallocate(hammhso)
          allocate(hammhso(ndimh,ndimh,3))
          call aughsoc(qp, ohsozz,ohsopm, ndimh, hammhso)
          iqbk=iq !q index for hammhso
          if(socmatrix) write(ifiv) iq
          if(socmatrix) write(ifiv) hammhso 
       endif
       if(lwvxc) then
          open(newunit=ifievec, file='evec'//trim(xt(iq))//trim(xt(isp)),form='unformatted')
          open(newunit=ifiv,    file= 'vxc'//trim(xt(iq))//trim(xt(isp)),form='unformatted')
          write(ifievec) ndimh, nlmto
          write(ifiv)    ndimh, nlmto !nlmo: MTO basis
       endif
       if(emptyrun) then !set dummy to avoid error exit
          evec=1d0 
          evl(:,iq,isp)=[(i*0.1,i=1,ndimhx)]
          nev=ndimh
          goto 1212
       endif
       GetHamiltonianAndDiagonalize: block
         integer:: iprint
         if(lso==1) then !L.S case nspc=2
            do ispc=1,nspc  ! nspc==2 ! LDA part of Hamiltonian and overlap matrices for this qp ---
               call hambl(ispc,qp,spotx,vconst,osig,otau,oppix, vxc(:,:),          ovlm(:,ispc,:,ispc))
               call hambl(ispc,qp,smpot,vconst,osig,otau,oppi, hamm(:,ispc,:,ispc),ovlm(:,ispc,:,ispc))
               vxc = hamm(:,ispc,:,ispc) - vxc ! vxc(LDA) part
               hamm(:,ispc,:,ispc)= hamm(:,ispc,:,ispc) + hammhso(:,:,ispc) !spin-diag SOC elements (1,1), (2,2) added
               if(lwvxc) write(ifiv) vxc
            enddo
            hamm(:,1,:,2)= hammhso(:,:,3)                    !spin-offdiagonal SOC elements (1,2) added
            hamm(:,2,:,1)= transpose(dconjg(hammhso(:,:,3)))
            if(sigmamode) then !Add  Vxc(QSGW)-Vxc 
               do ispc=1,nspc
                  call getsenex(qp, ispc, ndimh, ovlm(:,ispc,:,ispc)) !bugfix at 2024-4-24 obata: ispc was 1 when 2023-9-20
                  hamm(:,ispc,:,ispc) = hamm(:,ispc,:,ispc) + ham_scaledsigma*senex !sene= Vxc(QSGW)-Vxc(LDA)
                  call dsene()
               enddo
            endif
         else ! lso=0 (No SO) or lso=2(Lz.Sz)  Spin Diagonal case.spin diagonal) nspc=1 only
            call hambl(isp,qp,spotx,vconst,osig,otau,oppix, vxc,           ovlm(:,1,:,1)) !vxc=<F_i|H(LDA)-vxc(LDA)|F_j>
            call hambl(isp,qp,smpot,vconst,osig,otau,oppi,  hamm(:,1,:,1), ovlm(:,1,:,1)) !ham=<F_i|H(LDA)|F_j> and ovl=<F_i|F_j>
            if(lso==2) hamm(:,1,:,1) = hamm(:,1,:,1) + hammhso(:,:,isp) !diagonal part of SOC matrix added for Lz.Sz mode.
            vxc = hamm(:,1,:,1) - vxc ! vxc(LDA) part
            if(lwvxc) write(ifiv) vxc
            if(sigmamode) then !Add  Vxc(QSGW)-Vxc 
               call getsenex(qp,isp,ndimh,ovlm(:,1,:,1))
               hamm(:,1,:, 1) = hamm(:,1,:,1) + ham_scaledsigma*senex !senex= Vxc(QSGW)-Vxc(LDA)
               call dsene()
            endif
         endif
         if (mod(iq,10) /= 1) call pshpr(iprint()-6)  
         epsovl = ham_oveps
         evec=-1d99
         evl(:,iq,isp)=1d99
         if(magexist) then
            if(nspc==2) then
               hamm(:,1,:,1)= hamm(:,1,:,1) - vmag/2d0*ovlm(:,1,:,1)
               hamm(:,2,:,2)= hamm(:,2,:,2) + vmag/2d0*ovlm(:,2,:,2)
            else
               if(isp==1) hamm(:,1,:,1)= hamm(:,1,:,1) - vmag/2d0*ovlm(:,1,:,1)
               if(isp==2) hamm(:,1,:,1)= hamm(:,1,:,1) + vmag/2d0*ovlm(:,1,:,1)
            endif   
         endif
         call zhev_tk4(ndimhx,hamm,ovlm,ndimhx,nev,evl(1,iq,isp),evec,epsovl) !get evl and evec. 
       endblock GetHamiltonianAndDiagonalize
1212   continue
       if (lwvxc) write(ifievec) qp, evec(1:ndimhx,1:ndimhx),nev
       if(emptyrun) then
          allocate(pwz(ngp,ndimh)) !dummy
          goto 1214
       endif   
       write(stdo,"(' sugw:  kpt isp=',i8,i2,' of ',i8, ' k= ',3f9.5, ' ndimh= ',i5, &
            ' irank=',i4, ' lwvxc=',l,' nev=',i5)")  iq,isp,nqnum,qp,ndimh,procid,lwvxc,nev
       write(stdo,"(9f8.4)") (evl(i,iq,isp), i=1,nev)
       evl(1+nev:,iq,isp)=1d20 
       if(mod(iq,10) /= 1) call poppr
       if(debug) write(stdo,"(' sugw:procid iq isp lwvxc= ',3i3,' ',l)")procid, iq,isp,lwvxc
       nlmax = (lmxax+1)**2 
       gwcphi2: block
         !i   nlmax :leading dimension of aus
         !i   ndham :dimensions aus
         !i   nev   :number of eigenvectors to accumulate cphi
         !i   lmxax :dimensions nlindx
         !i   nlindx: offset that set index in cphi for element (ipqn,l,ib)
         !i   aus   : values of (phi,phidot,pz) at MT sphere boundary; see makusq
         !o Outputs
         !o   cphi :coefficients to phi,phidot,phiz 
         !o        :cphi(ichan,iv) :
         !o          ichan = orbital channel, i.e. one of phi,phidot,phiz (phiz is local orbital for a given site, l, and m; see nlindx.
         !o             iv = eigenvector
         !o   cphiw: diagonal matrix elements, one for each eigenvector. only for check
         !o        : cphiw(1,iv) = <cphi(iv) | overlap | cphi(iv)>
         use m_locpot,only: rotp
         use m_mkpot,only : sab_rv
         integer:: ilm,im,iv
         complex(8):: auasaz(3),aus_zv(nlmax,ndham,3,nsp,nbas)
         call makusq(nbas,[-999], nev,  isp, 1,qp,evec, aus_zv )
         do ispc=1,nspc
            if(lso==1) ispx=ispc
            if(lso/=1) ispx=isp
            !call gwcphi(ispx,nsp,nlmax,ndham,nev,nbas,lmxax,nlindx,ndima,aus_zv, cphi(1,1,ispx),cphiw(1,ispx ))
            cphiw=0d0
            do ib = 1, nbas
               do iv = 1, nev
                  ilm = 0
                  do  l = 0, lmxa(ispec(ib))
                     do im = 1, 2*l+1
                        ilm = ilm+1
                        auasaz=aus_zv(ilm,iv,1:3,ispx,ib) !coefficient for (u,s,gz)
                        cphiw(iv,ispx) = cphiw(iv,ispx) + sum(dconjg(auasaz)*matmul(sab_rv(:,:,l+1,ispx,ib),auasaz)) 
                        ! cphi corresponds to coefficients of augmented functions for {phi,phidot,gz(val=slo=0)}, which are not orthnormal. 
                        cphi(nlindx(1:2,l,ib)+im,iv,ispx)= matmul(auasaz(1:2),rotp(l,ispx,:,:,ib))
                        if (nlindx(3,l,ib) >= 0) cphi(nlindx(3,l,ib) + im,iv,ispx) = auasaz(3)
                     enddo
                  enddo
               enddo
            enddo
         enddo
       endblock gwcphi2
       if(ngp > 0) then !IPW expansion of eigenfunctions pwz 
          !  IPW orthogonalization is |IPWorth_G> = sum_G2 |IPW_G2> O^-1_{G2,G1}, where O_{G1,G2}=<IPW_G1|IPW_G2>
          !  ppovl: ppovl_G1,G2 = O_G1,G2 = <IPW_G1 IPW_G2>
          !  phovl: <IPW_G1 | basis function = smooth Hankel or APW >   
          !    pwz:  IPW expansion of eigen function    
          allocate(ppovl(ngp,ngp),pwz(ngp,ndimhx),phovl(ngp,ndimh))
          call pwmat(nbas,ndimh,napw,igv2x,qp,ngp,nlmax,ngvecp(1,1,iq),gmax, ppovl, phovl ) 
          pwz = matmul(phovl(1:ngp,1:ndimh),evec(1:ndimh,1:ndimhx))
          if(lso==1) pwz = matmul(phovl(1:ngp,1:ndimh),evec(ndimh+1:2*ndimh,1:ndimhx))
          deallocate(phovl)
          if (lchk >= 1) then
             allocate(pzovl,source=pwz)
             allocate(ppovld(ngp)) ! extract diagonal before ppovl overwritten
             forall(i = 1:ngp) ppovld(i) = ppovl(i,i)
          endif
          call matcinv(ngp,ppovl)! inversion of hermitian ppovl
          pwz = matmul(ppovl,pwz) !pwz= O^-1 *phovl * evec  ! IPW expansion of eigenfunction
          deallocate(ppovl)
          if (lchk >= 1) then
             allocate(testc(ndimh,ndimh),testcd(ndimh))
             testc=matmul(transpose(dconjg(pzovl)),pwz)
             deallocate(pzovl)
             testcd = [(sum(dconjg(pwz(:,i))*ppovld*pwz(:,i)),i=1,ndimhx)]
             deallocate(ppovld)
             ! xx(1) = sum over all augmentation w.f.  cphi+ ovl cphi
             ! xx(3) = IPW contribution to phi+ phi.   xx(1)+xx(3) should be close to unity.
             ! [xx(4) = IPW contribution to phi+ phi, using diagonal part only] 
             write(ifinormchk,"('# iq',i5,'   q',3f12.6:'  shortened q',3f12.6)") iq,qp
             do  i1 = 1, ndimhx
                xx(1) = sum(cphiw(i1,1:nspc))
                do  i2 = 1, ndimh
                   xx(3) = testc(i1,i2)
                   xx(4) = testcd(i1)      !if(i1==i2)write(ifinormchk,'(f12.5,5f14.6)')evl(i1,isp),xx(3),xx(4),xx(1),xx(1)+xx(3)
                   if(i1==i2)write(ifinormchk,'(f12.5,5f14.6)')evl(i1,iq,isp),xx(3),xx(4),xx(1),xx(1)+xx(3)
                enddo
             enddo
             deallocate(testc,testcd)
             write(ifinormchk,*)
          endif
       endif
       if(debug) write (stdo,"('q ndimh=',3f10.5,i10)") qp, ndimh
       !     Interstitial part of eigenfunction overlap:
       !     <psi_n| psi_n'> = sum_G1,G2 (pwz_G1,n|IPW_G1>)+  (pwz_G2,n'|IPW_G2>)
       !     = sum_G1,G2 (pwz_G1,n)+ ppovl_G1,G2 (PWZ_G2,n') = (PWZ)+ O (PWZ) = (PZOVL)+ (PWZ)  (old style)
       allocate(testc(ndimh,ndimh),source=matmul(transpose(dconjg(evec)),vxc))
       vxclda(:,iq,isp)=0d0
       forall(i1 = 1:ndimhx) vxclda(i1,iq,isp) = sum(testc(i1,1:ndimh) * evec(1:ndimh,i1))  !<i|Vxc^lda|i>
       deallocate(testc)
1214   continue
       
       writegwb: block
         use m_hamindex0,only: nindx,ibasindx
         integer::iband,ibas,iqqisp,ix,m,nm
!open(newunit=ifigwb_, file='gwb'//trim(xt(iqq))//trim(xt(isp)),form='unformatted')
!   evl(:,iqq,isp)=1d20
!   read(ifigwb_) evl(1:ndimh,iqq,isp),cphir(1:ndima,1:ndimh),geigr(1:ngp,1:ndimh,isp),vxclda(1:ndimh,iqq,isp),nev !nev is # of eigenfunctions.
         geigr(1:ngpmx,1:ndimh,isp)=0d0
         geigr(1:ngp,1:ndimh,isp)=pwz
         do ibas=1,nbas
            do ix = 1,ndima !nindx is for avoiding degeneracy. See zzpi.
               if(ibasindx(ix)==ibas) cphi(ix,1:nev,isp) = cphi(ix, 1:nev,isp)/sqrt(1d0+0.1d0*nindx(ix))
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
               cphix (iord(m,1:nm,l,ib),iband) = cphix (iord(m,1:nm,l,ib),iband) + zzpi(1:nm,n,l,ic,isp)*cphi(ix,iband,isp)
            enddo
         enddo
         iqqisp= isp + nsp*(iq-1)
         write(ifcphi,  rec=iqqisp)  cphix(1:ndima,1:nbandmx)
!         close(ifigwb_)
         iqqisp= isp + nsp*(iq-1)
         if(ngpmx/=0) write(ifgeig,  rec=iqqisp)  geigr(1:ngpmx,1:nbandmx,isp)
       endblock writegwb
       
       open(newunit=ifigwb, file='gwb' //trim(xt(iq))//trim(xt(isp)),form='unformatted')
       if(lso/=1) write(ifigwb) evl(1:ndimhx,iq,isp), vxclda(1:ndimhx,iq,isp) !,cphi(:,:,isp),   pwz,,nev
       if(lso==1) write(ifigwb) evl(1:ndimhx,iq,isp), vxclda(1:ndimhx,iq,isp) !,cphi(:,:,1:nspc),pwz,,nev
       close(ifigwb)
       
       deallocate(pwz)
       if (lwvxc) close(ifiv)
       if (lwvxc) close(ifievec)
       deallocate(hamm,ovlm,evec,vxc,cphi,cphiw)!,evl,vxclda)
1001 enddo iqisploop
    close(ifcphi)
    close(ifgeig)
    !
    call mpi_barrier(comm,ierr)
    if(master_mpi) then
 !       open(newunit=ifv,file='VXCFP',form='unformatted')
 !       write(ifv) nbandmx,nqirr
         irreducibleqloop: do 1200 iqq = 1,nqirr !irreducible points. At qirr, we calculated eigenfunctions.
            qirr(:,iqq)= qplist(:,iqq)
 !           write(ifv) qirr(1:3,iqq), vxclda(1:nbandmx,iqq,1:nsp) ! VXCFP
  1200   enddo irreducibleqloop
! Evalue
!       open(newunit=ifev,file='EValue',form='unformatted')
!       write(ifev) nbandmx, nqirr, nsp
!       write(ifev) qirr(1:3,1:nqirr) !qirr
!       write(ifev) evl(1:nbandmx, 1:nqirr, 1:nsp )
!       close(ifev)
    endif
    if(master_mpi) call rdata4gw()
    call tcx('m_sugw_init')
  end subroutine m_sugw_init
end module m_sugw


  ! subroutine gwcphi(isp,nsp,nlmax,ndham,nev,nbas,lmxax,nlindx,ndima,aus, cphi,cphiw)!- Project (phi,phidot) onto MT sphere
  !   use m_struc_def
  !   use m_lmfinit,only: ispec
  !   use m_locpot,only: rotp
  !   use m_mkpot,only: sab_rv
  !   !i   isp   :current spin channel (1 or 2)
  !   !i   nsp   :2 for spin-polarized case, otherwise 1
  !   !i   nlmax :leading dimension of aus
  !   !i   ndham :dimensions aus
  !   !i   nev   :number of eigenvectors to accumulate cphi
  !   !i   nbas  :number of sites in list
  !   !i   lmxax :dimensions nlindx
  !   !i   nlindx:offset that set index in cphi for element (ipqn,l,ib)
  !   !i   ndima :leading dimension of cphi
  !   !i   aus   :values of (phi,phidot) at MT sphere boundary; see makusq
  !   !o Outputs
  !   !o   cphi :coefficients to phi,phidot,phiz following Kotani conventions
  !   !o        :cphi(ichan,iv) :
  !   !o           ichan = orbital channel, i.e. one of (phi,phidot or phiz(=gz))
  !   !o           for a given site, l, and m; see nlindx.
  !   !o           iv = eigenvector
  !   !o   cphiw:diagonal matrix elements, one for each eigenvector.
  !   !o        :cphiw(1,iv) = <cphi(iv) | overlap | cphi(iv)>
  !   implicit none
  !   integer :: isp,nsp,nlmax,ndham,nbas,nev,lmxax,ndima, nlindx(3,0:lmxax,nbas),ichan,ib,is,igetss,iv,ilm,l,im,i
  !   real(8):: cphiw(nev),wgt
  !   complex(8):: au,as,az,sqrsz(3),auas(2),auasaz(3),aus(nlmax,ndham,3,nsp,*),cphi(ndima,nev)
  !   cphiw=0d0
  !   do ib = 1, nbas
  !      do iv = 1, nev
  !         ilm = 0
  !         do  l = 0, lmxa(ispec(ib))
  !            do im = 1, 2*l+1
  !               ilm = ilm+1
  !               auasaz=aus(ilm,iv,1:3,isp,ib) !coefficient for (u,s,gz)
  !               cphiw(iv) = cphiw(iv) + sum(dconjg(auasaz)*matmul(sab_rv(:,:,l+1,isp,ib),auasaz)) 
  !               !  cphi corresponds to coefficients of augmented functions for {phi,phidot,gz(val=slo=0)}, which are not orthnormal. 
  !               cphi(nlindx(1:2,l,ib)+im,iv)= matmul(auasaz(1:2),rotp(l,isp,:,:,ib))
  !               if (nlindx(3,l,ib) >= 0) cphi(nlindx(3,l,ib) + im,iv) = auasaz(3)
  !            enddo
  !         enddo
  !      enddo
  !   enddo
  ! end subroutine gwcphi
