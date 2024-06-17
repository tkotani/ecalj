module m_sugw
  use m_lgunit,only:stdo
  use m_lmfinit,only: z_i=>z,nr_i=>nr,lmxa,rmt_i=>rmt,spec_a
  real(8),allocatable,public::ecore(:,:,:),gcore(:,:,:,:),gval(:,:,:,:,:)
  integer,public::  ndima, ndham, lmxamx, ncoremx,nqirr,nqibz
  integer,allocatable,public:: konfig(:,:),ncores(:)
  ! ldim2,nbandmx,lmxamx, ncoremx,nrmx,plat,alat,nqirr
!  real(8),allocatable::zz(:)
  private
  public:: m_sugw_init
contains
  subroutine m_sugw_init (socmatrix,eferm,vmag,qval)
    use m_ext,only:   sname
    use m_suham,only: ndham=>ham_ndham !max dimension of hamiltonian +napwad (for so=0,2)
    use m_lmfinit,only: ham_pwmode,pwemax,ham_oveps,lrsig=>ham_lsig,nlmto,lso,nspx !nspx=nsp/nspc
    use m_lmfinit,only: ham_scaledsigma,alat=>lat_alat,nsp,nspc,ispec,nspec
    use m_lmfinit,only: nbas,n0,nppn,nkap0,slabl,nmcorex=>nmcore,iantiferro,lmxax,nrmx
    use m_lattic,only: plat=>lat_plat, qlat=>lat_qlat,bas=>rv_a_opos
    use m_supot,only: n1,n2,n3, gmax=>lat_gmax
    use m_rdsigm2,only: getsenex, senex,dsene
    use m_mkpot,only: smpot=>osmpot, vconst, vesrmt
    use m_mkpot,only: osig, otau, oppi, ohsozz,ohsopm, oppix,spotx
    use m_MPItk,only: numproc=>nsize,procid,master,master_mpi,comm
    use m_igv2x,only: napw,ndimh,ndimhx,igv2x,m_Igv2x_setiq,ndimhall
    use m_elocp,only: rsmlss=>rsml, ehlss=>ehl
    use m_qplist,only: qplist,ngplist,ngvecp,iqibzmax,niqisp,iqproc,isproc
    use m_hamindex0,only: Readhamindex0, nlindx,nclass,iclass=>iclasst
    use m_density,only: v0pot,pnuall,pnzall
    use m_augmbl,only: aughsoc
    use m_makusq,only: makusq
    use m_pwmat,only: pwmat
    use m_ftox
    use m_zhev,only: zhev_tk4
    use m_hambl,only: hambl
    implicit none
    intent(in)::          socmatrix,eferm,vmag,qval
    !! == Driver for fpgw (to prepare eigenfuncitons for fpgw) ==
    !! NOTE: following documents are not carefully examined. Not believe everything.
    ! i Inputs
    ! ! qval is passed to gwb
    ! i   osig,otau,oppi  augmentation matrices, s_rv1
    ! i   sham%rv_a_ohrs: real space Sigma_vxc
    ! i   nbas  :size of basis
    ! i   smpot :smooth potential on uniform mesh (mkpot.f)
    ! i   vconst:constant to be added to potential
    ! i   lcplxp:0 if ppi is real; 1 if ppi is complex
    ! i   ppn   :potential parameters, nmto style
    ! i   vesrmt  :electrostatic potential at MT boundaries
    ! i   jobgw :-999 prompt for and read jobgw from stdin
    ! i         :0 create files SYMOPS,LATTC,CLASS,NLAindx
    ! i         :1 create files gwb,gw1,gw2,gwa,vxc,evec,rhoMT.*,normchk
    ! o Outputs
    ! o   Files written according to jobgw
    !!
    ! o   gwb.head:  Information about eigenfunctions
    ! o   gwb.*:  Information about eigenfunctions
    !!        write(ifigwb) evl(1:ndimh,isp), cphi(:,:,isp), pwz, vxclda(1:ndimh), nev
    !!        pwz  (G vectors for psi,vcou; PW expansion of z)
    ! o   site data.
    ! o        *for each site:
    ! o           z, nr, a, b, rmt, lmxa, nsp, ncore
    ! o           konfig(0:lmxax) : note
    ! o           rofi: radial mesh
    ! o        *for each l, spin isp
    ! o             l, isp
    ! o             radial w.f. gval: phi
    ! o             radial w.f. gval: phidot
    ! o             radial w.f. gval: phiz    written if konfig(l)>10
    ! o             *for each l, spin, konf
    ! o                icore, l, isp, konf, ecore(icore)+vshft
    ! o                gcore(1:nr,1,icore)
    ! o   evec: eigenvectors.
    ! o         ndham, nsp, nnn, nqnum
    ! o        *for each q-point and spin:
    ! o           q, evec(1:ndimh,1:ndimh)
    ! o   vxc:  matrix elements of XC potential
    ! o         ndham, nsp, nnn
    ! o        *for each q-point and spin:
    ! o           q, vxc
    integer :: lchk=1,i,i1,i2,iat,ib,ibr,icore,ierr,ifeigen,&
         ifiqg,iflband(2),ifqeigen,ifsyml,igets,igetss,iix,iline, &
         im1,im2,ipb(nbas),ipqn,ipr,iprint,iq,is,isp,ispc,j,job,k1, &
         k2,k3,konf,l,ldim,loldpw, lsig,mx,mxint, ncore,nevl,nev,nglob,ngp,ngp_p, &
         ngpmx,nline,nlinemax,nlmax,nmx,nn1,nn2,nnn, &
         nphimx,npqn,nqbz,nqnum,nqnumx,nqtot,nr,iqibz,imx, & !,nqibz
         ifigwb,ifinormchk,ifigw1,ifildima,ifigwn,ifigwbhead, &
         ificlass,ifievec,ifievecx,ifigw2,ifiqbz,ifievv,idat
    integer:: ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0),norb,nnnx,ifiv
    integer:: ispSS,ispEE,ispx,iqbk,konfigk,konfz,icore2,icore2o,ic,nmcore
    real(8):: rsml(n0), ehl(n0) ,eferm,qval, vmag,vnow
    real(8):: QpGcut_psi,QpGcut_cou,dum,xx(5),a,z,rmt(nbas),vshft, qp(3),qpos,q_p(3), epsovl
    real(8),allocatable:: rofi(:),rwgt(:),evl(:,:),vxclda(:), cphiw(:,:) 
    real(8),pointer:: pnu(:,:),pnz(:,:)
    integer,allocatable :: konft(:,:,:),iiyf(:),ibidx(:,:),nqq(:)
    complex(8),allocatable :: aus_zv(:), hamm(:,:,:,:),ovlm(:,:,:,:),ovlmtoi(:,:),ovliovl(:,:) ,hammhso(:,:,:)
    complex(8),allocatable:: evec(:,:),evec0(:,:),vxc(:,:),&
         ppovl(:,:),phovl(:,:),pwh(:,:),pwz(:,:),pzovl(:,:), pwz0(:,:),&
         testc(:,:),testcd(:),ppovld(:),cphi(:,:,:),cphi0(:,:,:),cphi_p(:,:,:), &
         geig(:,:,:),geig_p(:,:,:),sene(:,:)
    logical :: lwvxc,cmdopt0, emptyrun, magexist, debug=.false., sigmamode 
    logical,optional:: socmatrix 
    character(8) :: xt
    character*256:: ext,sprocid,extn
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
       rmt(ib)=rmt_i(is)
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
    write(stdo,ftox)' ... Generate core wave functions ncoremx nphimx per spin=',ncoremx,nphimx
    allocate(gval(nrmx,0:lmxax,nphimx,nsp,nclass), ecore(ncoremx,nsp,nclass),gcore(nrmx,ncoremx,nsp,nclass),source=0d0) 
    do ib = 1, nbas
       is=ispec(ib)
       ic=iclass(ib)
       pnu=>pnuall(:,1:nsp,ib)
       pnz=>pnzall(:,1:nsp,ib)
       a= spec_a(is)
       nr=nr_i(is)
       z= z_i(is)
       nmcore= nmcorex(is)
       CreateAugmentedWaveFunctions: block
         real(8):: gvall(nr*2,0:lmxa(is),nphimx,nsp),gcorel(nr,2,ncoremx*nsp),ecorel(ncoremx*nsp)
         allocate(rofi(nr),rwgt(nr))
         call radmsh(rmt(ib),a,nr,rofi)
         call radwgt(rmt(ib),a,nr,rwgt)
         rsml= rsmlss(:,is)
         ehl = ehlss(:,is)
!         call atwf(a,lmxa(is),nr,nsp,pnu,pnz,rsml,ehl,rmt(ib),z,v0pot(ib)%v,nphimx,ncore,konfig(0,ib),ecorel, gcorel,gvall,nmcore)
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
                  gcore(1:nr,icore2,isp,ic) = gcorel(1:nr,1,icore)
               enddo
            enddo
         enddo
       endblock CreateAugmentedWaveFunctions
       deallocate(rofi,rwgt)
    enddo
    lmxamx=maxval(lmxa)
    ndima = 0
    do  ib = 1, nbas
       pnz=>pnzall(:,1:nsp,ib)
       if(sum(abs(pnz(:,1:nsp)-pnzall(:,1:nsp,ib)))>1d-9) call rx('sugw xxx1aaa')
       ndima = ndima + sum([((2*l+1)*merge(3,2,pnz(l+1,1)/=0), l=0,lmxa(ispec(ib)))])
    enddo
! gwb.head is virtually unused (only used in wannier plot mode. but easily removed...)
    if(master_mpi) then !only for wannier plot mode. read by lmf2gw. Irrelevant
       open(newunit=ifigwbhead,file='gwb.head',form='unformatted') 
       write(ifigwbhead)nbas,nsp,ndima,ndham,maxval(lmxa),ncoremx,nrmx,plat,alat,nqirr,nqibz
       write(ifigwbhead)bas,lmxa(ispec(1:nbas)),qplist,ngplist,ndimhall,qval
       close(ifigwbhead)
    endif
! IPW part ---------------------------
    ! For SO=1,  nsp=2, nspc=2, nspx=1,   ndimhx=ndimh*2
    ! For SO/=1,        nspc=1, nspx=nsp, ndimhx=ndimh (nsp=1 or 2)
    !! == GW setup loop over k-points ==
    !! --- Evecs and matrix elements of vxc for irr qp ---
    !!    Note: this routine should use only irr qp.
    !! == Main loop for eigenfunction generation ==
    open(newunit=ifiqg,file='QGpsi',form='unformatted')
    read(ifiqg ) nqnum, ngpmx ,QpGcut_psi,nqbz,nqirr,imx,nqibz
    if (lchk>=1 ) then
       open(newunit=ifinormchk,file='norm.'//'procid'//trim(xt(procid))//'.chk')
       write(ifinormchk,"(a)") '#     eval          IPW        IPW(diag)    Onsite(tot)      Total ! lmfgw'
    endif
    if(ham_scaledsigma/=1d0 .AND. sigmamode) write(stdo,*)' Scaled Sigma method: ScaledSigma=',ham_scaledsigma
    iqisploop: do 1001 idat=1,niqisp !iq = iqini,iqend ! iqini:iqend for this procid
       iq  = iqproc(idat)
       isp = isproc(idat) ! NOTE: isp=1:nspx=nsp/nspc
       qp  = qplist(:,iq) ! q vector
       ngp = ngplist(iq)  ! number of planewaves for PMT basis
       lwvxc = iq<=iqibzmax
       if (cmdopt0('--novxc')) lwvxc = .FALSE. 
       if (socmatrix) lwvxc = .TRUE. 
       call m_Igv2x_setiq(iq) ! Set napw ndimh ndimhx, and igv2x
       allocate(hamm(ndimh,nspc,ndimh,nspc),ovlm(ndimh,nspc,ndimh,nspc)) !Spin-offdiagonal block included since nspc=2 for lso=1.
       allocate(evec(ndimhx,ndimhx),vxc(ndimh,ndimh))
       allocate(cphi(ndima,ndimhx,nsp),cphiw(ndimhx,nsp))
       allocate(evl(ndimhx,nspx),vxclda(ndimhx))
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
       open(newunit=ifigwb, file='gwb' //trim(xt(iq))//trim(xt(isp)),form='unformatted')
       if(lwvxc) then
          open(newunit=ifievec, file='evec'//trim(xt(iq))//trim(xt(isp)),form='unformatted')
          open(newunit=ifiv,    file= 'vxc'//trim(xt(iq))//trim(xt(isp)),form='unformatted')
          write(ifievec) ndimh, nlmto
          write(ifiv)    ndimh, nlmto !nlmo: MTO basis
       endif
       if(emptyrun) then !set dummy to avoid error exit
          evec=1d0 
          evl(:,isp)=[(i*0.1,i=1,ndimhx)]
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
         evl(:,isp)=1d99
         if(magexist) then
            if(nspc==2) then
               hamm(:,1,:,1)= hamm(:,1,:,1) - vmag/2d0*ovlm(:,1,:,1)
               hamm(:,2,:,2)= hamm(:,2,:,2) + vmag/2d0*ovlm(:,2,:,2)
            else
               if(isp==1) hamm(:,1,:,1)= hamm(:,1,:,1) - vmag/2d0*ovlm(:,1,:,1)
               if(isp==2) hamm(:,1,:,1)= hamm(:,1,:,1) + vmag/2d0*ovlm(:,1,:,1)
            endif   
         endif
         call zhev_tk4(ndimhx,hamm,ovlm,ndimhx,nev,evl(1,isp),evec,epsovl) !get evl and evec. 
       endblock GetHamiltonianAndDiagonalize
1212   continue
       if (lwvxc) write(ifievec) qp, evec(1:ndimhx,1:ndimhx),nev
       if(emptyrun) then
          allocate(pwz(ngp,ndimh)) !dummy
          goto 1214
       endif   
       write(stdo,"(' sugw:  kpt isp=',i8,i2,' of ',i8, ' k= ',3f9.5, ' ndimh= ',i5, &
            ' irank=',i4, ' lwvxc=',l,' nev=',i5)")  iq,isp,nqnum,qp,ndimh,procid,lwvxc,nev
       write(stdo,"(9f8.4)") (evl(i,isp), i=1,nev)
       evl(1+nev:,isp)=1d20 ! See rdata4gw_v2
       if(mod(iq,10) /= 1) call poppr
       if(debug) write(stdo,"(' sugw:procid iq isp lwvxc= ',3i3,' ',l)")procid, iq,isp,lwvxc
!       if(magexist) then
!          if(nspc==2) call rx('not yet implemented magexist=T and SO=1')
!          if(isp==1) evl(1:ndimh,isp)=evl(1:ndimh,isp) - vmag/2d0
!          if(isp==2) evl(1:ndimh,isp)=evl(1:ndimh,isp) + vmag/2d0
!       endif
       nlmax = (lmxax+1)**2 
       allocate(aus_zv(nlmax*ndham*3*nsp*nbas))     ! Project wf into augmentation spheres, Kotani conventions ---
       aus_zv=0d0
       call makusq(nbas,[-999], nev,  isp, 1,qp,evec, aus_zv )
       do ispc=1,nspc
          if(lso==1) ispx=ispc
          if(lso/=1) ispx=isp
          call gwcphi(ispx,nsp,nlmax,ndham,nev,nbas,lmxax,nlindx,ndima,aus_zv, cphi(1,1,ispx),cphiw(1,ispx ))
       enddo
       deallocate(aus_zv)
       !     ! We keep note in the followings, but be careful (may contain bugs)...
       !     !  --- Overlap of IPWs, PW expansion of eigenfunctions pwz ---
       !     !      The IPW basis consisting of PWs with holes knocked out of spheres,
       !     !      IPWs must be orthogonalized: IPW -> IPWbar; see PRB 76, 165106.
       !     !         |IPWbar_G> = sum_G2 |IPW_G2> O^-1_G2,G1, where O_G1,G2=<IPW_G1|IPW_G2>
       !     ! == Definitions ==
       !     * ppovl = overlap of IPWs (Generated by pwmat and pwmat2)
       !     ppovl_G1,G2 = O_G1,G2 = <IPW_G1 IPW_G2>
       !     * pwh = PW expansion of basis function (Generated by pwmat2 only)
       !     basis_j> = sum_G2 PWH_G2,j |IPW_G2>
       !     * Matrix elements (overlap) of IPW and basis function (Generated by pwmat only)
       !     phovl_G1,j = sum_G2 ppovl_G1,G2 pwh'_G2,j  (matrix form PHOVL = O * PWH')
       !     Note: phovl is only used as an intermediate construction, old branch
       !     * Note: pwh' is expanded to the LMTO cutoff gmax while
       !     pwh  is expanded to the GW cutoff QpGcut_psi
       !     Thus O^-1 PHOVL will not identically recover PWH.
       !     The original branch (loldpw=0) uses effectively PWH'; the new one uses PWH.
       !     This is a major distinction between the two (see Remarks)
       !     * PW expansion of eigenfunction:
       !     |psi_n> = sum_j z_jn |basis_j>
       !     * Define pwz_G,n = sum_j PWH_G2,j z_jn  (in matrix form: PWZ = PWH Z)
       !     Then
       !     |psi_n> = sum_j,G2 z_jn PWH_G2,j |IPW_G2> = sum_G2 PWZ_G2,n |IPW_G2>
       !     Overlap of IPW and eigenfunction:
       !     PZOVL_G1,n = <IPW(G1) psi_n> = sum_G2 O_G1_G2 PWZ_G2,n
       !     PZOVL = O * PWZ (matrix form) <--- old
       !     Note: pzovl is only used as an intermediate construction, old branch
       if(ngp > 0) then
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
          pwz = matmul(ppovl,pwz)
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
                   if(i1==i2)write(ifinormchk,'(f12.5,5f14.6)')evl(i1,isp),xx(3),xx(4),xx(1),xx(1)+xx(3)
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
       forall(i1 = 1:ndimhx) vxclda(i1) = sum(testc(i1,1:ndimh) * evec(1:ndimh,i1))  !<i|Vxc^lda|i>
       deallocate(testc)
1214   continue
       if(lso/=1) write(ifigwb) evl(1:ndimhx,isp),cphi(:,:,isp),   pwz,vxclda(1:ndimhx),nev
       if(lso==1) write(ifigwb) evl(1:ndimhx,isp),cphi(:,:,1:nspc),pwz,vxclda(1:ndimhx),nev
       deallocate(pwz)
       close(ifigwb)
       if (lwvxc) close(ifiv)
       if (lwvxc) close(ifievec)
       deallocate(hamm,ovlm,evec,vxc,cphi,cphiw,evl,vxclda)
1001 enddo iqisploop
    close(ifiqg)
    call mpi_barrier(comm,ierr)
    if(master_mpi) call rdata4gw()
    call tcx('m_sugw_init')
  end subroutine m_sugw_init
  
  subroutine gwcphi(isp,nsp,nlmax,ndham,nev,nbas,lmxax,nlindx,ndima,aus, cphi,cphiw)!- Project (phi,phidot) onto MT sphere
    use m_struc_def
    use m_lmfinit,only: ispec
    use m_locpot,only: rotp
    use m_mkpot,only: sab_rv
    !i   isp   :current spin channel (1 or 2)
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   nlmax :leading dimension of aus
    !i   ndham :dimensions aus
    !i   nev   :number of eigenvectors to accumulate cphi
    !i   nbas  :number of sites in list
    !i   lmxax :dimensions nlindx
    !i   nlindx:offset that set index in cphi for element (ipqn,l,ib)
    !i   ndima :leading dimension of cphi
    !i   aus   :values of (phi,phidot) at MT sphere boundary; see makusq
    !o Outputs
    !o   cphi :coefficients to phi,phidot,phiz following Kotani conventions
    !o        :cphi(ichan,iv) :
    !o           ichan = orbital channel, i.e. one of (phi,phidot or phiz(=gz))
    !o           for a given site, l, and m; see nlindx.
    !o           iv = eigenvector
    !o   cphiw:diagonal matrix elements, one for each eigenvector.
    !o        :cphiw(1,iv) = <cphi(iv) | overlap | cphi(iv)>
    implicit none
    integer :: isp,nsp,nlmax,ndham,nbas,nev,lmxax,ndima, nlindx(3,0:lmxax,nbas),ichan,ib,is,igetss,iv,ilm,l,im,i
    real(8):: cphiw(nev),wgt
    complex(8):: au,as,az,sqrsz(3),auas(2),auasaz(3),aus(nlmax,ndham,3,nsp,*),cphi(ndima,nev)
    cphiw=0d0
    do ib = 1, nbas
       do iv = 1, nev
          ilm = 0
          do  l = 0, lmxa(ispec(ib))
             do im = 1, 2*l+1
                ilm = ilm+1
                auasaz=aus(ilm,iv,1:3,isp,ib) !coefficient for (u,s,gz)
                cphiw(iv) = cphiw(iv) + sum(dconjg(auasaz)*matmul(sab_rv(:,:,l+1,isp,ib),auasaz)) 
                !  cphi corresponds to coefficients of augmented functions for {phi,phidot,gz(val=slo=0)}, which are not orthnormal. 
                cphi(nlindx(1:2,l,ib)+im,iv)= matmul(auasaz(1:2),rotp(l,isp,:,:,ib))
                if (nlindx(3,l,ib) >= 0) cphi(nlindx(3,l,ib) + im,iv) = auasaz(3)
             enddo
          enddo
       enddo
    enddo
  end subroutine gwcphi
end module m_sugw
