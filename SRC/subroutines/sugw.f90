module m_sugw
  use m_lgunit,only:stdo
  use m_lmfinit,only: z_i=>z,nr_i=>nr,lmxa_i=>lmxa,rmt_i=>rmt,spec_a
  integer,public::  ndima, ndham, lmxamx, ncoremx,nrmx,nqirr,nqibz
  integer,allocatable,public:: konfig(:,:),lmxa(:),ncores(:)
  ! ldim2,nbandmx,lmxamx, ncoremx,nrmx,plat,alat,nqirr
!  real(8),allocatable::zz(:)
  private
  public:: m_sugw_init
contains
  subroutine m_sugw_init (socmatrix,eferm,vmag,qval)
    use m_ext,only:   sname
    use m_suham,only: ndham=>ham_ndham !max dimension of hamiltonian +napwad (for so=0,2)
    use m_lmfinit,only: ham_pwmode,pwemax,ham_oveps,lrsig=>ham_lsig,nlmto,lso
    use m_lmfinit,only: ham_scaledsigma,lat_alat,nsp,nspc,ispec,nspec
    use m_lmfinit,only: nbas,n0,nppn,nkap0,slabl,nmcorex=>nmcore,iantiferro,lmxax
    use m_lattic,only: lat_plat, lat_qlat,rv_a_opos
    use m_supot,only: n1,n2,n3, lat_gmax
    use m_rdsigm2,only: getsenex, senex,dsene
    use m_mkpot,only: smpot=>osmpot, vconst, vesrmt
    use m_mkpot,only: osig, otau, oppi, ohsozz,ohsopm, oppix,spotx
    use m_MPItk,only: numproc=>nsize,procid,master,master_mpi,comm
    use m_igv2x,only: napw,ndimh,ndimhx,igv2x,m_Igv2x_setiq,ndimhall
    use m_elocp,only: rsmlss=>rsml, ehlss=>ehl
    use m_qplist,only: qplist,ngplist,ngvecp,iqibzmax,niqisp,iqproc,isproc
    use m_hamindex0,only: Readhamindex0, nlindx,nclass
    use m_density,only: v0pot,pnuall,pnzall
    use m_augmbl,only: aughsoc
    use m_makusq,only: makusq
    use m_pwmat,only: pwmat
    use m_atwf,only: atwf,makrwf
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
    ! o   gwa:  site data.
    ! o        *for each site:
    ! o           z, nr, a, b, rmt, lmaxa, nsp, ncore
    ! o           konfig(0:lmaxa) : note
    ! o           rofi: radial mesh
    ! o          *for each l, spin isp
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
    !! memo: check shortned mechanism.
    !! ----------------------------------------------------------------------
    integer:: ham_lsig
    integer:: jobgw=1 , lh(n0)
    real(8):: rsml(n0), ehl(n0) ,eferm,qval
    logical :: lwvxc,cmdopt0!,FPMTmode=.false.
    !logical,optional:: FPMTmodein
    integer :: i,i1,i2,iat,ib,ibr,icore,ierr,ifeigen,&! & ifi,
         ifiqg,iflband(2),ifqeigen,ifsyml,igets,igetss,iix,iline, &
         im1,im2,ipb(nbas),ipqn,ipr,iprint,iq,is,isp,ispc,j,job,k1, &
         k2,k3,konf,l,lchk,ldim,loldpw, &
         lmaxa,lsig,mx,mxint,&!nat, &
         ncore,nevl,nev,nglob,ngp,ngp_p, &
         ngpmx,nline,nlinemax,nlmax,nmx,nn1,nn2,nnn, &
         nphimx,npqn,nqbz,nqnum,nqnumx,nqtot,nr,iqibz,imx, & !,nqibz
         ifigwb,ifigwa,ifinormchk,ifigw1,ifildima,ifigwn,ifigwbhead, &
         ificlass,ifievec,ifievecx,ifigw2,ifiqbz,ifievv,idat
    complex(8),allocatable :: aus_zv(:)
    real(8),allocatable :: ww_rv(:)
    real(8):: QpGcut_psi,QpGcut_cou,dum,xx(5),gmax,ecore(50),a,z,rmt(nbas),b,vshft, &
         alat,alfa,ef0,plat(3,3),qlat(3,3),qp(3),qpos,q_p(3), epsovl! pnu(n0,2),pnz(n0,2)
    real(8),pointer:: pnu(:,:),pnz(:,:)
    integer,allocatable:: ips(:),ipc(:),ipcx(:), ngvecp_p(:,:) 
    integer,allocatable :: konft(:,:,:),iiyf(:),ibidx(:,:),nqq(:)
    real(8),allocatable:: wk(:,:), &
         bas(:,:),rofi(:),rwgt(:),gcore(:,:,:),gval(:,:,:,:,:),evl(:,:),vxclda(:)
    real(8),allocatable::  cphiw(:,:) !!ovv(:,:),evl_p(:,:), qq1(:,:),qq2(:,:),
    complex(8),allocatable:: evec(:,:),evec0(:,:),vxc(:,:),&! & dipo(:,:,:),
         ppovl(:,:),phovl(:,:),pwh(:,:),pwz(:,:),pzovl(:,:), pwz0(:,:),&
         testc(:,:),testcd(:),ppovld(:),cphi(:,:,:),cphi0(:,:,:),cphi_p(:,:,:), &
         geig(:,:,:),geig_p(:,:,:),sene(:,:)
    integer::isize_ham(3)
    integer :: ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0),norb
    integer :: nspx,nnnx,ifiv
    character strn*120
    integer :: pwmode
    double precision :: pwgmin,pwgmax,eomin
    integer:: w(1) !! dummy
    real(8):: rdummy(1)
    real(8):: dnn(3),qlatinv(3,3),qout(3),qtarget(3),qrr(3),axx,bxx,qxxx(3),qxx1(3),qxx2(3)
    integer:: iqzz,nqzz,iiiii
    logical:: debug=.false.
    complex(8),allocatable:: evecout(:,:),evecr(:,:)
    real(8),allocatable:: qzz(:,:)
    logical:: l_dummy_isanrg,isanrg,oncewrite,newsigmasw,sigmamode !,noshorbz
    integer::   irr, nmcore !feb2012takao nk1,nk2,nk3,,nqirr
    integer:: nbalance,nrr,ifiproc,iproc,iadd,ldw
    integer:: iqq
    logical:: rank0
    character*256:: ext,sprocid,extn
    logical ::nexist,magexist  !! june2013 magfield is added for fsmom mode.
    real(8):: vnow
    integer:: ifimag
    integer,parameter:: noutmx=48
    integer:: iout,nout,nlatout(3,noutmx),iapw,ival
    real(8)::ppin(3)
    real(8):: rlatp(3,3),xmx2(3)
    real(8):: qqq(3),vmag
    integer:: mode,iwdummy,jx
    complex(8),allocatable:: hamm(:,:,:,:),ovlm(:,:,:,:),ovlmtoi(:,:),ovliovl(:,:) ,hammhso(:,:,:)
    integer:: lpdiag=0!,nrmx!,ncoremx!,ndham
    character spid*8
    character(8) :: xt
    logical,optional:: socmatrix !dipolematrix,
    integer:: ispSS,ispEE,ispx,iqbk
    logical:: emptyrun
    include "mpif.h"
    call tcn ('m_sugw_init')
    emptyrun=cmdopt0('--emptyrun')
    if(master_mpi) write(stdo,"(a)") 'm_sugw_init: start'
    sigmamode = mod(lrsig,10) .ne. 0
    call getpr(ipr)
    write(stdo,"('MagField added to Hailtonian -vmag/2 for isp=1, +vmag/2 for isp=2: vmag(Ry)=',d13.6)") vmag
    magexist= abs(vmag)>1d-6
    alat=lat_alat
    plat =lat_plat
    qlat =lat_qlat
    gmax =lat_gmax
    lchk = 1
    lwvxc = .not. cmdopt0('--novxc')
    ldim=nlmto
    pwmode= ham_pwmode
    ndima = 0
    do  ib = 1, nbas
       is=ispec(ib)
       lmaxa=lmxa_i(is)
       pnz=>pnzall(:,1:nsp,ib)
       if(sum(abs(pnz(:,1:nsp)-pnzall(:,1:nsp,ib)))>1d-9) call rx('sugw xxx1aaa')
!       lmxax = max(lmxax,lmaxa)
!       if (lmaxa > -1) then
       do  l = 0, lmaxa
          npqn = 2
          if (pnz(l+1,1) /= 0) npqn = 3
          ndima = ndima + npqn*(2*l+1)
       enddo
!       endif
    enddo
    call Readhamindex0() ! ==== Read file NLAindx ====
    open(newunit=ifiqg,file='QGpsi',form='unformatted')
    read(ifiqg ) nqnum, ngpmx ,QpGcut_psi,nqbz,nqirr,imx,nqibz
    !! ==== Write, or read past header information, file gwb ====
    sprocid='procid'//trim(xt(procid))
    if(master_mpi) open(newunit=ifigwa,file='gwa',form='unformatted')
    ef0 = 1d99        
    allocate(ips(nbas),lmxa(nbas),bas(3,nbas))
    bas(:,1:nbas)=rv_a_opos(:,1:nbas) 
    lmxa(1:nbas) = lmxa_i(ispec(1:nbas))
    print *,'iantiferro=',iantiferro
    nphimx = 0
    allocate(konfig(0:lmxax,nbas))
    do  ib = 1, nbas
       is=ispec(ib) 
       pnu=>pnuall(:,1:nsp,ib)
       pnz=>pnzall(:,1:nsp,ib)
       a=spec_a(is)
       nr=nr_i(is)
       z= z_i(is)
       rmt(ib)=rmt_i(is)
       lmaxa = lmxa_i(is)
       nmcore= nmcorex(is)
       if (lmaxa > -1) then
          call atwf(0,a,lmaxa,nr,nsp,pnu,pnz,rsml,ehl,rmt(ib),z,rdummy,i1,ncore,konfig(0,ib),ecore,rdummy,rdummy,nmcore)
          nphimx = max(nphimx,i1)
       endif
    enddo
    !! Atom data (gwa)
    write(stdo,*)' ... Generate core wave functions (file gwa)'
    ncoremx=0
    nrmx=0
    !    allocate(zz(nclass))
    allocate(ncores(nspec))
    do  ib = 1, nbas
       is=ispec(ib) 
       pnu=>pnuall(:,1:nsp,ib)
       pnz=>pnzall(:,1:nsp,ib)
       a= spec_a(is)
       nr=nr_i(is)
       z= z_i(is)
       rmt(ib)= rmt_i(is)
       lmaxa =  lmxa_i(is)
       spid=slabl(is)
       if (lmaxa > -1) then
          call atwf ( 0 , a , lmaxa , nr , nsp , pnu , pnz , rsml , ehl &
               , rmt ( ib ) , z , v0pot(ib)%v , i1 , ncore , konfig(0,ib) , ecore , rdummy , rdummy ,nmcore)
          allocate(rofi(nr),rwgt(nr),gcore(nr,2,ncore))
          allocate(gval(nr,2,0:lmaxa,nphimx,nsp))
          gval=0d0
          !     !         Create augmented wave functions for this atom
          rsml=rsmlss(:,is)
          ehl= ehlss(:,is)
          call atwf ( 03 , a , lmaxa , nr , nsp , pnu , pnz , rsml , ehl &
               , rmt ( ib ) , z , v0pot(ib)%v , nphimx , ncore , konfig(0,ib) , ecore , gcore , gval ,nmcore)
          if(nr     >nrmx   ) nrmx    = nr
          if(ncore  >ncoremx) ncoremx = ncore
          !     !         Header data for this atom
          b = rmt(ib)/(dexp(a*nr-a)-1d0)
          call radmsh(rmt(ib),a,nr,rofi)
          call radwgt(rmt(ib),a,nr,rwgt)
          !zz(is)=z
          !          ncore_d=
          ncores(is)=ncore/nsp
          if(master_mpi)write(ifigwa) z, nr, a, b, rmt(ib), lmaxa, nsp, ncore,spid
          if(master_mpi)write(ifigwa) konfig(0:lmaxa,ib)
          if(master_mpi)write(ifigwa) rofi
          !     !         Write orthonormalized valence wave functions for this atom
          do  l = 0, lmaxa
             do  i = 1, nsp
                if(master_mpi)write(ifigwa) l,i
                if(master_mpi)write(ifigwa) gval(1:nr,1,l,1,i)
                if(master_mpi)write(ifigwa) gval(1:nr,1,l,2,i)
                if (konfig(l,ib) >= 10 .AND. master_mpi) write(ifigwa) gval(1:nr,1,l,3,i)
             enddo
          enddo
          !     !         Core wave functions for this atom
          icore = 0
          vshft = vesrmt(ib)
          !     !        As of v6.11, shift is included in v0, passed in vval to
          !     !         locpot, in routine mkpot.f
          vshft = 0
          do  l = 0, lmaxa
             do  isp = 1, nsp
                do  konf = l+1, mod(konfig(l,ib),10)-1
                   icore = icore+1
                   if(master_mpi)write(ifigwa) icore, l, isp, konf, ecore(icore)+vshft
                   if(master_mpi)write(ifigwa) gcore(1:nr,1,icore)
                enddo
             enddo
          enddo
          !   radial integral test block radial functions when energy/pnu are changing.
          radint: block
            integer::ir,ie
            real(8):: enu,sum1,sum2,gfun(1:nr,2),gpfun(1:nr,2,4), rgfun(1:nr)
            real(8):: phi,dphi,phip,dphip,p,pnux(4)
            if(cmdopt0('--radialintegraltest')) then
               l = 3            !0, lmaxa
               do i = 1, nsp
                  do ie=-12,10
                     pnux= [0d0,0d0,0d0,4.65d0 + ie*0.05d0] !pnu is changing
                     if(ie==7) cycle
                     call makrwf(z,rmt(ib),l,v0pot(ib)%v(:,i),a,nr,rofi,pnux,2, gfun,gpfun,enu, phi,dphi,phip,dphip,p)
                     call gintxx(gfun(1:nr,1),gfun(1:nr,1),a,b,nr,sum1)
                     rgfun=[(rofi(ir)*gfun(ir,1),ir=1,nr)]
                     if(abs(sum1-1d0)>0.01d0) call rx('sugw: need to check normalization 111')
                     call gintxx(rgfun,gfun(1:nr,1),a,b,nr,sum2)
                     write(stdo,ftox)'\int phiphi : ib l isp=',ib,l,i,' pnu \int rphiphi=',ftof(pnux,3),ftof(sum2,3)
                  enddo
                  write(stdo,*)
               enddo
            endif
          endblock radint
          deallocate(rofi,rwgt,gcore,gval)
       endif
    enddo
    if(cmdopt0('--radialintegraltest')) call rx0('end of radial integral')
    if(master_mpi) then
       write(ifigwa) iantiferro(1:nbas) !iantiferro may2015
       close(ifigwa)
    endif
    lmxamx=maxval(lmxa(1:nbas))
    if(master_mpi) then
       open(newunit=ifigwbhead,file='gwb.head',form='unformatted')
       ncoremx=ncoremx/nsp
       write(ifigwbhead)nbas,nsp,ndima,ndham,maxval(lmxa(1:nbas)),ncoremx,nrmx,plat,alat,nqirr,nqibz
       write(ifigwbhead)bas,lmxa,qplist,ngplist,ndimhall,qval
       close(ifigwbhead)
    endif
    deallocate(ips,lmxa)
    nspx=nsp/nspc
    ! For SO=1,  nsp=2, nspc=2, nspx=1,   ndimhx=ndimh*2
    ! For SO/=1,        nspc=1, nspx=nsp, ndimhx=ndimh (nsp=1 or 2)
    !! == GW setup loop over k-points ==
    if (lchk>=1 ) then
       open(newunit=ifinormchk,file='norm.'//trim(sprocid)//'.chk')
       write(ifinormchk,"(a)") '#     eval          IPW        IPW(diag)    Onsite(tot)      Total ! lmfgw'
    endif
    !! --- Evecs and matrix elements of vxc for irr qp ---
    !!    Note: this routine should use only irr qp.
    !! == Main loop for eigenfunction generation ==
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
          write(ifievec) ndimh, ldim
          write(ifiv)    ndimh, ldim
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
       call makusq(nbas,[-999], nev,  isp, 1 , qp , evec ,  aus_zv )
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
             allocate(pzovl(ngp,ndimhx))
             pzovl = pwz
             allocate(ppovld(ngp)) ! extract diagonal before ppovl overwritten
             do  i = 1, ngp
                ppovld(i) = ppovl(i,i)
             enddo
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
!rdata4gw_v2
    rdata: block
!      use m_rdata4gw,only: rdata4gw
      call mpi_barrier(comm,ierr)
      if(master_mpi) then
         call rdata4gw()
      endif   
    endblock rdata
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
    integer :: isp,nsp,nlmax,ndham,nbas,nev,lmxax,ndima, nlindx(3,0:lmxax,nbas),& !,ipb(nbas)
         lmxa,ichan,ib,is,igetss,iv,ilm,l,im,i,ibas
    integer,parameter ::n0=10
    real(8):: cphiw(nev),wgt
    complex(8):: au,as,az,sqrsz(3),auas(2),auasaz(3),aus(nlmax,ndham,3,nsp,*),cphi(ndima,nev)
    cphiw=0d0
    do ib = 1, nbas
       do iv = 1, nev
          ilm = 0
          do  l = 0, lmxa_i(ispec(ib))
             do im = 1, 2*l+1
                ilm = ilm+1
                auasaz=aus(ilm,iv,1:3,isp,ib) !coefficient for (u,s,gz)
                cphiw(iv) = cphiw(iv) + sum(dconjg(auasaz)*matmul(sab_rv(:,:,l+1,isp,ib),auasaz)) 
                !  cphi corresponds to coefficients of augmented functions for
                !   {phi,phidot,gz(val=slo=0)}, which are not orthnormal. 
                cphi(nlindx(1:2,l,ib)+im,iv)= matmul(auasaz(1:2),rotp(l,isp,:,:,ib))
                if (nlindx(3,l,ib) >= 0) cphi(nlindx(3,l,ib) + im,iv) = auasaz(3)
             enddo
          enddo
       enddo
    enddo
  end subroutine gwcphi
end module m_sugw
