module m_sugw
  use m_lgunit,only:stdo
  private
  public:: m_sugw_init
contains
  subroutine m_sugw_init (socmatrix,eferm) !dipolematrix,
    use m_ext,only:   sname
    use m_struc_def,only: s_rv1, s_spec, s_site
    use m_suham,only: ndham=>ham_ndham !max dimension of hamiltonian +napwad (for so=0,2)
    use m_lmfinit, only: &
         ham_pwmode,pwemin,pwemax,ham_oveps,lrsig=>ham_lsig,nlmto,lso, &
         ham_scaledsigma,lat_alat,mxorb,nkaph,nsp,nspc,nl,mxorb,ssite=>v_ssite,sspec=>v_sspec,&
         nbas,n0,nppn,nkap0,slabl,nmcorex=>nmcore,iantiferro
    use m_lattic,only: lat_plat, lat_qlat,rv_a_opos
    use m_supot,only: lat_nabc, lat_gmax
    use m_rdsigm2,only: getsenex, senex,dsene
    use m_mksym,only: lat_npgrp,lat_nsgrp
    !      use m_shortn3,only: shortn3_initialize,shortn3
    use m_mkpot,only: &
         smpot=>osmpot, sv_p_osig, sv_p_otau, sv_p_oppi, vconst, vesrmt, &
         sv_p_osig, sv_p_otau, sv_p_oppi, ohsozz,ohsopm, ppn=>ppnl_rv, &
         sv_p_oppix,spotx!, spotxd, sv_p_oppixd
    use m_MPItk,only: numproc=>nsize,procid,master,master_mpi
    use m_igv2x,only: napw,ndimh,ndimhx,igv2x,m_Igv2x_setiq,ndimhall
    use m_elocp,only: rsmlss=>rsml, ehlss=>ehl
    use m_qplist,only: qplist,ngplist,ngvecp, iqini,iqend,ispini,ispend,iqibzmax      !for current rank iprocq,
    use m_hamindex0,only: Readhamindex0, nlindx
    use m_ftox
    use m_density,only: v0pot,pnuall,pnzall
   
    implicit none
    !! == Driver for fpgw (to prepare eigenfuncitons for fpgw) ==
    !! NOTE: following documents are not carefully examined. Not believe everything.
    ! i Inputs
    ! i   ssite,sspec,slat,sham :struct defined in m_struc_def.
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
    real(8):: rsml(n0), ehl(n0) ,eferm
    logical :: lwvxc,cmdopt0
    integer :: i,i1,i2,iat,ib,ibr,icore,ierr,ifeigen,&! & ifi,
         ifiqg,iflband(2),ifqeigen,ifsyml,igets,igetss,iix,iline, &
         im1,im2,ipb(nbas),ipqn,ipr,iprint,iq,is,isp,ispc,j,job,k1, &
         k2,k3,konf,konfig(0:n0),l,lchk,ldim,loldpw, &
         lmaxa,lmxax,lsig,mx,mxint,n1,n2,n3,nat, &
         ncore,ndima,nevl,nev,nglob,ngp,ngp_p, &
         ngpmx,nline,nlinemax,nlmax,nmx,nn1,nn2,nnn,npgrp, &
         nphimx,npqn,nqbz,nqibz,nqnum,nqnumx,nqtot,nr,nsgrp,iqibz,imx, &
         ifigwb,ifigwa,ifinormchk,ifigw1,ifildima,ifigwn,ifigwbhead, &
         ificlass,ifievec,ifievecx,ifigw2,ifiqbz,ifievv
    complex(8),allocatable :: aus_zv(:)
    real(8),allocatable :: ww_rv(:)
    real(8):: QpGcut_psi,QpGcut_cou,dum,dval, &
         xx(5),gmax,ecore(50),a,z,rmt(nbas),b,vshft, &
         alat,alfa,ef0,plat(3,3),qlat(3,3),qp(3),qpos,q_p(3), &! & ,qx(3)
         epsovl,dgets!, pnu(n0,2),pnz(n0,2)
    real(8),pointer:: pnu(:,:),pnz(:,:)
    integer,allocatable:: ips(:),ipc(:),ipcx(:),lmxa(:), &
         ngvecp_p(:,:) !takao feb2012 ,ngvecc(:,:) ,ngvecp(:,:)
    integer,allocatable :: konft(:,:,:),iiyf(:),ibidx(:,:),nqq(:)
    real(8),allocatable:: wk(:,:), &
         bas(:,:),rofi(:),rwgt(:),gcore(:,:,:),gval(:,:,:,:,:),evl(:,:),vxclda(:)
    real(8),allocatable::  cphin(:,:,:) !!ovv(:,:),evl_p(:,:), qq1(:,:),qq2(:,:),
    complex(8),allocatable:: ham(:,:),ovl(:,:),evec(:,:),vxc(:,:),&! & dipo(:,:,:),
         ppovl(:,:),phovl(:,:),pwh(:,:),pwz(:,:),pzovl(:,:), &
         testc(:,:),testcd(:),ppovld(:),cphi(:,:,:),cphi_p(:,:,:), &
         geig(:,:,:),geig_p(:,:,:),sene(:,:)
    integer::isize_ham(3)
    integer :: ltab(n0*nkap0),ktab(n0*nkap0),offl(n0*nkap0),norb
    integer :: nspx,nnnx,ifiv
    character strn*120
    integer :: pwmode
    double precision :: pwgmin,pwgmax,eomin

    ! ... for reading self-energy
    !      integer nqsig
    ! ... for band plotting
    !      real(8),allocatable :: ovvx(:,:,:)
    !      real(8) ::ovvpp
    !      integer idxdn(n0,nkap0)
    !      character lsym(0:n0-1)*1, lorb(3)*1, dig(9)*1, strn4*4
    !      data lsym /'s','p','d','f','g','5','6','7','8','9'/
    !      data lorb /'p','d','l'/
    !      data dig /'1','2','3','4','5','6','7','8','9'/
    integer:: w(1) !! dummy
    real(8):: dnn(3),qlatinv(3,3),qout(3),qtarget(3),qrr(3),axx,bxx,qxxx(3),qxx1(3),qxx2(3)
    integer:: inn(3),iqzz,nqzz,iiiii
    logical:: debug=.false.
    complex(8),allocatable:: evecout(:,:),evecr(:,:)
    real(8),allocatable:: qzz(:,:)
    !      integer,allocatable:: igv2(:,:)
    logical:: l_dummy_isanrg,isanrg,oncewrite,newsigmasw,sigmamode !,noshorbz
    integer::   irr,nqirr, nmcore !feb2012takao nk1,nk2,nk3,
    integer:: nbalance,nrr,ifiproc,iproc,iadd,ldw
    integer:: iqq
    logical:: rank0
    character*256:: ext,sprocid,extn
    logical ::nexist,magexist  !! june2013 magfield is added for fsmom mode.
    real(8):: vnow
    integer:: ifimag,ifile_handle !=9078
    integer,parameter:: noutmx=48
    integer:: iout,nout,nlatout(3,noutmx),iapw,ival
    real(8)::ppin(3)
    real(8):: rlatp(3,3),xmx2(3)
    real(8):: qqq(3)
    integer:: mode,iwdummy,jx
    complex(8),allocatable:: hamm(:,:,:),ovlm(:,:,:),ovlmtoi(:,:),ovliovl(:,:) ,hammhso(:,:,:)
!    integer,allocatable:: iantiferro(:)
    integer:: lpdiag=0,nrmx,ncoremx!,ndham
    character spid*8
    character(8) :: xt
    logical,optional:: socmatrix !dipolematrix,
    integer:: ispSS,ispEE
!#if MPI | MPIK
    include "mpif.h"
!#endif
    !! --- Setup ---
    call tcn ('m_sugw_init')
    !      stdo = lgunit(1)
    if(master_mpi) write(stdo,"(a)") 'm_sugw_init: start'
    sigmamode = mod(lrsig,10) .ne. 0
    !! june2013 magfield is added
    open(newunit=ifimag,file='MagField',status='old',err=112)
    read(ifimag,*,err=112) vnow
    close(ifimag)
    magexist=.true.
    write(6,"('Add mag.field to eval. -vnow/2 for isp=1, +vnow/2 for isp=2. vnow=',d13.6)")vnow
    goto 113
112 continue
    magexist=.false.
113 continue

    !!--------------------
    call getpr(ipr)
    alat=lat_alat
    n1 = lat_nabc(1)
    n2 = lat_nabc(2)
    n3 = lat_nabc(3)
    plat =lat_plat
    qlat =lat_qlat
    gmax =lat_gmax
    npgrp =lat_npgrp
    nsgrp =lat_nsgrp
    call fftz30(n1,n2,n3,k1,k2,k3)
    lchk = 1
    lwvxc = .not. cmdopt0('--novxc')
    !!  Count number of atoms : exclude floating orbitals
    nat = 0
    do  i = 1, nbas
       is = ssite(i)%spec
       lmaxa = int(sspec(is)%lmxa)
       if (lmaxa > -1) then
          nat = nat + 1
       endif
       ipb(i) = nat
    enddo
    if(nat/=nbas) call rx('sugw: gw is only for no floating orbital case') !2022jan

    !      lsig = 1
    ldim=nlmto
    pwmode= ham_pwmode
    !! ndima
    ndima = 0
    lmxax = -1
    do  ib = 1, nbas
       is=ssite(ib)%spec
       lmaxa=sspec(is)%lmxa
       !pnz  =ssite(ib)%pz
       pnz=>pnzall(:,1:nsp,ib)
       if(sum(abs(pnz(:,1:nsp)-pnzall(:,1:nsp,ib)))>1d-9) call rx('sugw xxx1aaa')
       lmxax = max(lmxax,lmaxa)
       if (lmaxa > -1) then
          do  l = 0, lmaxa
             npqn = 2
             if (pnz(l+1,1) /= 0) npqn = 3
             ndima = ndima + npqn*(2*l+1)
          enddo
       endif
    enddo
    !! ==== Read file NLAindx ====
    call Readhamindex0()
    open(newunit=ifiqg,file='QGpsi',form='unformatted')
    read(ifiqg ) nqnum, ngpmx ,QpGcut_psi,nqbz,nqirr,imx,nqibz
    !! ==== Write, or read past header information, file gwb ====
    sprocid='procid'//trim(xt(procid))
    if(master_mpi) open(newunit=ifigwa,file='gwa',form='unformatted')
    ef0 = 1d99        !dummy
    allocate(ips(nbas),lmxa(nat),bas(3,nat))!,iantiferro(nat))
    iat = 0
    do  i = 1, nbas
       lmaxa = sspec(int(ssite(i)%spec))%lmxa 
       if (lmaxa > -1) then
          iat = iat + 1
          if (iat > nat) call rx('bug in sugw')
          bas(:,iat)=rv_a_opos(:,i) !ssite(i)%pos
          lmxa(iat) = lmaxa
!          iantiferro(iat)=ssite(i)%iantiferro
       endif
    enddo
    print *,'iantiferro=',iantiferro
    !!   ... Determine nphimx
    nphimx = 0
    do  ib = 1, nbas
       is=ssite(ib)%spec
!       pnu=ssite(i)%pnu
!       pnz=ssite(i)%pz
!       if(sum(abs(pnz(:,1:nsp)-pnzall(:,1:nsp,i)))>1d-9) call rx('sugw xxx111')
!       if(sum(abs(pnu(:,1:nsp)-pnuall(:,1:nsp,i)))>1d-9) call rx('sugw xxx222')
       pnu=>pnuall(:,1:nsp,ib)
       pnz=>pnzall(:,1:nsp,ib)
       a=sspec(is)%a
       nr=sspec(is)%nr
       z=sspec(is)%z
       rmt(ib)=sspec(is)%rmt
       lmaxa = int(sspec(is)%lmxa)
       nmcore= nmcorex(is)
       if (lmaxa > -1) then
          call atwf(0,a,lmaxa,nr,nsp,pnu,pnz,rsml,ehl,rmt(ib),z,w,i1,ncore,konfig,ecore,w,w,nmcore)
          nphimx = max(nphimx,i1)
       endif
    enddo
    !! Atom data (gwa)
    write(stdo,*)' ... Generate core wave functions (file gwa)'
    ncoremx=0
    nrmx=0
    do  ib = 1, nbas
       is=ssite(ib)%spec
!       pnu=ssite(ib)%pnu
!       pnz=ssite(ib)%pz
       pnu=>pnuall(:,1:nsp,ib)
       pnz=>pnzall(:,1:nsp,ib)
       a=sspec(is)%a
       nr=sspec(is)%nr
       z=sspec(is)%z
       rmt(ib)=sspec(is)%rmt
       lmaxa = int(sspec(is)%lmxa)
       spid=slabl(is) !sspec(is)%name
       if (lmaxa > -1) then
          call atwf ( 0 , a , lmaxa , nr , nsp , pnu , pnz , rsml , ehl &
               , rmt ( ib ) , z , v0pot(ib)%v , i1 , ncore , konfig , ecore , w &
               , w ,nmcore)
          allocate(rofi(nr),rwgt(nr),gcore(nr,2,ncore))
          allocate(gval(nr,2,0:lmaxa,nphimx,nsp))
          gval=0d0
          !     !         Create augmented wave functions for this atom
          rsml=rsmlss(:,is)
          ehl= ehlss(:,is)
          call atwf ( 03 , a , lmaxa , nr , nsp , pnu , pnz , rsml , ehl &
               , rmt ( ib ) , z , v0pot(ib)%v , nphimx , ncore , konfig , ecore &
               , gcore , gval ,nmcore)
          if(nr     >nrmx   ) nrmx    = nr
          if(ncore  >ncoremx) ncoremx = ncore
          !     !         Header data for this atom
          b = rmt(ib)/(dexp(a*nr-a)-1d0)
          call radmsh(rmt(ib),a,nr,rofi)
          call radwgt(rmt(ib),a,nr,rwgt)
          if(master_mpi)write(ifigwa) z, nr, a, b, rmt(ib), lmaxa, nsp, ncore,spid
          if(master_mpi)write(ifigwa) konfig(0:lmaxa)
          if(master_mpi)write(ifigwa) rofi
          !     !         Write orthonormalized valence wave functions for this atom
          do  l = 0, lmaxa
             do  i = 1, nsp
                if(master_mpi)write(ifigwa) l,i
                if(master_mpi)write(ifigwa) gval(1:nr,1,l,1,i)
                if(master_mpi)write(ifigwa) gval(1:nr,1,l,2,i)
                if (konfig(l) >= 10 .AND. master_mpi) write(ifigwa) gval(1:nr,1,l,3,i)
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
                do  konf = l+1, mod(konfig(l),10)-1
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
                     call makrwf(0,z,rmt(ib),l,v0pot(ib)%v(1+(i-1)*nr:nr*i) &
                          ,a,nr,rofi,pnux,2,  gfun,gpfun,enu, phi,dphi,phip,dphip,p)
                     call gintxx(gfun(1:nr,1),gfun(1:nr,1),a,b,nr,sum1)
                     rgfun=[(rofi(ir)*gfun(ir,1),ir=1,nr)]
                     if(abs(sum1-1d0)>0.01d0) call rx('sugw: need to check normalization 111')
                     call gintxx(rgfun,gfun(1:nr,1),a,b,nr,sum2)
                     write(stdo,ftox)'\int phiphi : ib l isp=',ib,l,i, &
                          ' pnu \int rphiphi=',ftof(pnux,3),ftof(sum2,3)
                  enddo
                  write(stdo,*)
               enddo
               l = 3            !0, lmaxa
               do i = 1, nsp
                  do ie=-10,10
                     enu= eferm + 0.2d0*ie !enu is changing
                     call makrwf(1,z,rmt(ib),l,v0pot(ib)%v(1+(i-1)*nr:nr*i) &
                          ,a,nr,rofi,1d10,2,  gfun,gpfun,enu, phi,dphi,phip,dphip,p)
                     call gintxx(gfun(1:nr,1),gfun(1:nr,1),a,b,nr,sum1)
                     rgfun=[(rofi(ir)*gfun(ir,1),ir=1,nr)]
                     if(abs(sum1-1d0)>0.01d0) call rx('sugw: need to check normalization 222')
                     call gintxx(rgfun,gfun(1:nr,1),a,b,nr,sum2)
                     write(stdo,ftox)'\int phiphi : ib l isp=',ib,l,i, &
                          ' enu \int rphiphi=',ftof(enu,3),ftof(sum2,3)
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
       write(ifigwa) iantiferro(1:nat) !iantiferro may2015
       close(ifigwa)
    endif

    !     ndham= maxval(ndimhall) !some inconsistency if we assume this.
    !      See npwpad in suham.F I think npwpad=0 if essentially ok. but current version cause inconsistency.
    if(master_mpi) then
       open(newunit=ifigwbhead,file='gwb.head',form='unformatted')
       write(ifigwbhead) nat,nsp,ndima,ndham, maxval(lmxa(1:nat)),ncoremx/nsp,nrmx,plat,alat,nqirr,nqibz
       write(ifigwbhead) bas,lmxa,qplist,ngplist,ndimhall
       close(ifigwbhead)
    endif
    deallocate(ips,lmxa)
    allocate(evl(ndham,nsp),vxclda(ndham))

    !! == GW setup loop over k-points ==
    if (lchk>=1 ) then
       open(newunit=ifinormchk,file='norm.'//trim(sprocid)//'.chk')
       write(ifinormchk,"(a)") &
            '#     eval          IPW        IPW(diag)    Onsite(tot)   Onsite(phi)      Total ! lmfgw'
    endif
    !! --- Evecs and matrix elements of vxc for irr qp ---
    !!    Note: this routine should use only irr qp.
    !! == Main loop for eigenfunction generation ==
    if(ham_scaledsigma/=1d0 .AND. sigmamode) write(stdo,*)' Scaled Sigma method: ScaledSigma=',ham_scaledsigma
    do 1001 iq = iqini,iqend ! iqini:iqend for this procid
       !        if (dipolematrix.and.iq>nqbz) exit
       qp  = qplist(:,iq)     !  ... For this qp, G vectors for PW basis and hamiltonian dimension
       ngp = ngplist(iq)
       lwvxc = iq<=iqibzmax
       if (cmdopt0('--novxc')) lwvxc = .FALSE. 
       if (socmatrix) lwvxc = .TRUE. 
       call m_Igv2x_setiq(iq)! (qp)    ! Set napw ndimh ndimhx and igv2x
       allocate(ham(ndimh,ndimh),ovl(ndimh,ndimh),evec(ndimh,ndimh), vxc(ndimh,ndimh))
       allocate(cphi(ndima,ndimh,nsp),cphin(2,ndimh,nsp))
       !        if(dipolematrix) allocate(dipo(ndimh,ndimh,3))
       ispSS=1
       ispEE=nsp
       if(iq==iqini .AND. ispini==2) ispSS=2
       if(iq==iqend .AND. ispend==1) ispEE=1
       if(lso/=0 .OR. socmatrix) then
          allocate(hammhso(ndimh,ndimh,3))
          call aughsoc(qp, ohsozz,ohsopm, ndimh, hammhso)
       endif
       do 1002 isp = ispSS,ispEE
          open(newunit=ifigwb,file='gwb'//trim(xt(iq))//trim(xt(isp)),form='unformatted')
          if(lwvxc) then !lsig>0 .AND. ( .NOT. cmdopt0('--novxc'))) then
             open(newunit=ifievec,   file='evec'//trim(xt(iq))//trim(xt(isp)),form='unformatted')
             open(newunit=ifiv,      file='vxc'//trim(xt(iq))//trim(xt(isp)),form='unformatted')
             write(ifievec) ndimh, ldim
             write(ifiv)    ndimh, ldim
          endif
          !! note. This was intended for dipole but wrong since x,y,z are not periodic
          !          if(dipolematrix)then  !From spotd,and ppixd, dipo(=x,y,z)=<F_i|x,y,z|F_j> . Experimental.
          !            write(stdo,'("dipole matrix calculated")')    !is this is ambiguous for x+dx shifts but good for wannier?
          !            call hambl(isp,qp,spotxd(:,:,:,:,1),vconst,sv_p_osig,sv_p_otau,sv_p_oppixd(:,:,1),dipo(:,:,1),ovl)
          !            call hambl(isp,qp,spotxd(:,:,:,:,2),vconst,sv_p_osig,sv_p_otau,sv_p_oppixd(:,:,2),dipo(:,:,2),ovl)
          !            call hambl(isp,qp,spotxd(:,:,:,:,3),vconst,sv_p_osig,sv_p_otau,sv_p_oppixd(:,:,3),dipo(:,:,3),ovl)
          !          endif
          !! LDA Hamiltonian and overlap matrices for this qp ---
          call hambl(isp,qp,spotx,vconst,sv_p_osig, sv_p_otau,sv_p_oppix, vxc,ovl)!vxc=<F_i|H(LDA)-vxc(LDA)|F_j>
          call hambl(isp,qp,smpot,vconst,sv_p_osig, sv_p_otau, sv_p_oppi, ham,ovl)!ham=<F_i|H(LDA)|F_j> and ovl=<F_i|F_j>
          if(lso==2) ham(:,:) = ham(:,:) + hammhso(:,:,isp) !diagonal part of SOC matrix added for Lz.Sz mode.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !          qp  = qplist(:,iq)-[1d-6,2d-6,3d-6] !A trick to shift qp to avoid ambiguity of degeneracy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !          if(dipolematrix) then
          !             dipo(:,:,1)= dipo(:,:,1)-vxc ! dipo(i,j,1) = <Fi|x|Fj>
          !             dipo(:,:,2)= dipo(:,:,2)-vxc ! dipo(i,j,2) = <Fi|y|Fj>
          !             dipo(:,:,3)= dipo(:,:,3)-vxc ! dipo(i,j,3) = <Fi|z|Fj>
          !          endif
          vxc = ham - vxc ! vxc(LDA) part
          if (lwvxc) write(ifiv) vxc
          !          if (dipolematrix) write(ifiv) dipo,ovl
          if (socmatrix .AND. isp==1) write(ifiv) hammhso
          !! LDA + sigma Hamiltonian for this qp ---
          if(sigmamode) then ! See m_bandcal.F
             call getsenex(qp, isp, ndimh,ovl)
             ham(:,:) = ham(:,:) + ham_scaledsigma*senex
             call dsene() !delete senex
          endif
          !!   --- Branch jobgw = 1 : make cphi, matrix elements ---
          if (mod(iq,10) /= 1) call pshpr(iprint()-6)
          if (nspc == 2) call rx('diagonalization not ready for nspc=2')
          epsovl = ham_oveps
          !! we fined evl(i), z(:,i), i=1,nev !nov2015
          evec=-1d99
          evl(:,isp)=1d99
          call zhev_tk4(ndimh,ham,ovl,ndimh,nev,evl(1,isp),evec,epsovl)
          write(6,"(' sugw:  kpt isp=',i8,i2,' of ',i8, ' k= ',3f9.5, ' ndimh= ',i5, &
               ' irank=',i4, ' lwvxc=',l,' nev=',i5)")  iq,isp,nqnum,qp,ndimh,procid,lwvxc,nev
          call prtev(evec,ndimh,evl(1,isp),ndimh,nev)
          if (ndham>nev ) evl(1+nev:ndham,isp)=1d20 ! See rdata4gw_v2
          if(mod(iq,10) /= 1) call poppr
          if(debug) write(6,"(' sugw:procid iq isp lwvxc= ',3i3,' ',l)")procid, iq,isp,lwvxc
          if (lwvxc) write(ifievec) qp, evec(1:ndimh,1:ndimh),nev
          if(magexist) then
             if(isp==1) evl(1:ndimh,isp)=evl(1:ndimh,isp) - vnow/2d0
             if(isp==2) evl(1:ndimh,isp)=evl(1:ndimh,isp) + vnow/2d0
          endif
          !$$$  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          !$$$  !! wave function rotation test for si 5x5x5.
          !$$$  if(.false.) then
          !$$$  c       if(iq==8) then
          !$$$  ifievv=1012
          !$$$  open( ifievv, file='evecx' )
          !$$$  ifiqbz=1025
          !$$$  open(ifiqbz,file='QBZ')
          !$$$  read(ifiqbz,*) nqzz
          !$$$  allocate(qzz(3,nqzz))
          !$$$  do i=1,nqzz
          !$$$  read(ifiqbz,*) qzz(:,i)
          !$$$  enddo
          !$$$
          !$$$  allocate(evecout(ndimh,ndimh))
          !$$$  do iqzz = 1,nqzz
          !$$$  qtarget = qzz(:,iqzz)
          !$$$  c
          !$$$  c                call rotwv(q,qtarget,ndimh,napw,ndimh, plat,qlat,evec,evecout,ierr)
          !$$$  call rotwv(q,qtarget,ndimh,napw,ndimh, evec,evecout,ierr)
          !$$$  if(ierr/=0) cycle
          !$$$
          !$$$  do
          !$$$  read(ifievv,*,end=1019) iiiii
          !$$$  read(ifievv,*) qrr, ndimh
          !$$$  print *,'qrr ndimh=',qrr,ndimh
          !$$$  allocate(evecr(ndimh,ndimh))
          !$$$  print *,'qrr ndimh xxx=',qrr,ndimh
          !$$$  do j= 1,ndimh
          !$$$  do i= 1,ndimh
          !$$$  read(ifievv,*) axx,bxx
          !$$$  evecr(i,j)=dcmplx(axx,bxx)
          !$$$  enddo
          !$$$  enddo
          !$$$  if(sum(abs(qrr-qtarget))<1d-8) exit
          !$$$  deallocate(evecr)
          !$$$  enddo
          !$$$  rewind ifievv
          !$$$  c$$$                if(.not.noshorbz()) then
          !$$$  c$$$                  call shorbz(qtarget,qxxx,qlat,plat)
          !$$$  c$$$                else
          !$$$  c$$$                  qxxx=qtarget
          !$$$  c$$$                endif
          !$$$  qxxx=qtarget
          !$$$  write(1013,"(i10)") 11111
          !$$$  write(1013,"(3f8.3,i10)") qtarget,ndimh
          !$$$  c           write(1013,"(i10,3f8.3,i10)") 11111,qout,ndimh
          !$$$  do j=1,ndimh
          !$$$  do i=1,ndimh
          !$$$  if(abs(evecout(i,j))+abs(evecr(i,j))>1d-4) then
          !$$$  write(1013,"(2i4,2d13.5,2x,2d13.5,2x,2d12.4,2x,2d12.4,3f8.3)")
          !$$$  &                i,j,evecout(i,j), evecr(i,j), evecout(i,j)/evecr(i,j),
          !$$$  &                abs(evecout(i,j)), abs(evecr(i,j)),qxxx
          !$$$  endif
          !$$$  enddo
          !$$$  write(1013,*)
          !$$$  enddo
          !$$$
          !$$$  deallocate(evecr)
          !$$$  enddo
          !$$$  deallocate(evecout)
          !$$$  stop 'test end xxxxxxxx'
          !$$$  1019         continue
          !$$$  stop 'uuuuuuuuuuuuuuu'
          !$$$  endif
          !$$$  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

          !     !   --- Project wf into augmentation spheres, Kotani conventions ---
          !     nlmax = globalvariables%mxorb / globalvariables%nkaph
          nlmax = mxorb / nkaph
          allocate(aus_zv(nlmax*ndham*3*nsp*nbas))
          aus_zv(:)=0.0d0
          call makusq ( 1 ,  nbas,0, nev,  isp, 1 , qp , evec , aus_zv )
          call gwcphi (sspec , isp , nsp , nlmax , ndham , nev &
               , nbas , ipb , lmxax , nlindx , ndima , ppn , aus_zv , cphi &
               ( 1 , 1 , isp ) , cphin ( 1 , 1 , isp ) )
          if (allocated(aus_zv)) deallocate(aus_zv)
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
          if(debug) print *,'ppppppppppp ngp=',ngp
          if(ngp > 0) then
             allocate(ppovl(ngp,ngp),pwz(ngp,ndimh))
             allocate(phovl(ngp,ndimh))
             !     Pass qx to pwmat (or pwmat2):
             !     qx = (unshortened) q if 2s digit loldpw = 0
             !     qx = (shortened)  qp if 2s digit loldpw = 1
             !     Old convention: call pwmat
             !     if (mod(loldpw,2) .eq. 0) then

             !     !  We have  q+G(igvx; internal in pwmat) = qp + G(igv2)
             !     !  Thus, we have
             !     !           igv(internally in pwmat) = igv2 + qlatinv*(qp-q)
             !     !            inn = qlatinv*(qp-q)
!!!!  NOTE: from 20Sep2012, we set qp=q (shorbz is in hambl.F. In other words, q supplied to hambl
!!!!  do not need to be short enough).
             !             if(sum(abs(qp-q))>1d-8) stop 'sugw:qp/=q; qp=p sep2012'
             inn=0
             call pwmat (  ssite , sspec , nbas , ndimh , napw,&
                  igv2x, qp , ngp , nlmax , ngvecp(1,1,iq) , gmax , inn, ppovl, phovl )
             pwz=matmul(phovl,evec)
             ! all zgemm('N','N',ngp,ndimh,ndimh,(1d0,0d0),phovl,ngp, evec,ndimh,(0d0,0d0),pwz,ngp)
             if(debug) print *,'sss: pwz =',sum(abs(pwz))
             if(debug) print *,'sss:       '
             deallocate(phovl)
             if (lchk >= 1) then
                allocate(pzovl(ngp,ndimh))
                pzovl = pwz
                allocate(ppovld(ngp)) ! extract diagonal before ppovl overwritten
                do  i = 1, ngp
                   ppovld(i) = ppovl(i,i)
                enddo
             endif
             !     ! inversion of hermitian ppovl
             call matcinv(ngp,ppovl)
             pwz = matmul(ppovl,pwz)
             deallocate(ppovl)
          endif
          if(debug) write (6,"('q ndimh=',3f10.5,i10)") qp, ndimh
          !          write(ifigw1) qp, (evl(i,isp),i=1,ndimh)
          !  call prmx('e(H+sigma-vxc)',evl,ndimh,ndimh,1)
          !     ... Overlap checking.   Define:
          !     Interstitial part of eigenfunction overlap:
          !     <psi_n||psi_n'>
          !     = sum_G1,G2 (pwz_G1,n|IPW_G1>)+  (pwz_G2,n'|IPW_G2>)
          !     = sum_G1,G2 (pwz_G1,n)+ ppovl_G1,G2 (PWZ_G2,n')
          !     = (PWZ)+ O (PWZ) = (PZOVL)+ (PWZ)  (old style)
          if (lchk >= 1 .AND. ngp > 0) then
             allocate(testc(ndimh,ndimh),testcd(ndimh))
             testc=matmul(transpose(dconjg(pzovl)),pwz)
             !             call zgemm('C','N',ndimh,ndimh,ngp,(1d0,0d0),pzovl,ngp,pwz,ngp,(0d0,0d0),testc,ndimh)
             deallocate(pzovl)
             !     ! call zprm('(PWZ)+^-1 O_i^-1 (PWZ)',2,testc,ndimh,ndimh, ndimh)
             do  i = 1, ndimh
                testcd(i) = sum(dconjg(pwz(:,i))*ppovld*pwz(:,i))
             enddo
             deallocate(ppovld)
             !     xx(1) = sum over all augmentation w.f.  cphi+ ovl cphi
             !     xx(2) = sum over augmentation phi only.
             !     xx(3) = IPW contribution to phi+ phi
             !     xx(4) = IPW contribution to phi+ phi, using diagonal part only
             !             if (abs(sum(q-qp)) .gt. 1d-10) then
             !                write(ifinormchk,555) iq,qp,qp
             !             else
             write(ifinormchk,555) iq,qp
             !             endif
555          format('# iq',i5,'   q',3f12.6:'  shortened q',3f12.6)
             do  i1 = 1, ndimh
                xx(1) = cphin(1,i1,isp)
                xx(2) = cphin(2,i1,isp)
                do  i2 = 1, ndimh
                   !     xx(1) = sum(dconjg(cphi(1:ndima,i1,isp))*cphi(1:ndima,i2,isp))
                   !     xx(2) = sum(dconjg(cphi(1:nchan,i1,isp))*cphi(1:nchan,i2,isp))
                   xx(3) = testc(i1,i2)
                   xx(4) = testcd(i1)
                   if (i1==i2) write(ifinormchk,'(f12.5,5f14.6)') &
                        evl(i1,isp),xx(3),xx(4),xx(1),xx(2),xx(1)+xx(3)
                enddo
             enddo
             deallocate(testc,testcd)
             write(ifinormchk,*)
          endif                 !------------- end of lchk=1
          allocate(testc(ndimh,ndimh))
          !          call zgemm('C','N',ndimh,ndimh,ndimh,(1d0,0d0),evec,ndimh,vxc,ndimh,(0d0,0d0),testc,ndimh)
          testc=matmul(transpose(dconjg(evec)),vxc)
          do i1 = 1, ndimh
             vxclda(i1) =  sum(testc(i1,1:ndimh) * evec(1:ndimh,i1))  !<i|Vxc^lda|i>
          enddo
          write(ifigwb) evl(1:ndimh,isp),cphi(:,:,isp),pwz,vxclda(1:ndimh),nev
          deallocate(testc,pwz)
          close(ifigwb)
          if ( .NOT. cmdopt0('--novxc')) then
             close(ifiv)
             close(ifievec)
          endif
          continue             ! Loop over spins
1002   enddo
       if(allocated(hammhso)) deallocate(hammhso)
       !       if(allocated(dipo)) deallocate(dipo)
       deallocate(ham,ovl,evec,vxc,cphi,cphin)
       continue ! ===== Loop over qp ===============================
1001 enddo
    deallocate(evl)
    close(ifiqg)
    call tcx('m_sugw_init')
  end subroutine m_sugw_init
end module m_sugw

!$$$      subroutine ioaindx(npqn,lmxax,nbas,ndima,nlindx)
!$$$C-  read NlAindx
!$$$C ----------------------------------------------------------------------
!$$$Ci Inputs
!$$$Ci   npqn  :leading dimension of nlindx
!$$$Ci   lmxax :second dimension of nlindx
!$$$Ci   nbas  :size of basis
!$$$Ci   ndima :number of augmentation channels
!$$$Ci   nlindx:pointer to augmentation channels by site
!$$$      implicit none
!$$$      integer npqn,lmxax,nbas,ndima,ifi,nlindx(npqn,0:lmxax,nbas)
!$$$      character outs*80
!$$$      integer i,ipqn,l,ib,ii
!$$$      open(newunit=ifi,file='NLAindx')
!$$$      nlindx = -1
!$$$      read(ifi,'(a)') outs
!$$$      read(ifi,*) i
!$$$      if (ndima .gt. 0 .and. i .ne. ndima)
!$$$     .  call rx('ioaindx: file mismatch NLAindx')
!$$$      ndima = i
!$$$      do  i = 1, ndima
!$$$          read(ifi,'(a)',err=101,end=101) outs
!$$$          read(outs,*) ipqn,l,ib,ii
!$$$          nlindx(ipqn,l,ib) = ii
!$$$       enddo
!$$$ 101   continue
!$$$      close(ifi)
!$$$      end subroutine ioaindx
