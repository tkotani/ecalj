subroutine wmatqk_mpi(kount,irot,nrws1,nrws2,nrws,  tr, iatomp, &
     rws1,rws2, nsp,isp, &! & ifcphi jan2004,ifrb,ifcb,ifrhb,ifchb,
     ifrcw,ifrcwi, qbas,ginv, qibz,qbz,wk,nstbz,wik,nstar,irk, & ! & koun,,iindxk
     iclass,mdim,nlnmv,nlnmc,  icore,ncore,imdim, &
     ppb, freq_r,freqx,wx,expa,ua,dw,& ! & deltaw,freq
     nlmto,nqibz,nqbz,nctot, &
     nl,nnc,natom,natomx, &
     nlnmx,mdimx,nbloch,ngrp,nw_i,nw,nrw,niw,niwx,nq, &
     nblochpmx ,ngpmx,ngcmx, &! & ngveccBr,!Jan2004
     wgt0,wqt,nq0i,q0i,symope,alat, shtv,nband, &! ifvcfpout, &
     exchange, pomatr, qrr,nnr,nor,nnmx,nomx,nkpo, nwf,  rw_w,cw_w,rw_iw,cw_iw)
  use m_zmel_old,only: drvmelp3
  use m_ftox
  use m_readqg,only: readqg0
  use m_readeigen,only:readcphiw
  use m_keyvalue,only: getkeyvalue
  use m_read_bzdata,only: wklm
!
!  use rsmpi_rotkindex,only:nk_local_rotk,ik_index_rotk
  implicit none
  integer :: ntq, natom,nqbz,nqibz,ngrp,nq,nw_i,nw,niw, natomx,&
       nband,  nlmto, nq0i,nctot,mbytes,iwksize,nlmtobnd,nstate,nstatex, &
       irot,  iqisp,ikpisp,isp,nsp,  nlnmx, iq, ip, it,itp, it2, itp2,  iiclass,mdim(*), &
       ifrcw,ifrcwi, ndummy1,ndummy2,kx,kr,kr2,kr3,ngc,ngb,nbloch, &
       kp,nt0,nocc, nt0p,nt0m,irkp,i,nt0org,nmax,nt,ntp0, &
       nbmax,nl,nnc, nblochpmx,ix,nx,iw,iwp,ixs,ixsmx, &
       mdimx, nwx,niwx, iatomp(natom), &
       nstar(nqibz),irk(nqibz,ngrp),kount(nqibz),nwf !,iclose
  real(8) :: q(3),qbas(3*3),ginv(3*3),tr(3,natom), &
       wk(nqbz),wik(nqibz),qibz(3,nqibz),qbz(3,nqbz), &
       freqx(niw),wx(niw),expa(niw), &
       eq(nband), &
       !     &   ekq(nband), ekc(nctot+nband),
       tpi,ef,ef2,esmr,esmr2,efp,efm,wtx,wfac,wfacx,we,esmrx,ua, &
       dw,wtt,wexx,www,exx,exxq,weight
  integer:: ngpmx, ngcmx,  igc, nadd(3)
  real(8) :: wgt0(nq0i,ngrp),wqt(nq0i),qk(3), qbasinv(3,3),qdiff(3),add(3),symope(3,3), &
       qxx(3),q0i(1:3,1:nq0i),shtv(3),alat, & !,ecore(nctot), &
       ppb(1) !pdb(1),dpb(1),ddb(1)
  real(8),allocatable:: rmelt(:,:,:),cmelt(:,:,:), rmelt2(:,:,:),cmelt2(:,:,:), &
       rmelt3(:,:,:,:),cmelt3(:,:,:,:)
  complex(8),allocatable :: zz(:),zmel(:,:,:),zzmel(:,:,:), &
       zw (:,:), zwz(:,:,:), zwz0(:,:),zwzi(:,:),zwz00(:,:), &
       zmelt(:,:,:),zmelc(:,:,:,:),zmelcc(:,:,:,:)
  logical :: exchange,screen,cohtest,tote
  real(8),allocatable:: &
       w1p(:,:,:),w2p(:,:,:)
  complex(8),allocatable :: z1p(:,:,:),vcoul(:,:),vcoult(:,:)
  logical :: debug=.false.
  integer :: ibl,iii,ivsumxxx,ifexsp 
  integer,save::ifzwz=-999
  integer :: iwini, iwend, ia
  real(8) :: rw_w(nwf,nwf,nwf,nwf,nrws,0:nrw), &
       cw_w(nwf,nwf,nwf,nwf,nrws,0:nrw), &
       rw_iw(nwf,nwf,nwf,nwf,nrws,niw), &
       cw_iw(nwf,nwf,nwf,nwf,nrws,niw)
  complex(8),allocatable:: expikt(:)
  complex(8):: img=(0d0,1d0)
  complex(8):: cphiq(nlmto,nband), cphikq(nlmto,nband)   , cphiqtmp(nlmto,nband)
  integer :: nt_max, igb1,igb2,iigb, nw_w
  complex(8),allocatable:: zmel1(:)
  complex(8), allocatable :: zw_(:,:) !,zzmel(:,:)
  complex(8), allocatable :: zwz2(:,:),zw2(:,:,:,:) !0 variant
  complex(8) ::  zz2 ,zwz3(3)
  real(8) :: dd,omg_c,dw2
  real(8) :: freq_r(nw_i:nw)
  complex(8), allocatable :: zw3(:,:,:)
  real(8)::weavx,wfaccut=1d-10
  logical :: GaussSmear
  real(8) :: ddw !ebmx,
  integer:: nbmxe,nstatetot !nbmx,

  !      integer:: n_index_qbz
  !      integer:: index_qbz(n_index_qbz,n_index_qbz,n_index_qbz)

  integer::nlnmv(*),nlnmc(*),iclass(*),icore(*),ncore(*),imdim(*)

  integer::verbose,nstbz(nqbz),iqini,iqend !bzcase,
  real(8):: wgtq0p

  integer:: iqindx,nrec,kxx
  real(8)::quu(3),qibz_k(3),qbz_kr(3)
  logical :: onlyQ0P, onlyimagaxis ,noq0p !,noq0p,test_omitq0p,

  logical ::zwz3mode
  !      logical ::testimx=.false.

  real(8):: ua_,expa_(niw),ua2,freqw,freqw1,ratio,ua2_(niw)
  logical :: ua_auto
  integer:: icc=0
  real(8),allocatable:: uaa(:,:)

  !      logical ::testimx=.false.
  ! ccc zvz test cccccccccccccccccccccccccc
  integer:: ngbx
  !      complex(8):: vcoul(ngbx,ngbx)
  complex(8),allocatable:: vzz(:,:,:),aaa(:)
  complex(8):: zvz,zvz1
  integer:: ib1,ib2,ifix
  ! ccccccccccccccccccccccccccccccccc
  integer ::nbcut,nbcutc
  logical ::iww2=.true., oncew


  !...
  logical::smbasis
  integer:: nn,no,ifpomat,isx,iqx
  complex(8),allocatable:: pomat(:,:)
  real(8):: q_r(3)
  integer:: nnmx,nomx,nkpo, nnr(nkpo),nor(nkpo)
  complex(8):: pomatr(nnmx,nomx,nkpo)
  real(8):: qrr(3,nkpo)

  real(8):: elxx,ehxx,ekxx,efxx
  integer:: ixsmin,iwm,iir,nwxi
  real(8)   :: fffr(3)
  complex(8):: zwzz(3)
  integer :: nqbz2,nwf2,iko_ix,iko_fx,iqtmp,ifmlw,nko,iqk &
       ,ifi,in1,in2,imp,ilp,ii,jj,nrws,nrws1,nrws2 &
       ,ir1,ir2,ir3,ir,nrw
  real(8) :: norm2,qtmp(3),rws1(3,nrws1),rws2(3,nrws2),tmp
  complex(8) :: ztmp,expiqR1(nrws1),expiqR2
  complex(8),allocatable :: cnk(:,:,:),zmel2(:,:,:),zmel3(:,:,:)
  integer :: itq(nwf)
  complex(8) :: weightc(nrws1),zmeltt1
  complex(8),allocatable:: ppovl(:,:),ppovlz(:,:),zcousq(:,:),zmeltt(:,:,:)
  real(8),allocatable::vcoud(:),vcousq(:)
  real(8),parameter:: pi=4d0*atan(1d0),fpi=4d0*4d0*atan(1d0)
  integer:: ifvcoud,ivc,ngb0
  real(8)::qvv(3),vc
  character(10):: i2char
  complex(8)::w3p
  integer:: mrecl,nprecx,nwordr,il
  integer :: kx_local
  debug=.false.
  if(verbose()>=90) debug= .TRUE. 
   write(6,ftox)' nnnnnnnnnn wmatqk_mpi: nrws nrws1 nrws2       ',nrws,nrws1,nrws2
  call getkeyvalue("GWinput","nbcutlow_sig",nbcut, default=0 )
  nbcutc=nctot+nbcut
  tpi         = 8d0*datan(1.d0)
  nlmtobnd    = nlmto*nband
  nstatetot      = nctot + nband
  call minv33(qbas,qbasinv)
  if(debug) write(6,*) ' sxcf: 1'
  allocate(expikt(natom))
  if(abs(sum(qibz(:,1)**2))/=0d0) call rx( ' sxcf assumes 1st qibz/=0 ')
  if(abs(sum( qbz(:,1)**2))/=0d0) call rx( ' sxcf assumes 1st qbz /=0 ')
  do it = 1,nwf
     itq(it) = it
  enddo
  kx = 1  ! qibz(:,1)=0 contribution for kcount
  if(irk(kx,irot)/=0) kount(kx)= kount(kx) + 1
  iqini=2
  iqend=nqibz+nq0i
  iqini=1
  iqend=nqibz            !no sum for offset-Gamma points.
  call getkeyvalue("GWinput","TestOnlyQ0P",onlyq0p,default=.false.)
  call getkeyvalue("GWinput","TestNoQ0P",noq0p,default=.false.)
  if ( .NOT. noq0p) &
       call getkeyvalue("GWinput","NoQ0P",noq0p,default= .FALSE. )
  if(noq0p)write(*,*)'noq0p mode'
  if(noq0p) iqend=nqibz
  kxloop: do 1100 kx=1,nqibz !kx_local = 1,nk_local_rotk(irot)
    !     kx = ik_index_rotk(irot,kx_local)
    !     write(6,*) ' kkkkk sxcf: kkkkk goto loop kx=',kx
    if( kx <= nqibz ) then
      kr = irk(kx,irot) ! index for rotated k in the FBZ
      qibz_k= qibz(:,kx)
      if(kr/=0) qbz_kr= qbz (:,kr) !feb2006
    else
      kr=-99999 !for sanity check
      qibz_k= 0d0
      qbz_kr= 0d0
    endif
    call readqg0('QGcou',qibz_k,  quu,ngc)
    ngb = nbloch + ngc
    !! ===Readin diagonalized Coulomb interaction===
    !! note sep102012takao
    !!  Vcoud file is sequential file Vcoulomb matrix for qibz_k.
    !!  A possible choice for paralellization is "Vcoud.ID" files where ID=kx
    !!  Vould file is written in hvccfp0.m.F.
    !! For correlation, W-v is read instead of Vcoud file (ifrcw,ifrcwi for WVR and WVI)
    !! These can be also separeted into WVR.ID and WVI.ID files.
    qxx=qibz_k
    !           if(kx<=nqibz) qxx=qibz_k
    !           if(kx>nqibz ) qxx=q0i(:,kx-nqibz)
    open(newunit=ifvcoud,file=trim('Vcoud.'//i2char(kx)),form='unformatted')
    do
      read(ifvcoud) ngb0
      read(ifvcoud) qvv
      write(6,"('readin qvv ngb0=',3f9.4,5i5)")qvv,ngb0,ngc,ngb,nbloch
      !              write(6,"('readin qxx ngb0=',3f9.4,i5)")qxx
      if(allocated(vcoud)) deallocate(vcoud)
      allocate( zcousq(ngb0,ngb0),vcoud(ngb0) )
      read(ifvcoud) vcoud
      read(ifvcoud) zcousq
      if(sum(abs(qvv-qxx))<1d-6) goto 1133
    enddo
    if(sum(abs(qvv-qxx))>1d-6) then
      write(6,*)'qvv =',qvv
      write(6,*)'qxx=',qxx,kx
      call rx( 'wmatK: qvv/=qibz(:,kx) hvcc is not consistent')
    endif
1133 continue
    close(ifvcoud)! = iclose('Vcoud.'//i2char(kx))
    if( ngb0/=ngb ) then !sanity check
      write(6,*)' qxx ngb0 ngb=',qxx,ngb0,ngb
      call rx( 'hsfp0.m.f:ngb0/=ngb')
    endif
    !! <I|v|J>= \sum_mu ppovl*zcousq(:,mu) v^mu (Zcousq^*(:,mu) ppovl)
    !! zmel contains O^-1=<I|J>^-1 factor. zmel(phi phi J)= <phi_q,itp |phi_q-rk,it B_rk,I> O^-1_IJ
    !! ppovlz= O Zcousq
    !! (V_IJ - vcoud_mu O_IJ) Zcousq(J, mu)=0, where Z is normalized with O_IJ.
    if(allocated(ppovlz)) deallocate(ppovlz)
    allocate(ppovl(ngc,ngc),ppovlz(ngb,ngb))
    call readppovl0(qibz_k,ngc,ppovl)
    ppovlz(1:nbloch,:) = zcousq(1:nbloch,:)
    ppovlz(nbloch+1:nbloch+ngc,:) = matmul(ppovl,zcousq(nbloch+1:nbloch+ngc,:))
    deallocate(zcousq,ppovl)
    if (kr == 0) cycle
    if(OnlyQ0P .AND. kx<=nqibz) then
      if(exchange) deallocate(vcoul)
      cycle
    endif
    !! phase factor for off-site W
    do ir1=1,nrws1
      expiqR1(ir1) = exp(-img*tpi* sum(qbz_kr(:)*rws1(:,ir1)))
      !        write(6,*)'nnnnnnnniirrr ',qbz_kr(:),' vvvvv',rws1(:,ir1)
    enddo

    !! ===================================================================
    write(6,*)'nnnnnnnnnnnn at 735',ngb,nwf,nrws2,sum(abs(expiqR1))
    allocate( rmelt3(ngb,nwf,nwf,nrws2),cmelt3(ngb,nwf,nwf,nrws2))
    !     write(6,*)'nnnnnnnnnnnn at ssss'
    rmelt3 = 0d0
    cmelt3 = 0d0
    !! loop over FBZ
    do iq = 1,nqbz
      q(:) = qbz(:,iq)
      call readcphiW (qbz(:,iq), nlmto,isp, quu, cphiq)
      qk =  q - qbz_kr          ! qbz(:,kr)
      call  readcphiW(qk, nlmto,isp, quu, cphikq)
      do ia = 1,natom
        expikt(ia) = exp(img*tpi* sum(qibz_k*tr(:,ia)) ) !  write(6,'(" phase ",i3,2d12.4)')ia,expikt(ia)
      end do
      if(debug) write(6,*) ' sxcf: tr=',tr
      if(debug) write(6,*) ' sxcf: goto psicb2'
      nbmax = nwf
      nt   = nctot + nbmax      ! = nstate for the case of correlation
      ntp0 = nwf
      allocate( zzmel(nbloch,nt,ntp0)) ! rk,ibloch  q-rk,it  q,itp
      zzmel = 0d0
      call psi2b_v2 (nbmax, ntp0, iclass, &
           dreal(expikt(1:natom)),dimag(expikt(1:natom)), &
           cphikq,      &        ! & rbkq,cbkq,rhbkq,chbkq, !  q-rk nko
           cphiq,            &   ! & rbq,cbq,rhbq,chbq,     !  q    nko
           ppb,              &   ! & pdb,dpb,ddb,
           nlnmv,nlnmc,mdim,nctot, &
           imdim,iatomp, &
           mdimx,nlmto,nbloch,nlnmx, nband, nt,ntp0, &
           natom,natom, &
           zzmel)               ! rmel,cmel)
      if(debug) write(6,"('sum of zmel abszmel=',4d23.16)") &
           sum(zzmel),sum(abs(zzmel) )
      allocate( zmelt(ngb, nctot+nbmax, ntp0) ) ! rk,ibloch  q-rk,it  q,itp
      if(debug) write(6,*) ' sxcf_fal2: goto drvmelp xxxxx1',ngb,nctot,nbmax,ntp0
      call drvmelp3( q,   ntp0, &! &  q in FBZ  q,itp
           q-qbz_kr, nbmax,  &! &  q-rk,it
           qibz_k,           &! &  k in IBZ for e-product basis
           isp,ginv, &
           ngc,ngcmx,ngpmx,nband,itq, &
           symope, shtv, qbas, qbasinv,qibz,qbz,nqbz,nqibz, &
           dreal(zzmel), dimag(zzmel), nbloch, nt,nctot, &
           zmelt)
      deallocate(zzmel) !rmel,cmel)
      if(debug) write(6,*) ' sxcf: goto wtt'
      if(debug) write(6,"('sum of rmelt cmelt=',4d23.16)") sum(zmelt)
      do ir2=1,nrws2
        expiqR2 = exp( img*tpi* sum(q(:)*rws2(:,ir2)))
        rmelt3(:,:,:,ir2) = rmelt3(:,:,:,ir2) + wk(iq) * dreal(zmelt(:,:,:)*expiqR2)
        cmelt3(:,:,:,ir2) = cmelt3(:,:,:,ir2) + wk(iq) * dimag(zmelt(:,:,:)*expiqR2)
      enddo  ! ir2
      deallocate(zmelt)
    enddo
    !! ===================================================================



    if(kx<= nqibz) then
      wtt = wk(kr)           !         wtx = 1d0
    else
      wtt = wk(1)*wgt0(kx-nqibz,irot) ! wtx = wgt0(kx-nqibz,irot)
      if(abs(wk(1)-1d0/dble(nqbz))>1d-10) stop 'sxcf:wk(1) inconsistent'
    endif
    weight = wtt
    if(debug) then
      write(6,"(' kx wtt=',i4,f12.8)") kx,wtt
    endif
    do ir1=1,nrws1
      weightc(ir1) = weight*expiqR1(ir1)
    enddo
    !--------------------------------------------------------
    ! --- bare Coulomb section ---
    !--------------------------------------------------------
    if(exchange) then
      if (debug) write(*,*) 'bare coulomb section begins'
      allocate(zmel1(ngb))
      allocate(zmel(ngb, nwf, nwf))
      !        if( .NOT. newansisoW()) allocate(vcoult(1:ngb,1:ngb),z1p(ngb,nwf,nwf))
      do ir2=1,nrws2
        !!  (rmelt3,cmelt3)  !rk,ibloch  q-rk,it  q,itp
        zmel  = dcmplx (rmelt3(:,:,:,ir2),cmelt3(:,:,:,ir2)) !<psi_itp|psi_it B>
        ! based on E_I basis. See Christoph's paper
        allocate(zmeltt(nwf,nwf,ngb))
        do itp= 1,nwf
          do it = 1,nwf
            do ivc=1,ngb
              zmeltt(it,itp,ivc)=sum(zmel(:,it,itp)*ppovlz(:,ivc)) ! <psi_itp|psi_it E_I> (I=ivc)
            enddo
          enddo
        enddo
        do ir3=1,nrws2
          do itp2 = 1,nwf
            do it2  = 1,nwf
              do it   = 1,nwf
                do itp  = 1,nwf
                  zmel1(:)=dcmplx(rmelt3(:,it,itp,ir3),cmelt3(:,it,itp,ir3)) ! <psi_itp|psi_it B_I>
                  w3p=0d0
                  do ivc=1,ngb
                    zmeltt1 =  sum( zmel1(:)*ppovlz(:,ivc) ) !<psi_itp|psi_it E_I>
                    if(ivc==1 .AND. kx==iqini) then
                      vc= wklm(1)* fpi*sqrt(fpi) /wk(kx) !kx right?
                    else
                      vc= vcoud(ivc)
                    endif
                    w3p=w3p + zmeltt(it2,itp2,ivc)  *vc* dconjg(zmeltt1)
                    ! <psi_itp2|psi_it2 E_I> *vc* <E_I psi_it|psi_itp>
                  enddo
                  ztmp= w3p
                  do ir1=1,nrws1
                    ir = ir1 + (ir2-1 + (ir3-1)*nrws2)*nrws1
                    rw_w(itp2,it2,it,itp,ir,0) = rw_w(itp2,it2,it,itp,ir,0) + dreal(ztmp*weightc(ir1))
                    cw_w(itp2,it2,it,itp,ir,0) = cw_w(itp2,it2,it,itp,ir,0) + dimag(ztmp*weightc(ir1))
                  enddo ! ir1
                enddo
              enddo
            enddo
          enddo
        enddo ! ir3
        deallocate(zmeltt)
        !endif
      enddo ! ir2
      if(allocated(vcoul)) deallocate(vcoul,vcoult)
      if(allocated(z1p))   deallocate(z1p)
      if(allocated(rmelt3)) deallocate(rmelt3,cmelt3,zmel1, zmel)
      if (debug) write(*,*) 'bare coulomb section finished'
      !! --- End of bare-Coulomb section --------------
    else
      !--------------------------------------------------------------------------
      !--- screening effect section----------------------------------------------
      !--------------------------------------------------------------------------
      ! S[i,j] <psi(q,t) |psi(q-rk,n) B(rk,i)>
      !                Wc(k,0)(i,j) > <B(rk,j) psi(q-rk,n') |psi(q,t')>
      ! t-->itp   n -->it
      ! t'-->itp2 n'-->it2
      !--------------------------------------------------------------
      !!--- The matrix elements zmelc.
      !! zmelc  = < E(rk,j) psi(q-rk,it) | psi(q,itp) >
      !! E basis is ghe Christoph's basis diagonalize the Coulomb inteaction
      allocate( zmelc (ngb, nwf, nwf,nrws2),zmelcc (ngb, nwf, nwf,nrws2), &
           zw (nblochpmx,nblochpmx), &
           zw2(nwf,nwf,nwf,nwf) )
      zmelcc = dcmplx (rmelt3,-cmelt3) !zmelcc = < B(rk,j) psi(q-rk,it) | psi(q,itp) >
      do it =1,nwf
        do itp=1,nwf
          do ir2=1,nrws2
            zmelc(:,it,itp,ir2) =  matmul(zmelcc(:,it,itp,ir2),dconjg(ppovlz(:,:)))
          enddo
        enddo
      enddo
      deallocate(rmelt3,cmelt3,zmelcc)
      if(debug) write(6,*)' end of zmel'
      !====================================================================
      ! Wc(qt,w) along the imaginary axis
      !====================================================================
      !------------------------------------------------
      ! loop over w' = (1-x)/x, frequencies in Wc(k,w')
      ! {x} are gaussian points between (0,1)
      !------------------------------------------------
      nx  = niw
      nprecx=8
      mrecl  = nprecx*2*nblochpmx*nblochpmx/nwordr()
      open(newunit=ifrcwi,file='WVI.'//i2char(kx),action='read',form='unformatted',access='direct',recl=mrecl)
      do ix = 1,nx     ! imaginary frequency w'-loop
        ! ccccccccccccccccccccccccccccccc
        !          nrec=(kx-2)*niw+ix
        !          if(bzcase()==2) nrec= (kx-1)*niw+ix
        nrec=ix
        read(ifrcwi,rec=nrec) zw  ! Readin W-v on imag axis
        !          read(ifrcwi,rec=((kx-2)*niw+ix)) zw  ! Readin W-v on imag axis
        ! ccccccccccccccccccccccccccccccc

        ! zwz= S[i,j] <psi(q,t) |psi(q-rk,n) B(rk,i)>
        !                Wc(k,iw')(i,j) > <B(rk,j) psi(q-rk,n) |psi(q,t)>
        !        do itp = 1,ntp0
        !        do  it = 1,nstate
        !          zwz(ix,it,itp) = sum(
        !     &   dconjg(zmel(:,it,itp)),matmul(zw(1:ngb,1:ngb),zmel(:,it,itp)) )
        !        enddo
        !        enddo
        do ir3=1,nrws2
          do ir2=1,nrws2
            call matzwz3( zw(1:ngb,1:ngb), zmelc(:,:,:,ir2), zmelc(:,:,:,ir3), &
                 nwf,nwf,ngb, &
                 zw2)
            do ir1=1,nrws1
              ir = ir1 + (ir2-1 + (ir3-1)*nrws2)*nrws1
              rw_iw(:,:,:,:,ir,ix) = rw_iw(:,:,:,:,ir,ix) + dreal(zw2(:,:,:,:) * weightc(ir1))
              cw_iw(:,:,:,:,ir,ix) = cw_iw(:,:,:,:,ir,ix) + dimag(zw2(:,:,:,:) * weightc(ir1))
            enddo ! ir1
          enddo ! ir2
        enddo ! ir3
      enddo
      close(ifrcwi)! = iclose('WVI.'//i2char(kx))

      !====================================================================
      ! Wc(qt,w) along the real axis
      !====================================================================
      if(debug) write(6,*)' go to poles'
      nprecx=8
      mrecl  = nprecx*2*nblochpmx*nblochpmx/nwordr()
      open(newunit=ifrcw, file='WVR.'//i2char(kx),action='read',form='unformatted',access='direct',recl=mrecl)
      do      ix = 0,nrw                    ! real frequency w'-loop
        ! ccccccccccccccccccccccccc
        !          nrec=(kx-2)*(nw+1-nw_i)+ ix-nw_i+1
        !          if(bzcase()==2) nrec= (kx-1)*(nw+1-nw_i)+ ix-nw_i+1
        nrec=ix-nw_i+1
        read(ifrcw,rec=nrec) zw
        ! cccccccccccccccccccccccc
        ! zwz = S[i,j] <psi(q,t) |psi(q-rk,n) B(rk,i)> Wc(k,iw')(i,j) > <B(rk,j) psi(q-rk,n) |psi(q,t)>
        do ir3=1,nrws2
          do ir2=1,nrws2
            call matzwz3( zw(1:ngb,1:ngb), zmelc(:,:,:,ir2), zmelc(:,:,:,ir3), &
                 nwf,nwf,ngb, &
                 zw2)
            do ir1=1,nrws1
              ir = ir1 + (ir2-1 + (ir3-1)*nrws2)*nrws1
              rw_w(:,:,:,:,ir,ix)  = rw_w(:,:,:,:,ir,ix) + dreal(zw2(:,:,:,:) * weightc(ir1))
              cw_w(:,:,:,:,ir,ix)  = cw_w(:,:,:,:,ir,ix) + dimag(zw2(:,:,:,:) * weightc(ir1))
            enddo ! ir1
          enddo ! ir2
        enddo ! ir3
      enddo
      close(ifrcw)! = iclose('WVR.'//i2char(kx))
      deallocate(zmelc,zw,zw2)
      if(debug) write(6,*)' end of screen-if'
      ! end of if (exchange)
    endif
    continue  ! end of k-loop
1100 enddo kxloop
  return
 
contains
  !! psi2b_v2 and psicb_v2 are older versions of psi2b_v3 and psicb_v3 in m_zmel
  subroutine psi2b_v2(nt0,ntp0,iclass,coskt,sinkt, &
       cphik, cphikq,    ppb,  nlnmv,nlnmc,mdim,nctot,imdim,iatomp, &
       mdimx,nlmto,nbloch,nlnmx,noccxv,nt,ntp, &
       natom,natomx, &
       zpsi2b)
    ! originaly 92.03.17 by Ferdi.
    ! takao modified at Apr 2002
    ! calculates <psi(k',t') | psi(k,t) B(R,i)>
    ! for all R
    ! psi(k,t) = sum(RLn) b(RLn,k,t)*X(RLn,k)
    ! B(R,i)   = Bloch orthonormal product basis for atom R
    ! psi(k,t) is stored after nctot

    ! nt0        = no. t
    ! ntp0       = no. t'
    ! coskt,sinkt= exp(ik.T)
    ! cphik b(k)
    ! cphikq b(k')

    ! ppb        = <phi(RLn) phi(RL'n') B(R,i)>

    ! ddb        = <phidot(RLn) phidot(RL'n') B(R,i)>, s. ppbl.f
    ! nlnmv      = number of l,n,m for valence
    ! nlnmc      = number of n,l,m for core states
    ! mdim       = number of optimal product basis functions for each class
    ! nctot      = total no. allowed core states
    ! nbloch     = total no. optimal product basis
    ! nlnmx      = maximum number of l,n,m
    ! noccxv     = max. no. occupied valence states
    ! nt         = maximum number of occupied states
    ! ntp        = maximum number of unoccupied states

    ! zpsi2b     =  the matrix elements
    implicit real*8(a-h,o-z)
    implicit integer (i-n)
    integer::nt0,ntp0,natom,nlmto,natomx,nlnmc,mdim,nctot
    complex(8):: cphik(nlmto,noccxv),cphikq(nlmto,ntp0) &
         ,zpsi2b(nbloch,nt,ntp),phase
    dimension &
         ppb(nlnmx,nlnmx,mdimx,natom), &
         nlnmv(natom),nlnmc(natom),mdim(natom),iclass(natom), &
         coskt(natom),sinkt(natom),imdim(natom),iatomp(natom)

    integer,allocatable::iasx(:)
    integer :: nzppb1,nzppb2
    complex(8),allocatable :: zz(:,:), zppb(:,:)
    complex(8) :: alpha,beta

    ! zppb is used as work array for ppb(:,:,i,ic) and for zpsi2b(ib,:,:).
    nzppb1=max(nt0,nlnmx)
    nzppb2=max(ntp0,nlnmx)
    allocate( zz(nlnmx,ntp) )
    allocate( zppb(nzppb1,nzppb2) )
    beta=0d0  ; alpha=1d0

    ! check dimensions
    if (ntp0 > ntp) call rx( 'psi2bc: ntp exceeded')
    if (mdimx /= maxval(mdim)) call rx( 'psi2bc: wrong mdimx')
    if (nctot+nt0 > nt) call rx( 'psi2bc: nt exceeded')
    if (nt0 > noccxv) call rx( 'psi2bc: noccxv exceeded')
    if ( sum(mdim(iclass(1:natom)))/= nbloch ) call rx( 'psi2b_v2: wrong nbloch')
    allocate(iasx(natom))
    ias = 1
    do ia = 1,natom
       iasx(ia) = ias
       ias = ias + nlnmv(iclass(ia))
    enddo
    if(ias-1/=nlmto) call rx( ' psi2b_v2:sum(nlnmv)/= nlmto')
    ! loop over atoms
    do  ia = 1,natom
       ic   = iclass(ia)
       nc   = nlnmc(ic)
       nv   = nlnmv(ic)
       nc1  = nc + 1
       if (nc+ nlnmv(ic) > nlnmx) call rx( 'psi2b_v2: nlnmx exceeded')
       phase= dcmplx(coskt(ia),sinkt(ia))
       ias  = iasx(ia)
       iap  = iatomp(ia)
       icp  = iclass(iap)
       do   i = 1,mdim(icp) ! loop over optimal product basis
          zppb(1:nv,1:nv) = ppb(nc+1:nc+nv,nc+1:nc+nv,i,icp)
          call zgemm('T','N',nv,ntp0,nv, &
               alpha, zppb,nzppb1, cphikq(ias,1), nlmto,  beta, &
               zz,  nlnmx )
          do itp = 1,ntp0
             do jp = 1,nv
                zz(jp,itp)= dconjg(zz(jp,itp) )
             enddo
          enddo
          !----------------------------------------------------
          ! <psi(k+q,t') | psi(k,t) B(i)>
          call zgemm('T','N', nt0,ntp0,nv, &
               phase, cphik(ias,1),nlmto, zz,nlnmx, beta, &
               zppb, nzppb1 )
          ib = imdim(iap)-1+i
          zpsi2b(ib,nctot+1:nctot+nt0,1:ntp0)=zppb(1:nt0,1:ntp0)
          !------------------------------------------------------
       end do !end of optimal product basis-loop
    end do !end of atom-loop
    !      deallocate(rr,cc,iasx)
    deallocate(zz,zppb,iasx)
  end subroutine psi2b_v2
  subroutine psicb_v2 (icore,ncore,ntp0,iclass,coskt,sinkt, &
       cphikq, ppb,    nlnmv,nlnmc,mdim, &
       imdim,iatomp, &
       mdimx,nlmto,nbloch,nlnmx,nt,ntp,natom,natomx, &
       nl,nnc, &
       zpsi2b)! rpsi2b,cpsi2b)
    ! written by Ferdi  92.03.17
    ! takao modified at Apr 2002

    ! calculates <psi(k+q,t') | core(k,t) B(R,i)>
    ! for all R
    ! psi(k,t) = S[RLn] b(RLn,k,t)*X(RLn,k)
    !          = S[RLn] b(RLn,k,t)*Phi(RLn,k) + hb(RLn,k,t)*Phidot(RLn,k)
    ! core(k,t)= core states
    ! B(R,i)   = Bloch orthonormal product basis for atom R

    ! <psi(k+q,t') | core(k,t) B(R,i)>
    ! = S[RLn]  b(RLn,k+q,t')^* <Phi(RLn)    |core(k,t) B(R,i)>
    ! + S[RLn] hb(RLn,k+q,t')^* <Phidot(RLn) |core(k,t) B(R,i)>

    ! ntp0       = no. unoccupied states
    ! coskt,sinkt= exp(ik.T)
    ! cphikq  = real and imaginary part of b(k+q).
    !            coefficients of eigenfunctions for argumentationwaves in each MT

    ! icore      = index for core states
    ! ncore      = no. core states in each class
    ! ppb        = <Phi(RLn) Phi(RL'n') B(R,i)>

    ! nlnmv      = number of l,n,m for valence
    ! nlnmc      = number of l,n,m for core states
    ! mdim       = number of optimal product basis functions for each class
    ! nbloch     = total no. optimal product basis
    ! nlnmx      = maximum number of l,n,m
    ! nt         = maximum number of occupied states
    ! ntp        = maximum number of unoccupied states

    ! zpsi2b     =  the matrix elements

    implicit real*8(a-h,o-z)
    implicit integer (i-n)
    complex(8):: cphikq(nlmto,ntp0),zpsi2b(nbloch,nt,ntp),phase
    dimension &
         icore(nl*nl*nnc,natom),ncore(natom), &
         ppb(nlnmx,nlnmx,mdimx,natom), &
         nlnmv(natom),nlnmc(natom),mdim(natom),iclass(natom), &
         coskt(natom),sinkt(natom),imdim(natom),iatomp(natom)
    zpsi2b = 0d0
    ib         = 0
    ias        = 1
    ics        = 0
    do      ia = 1,natom
       ic         = iclass(ia)
       nc         = nlnmc(ic)
       nv         = nlnmv(ic)
       nc1        = nc + 1
       phase  =  dcmplx(coskt(ia),sinkt(ia))
       ! loop over optimal product basis
       iap        = iatomp(ia)
       icp        = iclass(iap)
       ib         = imdim(iap)-1
       do       i = 1,mdim(icp)
          ib         = ib + 1 
          ! S[Ln] bkq(Ln,t')^(*) * <Phi(Ln) core(L'n') B(i)> etc.
          ! for a given i, for all L'n' and t'
          ! bkq is complex but < > is real
          do     itp = 1,ntp0
             do      it = 1,ncore(ic)
                icr        = icore(it,ic)
                zpsi2b(ib,ics+it,itp) = phase* &
                     dconjg(sum(cphikq(ias:ias+nv-1,itp)*ppb(nc1:nc+nv,icr,i,icp)))

                ! end of t'(unoccupied)-loop
             end do
             ! end of t(occupied)-loop
          end do

          ! end of optimal product basis-loop
       end do

       ! end of atom-loop
       ias        = ias + nlnmv(ic)
       ics        = ics + ncore(ic)
    end do

    return
  end subroutine psicb_v2

  subroutine matzwz3(zw,zmel1,zmel2, ntp0,nstate,ngb, zwz)
    implicit none
    integer :: nstate,ntp0,itp,it,itp2,it2,ngb
    complex(8) :: zw(ngb,ngb),zmel1(ngb,nstate,ntp0), &
         zmel2(ngb,nstate,ntp0), &
         zwz(ntp0,nstate,nstate,ntp0)
    complex(8), allocatable :: CC(:,:,:)
    allocate(CC(ngb,nstate,ntp0) )
    call matm(zw,zmel2,cc, ngb, ngb, nstate*ntp0)
    do itp2 = 1,ntp0
       do  it2 = 1,nstate
          do  it  = 1,nstate
             do itp  = 1,ntp0
                zwz(itp,it,it2,itp2) &
                     = sum( dconjg(zmel1(1:ngb,it,itp))*CC(1:ngb,it2,itp2))
             enddo
          enddo
       enddo
    enddo
    deallocate(CC)
  end subroutine matzwz3

  subroutine readppovl0(q,ngc,ppovl)
    implicit none
    integer, intent(in) :: ngc
    complex(8), intent(out) :: ppovl(ngc,ngc)
    real(8), intent(in) :: q(3)
    integer:: ngc_r,ippovl0
    real(8):: qx(3),tolq=1d-8
    open(newunit=ippovl0,file='PPOVL0',form='unformatted')
    do
      read(ippovl0) qx,ngc_r
      if(sum(abs(qx-q))<tolq) then
        if(ngc_r/=ngc) call rx( 'readin ppovl: ngc_r/=ngc')
        read(ippovl0) ppovl
        exit
      endif
    enddo
    close(ippovl0)
  end subroutine readppovl0

  
end subroutine wmatqk_mpi

