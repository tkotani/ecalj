!>Write W-v. Gamma-cell averaged W-v
module m_llw
  use m_lgunit,only:stdo
  use m_ftox
  use m_rdpp,only: nblochpmx
  use m_genallcf_v3,only: natom,nspin,nl,nn, ndima,nlnmx, nctot,alat, deltaw, plat, pos, ecore, tpioa
  use m_freq,only: frhis,freq_r,freq_i, nwhis,nw_i,nw,npm,niw !output of getfreq
  use m_qbze,only: Setqbze, nqbze,nqibze,qbze,qibze
  use m_read_bzdata,only: Read_bzdata, ngrp2=>ngrp,nqbz,nqibz,n1,n2,n3,ginv,dq_,qbz,wbz,qibz,wibz, ntetf,idtetf,ib1bz
  use m_read_bzdata,only: qbzw,nqbzw, q0i,nq0i ,nq0iadd,ixyz
  use m_readVcoud,only: vcousq, ngb
  use m_rdpp,only: nbloch,mrecl
  use m_x0kf,only: zxq,zxqi
  use m_mpi, only: mpi__root_k, mpi__root_q, mpi__size_b,ipr
  use m_zmel, only: m2e_prod_basis
#ifdef __MP
  use m_mpi, only: MPI__GatherXqw => MPI__GatherXqw_kind4
#else
  use m_mpi, only: MPI__GatherXqw => MPI__GatherXqw
#endif
  use m_kind,only: kp => kindrcxq
  use m_stopwatch
  use m_blas, only: m_op_c, m_op_t
#ifdef __GPU
  use m_lapack, only: zminv => zminv_d
#else
  use m_lapack, only: zminv => zminv_h
#endif
#if defined(__MP) && defined(__GPU)
  use m_blas, only: gemm => cmm_d
#elif defined(__MP)
  use m_blas, only: gemm => cmm_h
#elif defined(__GPU)
  use m_blas, only: gemm => zmm_d
#else
  use m_blas, only: gemm => zmm_h
#endif
  implicit none
  public:: WVRllwR,WVIllwI,  MPI__sendllw,MPI__sendllw2
  complex(8),allocatable,protected,public:: llw(:,:), llwI(:,:)
  complex(8),allocatable,protected,public:: wmuk(:,:)
  logical,protected,public:: w4pmode
  integer,protected,public:: ngbq0
  private
  complex(8),allocatable:: zw0(:,:)
  complex(kind=kp), allocatable :: zw(:,:)
  complex(8),allocatable:: epstinv(:,:),epstilde(:,:)
  real(8),parameter:: pi=4d0*datan(1d0),fourpi = 4d0*pi
  logical :: tzminv=.true.
contains
  subroutine WVRllwR(q,iq,nmbas1,nmbas2,is_m_basis)
    use m_readqg,only: Readqg0
    intent(in)::       q,iq,    nmbas1,nmbas2 !zxq can be twiced when nspin=2
    logical, intent(in) :: is_m_basis
    integer:: iq,iq0,nwmax,nwmin,iw,imode,ix,igb1,igb2,ifllw
    integer:: nmbas1,nmbas2,ngc0,ifw4p,ifrcw,mreclx
    real(8):: frr,q(3),vcou1,quu(3),eee
    logical::  localfieldcorrectionllw,cmdopt0,emptyrun
    complex(kind=kp), allocatable :: zxqw(:,:), x_m2e(:,:)
    logical,save:: init=.true.
    type(stopwatch) :: t_sw_matinv, t_sw_x_gather, t_sw_x_m2e_xf
    integer :: istat
!    complex(8):: zxq(nmbas1,nmbas2,nw_i:nw)
    character(10):: i2char
    mreclx=mrecl
    emptyrun=.false. !cmdopt0('--emptyrun')
    if(init) then !initialization related to w4pmode, zw, tpioa...
       allocate( llw(nw_i:nw,nq0i),source=(1d99,0d0) )
       if(sum(ixyz)/=0) w4pmode= .TRUE. 
       if(w4pmode) allocate( wmuk(2:nblochpmx,3),source=(1d99,0d0))
       allocate( zw(nblochpmx,nblochpmx),source=(0d0,0d0) )
       init=.false.
    endif
    call stopwatch_init(t_sw_matinv, 'matinv')
    call stopwatch_init(t_sw_x_gather, 'gather')
    call stopwatch_init(t_sw_x_m2e_xf, 'xf chi: M2E')
    call readqg0('QGcou', (/0d0,0d0,0d0/),  quu,ngc0) ! ngb is q-dependent. released at the end of WVIllwi
    ngbq0 = nbloch+ngc0
    allocate( zw0(ngb,ngb), epstilde(ngb,ngb), epstinv(ngb,ngb))
    allocate (zxqw(ngb, ngb), x_m2e(ngb,ngb))
    !$acc enter data create(zxqw, epstilde, epstinv, zw, zw0, x_m2e) copyin(vcousq)
    if(nspin == 1) then
      !$acc kernels present(zxq)
      zxq(:,:,:) = 2d0*zxq(:,:,:)
      !$acc end kernels
    endif
    nwmax = nw
    nwmin = nw_i
    if(ipr)write(stdo,ftox)" === trace check for W-V === nqibz nwmin nwmax=",nqibz,nwmin,nwmax, 'iq q=',iq,ftof(q)
    if(ipr)write(stdo,ftox) 'size of zxq:',size(zxq,1), size(zxq,2), size(zxq,3)
    call flush(stdo)
    if(iq<=nqibz) then        !for mmmw
      if(mpi__root_q) then
        open(newunit=ifrcw, file='__WVR.'//i2char(iq),form='unformatted',access='direct',recl=mreclx)
      endif
      iwloop: do 1015 iw  = nwmin,nwmax
        if(emptyrun) exit
        frr= dsign(freq_r(abs(iw)),dble(iw))
        if(iq==1) then
          ix=1
          !$acc kernels
          zw0(:,1)=0d0
          zw0(1,:)=0d0
          !$acc end kernels
        else
          ix=0
        endif
        call stopwatch_start(t_sw_x_gather)
        if(mpi__size_b == 1) then
          !$acc kernels
          zxqw(:,:) = zxq(:,:,iw)
          !$acc end kernels
        else
          !$acc update host(zxq(1:nmbas1,1:nmbas2,iw))
          call MPI__GatherXqw(zxq(:,:,iw), zxqw, nmbas1, nmbas2)
          !$acc update device(zxqw)
        endif
        call stopwatch_pause(t_sw_x_gather)

        MToEBasisTransformation1: if(is_m_basis) then
          call stopwatch_start(t_sw_x_m2e_xf)
          !$acc host_data use_device(zxqw, m2e_prod_basis, x_m2e)
          istat = gemm(zxqw, m2e_prod_basis, x_m2e, ngb, ngb, ngb)
          istat = gemm(m2e_prod_basis, x_m2e, zxqw, ngb, ngb, ngb, opA=m_op_C)
          !$acc end host_data
          call stopwatch_pause(t_sw_x_m2e_xf)
        endif MToEBasisTransformation1

        rootq_if1:if(mpi__root_q) then
        !$acc kernels loop independent collapse(2)
        do igb1=ix+1,ngb !  Eqs.(37),(38) in PRB81 125102 (Friedlich)
          do igb2=ix+1,ngb
            epstilde(igb1,igb2)= -vcousq(igb1)*zxqw(igb1,igb2)*vcousq(igb2)
            if(igb1==igb2) epstilde(igb1,igb2)=1+epstilde(igb1,igb2)
          enddo
        enddo
        !$acc end kernels
        !$acc kernels 
        epstinv(ix+1:ngb,ix+1:ngb)=epstilde(ix+1:ngb,ix+1:ngb)
        !$acc end kernels
        call stopwatch_start(t_sw_matinv)
        if(tzminv) then
          !$acc data copy(epstinv)
          !$acc host_data use_device(epstinv)
          istat = zminv(epstinv(ix+1,ix+1), n=ngb-ix, lda=ngb)
          !$acc end host_data
          !$acc end data
        else
          !$acc update host(epstinv)
          call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
        endif
        call stopwatch_pause(t_sw_matinv)
        !  w4p writing eps
        if(iw==0 .AND. w4pmode) then ! static epstinv is saved. For q=0 epstilde (mu=1 skipped). For q/=0 full matrix inversion. ix=1 is set for q=0)
          !$acc update host(epstinv)
          open(newunit=ifw4p,file='__W4PHONON.'//i2char(iq),form='unformatted')
          write(ifw4p) iq,q,ngb,ix !ix=0, or ix=1 for q=0 (iq=1)
          write(ifw4p) epstinv(ix+1:ngb,ix+1:ngb)
          close(ifw4p)
        endif
        !$acc kernels loop independent collapse(2)
        do igb1=1+ix,ngb
          do igb2=1+ix,ngb
            zw0(igb1,igb2)= vcousq(igb1)*epstinv(igb1,igb2)*vcousq(igb2)
            if(igb1==igb2) zw0(igb1,igb2)= zw0(igb1,igb2)-vcousq(igb1)*vcousq(igb2)
          enddo
        enddo
        !$acc end kernels
        !$acc kernels
        zw(:,:) = (0_kp,0_kp)
        zw(1:ngb,1:ngb) = cmplx(zw0(1:ngb,1:ngb),kind=kp)
        !$acc end kernels
        !$acc update host(zw)
1012    continue
        write(ifrcw, rec= iw-nw_i+1 ) zw !  WP = vsc-v
        call tr_chkwrite("freq_r iq iw realomg trwv=", zw, iw, frr,nblochpmx, nbloch,ngb,iq)
        endif rootq_if1
1015  enddo iwloop
      if(mpi__root_q) close(ifrcw)
    else  ! llw, Wing elements of W. See PRB81 125102
      iq0 = iq - nqibz
      vcou1 = fourpi/sum(q**2*tpioa**2) ! --> vcousq(1)**2!  !fourpi/sum(q**2*tpioa**2-eee)
      do 1115 iw  = nwmin,nwmax
        if(emptyrun) exit
        frr= dsign(freq_r(abs(iw)),dble(iw))
        !! Full inversion to calculalte eps with LFC.
        !if(localfieldcorrectionllw()) then
        call stopwatch_start(t_sw_x_gather)
        if(mpi__size_b == 1) then
          !$acc kernels
          zxqw(:,:) = zxq(:,:,iw)
          !$acc end kernels
        else
          !$acc update host(zxq(1:nmbas1,1:nmbas2,iw))
          call MPI__GatherXqw(zxq(:,:,iw), zxqw, nmbas1, nmbas2)
          !$acc update device(zxqw)
        endif
        !for log output 
        !$acc update host(zxqw(1,1))
        call stopwatch_pause(t_sw_x_gather)
        MToEBasisTransformation2: if(is_m_basis) then
          call stopwatch_start(t_sw_x_m2e_xf)
          !$acc host_data use_device(zxqw, m2e_prod_basis, x_m2e)
          istat = gemm(zxqw, m2e_prod_basis, x_m2e, ngb, ngb, ngb)
          istat = gemm(m2e_prod_basis, x_m2e, zxqw, ngb, ngb, ngb, opA=m_op_C)
          !$acc end host_data
          call stopwatch_pause(t_sw_x_m2e_xf)
        endif MToEBasisTransformation2
        rootq_if2:if(mpi__root_q) then
        ix=0
        eee=0d0
        !$acc kernels loop independent collapse(2)
        do igb1=ix+1,ngb
          do igb2=ix+1,ngb
            if(igb1==1 .AND. igb2==1) then
              epstilde(igb1,igb2)= 1d0 - vcou1*zxqw(1,1)
              cycle
            endif
            epstilde(igb1,igb2)= -vcousq(igb1)*zxqw(igb1,igb2)*vcousq(igb2)
            if(igb1==igb2) then
              epstilde(igb1,igb2)=1d0 + epstilde(igb1,igb2)
            endif
          enddo
        enddo
        !$acc end kernels
        !$acc kernels 
        epstinv(ix+1:ngb,ix+1:ngb)=epstilde(ix+1:ngb,ix+1:ngb)
        !$acc end kernels
        call stopwatch_start(t_sw_matinv)
        if(tzminv) then
          !$acc data copy(epstinv)
          !$acc host_data use_device(epstinv)
          istat = zminv(epstinv(ix+1,ix+1), n=ngb-ix, lda=ngb)
          !$acc end host_data
          !$acc end data
        else
          !$acc update host(epstinv)
          call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
          !$acc update device(epstinv)
        endif
        !$acc update host(epstinv)
        call stopwatch_pause(t_sw_matinv)
        if(iq0<=nq0i) llw(iw,iq0)= 1d0/epstinv(1,1)
        !     ! Wing elements calculation july2016    ! We need check nqb is the same as that of q=0
        if(ixyz(iq0)/=0 .AND. iw==0) then
          if(ngb/=ngbq0) then
            if(ipr)write(6,*)q,iq0,ngb,ngbq0
            call rx('hx0p0_sc: ngb/=ngbq0')
          endif
          wmuk(2:ngb,ixyz(iq0))=epstinv(1,2:ngb)/epstinv(1,1) ! this is dot(q(:)*w_mu(:,igb)). See PRB125102(2016) eq.(36)
        endif
        !else
        !   if(iq0<=nq0i) llw(iw,iq0)= 1d0 - vcou1*zxq(1,1,iw)
        !endif
        if(iq0<=nq0i.and.ipr) write(6,"('epsWVR: iq iw_R omg(iw) eps(wFC) eps(woLFC) ', &
             2i5,x,10(d13.6,2x,d13.6,x,d13.6,2x,d13.6,x,d13.6))") &
             iq,iw,freq_r(iw),llw(iw,iq0),1d0-vcou1*zxqw(1,1)
             ! iq,iw,freq_r(iw),llw(iw,iq0),1d0-vcou1*zxq(1,1,iw)
        continue               !iw
        endif rootq_if2
1115  enddo
    endif
    !$acc exit data delete(zxqw, epstilde, epstinv, zw, zw0, vcousq, x_m2e)
    deallocate(zxqw)
    if(mpi__root_q) then 
      call stopwatch_show(t_sw_x_gather)
      call stopwatch_show(t_sw_matinv)
      if(is_m_basis) call stopwatch_show(t_sw_x_m2e_xf)
    endif
  end subroutine WVRllwR
  subroutine WVIllwi(q,iq,nmbas1,nmbas2,is_m_basis)
    intent(in)::       q,iq,     nmbas1,nmbas2 !zxqi can be twiced when nspin=2
    integer:: nmbas1,nmbas2,mreclx
    integer:: iq,iq0,nwmax,nwmin,iw,imode,ix,igb1,igb2,ifllwi,ifrcwi
    real(8):: frr,q(3),vcou1
    logical::  localfieldcorrectionllw,cmdopt0,emptyrun
    logical, intent(in) :: is_m_basis
    logical,save:: init=.true.
    complex(kind=kp), allocatable :: zxqw(:,:), x_m2e(:,:)
!    complex(8):: zxqi(nmbas1,nmbas2,niw)
    character(10):: i2char
    integer :: istat
    type(stopwatch) :: t_sw_matinv, t_sw_x_gather, t_sw_x_m2e_xf
    allocate (zxqw(ngb, ngb), x_m2e(ngb,ngb))
    !$acc enter data create(zxqw, x_m2e, epstilde, epstinv, zw, zw0) copyin(vcousq)
    mreclx=mrecl
    emptyrun=.false. !cmdopt0('--emptyrun')
    if(init) then
       allocate( llwI(niw,nq0i) )
       init=.false.
       llwI= 1d99
    endif
    call stopwatch_init(t_sw_matinv, 'matinv')
    call stopwatch_init(t_sw_x_gather, 'gather')
    call stopwatch_init(t_sw_x_m2e_xf, 'xf chi: M2E')
    if(ipr)write(6,*)'WVRllwI: init'
    if (nspin == 1) then
      !$acc kernels present(zxqi)
      zxqi(:,:,:) = 2d0*zxqi(:,:,:) ! if paramagnetic, multiply x0 by 2
      !$acc end kernels
    endif
    if( iq<=nqibz ) then
       if(mpi__root_q) then
         open(newunit=ifrcwi,file='__WVI.'//i2char(iq),form='unformatted',access='direct',recl=mreclx)
       endif
       do 1016 iw  = 1,niw
          if(emptyrun) exit
          !!  Eqs.(37),(38) in PRB81 125102
          if(iq==1) then
             ix=1
             !$acc kernels
             zw0(:,1)=0d0
             zw0(1,:)=0d0
             !$acc end kernels
          else
             ix=0
          endif
          call stopwatch_start(t_sw_x_gather)
          if(mpi__size_b == 1) then
            !$acc kernels
            zxqw(:,:) = zxqi(:,:,iw)
            !$acc end kernels
          else
            !$acc update host(zxqi(1:nmbas1,1:nmbas2,iw))
            call MPI__GatherXqw(zxqi(:,:,iw), zxqw, nmbas1, nmbas2)
            !$acc update device(zxqw)
          endif
          call stopwatch_pause(t_sw_x_gather)
          MToEBasisTransformation1: if(is_m_basis) then
            call stopwatch_start(t_sw_x_m2e_xf)
            !$acc host_data use_device(zxqw, m2e_prod_basis, x_m2e)
            istat = gemm(zxqw, m2e_prod_basis, x_m2e, ngb, ngb, ngb)
            istat = gemm(m2e_prod_basis, x_m2e, zxqw, ngb, ngb, ngb, opA=m_op_C)
            !$acc end host_data
            call stopwatch_pause(t_sw_x_m2e_xf)
          endif MToEBasisTransformation1
          if(mpi__root_q) then
          !$acc kernels loop independent collapse(2)
          do igb1=ix+1,ngb
             do igb2=ix+1,ngb
                epstilde(igb1,igb2)= -vcousq(igb1)*zxqw(igb1,igb2)*vcousq(igb2)
                if(igb1==igb2) epstilde(igb1,igb2)=1+epstilde(igb1,igb2)
             enddo
          enddo
          !$acc end kernels
          !$acc kernels 
          epstinv(:,:) = epstilde(:,:)
          !$acc end kernels
          call stopwatch_start(t_sw_matinv)
          if(tzminv) then
            !$acc data copy(epstinv)
            !$acc host_data use_device(epstinv)
            istat = zminv(epstinv(ix+1,ix+1), n=ngb-ix, lda=ngb)
            !$acc end host_data
            !$acc end data
          else
            !$acc update host(epstinv)
            call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
            !$acc update device(epstinv)
          endif
          call stopwatch_pause(t_sw_matinv)
          !$acc kernels loop independent collapse(2)
          do igb1=ix+1,ngb
             do igb2=ix+1,ngb
                zw0(igb1,igb2)= vcousq(igb1)*epstinv(igb1,igb2)*vcousq(igb2)
                if(igb1==igb2) zw0(igb1,igb2)= zw0(igb1,igb2)-vcousq(igb1)*vcousq(igb2)
             enddo
          enddo
          !$acc end kernels
1014      continue
          !$acc kernels
          zw(:,:) = (0_kp,0_kp)
          zw(1:ngb,1:ngb) = cmplx(zw0(1:ngb,1:ngb),kind=kp)
          !$acc end kernels
          !$acc update host(zw)
          write(ifrcwi, rec= iw)  zw !  WP = vsc-v
          call tr_chkwrite("freq_i iq iw imgomg trwv=",zw,iw,freq_i(iw),nblochpmx,nbloch,ngb,iq)
          endif
1016   enddo
       if(mpi__root_q) close(ifrcwi)
    else
       !! Full inversion to calculalte eps with LFC.
       iq0 = iq - nqibz
       vcou1 = fourpi/sum(q**2*tpioa**2) ! --> vcousq(1)**2!  !fourpi/sum(q**2*tpioa**2-eee)
       do 1116 iw  = 1,niw
          if(emptyrun) exit
          !if(localfieldcorrectionllw()) then
          call stopwatch_start(t_sw_x_gather)
          if(mpi__size_b == 1) then
            !$acc kernels
            zxqw(:,:) = zxqi(:,:,iw)
            !$acc end kernels
          else
            !$acc update host(zxqi(1:nmbas1,1:nmbas2,iw))
            call MPI__GatherXqw(zxqi(:,:,iw), zxqw, nmbas1, nmbas2)
            !$acc update device(zxqw)
          endif
          !$acc update host(zxqw(1,1))
          call stopwatch_pause(t_sw_x_gather)
          MToEBasisTransformation2: if(is_m_basis) then
            call stopwatch_start(t_sw_x_m2e_xf)
            !$acc host_data use_device(zxqw, m2e_prod_basis, x_m2e)
            istat = gemm(zxqw, m2e_prod_basis, x_m2e, ngb, ngb, ngb)
            istat = gemm(m2e_prod_basis, x_m2e, zxqw, ngb, ngb, ngb, opA=m_op_C)
            !$acc end host_data
          call stopwatch_pause(t_sw_x_m2e_xf)
          endif MToEBasisTransformation2
          if(mpi__root_q) then
             ix=0
             !$acc kernels loop independent collapse(2)
             do igb1=ix+1,ngb
                do igb2=ix+1,ngb
                   if(igb1==1 .AND. igb2==1) then
                      ! epstilde(igb1,igb2)= 1d0 - vcou1*zxqi(1,1,iw)
                      epstilde(igb1,igb2)= 1d0 - vcou1*zxqw(1,1)
                      cycle
                   endif
                   ! epstilde(igb1,igb2)= -vcousq(igb1)*zxqi(igb1,igb2,iw)*vcousq(igb2)
                   epstilde(igb1,igb2)= -vcousq(igb1)*zxqw(igb1,igb2)*vcousq(igb2)
                   if(igb1==igb2) then
                      epstilde(igb1,igb2)=1d0 + epstilde(igb1,igb2)
                   endif
                enddo
             enddo
             !$acc end kernels
             !$acc kernels 
             epstinv(ix+1:ngb,ix+1:ngb)=epstilde(ix+1:ngb,ix+1:ngb)
             !$acc end kernels
             call stopwatch_start(t_sw_matinv)
             if(tzminv) then
               !$acc data copy(epstinv)
               !$acc host_data use_device(epstinv)
               istat = zminv(epstinv(ix+1,ix+1), n=ngb-ix, lda=ngb)
               !$acc end host_data
               !$acc end data
             else
               !$acc update host(epstinv)
               call matcinv(ngb-ix,epstinv(ix+1:ngb,ix+1:ngb))
               !$acc update device(epstinv)
             endif
             !$acc update host(epstinv)
             call stopwatch_pause(t_sw_matinv)
             if(iq0<=nq0i) llwI(iw,iq0)= 1d0/epstinv(1,1)
          !else
          !   if(iq0<=nq0i) llwI(iw,iq0)=  1d0 -vcou1*zxqi(1,1,iw)
          !endif
          ! if(iq0<=nq0i) write(6,"('iq iw_img eps(wLFC) eps(noLFC)',i4,i4,2f10.4,2x,2f10.4)") &
          !      iq,iw,llwI(iw,iq0),1d0-vcou1*zxqi(1,1,iw)
             if(iq0<=nq0i) write(6,"('iq iw_img eps(wLFC) eps(noLFC)',i4,i4,2f10.4,2x,2f10.4)") &
                  iq,iw,llwI(iw,iq0),1d0-vcou1*zxqw(1,1)
          endif
1116   enddo
    endif
    !$acc exit data delete(zxqw, epstilde, epstinv, zw, zw0, vcousq, x_m2e) 
    deallocate(epstinv,epstilde,zw0)
    deallocate(zxqw)
    if(mpi__root_q) then 
      call stopwatch_show(t_sw_x_gather)
      call stopwatch_show(t_sw_matinv)
      if(is_m_basis) call stopwatch_show(t_sw_x_m2e_xf)
    endif
  end subroutine WVILLWI
  subroutine MPI__sendllw2(iqxend,MPI__ranktab) !for hx0fp0
    use m_mpi,only: MPI__root,MPI__DbleCOMPLEXsend,MPI__DbleCOMPLEXrecv,MPI__rank,MPI__size
    intent(in)::             iqxend
    integer:: iq0,dest,src,iq,iqxend,MPI__ranktab(:)
    !! === Recieve llw and llwI at node 0, where q=0(iq=1) is calculated. ===
    if(MPI__size==1) return
    do iq=nqibz+1,iqxend
      iq0 = iq - nqibz
      if(MPI__ranktab(iq)==0) cycle 
      if(MPI__ranktab(iq) == MPI__rank) then
        dest=0
        call MPI__DbleCOMPLEXsend(llw(nw_i,iq0),(nw-nw_i+1),dest)
        call MPI__DbleCOMPLEXsend(llwI(1,iq0),niw,dest)
      elseif(MPI__root) then
        src=MPI__ranktab(iq)
        call MPI__DbleCOMPLEXrecv(llw(nw_i,iq0),(nw-nw_i+1),src)
        call MPI__DbleCOMPLEXrecv(llwI(1,iq0),niw,src)
      endif
    enddo
  end subroutine MPI__sendllw2
  subroutine MPI__sendllw(iqxend,MPI__Qranktab) !for hx0fp0_sc
    use m_mpi,only: MPI__DbleCOMPLEXsendQ,MPI__DbleCOMPLEXrecvQ,MPI__size,MPI__rank,MPI__root
    ! === Recieve llw and llwI at node 0, where q=0(iq=1) is calculated. ===
    intent(in)::            iqxend
    integer:: iq0,dest,src,iq,iqxend,MPI__Qranktab(:)
    if(MPI__size==1) return
    do iq=nqibz+1,iqxend
      iq0 = iq - nqibz
      if(MPI__Qranktab(iq)==0) cycle
      if(MPI__Qranktab(iq) == MPI__rank) then
        dest=0
        if(iq0<=nq0i) then
          call MPI__DbleCOMPLEXsendQ(llw(nw_i,iq0),(nw-nw_i+1),dest)
          call MPI__DbleCOMPLEXsendQ(llwI(1,iq0),niw,dest)
        endif
        if(ixyz(iq0)/=0) call MPI__DbleCOMPLEXsendQ(wmuk(2:ngbq0,ixyz(iq0)),ngbq0-1,dest)
      elseif(MPI__root) then
        src=MPI__Qranktab(iq)
        if(iq0<=nq0i) then
          call MPI__DbleCOMPLEXrecvQ(llw(nw_i,iq0),(nw-nw_i+1),src)
          call MPI__DbleCOMPLEXrecvQ(llwI(1,iq0),niw,src)
        endif
        if(ixyz(iq0)/=0) call MPI__DbleCOMPLEXrecvQ(wmuk(2:ngbq0,ixyz(iq0)),ngbq0-1,src)
      endif
    enddo
  end subroutine MPI__sendllw
end module m_llw
!===================================================================
subroutine tr_chkwrite(tagname,zw,iw,freqq,nblochpmx,nbloch,ngb,iq)
  use m_kind,only: kp => kindrcxq
  use m_mpi,only:ipr
  implicit none
  integer:: nblochpmx,nbloch,ngb,iw,i,iq
  complex(kind=kp):: zw(nblochpmx,nblochpmx)
  complex(8):: trwv,trwv2
  real(8):: freqq
  character*(*)::tagname
  trwv=0d0
  do i = 1,nbloch
     trwv = trwv + zw(i,i)
  enddo
  trwv2 = 0d0
  do i = 1,ngb
     trwv2 = trwv2 + zw(i,i)
  enddo  !  write(6,'(" realomg trwv=",2i6,4d22.14)') iq,iw,trwv(iw),trwv2(iw)
  if(ipr) write(6,'(a,f10.4,2i5,4d22.14)')tagname,freqq,iq,iw,trwv,trwv2
  if(ipr) call flush(6)
end subroutine tr_chkwrite
