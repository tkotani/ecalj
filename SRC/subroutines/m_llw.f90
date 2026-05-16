!>Write W-v. Gamma-cell averaged W-v
module m_llw
  use m_lgunit,only:stdo
  use m_ftox
  use m_rdpp,only: nblochpmx
  use m_struct_from_lmf,only: natom,nspin,nl,alat,plat,pos,tpioa; use m_gw_product_basis,only: nn,ndima,nlnmx; use m_core_state,only: nctot,ecore; use m_gw_user_config,only: deltaw
  use m_freq,only: frhis,freq_r,freq_i, nwhis,nw_i,nw,npm,niw !output of getfreq
  use m_qbze,only: Setqbze, nqbze,nqibze,qbze,qibze
  use m_read_bzdata,only: Read_bzdata, ngrp2=>ngrp,nqbz,nqibz,n1,n2,n3,ginv,dq_,qbz,wbz,qibz,wibz, ntetf,idtetf,ib1bz
  use m_read_bzdata,only: qbzw,nqbzw, q0i,nq0i ,nq0iadd,ixyz
  use m_readVcoud,only: vcousq, ngb
  use m_rdpp,only: nbloch,mrecl
  use m_x0kf,only: zxq,zxqi
  use m_mpi, only: mpi__root_k => mpi__root_k_xq, mpi__root_q, &
                   mpi__size_b => mpi__size_b_xq, ipr, comm_root_k => comm_root_k_xq, &
                   mpi__rank_b => mpi__rank_b_xq, mpi__rank_root_k => mpi__rank_root_k_xq, &
                   MPI__AllreduceSum
  use mpi
  use m_zmel, only: m2e_prod_basis
#ifdef __MP
  use m_mpi, only: MPI__GatherXqw => MPI__GatherXqw_c
#else
  use m_mpi, only: MPI__GatherXqw => MPI__GatherXqw
#endif
  ! Step WA1/WB.3a/WB.3e: W data goes through the m_wv_storage singleton
  ! via wv_put_real / wv_put_imag (SHM backend: shm_wvr/shm_wvi windows).
  use m_wv_storage, only: &
       wv_open_iq_real_for_write, wv_open_iq_imag_for_write, &
       wv_put_real, wv_put_imag, wv_close_iq_for_write
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
  public:: WVRllwR,WVIllwI,  MPI__sendllw,MPI__sendllw2,MPI__sendllw_q, &
           MPI__irecvllw_q, MPI__isendllw_q, MPI__waitllw
  complex(8),allocatable,protected,public:: llw(:,:), llwI(:,:)
  complex(8),allocatable,protected,public:: wmuk(:,:)
  logical,protected,public:: w4pmode
  integer,protected,public:: ngbq0
  private
  real(8),parameter:: pi=4d0*datan(1d0),fourpi = 4d0*pi
  ! Non-blocking llw transfer: request list shared across irecv/isend/wait calls.
  integer, save :: llw_nreqs = 0
  integer, save :: llw_reqs(60)   ! generous upper bound (3 msgs × ~20 aux q-pts)
contains
  subroutine WVRllwR(q,iq,nmbas1,nmbas2,is_x0_m_basis,is_wc_m_basis)
    use m_readqg,only: Readqg0
    intent(in)::       q,iq,    nmbas1,nmbas2 !zxq can be twiced when nspin=2
    logical, intent(in) :: is_x0_m_basis, is_wc_m_basis
    integer:: iq,iq0,nwmax,nwmin,iw,imode,ix,igb1,igb2,ifllw
    integer:: nmbas1,nmbas2,ngc0,ifw4p
    real(8):: frr,q(3),vcou1,quu(3),eee
    logical::  localfieldcorrectionllw,cmdopt0
    logical,save:: init=.true.
    type(stopwatch) :: t_sw_matinv, t_sw_x_gather, t_sw_x_m2e_xf
    integer :: istat
    character(10):: i2char
    complex(kind=kp), allocatable :: zw(:,:), zxqw(:,:), x_m2e(:,:)
    complex(8), allocatable, target :: epstilde(:,:)
    complex(8), pointer :: epstinv(:,:) => null()
    integer :: iwblock, jw, irank, ierr
#ifdef __GPU
    attributes(device) :: epstinv, epstilde
#endif
    if(init) then !initialization related to w4pmode, zw, tpioa...
       allocate( llw(nw_i:nw,nq0i),source=(0d0,0d0) )
       if(sum(ixyz)/=0) w4pmode= .TRUE. 
       if(w4pmode) allocate( wmuk(2:nblochpmx,3),source=(0d0,0d0))
       init=.false.
    endif
    call stopwatch_init(t_sw_matinv, 'matinv')
    call stopwatch_init(t_sw_x_gather, 'gather')
    call stopwatch_init(t_sw_x_m2e_xf, 'xf chi/W: M2E/E2M')
    call readqg0('QGcou', (/0d0,0d0,0d0/),  quu,ngc0) ! ngb is q-dependent. released at the end of WVIllwi
    ngbq0 = nbloch+ngc0
    allocate(epstilde(ngb,ngb))
    allocate(zxqw(ngb, ngb))
    allocate(x_m2e(ngb,ngb))
    allocate(zw(nblochpmx,nblochpmx))

    !$acc enter data create(zxqw, zw, x_m2e) copyin(vcousq)
    if(nspin == 1) then
      ! SHM: only rank 0 scales shared memory; barrier ensures others see the result.
      if (mpi__rank_root_k == 0) zxq(:,:,:) = 2d0*zxq(:,:,:)
      call MPI_barrier(comm_root_k, ierr)
    endif
    nwmax = nw
    nwmin = nw_i
    if(ipr)write(stdo,ftox)" === trace check for W-V === nqibz nwmin nwmax=",nqibz,nwmin,nwmax, 'iq q=',iq,ftof(q)
    if(ipr)write(stdo,ftox) 'size of zxq:',size(zxq,1), size(zxq,2), size(zxq,3)
    call flush(stdo)
    if(iq<=nqibz) then        !for mmmw
      call wv_open_iq_real_for_write(iq, comm=comm_root_k)
      ix = merge(1, 0, iq == 1)
      iwloop: do 1015 iwblock = nwmin, nwmax, mpi__size_b
         iw = iwblock + mpi__rank_b
         if(iw > nwmax) cycle
        !$acc kernels
          zw(:,:) = (0_kp, 0_kp)
        !$acc end kernels
        call stopwatch_start(t_sw_x_gather)
        zxqw(:,:) = zxq(:,:,iw)    ! all root_k have full shm_wvr; direct copy
        !$acc update device(zxqw)
        call stopwatch_pause(t_sw_x_gather)
        MToEBasisTransformation1: if(is_x0_m_basis) then
          call stopwatch_start(t_sw_x_m2e_xf)
          !$acc host_data use_device(zxqw, m2e_prod_basis, x_m2e)
          istat = gemm(zxqw, m2e_prod_basis, x_m2e, ngb, ngb, ngb)
          istat = gemm(m2e_prod_basis, x_m2e, zxqw, ngb, ngb, ngb, opA=m_op_C)
          !$acc end host_data
          call stopwatch_pause(t_sw_x_m2e_xf)
        endif MToEBasisTransformation1
        !$acc kernels loop independent collapse(2)
        do igb2=ix+1,ngb
          do igb1=ix+1,ngb !  Eqs.(37),(38) in PRB81 125102 (Friedlich)
            epstilde(igb1,igb2)= -vcousq(igb1)*zxqw(igb1,igb2)*vcousq(igb2)
            if(igb1==igb2) epstilde(igb1,igb2)=1+epstilde(igb1,igb2)
          enddo
        enddo
        !$acc end kernels
        call stopwatch_start(t_sw_matinv)
        istat = zminv(epstilde(ix+1,ix+1), n=ngb-ix, lda=ngb)
        epstinv => epstilde
        call stopwatch_pause(t_sw_matinv)
        !  w4p writing eps
        if(iw==0 .AND. w4pmode) then ! static epstinv is saved. For q=0 epstilde (mu=1 skipped). For q/=0 full matrix inversion. ix=1 is set for q=0)
          block
          real(kind=kp) :: epstinv_h(ngb,ngb)
          epstinv_h(:,:) = epstinv(:,:) !copy from CPU to GPU
          open(newunit=ifw4p,file='__W4PHONON.'//i2char(iq),form='unformatted')
          write(ifw4p) iq,q,ngb,ix !ix=0, or ix=1 for q=0 (iq=1)
          write(ifw4p) epstinv_h(ix+1:ngb,ix+1:ngb)
          close(ifw4p)
          endblock
        endif
        !$acc kernels loop independent collapse(2)
        do igb2=1+ix,ngb
          do igb1=1+ix,ngb
            zw(igb1,igb2)= vcousq(igb1)*epstinv(igb1,igb2)*vcousq(igb2)
            if(igb1==igb2) zw(igb1,igb2)= zw(igb1,igb2)-vcousq(igb1)*vcousq(igb2)
          enddo
        enddo
        !$acc end kernels
        nullify(epstinv)
        ! W^[M] = AW^[E]A^\dagger
        EtoMBasisTransformation: if(is_wc_m_basis .and. iq /= 1) then
          call stopwatch_start(t_sw_x_m2e_xf)
          !$acc host_data use_device(zw, m2e_prod_basis, x_m2e)
          istat = gemm(zw, m2e_prod_basis, x_m2e, ngb, ngb, ngb, opB=m_op_C, ldA=nblochpmx)
          istat = gemm(m2e_prod_basis, x_m2e, zw, ngb, ngb, ngb, ldC=nblochpmx)
          !$acc end host_data
          call stopwatch_pause(t_sw_x_m2e_xf)
        endif EtoMBasisTransformation
        !$acc update host(zw)
        ! write(ifrcw, rec= iw-nw_i+1 ) zw !  WP = vsc-v
        call wv_put_real(iw, zw(1:nblochpmx, 1:nblochpmx))
        frr= dsign(freq_r(abs(iw)),dble(iw))
        call tr_chkwrite("freq_r iq iw realomg trwv=", zw, iw, frr,nblochpmx, nbloch,ngb,iq)
1015  enddo iwloop
      call wv_close_iq_for_write()
    else  ! llw, Wing elements of W. See PRB81 125102
      iq0 = iq - nqibz
      vcou1 = fourpi/sum(q**2*tpioa**2) ! --> vcousq(1)**2!  !fourpi/sum(q**2*tpioa**2-eee)
      do 1115 iwblock = nwmin, nwmax, mpi__size_b
        iw = iwblock + mpi__rank_b
        if(iw > nwmax) cycle
        call stopwatch_start(t_sw_x_gather)
        zxqw(:,:) = zxq(:,:,iw)    ! all root_k have full shm_wvr; direct copy
        !$acc update device(zxqw)
        !$acc update host(zxqw(1,1))
        call stopwatch_pause(t_sw_x_gather)
        MToEBasisTransformation2: if(is_x0_m_basis) then
          call stopwatch_start(t_sw_x_m2e_xf)
          !$acc host_data use_device(zxqw, m2e_prod_basis, x_m2e)
          istat = gemm(zxqw, m2e_prod_basis, x_m2e, ngb, ngb, ngb)
          istat = gemm(m2e_prod_basis, x_m2e, zxqw, ngb, ngb, ngb, opA=m_op_C)
          !$acc end host_data
          call stopwatch_pause(t_sw_x_m2e_xf)
        endif MToEBasisTransformation2
        ix=0
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
        call stopwatch_start(t_sw_matinv)
        istat = zminv(epstilde(ix+1,ix+1), n=ngb-ix, lda=ngb)
        epstinv => epstilde
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
        nullify(epstinv)
        if(iq0<=nq0i) write(stdo,"('epsWVR: iq iw_R omg(iw) eps(wFC) eps(woLFC) ', &
             2i5,x,10(d13.6,2x,d13.6,x,d13.6,2x,d13.6,x,d13.6))") &
             iq,iw,freq_r(iw),llw(iw,iq0),1d0-vcou1*zxqw(1,1)
        continue               !iw
1115  enddo
      if(iq0 <=nq0i) call MPI__AllreduceSum(llw(nwmin,iq0), nwmax-nwmin+1, communicator=comm_root_k)
      if(ixyz(iq0)/=0) call MPI__AllreduceSum(wmuk(2,ixyz(iq0)), ngb-1, communicator=comm_root_k)
    endif
    !$acc exit data delete(zxqw, zw, vcousq, x_m2e)
    deallocate(zw, zxqw, epstilde, x_m2e)
    if(mpi__root_q) then
      call stopwatch_show(t_sw_x_gather)
      call stopwatch_show(t_sw_matinv)
      if(is_x0_m_basis .or. is_wc_m_basis) call stopwatch_show(t_sw_x_m2e_xf)
    endif
  end subroutine WVRllwR
  subroutine WVIllwI(q,iq,nmbas1,nmbas2,is_x0_m_basis,is_wc_m_basis)
    intent(in)::       q,iq,     nmbas1,nmbas2 !zxqi can be twiced when nspin=2
    integer:: nmbas1,nmbas2
    integer:: iq,iq0,nwmax,nwmin,iw,imode,ix,igb1,igb2,ifllwi
    real(8):: frr,q(3),vcou1
    logical::  localfieldcorrectionllw,cmdopt0
    logical, intent(in) :: is_x0_m_basis, is_wc_m_basis
    logical,save:: init=.true.
!    complex(8):: zxqi(nmbas1,nmbas2,niw)
    character(10):: i2char
    integer :: istat
    type(stopwatch) :: t_sw_matinv, t_sw_x_gather, t_sw_x_m2e_xf
    complex(kind=kp), allocatable :: zw(:,:), zxqw(:,:), x_m2e(:,:)
    complex(8), allocatable, target :: epstilde(:,:)
    complex(8), pointer :: epstinv(:,:) => null()
    integer :: iwblock, jw, irank, ierr
#ifdef __GPU
    attributes(device) :: epstinv, epstilde
#endif
    allocate(epstilde(ngb,ngb))
    allocate(zxqw(ngb, ngb))
    allocate(x_m2e(ngb,ngb))
    allocate(zw(nblochpmx,nblochpmx))
    !$acc enter data create(zxqw, x_m2e, zw) copyin(vcousq)
    if(init) then
       allocate(llwI(niw,nq0i), source=(0d0,0d0))
       init=.false.
    endif
    call stopwatch_init(t_sw_matinv, 'matinv')
    call stopwatch_init(t_sw_x_gather, 'gather')
    call stopwatch_init(t_sw_x_m2e_xf, 'xf chi: M2E')
    if(ipr)write(6,*)'WVRllwI: init'
    if (nspin == 1) then
      ! SHM: only rank 0 scales shared memory; barrier ensures others see the result.
      if (mpi__rank_root_k == 0) zxqi(:,:,:) = 2d0*zxqi(:,:,:)
      call MPI_barrier(comm_root_k, ierr)
    endif
    if( iq<=nqibz ) then
       call wv_open_iq_imag_for_write(iq, comm=comm_root_k)
       ix = merge(1, 0, iq == 1)
       do 1016 iwblock = 1, niw, mpi__size_b
          iw = iwblock + mpi__rank_b
          if(iw > niw) cycle
          !!  Eqs.(37),(38) in PRB81 125102
          !$acc kernels
          zw(:,:) = (0_kp, 0_kp)
          !$acc end kernels
          call stopwatch_start(t_sw_x_gather)
          zxqw(:,:) = zxqi(:,:,iw)    ! all root_k have full shm_wvi; direct copy
          call stopwatch_pause(t_sw_x_gather)
          MToEBasisTransformation1: if(is_x0_m_basis) then
            call stopwatch_start(t_sw_x_m2e_xf)
            !$acc host_data use_device(zxqw, m2e_prod_basis, x_m2e)
            istat = gemm(zxqw, m2e_prod_basis, x_m2e, ngb, ngb, ngb)
            istat = gemm(m2e_prod_basis, x_m2e, zxqw, ngb, ngb, ngb, opA=m_op_C)
            !$acc end host_data
            call stopwatch_pause(t_sw_x_m2e_xf)
          endif MToEBasisTransformation1
          !$acc kernels loop independent collapse(2)
          do igb2=ix+1,ngb
             do igb1=ix+1,ngb
                epstilde(igb1,igb2)= -vcousq(igb1)*zxqw(igb1,igb2)*vcousq(igb2)
                if(igb1==igb2) epstilde(igb1,igb2)=1+epstilde(igb1,igb2)
             enddo
          enddo
          !$acc end kernels
          call stopwatch_start(t_sw_matinv)
          istat = zminv(epstilde(ix+1,ix+1), n=ngb-ix, lda=ngb)
          epstinv => epstilde
          call stopwatch_pause(t_sw_matinv)
          !$acc kernels loop independent collapse(2)
          do igb2=ix+1,ngb
             do igb1=ix+1,ngb
                zw(igb1,igb2)= vcousq(igb1)*epstinv(igb1,igb2)*vcousq(igb2)
                if(igb1==igb2) zw(igb1,igb2)= zw(igb1,igb2)-vcousq(igb1)*vcousq(igb2)
             enddo
          enddo
          !$acc end kernels
          nullify(epstinv)
          ! W^[M] = AW^[E]A^\dagger except iq==1
          EtoMBasisTransformation: if(is_wc_m_basis .and. iq /= 1) then
            call stopwatch_start(t_sw_x_m2e_xf)
            !$acc host_data use_device(zw, m2e_prod_basis, x_m2e)
            istat = gemm(zw, m2e_prod_basis, x_m2e, ngb, ngb, ngb, opB=m_op_C, ldA=nblochpmx)
            istat = gemm(m2e_prod_basis, x_m2e, zw, ngb, ngb, ngb, ldC=nblochpmx)
            !$acc end host_data
            call stopwatch_pause(t_sw_x_m2e_xf)
          endif EtoMBasisTransformation
          !$acc update host(zw)
          ! write(ifrcwi, rec= iw)  zw !  WP = vsc-v
          call wv_put_imag(iw, zw(1:nblochpmx, 1:nblochpmx))
          call tr_chkwrite("freq_i iq iw imgomg trwv=",zw,iw,freq_i(iw),nblochpmx,nbloch,ngb,iq)
1016   enddo
       call wv_close_iq_for_write()
    else
       !! Full inversion to calculalte eps with LFC.
       iq0 = iq - nqibz
       vcou1 = fourpi/sum(q**2*tpioa**2) ! --> vcousq(1)**2!  !fourpi/sum(q**2*tpioa**2-eee)
       do 1116 iwblock = 1, niw, mpi__size_b
          iw = iwblock + mpi__rank_b
          if(iw > niw) cycle
          !if(localfieldcorrectionllw()) then
          call stopwatch_start(t_sw_x_gather)
          zxqw(:,:) = zxqi(:,:,iw)    ! all root_k have full shm_wvi; direct copy
          !$acc update host(zxqw(1,1))
          call stopwatch_pause(t_sw_x_gather)
          MToEBasisTransformation2: if(is_x0_m_basis) then
            call stopwatch_start(t_sw_x_m2e_xf)
            !$acc host_data use_device(zxqw, m2e_prod_basis, x_m2e)
            istat = gemm(zxqw, m2e_prod_basis, x_m2e, ngb, ngb, ngb)
            istat = gemm(m2e_prod_basis, x_m2e, zxqw, ngb, ngb, ngb, opA=m_op_C)
            !$acc end host_data
          call stopwatch_pause(t_sw_x_m2e_xf)
          endif MToEBasisTransformation2
           ix=0
           !$acc kernels loop independent collapse(2)
           do igb2=ix+1,ngb
              do igb1=ix+1,ngb
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
           call stopwatch_start(t_sw_matinv)
           istat = zminv(epstilde(ix+1,ix+1), n=ngb-ix, lda=ngb)
           epstinv => epstilde
           call stopwatch_pause(t_sw_matinv)
           if(iq0<=nq0i) llwI(iw,iq0)= 1d0/epstinv(1,1) !copy to CPU
          !else
          !   if(iq0<=nq0i) llwI(iw,iq0)=  1d0 -vcou1*zxqi(1,1,iw)
          !endif
          ! if(iq0<=nq0i) write(6,"('iq iw_img eps(wLFC) eps(noLFC)',i4,i4,2f10.4,2x,2f10.4)") &
          !      iq,iw,llwI(iw,iq0),1d0-vcou1*zxqi(1,1,iw)
           if(iq0<=nq0i) write(stdo,"('iq iw_img eps(wLFC) eps(noLFC)',i4,i4,2f10.4,2x,2f10.4)") &
                iq,iw,llwI(iw,iq0),1d0-vcou1*zxqw(1,1)
           nullify(epstinv)
1116   enddo
       if(iq0 <=nq0i) call MPI__AllreduceSum(llwI(1,iq0), niw, communicator=comm_root_k)
    endif
    !$acc exit data delete(zxqw, zw, vcousq, x_m2e)
    deallocate(zxqw, zw, epstilde, x_m2e)
    if(mpi__root_q) then
      call stopwatch_show(t_sw_x_gather)
      call stopwatch_show(t_sw_matinv)
      if(is_x0_m_basis .or. is_wc_m_basis) call stopwatch_show(t_sw_x_m2e_xf)
    endif
  end subroutine WVIllWI
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

  subroutine llw_add_req(req)
    integer, intent(in) :: req
    llw_nreqs = llw_nreqs + 1
    if (llw_nreqs > size(llw_reqs)) call rx('MPI__llw: llw_reqs overflow — increase llw_reqs size')
    llw_reqs(llw_nreqs) = req
  end subroutine llw_add_req

  subroutine MPI__irecvllw_q(iq0, src, dest)
    ! Post non-blocking Irecv(s) for llw/llwI/wmuk of auxiliary q-point iq0.
    ! Only dest rank participates; all others return immediately.
    use m_mpi, only: MPI__rank, MPI__size, comm
    integer, intent(in) :: iq0, src, dest
    integer :: ierr, req, tag
    if (MPI__size == 1 .or. src == dest .or. MPI__rank /= dest) return
    tag = iq0 * 4
    if (iq0 <= nq0i) then
      call MPI_Irecv(llw(nw_i,iq0),  nw-nw_i+1, MPI_COMPLEX16, src, tag,   comm, req, ierr); call llw_add_req(req)
      call MPI_Irecv(llwI(1,iq0),    niw,        MPI_COMPLEX16, src, tag+1, comm, req, ierr); call llw_add_req(req)
    end if
    if (ixyz(iq0) /= 0) then
      call MPI_Irecv(wmuk(2,ixyz(iq0)), ngbq0-1, MPI_COMPLEX16, src, tag+2, comm, req, ierr); call llw_add_req(req)
    end if
  end subroutine MPI__irecvllw_q

  subroutine MPI__isendllw_q(iq0, src, dest)
    ! Post non-blocking Isend(s) for llw/llwI/wmuk of auxiliary q-point iq0.
    ! Only src rank participates; all others return immediately.
    use m_mpi, only: MPI__rank, MPI__size, comm
    integer, intent(in) :: iq0, src, dest
    integer :: ierr, req, tag
    if (MPI__size == 1 .or. src == dest .or. MPI__rank /= src) return
    tag = iq0 * 4
    if (iq0 <= nq0i) then
      call MPI_Isend(llw(nw_i,iq0),  nw-nw_i+1, MPI_COMPLEX16, dest, tag,   comm, req, ierr); call llw_add_req(req)
      call MPI_Isend(llwI(1,iq0),    niw,        MPI_COMPLEX16, dest, tag+1, comm, req, ierr); call llw_add_req(req)
    end if
    if (ixyz(iq0) /= 0) then
      call MPI_Isend(wmuk(2,ixyz(iq0)), ngbq0-1, MPI_COMPLEX16, dest, tag+2, comm, req, ierr); call llw_add_req(req)
    end if
  end subroutine MPI__isendllw_q

  subroutine MPI__waitllw()
    ! Wait for all pending non-blocking llw transfers to complete.
    integer :: ierr
    if (llw_nreqs == 0) return
    call MPI_Waitall(llw_nreqs, llw_reqs, MPI_STATUSES_IGNORE, ierr)
    llw_nreqs = 0
  end subroutine MPI__waitllw

  subroutine MPI__sendllw_q(iq0, src, dest)
    ! Send/recv llw/llwI/wmuk for one auxiliary q-point iq0 (1-based).
    ! src: rank that computed iq0; dest: rank to deliver to (typically 0).
    ! Only src and dest participate; all other ranks return immediately.
    use m_mpi,only: MPI__DbleCOMPLEXsendQ,MPI__DbleCOMPLEXrecvQ,MPI__size,MPI__rank
    integer, intent(in) :: iq0, src, dest
    if(MPI__size==1 .or. src==dest) return
    if(MPI__rank == src) then
      if(iq0 <= nq0i) then
        call MPI__DbleCOMPLEXsendQ(llw(nw_i,iq0),(nw-nw_i+1),dest)
        call MPI__DbleCOMPLEXsendQ(llwI(1,iq0),niw,dest)
      endif
      if(ixyz(iq0)/=0) call MPI__DbleCOMPLEXsendQ(wmuk(2:ngbq0,ixyz(iq0)),ngbq0-1,dest)
    elseif(MPI__rank == dest) then
      if(iq0 <= nq0i) then
        call MPI__DbleCOMPLEXrecvQ(llw(nw_i,iq0),(nw-nw_i+1),src)
        call MPI__DbleCOMPLEXrecvQ(llwI(1,iq0),niw,src)
      endif
      if(ixyz(iq0)/=0) call MPI__DbleCOMPLEXrecvQ(wmuk(2:ngbq0,ixyz(iq0)),ngbq0-1,src)
    endif
  end subroutine MPI__sendllw_q
end module m_llw
!===================================================================
subroutine tr_chkwrite(tagname,zw,iw,freqq,nblochpmx,nbloch,ngb,iq)
  use m_lgunit,only:stdo
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
  write(stdo,'(a,f10.4,2i5,4d22.14)')tagname,freqq,iq,iw,trwv,trwv2
  call flush(6)
end subroutine tr_chkwrite
