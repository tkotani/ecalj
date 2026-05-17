!> get zxq and zxqi for given q
module m_x0kf
  use m_lgunit,only: stdo
  use m_keyvalue,only : Getkeyvalue
  use m_GWinput, only: gwinput_init, gwinput_loaded, tg_zmel_max_size => zmel_max_size
  use m_pkm4crpa,only : Readpkm4crpa
  use m_zmel,only: build_zmel, zmel
  use m_freq,only: npm, nwhis
  use m_struct_from_lmf,only: nsp=>nspin, nband; use m_gw_product_basis,only: ndima; use m_core_state,only: nctot
  use m_read_bzdata,only:  nqbz,ginv,nqibz,  rk=>qbz,wk=>wbz
  use m_rdpp,only: nbloch
  use m_readqg,only: ngpmx,ngcmx
  use m_qbze,only: nqbze
  use m_hamindex,only: ngrp
  use m_tetwt,only:  gettetwt,tetdeallocate, whw,ihw,nhw,jhw,n1b,n2b,nbnb,nbnbx,nhwtot
  use m_ftox
  use m_readVcoud,only:   vcousq,zcousq,ngb,ngc
  use m_kind,only: kp => kindrcxq
  use m_mpi,only: ipr, mpi__root_k => mpi__root_k_xq
  use m_wv_storage, only: shm_wvr, shm_wvi, wv_ngb
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
  public:: x0kf_zxq, deallocatezxq, deallocatezxqi
  complex(kind=kp), public, pointer:: zxq(:,:,:) => null()
  complex(kind=kp), public, pointer, contiguous :: zxqi(:,:,:) => null()
  complex(kind=kp), allocatable :: rcxq(:,:,:)
#ifdef __GPU
  attributes(device) :: rcxq
#endif
  integer,public::npr
  private

  integer:: ncount,ncoun
  integer,allocatable:: nkmin(:), nkmax(:),nkqmin(:),nkqmax(:),kc(:)
  integer,allocatable:: icounkmin(:),icounkmax(:)

  real(8),public,allocatable:: whwc(:)
  integer,allocatable,public:: iwini(:),iwend(:),itc(:),itpc(:),jpmc(:),icouini(:)



  integer,public::icounkmink,icounkmaxk
  logical,external:: cmdopt0
  logical:: debug = .false.
contains

  function X0kf_v4hz_init(job, q, isp_k, isp_kq, iq, crpa, ikbz_in, fkbz_in) result(ierr)
    implicit none
    integer, intent(in) :: job, isp_k, isp_kq, iq
    real(8), intent(in) :: q(3)
    logical, intent(in) :: crpa
    integer, intent(in), optional :: ikbz_in, fkbz_in
    integer :: ierr, jpm, ibib, iw, it, itp, icount, ncc, icoun, k
    real(8) :: imagweight, wpw_k, wpw_kq
    integer :: ikbz, fkbz
    ikbz = 1
    fkbz = nqbz
    if (present(ikbz_in) .and. present(fkbz_in)) then
      ikbz = ikbz_in
      fkbz = fkbz_in
    endif
    if (ipr) write(stdo,'(" x0kf_v4hz_init: job q =",i3,3f8.4)') job, q
    ierr = -1
    ncc = merge(0, nctot, npm==1)
    if (job==0) then
      if (allocated(nkmin))  deallocate(nkmin)
      if (allocated(nkqmin)) deallocate(nkqmin)
      if (allocated(nkmax))  deallocate(nkmax)
      if (allocated(nkqmax)) deallocate(nkqmax)
      allocate(nkmin(ikbz:fkbz), nkqmin(ikbz:fkbz), source= 999999)
      allocate(nkmax(ikbz:fkbz), nkqmax(ikbz:fkbz), source=-999999)
    endif
    if (job==1) then
      if (allocated(whwc))     deallocate(whwc)
      if (allocated(kc))       deallocate(kc)
      if (allocated(iwini))    deallocate(iwini)
      if (allocated(iwend))    deallocate(iwend)
      if (allocated(itc))      deallocate(itc)
      if (allocated(itpc))     deallocate(itpc)
      if (allocated(jpmc))     deallocate(jpmc)
      if (allocated(icouini))  deallocate(icouini)
      allocate(whwc(ncount), kc(ncoun), iwini(ncoun), iwend(ncoun), &
               itc(ncoun), itpc(ncoun), jpmc(ncoun), icouini(ncoun))
      if (allocated(icounkmin)) deallocate(icounkmin)
      if (allocated(icounkmax)) deallocate(icounkmax)
      allocate(icounkmin(ikbz:fkbz), icounkmax(ikbz:fkbz))
    endif
    icount = 0
    icoun  = 0
    do k = 1, nqbz
      if (k < ikbz .OR. k > fkbz) cycle
      if (job==1) icounkmin(k) = icoun+1
      if (job==0) then
        do jpm = 1, npm
          do ibib = 1, nbnb(k,jpm)
            nkmin(k)  = min(n1b(ibib,k,jpm), nkmin(k))
            nkqmin(k) = min(n2b(ibib,k,jpm), nkqmin(k))
            if (n1b(ibib,k,jpm) <= nband) nkmax(k)  = max(n1b(ibib,k,jpm), nkmax(k))
            if (n2b(ibib,k,jpm) <= nband) nkqmax(k) = max(n2b(ibib,k,jpm), nkqmax(k))
          enddo
        enddo
      endif
      flush(stdo)
      if (npm==2 .AND. nkqmin(k)/=1) call rx(" When npm==2, nkqmin==1 should be.")
      do jpm = 1, npm
        do ibib = 1, nbnb(k,jpm)
          if (ihw(ibib,k,jpm)+nhw(ibib,k,jpm)-1 > nwhis) call rx("x0kf_v4hz: iw>nwhis")
          if (n1b(ibib,k,jpm) > nkmax(k) ) cycle
          if (n2b(ibib,k,jpm) > nkqmax(k)) cycle
          it  = merge(nctot+n1b(ibib,k,jpm),             n1b(ibib,k,jpm)-nband,             n1b(ibib,k,jpm)<=nband)
          itp = merge(ncc  +n2b(ibib,k,jpm)-nkqmin(k)+1, n2b(ibib,k,jpm)-nkqmin(k)+1-nband, n2b(ibib,k,jpm)<=nband)
          if (crpa) then
            wpw_k  = merge(readpkm4crpa(n1b(ibib,k,jpm),   rk(:,k), isp_k),  0d0, n1b(ibib,k,jpm)<=nband)
            wpw_kq = merge(readpkm4crpa(n2b(ibib,k,jpm), q+rk(:,k), isp_kq), 0d0, n2b(ibib,k,jpm)<=nband)
          endif
          icoun = icoun+1
          if (job==1) then
            kc    (icoun) = k
            itc   (icoun) = it
            itpc  (icoun) = itp
            jpmc  (icoun) = jpm
            iwini (icoun) = ihw(ibib,k,jpm)
            iwend (icoun) = ihw(ibib,k,jpm)+nhw(ibib,k,jpm)-1
            icouini(icoun) = icount+1
          endif
          do iw = ihw(ibib,k,jpm), ihw(ibib,k,jpm)+nhw(ibib,k,jpm)-1
            imagweight = whw(jhw(ibib,k,jpm)+iw-ihw(ibib,k,jpm))
            if (crpa) imagweight = imagweight*(1d0-wpw_k*wpw_kq)
            icount = icount+1
            if (job==1) whwc(icount) = imagweight
          enddo
        enddo
      enddo
      if (job==1) icounkmax(k) = icoun
    enddo
    ncount = icount
    ncoun  = icoun
    if (job==0 .and. ipr) write(stdo,"('x0kf_v4hz_init: job=0 ncount ncoun nqibz=',3i8)") ncount, ncoun, nqibz
    ierr = 0
  end function x0kf_v4hz_init

  subroutine x0kf_zxq(realomega, imagomega, q, iq, npr, schi, crpa, chipm, nolfco, q00, zzr, is_m_basis)
    use m_readgwinput,only: ecut, ecuts
    use m_dpsion,only: dpsion5, dpsion_init, &
                      dpsion_chiq        => dpsion_chiq_h, &
                      dpsion_setup_rcxq  => dpsion_setup_rcxq_h
    use m_freq,only: nw_i, nw_w=>nw, niwt=>niw
    use m_freq,only: nw, niw
    use m_readeigen,only:readeval
    use m_zmel,only: set_m2e_prod_basis, set_m2e_prod_basis_chipm
    use m_stopwatch
    use m_readVcoud, only: ReleaseZcousq
    use m_mpi,only: mpi__rank_b => mpi__rank_b_xq, mpi__size_b => mpi__size_b_xq, &
                    mpi__root_b => mpi__root_b_xq, comm_b => comm_b_xq, comm_q, &
                    mpi__rank_k => mpi__rank_k_xq, mpi__size_k => mpi__size_k_xq, &
                    mpi__root_k => mpi__root_k_xq, mpi__rank_root_k => mpi__rank_root_k_xq, &
                    comm_k => comm_k_xq, comm_root_k => comm_root_k_xq
#ifdef __MP
    use m_mpi,only: MPI__reduceSum => MPI__reduceSum_c
#else
    use m_mpi,only: MPI__reduceSum
#endif
    use m_gpu, only: use_gpu
    use mpi
    implicit none
    intent(in)::      realomega, imagomega, q, iq, npr, schi, crpa, chipm, nolfco, q00, zzr
    logical:: realomega, imagomega, crpa, chipm, nolfco, is_m_basis
    integer:: iq, isp_k, isp_kq, ix0, is, isf, kx, ierr, npr, k, k_lo, k_hi
    integer:: iw_lo, iw_hi, iw_chunk
    real(8),optional:: q00(3)
    complex(8),optional:: zzr(:,:)
    real(8):: q(3), schi, ekxx1(nband,nqbz), ekxx2(nband,nqbz)
    character(10) :: i2char
    logical :: tetwtk = .false.
    real(8) :: zmel_max_size
    type(stopwatch) :: t_sw_zmel, t_sw_x0, t_sw_dpsion

    ! Omega parallelism: split flat range (1-npm)*nwhis:nwhis across comm_b ranks.
    iw_chunk = (npm*nwhis + 1 + mpi__size_b - 1) / mpi__size_b
    iw_lo = (1-npm)*nwhis + mpi__rank_b * iw_chunk
    iw_hi = min(iw_lo + iw_chunk - 1, nwhis)
    ! k-point parallelism: split nqbz across comm_k ranks.
    k_lo = mpi__rank_k * ((nqbz + mpi__size_k - 1) / mpi__size_k) + 1
    k_hi = min(k_lo + (nqbz + mpi__size_k - 1) / mpi__size_k - 1, nqbz)
    write(stdo,*) 'x0kf_zxq: k_lo, k_hi, nwhis, nw_i, nw, iw_lo, iw_hi', k_lo, k_hi, nwhis, nw_i, nw, iw_lo, iw_hi
    if (npm /= 1)      call rx('x0kf_zxq: npm/=1 not supported')
    if (wv_ngb /= npr) call rx('x0kf_zxq: wv_ngb /= npr (shm_wvr size mismatch)')

    if (cmdopt0('--tetwtk')) tetwtk = .true.
    call gwinput_init()
    if (gwinput_loaded) then
      zmel_max_size = tg_zmel_max_size
    else
      call rx('m_GWinput: legacy GWinput reader is disabled. GWinput.toml is required.')
    endif
    if (zmel_max_size < 0.001d0) zmel_max_size = 1d0
    if (chipm .AND. nolfco) then
      call set_m2e_prod_basis_chipm(zzr, npr)
    else
      call set_m2e_prod_basis(npr=npr)
    endif
    call ReleaseZcousq()
    if (associated(zxq)) nullify(zxq)
    if (allocated(rcxq)) deallocate(rcxq)
    if (nw_w > nwhis) call rx('nwhis is smaller than nw_w')
    if (ipr) write(stdo,ftox)' size of rcxq:', npr, nwhis*npm+1
    call flush(stdo)
    allocate(rcxq(1:npr, 1:npr, iw_lo:iw_hi))
    !$acc kernels
    rcxq = (0_kp, 0_kp)
    !$acc end kernels
    debug = cmdopt0('--debugzmel')
    isloop: do isp_k = 1, nsp
      GETtetrahedronWeight: block
        isp_kq = merge(3-isp_k, isp_k, chipm)
        do kx = 1, nqbz
          ekxx1(1:nband,kx) = readeval(  rk(:,kx), isp_k )
          ekxx2(1:nband,kx) = readeval(q+rk(:,kx), isp_kq)
        enddo
        if (.not.tetwtk) then
          call gettetwt(q, iq, isp_k, isp_kq, ekxx1, ekxx2, nband=nband)
          ierr = x0kf_v4hz_init(0, q, isp_k, isp_kq, iq, crpa)
          ierr = x0kf_v4hz_init(1, q, isp_k, isp_kq, iq, crpa)
          call tetdeallocate()
        endif
      end block GETtetrahedronWeight
      x0kf_v4hz_block: block
        use m_mem,only: writemem
        integer:: k, jpm, ibib, iw, igb2, igb1, it, itp, nkmax1, nkqmax1, ib1, ib2, ngcx, ix, iy, igb
        integer:: izmel, nmtot, nqtot, iwmax, ifi0, icoucold, icoun, icount, kold
        real(8):: imagweight, wpw_k, wpw_kq, qa, q0a
        complex(8):: img=(0d0,1d0)
        call cputid(0)
        call stopwatch_init(t_sw_zmel, 'zmel_gemm')
        call stopwatch_init(t_sw_x0,   'x0_gemm')
        kloop: do k = k_lo, k_hi
          if (tetwtk) then
            call gettetwt(q, iq, isp_k, isp_kq, ekxx1, ekxx2, nband=nband, ikbz_in=k, fkbz_in=k)
            ierr = x0kf_v4hz_init(0, q, isp_k, isp_kq, iq, crpa, ikbz_in=k, fkbz_in=k)
            ierr = x0kf_v4hz_init(1, q, isp_k, isp_kq, iq, crpa, ikbz_in=k, fkbz_in=k)
            call tetdeallocate()
          endif
          icounkmink = icounkmin(k)
          icounkmaxk = icounkmax(k)
          if (debug.and.ipr) write(stdo,ftox) 'ggggggggg goto build_zmel', k, nkmin(k), nkmax(k), nctot
          NMBATCH: block
            integer :: nsize, nns, ibatch, nbatch, ns12, ns1, ns2
            integer, allocatable :: ns1lists(:), ns2lists(:)
            nsize = (nkqmax(k)-nkqmin(k))*npr
            nns   = (nkmax(k) - nkmin(k) + 1)
            nbatch = ceiling(dble(nns)*nsize*16/1000**3/zmel_max_size)
            allocate(ns1lists(nbatch), ns2lists(nbatch))
            ns1 = nkmin(k) + nctot
            do ibatch = 1, nbatch
              ns12 = (nns + ibatch - 1)/nbatch
              ns1lists(ibatch) = ns1
              ns2lists(ibatch) = ns1 + ns12 - 1
              ns1 = ns2lists(ibatch) + 1
            enddo
            do ibatch = 1, nbatch
              ns1  = ns1lists(ibatch)
              ns2  = ns2lists(ibatch)
              ns12 = ns2 - ns1 + 1
              if (ns12 == 0) cycle
              call stopwatch_start(t_sw_zmel)
              if (ipr) write(stdo,ftox) 'zmel_batch:', ibatch, ns1, ns2, nbatch
              if (debug) call writemem('xxxx start build_zmel')
              call build_zmel(q=q+rk(:,k), kvec=q, irot=1, rkvec=q, ns1=ns1, ns2=ns2, ispm=isp_k, &
                   nqini=nkqmin(k), nqmax=nkqmax(k), ispq=isp_kq, nctot=nctot, ncc=merge(0,nctot,npm==1), &
                   zmelconjg=.true., is_m_basis=is_m_basis, mpi_mode=.not.use_gpu, comm=comm_b)
              if (debug) call writemem('xxxx end build_zmel')
              call stopwatch_pause(t_sw_zmel)
              call stopwatch_start(t_sw_x0)
              if (debug) call writemem('xxxx start x0gemm')
              call x0gemm(rcxq, npr, nwhis, npm, ns1, ns2, iw_lo, iw_hi)
              if (debug) call writemem('xxxx end of x0gemm')
              call stopwatch_pause(t_sw_x0)
            enddo
            deallocate(ns1lists, ns2lists)
          end block NMBATCH
          if (ipr) write(stdo,ftox) 'end of k:', k, ' of:', nqbz, &
              'zmel:', ftof(stopwatch_lap_time(t_sw_zmel),4), '(sec)', &
              ' x0:', ftof(stopwatch_lap_time(t_sw_x0),4), '(sec)'
          call flush(6)
        enddo kloop
        call stopwatch_show(t_sw_zmel)
        call stopwatch_show(t_sw_x0)
        call cputid(0)
        if (debug.and.ipr) write(stdo,ftox)"--- x0kf_v4hz: end: sumcheck abs(rcxq)=", sum(abs(rcxq(:,:,:)))
      end block x0kf_v4hz_block
      deallocate(whwc, kc, iwini, iwend, itc, itpc, jpmc, icouini, nkmin, nkmax, nkqmin, nkqmax, icounkmin, icounkmax)
      HilbertTransformation: if (isp_k==nsp .OR. chipm) then
        !Get real part. When chipm=T, do dpsion5 for every isp_k; When =F, do dpsion5 after rxcq accumulated for spins
        mpi_k_accumulate: block
          ! Unified: n_kpara=1 uses single-rank comm_k (reduce is no-op); GPU implicit device→host.
          integer :: iw
          do iw = iw_lo, iw_hi
            if (iw == 0) then
              ! iw=0 is never accumulated (x0gemm skips it) but dpsion reads chi0(:,:,0).
              ! Zero shm_wvr to avoid stale W from the previous q-point.
              if (mpi__root_k) shm_wvr(:,:, iw - (1-npm)*nwhis + 1) = (0_kp, 0_kp)
              cycle
            end if
            block
              complex(kp) :: iw_slice(npr, npr)
              iw_slice = rcxq(:,:,iw)
              call MPI__reduceSum(0, iw_slice(1,1), npr*npr, communicator=comm_k)
              if (mpi__root_k) shm_wvr(:,:, iw - (1-npm)*nwhis + 1) = iw_slice
            end block
          enddo
          deallocate(rcxq)
          call MPI_barrier(comm_q, ierr)
        end block mpi_k_accumulate
        if (mpi__root_k) then
          block
            complex(kp), pointer, contiguous :: chi0(:,:,:)
            chi0(1:npr, 1:npr, (1-npm)*nwhis:nwhis) => shm_wvr
            if (mpi__rank_root_k == 0) then
              if (imagomega) zxqi(1:npr, 1:npr, 1:niw) => shm_wvi
              call stopwatch_init(t_sw_dpsion, 'dpsion')
              call stopwatch_start(t_sw_dpsion)
              call dpsion_init(realomega, imagomega, chipm)
              call dpsion_chiq(realomega, imagomega, chipm, chi0, zxqi, npr, npr, schi, isp_k, ecut)
              call stopwatch_pause(t_sw_dpsion)
              call stopwatch_show(t_sw_dpsion)
            endif
            call MPI_barrier(comm_root_k, ierr)
            if (realomega) zxq(1:,1:,nw_i:) => shm_wvr(1:npr, 1:npr, nw_i-(1-npm)*nwhis+1 : nw_w-(1-npm)*nwhis+1)
            if (chipm) then
              if (mpi__rank_root_k == 0) call dpsion_setup_rcxq(chi0, npr, npr, isp_k)
              call MPI_barrier(comm_root_k, ierr)
            endif
            if (imagomega) zxqi(1:npr, 1:npr, 1:niw) => shm_wvi
          end block
        endif
        if (chipm .and. isp_k /= nsp) then
          allocate(rcxq(1:npr, 1:npr, iw_lo:iw_hi))
          !$acc kernels
          rcxq = (0_kp, 0_kp)
          !$acc end kernels
        endif
      endif HilbertTransformation
    enddo isloop
  end subroutine x0kf_zxq

  subroutine deallocatezxq()
    nullify(zxq)
  end subroutine deallocatezxq

  subroutine deallocatezxqi()
    if (.not. associated(zxqi)) return
    nullify(zxqi)
  end subroutine deallocatezxqi

end module m_x0kf

!! === calculate chi0, or chi0_pm ===
!!
!! ppovl= <I|J> = O , V_IJ=<I|v|J>
!! (V_IJ - vcoud_mu O_IJ) Zcousq(J, mu)=0, where Z is normalized with O_IJ.
!! <I|v|J>= \sum_mu ppovl*zcousq(:,mu) v^mu (Zcousq^*(:,mu) ppovl)
!!
!! zmelt contains O^-1=<I|J>^-1 factor. Thus zmelt(phi phi J)= <phi |phi I> O^-1_IJ
!! ppovlz(I, mu) = \sum_J O_IJ Zcousq(J, mu)
!!
!!  rcxq (npr,npr,nwhis,npm): for given q,
!!       rcxq(I,J,iw,ipm) = Im (chi0(omega))= \sum_k <I_q psi_k|psi_(q+k)> <psi_(q+k)|psi_k> \delta(\omega- (e_i-ej))
!!        When npm=2 we calculate negative energy part. (time-reversal asymmetry)
!!
! note: zmel: matrix element <phi phi |M>
!! ndima   = total number of atomic basis functions within MT
!! nqbz    = number of k-points in the 1st BZ

! z1p = <M_ibg1 psi_it | psi_itp> < psi_itp | psi_it M_ibg2 >
!  zxq(iw,ibg1,igb2) = \sum_ibib imgw(iw,ibib)* z1p(ibib, igb1,igb2) !ibib means band pair (occ,unocc)
!
!  zzmel(1:nbloch, ib_k,ib_kq)
!      ib_k =[1:nctot]              core
!      ib_k =[nctot+nkmin:nctot+nkmax]  valence
!      ib_kq =[1:ncc]             core
!      ib_kq =[ncc+nkqmin:ncc+nkqmax]  valence range [nkqmin,nkqmax]
!   If jpm=1, ncc=0.  !   If jpm=2, ncc=ncore. nkqmin(k)=1 should be.
! NOTE:
!  q+rk n2b vec_kq  vec_kq_g geig_kq cphi_kq  ngp_kq ngvecp_kq  isp_kq
!    rk n1b vec_k   vec_k_g  geig_k  cphi_k   ngp_k  ngvecp_k   isp_k
!! -------------------------------------------------------------------------------
! note: for usual correlation mode, I think nctot=0
!!--- For dielectric funciton, we use irot=1 kvec=rkvec=q
!            < MPB      middle   |   end >
!!              q      rkvec     | q + rkvec
!                      nkmin:nt0 | nkqmin:ntp0
!                         occ    | unocc
!                      (nkmin=1)
!                      (cphi_k  | cphi_kq !in x0kf)
!!     rkvec= rk(:,k)   ! <phi(q+rk,nqmax)|phi(rk,nctot+nmmax)  MPB(q,ngb )>
!!     qbz_kr= rk(:,k)  !
!!     qibz_k= rk(:,k)  ! k
!! Get_zmelt in m_zmel gives the matrix element zmel,  ZO^-1 <MPB psi|psi> , where ZO is ppovlz
!! zmel(ngb, nctot+nt0,  ncc+ntp0) in m_zmel
!            nkmin:nt0, nkqmin:ntp0, where nt0=nkmax-nkmin+1  , ntp0=nkqmax-nkqmin+1
