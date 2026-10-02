!> get zxq and zxqi for given q
module m_x0kf
  use m_lgunit,only: stdo
  use m_keyvalue,only : Getkeyvalue
  use m_GWinput, only: gwinput_init, gwinput_loaded, tg_zmel_batch_gb => zmel_batch_gb
  use m_pkm4crpa,only : Readpkm4crpa
  use m_zmel,only: build_zmel, zmel
  use m_freq,only: npm, nwhis, frhis
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
use m_cmdopt_registry, only: c0_debugzmel, c0_tetwtk
  use iso_c_binding, only: c_int, c_char, c_null_char
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
  public:: x0kf_zxq, deallocatezxq, deallocatezxqi, x0kf_tetwt_write
  complex(kind=kp), public, pointer:: zxq(:,:,:) => null()
  complex(kind=kp), public, pointer, contiguous :: zxqi(:,:,:) => null()
  complex(kind=kp), allocatable :: rcxq(:,:,:)
#ifdef __GPU
  attributes(device) :: rcxq
#endif
  private
  integer :: ncount, ncoun
  integer, allocatable :: nkmin(:), nkmax(:),nkqmin(:),nkqmax(:),kc(:)
  integer, allocatable :: icounkmin(:), icounkmax(:)
  real(8), allocatable :: whwc(:)
  integer, allocatable :: iwini(:),iwend(:),itc(:),itpc(:),jpmc(:),icouini(:)
  logical, allocatable :: intrac(:)   ! pair is intraband (n1b == n2b): the Drude term, exempt from chi0_filterw
  real(8), allocatable :: fcw(:)      ! [gw] chi0_filterw: factor per histogram bin (1 = off)
  ! The tetrahedron weights in the form x0kf_zxq uses (the arrays x0kf_v4hz_init makes, all k) can come from a file
  ! __TETWT.<iq>.<isp> that another run of this program wrote on the CPU cores (hgw --tetwt_write, which gwsc starts
  ! next to hgw).  x0kf_zxq uses a file only when everything the weights depend on matches: the sizes, q, the band
  ! energies at every k and k+q (checksum), the histogram bins, E_F, kBT and the band cut (tetwt_pkey); otherwise
  ! (no file yet, stale, k split over ranks, cRPA, chi+-) it computes them.  Same code on the same input: bit-identical.
  ! (2026-09-27 04:24)
  integer, parameter :: tetwt_tag = 20260927   ! format tag of __TETWT.*: change it when the records of tetwt_save change
  interface
    integer(c_int) function c_rename(old, new) bind(C, name='rename')
      import :: c_int, c_char
      character(kind=c_char), intent(in) :: old(*), new(*)
    end function c_rename
  end interface
  logical :: debug = .false.
contains
  !> [gw] chi0_filterw = [wc, dw] (eV): factor 1/(1+exp((wc-|omega|)/dw)) of each histogram bin
  !! (bin center, frhis in Hartree), applied to the tetrahedron weights of the interband pairs
  !! (and of the intraband ones too when chi0_filterw_drude = false) in accumulate_chi0.
  !! Filtering the weights pair by pair is the same as filtering the Im chi0 histogram of those
  !! pairs; Re chi0(omega) and chi0(i omega) then follow from the same Im chi0 (dpsion5), so W stays
  !! consistent on both axes.  ecaljdoc/MD/research_log.md 2026-09-20 23:26.
  subroutine chi0_filterw_setup()
    use m_GWinput, only: chi0_filterw, chi0_filterw_set, chi0_filterw_drude
    use m_freq, only: frhis
    integer :: iw
    real(8) :: x
    if (allocated(fcw)) deallocate(fcw)
    allocate(fcw(nwhis), source = 1d0)
    call gwinput_init()
    if (.not. gwinput_loaded) return
    if (.not. chi0_filterw_set) return
    do iw = 1, nwhis
      x = (chi0_filterw(1) - 27.211386d0*.5d0*(frhis(iw)+frhis(iw+1)))/chi0_filterw(2)
      fcw(iw) = merge(0d0, merge(1d0, 1d0/(1d0+exp(x)), x < -40d0), x > 40d0)
    enddo
    if (ipr) write(stdo,'(" chi0_filterw (eV) wc dw =",2f8.3,"  drude kept =",l2, &
         "  : transitions below wc removed from chi0")') chi0_filterw, chi0_filterw_drude
  end subroutine chi0_filterw_setup

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
      if (allocated(intrac))   deallocate(intrac)
      allocate(whwc(ncount), kc(ncoun), iwini(ncoun), iwend(ncoun), &
               itc(ncoun), itpc(ncoun), jpmc(ncoun), icouini(ncoun), intrac(ncoun))
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
            intrac(icoun) = n1b(ibib,k,jpm) == n2b(ibib,k,jpm) .and. n1b(ibib,k,jpm) <= nband
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

  subroutine x0kf_zxq(realomega, imagomega, q, iq, npr, schi, crpa, chipm, nolfco, zzr, is_m_basis)
    use m_readgwinput,only: ecut, ecuts
    use m_dpsion,only: dpsion5, dpsion_init, &
                      dpsion_chiq        => dpsion_chiq_h, &
                      dpsion_setup_rcxq  => dpsion_setup_rcxq_h, &
                      dpsion_chiq_dev    => dpsion_chiq_d
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
    use m_GWinput, only: SmearX0
    use mpi
    implicit none
    intent(in)::      realomega, imagomega, q, iq, npr, schi, crpa, chipm, nolfco, zzr
    logical:: realomega, imagomega, crpa, chipm, nolfco, is_m_basis
    integer:: iq, isp_k, isp_kq, ix0, is, isf, kx, ierr, npr, k, k_lo, k_hi
    integer:: iw_lo, iw_hi, iw_chunk
    complex(8),optional:: zzr(:,:)
    real(8):: q(3), schi, ekxx1(nband,nqbz), ekxx2(nband,nqbz)
    character(10) :: i2char
    logical :: tetwtk = .false., hilbert_on_device, loaded
    real(8) :: zmel_batch_gb
    type(stopwatch) :: t_sw_zmel, t_sw_x0, t_sw_dpsion

    ! Omega parallelism: split flat range (1-npm)*nwhis:nwhis across comm_b ranks.
    iw_chunk = (npm*nwhis + 1 + mpi__size_b - 1) / mpi__size_b
    iw_lo = (1-npm)*nwhis + mpi__rank_b * iw_chunk
    iw_hi = min(iw_lo + iw_chunk - 1, nwhis)
    ! k-point parallelism: split nqbz across comm_k ranks.
    k_lo = mpi__rank_k * ((nqbz + mpi__size_k - 1) / mpi__size_k) + 1
    k_hi = min(k_lo + (nqbz + mpi__size_k - 1) / mpi__size_k - 1, nqbz)
    if (ipr) write(stdo,'(1X,A,7I6)') 'x0kf_zxq: k_lo k_hi nwhis nw_i nw iw_lo iw_hi =', &
                                       k_lo, k_hi, nwhis, nw_i, nw, iw_lo, iw_hi
    if (npm /= 1)      call rx('x0kf_zxq: npm/=1 not supported')
    if (wv_ngb /= npr) call rx('x0kf_zxq: wv_ngb /= npr (shm_wvr size mismatch)')

    if (c0_tetwtk) tetwtk = .true.
    call gwinput_init()
    if (gwinput_loaded) then
      zmel_batch_gb = tg_zmel_batch_gb
    else
      call rx('m_GWinput: legacy GWinput reader is disabled; ctrlg.<sname>.toml is required.')
    endif
    if (zmel_batch_gb < 0.001d0) zmel_batch_gb = 0.4d0
    if (chipm .AND. nolfco) then
      call set_m2e_prod_basis_chipm(zzr, npr)
    else
      call set_m2e_prod_basis(npr=npr)
    endif
    call ReleaseZcousq()
    if (associated(zxq)) nullify(zxq)
    if (allocated(rcxq)) deallocate(rcxq)
    if (nw_w > nwhis) call rx('nwhis is smaller than nw_w')
    if (ipr) write(stdo,'(1X,A,I0,A,I0,A,F7.3,A)') &
        'rcxq(npr,npr,niw): npr=', npr, '  niw=', iw_hi-iw_lo+1, &
        '  mem=', real(npr,8)**2 * real(iw_hi-iw_lo+1,8) * real(2*kp,8) / 1d9, ' GB'
    call flush(stdo)
    allocate(rcxq(1:npr, 1:npr, iw_lo:iw_hi))
    !$acc kernels
    rcxq = (0_kp, 0_kp)
    !$acc end kernels
    call chi0_filterw_setup()   ! [gw] chi0_filterw (normally off): fcw = 1
    debug = c0_debugzmel
    isloop: do isp_k = 1, nsp
      GETtetrahedronWeight: block
        isp_kq = merge(3-isp_k, isp_k, chipm)
        do kx = 1, nqbz
          ekxx1(1:nband,kx) = readeval(  rk(:,kx), isp_k )
          ekxx2(1:nband,kx) = readeval(q+rk(:,kx), isp_kq)
        enddo
        if (.not.tetwtk) then
          ! Tetrahedron weights only for this rank's k slice k_lo:k_hi (MPI k-parallel over comm_k):
          ! tetwt5 skips the tetrahedra with no vertex in the range, so the cost is ~(k_hi-k_lo+1)/nqbz
          ! of the full loop.  With one k rank this is the old all-k call.  (2026-09-22; replaces the
          ! OpenMP variant.)  The weights of a k point do not depend on the range (the degenerate-band
          ! symmetrization acts within a k point), so the result is independent of the split.
          ! All k on this rank: the weights may already be in __TETWT.<iq>.<isp> (hgw --tetwt_write).
          loaded = .false.
          if (k_lo == 1 .and. k_hi == nqbz .and. .not.chipm .and. .not.crpa) then
            loaded = tetwt_load(iq, isp_k, q, ekxx1, ekxx2)
          endif
          if (.not.loaded) then
            call gettetwt(q, iq, isp_k, isp_kq, ekxx1, ekxx2, nband=nband, ikbz_in=k_lo, fkbz_in=k_hi)
            ierr = x0kf_v4hz_init(0, q, isp_k, isp_kq, iq, crpa, ikbz_in=k_lo, fkbz_in=k_hi)
            ierr = x0kf_v4hz_init(1, q, isp_k, isp_kq, iq, crpa, ikbz_in=k_lo, fkbz_in=k_hi)
            call tetdeallocate()
          endif
        endif
      end block GETtetrahedronWeight
      x0kf_v4hz_block: block
        use m_mem,only: writemem
        integer:: k, jpm, ibib, iw, igb2, igb1, it, itp, nkmax1, nkqmax1, ib1, ib2, ngcx, ix, iy, igb
        integer:: izmel, nmtot, nqtot, iwmax, ifi0, icoucold, icoun, icount, kold
        integer:: ng, nsg, nqg, kbmax, nq12, iq12, is12
        integer, allocatable :: gk(:), gns1(:)
        complex(kind=kp), allocatable :: zg(:,:,:,:)
#ifdef __GPU
        attributes(device) :: zg
#endif
        real(8):: imagweight, wpw_k, wpw_kq, qa, q0a
        complex(8):: img=(0d0,1d0)
        call cputid(0)
        call stopwatch_init(t_sw_zmel, 'zmel_gemm')
        call stopwatch_init(t_sw_x0,   'x0_gemm')
        ! Group up to kbmax k points whose zmel is built in one piece and accumulate them together
        ! (accumulate_chi0_group); a k point split into NMBATCH pieces, tetwtk, or npm=2 goes the old way.
        ng = 0
        if (.not. tetwtk .and. npm == 1) then
          nsg = maxval(nkmax(k_lo:k_hi) - nkmin(k_lo:k_hi) + 1)
          nqg = maxval(nkqmax(k_lo:k_hi) - nkqmin(k_lo:k_hi) + 1) + merge(0, nctot, npm==1)
          kbmax = int(min(8d0, 1d9/(real(npr,8)*nsg*nqg*2*kp)))
          if (kbmax >= 2) then
            allocate(zg(npr, nsg, nqg, kbmax), gk(kbmax), gns1(kbmax))
          endif
        endif
        kloop: do k = k_lo, k_hi
          if (tetwtk) then
            call gettetwt(q, iq, isp_k, isp_kq, ekxx1, ekxx2, nband=nband, ikbz_in=k, fkbz_in=k)
            ierr = x0kf_v4hz_init(0, q, isp_k, isp_kq, iq, crpa, ikbz_in=k, fkbz_in=k)
            ierr = x0kf_v4hz_init(1, q, isp_k, isp_kq, iq, crpa, ikbz_in=k, fkbz_in=k)
            call tetdeallocate()
          endif
          if (debug.and.ipr) write(stdo,ftox) 'ggggggggg goto build_zmel', k, nkmin(k), nkmax(k), nctot
          NMBATCH: block
            integer :: nsize, nns, ibatch, nbatch, ns12, ns1, ns2
            integer, allocatable :: ns1lists(:), ns2lists(:)
            nsize = (nkqmax(k)-nkqmin(k))*npr
            nns   = (nkmax(k) - nkmin(k) + 1)
            nbatch = ceiling(dble(nns)*nsize*16/1000**3/zmel_batch_gb)
            allocate(ns1lists(nbatch), ns2lists(nbatch))
            ns1 = nkmin(k) + nctot
            do ibatch = 1, nbatch
              ns12 = (nns + ibatch - 1)/nbatch
              ns1lists(ibatch) = ns1
              ns2lists(ibatch) = ns1 + ns12 - 1
              ns1 = ns2lists(ibatch) + 1
            enddo
            if (nbatch == 1 .and. allocated(zg)) then       ! the whole zmel of k in one piece: into the group
              ns1 = ns1lists(1)
              ns2 = ns2lists(1)
              if (ns2 >= ns1) then
                call stopwatch_start(t_sw_zmel)
                call build_zmel(q=q+rk(:,k), kvec=q, irot=1, rkvec=q, ns1=ns1, ns2=ns2, ispm=isp_k, &
                     nqini=nkqmin(k), nqmax=nkqmax(k), ispq=isp_kq, nctot=nctot, ncc=merge(0,nctot,npm==1), &
                     zmelconjg=.true., is_m_basis=is_m_basis, mpi_mode=.not.use_gpu, comm=comm_b)
                call stopwatch_pause(t_sw_zmel)
                ng = ng + 1
                gk(ng) = k
                gns1(ng) = ns1
                nq12 = nkqmax(k) - nkqmin(k) + 1 + merge(0, nctot, npm==1)
                !$acc parallel loop collapse(3) present(zmel)
                do iq12 = 1, nq12
                  do is12 = 1, ns2-ns1+1
                    do igb = 1, npr
                      zg(igb,is12,iq12,ng) = zmel(igb,ns1+is12-1,iq12)
                    enddo
                  enddo
                enddo
              endif
              if (ng == kbmax .or. (k == k_hi .and. ng > 0)) then
                call stopwatch_start(t_sw_x0)
                call accumulate_chi0_group(ng, gk, gns1, zg, iw_lo, iw_hi, npr)
                call stopwatch_pause(t_sw_x0)
                ng = 0
              endif
              nbatch = 0                                   ! the loop below does nothing
            elseif (ng > 0) then                           ! this k is split: accumulate the group first
              call stopwatch_start(t_sw_x0)
              call accumulate_chi0_group(ng, gk, gns1, zg, iw_lo, iw_hi, npr)
              call stopwatch_pause(t_sw_x0)
              ng = 0
            endif
            do ibatch = 1, nbatch
              ns1  = ns1lists(ibatch)
              ns2  = ns2lists(ibatch)
              ns12 = ns2 - ns1 + 1
              if (ns12 == 0) cycle
              if (ipr) write(stdo,ftox) 'zmel_batch:', ibatch, ns1, ns2, nbatch

              if (debug) call writemem('xxxx start build_zmel')
              call stopwatch_start(t_sw_zmel)
              call build_zmel(q=q+rk(:,k), kvec=q, irot=1, rkvec=q, ns1=ns1, ns2=ns2, ispm=isp_k, &
                   nqini=nkqmin(k), nqmax=nkqmax(k), ispq=isp_kq, nctot=nctot, ncc=merge(0,nctot,npm==1), &
                   zmelconjg=.true., is_m_basis=is_m_basis, mpi_mode=.not.use_gpu, comm=comm_b)
              call stopwatch_pause(t_sw_zmel)
              if (debug) call writemem('xxxx end build_zmel')

              if (debug) call writemem('xxxx start accumulate_chi0')
              call stopwatch_start(t_sw_x0)
              call accumulate_chi0(ns1, ns2, iw_lo, iw_hi, npr, icounkmin(k), icounkmax(k))
              call stopwatch_pause(t_sw_x0)
              if (debug) call writemem('xxxx end of accumulate_chi0')
            enddo
            deallocate(ns1lists, ns2lists)
          end block NMBATCH
          if (ipr) write(stdo,ftox) 'end of k:', k, ' of:', nqbz, &
              'zmel:', ftof(stopwatch_lap_time(t_sw_zmel),4), '(sec)', &
              ' x0:', ftof(stopwatch_lap_time(t_sw_x0),4), '(sec)'
          if(ipr) call flush(stdo)
        enddo kloop
        if (ng > 0) then                                   ! a group left over (last k was split)
          call stopwatch_start(t_sw_x0)
          call accumulate_chi0_group(ng, gk, gns1, zg, iw_lo, iw_hi, npr)
          call stopwatch_pause(t_sw_x0)
        endif
        if (allocated(zg)) deallocate(zg, gk, gns1)
        call stopwatch_show(t_sw_zmel)
        call stopwatch_show(t_sw_x0)
        call cputid(0)
        if (debug.and.ipr) then
          block
            real(kp) :: sumcheck
            !$acc kernels
            sumcheck = sum(abs(rcxq(:,:,:)))
            !$acc end kernels
            write(stdo,ftox)"--- x0kf_v4hz: end: sumcheck abs(rcxq)=", sumcheck
          end block
        endif
      end block x0kf_v4hz_block
      deallocate(whwc, kc, iwini, iwend, itc, itpc, jpmc, icouini, intrac, nkmin, nkmax, nkqmin, nkqmax, icounkmin, icounkmax)
      HilbertTransformation: if (isp_k==nsp .OR. chipm) then
        !Get real part. When chipm=T, do dpsion5 for every isp_k; When =F, do dpsion5 after rxcq accumulated for spins
        ! One rank per q-group (no k or omega split) on a GPU holds all of chi0 in rcxq on the device: do the
        ! Hilbert transform there and copy the result into SHM.  The host path copies the histogram to SHM and
        ! transforms it on one CPU core while the GPU waits (12 s per q for npr~1050, nwhis~330; 2026-09-27).
        hilbert_on_device = use_gpu .and. mpi__size_k == 1 .and. mpi__size_b == 1 .and. .not.chipm
        if (hilbert_on_device) then
          block
            complex(kind=kp), allocatable :: zxqi_d(:,:,:)
#ifdef __GPU
            attributes(device) :: zxqi_d
#endif
            integer :: iw
            call stopwatch_init(t_sw_dpsion, 'dpsion(device)')
            call stopwatch_start(t_sw_dpsion)
            allocate(zxqi_d(npr, npr, niw))
            call dpsion_init(realomega, imagomega, chipm)
            call dpsion_chiq_dev(realomega, imagomega, chipm, rcxq, zxqi_d, npr, npr, schi, isp_k, ecut)
            do iw = iw_lo, iw_hi   ! all slices, as the host path leaves them (iw=0: chi0 at omega=0)
              shm_wvr(:,:, iw - (1-npm)*nwhis + 1) = rcxq(:,:,iw)
            enddo
            if (imagomega) shm_wvi(:,:,1:niw) = zxqi_d(:,:,1:niw)
            deallocate(zxqi_d)
            deallocate(rcxq)
            call stopwatch_pause(t_sw_dpsion)
            call stopwatch_show(t_sw_dpsion)
          end block
          call MPI_barrier(comm_q, ierr)
          if (realomega) zxq(1:,1:,nw_i:) => shm_wvr(1:npr, 1:npr, nw_i-(1-npm)*nwhis+1 : nw_w-(1-npm)*nwhis+1)
          if (imagomega) zxqi(1:npr, 1:npr, 1:niw) => shm_wvi
        else
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
        endif ! hilbert_on_device
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

  subroutine accumulate_chi0(ns1, ns2, iw_lo, iw_hi, npr, icounkmink, icounkmaxk)
    use m_blas, only: m_op_c
    use m_GWinput, only: chi0_filterw_drude
#ifdef __GPU
    use openacc
    use cudafor
#endif
    implicit none
    integer, intent(in) :: ns1, ns2, iw_lo, iw_hi, npr, icounkmink, icounkmaxk
    integer :: icoun, igb1, igb2, iw, jpm, iw_pos, it, itp, ittp, nttp_max, ierr
    integer :: pos_lo(2), pos_hi(2)
    integer, allocatable :: nttp(:,:), itw(:,:,:), itpw(:,:,:)
    complex(kind=kp), allocatable :: zw(:,:), wzw(:,:)
    complex(kind=kp), parameter :: CONE = (1_kp, 0_kp)
    real(8), allocatable :: hilbert_w(:,:,:)
#ifdef __GPU
    attributes(device) :: zw, wzw
#endif
    pos_lo(1) = max(iw_lo, 1);    pos_hi(1) = min(iw_hi, nwhis)
    pos_lo(2) = max(1, -iw_hi);   pos_hi(2) = min(nwhis, -iw_lo)

    allocate(nttp(nwhis,npm), source = 0)
    do icoun = icounkmink, icounkmaxk
      jpm = jpmc(icoun)
      do iw = max(iwini(icoun), pos_lo(jpm)), min(iwend(icoun), pos_hi(jpm))
        nttp(iw,jpm) = nttp(iw,jpm) + 1
      enddo
    enddo

    nttp_max = maxval(nttp(1:nwhis,1:npm))
    if(debug) write(stdo, ftox)'nttp_max = ', nttp_max
    allocate(itw(nttp_max,nwhis,npm), source = 0)
    allocate(itpw(nttp_max,nwhis,npm), source = 0)
    allocate(hilbert_w(nttp_max,nwhis,npm), source = 0d0)

    nttp(1:nwhis,1:npm) = 0
    do icoun = icounkmink, icounkmaxk
      jpm = jpmc(icoun)
      it  = itc(icoun)
      itp = itpc(icoun)
      if(it < ns1 .or. it > ns2) cycle
      do iw = max(iwini(icoun), pos_lo(jpm)), min(iwend(icoun), pos_hi(jpm))
        nttp(iw,jpm) = nttp(iw,jpm) + 1
        ittp = nttp(iw,jpm)
        itw(ittp,iw,jpm)       = it
        itpw(ittp,iw,jpm)      = itp
        hilbert_w(ittp,iw,jpm) = whwc(iw-iwini(icoun)+icouini(icoun))
        if (.not.(intrac(icoun) .and. chi0_filterw_drude)) hilbert_w(ittp,iw,jpm) = hilbert_w(ittp,iw,jpm)*fcw(iw)
      enddo
    enddo

    allocate(zw(nttp_max,npr), wzw(nttp_max,npr))
    !$acc data copyin(hilbert_w, itw, itpw) present(zmel)
    do iw = iw_lo, iw_hi
      if (iw == 0) cycle
      if (iw > 0) then; jpm = 1; iw_pos = iw
      else;             jpm = 2; iw_pos = -iw
      endif
      if (nttp(iw_pos,jpm) < 1) cycle
      !$acc kernels loop independent collapse(2)
      do ittp = 1, nttp(iw_pos,jpm)
        do igb1 = 1, npr
          it  = itw(ittp,iw_pos,jpm); itp = itpw(ittp,iw_pos,jpm)
          zw(ittp,igb1) = cmplx(zmel(igb1,it,itp),kind=kp)
        enddo
      enddo
      !$acc end kernels
      !$acc kernels loop independent collapse(2)
      do igb2 = 1, npr
        do ittp = 1, nttp(iw_pos,jpm)
          wzw(ittp,igb2) = cmplx(zw(ittp,igb2)*hilbert_w(ittp,iw_pos,jpm),kind=kp)
        enddo
      enddo
      !$acc end kernels
      ierr = gemm(zw, wzw, rcxq(1,1,iw), npr, npr, nttp(iw_pos,jpm), &
              &  opA = m_op_C, beta = CONE, ldA = nttp_max, ldB = nttp_max)
    enddo
    !$acc end data

    deallocate(itw, itpw, hilbert_w, wzw, zw, nttp)
  end subroutine accumulate_chi0

  subroutine accumulate_chi0_group(ng, gk, gns1, zg, iw_lo, iw_hi, npr)
    !> accumulate_chi0 for ng k points at once (the zmel of k point gk(i) is zg(:,:,:,i), its middle states from gns1(i)):
    !> one product per histogram bin with the pairs of all ng k points, instead of one per bin and k point.  A bin has
    !> a few tens of pairs per k point, so the per-k products are small and bound by launches.  Only for full zmel
    !> (no NMBATCH split), not tetwtk, and npm = 1: the pairs go to the positive bins only (jpmc is not looked at).
    !> npm = 2 (time reversal broken, x0kf_zxq stops for it now) must go through accumulate_chi0.
    !> (2026-09-27 17:36; the npm = 1 guard at the call 2026-09-27 22:00)
    use m_blas, only: m_op_c
    use m_GWinput, only: chi0_filterw_drude
    implicit none
    integer, intent(in) :: ng, gk(ng), gns1(ng), iw_lo, iw_hi, npr
    complex(kind=kp) :: zg(:,:,:,:)
#ifdef __GPU
    attributes(device) :: zg
#endif
    integer :: i, k, icoun, iw, it, itp, ittp, nttp_max, igb, ierr, pos_lo, pos_hi, ntw
    integer, allocatable :: nttp(:), gsl(:,:), git(:,:), gitp(:,:)
    real(8), allocatable :: gw(:,:)
    complex(kind=kp), allocatable :: zw(:,:), wzw(:,:)
    complex(kind=kp), parameter :: CONE = (1_kp, 0_kp)
#ifdef __GPU
    attributes(device) :: zw, wzw
#endif
    pos_lo = max(iw_lo, 1); pos_hi = min(iw_hi, nwhis)
    allocate(nttp(nwhis), source = 0)
    do i = 1, ng
      k = gk(i)
      do icoun = icounkmin(k), icounkmax(k)
        do iw = max(iwini(icoun), pos_lo), min(iwend(icoun), pos_hi)
          nttp(iw) = nttp(iw) + 1
        enddo
      enddo
    enddo
    nttp_max = maxval(nttp)
    if (nttp_max < 1) then
      deallocate(nttp)
      return
    endif
    allocate(gsl(nttp_max,nwhis), git(nttp_max,nwhis), gitp(nttp_max,nwhis), source = 1)
    allocate(gw(nttp_max,nwhis), source = 0d0)
    nttp = 0
    do i = 1, ng
      k = gk(i)
      do icoun = icounkmin(k), icounkmax(k)
        it  = itc(icoun)
        itp = itpc(icoun)
        do iw = max(iwini(icoun), pos_lo), min(iwend(icoun), pos_hi)
          nttp(iw) = nttp(iw) + 1
          ittp = nttp(iw)
          gsl(ittp,iw)  = i
          git(ittp,iw)  = it - gns1(i) + 1
          gitp(ittp,iw) = itp
          gw(ittp,iw)   = whwc(iw-iwini(icoun)+icouini(icoun))
          if (.not.(intrac(icoun) .and. chi0_filterw_drude)) gw(ittp,iw) = gw(ittp,iw)*fcw(iw)
        enddo
      enddo
    enddo
    allocate(zw(nttp_max,npr), wzw(nttp_max,npr))
    !$acc data copyin(gsl, git, gitp, gw)
    do iw = pos_lo, pos_hi
      ntw = nttp(iw)
      if (ntw < 1) cycle
      !$acc parallel loop collapse(2) present(gsl, git, gitp, gw)
      do igb = 1, npr
        do ittp = 1, ntw
          zw(ittp,igb)  = zg(igb, git(ittp,iw), gitp(ittp,iw), gsl(ittp,iw))
          wzw(ittp,igb) = cmplx(zw(ittp,igb)*gw(ittp,iw), kind=kp)
        enddo
      enddo
      ierr = gemm(zw, wzw, rcxq(1,1,iw), npr, npr, ntw, opA = m_op_C, beta = CONE, ldA = nttp_max, ldB = nttp_max)
    enddo
    !$acc end data
    deallocate(zw, wzw, nttp, gsl, git, gitp, gw)
  end subroutine accumulate_chi0_group


  subroutine x0kf_tetwt_write(q, iq)
    !> hgw --tetwt_write: the tetrahedron weights of q for all k, as x0kf_zxq uses them, to __TETWT.<iq>.<isp>.
    use m_readeigen,only: readeval
    real(8), intent(in) :: q(3)
    integer, intent(in) :: iq
    real(8) :: ekxx1(nband,nqbz), ekxx2(nband,nqbz)
    integer :: isp_k, kx, ierr
    do isp_k = 1, nsp
      do kx = 1, nqbz
        ekxx1(1:nband,kx) = readeval(  rk(:,kx), isp_k)
        ekxx2(1:nband,kx) = readeval(q+rk(:,kx), isp_k)
      enddo
      call gettetwt(q, iq, isp_k, isp_k, ekxx1, ekxx2, nband=nband, ikbz_in=1, fkbz_in=nqbz)
      ierr = x0kf_v4hz_init(0, q, isp_k, isp_k, iq, .false., ikbz_in=1, fkbz_in=nqbz)
      ierr = x0kf_v4hz_init(1, q, isp_k, isp_k, iq, .false., ikbz_in=1, fkbz_in=nqbz)
      call tetdeallocate()
      call tetwt_save(iq, isp_k, q, ekxx1, ekxx2)
    enddo
  end subroutine x0kf_tetwt_write

  character(64) function tetwt_file(iq, isp)
    integer, intent(in) :: iq, isp
    write(tetwt_file,'("__TETWT.",i0,".",i0)') iq, isp
  end function tetwt_file

  real(8) function tetwt_ekey(ekxx1, ekxx2) result(key)
    !> Checksum of the band energies at k and k+q (weights differ by band and k, so a swap also shows).
    real(8), intent(in) :: ekxx1(:,:), ekxx2(:,:)
    integer :: ib, kx
    key = 0d0
    do kx = 1, size(ekxx1, 2)
      do ib = 1, size(ekxx1, 1)
        key = key + ekxx1(ib,kx)*(1d0 + 1d-3*mod(7*ib + 13*kx, 101)) + ekxx2(ib,kx)*(1d0 + 1d-3*mod(11*ib + 17*kx, 103))
      enddo
    enddo
  end function tetwt_ekey

  real(8) function tetwt_pkey() result(key)
    !> Checksum of the rest the weights depend on: the histogram bins, E_F (and finite T, its shift), the band cut.
    use m_ReadEfermi, only: ef
    use m_readgwinput, only: ebmx, nbmx, mtet
    use m_GWinput, only: t_tetrakbt
    use m_cmdopt_registry, only: c2_EfermiShifteV, c2_EfermiShifteV_set
    integer :: i
    key = 0d0
    do i = 1, nwhis + 1
      key = key + frhis(i)*(1d0 + 1d-3*mod(i, 97))
    enddo
    key = key + 3.1d0*ef + 5.3d0*t_tetrakbt + 7.1d0*ebmx + 1.3d0*nbmx + 0.7d0*sum(mtet)
    if (c2_EfermiShifteV_set) key = key + 11.3d0*c2_EfermiShifteV
  end function tetwt_pkey

  subroutine tetwt_save(iq, isp, q, ekxx1, ekxx2)
    !> The arrays of x0kf_v4hz_init (all k) to __TETWT.<iq>.<isp>, through a temporary name and rename.
    integer, intent(in) :: iq, isp
    real(8), intent(in) :: q(3), ekxx1(:,:), ekxx2(:,:)
    character(64) :: fn
    integer :: ifi
    fn = tetwt_file(iq, isp)
    open(newunit=ifi, file=trim(fn)//'.tmp', access='stream', form='unformatted', status='replace')
    write(ifi) tetwt_tag, nqbz, nband, nctot, npm, nwhis, q, tetwt_ekey(ekxx1, ekxx2), tetwt_pkey(), ncount, ncoun
    write(ifi) nkmin, nkmax, nkqmin, nkqmax, icounkmin, icounkmax
    write(ifi) kc, iwini, iwend, itc, itpc, jpmc, icouini, intrac, whwc
    close(ifi)
    if (c_rename(trim(fn)//'.tmp'//c_null_char, trim(fn)//c_null_char) /= 0) call rx('tetwt_save: rename failed '//trim(fn))
    if (ipr) write(stdo,'(1x,a,2i8)') 'tetwt: wrote '//trim(fn)//'  ncount ncoun =', ncount, ncoun
  end subroutine tetwt_save

  logical function tetwt_load(iq, isp, q, ekxx1, ekxx2) result(ok)
    !> Read __TETWT.<iq>.<isp> into the arrays of x0kf_v4hz_init (all k) when it matches this q and these band
    !> energies; .false. (and nothing changed that x0kf_v4hz_init would not reset) otherwise.
    integer, intent(in) :: iq, isp
    real(8), intent(in) :: q(3), ekxx1(:,:), ekxx2(:,:)
    character(64) :: fn
    integer :: ifi, ios, tag, n1, n2, n3, n4, n5, nc, nco
    real(8) :: qf(3), ek, fk
    logical :: ex
    ok = .false.
    fn = tetwt_file(iq, isp)
    inquire(file=trim(fn), exist=ex)
    if (.not.ex) return
    open(newunit=ifi, file=trim(fn), access='stream', form='unformatted', status='old', action='read', iostat=ios)
    if (ios /= 0) return
    read(ifi, iostat=ios) tag, n1, n2, n3, n4, n5, qf, ek, fk, nc, nco
    if (ios /= 0 .or. tag /= tetwt_tag .or. n1 /= nqbz .or. n2 /= nband .or. n3 /= nctot .or. n4 /= npm .or. &
        n5 /= nwhis .or. any(qf /= q) .or. ek /= tetwt_ekey(ekxx1, ekxx2) .or. fk /= tetwt_pkey()) then
      close(ifi)
      if (ipr) write(stdo,'(1x,a)') 'tetwt: '//trim(fn)//' does not match this q or these bands; computing the weights'
      return
    endif
    ncount = nc
    ncoun  = nco
    if (allocated(nkmin))     deallocate(nkmin)
    if (allocated(nkmax))     deallocate(nkmax)
    if (allocated(nkqmin))    deallocate(nkqmin)
    if (allocated(nkqmax))    deallocate(nkqmax)
    if (allocated(icounkmin)) deallocate(icounkmin)
    if (allocated(icounkmax)) deallocate(icounkmax)
    if (allocated(whwc))      deallocate(whwc)
    if (allocated(kc))        deallocate(kc)
    if (allocated(iwini))     deallocate(iwini)
    if (allocated(iwend))     deallocate(iwend)
    if (allocated(itc))       deallocate(itc)
    if (allocated(itpc))      deallocate(itpc)
    if (allocated(jpmc))      deallocate(jpmc)
    if (allocated(icouini))   deallocate(icouini)
    if (allocated(intrac))    deallocate(intrac)
    allocate(nkmin(nqbz), nkmax(nqbz), nkqmin(nqbz), nkqmax(nqbz), icounkmin(nqbz), icounkmax(nqbz))
    allocate(whwc(ncount), kc(ncoun), iwini(ncoun), iwend(ncoun), itc(ncoun), itpc(ncoun), jpmc(ncoun), &
             icouini(ncoun), intrac(ncoun))
    read(ifi, iostat=ios) nkmin, nkmax, nkqmin, nkqmax, icounkmin, icounkmax
    if (ios == 0) read(ifi, iostat=ios) kc, iwini, iwend, itc, itpc, jpmc, icouini, intrac, whwc
    close(ifi)
    ok = ios == 0
    if (ipr) write(stdo,'(1x,a,l2)') 'tetwt: weights read from '//trim(fn)//':', ok
  end function tetwt_load
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
