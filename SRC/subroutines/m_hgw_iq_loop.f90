!> Per-iq production loop: compute chi0(iq) via x0kf_zxq, then derive W-v(iq)
!> via WVRllwR/WVIllwI on root_k. When streaming_consume=.true., interleave
!> per-iq sxcf consumption (Phase 1-C streaming) — for each iq <= nqibz, after
!> WV is in the m_wv_storage MEMORY_3D current buffer and synced across all
!> ranks, call sxcf_correlation_step_kx(kx=iq) immediately. iq=1/Gamma is
!> processed last so W0w0i can operate on the current buffer directly after
!> the loop returns (no separate iq=1 saved slot needed).
!>
!> MEMORY_3D buffer setup per-iq:
!>   Qtask ranks: x0kf_zxq calls wv_assoc_real_buf(rcxq) + wv_init_imag_buf.
!>   Non-Qtask ranks: ngb is broadcast via MPI_Allreduce(MAX) so they can call
!>   wv_alloc_zero_bufs(ngb, ...) before wv_sync_current.
module m_hgw_iq_loop
contains
  subroutine run_iq_loop(iqxini, iqxend, mpi__Qtask, realomega, imagomega, &
                          streaming_consume)
    use m_qbze,        only: qibze
    use m_x0kf,        only: x0kf_zxq, deallocatezxq, deallocatezxqi
    use m_llw,         only: WVRllwR, WVIllwI
    use m_readVcoud,   only: Readvcoud, ngb
    use m_read_bzdata, only: nqibz
    use m_mpi,         only: mpi__root_k, comm_k, comm_b, comm, mpi__rank, &
                              MPI__Setnpr_col, ipr
    use m_lgunit,      only: stdo
    use m_ftox
    use m_wv_storage,  only: wv_sync_current, wv_alloc_zero_bufs, &
                              wv_backend, WV_BACKEND_MEMORY_3D
    use m_sxcf_sc,     only: sxcf_correlation_step_kx
    use m_hsfp0_sc,    only: hs_ef, hs_esmr, hs_nspinmx
    use m_freq,        only: nw_i, nw, niw
    use mpi
    integer, intent(in) :: iqxini, iqxend
    logical, intent(in) :: mpi__Qtask(iqxini:iqxend)
    logical, intent(in) :: realomega, imagomega
    logical, intent(in), optional :: streaming_consume
    logical :: stream
    integer :: iq, npr, npr_col, ierr, ngb_cur, ierr2
    real(8) :: qp(3)
    real(8), parameter :: schi = -9999d0
    stream = .false.
    if (present(streaming_consume)) stream = streaming_consume

    ! iq=2,...,iqxend first (including auxiliary q0 points beyond nqibz).
    do iq = max(iqxini, 2), iqxend
      qp = qibze(:,iq)
      if (mpi__Qtask(iq)) then
        if (ipr) write(stdo,*) 'mpi_rank in IQ loop:', iq, mpi__rank
        call cputid(0)
        if (ipr) write(stdo,ftox) 'do 1001: iq q=', iq, ftof(qp,4), ' of nq=', iqxend
        call Readvcoud(qp, iq, NoVcou=.false.)
        npr = ngb
        call MPI__Setnpr_col(npr, npr_col)
        call x0kf_zxq(realomega, imagomega, qp, iq, npr, schi, &
                      crpa=.false., chipm=.false., nolfco=.false., is_m_basis=.true.)
        if (mpi__root_k) then
          call WVRllwR(qp, iq, npr, npr_col, is_x0_m_basis=.true., is_wc_m_basis=.true.)
          call deallocatezxq()
          call WVIllwI(qp, iq, npr, npr_col, is_x0_m_basis=.true., is_wc_m_basis=.true.)
          call deallocatezxqi()
        endif
        call mpi_barrier(comm_k, ierr)
        call mpi_barrier(comm_b, ierr)
      endif
      if (stream .and. iq <= nqibz) then
        ! Only root_k within the Qtask group wrote W-v; all other ranks need
        ! zero-filled buffers for AllreduceSum (non-root_k Qtask ranks hold
        ! partial chi0 in wv_real_buf => rcxq, which must not contribute).
        ngb_cur = merge(ngb, 0, mpi__Qtask(iq))
        call MPI_Allreduce(MPI_IN_PLACE, ngb_cur, 1, MPI_INTEGER, MPI_MAX, comm, ierr2)
        if (.not. (mpi__Qtask(iq) .and. mpi__root_k) .and. wv_backend == WV_BACKEND_MEMORY_3D) &
          call wv_alloc_zero_bufs(ngb_cur, niw, nw_i, nw)
        call wv_sync_current(comm)
        call sxcf_correlation_step_kx(iq, hs_ef, hs_esmr, hs_nspinmx)
      endif
    enddo

    ! iq=1 (Gamma) last: after this, current holds WV(iq=1) ready for
    ! W0w0i correction in the caller (no save/restore needed).
    if (iqxini <= 1 .and. 1 <= iqxend) then
      iq = 1; qp = qibze(:,1)
      if (mpi__Qtask(1)) then
        if (ipr) write(stdo,*) 'mpi_rank in IQ loop:', 1, mpi__rank
        call cputid(0)
        if (ipr) write(stdo,ftox) 'do 1001: iq q=', 1, ftof(qp,4), ' of nq=', iqxend
        call Readvcoud(qp, 1, NoVcou=.false.)
        npr = ngb
        call MPI__Setnpr_col(npr, npr_col)
        call x0kf_zxq(realomega, imagomega, qp, 1, npr, schi, &
                      crpa=.false., chipm=.false., nolfco=.false., is_m_basis=.true.)
        if (mpi__root_k) then
          call WVRllwR(qp, 1, npr, npr_col, is_x0_m_basis=.true., is_wc_m_basis=.true.)
          call deallocatezxq()
          call WVIllwI(qp, 1, npr, npr_col, is_x0_m_basis=.true., is_wc_m_basis=.true.)
          call deallocatezxqi()
        endif
        call mpi_barrier(comm_k, ierr)
        call mpi_barrier(comm_b, ierr)
      endif
      ! For iq=1 in Phase 3 all ranks are Qtask; only root_k wrote W-v via
      ! WVRllwR. Non-root_k ranks must contribute zero to AllreduceSum.
      if (stream) then
        ngb_cur = merge(ngb, 0, mpi__Qtask(1))
        call MPI_Allreduce(MPI_IN_PLACE, ngb_cur, 1, MPI_INTEGER, MPI_MAX, comm, ierr2)
        if (.not. (mpi__Qtask(1) .and. mpi__root_k) .and. wv_backend == WV_BACKEND_MEMORY_3D) &
          call wv_alloc_zero_bufs(ngb_cur, niw, nw_i, nw)
        call wv_sync_current(comm)
      endif
    endif

  end subroutine run_iq_loop
end module m_hgw_iq_loop
