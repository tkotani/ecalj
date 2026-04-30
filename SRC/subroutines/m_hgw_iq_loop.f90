!> Per-iq production loop: compute chi0(iq) via x0kf_zxq, then derive W-v(iq)
!> via WVRllwR/WVIllwI on root_k. When streaming_consume=.true., interleave
!> per-iq sxcf consumption (Phase 1-C streaming) — for each iq <= nqibz, after
!> WV is in the m_wv_storage MEMORY_3D current buffer and synced across all
!> ranks, either save it as the iq=1 slot (iq==1, deferred until W0w0i runs)
!> or call sxcf_correlation_step_kx(kx=iq) right away (iq>1).
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
    use m_wv_storage,  only: wv_zero_current, wv_sync_current, wv_save_iq1
    use m_sxcf_sc,     only: sxcf_correlation_step_kx
    use m_hsfp0_sc,    only: hs_ef, hs_esmr, hs_nspinmx
    integer, intent(in) :: iqxini, iqxend
    logical, intent(in) :: mpi__Qtask(iqxini:iqxend), realomega, imagomega
    logical, intent(in), optional :: streaming_consume
    logical :: stream
    integer :: iq, npr, npr_col, ierr
    real(8) :: qp(3)
    real(8), parameter :: schi = -9999d0
    stream = .false.
    if (present(streaming_consume)) stream = streaming_consume
    do iq = iqxini, iqxend
      qp = qibze(:,iq)
      ! In streaming mode, every rank zeros the current W buffer for iq <= nqibz
      ! so the subsequent Allreduce SUM yields just the owner's contribution.
      if (stream .and. iq <= nqibz) call wv_zero_current()
      OwnerWork: if (mpi__Qtask(iq)) then
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
      endif OwnerWork
      ! Streaming consume: distribute W(iq) to all ranks, then save (iq==1) or
      ! consume (iq>1, run sxcf step_kx). Skip for iq > nqibz (those iqs only
      ! contribute to llw / wmuk and have no W to consume).
      if (stream .and. iq <= nqibz) then
        call wv_sync_current(comm)
        if (iq == 1) then
          call wv_save_iq1()  ! defer iq=1 consume until after W0w0i correction
        else
          call sxcf_correlation_step_kx(iq, hs_ef, hs_esmr, hs_nspinmx)
        endif
      endif
    enddo
  end subroutine run_iq_loop
end module m_hgw_iq_loop
