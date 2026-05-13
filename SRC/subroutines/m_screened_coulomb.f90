!> Screened-Coulomb building primitive: compute W-v(iq) for a single q-point.
module m_screened_coulomb
contains

  !> Compute W-v for a single iq: Readvcoud -> x0kf_zxq -> WVRllwR/WVIllwI.
  !> Analog of sxcf_correlation_step_kx for the screened-Coulomb building phase.
  subroutine build_screened_coulomb_step_kx(iq, qp, realomega, imagomega)
    use m_x0kf,        only: x0kf_zxq, deallocatezxq, deallocatezxqi
    use m_llw,         only: WVRllwR, WVIllwI
    use m_readVcoud,   only: Readvcoud, ngb
    use m_mpi,         only: mpi__root_k, mpi__rank, comm_k, comm_b, comm_q, MPI__Setnpr_col, ipr
    use m_wv_storage,  only: wv_backend, WV_BACKEND_MEMORY_3D, WV_BACKEND_SHM, &
                             wv_init_shm, wv_alloc_recv_bufs
    use m_freq,        only: npm, nwhis, wv_niw => niw, nw_i, nw
    use m_lgunit,      only: stdo
    use m_ftox
    use mpi
    integer, intent(in) :: iq
    real(8), intent(in) :: qp(3)
    logical, intent(in) :: realomega, imagomega
    integer :: npr, npr_col, ierr
    real(8), parameter :: schi = -9999d0
    if (ipr) write(stdo,*) 'mpi_rank in IQ loop:', iq, mpi__rank
    call cputid(0)
    if (ipr) write(stdo,ftox) 'do 1001: iq q=', iq, ftof(qp,4)
    call Readvcoud(qp, iq, NoVcou=.false.)
    if (wv_backend == WV_BACKEND_SHM) &
      call wv_init_shm(ngb, wv_niw, (1-npm)*nwhis, nwhis, comm_q)
    if (mpi__root_k .and. wv_backend == WV_BACKEND_MEMORY_3D) &
      call wv_alloc_recv_bufs(ngb, wv_niw, nw_i, nw)
    npr = ngb
    call MPI__Setnpr_col(npr, npr_col)
    call x0kf_zxq(realomega, imagomega, qp, iq, npr, schi, &
                  crpa=.false., chipm=.false., nolfco=.false., is_m_basis=.true.)
    if (mpi__root_k) then
      call WVRllwR(qp, iq, npr, npr_col, is_x0_m_basis=.true., is_wc_m_basis=.true.)
      if (wv_backend /= WV_BACKEND_MEMORY_3D) call deallocatezxq()
      call WVIllwI(qp, iq, npr, npr_col, is_x0_m_basis=.true., is_wc_m_basis=.true.)
      if (wv_backend /= WV_BACKEND_MEMORY_3D) call deallocatezxqi()
    endif
    call mpi_barrier(comm_k, ierr)
    call mpi_barrier(comm_b, ierr)
  end subroutine build_screened_coulomb_step_kx

end module m_screened_coulomb
