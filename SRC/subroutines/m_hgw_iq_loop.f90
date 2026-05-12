!> Screened-Coulomb building primitive: compute W-v(iq) for a single q-point.
module m_hgw_iq_loop
contains

  !> Compute W-v for a single iq: Readvcoud -> x0kf_zxq -> WVRllwR/WVIllwI.
  !> Analog of sxcf_correlation_step_kx for the screened-Coulomb building phase.
  subroutine build_screened_coulomb_step_kx(iq, qp, realomega, imagomega)
    use m_x0kf,        only: x0kf_zxq, deallocatezxq, deallocatezxqi
    use m_llw,         only: WVRllwR, WVIllwI
    use m_readVcoud,   only: Readvcoud, ngb
    use m_mpi,         only: mpi__root_k, mpi__rank, comm_k, comm_b, MPI__Setnpr_col, ipr
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
  end subroutine build_screened_coulomb_step_kx

end module m_hgw_iq_loop
