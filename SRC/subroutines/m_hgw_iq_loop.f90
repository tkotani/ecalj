!> Per-iq production loop: compute chi0(iq) via x0kf_zxq, then derive W-v(iq)
!> via WVRllwR/WVIllwI on root_k. Encapsulates the iq loop that lives at the
!> top level of hrcxq so streaming callers can either call this verbatim
!> (FILE-mode = current behavior, W goes to __WVR.<iq> / __WVI.<iq>) or, in a
!> later step, interleave per-iq sxcf consumption (Phase 1-C streaming).
module m_hgw_iq_loop
contains
  subroutine run_iq_loop(iqxini, iqxend, mpi__Qtask, realomega, imagomega)
    use m_qbze,      only: qibze
    use m_x0kf,      only: x0kf_zxq, deallocatezxq, deallocatezxqi
    use m_llw,       only: WVRllwR, WVIllwI
    use m_readVcoud, only: Readvcoud, ngb
    use m_mpi,       only: mpi__root_k, comm_k, comm_b, mpi__rank, MPI__Setnpr_col, ipr
    use m_lgunit,    only: stdo
    use m_ftox
    integer, intent(in) :: iqxini, iqxend
    logical, intent(in) :: mpi__Qtask(iqxini:iqxend), realomega, imagomega
    integer :: iq, npr, npr_col, ierr
    real(8) :: qp(3)
    real(8), parameter :: schi = -9999d0
    do iq = iqxini, iqxend
      if (.not. mpi__Qtask(iq)) cycle
      if (ipr) write(stdo,*) 'mpi_rank in IQ loop:', iq, mpi__rank
      call cputid(0)
      qp = qibze(:,iq)
      if (ipr) write(stdo,ftox) 'do 1001: iq q=', iq, ftof(qp,4), ' of nq=', iqxend
      call Readvcoud(qp, iq, NoVcou=.false.)  ! Readin vcousq, zcousq, ngb, ngc for the Coulomb matrix
      npr = ngb
      call MPI__Setnpr_col(npr, npr_col)      ! split of npr (column) for MPI color_b
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
    enddo
  end subroutine run_iq_loop
end module m_hgw_iq_loop
