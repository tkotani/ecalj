!> GW main driver: build Im(chi0), Hilbert transform → W; stream W in-process
!> for exchange and correlation self-energy via hsfp0_sc kernels.
module m_hgw
  contains
subroutine hgw(do_correlation, do_exchange)
  !> Streaming flow (SHM backend):
  !>   Phase 1: iq=2..nqibz — build W per-iq, consume for Sc immediately.
  !>   Phase 2: auxiliary q0 points (iq > nqibz) — build W/llw only.
  !>   Phase 3: iq=1/Gamma — build W, apply W0w0i, consume for Sc.
  !> Exchange (do_exchange=.true.) runs before SplitXq using the global
  !> k-distribution, then hsfp0_sc_sc_writeout emits SECU/SEC2U.
  use m_ReadEfermi,only: Readefermi
  use m_readqg,only: Readngmx2
  use m_hamindex,only: Readhamindex, symgg=>symops, ngrp
  use m_readeigen,only: Init_readeigen,Init_readeigen2
  use m_read_bzdata,only:Read_bzdata, nq0i,nq0iadd,nqibz,q0i
  use m_genallcf_v3,only: Genallcf_v3
  use m_rdpp,only: mrecl,nblochpmx,nprecx
  use m_zmel,only: Mptauof_zmel
  use m_itq,only:  Setitq
  use m_freq,only: Getfreq2,freq_r,nw_i,nw,niw
  use m_w0w0i,only: W0w0i
  use m_readgwinput,only: ReadGwinputKeys
  use m_qbze,only:  Setqbze,qibze
  use m_llw,only: MPI__sendllw_q
  use m_mpi,only: MPI__Initialize,MPI__root,MPI__rank,MPI__size,MPI__consoleout,comm, &
                & MPI__SplitXq, MPI__FreeSplitXq, &
                & mpi__root_k, ipr, comm_q, mpi__rank_q
  use m_lgunit,only: m_lgunit_init,stdo
  use m_ftox
  use m_gpu,only: gpu_init
  use m_hsfp0_sc,only: hsfp0_sc, hsfp0_sc_setup, hsfp0_sc_writeout, &
                       hs_ef, hs_esmr, hs_nspinmx
  use m_screened_coulomb,only: build_screened_coulomb_step_kx
  use m_x0kf,only: deallocatezxq, deallocatezxqi
  use m_wv_storage,only: wv_dealloc
  use m_sxcf_sc,only: sxcf_correlation_init, sxcf_correlation_step_kx, &
                      sxcf_correlation_finalize
  use m_sxcf_count,only: sxcf_scz_count, mpi_assign_qtask_lpt
  use mpi
  implicit none
  logical, intent(in) :: do_correlation, do_exchange
  integer :: iq, iqxend, iw, ifwd, verbose, ifif, ierr
  real(8) :: ua=1d0, qp(3)
  logical :: debug=.false., realomega, imagomega
  logical :: hx0, iprintx=.false.
  logical, allocatable :: mpi__Qtask(:)
  integer :: worker_inQtask
  call MPI__Initialize()
  call gpu_init(comm)
  call M_lgunit_init()
  call MPI__consoleout('hgw')
  call cputid (0)
  if(verbose()>=100) debug= .TRUE.
  call Genallcf_v3(incwfx=-1)
  call Read_BZDATA(hx0)
  call Readefermi()
  call ReadGWinputKeys()
  call Readngmx2()
  call Setqbze()
  ! (symgg, ngrp) is a superset of identity-only setup; hgw only uses symop 1.
  call Readhamindex()
  call Mptauof_zmel(symgg, ngrp)
  call Setitq()
  call Init_readeigen()
  call Init_readeigen2()
  realomega = .true.
  imagomega = .true.
  call Getfreq2(.false.,realomega,imagomega,ua,iprintx)
  if(MPI__root) call writewvfreq()
  iqxend = nqibz + nq0i + nq0iadd
  if(ipr) write(stdo,'(1X,A,I5)') 'MPI: worker_inQtask (omega-parallel):', mpi__size

  ! Exchange runs before any SplitXq so mpi__size_k=mpi__size (global k-distribution).
  if (do_exchange) then
     if(ipr) write(stdo,ftox) ' hgw: starting in-process hsfp0_sc(--job=1) exchange phase'
     call hsfp0_sc(skip_init=.true., skip_rx0=.true., ixc_in=1)
  endif

  if(sum(qibze(:,1)**2)>1d-10) call rx(' hgw: sanity check. |q(iq=1)| /= 0')

  ! Phase 1: iq=2..nqibz — build W-v and consume immediately per-iq.
  ! SplitXq before hsfp0_sc_setup so rankdivider uses q-group-local mpi__rank_k/mpi__size_k.
  worker_inQtask = mpi__size
  call MPI__SplitXq(1, worker_inQtask)
  call hsfp0_sc_setup(skip_init=.true., ixc_in=2)
  call sxcf_correlation_init(hs_ef, hs_esmr, hs_nspinmx)
  allocate( mpi__Qtask(2:nqibz) )
  call mpi_assign_qtask_lpt(2, nqibz, hs_nspinmx, mpi__size/worker_inQtask, &
                             mpi__rank/worker_inQtask, worker_inQtask, mpi__Qtask)
  if(ipr) write(stdo,ftox) 'Phase1: mpi_rank',mpi__rank,'worker_inQtask',worker_inQtask,'mpi__Qtask=',mpi__Qtask
  call flush(stdo)
  do iq = 2, nqibz
    if (.not. mpi__Qtask(iq)) cycle
    qp = qibze(:,iq)
    call build_screened_coulomb_step_kx(iq, qp, realomega, imagomega)
    ! SHM data already visible via shared memory; barrier ensures all ranks
    ! see W before sxcf reads, and also synchronises before wv_init_shm (iq+1).
    call MPI_barrier(comm_q, ierr)
    call sxcf_correlation_step_kx(iq, hs_ef, hs_esmr, hs_nspinmx)
    ! Ensure all ranks finish sxcf(iq) before any rank enters wv_init_shm for iq+1
    ! (MPI_Win_free is collective on comm_q).
    call MPI_barrier(comm_q, ierr)
  enddo
  call MPI_barrier(comm, ierr)
  deallocate(mpi__Qtask)
  call MPI__FreeSplitXq()

  ! Phase 2: auxiliary q0 points (iq > nqibz) — W/llw only, no consumption.
  worker_inQtask = mpi__size
  call MPI__SplitXq(1, worker_inQtask)
  allocate( mpi__Qtask(nqibz+1:iqxend), source=[(mod(iq-nqibz-1,mpi__size/worker_inQtask)==mpi__rank/worker_inQtask,iq=nqibz+1,iqxend)])
  if(ipr) write(stdo,ftox) 'Phase2: mpi_rank',mpi__rank,'worker_inQtask',worker_inQtask,'mpi__Qtask=',mpi__Qtask
  call flush(stdo)
  if (nqibz < iqxend) then
    ! Pass 1: each q-group builds its assigned iq points.
    do iq = nqibz+1, iqxend
      if (.not. mpi__Qtask(iq)) cycle
      qp = qibze(:,iq)
      call build_screened_coulomb_step_kx(iq, qp, realomega, imagomega)
    enddo
    call MPI_barrier(comm, ierr)
    ! Pass 2: gather llw to rank 0 (blocking send/recv safe after barrier).
    ! src = first rank of the owning q-group (valid for n_bpara=1).
    do iq = nqibz+1, iqxend
      call MPI__sendllw_q(iq-nqibz, mod(iq-nqibz-1,mpi__size/worker_inQtask)*worker_inQtask, 0)
    enddo
  endif
  deallocate(mpi__Qtask)
  call MPI__FreeSplitXq()

  ! Phase 3: iq=1/Gamma — all ranks collaborate; build W-v, apply W0w0i, consume.
  worker_inQtask = mpi__size
  call MPI__SplitXq(1, worker_inQtask)
  if(ipr) write(stdo,ftox) 'Phase3: mpi_rank',mpi__rank
  call flush(stdo)
  call sxcf_scz_count(hs_ef, hs_esmr, .false., 2, hs_nspinmx)
  qp = qibze(:,1)
  call build_screened_coulomb_step_kx(1, qp, realomega, imagomega)
  call MPI_barrier(comm, ierr)
  if (MPI__rank == 0) call W0w0i(nw_i, nw, nq0i, niw, q0i, is_wc_m_basis=.true.)
  ! Barrier ensures W0w0i correction in shm_wvr is visible to all ranks before sxcf.
  call MPI_barrier(comm_q, ierr)
  call sxcf_correlation_step_kx(1, hs_ef, hs_esmr, hs_nspinmx)
  call sxcf_correlation_finalize()
  call wv_dealloc()
  call hsfp0_sc_writeout(skip_rx0=.true.)
  call MPI__FreeSplitXq()

  if(ipr) write(stdo,ftox) '--- end of hgw --- irank=',MPI__rank
  call cputid(0)
  call rx0( ' OK! hgw finished')

  contains
  subroutine writewvfreq() !writeonly
     open(newunit=ifwd, file='__WV.d')
     write(ifwd,"(1x,10i14)") nprecx, mrecl, nblochpmx, nw+1,niw, nqibz + nq0i-1, nw_i
     close(ifwd)
     open(newunit=ifif,file='freq_r')
     write(ifif,"(2i8,'  !(a.u.=2Ry)')") nw+1, nw_i
     do iw= nw_i,-1
        write(ifif,"(d23.15,2x,i6)") -freq_r(-iw),iw
     enddo
     do iw= 0,nw
        write(ifif,"(d23.15,2x,i6)") freq_r(iw),iw
     enddo
     close(ifif)
  end subroutine writewvfreq
end subroutine hgw
end module m_hgw
