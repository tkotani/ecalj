!> Calculate Im(chi0), do Hilbert transformation, and (in streaming mode)
!> consume the resulting W in-process for correlation self-energy.
module m_hrcxq
  contains
subroutine hrcxq(do_correlation, do_exchange)
  !> When do_correlation=.true. we run Phase 1-C streaming: per-iq W is held
  !> only in the m_wv_storage MEMORY_3D singleton, sxcf step_kx consumes it
  !> immediately (iq>1); iq=1/Gamma is processed last and consumed after W0w0i
  !> correction. SECU/SEC2U are written via hsfp0_sc_writeout. No __WVR/__WVI.
  !> When do_correlation is absent or .false., we run the legacy FILE flow:
  !> WVRllwR/WVIllwI emit __WVR.<iq>/__WVI.<iq>, W0w0i edits __WVR.1/__WVI.1
  !> in place, then (optionally) hsfp0_sc(--job=1) writes SEXU/SEX2U.
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
  use m_llw,only: MPI__sendllw, MPI__sendllw_q
  use m_mpi,only: MPI__Initialize,MPI__root,MPI__rank,MPI__size,MPI__consoleout,comm, &
                & MPI__SplitXq, MPI__FreeSplitXq, mpi__root_k, ipr, comm_q, mpi__rank_q
  use m_lgunit,only: m_lgunit_init,stdo
  use m_ftox
  use m_gpu,only: gpu_init
  use m_hsfp0_sc,only: hsfp0_sc, hsfp0_sc_setup, hsfp0_sc_writeout, &
                       hs_ef, hs_esmr, hs_nspinmx
  use m_hgw_iq_loop,only: build_screened_coulomb_step_kx
  use m_readVcoud,only: ngb
  use m_wv_storage,only: wv_init_file, wv_init_memory_3d, wv_dealloc, &
                         wv_bcast_current, wv_sync_current, wv_alloc_zero_bufs, &
                         wv_backend, WV_BACKEND_MEMORY_3D
  use m_sxcf_sc,only: sxcf_correlation_init, sxcf_correlation_step_kx, &
                      sxcf_correlation_finalize
  use m_sxcf_count,only: sxcf_scz_count, mpi_assign_qtask_lpt
  use mpi
  implicit none
  logical, intent(in), optional :: do_correlation, do_exchange
  integer :: iq, iqxini, iqxend, iw, ifwd, verbose, ifif, ierr
  real(8) :: ua=1d0, qp(3)
  logical :: debug=.false., realomega, imagomega
  logical :: hx0, iprintx=.false.
  logical :: cmdopt2
  logical :: streaming
  character(20) :: outs=''
  logical, allocatable :: mpi__Qtask(:)
  integer, allocatable :: mpi__Qrank(:)
  integer :: n_kpara = 1, n_bpara = 1, worker_inQtask
  call MPI__Initialize()
  call gpu_init(comm)
  call M_lgunit_init()
  call MPI__consoleout('hrcxq')
  call cputid (0)
  if(verbose()>=100) debug= .TRUE.
  call Genallcf_v3(incwfx=-1)
  call Read_BZDATA(hx0)
  call Readefermi()
  call ReadGWinputKeys()
  call Readngmx2()
  call Setqbze()
  ! (symgg, ngrp) is a superset of identity-only setup; hrcxq only uses symop 1.
  call Readhamindex()
  call Mptauof_zmel(symgg, ngrp)
  call Setitq()
  call Init_readeigen()
  call Init_readeigen2()
  realomega = .true.
  imagomega = .true.
  call Getfreq2(.false.,realomega,imagomega,ua,iprintx)
  if(MPI__root) call writewvfreq()
  iqxini = 1
  iqxend = nqibz + nq0i + nq0iadd
  n_kpara = mpi__size
  if(ipr) write(stdo,'(1X,A,2I5)') 'MPI: n_bpara, n_kpara:', n_bpara, n_kpara

  ! Exchange runs before any SplitXq so mpi__size_k=mpi__size (global k-distribution).
  if (present(do_exchange)) then
     if (do_exchange) then
        if(ipr) write(stdo,ftox) ' hrcxq: starting in-process hsfp0_sc(--job=1) exchange phase'
        call hsfp0_sc(skip_init=.true., skip_rx0=.true., ixc_in=1)
     endif
  endif

  if(sum(qibze(:,1)**2)>1d-10) call rx(' hx0fp0.sc: sanity check. |q(iqx)| /= 0')

  streaming = .false.
  if (present(do_correlation)) streaming = do_correlation

  StreamingOrFile: if (streaming) then
     ! Phase 1: iq=2..nqibz — build W-v and consume immediately per-iq.
     ! SplitXq before hsfp0_sc_setup so rankdivider uses q-group-local mpi__rank_k/mpi__size_k.
     call MPI__SplitXq(n_bpara, n_kpara)
     call hsfp0_sc_setup(skip_init=.true., ixc_in=2)
     call wv_init_memory_3d(nw_i, nw)
     call sxcf_correlation_init(hs_ef, hs_esmr, hs_nspinmx)
     worker_inQtask = n_bpara * n_kpara
     allocate( mpi__Qtask(2:nqibz), mpi__Qrank(2:nqibz) )
     call mpi_assign_qtask_lpt(2, nqibz, hs_nspinmx, mpi__size/worker_inQtask, &
                                mpi__rank/worker_inQtask, worker_inQtask, mpi__Qtask, mpi__Qrank)
     deallocate(mpi__Qrank)
     if(ipr) write(stdo,ftox) 'Phase1: mpi_rank',mpi__rank,'worker_inQtask',worker_inQtask,'mpi__Qtask=',mpi__Qtask
     call flush(stdo)
     do iq = 2, nqibz
       if (.not. mpi__Qtask(iq)) cycle
       qp = qibze(:,iq)
       call build_screened_coulomb_step_kx(iq, qp, realomega, imagomega)
       if (.not. mpi__root_k .and. wv_backend == WV_BACKEND_MEMORY_3D) call wv_alloc_zero_bufs(ngb, niw, nw_i, nw)
       call wv_sync_current(comm_q)
       call sxcf_correlation_step_kx(iq, hs_ef, hs_esmr, hs_nspinmx)
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
       ! src = global rank of mpi__root_k in the owning q-group = first rank of that group.
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
     if (.not. mpi__root_k .and. wv_backend == WV_BACKEND_MEMORY_3D) call wv_alloc_zero_bufs(ngb, niw, nw_i, nw)
     call wv_sync_current(comm_q)
     call MPI_barrier(comm, ierr)
     if (MPI__rank == 0) call W0w0i(nw_i, nw, nq0i, niw, q0i, is_wc_m_basis=.true.)
     ! Broadcast the corrected current buffer to all ranks, then consume.
     call wv_bcast_current(0, comm)
     call sxcf_correlation_step_kx(1, hs_ef, hs_esmr, hs_nspinmx)
     call sxcf_correlation_finalize()
     call wv_dealloc()
     call hsfp0_sc_writeout(skip_rx0=.true.)
     call MPI__FreeSplitXq()
     if(ipr) write(stdo,ftox) ' hrcxq+hsfp0_sc combined: finished (streaming)'
  else StreamingOrFile
     ! ---- Legacy FILE mode (no in-process correlation consume) ----
     ! Process iq=2,...,iqxend first, then iq=1 last (same ordering as streaming).
     call MPI__SplitXq(n_bpara, n_kpara)
     worker_inQtask = n_bpara * n_kpara
     allocate( mpi__Qrank(iqxini:iqxend), source=[(mod(iq-1,mpi__size/worker_inQtask)*worker_inQtask           ,iq=iqxini,iqxend)])
     allocate( mpi__Qtask(iqxini:iqxend), source=[(mod(iq-1,mpi__size/worker_inQtask)==mpi__rank/worker_inQtask,iq=iqxini,iqxend)])
     if(ipr) write(stdo,ftox) 'FileMode: mpi_rank',mpi__rank,'worker_inQtask',worker_inQtask,'mpi__Qtask=',mpi__Qtask
     call flush(stdo)
     call wv_init_file(mreclx=mrecl, nw_i=nw_i)
     do iq = 2, iqxend
       qp = qibze(:,iq)
       if (mpi__Qtask(iq)) call build_screened_coulomb_step_kx(iq, qp, realomega, imagomega)
     enddo
     qp = qibze(:,1)
     if (mpi__Qtask(1)) call build_screened_coulomb_step_kx(1, qp, realomega, imagomega)
     call MPI_barrier(comm, ierr)
     call MPI__sendllw(iqxend, MPI__Qrank)
     if (MPI__rank == 0) call W0w0i(nw_i, nw, nq0i, niw, q0i, is_wc_m_basis=.true.)
     call MPI__FreeSplitXq()
  endif StreamingOrFile

  if(ipr) write(stdo,ftox) '--- end of hrcxq --- irank=',MPI__rank
  call cputid(0)
  call rx0( ' OK! hrcxq WV generated')

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
end subroutine hrcxq
end module m_hrcxq
