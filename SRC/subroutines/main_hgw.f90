!> GW main driver: build Im(chi0), Hilbert transform → W; stream W in-process
!> for exchange and correlation self-energy via hsfp0_sc kernels.
module m_hgw
  contains
subroutine hgw(do_correlation, do_exchange)
  !> Streaming flow (SHM backend), unified loop iq = iqxend → 1:
  !>   iq > nqibz: auxiliary q0 — build W/llw only (no Sc).
  !>   iq <= nqibz: regular — build W, consume for Sc (W0w0i at iq=1).
  !> MPI layout: mpi__size = n_qgroup × worker_inQtask (q × per-q-group).
  !>   MPI__InitQgroups sets comm_q = intra-node communicator (ppn ranks).
  !>   MPI__SplitXq creates comm_k_xq for build_screened_coulomb (called just before main loop).
  !>   MPI__SplitSxc creates comm_k_sxc for exchange (n_bpara=1) and correlation (n_bpara>=1).
  !> Loop gate: both regular and aux iq assigned by unified LPT (mpi_assign_qtask_lpt);
  !>   aux weight ≈ avg_regular_weight/3, reflecting W-build:Sc ≈ 1:3 cost ratio.
  use m_ReadEfermi,only: Readefermi
  use m_readqg,only: Readngmx2
  use m_hamindex,only: Readhamindex, symgg=>symops, ngrp
  use m_readeigen,only: Init_readeigen,Init_readeigen2
  use m_read_bzdata,only:Read_bzdata, nq0i,nq0iadd,nqibz,q0i
  use m_genallcf_v3,only: Genallcf_v3
  use m_rdpp,only: mrecl,nblochpmx,nprecx
  use m_zmel,only: Mptauof_zmel
  use m_itq,only:  Setitq
  use m_freq,only: Getfreq2,freq_r,nw_i,nw,niw, nwhis_hgw=>nwhis, npm_hgw=>npm
  use m_w0w0i,only: W0w0i
  use m_readgwinput,only: ReadGwinputKeys
  use m_GWinput,only: mpi_worker_exch, mpi_worker_corr, zmel_batch_gb
  use m_qbze,only:  Setqbze,qibze
  use m_llw,only: MPI__irecvllw_q, MPI__isendllw_q, MPI__waitllw, MPI__llw_alloc_bufs
  use m_mpi,only: MPI__Initialize, MPI__InitQgroups, MPI__FreeQgroups, &
                & MPI__SplitXq, MPI__FreeXq, MPI__SplitSxc, MPI__FreeSxc, &
                & MPI__AutoSetup, MPI__PrintSummary, &
                & MPI__root, MPI__rank, MPI__size, MPI__consoleout, comm, &
                & ipr, comm_q, mpi__rank_q, mpi__root_q, &
                & worker_inQtask, n_qgroup, iq_qgroup, qgroup_root
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
  use m_sxcf_count,only: iq1_dest, reg_grp_assign, aux_grp_assign
  use m_stopwatch, only: stopwatch, stopwatch_begin, stopwatch_pause_and_show
  use mpi
  implicit none
  logical, intent(in) :: do_correlation, do_exchange
  integer :: iq, iqxend, iw, ifwd, verbose, ifif, ierr
  integer :: n_bpara_xq, n_kpara_xq, n_bpara_sxc, n_kpara_sxc, worker_auto, worker_exch
  integer :: src_group
  real(8) :: ua=1d0, qp(3)
  type(stopwatch) :: sw_q
  logical :: debug=.false., realomega, imagomega
  logical :: hx0, iprintx=.false.
  call MPI__Initialize()
  call gpu_init(comm, share_gpu=.true.)
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

  call MPI__AutoSetup(nblochpmx, nwhis_hgw, npm_hgw, niw, nqibz, &
                      n_bpara_sxc_hint=1, &
                      worker_out=worker_auto, worker_exch_out=worker_exch, &
                      n_bpara_xq_out=n_bpara_xq, n_kpara_xq_out=n_kpara_xq, &
                      n_bpara_sxc_out=n_bpara_sxc, n_kpara_sxc_out=n_kpara_sxc, &
                      zmel_per_rank_gb=zmel_batch_gb)
  if (mpi_worker_exch > 0) then
    if (mod(MPI__size, mpi_worker_exch) /= 0) &
      call rx('mpi_worker_exch must divide mpi__size')
    worker_exch = mpi_worker_exch
  endif
  if (mpi_worker_corr > 0) then
    if (mod(MPI__size, mpi_worker_corr) /= 0) &
      call rx('mpi_worker_corr must divide mpi__size')
    worker_auto = mpi_worker_corr
  endif

  if(MPI__root) call writewvfreq()
  iqxend = nqibz + nq0i + nq0iadd

  if (do_exchange) then
     call MPI__InitQgroups(worker_in=worker_exch)
     if(ipr) write(stdo,'(1X,A,3I5)') 'hgw(exch): worker_inQtask n_qgroup iqxend:', &
                                        worker_inQtask, n_qgroup, iqxend
     call MPI__SplitSxc(1, worker_inQtask)
     if(ipr) write(stdo,ftox) ' hgw: starting in-process hsfp0_sc(--job=1) exchange phase'
     call hsfp0_sc(skip_init=.true., skip_rx0=.true., ixc_in=1)
     call MPI__FreeSxc()
     call MPI__FreeQgroups()
  endif

  if(sum(qibze(:,1)**2)>1d-10) call rx(' hgw: sanity check. |q(iq=1)| /= 0')

  call MPI__InitQgroups(worker_in=worker_auto)
  if(ipr) write(stdo,'(1X,A,3I5)') 'hgw(corr): worker_inQtask n_qgroup iqxend:', &
                                     worker_inQtask, n_qgroup, iqxend
  call MPI__SplitXq(n_bpara_xq, n_kpara_xq)
  call MPI__SplitSxc(n_bpara_sxc, n_kpara_sxc)

  call hsfp0_sc_setup(skip_init=.true., ixc_in=2)
  call sxcf_correlation_init(hs_ef, hs_esmr, hs_nspinmx)
  ! Print after sxcf_correlation_init so q-point assignment (reg_grp_assign/aux_grp_assign)
  ! is available for display.  MPI__PrintSummary is read-only; reordering is safe.
  call MPI__PrintSummary(n_bpara_xq, n_kpara_xq, n_bpara_sxc, n_kpara_sxc, &
                          reg_grp=reg_grp_assign, aux_grp=aux_grp_assign)

  if(ipr) write(stdo,ftox) 'hgw: unified loop iqxend→1, mpi_rank=',MPI__rank
  call flush(stdo)

  ! Ensure llw/llwI/wmuk are allocated before pre-posting Irecvs.
  call MPI__llw_alloc_bufs()

  ! Pre-post Irecvs for auxiliary llw at iq1_dest (qroot of the group assigned iq=1 by LPT).
  ! Aux iq assigned by unified LPT (aux_grp_assign); src = qgroup_root of that group.
  ! MPI__irecvllw_q is a no-op on all ranks except dest (iq1_dest).
  do iq = nqibz+1, iqxend
    src_group = aux_grp_assign(iq-nqibz)
    call MPI__irecvllw_q(iq-nqibz, qgroup_root(src_group), iq1_dest)
  enddo

  do iq = iqxend, 1, -1
    if (iq > nqibz .and. aux_grp_assign(iq-nqibz) /= iq_qgroup) cycle
    if (iq <= nqibz .and.     reg_grp_assign(iq)  /= iq_qgroup) cycle
    qp = qibze(:,iq)

    if (mpi__root_q) write(stdo,'(1X,A,I0,A)') 'hgw iq=',iq,': W-build'
    if (mpi__root_q) call flush(stdo)
    call stopwatch_begin(sw_q, 'W-build(per-q-point)')
    call build_screened_coulomb_step_kx(iq, qp, realomega, imagomega)
    if (mpi__root_q) call stopwatch_pause_and_show(sw_q)

    if (iq > nqibz) then
      ! Auxiliary q-point: send llw to iq1_dest (qroot of group assigned iq=1).
      call MPI__isendllw_q(iq-nqibz, qgroup_root(iq_qgroup), iq1_dest)
    else
      if (iq == 1) then
        ! Group assigned iq=1 collects all aux llw then calls W0w0i.
        call MPI__waitllw()
        if (mpi__root_q) call W0w0i(nw_i, nw, nq0i, niw, q0i, is_wc_m_basis=.true.)
        call MPI_barrier(comm_q, ierr)
      end if
      if (mpi__root_q) write(stdo,'(1X,A,I0,A)') 'hgw iq=',iq,': Sc'
      if (mpi__root_q) call flush(stdo)
      call stopwatch_begin(sw_q, 'Sc(per-q-point)')
      call sxcf_correlation_step_kx(iq, hs_ef, hs_esmr, hs_nspinmx)
      if (mpi__root_q) call stopwatch_pause_and_show(sw_q)
    end if
    call wv_dealloc()
  enddo
  call MPI__waitllw()  ! ensure pending Isends complete

  call sxcf_correlation_finalize()
  call hsfp0_sc_writeout(skip_rx0=.true.)
  call MPI__FreeXq()
  call MPI__FreeSxc()
  call MPI__FreeQgroups()

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
