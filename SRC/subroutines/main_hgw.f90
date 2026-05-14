!> GW main driver: build Im(chi0), Hilbert transform → W; stream W in-process
!> for exchange and correlation self-energy via hsfp0_sc kernels.
module m_hgw
  contains
subroutine hgw(do_correlation, do_exchange)
  !> Streaming flow (SHM backend), unified loop iq = iqxend → 1:
  !>   iq > nqibz: auxiliary q0 — build W/llw; only q-group 0 processes these.
  !>   iq <= nqibz: regular — build W, consume for Sc (W0w0i at iq=1).
  !> MPI layout: mpi__size = n_qgroup × n_kpara × n_bpara (q × k × ω).
  !>   Exchange runs with SplitXq(1, mpi__size) for full k-parallel; then
  !>   FreeSplitXq/SplitXq(n_bpara, n_kpara) for the correlation q-loop.
  !> mpi__Qtask(iq) gates which q-group processes each iq; auxiliary iq's
  !>   are always assigned to q-group 0 so rank 0 (root_k) holds llw.
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
  use m_llw,only: MPI__irecvllw_q, MPI__isendllw_q, MPI__waitllw
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
  use mpi
  implicit none
  logical, intent(in) :: do_correlation, do_exchange
  integer :: iq, iqxend, iw, ifwd, verbose, ifif, ierr
  integer :: n_bpara, n_kpara, worker_inQtask, n_qgroup, iq_qgroup
  logical, allocatable :: mpi__Qtask(:)
  real(8) :: ua=1d0, qp(3)
  logical :: debug=.false., realomega, imagomega
  logical :: hx0, iprintx=.false.
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
  if(ipr) write(stdo,'(1X,A,I5)') 'MPI: nranks (omega-parallel):', mpi__size

  ! Exchange: SplitXq(1, mpi__size) gives full k-parallel (mpi__size_k=mpi__size).
  ! comm_b/comm_k must be initialized before hsfp0_sc; skip_init=.true. omits the
  ! internal SplitXq call, so we do it explicitly here.
  call MPI__SplitXq(1, mpi__size)
  if (do_exchange) then
     if(ipr) write(stdo,ftox) ' hgw: starting in-process hsfp0_sc(--job=1) exchange phase'
     call hsfp0_sc(skip_init=.true., skip_rx0=.true., ixc_in=1)
  endif
  call MPI__FreeSplitXq()

  if(sum(qibze(:,1)**2)>1d-10) call rx(' hgw: sanity check. |q(iq=1)| /= 0')

  ! Correlation: q × k × ω 3-way parallel.
  ! n_bpara (ω), n_kpara (k), n_qgroup = mpi__size/(n_bpara*n_kpara) (q).
  ! Default n_bpara=1, n_kpara=mpi__size → 1 q-group (backward compatible).
  n_bpara = 1
  n_kpara = mpi__size
  worker_inQtask = n_bpara * n_kpara
  n_qgroup = mpi__size / worker_inQtask
  iq_qgroup = mpi__rank / worker_inQtask
  if(ipr) write(stdo,'(1X,A,4I5)') 'hgw: n_bpara n_kpara n_qgroup iqxend:', &
                                     n_bpara, n_kpara, n_qgroup, iqxend
  call MPI__SplitXq(n_bpara, n_kpara)

  ! mpi__Qtask(iq): true if this rank's q-group should process iq.
  !   Regular iq (1..nqibz): round-robin among q-groups.
  !   Auxiliary iq (nqibz+1..iqxend): q-group 0 only; rank 0 (root_k) writes
  !   llw into WVRllwR/WVIllwI directly — no MPI transfer needed.
  allocate(mpi__Qtask(1:iqxend))
  mpi__Qtask(1:nqibz)       = [(mod(iq-1, n_qgroup) == iq_qgroup, iq=1,nqibz)]
  mpi__Qtask(nqibz+1:iqxend) = (iq_qgroup == 0)
  if(ipr) write(stdo,ftox) 'hgw: mpi_rank iq_qgroup mpi__Qtask=', &
                             mpi__rank, iq_qgroup, mpi__Qtask

  call hsfp0_sc_setup(skip_init=.true., ixc_in=2)
  call sxcf_correlation_init(hs_ef, hs_esmr, hs_nspinmx)
  if(ipr) write(stdo,ftox) 'hgw: unified loop iqxend→1, mpi_rank=',MPI__rank
  call flush(stdo)

  ! Pre-post Irecvs for auxiliary llw on rank 0.
  ! Auxiliary iq's processed by q-group 0; rank 0 is always root_k of q-group 0,
  ! so src=dest=0 → all calls are no-ops (llw already on rank 0).
  do iq = nqibz+1, iqxend
    call MPI__irecvllw_q(iq-nqibz, 0, 0)
  enddo

  do iq = iqxend, 1, -1
    if (.not. mpi__Qtask(iq)) cycle
    qp = qibze(:,iq)
    call build_screened_coulomb_step_kx(iq, qp, realomega, imagomega)
    if (iq > nqibz) then
      ! Auxiliary q-point: llw written by rank 0 directly; send is a no-op.
      call MPI__isendllw_q(iq-nqibz, 0, 0)
    else
      ! Regular q-point: barrier so all ranks in q-group see W in SHM, then consume.
      call MPI_barrier(comm_q, ierr)
      if (iq == 1) then
        ! iq=1 always in q-group 0 (round-robin: mod(0,n_qgroup)=0).
        ! Wait for all auxiliary llw, then apply W0w0i correction on rank 0.
        call MPI__waitllw()
        if (MPI__rank == 0) call W0w0i(nw_i, nw, nq0i, niw, q0i, is_wc_m_basis=.true.)
        call MPI_barrier(comm_q, ierr)
      end if
      call sxcf_correlation_step_kx(iq, hs_ef, hs_esmr, hs_nspinmx)
      call MPI_barrier(comm_q, ierr)
    end if
  enddo

  call sxcf_correlation_finalize()
  call wv_dealloc()
  call hsfp0_sc_writeout(skip_rx0=.true.)
  call MPI__FreeSplitXq()
  deallocate(mpi__Qtask)

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
