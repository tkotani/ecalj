!> Pre-compute k-point dependent data for rsibl using MPI shared memory.
!! All MPI ranks compute setup for their own k-points in parallel,
!! writing results into shared memory arrays (padded to ng_max).
!! Any rank can then read all k-points' data without MPI communication.
module m_rsibl_setup
  use m_lgunit,only: stdo
  use m_sharedmem
  implicit none
  public :: rsibl_setup_all, rsibl_setup_clean
  private

  integer, parameter :: nermx = 100
  integer, protected, public, save :: rsibl_net, rsibl_nrt, rsibl_ltop, rsibl_nlmtop
  real(8), public, save :: etab_s(nermx), rtab_s(nermx)
  integer, public, save :: ipet_s(10,5,nermx), iprt_s(10,5,nermx)
  logical, save :: tbhsi_done = .false.

  integer, protected, public, save :: rsibl_ng_max = 0
  integer, pointer, public, save :: rsibl_ng_all(:) => null()
  integer, pointer, public, save :: rsibl_nlmto_all(:) => null()
  integer, pointer, public, save :: rsibl_napw_all(:) => null()
  integer, pointer, public, save :: rsibl_ngmax_all(:) => null()
  real(8), pointer, public, save :: rsibl_he_all(:,:,:) => null()
  real(8), pointer, public, save :: rsibl_hr_all(:,:,:) => null()
  real(8), pointer, public, save :: rsibl_yl_all(:,:,:) => null()
  real(8), pointer, public, save :: rsibl_ogv_all(:,:,:) => null()
  real(8), pointer, public, save :: rsibl_wogq_all(:,:,:) => null()
  integer, pointer, public, save :: rsibl_iv_all(:,:) => null()
  integer, pointer, public, save :: rsibl_igv_all(:,:,:) => null()
  integer, pointer, public, save :: rsibl_ivp_all(:,:) => null()

  integer, save :: nwin = 0
  integer, save :: win_ids(20) = 0
  logical, public, save :: rsibl_setup_done = .false.

contains

  subroutine rsibl_setup_all(nkp, iqlist, isplist, nlocal)
    use m_lmfinit,only: alat=>lat_alat, nspec, n0, nkap0
    use m_lattic,only: qlat=>lat_qlat, plat=>lat_plat
    use m_supot,only: n1, n2, n3, lat_ng, gmax=>lat_gmax
    use m_sugcut,only: ngcut
    use m_hsibl,only: hsibl1
    use m_shortn3,only: gvlst2
    use m_nvfortran,only: findloc
    use m_igv2x,only: m_igv2x_getiq, t_igv2x_data, nbandmx
    use m_qplist,only: qplist
    use m_MPItk,only: procid_l=>procid, master_mpi, comm
    implicit none
    integer, intent(in) :: nkp, nlocal
    integer, intent(in) :: iqlist(nlocal), isplist(nlocal)
    type(t_igv2x_data) :: kdat
    integer :: iq, idat, ng, ig, jg, ixx(1), iprint, napw_max, ierr
    real(8) :: xx(1), q(3), q0(3)
    real(8), allocatable :: w_og2(:)

    if(rsibl_setup_done) return
    call tcn('rsibl_setup_all')

    ! k-independent tables
    if(.not. tbhsi_done) then
      call tbhsi(nspec, nermx, rsibl_net, etab_s, ipet_s, rsibl_nrt, rtab_s, iprt_s, rsibl_ltop)
      rsibl_nlmtop = (rsibl_ltop+1)**2
      tbhsi_done = .true.
    endif

    ! Determine ng_max across all k-points
    rsibl_ng_max = 0
    do idat = 1, nlocal
      iq = iqlist(idat); q = qplist(:, iq); ng = lat_ng
      call pshpr(iprint()-30)
      call gvlst2(alat, plat, q, n1, n2, n3, 0d0, gmax, [0], 000, 0, ng, ixx, xx, ixx)
      call poppr
      rsibl_ng_max = max(rsibl_ng_max, ng)
    enddo
    call mpibc2_int_max(rsibl_ng_max)
    napw_max = nbandmx

    ! Allocate shared memory (collective)
    call shm_init()
    call shm_alloc_i4_1d(rsibl_ng_all, nkp, win_ids(1)); nwin=1
    call shm_alloc_i4_1d(rsibl_nlmto_all, nkp, win_ids(2)); nwin=2
    call shm_alloc_i4_1d(rsibl_napw_all, nkp, win_ids(3)); nwin=3
    call shm_alloc_i4_1d(rsibl_ngmax_all, nkp, win_ids(4)); nwin=4
    call shm_alloc_r8_3d(rsibl_he_all, rsibl_ng_max, rsibl_net, nkp, win_ids(5)); nwin=5
    call shm_alloc_r8_3d(rsibl_hr_all, rsibl_ng_max, rsibl_nrt, nkp, win_ids(6)); nwin=6
    call shm_alloc_r8_3d(rsibl_yl_all, rsibl_ng_max, rsibl_nlmtop, nkp, win_ids(7)); nwin=7
    call shm_alloc_r8_3d(rsibl_ogv_all, rsibl_ng_max, 3, nkp, win_ids(8)); nwin=8
    call shm_alloc_r8_3d(rsibl_wogq_all, rsibl_ng_max, 3, nkp, win_ids(9)); nwin=9
    call shm_alloc_i4_2d(rsibl_iv_all, rsibl_ng_max*3, nkp, win_ids(10)); nwin=10
    call shm_alloc_i4_3d(rsibl_igv_all, rsibl_ng_max, 3, nkp, win_ids(11)); nwin=11
    call shm_alloc_i4_2d(rsibl_ivp_all, napw_max, nkp, win_ids(12)); nwin=12

    ! Zero-fill on shared rank 0
    if(shm_rank() == 0) then
      rsibl_ng_all=0; rsibl_nlmto_all=0; rsibl_napw_all=0; rsibl_ngmax_all=0
      rsibl_he_all=0d0; rsibl_hr_all=0d0; rsibl_yl_all=0d0
      rsibl_ogv_all=0d0; rsibl_wogq_all=0d0
      rsibl_iv_all=0; rsibl_igv_all=0; rsibl_ivp_all=0
    endif
    do idat=1,nwin; call shm_barrier(win_ids(idat)); enddo

    ! Each rank fills its k-points
    do idat = 1, nlocal
      iq = iqlist(idat); q = qplist(:, iq)
      call m_igv2x_getiq(iq, kdat)
      ng = lat_ng
      call pshpr(iprint()-30)
      call gvlst2(alat, plat, q, n1, n2, n3, 0d0, gmax, [0], 000, 0, ng, ixx, xx, ixx)
      rsibl_ng_all(iq) = ng
      rsibl_nlmto_all(iq) = kdat%ndimh - kdat%napw
      rsibl_napw_all(iq) = kdat%napw
      call gvlst2(alat, plat, q, n1, n2, n3, 0d0, gmax, [0], 509, ng, ng, &
           rsibl_iv_all(1,iq), rsibl_ogv_all(1,1,iq), rsibl_igv_all(1,1,iq))
      call poppr
      if(kdat%napw > 0) then
        do ig = 1, kdat%napw
          rsibl_ivp_all(ig, iq) = findloc( &
               [(sum(abs(rsibl_igv_all(jg,:,iq)-kdat%igv2x(:,ig)))==0, jg=1,ng)], value=.true., dim=1)
        enddo
      endif
      allocate(w_og2(ng)); q0=0d0
      if(rsibl_nlmto_all(iq)>0) &
        call hsibl1(rsibl_net, etab_s, rsibl_nrt, rtab_s, rsibl_ltop, alat, q0, ng, &
             rsibl_ogv_all(1,1,iq), rsibl_wogq_all(1,1,iq), w_og2, &
             rsibl_yl_all(1,1,iq), rsibl_he_all(1,1,iq), rsibl_hr_all(1,1,iq))
      deallocate(w_og2)
      rsibl_ngmax_all(iq) = min(maxval(ngcut), ng)
    enddo

    ! Barrier
    do idat=1,nwin; call shm_barrier(win_ids(idat)); enddo
    rsibl_setup_done = .true.
    if(master_mpi) write(stdo,'(a,i5,a,i6)') ' rsibl_setup_all: nkp=',nkp,' ng_max=',rsibl_ng_max
    call tcx('rsibl_setup_all')
  end subroutine

  subroutine rsibl_setup_clean()
    integer :: i
    do i=1,nwin; call shm_free(win_ids(i)); win_ids(i)=0; enddo
    nwin=0; rsibl_setup_done=.false.
    nullify(rsibl_ng_all,rsibl_nlmto_all,rsibl_napw_all,rsibl_ngmax_all)
    nullify(rsibl_he_all,rsibl_hr_all,rsibl_yl_all,rsibl_ogv_all,rsibl_wogq_all)
    nullify(rsibl_iv_all,rsibl_igv_all,rsibl_ivp_all)
  end subroutine

  subroutine mpibc2_int_max(val)
    use m_MPItk, only: comm
    use mpi, only: MPI_INTEGER, MPI_MAX
    implicit none
    integer, intent(inout) :: val
    integer :: recv, ierr
    call MPI_ALLREDUCE(val, recv, 1, MPI_INTEGER, MPI_MAX, comm, ierr)
    val = recv
  end subroutine

end module m_rsibl_setup
