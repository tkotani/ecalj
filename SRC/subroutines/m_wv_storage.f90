!> m_wv_storage — handle for the W-V (W minus v) data movement between phases.
!!
!! The W-V matrices are produced per-iq by m_llw:WVRllwR / WVIllwI, modified at
!! iq=1 by m_w0w0i:modifyWV0 (Gamma-cell W(0) correction), and consumed per-kx by
!! m_sxcf_sc (correlation self-energy). Three modules share access; the data may
!! live on disk (__WVR.<iq> / __WVI.<iq> files) or in memory (Phase 1-C 3D
!! streaming buffer; iq=1/Gamma is produced last so W0w0i can operate on current
!! buffer directly without a separate iq=1 saved slot).
!!
!! Per ecalj convention this module exposes a singleton: subroutines take inputs
!! by argument and operate on module-level state ("output" via module variables).
!! Backends:
!!   WV_BACKEND_FILE       openm/writem/read/close on __WVR.<iq>/__WVI.<iq>
!!   WV_BACKEND_MEMORY_3D  store/load via the module-level 3D buffers
!!
!! API conventions:
!!   - call wv_init_file(mreclx, nw_i)  or  wv_init_memory_3d(nw_lo, nw_hi)
!!   - MEMORY_3D per-iq setup (called from x0kf_zxq on Qtask/root_k ranks):
!!       wv_assoc_real_buf(rcxq)          — wv_real_buf => rcxq (no alloc)
!!       wv_assoc_imag_buf(zxqi)          — wv_imag_buf => zxqi (no alloc; W written in-place)
!!   - MEMORY_3D per-iq setup (called on non-root_k ranks, after ngb broadcast):
!!       wv_alloc_zero_bufs(ngbx, niwx)   — alloc zero-filled real+imag bufs for Bcast receive
!!   - Writer side per iq:
!!       wv_open_iq_real_for_write(iq, comm=...);   wv_open_iq_imag_for_write(iq, comm=...)
!!       wv_put_real(iw, zw);                       wv_put_imag(iw, zw)
!!       wv_close_iq_for_write()
!!   - Reader side per iq:
!!       wv_open_iq_for_read(iq, want_real, want_imag)
!!       wv_get_real(iw, zw_out);                   wv_get_imag(iw, zw_out)
!!       wv_close_iq_for_read()
!!   - Modify mode (rank-0 read+modify+write on iq=1):
!!       wv_open_iq_real_for_modify(iq);            wv_open_iq_imag_for_modify(iq)
!!       wv_modify_get_real(iw, zw);                wv_modify_put_real(iw, zw)
!!       wv_modify_get_imag(iw, zw);                wv_modify_put_imag(iw, zw)
!!       wv_close_iq_for_modify()
!!   - call wv_dealloc() at the end.
!!   - MEMORY_3D streaming helpers:
!!       wv_zero_current() — no-op (zeroing now handled by x0kf_zxq / wv_alloc_zero_bufs)
!!       wv_sync_current(comm); wv_bcast_current(root, comm)
module m_wv_storage
  use m_kind, only: kp => kindrcxq
  use m_mpiio, only: openm, closem
#ifdef __MP
  use m_mpiio, only: writem => writem_c
#else
  use m_mpiio, only: writem
#endif
  implicit none
  private

  integer, parameter, public :: WV_BACKEND_FILE      = 1
  integer, parameter, public :: WV_BACKEND_MEMORY_3D = 2  !> per-iq streaming buffer + iq=1 saved

  ! ---- module-level singleton state ----
  integer, protected, public :: wv_backend = WV_BACKEND_FILE
  integer, protected, public :: wv_cur_iq  = 0
  integer :: wv_mreclx = 0
  integer :: wv_nw_i   = 0
  integer :: wv_real_unit = -1
  integer :: wv_imag_unit = -1
  logical :: wv_in_modify_mode = .false.  ! FILE: use standard write (not MPI-IO) in wv_put_*
  ! MEMORY_3D backend buffers — per-iq, sized (ngb, ngb, ...) not (nblochpmx, nblochpmx, ...).
  !   wv_real_buf  => rcxq (Qtask ranks) or wv_real_zero (non-Qtask, zero-filled)
  !   wv_imag_buf  => zxqi (Qtask/root_k) or wv_imag_zero (non-root_k, zero-filled)
  ! Lifetime: wv_real_buf/wv_imag_buf valid from x0kf_zxq/wv_alloc_zero_bufs until caller deallocs.
  complex(kp), pointer     :: wv_real_buf(:,:,:) => null()
  complex(kp), allocatable, target :: wv_real_zero(:,:,:)  ! backing store for non-Qtask zero buf
  complex(kp), pointer     :: wv_imag_buf(:,:,:) => null()
  complex(kp), allocatable, target :: wv_imag_zero(:,:,:)  ! backing store for non-root_k zero buf

  public :: wv_init_file, wv_init_memory_3d, wv_dealloc
  public :: wv_assoc_real_buf, wv_assoc_imag_buf, wv_alloc_zero_bufs
  public :: wv_open_iq_real_for_write, wv_open_iq_imag_for_write
  public :: wv_close_iq_for_write
  public :: wv_open_iq_for_read, wv_close_iq_for_read
  public :: wv_put_real, wv_put_imag, wv_get_real, wv_get_imag
  public :: wv_open_iq_real_for_modify, wv_open_iq_imag_for_modify
  public :: wv_close_iq_for_modify
  public :: wv_zero_current
  public :: wv_sync_current, wv_bcast_current

contains

  subroutine wv_init_file(mreclx, nw_i)
    !> Configure singleton for FILE backend. No I/O yet — actual openm happens
    !> in wv_open_iq_*_for_*. comm is per-call (openm-coordination scope).
    integer, intent(in) :: mreclx, nw_i
    wv_backend   = WV_BACKEND_FILE
    wv_mreclx    = mreclx
    wv_nw_i      = nw_i
    wv_cur_iq    = 0
    wv_real_unit = -1
    wv_imag_unit = -1
  end subroutine wv_init_file

  subroutine wv_init_memory_3d(nw_lo, nw_hi)
    !> Configure singleton for MEMORY_3D streaming. Mode setter only — buffers
    !> are set up per-iq via wv_assoc_real_buf / wv_assoc_imag_buf / wv_alloc_zero_bufs.
    integer, intent(in) :: nw_lo, nw_hi
    wv_backend = WV_BACKEND_MEMORY_3D
    wv_nw_i    = nw_lo
    wv_cur_iq  = 0
    nullify(wv_real_buf)
    if (allocated(wv_real_zero)) deallocate(wv_real_zero)
    nullify(wv_imag_buf)
    if (allocated(wv_imag_zero)) deallocate(wv_imag_zero)
  end subroutine wv_init_memory_3d

  subroutine wv_assoc_real_buf(rcxq_target)
    !> Qtask ranks: associate wv_real_buf with rcxq (no allocation).
    !> rcxq covers (ngb, ngb, nw_i:nwhis*npm) — wv_put/get_real use iw in nw_i:nw.
    complex(kp), intent(in), target :: rcxq_target(:,:,:)
    wv_real_buf => rcxq_target
  end subroutine wv_assoc_real_buf

  subroutine wv_assoc_imag_buf(zxqi_target)
    !> Qtask root_k ranks: associate wv_imag_buf with zxqi (no allocation).
    !> WVIllwI writes W(iw) in-place into zxqi(:,:,iw) via wv_put_imag.
    complex(kp), intent(in), target :: zxqi_target(:,:,:)
    wv_imag_buf => zxqi_target
  end subroutine wv_assoc_imag_buf

  subroutine wv_alloc_zero_bufs(ngbx, niwx, nw_lo, nw_hi)
    !> Non-Qtask ranks: allocate zero-filled real+imag buffers so they
    !> contribute 0 to the AllreduceSum while sharing the same size as Qtask.
    !> 1-based 3rd dim matches the pointer-from-dummy lower bound on Qtask ranks.
    !> ngb is broadcast from the Qtask rank before calling this.
    integer, intent(in) :: ngbx, niwx, nw_lo, nw_hi
    if (allocated(wv_real_zero)) deallocate(wv_real_zero)
    allocate(wv_real_zero(ngbx, ngbx, nw_hi - nw_lo + 1), source=cmplx(0,0,kp))
    wv_real_buf => wv_real_zero
    if (allocated(wv_imag_zero)) deallocate(wv_imag_zero)
    allocate(wv_imag_zero(ngbx, ngbx, niwx),               source=cmplx(0,0,kp))
    wv_imag_buf => wv_imag_zero
  end subroutine wv_alloc_zero_bufs

  subroutine wv_dealloc()
    !> Close any lingering file units (FILE) and release 3D buffers (MEMORY_3D).
    integer :: istat
    if (wv_backend == WV_BACKEND_FILE) then
       if (wv_real_unit > 0) then; istat = closem(wv_real_unit); wv_real_unit = -1; endif
       if (wv_imag_unit > 0) then; istat = closem(wv_imag_unit); wv_imag_unit = -1; endif
    else if (wv_backend == WV_BACKEND_MEMORY_3D) then
       nullify(wv_real_buf)                                    ! owned by rcxq in m_x0kf
       if (allocated(wv_real_zero)) deallocate(wv_real_zero)  ! non-Qtask backing store
       nullify(wv_imag_buf)                                    ! owned by zxqi in m_x0kf
       if (allocated(wv_imag_zero)) deallocate(wv_imag_zero)  ! non-root_k backing store
    endif
  end subroutine wv_dealloc

  ! ---- write mode (MPI-coordinated openm via comm) ----
  subroutine wv_open_iq_real_for_write(iq, comm)
    integer, intent(in) :: iq
    integer, intent(in), optional :: comm
    integer :: istat, c
    character(10) :: i2char
    wv_cur_iq = iq
    if (wv_backend /= WV_BACKEND_FILE) return
    c = -1
    if (present(comm)) c = comm
    istat = openm(newunit=wv_real_unit, file='__WVR.'//i2char(iq), &
                  recl=wv_mreclx, comm=c)
  end subroutine wv_open_iq_real_for_write

  subroutine wv_open_iq_imag_for_write(iq, comm)
    integer, intent(in) :: iq
    integer, intent(in), optional :: comm
    integer :: istat, c
    character(10) :: i2char
    wv_cur_iq = iq
    if (wv_backend /= WV_BACKEND_FILE) return
    c = -1
    if (present(comm)) c = comm
    istat = openm(newunit=wv_imag_unit, file='__WVI.'//i2char(iq), &
                  recl=wv_mreclx, comm=c)
  end subroutine wv_open_iq_imag_for_write

  subroutine wv_close_iq_for_write()
    integer :: istat
    if (wv_backend /= WV_BACKEND_FILE) return
    if (wv_real_unit > 0) then; istat = closem(wv_real_unit); wv_real_unit = -1; endif
    if (wv_imag_unit > 0) then; istat = closem(wv_imag_unit); wv_imag_unit = -1; endif
  end subroutine wv_close_iq_for_write

  subroutine wv_put_real(iw, zw)
    integer,     intent(in) :: iw
    complex(kp), intent(in) :: zw(:,:)
    integer :: istat, n
    if (wv_backend == WV_BACKEND_FILE) then
       if (wv_in_modify_mode) then
          write(wv_real_unit, rec=iw - wv_nw_i + 1) zw  ! rank-0 only, standard Fortran write
       else
          istat = writem(wv_real_unit, rec=iw - wv_nw_i + 1, data=zw)
       endif
    else
       n = size(wv_real_buf, 1)
       wv_real_buf(1:n, 1:n, iw - wv_nw_i + 1) = zw(1:n, 1:n)
    endif
  end subroutine wv_put_real

  subroutine wv_put_imag(iw, zw)
    integer,     intent(in) :: iw
    complex(kp), intent(in) :: zw(:,:)
    integer :: istat, n
    if (wv_backend == WV_BACKEND_FILE) then
       if (wv_in_modify_mode) then
          write(wv_imag_unit, rec=iw) zw  ! rank-0 only, standard Fortran write
       else
          istat = writem(wv_imag_unit, rec=iw, data=zw)
       endif
    else
       n = size(wv_imag_buf, 1)
       wv_imag_buf(1:n, 1:n, iw) = zw(1:n, 1:n)
    endif
  end subroutine wv_put_imag

  ! ---- read mode (each rank opens independently for FILE) ----
  subroutine wv_open_iq_for_read(iq, want_real, want_imag)
    integer, intent(in) :: iq
    logical, intent(in) :: want_real, want_imag
    character(10) :: i2char
    wv_cur_iq = iq
    if (wv_backend /= WV_BACKEND_FILE) return
    if (want_real) then
       open(newunit=wv_real_unit, file='__WVR.'//i2char(iq), action='read', &
            form='unformatted', access='direct', recl=wv_mreclx)
    endif
    if (want_imag) then
       open(newunit=wv_imag_unit, file='__WVI.'//i2char(iq), action='read', &
            form='unformatted', access='direct', recl=wv_mreclx)
    endif
  end subroutine wv_open_iq_for_read

  subroutine wv_close_iq_for_read()
    if (wv_backend /= WV_BACKEND_FILE) return
    if (wv_real_unit > 0) then; close(wv_real_unit); wv_real_unit = -1; endif
    if (wv_imag_unit > 0) then; close(wv_imag_unit); wv_imag_unit = -1; endif
  end subroutine wv_close_iq_for_read

  subroutine wv_get_real(iw, zw)
    integer,     intent(in)  :: iw
    complex(kp), intent(out) :: zw(:,:)
    integer :: n
    if (wv_backend == WV_BACKEND_FILE) then
       read(wv_real_unit, rec=iw - wv_nw_i + 1) zw
    else
       n = size(wv_real_buf, 1)
       zw(1:n, 1:n) = wv_real_buf(1:n, 1:n, iw - wv_nw_i + 1)
    endif
  end subroutine wv_get_real

  subroutine wv_get_imag(iw, zw)
    integer,     intent(in)  :: iw
    complex(kp), intent(out) :: zw(:,:)
    integer :: n
    if (wv_backend == WV_BACKEND_FILE) then
       read(wv_imag_unit, rec=iw) zw
    else
       n = size(wv_imag_buf, 1)
       zw(1:n, 1:n) = wv_imag_buf(1:n, 1:n, iw)
    endif
  end subroutine wv_get_imag

  ! ---- modify mode (rank-0 read+write, no MPI coordination) ----
  ! FILE: read+write on direct-access unit. MEMORY_3D: read+write on iq=1 saved slot.
  subroutine wv_open_iq_real_for_modify(iq)
    integer, intent(in) :: iq
    character(10) :: i2char
    wv_cur_iq = iq
    wv_in_modify_mode = .true.
    if (wv_backend /= WV_BACKEND_FILE) return
    open(newunit=wv_real_unit, file='__WVR.'//i2char(iq), &
         form='unformatted', status='old', access='direct', recl=wv_mreclx)
  end subroutine wv_open_iq_real_for_modify

  subroutine wv_open_iq_imag_for_modify(iq)
    integer, intent(in) :: iq
    character(10) :: i2char
    wv_cur_iq = iq
    wv_in_modify_mode = .true.
    if (wv_backend /= WV_BACKEND_FILE) return
    open(newunit=wv_imag_unit, file='__WVI.'//i2char(iq), &
         form='unformatted', status='old', access='direct', recl=wv_mreclx)
  end subroutine wv_open_iq_imag_for_modify

  subroutine wv_close_iq_for_modify()
    wv_in_modify_mode = .false.
    if (wv_backend /= WV_BACKEND_FILE) return
    if (wv_real_unit > 0) then; close(wv_real_unit); wv_real_unit = -1; endif
    if (wv_imag_unit > 0) then; close(wv_imag_unit); wv_imag_unit = -1; endif
  end subroutine wv_close_iq_for_modify

  ! =============================================================
  ! MEMORY_3D streaming-specific helpers
  ! =============================================================
  subroutine wv_zero_current()
    !> No-op in MEMORY_3D mode: zeroing is now handled per-iq by
    !> x0kf_zxq (rcxq/zxqi) and wv_alloc_zero_bufs.
    if (wv_backend /= WV_BACKEND_MEMORY_3D) return
  end subroutine wv_zero_current

  subroutine wv_sync_current(comm)
    !> Bcast the current 3D buffers from mpi__root_k (rank 0 in comm) to all
    !> ranks. Non-root_k ranks must have buffers allocated (via wv_alloc_zero_bufs)
    !> before calling so MPI_Bcast has a valid receive buffer.
    use mpi
    integer, intent(in) :: comm
    integer :: ierr, mpi_type
    if (wv_backend /= WV_BACKEND_MEMORY_3D) return
#ifdef __MP
    mpi_type = MPI_COMPLEX
#else
    mpi_type = MPI_DOUBLE_COMPLEX
#endif
    if (associated(wv_real_buf)) &
      call MPI_Bcast(wv_real_buf, size(wv_real_buf), mpi_type, 0, comm, ierr)
    if (associated(wv_imag_buf)) &
      call MPI_Bcast(wv_imag_buf, size(wv_imag_buf), mpi_type, 0, comm, ierr)
  end subroutine wv_sync_current

  subroutine wv_bcast_current(root, comm)
    !> Broadcast the current buffer (W0w0i-corrected on root) to all ranks.
    use mpi
    integer, intent(in) :: root, comm
    integer :: ierr, mpi_type
    if (wv_backend /= WV_BACKEND_MEMORY_3D) return
#ifdef __MP
    mpi_type = MPI_COMPLEX
#else
    mpi_type = MPI_DOUBLE_COMPLEX
#endif
    if (associated(wv_real_buf)) &
      call MPI_Bcast(wv_real_buf, size(wv_real_buf), mpi_type, root, comm, ierr)
    if (associated(wv_imag_buf)) &
      call MPI_Bcast(wv_imag_buf, size(wv_imag_buf), mpi_type, root, comm, ierr)
  end subroutine wv_bcast_current

end module m_wv_storage
