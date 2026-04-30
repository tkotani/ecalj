!> m_wv_storage — handle for the W-V (W minus v) data movement between phases.
!!
!! The W-V matrices are produced per-iq by m_llw:WVRllwR / WVIllwI, modified at
!! iq=1 by m_w0w0i:modifyWV0 (Gamma-cell W(0) correction), and consumed per-kx by
!! m_sxcf_sc (correlation self-energy). Three modules share access; the data may
!! live on disk (__WVR.<iq> / __WVI.<iq> files) or in memory (Phase 1-C 3D
!! streaming buffer with current iq + saved iq=1 slot).
!!
!! Per ecalj convention this module exposes a singleton: subroutines take inputs
!! by argument and operate on module-level state ("output" via module variables).
!! Backends:
!!   WV_BACKEND_FILE       openm/writem/read/close on __WVR.<iq>/__WVI.<iq>
!!   WV_BACKEND_MEMORY_3D  store/load via the module-level 3D buffers
!!
!! API conventions:
!!   - call wv_init_file(mreclx, nw_i)  or  wv_init_memory_3d(nbpmx, nw_lo, nw_hi, niwx)
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
!!       wv_zero_current(); wv_save_iq1(); wv_restore_iq1()
!!       wv_sync_current(comm); wv_bcast_iq1(root, comm)
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
  ! MEMORY_3D backend buffers. Shape:
  !   wv_real_buf  (nblochpmx, nblochpmx, nw_i:nw)
  !   wv_imag_buf  (nblochpmx, nblochpmx, 1:niw)
  ! The "current" buffer is overwritten as iq advances; the iq1 slot is saved
  ! by wv_save_iq1 before W0w0i and restored after wv_bcast_iq1.
  complex(kp), allocatable :: wv_real_buf(:,:,:)
  complex(kp), allocatable :: wv_imag_buf(:,:,:)
  complex(kp), allocatable :: wv_real_buf_iq1(:,:,:)
  complex(kp), allocatable :: wv_imag_buf_iq1(:,:,:)

  public :: wv_init_file, wv_init_memory_3d, wv_dealloc
  public :: wv_open_iq_real_for_write, wv_open_iq_imag_for_write
  public :: wv_close_iq_for_write
  public :: wv_open_iq_for_read, wv_close_iq_for_read
  public :: wv_put_real, wv_put_imag, wv_get_real, wv_get_imag
  public :: wv_open_iq_real_for_modify, wv_open_iq_imag_for_modify
  public :: wv_modify_get_real, wv_modify_put_real
  public :: wv_modify_get_imag, wv_modify_put_imag
  public :: wv_close_iq_for_modify
  public :: wv_zero_current
  public :: wv_save_iq1, wv_restore_iq1
  public :: wv_sync_current, wv_bcast_iq1

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

  subroutine wv_init_memory_3d(nbpmx, nw_lo, nw_hi, niwx)
    !> Configure singleton for MEMORY_3D streaming. Allocates 3D current buffer
    !> and 3D iq=1 saved buffer, both initialized to zero.
    integer, intent(in) :: nbpmx, nw_lo, nw_hi, niwx
    wv_backend = WV_BACKEND_MEMORY_3D
    wv_nw_i    = nw_lo
    wv_cur_iq  = 0
    if (allocated(wv_real_buf))     deallocate(wv_real_buf)
    if (allocated(wv_imag_buf))     deallocate(wv_imag_buf)
    if (allocated(wv_real_buf_iq1)) deallocate(wv_real_buf_iq1)
    if (allocated(wv_imag_buf_iq1)) deallocate(wv_imag_buf_iq1)
    allocate(wv_real_buf(nbpmx, nbpmx, nw_lo:nw_hi),     source=cmplx(0,0,kp))
    allocate(wv_imag_buf(nbpmx, nbpmx, 1:niwx),          source=cmplx(0,0,kp))
    allocate(wv_real_buf_iq1(nbpmx, nbpmx, nw_lo:nw_hi), source=cmplx(0,0,kp))
    allocate(wv_imag_buf_iq1(nbpmx, nbpmx, 1:niwx),      source=cmplx(0,0,kp))
  end subroutine wv_init_memory_3d

  subroutine wv_dealloc()
    !> Close any lingering file units (FILE) and free 3D buffers (MEMORY_3D).
    integer :: istat
    if (wv_backend == WV_BACKEND_FILE) then
       if (wv_real_unit > 0) then; istat = closem(wv_real_unit); wv_real_unit = -1; endif
       if (wv_imag_unit > 0) then; istat = closem(wv_imag_unit); wv_imag_unit = -1; endif
    else if (wv_backend == WV_BACKEND_MEMORY_3D) then
       if (allocated(wv_real_buf))     deallocate(wv_real_buf)
       if (allocated(wv_imag_buf))     deallocate(wv_imag_buf)
       if (allocated(wv_real_buf_iq1)) deallocate(wv_real_buf_iq1)
       if (allocated(wv_imag_buf_iq1)) deallocate(wv_imag_buf_iq1)
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
       istat = writem(wv_real_unit, rec=iw - wv_nw_i + 1, data=zw)
    else
       n = size(wv_real_buf, 1)
       wv_real_buf(1:n, 1:n, iw) = zw(1:n, 1:n)
    endif
  end subroutine wv_put_real

  subroutine wv_put_imag(iw, zw)
    integer,     intent(in) :: iw
    complex(kp), intent(in) :: zw(:,:)
    integer :: istat, n
    if (wv_backend == WV_BACKEND_FILE) then
       istat = writem(wv_imag_unit, rec=iw, data=zw)
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
       zw(1:n, 1:n) = wv_real_buf(1:n, 1:n, iw)
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
    if (wv_backend /= WV_BACKEND_FILE) return
    open(newunit=wv_real_unit, file='__WVR.'//i2char(iq), &
         form='unformatted', status='old', access='direct', recl=wv_mreclx)
  end subroutine wv_open_iq_real_for_modify

  subroutine wv_open_iq_imag_for_modify(iq)
    integer, intent(in) :: iq
    character(10) :: i2char
    wv_cur_iq = iq
    if (wv_backend /= WV_BACKEND_FILE) return
    open(newunit=wv_imag_unit, file='__WVI.'//i2char(iq), &
         form='unformatted', status='old', access='direct', recl=wv_mreclx)
  end subroutine wv_open_iq_imag_for_modify

  subroutine wv_modify_get_real(iw, zw)
    integer,     intent(in)  :: iw
    complex(kp), intent(out) :: zw(:,:)
    integer :: n
    if (wv_backend == WV_BACKEND_FILE) then
       read(wv_real_unit, rec=iw - wv_nw_i + 1) zw
    else
       n = size(wv_real_buf_iq1, 1)
       zw(1:n, 1:n) = wv_real_buf_iq1(1:n, 1:n, iw)
    endif
  end subroutine wv_modify_get_real

  subroutine wv_modify_put_real(iw, zw)
    integer,     intent(in) :: iw
    complex(kp), intent(in) :: zw(:,:)
    integer :: n
    if (wv_backend == WV_BACKEND_FILE) then
       write(wv_real_unit, rec=iw - wv_nw_i + 1) zw
    else
       n = size(wv_real_buf_iq1, 1)
       wv_real_buf_iq1(1:n, 1:n, iw) = zw(1:n, 1:n)
    endif
  end subroutine wv_modify_put_real

  subroutine wv_modify_get_imag(iw, zw)
    integer,     intent(in)  :: iw
    complex(kp), intent(out) :: zw(:,:)
    integer :: n
    if (wv_backend == WV_BACKEND_FILE) then
       read(wv_imag_unit, rec=iw) zw
    else
       n = size(wv_imag_buf_iq1, 1)
       zw(1:n, 1:n) = wv_imag_buf_iq1(1:n, 1:n, iw)
    endif
  end subroutine wv_modify_get_imag

  subroutine wv_modify_put_imag(iw, zw)
    integer,     intent(in) :: iw
    complex(kp), intent(in) :: zw(:,:)
    integer :: n
    if (wv_backend == WV_BACKEND_FILE) then
       write(wv_imag_unit, rec=iw) zw
    else
       n = size(wv_imag_buf_iq1, 1)
       wv_imag_buf_iq1(1:n, 1:n, iw) = zw(1:n, 1:n)
    endif
  end subroutine wv_modify_put_imag

  subroutine wv_close_iq_for_modify()
    if (wv_backend /= WV_BACKEND_FILE) return
    if (wv_real_unit > 0) then; close(wv_real_unit); wv_real_unit = -1; endif
    if (wv_imag_unit > 0) then; close(wv_imag_unit); wv_imag_unit = -1; endif
  end subroutine wv_close_iq_for_modify

  ! =============================================================
  ! MEMORY_3D streaming-specific helpers
  ! =============================================================
  subroutine wv_zero_current()
    !> Zero the current (per-iq) buffers. Caller invokes before each iq's WV
    !> computation so non-owner ranks contribute 0 to the subsequent Allreduce.
    if (wv_backend /= WV_BACKEND_MEMORY_3D) return
    wv_real_buf = cmplx(0, 0, kp)
    wv_imag_buf = cmplx(0, 0, kp)
  end subroutine wv_zero_current

  subroutine wv_save_iq1()
    !> Copy the current buffer into the iq=1 saved slot. Called once after the
    !> iq=1 main-loop iteration, before W0w0i runs.
    if (wv_backend /= WV_BACKEND_MEMORY_3D) return
    wv_real_buf_iq1 = wv_real_buf
    wv_imag_buf_iq1 = wv_imag_buf
  end subroutine wv_save_iq1

  subroutine wv_restore_iq1()
    !> Copy iq=1 saved (W0w0i-corrected) data back into the current buffer so
    !> step_kx(kx=1) can read it via wv_get_real / wv_get_imag.
    if (wv_backend /= WV_BACKEND_MEMORY_3D) return
    wv_real_buf = wv_real_buf_iq1
    wv_imag_buf = wv_imag_buf_iq1
  end subroutine wv_restore_iq1

  subroutine wv_sync_current(comm)
    !> Allreduce(SUM) the current 3D buffers across all ranks. Owner has data,
    !> others have 0 (after zero_current), so SUM correctly distributes.
    use mpi
    integer, intent(in) :: comm
    integer :: ierr, mpi_type
    if (wv_backend /= WV_BACKEND_MEMORY_3D) return
#ifdef __MP
    mpi_type = MPI_COMPLEX
#else
    mpi_type = MPI_DOUBLE_COMPLEX
#endif
    call MPI_Allreduce(MPI_IN_PLACE, wv_real_buf, size(wv_real_buf), mpi_type, MPI_SUM, comm, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, wv_imag_buf, size(wv_imag_buf), mpi_type, MPI_SUM, comm, ierr)
  end subroutine wv_sync_current

  subroutine wv_bcast_iq1(root, comm)
    !> Broadcast the iq=1 saved slot (W0w0i-corrected on root) to all ranks.
    use mpi
    integer, intent(in) :: root, comm
    integer :: ierr, mpi_type
    if (wv_backend /= WV_BACKEND_MEMORY_3D) return
#ifdef __MP
    mpi_type = MPI_COMPLEX
#else
    mpi_type = MPI_DOUBLE_COMPLEX
#endif
    call MPI_Bcast(wv_real_buf_iq1, size(wv_real_buf_iq1), mpi_type, root, comm, ierr)
    call MPI_Bcast(wv_imag_buf_iq1, size(wv_imag_buf_iq1), mpi_type, root, comm, ierr)
  end subroutine wv_bcast_iq1

end module m_wv_storage
