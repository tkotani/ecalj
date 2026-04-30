!> m_wv_storage — handle for the W-V (W minus v) data movement between phases.
!!
!! The W-V matrices are produced per-iq by m_llw:WVRllwR / WVIllwI, modified at
!! iq=1 by m_w0w0i:modifyWV0 (Gamma-cell W(0) correction), and consumed per-kx by
!! m_sxcf_sc (correlation self-energy). Three modules share access; the data may
!! live on disk (__WVR.<iq> / __WVI.<iq> files) or in memory (Phase 1-B 4D buffer
!! / Phase 1-C 3D streaming buffer).
!!
!! This module wraps the storage as a derived type so callers do not branch on
!! backend. The abstraction is grown in steps:
!!   - First step (this commit): FILE backend only. Replaces direct calls to
!!     openm / writem / open / read / close in m_llw, m_w0w0i, m_sxcf_sc.
!!   - Later step: MEMORY_3D streaming backend (Phase 1-C). Adds save/restore
!!     iq=1 slot, MPI sync helpers.
!!
!! API conventions:
!!   - call wv_init_file(self, mreclx, comm, nw_i) before any access.
!!   - Writer side per iq:
!!       wv_open_iq_real_for_write(self, iq);   wv_open_iq_imag_for_write(self, iq)
!!       wv_put_real(self, iw, zw);             wv_put_imag(self, iw, zw)
!!       wv_close_iq_for_write(self)
!!   - Reader side per iq:
!!       wv_open_iq_for_read(self, iq)
!!       wv_get_real(self, iw, zw_out);          wv_get_imag(self, iw, zw_out)
!!       wv_close_iq_for_read(self)
!!   - call wv_dealloc(self) at the end (no-op for FILE backend).
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

  type, public :: wv_storage
    integer :: backend = WV_BACKEND_FILE
    !> FILE backend
    integer :: mreclx = 0
    integer :: comm   = -1
    integer :: nw_i   = 0
    integer :: cur_iq = 0
    integer :: ifrcw_unit  = -1
    integer :: ifrcwi_unit = -1
    !> MEMORY_3D backend: current iq buffer + saved iq=1 slot.
    !>   Shape: real_buf  (nblochpmx, nblochpmx, nw_i:nw)
    !>          imag_buf  (nblochpmx, nblochpmx, 1:niw)
    !> The current buffer is overwritten as iq advances; the iq1 slot is
    !> saved by save_iq1 before W0w0i and restored after Bcast.
    complex(kp), allocatable :: real_buf(:,:,:)
    complex(kp), allocatable :: imag_buf(:,:,:)
    complex(kp), allocatable :: real_buf_iq1(:,:,:)
    complex(kp), allocatable :: imag_buf_iq1(:,:,:)
  end type wv_storage

  public :: wv_init_file, wv_dealloc
  public :: wv_open_iq_real_for_write, wv_open_iq_imag_for_write
  public :: wv_close_iq_for_write
  public :: wv_open_iq_for_read, wv_close_iq_for_read
  public :: wv_put_real, wv_put_imag, wv_get_real, wv_get_imag
  public :: wv_open_iq_real_for_modify, wv_open_iq_imag_for_modify
  public :: wv_modify_get_real, wv_modify_put_real
  public :: wv_modify_get_imag, wv_modify_put_imag
  public :: wv_close_iq_for_modify
  ! MEMORY_3D streaming backend
  public :: wv_init_memory_3d
  public :: wv_zero_current
  public :: wv_save_iq1, wv_restore_iq1
  public :: wv_sync_current, wv_bcast_iq1

contains

  subroutine wv_init_file(self, mreclx, comm, nw_i)
    !> Configure storage for FILE backend. Caller-side state only; no I/O yet.
    type(wv_storage), intent(out) :: self
    integer, intent(in) :: mreclx, comm, nw_i
    self%backend = WV_BACKEND_FILE
    self%mreclx = mreclx
    self%comm = comm
    self%nw_i = nw_i
    self%cur_iq = 0
    self%ifrcw_unit = -1
    self%ifrcwi_unit = -1
  end subroutine wv_init_file

  subroutine wv_dealloc(self)
    !> Close any lingering file units (FILE) and free 3D buffers (MEMORY_3D).
    type(wv_storage), intent(inout) :: self
    integer :: istat
    if (self%backend == WV_BACKEND_FILE) then
       if (self%ifrcw_unit  > 0) then; istat = closem(self%ifrcw_unit);  self%ifrcw_unit  = -1; endif
       if (self%ifrcwi_unit > 0) then; istat = closem(self%ifrcwi_unit); self%ifrcwi_unit = -1; endif
    else if (self%backend == WV_BACKEND_MEMORY_3D) then
       if (allocated(self%real_buf))     deallocate(self%real_buf)
       if (allocated(self%imag_buf))     deallocate(self%imag_buf)
       if (allocated(self%real_buf_iq1)) deallocate(self%real_buf_iq1)
       if (allocated(self%imag_buf_iq1)) deallocate(self%imag_buf_iq1)
    endif
  end subroutine wv_dealloc

  subroutine wv_open_iq_real_for_write(self, iq)
    !> FILE: open WVR.<iq> via openm. MEMORY_3D: track current iq only.
    type(wv_storage), intent(inout) :: self
    integer, intent(in) :: iq
    integer :: istat
    character(10) :: i2char
    self%cur_iq = iq
    if (self%backend /= WV_BACKEND_FILE) return
    istat = openm(newunit=self%ifrcw_unit, file='__WVR.'//i2char(iq), &
                  recl=self%mreclx, comm=self%comm)
  end subroutine wv_open_iq_real_for_write

  subroutine wv_open_iq_imag_for_write(self, iq)
    !> FILE: open WVI.<iq> via openm. MEMORY_3D: track current iq only.
    type(wv_storage), intent(inout) :: self
    integer, intent(in) :: iq
    integer :: istat
    character(10) :: i2char
    self%cur_iq = iq
    if (self%backend /= WV_BACKEND_FILE) return
    istat = openm(newunit=self%ifrcwi_unit, file='__WVI.'//i2char(iq), &
                  recl=self%mreclx, comm=self%comm)
  end subroutine wv_open_iq_imag_for_write

  subroutine wv_close_iq_for_write(self)
    !> FILE: close any open units. MEMORY_3D: no-op (data remains in buffer).
    type(wv_storage), intent(inout) :: self
    integer :: istat
    if (self%backend /= WV_BACKEND_FILE) return
    if (self%ifrcw_unit  > 0) then; istat = closem(self%ifrcw_unit);  self%ifrcw_unit  = -1; endif
    if (self%ifrcwi_unit > 0) then; istat = closem(self%ifrcwi_unit); self%ifrcwi_unit = -1; endif
  end subroutine wv_close_iq_for_write

  subroutine wv_put_real(self, iw, zw)
    !> FILE: writem record. MEMORY_3D: store into current real buffer slot.
    type(wv_storage), intent(inout) :: self
    integer, intent(in) :: iw
    complex(kp), intent(in) :: zw(:,:)
    integer :: istat, n
    if (self%backend == WV_BACKEND_FILE) then
       istat = writem(self%ifrcw_unit, rec=iw - self%nw_i + 1, data=zw)
    else
       n = size(self%real_buf, 1)
       self%real_buf(1:n, 1:n, iw) = zw(1:n, 1:n)
    endif
  end subroutine wv_put_real

  subroutine wv_put_imag(self, iw, zw)
    !> FILE: writem record. MEMORY_3D: store into current imag buffer slot.
    type(wv_storage), intent(inout) :: self
    integer, intent(in) :: iw
    complex(kp), intent(in) :: zw(:,:)
    integer :: istat, n
    if (self%backend == WV_BACKEND_FILE) then
       istat = writem(self%ifrcwi_unit, rec=iw, data=zw)
    else
       n = size(self%imag_buf, 1)
       self%imag_buf(1:n, 1:n, iw) = zw(1:n, 1:n)
    endif
  end subroutine wv_put_imag

  subroutine wv_open_iq_for_read(self, iq, want_real, want_imag)
    !> FILE: open the WV files. MEMORY_3D: track current iq only — caller is
    !> responsible for ensuring the buffer contains the requested iq's data
    !> (via wv_zero_current + wv_put_* + wv_sync_current upstream).
    type(wv_storage), intent(inout) :: self
    integer, intent(in) :: iq
    logical, intent(in) :: want_real, want_imag
    character(10) :: i2char
    self%cur_iq = iq
    if (self%backend /= WV_BACKEND_FILE) return
    if (want_real) then
       open(newunit=self%ifrcw_unit, file='__WVR.'//i2char(iq), action='read', &
            form='unformatted', access='direct', recl=self%mreclx)
    endif
    if (want_imag) then
       open(newunit=self%ifrcwi_unit, file='__WVI.'//i2char(iq), action='read', &
            form='unformatted', access='direct', recl=self%mreclx)
    endif
  end subroutine wv_open_iq_for_read

  subroutine wv_close_iq_for_read(self)
    type(wv_storage), intent(inout) :: self
    if (self%backend /= WV_BACKEND_FILE) return
    if (self%ifrcw_unit  > 0) then; close(self%ifrcw_unit);  self%ifrcw_unit  = -1; endif
    if (self%ifrcwi_unit > 0) then; close(self%ifrcwi_unit); self%ifrcwi_unit = -1; endif
  end subroutine wv_close_iq_for_read

  subroutine wv_get_real(self, iw, zw)
    !> FILE: read record from WVR. MEMORY_3D: read from current real buffer.
    type(wv_storage), intent(in)  :: self
    integer,          intent(in)  :: iw
    complex(kp),      intent(out) :: zw(:,:)
    integer :: n
    if (self%backend == WV_BACKEND_FILE) then
       read(self%ifrcw_unit, rec=iw - self%nw_i + 1) zw
    else
       n = size(self%real_buf, 1)
       zw(1:n, 1:n) = self%real_buf(1:n, 1:n, iw)
    endif
  end subroutine wv_get_real

  subroutine wv_get_imag(self, iw, zw)
    !> FILE: read record from WVI. MEMORY_3D: read from current imag buffer.
    type(wv_storage), intent(in)  :: self
    integer,          intent(in)  :: iw
    complex(kp),      intent(out) :: zw(:,:)
    integer :: n
    if (self%backend == WV_BACKEND_FILE) then
       read(self%ifrcwi_unit, rec=iw) zw
    else
       n = size(self%imag_buf, 1)
       zw(1:n, 1:n) = self%imag_buf(1:n, 1:n, iw)
    endif
  end subroutine wv_get_imag

  ! ---- modify mode (rank-0 read+write, no MPI coordination) ----
  ! FILE: read+write on direct-access unit. MEMORY_3D: read+write on iq=1 saved slot.
  subroutine wv_open_iq_real_for_modify(self, iq)
    type(wv_storage), intent(inout) :: self
    integer, intent(in) :: iq
    character(10) :: i2char
    self%cur_iq = iq
    if (self%backend /= WV_BACKEND_FILE) return
    open(newunit=self%ifrcw_unit, file='__WVR.'//i2char(iq), &
         form='unformatted', status='old', access='direct', recl=self%mreclx)
  end subroutine wv_open_iq_real_for_modify

  subroutine wv_open_iq_imag_for_modify(self, iq)
    type(wv_storage), intent(inout) :: self
    integer, intent(in) :: iq
    character(10) :: i2char
    self%cur_iq = iq
    if (self%backend /= WV_BACKEND_FILE) return
    open(newunit=self%ifrcwi_unit, file='__WVI.'//i2char(iq), &
         form='unformatted', status='old', access='direct', recl=self%mreclx)
  end subroutine wv_open_iq_imag_for_modify

  subroutine wv_modify_get_real(self, iw, zw)
    type(wv_storage), intent(in) :: self
    integer, intent(in) :: iw
    complex(kp), intent(out) :: zw(:,:)
    integer :: n
    if (self%backend == WV_BACKEND_FILE) then
       read(self%ifrcw_unit, rec=iw - self%nw_i + 1) zw
    else
       n = size(self%real_buf_iq1, 1)
       zw(1:n, 1:n) = self%real_buf_iq1(1:n, 1:n, iw)
    endif
  end subroutine wv_modify_get_real

  subroutine wv_modify_put_real(self, iw, zw)
    type(wv_storage), intent(inout) :: self
    integer, intent(in) :: iw
    complex(kp), intent(in) :: zw(:,:)
    integer :: n
    if (self%backend == WV_BACKEND_FILE) then
       write(self%ifrcw_unit, rec=iw - self%nw_i + 1) zw
    else
       n = size(self%real_buf_iq1, 1)
       self%real_buf_iq1(1:n, 1:n, iw) = zw(1:n, 1:n)
    endif
  end subroutine wv_modify_put_real

  subroutine wv_modify_get_imag(self, iw, zw)
    type(wv_storage), intent(in) :: self
    integer, intent(in) :: iw
    complex(kp), intent(out) :: zw(:,:)
    integer :: n
    if (self%backend == WV_BACKEND_FILE) then
       read(self%ifrcwi_unit, rec=iw) zw
    else
       n = size(self%imag_buf_iq1, 1)
       zw(1:n, 1:n) = self%imag_buf_iq1(1:n, 1:n, iw)
    endif
  end subroutine wv_modify_get_imag

  subroutine wv_modify_put_imag(self, iw, zw)
    type(wv_storage), intent(inout) :: self
    integer, intent(in) :: iw
    complex(kp), intent(in) :: zw(:,:)
    integer :: n
    if (self%backend == WV_BACKEND_FILE) then
       write(self%ifrcwi_unit, rec=iw) zw
    else
       n = size(self%imag_buf_iq1, 1)
       self%imag_buf_iq1(1:n, 1:n, iw) = zw(1:n, 1:n)
    endif
  end subroutine wv_modify_put_imag

  subroutine wv_close_iq_for_modify(self)
    type(wv_storage), intent(inout) :: self
    if (self%backend /= WV_BACKEND_FILE) return
    if (self%ifrcw_unit  > 0) then; close(self%ifrcw_unit);  self%ifrcw_unit  = -1; endif
    if (self%ifrcwi_unit > 0) then; close(self%ifrcwi_unit); self%ifrcwi_unit = -1; endif
  end subroutine wv_close_iq_for_modify

  ! =============================================================
  ! MEMORY_3D streaming-specific methods
  ! =============================================================
  subroutine wv_init_memory_3d(self, nbpmx, nw_lo, nw_hi, niwx)
    !> Configure storage for MEMORY_3D streaming. Allocates 3D current buffer
    !> and 3D iq=1 saved buffer, both initialized to zero.
    type(wv_storage), intent(out) :: self
    integer, intent(in) :: nbpmx, nw_lo, nw_hi, niwx
    self%backend = WV_BACKEND_MEMORY_3D
    self%nw_i = nw_lo
    allocate(self%real_buf(nbpmx, nbpmx, nw_lo:nw_hi),     source=cmplx(0,0,kp))
    allocate(self%imag_buf(nbpmx, nbpmx, 1:niwx),          source=cmplx(0,0,kp))
    allocate(self%real_buf_iq1(nbpmx, nbpmx, nw_lo:nw_hi), source=cmplx(0,0,kp))
    allocate(self%imag_buf_iq1(nbpmx, nbpmx, 1:niwx),      source=cmplx(0,0,kp))
  end subroutine wv_init_memory_3d

  subroutine wv_zero_current(self)
    !> Zero the current (per-iq) buffers. Caller invokes before each iq's
    !> WV computation so non-owner ranks contribute 0 to the subsequent Allreduce.
    type(wv_storage), intent(inout) :: self
    if (self%backend /= WV_BACKEND_MEMORY_3D) return
    self%real_buf = cmplx(0, 0, kp)
    self%imag_buf = cmplx(0, 0, kp)
  end subroutine wv_zero_current

  subroutine wv_save_iq1(self)
    !> Copy the current buffer into the iq=1 saved slot. Called once after
    !> the iq=1 main-loop iteration, before W0w0i runs.
    type(wv_storage), intent(inout) :: self
    if (self%backend /= WV_BACKEND_MEMORY_3D) return
    self%real_buf_iq1 = self%real_buf
    self%imag_buf_iq1 = self%imag_buf
  end subroutine wv_save_iq1

  subroutine wv_restore_iq1(self)
    !> Copy iq=1 saved (W0w0i-corrected) data back into the current buffer
    !> so that step_kx(kx=1) can read it via wv_get_real / wv_get_imag.
    type(wv_storage), intent(inout) :: self
    if (self%backend /= WV_BACKEND_MEMORY_3D) return
    self%real_buf = self%real_buf_iq1
    self%imag_buf = self%imag_buf_iq1
  end subroutine wv_restore_iq1

  subroutine wv_sync_current(self, comm)
    !> Allreduce(SUM) the current 3D buffers across all ranks. Owner has data,
    !> others have 0 (after zero_current), so SUM correctly distributes.
    use mpi
    type(wv_storage), intent(inout) :: self
    integer, intent(in) :: comm
    integer :: ierr, mpi_type
    if (self%backend /= WV_BACKEND_MEMORY_3D) return
#ifdef __MP
    mpi_type = MPI_COMPLEX
#else
    mpi_type = MPI_DOUBLE_COMPLEX
#endif
    call MPI_Allreduce(MPI_IN_PLACE, self%real_buf, size(self%real_buf), mpi_type, MPI_SUM, comm, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, self%imag_buf, size(self%imag_buf), mpi_type, MPI_SUM, comm, ierr)
  end subroutine wv_sync_current

  subroutine wv_bcast_iq1(self, root, comm)
    !> Broadcast the iq=1 saved slot (W0w0i-corrected on root) to all ranks.
    use mpi
    type(wv_storage), intent(inout) :: self
    integer, intent(in) :: root, comm
    integer :: ierr, mpi_type
    if (self%backend /= WV_BACKEND_MEMORY_3D) return
#ifdef __MP
    mpi_type = MPI_COMPLEX
#else
    mpi_type = MPI_DOUBLE_COMPLEX
#endif
    call MPI_Bcast(self%real_buf_iq1, size(self%real_buf_iq1), mpi_type, root, comm, ierr)
    call MPI_Bcast(self%imag_buf_iq1, size(self%imag_buf_iq1), mpi_type, root, comm, ierr)
  end subroutine wv_bcast_iq1

end module m_wv_storage
