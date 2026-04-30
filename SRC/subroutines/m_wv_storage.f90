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

  integer, parameter, public :: WV_BACKEND_FILE = 1
  ! Reserved for the next step:
  ! integer, parameter, public :: WV_BACKEND_MEMORY_3D = 2

  type, public :: wv_storage
    integer :: backend = WV_BACKEND_FILE
    integer :: mreclx = 0     !> bytes-per-record for the underlying file
    integer :: comm   = -1    !> MPI communicator for openm coordination (writer side)
    integer :: nw_i   = 0     !> real-axis lower bound; record offset = iw - nw_i + 1
    !> per-iq state. For the writer the units are obtained from openm; for the
    !> reader they are plain Fortran direct-access units.
    integer :: cur_iq = 0
    integer :: ifrcw_unit  = -1  !> real-axis (WVR.<iq>) unit
    integer :: ifrcwi_unit = -1  !> imag-axis (WVI.<iq>) unit
  end type wv_storage

  public :: wv_init_file, wv_dealloc
  public :: wv_open_iq_real_for_write, wv_open_iq_imag_for_write
  public :: wv_close_iq_for_write
  public :: wv_open_iq_for_read, wv_close_iq_for_read
  public :: wv_put_real, wv_put_imag, wv_get_real, wv_get_imag
  ! Modify mode: rank 0 only (e.g. m_w0w0i:modifyWV0). Plain Fortran read/write
  ! on a status='old' unit, no MPI coordination.
  public :: wv_open_iq_real_for_modify, wv_open_iq_imag_for_modify
  public :: wv_modify_get_real, wv_modify_put_real
  public :: wv_modify_get_imag, wv_modify_put_imag
  public :: wv_close_iq_for_modify

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
    !> Close any lingering file units. For FILE backend usually a no-op (caller
    !> is expected to close per-iq), but defensive.
    type(wv_storage), intent(inout) :: self
    integer :: istat
    if (self%backend /= WV_BACKEND_FILE) return
    if (self%ifrcw_unit  > 0) then; istat = closem(self%ifrcw_unit);  self%ifrcw_unit  = -1; endif
    if (self%ifrcwi_unit > 0) then; istat = closem(self%ifrcwi_unit); self%ifrcwi_unit = -1; endif
  end subroutine wv_dealloc

  subroutine wv_open_iq_real_for_write(self, iq)
    !> Open the WVR.<iq> file for MPI-coordinated writes (openm).
    type(wv_storage), intent(inout) :: self
    integer, intent(in) :: iq
    integer :: istat
    character(10) :: i2char
    self%cur_iq = iq
    istat = openm(newunit=self%ifrcw_unit, file='__WVR.'//i2char(iq), &
                  recl=self%mreclx, comm=self%comm)
  end subroutine wv_open_iq_real_for_write

  subroutine wv_open_iq_imag_for_write(self, iq)
    !> Open the WVI.<iq> file for MPI-coordinated writes (openm).
    type(wv_storage), intent(inout) :: self
    integer, intent(in) :: iq
    integer :: istat
    character(10) :: i2char
    self%cur_iq = iq
    istat = openm(newunit=self%ifrcwi_unit, file='__WVI.'//i2char(iq), &
                  recl=self%mreclx, comm=self%comm)
  end subroutine wv_open_iq_imag_for_write

  subroutine wv_close_iq_for_write(self)
    !> Close any write-side units that were opened.
    type(wv_storage), intent(inout) :: self
    integer :: istat
    if (self%ifrcw_unit  > 0) then; istat = closem(self%ifrcw_unit);  self%ifrcw_unit  = -1; endif
    if (self%ifrcwi_unit > 0) then; istat = closem(self%ifrcwi_unit); self%ifrcwi_unit = -1; endif
  end subroutine wv_close_iq_for_write

  subroutine wv_put_real(self, iw, zw)
    !> Write zw at real-axis index iw (record number = iw - nw_i + 1).
    type(wv_storage), intent(inout) :: self
    integer, intent(in) :: iw
    complex(kp), intent(in) :: zw(:,:)
    integer :: istat
    istat = writem(self%ifrcw_unit, rec=iw - self%nw_i + 1, data=zw)
  end subroutine wv_put_real

  subroutine wv_put_imag(self, iw, zw)
    !> Write zw at imag-axis index iw (record number = iw, 1-based).
    type(wv_storage), intent(inout) :: self
    integer, intent(in) :: iw
    complex(kp), intent(in) :: zw(:,:)
    integer :: istat
    istat = writem(self%ifrcwi_unit, rec=iw, data=zw)
  end subroutine wv_put_imag

  subroutine wv_open_iq_for_read(self, iq, want_real, want_imag)
    !> Open the WVR.<iq> and/or WVI.<iq> files for direct-access reads.
    !> The reader path uses ordinary Fortran open (no MPI coordination).
    type(wv_storage), intent(inout) :: self
    integer, intent(in) :: iq
    logical, intent(in) :: want_real, want_imag
    character(10) :: i2char
    self%cur_iq = iq
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
    if (self%ifrcw_unit  > 0) then; close(self%ifrcw_unit);  self%ifrcw_unit  = -1; endif
    if (self%ifrcwi_unit > 0) then; close(self%ifrcwi_unit); self%ifrcwi_unit = -1; endif
  end subroutine wv_close_iq_for_read

  subroutine wv_get_real(self, iw, zw)
    type(wv_storage), intent(in)  :: self
    integer,          intent(in)  :: iw
    complex(kp),      intent(out) :: zw(:,:)
    read(self%ifrcw_unit, rec=iw - self%nw_i + 1) zw
  end subroutine wv_get_real

  subroutine wv_get_imag(self, iw, zw)
    type(wv_storage), intent(in)  :: self
    integer,          intent(in)  :: iw
    complex(kp),      intent(out) :: zw(:,:)
    read(self%ifrcwi_unit, rec=iw) zw
  end subroutine wv_get_imag

  ! ---- modify mode (rank-0 read+write, no MPI coordination) ----
  subroutine wv_open_iq_real_for_modify(self, iq)
    type(wv_storage), intent(inout) :: self
    integer, intent(in) :: iq
    character(10) :: i2char
    self%cur_iq = iq
    open(newunit=self%ifrcw_unit, file='__WVR.'//i2char(iq), &
         form='unformatted', status='old', access='direct', recl=self%mreclx)
  end subroutine wv_open_iq_real_for_modify

  subroutine wv_open_iq_imag_for_modify(self, iq)
    type(wv_storage), intent(inout) :: self
    integer, intent(in) :: iq
    character(10) :: i2char
    self%cur_iq = iq
    open(newunit=self%ifrcwi_unit, file='__WVI.'//i2char(iq), &
         form='unformatted', status='old', access='direct', recl=self%mreclx)
  end subroutine wv_open_iq_imag_for_modify

  subroutine wv_modify_get_real(self, iw, zw)
    type(wv_storage), intent(in) :: self
    integer, intent(in) :: iw
    complex(kp), intent(out) :: zw(:,:)
    read(self%ifrcw_unit, rec=iw - self%nw_i + 1) zw
  end subroutine wv_modify_get_real

  subroutine wv_modify_put_real(self, iw, zw)
    type(wv_storage), intent(in) :: self
    integer, intent(in) :: iw
    complex(kp), intent(in) :: zw(:,:)
    write(self%ifrcw_unit, rec=iw - self%nw_i + 1) zw
  end subroutine wv_modify_put_real

  subroutine wv_modify_get_imag(self, iw, zw)
    type(wv_storage), intent(in) :: self
    integer, intent(in) :: iw
    complex(kp), intent(out) :: zw(:,:)
    read(self%ifrcwi_unit, rec=iw) zw
  end subroutine wv_modify_get_imag

  subroutine wv_modify_put_imag(self, iw, zw)
    type(wv_storage), intent(in) :: self
    integer, intent(in) :: iw
    complex(kp), intent(in) :: zw(:,:)
    write(self%ifrcwi_unit, rec=iw) zw
  end subroutine wv_modify_put_imag

  subroutine wv_close_iq_for_modify(self)
    type(wv_storage), intent(inout) :: self
    if (self%ifrcw_unit  > 0) then; close(self%ifrcw_unit);  self%ifrcw_unit  = -1; endif
    if (self%ifrcwi_unit > 0) then; close(self%ifrcwi_unit); self%ifrcwi_unit = -1; endif
  end subroutine wv_close_iq_for_modify

end module m_wv_storage
