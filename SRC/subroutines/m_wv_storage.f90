!> m_wv_storage — handle for the W-V (W minus v) data movement between phases.
!!
!! The W-V matrices are produced per-iq by m_llw:WVRllwR / WVIllwI, modified at
!! iq=1 by m_w0w0i:modifyWV0 (Gamma-cell W(0) correction), and consumed per-kx by
!! m_sxcf_sc (correlation self-energy). Three modules share access; the data may
!! live on disk (__WVR.<iq> / __WVI.<iq> files) or in MPI shared-memory windows.
!!
!! Per ecalj convention this module exposes a singleton: subroutines take inputs
!! by argument and operate on module-level state ("output" via module variables).
!! Backends:
!!   WV_BACKEND_FILE  openm/writem/read/close on __WVR.<iq>/__WVI.<iq>
!!   WV_BACKEND_SHM   MPI shared memory window on comm_q node
!!
!! API conventions:
!!   - call wv_init_file(mreclx, nw_i)
!!   - SHM: wv_init_shm(ngbx, niwx, rcxq_lo, rcxq_hi, comm) called per-iq inside
!!       build_screened_coulomb_step_kx (after ngb is known from Readvcoud).
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

  integer, parameter :: WV_BACKEND_FILE = 1
  integer, parameter :: WV_BACKEND_SHM  = 2  !> MPI shared memory window on comm_q node

  ! ---- module-level singleton state ----
  integer, protected :: wv_backend = WV_BACKEND_FILE
  integer, protected :: wv_cur_iq  = 0
  integer :: wv_mreclx = 0
  integer :: wv_nw_i   = 0
  integer :: wv_real_unit = -1
  integer :: wv_imag_unit = -1
  logical :: wv_in_modify_mode = .false.  ! FILE: use standard write (not MPI-IO) in wv_put_*
  ! SHM backend: Wc-only buffers in shared memory (one copy per node via MPI_Win_allocate_shared).
  ! shm_wvr holds W-V real axis (Wc_real); shm_wvi holds W-V imag axis (Wc_imag).
  ! rcxq (chi0 spectral weight) stays private in x0kf_zxq; only the Wc output goes here.
  integer :: wv_shm_win_real = -1
  integer :: wv_shm_win_imag = -1
  integer, public :: wv_ngb = 0
  complex(kp), public, pointer, contiguous :: shm_wvr(:,:,:) => null()
  complex(kp), public, pointer, contiguous :: shm_wvi(:,:,:) => null()

  public :: wv_init_file, wv_init_shm, wv_dealloc, wv_dump_shm_to_file
  public :: wv_open_iq_real_for_write, wv_open_iq_imag_for_write
  public :: wv_close_iq_for_write
  public :: wv_open_iq_for_read, wv_close_iq_for_read
  public :: wv_put_real, wv_put_imag, wv_get_real, wv_get_imag
  public :: wv_open_iq_real_for_modify, wv_open_iq_imag_for_modify
  public :: wv_close_iq_for_modify

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

  subroutine wv_init_shm(ngbx, niwx, rcxq_lo, rcxq_hi, comm, mreclx)
    !> Configure singleton for SHM backend.
    !> shm_wvr(ngbx, ngbx, rcxq_lo:rcxq_hi) is the shared rcxq:
    !>   phase 1 — spectral bins (x0 accumulation, k-loop)
    !>   phase 2 — chi0_real (after KK, dpsion_chiq in-place)
    !>   phase 3 — Wc real   (after Dyson, WVRllwR in-place)
    !> shm_wvi(ngbx, ngbx, niwx) mirrors zxqi (chi0_imag → Wc_imag).
    !> rcxq_lo = (1-npm)*nwhis,  rcxq_hi = nwhis.
    !> Rank 0 in comm contributes the allocation; all ranks get a Fortran pointer
    !> to the same physical memory via MPI_Win_shared_query.
    use iso_c_binding, only: c_ptr, c_f_pointer
    use mpi
    integer, intent(in) :: ngbx, niwx, rcxq_lo, rcxq_hi, comm
    integer, intent(in), optional :: mreclx
    integer(MPI_ADDRESS_KIND) :: nbytes, sz
    integer :: disp_unit, ierr, my_rank, nrcxq
    type(c_ptr) :: baseptr
    complex(kp) :: tmp_c

    wv_backend = WV_BACKEND_SHM
    wv_nw_i    = rcxq_lo
    if (present(mreclx)) wv_mreclx = mreclx   ! lower bound of shm_wvr 3rd dim; wv_put/get_real uses iw-wv_nw_i+1
    wv_ngb     = ngbx
    wv_cur_iq  = 0
    nrcxq      = rcxq_hi - rcxq_lo + 1
    call MPI_Comm_rank(comm, my_rank, ierr)
    disp_unit = 1

    if (wv_shm_win_real /= -1) call MPI_Win_free(wv_shm_win_real, ierr)
    nbytes = merge(int(ngbx,MPI_ADDRESS_KIND)**2 * int(nrcxq,MPI_ADDRESS_KIND) &
                   * int(storage_size(tmp_c)/8,MPI_ADDRESS_KIND), &
                   0_MPI_ADDRESS_KIND, my_rank == 0)
    call MPI_Win_allocate_shared(nbytes, disp_unit, MPI_INFO_NULL, comm, baseptr, wv_shm_win_real, ierr)
    call MPI_Win_shared_query(wv_shm_win_real, 0, sz, disp_unit, baseptr, ierr)
    call c_f_pointer(baseptr, shm_wvr, [ngbx, ngbx, nrcxq])

    if (wv_shm_win_imag /= -1) call MPI_Win_free(wv_shm_win_imag, ierr)
    nbytes = merge(int(ngbx,MPI_ADDRESS_KIND)**2 * int(niwx,MPI_ADDRESS_KIND) &
                   * int(storage_size(tmp_c)/8,MPI_ADDRESS_KIND), &
                   0_MPI_ADDRESS_KIND, my_rank == 0)
    call MPI_Win_allocate_shared(nbytes, disp_unit, MPI_INFO_NULL, comm, baseptr, wv_shm_win_imag, ierr)
    call MPI_Win_shared_query(wv_shm_win_imag, 0, sz, disp_unit, baseptr, ierr)
    call c_f_pointer(baseptr, shm_wvi, [ngbx, ngbx, niwx])
  end subroutine wv_init_shm

  subroutine wv_dealloc()
    !> Close any lingering file units (FILE) or release SHM windows.
    use mpi
    integer :: istat, ierr
    if (wv_backend == WV_BACKEND_FILE) then
       if (wv_real_unit > 0) then; istat = closem(wv_real_unit); wv_real_unit = -1; endif
       if (wv_imag_unit > 0) then; istat = closem(wv_imag_unit); wv_imag_unit = -1; endif
    else if (wv_backend == WV_BACKEND_SHM) then
       nullify(shm_wvr)
       nullify(shm_wvi)
       if (wv_shm_win_real /= -1) then
         call MPI_Win_free(wv_shm_win_real, ierr); wv_shm_win_real = -1
       endif
       if (wv_shm_win_imag /= -1) then
         call MPI_Win_free(wv_shm_win_imag, ierr); wv_shm_win_imag = -1
       endif
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
    integer :: istat
    if (wv_backend == WV_BACKEND_SHM .and. associated(shm_wvr)) then
       shm_wvr(1:wv_ngb, 1:wv_ngb, iw - wv_nw_i + 1) = zw(1:wv_ngb, 1:wv_ngb)
    elseif (wv_in_modify_mode) then
       write(wv_real_unit, rec=iw - wv_nw_i + 1) zw  ! rank-0 only, standard Fortran write
    else
       istat = writem(wv_real_unit, rec=iw - wv_nw_i + 1, data=zw)
    endif
  end subroutine wv_put_real

  subroutine wv_put_imag(iw, zw)
    integer,     intent(in) :: iw
    complex(kp), intent(in) :: zw(:,:)
    integer :: istat
    if (wv_backend == WV_BACKEND_SHM .and. associated(shm_wvi)) then
       shm_wvi(1:wv_ngb, 1:wv_ngb, iw) = zw(1:wv_ngb, 1:wv_ngb)
    elseif (wv_in_modify_mode) then
       write(wv_imag_unit, rec=iw) zw  ! rank-0 only, standard Fortran write
    else
       istat = writem(wv_imag_unit, rec=iw, data=zw)
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
    if (wv_backend == WV_BACKEND_SHM .and. associated(shm_wvr)) then
       zw(1:wv_ngb, 1:wv_ngb) = shm_wvr(1:wv_ngb, 1:wv_ngb, iw - wv_nw_i + 1)
    else
       read(wv_real_unit, rec=iw - wv_nw_i + 1) zw
    endif
  end subroutine wv_get_real

  subroutine wv_get_imag(iw, zw)
    integer,     intent(in)  :: iw
    complex(kp), intent(out) :: zw(:,:)
    if (wv_backend == WV_BACKEND_SHM .and. associated(shm_wvi)) then
       zw(1:wv_ngb, 1:wv_ngb) = shm_wvi(1:wv_ngb, 1:wv_ngb, iw)
    else
       read(wv_imag_unit, rec=iw) zw
    endif
  end subroutine wv_get_imag

  ! ---- modify mode (rank-0 read+write, no MPI coordination) ----
  subroutine wv_open_iq_real_for_modify(iq)
    integer, intent(in) :: iq
    character(10) :: i2char
    wv_cur_iq = iq
    wv_in_modify_mode = .true.
    if (wv_backend == WV_BACKEND_SHM .and. associated(shm_wvr)) return  ! use SHM directly
    open(newunit=wv_real_unit, file='__WVR.'//i2char(iq), &
         form='unformatted', status='old', access='direct', recl=wv_mreclx)
  end subroutine wv_open_iq_real_for_modify

  subroutine wv_open_iq_imag_for_modify(iq)
    integer, intent(in) :: iq
    character(10) :: i2char
    wv_cur_iq = iq
    wv_in_modify_mode = .true.
    if (wv_backend == WV_BACKEND_SHM .and. associated(shm_wvi)) return  ! use SHM directly
    open(newunit=wv_imag_unit, file='__WVI.'//i2char(iq), &
         form='unformatted', status='old', access='direct', recl=wv_mreclx)
  end subroutine wv_open_iq_imag_for_modify

  subroutine wv_close_iq_for_modify()
    wv_in_modify_mode = .false.
    if (wv_real_unit > 0) then; close(wv_real_unit); wv_real_unit = -1; endif
    if (wv_imag_unit > 0) then; close(wv_imag_unit); wv_imag_unit = -1; endif
  end subroutine wv_close_iq_for_modify

  subroutine wv_dump_shm_to_file(iq, mreclx, nblochpmx_in, want_real, want_imag)
    !> Dump SHM backend buffers shm_wvr/shm_wvi to __WVR.<iq>/__WVI.<iq> files.
    !> Must be called by rank 0 only; records match FILE-mode wv_put_real format:
    !> nblochpmx_in × nblochpmx_in complex(kp) per record, ngb block in top-left corner.
    integer, intent(in) :: iq, mreclx, nblochpmx_in
    logical, intent(in) :: want_real, want_imag
    integer :: irec, fu
    character(10) :: i2char
    complex(kp) :: buf(nblochpmx_in, nblochpmx_in)
    if (want_real) then
      open(newunit=fu, file='__WVR.'//i2char(iq), form='unformatted', &
           status='replace', access='direct', recl=mreclx)
      buf = (0_kp, 0_kp)
      do irec = 1, size(shm_wvr, 3)
        buf(1:wv_ngb, 1:wv_ngb) = shm_wvr(1:wv_ngb, 1:wv_ngb, irec)
        write(fu, rec=irec) buf
      enddo
      close(fu)
    endif
    if (want_imag) then
      open(newunit=fu, file='__WVI.'//i2char(iq), form='unformatted', &
           status='replace', access='direct', recl=mreclx)
      buf = (0_kp, 0_kp)
      do irec = 1, size(shm_wvi, 3)
        buf(1:wv_ngb, 1:wv_ngb) = shm_wvi(1:wv_ngb, 1:wv_ngb, irec)
        write(fu, rec=irec) buf
      enddo
      close(fu)
    endif
  end subroutine wv_dump_shm_to_file

end module m_wv_storage
