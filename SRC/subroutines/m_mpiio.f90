module m_mpiio !MPI-IO. Fixed length recl
  use m_nvfortran
  use mpi
  implicit none
  ! mpiio_buf: byte buffer for packing (write) and unpacking (read) heterogeneous data.
  ! buf_put appends data to the buffer; buf_get extracts data sequentially.
  ! writem_buf writes the buffer to file; readm_buf reads a record into the buffer.
  ! This avoids c_loc/c_ptr, so TARGET attributes and contiguous arrays are not required.
  type :: mpiio_buf
    integer(1), allocatable :: bytes(:)
    integer :: pos = 1
  end type mpiio_buf
  public :: openm, writem, readm, closem, openedm
  public :: writem_c, writem_d, readm_d
  public :: mpiio_buf, buf_put, buf_get, buf_reset, writem_buf, readm_buf
  interface buf_put
    module procedure buf_put_real8_0d,    buf_put_real8_1d,    buf_put_real8_2d,    buf_put_real8_3d
    module procedure buf_put_int4_0d,     buf_put_int4_1d,     buf_put_int4_2d,     buf_put_int4_3d
    module procedure buf_put_complex8_0d, buf_put_complex8_1d, buf_put_complex8_2d, buf_put_complex8_3d
  end interface buf_put
  interface buf_get
    module procedure buf_get_real8_0d,    buf_get_real8_1d,    buf_get_real8_2d,    buf_get_real8_3d
    module procedure buf_get_int4_0d,     buf_get_int4_1d,     buf_get_int4_2d,     buf_get_int4_3d
    module procedure buf_get_complex8_0d, buf_get_complex8_1d, buf_get_complex8_2d, buf_get_complex8_3d
  end interface buf_get
  private

  ! Entry in the open-file table: MPI_FILE_NULL fh means the slot is free.
  type :: mpiio_entry
    integer                        :: fh    = MPI_FILE_NULL
    integer(kind=mpi_offset_kind)  :: recl  = 0
    character(256)                 :: fname = ''
  end type mpiio_entry

  integer, parameter        :: nfmax = 32   ! max simultaneously open MPI-IO files
  integer, parameter        :: nsize = 16   ! sizeof(complex(8)) in bytes
  type(mpiio_entry)         :: fh_table(nfmax)
  integer                   :: ierr

contains
  function openm(newunit, file, recl, comm) result(i)
    integer,      intent(out) :: newunit
    character(*), intent(in)  :: file
    integer,      intent(in)  :: recl
    integer,      intent(in), optional :: comm
    integer :: i, ifx, comm_in
    comm_in = MPI_COMM_WORLD
    if (present(comm)) comm_in = comm
    call mpi_file_open(comm_in, trim(file), mpi_mode_rdwr + mpi_mode_create, MPI_INFO_NULL, newunit, ierr)
    ifx = findloc(fh_table%fh == MPI_FILE_NULL, dim=1, value=.True.)
    if (ifx <= 0) call rx('m_mpiio:openm: too many simultaneously open files')
    fh_table(ifx)%fh    = newunit
    fh_table(ifx)%recl  = recl
    fh_table(ifx)%fname = trim(file)
    i = 0
  end function openm

  function writem(unit, rec, data) result(i)
    integer,     intent(in) :: unit, rec
    complex(8),  intent(in) :: data(1)
    integer :: i, ifx
    integer :: status(MPI_Status_size)
    integer(mpi_offset_kind) :: offset
    ifx = find_slot(unit)
    offset = (rec-1) * fh_table(ifx)%recl
    call mpi_file_write_at(fh_table(ifx)%fh, offset, data, int(fh_table(ifx)%recl/nsize), MPI_DOUBLE_COMPLEX, status, ierr)
    i = 0
  end function writem

  function writem_c(unit, rec, data) result(i)
    integer,    intent(in) :: unit, rec
    complex(4), intent(in) :: data(1)
    integer :: i, ifx
    integer :: status(MPI_Status_size)
    integer(mpi_offset_kind) :: offset
    ifx = find_slot(unit)
    offset = (rec-1) * fh_table(ifx)%recl
    call mpi_file_write_at(fh_table(ifx)%fh, offset, data, int(fh_table(ifx)%recl/8), MPI_COMPLEX, status, ierr)
    i = 0
  end function writem_c

  function writem_d(unit, rec, data) result(i)
    integer, intent(in) :: unit, rec
    real(8), intent(in) :: data(1)
    integer :: i, ifx
    integer :: status(MPI_Status_size)
    integer(mpi_offset_kind) :: offset
    ifx = find_slot(unit)
    offset = (rec-1) * fh_table(ifx)%recl
    call mpi_file_write_at(fh_table(ifx)%fh, offset, data, int(fh_table(ifx)%recl/8), MPI_DOUBLE_PRECISION, status, ierr)
    i = 0
  end function writem_d

  function readm(unit, rec, data) result(i)
    integer,    intent(in)  :: unit, rec
    complex(8), intent(out) :: data(1)
    integer :: i, ifx
    integer :: status(MPI_STATUS_SIZE)
    integer(mpi_offset_kind) :: offset
    ifx = find_slot(unit)
    offset = (rec-1) * fh_table(ifx)%recl
    call mpi_file_read_at(fh_table(ifx)%fh, offset, data, int(fh_table(ifx)%recl/nsize), MPI_DOUBLE_COMPLEX, status, ierr)
    i = 0
  end function readm

  function readm_d(unit, rec, data) result(i)
    integer, intent(in)  :: unit, rec
    real(8), intent(out) :: data(1)
    integer :: i, ifx
    integer :: status(MPI_STATUS_SIZE)
    integer(mpi_offset_kind) :: offset
    ifx = find_slot(unit)
    offset = (rec-1) * fh_table(ifx)%recl
    call mpi_file_read_at(fh_table(ifx)%fh, offset, data, int(fh_table(ifx)%recl/8), MPI_DOUBLE_PRECISION, status, ierr)
    i = 0
  end function readm_d

  function closem(unit) result(i)
    integer, intent(in) :: unit
    integer :: i, ifx
    ifx = find_slot(unit)
    call mpi_file_close(fh_table(ifx)%fh, ierr)
    fh_table(ifx)%fh    = MPI_FILE_NULL
    fh_table(ifx)%fname = ''
    i = 0
  end function closem

  function openedm(unit) result(is_open)
    integer, intent(in) :: unit
    logical :: is_open
    is_open = findloc(fh_table%fh == unit, dim=1, value=.True.) > 0
  end function openedm

  integer function find_slot(unit) result(ifx)
    integer, intent(in) :: unit
    ifx = findloc(fh_table%fh == unit, dim=1, value=.True.)
    if (ifx <= 0) call rx('m_mpiio: MPI file handle not in table')
  end function find_slot

  subroutine buf_reset(buf)
    type(mpiio_buf), intent(inout) :: buf
    buf%pos = 1
  end subroutine buf_reset

  subroutine buf_grow(buf, needed)
    type(mpiio_buf), intent(inout) :: buf
    integer, intent(in) :: needed
    integer :: new_size
    integer(1), allocatable :: tmp(:)
    if (.not. allocated(buf%bytes)) then
      allocate(buf%bytes(max(needed, 1024)))
      return
    end if
    if (buf%pos + needed - 1 > size(buf%bytes)) then
      new_size = max(size(buf%bytes)*2, buf%pos + needed - 1)
      allocate(tmp(new_size))
      tmp(1:buf%pos-1) = buf%bytes(1:buf%pos-1)
      call move_alloc(tmp, buf%bytes)
    end if
  end subroutine buf_grow

  integer function writem_buf(unit, rec, buf) result(i)
    integer,         intent(in) :: unit, rec
    type(mpiio_buf), intent(in) :: buf
    integer :: ifx, nbytes
    integer(mpi_offset_kind) :: offset
    integer :: status(MPI_STATUS_SIZE)
    ifx = find_slot(unit)
    nbytes = buf%pos - 1
    if (nbytes /= fh_table(ifx)%recl) then
      write(6,*) 'm_mpiio:writem_buf: buffer size mismatch (file='//trim(fh_table(ifx)%fname)//'):', nbytes, fh_table(ifx)%recl
      call rx('m_mpiio:writem_buf: buffer size /= recl (file='//trim(fh_table(ifx)%fname)//')')
    end if
    offset = (rec-1) * fh_table(ifx)%recl
    call MPI_File_write_at(fh_table(ifx)%fh, offset, buf%bytes(1), nbytes, MPI_BYTE, status, ierr)
    i = ierr
  end function writem_buf

  integer function readm_buf(unit, rec, buf) result(i)
    integer,         intent(in)    :: unit, rec
    type(mpiio_buf), intent(inout) :: buf
    integer :: ifx, nbytes
    integer(mpi_offset_kind) :: offset
    integer :: status(MPI_STATUS_SIZE)
    ifx = find_slot(unit)
    nbytes = int(fh_table(ifx)%recl)
    if (.not. allocated(buf%bytes) .or. size(buf%bytes) /= nbytes) then
      if (allocated(buf%bytes)) deallocate(buf%bytes)
      allocate(buf%bytes(nbytes))
    end if
    buf%pos = 1
    offset = (rec-1) * fh_table(ifx)%recl
    call MPI_File_read_at(fh_table(ifx)%fh, offset, buf%bytes(1), nbytes, MPI_BYTE, status, ierr)
    i = ierr
  end function readm_buf

  ! ===== buf_put =====
  subroutine buf_put_real8_0d(buf, var)
    type(mpiio_buf), intent(inout) :: buf; real(8), intent(in) :: var
    integer, parameter :: nb = 8
    call buf_grow(buf, nb)
    buf%bytes(buf%pos:buf%pos+nb-1) = transfer(var, buf%bytes(buf%pos:buf%pos+nb-1))
    buf%pos = buf%pos + nb
  end subroutine
  subroutine buf_put_real8_1d(buf, var)
    type(mpiio_buf), intent(inout) :: buf; real(8), intent(in) :: var(:)
    integer :: nb
    nb = size(var)*8
    call buf_grow(buf, nb)
    buf%bytes(buf%pos:buf%pos+nb-1) = transfer(var, buf%bytes(buf%pos:buf%pos+nb-1))
    buf%pos = buf%pos + nb
  end subroutine
  subroutine buf_put_real8_2d(buf, var)
    type(mpiio_buf), intent(inout) :: buf; real(8), intent(in) :: var(:,:)
    integer :: nb
    nb = size(var)*8
    call buf_grow(buf, nb)
    buf%bytes(buf%pos:buf%pos+nb-1) = transfer(var, buf%bytes(buf%pos:buf%pos+nb-1))
    buf%pos = buf%pos + nb
  end subroutine
  subroutine buf_put_real8_3d(buf, var)
    type(mpiio_buf), intent(inout) :: buf; real(8), intent(in) :: var(:,:,:)
    integer :: nb
    nb = size(var)*8
    call buf_grow(buf, nb)
    buf%bytes(buf%pos:buf%pos+nb-1) = transfer(var, buf%bytes(buf%pos:buf%pos+nb-1))
    buf%pos = buf%pos + nb
  end subroutine
  subroutine buf_put_int4_0d(buf, var)
    type(mpiio_buf), intent(inout) :: buf; integer(4), intent(in) :: var
    integer, parameter :: nb = 4
    call buf_grow(buf, nb)
    buf%bytes(buf%pos:buf%pos+nb-1) = transfer(var, buf%bytes(buf%pos:buf%pos+nb-1))
    buf%pos = buf%pos + nb
  end subroutine
  subroutine buf_put_int4_1d(buf, var)
    type(mpiio_buf), intent(inout) :: buf; integer(4), intent(in) :: var(:)
    integer :: nb
    nb = size(var)*4
    call buf_grow(buf, nb)
    buf%bytes(buf%pos:buf%pos+nb-1) = transfer(var, buf%bytes(buf%pos:buf%pos+nb-1))
    buf%pos = buf%pos + nb
  end subroutine
  subroutine buf_put_int4_2d(buf, var)
    type(mpiio_buf), intent(inout) :: buf; integer(4), intent(in) :: var(:,:)
    integer :: nb
    nb = size(var)*4
    call buf_grow(buf, nb)
    buf%bytes(buf%pos:buf%pos+nb-1) = transfer(var, buf%bytes(buf%pos:buf%pos+nb-1))
    buf%pos = buf%pos + nb
  end subroutine
  subroutine buf_put_int4_3d(buf, var)
    type(mpiio_buf), intent(inout) :: buf; integer(4), intent(in) :: var(:,:,:)
    integer :: nb
    nb = size(var)*4
    call buf_grow(buf, nb)
    buf%bytes(buf%pos:buf%pos+nb-1) = transfer(var, buf%bytes(buf%pos:buf%pos+nb-1))
    buf%pos = buf%pos + nb
  end subroutine
  subroutine buf_put_complex8_0d(buf, var)
    type(mpiio_buf), intent(inout) :: buf; complex(8), intent(in) :: var
    integer, parameter :: nb = 16
    call buf_grow(buf, nb)
    buf%bytes(buf%pos:buf%pos+nb-1) = transfer(var, buf%bytes(buf%pos:buf%pos+nb-1))
    buf%pos = buf%pos + nb
  end subroutine
  subroutine buf_put_complex8_1d(buf, var)
    type(mpiio_buf), intent(inout) :: buf; complex(8), intent(in) :: var(:)
    integer :: nb
    nb = size(var)*16
    call buf_grow(buf, nb)
    buf%bytes(buf%pos:buf%pos+nb-1) = transfer(var, buf%bytes(buf%pos:buf%pos+nb-1))
    buf%pos = buf%pos + nb
  end subroutine
  subroutine buf_put_complex8_2d(buf, var)
    type(mpiio_buf), intent(inout) :: buf; complex(8), intent(in) :: var(:,:)
    integer :: nb
    nb = size(var)*16
    call buf_grow(buf, nb)
    buf%bytes(buf%pos:buf%pos+nb-1) = transfer(var, buf%bytes(buf%pos:buf%pos+nb-1))
    buf%pos = buf%pos + nb
  end subroutine
  subroutine buf_put_complex8_3d(buf, var)
    type(mpiio_buf), intent(inout) :: buf; complex(8), intent(in) :: var(:,:,:)
    integer :: nb
    nb = size(var)*16
    call buf_grow(buf, nb)
    buf%bytes(buf%pos:buf%pos+nb-1) = transfer(var, buf%bytes(buf%pos:buf%pos+nb-1))
    buf%pos = buf%pos + nb
  end subroutine

  ! ===== buf_get =====
  subroutine buf_get_real8_0d(buf, var)
    type(mpiio_buf), intent(inout) :: buf; real(8), intent(out) :: var
    integer, parameter :: nb = 8
    var = transfer(buf%bytes(buf%pos:buf%pos+nb-1), var)
    buf%pos = buf%pos + nb
  end subroutine
  subroutine buf_get_real8_1d(buf, var)
    type(mpiio_buf), intent(inout) :: buf; real(8), intent(out) :: var(:)
    integer :: nb
    nb = size(var)*8
    var = transfer(buf%bytes(buf%pos:buf%pos+nb-1), var)
    buf%pos = buf%pos + nb
  end subroutine
  subroutine buf_get_real8_2d(buf, var)
    type(mpiio_buf), intent(inout) :: buf; real(8), intent(out) :: var(:,:)
    integer :: nb
    nb = size(var)*8
    var = reshape(transfer(buf%bytes(buf%pos:buf%pos+nb-1), [0d0]), shape(var))
    buf%pos = buf%pos + nb
  end subroutine
  subroutine buf_get_real8_3d(buf, var)
    type(mpiio_buf), intent(inout) :: buf; real(8), intent(out) :: var(:,:,:)
    integer :: nb
    nb = size(var)*8
    var = reshape(transfer(buf%bytes(buf%pos:buf%pos+nb-1), [0d0]), shape(var))
    buf%pos = buf%pos + nb
  end subroutine
  subroutine buf_get_int4_0d(buf, var)
    type(mpiio_buf), intent(inout) :: buf; integer(4), intent(out) :: var
    integer, parameter :: nb = 4
    var = transfer(buf%bytes(buf%pos:buf%pos+nb-1), var)
    buf%pos = buf%pos + nb
  end subroutine
  subroutine buf_get_int4_1d(buf, var)
    type(mpiio_buf), intent(inout) :: buf; integer(4), intent(out) :: var(:)
    integer :: nb
    nb = size(var)*4
    var = transfer(buf%bytes(buf%pos:buf%pos+nb-1), var)
    buf%pos = buf%pos + nb
  end subroutine
  subroutine buf_get_int4_2d(buf, var)
    type(mpiio_buf), intent(inout) :: buf; integer(4), intent(out) :: var(:,:)
    integer :: nb
    nb = size(var)*4
    var = reshape(transfer(buf%bytes(buf%pos:buf%pos+nb-1), [0_4]), shape(var))
    buf%pos = buf%pos + nb
  end subroutine
  subroutine buf_get_int4_3d(buf, var)
    type(mpiio_buf), intent(inout) :: buf; integer(4), intent(out) :: var(:,:,:)
    integer :: nb
    nb = size(var)*4
    var = reshape(transfer(buf%bytes(buf%pos:buf%pos+nb-1), [0_4]), shape(var))
    buf%pos = buf%pos + nb
  end subroutine
  subroutine buf_get_complex8_0d(buf, var)
    type(mpiio_buf), intent(inout) :: buf; complex(8), intent(out) :: var
    integer, parameter :: nb = 16
    var = transfer(buf%bytes(buf%pos:buf%pos+nb-1), var)
    buf%pos = buf%pos + nb
  end subroutine
  subroutine buf_get_complex8_1d(buf, var)
    type(mpiio_buf), intent(inout) :: buf; complex(8), intent(out) :: var(:)
    integer :: nb
    nb = size(var)*16
    var = transfer(buf%bytes(buf%pos:buf%pos+nb-1), var)
    buf%pos = buf%pos + nb
  end subroutine
  subroutine buf_get_complex8_2d(buf, var)
    type(mpiio_buf), intent(inout) :: buf; complex(8), intent(out) :: var(:,:)
    integer :: nb
    nb = size(var)*16
    var = reshape(transfer(buf%bytes(buf%pos:buf%pos+nb-1), [(cmplx(0d0,0d0,8))]), shape(var))
    buf%pos = buf%pos + nb
  end subroutine
  subroutine buf_get_complex8_3d(buf, var)
    type(mpiio_buf), intent(inout) :: buf; complex(8), intent(out) :: var(:,:,:)
    integer :: nb
    nb = size(var)*16
    var = reshape(transfer(buf%bytes(buf%pos:buf%pos+nb-1), [(cmplx(0d0,0d0,8))]), shape(var))
    buf%pos = buf%pos + nb
  end subroutine
end module m_mpiio
