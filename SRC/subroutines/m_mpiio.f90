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
  integer, parameter :: nfmax=1000, nsize=16 !maxsize of opened file by openm
  integer :: ierr, fhl(nfmax)=-9999, iff=0   ! -9999 is used as a missing value indicator
  integer(kind=mpi_offset_kind) :: recll(nfmax)
contains
  function openm(newunit, file, recl, comm) result(i) !recl=16*size
    integer ::      newunit,      recl, info, comm_in
    character(*) ::       file
    integer, intent(in), optional :: comm
    integer :: i
    info = mpi_info_null
    comm_in = MPI_COMM_WORLD
    if (present(comm)) comm_in = comm
    call mpi_file_open(comm_in, trim(file), mpi_mode_rdwr + mpi_mode_create, MPI_INFO_NULL, newunit, ierr)
    iff = iff + 1
    if (iff > nfmax) call rx('m_mpiio:iff>nfmax')
    fhl(iff)   = newunit
    recll(iff) = recl !in byte
    i = 0
  end function openm
  function writem(unit, rec, data) result(i)
    integer :: unit
    integer(mpi_offset_kind) :: offset
    integer :: rec, count
    complex(8) :: data(1)
    integer :: i, ifx
    integer :: status(MPI_Status_size)
    ifx = findloc(unit==fhl(1:iff), dim=1, value=.True.)
    offset = (rec-1)*recll(ifx)
    count  = recll(ifx)/nsize
    call mpi_file_write_at(fhl(ifx), offset, data, count, MPI_DOUBLE_COMPLEX, status, ierr)
    i = 0
  end function writem
  function writem_c(unit, rec, data) result(i)
    integer :: unit
    integer(mpi_offset_kind) :: offset
    integer :: rec, count
    complex(4) :: data(1)
    integer :: i, ifx
    integer :: status(MPI_Status_size)
    ifx = findloc(unit==fhl(1:iff), dim=1, value=.True.)
    offset = (rec-1)*recll(ifx)
    count  = recll(ifx)/8
    call mpi_file_write_at(fhl(ifx), offset, data, count, MPI_COMPLEX, status, ierr)
    i = 0
  end function writem_c
  function writem_d(unit, rec, data) result(i)
    integer :: unit
    integer(mpi_offset_kind) :: offset
    integer :: rec, count
    real(8) :: data(1)
    integer :: i, ifx
    integer :: status(MPI_Status_size)
    ifx = findloc(unit==fhl(1:iff), dim=1, value=.True.)
    offset = (rec-1)*recll(ifx)
    count  = recll(ifx)/8
    call mpi_file_write_at(fhl(ifx), offset, data, count, MPI_DOUBLE_PRECISION, status, ierr)
    i = 0
  end function writem_d
  function readm(unit, rec, data) result(i)
    integer :: unit
    integer :: rec, count
    integer(mpi_offset_kind) :: offset
    integer :: status(MPI_STATUS_SIZE)
    complex(8) :: data(1)
    integer :: i, ifx
    ifx = findloc(unit==fhl, dim=1, value=.True.)
    offset = (rec-1)*recll(ifx)
    count  = recll(ifx)/nsize
    call mpi_file_read_at(fhl(ifx), offset, data, count, MPI_DOUBLE_COMPLEX, status, ierr)
    i = 0
  end function readm
  function readm_d(unit, rec, data) result(i)
    integer :: unit
    integer :: rec, count
    integer(mpi_offset_kind) :: offset
    integer :: status(MPI_STATUS_SIZE)
    real(8) :: data(1)
    integer :: i, ifx
    ifx = findloc(unit==fhl, dim=1, value=.True.)
    offset = (rec-1)*recll(ifx)
    count  = recll(ifx)/8
    call mpi_file_read_at(fhl(ifx), offset, data, count, MPI_DOUBLE_PRECISION, status, ierr)
    i = 0
  end function readm_d
  function closem(unit) result(i)
    integer :: unit
    integer :: i, ifx
    ifx = findloc(unit==fhl(1:iff), dim=1, value=.True.)
    fhl(ifx) = -9999
    call mpi_file_close(unit, ierr)
    i = 0
  end function closem
  function openedm(unit) result(is_open)
    integer :: unit
    logical :: is_open
    integer :: ifx
    ifx = findloc(unit==fhl(1:iff), dim=1, value=.True.)
    is_open = ifx > 0
  end function openedm

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
    integer, intent(in) :: unit, rec
    type(mpiio_buf), intent(in) :: buf
    integer :: ifx, nbytes
    integer(mpi_offset_kind) :: offset
    integer :: status(MPI_STATUS_SIZE)
    ifx = findloc(unit == fhl(1:iff), dim=1, value=.true.)
    if (ifx <= 0) call rx('m_mpiio:writem_buf: unit not opened')
    nbytes = buf%pos - 1
    if (nbytes /= recll(ifx)) then
      write(6,*) 'm_mpiio:writem_buf: buffer size mismatch:', nbytes, recll(ifx)
      call rx('m_mpiio:writem_buf: buffer size /= recll')
    end if
    offset = (rec-1)*recll(ifx)
    call MPI_File_write_at(fhl(ifx), offset, buf%bytes(1), nbytes, MPI_BYTE, status, ierr)
    i = ierr
  end function writem_buf

  integer function readm_buf(unit, rec, buf) result(i)
    integer, intent(in) :: unit, rec
    type(mpiio_buf), intent(inout) :: buf
    integer :: ifx, nbytes
    integer(mpi_offset_kind) :: offset
    integer :: status(MPI_STATUS_SIZE)
    ifx = findloc(unit == fhl(1:iff), dim=1, value=.true.)
    if (ifx <= 0) call rx('m_mpiio:readm_buf: unit not opened')
    nbytes = int(recll(ifx))
    if (.not. allocated(buf%bytes) .or. size(buf%bytes) /= nbytes) then
      if (allocated(buf%bytes)) deallocate(buf%bytes)
      allocate(buf%bytes(nbytes))
    end if
    buf%pos = 1
    offset = (rec-1)*recll(ifx)
    call MPI_File_read_at(fhl(ifx), offset, buf%bytes(1), nbytes, MPI_BYTE, status, ierr)
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
