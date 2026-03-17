module m_mpiio !MPI-IO only for complex(8). Fixed length recl
  use iso_c_binding
  use m_nvfortran
  use mpi
  implicit none
  type :: record_item
    type(c_ptr) :: addr
    integer     :: count
    integer     :: mpi_type
  endtype record_item
  public:: openm,writem,readm,closem, openedm
  public:: writem_c, writem_d, readm_d
  public:: record_item, record_item_from, writem_struct, readm_struct
  interface record_item_from
    module procedure record_item_from_real8_0d
    module procedure record_item_from_real8_1d
    module procedure record_item_from_real8_2d
    module procedure record_item_from_real8_3d
    module procedure record_item_from_int4_0d
    module procedure record_item_from_int4_1d
    module procedure record_item_from_int4_2d
    module procedure record_item_from_int4_3d
    module procedure record_item_from_complex8_0d
    module procedure record_item_from_complex8_1d
    module procedure record_item_from_complex8_2d
    module procedure record_item_from_complex8_3d
  endinterface
  private
  integer,parameter::nfmax=1000, nsize=16 !maxsize of opened file by openm
  integer :: ierr,fhl(nfmax)=-9999,iff=0  ! -9999 is used as a missing value indicator (assumed not to occur as a valid value)
  integer(kind=mpi_offset_kind)::recll(nfmax)
contains
  function openm(newunit,file,recl,comm) result(i) !recl=16*size
    integer::    newunit,     recl,info,amode,comm_in
    character(*)::     file
    integer, intent(in), optional :: comm
    integer:: i
    info = mpi_info_null
    comm_in = MPI_COMM_WORLD
    if(present(comm)) comm_in=comm
    call mpi_file_open(comm_in, trim(file), mpi_mode_rdwr + mpi_mode_create,MPI_INFO_NULL, newunit,ierr)
    iff=iff+1
    if(iff>nfmax) call rx('m_mpiio:iff>nfmax')
    fhl(iff)   = newunit
    recll(iff) = recl !in byte
    i=0
  end function openm
  function writem(unit,rec,data) result(i)
    integer::unit
    integer(mpi_offset_kind) :: offset
    integer::rec,count
    complex(8):: data(1)
    integer:: i,ifx
    integer:: status(MPI_Status_size)
    ifx = findloc(unit==fhl(1:iff),dim=1,value=.True.)
    offset= (rec-1)*recll(ifx)
    count = recll(ifx)/nsize     !    write(6,*)'writemmmmm',ifx,offset,count 
    call mpi_file_write_at(fhl(ifx), offset, data, count, MPI_DOUBLE_COMPLEX, status, ierr)
    i=0
  end function writem
  function writem_c(unit,rec,data) result(i)
    integer::unit
    integer(mpi_offset_kind) :: offset
    integer::rec,count
    complex(4):: data(1)
    integer:: i,ifx
    integer:: status(MPI_Status_size)
    ifx = findloc(unit==fhl(1:iff),dim=1,value=.True.)
    offset= (rec-1)*recll(ifx)
    count = recll(ifx)/8 !    write(6,*)'writemmmmm',ifx,offset,count 
    call mpi_file_write_at(fhl(ifx), offset, data, count, MPI_COMPLEX, status, ierr)
    i=0
  end function writem_c
  function writem_d(unit,rec,data) result(i)
    integer::unit
    integer(mpi_offset_kind) :: offset
    integer::rec,count
    real(8):: data(1)
    integer:: i,ifx
    integer:: status(MPI_Status_size)
    ifx = findloc(unit==fhl(1:iff),dim=1,value=.True.)
    offset= (rec-1)*recll(ifx)
    count = recll(ifx)/8 !    write(6,*)'writemmmmm',ifx,offset,count 
    call mpi_file_write_at(fhl(ifx), offset, data, count, MPI_DOUBLE_PRECISION, status, ierr)
    i=0
  end function writem_d
  function readm(unit,rec,data) result(i)
    integer::unit
    integer::rec,count
    integer(mpi_offset_kind) :: offset
    integer::  status(MPI_STATUS_SIZE)
    complex(8):: data(1)
    integer:: i,ifx
    ifx = findloc(unit==fhl,dim=1,value=.True.)
    offset=(rec-1)*recll(ifx)
    count=recll(ifx)/nsize
    call mpi_file_read_at(fhl(ifx), offset, data, count, MPI_DOUBLE_COMPLEX, status, ierr)
    i=0
  end function readm
  function readm_d(unit,rec,data) result(i)
    integer::unit
    integer::rec,count
    integer(mpi_offset_kind) :: offset
    integer::  status(MPI_STATUS_SIZE)
    real(8):: data(1)
    integer:: i,ifx
    ifx = findloc(unit==fhl,dim=1,value=.True.)
    offset=(rec-1)*recll(ifx)
    count=recll(ifx)/8
    call mpi_file_read_at(fhl(ifx), offset, data, count, MPI_DOUBLE_PRECISION, status, ierr)
    i=0
  end function readm_d
  function closem(unit) result(i)
    integer::unit
    integer:: i, ifx
    ifx = findloc(unit==fhl(1:iff),dim=1,value=.True.)
    fhl(ifx)=-9999
    call mpi_file_close(unit, ierr)
    i=0
  end function closem
  function openedm(unit) result(is_open)
    integer::unit
    logical:: is_open
    integer:: ifx
    ifx = findloc(unit==fhl(1:iff),dim=1,value=.True.)
    if(ifx>0) then
      is_open = .true.
    else
      is_open = .false.
    endif
  end function openedm

  subroutine build_struct_type(items, filetype, record_bytes)
    type(record_item), intent(in) :: items(:)
    integer, intent(out) :: filetype
    integer(kind=MPI_ADDRESS_KIND), intent(out) :: record_bytes
    integer :: n, i, ierr_local
    integer, allocatable :: blocklen(:), types(:)
    integer(kind=MPI_ADDRESS_KIND), allocatable  :: disp(:)
    n = size(items)
    allocate(blocklen(n), types(n), disp(n))
    do i = 1, n
      blocklen(i) = items(i)%count
      types(i)    = items(i)%mpi_type
      disp(i) = transfer(items(i)%addr, 0_MPI_ADDRESS_KIND) ! extract address stored in c_ptr
    enddo
    call MPI_Type_create_struct(n, blocklen, disp, types, filetype, ierr_local)
    call MPI_Type_commit(filetype, ierr_local)
    record_bytes = 0
    do i = 1, n
      select case(types(i))
      case(MPI_DOUBLE_PRECISION)
        record_bytes = record_bytes + 8_MPI_ADDRESS_KIND * blocklen(i)
      case(MPI_INTEGER)
        record_bytes = record_bytes + 4_MPI_ADDRESS_KIND * blocklen(i)
      case(MPI_DOUBLE_COMPLEX)
        record_bytes = record_bytes + 16_MPI_ADDRESS_KIND * blocklen(i)
      case default
        call rx('m_mpiio:build_struct_type: unsupported MPI type')
      endselect
    enddo
    deallocate(blocklen, types, disp)
  end subroutine build_struct_type

  integer function writem_struct(unit, rec, items) result(i)
    integer, intent(in) :: unit, rec
    type(record_item), intent(in) :: items(:)
    integer :: ifx, filetype
    integer(kind=MPI_ADDRESS_KIND) :: record_bytes
    integer(kind=MPI_OFFSET_KIND) :: offset
    integer :: status(MPI_STATUS_SIZE)
    ifx = findloc(unit == fhl(1:iff), dim=1, value=.true.)
    if (ifx <= 0) call rx('m_mpiio:writem_struct: unit not opened')
    call build_struct_type(items, filetype, record_bytes)
    if (record_bytes /= recll(ifx)) then
      write(6,*) 'm_mpiio:writem_struct: record_bytes mismatch:', record_bytes, recll(ifx)
      call rx('m_mpiio:writem_struct: record_bytes /= recll')
    endif
    offset = (rec-1)*recll(ifx)
    call MPI_File_write_at(fhl(ifx), offset, MPI_BOTTOM, 1, filetype, status, ierr)
    call MPI_Type_free(filetype, ierr)
    i = 0
  end function writem_struct

  integer function readm_struct(unit, rec, items) result(i)
    integer, intent(in) :: unit, rec
    type(record_item), intent(in) :: items(:)
    integer :: ifx, filetype
    integer(kind=MPI_ADDRESS_KIND) :: record_bytes
    integer(kind=MPI_OFFSET_KIND) :: offset
    integer :: status(MPI_STATUS_SIZE)
    ifx = findloc(unit == fhl(1:iff), dim=1, value=.true.)
    if (ifx <= 0) call rx('m_mpiio:readm_struct: unit not opened')
    call build_struct_type(items, filetype, record_bytes)
    if (record_bytes /= recll(ifx)) then
      write(6,*) 'm_mpiio:readm_struct: record_bytes mismatch:', record_bytes, recll(ifx)
      call rx('m_mpiio:readm_struct: record_bytes /= recll')
    endif
    offset = (rec-1)*recll(ifx)
    call MPI_File_read_at(fhl(ifx), offset, MPI_BOTTOM, 1, filetype, status, ierr)
    call MPI_Type_free(filetype, ierr)
    i = 0
  end function readm_struct

  function record_item_from_real8_0d(var) result(item)
    type(record_item) :: item
    real(8), intent(in), target :: var
    item%addr     = c_loc(var)
    item%count    = 1
    item%mpi_type = MPI_DOUBLE_PRECISION
  end function
  function record_item_from_real8_1d(var) result(item)
    type(record_item) :: item
    real(8), intent(in), target :: var(:)
    item%addr     = c_loc(var(1))
    item%count    = size(var)
    item%mpi_type = MPI_DOUBLE_PRECISION
  end function
  function record_item_from_real8_2d(var) result(item)
    type(record_item) :: item
    real(8), intent(in), target :: var(:, :)
    item%addr     = c_loc(var(1,1))
    item%count    = size(var)
    item%mpi_type = MPI_DOUBLE_PRECISION
  end function
  function record_item_from_real8_3d(var) result(item)
    type(record_item) :: item
    real(8), intent(in), target :: var(:, :, :)
    item%addr     = c_loc(var(1,1,1))
    item%count    = size(var)
    item%mpi_type = MPI_DOUBLE_PRECISION
  end function
  function record_item_from_int4_0d(var) result(item)
    type(record_item) :: item
    integer(4), intent(in), target :: var
    item%addr     = c_loc(var)
    item%count    = 1
    item%mpi_type = MPI_INTEGER
  end function
  function record_item_from_int4_1d(var) result(item)
    type(record_item) :: item
    integer(4), intent(in), target :: var(:)
    item%addr     = c_loc(var(1))
    item%count    = size(var)
    item%mpi_type = MPI_INTEGER
  end function
  function record_item_from_int4_2d(var) result(item)
    type(record_item) :: item
    integer(4), intent(in), target :: var(:, :)
    item%addr     = c_loc(var(1,1))
    item%count    = size(var)
    item%mpi_type = MPI_INTEGER
  end function
  function record_item_from_int4_3d(var) result(item)
    type(record_item) :: item
    integer(4), intent(in), target :: var(:, :, :)
    item%addr     = c_loc(var(1,1,1))
    item%count    = size(var)
    item%mpi_type = MPI_INTEGER
  end function
  function record_item_from_complex8_0d(var) result(item)
    type(record_item) :: item
    complex(8), intent(in), target :: var
    item%addr     = c_loc(var)
    item%count    = 1
    item%mpi_type = MPI_DOUBLE_COMPLEX
  end function
  function record_item_from_complex8_1d(var) result(item)
    type(record_item) :: item
    complex(8), intent(in), target :: var(:)
    item%addr     = c_loc(var(1))
    item%count    = size(var)
    item%mpi_type = MPI_DOUBLE_COMPLEX
  end function
  function record_item_from_complex8_2d(var) result(item)
    type(record_item) :: item
    complex(8), intent(in), target :: var(:, :)
    item%addr     = c_loc(var(1,1))
    item%count    = size(var)
    item%mpi_type = MPI_DOUBLE_COMPLEX
  end function
  function record_item_from_complex8_3d(var) result(item)
    type(record_item) :: item
    complex(8), intent(in), target :: var(:, :, :)
    item%addr     = c_loc(var(1,1,1))
    item%count    = size(var)
    item%mpi_type = MPI_DOUBLE_COMPLEX
  end function
end module m_mpiio
