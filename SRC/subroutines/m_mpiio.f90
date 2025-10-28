module m_mpiio !MPI-IO only for complex(8). Fixed length recl
  use m_nvfortran
  use mpi
  implicit none
  public:: openm,writem,readm,closem
  public:: writem_kind4
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
  function writem_kind4(unit,rec,data) result(i)
    integer::unit
    integer(mpi_offset_kind) :: offset
    integer::rec,count
    complex(4):: data(1)
    integer:: i,ifx
    integer:: status(MPI_Status_size)
    ifx = findloc(unit==fhl(1:iff),dim=1,value=.True.)
    offset= (rec-1)*recll(ifx)
    count = recll(ifx)/nsize     !    write(6,*)'writemmmmm',ifx,offset,count 
    call mpi_file_write_at(fhl(ifx), offset, data, count, MPI_COMPLEX, status, ierr)
    i=0
  end function writem_kind4
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
  function closem(unit) result(i)
    integer::unit
    integer:: i, ifx
    ifx = findloc(unit==fhl(1:iff),dim=1,value=.True.)
    fhl(ifx)=-9999
    call mpi_file_close(unit, ierr)
    i=0
  end function closem
end module m_mpiio
