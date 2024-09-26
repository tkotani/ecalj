module m_mpiio !MPI-IO only for complex(8). Fixed length recl
  use m_nvfortran
  use mpi
  implicit none
  public:: openm,writem,readm,closem
  private
  integer,parameter::nfmax=100 ,nsize=16
  integer::ierr,recll(nfmax),fhl(nfmax)=-9999,iff=0
contains
  function openm(newunit,file,recl) result(i) !recl=16*size
    integer::    newunit,     recl,info,amode
    character(*)::     file
    integer:: i
    i=0
    info = mpi_info_null
    call mpi_file_open(MPI_COMM_WORLD, trim(file), mpi_mode_rdwr + mpi_mode_create,MPI_INFO_NULL, newunit,ierr)
    iff=iff+1
    if(iff>nfmax) call rx('m_mpiio:iff>nfmax')
    fhl(iff)   = newunit
    recll(iff) = recl !in byte
  end function openm
  function writem(fh,rec,data) result(i)
    integer::fh
    integer(mpi_offset_kind) :: offset
    integer::rec,count
    complex(8):: data(1)
    integer:: i,ifx
    integer:: status(MPI_Status_size)
    i=0
    ifx = fhl(findloc(fh==fhl,dim=1,value=.True.))
    offset= (rec-1)*recll(ifx)
    count = recll(ifx)/nsize
!    write(6,*)'fffffffffffff222',ifx,offset,count !,recll(ifx),rec
!    call mpi_file_set_view(ifx,offset,MPI_DOUBLE_COMPLEX,MPI_DOUBLE_COMPLEX,'native',MPI_INFO_NULL,ierr)
!    call mpi_file_write(ifx,data,count,MPI_DOUBLE_COMPLEX,mpi_status_ignore,ierr)
    call mpi_file_write_at(ifx, offset, data, count, MPI_DOUBLE_COMPLEX, status, ierr)
!    write(6,*)'fffffffffffff222xxx',ifx,offset,count,ierr
!    call rx('xxxxxxxxx')
  end function writem
  function readm(fh,rec,data)result(i)
    integer::fh
    integer::rec,count
    integer(mpi_offset_kind) :: offset
    integer::  status(MPI_STATUS_SIZE)
    complex(8):: data(1)
    integer:: i,ifx
    i=0
    ifx = fhl(findloc(fh==fhl,dim=1,value=.True.))
    offset=(rec-1)*recll(ifx)
    count=recll(ifx)/nsize
    call mpi_file_read_at(ifx, offset, data, count, MPI_DOUBLE_COMPLEX, status, ierr)
  end function readm
  function closem(fh) result(i)
    integer::fh
    integer:: i
    i=0
    call mpi_file_close(fh, ierr)
  end function closem
end module m_mpiio
