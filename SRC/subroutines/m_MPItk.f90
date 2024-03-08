!> MPI utility routines and variablis by TK
module m_MPItk
  use m_lgunit, only: stml, stdl
  public :: m_MPItk_init, procid, strprocid,master, nsize, master_mpi, xmpbnd2
  private
  integer:: procid, master = 0, nsize !,nproc
  include "mpif.h"
  logical:: master_mpi
  character(8) :: strprocid
contains
  subroutine m_MPItk_init()
    logical:: cmdopt0
    integer :: fext,ierr
    character(10):: i2char
    call mpi_comm_size(MPI_COMM_WORLD, nsize, ierr)
    call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
    strprocid = trim(i2char(procid))
    master_mpi = .false.
    if (procid == master) master_mpi = .TRUE.
  end subroutine m_MPItk_init
  subroutine xmpbnd2(kpproc, ndham, ndat, eb)       !- Collect eb from various processors (MPI)
    !i Inputs
    !i   kpproc
    !i   ndham :leading dimension of eb
    !i   nkp   :number of irreducible k-points (bzmesh.f)
    !i   nsp   :2 for spin-polarized case, otherwise 1
    !i   eb    :energy bands; alias eband
    !o Outputs
    implicit none
    integer:: kpproc(0:*), ndham, ndat
    double precision :: eb(ndham, ndat)
    integer :: i, ista, iend, ierr
    integer, dimension(:), allocatable :: offset, length
    real(8), allocatable :: buf_rv(:, :)
    allocate (offset(0:nsize), stat=ierr)
    allocate (length(0:nsize), stat=ierr)
    offset(0) = 0
    do i = 0, nsize - 1 !NOTE: ndat is divided into kpproc
      ista = kpproc(i)
      iend = kpproc(i + 1) - 1
      length(i) = (iend - ista + 1)*ndham !nsp*ndham
      offset(i + 1) = offset(i) + length(i)
    end do
    ista = kpproc(procid)
    allocate (buf_rv(ndham, ndat))
    call mpi_allgatherv(eb(1, ista), length(procid), mpi_double_precision, buf_rv, length, offset, mpi_double_precision,&
      mpi_comm_world, ierr)
    eb = buf_rv
    deallocate (buf_rv)
    deallocate (offset, stat=ierr)
    deallocate (length, stat=ierr)
  end subroutine xmpbnd2
end module m_MPItk
