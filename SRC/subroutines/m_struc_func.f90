!> mpi broadcase of spec
module m_struc_func
  use m_ftox,only: ftox
  use m_lgunit,only:stml
  public  mpibc1_s_spec
contains
  subroutine mpibc1_s_spec(struc,label)
    use m_struc_def, only: s_spec
    use m_MPItk,only: mlog
    implicit none
    type(s_spec):: struc
    character label*(*)
    include 'mpif.h'
    integer :: numprocs, ierr
    integer :: resultlen
    character(26) :: datim
    character(256) :: strn
    integer :: procid,master,i_data_size
    integer:: n=0
    master=0
    call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
    call mpi_bcast(struc%z,     1,MPI_REAL8   , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%rmt,   1,MPI_REAL8 , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%rg,    1,MPI_REAL8  , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%lmxa,  1,MPI_INTEGER, master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%lmxl,  1,MPI_INTEGER, master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%kmxt,  1,MPI_INTEGER, master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%rsmv,  1,MPI_REAL8  , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%a,     1,MPI_REAL8     , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%nr,    1,MPI_INTEGER  , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%lfoca, 1,MPI_INTEGER,master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%ctail, 1,MPI_REAL8 , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%etail, 1,MPI_REAL8 , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%stc,   1,MPI_REAL8   , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%lmxb,  1,MPI_INTEGER, master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%rfoca, 1,MPI_REAL8 , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%nxi,   1,MPI_INTEGER , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%qc,    1,MPI_REAL8    , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%exi, size(struc%exi),MPI_REAL8 , master, MPI_COMM_WORLD,ierr) 
    call mpi_bcast(struc%chfa,size(struc%chfa),MPI_REAL8, master, MPI_COMM_WORLD,ierr)
    if (mlog) then
       call gettime(datim)
       write(stml,ftox)'  '//datim//' Process',procid,'of',numprocs,' '//label
    endif
  end subroutine mpibc1_s_spec
end module m_struc_func
