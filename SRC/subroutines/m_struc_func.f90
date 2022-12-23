Module m_struc_func
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
    integer :: MAX_PROCS
    parameter (MAX_PROCS = 100)
    integer :: resultlen
    character*(MPI_MAX_PROCESSOR_NAME) name
    character(10) :: shortname(0:MAX_PROCS-1)
    character(26) :: datim
    integer :: namelen(0:MAX_PROCS-1)
    character(256) :: strn
    integer :: procid,master,i_data_size
    integer:: n=0
    master=0
    call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
    call mpi_bcast(struc%z, 1,MPI_REAL8   , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%rmt, 1,MPI_REAL8 , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%rg, 1,MPI_REAL8  , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%lmxa, 1,MPI_INTEGER, master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%lmxl, 1,MPI_INTEGER, master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%kmxt, 1,MPI_INTEGER, master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%rsmv, 1,MPI_REAL8  , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%a, 1,MPI_REAL8     , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%nr, 1,MPI_INTEGER  , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%lfoca, 1,MPI_INTEGER,master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%ctail, 1,MPI_REAL8 , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%etail, 1,MPI_REAL8 , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%stc, 1,MPI_REAL8   , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%lmxb, 1,MPI_INTEGER, master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%rfoca, 1,MPI_REAL8 , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%nxi, 1,MPI_INTEGER , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%qc, 1,MPI_REAL8    , master, MPI_COMM_WORLD,ierr)
    call mpi_bcast(struc%exi, size(struc%exi),MPI_REAL8 , master, MPI_COMM_WORLD,ierr) 
    call mpi_bcast(struc%chfa,size(struc%chfa),MPI_REAL8, master, MPI_COMM_WORLD,ierr)
    if (mlog) then
       call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
       shortname(procid)= name !call strcop(shortname(procid),name,10,'.',ierr)
       namelen(procid) = ierr-1
       call gettime(datim)
       write(stml,ftox)'  '//datim//' Process ',procid,' of ',numprocs,' on ' &
            //shortname(procid)(1:namelen(procid))//' '//label
    endif
  end subroutine mpibc1_s_spec
end module m_struc_func
