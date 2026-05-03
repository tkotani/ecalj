program main
  use mpi
  use m_cmdpath,only: setcmdpath
  use m_args,only:    m_setargs
  use m_ext,only:     m_ext_init
  use m_lmchk,only:lmchk
  implicit none
  integer:: ierr,comm,procid
  call mpi_init(ierr)
  comm = MPI_COMM_WORLD
  call mpi_comm_rank(comm, procid, ierr )
  call setcmdpath()
  call m_setargs()
  call m_ext_init()    ! Get sname, e.g. trim(sname)=si of ctrl.si
  call lmchk()
  call mpi_finalize(ierr)
  call exit(0)
end program main
