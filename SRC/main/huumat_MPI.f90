program main !wannier
  use mpi
  use m_cmdpath,only: setcmdpath, cmdpath
  use m_args,only:    m_setargs, argall
  use m_ext,only:     m_ext_init, sname
  use m_mpi, only: MPI__Initialize, mpi__root, comm, setipr, mpi__rank
  use m_uumat, only: uumatrix
  integer :: ierr

  call MPI__Initialize()
  call setipr(comm)
  call setcmdpath()
  call m_setargs()
  call m_ext_init()    ! Get sname, e.g. trim(sname)=si of ctrl.si
  call mpi_barrier(comm, ierr) !wait finishing of ctrl2ctrlp
  call uumatrix()
end program
