program main
  use m_cmdpath,only: setcmdpath
  use m_args,only: m_setargs
  use m_lmfham1,only: lmfham1
  integer:: ierr
  call mpi_init(ierr)
  call setcmdpath()
  call m_setargs()
  call lmfham1()
  call mpi_finalize(ierr)
  call exit(0)
end program main
