program main
  use m_cmdpath,only: setcmdpath
  use m_args,only: m_setargs
  integer:: ierr
  interface
     subroutine lmf(commin) bind(C)
       integer,optional:: commin
     end subroutine lmf
  end interface
  call mpi_init(ierr)
  call setcmdpath()
  call m_setargs()
  call lmf()
  call mpi_finalize(ierr)
  call exit(0)
end program main
