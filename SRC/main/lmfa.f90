program main
  use m_cmdpath,only: setcmdpath
  use m_args,only: m_setargs
  implicit none
  integer:: ierr
  include "mpif.h" 
  interface
     subroutine lmfa(commin) bind(C)  
       integer,optional:: commin
     end subroutine lmfa
  end interface
  call mpi_init(ierr)
  call setcmdpath()
  call m_setargs()
  call lmfa()
  call mpi_finalize(ierr)
  call exit(0)
end program main
