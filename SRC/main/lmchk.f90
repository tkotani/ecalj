program main
  use m_cmdpath,only: setcmdpath
  use m_args,only: m_setargs
  use m_lmchk,only:lmchk
  implicit none
  integer:: ierr
  include "mpif.h" 
  call mpi_init(ierr)
  call setcmdpath()
  call m_setargs()
  call lmchk()
  call mpi_finalize(ierr)
  call exit(0)
end program main
