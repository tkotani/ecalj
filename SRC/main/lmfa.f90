program main
  integer:: ierr
  interface
     subroutine lmfa() bind(C)  
     end subroutine lmfa
  end interface
  call mpi_init(ierr)
  call lmfa()
  call mpi_finalize(ierr)
  call exit(0)
end program main
