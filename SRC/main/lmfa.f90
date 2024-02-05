program main
  integer:: ierr
  interface
      subroutine lmfa() bind(C)  
      end subroutine lmfa
  end interface
  call mpi_init(ierr)
  call lmfa()
end program main
