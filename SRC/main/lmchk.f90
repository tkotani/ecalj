program main
  integer:: ierr
  interface
     subroutine lmchk() bind(C)
     end subroutine lmchk
  end interface
  call mpi_init(ierr)
  call lmchk()
endprogram main
