program main
   integer:: ierr
   interface
      subroutine lmf() bind(C)
      end subroutine lmf
   end interface
   call mpi_init(ierr)
   call lmf()
end program main