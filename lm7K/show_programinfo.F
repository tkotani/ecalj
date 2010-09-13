      subroutine show_programinfo(io)
      implicit none
      integer:: io
      character(17),parameter:: info='HOST_INFORMATION '
#if MPI|MPIK
      include 'mpif.h'
      integer:: procid, ierr
      call mpi_barrier( MPI_COMM_WORLD, ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
      if (procid.eq.0) then
#endif
       write(io,'(a,a,a)') info,'platform: ',___PLATFORM___
       write(io,'(a,a,a)') info,'compiler version: ',___FC_VERSION___
       write(io,'(a,a,a)') info,'FFLAGS (<=120): ',
     .___FFLAGS___
       write(io,'(a,a,a)') info,'LIBLOC (<=120): ',
     .___LIBLOC___
       write(io,'(a,a,a)') info,'uname -a (<=120): ',
     .___UNAME_A___
       write(io,'(a,a,a)') info,'/etc/issue: ',___ETC_ISSUE___
       write(io,'(a,a,a)') info,'git branch: ',___GIT_BRANCH___
       write(io,'(a,a,a)') info,'git commit: ',___GIT_COMMIT___
       write(io,'(a,a,a)') info,'linked at: ',___LINK_TIME___
#if MPI|MPIK
       endif
#endif
       end subroutine show_programinfo

C       program test
C       call show_programinfo(6)
C       end 
