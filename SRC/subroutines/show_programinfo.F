c$$$      subroutine show_programinfo(io)
c$$$      implicit none
c$$$      integer:: io
c$$$      character(17),parameter:: info='HOST_INFORMATION '
c$$$#if MPI|MPIK
c$$$      include 'mpif.h'
c$$$      integer:: procid, ierr
c$$$      call mpi_barrier( MPI_COMM_WORLD, ierr )
c$$$      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
c$$$      if (procid.eq.0) then
c$$$#endif
c$$$        write(io,'(a,a,a)') info,'platform: ',___PLATFORM___
c$$$        write(io,'(a,a,a)') info,'compiler version: ',___FC_VERSION___
c$$$        write(io,'(a,a,a)') info,'FFLAGS (<=120): ',
c$$$     .___FFLAGS___
c$$$        write(io,'(a,a,a)') info,'LIBLOC (<=120): ',
c$$$     .___LIBLOC___
c$$$        write(io,'(a,a,a)') info,'uname -a (<=120): ',
c$$$     .___UNAME_A___
c$$$        write(io,'(a,a,a)') info,'/etc/issue: ',___ETC_ISSUE___
c$$$        write(io,'(a,a,a)') info,'git branch: ',___GIT_BRANCH___
c$$$        write(io,'(a,a,a)') info,'git commit: ',___GIT_COMMIT___
c$$$        write(io,'(a,a,a)') info,'linked at: ',___LINK_TIME___
c$$$#if MPI|MPIK
c$$$      endif
c$$$#endif
c$$$      end subroutine show_programinfo
c$$$
c$$$C       program test
c$$$C       call show_programinfo(6)
c$$$C       end
