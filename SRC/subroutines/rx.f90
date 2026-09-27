! All exit routines with errors 2023 feb. Too much varieties is historical reasons.
! 2026-09-27 22:50: every error exit goes through rx_stop.  The message goes to stdout (the log of rank 0; the other
! ranks write stdout to /dev/null unless --fullstdo, see MPI__consoleout) and to stderr with the rank; then MPI_Abort
! stops all ranks, or error stop when MPI is not running.  Bugs fixed then: rx wrote to stdout only, rxs cut the
! message at 120 characters, and the normal exits rx0s (no MPI_Finalize: mpirun waited or reported an abnormal exit)
! and rx0 (m_mpi's comm, unset in programs that call mpi_init themselves) did not handle MPI properly.
subroutine rx_stop(msg, code) !error exit of all rx* below: message to stdout and stderr, then stop every rank
  use, intrinsic :: iso_fortran_env, only: error_unit
  use m_lgunit,only: stdo
  use m_mpi,only: procid
  use mpi
  implicit none
  character(*), intent(in) :: msg
  integer, intent(in) :: code
  logical :: inited, finalized
  integer :: ierr
  write(stdo,*) trim(msg)
  flush(stdo)
  write(error_unit,'(a,i0,a)') ' rank ',procid,': '//trim(msg)
  flush(error_unit)
  call MPI_Initialized(inited, ierr)
  call MPI_Finalized(finalized, ierr)
  if (inited .and. .not. finalized) call MPI_Abort(MPI_COMM_WORLD, code, ierr)
  call exit(code)
end subroutine rx_stop
subroutine rx(string) !error exit
  character*(*) string
  call rx_stop(string, 11)
end subroutine rx
subroutine rxi(string,iarg) ! Error exit, with a single integer at end
  character*(*) string
  integer:: iarg
  character(10):: i2char
  call rx_stop(' Exit -1 '//string//' '//trim(i2char(iarg)), 21)
end subroutine rxi
subroutine rxii(string,iarg,iarg2)
  character*(*) string
  integer:: iarg,iarg2
  character(10):: i2char
  call rx_stop(' Exit -1 '//string//' '//trim(i2char(iarg))//' '//trim(i2char(iarg2)), 31)
end subroutine rxii
subroutine rxiii(string,iarg,iarg2,iarg3)
  character*(*) string
  integer:: iarg,iarg2,iarg3
  character(10):: i2char
  call rx_stop(' Exit -1 '//string//' '//trim(i2char(iarg))//' '//trim(i2char(iarg2))//' '//trim(i2char(iarg3)), 41)
end subroutine rxiii
subroutine rx1(string,arg) ! Error exit, with a single argument
  use m_ftox
  character*(*) string
  double precision :: arg
  call rx_stop(' Exit -1 '//string//trim(ftof(arg)), 51)
end subroutine rx1
subroutine rx2(string,arg1,arg2) ! Error exit, with two arguments
  use m_ftox
  character*(*) string
  double precision :: arg1,arg2
  call rx_stop(' Exit -1 '//string//trim(ftof(arg1))//' '//trim(ftof(arg2)), 61)
end subroutine rx2
subroutine rx3(string,arg1,arg2,arg3) ! Error exit, with three arguments
  use m_ftox
  character*(*) string
  real(8):: arg1,arg2,arg3
  call rx_stop(' Exit -1 '//string//trim(ftof(arg1))//' '//trim(ftof(arg2))//' '//trim(ftof(arg3)), 71)
end subroutine rx3
subroutine rxs(string,msg) ! Error exit with extra string message
  character*(*) string,msg
  call rx(string//msg)
end subroutine rxs
subroutine rxx(test,string)
  logical :: test
  character*(*) string
  if (test) call rx(string)
end subroutine rxx
!----------- Normal exit
subroutine rx0s(string) !normal exit without the timing report: message by rank 0, then MPI_Finalize (all ranks call it)
  use mpi
  implicit none
  character*(*) string
  logical :: inited, finalized
  integer :: ierr, rank
  call MPI_Initialized(inited, ierr)
  call MPI_Finalized(finalized, ierr)
  rank = 0
  if (inited .and. .not. finalized) call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  if (rank == 0) write(6,"(/,a)") trim(string)//' ======================'
  if (inited .and. .not. finalized) call MPI_Finalize(ierr)
  call exit(0)
end subroutine rx0s
subroutine rx0(strng)! Normal exit: CPU time and the timing report by rank 0, then MPI_Finalize (all ranks call it)
  use m_mpi,only: comm, readtk
  use m_lgunit,only: stdo
  use mpi
  implicit none
  character*(*) strng
  double precision :: cpusec, tnew
  character(1) :: timeu
  character :: datim*24
  logical :: inited, finalized
  integer :: ierr, procid, comm_
  call MPI_Initialized(inited, ierr)
  call MPI_Finalized(finalized, ierr)
  procid = 0
  if (inited .and. .not. finalized) then
    comm_ = merge(comm, MPI_COMM_WORLD, readtk)   ! m_mpi's comm when MPI__Initialize set it
    call MPI_Comm_rank(comm_, procid, ierr)
    call MPI_Barrier(comm_, ierr)
  endif
  if (procid == 0) then
    if (cpusec() /= 0) then
      timeu = 's'
      tnew = cpusec()
      if (tnew > 3600) then
        timeu = 'm'
        tnew = tnew/60
        if (tnew > 200) then
          timeu = 'h'
          tnew = tnew/60
        endif
      endif
      call ftime(datim)
      write(stdo,"('CPU time:', f9.3,a1,5x,a,' on process=',i0)") tnew,timeu,datim,procid
    endif
    call tcprt(stdo)
    write(stdo,"(a,i0,a)") 'Exit 0 procid= ',procid,' '//trim(strng)
  endif
  if (inited .and. .not. finalized) call MPI_Finalize(ierr)
  call exit(0)
end subroutine rx0
