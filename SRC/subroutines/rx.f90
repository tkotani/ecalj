! All exit routines with errors 2023 feb. Too much varieties is historical reasons.
subroutine rx(string) !error exit
  use, intrinsic :: iso_fortran_env, only: error_unit
  use m_lgunit,only:stdo,stdl
  use m_mpi,only:   comm
  integer::ierr
  character*(*) string
  write(stdo,*) trim(string)
!  write(error_unit,*) trim(string)
  flush(stdo)
  call MPI_Abort(COMM, 11,ierr)
end subroutine rx
subroutine rxi(string,iarg) ! Error exit, with a single integer at end
  use, intrinsic :: iso_fortran_env, only: error_unit
  use m_lgunit,only:stdo,stdl
  use m_mpi,only:   comm
  character*(*) string
  integer:: iarg,ierr
  character(10):: i2char
  write(stdo,*)trim(' Exit -1 '//string//' '//trim(i2char(iarg)))
  write(error_unit,*)trim(' Exit -1 '//string//' '//trim(i2char(iarg)))
  call MPI_Abort(COMM, 21,ierr)
end subroutine rxi
subroutine rxii(string,iarg,iarg2) 
  use, intrinsic :: iso_fortran_env, only: error_unit
  use m_lgunit,only:stdo,stdl
  use m_mpi,only:   comm
  character*(*) string
  integer:: iarg,iarg2,ierr
  character(10):: i2char
  write(stdo,*)trim(' Exit -1 '//string//' '//trim(i2char(iarg))//' '//trim(i2char(iarg2)))
  write(error_unit,*)trim(' Exit -1 '//string//' '//trim(i2char(iarg))//' '//trim(i2char(iarg2)))
  call MPI_Abort(COMM, 31,ierr)
end subroutine rxii
subroutine rxiii(string,iarg,iarg2,iarg3) 
  use, intrinsic :: iso_fortran_env, only: error_unit
  use m_lgunit,only:stdo,stdl
  use m_mpi,only:   comm
  character*(*) string
  integer:: iarg,iarg2,iarg3,ierr
  character(10):: i2char
  write(stdo,*)trim(' Exit -1 '//string//' '//trim(i2char(iarg))//' '//trim(i2char(iarg2))//' '//trim(i2char(iarg3)))
  write(error_unit,*)trim(' Exit -1 '//string//' '//trim(i2char(iarg))//' '//trim(i2char(iarg2))//' '//trim(i2char(iarg3)))
  call MPI_Abort(COMM, 41,ierr)
end subroutine rxiii
subroutine rx1(string,arg) ! Error exit, with a single argument
  use, intrinsic :: iso_fortran_env, only: error_unit
  use m_ftox
  use m_lgunit,only:stdo,stdl
  use m_mpi,only:   comm
  integer:: ierr
  character(15):: f2a
  character*(*) string
  double precision :: arg
  write(stdo,*)trim(' Exit -1 '//string//trim(ftof(arg)))
  write(error_unit,*)trim(' Exit -1 '//string//trim(ftof(arg)))
  call MPI_Abort(COMM, 51,ierr)
end subroutine rx1
subroutine rx2(string,arg1,arg2) ! Error exit, with two arguments
  use, intrinsic :: iso_fortran_env, only: error_unit
  use m_lgunit,only:stdo,stdl
  use m_ftox
  use m_mpi,only:   comm
  integer::ierr
  character(15):: f2a
  character*(*) string
  double precision :: arg1,arg2
  write(stdo,*)trim(' Exit -1 '//string//trim(ftof(arg1))//' '//trim(ftof(arg2)))
  write(error_unit,*)trim(' Exit -1 '//string//trim(ftof(arg1))//' '//trim(ftof(arg2)))
  call MPI_Abort(COMM, 61,ierr)
end subroutine rx2
subroutine rx3(string,arg1,arg2,arg3) ! Error exit, with two arguments
  use, intrinsic :: iso_fortran_env, only: error_unit
  use m_lgunit,only:stdo,stdl
  use m_ftox
  use m_mpi,only:   comm
  integer:: ierr
  character(15):: f2a
  character*(*) string
  real(8):: arg1,arg2,arg3
  write(stdo,*)trim(' Exit -1 '//string//trim(ftof(arg1))//' '//trim(ftof(arg2))//' '//trim(ftof(arg3)))
  write(error_unit,*)trim(' Exit -1 '//string//trim(ftof(arg1))//' '//trim(ftof(arg2))//' '//trim(ftof(arg3)))
  call MPI_Abort(COMM, 71,ierr)
end subroutine rx3
subroutine rxs(string,msg) ! Error exit with extra string message
  character*(*) string,msg
  character(120) :: outs
  integer :: i
  outs = string // msg
  call rx(trim(outs))
end subroutine rxs
subroutine rxx(test,string)
  logical :: test
  character*(*) string
  if (test) call rx(string)
end subroutine rxx
!----------- Normal exit
subroutine rx0s(string) !normal exit for master_mpi
  use m_MPItk,only: master_mpi
  character*(*) string
  if(master_mpi) write(6,"(/,a)") trim(string)//' ======================'
  call exit(0)
end subroutine rx0s
subroutine rx0(strng)! Normal exit
  use m_mpi,only:   comm1=>comm
  use m_MPItk,only: comm2=>comm,readtk
  use m_lgunit,only:stdo,stdl
  implicit none
  integer :: iopt,abret
  character*(*) strng
  double precision :: args,arg2,arg3,argss(3)
  integer :: iprint,i,i2,scrwid,ix
  double precision :: cpusec,tnew
  character(1) :: timeu
  character(256) :: strn
  character :: datim*24,hostnm*20
  character(9):: ftoa9
  logical :: isopen,master_mpi
  integer :: master,ierr,comm,procid,mpi__info
  parameter (master = 0)
  include "mpif.h"
  comm=merge(comm2,comm1,readtk)
  call MPI_Comm_rank( comm, procid, mpi__info )
  master_mpi= procid==0
  call mpi_barrier(comm,ierr)
  if (master_mpi) then 
     if(cpusec() /= 0) then
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
  endif
  if(procid==master) write(stdo,"(a,i0,a,3d15.8)") 'Exit 0 procid= ',procid,' '//trim(strng)
  call MPI_FINALIZE(ierr)
  call exit(0)
end subroutine rx0
