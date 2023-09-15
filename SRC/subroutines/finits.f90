! All exit routines with errors 2023 feb. Too much varieties is historical reasons.
subroutine rx0s(string) !normal exit for master_mpi
  use m_MPItk,only: master_mpi
  character*(*) string
  if(master_mpi) write(6,"(/,a)") trim(string)//' ======================'
  call exit(0)
end subroutine rx0s
subroutine rx0(string) !normal exit
  use m_lgunit,only:stdo,stdl
  character*(*) string
  call fexit0(0,string)
  call rx0s(string)         !for single core exit
end subroutine rx0

subroutine rx(string) !error exit
  use m_lgunit,only:stdo,stdl
  character*(*) string
  call fexit0(-1,string)
  call rx0s(string)         !for single core exit
end subroutine rx
subroutine rxi(string,iarg) ! Error exit, with a single integer at end
  character*(*) string
  integer:: iarg
  character(10):: i2char
  call fexit0(-1,' Exit -1 '//string//' '//trim(i2char(iarg)))
end subroutine rxi
subroutine rxii(string,iarg,iarg2) 
  character*(*) string
  integer:: iarg,iarg2
  character(10):: i2char
  call fexit0(-1,' Exit -1 '//string//' '//trim(i2char(iarg))//' '//trim(i2char(iarg2)))
end subroutine rxii
subroutine rxiii(string,iarg,iarg2,iarg3) 
  character*(*) string
  integer:: iarg,iarg2,iarg3
  character(10):: i2char
  call fexit0(-1,' Exit -1 '//string//' '//trim(i2char(iarg))//' '//trim(i2char(iarg2))//' '//trim(i2char(iarg3)))
end subroutine rxiii

subroutine rx1(string,arg) ! Error exit, with a single argument
  use m_ftox
  use m_lgunit,only:stdo,stdl
  character(15):: f2a
  character*(*) string
  double precision :: arg
  call fexit0(-1,' Exit -1 '//string//trim(ftof(arg)))
end subroutine rx1
subroutine rx2(string,arg1,arg2) ! Error exit, with two arguments
  use m_ftox
  character(15):: f2a
  character*(*) string
  double precision :: arg1,arg2
  call fexit0(-1,' Exit -1 '//string//trim(ftof(arg1))//' '//trim(ftof(arg2)))
end subroutine rx2
subroutine rx3(string,arg1,arg2,arg3) ! Error exit, with two arguments
  use m_ftox
  character(15):: f2a
  character*(*) string
  real(8):: arg1,arg2,arg3
  call fexit0(-1,' Exit -1 '//string//trim(ftof(arg1))//' '//trim(ftof(arg2))//' '//trim(ftof(arg3)))
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

subroutine fexit0(retval,strng)! retval:  return value passed to operating system. /=0 for error exit
  use m_MPItk,only:  m_MPItk_finalize
  use m_lgunit,only:stdo,stdl
  implicit none
  integer :: retval,iopt,abret
  character*(*) strng
  double precision :: args,arg2,arg3,argss(3)
  integer :: iprint,i,i2,scrwid,mpipid,ix
  double precision :: cpusec,tnew
  character(1) :: timeu
  character(256) :: strn
  character :: datim*24,hostnm*20
  character(9):: ftoa9
  logical :: isopen
  integer :: master,procid,ierr
  parameter (master = 0)
  include "mpif.h"
  procid = mpipid(1)
  if (procid == master) then
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
  endif
  if(procid==master) call tcprt(stdo)
  if(retval/=0) then
     write(stdo,"(a,i0,a,i0,a,3d15.8)")"ERROR Exit ",retval,' procid= ',procid,' '//trim(strng)
     call MPI_Abort(MPI_COMM_WORLD, 911)
  else
     if(procid==master) write(stdo,"(a,i0,a,3d15.8)") 'Exit 0 procid= ',procid,' '//trim(strng)
  endif
  call MPI_FINALIZE(ierr)
  call exit(retval)
end subroutine fexit0
