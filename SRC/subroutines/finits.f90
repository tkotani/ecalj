subroutine fexit0(retval,strng)
  use m_MPItk,only:  m_MPItk_finalize
  use m_lgunit,only:stdo,stdl
  implicit none
  !i Inputs
  !i   retval:  return value passed to operating system
  !i   iopt decomposed into 3 one-digit numbers.
  !i   digit
  !i     1:  0: do not print string on exit;
  !i         9: print strng as Exit(retval): 'strng'
  !i      else: exit, using strn as a format statement and args a vector
  !i            of  c  double precision arguments
  !i   100:   0: do not print work array usage, else do
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
  integer :: master,procid,ierr,ia
  parameter (master = 0)
  include "mpif.h"
  call flush()
  ia=0
  goto 5
  entry fexit(retval,iopt,strng,args)
  ia=1
  goto 5
  entry fexit2(retval,iopt,strng,args,arg2)
  ia=2
  goto 5
  entry fexit3(retval,iopt,strng,args,arg2,arg3)
  ia=3
5 continue
  if(ia>=1) argss(1)=args
  if(ia>=2) argss(2)=arg2
  if(ia>=3) argss(3)=arg3
  procid = mpipid(1)
  if (procid == master) then
     if ( cpusec() /= 0) then
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
        !datim = ' '
        call ftime(datim)
        !          hostnm = ' '
        !          call get_environment_variable('HOST',hostnm)
        !write(stdo,*)'fexit0: exiting process=',procid
        write(stdo,"('CPU time:', f9.3,a1,5x,a,' on process=',i0)") &
             tnew,timeu,datim,procid!trim(adjustl(hostnm))
        !          write(stdl,10) tnew,timeu,datim,procid!trim(adjustl(hostnm))
     endif
  endif   !write(stdo,*)'finalizing procid=',procid
  if(procid==master) call tcprt(stdo)
  if(retval/=0) then
     write(stdo,"(a,i0,a,i0,a,3d15.8)") &
          "ERROR Exit ",retval,' procid= ',procid,' '//trim(strng),(argss(ix),ix=1,ia)
     !call MPI_abort(MPI_comm_world,retval,ierr)
  else
     if(procid==master)write(stdo,"(a,i0,a,3d15.8)")'Exit 0 procid= ',procid,' '//trim(strng),(argss(ix),ix=1,ia)
  endif
  call MPI_FINALIZE(ierr)
  call exit(retval)
end subroutine fexit0
subroutine rx0s(string) !normal exit for master_mpi
  use m_MPItk,only: master_mpi
  character*(*) string
  if(master_mpi) write(6,"(/,a)") trim(string)//' ======================'
  call exit(0)
end subroutine rx0s
subroutine rx(string) !error exit
  character*(*) string
  call flush()
  call fexit0(-1,string)
  call rx0s(string)         !for single core exit
end subroutine rx
subroutine rx0(string) !normal exit
  character*(*) string
  call flush()
  call fexit0(0,string)
  call rx0s(string)         !for single core exit
end subroutine rx0
subroutine rx1(string,arg) ! Error exit, with a single argument
  character(15):: f2a
  character*(*) string
  double precision :: arg
  call flush()
  call fexit0(-1,trim(' Exit -1 '//string//trim(f2a(arg))))
end subroutine rx1
subroutine rx2(string,arg1,arg2) ! Error exit, with two arguments
  character(15):: f2a
  character*(*) string
  double precision :: arg1,arg2
  call fexit0(-1,' Exit -1 '//string//trim(f2a(arg1))//' '//trim(f2a(arg2)))
end subroutine rx2
subroutine rxi(string,iarg) ! Error exit, with a single integer at end
  character*(*) string
  integer:: iarg
  character(10):: i2char
  call fexit0(-1,' Exit -1 '//string//' '//trim(i2char(iarg)))
end subroutine rxi
subroutine rxs(string,msg) ! Error exit with extra string message
  character*(*) string,msg
  character(120) :: outs
  integer :: i
  outs = string // msg
  call rx(trim(outs))
end subroutine rxs
subroutine rxs2(string,msg,msg2) ! Error exit with extra string messages
  character*(*) string,msg,msg2
  character(120) :: outs
  integer :: i
  outs = string // msg // msg2
  call rx(trim(outs))
end subroutine rxs2
subroutine rxs4(string,msg,msg2,msg3,msg4) ! Error exit with extra string messages
  character*(*) string,msg,msg2,msg3,msg4
  character(120) :: outs
  integer :: i
  outs = string // msg // msg2 // msg3 // msg4
  call skpblb(outs,len(outs),i)
  call rx(outs(1:i+1))
end subroutine rxs4
subroutine rxx(test,string)
  logical :: test
  character*(*) string
  if (test) call rx(string)
end subroutine rxx
subroutine rx_(string) !error exit routine
  character*(*) string
  write(6,"(/' ---- ',a)") string
  write(6,"(' ---- Error exit')")
  write(71,"('++ ',a)") string
  call exit(-1)
end subroutine rx_
subroutine cexit(pv,ps)
  implicit none
  include 'mpif.h'
  integer:: pv,ps,i
  integer:: status,ierr
  if (ps /= 0) then
     if (pv == 0) then
        call flush()
        call MPI_finalized(status,ierr)
        if (status == 0) then
           call MPI_finalize(ierr)
        endif
     endif
     call exit(pv)
  endif
end subroutine cexit

