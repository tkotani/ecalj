real(8) function memused() !in GB
  use iso_c_binding
  implicit none
  type, bind(C) :: c_timeval
     integer(c_long) :: tv_sec
     integer(c_long) :: tv_usec
  endtype c_timeval
  type, bind(C) :: c_rusage
     type(c_timeval) :: ru_utime
     type(c_timeval) :: ru_stime
     integer(c_long) :: ru_maxrss
     integer(c_long) :: ru_ixrss
     integer(c_long) :: ru_idrss
     integer(c_long) :: ru_isrss
     integer(c_long) :: ru_minflt
     integer(c_long) :: ru_majflt
     integer(c_long) :: ru_nswap
     integer(c_long) :: ru_inblock
     integer(c_long) :: ru_oublock
     integer(c_long) :: ru_msgsnd
     integer(c_long) :: ru_msgrcv
     integer(c_long) :: ru_nsignals
     integer(c_long) :: ru_nvcsw
     integer(c_long) :: ru_nivcsw
  end type c_rusage
  interface
     function getrusage(what,usage) bind(C, name="getrusage")
       import :: c_int, c_long, c_rusage
       integer(c_int) :: getrusage
       integer(c_int), value :: what
       type(c_rusage) :: usage
     end function getrusage
  end interface
  type(c_rusage) :: usage
  integer(c_int) :: ret
  integer :: mpi__info
  real(8)::k=1000
  ret = getrusage(0,usage)
  memused = usage%ru_maxrss/k**2 ! in GB
end function memused


character(18) function datetime()
  character(8)  :: date
  character(10) :: time
  character(5)  :: zone
  integer,dimension(8) :: values
  call date_and_time(date,time,zone,values)
  call date_and_time(DATE=date,ZONE=zone)
  call date_and_time(TIME=time)
  call date_and_time(VALUES=values)
  datetime=trim(date)//' '//trim(time)
  !print '(a,2x,a,2x,a)', date, time, zone
  !print '(8i5)', values
end function datetime

! real(8) function totalram() !GB
!   use iso_c_binding
!   implicit none
!   type, bind(C) :: sysinfo
!   integer(c_long) :: uptime
!   integer(c_long) :: loads(3)
!   integer(c_long) :: totalram
!   integer(c_long) :: freeram
!   integer(c_long) :: sharedram
!   integer(c_long) :: bufferram
!   integer(c_long) :: totalswap
!   integer(c_long) :: freeswap
!   integer(c_short) :: procs
!   integer(c_long) :: totalhigh
!   integer(c_long) :: freehigh
!   integer(c_int) :: mem_unit
!   character(c_char)::f(20-2*sizeof(c_long)-sizeof(c_int))
!   end type sysinfo
!   interface
!      function fsysinfo(info) bind(C, name="sysinfo")
!        import :: sysinfo
!        type(sysinfo), intent(out) :: info
!      end function fsysinfo
!   end interface
!   type(sysinfo) :: info
!   integer(c_int) :: ret
!   real(8)::k=1000
!   ret = fsysinfo(info)
!   totalram=info%totalram/k**3
!   ! print *, 'ret=', ret
!   ! print *, 'uptime=', info%uptime
!   ! print *, 'loads=', info%loads
!   ! print *, 'totalram=', info%totalram
!   ! print *, 'freeram=', info%freeram
!   ! print *, 'sharedram=', info%sharedram
!   ! print *, 'bufferram=', info%bufferram
!   ! print *, 'totalswap=', info%totalswap
!   ! print *, 'freeswap=', info%freeswap
!   ! print *, 'procs=', info%procs
!   ! print *, 'totalhigh=', info%totalhigh
!   ! print *, 'freehigh=', info%freehigh
!   ! print *, 'mem_unit=', info%mem_unit
! end function totalram

! program memory_usage !GB
!     implicit none
!     include "mpif.h"
!     integer:: k=1000,ierr,comm,procid,i
!     real(8), allocatable :: aaa(:),bbb(:),ccc(:)
!     real(8) :: memused,totalram,x,timei,timee,diff
!     real:: start,end
!     integer date_timei(8),date_timee(8),ddd(8)
!     character*10 b(3)
!     call MPI_INIT(ierr)
!     comm=MPI_COMM_WORLD
!     !call mpi_comm_size(comm, nsize, ierr)
!     call MPI_COMM_RANK(comm, procid, ierr )
!     do i=1,38!       call cpu_time(start)
!        call date_and_time(b(1), b(2), b(3), date_timei)
! !       timei=time()
!        allocate(aaa(i*100*k**2),source=1d0)
!        x= sum(aaa)
! !       timee=time()
! !       call cpu_time(end)
!        call date_and_time(b(1), b(2), b(3), date_timee)
!        ddd=date_timee-date_timei
!        diff= merge(ddd(7),60+ddd(7),ddd(7)>0)+ ddd(6)*60
!        write(6,"('procid imem',2i4,' Memuse1',f10.3,' GB RAM=',f10.3,' GB time=',f8.2)") procid,i,memused(),totalram(),diff
!        flush(6)
!        deallocate(aaa)
!     enddo
!     call MPI_FINALIZE(ierr)
! end program memory_usage
