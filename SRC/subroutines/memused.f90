! m_mem_node: system memory queries with no m_mpi dependency.
! m_mpi can freely use this module without a circular dependency.
module m_mem_node
  public mem_avail_node_gb, freeram, totalram
  private
contains
  ! Returns available node memory in GB: min(MemAvailable, cgroup limit).
  ! MemAvailable: includes reclaimable page cache, excludes swap.
  ! cgroup limit: enforced by SLURM/PBS/LSF job schedulers on HPC systems.
  real(8) function mem_avail_node_gb()
    integer :: funit, ios, c1, c2
    character(len=256) :: line, limitfile
    character(len=20)  :: key
    integer(8) :: val_kb, lim_bytes
    real(8) :: mem_avail, cg_lim

    mem_avail = 0d0
    open(newunit=funit, file='/proc/meminfo', status='old', iostat=ios)
    if (ios == 0) then
      do
        read(funit, '(A)', iostat=ios) line
        if (ios /= 0) exit
        read(line, *, iostat=ios) key, val_kb
        if (ios /= 0) cycle
        if (trim(key) == 'MemAvailable:') then
          mem_avail = real(val_kb, 8) / 1d6  ! kB → GB
          exit
        endif
      enddo
      close(funit)
    endif

    ! Find cgroup memory limit from /proc/self/cgroup
    cg_lim   = huge(1d0)
    limitfile = ''
    open(newunit=funit, file='/proc/self/cgroup', status='old', iostat=ios)
    if (ios == 0) then
      do
        read(funit, '(A)', iostat=ios) line
        if (ios /= 0) exit
        if (line(1:3) == '0::') then
          ! cgroups v2: "0::/<path>" → /sys/fs/cgroup/<path>/memory.max
          limitfile = '/sys/fs/cgroup' // trim(line(4:)) // '/memory.max'
          exit
        endif
        ! cgroups v1: find line whose subsystem list contains "memory"
        c1 = index(line, ':')
        c2 = index(line(c1+1:), ':') + c1
        if (c1 > 0 .and. c2 > c1 .and. index(line(c1+1:c2-1), 'memory') > 0) then
          limitfile = '/sys/fs/cgroup/memory' // trim(line(c2+1:)) &
                    // '/memory.limit_in_bytes'
          exit
        endif
      enddo
      close(funit)
    endif
    if (len_trim(limitfile) > 0) then
      open(newunit=funit, file=trim(limitfile), status='old', iostat=ios)
      if (ios == 0) then
        read(funit, '(A)', iostat=ios) line
        close(funit)
        line = adjustl(line)
        if (trim(line) /= 'max') then   ! 'max' = no cgroup limit (v2)
          read(line, *, iostat=ios) lim_bytes
          ! cgroups v1 "no limit" sentinel: value near INT64_MAX
          if (ios == 0 .and. lim_bytes > 0 .and. lim_bytes < 2_8**60) &
            cg_lim = real(lim_bytes, 8) / 1d9  ! bytes → GB
        endif
      endif
    endif

    mem_avail_node_gb = min(mem_avail, cg_lim)
  end function mem_avail_node_gb

  real(8) function freeram() !GB
    freeram = mem_avail_node_gb()
  end function freeram

  real(8) function totalram() !GB
    use iso_c_binding
    implicit none
    type, bind(C) :: t_sysinfo
       integer(c_long) :: uptime, loads(3), totalram, freeram, sharedram, bufferram
       integer(c_long) :: totalswap, freeswap
       integer(c_short) :: procs
       integer(c_long) :: totalhigh, freehigh
       integer(c_int) :: mem_unit
       character(c_char) :: f(20-2*sizeof(c_long)-sizeof(c_int))
    end type t_sysinfo
    interface
       integer function fsysinfo(info) bind(C, name="sysinfo")
         import :: t_sysinfo
         type(t_sysinfo), intent(out) :: info
       end function fsysinfo
    end interface
    type(t_sysinfo) :: info
    integer :: ret
    ret = fsysinfo(info)
    totalram = real(info%totalram, 8) * real(info%mem_unit, 8) / 1d9
  end function totalram
end module m_mem_node

module m_mem
  use m_mem_node, only: freeram, totalram, mem_avail_node_gb
  public writemem, memused, datetime, totalram, freeram, mem_avail_node_gb
  private
  real(8) :: mempeak=0d0
contains
  subroutine writemem(message)
    use m_ftox
    use m_lgunit,only:stdo
    use m_mpi,only: MPI__rank
    use m_gpu, only: mydev
#ifdef __GPU
    use openacc
    use cudafor
#endif
    character(*)  :: message
    character(len=128) :: memuse_gpu = ''
#ifdef __GPU
    real(8),parameter:: kk=1024,GG=kk**3 !afac for workspace of zhgv
    integer(8):: total_mem,free_mem, used_mem
    integer :: istat
    istat = cudaDeviceSynchronize()
    total_mem = acc_get_property(mydev, acc_device_nvidia, acc_property_memory)
    free_mem  = acc_get_property(mydev, acc_device_nvidia, acc_property_free_memory)
    used_mem  = total_mem - free_mem
    write(memuse_gpu,ftox,advance="no") ' (GPU)',ftof(dble(used_mem/GG),3),'GB'
#endif
    write(stdo,ftox)trim(message)//repeat(' ',mod(1000-len_trim(message),55))//&
         ' rank=',MPI__rank,'Memused (CPU)',ftof(memused(),3),'GB'//trim(memuse_gpu),datetime()

    flush(stdo)
    if(mempeak>memused()) mempeak=memused()
  end subroutine writemem

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

  character(23) function datetime()
    character(8)  :: date
    character(10) :: time
    character(5)  :: zone
    character(1)::sep='-',sp='-'
    integer,dimension(8) :: values
    call date_and_time(date,time,zone,values)
    call date_and_time(DATE=date,ZONE=zone)
    call date_and_time(TIME=time)
    datetime=date(1:4)//sp//date(5:6)//sp//date(7:8)//'T'//time(1:2)//':'//time(3:4)//':'//trim(time(5:10))
  end function datetime
end module m_mem
