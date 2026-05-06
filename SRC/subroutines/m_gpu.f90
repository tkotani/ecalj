module m_gpu
#ifdef __GPU
  use openacc
   use m_nvfortran, only: findloc
  use cudafor
#endif
  implicit none
  public :: gpu_init, check_memory_gpu, use_gpu, gpu_finalize, ngpu_ranks
  integer,public :: mydev
  logical, protected :: use_gpu = .false.
  integer, protected :: ngpu_ranks = 0  ! number of GPU ranks (= ndevs)
  private
  integer :: procid, nsize
  contains
  
  subroutine gpu_init(comm)
    use iso_c_binding
    use mpi
    implicit none
    integer, intent(in) :: comm
    integer :: status(mpi_status_size)
    integer :: ierr, ndevs, ndevs_tmp, mydev_tmp, hostid_tmp, i, hostid, nlocal_procs, ilocal_rank
    integer, allocatable :: hostids(:), rankids(:)
    logical :: cmdopt0
    interface
      function gethostid() bind(c)
        use iso_c_binding
        integer (c_int) :: gethostid
      end function gethostid
    end interface

#ifdef __GPU
    call mpi_comm_rank(comm, procid, ierr)
    call mpi_comm_size(comm, nsize, ierr)
    allocate(hostids(nsize), source = 0)
    allocate(rankids(nsize), source = 0)
    hostid = gethostid()
    call mpi_allgather(hostid, 1, mpi_integer, hostids, 1, mpi_integer, comm, ierr)
    call mpi_allgather(procid, 1, mpi_integer, rankids, 1, mpi_integer, comm, ierr)

    nlocal_procs = size(pack(hostids, hostids == hostid))
    ilocal_rank = findloc(pack(rankids, hostids == hostid), procid, dim=1) - 1

    ndevs = acc_get_num_devices(acc_device_nvidia)
    if(ndevs == 0) return
    mydev = mod(ilocal_rank, ndevs)

    ! Automatic GPU assignment: only first ndevs local ranks use GPU
    ! e.g. 2 GPUs → rank 0=GPU0, rank 1=GPU1, rank 2+=CPU only
    ngpu_ranks = min(ndevs, nlocal_procs)
    if(ilocal_rank >= ndevs) then
      use_gpu = .false.
    else
      use_gpu = .true.
    endif

    if(use_gpu) then
      call acc_set_device_num(mydev, acc_device_nvidia)
      call acc_init(acc_device_nvidia)
    endif
    ! call check_memory_gpu("gpu_init")

    if (procid == 0) then
      write(06,'(a,i6,x,2(a,i3),a,i12)') "i_procs:", procid, "gpuid:", mydev, "/", ndevs, " hostid:", hostid
    endif
    do i = 1, nsize-1
      if(procid == 0) then
        call mpi_recv(mydev_tmp, 1, mpi_integer, i, i, comm, status, ierr)
        call mpi_recv(ndevs_tmp, 1, mpi_integer, i, i, comm, status, ierr)
        call mpi_recv(hostid_tmp, 1, mpi_integer, i, i, comm, status, ierr)
        write(06,'(a,i6,x,2(a,i3),a,i12)') "i_procs:", i, "gpuid:", mydev_tmp, "/", ndevs_tmp, " hostid:", hostid_tmp
      elseif(procid == i) then
        call mpi_send(mydev, 1, mpi_integer, 0, procid, comm, ierr)
        call mpi_send(ndevs, 1, mpi_integer, 0, procid, comm, ierr)
        call mpi_send(hostid, 1, mpi_integer, 0, procid, comm, ierr)
      endif
      call mpi_barrier(comm, ierr)
    enddo
#endif
  end subroutine

  subroutine check_memory_gpu(keyword)
    character(len=*), intent(in) :: keyword
    character(len=1024) :: cmd
    integer :: ierr
#ifdef __GPU
    cmd = 'echo gpu_memory_check ' // trim(keyword)
    call execute_command_line(cmd)
    ierr = cudadevicesynchronize() 
    cmd = 'nvidia-smi --query-gpu=name,utilization.memory,memory.total,memory.free,memory.used --format=csv'
    call execute_command_line(cmd)
#endif
  end subroutine

  subroutine gpu_finalize()
#ifdef __GPU
    use m_blas, only: cublas_finalize
    use m_lapack, only: cusolver_finalize
    integer :: istat
    istat = cusolver_finalize()
    istat = cublas_finalize()
    call acc_shutdown(acc_device_nvidia)
#endif
    use_gpu = .false.
  end subroutine gpu_finalize
 !  real(8) function Amem_gpu() !Available memory in GPU 
! #ifdef __GPU
!     use cudafor
! #endif
!     character(len=1024) :: cmd
!     integer :: ierr,ifi,ndev,i
!     character(256):: temp
!     character(8),external:: xt
!     real(8)::total_mem,free_mem
!     Amem_GPU =1d10
! #ifdef __GPU
!     ierr = cudadevicesynchronize()
!     ierr = cudaMemGetInfo(Amem_gpu, total_mem)
! #endif
!   end function Amem_gpu
 
end module m_gpu
