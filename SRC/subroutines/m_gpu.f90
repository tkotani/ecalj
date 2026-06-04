module m_gpu
#ifdef __GPU
  use openacc
   use m_nvfortran, only: findloc
  use cudafor
#endif
  implicit none
  public :: gpu_init, check_memory_gpu, use_gpu, gpu_finalize, ngpu_ranks, ranks_per_gpu
  integer,public :: mydev
  logical, protected :: use_gpu = .false.
  integer, protected :: ngpu_ranks = 0    ! number of GPU ranks (= ndevs)
  integer, protected :: ranks_per_gpu = 1  ! MPI ranks sharing one GPU
  private
  integer :: procid, nsize
  contains
  
  subroutine gpu_init(comm, share_gpu)
    use iso_c_binding
    use mpi
    implicit none
    integer, intent(in) :: comm
    !> share_gpu=.true.: all local ranks share GPUs (round-robin); default .false. = only first ndevs ranks use GPU.
    logical, intent(in), optional :: share_gpu
    logical :: share_gpu_
    integer :: status(mpi_status_size)
    integer :: ierr, ndevs, ndevs_tmp, mydev_tmp, hostid_tmp, i, hostid, nlocal_procs, ilocal_rank
    integer, allocatable :: hostids(:), rankids(:)
    interface
      function gethostid() bind(c)
        use iso_c_binding
        integer (c_int) :: gethostid
      end function gethostid
    end interface

    share_gpu_ = .false.
    if (present(share_gpu)) share_gpu_ = share_gpu

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

    if (share_gpu_) then
      ! All local ranks share GPUs (round-robin); each initializes its assigned device.
      ngpu_ranks = nlocal_procs
      ranks_per_gpu = (nlocal_procs + ndevs - 1) / ndevs  ! ceiling(nlocal_procs/ndevs)
      use_gpu = .true.
    else
      ! Only first ndevs local ranks use GPU (one rank per GPU).
      ngpu_ranks = min(ndevs, nlocal_procs)
      ranks_per_gpu = 1
      use_gpu = (ilocal_rank < ndevs)
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
