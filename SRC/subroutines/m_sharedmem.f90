!> MPI shared memory manager for intra-node data sharing.
!! All ranks on the same node share the same physical memory via MPI-3 shared memory windows.
!! Rank 0 in the shared comm allocates; others attach via MPI_Win_shared_query.
module m_sharedmem
  use m_mpi, only: procid, numprocs=>nsize
  implicit none
  public :: shm_init, shm_finalize, shm_barrier
  public :: shm_alloc_r8_1d, shm_alloc_r8_2d, shm_alloc_r8_3d
  public :: shm_alloc_i4_1d, shm_alloc_i4_2d, shm_alloc_i4_3d
  public :: shm_alloc_c8_4d
  public :: shm_free
  public :: shm_comm, shm_rank, shm_nprocs
  private

  integer, parameter :: SHM_MAX_WIN = 32
  integer, save :: comm_shared = -1
  integer, save :: rank_shared = -1
  integer, save :: nprocs_shared = -1
  logical, save :: initialized = .false.
  integer, save :: win_handles(SHM_MAX_WIN) = -1
  logical, save :: win_used(SHM_MAX_WIN) = .false.

contains

  subroutine shm_init()
    use mpi
    implicit none
    integer :: ierr
    if(initialized) return
    call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, comm_shared, ierr)
    call MPI_Comm_rank(comm_shared, rank_shared, ierr)
    call MPI_Comm_size(comm_shared, nprocs_shared, ierr)
    initialized = .true.
  end subroutine

  function shm_comm() result(c); integer :: c; c = comm_shared; end function
  function shm_rank() result(r); integer :: r; r = rank_shared; end function
  function shm_nprocs() result(n); integer :: n; n = nprocs_shared; end function

  subroutine shm_alloc_r8_1d(ptr, n1, id)
    use mpi; use, intrinsic :: iso_c_binding
    implicit none
    real(8), pointer, intent(out) :: ptr(:)
    integer, intent(in) :: n1; integer, intent(out) :: id
    integer(MPI_ADDRESS_KIND) :: sz, szq; integer :: du,ierr,win; type(c_ptr) :: bp
    call shm_get_slot(id); du=8; sz=0
    if(rank_shared==0) sz=int(n1,MPI_ADDRESS_KIND)*8
    call MPI_Win_allocate_shared(sz,du,MPI_INFO_NULL,comm_shared,bp,win,ierr)
    if(rank_shared/=0) call MPI_Win_shared_query(win,0,szq,du,bp,ierr)
    call c_f_pointer(bp,ptr,[n1]); win_handles(id)=win
  end subroutine

  subroutine shm_alloc_r8_2d(ptr, n1, n2, id)
    use mpi; use, intrinsic :: iso_c_binding
    implicit none
    real(8), pointer, intent(out) :: ptr(:,:)
    integer, intent(in) :: n1,n2; integer, intent(out) :: id
    integer(MPI_ADDRESS_KIND) :: sz, szq; integer :: du,ierr,win; type(c_ptr) :: bp
    call shm_get_slot(id); du=8; sz=0
    if(rank_shared==0) sz=int(n1,MPI_ADDRESS_KIND)*int(n2,MPI_ADDRESS_KIND)*8
    call MPI_Win_allocate_shared(sz,du,MPI_INFO_NULL,comm_shared,bp,win,ierr)
    if(rank_shared/=0) call MPI_Win_shared_query(win,0,szq,du,bp,ierr)
    call c_f_pointer(bp,ptr,[n1,n2]); win_handles(id)=win
  end subroutine

  subroutine shm_alloc_r8_3d(ptr, n1, n2, n3, id)
    use mpi; use, intrinsic :: iso_c_binding
    implicit none
    real(8), pointer, intent(out) :: ptr(:,:,:)
    integer, intent(in) :: n1,n2,n3; integer, intent(out) :: id
    integer(MPI_ADDRESS_KIND) :: sz, szq; integer :: du,ierr,win; type(c_ptr) :: bp
    call shm_get_slot(id); du=8; sz=0
    if(rank_shared==0) sz=int(n1,MPI_ADDRESS_KIND)*int(n2,MPI_ADDRESS_KIND)*int(n3,MPI_ADDRESS_KIND)*8
    call MPI_Win_allocate_shared(sz,du,MPI_INFO_NULL,comm_shared,bp,win,ierr)
    if(rank_shared/=0) call MPI_Win_shared_query(win,0,szq,du,bp,ierr)
    call c_f_pointer(bp,ptr,[n1,n2,n3]); win_handles(id)=win
  end subroutine

  subroutine shm_alloc_i4_1d(ptr, n1, id)
    use mpi; use, intrinsic :: iso_c_binding
    implicit none
    integer, pointer, intent(out) :: ptr(:)
    integer, intent(in) :: n1; integer, intent(out) :: id
    integer(MPI_ADDRESS_KIND) :: sz, szq; integer :: du,ierr,win; type(c_ptr) :: bp
    call shm_get_slot(id); du=4; sz=0
    if(rank_shared==0) sz=int(n1,MPI_ADDRESS_KIND)*4
    call MPI_Win_allocate_shared(sz,du,MPI_INFO_NULL,comm_shared,bp,win,ierr)
    if(rank_shared/=0) call MPI_Win_shared_query(win,0,szq,du,bp,ierr)
    call c_f_pointer(bp,ptr,[n1]); win_handles(id)=win
  end subroutine

  subroutine shm_alloc_i4_2d(ptr, n1, n2, id)
    use mpi; use, intrinsic :: iso_c_binding
    implicit none
    integer, pointer, intent(out) :: ptr(:,:)
    integer, intent(in) :: n1,n2; integer, intent(out) :: id
    integer(MPI_ADDRESS_KIND) :: sz, szq; integer :: du,ierr,win; type(c_ptr) :: bp
    call shm_get_slot(id); du=4; sz=0
    if(rank_shared==0) sz=int(n1,MPI_ADDRESS_KIND)*int(n2,MPI_ADDRESS_KIND)*4
    call MPI_Win_allocate_shared(sz,du,MPI_INFO_NULL,comm_shared,bp,win,ierr)
    if(rank_shared/=0) call MPI_Win_shared_query(win,0,szq,du,bp,ierr)
    call c_f_pointer(bp,ptr,[n1,n2]); win_handles(id)=win
  end subroutine

  subroutine shm_alloc_i4_3d(ptr, n1, n2, n3, id)
    use mpi; use, intrinsic :: iso_c_binding
    implicit none
    integer, pointer, intent(out) :: ptr(:,:,:)
    integer, intent(in) :: n1,n2,n3; integer, intent(out) :: id
    integer(MPI_ADDRESS_KIND) :: sz, szq; integer :: du,ierr,win; type(c_ptr) :: bp
    call shm_get_slot(id); du=4; sz=0
    if(rank_shared==0) sz=int(n1,MPI_ADDRESS_KIND)*int(n2,MPI_ADDRESS_KIND)*int(n3,MPI_ADDRESS_KIND)*4
    call MPI_Win_allocate_shared(sz,du,MPI_INFO_NULL,comm_shared,bp,win,ierr)
    if(rank_shared/=0) call MPI_Win_shared_query(win,0,szq,du,bp,ierr)
    call c_f_pointer(bp,ptr,[n1,n2,n3]); win_handles(id)=win
  end subroutine

  subroutine shm_alloc_c8_4d(ptr, n1, n2, n3, n4, id)
    use mpi; use, intrinsic :: iso_c_binding
    implicit none
    complex(8), pointer, intent(out) :: ptr(:,:,:,:)
    integer, intent(in) :: n1,n2,n3,n4; integer, intent(out) :: id
    integer(MPI_ADDRESS_KIND) :: sz, szq; integer :: du,ierr,win; type(c_ptr) :: bp
    call shm_get_slot(id); du=16; sz=0
    if(rank_shared==0) sz=int(n1,MPI_ADDRESS_KIND)*int(n2,MPI_ADDRESS_KIND) &
                          *int(n3,MPI_ADDRESS_KIND)*int(n4,MPI_ADDRESS_KIND)*16
    call MPI_Win_allocate_shared(sz,du,MPI_INFO_NULL,comm_shared,bp,win,ierr)
    if(rank_shared/=0) call MPI_Win_shared_query(win,0,szq,du,bp,ierr)
    call c_f_pointer(bp,ptr,[n1,n2,n3,n4]); win_handles(id)=win
  end subroutine

  subroutine shm_barrier(id)
    use mpi
    implicit none
    integer, intent(in) :: id; integer :: ierr
    call MPI_Win_fence(0, win_handles(id), ierr)
  end subroutine

  subroutine shm_free(id)
    use mpi
    implicit none
    integer, intent(in) :: id; integer :: ierr
    if(id<1.or.id>SHM_MAX_WIN) return
    if(.not.win_used(id)) return
    call MPI_Win_free(win_handles(id), ierr)
    win_handles(id)=-1; win_used(id)=.false.
  end subroutine

  subroutine shm_finalize()
    use mpi
    implicit none
    integer :: i, ierr
    do i=1,SHM_MAX_WIN; if(win_used(i)) call shm_free(i); enddo
    if(comm_shared/=-1) call MPI_Comm_free(comm_shared, ierr)
    comm_shared=-1; initialized=.false.
  end subroutine

  subroutine shm_get_slot(id)
    implicit none
    integer, intent(out) :: id; integer :: i
    do i=1,SHM_MAX_WIN
      if(.not.win_used(i)) then; id=i; win_used(i)=.true.; return; endif
    enddo
    call rx('m_sharedmem: no free slots')
  end subroutine

end module m_sharedmem
