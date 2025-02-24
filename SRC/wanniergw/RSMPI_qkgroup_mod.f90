module rsmpi_qkgroup
  use rsmpi,only:mpi_comm_world,myrank_rsmpi,nproc_rsmpi
  implicit none
  public rsmpi_qkgroup_init
  integer,public :: nk_local_qkgroup ! number of k-points treated in the process
  integer,allocatable,public :: ik_index_qkgroup(:)

  private
  integer :: nqkgroup ! number of subgroups
  ! it is equivalent to the number of q-points       calculated at the same time
!  integer :: iqkgroup ! the subgroup the current process belongs
  integer :: nproc_qkgroup    ! number of processes in the subgroup
  integer :: myrank_qkgroup   ! rank of the calling process in the subgroup
!  logical :: file_io_qkgroup  ! The processes with this value .true.
  integer :: nq_local_qkgroup ! number of q-points treated in the subgroup
  integer,allocatable :: iq_index_qkgroup(:) ! index
  integer :: ierror_qkgroup   ! error check
  integer :: comm_qkgroup ! new communicator
  character(7) :: qkgroup_id   !
  integer :: ifile_qkgroup    ! file id
contains
  subroutine RSMPI_qkgroup_Init(Nq,Nk)
    implicit none
    integer,intent(in) :: Nq,Nk
    integer:: iqkgroup
    nqkgroup = merge(nq,nproc_rsmpi,Nq <= nproc_rsmpi)
    iqkgroup = get_my_qkgroup(myrank_rsmpi,nproc_rsmpi)
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,iqkgroup, myrank_rsmpi,comm_qkgroup,ierror_qkgroup)
    call MPI_COMM_SIZE(comm_qkgroup,nproc_qkgroup,ierror_qkgroup)
    call MPI_COMM_RANK(comm_qkgroup,myrank_qkgroup,ierror_qkgroup)
!    call set_q_local_qkgroup(Nq)
    call set_k_local_qkgroup(Nk)
  end subroutine RSMPI_qkgroup_Init
  subroutine set_k_local_qkgroup(Nk)
    implicit none
    integer,intent(in) :: Nk ! total number of k-points in full BZ
    integer,allocatable :: nk_local_all(:),ik_index_all(:,:)
    allocate(nk_local_all(nproc_qkgroup),      ik_index_all(nproc_qkgroup,Nk/nproc_qkgroup+1))
    call set_index_rsmpi(Nk,nproc_qkgroup,nk_local_all,ik_index_all)
    nk_local_qkgroup = nk_local_all(myrank_qkgroup+1) ! rank is 0,1,..,N-1
    allocate(ik_index_qkgroup(nk_local_qkgroup))
    ik_index_qkgroup(1:nk_local_qkgroup) =         ik_index_all(myrank_qkgroup+1,1:nk_local_qkgroup)
    deallocate(nk_local_all,ik_index_all)
  end subroutine set_k_local_qkgroup
  integer function get_my_qkgroup(myrank_world,nproc_world)
    implicit none
    integer,intent(in) :: myrank_world,nproc_world
    integer :: iqkg,red,iproc_local,iproc_world
    iproc_world=0
    do iqkg=1,nqkgroup
       if (iqkg <= mod(nproc_world,nqkgroup)) then
          red = 1
       else
          red = 0
       endif
       do iproc_local=1,nproc_world/nqkgroup+red
          if (myrank_world == iproc_world) then
             get_my_qkgroup = iqkg
             return
          endif
          iproc_world = iproc_world + 1
       enddo
    enddo
  end function get_my_qkgroup
  ! subroutine set_q_local_qkgroup(Nq)
  !   implicit none
  !   integer,intent(in) :: Nq ! total number of q-points in irr.BZ
  !   integer,allocatable :: nq_local_all(:),iq_index_all(:,:)
  !   allocate(nq_local_all(nqkgroup),       iq_index_all(nqkgroup,Nq/nqkgroup+1))
  !   call set_index_rsmpi(Nq,nqkgroup,nq_local_all,iq_index_all)
  !   nq_local_qkgroup = nq_local_all(iqkgroup)
  !   allocate(iq_index_qkgroup(nq_local_qkgroup))
  !   iq_index_qkgroup(1:nq_local_qkgroup) =   iq_index_all(iqkgroup,1:nq_local_qkgroup)
  !   deallocate(nq_local_all,iq_index_all)
  ! end subroutine set_q_local_qkgroup
end module RSMPI_qkgroup
