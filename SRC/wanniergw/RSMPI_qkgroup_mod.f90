!     MPI version of FPLMTO-GW code
!     R. Sakuma 2007

!     This module contains variables needed
!        for q- and k-point (double) parallelization,
!     and is used for the calculations like
!        F(q) = \sigma_{k}^{full B.Z.} f(q,k)
!     Here q lies in the irreducible B.Z.

!     All processes are divided into several subgroups,
!     and each q-point (which lies in irreducible B.Z.) is
!     treated by one subgroup.
!     For each q-point, the summation with respect to k-points (which
!     lies in FULL B.Z.) is performed by processes in the subgroup.

!     For this purpose, a local communicator "comm_qkgroup" is created
!                                            by using "MPI_Comm_Split".

module rsmpi_qkgroup
  use rsmpi
  implicit none
  ! S:

  ! MPI parameters
  integer :: nqkgroup ! number of subgroups
  ! it is equivalent to the number of q-points
  !                calculated at the same time

  integer :: iqkgroup ! the subgroup the current process belongs


  integer :: nproc_qkgroup    ! number of processes in the subgroup
  integer :: myrank_qkgroup   ! rank of the calling process in the subgroup


  logical :: file_io_qkgroup  ! The processes with this value .true.
  ! dump output file for each q (ex.W(q))
  ! (true only myrank_qkgroup==0)


  integer :: nq_local_qkgroup ! number of q-points treated in the subgroup
  integer,allocatable :: iq_index_qkgroup(:) ! index

  integer :: nk_local_qkgroup ! number of k-points treated in the process
  integer,allocatable :: ik_index_qkgroup(:)

  integer :: ierror_qkgroup   ! error check

  integer :: comm_qkgroup ! new communicator

  ! ID of the subgroup (ex. 0000012, 0123456)
  ! used mainly for output filename (ex. WVR.RSMPI0000012 )
  ! if total number of subgroups >= 10^8, increase the size of the array
  character(7) :: qkgroup_id   !
  integer :: ifile_qkgroup    ! file id

contains
  !------------------------------------------------------
  subroutine RSMPI_qkgroup_Init(Nq,Nk)
    ! INPUT
    !   Nq: total number of q-points
    !   Nk: total number of k-points

    ! VARIABLES TO BE SET
    !   nqgroup
    !   iqgroup
    !   file_io_qkgroup
    !   nq_local_qkgroup
    !   iq_index_qkgroup
    !   nk_local_qkgroup
    !   ik_index_qkgroup
    !   comm_qkgroup
    !   qkgroup_id

    implicit none
    integer,intent(in) :: Nq,Nk


    !   nproc_rsmpi: total number of processes (in MPI_COMM_WORLD)
    if (Nq <= nproc_rsmpi) then
       nqkgroup = Nq
    else
       nqkgroup = nproc_rsmpi
    endif

    if (Is_IO_Root_RSMPI()) then
       write(6,*) "RS:   --- RSMPI_qkgroup_Init ---"
       write(6,*) "RS: Number of qk groups = ",nqkgroup
    endif

    iqkgroup = get_my_qkgroup(myrank_rsmpi,nproc_rsmpi)

    ! create new (local) communicator..
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,iqkgroup, &
         myrank_rsmpi,comm_qkgroup,ierror_qkgroup)
    call MPI_COMM_SIZE(comm_qkgroup,nproc_qkgroup,ierror_qkgroup)
    call MPI_COMM_RANK(comm_qkgroup,myrank_qkgroup,ierror_qkgroup)

    ! For each subgroup, file I/O is performed by the process with rank=0
    if (myrank_qkgroup == 0) then
       file_io_qkgroup = .true.
    else
       file_io_qkgroup = .false.
    endif

    call set_q_local_qkgroup(Nq)
    call set_k_local_qkgroup(Nk)

    call set_group_id_qkgroup()


    ! print group
    call print_info_qkgroup(Nq,Nk)

    if (Is_IO_Root_RSMPI()) then
       write(6,*) "RS:   --- end of RSMPI_qkgroup_Init ---"
    endif

  end subroutine RSMPI_qkgroup_Init

  !--------------------------------------------------------
  integer function get_my_qkgroup(myrank_world,nproc_world)
    implicit none
    integer,intent(in) :: myrank_world,nproc_world
    ! local
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
  !--------------------------------------------------------
  subroutine set_q_local_qkgroup(Nq)
    implicit none
    integer,intent(in) :: Nq ! total number of q-points in irr.BZ

    integer,allocatable :: nq_local_all(:),iq_index_all(:,:)
    allocate(nq_local_all(nqkgroup), &
         iq_index_all(nqkgroup,Nq/nqkgroup+1))

    call set_index_rsmpi(Nq,nqkgroup,nq_local_all,iq_index_all)

    nq_local_qkgroup = nq_local_all(iqkgroup)
    allocate(iq_index_qkgroup(nq_local_qkgroup))

    iq_index_qkgroup(1:nq_local_qkgroup) = &
         iq_index_all(iqkgroup,1:nq_local_qkgroup)

    deallocate(nq_local_all,iq_index_all)
  end subroutine set_q_local_qkgroup
  !--------------------------------------------------------
  ! same as the above subroutine..
  ! just q is replaced to k
  subroutine set_k_local_qkgroup(Nk)
    implicit none
    integer,intent(in) :: Nk ! total number of k-points in full BZ

    integer,allocatable :: nk_local_all(:),ik_index_all(:,:)
    allocate(nk_local_all(nproc_qkgroup), &
         ik_index_all(nproc_qkgroup,Nk/nproc_qkgroup+1))

    call set_index_rsmpi(Nk,nproc_qkgroup,nk_local_all,ik_index_all)


    nk_local_qkgroup = nk_local_all(myrank_qkgroup+1) ! rank is 0,1,..,N-1
    allocate(ik_index_qkgroup(nk_local_qkgroup))

    ik_index_qkgroup(1:nk_local_qkgroup) = &
         ik_index_all(myrank_qkgroup+1,1:nk_local_qkgroup)

    deallocate(nk_local_all,ik_index_all)
  end subroutine set_k_local_qkgroup
  !--------------------------------------------------------
  subroutine set_group_id_qkgroup()
    ! set group id...
    implicit none

    if (iqkgroup/1000000 >= 10) then
       ! RS: I don't expect this will happen..
       write(6,*) &
            "RS:number of subgroups exceeds 10^6: Modify RSMPI_qkgroup_mod.F!"
       call MPI_ABORT(MPI_COMM_WORLD,99,ierror_qkgroup)
    endif
    write(qkgroup_id,'(i7.7)') iqkgroup

  end subroutine set_group_id_qkgroup
  !--------------------------------------------------------
  subroutine print_info_qkgroup(Nq,Nk)
    implicit none
    integer,intent(in) :: Nq,Nk
    ! First, print how processes are divided
    if (Is_IO_Root_RSMPI()) then
       write(6,*) "RS: ", nproc_rsmpi,"processes are divided into", &
            nqkgroup, " subgroups"
       write(6,*) "# rank_global group_id nproc_local rank_local"
    endif
    write(buf_rsmpi,*) myrank_id_rsmpi,"      ",qkgroup_id, &
         nproc_qkgroup,myrank_qkgroup
    call RSMPI_Write(6)

    ! q-points
    if (Is_IO_Root_RSMPI()) then
       write(6,*)
       write(6,*) "RS: ", Nq," q-points are divided as follows"
       write(6,*) "# iqkgroup nq_local  iq_index"
    endif
    buf_rsmpi = ""
    if (file_io_qkgroup) then
       write(buf_rsmpi,*) "> ",qkgroup_id,nq_local_qkgroup, &
            " : ", iq_index_qkgroup(1), " ... ,", &
            iq_index_qkgroup(nq_local_qkgroup)
    endif
    call RSMPI_Write(6)

    ! k-points
    if (Is_IO_Root_RSMPI()) then
       write(6,*)
       write(6,*) "RS: ", Nk," k-points are divided as follows"
       write(6,*) "# iqkgroup rank_local, nk_local, ikindex"
    endif
    if (nk_local_qkgroup > 0) then
       write(buf_rsmpi,*) "> ",qkgroup_id,myrank_qkgroup, nk_local_qkgroup, &
            " :",ik_index_qkgroup(1), " ... ,", &
            ik_index_qkgroup(nk_local_qkgroup)
    else
       write(buf_rsmpi,*) "> ",qkgroup_id,myrank_qkgroup, nk_local_qkgroup, &
            " :"
    endif
    call RSMPI_Write(6)

  end subroutine print_info_qkgroup
  !--------------------------------------------------------
end module RSMPI_qkgroup
