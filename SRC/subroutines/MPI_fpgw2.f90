module m_mpi !MPI utility for fpgw
  implicit none
  include "mpif.h"
  integer :: mpi__size
  integer :: mpi__rank
  logical :: mpi__root
  integer :: comm
  integer :: mpi__sizeMG 
  integer :: mpi__rankMG
!MPI for hrcxq
  integer :: comm_q, mpi__rank_q, mpi__size_q
  integer :: comm_k, mpi__rank_k, mpi__size_k
  integer :: comm_b, mpi__rank_b, mpi__size_b
  integer :: comm_root_k
  logical :: mpi__root_q, mpi__root_k
  integer, allocatable :: mpi__npr_col(:), mpi__ipr_col(:)
!MPI for Sc (hsfp0_sc --job=2)
  integer :: comm_w, mpi__rank_w, mpi__size_w
  logical :: mpi__root_w
  integer :: worker_intask = 1 !default used in ixc /= 2

  integer,private :: mpi__info
  integer,private:: ista(MPI_STATUS_SIZE )
contains
  subroutine MPI__Initialize(commin)
    implicit none
    character(1024*4) :: cwd, stdout
    integer,optional:: commin
    comm=MPI_COMM_WORLD
    if(present(commin)) comm= commin 
    !merge(commin,MPI_COMM_WORLD,present(commin))
    call getcwd(cwd)           ! get current working directory
    call MPI_Init( mpi__info ) ! current working directory is changed if mpirun is not used
    call MPI_Comm_rank( comm, mpi__rank, mpi__info )
    call MPI_Comm_size( comm, mpi__size, mpi__info )
    mpi__root= mpi__rank==0
    if( mpi__root ) call chdir(cwd)        ! recover current working directory
  end subroutine MPI__Initialize

!  MPI__SplitXq is only used in hrcxq for q-points, k-points, and MPB parallel.
! example in case of n_bpara = 2 and n_kpara  = 3
! mpi__rank                           : 0,1,2,3,4,5, 6,7,8,9,10,11
! color = mpi__rank/(n_bpara*n_kpara) : 0,0,0,0,0,0, 1,1,1,1, 1, 1,
! mpi__rank_q                         : 0,1,2,3,4,5, 0,1,2,3, 4, 5
! mpi__root_q                         : T,F,F,F,F,F, T,F,F,F, F, F
! color = mpi__rank_q/n_bpara         : 0,0,1,1,2,2  0,0,1,1, 2, 2
! color = mod(mpi__rank_q,n_bpara)    : 0,1,0,1,0,1  0,1,0,1, 0, 1

  subroutine MPI__SplitXq(n_bpara, n_kpara)
    implicit none
    integer, intent(in) :: n_bpara, n_kpara
    integer :: color

    color = mpi__rank/(n_bpara*n_kpara)
    call mpi_comm_split(comm, color, mpi__rank, comm_q, mpi__info)
    call mpi_comm_rank(comm_q, mpi__rank_q, mpi__info)
    call mpi_comm_size(comm_q, mpi__size_q, mpi__info)
    mpi__root_q = mpi__rank_q == 0

    color = mpi__rank_q/n_bpara
    call mpi_comm_split(comm_q, color, mpi__rank, comm_b, mpi__info)
    call mpi_comm_rank(comm_b, mpi__rank_b, mpi__info)
    call mpi_comm_size(comm_b, mpi__size_b, mpi__info)

    color = mod(mpi__rank_q,n_bpara)
    call mpi_comm_split(comm_q, color, mpi__rank, comm_k, mpi__info)
    call mpi_comm_rank(comm_k, mpi__rank_k, mpi__info)
    call mpi_comm_size(comm_k, mpi__size_k, mpi__info)

    color = merge(0, MPI_UNDEFINED, mpi__rank_k == 0)
    mpi__root_k = mpi__rank_k == 0
    call mpi_comm_split(comm_q, color, mpi__rank, comm_root_k, mpi__info)

    write(06,'(X,A,4I5,2L2)') "MPI: rank, rank_q, rank_k, rank_b, root_q, root_k ", &
                mpi__rank, mpi__rank_q, mpi__rank_k, mpi__rank_b, mpi__root_q, mpi__root_k

  end subroutine MPI__SplitXq

  subroutine MPI__Setnpr_col(npr, npr_col)
    integer, intent(in) :: npr
    integer, intent(out) :: npr_col
    integer :: irank_b, ipr_col
    if(.not.allocated(mpi__npr_col)) allocate(mpi__npr_col(0:mpi__size_b-1))
    if(.not.allocated(mpi__ipr_col)) allocate(mpi__ipr_col(0:mpi__size_b-1))
    ipr_col = 1
    do irank_b = 0, mpi__size_b - 1
      npr_col = (npr + irank_b)/mpi__size_b
      mpi__npr_col(irank_b) = npr_col
      mpi__ipr_col(irank_b) = ipr_col
      ipr_col = ipr_col + npr_col
      if (npr_col == 0) call rx("MPI__Setnpr_col: use small parallelization")
    enddo
    npr_col = mpi__npr_col(mpi__rank_b)
    ipr_col = mpi__ipr_col(mpi__rank_b)
    write(06,'(X,A,I5,2I7)') "mpi__rank_b, ipr_col, npr_col=", mpi__rank_b, ipr_col, npr_col
  end subroutine MPI__Setnpr_col

  subroutine MPI__GatherXqw(xqw, xqw_all, npr, npr_col)
    integer, intent(in) :: npr, npr_col
    complex(8), intent(in) :: xqw(npr,npr_col)
    complex(8), intent(out) :: xqw_all(npr,npr)  ! we suppose only column was split
    integer, allocatable :: data_disp(:), data_size(:)
    integer :: irank_b, rank
    xqw_all(1:npr,1:npr) = (0D0, 0D0)
    if(mpi__size_b == 1 .and. npr == npr_col) then
      xqw_all = xqw
      return
    endif
    allocate(data_size(0:mpi__size_b-1), data_disp(0:mpi__size_b-1))
    do irank_b = 0, mpi__size_b -1
      data_size(irank_b) = npr*mpi__npr_col(irank_b)
      data_disp(irank_b) = npr*(mpi__ipr_col(irank_b)-1)
    enddo
    call mpi_allgatherv(xqw, npr*npr_col, mpi_complex16, xqw_all, data_size, data_disp, &
                  &  mpi_complex16, comm_root_k, mpi__info)
    ! call mpi_gatherv(xqw, npr*npr_col, mpi_complex16, xqw_all, data_size, data_disp, &
    !               &  mpi_complex16, 0, comm_root_k, mpi__info)
    deallocate(data_size, data_disp)
  end subroutine MPI__GatherXqw
  subroutine MPI__GatherXqw_kind4(xqw, xqw_all, npr, npr_col)
    integer, intent(in) :: npr, npr_col
    complex(4), intent(in) :: xqw(npr,npr_col)
    complex(4), intent(out) :: xqw_all(npr,npr)  ! we suppose only column was split
    integer, allocatable :: data_disp(:), data_size(:)
    integer :: irank_b, rank
    xqw_all(1:npr,1:npr) = (0.0, 0.0)
    if(mpi__size_b == 1 .and. npr == npr_col) then
      xqw_all = xqw
      return
    endif
    allocate(data_size(0:mpi__size_b-1), data_disp(0:mpi__size_b-1))
    do irank_b = 0, mpi__size_b -1
      data_size(irank_b) = npr*mpi__npr_col(irank_b)
      data_disp(irank_b) = npr*(mpi__ipr_col(irank_b)-1)
    enddo
    call mpi_allgatherv(xqw, npr*npr_col, mpi_complex, xqw_all, data_size, data_disp, &
                  &  mpi_complex, comm_root_k, mpi__info)
    ! call mpi_gatherv(xqw, npr*npr_col, mpi_complex16, xqw_all, data_size, data_disp, &
    !               &  mpi_complex16, 0, comm_root_k, mpi__info)
    deallocate(data_size, data_disp)
  end subroutine MPI__GatherXqw_kind4
  integer function get_mpi_size(communicator) result(mpi_size)
    implicit none
    integer, intent(in), optional :: communicator
    integer :: comm_in, ierr
    comm_in = comm
    if(present(communicator)) comm_in = communicator
    call MPI_Comm_size(comm_in, mpi_size, ierr)
  end function get_mpi_size
  subroutine MPI__SplitSc(n_wpara)
    implicit none
    integer, intent(in) :: n_wpara
    integer :: color
    if(n_wpara < 1) call rx("MPI__SplitSc: n_wpara < 1")
    if(n_wpara > mpi__size) call rx("MPI__SplitSc: n_wpara > mpi__size")
    color = mpi__rank/n_wpara
    call mpi_comm_split(comm, color, mpi__rank, comm_w, mpi__info)
    call mpi_comm_rank(comm_w, mpi__rank_w, mpi__info)
    call mpi_comm_size(comm_w, mpi__size_w, mpi__info)
    mpi__root_w = mpi__rank_w == 0
  end subroutine MPI__SplitSc
  subroutine MPI__Initialize_magnon(commin)
    implicit none
    character(1024*4) :: cwd, stdout
    integer,optional:: commin
    comm=MPI_COMM_WORLD
    if(present(commin)) comm= commin
    !comm= merge(commin,MPI_COMM_WORLD,present(commin))
    call getcwd(cwd)          ! get current working directory
    call MPI_Init( mpi__info ) ! current working directory is changed if mpirun is not used
    call MPI_Comm_rank( comm, mpi__rankMG, mpi__info )
    call MPI_Comm_size( comm, mpi__sizeMG, mpi__info )
    mpi__root=  mpi__rankMG == 0 
    if( mpi__root ) call chdir(cwd)        ! recover current working directory
  end subroutine MPI__Initialize_magnon
  subroutine MPI__consoleout(idn)
    use m_lgunit,only:stdo,stdl
    implicit none
    character(1024*4) :: cwd, stdout
    character*(*):: idn
    if( mpi__size == 1 ) return
    if( mpi__root ) then
      write(6,"(' MPI outputs in each rank are in stdout.{RankId}.',a)")idn
      call flush(stdo)
    end if
    write(stdout,"('stdout.',i4.4,'.',a)") mpi__rank,idn
    open(unit=6,file=trim(stdout))
    write(6,"(a,i3)")" ### console output for rank=",mpi__rank
  end subroutine MPI__consoleout
  subroutine MPI__consoleout_magnon(idn,size_lim)
    use m_lgunit,only:stdo,stdl
    implicit none
    character(1024*4) :: cwd, stdout
    character*(*):: idn
    integer , intent(in) :: size_lim
    if(mpi__sizeMG > size_lim) mpi__sizeMG=size_lim !Reduce size for magnon (avoid memory leak)      
    if( mpi__sizeMG == 1 ) return
    if( mpi__root ) then
      write(6,"(' MPI outputs in each rank are in stdout.{RankId}.',a)")idn
      call flush(stdo)
    end if
    write(stdout,"('stdout.',i4.4,'.',a)") mpi__rankMG,idn
    open(unit=6,file=trim(stdout))
    write(6,"(a,i3)")" ### console output for rank=",mpi__rankMG
  end subroutine MPI__consoleout_magnon
  subroutine MPI__Broadcast( data )
    implicit none
    integer, intent(inout) :: data
    call MPI_Bcast( data, 1, MPI_INTEGER, 0, comm, mpi__info )
  end subroutine MPI__Broadcast
  subroutine MPI__REAL8send(data,n,dest)
    implicit none
    real(8):: data(n)
    integer :: n,dest,ierr
    call MPI_Send(data,n,MPI_REAL8,dest,mpi__rank, comm,ierr)
  end subroutine MPI__REAL8send
  subroutine MPI__REAL8recv(data,n,src)
    implicit none
    real(8):: data(n)
    integer :: n,src,ierr
    call MPI_Recv(data,n,MPI_REAL8,src,src, comm,ista,ierr)
  end subroutine MPI__REAL8recv
  subroutine MPI__DbleCOMPLEXsend(data,n,dest)
    implicit none
    complex(8):: data(n)
    integer :: n,dest,ierr
    call MPI_Send(data,n,MPI_COMPLEX16,dest,mpi__rank, comm,ierr)
  end subroutine MPI__DbleCOMPLEXsend
  subroutine MPI__DbleCOMPLEXrecv(data,n,src)
    implicit none
    complex(8):: data(n)
    integer :: n,src,ierr
    call MPI_Recv(data,n,MPI_COMPLEX16,src,src, comm,ista,ierr)
  end subroutine MPI__DbleCOMPLEXrecv
  subroutine MPI__DbleCOMPLEXsendQ(data,n,destQ)
    implicit none
    complex(8):: data(n)
    integer :: n,destQ,ierr
    call MPI_Send(data,n,MPI_COMPLEX16,destQ,0,comm,ierr)
  end subroutine MPI__DbleCOMPLEXsendQ
  subroutine MPI__DbleCOMPLEXrecvQ(data,n,srcQ)
    implicit none
    complex(8):: data(n)
    integer :: n,srcQ,ierr
    call MPI_Recv(data,n,MPI_COMPLEX16,srcQ,0,comm,ista,ierr)
  end subroutine MPI__DbleCOMPLEXrecvQ
  subroutine MPI__AllreduceSum( data, sizex, communicator)
    implicit none
    integer, intent(in) :: sizex
    complex(8), intent(inout) :: data(sizex)
    complex(8), allocatable   :: mpi__data(:) 
    integer, intent(in), optional :: communicator
    integer :: comm_in, mpi_size_comm_in, ierr
    if(mpi__size == 1) return
    comm_in = comm
    if(present(communicator)) comm_in = communicator
    mpi_size_comm_in = get_mpi_size(comm_in)
    if(mpi_size_comm_in == 1) return
    allocate(mpi__data(sizex))
    mpi__data = data
    call MPI_Allreduce( mpi__data, data, sizex, MPI_DOUBLE_COMPLEX, MPI_SUM, comm_in, mpi__info )
    deallocate( mpi__data )
  end subroutine MPI__AllreduceSum
  subroutine MPI__reduceSum( root, data, sizex, communicator)
    implicit none
    integer, intent(in) :: sizex,root
    complex(8), intent(inout) :: data(sizex)
    complex(8), allocatable   :: mpi__data(:) 
    integer, intent(in), optional :: communicator
    integer :: comm_in, mpi_size_comm_in, ierr
    if(mpi__size == 1) return
    comm_in = comm
    if(present(communicator)) comm_in = communicator
    mpi_size_comm_in = get_mpi_size(comm_in)
    if(mpi_size_comm_in == 1) return
    allocate(mpi__data(sizex))
    mpi__data = data
    call MPI_reduce( mpi__data, data, sizex, MPI_DOUBLE_COMPLEX, MPI_SUM, root, comm_in, mpi__info )
    deallocate( mpi__data )
    return
  end subroutine MPI__reduceSum
  subroutine MPI__reduceSum_kind4( root, data, sizex, communicator)
    implicit none
    integer, intent(in) :: sizex,root
    complex(4), intent(inout) :: data(sizex)
    complex(4), allocatable   :: mpi__data(:) 
    integer, intent(in), optional :: communicator
    integer :: comm_in, mpi_size_comm_in, ierr
    if( mpi__size == 1 ) return
    comm_in = comm
    if(present(communicator)) comm_in = communicator
    call MPI_Comm_size(comm_in, mpi_size_comm_in, ierr)
    if( mpi_size_comm_in == 1 ) return
    allocate(mpi__data(sizex))
    mpi__data = data
    call MPI_reduce( mpi__data, data, sizex, MPI_COMPLEX, MPI_SUM, root, comm_in, mpi__info )
    deallocate( mpi__data )
    return
  end subroutine MPI__reduceSum_kind4
  subroutine MPI__AllreduceMax( data, sizex )
    implicit none
    integer, intent(in) :: sizex
    integer, intent(inout) :: data(sizex)
    integer, allocatable   :: mpi__data(:) 
    if( mpi__size == 1 ) return
    allocate(mpi__data(sizex))
    mpi__data = data
    call MPI_Allreduce( mpi__data, data, sizex, MPI_INTEGER, MPI_MAX, comm, mpi__info )
    deallocate( mpi__data )
  end subroutine MPI__AllreduceMax
!MO Following subroutines are for CPU(host) and GPU(device) implementations 2024/12/27
  subroutine MPI__zBcast_h(data, sizex, communicator, sender)
    implicit none
    integer, intent(in) :: sizex
    complex(8), intent(inout) :: data(sizex)
    integer, intent(in), optional :: communicator, sender
    integer :: comm_in, sender_in, mpi_size_comm_in
    comm_in = comm
    sender_in = 0
    if(present(communicator)) comm_in = communicator
    if(present(sender)) sender_in = sender
    mpi_size_comm_in = get_mpi_size(comm_in)
    if(mpi_size_comm_in == 1) return
    call MPI_Bcast(data, sizex, MPI_DOUBLE_COMPLEX, sender_in, comm_in, mpi__info)
  end subroutine MPI__zBcast_h
#ifdef __GPU
  subroutine MPI__zBcast_d(data_d, sizex, communicator, sender)
    use cudafor
    implicit none
    integer, intent(in) :: sizex
    complex(8), intent(inout), device :: data_d(sizex)
    complex(8) :: data_h(sizex) !Host data for MPI communication
    integer, intent(in), optional :: communicator, sender
    integer :: comm_in, sender_in, mpi_size_comm_in
    comm_in = comm
    sender_in = 0
    if(present(communicator)) comm_in = communicator
    if(present(sender)) sender_in = sender
    mpi_size_comm_in = get_mpi_size(comm_in)
    if(mpi_size_comm_in == 1) return
    data_h(:) = data_d(:) ! copy to host
    call MPI_Bcast(data_h, sizex, MPI_DOUBLE_COMPLEX, sender_in, comm_in, mpi__info)
    data_d(:) = data_h(:) !copy to device
  end subroutine MPI__zBcast_d
#endif
end module m_mpi

subroutine MPI__sxcf_rankdivider(irkip_all,nspinmx,nqibz,ngrp,nq,irkip)
  use m_mpi,only: mpi__rank,mpi__size
  use m_mpi, only: worker_intask !set as 1 in the case of without omega parallelization
  implicit none
  integer, intent(out) :: irkip    (nspinmx,nqibz,ngrp,nq)
  integer, intent(in)  :: irkip_all(nspinmx,nqibz,ngrp,nq)
  integer, intent(in)  :: nspinmx,nqibz,ngrp,nq
  integer :: ispinmx,iqibz,igrp,iq
  integer :: total
  integer, allocatable :: vtotal(:)
  integer :: indexi, indexe
  integer :: p, ngroup
  if( mpi__size == 1 ) then
     irkip = irkip_all
     return
  end if
  total = count(irkip_all>0)
  ngroup = mpi__size/worker_intask
  write(6,"('MPI__sxcf_rankdivider:$')")
  write(6,"('nspinmx,nqibz,ngrp,nq,total=',5i6)") nspinmx,nqibz,ngrp,nq,total
  write(6,'(A,2I5)') 'MPI: Worker in Task, # of groups', worker_intask, ngroup
  allocate( vtotal(0:mpi__size-1) )
  ! vtotal(:) = total/mpi__size
  vtotal(:) = total/ngroup
  do p=1, mod(total, ngroup)
     vtotal(p-1) = vtotal(p-1) + 1
  end do
  indexe=0
  indexi=-999999
  do p=0, mpi__rank/worker_intask !same definition with color in SplitSc
     indexi = indexe+1
     indexe = indexi+vtotal(p)-1
  end do
  deallocate(vtotal)
  total = 0
  irkip(:,:,:,:) = 0
  do iq=1, nq
     do ispinmx=1, nspinmx
        do iqibz=1, nqibz
           do igrp=1, ngrp
              if( irkip_all(ispinmx,iqibz,igrp,iq) >0 ) then
                 total = total + 1
                 if( indexi<=total .and. total<=indexe ) then
                    irkip(ispinmx,iqibz,igrp,iq) = irkip_all(ispinmx,iqibz,igrp,iq)
                 endif
              endif
           enddo
        enddo
     enddo
  enddo
end subroutine MPI__sxcf_rankdivider
