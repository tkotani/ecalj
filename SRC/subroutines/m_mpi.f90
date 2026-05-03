module m_mpi !MPI utility (unified from m_mpi + m_MPItk)
  use mpi
  use m_lgunit, only: stdo, stdl
  implicit none
  integer :: mpi__size
  integer :: mpi__rank
  logical :: mpi__root
  integer :: comm
!-- m_MPItk compatible variables
  integer, protected :: procid, master = 0, nsize
  logical, protected :: master_mpi, readtk = .false.
  character(8), protected :: strprocid
!MPI for hrcxq
  integer :: comm_q, mpi__rank_q, mpi__size_q
  integer :: comm_k, mpi__rank_k, mpi__size_k
  integer :: comm_b, mpi__rank_b, mpi__size_b
  integer :: comm_root_k, mpi__rank_root_k, mpi__size_root_k
  logical :: mpi__root_q, mpi__root_k
  integer, allocatable :: mpi__npr_col(:), mpi__ipr_col(:)
!MPI for Sc (hsfp0_sc --job=2)
  integer :: comm_w, mpi__rank_w, mpi__size_w
  logical :: mpi__root_w,ipr=.true.
  integer :: worker_intask = 1 !default used in ixc /= 2
!Simple split of MPI communicator
  integer :: comm_s, mpi__rank_s, mpi__size_s, comm_s_idx
  logical :: mpi__root_s

  integer,private :: mpi__info
  integer,private:: ista(MPI_STATUS_SIZE )
contains
  subroutine setipr(comm)
    integer:: comm
    logical,external:: cmdopt0
    call MPI_Comm_rank( comm, mpi__rank, mpi__info )
    mpi__root= mpi__rank==0
    ipr=mpi__root
    if(cmdopt0('--fullstdo')) ipr=.true.
  end subroutine setipr
  subroutine MPI__Initialize(commin)
    implicit none
    character(1024*4) :: cwd, stdout
    character(10):: i2char
    integer,optional:: commin
    logical,external:: cmdopt0
    logical :: initialized
    comm=MPI_COMM_WORLD
    if(present(commin)) comm= commin 
    !merge(commin,MPI_COMM_WORLD,present(commin))
    call getcwd(cwd)           ! get current working directory
    call MPI_Initialized(initialized, mpi__info)
    if(.not. initialized) call MPI_Init( mpi__info ) ! current working directory is changed if mpirun is not used
    call MPI_Comm_rank( comm, mpi__rank, mpi__info )
    call MPI_Comm_size( comm, mpi__size, mpi__info )
    mpi__root= mpi__rank==0
    if( mpi__root ) call chdir(cwd)        ! recover current working directory
    ipr=mpi__root
    if(cmdopt0('--fullstdo')) ipr=.true.
    !-- m_MPItk compatible
    procid = mpi__rank
    nsize = mpi__size
    master_mpi = mpi__root
    readtk = .true.
    strprocid = trim(i2char(procid))
    call m_setargs_init()  ! ensures arglist / sname are populated for any Fortran binary
  end subroutine MPI__Initialize

  subroutine m_setargs_init()
    use m_args, only: m_setargs
    use m_ext,  only: m_ext_init
    call m_setargs()
    call m_ext_init()
  end subroutine m_setargs_init

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
    if (comm_root_k /= MPI_COMM_NULL) then
      call mpi_comm_rank(comm_root_k, mpi__rank_root_k, mpi__info)
      call mpi_comm_size(comm_root_k, mpi__size_root_k, mpi__info)
    endif
    if(ipr)write(06,'(X,A,4I5,2L2)') "MPI: rank, rank_q, rank_k, rank_b, root_q, root_k ", &
                mpi__rank, mpi__rank_q, mpi__rank_k, mpi__rank_b, mpi__root_q, mpi__root_k

  end subroutine MPI__SplitXq
  
  subroutine MPI__Split(n_split)
    implicit none
    integer, intent(in) :: n_split
    integer :: color
    color = mod(mpi__rank, n_split)
    comm_s_idx = color
    call mpi_comm_split(comm, color, mpi__rank, comm_s, mpi__info)
    call mpi_comm_rank(comm_s, mpi__rank_s, mpi__info)
    call mpi_comm_size(comm_s, mpi__size_s, mpi__info)
    mpi__root_s = mpi__rank_s == 0
    write(06,'(X,A,2I5,L2,I5)') "MPI: rank, rank_s, root_s, comm_s_idx ", mpi__rank, mpi__rank_s, mpi__root_s, comm_s_idx
  end subroutine MPI__Split

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
    if(ipr)write(06,'(X,A,I5,2I7)') "mpi__rank_b, ipr_col, npr_col=", mpi__rank_b, ipr_col, npr_col
  end subroutine MPI__Setnpr_col

  subroutine MPI__GatherXqw(xqw, xqw_all, npr, npr_col, collector_rank)
    integer, intent(in) :: npr, npr_col
    integer, intent(in), optional :: collector_rank
    complex(8), intent(in) :: xqw(npr,npr_col)
    complex(8), intent(inout) :: xqw_all(npr,npr)  ! we suppose only column was split
    integer, allocatable :: data_disp(:), data_size(:)
    integer :: irank_b, collector_rank_in
    if(mpi__size_b == 1 .and. npr == npr_col) then
      xqw_all(:,:) = xqw(:,:)
      return
    endif
    allocate(data_size(0:mpi__size_b-1), data_disp(0:mpi__size_b-1))
    do irank_b = 0, mpi__size_b -1
      data_size(irank_b) = npr*mpi__npr_col(irank_b)
      data_disp(irank_b) = npr*(mpi__ipr_col(irank_b)-1)
    enddo
    ! call mpi_allgatherv(xqw, npr*npr_col, mpi_complex16, xqw_all, data_size, data_disp, &
    !               &  mpi_complex16, comm_root_k, mpi__info)
    collector_rank_in = 0
    if(present(collector_rank)) collector_rank_in = collector_rank
    call mpi_gatherv(xqw, npr*npr_col, mpi_complex16, xqw_all, data_size, data_disp, &
                  &  mpi_complex16, collector_rank_in, comm_root_k, mpi__info)
    deallocate(data_size, data_disp)
  end subroutine MPI__GatherXqw
  subroutine MPI__GatherXqw_c(xqw, xqw_all, npr, npr_col, collector_rank)
    integer, intent(in) :: npr, npr_col
    integer, intent(in), optional :: collector_rank
    complex(4), intent(in) :: xqw(npr,npr_col)
    complex(4), intent(out) :: xqw_all(npr,npr)  ! we suppose only column was split
    integer, allocatable :: data_disp(:), data_size(:)
    integer :: irank_b, collector_rank_in
    if(mpi__size_b == 1 .and. npr == npr_col) then
      xqw_all(:,:) = xqw(:,:)
      return
    endif
    allocate(data_size(0:mpi__size_b-1), data_disp(0:mpi__size_b-1))
    do irank_b = 0, mpi__size_b -1
      data_size(irank_b) = npr*mpi__npr_col(irank_b)
      data_disp(irank_b) = npr*(mpi__ipr_col(irank_b)-1)
    enddo
    ! call mpi_allgatherv(xqw, npr*npr_col, mpi_complex, xqw_all, data_size, data_disp, &
    !               &  mpi_complex, comm_root_k, mpi__info)
    collector_rank_in = 0
    if(present(collector_rank)) collector_rank_in = collector_rank
    call mpi_gatherv(xqw, npr*npr_col, mpi_complex, xqw_all, data_size, data_disp, &
                  &  mpi_complex, collector_rank_in, comm_root_k, mpi__info)
    deallocate(data_size, data_disp)
  end subroutine MPI__GatherXqw_c
  integer function get_mpi_size(communicator) result(mpi_size)
    implicit none
    integer, intent(in), optional :: communicator
    integer :: comm_in, ierr
    comm_in = comm
    if(present(communicator)) comm_in = communicator
    call MPI_Comm_size(comm_in, mpi_size, ierr)
  end function get_mpi_size
  logical function get_mpi_master(communicator) result(mpi_master)
    integer, intent(in), optional :: communicator
    integer :: comm_in, ierr, mpi_rank
    comm_in = comm
    if(present(communicator)) comm_in = communicator
    call mpi_comm_rank(comm_in, mpi_rank, ierr)
    mpi_master = (mpi_rank == 0)
  end function get_mpi_master
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
  subroutine MPI__reduceSum_c( root, data, sizex, communicator)
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
  end subroutine MPI__reduceSum_c
!  subroutine MPI__AllreduceMax( data, sizex ) !currently unused
!    implicit none
!    integer, intent(in) :: sizex
!    integer, intent(inout) :: data(sizex)
!    integer, allocatable   :: mpi__data(:)
!    if( mpi__size == 1 ) return
!    allocate(mpi__data(sizex))
!    mpi__data = data
!    call MPI_Allreduce( mpi__data, data, sizex, MPI_INTEGER, MPI_MAX, comm, mpi__info )
!    deallocate( mpi__data )
!  end subroutine MPI__AllreduceMax
!MO Addtional subroutines for MPI 2024/12/28
  subroutine MPI__AllreduceAND(data, communicator)
    implicit none
    logical, intent(inout) :: data
    integer, intent(in), optional :: communicator
    logical :: mpi__data
    integer :: comm_in
    comm_in = comm
    if(present(communicator)) comm_in = communicator
    if(get_mpi_size(comm_in) == 1) return
    mpi__data =  data
    call MPI_Allreduce(mpi__data, data, 1, MPI_LOGICAL, MPI_LAND, comm_in, mpi__info)
  end subroutine MPI__AllreduceAND
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

  subroutine MPI__AllreduceSumReal( data, sizex, communicator)
    implicit none
    integer, intent(in) :: sizex
    real(8), intent(inout) :: data(sizex)
    real(8), allocatable   :: mpi__data(:) 
    integer, intent(in), optional :: communicator
    integer :: comm_in
    if( mpi__size == 1 ) return
    allocate(mpi__data(sizex))
    mpi__data = data
    comm_in = comm
    if(present(communicator)) comm_in = communicator
    call MPI_Allreduce( mpi__data, data, sizex, MPI_DOUBLE_PRECISION, MPI_SUM, comm_in, mpi__info )
    deallocate( mpi__data )
  end subroutine MPI__AllreduceSumReal
  subroutine MPI__AllreduceSumRealSca( data, communicator)
    implicit none
    real(8), intent(inout) :: data
    real(8) :: mpi__data
    integer, intent(in), optional :: communicator
    integer :: comm_in
    if( mpi__size == 1 ) return
    mpi__data = data
    comm_in = comm
    if(present(communicator)) comm_in = communicator
    call MPI_Allreduce( mpi__data, data, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_in, mpi__info )
  end subroutine MPI__AllreduceSumRealSca

!#ifdef __GPU
!  subroutine MPI__zBcast_d(data_d, sizex, communicator, sender) !currently unused
!    use cudafor
!    implicit none
!    integer, intent(in) :: sizex
!    complex(8), intent(inout), device :: data_d(sizex)
!    complex(8) :: data_h(sizex) !Host data for MPI communication
!    integer, intent(in), optional :: communicator, sender
!    integer :: comm_in, sender_in, mpi_size_comm_in
!    comm_in = comm
!    sender_in = 0
!    if(present(communicator)) comm_in = communicator
!    if(present(sender)) sender_in = sender
!    mpi_size_comm_in = get_mpi_size(comm_in)
!    if(mpi_size_comm_in == 1) return
!    data_h(:) = data_d(:) ! copy to host
!    call MPI_Bcast(data_h, sizex, MPI_DOUBLE_COMPLEX, sender_in, comm_in, mpi__info)
!    data_d(:) = data_h(:) !copy to device
!  end subroutine MPI__zBcast_d
!#endif

  subroutine xmpbnd2(kpproc, ndham, ndat, eb)  !- Collect eb from various processors (MPI)
    implicit none
    integer:: kpproc(0:*), ndham, ndat
    double precision :: eb(ndham, ndat)
    integer :: i, ista, iend, ierr
    integer, dimension(:), allocatable :: offset, length
    real(8), allocatable :: buf_rv(:, :)
    allocate (offset(0:nsize), length(0:nsize))
    offset(0) = 0
    do i = 0, nsize - 1
      ista = kpproc(i)
      iend = kpproc(i + 1) - 1
      length(i) = (iend - ista + 1)*ndham
      offset(i + 1) = offset(i) + length(i)
    end do
    ista = kpproc(procid)
    iend = kpproc(procid + 1) - 1
    allocate (buf_rv(ndham, ndat))
    call mpi_allgatherv(eb(1:ndham, ista:iend), length(procid), mpi_double_precision, buf_rv, length, offset, mpi_double_precision,&
      comm, ierr)
    eb = buf_rv
    deallocate (buf_rv, offset, length)
  end subroutine xmpbnd2

end module m_mpi

subroutine MPI__sxcf_rankdivider(irkip_all,nspinmx,nqibz,ngrp,nq,irkip)
  use m_mpi,only: mpi__rank,mpi__size
  use m_mpi, only: ipr,worker_intask !set as 1 in the case of without omega parallelization
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
  if(ipr)write(6,"('MPI__sxcf_rankdivider:$')")
  if(ipr)write(6,"('nspinmx,nqibz,ngrp,nq,total=',5i6)") nspinmx,nqibz,ngrp,nq,total
  if(ipr)write(6,'(A,2I5)') 'MPI: Worker in Task, # of groups', worker_intask, ngroup
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
