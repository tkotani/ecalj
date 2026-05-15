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
!MPI communicator hierarchy: q-group layer (comm_q) + two persistent intra-group splits
  integer :: comm_q, mpi__rank_q, mpi__size_q
  logical, protected :: mpi__root_q
  integer, protected :: iq_qgroup=0, n_qgroup=1, worker_inQtask=2

!-- Xq split (k-priority): used by build_screened_coulomb and exchange
!   comm_k_xq: all worker ranks (k-parallel);  comm_b_xq: ω-parallel (size=n_bpara)
  integer, protected :: comm_k_xq, mpi__rank_k_xq, mpi__size_k_xq
  integer, protected :: comm_b_xq, mpi__rank_b_xq, mpi__size_b_xq
  logical, protected :: mpi__root_k_xq, mpi__root_b_xq
  integer, protected :: comm_root_k_xq, mpi__rank_root_k_xq, mpi__size_root_k_xq

!-- Sxc split (k-priority): used by sxcf_correlation
!   comm_k_sxc: all worker ranks (k-parallel);  comm_b_sxc: ω-parallel (size=n_bpara)
  integer, protected :: comm_k_sxc, mpi__rank_k_sxc, mpi__size_k_sxc
  integer, protected :: comm_b_sxc, mpi__rank_b_sxc, mpi__size_b_sxc
  logical, protected :: mpi__root_k_sxc, mpi__root_b_sxc

  integer, allocatable :: mpi__npr_col(:), mpi__ipr_col(:)
  logical :: ipr=.true.
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
    call getcwd(cwd)
    call MPI_Initialized(initialized, mpi__info)
    if(.not. initialized) call MPI_Init( mpi__info )
    call MPI_Comm_rank( comm, mpi__rank, mpi__info )
    call MPI_Comm_size( comm, mpi__size, mpi__info )
    mpi__root= mpi__rank==0
    if( mpi__root ) call chdir(cwd)
    ipr=mpi__root
    if(cmdopt0('--fullstdo')) ipr=.true.
    !-- m_MPItk compatible
    procid = mpi__rank
    nsize = mpi__size
    master_mpi = mpi__root
    readtk = .true.
    strprocid = trim(i2char(procid))
    call m_setargs_init()
  end subroutine MPI__Initialize

  subroutine m_setargs_init()
    use m_args, only: m_setargs
    use m_ext,  only: m_ext_init
    call m_setargs()
    call m_ext_init()
  end subroutine m_setargs_init

  subroutine MPI__InitQgroups(worker_in)
    !> Set up intra-group communicator comm_q and derived variables.
    !> No arg: detect node topology via MPI_COMM_TYPE_SHARED.
    !> With worker_in: count-based split using the given worker_inQtask value.
    !> Call after MPI__Initialize; pair with MPI__FreeQgroups.
    integer, intent(in), optional :: worker_in
    integer :: color
    if (present(worker_in)) then
      worker_inQtask = worker_in
      color = mpi__rank / worker_inQtask
      call mpi_comm_split(comm, color, mpi__rank, comm_q, mpi__info)
    else
      call MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, mpi__rank, MPI_INFO_NULL, comm_q, mpi__info)
      call MPI_Comm_size(comm_q, worker_inQtask, mpi__info)
    endif
    call MPI_Comm_rank(comm_q, mpi__rank_q, mpi__info)
    mpi__size_q = worker_inQtask
    mpi__root_q = mpi__rank_q == 0
    n_qgroup  = mpi__size / worker_inQtask
    iq_qgroup = mpi__rank / worker_inQtask
  end subroutine MPI__InitQgroups

  subroutine MPI__SplitXq(n_bpara, n_kpara)
    !> k-priority split of comm_q into comm_k_xq (size=n_kpara) and comm_b_xq (size=n_bpara).
    !> Layout: rank_q = rank_b*n_kpara + rank_k → color_k=rank_q/n_kpara, color_b=mod(rank_q,n_kpara).
    !> Special cases: n_bpara=1 → all-k (full k-parallel); n_kpara=1 → all-b (full ω-parallel).
    !> Sets backward-compat aliases comm_k/b, rank_k/b, size_k/b, root_k/b.
    !> Pair with MPI__FreeXq.
    integer, intent(in) :: n_bpara, n_kpara
    integer :: color
    ! comm_b_xq: ranks sharing the same rank_k index → same b-group (size = n_bpara)
    color = mod(mpi__rank_q, n_kpara)
    call mpi_comm_split(comm_q, color, mpi__rank_q, comm_b_xq, mpi__info)
    call mpi_comm_rank(comm_b_xq, mpi__rank_b_xq, mpi__info)
    call mpi_comm_size(comm_b_xq, mpi__size_b_xq, mpi__info)
    mpi__root_b_xq = mpi__rank_b_xq == 0
    ! comm_k_xq: ranks sharing the same rank_b index → same k-group (size = n_kpara)
    color = mpi__rank_q / n_kpara
    call mpi_comm_split(comm_q, color, mpi__rank_q, comm_k_xq, mpi__info)
    call mpi_comm_rank(comm_k_xq, mpi__rank_k_xq, mpi__info)
    call mpi_comm_size(comm_k_xq, mpi__size_k_xq, mpi__info)
    mpi__root_k_xq = mpi__rank_k_xq == 0
    ! comm_root_k_xq: k-roots across b-groups (for GatherXqw)
    color = merge(0, MPI_UNDEFINED, mpi__rank_k_xq == 0)
    call mpi_comm_split(comm_q, color, mpi__rank_q, comm_root_k_xq, mpi__info)
    if (comm_root_k_xq /= MPI_COMM_NULL) then
      call mpi_comm_rank(comm_root_k_xq, mpi__rank_root_k_xq, mpi__info)
      call mpi_comm_size(comm_root_k_xq, mpi__size_root_k_xq, mpi__info)
    endif
    if(ipr) write(06,'(X,A,5I5,3L2)') &
      "MPI(Xq): rank rank_q rank_k_xq rank_b_xq n_bpara n_kpara root_q root_k root_b", &
      mpi__rank, mpi__rank_q, mpi__rank_k_xq, mpi__rank_b_xq, n_bpara, n_kpara, &
      mpi__root_q, mpi__root_k_xq, mpi__root_b_xq
  end subroutine MPI__SplitXq

  subroutine MPI__FreeXq()
    !> Free comm_k_xq, comm_b_xq, comm_root_k_xq created by MPI__SplitXq.
    implicit none
    call mpi_comm_free(comm_k_xq, mpi__info)
    call mpi_comm_free(comm_b_xq, mpi__info)
    if (comm_root_k_xq /= MPI_COMM_NULL) call mpi_comm_free(comm_root_k_xq, mpi__info)
    if (allocated(mpi__npr_col)) deallocate(mpi__npr_col)
    if (allocated(mpi__ipr_col)) deallocate(mpi__ipr_col)
  end subroutine MPI__FreeXq

  subroutine MPI__SplitSxc(n_bpara, n_kpara)
    !> k-priority split of comm_q into comm_k_sxc (size=n_kpara) and comm_b_sxc (size=n_bpara).
    !> Same layout as SplitXq: rank_q = rank_b*n_kpara + rank_k.
    !> Special cases: n_bpara=1 → all-k (full k-parallel); n_kpara=1 → all-b.
    !> Sets backward-compat alias comm_k/rank_k/size_k/root_k to Sxc k-versions.
    !> Pair with MPI__FreeSxc.
    integer, intent(in) :: n_bpara, n_kpara
    integer :: color
    ! comm_b_sxc: ranks sharing the same rank_k → same b-group (size = n_bpara)
    color = mod(mpi__rank_q, n_kpara)
    call mpi_comm_split(comm_q, color, mpi__rank_q, comm_b_sxc, mpi__info)
    call mpi_comm_rank(comm_b_sxc, mpi__rank_b_sxc, mpi__info)
    call mpi_comm_size(comm_b_sxc, mpi__size_b_sxc, mpi__info)
    mpi__root_b_sxc = mpi__rank_b_sxc == 0
    ! comm_k_sxc: ranks sharing the same rank_b → same k-group (size = n_kpara)
    color = mpi__rank_q / n_kpara
    call mpi_comm_split(comm_q, color, mpi__rank_q, comm_k_sxc, mpi__info)
    call mpi_comm_rank(comm_k_sxc, mpi__rank_k_sxc, mpi__info)
    call mpi_comm_size(comm_k_sxc, mpi__size_k_sxc, mpi__info)
    mpi__root_k_sxc = mpi__rank_k_sxc == 0
    if(ipr) write(06,'(X,A,5I5,3L2)') &
      "MPI(Sxc): rank rank_q rank_k_sxc rank_b_sxc n_bpara n_kpara root_q root_k root_b", &
      mpi__rank, mpi__rank_q, mpi__rank_k_sxc, mpi__rank_b_sxc, n_bpara, n_kpara, &
      mpi__root_q, mpi__root_k_sxc, mpi__root_b_sxc
  end subroutine MPI__SplitSxc

  subroutine MPI__FreeSxc()
    !> Free comm_k_sxc, comm_b_sxc created by MPI__SplitSxc.
    implicit none
    call mpi_comm_free(comm_k_sxc, mpi__info)
    call mpi_comm_free(comm_b_sxc, mpi__info)
  end subroutine MPI__FreeSxc

  subroutine MPI__FreeQgroups()
    !> Free comm_q created by MPI__InitQgroups. Call once at end of hgw.
    implicit none
    call mpi_comm_free(comm_q, mpi__info)
  end subroutine MPI__FreeQgroups

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
    if(.not.allocated(mpi__npr_col)) allocate(mpi__npr_col(0:mpi__size_b_xq-1))
    if(.not.allocated(mpi__ipr_col)) allocate(mpi__ipr_col(0:mpi__size_b_xq-1))
    ipr_col = 1
    do irank_b = 0, mpi__size_b_xq - 1
      npr_col = (npr + irank_b)/mpi__size_b_xq
      mpi__npr_col(irank_b) = npr_col
      mpi__ipr_col(irank_b) = ipr_col
      ipr_col = ipr_col + npr_col
      if (npr_col == 0) call rx("MPI__Setnpr_col: use small parallelization")
    enddo
    npr_col = mpi__npr_col(mpi__rank_b_xq)
    ipr_col = mpi__ipr_col(mpi__rank_b_xq)
    if(ipr)write(06,'(X,A,I5,2I7)') "mpi__rank_b_xq, ipr_col, npr_col=", mpi__rank_b_xq, ipr_col, npr_col
  end subroutine MPI__Setnpr_col

  subroutine MPI__GatherXqw(xqw, xqw_all, npr, npr_col, collector_rank)
    integer, intent(in) :: npr, npr_col
    integer, intent(in), optional :: collector_rank
    complex(8), intent(in) :: xqw(npr,npr_col)
    complex(8), intent(inout) :: xqw_all(npr,npr)
    integer, allocatable :: data_disp(:), data_size(:)
    integer :: irank_b, collector_rank_in
    if(mpi__size_b_xq == 1 .and. npr == npr_col) then
      xqw_all(:,:) = xqw(:,:)
      return
    endif
    allocate(data_size(0:mpi__size_b_xq-1), data_disp(0:mpi__size_b_xq-1))
    do irank_b = 0, mpi__size_b_xq -1
      data_size(irank_b) = npr*mpi__npr_col(irank_b)
      data_disp(irank_b) = npr*(mpi__ipr_col(irank_b)-1)
    enddo
    collector_rank_in = 0
    if(present(collector_rank)) collector_rank_in = collector_rank
    call mpi_gatherv(xqw, npr*npr_col, mpi_complex16, xqw_all, data_size, data_disp, &
                  &  mpi_complex16, collector_rank_in, comm_root_k_xq, mpi__info)
    deallocate(data_size, data_disp)
  end subroutine MPI__GatherXqw
  subroutine MPI__GatherXqw_c(xqw, xqw_all, npr, npr_col, collector_rank)
    integer, intent(in) :: npr, npr_col
    integer, intent(in), optional :: collector_rank
    complex(4), intent(in) :: xqw(npr,npr_col)
    complex(4), intent(out) :: xqw_all(npr,npr)
    integer, allocatable :: data_disp(:), data_size(:)
    integer :: irank_b, collector_rank_in
    if(mpi__size_b_xq == 1 .and. npr == npr_col) then
      xqw_all(:,:) = xqw(:,:)
      return
    endif
    allocate(data_size(0:mpi__size_b_xq-1), data_disp(0:mpi__size_b_xq-1))
    do irank_b = 0, mpi__size_b_xq -1
      data_size(irank_b) = npr*mpi__npr_col(irank_b)
      data_disp(irank_b) = npr*(mpi__ipr_col(irank_b)-1)
    enddo
    collector_rank_in = 0
    if(present(collector_rank)) collector_rank_in = collector_rank
    call mpi_gatherv(xqw, npr*npr_col, mpi_complex, xqw_all, data_size, data_disp, &
                  &  mpi_complex, collector_rank_in, comm_root_k_xq, mpi__info)
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
  subroutine MPI__consoleout(idn)
    use m_lgunit,only:stdo,stdl
    logical, save :: init = .false.
    character(1024*4) :: cwd, stdout
    character*(*):: idn
    if(init) return
    init = .true.
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

  subroutine xmpbnd2(kpproc, ndham, ndat, eb)
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
