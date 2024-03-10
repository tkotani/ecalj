module m_mpi !MPI utility for fpgw
  implicit none
  include "mpif.h"
  integer :: mpi__size
  integer :: mpi__rank
  logical :: mpi__root
  integer :: comm
  integer :: mpi__sizeMG 
  integer :: mpi__rankMG
  integer,private :: mpi__info
  integer,private:: ista(MPI_STATUS_SIZE )
contains
  subroutine MPI__Initialize(commin)
    implicit none
    character(1024*4) :: cwd, stdout
    integer,optional:: commin
    comm= merge(commin,MPI_COMM_WORLD,present(commin))
    call getcwd(cwd)           ! get current working directory
    call MPI_Init( mpi__info ) ! current working directory is changed if mpirun is not used
    call MPI_Comm_rank( comm, mpi__rank, mpi__info )
    call MPI_Comm_size( comm, mpi__size, mpi__info )
    mpi__root= mpi__rank==0
    if( mpi__root ) call chdir(cwd)        ! recover current working directory
  end subroutine MPI__Initialize
  subroutine MPI__Initialize_magnon(commin)
    implicit none
    character(1024*4) :: cwd, stdout
    integer,optional:: commin
    comm= merge(commin,MPI_COMM_WORLD,present(commin))
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
  subroutine MPI__AllreduceSum( data, sizex )
    implicit none
    integer, intent(in) :: sizex
    complex(8), intent(inout) :: data(sizex)
    complex(8), allocatable   :: mpi__data(:) 
    if( mpi__size == 1 ) return
    allocate(mpi__data(sizex))
    mpi__data = data
    call MPI_Allreduce( mpi__data, data, sizex, MPI_DOUBLE_COMPLEX, MPI_SUM, comm, mpi__info )
    deallocate( mpi__data )
  end subroutine MPI__AllreduceSum
  subroutine MPI__reduceSum( root, data, sizex )
    implicit none
    integer, intent(in) :: sizex,root
    complex(8), intent(inout) :: data(sizex)
    complex(8), allocatable   :: mpi__data(:) 
    if( mpi__size == 1 ) return
    allocate(mpi__data(sizex))
    mpi__data = data
    call MPI_reduce( mpi__data, data, sizex, MPI_DOUBLE_COMPLEX, MPI_SUM, root, comm, mpi__info )
    deallocate( mpi__data )
    return
  end subroutine MPI__reduceSum
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
end module m_mpi

subroutine MPI__sxcf_rankdivider(irkip_all,nspinmx,nqibz,ngrp,nq,irkip)
  use m_mpi,only: mpi__rank,mpi__size
  implicit none
  integer, intent(out) :: irkip    (nspinmx,nqibz,ngrp,nq)
  integer, intent(in)  :: irkip_all(nspinmx,nqibz,ngrp,nq)
  integer, intent(in)  :: nspinmx,nqibz,ngrp,nq
  integer :: ispinmx,iqibz,igrp,iq
  integer :: total
  integer, allocatable :: vtotal(:)
  integer :: indexi, indexe
  integer :: p
  if( mpi__size == 1 ) then
     irkip = irkip_all
     return
  end if
  total = count(irkip_all>0)
  write(6,"('MPI__sxcf_rankdivider:$')")
  write(6,"('nspinmx,nqibz,ngrp,nq,total=',5i6)") nspinmx,nqibz,ngrp,nq,total 
  allocate( vtotal(0:mpi__size-1) )
  vtotal(:) = total/mpi__size
  do p=1, mod(total,mpi__size)
     vtotal(p-1) = vtotal(p-1) + 1
  end do
  indexe=0
  indexi=-999999
  do p=0, mpi__rank
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
