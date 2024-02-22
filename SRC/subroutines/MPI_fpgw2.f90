module m_mpi !MPI utility for fpgw
  implicit none
  include "mpif.h"
  integer :: mpi__info
  integer :: mpi__size
  integer :: mpi__sizeMG = 1
  integer :: mpi__rank
  integer :: mpi__rankMG
  integer :: mpi__rankQ
  logical :: mpi__root
  logical :: mpi__rootQ
  integer :: mpi__comm = MPI_COMM_WORLD
  integer :: mpi__commQ
  integer:: ista(MPI_STATUS_SIZE )
  integer :: mpi__iini, mpi__iend
  logical, allocatable :: mpi__task(:)
  integer, allocatable :: mpi__ranktab(:)
  integer  :: mpi__MEq ! for magnon E(q) parallelization
contains
  subroutine MPI__Initialize
    implicit none
    character(1024*4) :: cwd, stdout
    call getcwd(cwd)          ! get current working directory
    call MPI_Init( mpi__info ) ! current working directory is changed if mpirun is not used
    call MPI_Comm_rank( MPI_COMM_WORLD, mpi__rank, mpi__info )
    call MPI_Comm_size( MPI_COMM_WORLD, mpi__size, mpi__info )
    if( mpi__rank == 0 ) then
       mpi__root = .true.
    else
       mpi__root = .false.
    end if
    if( mpi__root ) then
       call chdir(cwd)        ! recover current working directory
    endif
  end subroutine MPI__Initialize
  !     !======================================================
  subroutine MPI__Initialize_magnon()
    implicit none
    character(1024*4) :: cwd, stdout
    call getcwd(cwd)          ! get current working directory
    call MPI_Init( mpi__info ) ! current working directory is changed if mpirun is not used
    call MPI_Comm_rank( MPI_COMM_WORLD, mpi__rankMG, mpi__info )
    call MPI_Comm_size( MPI_COMM_WORLD, mpi__sizeMG, mpi__info )
    if( mpi__rankMG == 0 ) then
       mpi__root = .true.
    else
       mpi__root = .false.
    end if
    if( mpi__root ) then
       call chdir(cwd)        ! recover current working directory
    endif
    return
  end subroutine MPI__Initialize_magnon
  subroutine MPI__consoleout(idn)
    use m_lgunit,only:stdo,stdl
    implicit none
    character(1024*4) :: cwd, stdout
    character*(*):: idn
    ! console-output from different nodes to different files
    if( mpi__size > 1 ) then
       if( mpi__root ) then
          write(6,"(' MPI outputs in each rank are in stdout.{RankId}.',a)")idn
          call flush(stdo)
       end if
       !        write(6,"('   stdout.',i4.4,'.',a)") mpi__rank,idn
       write(stdout,"('stdout.',i4.4,'.',a)") mpi__rank,idn
       open(unit=6,file=trim(stdout))
       write(6,"(a,i3)")" ### console output for rank=",mpi__rank
    endif
    return
  end subroutine MPI__consoleout
  subroutine MPI__consoleout_magnon(idn,size_lim)
    use m_lgunit,only:stdo,stdl
    implicit none
    character(1024*4) :: cwd, stdout
    character*(*):: idn
    integer , intent(in) :: size_lim
    !C Reduce size for magnon (avoid memory leak)      
    if (mpi__sizeMG > size_lim) then
       mpi__sizeMG=size_lim
    endif
    ! console-output from different nodes to different files
    if( mpi__sizeMG > 1 ) then
       if( mpi__root ) then
          write(6,"(' MPI outputs in each rank are in stdout.{RankId}.',a)")idn
          call flush(stdo)
       end if
       !        write(6,"('   stdout.',i4.4,'.',a)") mpi__rank,idn
       write(stdout,"('stdout.',i4.4,'.',a)") mpi__rankMG,idn
       open(unit=6,file=trim(stdout))
       write(6,"(a,i3)")" ### console output for rank=",mpi__rankMG
    endif
    return
  end subroutine MPI__consoleout_magnon
  subroutine MPI__Barrier
    implicit none
    call MPI_Barrier( MPI_COMM_WORLD, mpi__info )
  end subroutine MPI__Barrier
  subroutine MPI__getRange( mpi__indexi, mpi__indexe, indexi, indexe )
    implicit none
    integer, intent(out) :: mpi__indexi, mpi__indexe
    integer, intent(in)  :: indexi, indexe
    integer, allocatable :: mpi__total(:)
    integer              :: total
    integer :: p
    allocate( mpi__total(0:mpi__size-1) )
    total = indexe-indexi+1
    mpi__total(:) = total/mpi__size
    do p=1, mod(total,mpi__size)
       mpi__total(p-1) = mpi__total(p-1) + 1
    end do
    mpi__indexe=indexi-1
    do p=0, mpi__rank
       mpi__indexi = mpi__indexe+1
       mpi__indexe = mpi__indexi+mpi__total(p)-1
    end do
    deallocate(mpi__total)
  end subroutine MPI__getRange
  subroutine MPI__Broadcast( data )
    implicit none
    integer, intent(inout) :: data
    call MPI_Bcast( data, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi__info )
  end subroutine MPI__Broadcast
  subroutine MPI__REAL8send(data,n,dest)
    implicit none
    real(8):: data(n)
    integer :: n,dest,ierr
    call MPI_Send(data,n,MPI_REAL8,dest,mpi__rank, MPI_COMM_WORLD,ierr)
  end subroutine MPI__REAL8send
  subroutine MPI__REAL8recv(data,n,src)
    implicit none
    real(8):: data(n)
    integer :: n,src,ierr
    call MPI_Recv(data,n,MPI_REAL8,src,src, MPI_COMM_WORLD,ista,ierr)
  end subroutine MPI__REAL8recv
  subroutine MPI__DbleCOMPLEXsend(data,n,dest)
    implicit none
    complex(8):: data(n)
    integer :: n,dest,ierr
    call MPI_Send(data,n,MPI_COMPLEX16,dest,mpi__rank, MPI_COMM_WORLD,ierr)
  end subroutine MPI__DbleCOMPLEXsend
  subroutine MPI__DbleCOMPLEXrecv(data,n,src)
    implicit none
    complex(8):: data(n)
    integer :: n,src,ierr
    call MPI_Recv(data,n,MPI_COMPLEX16,src,src, MPI_COMM_WORLD,ista,ierr)
  end subroutine MPI__DbleCOMPLEXrecv
  subroutine MPI__DbleCOMPLEXsendQ(data,n,destQ)
    implicit none
    complex(8):: data(n)
    integer :: n,destQ,ierr
    call MPI_Send(data,n,MPI_COMPLEX16,destQ,0,mpi__commQ,ierr)
  end subroutine MPI__DbleCOMPLEXsendQ
  subroutine MPI__DbleCOMPLEXrecvQ(data,n,srcQ)
    implicit none
    complex(8):: data(n)
    integer :: n,srcQ,ierr
    call MPI_Recv(data,n,MPI_COMPLEX16,srcQ,0,mpi__commQ,ista,ierr)
  end subroutine MPI__DbleCOMPLEXrecvQ
  subroutine MPI__AllreduceSum( data, sizex )
    implicit none
    integer, intent(in) :: sizex
    complex(8), intent(inout) :: data(sizex)
    complex(8), allocatable   :: mpi__data(:) 
    if( mpi__size == 1 ) return
    allocate(mpi__data(sizex))
    mpi__data = data
    call MPI_Allreduce( mpi__data, data, sizex,&
         MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, mpi__info )
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
    call MPI_reduce( mpi__data, data, sizex,&
         MPI_DOUBLE_COMPLEX, MPI_SUM, root, MPI_COMM_WORLD, mpi__info )
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
    call MPI_Allreduce( mpi__data, data, sizex,&
         MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, mpi__info )
    deallocate( mpi__data )
  end subroutine MPI__AllreduceMax
  !=========================================================================
  subroutine MPI__sxcf_rankdivider(irkip_all,nspinmx,nqibz,ngrp,nq,irkip)
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
  subroutine MPI__hmagnon_rankdivider(nqbz)
    implicit none
    integer, intent(in) :: nqbz
    integer :: iq,i
    allocate( mpi__task(1:nqbz),mpi__ranktab(1:nqbz) )
    mpi__task(:) = .false.
    mpi__ranktab(1:nqbz)=999999
    write(6,*) "mpi__sizeMG:",mpi__sizeMG
    mpi__MEq=1+nqbz/mpi__sizeMG
    write(6,*) "mpi__sizeMG:",mpi__MEq
    if( mpi__sizeMG == 1 ) then
       mpi__task(:) = .true.
       mpi__ranktab(:) = mpi__rankMG
       return
    endif
    if(mpi__rankMG==0) write(6,*) "MPI_hmagnon_rankdivider:"
    do iq=1,nqbz
       mpi__ranktab(iq) = mod(iq-1,mpi__sizeMG)  !rank_table for given iq. iq=1 must give rank=0
       if( mpi__ranktab(iq) == mpi__rankMG) then
          mpi__task(iq) = .true.               !mpi__task is nodeID-dependent.
       endif
       if(mpi__rankMG==0) then
          write(6,"('  iq irank=',2i5)")iq,mpi__ranktab(iq)
       endif
    enddo
    return
  end subroutine MPI__hmagnon_rankdivider
end module m_mpi

