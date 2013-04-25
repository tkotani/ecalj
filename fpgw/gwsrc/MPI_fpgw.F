      module m_mpi
      implicit none
      include "mpif.h"

      integer :: mpi__info
      integer :: mpi__size
      integer :: mpi__rank
      logical :: mpi__root
      integer :: mpi__iini, mpi__iend
      logical, allocatable ::  mpi__task(:)
      integer,allocatable :: mpi__ranktab(:)
      integer:: ista(MPI_STATUS_SIZE )

      interface MPI__Send 
      module procedure
     &  MPI__Send_i,  MPI__Send_iv, 
     &  MPI__Send_d,  MPI__Send_dv
      end interface

      interface MPI__Recv
      module procedure
     &  MPI__Recv_i,  MPI__Recv_iv, 
     &  MPI__Recv_d,  MPI__Recv_dv 
      end interface

      contains

      
      subroutine MPI__Initialize()
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
c! console-output from different nodes to different files
c      if( mpi__size > 1 ) then
c        if(mpi__root )write(*,*)'MPI console outputs to following files.'
c        write(*,"('   stdout.',i4.4,'.',a)") mpi__rank,idn
c        write(stdout,"('stdout.',i4.4,'.',a)") mpi__rank,idn
c        open(unit=6,file=trim(stdout))
c        write(*,"(a,i3)")" ### console output for rank=",mpi__rank
c      endif
c      if(mpi__root ) then
c        close(unit=6)
c      endif
      return
      end subroutine MPI__Initialize

      subroutine MPI__consoleout(idn)
      implicit none
      character(1024*4) :: cwd, stdout
      character*(*):: idn
! console-output from different nodes to different files
      if( mpi__size > 1 ) then
        if(mpi__root )write(*,*)'MPI console outputs to following files.'
        write(*,"('   stdout.',i4.4,'.',a)") mpi__rank,idn
        write(stdout,"('stdout.',i4.4,'.',a)") mpi__rank,idn
        open(unit=6,file=trim(stdout))
        write(*,"(a,i3)")" ### console output for rank=",mpi__rank
      endif
      return
      end subroutine MPI__consoleout


      subroutine MPI__Barrier
      implicit none

      call MPI_Barrier( MPI_COMM_WORLD, mpi__info )

      end subroutine MPI__Barrier

      subroutine MPI__Finalize
      implicit none

      call MPI_Finalize ( mpi__info )

      end subroutine MPI__Finalize



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

      return
      end subroutine MPI__getRange


      subroutine MPI__Broadcast( data )
      implicit none
      integer, intent(inout) :: data

      call MPI_Bcast( data, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi__info )

      return
      end subroutine MPI__Broadcast


      subroutine MPI__send_d(data,dest)
      implicit none
      real(8):: data
      integer :: n,dest,ierr
      n=1
      call MPI_Send(data,n,MPI_REAL8,dest,mpi__rank, MPI_COMM_WORLD,ierr)
      end subroutine MPI__send_d

      subroutine MPI__recv_d(data,src)
      implicit none
      real(8):: data
      integer :: n,src,ierr
      n=1
      call MPI_Recv(data,n,MPI_REAL8,src,src, MPI_COMM_WORLD,ista,ierr)
      end subroutine MPI__recv_d

      subroutine MPI__send_dv(data,dest)
      implicit none
      real(8):: data(:)
      integer :: n,dest,ierr
      n=size(data)
      call MPI_Send(data,n,MPI_REAL8,dest,mpi__rank, MPI_COMM_WORLD,ierr)
      end subroutine MPI__send_dv

      subroutine MPI__recv_dv(data,src)
      implicit none
      real(8):: data(:)
      integer :: n,src,ierr
      n=size(data)
      call MPI_Recv(data,n,MPI_REAL8,src,src, MPI_COMM_WORLD,ista,ierr)
      end subroutine MPI__recv_dv

      subroutine MPI__send_i(data,dest)
      implicit none
      integer:: data
      integer :: n,dest,ierr
      n=1
      call MPI_Send(data,n,MPI_INTEGER,dest,mpi__rank, MPI_COMM_WORLD,ierr)
      end subroutine MPI__send_i

      subroutine MPI__recv_i(data,src)
      implicit none
      integer:: data
      integer :: n,src,ierr
      n=1
      call MPI_Recv(data,n,MPI_INTEGER,src,src, MPI_COMM_WORLD,ista,ierr)
      end subroutine MPI__recv_i

      subroutine MPI__send_iv(data,dest)
      implicit none
      integer:: data(:)
      integer :: n,dest,ierr
      n=size(data)
      call MPI_Send(data,n,MPI_INTEGER,dest,mpi__rank, MPI_COMM_WORLD,ierr)
      end subroutine MPI__send_iv

      subroutine MPI__recv_iv(data,src)
      implicit none
      integer:: data(:)
      integer :: n,src,ierr
      n=size(data)
      call MPI_Recv(data,n,MPI_INTEGER,src,src, MPI_COMM_WORLD,ista,ierr)
      end subroutine MPI__recv_iv


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

      subroutine MPI__AllreduceSum( data, sizex )
      implicit none
      integer, intent(in) :: sizex
      complex(8), intent(inout) :: data(sizex)
      complex(8), allocatable   :: mpi__data(:) 

      if( mpi__size == 1 ) return

      allocate(mpi__data(sizex))
      mpi__data = data

      call MPI_Allreduce( mpi__data, data, sizex,
     &     MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, mpi__info )

      deallocate( mpi__data )

      return
      end subroutine MPI__AllreduceSum


      subroutine MPI__AllreduceMax( data, sizex )
      implicit none
      integer, intent(in) :: sizex
      integer, intent(inout) :: data(sizex)
      integer, allocatable   :: mpi__data(:) 

      if( mpi__size == 1 ) return

      allocate(mpi__data(sizex))
      mpi__data = data

      call MPI_Allreduce( mpi__data, data, sizex,
     &     MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, mpi__info )

      deallocate( mpi__data )

      return
      end subroutine MPI__AllreduceMax


      subroutine MPI__sxcf_rankdivider(irkip,irkip_all,nspinmx,nqibz,ngrp,nq)
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

c$$$      write(*,*) "irkip_all", mpi__rank
      total = 0
      do ispinmx=1, nspinmx
         do iq=1, nq
            do iqibz=1, nqibz
               do igrp=1, ngrp
                  if( irkip_all(ispinmx,iqibz,igrp,iq) /= 0 ) then
                     total = total + 1
c$$$                     write(*,*) ispinmx, iq, iqibz, igrp, irkip_all(ispinmx,iqibz,igrp,iq)
                  end if
               end do
            end do
         end do
      end do
      write(6,"('MPI__sxcf_rankdivider: nspinmx,nqibz,ngrp,nq,total=',5i6)") nspinmx,nqibz,ngrp,nq,total 
      
      allocate( vtotal(0:mpi__size-1) )
      vtotal(:) = total/mpi__size

      do p=1, mod(total,mpi__size)
         vtotal(p-1) = vtotal(p-1) + 1
      end do

      indexe=0
      do p=0, mpi__rank
         indexi = indexe+1
         indexe = indexi+vtotal(p)-1
      end do
      deallocate(vtotal)

      total = 0
      irkip(:,:,:,:) = 0

c$$$      write(*,*) "irkip", mpi__rank
      do ispinmx=1, nspinmx
         do iq=1, nq
            do iqibz=1, nqibz
               do igrp=1, ngrp
                  if( irkip_all(ispinmx,iqibz,igrp,iq) /= 0 ) then
                     total = total + 1
                     if( indexi<=total .and. total<=indexe ) then
                        irkip(ispinmx,iqibz,igrp,iq) = irkip_all(ispinmx,iqibz,igrp,iq)
c$$$                        write(*,*) ispinmx, iq, iqibz, igrp, irkip(ispinmx,iqibz,igrp,iq)
                     end if
                  end if
               end do
            end do
         end do
      end do

      return
      end subroutine MPI__sxcf_rankdivider

!!
      subroutine MPI__hx0fp0_rankdivider(iqxini,iqxend,nqibz)
      implicit none
      integer, intent(in) :: iqxini, iqxend, nqibz

      integer :: iq
      allocate( mpi__task(1:iqxend),mpi__ranktab(1:iqxend) )
      
      if( mpi__size == 1 ) then
         mpi__task(:) = .true.
         mpi__ranktab(iq) = mpi__rank
         return
      end if
!!
      mpi__task(:) = .false.
      do iq=iqxini, iqxend
        if(iq==1.or. iq>nqibz) then
           mpi__ranktab(iq) = 0
        else
           mpi__ranktab(iq) = mod(iq,mpi__size-1)+1
        endif  
!!
        if( mpi__rank == 0 ) then
            if( iq == 1 .or. iq>nqibz ) then
               mpi__task(iq) = .true.
            else
               mpi__task(iq) = .false.
            end if
        else
            if( iq == 1 .or. iq>nqibz ) then
               mpi__task(iq) = .false.
            else
               if( mpi__rank == mod(iq,mpi__size-1)+1 ) then
                  mpi__task(iq) = .true.
               else
                  mpi__task(iq) = .false.
               end if
            end if
        end if
      end do

      return
      end subroutine MPI__hx0fp0_rankdivider


      subroutine MPI__hx0fp0_rankdivider2(iqxini,iqxend)
      implicit none
      integer, intent(in) :: iqxini, iqxend
      integer :: iq,i
      allocate( mpi__task(1:iqxend),mpi__ranktab(1:iqxend) )
      mpi__task(:) = .false.
      mpi__ranktab(1:iqxend)=999999
      if( mpi__size == 1 ) then
         mpi__task(:) = .true.
         mpi__ranktab(iq) = mpi__rank
         return
      end if
      if(mpi__rank==0) write(6,*) "MPI_hx0fp0_rankdivider2:"
      do iq=iqxini, iqxend
         mpi__ranktab(iq) = mod(iq-1,mpi__size)  !rank_table for given iq. iq=1 must give rank=0
         if( mpi__ranktab(iq) == mpi__rank) then
            mpi__task(iq) = .true.               !mpi__task is nodeID-dependent.
         endif
         if(mpi__rank==0) then
           write(6,"('  iq irank=',2i5)")iq,mpi__ranktab(iq)
         endif
      enddo   
      return
      end subroutine MPI__hx0fp0_rankdivider2

      end module m_mpi
