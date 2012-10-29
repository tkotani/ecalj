      module m_mpi
      implicit none
      include "mpif.h"

      integer :: mpi__info
      integer :: mpi__size
      integer :: mpi__rank
      logical :: mpi__root
      integer :: mpi__iini, mpi__iend
      logical, allocatable ::  mpi__task(:)

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
      end if

c$$$      if( mpi__size > 1 ) then
c$$$         write(stdout,"('output-console.s',i3.3)") mpi__rank
c$$$         open(unit=6,file=stdout)
c$$$      end if
c$$$      if( .not. mpi__root ) then
c$$$         close(unit=6)
c$$$      end if


      return
      end subroutine MPI__Initialize


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


      subroutine MPI__AllreduceSum( data, size )
      implicit none
      complex(8), intent(inout) :: data(size)
      integer, intent(in) :: size
      complex(8), allocatable   :: mpi__data(:) 

      if( mpi__size == 1 ) return

      allocate(mpi__data(size))
      mpi__data = data

      call MPI_Allreduce( mpi__data, data, size,
     &     MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, mpi__info )

      deallocate( mpi__data )

      return
      end subroutine MPI__AllreduceSum


      subroutine MPI__AllreduceMax( data, size )
      implicit none
      integer, intent(inout) :: data(size)
      integer, intent(in) :: size
      integer, allocatable   :: mpi__data(:) 

      if( mpi__size == 1 ) return

      allocate(mpi__data(size))
      mpi__data = data

      call MPI_Allreduce( mpi__data, data, size,
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


      subroutine MPI__hx0fp0_rankdivider(iqxini,iqxend,nqibz)
      implicit none
      integer, intent(in) :: iqxini, iqxend, nqibz

      integer :: iq
      allocate( mpi__task(1:iqxend) )
      
      if( mpi__size == 1 ) then
         mpi__task(:) = .true.
         return
      end if

      mpi__task(:) = .false.
      do iq=iqxini, iqxend
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

      end module m_mpi
