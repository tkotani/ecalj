c#undef SCALAPACK
      module m_mpi
      implicit none
      include "mpif.h"

      integer :: mpi__info
      integer :: mpi__size
      integer :: mpi__sizeMG = 1
      integer :: mpi__sizeQ = 1
      integer :: mpi__sizeS = 1
      integer :: mpi__sizeP = 1
      integer :: mpi__sizeB = 1
      integer :: mpi__sizeM = 1
      integer :: mpi__sizeW = 1

      integer :: mpi__rank
      integer :: mpi__rankMG
      integer :: mpi__rankQ
c      integer :: mpi__rankS
c      integer :: mpi__rankP
c      integer :: mpi__rankB
c      integer :: mpi__rankM
c      integer :: mpi__rankW

      logical :: mpi__root
      logical :: mpi__rootQ
c      logical :: mpi__rootS
c      logical :: mpi__rootP
c      logical :: mpi__rootB
c      logical :: mpi__rootM

      integer :: mpi__comm = MPI_COMM_WORLD
      integer :: mpi__commQ
c      integer :: mpi__commS
c      integer :: mpi__commP
c      integer :: mpi__commB
c      integer :: mpi__commM

      integer:: ista(MPI_STATUS_SIZE )

      integer :: mpi__iini, mpi__iend

      logical, allocatable :: mpi__task(:)
      integer, allocatable :: mpi__ranktab(:)

C     ! for Q(nq) parallelization
      logical, allocatable :: mpi__Qtask(:)
      integer, allocatable :: mpi__Qranktab(:)

C     ! for S(nspin) parallelization
      logical, allocatable :: mpi__Stask(:)
      integer, allocatable :: mpi__Sranktab(:)
      integer :: mpi__Snall
      integer :: mpi__Sn, mpi__Ss, mpi__Se
      integer, allocatable :: mpi__Svn(:), mpi__Svs(:), mpi__Sve(:)


C     ! for P(npm) parallelization
      integer :: mpi__Pnall
      integer :: mpi__Pn, mpi__Ps, mpi__Pe
      integer, allocatable :: mpi__Pvn(:), mpi__Pvs(:), mpi__Pve(:)


C     ! for B(Kpoint&nbnb) parallelization
      integer :: mpi__Bnall
      integer :: mpi__Bn, mpi__Bs, mpi__Be
      integer, allocatable :: mpi__Bvn(:), mpi__Bvs(:), mpi__Bve(:)
      logical, allocatable :: mpi__Btask1(:)
      logical, allocatable :: mpi__Btask2(:,:)
      logical, allocatable :: mpi__Btask3(:,:,:)

C     ! for M(matrix) parallelization
      integer :: mpi__Mnall
      integer :: mpi__Mn, mpi__Ms, mpi__Me, mpi__Mf
      integer, allocatable :: mpi__Mvn(:), mpi__Mvs(:), mpi__Mve(:)
      integer :: mpi__Mhandle=0
      integer :: mpi__Mdescv(9)
      integer :: mpi__Mdescm(9)

C     ! for W(iwt) parallelization
      integer :: mpi__Wnall
      integer :: mpi__Wn, mpi__Ws, mpi__We
      integer, allocatable :: mpi__Wvn(:), mpi__Wvs(:), mpi__Wve(:)

C     ! for magnon E(q) parallelization
      integer  :: mpi__MEq


      interface MPI__Send 
      module procedure
     &     MPI__Send_i,  MPI__Send_iv, 
     &     MPI__Send_d,  MPI__Send_dv
      end interface

      interface MPI__Recv
      module procedure
     &     MPI__Recv_i,  MPI__Recv_iv, 
     &     MPI__Recv_d,  MPI__Recv_dv 
      end interface

      contains

C     !======================================================
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
c! console-output from different nodes to different files
c      if( mpi__size > 1 ) then
c        if(mpi__root )write(6,*)'MPI console outputs to following files.'
c        write(6,"('   stdout.',i4.4,'.',a)") mpi__rank,idn
c        write(stdout,"('stdout.',i4.4,'.',a)") mpi__rank,idn
c        open(unit=6,file=trim(stdout))
c        write(6,"(a,i3)")" ### console output for rank=",mpi__rank
c      endif
c      if(mpi__root ) then
c        close(unit=6)
c      endif
      return
      end subroutine MPI__Initialize

C     !======================================================
      subroutine MPI__Initialize_magnon()
      implicit none
      character(1024*4) :: cwd, stdout
      call getcwd(cwd)          ! get current working directory
c$$$CC Reduce size for magnon (avoid memory leak)      
c$$$      if (mpi__sizeMG > size_lim) then
c$$$         mpi__sizeMG=size_lim
c$$$      endif
c$$$      if (mpi__sizeMG > 10) then
c$$$         mpi__sizeMG=10
c$$$      endif
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
c! console-output from different nodes to different files
c      if( mpi__size > 1 ) then
c        if(mpi__root )write(6,*)'MPI console outputs to following files.'
c        write(6,"('   stdout.',i4.4,'.',a)") mpi__rank,idn
c        write(stdout,"('stdout.',i4.4,'.',a)") mpi__rank,idn
c        open(unit=6,file=trim(stdout))
c        write(6,"(a,i3)")" ### console output for rank=",mpi__rank
c      endif
c      if(mpi__root ) then
c        close(unit=6)
c      endif
      return
      end subroutine MPI__Initialize_magnon
      
C     !======================================================
      subroutine MPI__InitializeQ !MPI__InitializeQSPBM
      implicit none
      character(1024*4) :: cwd, stdout
      integer :: narg, iargc
      character(len=1024) :: arg
c      integer :: mpi__group
c      integer :: mpi__groupQ
c$$$      integer :: mpi__groupS
c$$$      integer :: mpi__groupP
c$$$      integer :: mpi__groupB
c$$$      integer :: mpi__groupM

      integer, allocatable :: ranklistQ(:)
c      integer, allocatable :: ranklistS(:)
c      integer, allocatable :: ranklistP(:)
c      integer, allocatable :: ranklistB(:)
c      integer, allocatable :: ranklistM(:)

c$$$      integer, allocatable :: vsizeQ(:)
c$$$      integer, allocatable :: vsizeS(:)
c$$$      integer, allocatable :: vsizeP(:)
c$$$      integer, allocatable :: vsizeB(:)
c$$$      integer, allocatable :: vsizeM(:)

      integer, allocatable :: vrankQ(:)
c      integer, allocatable :: vrankS(:)
c      integer, allocatable :: vrankP(:)
c      integer, allocatable :: vrankB(:)
c      integer, allocatable :: vrankM(:)
      integer :: i, j, n
c      integer, allocatable :: imap(:)

C     ! intialize MPI and
      call getcwd(cwd)          ! get current working directory

      call MPI_Init( mpi__info ) ! current working directory is changed if mpirun is not used
      call MPI_Comm_rank ( MPI_COMM_WORLD, mpi__rank,  mpi__info )
      call MPI_Comm_size ( MPI_COMM_WORLD, mpi__size,  mpi__info )

      mpi__root = ( mpi__rank == 0 )
      if( mpi__root ) then
         call chdir(cwd)        ! recover current working directory
      endif
      

c$$$C     ! setup S-parallel
c$$$C     !  commS for among other S, same Q,P,B,M
c$$$      allocate( ranklistS(0:mpi__sizeS-1) )
c$$$      i=0
c$$$      do j=0, mpi__size-1
c$$$         if(  vrankQ(j) == mpi__rankQ .and.
c$$$     .        vrankP(j) == mpi__rankP .and.
c$$$     .        vrankB(j) == mpi__rankB .and.
c$$$     .        vrankM(j) == mpi__rankM ) then
c$$$            ranklistS(i) = j
c$$$            i=i+1
c$$$         end if
c$$$      end do
c$$$      call MPI_Group_incl( mpi__group,
c$$$     .     mpi__sizeS, ranklistS, mpi__groupS, mpi__info )
c$$$      call MPI_Comm_create( MPI_COMM_WORLD,
c$$$     .     mpi__groupS, mpi__commS, mpi__info )
c$$$      deallocate( ranklistS )


c$$$C     ! setup P-parallel
c$$$C     !  commP for among other P, same Q,B,M
c$$$      allocate( ranklistP(0:mpi__sizeP-1) )
c$$$      i=0
c$$$      do j=0, mpi__size-1
c$$$         if(  vrankQ(j) == mpi__rankQ .and.
c$$$     .        vrankS(j) == mpi__rankS .and.
c$$$     .        vrankB(j) == mpi__rankB .and.
c$$$     .        vrankM(j) == mpi__rankM ) then
c$$$            ranklistP(i) = j
c$$$            i=i+1
c$$$         end if
c$$$      end do
c$$$      call MPI_Group_incl( mpi__group,
c$$$     .     mpi__sizeP, ranklistP, mpi__groupP, mpi__info )
c$$$      call MPI_Comm_create( MPI_COMM_WORLD,
c$$$     .     mpi__groupP, mpi__commP, mpi__info )
c$$$      deallocate( ranklistP )
      
c$$$C     ! setup B-parallel
c$$$C     !  commB for among other B, same Q,P,M
c$$$      allocate( ranklistB(0:mpi__sizeB-1) )
c$$$      i=0
c$$$      do j=0, mpi__size-1
c$$$         if(  vrankQ(j) == mpi__rankQ .and.
c$$$     .        vrankS(j) == mpi__rankS .and.
c$$$     .        vrankP(j) == mpi__rankP .and.
c$$$     .        vrankM(j) == mpi__rankM ) then
c$$$            ranklistB(i) = j
c$$$            i=i+1
c$$$         end if
c$$$      end do
c$$$      call MPI_Group_incl( mpi__group,
c$$$     .     mpi__sizeB, ranklistB, mpi__groupB, mpi__info )
c$$$      call MPI_Comm_create( MPI_COMM_WORLD,
c$$$     .     mpi__groupB, mpi__commB, mpi__info )
c$$$      deallocate( ranklistB )
c$$$
c$$$C     ! setup M-parallel      
c$$$C     !  commM for among other M, same Q,P,B
c$$$      allocate( ranklistM(0:mpi__sizeM-1) )
c$$$      i=0
c$$$      do j=0, mpi__size-1
c$$$         if(  vrankQ(j) == mpi__rankQ .and.
c$$$     .        vrankS(j) == mpi__rankS .and.
c$$$     .        vrankP(j) == mpi__rankP .and.
c$$$     .        vrankB(j) == mpi__rankB ) then
c$$$            ranklistM(i) = j
c$$$            i=i+1
c$$$         end if
c$$$      end do
c$$$      call MPI_Group_incl( mpi__group,
c$$$     .     mpi__sizeM, ranklistM, mpi__groupM, mpi__info )
c$$$      call MPI_Comm_create( MPI_COMM_WORLD,
c$$$     .     mpi__groupM, mpi__commM, mpi__info )
c$$$      deallocate( ranklistM )

c$$$#if SCALAPACK
c$$$C     ! setup ScaLapack process list
c$$$      allocate( imap(0:mpi__sizeM-1) )
c$$$      
c$$$      do i=0, mpi__sizeM-1
c$$$         imap(i) = mpi__rankQ * (mpi__sizeS*mpi__sizeP*mpi__sizeB*mpi__sizeM)
c$$$     .        +    mpi__rankS * (mpi__sizeP*mpi__sizeB*mpi__sizeM)
c$$$     .        +    mpi__rankP * (mpi__sizeB*mpi__sizeM)
c$$$     .        +    mpi__rankB * (mpi__sizeM) + i
c$$$      end do
c$$$
c$$$      call BLACS_GET( 0, 0, mpi__Mhandle )
c$$$      call BLACS_GRIDMAP( mpi__Mhandle, imap, 1, 1, mpi__sizeM )
c$$$      deallocate(imap)
c$$$#endif
      return
      end subroutine MPI__InitializeQ!SPBM



C     !======================================================
      subroutine MPI__consoleout(idn)
      implicit none
      character(1024*4) :: cwd, stdout
      character*(*):: idn
! console-output from different nodes to different files
      if( mpi__size > 1 ) then
        if( mpi__root ) then
           write(6,"(' MPI outputs in each rank are in stdout.{RankId}.',a)")idn
           call flush(6)
        end if
c        write(6,"('   stdout.',i4.4,'.',a)") mpi__rank,idn
        write(stdout,"('stdout.',i4.4,'.',a)") mpi__rank,idn
        open(unit=6,file=trim(stdout))
        write(6,"(a,i3)")" ### console output for rank=",mpi__rank
      endif
      return
      end subroutine MPI__consoleout

C     !======================================================
      subroutine MPI__consoleout_magnon(idn,size_lim)
      implicit none
      character(1024*4) :: cwd, stdout
      character*(*):: idn
      integer , intent(in) :: size_lim

CC Reduce size for magnon (avoid memory leak)      
      if (mpi__sizeMG > size_lim) then
         mpi__sizeMG=size_lim
      endif

! console-output from different nodes to different files
      if( mpi__sizeMG > 1 ) then
        if( mpi__root ) then
           write(6,"(' MPI outputs in each rank are in stdout.{RankId}.',a)")idn
           call flush(6)
        end if
c        write(6,"('   stdout.',i4.4,'.',a)") mpi__rank,idn
        write(stdout,"('stdout.',i4.4,'.',a)") mpi__rankMG,idn
        open(unit=6,file=trim(stdout))
        write(6,"(a,i3)")" ### console output for rank=",mpi__rankMG
      endif
      return
      end subroutine MPI__consoleout_magnon


C     !======================================================
      subroutine MPI__Barrier
      implicit none

      call MPI_Barrier( MPI_COMM_WORLD, mpi__info )

      end subroutine MPI__Barrier



C     !======================================================
      subroutine MPI__Finalize
      implicit none
c#if SCALAPACK
c      if( mpi__Mhandle /= 0 ) then
c         call BLACS_GRIDEXIT( mpi__Mhandle )
c      end if
c#endif
      call MPI_Finalize ( mpi__info )

      end subroutine MPI__Finalize



C     !======================================================
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


C     !======================================================
      subroutine MPI__Broadcast( data )
      implicit none
      integer, intent(inout) :: data

      call MPI_Bcast( data, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi__info )

      return
      end subroutine MPI__Broadcast


C     !======================================================
      subroutine MPI__send_d(data,dest)
      implicit none
      real(8):: data
      integer :: n,dest,ierr
      n=1
      call MPI_Send(data,n,MPI_REAL8,dest,mpi__rank, MPI_COMM_WORLD,ierr)
      end subroutine MPI__send_d

C     !======================================================
      subroutine MPI__recv_d(data,src)
      implicit none
      real(8):: data
      integer :: n,src,ierr
      n=1
      call MPI_Recv(data,n,MPI_REAL8,src,src, MPI_COMM_WORLD,ista,ierr)
      end subroutine MPI__recv_d

C     !======================================================
      subroutine MPI__send_dv(data,dest)
      implicit none
      real(8):: data(:)
      integer :: n,dest,ierr
      n=size(data)
      call MPI_Send(data,n,MPI_REAL8,dest,mpi__rank, MPI_COMM_WORLD,ierr)
      end subroutine MPI__send_dv

C     !======================================================
      subroutine MPI__recv_dv(data,src)
      implicit none
      real(8):: data(:)
      integer :: n,src,ierr
      n=size(data)
      call MPI_Recv(data,n,MPI_REAL8,src,src, MPI_COMM_WORLD,ista,ierr)
      end subroutine MPI__recv_dv

C     !======================================================
      subroutine MPI__send_i(data,dest)
      implicit none
      integer:: data
      integer :: n,dest,ierr
      n=1
      call MPI_Send(data,n,MPI_INTEGER,dest,mpi__rank, MPI_COMM_WORLD,ierr)
      end subroutine MPI__send_i

C     !======================================================
      subroutine MPI__recv_i(data,src)
      implicit none
      integer:: data
      integer :: n,src,ierr
      n=1
      call MPI_Recv(data,n,MPI_INTEGER,src,src, MPI_COMM_WORLD,ista,ierr)
      end subroutine MPI__recv_i

C     !======================================================
      subroutine MPI__send_iv(data,dest)
      implicit none
      integer:: data(:)
      integer :: n,dest,ierr
      n=size(data)
      call MPI_Send(data,n,MPI_INTEGER,dest,mpi__rank, MPI_COMM_WORLD,ierr)
      end subroutine MPI__send_iv

C     !======================================================
      subroutine MPI__recv_iv(data,src)
      implicit none
      integer:: data(:)
      integer :: n,src,ierr
      n=size(data)
      call MPI_Recv(data,n,MPI_INTEGER,src,src, MPI_COMM_WORLD,ista,ierr)
      end subroutine MPI__recv_iv

C     !======================================================
      subroutine MPI__REAL8send(data,n,dest)
      implicit none
      real(8):: data(n)
      integer :: n,dest,ierr
      call MPI_Send(data,n,MPI_REAL8,dest,mpi__rank, MPI_COMM_WORLD,ierr)
      end subroutine MPI__REAL8send

C     !======================================================
      subroutine MPI__REAL8recv(data,n,src)
      implicit none
      real(8):: data(n)
      integer :: n,src,ierr
      call MPI_Recv(data,n,MPI_REAL8,src,src, MPI_COMM_WORLD,ista,ierr)
      end subroutine MPI__REAL8recv

C     !======================================================
      subroutine MPI__DbleCOMPLEXsend(data,n,dest)
      implicit none
      complex(8):: data(n)
      integer :: n,dest,ierr
      call MPI_Send(data,n,MPI_COMPLEX16,dest,mpi__rank, MPI_COMM_WORLD,ierr)
      end subroutine MPI__DbleCOMPLEXsend

C     !======================================================
      subroutine MPI__DbleCOMPLEXrecv(data,n,src)
      implicit none
      complex(8):: data(n)
      integer :: n,src,ierr
      call MPI_Recv(data,n,MPI_COMPLEX16,src,src, MPI_COMM_WORLD,ista,ierr)
      end subroutine MPI__DbleCOMPLEXrecv

C     !======================================================
      subroutine MPI__DbleCOMPLEXsendQ(data,n,destQ)
      implicit none
      complex(8):: data(n)
      integer :: n,destQ,ierr
      call MPI_Send(data,n,MPI_COMPLEX16,destQ,0,mpi__commQ,ierr)
      end subroutine MPI__DbleCOMPLEXsendQ

C     !======================================================
      subroutine MPI__DbleCOMPLEXrecvQ(data,n,srcQ)
      implicit none
      complex(8):: data(n)
      integer :: n,srcQ,ierr
      call MPI_Recv(data,n,MPI_COMPLEX16,srcQ,0,mpi__commQ,ista,ierr)
      end subroutine MPI__DbleCOMPLEXrecvQ

C     !======================================================
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

C     !======================================================
      subroutine MPI__reduceSum( root, data, sizex )
      implicit none
      integer, intent(in) :: sizex,root
      complex(8), intent(inout) :: data(sizex)
      complex(8), allocatable   :: mpi__data(:) 

      if( mpi__size == 1 ) return

      allocate(mpi__data(sizex))
      mpi__data = data

      call MPI_reduce( mpi__data, data, sizex,
     &     MPI_DOUBLE_COMPLEX, MPI_SUM, root, MPI_COMM_WORLD, mpi__info )

      deallocate( mpi__data )

      return
      end subroutine MPI__reduceSum



C     !======================================================
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


c$$$C     !======================================================
c$$$      subroutine MPI__AllreduceSumS( data, sizex )
c$$$      implicit none
c$$$      integer, intent(in) :: sizex
c$$$      complex(8), intent(inout) :: data(sizex)
c$$$      complex(8), allocatable   :: mpi__data(:)
c$$$      if( mpi__sizeS == 1 ) return
c$$$      allocate(mpi__data(sizex))
c$$$      mpi__data = data
c$$$      call MPI_Allreduce( mpi__data, data, sizex,
c$$$     &     MPI_DOUBLE_COMPLEX, MPI_SUM, mpi__commS, mpi__info )
c$$$      deallocate( mpi__data )
c$$$      return
c$$$      end subroutine MPI__AllreduceSumS


C     !======================================================
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

c$$$      write(6,*) "irkip_all", mpi__rank
      total = 0
      do ispinmx=1, nspinmx
         do iq=1, nq
            do iqibz=1, nqibz
               do igrp=1, ngrp
                  if( irkip_all(ispinmx,iqibz,igrp,iq) /= 0 ) then
                     total = total + 1
c$$$                     write(6,*) ispinmx, iq, iqibz, igrp, irkip_all(ispinmx,iqibz,igrp,iq)
                  end if
               end do
            end do
         end do
      end do
      write(6,"('MPI__sxcf_rankdivider:$')")
      write(6,"('nspinmx,nqibz,ngrp,nq,total=',5i6)")
     &     nspinmx,nqibz,ngrp,nq,total 
      
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

c$$$      write(6,*) "irkip", mpi__rank
      do ispinmx=1, nspinmx
         do iq=1, nq
            do iqibz=1, nqibz
               do igrp=1, ngrp
                  if( irkip_all(ispinmx,iqibz,igrp,iq) /= 0 ) then
                     total = total + 1
                     if( indexi<=total .and. total<=indexe ) then
                        irkip(ispinmx,iqibz,igrp,iq) = irkip_all(ispinmx,iqibz,igrp,iq)
c$$$                        write(6,*) ispinmx, iq, iqibz, igrp, irkip(ispinmx,iqibz,igrp,iq)
                     end if
                  end if
               end do
            end do
         end do
      end do

      return
      end subroutine MPI__sxcf_rankdivider

!!
C     !======================================================
      subroutine MPI__hx0fp0_rankdivider(iqxini,iqxend,nqibz)
      implicit none
      integer, intent(in) :: iqxini, iqxend, nqibz

      integer :: iq
      allocate( mpi__task(1:iqxend),mpi__ranktab(1:iqxend) )
      
      if( mpi__size == 1 ) then
         mpi__task(:) = .true.
         mpi__ranktab(:) = mpi__rank
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


C     !======================================================
      subroutine MPI__hx0fp0_rankdivider2(iqxini,iqxend)
      implicit none
      integer, intent(in) :: iqxini, iqxend
      integer :: iq,i
      allocate( mpi__task(1:iqxend),mpi__ranktab(1:iqxend) )
      mpi__task(:) = .false.
      mpi__ranktab(1:iqxend)=999999
      if( mpi__size == 1 ) then
         mpi__task(:) = .true.
         mpi__ranktab(:) = mpi__rank
         write(6,*) "rankdivider"
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

C     !======================================================
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

C     !======================================================
      subroutine MPI__hx0fp0_rankdivider2Q(iqxini,iqxend)
!! not mpi_commq is used in m_llw      
      implicit none
      integer, intent(in) :: iqxini, iqxend
      integer :: iq,i
      integer,allocatable:: ranklistQ(:)
      integer :: mpi__group
      integer :: mpi__groupQ
      
      mpi__sizeQ = mpi__size
      mpi__rankQ = mod(mpi__rank,mpi__sizeQ) !vrankQ(mpi__rank)
      mpi__rootQ = ( mpi__rankQ == 0 )
      mpi__commQ = MPI_COMM_WORLD !! mpi__commQ is used in m_llw
c      
C     ! setup Q-parallel
C     !  commQ for among other Q, same S,P,B,M
c      allocate( ranklistQ(0:mpi__sizeQ-1) )
c      i=0
c      do j=0, mpi__size-1
cc     if(  vrankS(j) == mpi__rankS .and.
cc     .        vrankP(j) == mpi__rankP .and.
cc     .        vrankB(j) == mpi__rankB .and.
cc     .        vrankM(j) == mpi__rankM ) then
c         ranklistQ(i) = j
c         i=i+1
cc     end if
c      end do
c      call MPI_Comm_group( MPI_COMM_WORLD, mpi__group, mpi__info )
c      call MPI_Group_incl(  mpi__group, mpi__sizeQ, ranklistQ, mpi__groupQ, mpi__info )
c      call MPI_Comm_create( MPI_COMM_WORLD, mpi__groupQ, mpi__commQ, mpi__info ) 
c      deallocate( ranklistQ )
c      
      allocate( mpi__Qtask(1:iqxend), mpi__Qranktab(1:iqxend) )
      mpi__Qtask(:) = .false.
      mpi__Qranktab(1:iqxend)=999999
      if( mpi__sizeQ == 1 ) then
         mpi__Qtask(:) = .true.
         mpi__Qranktab(:) = mpi__rankQ
         return
      end if
      if(mpi__root) then
         write(6,*) "MPI_hx0fp0_rankdivider2Q:"
         write(6,'(a,$)')'mpi__Qtask='
         write(6,'(10L2)')mpi__Qtask
      endif
      do iq=iqxini, iqxend
         mpi__Qranktab(iq) = mod(iq-1,mpi__sizeQ)  ! rank_table for given iq. iq=1 must give rank=0
         if( mpi__Qranktab(iq) == mpi__rankQ) then
            mpi__Qtask(iq) = .true.                ! mpi__task is nodeID-dependent.
         endif
         if(mpi__root) write(6,"('  iq irank=',2i5)")iq,mpi__Qranktab(iq)
      enddo   
      end subroutine MPI__hx0fp0_rankdivider2Q
      
      end module m_mpi
      
      
c$$$C     !======================================================
c$$$      subroutine MPI__hx0fp0_rankdivider2S( nspinmx )
c$$$      implicit none
c$$$      integer, intent(in) :: nspinmx
c$$$      integer :: k, jpm, ibib
c$$$      integer :: p, n
c$$$      integer :: blocksize
c$$$      integer :: max_nbnb
c$$$C     ! for S(npm) parallelization
c$$$      mpi__Snall = nspinmx
c$$$      if( allocated(mpi__Svn) ) deallocate(mpi__Svn)
c$$$      if( allocated(mpi__Svs) ) deallocate(mpi__Svs)
c$$$      if( allocated(mpi__Sve) ) deallocate(mpi__Sve)
c$$$      allocate( mpi__Svn(0:mpi__sizeS-1) )
c$$$      allocate( mpi__Svs(0:mpi__sizeS-1) )
c$$$      allocate( mpi__Sve(0:mpi__sizeS-1) )
c$$$
c$$$      mpi__Svn(:) = mpi__Snall/mpi__sizeS
c$$$      if( mpi__Svn(0) == 0 ) then
c$$$         write(6,*) "Warning: too many processes for S parallelization"
c$$$         write(6,*) "  maximum -np ", mpi__Snall
c$$$         call rx('Warning: too many processes for S parallelization')
c$$$      end if
c$$$
c$$$      do p=0, mpi__sizeS-1
c$$$         if( sum(mpi__Svn(:)) == mpi__Snall ) exit
c$$$         mpi__Svn(p) = mpi__Svn(p) + 1
c$$$      end do
c$$$
c$$$      n=1
c$$$      do p=0, mpi__sizeS-1
c$$$         mpi__Svs(p) = n
c$$$         mpi__Sve(p) = mpi__Svs(p) + mpi__Svn(p) - 1
c$$$         n = n + mpi__Svn(p)
c$$$      end do
c$$$      mpi__Sn = mpi__Svn(mpi__rankS)
c$$$      mpi__Ss = mpi__Svs(mpi__rankS)
c$$$      mpi__Se = mpi__Sve(mpi__rankS)
c$$$
c$$$      if (mpi__root) then
c$$$      write(6,*)mpi__rank,'mpi__Svn',mpi__Svn
c$$$      write(6,*)mpi__rank,'mpi__Svs',mpi__Svs
c$$$      write(6,*)mpi__rank,'mpi__Sve',mpi__Sve
c$$$      endif
c$$$#if 0
c$$$      write(6,*) "DEBUG: mpi__Qnall=", size(mpi__Qtask)
c$$$      write(6,*) "DEBUG: mpi__Snall=", mpi__Snall
c$$$
c$$$      write(6,*) "DEBUG: mpi__Qn=", count(mpi__Qtask)
c$$$      write(6,*) "DEBUG: mpi__Sn=", mpi__Sn
c$$$#endif
c$$$      return
c$$$      end subroutine MPI__hx0fp0_rankdivider2S


c$$$C     !======================================================
c$$$      subroutine MPI__x0kf_rankdivider( nbnb, nqbz, npm, ngb, nwt )
c$$$      implicit none
c$$$
c$$$      integer, intent(in) :: nbnb(nqbz,npm), nqbz,npm
c$$$      integer, intent(in) :: ngb,nwt
c$$$      
c$$$      integer :: k, jpm, ibib
c$$$      integer :: p, n
c$$$      integer :: blocksize
c$$$      integer :: max_nbnb
c$$$
c$$$C     ! for P(npm) parallelization
c$$$      mpi__Pnall = npm
c$$$
c$$$      if( allocated(mpi__Pvn) ) deallocate(mpi__Pvn)
c$$$      if( allocated(mpi__Pvs) ) deallocate(mpi__Pvs)
c$$$      if( allocated(mpi__Pve) ) deallocate(mpi__Pve)
c$$$      allocate( mpi__Pvn(0:mpi__sizeP-1) )
c$$$      allocate( mpi__Pvs(0:mpi__sizeP-1) )
c$$$      allocate( mpi__Pve(0:mpi__sizeP-1) )
c$$$
c$$$      mpi__Pvn(:) = mpi__Pnall/mpi__sizeP
c$$$      if( mpi__Pvn(0) == 0 ) then
c$$$         write(6,*) "Warning: too many processes for P parallelization"
c$$$         write(6,*) "  maximum -np ", mpi__Pnall
c$$$         call rx("Warning: too many processes for P parallelization")
c$$$      end if
c$$$
c$$$      do p=0, mpi__sizeP-1
c$$$         if( sum(mpi__Pvn(:)) == mpi__Pnall ) exit
c$$$         mpi__Pvn(p) = mpi__Pvn(p) + 1
c$$$      end do
c$$$
c$$$      n=1
c$$$      do p=0, mpi__sizeP-1
c$$$         mpi__Pvs(p) = n
c$$$         mpi__Pve(p) = mpi__Pvs(p) + mpi__Pvn(p) - 1
c$$$         n = n + mpi__Pvn(p)
c$$$      end do
c$$$      mpi__Pn = mpi__Pvn(mpi__rankP)
c$$$      mpi__Ps = mpi__Pvs(mpi__rankP)
c$$$      mpi__Pe = mpi__Pve(mpi__rankP)
c$$$
c$$$
c$$$C     ! for B(Kpoint&nbnb) parallelization
c$$$      mpi__Bnall=sum(nbnb(1:nqbz,mpi__Ps:mpi__Pe))
c$$$
c$$$      
c$$$      if( allocated(mpi__Bvn) ) deallocate(mpi__Bvn)
c$$$      if( allocated(mpi__Bvs) ) deallocate(mpi__Bvs)
c$$$      if( allocated(mpi__Bve) ) deallocate(mpi__Bve)
c$$$      allocate( mpi__Bvn(0:mpi__sizeB-1) )
c$$$      allocate( mpi__Bvs(0:mpi__sizeB-1) )
c$$$      allocate( mpi__Bve(0:mpi__sizeB-1) )
c$$$
c$$$      mpi__Bvn(:) = mpi__Bnall/mpi__sizeB
c$$$      if( mpi__Bvn(0) == 0 ) then
c$$$         write(6,*) "Warning: too many processes for B parallelization"
c$$$         write(6,*) "  maximum -nb ", mpi__Bnall
c$$$         call rx("Warning: too many processes for B parallelization")
c$$$      end if
c$$$      do p=0, mpi__sizeB-1
c$$$         if( sum(mpi__Bvn(:)) == mpi__Bnall ) exit
c$$$         mpi__Bvn(p) = mpi__Bvn(p) + 1
c$$$      end do
c$$$
c$$$      n=1
c$$$      do p=0, mpi__sizeB-1
c$$$         mpi__Bvs(p) = n
c$$$         mpi__Bve(p) = mpi__Bvs(p) + mpi__Bvn(p) - 1
c$$$         n = n + mpi__Bvn(p)
c$$$      end do
c$$$      mpi__Bn = mpi__Bvn(mpi__rankB)
c$$$      mpi__Bs = mpi__Bvs(mpi__rankB)
c$$$      mpi__Be = mpi__Bve(mpi__rankB)
c$$$
c$$$      max_nbnb=maxval(nbnb(1:nqbz,mpi__Ps:mpi__Pe))      
c$$$
c$$$      if( allocated(mpi__Btask1) ) deallocate(mpi__Btask1)
c$$$      if( allocated(mpi__Btask2) ) deallocate(mpi__Btask2)
c$$$      if( allocated(mpi__Btask3) ) deallocate(mpi__Btask3)
c$$$      allocate( mpi__Btask1(nqbz) )
c$$$      allocate( mpi__Btask2(nqbz,mpi__Ps:mpi__Pe) )
c$$$      allocate( mpi__Btask3(nqbz,mpi__Ps:mpi__Pe,max_nbnb) )
c$$$
c$$$      mpi__Btask1(:)     = .false.
c$$$      mpi__Btask2(:,:)   = .false.
c$$$      mpi__Btask3(:,:,:) = .false.
c$$$
c$$$      n=1
c$$$      do k=1, nqbz
c$$$         do jpm=mpi__Ps, mpi__Pe
c$$$            do ibib=1, nbnb(k,jpm)
c$$$               if( mpi__Bs <= n .and. n<=mpi__Be ) then
c$$$                  mpi__Btask1(k) = .true.
c$$$                  mpi__Btask2(k,jpm) = .true.
c$$$                  mpi__Btask3(k,jpm,ibib) = .true.
c$$$               end if
c$$$               n=n+1
c$$$            end do
c$$$         end do
c$$$      end do
c$$$
c$$$
c$$$C     ! for M(matrix) parallelization
c$$$      if( allocated(mpi__Mvn) ) then
c$$$         deallocate( mpi__Mvn )
c$$$         deallocate( mpi__Mvs )
c$$$         deallocate( mpi__Mve )
c$$$      end if
c$$$      allocate( mpi__Mvn(0:mpi__sizeM-1) )
c$$$      allocate( mpi__Mvs(0:mpi__sizeM-1) )
c$$$      allocate( mpi__Mve(0:mpi__sizeM-1) )
c$$$    
c$$$      mpi__Mnall = ngb
c$$$      
c$$$C     !! the mean is larger
c$$$      mpi__Mvn(:) = int(ceiling(dble(mpi__Mnall)/dble(mpi__sizeM)))
c$$$C     !! the last block is smaller 
c$$$      mpi__Mvn(mpi__sizeM-1) = mpi__Mvn(mpi__sizeM-1)
c$$$     .     - (sum(mpi__Mvn(:))-mpi__Mnall)
c$$$      blocksize = maxval(mpi__Mvn(:))
c$$$c#if SCALAPACK           
c$$$c      call DescInit( mpi__Mdescv, 1, mpi__Mnall,
c$$$c     .     1, blocksize, 0, 0,
c$$$c     .     mpi__Mhandle, 1, mpi__info )
c$$$c      
c$$$c      call DescInit( mpi__Mdescm, mpi__Mnall, mpi__Mnall,
c$$$c     .     blocksize, blocksize, 0, 0,
c$$$c     .     mpi__Mhandle, mpi__Mnall, mpi__info )
c$$$c#endif
c$$$      n=1
c$$$      do p=0, mpi__sizeM-1
c$$$         mpi__Mvs(p) = n
c$$$         mpi__Mve(p) = mpi__Mvs(p) + mpi__Mvn(p) - 1
c$$$         n = n + mpi__Mvn(p)
c$$$      end do
c$$$      mpi__Ms = mpi__Mvs(mpi__rankM)
c$$$      mpi__Me = mpi__Mve(mpi__rankM)
c$$$      mpi__Mn = mpi__Mvn(mpi__rankM)
c$$$      mpi__Mf = mpi__Ms+blocksize-1
c$$$
c$$$C     ! for W(iw) parallelization
c$$$      mpi__Wnall = nwt
c$$$
c$$$      if( allocated(mpi__Wvn) ) deallocate(mpi__Wvn)
c$$$      if( allocated(mpi__Wvs) ) deallocate(mpi__Wvs)
c$$$      if( allocated(mpi__Wve) ) deallocate(mpi__Wve)
c$$$      allocate( mpi__Wvn(0:mpi__sizeW-1) )
c$$$      allocate( mpi__Wvs(0:mpi__sizeW-1) )
c$$$      allocate( mpi__Wve(0:mpi__sizeW-1) )
c$$$
c$$$      mpi__Wvn(:) = mpi__Wnall/mpi__sizeB
c$$$      if( mpi__Wvn(0) == 0 ) then
c$$$         write(6,*) "Warning: too many processes for B parallelization"
c$$$         write(6,*) "  maximum -np ", mpi__Wnall
c$$$         call rx("Warning: too many processes for B parallelization")
c$$$      end if
c$$$
c$$$      do p=0, mpi__sizeW-1
c$$$         if( sum(mpi__Wvn(:)) == mpi__Wnall ) exit
c$$$         mpi__Wvn(p) = mpi__Wvn(p) + 1
c$$$      end do
c$$$
c$$$      n=1
c$$$      do p=0, mpi__sizeW-1
c$$$         mpi__Wvs(p) = n
c$$$         mpi__Wve(p) = mpi__Wvs(p) + mpi__Wvn(p) - 1
c$$$         n = n + mpi__Wvn(p)
c$$$      end do
c$$$      mpi__Wn = mpi__Wvn(mpi__rankW)
c$$$      mpi__Ws = mpi__Wvs(mpi__rankW)
c$$$      mpi__We = mpi__Wve(mpi__rankW)
c$$$
c$$$      write(6,*) "DEBUG: mpi__Qnall=", size(mpi__Qtask)
c$$$      write(6,*) "DEBUG: mpi__Pnall=", mpi__Pnall
c$$$      write(6,*) "DEBUG: mpi__Bnall=", mpi__Bnall
c$$$      write(6,*) "DEBUG: mpi__Mnall=", mpi__Mnall
c$$$      write(6,*) "DEBUG: mpi__Wnall=", mpi__Wnall
c$$$
c$$$      write(6,*) "DEBUG: mpi__Qn=", count(mpi__Qtask)
c$$$      write(6,*) "DEBUG: mpi__Pn=", mpi__Pn
c$$$      write(6,*) "DEBUG: mpi__Bn=", mpi__Bn
c$$$      write(6,*) "DEBUG: mpi__Mn=", mpi__Mn
c$$$      write(6,*) "DEBUG: mpi__Wn=", mpi__Wn
c$$$
c$$$      return
c$$$      end subroutine MPI__x0kf_rankdivider

c$$$#if SCALAPACK
c$$$C     !======================================================
c$$$      subroutine MPI__dummy
c$$$      implicit none
c$$$
c$$$      call PZHEGVX              ! dummy call, to poll scalapack library
c$$$
c$$$      return
c$$$      end subroutine MPI__dummy
c$$$      
c$$$C     !======================================================
c$$$C     !!  y = A * x
c$$$      subroutine MPI__ZGEMV( matrixA, vectorX, vectorY )
c$$$      implicit none
c$$$      
c$$$      complex(8), intent(in)  :: matrixA(mpi__Mnall,mpi__Ms:mpi__Mf)
c$$$      complex(8), intent(in)  :: vectorX(mpi__Ms:mpi__Mf)
c$$$      complex(8), intent(out) :: vectorY(mpi__Ms:mpi__Mf)
c$$$      character,  parameter   :: transa = 'N'
c$$$      complex(8), parameter   :: alpha = (1.0d0,0.0d0)
c$$$      complex(8), parameter   :: beta  = (0.0d0,0.0d0)
c$$$      integer :: i, j
c$$$
c$$$      call PZGEMV( transa, mpi__Mnall, mpi__Mnall, alpha,
c$$$     .     matrixA, 1, 1, mpi__Mdescm,
c$$$     .     vectorX, 1, 1, mpi__Mdescv, 1, beta,
c$$$     .     vectorY, 1, 1, mpi__Mdescv, 1 )
c$$$
c$$$      return
c$$$      end subroutine MPI__ZGEMV
c$$$      
c$$$C     !======================================================
c$$$C     !! C = A*B
c$$$      subroutine MPI__ZGEMM( matrixA, matrixB, matrixC )
c$$$      implicit none
c$$$
c$$$      complex(8), intent(in)  :: matrixA(mpi__Mnall,mpi__Ms:mpi__Mf)
c$$$      complex(8), intent(in)  :: matrixB(mpi__Mnall,mpi__Ms:mpi__Mf)
c$$$      complex(8), intent(out) :: matrixC(mpi__Mnall,mpi__Ms:mpi__Mf)
c$$$      character,  parameter   :: transa = 'N'
c$$$      character,  parameter   :: transb = 'N'
c$$$      complex(8), parameter   :: alpha = (1.0d0,0.0d0)
c$$$      complex(8), parameter   :: beta  = (0.0d0,0.0d0)
c$$$      integer :: i, j, k
c$$$
c$$$      call PZGEMM( transa, transb,
c$$$     .     mpi__Mnall, mpi__Mnall, mpi__Mnall,
c$$$     .     alpha,
c$$$     .     matrixA, 1, 1, mpi__Mdescm,
c$$$     .     matrixB, 1, 1, mpi__Mdescm,
c$$$     .     beta,
c$$$     .     matrixC, 1, 1, mpi__Mdescm )
c$$$
c$$$      return
c$$$      end subroutine MPI__ZGEMM
c$$$
c$$$C     !======================================================
c$$$C     !!  C = x^H * y
c$$$      subroutine MPI__ZGERC( vectorX, vectorY, matrixC )
c$$$      implicit none
c$$$
c$$$      complex(8), intent(in)  :: vectorX(mpi__Ms:mpi__Mf)
c$$$      complex(8), intent(in)  :: vectorY(mpi__Ms:mpi__Mf)
c$$$      complex(8), intent(out) :: matrixC(mpi__Mnall,mpi__Ms:mpi__Mf)
c$$$      character,  parameter   :: transa = 'C', transb = 'N'
c$$$      complex(8), parameter   :: alpha = (1.0d0,0.0d0)
c$$$      complex(8), parameter   :: beta  = (0.0d0,0.0d0)
c$$$      integer :: i, j
c$$$
c$$$      call PZGEMM( transa, transb, mpi__Mnall, mpi__Mnall, 1, alpha,
c$$$     .     vectorX, 1, 1, mpi__Mdescv,
c$$$     .     vectorY, 1, 1, mpi__Mdescv, beta,
c$$$     .     matrixC, 1, 1, mpi__Mdescm )
c$$$
c$$$      return
c$$$      end subroutine MPI__ZGERC
c$$$
c$$$C     !======================================================
c$$$C     !!  C = a*A + C
c$$$      subroutine MPI__ZAMPM( alpha, matrixA, matrixC )
c$$$      implicit none
c$$$
c$$$      real(8), intent(in)       :: alpha
c$$$      complex(8), intent(in)    :: matrixA(mpi__Mnall,mpi__Ms:mpi__Mf)
c$$$      complex(8), intent(inout) :: matrixC(mpi__Mnall,mpi__Ms:mpi__Mf)
c$$$      integer :: i, j
c$$$
c$$$!$OMP PARALLEL DO PRIVATE(j) FIRSTPRIVATE(mpi__Ms,mpi__Mf)
c$$$      do j=mpi__Ms,mpi__Mf
c$$$         matrixC(:,j) = alpha*matrixA(:,j) + matrixC(:,j)
c$$$      end do
c$$$!$OMP END PARALLEL
c$$$
c$$$      return
c$$$      end subroutine MPI__ZAMPM
c$$$
c$$$C     !======================================================
c$$$C     !! B = reduce A
c$$$      subroutine MPI__AllreduceMatrixB( matrixA )
c$$$      implicit none
c$$$      complex(8), intent(inout) :: matrixA(mpi__Mnall,mpi__Ms:mpi__Me)
c$$$      complex(8), allocatable   :: matrixA_sum(:,:)
c$$$
c$$$      if( mpi__sizeB == 1 ) return
c$$$
c$$$      allocate( matrixA_sum(mpi__Mnall,mpi__Ms:mpi__Me) )
c$$$
c$$$      call MPI_Allreduce( matrixA, matrixA_sum,
c$$$     .     mpi__Mn*mpi__Mnall,
c$$$     .     MPI_DOUBLE_COMPLEX, MPI_SUM, mpi__commB, mpi__info )
c$$$      
c$$$      matrixA = matrixA_sum
c$$$      deallocate( matrixA_sum )
c$$$
c$$$      return
c$$$      end subroutine MPI__AllreduceMatrixB
c$$$
c$$$C     !======================================================
c$$$C     !! y = gather x
c$$$      subroutine MPI__AllgatherVectorM( vectorX, vectorY )
c$$$      implicit none
c$$$      complex(8), intent(in)  :: vectorX(mpi__Ms:mpi__Me)
c$$$      complex(8), intent(out) :: vectorY(mpi__Mnall)
c$$$
c$$$      if( mpi__sizeM == 1 ) then
c$$$         vectorY = vectorX
c$$$         return
c$$$      end if
c$$$
c$$$      call MPI_Allgatherv( vectorX(mpi__Ms), mpi__Mn,
c$$$     .     MPI_DOUBLE_COMPLEX,
c$$$     .     vectorY, (mpi__Mvn(:)+0), (mpi__Mvs(:)-1),
c$$$     .     MPI_DOUBLE_COMPLEX, mpi__commM, mpi__info )
c$$$
c$$$      return
c$$$      end subroutine MPI__AllgatherVectorM
c$$$
c$$$C     !======================================================
c$$$C     !! Aall = gather A
c$$$      subroutine MPI__AllgatherMatrixM( matrixA, matrixAall )
c$$$      implicit none
c$$$      complex(8), intent(in)   :: matrixA   (mpi__Mnall,mpi__Ms:mpi__Me)
c$$$      complex(8), intent(out)  :: matrixAall(mpi__Mnall,mpi__Mnall)
c$$$
c$$$      if( mpi__sizeM == 1 ) then
c$$$         matrixAall = matrixA
c$$$         return
c$$$      end if
c$$$
c$$$      call MPI_Allgatherv( matrixA(1,mpi__Ms), (mpi__Mn)*mpi__Mnall,
c$$$     .     MPI_DOUBLE_COMPLEX, matrixAall,
c$$$     .     (mpi__Mvn(:)+0)*mpi__Mnall, (mpi__Mvs(:)-1)*mpi__Mnall,
c$$$     .     MPI_DOUBLE_COMPLEX, mpi__commM, mpi__info )
c$$$
c$$$      return
c$$$      end subroutine MPI__AllgatherMatrixM
c$$$
c$$$
c$$$C     !======================================================
c$$$C     !! Aall = gather A
c$$$      subroutine MPI__AllgatherArrayP( arrayA, arrayAall )
c$$$      implicit none
c$$$      complex(8), intent(in)   :: arrayA   
c$$$     .     (mpi__Mnall,mpi__Mnall,mpi__Wnall,mpi__Ps:mpi__Pe)
c$$$      complex(8), intent(out)  :: arrayAall
c$$$     .     (mpi__Mnall,mpi__Mnall,mpi__Wnall,mpi__Pnall)
c$$$
c$$$      if( mpi__sizeP == 1 ) then
c$$$         arrayAall = arrayA
c$$$         return
c$$$      end if
c$$$
c$$$      call MPI_Allgatherv( arrayA(1,1,1,mpi__Ps),
c$$$     .     mpi__Mnall*mpi__Mnall*mpi__Wnall*(mpi__Pn),
c$$$     .     MPI_DOUBLE_COMPLEX, arrayAall,
c$$$     .     mpi__Mnall*mpi__Mnall*mpi__Wnall*(mpi__Pvn(:)+0),
c$$$     .     mpi__Mnall*mpi__Mnall*mpi__Wnall*(mpi__Pvs(:)-1),
c$$$     .     MPI_DOUBLE_COMPLEX, mpi__commP, mpi__info )
c$$$
c$$$      return
c$$$      end subroutine MPI__AllgatherArrayP
c$$$
c$$$
c$$$C     !======================================================
c$$$C     !! matrixB(:,s:e) = matrixA(s:e,:)
c$$$      subroutine MPI__AllMPtoAllMW( matrixA, matrixB )
c$$$      implicit none
c$$$      complex(8), intent(in)  :: matrixA
c$$$     .     (mpi__Mnall,mpi__Ms:mpi__Me,mpi__Wnall,mpi__Ps:mpi__Pe)
c$$$      complex(8), intent(out) :: matrixB
c$$$     .     (mpi__Mnall,mpi__Ms:mpi__Me,mpi__Ws:mpi__We,mpi__Pnall)
c$$$
c$$$      complex(8), allocatable :: matrix_send(:,:,:)
c$$$      complex(8), allocatable :: matrix_recv(:,:,:)
c$$$      integer, allocatable :: vcount_send(:), vdispl_send(:)
c$$$      integer, allocatable :: vcount_recv(:), vdispl_recv(:)
c$$$      integer :: i, j, p, index
c$$$
c$$$      if( mpi__sizeP == 1 ) then
c$$$         matrixB = matrixA
c$$$         return
c$$$      end if
c$$$
c$$$      allocate( matrix_send
c$$$     .     ( mpi__Mnall, mpi__Ms:mpi__Me, mpi__Wnall*mpi__Pn) )
c$$$      allocate( matrix_recv
c$$$     .     ( mpi__Mnall, mpi__Ms:mpi__Me, mpi__Wn*mpi__Pnall) )
c$$$      allocate( vcount_send(0:mpi__sizeP-1) )
c$$$      allocate( vdispl_send(0:mpi__sizeP-1) )
c$$$      allocate( vcount_recv(0:mpi__sizeP-1) )
c$$$      allocate( vdispl_recv(0:mpi__sizeP-1) )
c$$$
c$$$      index=0
c$$$      do p=0, mpi__sizeP-1
c$$$         do i=mpi__Ps, mpi__Pe
c$$$            do j=mpi__Wvs(p), mpi__Wve(p)
c$$$               index=index+1
c$$$               matrix_send(:,:,index) = matrixA(:,:,j,i)
c$$$            end do
c$$$         end do
c$$$
c$$$         vcount_send(p) = mpi__Mnall*mpi__Mn*(mpi__Wvn(p)-0)*mpi__Pn
c$$$         vdispl_send(p) = mpi__Mnall*mpi__Mn*(mpi__Wvs(p)-1)*mpi__Pn
c$$$         vcount_recv(p) = mpi__Mnall*mpi__Mn*mpi__Wn*(mpi__Pvn(p)-0)
c$$$         vdispl_recv(p) = mpi__Mnall*mpi__Mn*mpi__Wn*(mpi__Pvs(p)-1)
c$$$      end do
c$$$
c$$$      call MPI_AlltoAllv(
c$$$     .     matrix_send, vcount_send, vdispl_send, MPI_DOUBLE_COMPLEX,
c$$$     .     matrix_recv, vcount_recv, vdispl_recv, MPI_DOUBLE_COMPLEX,
c$$$     .     mpi__commP, mpi__info )
c$$$
c$$$      index=0
c$$$      do p=0, mpi__sizeP-1
c$$$         do i=mpi__Pvs(p), mpi__Pve(p)
c$$$            do j=mpi__Ws, mpi__We
c$$$               index=index+1
c$$$               matrixB(:,:,j,i) = matrix_recv(:,:,index)
c$$$            end do
c$$$         end do
c$$$      end do
c$$$      
c$$$      deallocate( matrix_send, matrix_recv )
c$$$      deallocate( vcount_send, vcount_recv )
c$$$      deallocate( vdispl_send, vdispl_recv )
c$$$      
c$$$      return
c$$$      end subroutine MPI__AllMPtoAllMW
c$$$#endif 
c$$$      end module m_mpi
