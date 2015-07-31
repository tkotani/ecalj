#undef SCALAPACK
      module m_mpi
      implicit none
      include "mpif.h"

      integer :: mpi__info
      integer :: mpi__size
      integer :: mpi__sizeQ = 1
      integer :: mpi__sizeS = 1
      integer :: mpi__sizeP = 1
      integer :: mpi__sizeB = 1
      integer :: mpi__sizeM = 1
      integer :: mpi__sizeW = 1
      integer :: mpi__sizeT = 1 ! for SXCFFAL2

      integer :: mpi__rank
      integer :: mpi__rankQ
      integer :: mpi__rankS
      integer :: mpi__rankP
      integer :: mpi__rankB
      integer :: mpi__rankM
      integer :: mpi__rankW
      integer :: mpi__rankT ! for SXCFFAL2

      logical :: mpi__root
      logical :: mpi__rootQ
      logical :: mpi__rootS
      logical :: mpi__rootP
      logical :: mpi__rootB
      logical :: mpi__rootM
      logical :: mpi__rootT ! for SXCFFAL2

      integer :: mpi__comm = MPI_COMM_WORLD
      integer :: mpi__commQ
      integer :: mpi__commS
      integer :: mpi__commP
      integer :: mpi__commB
      integer :: mpi__commM
      integer :: mpi__commT ! for SXCFFAL2

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

C     ! for T(nt) parallelization, for SXCFFAL2
      integer :: mpi__Tnall
      integer :: mpi__Tn, mpi__Ts, mpi__Te
      integer, allocatable :: mpi__Tvn(:), mpi__Tvs(:), mpi__Tve(:)
      integer :: mpi__T0n, mpi__T0s, mpi__T0e

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

      real*8 :: timer_s(20), timer_e(20)

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
      subroutine MPI__hx0fp0_initialize
      implicit none

      character(1024*4) :: cwd, stdout
      integer :: narg, iargc
      character(len=1024) :: arg

      integer :: mpi__group
      integer :: mpi__groupQ
      integer :: mpi__groupS
      integer :: mpi__groupP
      integer :: mpi__groupB
      integer :: mpi__groupM

      integer, allocatable :: ranklistQ(:)
      integer, allocatable :: ranklistS(:)
      integer, allocatable :: ranklistP(:)
      integer, allocatable :: ranklistB(:)
      integer, allocatable :: ranklistM(:)

      integer, allocatable :: vsizeQ(:)
      integer, allocatable :: vsizeS(:)
      integer, allocatable :: vsizeP(:)
      integer, allocatable :: vsizeB(:)
      integer, allocatable :: vsizeM(:)

      integer, allocatable :: vrankQ(:)
      integer, allocatable :: vrankS(:)
      integer, allocatable :: vrankP(:)
      integer, allocatable :: vrankB(:)
      integer, allocatable :: vrankM(:)

      integer :: i, j, n
      integer, allocatable :: imap(:)

C     ! intialize MPI and
      call getcwd(cwd)          ! get current working directory

      call MPI_Init( mpi__info ) ! current working directory is changed if mpirun is not used
      call MPI_Comm_rank ( MPI_COMM_WORLD, mpi__rank,  mpi__info )
      call MPI_Comm_size ( MPI_COMM_WORLD, mpi__size,  mpi__info )
      call MPI_Comm_group( MPI_COMM_WORLD, mpi__group, mpi__info )

      mpi__root = ( mpi__rank == 0 )

      if( mpi__root ) then
         call chdir(cwd)        ! recover current working directory
      endif

C     ! get arguments
      narg = iargc()
      do i=1, narg
         call getarg(i,arg)
           if( .false. ) then
         else if( arg == "-nq" ) then
            call getarg(i+1,arg)
            read(arg,*) mpi__sizeQ
         else if( arg == "-ns" ) then
            call getarg(i+1,arg)
            read(arg,*) mpi__sizeS
         else if( arg == "-np" ) then
            call getarg(i+1,arg)
            read(arg,*) mpi__sizeP
         else if( arg == "-nb" ) then
            call getarg(i+1,arg)
            read(arg,*) mpi__sizeB
         else if( arg == "-nm" ) then
            call getarg(i+1,arg)
            read(arg,*) mpi__sizeM
         end if
      end do
      mpi__sizeW = mpi__sizeP

C     ! in case no arguments, set as default
      if(  mpi__sizeQ == 1 .and. 
     .     mpi__sizeS == 1 .and. 
     .     mpi__sizeP == 1 .and. 
     .     mpi__sizeB == 1 .and. 
     .     mpi__sizeM == 1 ) then
         mpi__sizeQ = mpi__size
         mpi__sizeS = 1
         mpi__sizeP = 1
         mpi__sizeB = 1
         mpi__sizeM = 1
      end if

      if (mpi__root) then
      write(6,'(a,10I5)')'mpi_size,Q,S,P,B,M=', mpi__size,
     .   mpi__sizeQ,mpi__sizeS,mpi__sizeP,mpi__sizeB,mpi__sizeM
      endif

C     ! check inconsistent settings
      if( mpi__sizeQ*mpi__sizeS*mpi__sizeP*mpi__sizeB*mpi__sizeM /= mpi__size ) then
         write(6,*) "MPI size mismatch, mpirun -np size program -nq NQ -ns NS -np NP -nb NB -nm NM..."
         write(6,*) "size must be NQ*NS*NP*NB*NM"
         call MPI_Finalize ( mpi__info )
         call rx( ' MPI size mismatch' )
      end if

C     ! determine ranks
      allocate( vrankQ(0:mpi__size-1) )
      allocate( vrankS(0:mpi__size-1) )
      allocate( vrankP(0:mpi__size-1) )
      allocate( vrankB(0:mpi__size-1) )
      allocate( vrankM(0:mpi__size-1) )
      
      do i=0, mpi__size-1
         n = i
         vrankM(i) = mod(n,mpi__sizeM)
         n = n / mpi__sizeM
         vrankB(i) = mod(n,mpi__sizeB)
         n = n / mpi__sizeB
         vrankP(i) = mod(n,mpi__sizeP)
         n = n / mpi__sizeP
         vrankS(i) = mod(n,mpi__sizeS)
         n = n / mpi__sizeS
         vrankQ(i) = mod(n,mpi__sizeQ)
      end do

      if (mpi__root) then
      write(6,*)'MPI rank information'
      write(6,'(a,$)') 'vrankQ='; write(6,'(10i5)') vrankQ
      write(6,'(a,$)') 'vrankS='; write(6,'(10i5)') vrankS
      write(6,'(a,$)') 'vrankP='; write(6,'(10i5)') vrankP
      write(6,'(a,$)') 'vrankB='; write(6,'(10i5)') vrankB
      write(6,'(a,$)') 'vrankM='; write(6,'(10i5)') vrankM
      endif

C     ! rank of this process
      mpi__rankQ = vrankQ(mpi__rank)
      mpi__rankS = vrankS(mpi__rank)
      mpi__rankP = vrankP(mpi__rank)
      mpi__rankB = vrankB(mpi__rank)
      mpi__rankM = vrankM(mpi__rank)
      mpi__rankW = mpi__rankP


C     ! root or not of this process
      mpi__rootQ = ( mpi__rankQ == 0 )
      mpi__rootS = ( mpi__rankS == 0 )
      mpi__rootP = ( mpi__rankP == 0 )
      mpi__rootB = ( mpi__rankB == 0 )
      mpi__rootM = ( mpi__rankM == 0 )
      
      
C     ! setup Q-parallel
C     !  commQ for among other Q, same S,P,B,M
      allocate( ranklistQ(0:mpi__sizeQ-1) )
      i=0
      do j=0, mpi__size-1
         if(  vrankS(j) == mpi__rankS .and.
     .        vrankP(j) == mpi__rankP .and.
     .        vrankB(j) == mpi__rankB .and.
     .        vrankM(j) == mpi__rankM ) then
            ranklistQ(i) = j
            i=i+1
         end if
      end do
      call MPI_Group_incl( mpi__group,
     .     mpi__sizeQ, ranklistQ, mpi__groupQ, mpi__info )
      call MPI_Comm_create( MPI_COMM_WORLD,
     .     mpi__groupQ, mpi__commQ, mpi__info )
      deallocate( ranklistQ )


C     ! setup S-parallel
C     !  commS for among other S, same Q,P,B,M
      allocate( ranklistS(0:mpi__sizeS-1) )
      i=0
      do j=0, mpi__size-1
         if(  vrankQ(j) == mpi__rankQ .and.
     .        vrankP(j) == mpi__rankP .and.
     .        vrankB(j) == mpi__rankB .and.
     .        vrankM(j) == mpi__rankM ) then
            ranklistS(i) = j
            i=i+1
         end if
      end do
      call MPI_Group_incl( mpi__group,
     .     mpi__sizeS, ranklistS, mpi__groupS, mpi__info )
      call MPI_Comm_create( MPI_COMM_WORLD,
     .     mpi__groupS, mpi__commS, mpi__info )
      deallocate( ranklistS )


C     ! setup P-parallel
C     !  commP for among other P, same Q,B,M
      allocate( ranklistP(0:mpi__sizeP-1) )
      i=0
      do j=0, mpi__size-1
         if(  vrankQ(j) == mpi__rankQ .and.
     .        vrankS(j) == mpi__rankS .and.
     .        vrankB(j) == mpi__rankB .and.
     .        vrankM(j) == mpi__rankM ) then
            ranklistP(i) = j
            i=i+1
         end if
      end do
      call MPI_Group_incl( mpi__group,
     .     mpi__sizeP, ranklistP, mpi__groupP, mpi__info )
      call MPI_Comm_create( MPI_COMM_WORLD,
     .     mpi__groupP, mpi__commP, mpi__info )
      deallocate( ranklistP )
      
C     ! setup B-parallel
C     !  commB for among other B, same Q,P,M
      allocate( ranklistB(0:mpi__sizeB-1) )
      i=0
      do j=0, mpi__size-1
         if(  vrankQ(j) == mpi__rankQ .and.
     .        vrankS(j) == mpi__rankS .and.
     .        vrankP(j) == mpi__rankP .and.
     .        vrankM(j) == mpi__rankM ) then
            ranklistB(i) = j
            i=i+1
         end if
      end do
      call MPI_Group_incl( mpi__group,
     .     mpi__sizeB, ranklistB, mpi__groupB, mpi__info )
      call MPI_Comm_create( MPI_COMM_WORLD,
     .     mpi__groupB, mpi__commB, mpi__info )
      deallocate( ranklistB )

C     ! setup M-parallel      
C     !  commM for among other M, same Q,P,B
      allocate( ranklistM(0:mpi__sizeM-1) )
      i=0
      do j=0, mpi__size-1
         if(  vrankQ(j) == mpi__rankQ .and.
     .        vrankS(j) == mpi__rankS .and.
     .        vrankP(j) == mpi__rankP .and.
     .        vrankB(j) == mpi__rankB ) then
            ranklistM(i) = j
            i=i+1
         end if
      end do
      call MPI_Group_incl( mpi__group,
     .     mpi__sizeM, ranklistM, mpi__groupM, mpi__info )
      call MPI_Comm_create( MPI_COMM_WORLD,
     .     mpi__groupM, mpi__commM, mpi__info )
      deallocate( ranklistM )

#if SCALAPACK
C     ! setup ScaLapack process list
      allocate( imap(0:mpi__sizeM-1) )
      
      do i=0, mpi__sizeM-1
         imap(i) = mpi__rankQ * (mpi__sizeS*mpi__sizeP*mpi__sizeB*mpi__sizeM)
     .        +    mpi__rankS * (mpi__sizeP*mpi__sizeB*mpi__sizeM)
     .        +    mpi__rankP * (mpi__sizeB*mpi__sizeM)
     .        +    mpi__rankB * (mpi__sizeM) + i
      end do

      call BLACS_GET( 0, 0, mpi__Mhandle )
      call BLACS_GRIDMAP( mpi__Mhandle, imap, 1, 1, mpi__sizeM )
      deallocate(imap)
#endif
      return
      end subroutine MPI__hx0fp0_initialize

      
C     !======================================================
      subroutine MPI__sxcf_initialize
      implicit none

      character(1024*4) :: cwd, stdout
      integer :: narg, iargc
      character(len=1024) :: arg

      integer :: mpi__group
      integer :: mpi__groupQ
      integer :: mpi__groupP
      integer :: mpi__groupT

      integer, allocatable :: ranklistQ(:)
      integer, allocatable :: ranklistP(:)
      integer, allocatable :: ranklistT(:)

      integer, allocatable :: vsizeQ(:)
      integer, allocatable :: vsizeP(:)
      integer, allocatable :: vsizeT(:)

      integer, allocatable :: vrankQ(:)
      integer, allocatable :: vrankP(:)
      integer, allocatable :: vrankT(:)

      integer :: i, j, n
      integer, allocatable :: imap(:)

C     ! intialize MPI and
      call getcwd(cwd)          ! get current working directory

      call MPI_Init( mpi__info ) ! current working directory is changed if mpirun is not used
      call MPI_Comm_rank ( MPI_COMM_WORLD, mpi__rank,  mpi__info )
      call MPI_Comm_size ( MPI_COMM_WORLD, mpi__size,  mpi__info )
      call MPI_Comm_group( MPI_COMM_WORLD, mpi__group, mpi__info )

      mpi__root = ( mpi__rank == 0 )

      if( mpi__root ) then
         call chdir(cwd)        ! recover current working directory
      endif

C     ! get arguments
      narg = iargc()
      do i=1, narg
         call getarg(i,arg)
           if( .false. ) then
         else if( arg == "-nq" ) then
            call getarg(i+1,arg)
            read(arg,*) mpi__sizeQ
         else if( arg == "-np" ) then
            call getarg(i+1,arg)
            read(arg,*) mpi__sizeP
         else if( arg == "-nt" ) then
            call getarg(i+1,arg)
            read(arg,*) mpi__sizeT
         end if
      end do
      mpi__sizeW = mpi__sizeP

C     ! in case no arguments, set as default
      if(  mpi__sizeQ == 1 .and. 
     .     mpi__sizeP == 1 .and. 
     .     mpi__sizeT == 1 ) then
         mpi__sizeQ = mpi__size
         mpi__sizeP = 1
         mpi__sizeT = 1
         mpi__sizeM = 1
      end if

      if (mpi__root) then
      write(6,'(a,10I5)')'mpi_size,Q,P,T=', mpi__size,
     .   mpi__sizeQ,mpi__sizeP,mpi__sizeT
      endif

C     ! check inconsistent settings
      if( mpi__sizeQ*mpi__sizeP*mpi__sizeT /= mpi__size ) then
         write(6,*) "MPI size mismatch, mpirun -np size program -nq NQ -np NP -nt NT"
         write(6,*) "size must be NQ*NP*NT"
         call MPI_Finalize ( mpi__info )
         call rx( ' MPI size mismatch' )
      end if

C     ! determine ranks
      allocate( vrankQ(0:mpi__size-1) )
      allocate( vrankP(0:mpi__size-1) )
      allocate( vrankT(0:mpi__size-1) )

      do i=0, mpi__size-1
         n = i
         vrankT(i) = mod(n,mpi__sizeT)
         n = n / mpi__sizeT
         vrankP(i) = mod(n,mpi__sizeP)
         n = n / mpi__sizeP
         vrankQ(i) = mod(n,mpi__sizeQ)
      end do

      if (mpi__root) then
      write(6,*)'MPI rank information'
      write(6,'(a,$)') 'vrankQ='; write(6,'(10i5)') vrankQ
      write(6,'(a,$)') 'vrankP='; write(6,'(10i5)') vrankP
      write(6,'(a,$)') 'vrankT='; write(6,'(10i5)') vrankT
      endif

C     ! rank of this process
      mpi__rankQ = vrankQ(mpi__rank)
      mpi__rankP = vrankP(mpi__rank)
      mpi__rankT = vrankT(mpi__rank)


C     ! root or not of this process
      mpi__rootQ = ( mpi__rankQ == 0 )
      mpi__rootP = ( mpi__rankP == 0 )
      mpi__rootT = ( mpi__rankT == 0 )
      
      
C     ! setup Q-parallel
C     !  commQ for among other Q, same P,T
      allocate( ranklistQ(0:mpi__sizeQ-1) )
      i=0
      do j=0, mpi__size-1
         if(  vrankP(j) == mpi__rankP .and.
     .        vrankT(j) == mpi__rankT ) then
            ranklistQ(i) = j
            i=i+1
         end if
      end do
      call MPI_Group_incl( mpi__group,
     .     mpi__sizeQ, ranklistQ, mpi__groupQ, mpi__info )
      call MPI_Comm_create( MPI_COMM_WORLD,
     .     mpi__groupQ, mpi__commQ, mpi__info )
      deallocate( ranklistQ )


C     ! setup P-parallel
C     !  commP for among other P, same Q,T
      allocate( ranklistP(0:mpi__sizeP-1) )
      i=0
      do j=0, mpi__size-1
         if(  vrankQ(j) == mpi__rankQ .and.
     .        vrankT(j) == mpi__rankT ) then
            ranklistP(i) = j
            i=i+1
         end if
      end do
      call MPI_Group_incl( mpi__group,
     .     mpi__sizeP, ranklistP, mpi__groupP, mpi__info )
      call MPI_Comm_create( MPI_COMM_WORLD,
     .     mpi__groupP, mpi__commP, mpi__info )
      deallocate( ranklistP )
      

C     ! setup T-parallel      
C     !  commT for among other T, same Q,P
      allocate( ranklistT(0:mpi__sizeT-1) )
      i=0
      do j=0, mpi__size-1
         if(  vrankQ(j) == mpi__rankQ .and.
     .        vrankP(j) == mpi__rankP ) then
            ranklistT(i) = j
            i=i+1
         end if
      end do
      call MPI_Group_incl( mpi__group,
     .     mpi__sizeT, ranklistT, mpi__groupT, mpi__info )
      call MPI_Comm_create( MPI_COMM_WORLD,
     .     mpi__groupT, mpi__commT, mpi__info )
      deallocate( ranklistT )

      return
      end subroutine MPI__sxcf_initialize



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
      subroutine MPI__Barrier
      implicit none

      call MPI_Barrier( MPI_COMM_WORLD, mpi__info )

      end subroutine MPI__Barrier



C     !======================================================
      subroutine MPI__Finalize
      implicit none
#if SCALAPACK
      if( mpi__Mhandle /= 0 ) then
         call BLACS_GRIDEXIT( mpi__Mhandle )
      end if
#endif
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
      subroutine MPI__AllreduceSumQP( data, sizex )
      implicit none
      integer, intent(in) :: sizex
      complex(8), intent(inout) :: data(sizex)
      complex(8), allocatable   :: mpi__data(:) 

      if( mpi__sizeQ == 1 .and. mpi__sizeP == 1 ) return

      allocate(mpi__data(sizex))
      mpi__data = data
      call MPI_Allreduce( mpi__data, data, sizex,
     &     MPI_DOUBLE_COMPLEX, MPI_SUM, mpi__commP, mpi__info )
      mpi__data = data
      call MPI_Allreduce( mpi__data, data, sizex,
     &     MPI_DOUBLE_COMPLEX, MPI_SUM, mpi__commQ, mpi__info )

      deallocate( mpi__data )

      return
      end subroutine MPI__AllreduceSumQP

C     !======================================================
      subroutine MPI__AllreduceSumT( data, sizex )
      implicit none
      integer, intent(in) :: sizex
      complex(8), intent(inout) :: data(sizex)
      complex(8), allocatable   :: mpi__data(:) 

      if( mpi__sizeT == 1 ) return

      allocate(mpi__data(sizex))
      mpi__data = data

      call MPI_Allreduce( mpi__data, data, sizex,
     &     MPI_DOUBLE_COMPLEX, MPI_SUM, mpi__commT, mpi__info )

      deallocate( mpi__data )

      return
      end subroutine MPI__AllreduceSumT

C     !======================================================
      subroutine MPI__AllreduceSumP( data, sizex )
      implicit none
      integer, intent(in) :: sizex
      complex(8), intent(inout) :: data(sizex)
      complex(8), allocatable   :: mpi__data(:) 

      if( mpi__sizeP == 1 ) return

      allocate(mpi__data(sizex))
      mpi__data = data

      call MPI_Allreduce( mpi__data, data, sizex,
     &     MPI_DOUBLE_COMPLEX, MPI_SUM, mpi__commP, mpi__info )

      deallocate( mpi__data )

      return
      end subroutine MPI__AllreduceSumP


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


C     !======================================================
      subroutine MPI__AllreduceSumS( data, sizex )
      implicit none
      integer, intent(in) :: sizex
      complex(8), intent(inout) :: data(sizex)
      complex(8), allocatable   :: mpi__data(:)

      if( mpi__sizeS == 1 ) return

      allocate(mpi__data(sizex))
      mpi__data = data

      call MPI_Allreduce( mpi__data, data, sizex,
     &     MPI_DOUBLE_COMPLEX, MPI_SUM, mpi__commS, mpi__info )

      deallocate( mpi__data )

      return
      end subroutine MPI__AllreduceSumS

C     !======================================================
C     !! Aall = gather A
      subroutine MPI__AllgatherZmelP( zmel, ngb )
      implicit none
      integer ,intent(in) :: ngb
      complex(8), intent(inout) :: zmel(ngb,mpi__Tnall,mpi__Pnall)
      complex(8), allocatable   :: zmel_mpi(:,:,:)

      if( mpi__sizeP == 1 ) then
         return
      end if

      ! another buffer is required by some MPI versions
      allocate( zmel_mpi(ngb,mpi__Ts:mpi__Te,mpi__Pnall) )
      zmel_mpi(:,mpi__Ts:mpi__Te,mpi__Ps:mpi__Pe) = zmel(:,mpi__Ts:mpi__Te,mpi__Ps:mpi__Pe)
c$$$      allocate( zmel_mpi(ngb,mpi__Tnall,mpi__Ps:mpi__Pe) )
c$$$      zmel_mpi(:,:,mpi__Ps:mpi__Pe) = zmel(:,:,mpi__Ps:mpi__Pe)

      call MPI_Allgatherv( 
     .     zmel_mpi(1,mpi__Ts,mpi__Ps),
     .     ngb*mpi__Tn*mpi__Pn,
     .     MPI_DOUBLE_COMPLEX, zmel_mpi,
     .     ngb*mpi__Tn*mpi__Pvn(:),
     .     ngb*mpi__Tn*(mpi__Pvs(:)-1),
     .     MPI_DOUBLE_COMPLEX, mpi__commP, mpi__info )
c$$$      call MPI_Allgatherv( 
c$$$     .     zmel_mpi(1,1,mpi__Ps),
c$$$     .     ngb*mpi__Tnall*mpi__Pn,
c$$$     .     MPI_DOUBLE_COMPLEX, zmel,
c$$$     .     ngb*mpi__Tnall*(mpi__Pvn(:)+0),
c$$$     .     ngb*mpi__Tnall*(mpi__Pvs(:)-1),
c$$$     .     MPI_DOUBLE_COMPLEX, mpi__commP, mpi__info )
      zmel(:,mpi__Ts:mpi__Te,1:mpi__Pnall) = zmel_mpi(:,mpi__Ts:mpi__Te,1:mpi__Pnall)

      deallocate( zmel_mpi )

      return
      end subroutine MPI__AllgatherZmelP

C     !======================================================
      subroutine MPI__hsfp0_rankdivider(irkip,irkip_all,nspinmx,nqibz,ngrp,nq)
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
      write(6,"('MPI__hsfp0_rankdivider:$')")
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
      end subroutine MPI__hsfp0_rankdivider


C     !======================================================
      subroutine MPI__hsfp0_rankdividerQ(irkip,irkip_all,nspinmx,nqibz,ngrp,nq)
      implicit none

      integer, intent(out) :: irkip    (nspinmx,nqibz,ngrp,nq)
      integer, intent(in)  :: irkip_all(nspinmx,nqibz,ngrp,nq)
      integer, intent(in)  :: nspinmx,nqibz,ngrp,nq
      
      integer :: ispinmx,iqibz,igrp,iq
      integer :: total
      integer, allocatable :: vtotal(:)
      integer :: indexi, indexe
      integer :: p

      if( mpi__sizeQ == 1 ) then
         irkip = irkip_all
         return
      end if

c$$$      write(6,*) "irkip_all", mpi__rankQ
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
      write(6,"('MPI__hsfp0_rankdividerQ:$')")
      write(6,"('nspinmx,nqibz,ngrp,nq,total=',5i6)")
     &     nspinmx,nqibz,ngrp,nq,total 
      
      allocate( vtotal(0:mpi__sizeQ-1) )
      vtotal(:) = total/mpi__sizeQ

      do p=1, mod(total,mpi__sizeQ)
         vtotal(p-1) = vtotal(p-1) + 1
      end do

#if 1
      write(6,*) "DEBUG: mpi__Qnall=", sum(vtotal)
      write(6,*) "DEBUG: mpi__Qn=", vtotal(mpi__rankQ)
#endif

      indexe=0
      do p=0, mpi__rankQ
         indexi = indexe+1
         indexe = indexi+vtotal(p)-1
      end do
      deallocate(vtotal)

      total = 0
      irkip(:,:,:,:) = 0

      write(6,*) "irkip", mpi__rankQ
      do ispinmx=1, nspinmx
         do iq=1, nq
            do iqibz=1, nqibz
               do igrp=1, ngrp
                  if( irkip_all(ispinmx,iqibz,igrp,iq) /= 0 ) then
                     total = total + 1
                     if( indexi<=total .and. total<=indexe ) then
                        irkip(ispinmx,iqibz,igrp,iq) = irkip_all(ispinmx,iqibz,igrp,iq)
                        write(6,*) ispinmx, iq, iqibz, igrp, irkip(ispinmx,iqibz,igrp,iq)
                     end if
                  end if
               end do
            end do
         end do
      end do

      return
      end subroutine MPI__hsfp0_rankdividerQ

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
      subroutine MPI__hx0fp0_rankdivider2Q(iqxini,iqxend)
      implicit none
      integer, intent(in) :: iqxini, iqxend
      integer :: iq,i
      allocate( mpi__Qtask(1:iqxend), mpi__Qranktab(1:iqxend) )

      mpi__Qtask(:) = .false.
      mpi__Qranktab(1:iqxend)=999999
      if( mpi__sizeQ == 1 ) then
         mpi__Qtask(:) = .true.
         mpi__Qranktab(:) = mpi__rankQ
         return
      end if
      if(mpi__root) write(6,*) "MPI_hx0fp0_rankdivider2Q:"
      if(mpi__root) then
         write(6,'(a,$)')'mpi__Qtask='
         write(6,'(10L2)')mpi__Qtask
      endif
      do iq=iqxini, iqxend
         mpi__Qranktab(iq) = mod(iq-1,mpi__sizeQ)  !rank_table for given iq. iq=1 must give rank=0
         if( mpi__Qranktab(iq) == mpi__rankQ) then
            mpi__Qtask(iq) = .true.               !mpi__task is nodeID-dependent.
         endif
         if(mpi__root) then
           write(6,"('  iq irank=',2i5)")iq,mpi__Qranktab(iq)
         endif
      enddo   

      return
      end subroutine MPI__hx0fp0_rankdivider2Q


C     !======================================================
      subroutine MPI__hx0fp0_rankdivider2S( nspinmx )
      implicit none

      integer, intent(in) :: nspinmx
      
      integer :: k, jpm, ibib
      integer :: p, n
      integer :: blocksize
      integer :: max_nbnb

C     ! for S(npm) parallelization
      mpi__Snall = nspinmx

      if( allocated(mpi__Svn) ) deallocate(mpi__Svn)
      if( allocated(mpi__Svs) ) deallocate(mpi__Svs)
      if( allocated(mpi__Sve) ) deallocate(mpi__Sve)
      allocate( mpi__Svn(0:mpi__sizeS-1) )
      allocate( mpi__Svs(0:mpi__sizeS-1) )
      allocate( mpi__Sve(0:mpi__sizeS-1) )

      mpi__Svn(:) = mpi__Snall/mpi__sizeS
      if( mpi__Svn(0) == 0 ) then
         write(6,*) "Warning: too many processes for S parallelization"
         write(6,*) "  maximum -np ", mpi__Snall
         call rx('Warning: too many processes for S parallelization')
      end if

      do p=0, mpi__sizeS-1
         if( sum(mpi__Svn(:)) == mpi__Snall ) exit
         mpi__Svn(p) = mpi__Svn(p) + 1
      end do

      n=1
      do p=0, mpi__sizeS-1
         mpi__Svs(p) = n
         mpi__Sve(p) = mpi__Svs(p) + mpi__Svn(p) - 1
         n = n + mpi__Svn(p)
      end do
      mpi__Sn = mpi__Svn(mpi__rankS)
      mpi__Ss = mpi__Svs(mpi__rankS)
      mpi__Se = mpi__Sve(mpi__rankS)

      if (mpi__root) then
      write(6,*)mpi__rank,'mpi__Svn',mpi__Svn
      write(6,*)mpi__rank,'mpi__Svs',mpi__Svs
      write(6,*)mpi__rank,'mpi__Sve',mpi__Sve
      endif
#if 0
      write(6,*) "DEBUG: mpi__Qnall=", size(mpi__Qtask)
      write(6,*) "DEBUG: mpi__Snall=", mpi__Snall

      write(6,*) "DEBUG: mpi__Qn=", count(mpi__Qtask)
      write(6,*) "DEBUG: mpi__Sn=", mpi__Sn
#endif
      return
      end subroutine MPI__hx0fp0_rankdivider2S

C     !======================================================
      subroutine MPI__x0kf_rankdividerPBM( nbnb, nqbz, npm, ngb, nwt )
      implicit none

      integer, intent(in) :: nbnb(nqbz,npm), nqbz,npm
      integer, intent(in) :: ngb,nwt
      
      integer :: k, jpm, ibib
      integer :: p, n
      integer :: blocksize
      integer :: max_nbnb

C     ! for P(npm) parallelization
      mpi__Pnall = npm

      if( allocated(mpi__Pvn) ) deallocate(mpi__Pvn)
      if( allocated(mpi__Pvs) ) deallocate(mpi__Pvs)
      if( allocated(mpi__Pve) ) deallocate(mpi__Pve)
      allocate( mpi__Pvn(0:mpi__sizeP-1) )
      allocate( mpi__Pvs(0:mpi__sizeP-1) )
      allocate( mpi__Pve(0:mpi__sizeP-1) )

      mpi__Pvn(:) = mpi__Pnall/mpi__sizeP
      if( mpi__Pvn(0) == 0 ) then
         write(6,*) "Warning: too many processes for P parallelization"
         write(6,*) "  maximum -np ", mpi__Pnall
         call rx("Warning: too many processes for P parallelization")
      end if

      do p=0, mpi__sizeP-1
         if( sum(mpi__Pvn(:)) == mpi__Pnall ) exit
         mpi__Pvn(p) = mpi__Pvn(p) + 1
      end do

      n=1
      do p=0, mpi__sizeP-1
         mpi__Pvs(p) = n
         mpi__Pve(p) = mpi__Pvs(p) + mpi__Pvn(p) - 1
         n = n + mpi__Pvn(p)
      end do
      mpi__Pn = mpi__Pvn(mpi__rankP)
      mpi__Ps = mpi__Pvs(mpi__rankP)
      mpi__Pe = mpi__Pve(mpi__rankP)


C     ! for B(Kpoint&nbnb) parallelization
      mpi__Bnall=sum(nbnb(1:nqbz,mpi__Ps:mpi__Pe))

      
      if( allocated(mpi__Bvn) ) deallocate(mpi__Bvn)
      if( allocated(mpi__Bvs) ) deallocate(mpi__Bvs)
      if( allocated(mpi__Bve) ) deallocate(mpi__Bve)
      allocate( mpi__Bvn(0:mpi__sizeB-1) )
      allocate( mpi__Bvs(0:mpi__sizeB-1) )
      allocate( mpi__Bve(0:mpi__sizeB-1) )

      mpi__Bvn(:) = mpi__Bnall/mpi__sizeB
      if( mpi__Bvn(0) == 0 ) then
         write(6,*) "Warning: too many processes for B parallelization"
         write(6,*) "  maximum -nb ", mpi__Bnall
         call rx("Warning: too many processes for B parallelization")
      end if
      do p=0, mpi__sizeB-1
         if( sum(mpi__Bvn(:)) == mpi__Bnall ) exit
         mpi__Bvn(p) = mpi__Bvn(p) + 1
      end do

      n=1
      do p=0, mpi__sizeB-1
         mpi__Bvs(p) = n
         mpi__Bve(p) = mpi__Bvs(p) + mpi__Bvn(p) - 1
         n = n + mpi__Bvn(p)
      end do
      mpi__Bn = mpi__Bvn(mpi__rankB)
      mpi__Bs = mpi__Bvs(mpi__rankB)
      mpi__Be = mpi__Bve(mpi__rankB)

      max_nbnb=maxval(nbnb(1:nqbz,mpi__Ps:mpi__Pe))      

      if( allocated(mpi__Btask1) ) deallocate(mpi__Btask1)
      if( allocated(mpi__Btask2) ) deallocate(mpi__Btask2)
      if( allocated(mpi__Btask3) ) deallocate(mpi__Btask3)
      allocate( mpi__Btask1(nqbz) )
      allocate( mpi__Btask2(nqbz,mpi__Ps:mpi__Pe) )
      allocate( mpi__Btask3(nqbz,mpi__Ps:mpi__Pe,max_nbnb) )

      mpi__Btask1(:)     = .false.
      mpi__Btask2(:,:)   = .false.
      mpi__Btask3(:,:,:) = .false.

      n=1
      do k=1, nqbz
         do jpm=mpi__Ps, mpi__Pe
            do ibib=1, nbnb(k,jpm)
               if( mpi__Bs <= n .and. n<=mpi__Be ) then
                  mpi__Btask1(k) = .true.
                  mpi__Btask2(k,jpm) = .true.
                  mpi__Btask3(k,jpm,ibib) = .true.
               end if
               n=n+1
            end do
         end do
      end do


C     ! for M(matrix) parallelization
      if( allocated(mpi__Mvn) ) then
         deallocate( mpi__Mvn )
         deallocate( mpi__Mvs )
         deallocate( mpi__Mve )
      end if
      allocate( mpi__Mvn(0:mpi__sizeM-1) )
      allocate( mpi__Mvs(0:mpi__sizeM-1) )
      allocate( mpi__Mve(0:mpi__sizeM-1) )
    
      mpi__Mnall = ngb
      
C     !! the mean is larger
      mpi__Mvn(:) = int(ceiling(dble(mpi__Mnall)/dble(mpi__sizeM)))
C     !! the last block is smaller 
      mpi__Mvn(mpi__sizeM-1) = mpi__Mvn(mpi__sizeM-1)
     .     - (sum(mpi__Mvn(:))-mpi__Mnall)
      blocksize = maxval(mpi__Mvn(:))
#if SCALAPACK           
      call DescInit( mpi__Mdescv, 1, mpi__Mnall,
     .     1, blocksize, 0, 0,
     .     mpi__Mhandle, 1, mpi__info )
      
      call DescInit( mpi__Mdescm, mpi__Mnall, mpi__Mnall,
     .     blocksize, blocksize, 0, 0,
     .     mpi__Mhandle, mpi__Mnall, mpi__info )
#endif
      n=1
      do p=0, mpi__sizeM-1
         mpi__Mvs(p) = n
         mpi__Mve(p) = mpi__Mvs(p) + mpi__Mvn(p) - 1
         n = n + mpi__Mvn(p)
      end do
      mpi__Ms = mpi__Mvs(mpi__rankM)
      mpi__Me = mpi__Mve(mpi__rankM)
      mpi__Mn = mpi__Mvn(mpi__rankM)
      mpi__Mf = mpi__Ms+blocksize-1 ! equivalent end among all processes

C     ! for W(iw) parallelization
      mpi__Wnall = nwt

      if( allocated(mpi__Wvn) ) deallocate(mpi__Wvn)
      if( allocated(mpi__Wvs) ) deallocate(mpi__Wvs)
      if( allocated(mpi__Wve) ) deallocate(mpi__Wve)
      allocate( mpi__Wvn(0:mpi__sizeW-1) )
      allocate( mpi__Wvs(0:mpi__sizeW-1) )
      allocate( mpi__Wve(0:mpi__sizeW-1) )

      mpi__Wvn(:) = mpi__Wnall/mpi__sizeB
      if( mpi__Wvn(0) == 0 ) then
         write(6,*) "Warning: too many processes for B parallelization"
         write(6,*) "  maximum -np ", mpi__Wnall
         call rx("Warning: too many processes for B parallelization")
      end if

      do p=0, mpi__sizeW-1
         if( sum(mpi__Wvn(:)) == mpi__Wnall ) exit
         mpi__Wvn(p) = mpi__Wvn(p) + 1
      end do

      n=1
      do p=0, mpi__sizeW-1
         mpi__Wvs(p) = n
         mpi__Wve(p) = mpi__Wvs(p) + mpi__Wvn(p) - 1
         n = n + mpi__Wvn(p)
      end do
      mpi__Wn = mpi__Wvn(mpi__rankW)
      mpi__Ws = mpi__Wvs(mpi__rankW)
      mpi__We = mpi__Wve(mpi__rankW)
#if 1
      write(6,*) "DEBUG: mpi__Qnall=", size(mpi__Qtask)
      write(6,*) "DEBUG: mpi__Pnall=", mpi__Pnall
      write(6,*) "DEBUG: mpi__Bnall=", mpi__Bnall
      write(6,*) "DEBUG: mpi__Mnall=", mpi__Mnall
      write(6,*) "DEBUG: mpi__Wnall=", mpi__Wnall

      write(6,*) "DEBUG: mpi__Qn=", count(mpi__Qtask)
      write(6,*) "DEBUG: mpi__Pn=", mpi__Pn
      write(6,*) "DEBUG: mpi__Bn=", mpi__Bn
      write(6,*) "DEBUG: mpi__Mn=", mpi__Mn
      write(6,*) "DEBUG: mpi__Wn=", mpi__Wn
#endif
      return
      end subroutine MPI__x0kf_rankdividerPBM



C     !======================================================
      subroutine MPI__sxcf_rankdividerPT( ntqxx, nctot, nbmax )
      implicit none

      integer, intent(in) :: ntqxx, nctot, nbmax
      
      integer :: k, jpm, ibib
      integer :: p, n
      integer :: blocksize
      integer :: max_nbnb


C     ! for P(np) parallelization
      mpi__Pnall = ntqxx

      if( allocated(mpi__Pvn) ) deallocate(mpi__Pvn)
      if( allocated(mpi__Pvs) ) deallocate(mpi__Pvs)
      if( allocated(mpi__Pve) ) deallocate(mpi__Pve)
      allocate( mpi__Pvn(0:mpi__sizeP-1) )
      allocate( mpi__Pvs(0:mpi__sizeP-1) )
      allocate( mpi__Pve(0:mpi__sizeP-1) )

      mpi__Pvn(:) = mpi__Pnall/mpi__sizeP
      if( mpi__Pvn(0) == 0 ) then
         write(6,*) "Warning: too many processes for P parallelization"
         write(6,*) "  maximum -np ", mpi__Pnall
         call rx("Warning: too many processes for P parallelization")
      end if

      do p=0, mpi__sizeP-1
         if( sum(mpi__Pvn(:)) == mpi__Pnall ) exit
         mpi__Pvn(p) = mpi__Pvn(p) + 1
      end do

      n=1
      do p=0, mpi__sizeP-1
         mpi__Pvs(p) = n
         mpi__Pve(p) = mpi__Pvs(p) + mpi__Pvn(p) - 1
         n = n + mpi__Pvn(p)
      end do
      mpi__Pn = mpi__Pvn(mpi__rankP)
      mpi__Ps = mpi__Pvs(mpi__rankP)
      mpi__Pe = mpi__Pve(mpi__rankP)

C     ! for T(nt) parallelization
      mpi__Tnall = nctot+nbmax

      if( allocated(mpi__Tvn) ) deallocate(mpi__Tvn)
      if( allocated(mpi__Tvs) ) deallocate(mpi__Tvs)
      if( allocated(mpi__Tve) ) deallocate(mpi__Tve)
      allocate( mpi__Tvn(0:mpi__sizeT-1) )
      allocate( mpi__Tvs(0:mpi__sizeT-1) )
      allocate( mpi__Tve(0:mpi__sizeT-1) )

      mpi__Tvn(:) = mpi__Tnall/mpi__sizeT
      if( mpi__Tvn(0) == 0 ) then
         write(6,*) "Warning: too many processes for T parallelization"
         write(6,*) "  maximum -nt ", mpi__Tnall
         call rx("Warning: too many processes for T parallelization")
      end if

      do p=0, mpi__sizeT-1
         if( sum(mpi__Tvn(:)) == mpi__Tnall ) exit
         mpi__Tvn(p) = mpi__Tvn(p) + 1
      end do

      n=1
      do p=0, mpi__sizeT-1
         mpi__Tvs(p) = n
         mpi__Tve(p) = mpi__Tvs(p) + mpi__Tvn(p) - 1
         n = n + mpi__Tvn(p)
      end do
      mpi__Tn = mpi__Tvn(mpi__rankT)
      mpi__Ts = mpi__Tvs(mpi__rankT)
      mpi__Te = mpi__Tve(mpi__rankT)

      if( mpi__Te<nctot+1 ) then
         mpi__T0s = 0
         mpi__T0e = -1
      else if( mpi__Ts<nctot+1 .and. nctot+1<=mpi__Te ) then
         mpi__T0s = 1
         mpi__T0e = mpi__Te - nctot
      else if( nctot+1<=mpi__Ts ) then
         mpi__T0s = mpi__Ts - nctot
         mpi__T0e = mpi__Te - nctot
      end if
      mpi__T0n = mpi__T0e - mpi__T0s + 1
#if 1
      write(6,*) "DEBUG: mpi__Pnall=", mpi__Pnall
      write(6,*) "DEBUG: mpi__Pn   =", mpi__Pn
      write(6,*) "DEBUG: mpi__Tnall=", mpi__Tnall
      write(6,*) "DEBUG: mpi__Tn   =", mpi__Tn
      write(6,*) "DEBUG: mpi__T0n  =", mpi__T0n
#endif
      end subroutine MPI__sxcf_rankdividerPT


#if SCALAPACK
C     !======================================================
      subroutine MPI__dummy
      implicit none

      call PZHEGVX              ! dummy call, to poll scalapack library

      return
      end subroutine MPI__dummy
      
      
C     !======================================================
C     !!  a = x*y
      subroutine MPI__ZDOTU( a, vectorX, vectorY )
      implicit none

      complex(8), intent(out) :: a
      complex(8), intent(in)  :: vectorX(mpi__Mnall)
      complex(8), intent(in)  :: vectorY(mpi__Mnall)

      call PZDOTU( mpi__Mnall, a,
     .     vectorX(mpi__Ms), 1, 1, mpi__Mdescv, 1,
     .     vectorY(mpi__Ms), 1, 1, mpi__Mdescv, 1 )

      return
      end subroutine MPI__ZDOTU

      
C     !======================================================
C     !!  a = conjg(x)*y
      subroutine MPI__ZDOTC( a, vectorX, vectorY )
      implicit none

      complex(8), intent(out) :: a
      complex(8), intent(in)  :: vectorX(mpi__Ms:mpi__Mf)
      complex(8), intent(in)  :: vectorY(mpi__Ms:mpi__Mf)

      call PZDOTC( mpi__Mnall, a,
     .     vectorX, 1, 1, mpi__Mdescv, 1,
     .     vectorY, 1, 1, mpi__Mdescv, 1 )

      return
      end subroutine MPI__ZDOTC

      
C     !======================================================
C     !!  y = A * x
      subroutine MPI__ZGEMV( matrixA, vectorX, vectorY )
      implicit none
      
      complex(8), intent(in)  :: matrixA(mpi__Mnall,mpi__Ms:mpi__Mf)
      complex(8), intent(in)  :: vectorX(mpi__Ms:mpi__Mf)
      complex(8), intent(out) :: vectorY(mpi__Ms:mpi__Mf)
      character,  parameter   :: transa = 'N'
      complex(8), parameter   :: alpha = (1.0d0,0.0d0)
      complex(8), parameter   :: beta  = (0.0d0,0.0d0)
      integer :: i, j

      call PZGEMV( transa, mpi__Mnall, mpi__Mnall, alpha,
     .     matrixA, 1, 1, mpi__Mdescm,
     .     vectorX, 1, 1, mpi__Mdescv, 1, beta,
     .     vectorY, 1, 1, mpi__Mdescv, 1 )

      return
      end subroutine MPI__ZGEMV
      
C     !======================================================
C     !! C = A*B
      subroutine MPI__ZGEMM( matrixA, matrixB, matrixC )
      implicit none

      complex(8), intent(in)  :: matrixA(mpi__Mnall,mpi__Ms:mpi__Mf)
      complex(8), intent(in)  :: matrixB(mpi__Mnall,mpi__Ms:mpi__Mf)
      complex(8), intent(out) :: matrixC(mpi__Mnall,mpi__Ms:mpi__Mf)
      character,  parameter   :: transa = 'N'
      character,  parameter   :: transb = 'N'
      complex(8), parameter   :: alpha = (1.0d0,0.0d0)
      complex(8), parameter   :: beta  = (0.0d0,0.0d0)
      integer :: i, j, k

      call PZGEMM( transa, transb,
     .     mpi__Mnall, mpi__Mnall, mpi__Mnall,
     .     alpha,
     .     matrixA, 1, 1, mpi__Mdescm,
     .     matrixB, 1, 1, mpi__Mdescm,
     .     beta,
     .     matrixC, 1, 1, mpi__Mdescm )

      return
      end subroutine MPI__ZGEMM

C     !======================================================
C     !!  C = x^H * y
      subroutine MPI__ZGERC( vectorX, vectorY, matrixC )
      implicit none

      complex(8), intent(in)  :: vectorX(mpi__Ms:mpi__Mf)
      complex(8), intent(in)  :: vectorY(mpi__Ms:mpi__Mf)
      complex(8), intent(out) :: matrixC(mpi__Mnall,mpi__Ms:mpi__Mf)
      character,  parameter   :: transa = 'C', transb = 'N'
      complex(8), parameter   :: alpha = (1.0d0,0.0d0)
      complex(8), parameter   :: beta  = (0.0d0,0.0d0)
      integer :: i, j

      call PZGEMM( transa, transb, mpi__Mnall, mpi__Mnall, 1, alpha,
     .     vectorX, 1, 1, mpi__Mdescv,
     .     vectorY, 1, 1, mpi__Mdescv, beta,
     .     matrixC, 1, 1, mpi__Mdescm )

      return
      end subroutine MPI__ZGERC

C     !======================================================
C     !!  C = a*A + C
      subroutine MPI__ZAMPM( alpha, matrixA, matrixC )
      implicit none

      real(8), intent(in)       :: alpha
      complex(8), intent(in)    :: matrixA(mpi__Mnall,mpi__Ms:mpi__Mf)
      complex(8), intent(inout) :: matrixC(mpi__Mnall,mpi__Ms:mpi__Mf)
      integer :: i, j

!$OMP PARALLEL DO PRIVATE(j) FIRSTPRIVATE(mpi__Ms,mpi__Mf)
      do j=mpi__Ms,mpi__Mf
         matrixC(:,j) = alpha*matrixA(:,j) + matrixC(:,j)
      end do
!$OMP END PARALLEL

      return
      end subroutine MPI__ZAMPM


C     !======================================================
C     !! B = reduce A
      subroutine MPI__AllreduceMatrixB( matrixA )
      implicit none
      complex(8), intent(inout) :: matrixA(mpi__Mnall,mpi__Ms:mpi__Me)
      complex(8), allocatable   :: matrixA_sum(:,:)

      if( mpi__sizeB == 1 ) return

      allocate( matrixA_sum(mpi__Mnall,mpi__Ms:mpi__Me) )

      call MPI_Allreduce( matrixA, matrixA_sum,
     .     mpi__Mn*mpi__Mnall,
     .     MPI_DOUBLE_COMPLEX, MPI_SUM, mpi__commB, mpi__info )
      
      matrixA = matrixA_sum
      deallocate( matrixA_sum )

      return
      end subroutine MPI__AllreduceMatrixB

C     !======================================================
C     !! y = gather x
      subroutine MPI__AllgatherVectorM( vectorX, vectorY )
      implicit none
      complex(8), intent(in)  :: vectorX(mpi__Ms:mpi__Me)
      complex(8), intent(out) :: vectorY(mpi__Mnall)

      if( mpi__sizeM == 1 ) then
         vectorY = vectorX
         return
      end if

      call MPI_Allgatherv( vectorX(mpi__Ms), mpi__Mn,
     .     MPI_DOUBLE_COMPLEX,
     .     vectorY, (mpi__Mvn(:)+0), (mpi__Mvs(:)-1),
     .     MPI_DOUBLE_COMPLEX, mpi__commM, mpi__info )

      return
      end subroutine MPI__AllgatherVectorM

C     !======================================================
C     !! Aall = gather A
      subroutine MPI__AllgatherMatrixM( matrixA, matrixAall )
      implicit none
      complex(8), intent(in)   :: matrixA   (mpi__Mnall,mpi__Ms:mpi__Me)
      complex(8), intent(out)  :: matrixAall(mpi__Mnall,mpi__Mnall)

      if( mpi__sizeM == 1 ) then
         matrixAall = matrixA
         return
      end if

      call MPI_Allgatherv( matrixA(1,mpi__Ms), (mpi__Mn)*mpi__Mnall,
     .     MPI_DOUBLE_COMPLEX, matrixAall,
     .     (mpi__Mvn(:)+0)*mpi__Mnall, (mpi__Mvs(:)-1)*mpi__Mnall,
     .     MPI_DOUBLE_COMPLEX, mpi__commM, mpi__info )

      return
      end subroutine MPI__AllgatherMatrixM


C     !======================================================
C     !! Aall = gather A
      subroutine MPI__AllgatherArrayP( arrayA, arrayAall )
      implicit none
      complex(8), intent(in)   :: arrayA   
     .     (mpi__Mnall,mpi__Mnall,mpi__Wnall,mpi__Ps:mpi__Pe)
      complex(8), intent(out)  :: arrayAall
     .     (mpi__Mnall,mpi__Mnall,mpi__Wnall,mpi__Pnall)

      if( mpi__sizeP == 1 ) then
         arrayAall = arrayA
         return
      end if

      call MPI_Allgatherv( arrayA(1,1,1,mpi__Ps),
     .     mpi__Mnall*mpi__Mnall*mpi__Wnall*(mpi__Pn),
     .     MPI_DOUBLE_COMPLEX, arrayAall,
     .     mpi__Mnall*mpi__Mnall*mpi__Wnall*(mpi__Pvn(:)+0),
     .     mpi__Mnall*mpi__Mnall*mpi__Wnall*(mpi__Pvs(:)-1),
     .     MPI_DOUBLE_COMPLEX, mpi__commP, mpi__info )

      return
      end subroutine MPI__AllgatherArrayP


C     !======================================================
C     !! matrixB(:,s:e) = matrixA(s:e,:)
      subroutine MPI__AllMPtoAllMW( matrixA, matrixB )
      implicit none
      complex(8), intent(in)  :: matrixA
     .     (mpi__Mnall,mpi__Ms:mpi__Me,mpi__Wnall,mpi__Ps:mpi__Pe)
      complex(8), intent(out) :: matrixB
     .     (mpi__Mnall,mpi__Ms:mpi__Me,mpi__Ws:mpi__We,mpi__Pnall)

      complex(8), allocatable :: matrix_send(:,:,:)
      complex(8), allocatable :: matrix_recv(:,:,:)
      integer, allocatable :: vcount_send(:), vdispl_send(:)
      integer, allocatable :: vcount_recv(:), vdispl_recv(:)
      integer :: i, j, p, index

      if( mpi__sizeP == 1 ) then
         matrixB = matrixA
         return
      end if

      allocate( matrix_send
     .     ( mpi__Mnall, mpi__Ms:mpi__Me, mpi__Wnall*mpi__Pn) )
      allocate( matrix_recv
     .     ( mpi__Mnall, mpi__Ms:mpi__Me, mpi__Wn*mpi__Pnall) )
      allocate( vcount_send(0:mpi__sizeP-1) )
      allocate( vdispl_send(0:mpi__sizeP-1) )
      allocate( vcount_recv(0:mpi__sizeP-1) )
      allocate( vdispl_recv(0:mpi__sizeP-1) )

      index=0
      do p=0, mpi__sizeP-1
         do i=mpi__Ps, mpi__Pe
            do j=mpi__Wvs(p), mpi__Wve(p)
               index=index+1
               matrix_send(:,:,index) = matrixA(:,:,j,i)
            end do
         end do

         vcount_send(p) = mpi__Mnall*mpi__Mn*(mpi__Wvn(p)-0)*mpi__Pn
         vdispl_send(p) = mpi__Mnall*mpi__Mn*(mpi__Wvs(p)-1)*mpi__Pn
         vcount_recv(p) = mpi__Mnall*mpi__Mn*mpi__Wn*(mpi__Pvn(p)-0)
         vdispl_recv(p) = mpi__Mnall*mpi__Mn*mpi__Wn*(mpi__Pvs(p)-1)
      end do

      call MPI_AlltoAllv(
     .     matrix_send, vcount_send, vdispl_send, MPI_DOUBLE_COMPLEX,
     .     matrix_recv, vcount_recv, vdispl_recv, MPI_DOUBLE_COMPLEX,
     .     mpi__commP, mpi__info )

      index=0
      do p=0, mpi__sizeP-1
         do i=mpi__Pvs(p), mpi__Pve(p)
            do j=mpi__Ws, mpi__We
               index=index+1
               matrixB(:,:,j,i) = matrix_recv(:,:,index)
            end do
         end do
      end do
      
      deallocate( matrix_send, matrix_recv )
      deallocate( vcount_send, vcount_recv )
      deallocate( vdispl_send, vdispl_recv )
      
      return
      end subroutine MPI__AllMPtoAllMW
#endif 

      end module m_mpi
