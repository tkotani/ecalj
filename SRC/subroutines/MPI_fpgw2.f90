module m_mpi !MPI utility for fpgw
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
  !      integer :: mpi__rankS
  !      integer :: mpi__rankP
  !      integer :: mpi__rankB
  !      integer :: mpi__rankM
  !      integer :: mpi__rankW

  logical :: mpi__root
  logical :: mpi__rootQ
  !      logical :: mpi__rootS
  !      logical :: mpi__rootP
  !      logical :: mpi__rootB
  !      logical :: mpi__rootM

  integer :: mpi__comm = MPI_COMM_WORLD
  integer :: mpi__commQ
  !      integer :: mpi__commS
  !      integer :: mpi__commP
  !      integer :: mpi__commB
  !      integer :: mpi__commM

  integer:: ista(MPI_STATUS_SIZE )

  integer :: mpi__iini, mpi__iend

  logical, allocatable :: mpi__task(:)
  integer, allocatable :: mpi__ranktab(:)

  !     ! for Q(nq) parallelization
  logical, allocatable :: mpi__Qtask(:)
  integer, allocatable :: mpi__Qranktab(:)

  !     ! for S(nspin) parallelization
  logical, allocatable :: mpi__Stask(:)
  integer, allocatable :: mpi__Sranktab(:)
  integer :: mpi__Snall
  integer :: mpi__Sn, mpi__Ss, mpi__Se
  integer, allocatable :: mpi__Svn(:), mpi__Svs(:), mpi__Sve(:)


  !     ! for P(npm) parallelization
  integer :: mpi__Pnall
  integer :: mpi__Pn, mpi__Ps, mpi__Pe
  integer, allocatable :: mpi__Pvn(:), mpi__Pvs(:), mpi__Pve(:)


  !     ! for B(Kpoint&nbnb) parallelization
  integer :: mpi__Bnall
  integer :: mpi__Bn, mpi__Bs, mpi__Be
  integer, allocatable :: mpi__Bvn(:), mpi__Bvs(:), mpi__Bve(:)
  logical, allocatable :: mpi__Btask1(:)
  logical, allocatable :: mpi__Btask2(:,:)
  logical, allocatable :: mpi__Btask3(:,:,:)

  !     ! for M(matrix) parallelization
  integer :: mpi__Mnall
  integer :: mpi__Mn, mpi__Ms, mpi__Me, mpi__Mf
  integer, allocatable :: mpi__Mvn(:), mpi__Mvs(:), mpi__Mve(:)
  integer :: mpi__Mhandle=0
  integer :: mpi__Mdescv(9)
  integer :: mpi__Mdescm(9)

  !     ! for W(iwt) parallelization
  integer :: mpi__Wnall
  integer :: mpi__Wn, mpi__Ws, mpi__We
  integer, allocatable :: mpi__Wvn(:), mpi__Wvs(:), mpi__Wve(:)

  !     ! for magnon E(q) parallelization
  integer  :: mpi__MEq


  interface MPI__Send 
     module procedure &
          MPI__Send_i,  MPI__Send_iv, &
          MPI__Send_d,  MPI__Send_dv
  end interface MPI__Send

  interface MPI__Recv
     module procedure &
          MPI__Recv_i,  MPI__Recv_iv, &
          MPI__Recv_d,  MPI__Recv_dv 
  end interface MPI__Recv

contains

  !     !======================================================
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
    !! console-output from different nodes to different files
    !      if( mpi__size > 1 ) then
    !        if(mpi__root )write(6,*)'MPI console outputs to following files.'
    !        write(6,"('   stdout.',i4.4,'.',a)") mpi__rank,idn
    !        write(stdout,"('stdout.',i4.4,'.',a)") mpi__rank,idn
    !        open(unit=6,file=trim(stdout))
    !        write(6,"(a,i3)")" ### console output for rank=",mpi__rank
    !      endif
    !      if(mpi__root ) then
    !        close(unit=6)
    !      endif
    return
  end subroutine MPI__Initialize

  !     !======================================================
  subroutine MPI__Initialize_magnon()
    implicit none
    character(1024*4) :: cwd, stdout
    call getcwd(cwd)          ! get current working directory
    !$$$CC Reduce size for magnon (avoid memory leak)      
    !$$$      if (mpi__sizeMG > size_lim) then
    !$$$         mpi__sizeMG=size_lim
    !$$$      endif
    !$$$      if (mpi__sizeMG > 10) then
    !$$$         mpi__sizeMG=10
    !$$$      endif
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
    !! console-output from different nodes to different files
    !      if( mpi__size > 1 ) then
    !        if(mpi__root )write(6,*)'MPI console outputs to following files.'
    !        write(6,"('   stdout.',i4.4,'.',a)") mpi__rank,idn
    !        write(stdout,"('stdout.',i4.4,'.',a)") mpi__rank,idn
    !        open(unit=6,file=trim(stdout))
    !        write(6,"(a,i3)")" ### console output for rank=",mpi__rank
    !      endif
    !      if(mpi__root ) then
    !        close(unit=6)
    !      endif
    return
  end subroutine MPI__Initialize_magnon

  !     !======================================================
  subroutine MPI__InitializeQ !MPI__InitializeQSPBM not so much used now...
    implicit none
    character(1024*4) :: cwd, stdout
    integer :: narg, iargc
    character(len=1024) :: arg
    integer, allocatable :: ranklistQ(:)
    integer, allocatable :: vrankQ(:)
    integer :: i, j, n
    call getcwd(cwd)          ! get current working directory
    call MPI_Init( mpi__info ) ! current working directory is changed if mpirun is not used
    call MPI_Comm_rank ( MPI_COMM_WORLD, mpi__rank,  mpi__info )
    call MPI_Comm_size ( MPI_COMM_WORLD, mpi__size,  mpi__info )
    mpi__root = ( mpi__rank == 0 )
    if( mpi__root ) then
       call chdir(cwd)        ! recover current working directory
    endif
    return
  end subroutine MPI__InitializeQ!SPBM

  !     !======================================================
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

  !     !======================================================
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


  !     !======================================================
  subroutine MPI__Barrier
    implicit none

    call MPI_Barrier( MPI_COMM_WORLD, mpi__info )

  end subroutine MPI__Barrier



  !     !======================================================
  subroutine MPI__Finalize
    implicit none
    !#if SCALAPACK
    !      if( mpi__Mhandle /= 0 ) then
    !         call BLACS_GRIDEXIT( mpi__Mhandle )
    !      end if
    !#endif
    call MPI_Finalize ( mpi__info )

  end subroutine MPI__Finalize



  !     !======================================================
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


  !     !======================================================
  subroutine MPI__Broadcast( data )
    implicit none
    integer, intent(inout) :: data

    call MPI_Bcast( data, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi__info )

    return
  end subroutine MPI__Broadcast


  !     !======================================================
  subroutine MPI__send_d(data,dest)
    implicit none
    real(8):: data
    integer :: n,dest,ierr
    n=1
    call MPI_Send(data,n,MPI_REAL8,dest,mpi__rank, MPI_COMM_WORLD,ierr)
  end subroutine MPI__send_d

  !     !======================================================
  subroutine MPI__recv_d(data,src)
    implicit none
    real(8):: data
    integer :: n,src,ierr
    n=1
    call MPI_Recv(data,n,MPI_REAL8,src,src, MPI_COMM_WORLD,ista,ierr)
  end subroutine MPI__recv_d

  !     !======================================================
  subroutine MPI__send_dv(data,dest)
    implicit none
    real(8):: data(:)
    integer :: n,dest,ierr
    n=size(data)
    call MPI_Send(data,n,MPI_REAL8,dest,mpi__rank, MPI_COMM_WORLD,ierr)
  end subroutine MPI__send_dv

  !     !======================================================
  subroutine MPI__recv_dv(data,src)
    implicit none
    real(8):: data(:)
    integer :: n,src,ierr
    n=size(data)
    call MPI_Recv(data,n,MPI_REAL8,src,src, MPI_COMM_WORLD,ista,ierr)
  end subroutine MPI__recv_dv

  !     !======================================================
  subroutine MPI__send_i(data,dest)
    implicit none
    integer:: data
    integer :: n,dest,ierr
    n=1
    call MPI_Send(data,n,MPI_INTEGER,dest,mpi__rank, MPI_COMM_WORLD,ierr)
  end subroutine MPI__send_i

  !     !======================================================
  subroutine MPI__recv_i(data,src)
    implicit none
    integer:: data
    integer :: n,src,ierr
    n=1
    call MPI_Recv(data,n,MPI_INTEGER,src,src, MPI_COMM_WORLD,ista,ierr)
  end subroutine MPI__recv_i

  !     !======================================================
  subroutine MPI__send_iv(data,dest)
    implicit none
    integer:: data(:)
    integer :: n,dest,ierr
    n=size(data)
    call MPI_Send(data,n,MPI_INTEGER,dest,mpi__rank, MPI_COMM_WORLD,ierr)
  end subroutine MPI__send_iv

  !     !======================================================
  subroutine MPI__recv_iv(data,src)
    implicit none
    integer:: data(:)
    integer :: n,src,ierr
    n=size(data)
    call MPI_Recv(data,n,MPI_INTEGER,src,src, MPI_COMM_WORLD,ista,ierr)
  end subroutine MPI__recv_iv

  !     !======================================================
  subroutine MPI__REAL8send(data,n,dest)
    implicit none
    real(8):: data(n)
    integer :: n,dest,ierr
    call MPI_Send(data,n,MPI_REAL8,dest,mpi__rank, MPI_COMM_WORLD,ierr)
  end subroutine MPI__REAL8send

  !     !======================================================
  subroutine MPI__REAL8recv(data,n,src)
    implicit none
    real(8):: data(n)
    integer :: n,src,ierr
    call MPI_Recv(data,n,MPI_REAL8,src,src, MPI_COMM_WORLD,ista,ierr)
  end subroutine MPI__REAL8recv

  !     !======================================================
  subroutine MPI__DbleCOMPLEXsend(data,n,dest)
    implicit none
    complex(8):: data(n)
    integer :: n,dest,ierr
    call MPI_Send(data,n,MPI_COMPLEX16,dest,mpi__rank, MPI_COMM_WORLD,ierr)
  end subroutine MPI__DbleCOMPLEXsend

  !     !======================================================
  subroutine MPI__DbleCOMPLEXrecv(data,n,src)
    implicit none
    complex(8):: data(n)
    integer :: n,src,ierr
    call MPI_Recv(data,n,MPI_COMPLEX16,src,src, MPI_COMM_WORLD,ista,ierr)
  end subroutine MPI__DbleCOMPLEXrecv

  !     !======================================================
  subroutine MPI__DbleCOMPLEXsendQ(data,n,destQ)
    implicit none
    complex(8):: data(n)
    integer :: n,destQ,ierr
    call MPI_Send(data,n,MPI_COMPLEX16,destQ,0,mpi__commQ,ierr)
  end subroutine MPI__DbleCOMPLEXsendQ

  !     !======================================================
  subroutine MPI__DbleCOMPLEXrecvQ(data,n,srcQ)
    implicit none
    complex(8):: data(n)
    integer :: n,srcQ,ierr
    call MPI_Recv(data,n,MPI_COMPLEX16,srcQ,0,mpi__commQ,ista,ierr)
  end subroutine MPI__DbleCOMPLEXrecvQ

  !     !======================================================
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

    return
  end subroutine MPI__AllreduceSum

  !     !======================================================
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

  !     !======================================================
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

    return
  end subroutine MPI__AllreduceMax

  !     !======================================================
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

    return
  end subroutine MPI__sxcf_rankdivider

  !!
  !     !======================================================
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


  !     !======================================================
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

  !     !======================================================
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

  !     !======================================================
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
    !      
    !     ! setup Q-parallel
    !     !  commQ for among other Q, same S,P,B,M
    !      allocate( ranklistQ(0:mpi__sizeQ-1) )
    !      i=0
    !      do j=0, mpi__size-1
    !c     if(  vrankS(j) == mpi__rankS .and.
    !c     .        vrankP(j) == mpi__rankP .and.
    !c     .        vrankB(j) == mpi__rankB .and.
    !c     .        vrankM(j) == mpi__rankM ) then
    !         ranklistQ(i) = j
    !         i=i+1
    !c     end if
    !      end do
    !      call MPI_Comm_group( MPI_COMM_WORLD, mpi__group, mpi__info )
    !      call MPI_Group_incl(  mpi__group, mpi__sizeQ, ranklistQ, mpi__groupQ, mpi__info )
    !      call MPI_Comm_create( MPI_COMM_WORLD, mpi__groupQ, mpi__commQ, mpi__info ) 
    !      deallocate( ranklistQ )
    !      
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

