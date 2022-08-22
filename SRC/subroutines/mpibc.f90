! all reduce
subroutine mpibc2_int(vec,nnn,label)
  use m_MPItk,only: mlog
  character funnam*(1), label*(*)
  integer::nnn,cast
  integer:: vec(*)
  !      print *,trim(label),nnn,vec(1)
  cast=2
  funnam=''
  call mpibc2(vec,nnn,cast,mlog,funnam,label)
end subroutine mpibc2_int
subroutine mpibc2_real(vec,nnn,label)
  use m_MPItk,only: mlog
  character funnam*(1), label*(*)
  integer::nnn,cast
  real(8):: vec(*)
  !      print *,trim(label),nnn,vec(1)
  cast=4
  funnam=''
  call mpibc2(vec,nnn,cast,mlog,funnam,label)
end subroutine mpibc2_real
subroutine mpibc2_complex(vec,nnn,label)
  use m_MPItk,only: mlog
  character funnam*(1), label*(*)
  integer::nnn,cast
  complex(8):: vec(*)
  !      print *,trim(label),nnn
  cast=6
  funnam=''
  call mpibc2(vec,nnn,cast,mlog,funnam,label)
end subroutine mpibc2_complex
! master to world
subroutine mpibc1_logical(vec,nnn,label)
  use m_MPItk,only: mlog
  character funnam*(1), label*(*)
  integer::nnn,cast
  real(8):: vec(*)
  cast=1
  funnam=''
  call mpibc1(vec,nnn,cast,mlog,funnam,label)
end subroutine mpibc1_logical
subroutine mpibc1_int(vec,nnn,label)
  use m_MPItk,only: mlog
  character funnam*(1), label*(*)
  integer::nnn,cast
  real(8):: vec(*)
  cast=2
  funnam=''
  call mpibc1(vec,nnn,cast,mlog,funnam,label)
end subroutine mpibc1_int
subroutine mpibc1_real(vec,nnn,label)
  use m_MPItk,only: mlog
  character funnam*(1), label*(*)
  integer::nnn,cast
  real(8):: vec(*)
  cast=4
  funnam=''
  call mpibc1(vec,nnn,cast,mlog,funnam,label)
end subroutine mpibc1_real
subroutine mpibc1_complex(vec,nnn,label)
  use m_MPItk,only: mlog
  character funnam*(1), label*(*)
  integer::nnn,cast
  complex(8):: vec(*)
  !      print *,trim(label),nnn
  cast=6
  funnam=''
  call mpibc1(vec,nnn,cast,mlog,funnam,label)
end subroutine mpibc1_complex


subroutine mpibc1(vec,n,cast,mlog,funnam,label)
  !- Broadcasts a vector from master node to the world (MPI)
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   vec   :vector to broadcast
  !i   n     :length of vector
  !i   cast  :cast of vector:
  !i         : 1 logical
  !i         : 2 int
  !i         : 4 double
  !i         : 6 double complex
  !i   mlog  : T write message to mlog file
  !i   funnam:string used in writing message (function name)
  !i   label :string used in writing message (variable name)
  !o Outputs
  !l Local variables
  !r Remarks
  !r
  !u Updates
  !u   09 Jul 07 Can broadcast logical vectors
  !u   14 Apr 03 First created
  ! ----------------------------------------------------------------------
  !     implicit none
  include "mpif.h"
  integer :: numprocs, ierr
  integer :: MAX_PROCS
  parameter (MAX_PROCS = 100)
  integer :: resultlen
  character*(MPI_MAX_PROCESSOR_NAME) name
  character(10) :: shortname(0:MAX_PROCS-1)
  character(26) :: datim
  integer :: namelen(0:MAX_PROCS-1)
  character(256) :: strn
  !      logical lgunit
  integer :: procid,master
  ! ... Passed parameters
  logical :: mlog
  integer :: n,cast
  double precision :: vec(n)
  character funnam*(*), label*(*)
  ! ... Local parameters
  if (n <= 0) return
  master = 0
  call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )

  if (cast == 1) then
     call MPI_BCAST(vec,n,MPI_LOGICAL, &
          master,MPI_COMM_WORLD,ierr)
  elseif (cast == 2) then
     call MPI_BCAST(vec,n,MPI_INTEGER, &
          master,MPI_COMM_WORLD,ierr)

  elseif (cast == 4) then
     call MPI_BCAST(vec,n,MPI_DOUBLE_PRECISION, &
          master,MPI_COMM_WORLD,ierr)

  elseif (cast == 6) then
     call MPI_BCAST(vec,2*n,MPI_DOUBLE_PRECISION, &
          master,MPI_COMM_WORLD,ierr)

  else
     call rxi('mpibc1: cast not implemented',cast)

  endif

  if (mlog) then
     call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
     !call strcop(shortname(procid),name,10,'.',ierr)
     shortname(procid)=name
     namelen(procid) = ierr-1
     call gettime(datim)
     !        strn = ' '//funnam//' '//datim//' Process %i of %i on '
     !     .  //shortname(procid)(1:namelen(procid))//' bcast '//label//
     !     .  ' (%i %?#n==2#int##%?#n==4#d.p.##%?#n==6#d.c.##)'
     !        call awrit6(strn,' ',-256,lgunit(3),procid,numprocs,n,cast,cast,
     !     .  cast)
  endif

end subroutine mpibc1
!$$$      subroutine mpibc3(vec,n,cast,pid,mlog,funnam,label) !   this is removed at 19jun2021

subroutine mpibc2(vec,n,cast,mlog,funnam,label)
  use m_lgunit,only:stml
  !- Performs MPI_ALLREDUCE on a vector (MPI)
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   vec   :vector to broadcast
  !i   n     :length of vector
  !i   cast  :cast of vector:
  !i         : 2 int
  !i         : 4 double
  !i         : 6 double complex
  !i   mlog  : T write message to mlog file
  !i   funnam:string used in writing message (function name)
  !i   label :string used in writing message (variable name)
  !o Outputs
  !l Local variables
  !r Remarks
  !r   ALLREDUCE sums the contributions from all the individual threads
  !u Updates
  !u   14 Apr 03 First created
  ! ----------------------------------------------------------------------
  !     implicit none
  include "mpif.h"
  integer :: numprocs, ierr
  integer :: MAX_PROCS
  parameter (MAX_PROCS = 100)
  integer :: resultlen
  character*(MPI_MAX_PROCESSOR_NAME) name
  character(10) :: shortname(0:MAX_PROCS-1)
  character(26) :: datim
  integer :: namelen(0:MAX_PROCS-1)
  character(256) :: strn
  !      integer:: lgunit
  integer :: procid,master
  ! ... Passed parameters
  logical :: mlog
  integer :: n,cast
  double precision :: vec(n)
  character funnam*(*), label*(*)
  ! ... Local parameters
  integer, allocatable :: ibuf(:)
  real(8) ,allocatable :: dbuf(:)
  integer :: obuf

  if (n <= 0) return
  master = 0
  call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )

  if (cast == 2) then
     allocate(ibuf(n), stat=ierr)
     call MPI_ALLREDUCE(vec,ibuf,n, &
          MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
     call icopy(n,ibuf,1,vec,1)
     deallocate(ibuf, stat=ierr)
  elseif (cast == 4) then
     allocate(dbuf(n), stat=ierr)
     call MPI_ALLREDUCE(vec,dbuf,n, &
          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
     call dcopy(n,dbuf,1,vec,1)
     deallocate(dbuf, stat=ierr)
  elseif (cast == 6) then
     allocate(dbuf(2*n), stat=ierr)
     call MPI_ALLREDUCE(vec,dbuf,2*n, &
          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
     call dcopy(2*n,dbuf,1,vec,1)
     deallocate(dbuf, stat=ierr)
  else
     call rxi('mpibc2: cast not implemented',cast)
  endif

  if (mlog) then
     call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
     !call strcop(shortname(procid),name,10,'.',ierr)
     shortname(procid)=name
     namelen(procid) = ierr-1
     call gettime(datim)
     !        call awrit2(strn,' ',-256,lgunit(3),procid,numprocs)
     write(stml,"(a,i0,' of ',i0,' on ',a)") &
          ' '//funnam//' '//datim//' Process ',procid,numprocs, &
          shortname(procid)(1:namelen(procid))//' allreduce '//label
  endif
end subroutine mpibc2

integer function mpipid(mode)
  !- Returns MPI procid
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   mode  :0 return number of processors
  !i         :1 return procid
  !i         :2 calls MPI_BARRIER; returns ierr
  !i         :Otherwise, return 0
  !o Outputs
  !o   mpipid:procid or number of processors (see mode)
  !r Remarks
  !u Updates
  !u   24 Nov 05 Added mode 2
  !u   14 Apr 03 First created
  ! ----------------------------------------------------------------------
  !     implicit none
  include "mpif.h"
  integer :: numprocs, ierr, procid
  integer :: mode

  if (mode == 0) then
     call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
     mpipid = numprocs
  else if (mode == 1) then
     call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
     mpipid = procid
  else if (mode == 2) then
     call MPI_BARRIER( MPI_COMM_WORLD, ierr )
     mpipid = ierr
  else
     mpipid = 0
  endif
end function mpipid
!      subroutine fmain
!      integer n
!      n = mpipid(0)
!      print *, 'mpipid for number of processors:', n
!      n = mpipid(1)
!      print *, 'mpipid for processor id:', n
!      call MPI_FINALIZE(n)
!      end

