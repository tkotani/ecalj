!> MPI all reduce
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
subroutine mpibc1(vec,n,cast,mlog,funnam,label)  !- Broadcasts a vector from master node to the world (MPI)
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
  implicit none
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
  integer :: procid,master
  logical :: mlog
  integer :: n,cast
  double precision :: vec(n)
  character funnam*(*), label*(*)
  if (n <= 0) return
  master = 0
  call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
  if (cast == 1) then
     call MPI_BCAST(vec,n,MPI_LOGICAL,   master,MPI_COMM_WORLD,ierr)
  elseif (cast == 2) then
     call MPI_BCAST(vec,n,MPI_INTEGER,   master,MPI_COMM_WORLD,ierr)
  elseif (cast == 4) then
     call MPI_BCAST(vec,n,MPI_DOUBLE_PRECISION, master,MPI_COMM_WORLD,ierr)
  elseif (cast == 6) then
     call MPI_BCAST(vec,2*n,MPI_DOUBLE_PRECISION, master,MPI_COMM_WORLD,ierr)
  else
     call rxi('mpibc1: cast not implemented',cast)
  endif
  if (mlog) then
     call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
     shortname(procid)=name
     namelen(procid) = ierr-1
     call gettime(datim)
  endif
end subroutine mpibc1
subroutine mpibc2(vec,n,cast,mlog,funnam,label) !Performs MPI_ALLREDUCE on a vector (MPI)
  use m_lgunit,only:stml
  !i   vec   :vector to broadcast
  !i   n     :length of vector
  !i   cast  :cast of vector:
  !i         : 2 int
  !i         : 4 double
  !i         : 6 double complex
  !i   mlog  : T write message to mlog file
  !i   funnam:string used in writing message (function name)
  !i   label :string used in writing message (variable name)
  !r Remarks
  !r   ALLREDUCE sums the contributions from all the individual threads
  implicit none
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
  integer :: procid,master
  logical :: mlog
  integer :: n,cast
  double precision :: vec(n)
  character funnam*(*), label*(*)
  integer, allocatable :: ibuf(:)
  real(8) ,allocatable :: dbuf(:)
  integer :: obuf
  if (n <= 0) return
  master = 0
  call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
  call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
  if (cast == 2) then
     allocate(ibuf(n), stat=ierr)
     call MPI_ALLREDUCE(vec,ibuf,n,  MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
     call icopy(n,ibuf,1,vec,1)
     deallocate(ibuf, stat=ierr)
  elseif (cast == 4) then
     allocate(dbuf(n), stat=ierr)
     call MPI_ALLREDUCE(vec,dbuf,n,  MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
     call dcopy(n,dbuf,1,vec,1)
     deallocate(dbuf, stat=ierr)
  elseif (cast == 6) then
     allocate(dbuf(2*n), stat=ierr)
     call MPI_ALLREDUCE(vec,dbuf,2*n,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
     call dcopy(2*n,dbuf,1,vec,1)
     deallocate(dbuf, stat=ierr)
  else
     call rxi('mpibc2: cast not implemented',cast)
  endif
  if (mlog) then
     call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
     shortname(procid)=name
     namelen(procid) = ierr-1
     call gettime(datim)
     write(stml,"(a,i0,' of ',i0,' on ',a)") &
          ' '//funnam//' '//datim//' Process ',procid,numprocs, &
          shortname(procid)(1:namelen(procid))//' allreduce '//label
  endif
end subroutine mpibc2
subroutine icopy(n,dx,incx,dy,incy)   !     copies a vector, x, to a vector, y.  Adapted from:
  integer :: dx(1),dy(1)
  integer :: i,incx,incy,ix,iy,n
  ix = 1
  iy = 1
  if (incx < 0) ix = (1-n)*incx + 1
  if (incy < 0) iy = (1-n)*incy + 1
  do  10  i = 1, n
     dy(iy) = dx(ix)
     ix = ix + incx
     iy = iy + incy
10 enddo
end subroutine icopy
integer function mpipid(mode)  !- Returns MPI procid
  !i   mode  :0 return number of processors
  !i         :1 return procid
  !i         :2 calls MPI_BARRIER; returns ierr
  !i         :Otherwise, return 0
  !o   mpipid:procid or number of processors (see mode)
  implicit none
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
