!> MPI all reduce
subroutine mpibc2_int(vec,nnn,label)
  logical:: mlog
  character funnam*(1), label*(*)
  integer::nnn,cast
  integer:: vec(nnn)
  cast=2
  funnam=''
  call mpibc2(vec,nnn,cast,mlog,funnam,label)
end subroutine mpibc2_int
subroutine mpibc2_real(vec,nnn,label)
  logical:: mlog
  character funnam*(1), label*(*)
  integer::nnn,cast
  real(8):: vec(nnn)
  cast=4
  funnam=''
  call mpibc2(vec,nnn,cast,mlog,funnam,label)
end subroutine mpibc2_real
subroutine mpibc2_complex(vec,nnn,label)
  logical:: mlog
  character funnam*(1), label*(*)
  integer::nnn,cast
  complex(8):: vec(nnn)
  cast=6
  funnam=''
  call mpibc2(vec,nnn,cast,mlog,funnam,label)
end subroutine mpibc2_complex
subroutine mpibc1_logical(vec,nnn,label)
  logical:: mlog
  character funnam*(1), label*(*)
  integer::nnn,cast
  real(8):: vec(nnn)
  cast=1
  funnam=''
  call mpibc1(vec,nnn,cast,mlog,funnam,label)
end subroutine mpibc1_logical
subroutine mpibc1_int(vec,nnn,label)
  logical:: mlog
  character funnam*(1), label*(*)
  integer::nnn,cast
  real(8):: vec(nnn)
  cast=2
  funnam=''
  call mpibc1(vec,nnn,cast,mlog,funnam,label)
end subroutine mpibc1_int
subroutine mpibc1_real(vec,nnn,label)
  logical:: mlog
  character funnam*(1), label*(*)
  integer::nnn,cast
  real(8):: vec(nnn)
  cast=4
  funnam=''
  call mpibc1(vec,nnn,cast,mlog,funnam,label)
end subroutine mpibc1_real
subroutine mpibc1_complex(vec,nnn,label)
  logical:: mlog
  character funnam*(1), label*(*)
  integer::nnn,cast
  complex(8):: vec(nnn)
  cast=6
  funnam=''
  call mpibc1(vec,nnn,cast,mlog,funnam,label)
end subroutine mpibc1_complex
subroutine mpibc1(vec,n,cast,mlog,funnam,label)  !- Broadcasts a vector from master node to the world (MPI)
   use m_MPItk,only: procid, numprocs=>nsize,comm
  !i Inputs
  !i   vec   :vector to broadcast
  !i   n     :length of vector
  !i   cast  :cast of vector:
  !i         : 1 logical
  !i         : 2 int
  !i         : 4 double
  !i         : 6 double complex
  !i   mlog  : dummy
  !i   funnam:string used in writing message (function name)
  !i   label :string used in writing message (variable name)
  implicit none
  include "mpif.h"
  integer :: ierr
  integer :: MAX_PROCS
  parameter (MAX_PROCS = 100)
  integer :: resultlen
  character*(MPI_MAX_PROCESSOR_NAME) name
  character(10) :: shortname(0:MAX_PROCS-1)
  character(26) :: datim
  integer :: namelen(0:MAX_PROCS-1)
  character(256) :: strn
  integer :: master
  logical :: mlog
  integer :: n,cast
  double precision :: vec(n)
  character funnam*(*), label*(*)
  if (n <= 0) return
  master = 0
  if (cast == 1) then
     call MPI_BCAST(vec,n,MPI_LOGICAL,   master,comm,ierr)
  elseif (cast == 2) then
     call MPI_BCAST(vec,n,MPI_INTEGER,   master,comm,ierr)
  elseif (cast == 4) then
     call MPI_BCAST(vec,n,MPI_DOUBLE_PRECISION, master,comm,ierr)
  elseif (cast == 6) then
     call MPI_BCAST(vec,2*n,MPI_DOUBLE_PRECISION, master,comm,ierr)
  else
     call rxi('mpibc1: cast not implemented',cast)
  endif
end subroutine mpibc1
subroutine mpibc2(vec,n,cast,mlog,funnam,label) !Performs MPI_ALLREDUCE on a vector (MPI)
  use m_MPItk, only: procid, numprocs=>nsize,comm
  use m_lgunit,only:stml
  use m_ftox
  !i   vec   :vector to broadcast
  !i   n     :length of vector
  !i   cast  :cast of vector:
  !i         : 2 int
  !i         : 4 double
  !i         : 6 double complex
  !i   mlog  : dummy 
  !i   funnam:string used in writing message (function name)
  !i   label :string used in writing message (variable name)
  !r Remarks
  !r   ALLREDUCE sums the contributions from all the individual threads
  implicit none
  include "mpif.h"
  integer ::ierr
  integer :: MAX_PROCS
  parameter (MAX_PROCS = 100)
  integer :: resultlen
  character*(MPI_MAX_PROCESSOR_NAME) name
  character(10) :: shortname(0:MAX_PROCS-1)
  character(26) :: datim
  integer :: namelen(0:MAX_PROCS-1)
  character(256) :: strn
  integer :: master
  logical :: mlog
  integer :: n,cast
  double precision :: vec(n)
  character funnam*(*), label*(*)
  integer, allocatable :: ibuf(:)
  real(8) ,allocatable :: dbuf(:)
  integer :: obuf
  if (n <= 0) return
  master = 0
  if (cast == 2) then
     allocate(ibuf(n), stat=ierr)
     call MPI_ALLREDUCE(vec,ibuf,n,  MPI_INTEGER,MPI_SUM,comm,ierr)
     call icopy(n,ibuf,1,vec,1)     !vec=transfer(ibuf,vec) !this did not work in ifort ver2018 in ucgw
     deallocate(ibuf, stat=ierr)
  elseif (cast == 4) then
     allocate(dbuf(n), stat=ierr) 
     call MPI_ALLREDUCE(vec,dbuf,n,  MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
     !vec=transfer(dbuf,vec) 
     call dcopy(n,dbuf,1,vec,1)
     deallocate(dbuf, stat=ierr)
  elseif (cast == 6) then
     allocate(dbuf(2*n), stat=ierr)
     call MPI_ALLREDUCE(vec,dbuf,2*n,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
     call dcopy(2*n,dbuf,1,vec,1)
     !vec=transfer(dbuf,vec)
     deallocate(dbuf, stat=ierr)
  else
     call rxi('mpibc2: cast not implemented',cast)
  endif
end subroutine mpibc2
subroutine icopy(n,dx,incx,dy,incy)
  !     copies a vector, x, to a vector, y.  Adapted from:
  !     jack dongarra, linpack, 3/11/78.
  integer :: dx(*),dy(*)
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
