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
   use m_mpi,only: procid, numprocs=>nsize,comm
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
  use mpi
  implicit none
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
subroutine mpibc2(vec,n,cast,mlog,funnam,label) !MPI_ALLREDUCE SUM on `vec`, in-place on every rank
  use m_mpi, only: comm
  use mpi
  implicit none
  !i   vec   : data to be summed across all ranks; result is delivered
  !i           back into `vec` on every rank (MPI_IN_PLACE)
  !i   n     : length of vec
  !i   cast  : 2 = integer, 4 = double, 6 = double complex
  !i   mlog, funnam, label : retained for caller-API compatibility (unused)
  logical :: mlog
  integer :: n, cast, ierr
  double precision :: vec(n)
  character funnam*(*), label*(*)
  if (n <= 0) return
  select case (cast)
  case (2); call MPI_ALLREDUCE(MPI_IN_PLACE, vec, n,   MPI_INTEGER,          MPI_SUM, comm, ierr)
  case (4); call MPI_ALLREDUCE(MPI_IN_PLACE, vec, n,   MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
  case (6); call MPI_ALLREDUCE(MPI_IN_PLACE, vec, 2*n, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
  case default; call rxi('mpibc2: cast not implemented', cast)
  end select
end subroutine mpibc2
