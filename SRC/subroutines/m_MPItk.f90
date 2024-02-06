!> MPI utility routines and variablis by TK
function convcchar(instr) result(outstr) !convert char(1) to char(1024)
   use iso_c_binding
   integer:: i,nend
   character(1):: instr(1024)
   character(1024):: outstr,instr2
   forall(i=1:1024) instr2(i:i) = instr(i)
   nend = index(instr2, c_null_char)-1
   outstr=''
   outstr(1:nend)= instr2(1:nend)
end function convcchar
module m_prgnam
   character(32):: prgnamx = ''
contains
   subroutine set_prgnam(prgnam)
      character(*):: prgnam
      prgnamx = prgnam
   end subroutine set_prgnam
   subroutine set_prgnamc(prgnamc) bind(C)
      character(1024):: convcchar
      character(1):: prgnamc(*)
      prgnamx = trim(convcchar(prgnamc))
      write (*, *) 'prgnamx=', trim(prgnamx)
    end subroutine set_prgnamc
end module m_prgnam

module m_MPItk
   use m_ext, only: sname
   use m_lgunit, only: stml, stdl
   public :: m_MPItk_init, m_MPItk_finalize, procid, strprocid, master, nsize, master_mpi, mlog, mlog_MPIiq, xmpbnd2
   private
   integer:: procid, master = 0, nsize
   include "mpif.h"
   logical:: mlog, master_mpi
   integer :: numprocs, ierr, status(MPI_STATUS_SIZE)
   character*(MPI_MAX_PROCESSOR_NAME) name
   character*10, allocatable:: shortname(:)
   character :: datim(24)
   character(8) :: strprocid
   double precision :: starttime, endtime
   integer ::  resultlen, id, nproc
   character::  prgnam*32, ext*100
contains
   subroutine m_MPItk_init()
      use m_prgnam, only: prgnamx
      !character::  prgnamx*(*)
!      integer:: mpipid
      logical:: cmdopt0
      integer :: fext
      character(10):: i2char
      prgnam = prgnamx
      !    call mpi_init(ierr)
      call mpi_comm_size(MPI_COMM_WORLD, nsize, ierr)
      call mpi_comm_rank(MPI_COMM_WORLD, id, ierr)
      allocate (shortname(0:nsize - 1))
      call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, nproc, ierr )
      shortname(procid) = trim(name)
      strprocid = trim(i2char(procid))
      call Gettime(datim)
      mlog = cmdopt0('--mlog') !! set log for --mlog (not maintained well)
      if (mlog) write (stml, "(a)") ' lmf '//datim//' Process ' &
         //trim(i2char(procid))//' of '//trim(i2char(nproc - 1))//' on '//trim(shortname(procid))
!      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if (procid == master) ext = ''
      if (procid /= master) ext = '_'//trim(i2char(procid))
      master_mpi = .false.
      if (procid == master) master_mpi = .TRUE.
   end subroutine m_MPItk_init
   subroutine m_MPItk_finalize()
      character(256) :: strn
      character(26):: datim
      character(20):: hostnm
      real(8):: cpusec
!      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if (master_mpi) then
         call fdate(datim)
         datim = datim(1:26) !takao. when datim(1:24) write(6,*) do not gives CR at the ene of datim*26.
         call ftime(datim)
         hostnm = ' '
         call get_environment_variable('HOST', hostnm)
         write (stdl, 10) cpusec(), trim(datim), trim(adjustl(hostnm))
10       format(' CPU time:', f9.3, ' sec  ', a26, ' on ', a)
         close (stdl)
      end if
      call Mpi_finalize(ierr)
   end subroutine m_MPItk_finalize

   subroutine mlog_MPIiq(iq, iqini, iqend)
      integer:: iq, numprocs, iqini, iqend
      character(26) :: datim
      character(512):: aaachar
      character(10):: i2char
      numprocs = nsize
      if (mlog) then
         call gettime(datim)
         aaachar = ' bndfp '//datim//' Process '// &
                   trim(i2char(procid))//' of '//trim(i2char(numprocs))//' on '// &
                   trim(i2char(procid))//' starting k-points '// &
                   trim(i2char(iqini))//' to '//trim(i2char(iqend))
         write (stml, "(a)") trim(aaachar)
      end if
   end subroutine mlog_MPIiq

   ! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
   subroutine xmpbnd2(kpproc, ndham, ndat, eb)
      !- Collect eb from various processors (MPI)
      ! ----------------------------------------------------------------------
      !i Inputs
      !i   kpproc
      !i   ndham :leading dimension of eb
      !i   nkp   :number of irreducible k-points (bzmesh.f)
      !i   nsp   :2 for spin-polarized case, otherwise 1
      !i   eb    :energy bands; alias eband
      !o Outputs
      ! ----------------------------------------------------------------------
      implicit none
      integer:: kpproc(0:*), ndham, ndat
      double precision :: eb(ndham, ndat)
      integer :: i, ista, iend, ierr
      integer, dimension(:), allocatable :: offset, length
      real(8), allocatable :: buf_rv(:, :)
      allocate (offset(0:nsize), stat=ierr)
      allocate (length(0:nsize), stat=ierr)
      offset(0) = 0
      do i = 0, nsize - 1 !NOTE: ndat is divided into kpproc
         ista = kpproc(i)
         iend = kpproc(i + 1) - 1
         length(i) = (iend - ista + 1)*ndham !nsp*ndham
         offset(i + 1) = offset(i) + length(i)
      end do
      ista = kpproc(procid)
      allocate (buf_rv(ndham, ndat))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !    write(6,*)'ppppp',procid,ista,ndham,'nkp nspx=',nkp,nspx,'len=',length(procid)
      call mpi_allgatherv(eb(1, ista), length(procid), mpi_double_precision, buf_rv, length, offset &
                          , mpi_double_precision, mpi_comm_world, ierr)
      eb = buf_rv
      deallocate (buf_rv)
      deallocate (offset, stat=ierr)
      deallocate (length, stat=ierr)
   end subroutine xmpbnd2
end module m_MPItk
