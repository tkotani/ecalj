module m_MPItk
  use m_ext,only:sname
  use m_lgunit,only: stml,stdl
  public :: &
       m_MPItk_init, m_MPItk_finalize, procid, strprocid,master,nsize,master_mpi,mlog,mlog_MPIiq,&
       xmpbnd2

  private
  integer:: procid,master=0,nsize
  include "mpif.h"
  logical:: mlog,master_mpi
  integer :: numprocs, ierr, status(MPI_STATUS_SIZE)
  character*(MPI_MAX_PROCESSOR_NAME) name
  character*10,allocatable:: shortname(:)
  character :: datim(24)
  character(8) :: strprocid
  double precision :: starttime, endtime
  integer ::  resultlen, id,nproc
  character::  prgnam*32, ext*100
contains

  subroutine m_MPItk_init(prgnamx)
    character::  prgnamx*(*)
    integer:: mpipid
    logical:: cmdopt0
    integer :: fext
    character(10):: i2char
    prgnam=prgnamx
    call mpi_init(ierr)
    call mpi_comm_size(MPI_COMM_WORLD, nsize,ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, id,ierr)
    allocate(shortname(0:nsize-1))
    call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
    procid = mpipid(1) !
    nproc  = mpipid(0) ! num of processors
    shortname(procid) = trim(name)
    strprocid=trim(i2char(procid))
    call Gettime(datim)
    mlog = cmdopt0('--mlog') !! set log for --mlog (not maintained well)
    if(mlog) write(stml,"(a)")' lmf '//datim//' Process ' &
         //trim(i2char(procid))//' of '//trim(i2char(nproc-1))//' on '//trim(shortname(procid))
    call MPI_BARRIER( MPI_COMM_WORLD, ierr )
    if(procid == master) ext=''
    if(procid /= master) ext='_'//trim(i2char(procid))
    master_mpi=.false.
    if(procid==master) master_mpi= .TRUE. 
  end subroutine m_MPItk_init

  subroutine m_MPItk_finalize()
    character(256) :: strn 
    character(26):: datim
    character(20):: hostnm
    real(8):: cpusec
    call MPI_BARRIER( MPI_COMM_WORLD, ierr )
    if(master_mpi) then
       call fdate(datim)
       datim=datim(1:26) !takao. when datim(1:24) write(6,*) do not gives CR at the ene of datim*26.
       call ftime(datim)
       hostnm = ' '
       call get_environment_variable('HOST',hostnm)
       write(stdl,10) cpusec(), trim(datim),trim(adjustl(hostnm))
10     format(' CPU time:', f9.3,' sec  ',a26,' on ',a)
       close(stdl)
    endif
    call Mpi_finalize(ierr)
  end subroutine m_MPItk_finalize

  subroutine mlog_MPIiq(iq,iqini,iqend)
    integer:: iq,numprocs,iqini,iqend
    character(26) :: datim
    character(512):: aaachar
    character(10):: i2char
    numprocs=nsize
    if(mlog) then
       call gettime(datim)
       aaachar=' bndfp '//datim//' Process '// &
            trim(i2char(procid))//' of '//trim(i2char(numprocs))//' on '// &
            trim(i2char(procid))//' starting k-points '// &
            trim(i2char( iqini))//' to '//trim(i2char(iqend))
       write(stml,"(a)")trim(aaachar)
    endif
  end subroutine mlog_MPIiq

  ! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  subroutine xmpbnd2(kpproc,ndham,nkp,nspx,eb) !use m_MPItk,only: mlog, master_mpi, strprocid, numprocs=>nsize, mlog_MPIiq
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
    integer:: kpproc(0:*),ndham,nkp,nspx
    double precision :: eb(ndham,nspx*nkp)
    integer :: i,ista,iend, ierr
    integer, dimension(:),allocatable :: offset,length
    real(8) ,allocatable :: buf_rv(:,:)
    allocate (offset(0:nsize), stat=ierr)
    allocate (length(0:nsize), stat=ierr)
    offset(0) = 0
    do  i = 0, nsize-1
       ista = kpproc(i)
       iend = kpproc(i+1)-1
       length(i) = (iend - ista + 1)*ndham !nsp*ndham
       offset(i+1) = offset(i) + length(i)
    enddo
    ista = kpproc(procid)
    allocate(buf_rv(ndham,nkp*nspx))
    call mpi_allgatherv(eb(1,ista),length(procid),mpi_double_precision,buf_rv,length,offset &
         , mpi_double_precision , mpi_comm_world , ierr )
    eb= buf_rv 
    deallocate(buf_rv)
    deallocate(offset, stat=ierr)
    deallocate(length, stat=ierr)
  end subroutine xmpbnd2

end module m_MPItk
