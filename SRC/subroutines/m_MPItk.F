      module m_MPItk
      use m_ext,only:sname
      public ::
     &     m_MPItk_init, m_MPItk_finalize, procid, strprocid,master,nsize,master_mpi,mlog,mlog_MPIiq

      private
      integer:: procid,master=0,nsize
      include "mpif.h"
      logical:: mlog,master_mpi
      integer numprocs, ierr, status(MPI_STATUS_SIZE)
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10,allocatable:: shortname(:) 
      character*26 datim
      character*8 strprocid
      double precision starttime, endtime
      integer  resultlen, id,nproc
      character::  prgnam*32, ext*100
      contains
      
      subroutine m_MPItk_init(prgnamx)
#if MPIK      
      character::  prgnamx*(*)
      integer:: mpipid
      logical:: cmdopt0
      character(5) :: charnum5
      character*8 charext
      integer fext,lgunit
      prgnam=prgnamx
      call mpi_init(ierr)
      call mpi_comm_size(MPI_COMM_WORLD, nsize,ierr)
      call mpi_comm_rank(MPI_COMM_WORLD, id,ierr)
      allocate(shortname(0:nsize-1))
      call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
      procid = mpipid(1) !
      nproc  = mpipid(0) ! num of processors
      shortname(procid) = trim(name)
      strprocid=charext(procid)
      call Gettime(datim)
      call Finits() !read and set arguments addsyv in symvar.F
      mlog = cmdopt0('--mlog') !! set log for --mlog (not maintained well)
      if(mlog) write(lgunit(3),"(a)")' lmf '//datim//' Process '
     &  //trim(charnum5(procid))//' of '//trim(charnum5(nproc-1))//' on '//trim(shortname(procid)) 
      call MPI_BARRIER( MPI_COMM_WORLD, ierr )
#else
      call Finits() !read and set arguments addsyv in symvar.F
      procid=0 
#endif      
      if(procid == master) ext=''
      if(procid /= master) ext='_'//trim(charext(procid))
      ierr = fext('.'//trim(sname)//trim(ext)) !set ext for fopen
      master_mpi=.false.
      if(procid==master) master_mpi=.true.
      end

      subroutine m_MPItk_finalize()
      character*256 strn, datim*26, hostnm*20
      integer:: lgunit,stdl
      real(8):: cpusec
#if MPIK
      call MPI_BARRIER( MPI_COMM_WORLD, ierr )
c     if( master_mpi ) call rx0(prgnam//' on '//trim(shortname(procid)))
      if(master_mpi) then
         call fdate(datim)
         datim=datim(1:24)      !takao. If this is not, write(6,*) gives CR at the ene of datim*26.
         call ftime(datim)
         hostnm = ' '
         call get_environment_variable('HOST',hostnm)
         stdl = lgunit(2)
         write(stdl,10) cpusec(), trim(datim),trim(adjustl(hostnm))
 10      format(' CPU time:', f9.3,' sec  ',a26,' on ',a)
         close(stdl)
      endif
      call Mpi_finalize(ierr)
#else
      close(stdl)
      call Rx0s(prgnam) !for single core exit
#endif
      end
      
      subroutine mlog_MPIiq(iq,iqini,iqend)
      integer:: iq,numprocs,lgunit,iqini,iqend
      character*26 datim
      character(512):: aaachar
      character(10):: i2char
      character*8:: charext
      numprocs=nsize
c      if(iq==iqini) then
         if(mlog) then
               call gettime(datim)
               aaachar=' bndfp '//datim//' Process '//
     &              trim(i2char(procid))//' of '//trim(i2char(numprocs))//' on '//
     &              trim(charext(procid))//' starting k-points '//
     &              trim(i2char( iqini))//' to '//trim(i2char(iqend))
               write(lgunit(3),"(a)")trim(aaachar)
         endif
c      endif   
      end
      
      end module
