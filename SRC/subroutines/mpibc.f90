! all reduce
      subroutine mpibc2_int(vec,nnn,label)
      use m_MPItk,only: mlog
      character funnam*(1), label*(*)
      integer::nnn,cast
      integer:: vec(*) 
c      print *,trim(label),nnn,vec(1)
      cast=2
      funnam=''
      call mpibc2(vec,nnn,cast,mlog,funnam,label)
      end
      subroutine mpibc2_real(vec,nnn,label)
      use m_MPItk,only: mlog
      character funnam*(1), label*(*)
      integer::nnn,cast
      real(8):: vec(*) 
c      print *,trim(label),nnn,vec(1)
      cast=4
      funnam=''
      call mpibc2(vec,nnn,cast,mlog,funnam,label)
      end
      subroutine mpibc2_complex(vec,nnn,label)
      use m_MPItk,only: mlog
      character funnam*(1), label*(*)
      integer::nnn,cast
      complex(8):: vec(*) 
c      print *,trim(label),nnn
      cast=6
      funnam=''
      call mpibc2(vec,nnn,cast,mlog,funnam,label)
      end
! master to world      
      subroutine mpibc1_logical(vec,nnn,label)
      use m_MPItk,only: mlog
      character funnam*(1), label*(*)
      integer::nnn,cast
      real(8):: vec(*) 
      cast=1
      funnam=''
      call mpibc1(vec,nnn,cast,mlog,funnam,label)
      end
      subroutine mpibc1_int(vec,nnn,label)
      use m_MPItk,only: mlog
      character funnam*(1), label*(*)
      integer::nnn,cast
      real(8):: vec(*) 
      cast=2
      funnam=''
      call mpibc1(vec,nnn,cast,mlog,funnam,label)
      end
      subroutine mpibc1_real(vec,nnn,label)
      use m_MPItk,only: mlog
      character funnam*(1), label*(*)
      integer::nnn,cast
      real(8):: vec(*) 
      cast=4
      funnam=''
      call mpibc1(vec,nnn,cast,mlog,funnam,label)
      end
      subroutine mpibc1_complex(vec,nnn,label)
      use m_MPItk,only: mlog
      character funnam*(1), label*(*)
      integer::nnn,cast
      complex(8):: vec(*) 
c      print *,trim(label),nnn
      cast=6
      funnam=''
      call mpibc1(vec,nnn,cast,mlog,funnam,label)
      end

      
      subroutine mpibc1(vec,n,cast,mlog,funnam,label)
C- Broadcasts a vector from master node to the world (MPI)
C ----------------------------------------------------------------------
Ci Inputs
Ci   vec   :vector to broadcast
Ci   n     :length of vector
Ci   cast  :cast of vector:
Ci         : 1 logical
Ci         : 2 int
Ci         : 4 double
Ci         : 6 double complex
Ci   mlog  : T write message to mlog file
Ci   funnam:string used in writing message (function name)
Ci   label :string used in writing message (variable name)
Co Outputs
Cl Local variables
Cr Remarks
Cr
Cu Updates
Cu   09 Jul 07 Can broadcast logical vectors
Cu   14 Apr 03 First created
C ----------------------------------------------------------------------
C     implicit none
#if MPI|MPIK
      include "mpif.h"
      integer numprocs, ierr
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
      character*256 strn
c      logical lgunit
      integer procid,master
#endif
C ... Passed parameters
      logical mlog
      integer n,cast
      double precision vec(n)
      character funnam*(*), label*(*)
C ... Local parameters

#if MPI|MPIK
      if (n .le. 0) return
      master = 0
      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )

      if (cast .eq. 1) then
        call MPI_BCAST(vec,n,MPI_LOGICAL,
     .  master,MPI_COMM_WORLD,ierr)
      elseif (cast .eq. 2) then
        call MPI_BCAST(vec,n,MPI_INTEGER,
     .  master,MPI_COMM_WORLD,ierr)

      elseif (cast .eq. 4) then
        call MPI_BCAST(vec,n,MPI_DOUBLE_PRECISION,
     .  master,MPI_COMM_WORLD,ierr)

      elseif (cast .eq. 6) then
        call MPI_BCAST(vec,2*n,MPI_DOUBLE_PRECISION,
     .  master,MPI_COMM_WORLD,ierr)

      else
        call rxi('mpibc1: cast not implemented',cast)

      endif

      if (mlog) then
        call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
        call strcop(shortname(procid),name,10,'.',ierr)
        namelen(procid) = ierr-1
        call gettime(datim)
C        strn = ' '//funnam//' '//datim//' Process %i of %i on '
C     .  //shortname(procid)(1:namelen(procid))//' bcast '//label//
C     .  ' (%i %?#n==2#int##%?#n==4#d.p.##%?#n==6#d.c.##)'
C        call awrit6(strn,' ',-256,lgunit(3),procid,numprocs,n,cast,cast,
C     .  cast)
      endif
#endif

      end
c$$$      subroutine mpibc3(vec,n,cast,pid,mlog,funnam,label) !   this is removed at 19jun2021
      
      subroutine mpibc2(vec,n,cast,mlog,funnam,label)
      use m_lgunit,only:stml
C- Performs MPI_ALLREDUCE on a vector (MPI)
C ----------------------------------------------------------------------
Ci Inputs
Ci   vec   :vector to broadcast
Ci   n     :length of vector
Ci   cast  :cast of vector:
Ci         : 2 int
Ci         : 4 double
Ci         : 6 double complex
Ci   mlog  : T write message to mlog file
Ci   funnam:string used in writing message (function name)
Ci   label :string used in writing message (variable name)
Co Outputs
Cl Local variables
Cr Remarks
Cr   ALLREDUCE sums the contributions from all the individual threads
Cu Updates
Cu   14 Apr 03 First created
C ----------------------------------------------------------------------
C     implicit none
#if MPI|MPIK
      include "mpif.h"
      integer numprocs, ierr
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
      character*256 strn
c      integer:: lgunit
      integer procid,master
#endif
C ... Passed parameters
      logical mlog
      integer n,cast
      double precision vec(n)
      character funnam*(*), label*(*)
#if MPI|MPIK
C ... Local parameters
      integer, allocatable :: ibuf(:)
      real(8) ,allocatable :: dbuf(:)
      integer obuf

      if (n .le. 0) return
      master = 0
      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )

      if (cast .eq. 2) then
        allocate(ibuf(n), stat=ierr)
        call MPI_ALLREDUCE(vec,ibuf,n,
     .  MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
        call icopy(n,ibuf,1,vec,1)
        deallocate(ibuf, stat=ierr)
      elseif (cast .eq. 4) then
        allocate(dbuf(n), stat=ierr)
        call MPI_ALLREDUCE(vec,dbuf,n,
     .  MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call dcopy(n,dbuf,1,vec,1)
        deallocate(dbuf, stat=ierr)
      elseif (cast .eq. 6) then
        allocate(dbuf(2*n), stat=ierr)
        call MPI_ALLREDUCE(vec,dbuf,2*n,
     .  MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
        call dcopy(2*n,dbuf,1,vec,1)
        deallocate(dbuf, stat=ierr)
      else
        call rxi('mpibc2: cast not implemented',cast)
      endif

      if (mlog) then
        call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
        call strcop(shortname(procid),name,10,'.',ierr)
        namelen(procid) = ierr-1
        call gettime(datim)
c        call awrit2(strn,' ',-256,lgunit(3),procid,numprocs)
        write(stml,"(a,i0,' of ',i0,' on ',a)")
     &       ' '//funnam//' '//datim//' Process ',procid,numprocs,
     &       shortname(procid)(1:namelen(procid))//' allreduce '//label
      endif
#endif
      end

      integer function mpipid(mode)
C- Returns MPI procid
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 return number of processors
Ci         :1 return procid
Ci         :2 calls MPI_BARRIER; returns ierr
Ci         :Otherwise, return 0
Co Outputs
Co   mpipid:procid or number of processors (see mode)
Cr Remarks
Cu Updates
Cu   24 Nov 05 Added mode 2
Cu   14 Apr 03 First created
C ----------------------------------------------------------------------
C     implicit none
#if MPI|MPIK
      include "mpif.h"
      integer numprocs, ierr, procid
#endif
      integer mode

#if MPI|MPIK
      if (mode .eq. 0) then
        call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
        mpipid = numprocs
      else if (mode .eq. 1) then
        call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
        mpipid = procid
      else if (mode .eq. 2) then
        call MPI_BARRIER( MPI_COMM_WORLD, ierr )
        mpipid = ierr
      else
        mpipid = 0
      endif
#else
      mpipid = 0
#endif
      end
C      subroutine fmain
C      integer n
C      n = mpipid(0)
C      print *, 'mpipid for number of processors:', n
C      n = mpipid(1)
C      print *, 'mpipid for processor id:', n
C      call MPI_FINALIZE(n)
C      end

