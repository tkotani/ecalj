C-----------------------------------------------------------------------
C tpzhev.f test pzhev.f
C-----------------------------------------------------------------------
      subroutine fmain
      implicit none
      integer i,j,n,nevmx,nev,i1mach,ndmx,lgunit,fopn,fext,fextg
      double precision eevmx,d1mach
      logical past, test, twice, cmdopt, a2bin, l, lov
      character*72 outs
      character*256 strn
      
C MPI
      include "mpif.h"
C#ifdefC MPE
C      include "mpef.h"
C      integer MPE_LOG_GET_EVENT_NUMBER,MPE_DESCRIBE_STATE,
C     .        MPE_LOG_EVENT,MPE_INIT_LOG,MPE_FINISH_LOG
C#endif 
      integer procid, master, numprocs, ierr, status(MPI_STATUS_SIZE)
      integer MAX_PROCS, dims(2)
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*20 ext
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
      double precision starttime, endtime
      logical mlog

C Heap allocation
      integer wksize
      parameter(wksize= 250 000 000)
      integer w(wksize)
C     Next two lines guarantee w is aligned along a d.p. boundary
      double precision ws
      equivalence (ws,w(1))
      common /w/ w

C Heap pointers
      integer oe,oh,os,ot,ochk1,ochk2

C MPI process configuration
      integer nb,nprow,npcol

C Initialise command line and heap
      call finits(2,0,0,ierr)
      call pshpr(0)
      call wkinit(wksize)
      call poppr

      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
      call strcop(shortname(procid),name,10,'.',i)
      namelen(procid) = i-1
      master = 0
      if (procid .ne. master) then
        call pshpr(0)
        call pshpr(0)
        do  i = 1, 4
          call sprt(i,0)
        enddo
      else
        call pshpr(30)
      endif
C#ifdefC MPE
C      ierr = MPE_INIT_LOG()
C#endif 
      mlog = cmdopt('--mlog',6,0,strn)
      i = fextg(ext)
      call MPI_BCAST(ext,20,MPI_CHARACTER,master,MPI_COMM_WORLD,ierr)
      if (procid .eq. master) then
        call gettime(datim)
        if (mlog) i = fopn('MLOG')
        if (mlog) then
          call awrit2(' tpzhev '//datim//' Process %i of %i on '
     .      //shortname(procid)(1:namelen(procid))//
     .      ' is master',' ',256,
     .      lgunit(3),procid,numprocs)
        endif
      else
        call strcat(ext,20,' ','_',1,' ',i)
        call bin2a(' ',0,0,procid,2,0,20,ext,i)
        ierr = fext(ext(1:i+1))
        if (mlog) ierr = fopn('MLOG')
        ierr = fextg(ext)
        call gettime(datim)
        if (mlog) then
          call awrit2(' tpzhev '//datim//' Process %i of %i on '
     .      //shortname(procid)(1:namelen(procid))//
     .      ' file extension is '//ext(2:i+1),' ',
     .      256,lgunit(3),procid,numprocs)
        endif
      endif
      call MPI_BARRIER( MPI_COMM_WORLD, ierr )

C MPI process configuration
      nb = 16
      nprow = numprocs
      npcol = -1

C Get dimension
      j = 3
      if (cmdopt('-n=',j,0,outs)) then
        if (a2bin(outs,n,2,0,' ',j,72)) goto 2
      endif
      if (procid .eq. master) print *,
     .  'Usage: tpzhev --g --test --twice -n=# -nb=# -nprow=# -npcol=#'
      goto 3
    2 continue 
      test = .false.
      if (cmdopt('--test',6,0,outs))  test = .true.
      twice = .false.
      if (cmdopt('--twice',7,0,outs)) twice = .true.
      lov = .false.
      if (cmdopt('--g',3,0,outs))     lov = .true.
      j = 4
      if (cmdopt('-nb=',j,0,outs))    l = a2bin(outs,nb,2,0,' ',j,72)
      j = 7
C MPI process configuration
      if (cmdopt('-nprow=',j,0,outs)) l = a2bin(outs,nprow,2,0,' ',j,72)
      j = 7
      if (cmdopt('-npcol=',j,0,outs)) l = a2bin(outs,npcol,2,0,' ',j,72)
      if (npcol .eq. -1) then
        dims(1) = 0
        dims(2) = 0
        call MPI_DIMS_CREATE(numprocs,2,dims,ierr)
        npcol = dims(1)
        nprow = dims(2)
        if (procid .eq. master) then
          call awrit2(
     .      ' MPI creating process configuration .. nprow=%i npcol=%i',
     .      ' ',256,lgunit(1),nprow,npcol)
        endif
      endif
      if (procid .eq. master) then
        call awrit4(' Testing pzhev, n=%i nb=%i nprow=%i npcol=%i',' ',
     .              128,lgunit(1),n,nb,nprow,npcol)
      endif
      nevmx = n
      eevmx = 1d12
    4 continue 
      call defrr(oe, n)
      call defcc(oh, n*n)
      if (lov) then
        call defcc(os, n*n)
      else
        os = 5
      endif
      call mkhs(lov,n,w(oh),w(os))
      if (procid .eq. master) starttime = MPI_WTIME()

      call pzhev(lov,n,oh,os,nb,nprow,npcol,eevmx,nevmx,nev,w(oe),ot)

      if (procid .eq. master) then
        endtime = MPI_WTIME()
        call awrit2('%N MPI Wall-clock time: %;3ds. Resolution %;9ds%N',
     .              ' ',128,lgunit(1),endtime-starttime,MPI_WTICK())
      endif

C Test the result
      if (test) then
        call mkhs(lov,n,w(oh),w(os))
        call defcc(ochk1, n*n)
        call zmpy(w(oh),2*n,2,1,w(ot),2*n,2,1,w(ochk1),2*n,2,1,n,n,n)
        if (lov) then
          call defcc(ochk2, n*n)
          call zmpy(w(os),2*n,2,1,w(ot),2*n,2,1,w(ochk2),2*n,2,1,n,n,n)
        else
          ochk2 = ot
        endif
        call check(n,w(ochk1),w(ochk2),w(oe),past)
        if (past .and. procid .eq. master) print *, 'Test passed.'
        if (.not. past .and. procid .eq. master) print *, 'Test failed.'
      endif
      call rlse(oe)
      if (twice) then
        twice = .false.
        goto 4
      endif
    3 continue 
C#ifdefC MPE
C      ierr = MPE_FINISH_LOG('tpzhev')
C#endif 
      if ( procid .eq. master ) then
        call fexit(i,011,
     .   'EXIT tpzhev on '//shortname(procid)(1:namelen(procid)),0)
      else
        call fexit(i,0,' ',0)
      endif
      end

      subroutine mkhs(lov,n,h,s)
      implicit none
      logical lov
      integer n,i,j
      double complex h(n,n),s(n,n)
      do i = 1, n
        do j = 1, n
          if (lov) s(i,j) = dcmplx(n-abs(i-j))
          if (i .eq. j) then
            h(i,j) = dcmplx(1d0/(dble(i+j)-1d0))+dcmplx(1d0)
          else
            h(i,j) = dcmplx(1d0/(dble(i+j)-1d0),dble(j-i))
          endif
        enddo
      enddo
      end


      subroutine check(n,c,zz,e,past)
      implicit none
      integer n,i,j,i1mach,iprint
      complex*16 c(n,n),zz(n,n)
      double precision e(1),tol
      logical past

      tol = 1d-10
      past = .true.
      do  10  i = 1, n
      do  10  j = 1, n
        if (abs(c(i,j) - e(j)*zz(i,j)) .gt. tol) then
          if (iprint() .gt. 60) then
            call awrit3('H Z  -  E O Z != 0: i=%i, j=%i,'//
     .        ' c(i,j) - e(j)*zz(i,j)=%d',' ',256,i1mach(2),i,j,
     .        c(i,j) - e(j)*zz(i,j))
          endif
          past = .false.
        endif
   10 continue 
      end
