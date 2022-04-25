      subroutine finits()       !job,fcn,fcargs,iarg)
      use m_ext,only:sname
      implicit none
!! Set the numbr of iprint(), extension, and -vnam=val
      integer iarg
      double precision fcargs(1)
      logical lsequ,lext
      integer i,fext,nargf,n,it(5),iv(5),a2vec,k
      logical:: cmdopt
      character strn*256
      character*100 extns
!! Command line arguments and extension ---
      lext = .false.
      do iarg = 1,nargf()-1
        call getarg(iarg,strn)
        extns = strn
!  ... v encountered ... parse variables
        if (lsequ(strn,'-v',2,' ',n)) then
           i = 2
           call parsyv(strn,len(strn),999,0,i)
        endif
        if (lsequ(strn,'-',1,' ',n)) cycle 
      enddo
      end

      subroutine fexit0(retval,iopt,strng)
      implicit none
C- Machine and compiler-dependent program termination
C ----------------------------------------------------------------------
Ci Inputs
Ci   retval:  return value passed to operating system
Ci   iopt decomposed into 3 one-digit numbers.
Ci   digit
Ci     1:  0: do not print string on exit; 
Ci         9: print strng as Exit(retval): 'strng'
Ci      else: exit, using strn as a format statement and args a vector
Ci            of  c  double precision arguments
Ci    10:   0: do not print cpu time, else do
Ci   100:   0: do not print work array usage, else do
Co Outputs
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
C Passed parameters 
      integer retval,iopt
      character*(*) strng
      double precision args,arg2,arg3,argss(3)
      integer fhndl,iprint,i,i2,getdig,scrwid,mpipid,ix
      double precision cpusec,tnew
      character*1 timeu
      character*256 strn, datim*26, hostnm*20
c      character(8):: charext
      character(9):: ftoa9
      logical isopen
      integer master,procid,ierr,stdo,stdl,lgunit,ia
      parameter (master = 0)
      include "mpif.h"
      ia=0
      goto 5
      
      entry fexit(retval,iopt,strng,args)
      ia=1
      goto 5
      
      entry fexit2(retval,iopt,strng,args,arg2)
      ia=2
      goto 5
      
      entry fexit3(retval,iopt,strng,args,arg2,arg3)
      ia=3
    5 continue
      if(ia>=1) argss(1)=args
      if(ia>=2) argss(2)=arg2
      if(ia>=3) argss(3)=arg3
      stdo=lgunit(1)
      stdl=lgunit(2)
      procid = mpipid(1)
      if (procid .eq. master.or. retval/=0) then
        i = mod(iopt,10) !getdig(iopt,0,10) 
        if (i .ne. 0) then
           if (i .eq. 9) ia=0
           write(stdo,"(a,i0,a,3d15.8)")"Exit ",retval,' '//trim(strng),(argss(ix),ix=1,ia)
        endif
        i = mod(iopt/10,10) !getdig(iopt,1,10)
        if (i .ne. 0 .and. cpusec() .ne. 0) then
          timeu = 's'
          tnew = cpusec()
          if (tnew .gt. 3600) then
            timeu = 'm'
            tnew = tnew/60
            if (tnew .gt. 200) then
              timeu = 'h'
              tnew = tnew/60
            endif
          endif
          datim = ' '
          call ftime(datim)
          hostnm = ' '
          call get_environment_variable('HOST',hostnm)
          write(stdo,10) tnew,timeu,datim,trim(adjustl(hostnm))
          write(stdl,10) tnew,timeu,datim,trim(adjustl(hostnm))
 10       format('CPU time:', f9.3,a1,5x,a26,' on ',a)
        endif
      endif
      call tcprt(stdo) 
      if(retval/=0) then
         write(6,"(a,i5)")'ERROR: fexit, error code=',retval
         write(6,*)
         call MPI_abort(MPI_comm_world,retval,ierr)   
         call MPI_FINALIZE(i)
         stop
      endif
      call cexit(retval,1)
      stop
      end

      subroutine rx0s(string) !normal exit for master_mpi
      use m_MPItk,only: master_mpi
      character*(*) string
      if(master_mpi) write(6,892) trim(string)//' ======================'
  892 format(/,a)
      call exit(0)
      end

      subroutine rx(string) !error exit
      character*(*) string
      call fexit0(-1,119,string)
      end

      subroutine rx0(string) !normal exit
      character*(*) string
      call fexit0(0,119,string)
      call rx0s(string)         !for single core exit
      end

      subroutine rx1(string,arg) ! Error exit, with a single argument
      character*(*) string
      double precision arg
      character*120 outs
      outs = ' Exit -1 '//string
      call fexit(-1,111,outs,arg)
      end
      
      subroutine rx2(string,arg1,arg2) ! Error exit, with two arguments
      character*(*) string
      double precision arg1,arg2
      character*120 outs
      outs = ' Exit -1 '//string
      call fexit2(-1,111,outs,arg1,arg2)
      end
      
      subroutine rxi(string,arg) ! Error exit, with a single integer at end
      character*(*) string
      integer:: arg
      character(10):: i2char
      character*120 outs
      outs = ' Exit -1 '//string//' '//trim(i2char(arg))
      call fexit0(-1,119,outs)
      end

      subroutine rxs(string,msg) ! Error exit with extra string message
      character*(*) string,msg
      character*120 outs
      integer i
      outs = string // msg
      call skpblb(outs,len(outs),i)
      call rx(outs(1:i+1))
      end
      
      subroutine rxs2(string,msg,msg2) ! Error exit with extra string messages
      character*(*) string,msg,msg2
      character*120 outs
      integer i
      outs = string // msg // msg2
      call skpblb(outs,len(outs),i)
      call rx(outs(1:i+1))
      end
      
      subroutine rxs4(string,msg,msg2,msg3,msg4) ! Error exit with extra string messages
      character*(*) string,msg,msg2,msg3,msg4
      character*120 outs
      integer i
      outs = string // msg // msg2 // msg3 // msg4
      call skpblb(outs,len(outs),i)
      call rx(outs(1:i+1))
      end
      
      subroutine rxx(test,string)
      logical test
      character*(*) string
      if (test) call rx(string)
      end
      
      subroutine rx_(string) !error exit routine
      character*(*) string
      write(6,892) string
      write(6,890)
  890 format(' ---- Error exit')
  892 format(/' ---- ',a)
      write(71,710) string
  710 format('++ ',a)
      call exit(-1)
      end
