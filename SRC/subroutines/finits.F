      subroutine finits(job,fcn,fcargs,iarg)
C- Machine and compiler-dependent inits for standard FORTRAN startup
C ----------------------------------------------------------------
Ci Inputs
Ci   job:  0, no command-line arguments; 
Ci         1, no switches (extens only);
Ci         2, [-vnam=val ...] [-pr#] [switches] extens;
Ci         3, call fcn for switches first
Ci   fcn,fcargs:  see Remarks
Co Outputs
Co   iarg: argument to extens, if found (0 if not)
Cr Remarks
Cr   finits parses command line arguments, ignoring args that
Cr   begin with "-", to find file extension.  When job=1,
Cr   switches -vnum are taken to be variable defs.
Cr   Calling prog can set own switches through external
Cr   fcn(iarg,fcargs).  fcn should return iarg as last arg parsed.
Cu Updates
Cu   16 Jun 06 Added MPI Wall clock time printout
Cu    3 Aug 04 Changed call to nargc with call to nargf
C ----------------------------------------------------------------
C     implicit none
C Passed parameters      
      integer job,iarg
      double precision fcargs(1)
      logical lsequ,lext
      integer i,fext,nargf,n,it(5),iv(5),a2vec
      character strn*256
      character*100 extns
#if MPI | MPIK
      integer master,procid
      parameter (master = 0)
      include "mpif.h"
      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, i )
#endif
      call awrit0('%u',' ',80,0)
      call initqu(.false.)

C --- Command line arguments and extension ---
      if (job .eq. 0) return
      iarg = 0
      if (job .eq. 3) call fcn(fcargs,iarg)
      lext = .false.
   10 iarg = iarg+1
      if (nargf() .gt. iarg) then
        call getarf(iarg,strn)
        extns = strn
        if (job .ge. 2) then
C     ... -pr encountered ... set verbosity
          if (lsequ(strn,'-pr',2,' ',n)) then
            i = 3
            n = a2vec(strn,len(strn),i,2,', ',2,2,5,it,iv)
            do  12  i = 1, n
   12       call sprt(i-1,iv(i))
          endif
C     ... v encountered ... parse variables
          if (lsequ(strn,'-v',2,' ',n)) then
            i = 2
            call parsyv(strn,len(strn),999,0,i)
c            call shosyv(0,0,0,i1mach(2))
c            stop
          endif
        endif
        if (lsequ(strn,'-',1,' ',n)) goto 10
        if (.not. lext) then
          i = fext('.'//extns)
          lext = .true.
          goto 10
        endif
      endif
      iarg = iarg-1
      end

      subroutine fexit(retval,iopt,strng,args)
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
Cu   03 Jul 03 Open file tmp with fopnT, for local directories (MPI)
C ----------------------------------------------------------------------
      implicit none
C Passed parameters 
      integer retval,iopt
      character*(*) strng
      double precision args(1),arg2(1),arg3(1)
C Local parameters 
      integer fopn,fhndl,iprint,i,i2,getdig,i1mach,scrwid,fopnT,mpipid
c      parameter (scrwid=300)
      double precision cpusec,tnew
      character*1 timeu
      character*256 strn, datim*26, hostnm*20
      logical isopen
      integer master,procid,ierr
      parameter (master = 0)
#if MPI | MPIK
      include "mpif.h"
#endif
      goto 5
      entry fexit3(retval,iopt,strng,args,arg2,arg3)
      entry fexit2(retval,iopt,strng,args,arg2)
    5 continue
      procid = mpipid(1)
      if (procid .eq. master.or. retval/=0) then
        i = getdig(iopt,0,10) 
        if (i .ne. 0) then
           if (i .eq. 9) then
              strn = ' Exit %i '//strng
c              call strip(strn,i,i2)
              call awrit1(strn,' ',len_trim(strn),i1mach(2),retval)
           else
             call awrit3(strng,strn,-len_trim(strn),i1mach(2),args,arg2,arg3)
           endif
        endif
        i = getdig(iopt,1,10)
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
c          call word(hostnm,1,i,i2)
c          i2 = max(i,i2)
          write(i1mach(2),10) tnew,timeu,datim,trim(adjustl(hostnm))
          if (fhndl('LOG') .ge. 0)  then
            if (isopen(fhndl('LOG'),.false.))
     .      write (fhndl('LOG'),10) tnew,timeu,datim,trim(adjustl(hostnm))
   10       format(' CPU time:', f9.3,a1,5x,a26,' on ',a)
          endif
        endif
      endif
      call tcprt(i1mach(2))
#if MPI | MPIK
      if(retval/=0) then
         print *,'fexit,fexit2,fexit3 error retval=',retval
         print *
         call MPI_abort(MPI_comm_world,retval,ierr)   
         call MPI_FINALIZE(i)
         stop
      endif
#endif 
      call cexit(retval,1)
      stop
      end

      subroutine rx0s(string)
c  write extra information before error exit in standard format
      character*(*) string
      write(6,892) string
  892 format(/' ---- ',a,a)
      call exit(0)
      end

      subroutine rx(string)
C- Error exit
C     implicit none
      character*(*) string
      call fexit(-1,119,string,0d0)
      end

      subroutine rx0(string)
C- Normal exit
C     implicit none
      character*(*) string
      call fexit(0,119,string,0d0)
      end

      subroutine rx1(string,arg)
C- Error exit, with a single argument
C     implicit none
      character*(*) string
      double precision arg
      character*120 outs
      outs = '%N Exit -1 '//string
      call fexit(-1,111,outs,arg)
      end
      subroutine rx2(string,arg1,arg2)
C- Error exit, with two arguments
C     implicit none
      character*(*) string
      double precision arg1,arg2
      character*120 outs
      outs = '%N Exit -1 '//string
      call fexit2(-1,111,outs,arg1,arg2)
      end
      subroutine rxi(string,arg)
C- Error exit, with a single integer at end
C     implicit none
      character*(*) string
      double precision arg
      character*120 outs
      outs = '%N Exit -1 '//string//' %i'
      call fexit(-1,111,outs,arg)
      end

      subroutine rxs(string,msg)
C- Error exit with extra string message
C     implicit none
      character*(*) string,msg
      character*120 outs
      integer i
      outs = string // msg
      call skpblb(outs,len(outs),i)
      call rx(outs(1:i+1))
      end
      subroutine rxs2(string,msg,msg2)
C- Error exit with extra string messages
C     implicit none
      character*(*) string,msg,msg2
      character*120 outs
      integer i
      outs = string // msg // msg2
      call skpblb(outs,len(outs),i)
      call rx(outs(1:i+1))
      end
      subroutine rxs4(string,msg,msg2,msg3,msg4)
C- Error exit with extra string messages
C     implicit none
      character*(*) string,msg,msg2,msg3,msg4
      character*120 outs
      integer i
      outs = string // msg // msg2 // msg3 // msg4
      call skpblb(outs,len(outs),i)
      call rx(outs(1:i+1))
      end
      subroutine rxx(test,string)
C- Test for error exit
C     implicit none
      logical test
      character*(*) string

      if (test) call rx(string)
      end
