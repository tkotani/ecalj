      logical function cmdopt(argstr,strln,nargs,outstr)
C- Determines whether a command-line argument supplied, and its argument
C ----------------------------------------------------------------
Ci Inputs
Ci   argstr,strln: command-line string to search; search to strln chars
Ci   nargs:        number of arguments associated with argstr
Co Outputs
Co   cmdopt: true if argument found, else false
Co   outstr (nargs>0 only) nth string after string argstr
Cr Remarks
Cu    3 Aug 04 Changed call to nargc with call to nargf
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      character*(*) argstr,outstr
      integer nargs,strln
C Local parameters
      logical lsequ
      integer iarg,nargf,idum,nxarg
      character*120 strn

      cmdopt = .false.
      iarg = 0
   10 iarg = iarg+1
C ... A usual command-line argument
      if (nargf() .gt. iarg) then
        call getarf(iarg,strn)
C ... If not on the command-line, try 'extra' arguments
      else
        call ncmdop(nxarg)
*       print *, nxarg,iarg-nargf()
        if (nxarg .le. iarg-nargf()) return
        call gcmdop(iarg-nargf()+1,strn)
      endif
      if (.not. lsequ(strn,argstr,strln,' ',idum)) goto 10
      cmdopt = .true.
      outstr = ' '
      if (nargf() .gt. iarg+nargs) then
        call getarf(iarg+nargs,outstr)
      elseif (nargf() .gt. iarg+nargs+nxarg) then
        call rx('bug in CMDOPT')
      else
        call gcmdop(iarg-nargf()+1,outstr)
      endif

      end

      subroutine acmdop(strn,lstr,opt)
C- Append strn to 'extra' command options
C ----------------------------------------------------------------
Ci Inputs
Ci   strn:  (acmdop,opt=0) is appended to internal cmdarg
Ci   lstr:  length of input string
Ci   opt:   0 append strn
Ci          1 print out table
Ci   iarg   (gcmdop) retrieve argument n
Co Outputs
Co   n       (ncmdop) number of arguments in list
Co   n       (gcmdop) number of arguments in list
Cr Remarks
C ----------------------------------------------------------------
C     implicit none
      integer lstr,opt
      character*1 strn(1), sout*(*)
C Local variables
      integer mxarg,lcmd
#if SUN_ULTRA
      parameter (mxarg=2000,lcmd=2048)
#else
      parameter (mxarg=2000,lcmd=20000)
#endif
      integer marker(0:mxarg),nxarg,i1,i2,it,ia,i,lgunit,n,iarg
      character*(lcmd) cmdarg, ch*3
      save nxarg,marker
      data cmdarg /' '/ ch /' "'''/, nxarg /0/

*     print *, (strn(i), i=1,lstr)
      marker(0) = 1
      if (opt .eq. 1) goto 100

      i2 = -1
   10 continue
      i1 = i2+1
C --- Continue until all command arguments exhausted ---
      call skipbl(strn,lstr,i1)
      if (i1 .lt. lstr) then
        ia = marker(nxarg)
C   ... Find i2 : points to past last char of the argument
        i2 = i1
   12   i1 = i2
        call chrps2(strn,ch,3,lstr,i2,it)
        call strcop(cmdarg(ia:),strn(i1+1),i2-i1,' ',i)
        ia = ia+i
        if (ia .ge. lcmd) call rx('acmdop: increase lcmd')
*       print *, cmdarg(1:ia)
C   ... A quote encountered ... continue copying string
        if (it .gt. 1) then
          i2 = i2+1
          i1 = i2
          call chrpos(strn,ch(it:it),lstr,i2)
          call strncp(cmdarg,strn,ia,i1+1,i2-i1)
          ia = ia+i2-i1
          if (ia .ge. lcmd) call rx('acmdop: increase lcmd')
*          print *, cmdarg(1:ia)
          i2 = i2+1
          if (i2 .le. lstr) goto 12
        endif
C   ... End of this argument ... start on another
        nxarg = nxarg+1
        if (nxarg .gt. mxarg) call rx('acmdop: increase mxarg')
        marker(nxarg) = ia
        goto 10
      endif
      if (opt .eq. 0) return

C  --- Printout of command line arguments ---      
  100 continue
      call awrit1(' acmdop:  %i file command switches:',
     .  ' ',80,lgunit(1),nxarg)
      do  20  i = 1, nxarg
   20 print 333, i, cmdarg(marker(i-1):marker(i)-1)
  333 format(i4,2x,'"',a,'"')
      return

      entry ncmdop(n)
      n = nxarg
      return

      entry gcmdop(iarg,sout)
      if (iarg .gt. nxarg) then
        sout = ' '
        return
      endif
      sout = cmdarg(marker(iarg-1):marker(iarg)-1)
      return

      end

