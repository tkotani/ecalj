C  error exit routines.
C  some routines print out a number before exiting.

      subroutine rxx(test,string)
C- Test for error exit
C     implicit none
      logical test
      character*(*) string

      if (test) call rx(string)
      end

      subroutine rx(string)
c  routine for error exit
      character*(*) string
      write(6,892) string
      write(6,890)
  890 format(' ---- ERROR EXIT! ')
  892 format(/' ---- ',a)
c ... write to log file also
      write(0,710) string
  710 format('ERROR EXIT! ',a)
      call exit(-1)
      end

      subroutine rxs(string,msg)
c  error exit with extra string message
      character*(*) string
      character*(*) msg
      call strip(msg,i1,i2)
      write(6,892) string,msg(i1:i2)
      write(6,890)
  890 format(' ERROR EXIT! ')
  892 format(/' ---- ',a,a)
c ... write to log file also
      write(0,710) string,msg(i1:i2)
  710 format('ERROR EXIT! ',a,a)
      call exit(-1)
      end

      subroutine rxi(string,i)
c  error exit with printout of one integer
      character*(*) string
      write(6,892) string,i
      write(6,890)
  892 format(/' ---- ',a,i8)
  890 format(' ERROR EXIT! ')
      write(0,710) string,i
  710 format('ERROR EXIT! ',a,i8)
      call exit(-1)
      end

      subroutine rxf(string,a)
c  error exit with printout of one float
      real*8 a
      character*(*) string
      write(6,892) string,a
      write(6,890)
  892 format(/' ---- ',a,f15.8)
  890 format(' ---- Error exit')
      write(0,710) string,a
  710 format('ERROR EXIT! ',a,f15.8)
      call exit(-1)
      end

      subroutine rx0(string,msg)
c  write extra information before error exit in standard format
      character*(*) string
      character*(*) msg
      call strip(msg,i1,i2)
      write(6,892) string,msg(i1:i2)
  892 format(/' ---- ',a,a)
      write(0,710) string,msg(i1:i2)
  710 format('ERROR EXIT! ',a,a)
      end

      subroutine rx00(string,msg)
c  write extra information to log but not to listing
      character*(*) string
      character*(*) msg
      call strip(msg,i1,i2)
      write(0,710) string,msg(i1:i2)
  710 format('ERROR EXIT! ',a,a)
      end

c --- subs to write warning messages
      subroutine wg(string)
C- write warning to listing and log
      character*(*) string
      write(6,892) string
  892 format(' ---- warning: ',a)
      write(0,710) string
  710 format('== ',a)
      end

      subroutine wgs(string,msg)
C- write warning to listing and log
      character*(*) string,msg
      write(6,892) string,msg
  892 format(' ---- warning: ',a,a)
      write(0,710) string,msg
  710 format('== ',a,a)
      end

      subroutine wgi(string,i)
C- write warning to listing and log with one integer
      character*(*) string
      write(6,892) string,i
  892 format(' ---- warning: ',a,i8)
      write(0,710) string,i
  710 format('== ',a,i8)
      end


      subroutine wgf(string,a)
C- write warning to listing and log with one float
      character*(*) string
      real*8 a
      write(6,892) string,a
  892 format(' ---- warning: ',a,f15.8)
      write(0,710) string,a
  710 format('== ',a,f15.8)
      end

