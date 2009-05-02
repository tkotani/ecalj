c Time routines -------------------------------------------------------------
c This dummy time routine
#ifdef ETIME_DUMMY
      real*4 function etime(ar)
      real*4 ar(2)
      ar = 0.0
      end
#endif
c
c       subroutine fdate(s)
c       character(*) s
c       s = ' '
c       end
c----------------------------------------------------------------------------





c----------------------------------------------------------------------
      double precision function cpumin(ifile,cpu0)
      implicit real*8(a-h,o-z)

c cpu time in minutes for RISC: mclock()/600.d0
c cpusek is the total cputime - cpu0

c     cpumin   = mclock()/600.d0 - cpu0
      cpusec   = (secnds(0.0) - cpu0)
      cpumin   = cpusec/60.d0
c     cpusec=0d0
c     cpumin=0d0
c----
      if(ifile .le. 0)then
        write(*,6400)cpusec,cpumin
      else
        write(ifile,6400)cpusec,cpumin
      end if
 6400 format(/' cpu time = ',f15.4,' secs ',f15.4,' mins ')

      return
      end
c----------------------------------------------------------------------
c     double precision function cpusec(cpu0)
c     implicit real*8(a-h,o-z)
c
c cpu time in seconds for RISC
c cpusek is the total cputime - cpu0
c
c     cpusec   = mclock()/100.d0 - cpu0
c     return
c     end
c----------------------------------------------------------------------
      subroutine cputid(ifile)
      implicit real*8(a-h,o-z)

c cpu time in seconds and minutes for RISC IBM: mclock()/100 (secs)
c                     for HP: secnds(0.0)
c                        DEC: etime(etw)
c ifile = file number where the cpu time information is
c         to be printed
c ifile = 0 ==> screen

      real*4 cpuetime, etw(2),cpulast, etime
      data cpu0 /-1.d0/
      save cpu0

      if(ifile .lt. 0)stop 'cputid: negative unit number'
c     if (cpu0 .lt. 0.d0) cpu0 = mclock()/100.d0
c     cpusec   = mclock()/100.d0

c     if (cpu0 .lt. 0.d0) cpu0 = secnds(0.0)
c     cpusec   = secnds(0.0) - cpu0
c     cpusec   = secnds(0.0)

      cpuetime = etime(etw)
      if (cpu0 .lt. 0.d0) cpu0 = cpuetime
      cpusec   = cpuetime - cpu0
      cpumin   = cpusec/60.d0

      if(ifile .eq. 0)then
        write(*,6400)cpusec,cpumin
      else
        write(ifile,6400)cpusec,cpumin
      end if

c 6400 format(/' cpu time = ',f15.4,' secs ',f15.4,' mins ')
 6400 format(' CPU= ',f12.4,' secs',f7.1,' mins')
      return
      end
