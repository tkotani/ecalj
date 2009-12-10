      subroutine cpudel(unit,strn,delt)
C- incremental cup time, in seconds
C ----------------------------------------------------------------------
Ci Inputs:
Ci   unit>=0: printout of time to unit; else no printing
Ci   unit<-1: delt not calculated
Co Outputs:
Co   delt: incremental cpu time since last call, seconds
Cr Remarks
Cr   Uses cpusec
C ----------------------------------------------------------------------
C Passed parameters
      character*(*) strn
      integer unit
      double precision delt
C Local parameters
      character*1 timeu, outs*80
      double precision cpusec,told,tnew
      save told
      data told /0d0/

      if (unit .lt. -1) return
      tnew = cpusec()
      delt = tnew - told
      told = tnew
      timeu = 's'
      if (tnew .gt. 60) then
        timeu = 'm'
        tnew = tnew/60
        if (tnew .gt. 60) then
        timeu = 'h'
        tnew = tnew/60
        endif
      endif
      if (unit .ge. 0 .and. tnew .gt. 0d0) then
        outs = ' '
        write(outs,333) strn
        call awrit2('%a  %1;3,3g%53ptotal:  %1,3;3g'//timeu,
     .    outs,len(outs),-unit,delt,tnew)
      endif
C      if (unit .ge. 0 .and. tnew .gt. 0d0)
C     .  write(unit,333) strn,delt,tnew,timeu
  333 format(' cpudel',a25,'  time(s):',g10.3,'  total:',f8.3,a1)
      end
      double precision function cpusec()
C- process cputime, in seconds
C ----------------------------------------------------------------------
Ci Inputs:
Ci   none
Co Outputs:
Co   returns cpu time, in seconds
Cr  Remarks
Cr    On the Apollo: time (in 4 microsecond units) is
Cr    (time(1)*65536 + time(2))*65536 + time(3)
C ----------------------------------------------------------------------
#if APOLLO
      integer time(3)
      integer*2 t2(3)
      call proc1_$get_cput(t2)
      time(1) = t2(1)
      time(2) = t2(2)
      time(3) = t2(3)
      if (time(1) .lt. 0) time(1) = time(1) + 65536
      if (time(2) .lt. 0) time(2) = time(2) + 65536
      if (time(3) .lt. 0) time(3) = time(3) + 65536
      cpusec = 4d-6 * 
     .  ((dble(time(1))*65536 + dble(time(2)))*65536 + dble(time(3)))
#elif AIX
      integer mclock
      cpusec = dble(mclock())/100
#elif CRAY | HASSECOND
      real(4):: second 
c           ---- not real::second
      cpusec = second()
#elif VANILLA
      cpusec = 0
#else
      real(4):: tarray(2),timet
c            --------- not real:: tarray(2),timet
      real etime
      real etime_
      tarray(2)=0.0
#if POWER_PC
      timet = etime_(tarray)
#elif HASETIME
      timet = etime(tarray)
#endif
      timet = tarray(1)
      cpusec = dble(timet)
#endif
      end

