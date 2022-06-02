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
c$$$      if (unit .ge. 0 .and. tnew .gt. 0d0) then
c$$$        outs = ' '
c$$$        write(outs,333) strn
c$$$        call awrit2('%a  %1;3,3g%53ptotal:  %1,3;3g'//timeu,
c$$$     .       outs,len(outs),-unit,delt,tnew)
c$$$      endif
      if (unit .ge. 0 .and. tnew .gt. 0d0)
     .  write(unit,333) strn,delt,tnew,timeu
  333 format(' cpudel',a25,'  time(s):',g10.3,'  total:',f8.3,a1)
      end
      double precision function cpusec()
      implicit none
      logical,save::firsttime=.true.
      integer,save::i1
      integer:: i2,irate,imax
      real(4)::diff

      if (firsttime) then 
         call system_clock(i1)
         firsttime=.false.
         cpusec=0.0d0
        ! cpumin=0.0d0
      else
         call system_clock(i2,irate,imax)
         diff=i2-i1
         if (diff<0) diff=imax-i1+i2
         diff=diff/dble(irate)
         cpusec=diff
        ! cpumin=cpusec/60.d0
      endif
      end function cpusec
