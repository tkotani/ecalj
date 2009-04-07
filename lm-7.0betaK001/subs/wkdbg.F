      subroutine wkdbg()
      call wkinfo
      call wkchk('wkdbg')
      call wkprnt(1)
C      call tcget(ltc)
C      if (ltc .eq. 0) then
C        call tc('on')
C      else
C        call tc('from wkdbg')
C      endif
      end

      subroutine wkdbg2()
      integer mfree,mnow,mmax
      call wquery(mfree,mnow,mmax)
      print '(a,i9,a,i9)',' mem used now = ',mnow,
     .  ' max mem used = ',mmax
      call wkprnt(0)
      end

