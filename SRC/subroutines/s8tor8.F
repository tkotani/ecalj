c$$$      subroutine s8tor8(str,dbl)
c$$$      character(8)::dbl,str
c$$$      dbl=str
c$$$C- Store character*(*) string in a double-precision variable.
c$$$C  Appears to conform to ansi standard.  Machines using
c$$$C  ansi arithmetic can store up to 8 characters.
c$$$C     implicit none
c$$$c$$$      double precision dbl
c$$$c$$$      character*(*) str
c$$$c$$$      character*8 strn
c$$$c$$$      strn = str
c$$$c$$$      read(strn,100) dbl
c$$$c$$$  100 format(a)
c$$$      end

c$$$      subroutine r8tos8(dbl,str)
c$$$C- Retrieve character*8 string from a real*8 variable
c$$$C     implicit none
c$$$      character(8)::dbl,str
c$$$      str=dbl
c$$$c$$$      double precision dbl
c$$$c$$$      character*(*) str
c$$$c$$$      character*8 strn
c$$$c$$$      write(strn,100) dbl
c$$$c$$$  100 format(a)
c$$$c$$$      str = strn
c$$$      end

