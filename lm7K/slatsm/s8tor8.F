      subroutine s8tor8(str,dbl)
      character(8)::dbl,str
      dbl=str
C- Store character*(*) string in a double-precision variable.
C  Appears to conform to ansi standard.  Machines using
C  ansi arithmetic can store up to 8 characters.
C     implicit none
c$$$      double precision dbl
c$$$      character*(*) str
c$$$      character*8 strn
c$$$      strn = str
c$$$      read(strn,100) dbl
c$$$  100 format(a)
      end
      subroutine r8tos8(dbl,str)
C- Retrieve character*8 string from a real*8 variable
C     implicit none
      character(8)::dbl,str
      str=dbl
c$$$      double precision dbl
c$$$      character*(*) str
c$$$      character*8 strn
c$$$      write(strn,100) dbl
c$$$  100 format(a)
c$$$      str = strn
      end

