      subroutine s8tor8(str,dbl)
C- Store character*(*) string in a double-precision variable.
C  Appears to conform to ansi standard.  Machines using
C  ansi arithmetic can store up to 8 characters.
C     implicit none
      double precision dbl
      character*(*) str
      character*8 strn
      strn = str
      read(strn,100) dbl
  100 format(a)
      end
      subroutine r8tos8(dbl,str)
C- Retrieve character*8 string from a real*8 variable
C     implicit none
      double precision dbl
      character*(*) str
      character*8 strn
      write(strn,100) dbl
  100 format(a)
      str = strn
      end

