      subroutine strip(str,i1,i2)
C- Returns indices to first and last nonblank characters in a string
C     implicit none
      integer i1,i2
      character*(*) str
      integer i
      i1 = 0
      do  1  i = 1, len(str)
        if(str(i:i) .ne. ' ') then
          i1 = i
          goto 2
        endif
    1 continue
      i1 = 1
      i2 = 0
      return
    2 i2 = len(str) + 1
      do  3  i = len(str), 1, -1
        if(str(i:i) .ne. ' ') then
          i2 = i
          goto 4
        endif
    3 continue
    4 continue
      end

