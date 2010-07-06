      subroutine extstr(strn,lstr,i1,i2)
C- Extracts string that may be in quotes '...' or "..."
C Outputs
C   i1,i2:  strn(i1..i2) is resultant
C     implicit none
      integer lstr,i1,i2
      character*1 strn(lstr)
      character*81 ch(3)
      integer k
      data ch /'''','"',' '/

      k = 3
      if (strn(1) .eq. ch(1)) k=1
      if (strn(1) .eq. ch(2)) k=2
      i1 = 1
      if (k .ne. 3) i1 = 2
      i2 = 0
      call chrpos(strn(i1),ch(k),lstr-i1+1,i2)
      i2 = i2+i1
      if (strn(i2) .eq. ch(k)) i2 = i2-1

      end


