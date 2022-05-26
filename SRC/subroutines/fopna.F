      subroutine namcat(dest,fnam1,fnam2,lcase)
C- Concatenates two names, optionally lowering the case
      implicit none
      logical lcase
      character*(*) fnam1,fnam2,dest
      integer iprint,i1,i2,j1,j2
      call word(fnam1,1,i1,i2)
      call word(fnam2,1,j1,j2)
      i2 = max(i1-1,i2)
      j2 = max(j1,j2)
      if (i2 .ge. i1) then
        dest = fnam1(i1:i2)//fnam2(j1:j2)
      else
        dest = fnam2(j1:j2)
      endif
      if (i2-i1+j2-j1+2 .gt. len(dest) .and. iprint() .gt. 0) then
        print *, 'namcat: string "',fnam1(i1:i2),fnam2(j1:j2),
     .  '" was truncated to "',dest,'"'
      endif
      if (lcase) call locase(dest)
      end
      
