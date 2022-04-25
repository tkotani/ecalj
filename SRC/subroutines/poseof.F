      subroutine poseof(iunit)
C- Positions file handle at end-of-file
      integer iunit
      integer i,nrec
      nrec = 0
      rewind iunit
      do  10  i = 1, 100000000
         read(iunit,"(a1)",end=90,err=91)
         nrec = i
   10 continue
      write(*,"(' POSEOF: no EOF found for file',i3)") iunit
      return
   90 continue
      backspace iunit
   91 continue
      end

