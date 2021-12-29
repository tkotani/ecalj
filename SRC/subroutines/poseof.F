      subroutine headl3(name,vrsn,wksiz,ifi)
C-  Puts a heading line into file, unit ifi
C     implicit none
      integer ifi,wksiz
      double precision vrsn
      character*8 name,outs*80

      outs = ' '
      if (vrsn .eq. 0) call headl2(name,wksiz,ifi)

      if (wksiz .eq. 0)
     .call awrit1(' -----------------------  START '
     .//name//'%a (%d)  -----------------------',
     .outs,-80,-ifi,vrsn)
      if (wksiz .ne. 0)
     .call awrit2(' -----------------------  START '
     .//name//'%a (v%d, %iK)  -----------------------',
     .outs,-80,-ifi,vrsn,(wksiz+499)/1000)
      end

      subroutine headl2(name,wksiz,ifi)
C-  Puts a heading line into file, unit ifi
C     implicit none
      integer ifi,wksiz
      character*8 name,outs*80
      outs = ' '
      if (wksiz .eq. 0) write(ifi,300) name
  300 format(' -----------------------  START ',a8,
     .' -----------------------')
      if (wksiz .ne. 0) call awrit1(' -----------------------  START '
     .//name//'%a (%iK)  -----------------------',
     .outs,-80,-ifi,(wksiz+499)/1000)
      end
      
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

