      subroutine headl3(name,vrsn,wksiz,ifi)
C-  Puts a heading line into file, unit ifi
C     implicit none
      integer ifi,wksiz
      double precision vrsn
      character*8 name,outs*80

      outs = ' '
      if (vrsn .eq. 0) call headl2(name,wksiz,ifi)

      if (wksiz .eq. 0)
     .  call awrit1(' -----------------------  START '
     .  //name//'%a (%d)  -----------------------',
     .  outs,-80,-ifi,vrsn)
      if (wksiz .ne. 0)
     .  call awrit2(' -----------------------  START '
     .  //name//'%a (v%d, %iK)  -----------------------',
     .  outs,-80,-ifi,vrsn,(wksiz+499)/1000)
      end

      subroutine headl2(name,wksiz,ifi)
C-  Puts a heading line into file, unit ifi
C     implicit none
      integer ifi,wksiz
      character*8 name,outs*80

      outs = ' '
      if (wksiz .eq. 0) write(ifi,300) name
  300 format(' -----------------------  START ',a8,
     .  ' -----------------------')
      if (wksiz .ne. 0) call awrit1(' -----------------------  START '
     .  //name//'%a (%iK)  -----------------------',
     .  outs,-80,-ifi,(wksiz+499)/1000)
      end
      subroutine poseof(iunit)
C- Positions file handle at end-of-file
C     implicit none
C Passed parameters 
      integer iunit
C Local parameters 
      integer i,nrec
#if SVS | DEC | DECA
      double precision dummy
#endif

      nrec = 0
      rewind iunit
      do  10  i = 1, 100000000
C ... Avoid read bug ...
#if SVS | DEC | DECA     
        read(iunit,100,end=90,err=91) dummy
#endif                
C ... Avoid AIX bug
#if AIX
        read(iunit,100,end=90,err=90)
#else
        read(iunit,100,end=90,err=91)
#endif
        nrec = i
  100   format(a1)
   10 continue
      write(*,200) iunit
  200 format(' POSEOF: no EOF found for file',i3)
      return
   90 continue
C ... Avoid AIX bug in xlf90 compiler
#ifndef AIX_xlf90
      backspace iunit
#endif
C ... If backspace doesn't work, do it this way ...
C     rewind iunit
C     do 11 i=1,nrec
C  11 read(iunit,100)
   91 continue
      end

