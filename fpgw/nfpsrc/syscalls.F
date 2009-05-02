c ----- sysinf -- return system information as a string ------
      subroutine sysinf(job,nt,mt,text)
c  Returns various items of system information as a string in text.
c  Returns also: nt, number of characters set
c                mt, position of last non-blank character
c  Make sure str is dimensioned large enough in the calling program.
c     job=1    returns date and time
c     job=2    returns hostname
c     job=3    returns current user name

      character* (*) text
      character*50 str

C     dummy - returns blank string but right length
      if (job.eq.1) then
        nt=24
      else if (job.eq.2) then
        nt=12
      else if (job.eq.3) then
        nt=8
      else
        nt=0
      endif
      if (nt .gt. 0) text(1:nt)=" "

c  find last non-blank char
      mt=0
      do i=1,nt
        if (text(i:i) .ne. ' ') mt=i
      enddo

c   fill up with blanks
      text(mt+1:)=' '

      end

c -------- klock() --------
      integer function klock()
c  Returns the used cpu time in units of 1/100th of a second.
      klock= 0
      end









