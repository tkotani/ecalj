      subroutine dfdump(array,length,ifile)
C- ASCII I/O of an array
C     implicit none
      integer ifile,length
      double precision array(length)
      if (ifile .gt. 0) read(ifile,333) array
      if (ifile .lt. 0) write(-ifile,333) array
  333 format(1p,4e20.13)
      end
      subroutine dfsdmp(array,n1,n2,ifile)
C- ASCII I/O of an array segment
      integer n1,n2,ifile,length
      double precision array(n2)
      length = n2-n1+1
      if (length .gt. 0) call dfdump(array(n1),length,ifile)
      end
      logical function lfdmp(array,length,ifile)
C- ASCII I/O of an array, returning T if I/O without error or EOF
C     implicit none
      integer length,ifile
      double precision array(length)

      lfdmp = .true.
      if (ifile .gt. 0) read(ifile,333,end=90,err=90) array
      if (ifile .lt. 0) write(-ifile,333) array
  333 format(1p,4e20.13)
      return

   90 lfdmp = .false.
      end

