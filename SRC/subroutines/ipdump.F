      subroutine ipdump(array,length,ifile)
C- Binary I/O of an integer array
      integer array(length)
      if (ifile .gt. 0) read(ifile) array
      if (ifile .lt. 0) write(-ifile) array
      end
      subroutine ipsdmp(array,n1,n2,ifile)
C- Binary I/O of an integer array segment
      integer n1,n2,ifile,length
      integer array(n2)
      length = n2-n1+1
      if (length .gt. 0) call ipdump(array(n1),length,ifile)
      end
      logical function lidump(array,length,ifile)
C- Binary I/O of an array, returning T if I/O without error or EOF
C     implicit none
      integer length,ifile
      integer array(length),xx,yy
      lidump = .true.
      if (ifile .gt. 0) then
        yy = array(length)
C       (some random number)
        xx = -142567212
        array(length) = xx
        read(ifile,end=90,err=91) array
        if (xx .ne. array(length)) return
        array(length) = yy
        goto 90
   90   continue
   91   continue
        lidump = .false.
      else
        write(-ifile) array 
      endif
      end

