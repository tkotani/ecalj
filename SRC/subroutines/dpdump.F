      subroutine dpdump(array,length,ifile)
C     - Binary I/O of an array
      integer:: length,ifile
      double precision array(length)
      if (ifile .gt. 0) read(ifile) array
      if (ifile .lt. 0) write(-ifile) array
      end
      subroutine dpsdmp(array,n1,n2,ifile)
C- Binary I/O of an array segment
      integer n1,n2,ifile,length
      double precision array(n2)
      length = n2-n1+1
      if (length .gt. 0) call dpdump(array(n1),length,ifile)
      end
      logical function lddump(array,length,ifile)
C- Binary I/O of an array, returning T if I/O without error or EOF
C     implicit none
      integer length,ifile
      double precision array(length),xx,yy
      lddump = .true.
      if (ifile .gt. 0) then
        yy = array(length)
C       (some random number)
        xx = -1.9283746d0*datan(1d0)
        array(length) = xx
        read(ifile,end=90,err=91) array
        if (xx .ne. array(length)) return
        array(length) = yy
        goto 90
   90   continue
   91   continue
        lddump = .false.
      else
        write(-ifile) array 
      endif
      end
C test lddump
C      subroutine fmain
C      implicit none
C      double precision x(5),y(5)
C      integer ifi,fopna,lddump
C      logical lx
C      data x /1d0,2d0,3d0,4d0,5d0/ y /-1d0,-2d0,-3d0,-4d0,-5d0/
C
C      ifi = fopna('out',-1,4)
C      rewind ifi
C      call dpdump(x,5,-ifi)
C      call dpdump(y,5,-ifi)
C      call dpdump(x,4,-ifi)
C      rewind ifi
C
C      lx = lddump(x,5,ifi)
C      print *, lx
C      lx = lddump(y,5,ifi)
C      print *, lx
C      lx = lddump(y,5,ifi)
C      print *, lx
C      call awrit1('%5:1d',' ',80,6,y)
C      end

