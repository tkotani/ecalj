      subroutine bsortl(vec,nkd)
C     implicit none
C- Sorts a table of vectors by increasing length
      integer nkd
      double precision vec(3,1)
      double precision dsqr,xx
      integer i,is,j
      logical lswap

      dsqr(i) = vec(1,i)**2 + vec(2,i)**2 + vec(3,i)**2

   10 lswap = .false.
      do  20  is = 2, nkd
        if (dsqr(is) .ge. dsqr(is-1)) goto 20
        lswap = .true.
        do  30  j = 1, 3
          xx = vec(j,is)
          vec(j,is) = vec(j,is-1)
          vec(j,is-1) = xx
   30   continue
   20 continue
      if (lswap) goto 10

c      print *, nkd, ' lattice vectors'
c      print 300, (dsqr(is), vec(1,is), vec(2,is), vec(3,is), is=1,nkd)
c  300 format(4f12.5)
c      stop
      end

