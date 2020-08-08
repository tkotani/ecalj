      subroutine sortea(ea,ieaord,n,isig)
c mini-sort routine.
      implicit real*8(a-h,o-z)
      real(8)::        ea(n)
      integer(4):: ieaord(n)
c sorting of ea
      isig = 1
*poption noparallel
      do i = 1,n
        ieaord(i) = i
      enddo
*poption noparallel
      do ix= 2,n
*poption noparallel
        do i=ix,2,-1
          if( ea(ieaord(i-1)) >ea(ieaord(i) ) ) then
#ifdef EXPAND_ISWAP
            itmp = ieaord(i-1)
            ieaord(i-1) = ieaord(i)
            ieaord(i) = itmp

#else
            call iswap (ieaord(i-1),ieaord(i))
#endif
            isig= -isig
            cycle
          endif
          exit
        enddo
      enddo
      end
#ifndef EXPAND_ISWAP
      subroutine iswap(i,j)
      iwork = j
      j = i
      i = iwork
      end
#endif








