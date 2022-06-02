subroutine dpcopy(afrom,ato,n1,n2,fac)
  !- Copy and scale a portion of a vector
  !     implicit none
  integer :: n1,n2,i
  double precision :: afrom(1),ato(1),fac
  if (fac /= 1d0) goto 100
  call dcopy(n2-n1+1,afrom(n1),1,ato(n1),1)
  return
  ! --- fac not unity ---
100 continue
  do    i = n1, n2
     ato(i) = fac*afrom(i)
  enddo
end subroutine dpcopy
!      program test
!      implicit none
!      double precision from(10), to(10)
!      integer i

!      do  10  i = 1, 10
!      to(i) = 0
!   10 from(i) = i

!      call dpcopy(from,to,3,7,2d0)
!      print 333, to
!      call dpcopy(from,to,3,7,1d0)
!      print 333, to
!  333 format(10f7.3)
!      end

