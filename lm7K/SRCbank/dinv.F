c takao modified this based on delw_2.2.
      integer function dinv(cs,n,lda,a)
C- Inversion of a real matrix
C     implicit none
      character*1 cs
      integer n,lda
      double precision a(lda,lda)
      integer:: ldw , i 
      real(8) ,allocatable :: work(:)
      allocate(work(ldw*(n+1)))
      call dqinv ( cs , a , lda , 2 , n , work, ldw , i ) 
      deallocate(work)
      dinv = i
      end

