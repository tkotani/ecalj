      integer function dinv(cs,n,lda,a)
C- Inversion of a real matrix
C     implicit none
      character*1 cs
      integer n,lda
      double precision a(lda,lda)
      integer ldw,ow,i
      integer w(1)
      common /w/ w

      ldw = n
      call defrr(ow,  ldw*(n+1))
      call dqinv(cs,a,lda,2,n,w(ow),ldw,i)
      call rlse(ow)
      dinv = i
      end

