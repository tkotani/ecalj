      subroutine fmain
      implicit none
      integer iset,isig
      integer n1,n2,n3
      integer k1,k2,k3
      double precision simple(2,4,4),ans(2,4,4)
      logical lok
      integer w(1000)
      common /w/ w

      call pshpr(0)
      call wkinit(1000)
      call poppr

      call dvset(simple,1,32,.0625d0)
      call dpzero(ans,32)
      ans(1,1,1) = 1d0
      ans(2,1,1) = 1d0
      n1 = 2
      n2 = 2
      n3 = 4
      iset = 0
      isig = 1

      call fftz30(n1,n2,n3,k1,k2,k3)
      call fftz3(simple,n1,n2,n3,k1,k2,k3,1,iset,isig)

      
      print '(a)', ' tfft: check fft for simple 2x2x4 matrix'
      CALL diff('compare fftz3 ', simple, ans, 2*n1, 2*n1, n3, lok )

      end

      SUBROUTINE diff (strn, sc, pc, ldc, n, m, lok )
C- Compare the two arrays for differences
      implicit none
      character*(*) strn
      INTEGER ldc, n, m, i, j
      logical lok
      double precision sc( ldc, * )
      double precision pc( ldc, * )
      
      DO i = 1, n
         DO j = 1, m
            IF ( dABS( sc(i,j) - pc(i,j) ) .GT. 1.0d-6 ) THEN
               WRITE (6,100) strn
               lok = .false.
               RETURN
               END IF
            END DO
        END DO
      WRITE (6,110) strn
      lok = .true.
      RETURN
  100 FORMAT(1X,a,'*** ERROR ***   Arrays Have Different Results')
  110 FORMAT(1X,a,'... arrays have the same results')
      END
