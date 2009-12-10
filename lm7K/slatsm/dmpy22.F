       subroutine dmpy22(a,b,c)
C- Multiplies 2x2 matrices a*b
C ----------------------------------------------------------------------
Ci Inputs
Ci   a     :left matrix
Ci   b     :right matrix
Co Outputs
Co   c     :Product a*b
Cr Remarks
Cr   It is permissible for any of a,b,c to use the same address space
Cu Updates
Cu   17 Mar 03 First created (from A. Chantis)
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
       double precision a(2,2), b(2,2), c(2,2)
C ... Local parameters
       double precision cloc(2,2)

       cloc(1,1) = a(1,1)*b(1,1) + a(1,2)*b(2,1)
       cloc(1,2) = a(1,1)*b(1,2) + a(1,2)*b(2,2)
       cloc(2,1) = a(2,1)*b(1,1) + a(2,2)*b(2,1)
       cloc(2,2) = a(2,1)*b(1,2) + a(2,2)*b(2,2)

       c(1,1) = cloc(1,1)
       c(1,2) = cloc(1,2)
       c(2,1) = cloc(2,1)
       c(2,2) = cloc(2,2)

       end

       subroutine dinv22(a,ainv)
C- Inverse of a real 2x2 matrix
C ----------------------------------------------------------------------
Ci Inputs
Ci   a     :the matrix to be inverted
Co Outputs
Co   ainv  :Inverse of a
Cr Remarks
Cr   It is permissible for a and ainv to use the same address space
Cu Updates
Cu   17 Mar 03 First created (from A. Chantis)
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      double precision a(2,2), ainv(2,2)
C ... Local parameters
      double precision det,aloc(2,2)

      det = a(1,1)*a(2,2) - a(1,2)*a(2,1)
      if (det .eq. 0d0) call rx('INV22: vanishing determinant')

      aloc(1,1) = a(2,2)/det
      aloc(2,2) = a(1,1)/det
      aloc(1,2) = -a(1,2)/det
      aloc(2,1) = -a(2,1)/det

      ainv(1,1) = aloc(1,1)
      ainv(2,2) = aloc(2,2)
      ainv(1,2) = aloc(1,2)
      ainv(2,1) = aloc(2,1)
      end

C     Test
C      subroutine fmain
C
C      double precision a(2,2), ainv(2,2)
C
C      a(1,1) = 1d0
C      a(1,2) = 2d0
C      a(2,1) = 4d0
C      a(2,2) = 3d0
C
C      call dinv22(a,ainv)
C      write(*,'(''% real'')')
C      print 333, dble(a)
C      print 333, dble(ainv)
C  333 format(2f15.10)
C
C      end

