      subroutine zmpy22(a,b,c)
C- Multiplies 2x2 complex matrices a*b
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
       double complex a(2,2), b(2,2), c(2,2)
C ... Local parameters
       double complex cloc(2,2)

       cloc(1,1) = a(1,1)*b(1,1) + a(1,2)*b(2,1)
       cloc(1,2) = a(1,1)*b(1,2) + a(1,2)*b(2,2)
       cloc(2,1) = a(2,1)*b(1,1) + a(2,2)*b(2,1)
       cloc(2,2) = a(2,1)*b(1,2) + a(2,2)*b(2,2)

       c(1,1) = cloc(1,1)
       c(1,2) = cloc(1,2)
       c(2,1) = cloc(2,1)
       c(2,2) = cloc(2,2)

       end

      subroutine zinv22(a,ainv)
C- Inverse of a double complex 2x2 matrix
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
      double complex a(2,2), ainv(2,2)
C ... Local parameters
      double complex det,aloc(2,2)

      det = a(1,1)*a(2,2) - a(1,2)*a(2,1)
      if (dble(det) .eq. 0d0 .and. dimag(det) .eq. 0d0)
     .  call rx('ZINV22: vanishing determinant')

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
C      double complex a(2,2), ainv(2,2)
C
C      a(1,1) = (1d0,2d0)
C      a(1,2) = (2d0,4d0)
C      a(2,1) = (4d0,2d0)
C      a(2,2) = (3d0,5d0)
C
C      call zinv22(a,ainv)
C      write(*,'(''% complex'')')
CC      print 333, dble(a)
CC      print 333, dimag(a)
C      print 333, dble(ainv)
C      print 333, dimag(ainv)
C  333 format(2f15.10)
C
C      end

