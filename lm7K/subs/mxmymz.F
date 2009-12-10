      SUBROUTINE MXMYMZ(KIN,K,LX,LY,LZ)
C- Do mirrors in x,y,z if lx,ly,lz=1, respectively
C ----------------------------------------------------------------
Ci Inputs
Ci
Co Outputs
Co
Cr Remarks
Cr
C ----------------------------------------------------------------
C     implicit none
      integer KIN(3),K(3),lx,ly,lz
      K(1) = KIN(1)
      K(2) = KIN(2)
      K(3) = KIN(3)
      IF (LX .EQ. 1) K(1) = 1-K(1)
      IF (LY .EQ. 1) K(2) = 1-K(2)
      IF (LZ .EQ. 1) K(3) = 1-K(3)
      RETURN
      END

