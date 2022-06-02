SUBROUTINE MXMYMZ(KIN,K,LX,LY,LZ)
  !- Do mirrors in x,y,z if lx,ly,lz=1, respectively
  ! ----------------------------------------------------------------
  !i Inputs
  !i
  !o Outputs
  !o
  !r Remarks
  !r
  ! ----------------------------------------------------------------
  !     implicit none
  integer :: KIN(3),K(3),lx,ly,lz
  K(1) = KIN(1)
  K(2) = KIN(2)
  K(3) = KIN(3)
  IF (LX == 1) K(1) = 1-K(1)
  IF (LY == 1) K(2) = 1-K(2)
  IF (LZ == 1) K(3) = 1-K(3)
  RETURN
END SUBROUTINE MXMYMZ

