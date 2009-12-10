      subroutine dinv33(matrix,iopt,invrse,det)
C- Inverts 3x3 matrix
C ----------------------------------------------------------------
Ci Inputs
Ci   matrix:  matrix to be inverted
Ci   iopt:  if 0, usual inverse
Ci             1, transpose of inverse
Ci             2, 2*pi*inverse
Ci             3, 2*pi*transpose of inverse
Co Outputs
Co   invrse    see iopt
Co   det:      determinant, or det/2*pi (sign ok ??)
Cr Remarks
Cr   To generate reciprocal lattice vectors, call dinv33(plat,3,rlat)
C ----------------------------------------------------------------
C     implicit none
      integer iopt,i,j
      double precision matrix(3,3),invrse(3,3),det,ddot
      double precision xx
      call cross(matrix(1,2),matrix(1,3),invrse     )
      call cross(matrix(1,3),matrix     ,invrse(1,2))
      call cross(matrix     ,matrix(1,2),invrse(1,3))
      det = ddot(3,matrix,1,invrse,1)
      if (det .eq. 0d0) call rx('INV33: vanishing determinant')
      if (iopt .ge. 2) det = det/(8*datan(1d0))
      if (mod(iopt,2) .eq. 0) then
        do  10  i = 1, 3
          do  10  j = i+1, 3
          xx = invrse(i,j)
          invrse(i,j) = invrse(j,i)
          invrse(j,i) = xx
   10   continue
      endif
      call dscal(9,1/det,invrse,1)
      end

