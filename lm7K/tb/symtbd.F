      subroutine symtbd(ng,g,gd)
C- Make transformation table for use in symmetrizing d-orbitals
C ----------------------------------------------------------------------
Ci Inputs
Ci   ng:  number of group operations
Ci   g:  3x3 rotation matrices (Cartesian)
Co Outputs
Co   gd:  transformation table for d-orbitals
Cr Remarks
Cr   Under symmetry operation ig
Cr     d_new(m) = \sum_m' gd(m,m',ig) d(m')
Cr   where the m indices are:
Cr      1    2    3      4          5
Cr     xy   yz   zx   x^2-y^2   3z^2-r^2 
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer ng
      double precision g(3,3,ng),gd(5,5,ng)
C Local variables
      integer ig,m,i,j
      double precision r3,byr3

      r3 = dsqrt(3d0)
      byr3 = 1d0 / r3

      do  20  ig = 1, ng
        do  10  m = 1, 5
          if (m .eq. 1 .or. m .eq. 2 .or. m .eq. 3) then
C --- xy, yz, zx ---
            if (m .eq. 1) then
              i = 1
              j = 2
            elseif (m .eq. 2) then
              i = 2
              j = 3
            elseif (m .eq. 3) then
              i = 3
              j = 1
            endif
            gd(m,1,ig) = g(1,i,ig)*g(2,j,ig) + g(2,i,ig)*g(1,j,ig)
            gd(m,2,ig) = g(2,i,ig)*g(3,j,ig) + g(3,i,ig)*g(2,j,ig)
            gd(m,3,ig) = g(3,i,ig)*g(1,j,ig) + g(1,i,ig)*g(3,j,ig)
            gd(m,4,ig) = g(1,i,ig)*g(1,j,ig) - g(2,i,ig)*g(2,j,ig)
            gd(m,5,ig) = g(3,i,ig)*g(3,j,ig)*r3
          elseif (m .eq. 4) then
C --- x^2-y^2 ---
            gd(m,1,ig) = g(1,1,ig)*g(2,1,ig) - g(2,2,ig)*g(1,2,ig)
            gd(m,2,ig) = g(2,1,ig)*g(3,1,ig) - g(3,2,ig)*g(2,2,ig)
            gd(m,3,ig) = g(3,1,ig)*g(1,1,ig) - g(1,2,ig)*g(3,2,ig)
            gd(m,4,ig) = (g(1,1,ig)*g(1,1,ig) - g(2,1,ig)*g(2,1,ig) 
     .        - g(1,2,ig)*g(1,2,ig) + g(2,2,ig)*g(2,2,ig))*0.5d0
            gd(m,5,ig) = (g(3,1,ig)*g(3,1,ig)
     .                  - g(3,2,ig)*g(3,2,ig))*r3*0.5d0
          elseif (m .eq. 5) then
C --- 3z^2-r^2 ---
            gd(m,1,ig) = (2*g(1,3,ig)*g(2,3,ig) - g(1,1,ig)*g(2,1,ig)
     .                 -  g(2,2,ig)*g(1,2,ig))*byr3
            gd(m,2,ig) = (2*g(2,3,ig)*g(3,3,ig) - g(2,1,ig)*g(3,1,ig)
     .                 -  g(3,2,ig)*g(2,2,ig))*byr3
            gd(m,3,ig) = (2*g(3,3,ig)*g(1,3,ig) - g(3,1,ig)*g(1,1,ig)
     .                 -  g(1,2,ig)*g(3,2,ig))*byr3
            gd(m,4,ig) = (2*g(1,3,ig)*g(1,3,ig) - 2*g(2,3,ig)*g(2,3,ig)
     .        - g(1,1,ig)*g(1,1,ig) + g(2,1,ig)*g(2,1,ig)
     .        - g(1,2,ig)*g(1,2,ig) + g(2,2,ig)*g(2,2,ig))*byr3*0.5d0
            gd(m,5,ig) = (2*g(3,3,ig)*g(3,3,ig) - g(3,1,ig)*g(3,1,ig)
     .                 -    g(3,2,ig)*g(3,2,ig))*0.5d0
          endif
   10   continue
   20 continue

      end
