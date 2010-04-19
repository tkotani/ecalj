      subroutine skhso(nl,h)
C- Slater-Koster spin-orbit couplings
C ----------------------------------------------------------------------
Ci Inputs
Ci   nl: 1, 2, or 3 for s, sp, spd
Co Outputs
Co   h: slater-koster spin-orbit hamiltonian for one site
Cr Remarks
Cr  Calculates the spin-orbit matrix h_so
Cr         1     2     3     4     5     6     7     8         9
Cr         s     x     y     z    xy    yz    zx  x^2-y^2   3z^2-r^2
Cr  h_so is zero for s-states
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nl
      double precision h(nl**2,nl**2,4,2)
C Local variables
      double precision w3
      integer nl2,j,k

      if (nl .lt. 2) return

      w3 = dsqrt(3d0)

C     --- p block ---
C     spin: + +
      h(2,3,1,2) = -1d0
C     spin: + -
      h(2,4,2,1) =  1d0
      h(3,4,2,2) = -1d0

      if (nl .lt. 3) goto 10
      
C     --- d block ---
C     spin: + +
      h(5,8,1,2) =  2d0
      h(6,7,1,2) =  1d0
C     spin: + -
      h(5,6,2,1) =  1d0
      h(5,7,2,2) = -1d0
      h(6,8,2,2) = -1d0
      h(6,9,2,2) = -w3
      h(7,8,2,1) = -1d0
      h(7,9,2,1) =  w3

   10 continue

C     spin: lower half of + + and + -; all of - + and - -
      nl2 = nl**2
      do  20  k = 2, nl2
        do  20  j = 2, nl2
          if (j .gt. k) then
            h(j,k,1,2) = -h(k,j,1,2)
            h(j,k,2,1) = -h(k,j,2,1)
            h(j,k,2,2) = -h(k,j,2,2)
          endif
          h(j,k,3,1) =  h(k,j,2,1)
          h(j,k,3,2) = -h(k,j,2,2)
          h(j,k,4,2) = -h(j,k,1,2)
   20 continue

      end

