      subroutine cpvprm(mode,i1,i2,iprm,a,aperm)
C- Copy one vector to another in permuted order
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :not referenced for now.
Ci   i1    :copy elements from a into aperm(i1:i2)
Ci   i2    :  "     "      "    "      "   "
Ci   iprm  :permutation table
Ci   a     :source array 
Co Outputs
Co   aperm :denstination array, dimensioned i1:i2
Cr Remarks
Cr   
Cu Updates
Cu   26 Nov 02 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,i1,i2,iprm(i2)
      double precision a(*),aperm(i1:i2)
C ... Local variables
      integer i

      if (mode .eq. 0) then
        do  i = i1, i2
          aperm(i) = a(iprm(i))
        enddo
      else
        do  i = i1, i2
          aperm(iprm(i)) = a(i)
        enddo
      endif

      end

