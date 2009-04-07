      double precision function volsph(nspec,nrspec,wsr)
C- Sum-of-sphere volumes
C ----------------------------------------------------------------------
Ci Inputs:
Ci  nspec :number of species
Ci  nrspec:number of atoms in the ith species
Ci  wsr   :Wigner-Seitz sphere radius
Co Outputs:
Co  volsph:sum-of-sphere volumes, in units of wsr**3
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer nspec,nrspec(*)
      double precision wsr(*)
C Local variables:
      integer ic
      double precision fpi3
      parameter(fpi3=4.18879020478639053d0)

      volsph = 0d0
      do  ic = 1, nspec
        volsph = volsph + fpi3*nrspec(ic)*wsr(ic)**3
      enddo
      end

