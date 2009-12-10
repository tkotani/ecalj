      subroutine gskham(dl,dm,dn,V,nl,ndimL,h)
C- Extension to Slater-Koster couplings from Srini's extra term
C ----------------------------------------------------------------------
Ci Inputs
Ci   dl,dm,dn: direction cosines
Ci   V: sss,sps,pps,ppp,usp,uxy,upps,uppp
Ci   nl: 1,2 or 3 for s, sp, spd; ndimL: first dimension of h
Co Outputs
Co   h: Addition to Slater-Koster hamiltonian
Cr Remarks
Cr  Calculates the structure constant matrix h using the s-s-sigma
Cr  etc. constants and the direction cosines following Koster and Slater
Cr  Phys. Rev. 94, 1498 table 1
Cr         1     2     3     4     5     6     7     8     9
Cr         s     x     y     z    xy    yz    zx  x^2-y^2   3x^2-1 
Cr  V has only 10 elements, which assumes symmetry in sp ps, sd ds etc..
Cr  If these matrix elements are different the hamiltonian is 
Cr  symmetrised later. See swapV
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nl,ndimL
      double precision dl,dm,dn,h(ndimL,9),V(10)
C Local variables
      double precision dlm,dmn,dnl

      dlm = dl*dm
      dmn = dm*dn
      dnl = dl*dn

      h(1,2) = h(1,2) + dmn*V(5)
      h(1,3) = h(1,3) + dnl*V(5)
      h(1,4) = h(1,4) + dlm*V(5)
      h(2,1) = h(2,1) - dmn*V(5)
      h(3,1) = h(3,1) - dnl*V(5)
      h(4,1) = h(4,1) - dlm*V(5)

      h(2,3) = h(2,3) + dlm*V(6)
      h(2,4) = h(2,4) + dnl*V(6)
      h(3,4) = h(3,4) + dmn*V(6)
      h(3,2) = h(3,2) + dlm*V(6)
      h(4,2) = h(4,2) + dnl*V(6)
      h(4,3) = h(4,3) + dmn*V(6)

      h(2,2) = h(2,2) + dl*dl*(V(7)-V(8))
      h(3,3) = h(3,3) + dm*dm*(V(7)-V(8))
      h(4,4) = h(4,4) + dn*dn*(V(7)-V(8))

      end

