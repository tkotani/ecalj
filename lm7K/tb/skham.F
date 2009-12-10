      subroutine skham(dl,dm,dn,V,nl,ndimL,h)
C- Slater-Koster couplings from matrix elements and direction cosines
C ----------------------------------------------------------------------
Ci Inputs
Ci   dl,dm,dn: direction cosines
Ci   V: sss,sps,pps,ppp,sds,pds,pdp,dds,ddp,ddd (matrix elements)
Ci   nl: 1,2 or 3 for s, sp, spd; ndimL: first dimension of h
Co Outputs
Co   h: slater-koster hamiltonian
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
      double precision dl2,dm2,dn2,dlm,dmn,dnl,w3
      integer il(9),nl2,j,k
      data il/1,-1,-1,-1,1,1,1,1,1/

      w3 = dsqrt(3d0)

      dl2 = dl**2
      dm2 = dm**2
      dn2 = dn**2
      dlm = dl*dm
      dmn = dm*dn
      dnl = dn*dl

C     --- ss block ---
      h(1,1) = V(1)

      if (nl .lt. 2) goto 10

C     --- sp block ---
      h(1,2) = dl*V(2)
      h(1,3) = dm*V(2)
      h(1,4) = dn*V(2)

C     --- pp block ---
C     x,x  y,y  z,z ...
      h(2,2) = dl2*V(3) + (1-dl2)*V(4)
      h(3,3) = dm2*V(3) + (1-dm2)*V(4)
      h(4,4) = dn2*V(3) + (1-dn2)*V(4)
C     x,y  y,z  z,x ...
      h(2,3) = dlm*(V(3) - V(4))
      h(3,4) = dmn*(V(3) - V(4))
      h(2,4) = dnl*(V(3) - V(4))

      if (nl .lt. 3) goto 10
      
C     --- sd block ---
      h(1,5) = w3*dlm*V(5)
      h(1,6) = w3*dmn*V(5)
      h(1,7) = w3*dnl*V(5)
      h(1,8) = .5d0*w3*(dl2 - dm2)*V(5)
      h(1,9) = (dn2-.5d0*(dl2 + dm2))*V(5)

C     --- pd block ---
      h(2,6) = dlm*dn*(w3*V(6) - 2*V(7))
      h(3,7) = h(2,6)
      h(4,5) = h(2,6)

      h(2,7) = w3*dl2*V(6) + (1 - 2*dl2)*V(7)
      h(2,5) = h(2,7)*dm
      h(2,7) = h(2,7)*dn
      h(3,5) = w3*dm2*V(6) + (1 - 2*dm2)*V(7)
      h(3,6) = h(3,5)*dn
      h(3,5) = h(3,5)*dl
      h(4,6) = w3*dn2*V(6) + (1 - 2*dn2)*V(7)
      h(4,7) = h(4,6)*dl
      h(4,6) = h(4,6)*dm

      h(2,9) = (dn2-.5d0*(dl2 + dm2))*V(6)
      h(3,9) = h(2,9) - w3*dn2*V(7)
      h(4,9) = h(2,9) + w3*(dl2 + dm2)*V(7)
      h(2,9) = h(3,9)*dl
      h(3,9) = h(3,9)*dm
      h(4,9) = h(4,9)*dn

      h(2,8)= .5d0*w3*(dl2 - dm2)*V(6) - (dl2-dm2)*V(7)
      h(3,8) = (h(2,8) - V(7))*dm
      h(4,8) = h(2,8)*dn
      h(2,8) = (h(2,8) + V(7))*dl

C     --- dd block ---
C     xy,xy  yz,yz  zx,zx ...
      h(5,5) = 3*dl2*dm2*V(8) +
     .  (dl2 + dm2 - 4*dl2*dm2)*V(9) + (dn2 + dl2*dm2)*V(10)
      h(6,6) = 3*dm2*dn2*V(8) +
     .  (dm2 + dn2 - 4*dm2*dn2)*V(9) + (dl2 + dm2*dn2)*V(10)
      h(7,7) = 3*dn2*dl2*V(8) +
     .  (dn2 + dl2 - 4*dn2*dl2)*V(9) + (dm2 + dn2*dl2)*V(10)
C     xy,yz  yz,zx  zx,xy ...
      h(5,6) = (3*dm2*V(8) + (1 - 4*dm2)*V(9) + (dm2 - 1)*V(10))*dnl
      h(6,7) = (3*dn2*V(8) + (1 - 4*dn2)*V(9) + (dn2 - 1)*V(10))*dlm
      h(5,7) = (3*dl2*V(8) + (1 - 4*dl2)*V(9) + (dl2 - 1)*V(10))*dmn
C     xy,3z^2-1  yz,3z^2-1  zx,3z^2-1 ...
      h(5,9) = (dn2 - .5d0*(dl2 + dm2))*V(8)
      h(6,9) = h(5,9)
      h(5,9) = h(5,9) - 2.d0*dn2*V(9) + .5d0*(1 + dn2)*V(10)
      h(5,9) = h(5,9)*dlm*w3
      h(6,9) = h(6,9) + (dl2 + dm2 - dn2)*V(9) - .5d0*(dl2 + dm2)*V(10)
      h(7,9) = h(6,9)*dnl*w3
      h(6,9) = h(6,9)*dmn*w3
C     xy,x^2-y^2  yz,x^2-y^2   zx,x^2-y^2 ...
      h(5,8) = (dl2 - dm2)*(1.5d0*V(8) - 2.d0*V(9) + 0.5d0*V(10))
      h(6,8) = h(5,8) - V(9) + V(10)
      h(7,8) = h(5,8) + V(9) - V(10)
      h(5,8) = h(5,8)*dlm
      h(6,8) = h(6,8)*dmn
      h(7,8) = h(7,8)*dnl
C     3x^2-1 ...
      h(9,9) = (dn2 - 0.5d0*(dl2 + dm2))**2*V(8) +
     .  3*dn2*(dl2 + dm2)*V(9) + .75d0*(dl2 + dm2)**2*V(10)
C     x^2-y^2 ...
      h(8,8) = .75d0*(dl2 - dm2)**2*V(8) +
     .  (dl2 + dm2 - (dl2 - dm2)**2)*V(9) +
     .  (dn2 + .25d0*(dl2 - dm2)**2)*V(10)
C     x^2-y^2  3x^2-1 ...
      h(8,9) = .5d0*w3*(dl2 - dm2)*(dn2 - .5d0*(dl2 + dm2))*V(8) +
     .  w3*dn2*(dm2 - dl2)*V(9) +
     .  .25d0*w3*(1 + dn2)*(dl2 - dm2)*V(10)

   10 continue
      nl2 = nl**2
      do  20  k = 1, nl2-1
        do  20  j = k+1, nl2
        h(j,k) = h(k,j)*il(k)*il(j)
   20 continue

      end
