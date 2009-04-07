      subroutine dskham(l,m,n,r,V,dV,nl,ndimL,cryf,dhx,dhy,dhz,dhr)
C- space derivatives of Slater-Koster hamiltonian
C ----------------------------------------------------------------------
Ci Inputs
Ci   l,m,n: direction cosines; r: bond length
Ci   V: sss,sps,pps,ppp,sds,pds,pdp,dds,ddp,ddd (matrix elements)
Ci  dV: radial derivatives of V. V and dV evaluated at r
Ci  nl: 1,2 or 3 for s, sp, spd; ndimL: first dimension of dhx,y,z
Ci  cryf: if true then crystal field derivatives
Co Outputs
Co   dh(x,y,z): slater-koster hamiltonian derivative; dhr: radial deriv.
Cr         1     2     3     4     5     6     7     8       9
Cr         s     x     y     z    xy    yz    zx  x^2-y^2   3x^2-1 
Cr  V  and dV have only 10 elements, which assumes symmetry in 
Cr  sp ps, sd ds etc.. If these matrix elements are different the 
Cr  hamiltonian gradient is symmetrised later. See swapdV
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nl,ndimL
      double precision l,m,n,r,V(10),dV(10),dhx(ndimL,9),dhy(ndimL,9),
     .  dhz(ndimL,9),dhr(ndimL,9)
      logical cryf
C Local
      double precision r3,wk(9,9),l2,m2,n2,lmn,l2m2,n2m2,n2l2,l2m2mi,
     .  l2m2pl,l2n2pl,n2m2pl
      integer il(9),i,j,nl2
      data il/1,-1,-1,-1,1,1,1,1,1/

      nl2 = nl*nl
      l2=l*l
      m2=m*m
      n2=n*n
      lmn=l*m*n
      l2m2=l2*m2
      n2m2=n2*m2
      n2l2=n2*l2
      l2m2mi=l2-m2
      l2m2pl=l2+m2
      l2n2pl=l2+n2
      n2m2pl=n2+m2

      r3 = dsqrt(3d0)

C ---- Angular derivatives ----

C     --- ss block ---
      dhx(1,1) = 0d0
      dhy(1,1) = 0d0
      dhz(1,1) = 0d0

      if (nl .lt. 2) goto 10000

C     --- sp block ---
      dhx(1,2) = V(2)*(l*l - 1d0)
      dhy(1,2) = V(2)*l*m
      dhz(1,2) = V(2)*l*n
      dhx(1,3) = V(2)*m*l
      dhy(1,3) = V(2)*(m*m - 1d0)
      dhz(1,3) = V(2)*m*n
      dhx(1,4) = V(2)*n*l
      dhy(1,4) = V(2)*n*m
      dhz(1,4) = V(2)*(n*n - 1d0)

C     --- pp block ---
      dhx(2,2) = (V(3) - V(4))*2d0*l*(l*l - 1d0)
      dhy(2,2) = (V(3) - V(4))*2d0*l*l*m
      dhz(2,2) = (V(3) - V(4))*2d0*l*l*n
      dhx(3,3) = (V(3) - V(4))*2d0*m*m*l
      dhy(3,3) = (V(3) - V(4))*2d0*m*(m*m - 1d0)
      dhz(3,3) = (V(3) - V(4))*2d0*m*m*n
      dhx(4,4) = (V(3) - V(4))*2d0*n*n*l
      dhy(4,4) = (V(3) - V(4))*2d0*n*n*m
      dhz(4,4) = (V(3) - V(4))*2d0*n*(n*n - 1d0)
      
      dhx(2,3) = (V(3) - V(4))*(2d0*l*m*l - m)
      dhy(2,3) = (V(3) - V(4))*(2d0*l*m*m - l)
      dhz(2,3) = (V(3) - V(4))*2d0*l*m*n
      dhx(3,4) = (V(3) - V(4))*2d0*m*n*l
      dhy(3,4) = (V(3) - V(4))*(2d0*m*n*m - n)
      dhz(3,4) = (V(3) - V(4))*(2d0*m*n*n - m)
      dhx(2,4) = (V(3) - V(4))*(2d0*l*n*l - n)
      dhy(2,4) = (V(3) - V(4))*2d0*n*m*l
      dhz(2,4) = (V(3) - V(4))*(2d0*l*n*n - l)

      if (nl .lt. 3) goto 10000

C     --- sd block ---
      dhx(1,5) = r3*V(5)*(2d0*l*l*m - m)
      dhy(1,5) = r3*V(5)*(2d0*m*l*m - l)
      dhz(1,5) = r3*V(5)*2d0*n*l*m
      dhx(1,6) = dhz(1,5)
      dhy(1,6) = r3*V(5)*(2d0*m*m*n - n)
      dhz(1,6) = r3*V(5)*(2d0*n*m*n - m)
      dhx(1,7) = r3*V(5)*(2d0*l*n*l - n)
      dhy(1,7) = dhz(1,5)
      dhz(1,7) = r3*V(5)*(2d0*n*n*l - l)
      dhx(1,8) = r3*V(5)*(l*(l*l-m*m) - l)
      dhy(1,8) = r3*V(5)*(m*(l*l-m*m) + m)
      dhz(1,8) = r3*V(5)*n*(l*l-m*m)
      dhx(1,9) = V(5)*(2d0*l*(n*n - (l*l + m*m)/2d0) + l)
      dhy(1,9) = V(5)*(2d0*m*(n*n - (l*l + m*m)/2d0) + m)
      dhz(1,9) = V(5)*(2d0*(n*(n*n - (l*l + m*m)/2d0) - n))

C     --- pd block ---
C x,xy :      
      dhx(2,5) = r3*(l*3d0*l*l*m - 2d0*l*m)*V(6)
     .          + (l*m*(1d0 - 6d0*l*l) + 4d0*l*m)*V(7)
      dhy(2,5) = r3*(m*3d0*l*l*m - l*l)*V(6) 
     .          + (m*m*(1d0 - 6d0*l*l) - (1d0 - 2d0*l*l))*V(7)
      dhz(2,5) = r3*n*3d0*l*l*m*V(6) + n*m*(1d0 - 6d0*l*l)*V(7)
C y,yz :
      dhx(3,6) = r3*l*3d0*m*m*n*V(6) + l*n*(1d0 - 6d0*m*m)*V(7)
      dhy(3,6) = r3*(m*3d0*m*m*n - 2d0*m*n)*V(6)
     .          + (m*n*(1d0 - 6d0*m*m) + 4d0*m*n)*V(7)
      dhz(3,6) = r3*(n*3d0*m*m*n - m*m)*V(6)
     .          + (n*n*(1d0 - 6d0*m*m) - (1d0 - 2d0*m*m))*V(7)
C z,zx :
      dhx(4,7) = r3*(l*3d0*n*n*l - n*n)*V(6)
     .          + (l*l*(1d0 - 6d0*n*n) - (1d0 - 2d0*n*n))*V(7)
      dhy(4,7) = r3*m*3d0*n*n*l*V(6) + m*l*(1d0 - 6d0*n*n)*V(7)
      dhz(4,7) = r3*(n*3d0*n*n*l - 2d0*n*l)*V(6)
     .          + (n*l*(1d0 - 6d0*n*n) + 4d0*n*l)*V(7)
C x,zx :
      dhx(2,7) = r3*(l*3d0*l*l*n - 2d0*l*n)*V(6)
     .          + (l*n*(1d0 - 6d0*l*l) + 4d0*l*n)*V(7)
      dhy(2,7) = r3*m*3d0*l*l*n*V(6) + m*n*(1d0 - 6d0*l*l)*V(7)
      dhz(2,7) = r3*(n*3d0*l*l*n - l*l)*V(6) 
     .          + (n*n*(1d0 - 6d0*l*l) - (1d0 - 2d0*l*l))*V(7)
C y,xy :
      dhx(3,5) = r3*(l*3d0*m*m*l - m*m)*V(6)
     .          + (l*l*(1d0 - 6d0*m*m) - (1d0 - 2d0*m*m))*V(7)
      dhy(3,5) = r3*(m*3d0*m*m*l - 2d0*m*l)*V(6)
     .          + (m*l*(1d0 - 6d0*m*m) + 4d0*m*l)*V(7)
      dhz(3,5) = r3*n*3d0*m*m*l*V(6) + n*l*(1d0 - 6d0*m*m)*V(7)
C z,yz :
      dhx(4,6) = r3*l*3d0*n*n*m*V(6) + l*m*(1d0 - 6d0*n*n)*V(7)
      dhy(4,6) = r3*(m*3d0*n*n*m - n*n)*V(6)
     .          + (m*m*(1d0 - 6d0*n*n) - (1d0 - 2d0*n*n))*V(7)
      dhz(4,6) = r3*(n*3d0*n*n*m - 2d0*n*m)*V(6)
     .          + (n*m*(1d0 - 6d0*n*n) + 4d0*n*m)*V(7)
C x,yz y,zx z,xy :
      dhx(2,6) = (l*3d0*l*m*n - m*n)*(r3*V(6) - 2d0*V(7))
      dhy(2,6) = (m*3d0*l*m*n - l*n)*(r3*V(6) - 2d0*V(7))
      dhz(2,6) = (n*3d0*l*m*n - l*m)*(r3*V(6) - 2d0*V(7))
      dhx(3,7) = dhx(2,6)
      dhy(3,7) = dhy(2,6)
      dhz(3,7) = dhz(2,6)
      dhx(4,5) = dhx(2,6)
      dhy(4,5) = dhy(2,6)
      dhz(4,5) = dhz(2,6)
C x,x^2-y^2 :
      dhx(2,8) = (l*3d0*l*(m*m - l*l) - m*m + 3d0*l*l)
     .           *(V(7) - r3*V(6)/2d0) + (l*l - 1d0)*V(7)
      dhy(2,8) = (m*3d0*l*(m*m - l*l) - 2d0*l*m)*(V(7) - r3*V(6)/2d0)
     .          + m*l*V(7)
      dhz(2,8) = n*3d0*l*(m*m - l*l)*(V(7) - r3*V(6)/2d0) + n*l*V(7)
C y,x^2-y^2 :
      dhx(3,8) = (l*3d0*m*(l*l - m*m) - 2d0*m*l)*(r3*V(6)/2d0 - V(7))
     .          - l*m*V(7)
      dhy(3,8) = (m*3d0*m*(l*l - m*m) - l*l + 3d0*m*m)
     .           *(r3*V(6)/2d0 - V(7)) - (m*m - 1d0)*V(7)
      dhz(3,8) = n*3d0*m*(l*l - m*m)*(r3*V(6)/2d0 - V(7)) - n*m*V(7)
C z,x^2-y^2
      dhx(4,8) = (l*3d0*n*(l*l - m*m) - 2d0*l*n)*(r3*V(6)/2d0 - V(7))
      dhy(4,8) = (m*3d0*n*(l*l - m*m) + 2d0*m*n)*(r3*V(6)/2d0 - V(7))
      dhz(4,8) = (n*3d0*n*(l*l - m*m)
     .                             - (l*l - m*m))*(r3*V(6)/2d0 - V(7))
C x,3z^2-1
      dhx(2,9) = (l*3d0*l*n*n - n*n)*(V(6) - r3*V(7))
     .          - (l*3d0*l*(l*l + m*m)/2d0 - (3d0*l*l + m*m)/2d0)*V(6)
      dhy(2,9) = m*3d0*l*n*n*(V(6) - r3*V(7))
     .          - (m*3d0*l*(l*l + m*m)/2d0 - l*m)*V(6)
      dhz(2,9) = (n*3d0*l*n*n - 2d0*l*n)*(V(6) - r3*V(7))
     .          - n*3d0*l*(l*l + m*m)/2d0*v(6)
C y,3z^2-1
      dhx(3,9) = dhy(2,9)
      dhy(3,9) = (m*3d0*m*n*n - n*n)*(V(6) - r3*V(7))
     .          - (m*3d0*m*(l*l + m*m)/2d0 - (l*l + 3d0*m*m)/2d0)*V(6)
      dhz(3,9) = (n*3d0*m*n*n - 2d0*m*n)*(V(6) - r3*V(7))
     .          - n*3d0*m*(l*l + m*m)/2d0*V(6)
C z,3z^2-1
      dhx(4,9) = l*3d0*n*n*n*V(6)
     .            - (l*3d0*n*(l*l + m*m) - 2d0*n*l)*(V(6)/2d0 - r3*V(7))
      dhy(4,9) = m*3d0*n*n*n*V(6)
     .            - (m*3d0*n*(l*l + m*m) - 2d0*n*m)*(V(6)/2d0 - r3*V(7))
      dhz(4,9) = (n*3d0*n*n*n - 3d0*n*n)*V(6)
     .        - (n*3d0*n*(l*l + m*m) - (l*l + m*m))*(V(6)/2d0 - r3*V(7))

C     --- dd block (taken from Mike Finnis) ---
      dhx(5,5) = l*(m2*(6d0 - 12d0*l2)*V(8)
     .   + (2d0 - 2d0*l2 - 10d0*m2 + 16d0*l2m2)*V(9)
     .   - (2d0 - 2d0*l2 - 4d0*m2 + 4d0*l2m2)*V(10))
      dhy(5,5) = m*(l2*(6d0 - 12d0*m2)*V(8)
     .   + (2d0 - 2d0*m2 - 10d0*l2 + 16d0*l2m2)*V(9)
     .   - (2d0 - 2d0*m2 - 4d0*l2 + 4d0*l2m2)*V(10))
      dhz(5,5) = n*( - 12d0*l2m2*V(8) - (2d0*l2m2pl - 16d0*l2m2)*V(9)
     .   - (4d0*l2m2 - 2d0*l2m2pl)*V(10))
      dhx(5,6) = n*(m2*(3d0 - 12d0*l2)*V(8)
     .   + (1d0 - 2d0*l2 - 4d0*m2 + 16d0*l2m2)*V(9)
     .   - (1d0 - 2d0*l2 - m2 + 4d0*l2m2)*V(10))
      dhy(5,6) = lmn*((6d0 - 12d0*m2)*V(8) - (10d0 - 16d0*m2)*V(9)
     .   + (4d0 - 4d0*m2)*V(10))
      dhz(5,6) = l*(m2*(3d0 - 12d0*n2)*V(8)
     .   + (1d0 - 2d0*n2 - 4d0*m2 + 16d0*n2m2)*V(9)
     .   - (1d0 - 2d0*n2 - m2 + 4d0*n2m2)*V(10))
      dhx(5,7) = lmn*((6d0 - 12d0*l2)*V(8)
     .   - (10d0 - 16d0*l2)*V(9) + (4d0 - 4d0*l2)*V(10))
      dhy(5,7) = n*(l2*(3d0 - 12d0*m2)*V(8)
     .   + (1d0 - 2d0*m2 - 4d0*l2 + 16d0*l2m2)*V(9)
     .   - (1d0 - 2d0*m2 - l2 + 4d0*l2m2)*V(10))
      dhz(5,7) = m*(l2*(3d0 - 12d0*n2)*V(8)
     .   + (1d0 - 2d0*n2 - 4d0*l2 + 16d0*n2l2)*V(9)
     .   - (1d0 - 2d0*n2 - l2 + 4d0*n2l2)*V(10))
      dhx(5,8) = m*((1.5d0*(3d0*l2 - m2) - 6d0*l2*l2m2mi)*V(8)
     .   + (2d0*m2 - 6d0*l2 + 8d0*l2*l2m2mi)*V(9)
     .   + (1.5d0*l2 - 0.5d0*m2 - 2d0*l2*l2m2mi)*V(10))
      dhy(5,8) =  - l*((1.5d0*(3d0*m2 - l2) + 6d0*m2*l2m2mi)*V(8)
     .   + (2d0*l2 - 6d0*m2 - 8d0*m2*l2m2mi)*V(9)
     .   + (1.5d0*m2 - 0.5d0*l2 + 2d0*m2*l2m2mi)*V(10))
      dhz(5,8) =  - l2m2mi*lmn*(6d0*V(8) - 8d0*V(9) + 2d0*V(10))
      dhx(5,9) = r3*m*((1d0 - 1.5d0*m2 + l2*(6d0*l2m2pl - 6.5d0))*V(8)
     .   - 2d0*(1d0 - m2 - l2*(5d0 - 4d0*l2m2pl))*V(9)
     .   + (1d0 - 0.5d0*m2 - l2*(1.5d0 + 2d0*n2))*V(10))
      dhy(5,9) = r3*l*((1d0 - 1.5d0*l2 + m2*(6d0*l2m2pl - 6.5d0))*V(8)
     .   - 2d0*(1d0 - l2 - m2*(5d0 - 4d0*l2m2pl))*V(9)
     .   + (1d0 - 0.5d0*l2 - m2*(1.5d0 + 2d0*n2))*V(10))
      dhz(5,9) =  - lmn*r3*((2d0 - 6d0*l2m2pl)*V(8)
     .   - (4d0 - 8d0*l2m2pl)*V(9) + 2*n2*V(10))
      dhx(6,6) = l*( - 12d0*n2m2*V(8) - (2d0*n2m2pl - 16d0*n2m2)*V(9)
     .   - (4d0*n2m2 - 2d0*n2m2pl)*V(10))
      dhy(6,6) = m*(n2*(6d0 - 12d0*m2)*V(8)
     .   + (2d0 - 2d0*m2 - 10d0*n2 + 16d0*n2m2)*V(9)
     .   - (2d0 - 2d0*m2 - 4d0*n2 + 4d0*n2m2)*V(10))
      dhz(6,6) = n*(m2*(6d0 - 12d0*n2)*V(8)
     .   + (2d0 - 2d0*n2 - 10d0*m2 + 16d0*n2m2)*V(9)
     .   - (2d0 - 2d0*n2 - 4d0*m2 + 4d0*n2m2)*V(10))
      dhx(6,7) = m*(n2*(3d0 - 12d0*l2)*V(8)
     .   + (1d0 - 2d0*l2 - 4d0*n2 + 16d0*n2l2)*V(9)
     .   - (1d0 - 2d0*l2 - n2 + 4d0*n2l2)*V(10))
      dhy(6,7) = l*(n2*(3d0 - 12d0*m2)*V(8)
     .   + (1d0 - 2d0*m2 - 4d0*n2 + 16d0*n2m2)*V(9)
     .   - (1d0 - 2d0*m2 - n2 + 4d0*n2m2)*V(10))
      dhz(6,7) = lmn*((6d0 - 12d0*n2)*V(8) - (10d0 - 16d0*n2)*V(9)
     .   + (4d0 - 4d0*n2)*V(10))
      dhx(6,8) = lmn*((3d0 - 6d0*l2m2mi)*V(8)
     .   - (2d0 - 8d0*l2m2mi)*V(9) - (1d0 + 2d0*l2m2mi)*V(10))
      dhy(6,8) = n*((1.5d0*l2 - m2*(4.5d0 + 6d0*l2m2mi))*V(8)
     .   - (1d0 + 2d0*l2 - 8d0*m2*(1d0 + l2m2mi))*V(9)
     .   + (1d0 + 0.5d0*l2 - m2*(3.5d0 + 2d0*l2m2mi))*V(10))
      dhz(6,8) = m*( - l2m2mi*(4.5d0 - 6d0*l2m2pl)*V(8)
     .   + (4d0*n2*l2m2mi + (1d0 + 2d0*l2m2mi)*(1d0 - 2d0*l2m2pl))*V(9)
     .   - (n2*l2m2mi + (1d0 + 0.5d0*l2m2mi)*(1d0 - 2d0*l2m2pl))*V(10))
      dhx(6,9) =  - lmn*r3*((5d0 - 6d0*l2m2pl)*V(8)
     .   - (6d0 - 8d0*l2m2pl)*V(9) + (1d0 - 2d0*l2m2pl)*V(10))
      dhy(6,9) = n*r3*((1d0 - 1.5d0*l2 + m2*(6d0*l2m2pl - 6.5d0))*V(8)
     .   + (4d0*n2m2 + (1d0 - 2d0*m2)*(2d0*l2m2pl - 1d0))*V(9)
     .   - (n2m2 + 0.5d0*l2m2pl*(1d0 - 2d0*m2))*V(10))
      dhz(6,9) = m*r3*((3d0*n2*l2m2pl - (1d0 - 1.5d0*l2m2pl)
     .  *(1d0 - 2d0*l2m2pl))*V(8) - (4d0*n2*l2m2pl
     .   + (2d0*l2m2pl - 1d0)*(1d0 - 2d0*l2m2pl))*V(9)
     .   + l2m2pl*0.5d0*(3d0 - 4d0*l2m2pl)*V(10))
      dhx(7,7) = l*(n2*(6d0 - 12d0*l2)*V(8)
     .   + (2d0 - 2d0*l2 - 10d0*n2 + 16d0*n2l2)*V(9)
     .   - (2d0 - 2d0*l2 - 4d0*n2 + 4d0*n2l2)*V(10))
      dhy(7,7) = m*( - 12d0*n2l2*V(8) - (2d0*l2n2pl - 16d0*n2l2)*V(9)
     .   - (4d0*n2l2 - 2d0*l2n2pl)*V(10))
      dhz(7,7) = n*(l2*(6d0 - 12d0*n2)*V(8)
     .   + (2d0 - 2d0*n2 - 10d0*l2 + 16d0*n2l2)*V(9)
     .   - (2d0 - 2d0*n2 - 4d0*l2 + 4d0*n2l2)*V(10))
      dhx(7,8) =  - n*((1.5d0*m2 - l2*(4.5d0 - 6d0*l2m2mi))*V(8)
     .   - (1d0 + 2d0*m2 - 8d0*l2*(1d0 - l2m2mi))*V(9)
     .   + (1d0 + 0.5d0*m2 - l2*(3.5d0 - 2d0*l2m2mi))*V(10))
      dhy(7,8) =  - lmn*((3d0 + 6d0*l2m2mi)*V(8)
     .   - (2d0 + 8d0*l2m2mi)*V(9) - (1d0 - 2d0*l2m2mi)*V(10))
      dhz(7,8) =  - l*(l2m2mi*(4.5d0 - 6d0*l2m2pl)*V(8)
     .   + ( - 4d0*n2*l2m2mi + (1d0 - 2d0*l2m2mi)*(1d0 - 2d0*l2m2pl))
     .  *V(9) - ( - n2*l2m2mi + (1d0 - 0.5d0*l2m2mi)*(1d0 - 2d0*l2m2pl))
     .  *V(10))
      dhx(7,9) = n*r3*((1d0 - 1.5d0*m2 + l2*(6d0*l2m2pl - 6.5d0))*V(8)
     .   + (4d0*n2l2 + (1d0 - 2d0*l2)*(2d0*l2m2pl - 1d0))*V(9)
     .   - (n2l2 + 0.5d0*l2m2pl*(1d0 - 2d0*l2))*V(10))
      dhy(7,9) = dhx(6,9)
      dhz(7,9) = l*r3*((3d0*n2*l2m2pl - (1d0 - 1.5d0*l2m2pl)
     .  *(1d0 - 2d0*l2m2pl))*V(8) - (4d0*n2*l2m2pl + (2d0
     .  *l2m2pl - 1d0)*(1d0 - 2d0*l2m2pl))*V(9)
     .   + l2m2pl*0.5d0*(3d0 - 4d0*l2m2pl)*V(10))
      dhx(8,8) = 3d0*l2m2mi*l*(1d0 - l2m2mi)*V(8)
     .   + 2d0*l*(n2 - 2d0*(1d0 - l2m2mi)*l2m2mi)*V(9)
     .   - 2d0*l*(n2 - 0.5d0*(1d0 - l2m2mi)*l2m2mi)*V(10)
      dhy(8,8) = ( - 3d0*l2m2mi*m*(1d0 + l2m2mi))*V(8)
     .   + 2d0*m*(n2 + 2d0*(1d0 + l2m2mi)*l2m2mi)*V(9)
     .   - 2d0*m*(n2 + 0.5d0*(1d0 + l2m2mi)*l2m2mi)*V(10)
      dhz(8,8) = n*( - 3d0*l2m2mi*l2m2mi*V(8)
     .   - 2d0*(l2m2pl - 2d0*l2m2mi*l2m2mi)*V(9)
     .   + (2d0*l2m2pl - l2m2mi*l2m2mi)*V(10))
      dhx(8,9) = r3*l*(((1d0 - 3d0*l2)*(1d0 - l2) + m2*(1d0 - 3d0
     .  *m2))*V(8) - 2d0*((1d0 - 2d0*l2)*(1d0 - l2) + m2*(1 - 2d0*m2))
     .  *V(9) + ((1d0 - l2)*(1d0 - l2) + m2*(1d0 - m2))*V(10))
      dhy(8,9) =  - r3*m*(((1d0 - 3d0*m2)*(1d0 - m2) + l2*(1d0 - 3d0
     .  *l2))*V(8) - 2d0*((1d0 - 2d0*m2)*(1d0 - m2) + l2*(1 - 2d0
     .  *l2))*V(9) + ((1d0 - m2)*(1d0 - m2) + l2*(1d0 - l2))*V(10))
      dhz(8,9) =  - r3*n*((l2*(1d0 - 3d0*l2) - m2*(1d0 - 3d0*m2))*V(8)
     .   - 2d0*(l2*(1d0 - 2d0*l2) - m2*(1d0 - 2d0*m2))*V(9)
     .   + (l2*(1d0 - l2) - m2*(1d0 - m2))*V(10))
      dhx(9,9) = l*n2*((9d0*l2m2pl - 6d0)*V(8)
     .   + (6d0 - 12d0*l2m2pl)*V(9) + 3d0*l2m2pl*V(10))
      dhy(9,9) = m*n2*((9d0*l2m2pl - 6d0)*V(8)
     .   + (6d0 - 12d0*l2m2pl)*V(9) + 3d0*l2m2pl*V(10))
      dhz(9,9) = n*l2m2pl*((6d0 - 9d0*l2m2pl)*V(8)
     .   - (6d0 - 12d0*l2m2pl)*V(9) - 3d0*l2m2pl*V(10))
      do  1  j = 5, 9
      do  1  i = 5, j
        dhx(i,j) = -dhx(i,j)
        dhy(i,j) = -dhy(i,j)
        dhz(i,j) = -dhz(i,j)
    1 continue

10000 continue
C  ---- Radial derivatives ---- 
      call skham(l,m,n,dV,nl,9,wk)

C  ---- Angular + Radial and make lower triangle ---- 
      do  2  j = 1, nl2
      do  2  i = 1, j
        dhx(i,j) = dhx(i,j)/r  - wk(i,j)*l
        dhy(i,j) = dhy(i,j)/r  - wk(i,j)*m
        dhz(i,j) = dhz(i,j)/r  - wk(i,j)*n
        dhr(i,j) = wk(i,j)*r
    2 continue
      if (.not. cryf) then
        do  3  j = 1, nl2-1
        do  3  i = j+1, nl2
          dhx(i,j) = dhx(j,i)*il(i)*il(j)
          dhy(i,j) = dhy(j,i)*il(i)*il(j)
          dhz(i,j) = dhz(j,i)*il(i)*il(j)
          dhr(i,j) = dhr(j,i)*il(i)*il(j)
    3   continue
      else
        do  4  j = 1, nl2-1
        do  4  i = j+1, nl2
          dhx(i,j) = dhx(j,i)
          dhy(i,j) = dhy(j,i)
          dhz(i,j) = dhz(j,i)
          dhr(i,j) = dhr(j,i)
    4   continue
      endif
      end
