      double precision function drr2(plat,tau1,tau2,i,j,k,dr)

C- Calculates the vector connecting two sites in a solid
C ----------------------------------------------------------------
Ci Inputs
Ci   plat: primitive lattice vectors
Ci   tau1,tau2: basis vectors of the two sites
Ci   i,j,k:the number of primitive lattice vectors separating sites
Co Outputs
Co   dr:   connecting vector tau2 - tau1
Co   drr2: square of the length of this vector
Cr Remarks
Cr   Using the TB package and a table of indices iax, the connecting
Cr   vector and the square of the distance is obtained by
Cr      rsqr = drr2(plat,bas(1,iax(1)),bas(1,iax(2)),
Cr     .            iax(3),iax(4),iax(5),dr)
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer i,j,k
      double precision dr(3)
      double precision plat(3,3),tau1(3),tau2(3)
C Local parameters
      integer ix

      drr2 = 0.d0
      do  10  ix = 1, 3
        dr(ix) = tau2(ix) - tau1(ix) +
     .           plat(ix,1)*i + plat(ix,2)*j + plat(ix,3)*k
        drr2 = drr2 + dr(ix)**2
   10 continue
      end function drr2


