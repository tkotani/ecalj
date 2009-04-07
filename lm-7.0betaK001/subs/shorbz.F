      subroutine shorbz(p,pout,plat,qlat)
C- Shortens vector to equivalent in first Brillouin zone.
C ----------------------------------------------------------------
Ci Inputs:  
Ci   plat,qlat lattice vectors and inverse
Ci   p         vector to shorten
Co Outputs:
Co   pout      shortened p.   pout and p can point to the same address.
Cr Remarks
Cr   Switch around plat,qlat to shorten reciprocal space vectors.
Cr   Jan 1997 Adapted from shorps to fix bug:  Example:
Cr   plat=  -0.5  0.5  1.7517  0.5  -0.5  1.7517  0.5  0.5  -1.7517
Cr   p= 0.0 -0.5 -1.26384
Cr   Should get pout -> 0.5 0.0 0.48786, not -0.5 1.0 0.48786.
C ----------------------------------------------------------------
C     implicit none
      double precision p(3),pout(3),plat(3,3),qlat(3,3),x(3),x0,xx,a2,ap
      double precision tol
      parameter (tol=-1d-10)
      integer i,j,m,j2min,j3min,j1,j2,j3

C --- Reduce to unit cell centered at origin ---
      do  1  i = 1, 3
C ... x is projection of pin along plat(i), with multiples of p removed
      x0 = p(1)*qlat(1,i)+p(2)*qlat(2,i)+p(3)*qlat(3,i)
      xx = idnint(x0)
    1 x(i) = x0-xx
C ... pout is x rotated back to Cartesian coordinates
      do  2  m = 1, 3
    2 pout(m) = x(1)*plat(m,1)+x(2)*plat(m,2)+x(3)*plat(m,3)

C --- Try shortening by adding +/- basis vectors ---
   15   continue
        do  10  j1 =  0, 1
        j2min = -1
        if (j1 .eq. 0) j2min = 0
        do  10  j2 = j2min, 1
        j3min = -1
        if (j1 .eq. 0 .and. j2 .eq. 0) j3min = 0
        do  10  j3 = j3min, 1

C     ... (-1,0,1) (plat(1) + (-1,0,1) plat(2)) + (-1,0,1) plat(3))
          do  17  i = 1, 3
   17     x(i) = plat(i,1)*j1 + plat(i,2)*j2 + plat(i,3)*j3
          a2 = x(1)*x(1) + x(2)*x(2) + x(3)*x(3)
          ap = pout(1)*x(1) + pout(2)*x(2) + pout(3)*x(3)
          j = 0
          if (a2 + 2*ap .lt. tol) j = 1
          if (a2 - 2*ap .lt. tol) j = -1
          if (j .ne. 0) then
            pout(1) = pout(1) + j*x(1)
            pout(2) = pout(2) + j*x(2)
            pout(3) = pout(3) + j*x(3)
            goto 15
          endif
   10   continue

      end
C      subroutine fmain
C      implicit none
C      double precision plat(9),qlat(9),p(3),p1(3),xx
C      integer mode(3)
C
C      data plat /-0.5d0,0.5d0,1.7517d0,
C     .            0.5d0,-.5d0,1.7517d0,
C     .            0.5d0,0.5d0,-1.7517d0/
C
C      integer w(10000)
C      common /w/ w
C
C      data p /0.0d0,-0.5d0,-1.2638400000000001d0/
C
C      call wkinit(10000)
C
CC ... qlat = (plat^-1)^T so that qlat^T . plat = 1
C      call mkqlat(plat,qlat,xx)
C
C      call shorbz(p,p1,plat,qlat)
C      call prmx('p1 from shorbz',p1,1,1,3)
C
C      mode(1) = 2
C      mode(2) = 2
C      mode(3) = 3
C      call shorps(1,plat,mode,p,p1)
C      call prmx('p1 from shorps',p1,1,1,3)
C      end

