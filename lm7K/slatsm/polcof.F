      subroutine polcof(lopt,x0,x,y,n,cof)
C- Coefficients to interpolating polynomial given a tabulated function
C  Adapted from Numerical Recipes
C  lopt nonzero => evaluate intermediate results in quad precision
C  If compiler does not allow it, quad and double are the same.
C     implicit none
      integer n,lopt
      double precision x0,x(n),y(n),cof(n)
      integer i,j,k,nmax
      parameter (nmax=50)
      double precision s(nmax),phi,ff,b,dxi,dxj
#if NOQUAD
      double precision zs(nmax),zphi,zff,zb,zdxi,zdxj
#else
      real*16 zs(nmax),zphi,zff,zb,zdxi,zdxj
#endif

C ... double precision
      if (lopt .eq. 0) then
        do  11  i = 1, n
          s(i) = 0
          cof(i) = 0
   11   continue
        s(n) = -(x(1)-x0)
        do  13  i = 2, n
          dxi = x(i)-x0
          do  12  j = n+1-i, n-1
   12     s(j) = s(j) - dxi*s(j+1)
          s(n) = s(n) - dxi
   13   continue
        do  16  j = 1, n
          dxj = x(j)-x0
          phi = n
          do  14  k = n-1, 1, -1
   14     phi = k*s(k+1) + dxj*phi
          ff = y(j)/phi
          b = 1
          do  15  k = n, 1,-1
            cof(k) = cof(k) + b*ff
            b = s(k) + dxj*b
   15     continue
   16   continue

C ... quad precision
      else
        do  21  i = 1, n
          zs(i) = 0
          cof(i) = 0
   21   continue
        zs(n) = -(x(1)-x0)
        do  23  i = 2, n
          zdxi = x(i)-x0
          do  22  j = n+1-i, n-1
   22     zs(j) = zs(j) - zdxi*zs(j+1)
          zs(n) = zs(n) - zdxi
   23   continue
        do  26  j = 1, n
          zdxj = x(j)-x0
          zphi = n
          do  24  k = n-1, 1, -1
   24     zphi = k*zs(k+1) + zdxj*zphi
          zff = y(j)/zphi
          zb = 1
          do  25  k = n, 1,-1
            cof(k) = cof(k) + zb*zff
            zb = zs(k) + zdxj*zb
   25     continue
   26   continue
      endif

      end

