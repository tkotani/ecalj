      subroutine mklagw(n,iopt,alfa,z0,z,w,ipr)
C- Quadrature weights for Laguerre Gaussian quadrature
C ----------------------------------------------------------------------
Ci Inputs
Ci   n,z0: numerical integration in interval (z0,infinity) with n points
Ci   alfa: Estimate of exponent to fn's asymtotic behavior (see Remarks)
Ci   iopt: 0, integration with Jacobian = exp(-alfa z) (see Remarks)
Ci         1, integration with Jacobian = 1
Co Outputs
Co   z,w
Cr Remarks
Cr   If iopt=0, integrates any linear combination of functions z**i
Cr     exactly for 0 <= i < 2n, where the Jacobian is exp(-alfa z) dz
Cr   If iopt=1, integrates any l.c. of functions z**i exp(-alfa z)
Cr     exactly for 0 <= i < 2n, where the Jacobian is dz
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer iopt,n,ipr
      double precision w(1),z(1),alfa,z0
C Local parameters
      integer in,iroot
      double precision delta,lagn,lnolnp,smachp,root,d1mach
      logical last

      smachp = dsqrt(d1mach(3))/100

C --- Find all the roots of L_n (L_1=1-x) ---
      z(1) = 1
      do  20  in = 2, n
C ... First guess for z
        z(in) = (z(in-1) + 3*z(in-1)/in)/(1.d0 - 1.d0/in)
        do  18  iroot = 1, in
          root = z(iroot) - z(iroot)/in
          if (ipr .ge. 110)
     .      print '('' mklagw: starting guess is'',i3,f9.5)', iroot,root
          delta = 2*smachp
   10     continue
          last = dabs(delta) .lt. smachp
          delta = -lnolnp(in,root)
          root = delta+root
          if (ipr .gt. 120) write(*,*) 'delta is ',delta,root,last
          if (.not. last) goto 10
          z(iroot) = root
          if (ipr .ge. 120) print 333, in,iroot,root,delta
  333     format('  n=',i2,': found root',i3,f10.5,': delta=',1pe10.1)
   18   continue
   20 continue

C --- Debugging: check for identical roots ---
      do  30  in = 1, n-1
   30 if (dabs(z(in)-z(in+1)).lt.smachp) stop 'mklagw: identical roots'

C --- Weights for alfa=1, z0=0 ---
      do  40  in = 1, n
   40 w(in) = z(in)/(lagn(n+1,z(in))*(n+1))**2

C --- Shift points and weights for alfa and z0 ---
      do  50  in = 1, n
        z(in) = z(in) / alfa + z0
        w(in) = w(in) * dexp(-alfa*z0) / alfa
   50 continue

C --- Scale weights for opt=1 ---
      if (iopt .eq. 1) then
        do  60  in = 1, n
   60   w(in) = w(in)*dexp(alfa*z(in))
      endif

C --- Printout ---
      if (ipr .lt. 100) return
      print *, 'mklagw: points and weights'
      do  70  in = 1, n
        print 71, z(in),w(in)
   71   format(1x,3(1pe18.10))
   70 continue

      end

      double precision function lagn(n,x)
C- Calculates Laguerre polynomical using a recursion relation
C ----------------------------------------------------------------------
Ci Inputs
Ci   n,x
Co Outputs
Co   lagn: L_n(x)
Cr Remarks
Cr   Recursion relation is
Cr   L_n = [(2*n-1-x)*L_(n-1) - (n-1)*L_(n-2)]/n
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer n
      double precision x
C Local parameters
      integer j
      double precision cj,jljm1,ljp1

C jljm1 is j*L_(j-1); cj is 2*j - 1 - x;  ljp1 is L_j + 1
      jljm1 = 0
      lagn = 1
      cj = 1-x
      do  100  j = 1, n
        ljp1 = (cj*lagn - jljm1)/j
        jljm1 = j*lagn
        cj = cj + 2
        lagn = ljp1
  100 continue
      end
      double precision function lnolnp(n,x)
C- Calculates derivative of Laguerre polynomical with recursion relation
C ----------------------------------------------------------------------
Ci Inputs
Ci   n,x
Co Outputs
Co   lagn: L_n(x)
Cr Remarks
Cr   Recursion relation for L is found in lagn; also
Cr   L'_n(x) = n/x [ L_n - L_(n-1) ]
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer n
      double precision x
C Local parameters
      integer j
      double precision cj,jljm1,ljp1,ln

C --- jljm1 is j*L_(j-1);  cj is 2*j - 1 - x; ljp1 is L_j + 1 ---
      jljm1 = 0
      ln = 1
      cj = 1-x
      do  100  j = 1, n
        ljp1 = (cj*ln - jljm1)/j
        jljm1 = j*ln
        cj = cj + 2
        ln = ljp1
  100 continue
      lnolnp = ln*x/(n*ln - jljm1)
      end

c$$$#if TEST
c$$$C shows z and weights integrate x**i from i=0..2n-1 exactly
c$$$      implicit none
c$$$      double precision z(200),w(200),alfa,z0
c$$$      double precision zz(200),ww(200),r1,r2,r3,resg
c$$$      integer n,iopt,ipr,i,j,nn
c$$$      double precision facti,dsum,res
c$$$
c$$$      iopt = 0
c$$$      ipr = 80
c$$$      print *, 'n,alfa,z0,iopt='
c$$$      read(*,*) n,alfa,z0,iopt
c$$$      call mklagw(n,iopt,alfa,z0,z,w,ipr)
c$$$      nn = 200
c$$$      call mklegw(nn,zz,ww,0)
c$$$
c$$$      facti = 1
c$$$      print *, 'sum z,w is', dsum(n,z,1), dsum(n,w,1)
c$$$      print *,
c$$$     .  '  i          int           int/int(0,infty)          err'
c$$$      do  10  i = 1, 2*n+1
c$$$        facti = facti * i
c$$$        res = 0
c$$$        if (iopt .eq. 0) then
c$$$          do  20  j = 1, n
c$$$   20     res = res + w(j)*z(j)**i
c$$$        else
c$$$          do  22  j = 1, n
c$$$   22     res = res + w(j)*z(j)**i*exp(-alfa*z(j))
c$$$        endif
c$$$
c$$$
c$$$C ...  Do integral from 0 to z0 using mlegw
c$$$        r1 = (z0-0)/2
c$$$        r2 = (z0+0)/2
c$$$        resg = 0
c$$$        do  40  j = 1, nn
c$$$
c$$$        r3 = r1*zz(nn-j+1) + r2
c$$$   40   resg = resg + ww(j-1+1)*r1*r3**i*exp(-alfa*r3)
c$$$        print 333, i, res, res/facti*alfa**(i+1),
c$$$     .    (res+resg)/facti*alfa**(i+1)-1
c$$$  333   format(i4,1x,3f20.15)
c$$$        if (i .eq. 2*n-1) print *, ' ------ exact ends here ---- '
c$$$   10 continue
c$$$      end
c$$$#endif

