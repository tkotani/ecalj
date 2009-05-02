      subroutine mklegw(n,z,w,ipr)
c- quadrature weights for legendre gaussian quadrature
c ----------------------------------------------------------------------
ci inputs
ci   n:   number of mesh points for numerical integration
ci        (on interval -1, 1)
ci   ipr: verbosity
co outputs
co   z,w
cr remarks
cr   integrates 2(n-1) order polynomial exactly on interval (-1,1)
c ----------------------------------------------------------------------
C     implicit none
c passed parameters
      integer n
      double precision z(1),w(1)
      double precision plegn,pnp
c local parameters
      integer iroot,in,ipr
      double precision root,delta,machep,d1mach,pi

      pi = 4 * datan(1d0)
c|    machep = 10*d1mach(3)
      machep=1d-14

c --- find all the roots of p_n ---
      do  200  iroot = 1, n
        z(iroot) = dcos(pi*(2*iroot-.5d0)/(2*n+1))
        root = z(iroot)
        if (ipr .ge. 110)
     .    print *, 'mklegw: ',iroot,' starting guess is ',root
  100   continue
        delta = -plegn(n,root)/pnp(n,root)
        if (ipr .ge. 110) write(*,*) 'delta is ',delta
        root = root + delta
        if (dabs(delta) .gt. dabs(machep*root)) goto 100
        z(iroot) = root
  200 continue

c --- debugging:  check for identical roots ---
      do  300  in = 1, n-1
        if (dabs(z(in)-z(in+1)) .lt. machep)
     .    stop 'mklegw: identical roots'
  300 continue

c --- make the weights ---
      do  700  in = 1, n
        w(in) = 2/((1-z(in)**2)*(pnp(n,z(in))**2))
  700 continue

c --- printout ---
      if (ipr .lt. 50) return
      print *
      print *, 'mklegw: roots and weighting factors'
      do  400  in = 1, n
        write(*,410) z(in),w(in)
  410   format(1x,2(1pe26.15))
  400 continue
      end

      double precision function plegn(n,x)
c- calculates legendre polynomical using a recursion relation
c ----------------------------------------------------------------------
ci inputs
ci   n,x
co outputs
co   plegn: p_n(x)
cr remarks
cr   recursion relation is
cr   p_n = ((2*n-1)*x*p_(n-1) - (n-1)*p_(n-2))/n
c ----------------------------------------------------------------------
C     implicit none
c passed parameters
      integer n
      double precision x
c local parameters
      double precision jpjm1,cj,pjp1
      integer j

c jpjm1 is j*p_(j-1);  cj is 2*j - 1;  pjp1 is p_(j+1)
      jpjm1 = 0
      plegn = 1
      cj = 1
      do  100  j = 1, n
        pjp1 = (cj*x*plegn - jpjm1)/j
        jpjm1 = j*plegn
        cj = cj + 2
        plegn = pjp1
  100 continue
      end

      double precision function pnp(n,x)
c- calculates derivative of legendre polynomical from recursion relation
c ----------------------------------------------------------------------
ci inputs
ci   n,x
co outputs
co   pnp: p'_n(x)
cr remarks
cr   recursion relations for p and p' are
cr   p_n (x) = ((2*n-1)*x*p_(n-1) - (n-1)*p_(n-2))/n
cr   p'_n(x) = n/(1-x^2) (-x*p_n + p_(n-1) )
c ----------------------------------------------------------------------
C     implicit none
c passed parameters
      integer n
      double precision x
c local parameters
      double precision jpjm1,cj,pjp1,ln
      integer j

c jpjm1 is j*p_(j-1);  cj is 2*j - 1;  pjp1 is p_(j+1)
      jpjm1 = 0
      ln = 1
      cj = 1
      do  100  j = 1, n
        pjp1 = (cj*x*ln - jpjm1)/j
        jpjm1 = j*ln
        cj = cj + 2
        ln = pjp1
  100 continue
      pnp = 1/(1-x**2)*(-n*x*ln + jpjm1)
      end
