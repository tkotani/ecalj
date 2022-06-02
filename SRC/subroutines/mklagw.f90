subroutine mklagw(n,iopt,alfa,z0,z,w,ipr)
  !- Quadrature weights for Laguerre Gaussian quadrature
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   n,z0: numerical integration in interval (z0,infinity) with n points
  !i   alfa: Estimate of exponent to fn's asymtotic behavior (see Remarks)
  !i   iopt: 0, integration with Jacobian = exp(-alfa z) (see Remarks)
  !i         1, integration with Jacobian = 1
  !o Outputs
  !o   z,w
  !r Remarks
  !r   If iopt=0, integrates any linear combination of functions z**i
  !r     exactly for 0 <= i < 2n, where the Jacobian is exp(-alfa z) dz
  !r   If iopt=1, integrates any l.c. of functions z**i exp(-alfa z)
  !r     exactly for 0 <= i < 2n, where the Jacobian is dz
  ! ----------------------------------------------------------------------
  !     implicit none
  ! Passed parameters
  integer :: iopt,n,ipr
  double precision :: w(1),z(*),alfa,z0
  ! Local parameters
  integer :: in,iroot
  double precision :: delta,lagn,lnolnp,smachp,root,d1mach
  logical :: last

  smachp = dsqrt(d1mach(3))/100

  ! --- Find all the roots of L_n (L_1=1-x) ---
  z(1) = 1
  do  20  in = 2, n
     ! ... First guess for z
     z(in) = (z(in-1) + 3*z(in-1)/in)/(1.d0 - 1.d0/in)
     do  18  iroot = 1, in
        root = z(iroot) - z(iroot)/in
        if (ipr >= 110) &
             print '('' mklagw: starting guess is'',i3,f9.5)', iroot,root
        delta = 2*smachp
10      continue
        last = dabs(delta) .lt. smachp
        delta = -lnolnp(in,root)
        root = delta+root
        if (ipr > 120) write(*,*) 'delta is ',delta,root,last
        if ( .NOT. last) goto 10
        z(iroot) = root
        if (ipr >= 120) print 333, in,iroot,root,delta
333     format('  n=',i2,': found root',i3,f10.5,': delta=',1pe10.1)
18   enddo
20 enddo

  ! --- Debugging: check for identical roots ---
  do    in = 1, n-1
     if (dabs(z(in)-z(in+1)) < smachp) stop 'mklagw: identical roots'
  enddo

  ! --- Weights for alfa=1, z0=0 ---
  do    in = 1, n
     w(in) = z(in)/(lagn(n+1,z(in))*(n+1))**2
  enddo

  ! --- Shift points and weights for alfa and z0 ---
  do  50  in = 1, n
     z(in) = z(in) / alfa + z0
     w(in) = w(in) * dexp(-alfa*z0) / alfa
50 enddo

  ! --- Scale weights for opt=1 ---
  if (iopt == 1) then
     do    in = 1, n
        w(in) = w(in)*dexp(alfa*z(in))
     enddo
  endif

  ! --- Printout ---
  if (ipr < 100) return
  print *, 'mklagw: points and weights'
  do  70  in = 1, n
     print 71, z(in),w(in)
71   format(1x,3(1pe18.10))
70 enddo

end subroutine mklagw

double precision function lagn(n,x)
  !- Calculates Laguerre polynomical using a recursion relation
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   n,x
  !o Outputs
  !o   lagn: L_n(x)
  !r Remarks
  !r   Recursion relation is
  !r   L_n = [(2*n-1-x)*L_(n-1) - (n-1)*L_(n-2)]/n
  ! ----------------------------------------------------------------------
  !     implicit none
  ! Passed parameters
  integer :: n
  double precision :: x
  ! Local parameters
  integer :: j
  double precision :: cj,jljm1,ljp1

  ! jljm1 is j*L_(j-1); cj is 2*j - 1 - x;  ljp1 is L_j + 1
  jljm1 = 0
  lagn = 1
  cj = 1-x
  do  100  j = 1, n
     ljp1 = (cj*lagn - jljm1)/j
     jljm1 = j*lagn
     cj = cj + 2
     lagn = ljp1
100 enddo
END function lagn
double precision function lnolnp(n,x)
  !- Calculates derivative of Laguerre polynomical with recursion relation
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   n,x
  !o Outputs
  !o   lagn: L_n(x)
  !r Remarks
  !r   Recursion relation for L is found in lagn; also
  !r   L'_n(x) = n/x [ L_n - L_(n-1) ]
  ! ----------------------------------------------------------------------
  !     implicit none
  ! Passed parameters
  integer :: n
  double precision :: x
  ! Local parameters
  integer :: j
  double precision :: cj,jljm1,ljp1,ln

  ! --- jljm1 is j*L_(j-1);  cj is 2*j - 1 - x; ljp1 is L_j + 1 ---
  jljm1 = 0
  ln = 1
  cj = 1-x
  do  100  j = 1, n
     ljp1 = (cj*ln - jljm1)/j
     jljm1 = j*ln
     cj = cj + 2
     ln = ljp1
100 enddo
  lnolnp = ln*x/(n*ln - jljm1)
END function lnolnp

!$$$#if TEST
!$$$C shows z and weights integrate x**i from i=0..2n-1 exactly
!$$$      implicit none
!$$$      double precision z(200),w(200),alfa,z0
!$$$      double precision zz(200),ww(200),r1,r2,r3,resg
!$$$      integer n,iopt,ipr,i,j,nn
!$$$      double precision facti,dsum,res
!$$$
!$$$      iopt = 0
!$$$      ipr = 80
!$$$      print *, 'n,alfa,z0,iopt='
!$$$      read(*,*) n,alfa,z0,iopt
!$$$      call mklagw(n,iopt,alfa,z0,z,w,ipr)
!$$$      nn = 200
!$$$      call mklegw(nn,zz,ww,0)
!$$$
!$$$      facti = 1
!$$$      print *, 'sum z,w is', dsum(n,z,1), dsum(n,w,1)
!$$$      print *,
!$$$     .  '  i          int           int/int(0,infty)          err'
!$$$      do  10  i = 1, 2*n+1
!$$$        facti = facti * i
!$$$        res = 0
!$$$        if (iopt .eq. 0) then
!$$$          do  20  j = 1, n
!$$$   20     res = res + w(j)*z(j)**i
!$$$        else
!$$$          do  22  j = 1, n
!$$$   22     res = res + w(j)*z(j)**i*exp(-alfa*z(j))
!$$$        endif
!$$$
!$$$
!$$$C ...  Do integral from 0 to z0 using mlegw
!$$$        r1 = (z0-0)/2
!$$$        r2 = (z0+0)/2
!$$$        resg = 0
!$$$        do  40  j = 1, nn
!$$$
!$$$        r3 = r1*zz(nn-j+1) + r2
!$$$   40   resg = resg + ww(j-1+1)*r1*r3**i*exp(-alfa*r3)
!$$$        print 333, i, res, res/facti*alfa**(i+1),
!$$$     .    (res+resg)/facti*alfa**(i+1)-1
!$$$  333   format(i4,1x,3f20.15)
!$$$        if (i .eq. 2*n-1) print *, ' ------ exact ends here ---- '
!$$$   10 continue
!$$$      end
!$$$#endif

