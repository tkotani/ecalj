subroutine mklegw(n,z,w,ipr)
  !- Quadrature weights for Legendre Gaussian quadrature
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   n:   Number of mesh points for numerical integration
  !i        (on interval -1, 1)
  !i   ipr: verbosity
  !o Outputs
  !o   z,w
  !r Remarks
  !r   Integrates 2(n-1) order polynomial exactly on interval (-1,1)
  ! ----------------------------------------------------------------------
  implicit none
  integer :: n
  double precision :: z(1),w(1)
  double precision :: plegn,pnp
  integer :: iroot,in,ipr
  double precision :: root,delta,machep,d1mach,pi
  pi = 4 * datan(1d0)
  machep=1d-14
  ! --- Find all the roots of P_n ---
  do  200  iroot = 1, n
     z(iroot) = dcos(pi*(2*iroot-.5d0)/(2*n+1))
     root = z(iroot)
     if (ipr >= 110) &
          print *, 'mklegw: ',iroot,' starting guess is ',root
100  continue
     delta = -plegn(n,root)/pnp(n,root)
     if (ipr >= 110) write(*,*) 'delta is ',delta
     root = root + delta
     if (dabs(delta) > dabs(machep*root)) goto 100
     z(iroot) = root
200 enddo

  ! --- debugging:  check for identical roots ---
  do  300  in = 1, n-1
     if (dabs(z(in)-z(in+1)) < machep) &
                                ! top2rx 2013.08.09 kino     .    stop 'mklegw: identical roots'
          call rx( 'mklegw: identical roots')
300 enddo
  ! --- Make the weights ---
  do  700  in = 1, n
     w(in) = 2/((1-z(in)**2)*(pnp(n,z(in))**2))
700 enddo
  ! --- Printout ---
  if (ipr < 80) return
  print *
  print *, 'mklegw: roots and weighting factors'
  do  400  in = 1, n
     write(*,410) z(in),w(in)
410  format(1x,2(1pe26.15))
400 enddo
end subroutine mklegw
subroutine gausq(n,x1,x2,xp,wp,mode,ipr)
  !- Quadrature points and weights for Legendre Gaussian quadrature
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   n,x1,x2: Number of mesh points and interval (x1,x2)
  !i            NB: x2 has different meaning when 10s digit of mode nonzero
  !i   mode: 0, straight gaussian quadrature
  !i         1s digit: n, scale weights by x**n
  !i        10s digit: 1, Quadrature for int_x1^infty.  Internal
  !i                      transformation and integration over y=exp(-x/x2).
  !i   ipr: verbosity
  !o Outputs
  !o   xp,wp
  ! ----------------------------------------------------------------------
  !     implicit none
  integer :: n,ipr,i,nn,mode
  double precision :: x1,x2,xp(n),wp(n),xmin,xmax,y,d1mach
  character(20) :: fmt

  ! --- xmin and xmax from x1,x2, depending on mode ---
  if (mod(mode,100)/10 == 0) then
     xmin = x1
     xmax = x2
  else
     xmin = 0
     xmax = exp(-x1/x2)
  endif

  ! --- Straight Gaussian quadrature for interval (xmin,xmax) ---
  call mklegw(n,xp,wp,0)
  call dscal(n,-(xmax-xmin)/2,xp,1)
  call dscal(n,(xmax-xmin)/2,wp,1)
  do   i = 1, n
     xp(i) = xp(i) + (xmax+xmin)/2
  enddo

  ! --- 10s digit of mode nonzero => shift points y to x ---
  if (mod(mode,100)/10 == 1) then
     do  20  i = 1, n
        y = xp(i)
        xp(i) = -x2*dlog(y)
        wp(i) = x2*wp(i)/y
20   enddo
  endif

  ! --- 1s digit of mode = n => weight points by x**n
  nn = mod(mode,10)
  if (nn > 0) then
     do  i = 1, n
        wp(i) = wp(i)*xp(i)**nn
     enddo
  endif

  ! --- Printout ---
  if (ipr < 80) return
  print 333, mode
333 format(' gausq: points and weights, mode=',i3)
  !      i = -dlog(d1mach(3))/dlog(10d0)
  !      call awrit2('%x(1x,2(1pe%i.%i))',fmt,len(fmt),0,i+10,i)
  do  400  i = 1, n
     write(*,fmt) xp(i),wp(i)
410  format(1x,2(1pe26.15))
400 enddo

end subroutine gausq
double precision function plegn(n,x)
  !- Calculates Legendre polynomical using a recursion relation
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   n,x
  !o Outputs
  !o   plegn: P_n(x)
  !r Remarks
  !r   Recursion relation is
  !r   P_n = [(2*n-1)*x*P_(n-1) - (n-1)*P_(n-2)]/n
  ! ----------------------------------------------------------------------
  !     implicit none
  ! Passed parameters
  integer :: n
  double precision :: x
  ! Local parameters
  double precision :: jpjm1,cj,pjp1
  integer :: j

  ! jpjm1 is j*p_(j-1);  cj is 2*j - 1;  pjp1 is p_(j+1)
  jpjm1 = 0
  plegn = 1
  cj = 1
  do  100  j = 1, n
     pjp1 = (cj*x*plegn - jpjm1)/j
     jpjm1 = j*plegn
     cj = cj + 2
     plegn = pjp1
100 enddo
END function plegn

double precision function pnp(n,x)
  !- Calculates derivative of Legendre polynomical from recursion relation
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   n,x
  !o Outputs
  !o   pnp: P'_n(x)
  !r Remarks
  !r   Recursion relations for P and P' are
  !r   P_n (x) = [(2*n-1)*x*P_(n-1) - (n-1)*P_(n-2)]/n
  !r   P'_n(x) = n/(1-x^2) [-x*P_n + P_(n-1) ]
  ! ----------------------------------------------------------------------
  !     implicit none
  ! Passed parameters
  integer :: n
  double precision :: x
  ! Local parameters
  double precision :: jpjm1,cj,pjp1,ln
  integer :: j

  ! jpjm1 is j*p_(j-1);  cj is 2*j - 1;  pjp1 is p_(j+1)
  jpjm1 = 0
  ln = 1
  cj = 1
  do  100  j = 1, n
     pjp1 = (cj*x*ln - jpjm1)/j
     jpjm1 = j*ln
     cj = cj + 2
     ln = pjp1
100 enddo
  pnp = 1/(1-x**2)*(-n*x*ln + jpjm1)
END function pnp

