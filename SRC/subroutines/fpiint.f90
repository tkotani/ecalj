module m_fpiint !- Points and weights for integration on a sphere surface
  public fpiint
  private
  contains
subroutine fpiint(nx,np,nxp,x,w) !- Points and weights for integration on a sphere surface
  ! takao. Maybe it will be better to replace this with spherical design
  ! see http://www2.research.att.com/~njas/sphdesigns/
  !i Inputs
  !i   nx    :number of points in polar angle (Legendre integration)
  !i          Use nx<0 for special points; see Remarks.
  !i   np    :number of points in phi (uniform points)
  !i         :np=0 => np depends on nx making dphi approx constant.
  !i          for nx<0, np is not used.
  !o Outputs
  !o   nxp   :total number of number of points in quadrature
  !o   x     :cartesian coordinates of points on unit sphere
  !o   w     :weights corresponding to x
  !r Remarks
  !r   fpiint generates a mesh of points on a unit sphere for angular
  !r   integration, either using a set of special generated from the
  !r   Platonic solids, or by integrating the polar angle with Legendre
  !r   gaussian quadrature and the azimuthal angle with a set of evenly
  !r   spaced points on a circle.
  !r   For special points, invoke fpiint with one of the following:
  !     nx= -4 integrates any ylm<=2 exactly (tetrahedron)
  !     nx= -6 integrates any ylm<=3 exactly (faces of cube)
  !     nx= -8 integrates any ylm<=3 exactly (cube)
  !        -12 integrates any ylm<=5 exactly (icosahedron)
  !        -20 integrates any ylm<=5 exactly (faces of icosahedron)
  !        -30 integrates any ylm<=5 exactly (sides of icosahedron)
  !        -60 integrates any ylm<=5 exactly (buckeyball)
  !        -32 integrates any ylm<=9 exactly  (combination of 12,20)
  !        -62 integrates any ylm<=11 exactly (combination of 12,20,30)
  !        -92 integrates any ylm<=11 exactly (combination of 12,20,60)
  !       -122 integrates any ylm<=15 exactly (combination of 12,20,30,60)
  ! ----------------------------------------------------------------------
  !     implicit none
  integer :: nx,np,nxp
  double precision :: x(3,*),w(*),fpi
  integer :: i,iprint

  if (nx >= 0) then
     call nintsp(nx,np,nxp,x,w)
  else
     fpi = 16*datan(1d0)
     nxp = -nx
     if (nx == -32) then
        call platsl(x,12)
        call platsl(x(1,13),20)
        do   i = 1, 12
           w(i) = 5d0*fpi/(14*12)
        enddo
        do    i = 13, 32
           w(i) = 9d0*fpi/(14*20)
        enddo
     elseif (nx == -62) then
        call platsl(x,12)
        call platsl(x(1,13),20)
        call platsl(x(1,33),30)
        do  i = 1, 12
           w(i) = 125d0*fpi/(14*12*33)
        enddo
        do   i = 13, 32
           w(i) = 81d0*fpi/(14*20*33)
        enddo
        do   i = 33, 62
           w(i) = 256d0*fpi/(14*30*33)
        enddo
     elseif (nx == -92) then
        call platsl(x,12)
        call platsl(x(1,13),20)
        call platsl(x(1,33),60)
        do    i = 1, 12
           w(i) = 1/12.34817490904537d0*fpi/12
        enddo
        do   i = 13, 32
           w(i) = 2.986997567806883d0/12.34817490904537d0*fpi/20
        enddo
        do    i = 33, 92
           w(i) = 8.361177341238484d0/12.34817490904537d0*fpi/60
        enddo
     elseif (nx == -122) then
        call platsl(x,12)
        call platsl(x(1,13),20)
        call platsl(x(1,33),30)
        call platsl(x(1,63),60)
        do   i = 1, 12
           w(i) = (0.0939463041645901d0)*fpi/12
        enddo
        do    i = 13, 32
           w(i) = (0.2373458837681504d0)*fpi/20
        enddo
        do    i = 33, 92
           w(i) = (0.0378880378880377d0)*fpi/30
        enddo
        do    i = 63, 122
           w(i) = (0.6308197741792218d0)*fpi/60
        enddo
     else
        call platsl(x,nxp)
        do   i = 1, nxp
           w(i) = fpi/nxp
        enddo
     endif
  endif
  ! --- Printout ---
  if (iprint() < 80) return
  print '(/'' fpiint:'',i5,'' points generated:'')', nxp
  do  i = 1, nxp
     print 333, i, x(1,i), x(2,i), x(3,i), w(i)
  enddo
333 format(i3,4f20.15)
end subroutine fpiint

subroutine nintsp(nx,np,nxp,x,w)
  !- Points and weights for Legendre integration on a sphere surface
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   nx    :number of points in polar angle (Legendre integration)
  !i   np    :number of points in phi (uniform points)
  !i         :np=0 => np depends on nx making dphi approx constant.
  !o Outputs
  !o   nxp   :total number of number of points in quadrature
  !o   x     :cartesian coordinates of points on unit sphere
  !o   w     :weights corresponding to x
  !r Remarks
  !r   nintsp generates a mesh of points on a unit sphere for angular
  !r   integration, using Legendre gaussian quadrature for the polar
  !r   angle, and evenly spaced points on a circle for phi integration
  !r   See fpiint as an alternative routine that uses special points.
  ! ----------------------------------------------------------------------
  !     implicit none
  integer :: nx,np,nxp
  double precision :: x(3,1),w(1)
  integer :: ix,ip,i,npi
  double precision :: xtab(100),wx(100),sint,cost,phi,pi2

  pi2 = 8*datan(1d0)
  call mklegw(nx,xtab,wx,0)
  i = 0
  do    ip = 1, np
     do  ix = 1, nx
        i = i+1
        cost = xtab(ix)
        phi =  (pi2*ip)/np
        sint = dsqrt(1-cost**2)
        x(1,i)=dcos(phi)*sint
        x(2,i)=dsin(phi)*sint
        x(3,i)=cost
        w(i)  =wx(ix)*(pi2/np)
     enddo
  enddo
  nxp = nx*np
  if (np > 0) return

  ! --- Let np depend on x ---
  i = 0
  do  20  ix = 1, nx
     cost = xtab(ix)
     sint = dsqrt(1-cost**2)
     npi  = max(nint(sint*2*nx),1)
     do  30  ip = 1, npi
        i = i+1
        phi =  (pi2*ip)/npi
        x(1,i)=dcos(phi)*sint
        x(2,i)=dsin(phi)*sint
        x(3,i)=cost
        w(i)  =wx(ix)*(pi2/npi)
30   enddo
20 enddo
  nxp = i
end subroutine nintsp
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
  double precision :: plegn!,pnp
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

subroutine platsl(r,npt)
  !- Generates points corresponding platonic and related solids.
  !  npt (number of points) must be 4, 6, 8, 12, 20, 30 or 60
  !  30 generates points midway between the sides of the icosahedron.
  !  60 generates points for the buckeyball.
  !     implicit none
  integer :: npt
  double precision :: r(3,npt)
  double precision :: d,cp,cx,pi,sx,sp,ddot,x(3,12),r2(3),wb,wi,rm(3,3)
  integer :: i
  ! --- npt 4 or 8 ---
  if (npt == 4 .OR. npt == 8) then
     cx = 1/dsqrt(3d0)
     r(1:3,1:4) = cx
     do  i = 1, 3
        r(i,i) = -cx
        r(i,4) = -cx
     enddo
     if (npt == 8) then
        r(1:3,5:8) = -r(1:3,1:4)
     endif
     return
  endif
  ! --- npt 6 ---
  if (npt == 6) then
     r(1:3,1:6) = 0d0
     do  i = 1, 3
        r(i,i) = 1
        r(i,i+3) = -1
     enddo
     return
  endif

  if (npt /= 12 .AND. npt /= 20 .AND. npt /= 30 .AND. &
       npt /= 60) call rx('platsl: bad npt')
  pi = 4*datan(1d0)
  cp = dcos(2*pi/5)
  sp = dsqrt(1-cp**2)
  cx = 1/dsqrt(5d0)
  sx = dsqrt(1-cx**2)

  ! --- Icosahedron ---
  ! ... Top 6 points
  r(1,6) = 0d0
  r(2,6) = 0d0
  r(3,6) = 1d0
  r(1,1) = sx*r(3,6)
  r(2,1) = 0d0
  r(3,1) = cx*r(3,6)
  do  1210  i = 2, 5
     r(1,i) = cp*r(1,i-1) - sp*r(2,i-1)
     r(2,i) = sp*r(1,i-1) + cp*r(2,i-1)
     r(3,i) = r(3,1)
1210 enddo

  ! ... Bottom 6 points of icosahedron
  r(1,12) = 0d0
  r(2,12) = 0d0
  r(3,12) =-1d0
  r(1,7) = sx*r(3,12)
  r(2,7) = 0d0
  r(3,7) = cx*r(3,12)
  do  1220  i = 8,11
     r(1,i) = cp*r(1,i-1) - sp*r(2,i-1)
     r(2,i) = sp*r(1,i-1) + cp*r(2,i-1)
     r(3,i) = r(3,7)
1220 enddo

  ! --- Points on the twenty faces of the icosahedron ---
  if (npt == 20) then
     do i = 1, 3
        r2(i) = r(i,1)+r(i,2)+r(i,6)
     enddo
     d = dsqrt(ddot(3,r2,1,r2,1))
     call dpscop(r,x,36,1,1,1/d)
     do  2040  i = 1, 5
        r(1,i) = x(1,i) + x(1,1+mod(i,5)) + x(1,6)
        r(2,i) = x(2,i) + x(2,1+mod(i,5)) + x(2,6)
        r(3,i) = x(3,i) + x(3,1+mod(i,5)) + x(3,6)

        r(1,5+i) = x(1,6+i) + x(1,7+mod(i,5)) + x(1,12)
        r(2,5+i) = x(2,6+i) + x(2,7+mod(i,5)) + x(2,12)
        r(3,5+i) = x(3,6+i) + x(3,7+mod(i,5)) + x(3,12)

        r(1,10+i) = x(1,i) + x(1,1+mod(i,5)) + x(1,7+mod(i+2,5))
        r(2,10+i) = x(2,i) + x(2,1+mod(i,5)) + x(2,7+mod(i+2,5))
        r(3,10+i) = x(3,i) + x(3,1+mod(i,5)) + x(3,7+mod(i+2,5))

        r(1,15+i) = x(1,6+i) + x(1,7+mod(i,5)) + x(1,1+mod(i+2,5))
        r(2,15+i) = x(2,6+i) + x(2,7+mod(i,5)) + x(2,1+mod(i+2,5))
        r(3,15+i) = x(3,6+i) + x(3,7+mod(i,5)) + x(3,1+mod(i+2,5))

2040 enddo
  endif

  ! --- Points on the thirty edges of the icosahedron ---
  if (npt == 30) then
     do   i = 1, 3
        r2(i) = r(i,1)+r(i,2)
     enddo
     d = dsqrt(ddot(3,r2,1,r2,1))
     call dpscop(r,x,36,1,1,1/d)
     do  3040  i = 1, 5
        r(1,i) = x(1,i) + x(1,6)
        r(2,i) = x(2,i) + x(2,6)
        r(3,i) = x(3,i) + x(3,6)

        r(1,5+i) = x(1,6+i) + x(1,12)
        r(2,5+i) = x(2,6+i) + x(2,12)
        r(3,5+i) = x(3,6+i) + x(3,12)

        r(1,10+i) = x(1,i) + x(1,7+mod(i+2,5))
        r(2,10+i) = x(2,i) + x(2,7+mod(i+2,5))
        r(3,10+i) = x(3,i) + x(3,7+mod(i+2,5))

        r(1,15+i) = x(1,6+i) + x(1,1+mod(i+2,5))
        r(2,15+i) = x(2,6+i) + x(2,1+mod(i+2,5))
        r(3,15+i) = x(3,6+i) + x(3,1+mod(i+2,5))

        r(1,20+i) = x(1,i) + x(1,1+mod(i,5))
        r(2,20+i) = x(2,i) + x(2,1+mod(i,5))
        r(3,20+i) = x(3,i) + x(3,1+mod(i,5))

        r(1,25+i) = x(1,6+i) + x(1,7+mod(i,5))
        r(2,25+i) = x(2,6+i) + x(2,7+mod(i,5))
        r(3,25+i) = x(3,6+i) + x(3,7+mod(i,5))

3040 enddo
  endif

  ! --- Buckeyball ---
  ! The 30 edges of the icosahedron lie at the midpoints of the arcs
  ! connecting the 60 buckeyball points.  Determine arc length wb from
  ! sol'n of sin^2(wi-wb)cos(2pi/5) + cos^2(wi-wb) = cos(2*wb)
  ! where wi is half the icosahedral angle.  Orient making a
  ! pentagon at the top, with symmetry r5z.  This is replicated 5 times
  ! by rotation by euler angles z:i*2pi/5,y:2*wi,z:pi for i=0..4 to
  ! generate the 30 points for the top half of the buckeyball.
  if (npt == 60) then
     wi = datan(2/(1+sqrt(5d0)))
     wb = 0.2031689460497291d0
     ! ...   Five points of top pentagon
     r(1,1) = dsin(wi-wb)
     r(2,1) = 0
     r(3,1) = dcos(wi-wb)
     do  6010  i = 2, 5
        r(1,i) = cp*r(1,i-1) - sp*r(2,i-1)
        r(2,i) = sp*r(1,i-1) + cp*r(2,i-1)
        r(3,i) = r(3,1)
6010 enddo
     ! ...   Replicate these 5 points by rotations to other pentagons
     do  6020  i = 0, 4
        call eulerm(i*2*pi/5,2*wi,pi,rm)
        call dgemm('T','N',3,5,3,1d0,rm,3,r,3,0d0,r(1,6+5*i),3)
6020 enddo
     ! ...   Bottom 60 points are inversions of top 60
     do  6030  i = 1, 30
        r(1,61-i) = -r(1,i)
        r(2,61-i) = -r(2,i)
        r(3,61-i) = -r(3,i)
6030 enddo
  endif

  !      print '(/'' platsl:'',i5,'' points generated:'')', npt
  !      do  90  i = 1, npt
  !   90 print 333, i, r(1,i), r(2,i), r(3,i)
  !  333 format(i3,3f20.15)

  !  999 continue
  !      do  91  i = 1, npt
  !      do  91  j = 1, npt
  !   91 print 334, i,j,ddot(3,r(1,i),1,r(1,j),1)
  !  334 format(2i3,f20.10)

end subroutine platsl
subroutine eulerm(alpha,beta,gamma,r)
  !- Generate the rotation matrix corresponding to Euler angles
  !     implicit none
  double precision :: r(3,3),alpha,beta,gamma
  double precision :: ca,cb,cg,sa,sb,sg

  ca = dcos(alpha)
  sa = dsin(alpha)
  cb = dcos(beta)
  sb = dsin(beta)
  cg = dcos(gamma)
  sg = dsin(gamma)

  r(1,1) =  ca*cb*cg - sa*sg
  r(2,1) = -ca*cb*sg - sa*cg
  r(3,1) =  ca*sb
  r(1,2) =  sa*cb*cg + ca*sg
  r(2,2) = -sa*cb*sg + ca*cg
  r(3,2) =  sa*sb
  r(1,3) = -sb*cg
  r(2,3) =  sb*sg
  r(3,3) =  cb

  !      print 335, ((r(i,j),j=1,3),i=1,3)
  !  335 format((3f15.9))
end subroutine eulerm

endmodule m_fpiint
