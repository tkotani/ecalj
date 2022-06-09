
subroutine tqlrat(n,d,e2,ierr)
  integer :: i,j,l,m,n,ii,l1,mml,ierr
  double precision :: d(n),e2(n)
  double precision :: b=1d99,c=1d99,f,g=1d99,h,p,r,s,t,epslon,pythag
  !     this subroutine is a translation of the algol procedure tqlrat,
  !     algorithm 464, comm. acm 16, 689(1973) by reinsch.

  !     this subroutine finds the eigenvalues of a symmetric
  !     tridiagonal matrix by the rational ql method.

  !     on input

  !        n is the order of the matrix.

  !        d contains the diagonal elements of the input matrix.

  !        e2 contains the squares of the subdiagonal elements of the
  !          input matrix in its last n-1 positions.  e2(1) is arbitrary.

  !      on output

  !        d contains the eigenvalues in ascending order.  if an
  !          error exit is made, the eigenvalues are correct and
  !          ordered for indices 1,2,...ierr-1, but may not be
  !          the smallest eigenvalues.

  !        e2 has been destroyed.

  !        ierr is set to
  !          zero       for normal return,
  !          j          if the j-th eigenvalue has not been
  !                     determined after 30 iterations.

  !     calls pythag for  dsqrt(a*a + b*b) .

  !     questions and comments should be directed to burton s. garbow,
  !     mathematics and computer science div, argonne national laboratory

  !     this version dated august 1983.

  !     ------------------------------------------------------------------

  ierr = 0
  if (n == 1) go to 1001

  do  i = 2, n
     e2(i-1) = e2(i)
  enddo

  f = 0.0d0
  t = 0.0d0
  e2(n) = 0.0d0

  do 290 l = 1, n
     j = 0
     h = dabs(d(l)) + dsqrt(e2(l))
     if (t > h) go to 105
     t = h
     b = epslon(t)
     c = b * b
     !     .......... look for small squared sub-diagonal element ..........
105  do 110 m = l, n
        if (e2(m) <= c) go to 120
        !     .......... e2(n) is always zero, so there is no exit
        !                through the bottom of the loop ..........
110  enddo

120  if (m == l) go to 210
130  if (j == 30) go to 1000
     j = j + 1
     !     .......... form shift ..........
     l1 = l + 1
     s = dsqrt(e2(l))
     g = d(l)
     p = (d(l1) - g) / (2.0d0 * s)
     r = pythag(p,1.0d0)
     d(l) = s / (p + dsign(r,p))
     h = g - d(l)

     do i = l1, n
        d(i) = d(i) - h
     enddo

     f = f + h
     !     .......... rational ql transformation ..........
     g = d(m)
     if (g == 0.0d0) g = b
     h = g
     s = 0.0d0
     mml = m - l
     !     .......... for i=m-1 step -1 until l do -- ..........
     do 200 ii = 1, mml
        i = m - ii
        p = g * h
        r = p + e2(i)
        e2(i+1) = s * r
        s = e2(i) / r
        d(i+1) = h + s * (h + d(i))
        g = d(i) - e2(i) / g
        if (g == 0.0d0) g = b
        h = g * p / r
200  enddo

     e2(l) = s * g
     d(l) = h
     !     .......... guard against underflow in convergence test ..........
     if (h == 0.0d0) go to 210
     if (dabs(e2(l)) <= dabs(c/h)) go to 210
     e2(l) = h * e2(l)
     if (e2(l) /= 0.0d0) go to 130
210  p = d(l) + f
     !     .......... order eigenvalues ..........
     if (l == 1) go to 250
     !     .......... for i=l step -1 until 2 do -- ..........
     do 230 ii = 2, l
        i = l + 2 - ii
        if (p >= d(i-1)) go to 270
        d(i) = d(i-1)
230  enddo

250  i = 1
270  d(i) = p
290 enddo
  go to 1001
  !     .......... set error -- no convergence to an
  !                eigenvalue after 30 iterations ..........
1000 ierr = l
1001 return
end subroutine tqlrat

subroutine tinvit(nm,n,d,e,e2,m,w,ind,z, &
     ierr,rv1,rv2,rv3,rv4,rv6)
  integer :: i,j,m,n,p,q,r,s,ii,ip,jj,nm,its,tag,ierr,group
  double precision :: d(n),e(n),e2(n),w(m),z(nm,m), &
       rv1(n),rv2(n),rv3(n),rv4(n),rv6(n)
  double precision :: u,v,uk,xu,x0,x1,eps2,eps3,eps4,norm,order,epslon, &
       pythag
  integer :: ind(m)

  !     this subroutine is a translation of the inverse iteration tech-
  !     nique in the algol procedure tristurm by peters and wilkinson.
  !     handbook for auto. comp., vol.ii-linear algebra, 418-439(1971).

  !     this subroutine finds those eigenvectors of a tridiagonal
  !     symmetric matrix corresponding to specified eigenvalues,
  !     using inverse iteration.

  !     on input

  !        nm must be set to the row dimension of two-dimensional
  !          array parameters as declared in the calling program
  !          dimension statement.

  !        n is the order of the matrix.

  !        d contains the diagonal elements of the input matrix.

  !        e contains the subdiagonal elements of the input matrix
  !          in its last n-1 positions.  e(1) is arbitrary.

  !        e2 contains the squares of the corresponding elements of e,
  !          with zeros corresponding to negligible elements of e.
  !          e(i) is considered negligible if it is not larger than
  !          the product of the relative machine precision and the sum
  !          of the magnitudes of d(i) and d(i-1).  e2(1) must contain
  !          0.0d0 if the eigenvalues are in ascending order, or 2.0d0
  !          if the eigenvalues are in descending order.  if  bisect,
  !          tridib, or  imtqlv  has been used to find the eigenvalues,
  !          their output e2 array is exactly what is expected here.

  !        m is the number of specified eigenvalues.

  !        w contains the m eigenvalues in ascending or descending order.

  !        ind contains in its first m positions the submatrix indices
  !          associated with the corresponding eigenvalues in w --
  !          1 for eigenvalues belonging to the first submatrix from
  !          the top, 2 for those belonging to the second submatrix, etc.

  !     on output

  !        all input arrays are unaltered.

  !        z contains the associated set of orthonormal eigenvectors.
  !          any vector which fails to converge is set to zero.

  !        ierr is set to
  !          zero       for normal return,
  !          -r         if the eigenvector corresponding to the r-th
  !                     eigenvalue fails to converge in 5 iterations.

  !        rv1, rv2, rv3, rv4, and rv6 are temporary storage arrays.

  !     calls pythag for  dsqrt(a*a + b*b) .

  !     questions and comments should be directed to burton s. garbow,
  !     mathematics and computer science div, argonne national laboratory

  !     this version dated august 1983.

  !     ------------------------------------------------------------------

  call tcn('tinvit')

  ierr = 0
  if (m == 0) go to 1001
  tag = 0
  order = 1.0d0 - e2(1)
  q = 0
  !     .......... establish and process next submatrix ..........
100 p = q + 1

  do 120 q = p, n
     if (q == n) go to 140
     if (e2(q+1) == 0.0d0) go to 140
120 enddo
  !     .......... find vectors by inverse iteration ..........
140 tag = tag + 1
  s = 0

  do 920 r = 1, m
     if (ind(r) /= tag) go to 920
     its = 1
     x1 = w(r)
     if (s /= 0) go to 510
     !     .......... check for isolated root ..........
     xu = 1.0d0
     if (p /= q) go to 490
     rv6(p) = 1.0d0
     go to 870
490  continue
     norm = dabs(d(p))
     ip = p + 1

     do  i = ip, q
        norm = dmax1(norm, dabs(d(i))+dabs(e(i)))
     enddo
     !     .......... eps2 is the criterion for grouping,
     !                eps3 replaces zero pivots and equal
     !                roots are modified by eps3,
     !                eps4 is taken very small to avoid overflow ..........
     eps2 = 1.0d-3 * norm
     eps3 = epslon(norm)
     uk = q - p + 1
     eps4 = uk * eps3
     uk = eps4 / dsqrt(uk)
     s = p
505  group = 0
     go to 520
     !     .......... look for close or coincident roots ..........
510  if (dabs(x1-x0) >= eps2) go to 505
     group = group + 1
     if (order * (x1 - x0) <= 0.0d0) x1 = x0 + order * eps3
     !     .......... elimination with interchanges and
     !                initialization of vector ..........
520  v = 0.0d0

     do 580 i = p, q
        rv6(i) = uk
        if (i == p) go to 560
        if (dabs(e(i)) < dabs(u)) go to 540
        !     .......... warning -- a divide check may occur here if
        !                e2 array has not been specified correctly ..........
        xu = u / e(i)
        rv4(i) = xu
        rv1(i-1) = e(i)
        rv2(i-1) = d(i) - x1
        rv3(i-1) = 0.0d0
        if (i /= q) rv3(i-1) = e(i+1)
        u = v - xu * rv2(i-1)
        v = -xu * rv3(i-1)
        go to 580
540     xu = e(i) / u
        rv4(i) = xu
        rv1(i-1) = u
        rv2(i-1) = v
        rv3(i-1) = 0.0d0
560     u = d(i) - x1 - xu * v
        if (i /= q) v = e(i+1)
580  enddo

     if (u == 0.0d0) u = eps3
     rv1(q) = u
     rv2(q) = 0.0d0
     rv3(q) = 0.0d0
     !     .......... back substitution
     !                for i=q step -1 until p do -- ..........
600  do 620 ii = p, q
        i = p + q - ii
        rv6(i) = (rv6(i) - u * rv2(i) - v * rv3(i)) / rv1(i)
        v = u
        u = rv6(i)
620  enddo
     !     .......... orthogonalize with respect to previous
     !                members of group ..........
     if (group == 0) go to 700
     j = r

     do 680 jj = 1, group
630     j = j - 1
        if (ind(j) /= tag) go to 630
        xu = 0.0d0

        do i = p, q
           xu = xu + rv6(i) * z(i,j)
        enddo

        do i = p, q
           rv6(i) = rv6(i) - xu * z(i,j)
        enddo
680  enddo

700  norm = 0.0d0

     do i = p, q
        norm = norm + dabs(rv6(i))
     enddo

     if (norm >= 1.0d0) go to 840
     !     .......... forward substitution ..........
     if (its == 5) go to 830
     if (norm /= 0.0d0) go to 740
     rv6(s) = eps4
     s = s + 1
     if (s > q) s = p
     go to 780
740  continue
     xu = eps4 / norm

     do i = p, q
        rv6(i) = rv6(i) * xu
     enddo
     !     .......... elimination operations on next vector
     !                iterate ..........
780  continue
     do 820 i = ip, q
        u = rv6(i)
        !     .......... if rv1(i-1) .eq. e(i), a row interchange
        !                was performed earlier in the
        !                triangularization process ..........
        if (rv1(i-1) /= e(i)) go to 800
        u = rv6(i-1)
        rv6(i-1) = rv6(i)
800     continue
        rv6(i) = u - rv4(i) * rv6(i-1)
820  enddo

     its = its + 1
     go to 600
     !     .......... set error -- non-converged eigenvector ..........
830  continue
     ierr = -r
     xu = 0.0d0
     go to 870
     !     .......... normalize so that sum of squares is
     !                1 and expand to full order ..........
840  u = 0.0d0

     do i = p, q
        u = pythag(u,rv6(i))
     enddo

     xu = 1.0d0 / u

870  continue
     do i = 1, n
        z(i,r) = 0.0d0
     enddo

     do i = p, q
        z(i,r) = rv6(i) * xu
     enddo

     x0 = x1
920 enddo

  if (q < n) go to 100
1001 continue
  call tcx('tinvit')
end subroutine tinvit

real(8) function epslon (x)
  double precision :: x
  double precision :: a,b,c,eps

  !     this program should function properly on all systems
  !     satisfying the following two assumptions,
  !        1.  the base used in representing floating point
  !            numbers is not a power of three.
  !        2.  the quantity  a  in statement 10 is represented to
  !            the accuracy used in floating point variables
  !            that are stored in memory.
  !     the statement number 10 and the go to 10 are intended to
  !     force optimizing compilers to generate code satisfying
  !     assumption 2.
  !     under these assumptions, it should be true that,
  !            a  is not exactly equal to four-thirds,
  !            b  has a zero for its last bit or digit,
  !            c  is not exactly equal to one,
  !            eps  measures the separation of 1.0 from
  !                 the next larger floating point number.
  !     the developers of eispack would appreciate being informed
  !     about any systems where these assumptions do not hold.

  !     this version dated 4/6/83.

  a = 4.0d0/3.0d0
10 b = a - 1.0d0
  c = b + b + b
  eps = dabs(c-1.0d0)
  if (eps == 0.0d0) go to 10
  epslon = eps*dabs(x)
  return
END function epslon

double precision function pythag(a,b)
  !     implicit none
  double precision :: a,b
  !     finds dsqrt(a**2+b**2) without overflow or destructive underflow
  double precision :: p,r,s,t,u
  p = dmax1(dabs(a),dabs(b))
  if (p == 0.0d0) go to 20
  r = (dmin1(dabs(a),dabs(b))/p)**2
10 continue
  t = 4.0d0 + r
  if (t == 4.0d0) go to 20
  s = r/t
  u = 1.0d0 + 2.0d0*s
  p = u*p
  r = (s/u)**2 * r
  go to 10
20 pythag = p
  return
END function pythag

