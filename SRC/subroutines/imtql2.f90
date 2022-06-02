subroutine imtql2(nm,n,d,e,z,ierr)

  integer :: i,j,k,l,m,n,ii,nm,mml,ierr
  double precision :: d(n),e(n),z(nm,n)
  double precision :: b,c,f,g,p,r,s,tst1,tst2,pythag,d1mach,tol

  !     this subroutine is a translation of the algol procedure imtql2,
  !     num. math. 12, 377-383(1968) by martin and wilkinson,
  !     as modified in num. math. 15, 450(1970) by dubrulle.
  !     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).

  !     this subroutine finds the eigenvalues and eigenvectors
  !     of a symmetric tridiagonal matrix by the implicit ql method.
  !     the eigenvectors of a full symmetric matrix can also
  !     be found if  tred2  has been used to reduce this
  !     full matrix to tridiagonal form.

  !     on input

  !        nm must be set to the row dimension of two-dimensional
  !          array parameters as declared in the calling program
  !          dimension statement.

  !        n is the order of the matrix.

  !        d contains the diagonal elements of the input matrix.

  !        e contains the subdiagonal elements of the input matrix
  !          in its last n-1 positions.  e(1) is arbitrary.

  !        z contains the transformation matrix produced in the
  !          reduction by  tred2, if performed.  if the eigenvectors
  !          of the tridiagonal matrix are desired, z must contain
  !          the identity matrix.

  !      on output

  !        d contains the eigenvalues in ascending order.  if an
  !          error exit is made, the eigenvalues are correct but
  !          unordered for indices 1,2,...,ierr-1.

  !        e has been destroyed.

  !        z contains orthonormal eigenvectors of the symmetric
  !          tridiagonal (or full) matrix.  if an error exit is made,
  !          z contains the eigenvectors associated with the stored
  !          eigenvalues.

  !        ierr is set to
  !          zero       for normal return,
  !          j          if the j-th eigenvalue has not been
  !                     determined after 30 iterations.

  !     calls pythag for  dsqrt(a*a + b*b) .

  !     questions and comments should be directed to burton s. garbow,
  !     mathematics and computer science div, argonne national laboratory

  !     this version dated august 1983.

  !     ------------------------------------------------------------------

  call tcn('imtql2')

  ierr = 0
  if (n == 1) go to 1001
  tol = 2*d1mach(3)

  do  i = 2, n
     e(i-1) = e(i)
  enddo

  e(n) = 0.0d0

  do 240 l = 1, n
     j = 0
     !     .......... look for small sub-diagonal element ..........
105  do 110 m = l, n
        if (m == n) go to 120
        tst1 = dabs(d(m)) + dabs(d(m+1))
        tst2 = tst1 + dabs(e(m))
        if (dabs(tst2-tst1) < tol*dabs(tst2)) go to 120
110  enddo

120  p = d(l)
     if (m == l) go to 240
     if (j == 30) go to 1000
     j = j + 1
     !     .......... form shift ..........
     g = (d(l+1) - p) / (2.0d0 * e(l))
     r = pythag(g,1.0d0)
     g = d(m) - p + e(l) / (g + dsign(r,g))
     s = 1.0d0
     c = 1.0d0
     p = 0.0d0
     mml = m - l
     !     .......... for i=m-1 step -1 until l do -- ..........
     do 200 ii = 1, mml
        i = m - ii
        f = s * e(i)
        b = c * e(i)
        r = pythag(f,g)
        e(i+1) = r
        if (r == 0.0d0) go to 210
        s = f / r
        c = g / r
        g = d(i+1) - p
        r = (d(i) - g) * s + 2.0d0 * c * b
        p = s * r
        d(i+1) = g + p
        g = c * r - b
        !     .......... form vector ..........
        do 180 k = 1, n
           f = z(k,i+1)
           z(k,i+1) = s * z(k,i) + c * f
           z(k,i) = c * z(k,i) - s * f
180     enddo

200  enddo

     d(l) = d(l) - p
     e(l) = g
     e(m) = 0.0d0
     go to 105
     !     .......... recover from underflow ..........
210  d(i+1) = d(i+1) - p
     e(m) = 0.0d0
     go to 105
240 enddo
  !     .......... order eigenvalues and eigenvectors ..........
  do 300 ii = 2, n
     i = ii - 1
     k = i
     p = d(i)

     do 260 j = ii, n
        if (d(j) >= p) go to 260
        k = j
        p = d(j)
260  enddo

     if (k == i) go to 300
     d(k) = d(i)
     d(i) = p

     do 280 j = 1, n
        p = z(j,i)
        z(j,i) = z(j,k)
        z(j,k) = p
280  enddo

300 enddo

  go to 1001
  !     .......... set error -- no convergence to an
  !                eigenvalue after 30 iterations ..........
1000 ierr = l
1001 continue
  call tcx('imtql2')
end subroutine imtql2

