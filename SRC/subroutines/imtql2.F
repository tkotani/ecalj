      subroutine imtql2(nm,n,d,e,z,ierr)
c
      integer i,j,k,l,m,n,ii,nm,mml,ierr
      double precision d(n),e(n),z(nm,n)
      double precision b,c,f,g,p,r,s,tst1,tst2,pythag,d1mach,tol
c
c     this subroutine is a translation of the algol procedure imtql2,
c     num. math. 12, 377-383(1968) by martin and wilkinson,
c     as modified in num. math. 15, 450(1970) by dubrulle.
c     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
c
c     this subroutine finds the eigenvalues and eigenvectors
c     of a symmetric tridiagonal matrix by the implicit ql method.
c     the eigenvectors of a full symmetric matrix can also
c     be found if  tred2  has been used to reduce this
c     full matrix to tridiagonal form.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c        z contains the transformation matrix produced in the
c          reduction by  tred2, if performed.  if the eigenvectors
c          of the tridiagonal matrix are desired, z must contain
c          the identity matrix.
c
c      on output
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct but
c          unordered for indices 1,2,...,ierr-1.
c
c        e has been destroyed.
c
c        z contains orthonormal eigenvectors of the symmetric
c          tridiagonal (or full) matrix.  if an error exit is made,
c          z contains the eigenvectors associated with the stored
c          eigenvalues.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      call tcn('imtql2')

      ierr = 0
      if (n .eq. 1) go to 1001
      tol = 2*d1mach(3)
c
      do  i = 2, n
         e(i-1) = e(i)
      enddo   
c
      e(n) = 0.0d0
c
      do 240 l = 1, n
        j = 0
c     .......... look for small sub-diagonal element ..........
  105   do 110 m = l, n
          if (m .eq. n) go to 120
          tst1 = dabs(d(m)) + dabs(d(m+1))
          tst2 = tst1 + dabs(e(m))
          if (dabs(tst2-tst1) .lt. tol*dabs(tst2)) go to 120
  110   continue
c
  120   p = d(l)
        if (m .eq. l) go to 240
        if (j .eq. 30) go to 1000
        j = j + 1
c     .......... form shift ..........
        g = (d(l+1) - p) / (2.0d0 * e(l))
        r = pythag(g,1.0d0)
        g = d(m) - p + e(l) / (g + dsign(r,g))
        s = 1.0d0
        c = 1.0d0
        p = 0.0d0
        mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
        do 200 ii = 1, mml
          i = m - ii
          f = s * e(i)
          b = c * e(i)
          r = pythag(f,g)
          e(i+1) = r
          if (r .eq. 0.0d0) go to 210
          s = f / r
          c = g / r
          g = d(i+1) - p
          r = (d(i) - g) * s + 2.0d0 * c * b
          p = s * r
          d(i+1) = g + p
          g = c * r - b
c     .......... form vector ..........
          do 180 k = 1, n
            f = z(k,i+1)
            z(k,i+1) = s * z(k,i) + c * f
            z(k,i) = c * z(k,i) - s * f
  180     continue
c
  200   continue
c
        d(l) = d(l) - p
        e(l) = g
        e(m) = 0.0d0
        go to 105
c     .......... recover from underflow ..........
  210   d(i+1) = d(i+1) - p
        e(m) = 0.0d0
        go to 105
  240 continue
c     .......... order eigenvalues and eigenvectors ..........
      do 300 ii = 2, n
        i = ii - 1
        k = i
        p = d(i)
c
        do 260 j = ii, n
          if (d(j) .ge. p) go to 260
          k = j
          p = d(j)
  260   continue
c
        if (k .eq. i) go to 300
        d(k) = d(i)
        d(i) = p
c
        do 280 j = 1, n
          p = z(j,i)
          z(j,i) = z(j,k)
          z(j,k) = p
  280   continue
c
  300 continue
c
      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 continue
      call tcx('imtql2')
      end

