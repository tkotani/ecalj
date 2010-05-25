      subroutine imtqlv(n,d,e,e2,w,ind,ierr,rv1)
c
      integer i,j,k,l,m,n,ii,mml,tag,ierr
      double precision d(n),e(n),e2(n),w(n),rv1(n)
      double precision b,c,f,g,p,r,s,tst1,tst2,pythag,d1mach,tol
      integer ind(n)
c
c     this subroutine is a variant of  imtql1  which is a translation of
c     algol procedure imtql1, num. math. 12, 377-383(1968) by martin and
c     wilkinson, as modified in num. math. 15, 450(1970) by dubrulle.
c     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
c
c     this subroutine finds the eigenvalues of a symmetric tridiagonal
c     matrix by the implicit ql method and associates with them
c     their corresponding submatrix indices.
c
c     on input
c
c        n is the order of the matrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c        e2 contains the squares of the corresponding elements of e.
c          e2(1) is arbitrary.
c
c     on output
c
c        d and e are unaltered.
c
c        elements of e2, corresponding to elements of e regarded
c          as negligible, have been replaced by zero causing the
c          matrix to split into a direct sum of submatrices.
c          e2(1) is also set to zero.
c
c        w contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct and
c          ordered for indices 1,2,...ierr-1, but may not be
c          the smallest eigenvalues.
c
c        ind contains the submatrix indices associated with the
c          corresponding eigenvalues in w -- 1 for eigenvalues
c          belonging to the first submatrix from the top,
c          2 for those belonging to the second submatrix, etc..
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c        rv1 is a temporary storage array.
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
      call tcn('imtqlv')
      ierr = 0
      k = 0
      tag = 0
      tol = 2*d1mach(3)
c
      do 100 i = 1, n
        w(i) = d(i)
        if (i .ne. 1) rv1(i-1) = e(i)
  100 continue
c
      e2(1) = 0.0d0
      rv1(n) = 0.0d0
c
      do 290 l = 1, n
        j = 0
c     .......... look for small sub-diagonal element ..........
  105   do 110 m = l, n
          if (m .eq. n) go to 120
          tst1 = dabs(w(m)) + dabs(w(m+1))
          tst2 = tst1 + dabs(rv1(m))
          if (dabs(tst2-tst1) .lt. tol*dabs(tst2)) go to 120
c     .......... guard against underflowed element of e2 ..........
          if (e2(m+1) .eq. 0.0d0) go to 125
  110   continue
c
  120   if (m .le. k) go to 130
        if (m .ne. n) e2(m+1) = 0.0d0
  125   k = m
        tag = tag + 1
  130   p = w(l)
        if (m .eq. l) go to 215
        if (j .eq. 30) go to 1000
        j = j + 1
c     .......... form shift ..........
        g = (w(l+1) - p) / (2.0d0 * rv1(l))
        r = pythag(g,1.0d0)
        g = w(m) - p + rv1(l) / (g + dsign(r,g))
        s = 1.0d0
        c = 1.0d0
        p = 0.0d0
        mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
        do 200 ii = 1, mml
          i = m - ii
          f = s * rv1(i)
          b = c * rv1(i)
          r = pythag(f,g)
          rv1(i+1) = r
          if (r .eq. 0.0d0) go to 210
          s = f / r
          c = g / r
          g = w(i+1) - p
          r = (w(i) - g) * s + 2.0d0 * c * b
          p = s * r
          w(i+1) = g + p
          g = c * r - b
  200   continue
c
        w(l) = w(l) - p
        rv1(l) = g
        rv1(m) = 0.0d0
        go to 105
c     .......... recover from underflow ..........
  210   w(i+1) = w(i+1) - p
        rv1(m) = 0.0d0
        go to 105
c     .......... order eigenvalues ..........
  215   if (l .eq. 1) go to 250
c     .......... for i=l step -1 until 2 do -- ..........
        do 230 ii = 2, l
          i = l + 2 - ii
          if (p .ge. w(i-1)) go to 270
          w(i) = w(i-1)
          ind(i) = ind(i-1)
  230   continue
c
  250   i = 1
  270   w(i) = p
        ind(i) = tag
  290 continue
c
      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 continue
      call tcx('imtqlv')
      end

