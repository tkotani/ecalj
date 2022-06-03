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

