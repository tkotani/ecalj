subroutine tred2(nm,n,a,d,e,z)

  integer :: i,j,k,l,n,ii,nm,jp1
  double precision :: a(nm,n),d(n),e(n),z(nm,n)
  double precision :: f,g,h,hh,scale

  !     this subroutine is a translation of the algol procedure tred2,
  !     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
  !     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).

  !     this subroutine reduces a real symmetric matrix to a
  !     symmetric tridiagonal matrix using and accumulating
  !     orthogonal similarity transformations.

  !     on input

  !        nm must be set to the row dimension of two-dimensional
  !          array parameters as declared in the calling program
  !          dimension statement.

  !        n is the order of the matrix.

  !        a contains the real symmetric input matrix.  only the
  !          lower triangle of the matrix need be supplied.

  !     on output

  !        d contains the diagonal elements of the tridiagonal matrix.

  !        e contains the subdiagonal elements of the tridiagonal
  !          matrix in its last n-1 positions.  e(1) is set to zero.

  !        z contains the orthogonal transformation matrix
  !          produced in the reduction.

  !        a and z may coincide.  if distinct, a is unaltered.

  !     questions and comments should be directed to burton s. garbow,
  !     mathematics and computer science div, argonne national laboratory

  !     this version dated august 1983.

  !     ------------------------------------------------------------------

  do 100 i = 1, n

     do j = i, n
        z(j,i) = a(j,i)
     enddo

     d(i) = a(n,i)
100 enddo

  if (n == 1) go to 510
  !     .......... for i=n step -1 until 2 do -- ..........
  do 300 ii = 2, n
     i = n + 2 - ii
     l = i - 1
     h = 0.0d0
     scale = 0.0d0
     if (l < 2) go to 130
     !     .......... scale row (algol tol then not needed) ..........
     do k = 1, l
        scale = scale + dabs(d(k))
     enddo

     if (scale /= 0.0d0) go to 140
130  e(i) = d(l)

     do 135 j = 1, l
        d(j) = z(l,j)
        z(i,j) = 0.0d0
        z(j,i) = 0.0d0
135  enddo

     go to 290

140  do 150 k = 1, l
        d(k) = d(k) / scale
        h = h + d(k) * d(k)
150  enddo

     f = d(l)
     g = -dsign(dsqrt(h),f)
     e(i) = scale * g
     h = h - f * g
     d(l) = f - g
     !     .......... form a*u ..........
     do j = 1, l
        e(j) = 0.0d0
     enddo

     do 240 j = 1, l
        f = d(j)
        z(j,i) = f
        g = e(j) + z(j,j) * f
        jp1 = j + 1
        if (l < jp1) go to 220

        do 200 k = jp1, l
           g = g + z(k,j) * d(k)
           e(k) = e(k) + z(k,j) * f
200     enddo

220     e(j) = g
240  enddo
     !     .......... form p ..........
     f = 0.0d0

     do 245 j = 1, l
        e(j) = e(j) / h
        f = f + e(j) * d(j)
245  enddo

     hh = f / (h + h)
     !     .......... form q ..........
     do  j = 1, l
        e(j) = e(j) - hh * d(j)
     enddo
     !     .......... form reduced a ..........
     do 280 j = 1, l
        f = d(j)
        g = e(j)

        do  k = j, l
           z(k,j) = z(k,j) - f * e(k) - g * d(k)
        enddo

        d(j) = z(l,j)
        z(i,j) = 0.0d0
280  enddo

290  d(i) = h
300 enddo
  !     .......... accumulation of transformation matrices ..........
  do 500 i = 2, n
     l = i - 1
     z(n,l) = z(l,l)
     z(l,l) = 1.0d0
     h = d(i)
     if (h == 0.0d0) go to 380

     do  k = 1, l
        d(k) = z(k,i) / h
     enddo

     do 361 j = 1, l
        g = 0.0d0

        do k = 1, l
           g = g + z(k,i) * z(k,j)
        enddo

        do 360 k = 1, l
           z(k,j) = z(k,j) - g * d(k)
360     enddo
361  enddo

380  continue
     do k = 1, l
        z(k,i) = 0.0d0
     enddo

500 enddo

510 continue
  do 520 i = 1, n
     d(i) = z(n,i)
     z(n,i) = 0.0d0
520 enddo

  z(n,n) = 1.0d0
  e(1) = 0.0d0
  return
end subroutine tred2
