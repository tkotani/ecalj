subroutine minfit(nm,m,n,a,w,ip,b,ierr,rv1)
  integer :: i,j,k,l,m,n,ii,ip,i1,kk,k1,ll,l1=9999,m1,nm,its,ierr
  double precision :: a(nm,n),w(n),b(nm,ip),rv1(n)
  double precision :: c,f,g,h,s,x,y,z,tst1,tst2,scale,pythag

  !     this subroutine is a translation of the algol procedure minfit,
  !     num. math. 14, 403-420(1970) by golub and reinsch.
  !     handbook for auto. comp., vol ii-linear algebra, 134-151(1971).

  !     this subroutine determines, towards the solution of the linear
  !                                                        t
  !     system ax=b, the singular value decomposition a=usv  of a real
  !                                         t
  !     m by n rectangular matrix, forming u b rather than u.  householder
  !     bidiagonalization and a variant of the qr algorithm are used.

  !     on input

  !        nm must be set to the row dimension of two-dimensional
  !          array parameters as declared in the calling program
  !          dimension statement.  note that nm must be at least
  !          as large as the maximum of m and n.

  !        m is the number of rows of a and b.

  !        n is the number of columns of a and the order of v.

  !        a contains the rectangular coefficient matrix of the system.

  !        ip is the number of columns of b.  ip can be zero.

  !        b contains the constant column matrix of the system
  !          if ip is not zero.  otherwise b is not referenced.

  !     on output

  !        a has been overwritten by the matrix v (orthogonal) of the
  !          decomposition in its first n rows and columns.  if an
  !          error exit is made, the columns of v corresponding to
  !          indices of correct singular values should be correct.

  !        w contains the n (non-negative) singular values of a (the
  !          diagonal elements of s).  they are unordered.  if an
  !          error exit is made, the singular values should be correct
  !          for indices ierr+1,ierr+2,...,n.

  !                                   t
  !        b has been overwritten by u b.  if an error exit is made,
  !                       t
  !          the rows of u b corresponding to indices of correct
  !          singular values should be correct.

  !        ierr is set to
  !          zero       for normal return,
  !          k          if the k-th singular value has not been
  !                     determined after 30 iterations.

  !        rv1 is a temporary storage array.

  !     calls pythag for  dsqrt(a*a + b*b) .

  !     questions and comments should be directed to burton s. garbow,
  !     mathematics and computer science div, argonne national laboratory

  !     this version dated august 1983.

  !     ------------------------------------------------------------------

  ierr = 0
  !     .......... householder reduction to bidiagonal form ..........
  g = 0.0d0
  scale = 0.0d0
  x = 0.0d0

  do 300 i = 1, n
     l = i + 1
     rv1(i) = scale * g
     g = 0.0d0
     s = 0.0d0
     scale = 0.0d0
     if (i > m) go to 210

     do k = i, m
        scale = scale + dabs(a(k,i))
     enddo

     if (scale == 0.0d0) go to 210

     do 130 k = i, m
        a(k,i) = a(k,i) / scale
        s = s + a(k,i)**2
130  enddo

     f = a(i,i)
     g = -dsign(dsqrt(s),f)
     h = f * g - s
     a(i,i) = f - g
     if (i == n) go to 160

     do 151 j = l, n
        s = 0.0d0

        do  k = i, m
           s = s + a(k,i) * a(k,j)
        enddo

        f = s / h

        do 150 k = i, m
           a(k,j) = a(k,j) + f * a(k,i)
150     enddo
151  enddo

160  if (ip == 0) go to 190

     do 181 j = 1, ip
        s = 0.0d0

        do  k = i, m
           s = s + a(k,i) * b(k,j)
        enddo

        f = s / h

        do 180 k = i, m
           b(k,j) = b(k,j) + f * a(k,i)
180     enddo
181  enddo

190  continue
     do k = i, m
        a(k,i) = scale * a(k,i)
     enddo

210  w(i) = scale * g
     g = 0.0d0
     s = 0.0d0
     scale = 0.0d0
     if (i > m .OR. i == n) go to 290

     do k = l, n
        scale = scale + dabs(a(i,k))
     enddo

     if (scale == 0.0d0) go to 290

     do 230 k = l, n
        a(i,k) = a(i,k) / scale
        s = s + a(i,k)**2
230  enddo

     f = a(i,l)
     g = -dsign(dsqrt(s),f)
     h = f * g - s
     a(i,l) = f - g

     do  k = l, n
        rv1(k) = a(i,k) / h
     enddo

     if (i == m) go to 270

     do 261 j = l, m
        s = 0.0d0

        do k = l, n
           s = s + a(j,k) * a(i,k)
        enddo

        do 260 k = l, n
           a(j,k) = a(j,k) + s * rv1(k)
260     enddo
261  enddo

270  continue
     do  k = l, n
        a(i,k) = scale * a(i,k)
     enddo
290  continue
     x = dmax1(x,dabs(w(i))+dabs(rv1(i)))
300 enddo
  !     .......... accumulation of right-hand transformations.
  !                for i=n step -1 until 1 do -- ..........
  do 400 ii = 1, n
     i = n + 1 - ii
     if (i == n) go to 390
     if (g == 0.0d0) go to 360

     do j = l, n
        !     .......... double division avoids possible underflow ..........
        a(j,i) = (a(i,j) / a(i,l)) / g
     enddo

     do 350 j = l, n
        s = 0.0d0

        do k = l, n
           s = s + a(i,k) * a(k,j)
        enddo

        do k = l, n
           a(k,j) = a(k,j) + s * a(k,i)
        enddo
350  enddo

360  continue
     do j = l, n
        a(i,j) = 0.0d0
        a(j,i) = 0.0d0
     enddo

390  continue
     a(i,i) = 1.0d0
     g = rv1(i)
     l = i
400 enddo

  if (m >= n .OR. ip == 0) go to 510
  m1 = m + 1

  do  i = m1, n
     do j = 1, ip
        b(i,j) = 0.0d0
     enddo
  enddo
  !     .......... diagonalization of the bidiagonal form ..........
510 continue
  tst1 = x
  !     .......... for k=n step -1 until 1 do -- ..........
  do 700 kk = 1, n
     k1 = n - kk
     k = k1 + 1
     its = 0
     !     .......... test for splitting.
     !                for l=k step -1 until 1 do -- ..........
520  continue
     do 530 ll = 1, k
        l1 = k - ll
        l = l1 + 1
        tst2 = tst1 + dabs(rv1(l))
        if (tst2 == tst1) go to 565
        !     .......... rv1(1) is always zero, so there is no exit
        !                through the bottom of the loop ..........
        tst2 = tst1 + dabs(w(l1))
        if (tst2 == tst1) go to 540
530  enddo
     !     .......... cancellation of rv1(l) if l greater than 1 ..........
540  continue
     c = 0.0d0
     s = 1.0d0

     do 560 i = l, k
        f = s * rv1(i)
        rv1(i) = c * rv1(i)
        tst2 = tst1 + dabs(f)
        if (tst2 == tst1) go to 565
        g = w(i)
        h = pythag(f,g)
        w(i) = h
        c = g / h
        s = -f / h
        if (ip == 0) go to 560

        do 550 j = 1, ip
           y = b(l1,j)
           z = b(i,j)
           b(l1,j) = y * c + z * s
           b(i,j) = -y * s + z * c
550     enddo

560  enddo
     !     .......... test for convergence ..........
565  continue
     z = w(k)
     if (l == k) go to 650
     !     .......... shift from bottom 2 by 2 minor ..........
     if (its == 30) go to 1000
     its = its + 1
     x = w(l)
     y = w(k1)
     g = rv1(k1)
     h = rv1(k)
     f = 0.5d0 * (((g + z) / h) * ((g - z) / y) + y / h - h / y)
     g = pythag(f,1.0d0)
     f = x - (z / x) * z + (h / x) * (y / (f + dsign(g,f)) - h)
     !     .......... next qr transformation ..........
     c = 1.0d0
     s = 1.0d0

     do 600 i1 = l, k1
        i = i1 + 1
        g = rv1(i)
        y = w(i)
        h = s * g
        g = c * g
        z = pythag(f,h)
        rv1(i1) = z
        c = f / z
        s = h / z
        f = x * c + g * s
        g = -x * s + g * c
        h = y * s
        y = y * c

        do 570 j = 1, n
           x = a(j,i1)
           z = a(j,i)
           a(j,i1) = x * c + z * s
           a(j,i) = -x * s + z * c
570     enddo

        z = pythag(f,h)
        w(i1) = z
        !     .......... rotation can be arbitrary if z is zero ..........
        if (z == 0.0d0) go to 580
        c = f / z
        s = h / z
580     continue
        f = c * g + s * y
        x = -s * g + c * y
        if (ip == 0) go to 600

        do 590 j = 1, ip
           y = b(i1,j)
           z = b(i,j)
           b(i1,j) = y * c + z * s
           b(i,j) = -y * s + z * c
590     enddo

600  enddo

     rv1(l) = 0.0d0
     rv1(k) = f
     w(k) = x
     go to 520
     !     .......... convergence ..........
650  continue
     if (z >= 0.0d0) go to 700
     !     .......... w(k) is made non-negative ..........
     w(k) = -z

     do j = 1, n
        a(j,k) = -a(j,k)
     enddo

700 enddo

  go to 1001
  !     .......... set error -- no convergence to a
  !                singular value after 30 iterations ..........
1000 ierr = k
1001 return
end subroutine minfit

