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

