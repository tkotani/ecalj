      subroutine dtridx(nm,n,ar,d,e,e2,tau)
c
C     implicit none
C     implicit real*8 (v)
      integer i,j,k,l,n,nm
      double precision ar(nm,n),d(n),e(n),e2(n),tau(n,2)
      double precision f,g,h,hh,scale
c
c     this subroutine is a translation of a complex analogue of
c     the algol procedure tred1, num. math. 11, 181-195(1968)
c     by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine reduces a complex hermitian matrix
c     to a real symmetric tridiagonal matrix using
c     unitary similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        ar and ai contain the real and imaginary parts,
c          respectively, of the complex hermitian input matrix.
c          only the lower triangle of the matrix need be supplied.
c
c     on output
c
c        ar and ai contain information about the unitary trans-
c          formations used in the reduction in their full lower
c          triangles.  their strict upper triangles and the
c          diagonal of ar are unaltered.
c
c        d contains the diagonal elements of the the tridiagonal matrix.
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is set to zero.
c
c        e2 contains the squares of the corresponding elements of e.
c          e2 may coincide with e if the squares are not needed.
c
c        tau contains further information about the transformations.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version adapted from august 1983 htridi.  Differences
C     with htridi are that indices of tau and a are permuted (uses
C     upper triangle of a)
c     ------------------------------------------------------------------
c
      call tcn('dtridx')
      tau(n,1) = 1.0d0
c
      do  100  i = 1, n
  100 d(i) = ar(i,i)
      do  300  i = n, 1, -1
         l = i - 1
         h = 0.0d0
         scale = 0.0d0
         if (l .lt. 1) go to 130
c     .......... scale row (algol tol then not needed) ..........
         do  120  k = 1, l
  120    scale = scale + dabs(ar(k,i))
         if (scale .ne. 0.0d0) go to 140
         tau(l,1) = 1.0d0
  130    e(i) = 0.0d0
         e2(i) = 0.0d0
         go to 290
c
  140    do  150  k = 1, l
            ar(k,i) = ar(k,i) / scale
            h = h + ar(k,i)**2
  150    continue
c
         e2(i) = scale * scale * h
         g = dsqrt(h)
         e(i) = scale * g
         f = dabs(ar(l,i))
c     .......... form next diagonal element of matrix t ..........
         if (f .eq. 0.0d0) go to 160
         tau(l,1) = ( - ar(l,i) * tau(i,1)) / f
         h = h + f * g
         g = 1.0d0 + g / f
         ar(l,i) = g * ar(l,i)
         if (l .eq. 1) go to 270
         go to 170
  160    tau(l,1) = -tau(i,1)
         ar(l,i) = g
  170    f = 0.0d0
         do  240  j = 1, l
c     .......... form element of a*u ..........
            g = 0.d0
            do  180  k = 1, j
               g = g + ar(k,j) * ar(k,i)
  180       continue
            do  200  k = j+1, l
               g = g + ar(j,k) * ar(k,i)
  200       continue
c     .......... form element of p ..........
  220       e(j) = g / h
            f = f + e(j) * ar(j,i)
  240    continue
c
         hh = f / (h + h)
c     .......... form reduced a ..........
         do  260  j = 1, l
            f = ar(j,i)
            g = e(j) - hh * f
            e(j) = g
c
            do  260  k = 1, j
               ar(k,j) = ar(k,j) - f * e(k) - g * ar(k,i)
  260    continue
  270    do  280  k = 1, l
            ar(k,i) = scale * ar(k,i)
  280    continue
  290    hh = d(i)
         d(i) = ar(i,i)
         tau(i,2) = -scale * dsqrt(h)
  300 continue
      call tcx('dtridx')
      end

