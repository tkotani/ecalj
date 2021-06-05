      subroutine htridi(nm,n,ar,ai,d,e,e2,tau)
c
c     implicit none
      integer i,j,k,l,n,nm
      double precision ar(nm,n),ai(nm,n),d(n),e(n),e2(n),tau(2,n)
      double precision f,g,h,fi,gi,hh,si,scale,pythag
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
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      tau(1,n) = 1.0d0
      tau(2,n) = 0.0d0
c
      do   i = 1, n
         d(i) = ar(i,i)
      enddo   
      do  300  i = n, 1, -1
        l = i - 1
        h = 0.0d0
        scale = 0.0d0
        if (l .lt. 1) go to 130
c     .......... scale row (algol tol then not needed) ..........
        do   k = 1, l
           scale = scale + dabs(ar(i,k)) + dabs(ai(i,k))
        enddo   
c
        if (scale .ne. 0.0d0) go to 140
        tau(1,l) = 1.0d0
        tau(2,l) = 0.0d0
 130    continue
        e(i) = 0.0d0
        e2(i) = 0.0d0
        go to 290
c
 140    continue
        do 150  k = 1, l
          ar(i,k) = ar(i,k) / scale
          ai(i,k) = ai(i,k) / scale
          h = h + ar(i,k) * ar(i,k) + ai(i,k) * ai(i,k)
  150   continue
c
        e2(i) = scale * scale * h
        g = dsqrt(h)
        e(i) = scale * g
        f = pythag(ar(i,l),ai(i,l))
c     .......... form next diagonal element of matrix t ..........
        if (f .eq. 0.0d0) go to 160
        tau(1,l) = (ai(i,l) * tau(2,i) - ar(i,l) * tau(1,i)) / f
        si = (ar(i,l) * tau(2,i) + ai(i,l) * tau(1,i)) / f
        h = h + f * g
        g = 1.0d0 + g / f
        ar(i,l) = g * ar(i,l)
        ai(i,l) = g * ai(i,l)
        if (l .eq. 1) go to 270
        go to 170
 160    continue
        tau(1,l) = -tau(1,i)
        si = tau(2,i)
        ar(i,l) = g
 170    continue
        f = 0.0d0
        g = tau(1,l)
C ... use tau1,1:l) and tau(2,1:l) as work array ...
        do  180  j = 1, l
          tau(1,j) = 0
          tau(2,j) = 0
  180   continue
        do  220  k = 1, l
          call yyxcpy(k-1,ar(i,k),-ai(i,k),ar(k,1),ai(k,1),nm,
     .    tau,tau(2,1),2,.true.)
          call yyaxpy(l-k+1,ar(i,k),-ai(i,k),ar(k,k),ai(k,k),1,
     .    tau(1,k),tau(2,k),2,.true.)
  220   continue
        call dcopy(l,tau,2,e,1)
        tau(1,l) = g
        do  240  j = 1, l
          e(j) = e(j)/h
          tau(2,j) = tau(2,j)/h
          f = f + e(j) * ar(i,j) - tau(2,j) * ai(i,j)         
  240   continue
        hh = f / (h + h)
c     .......... form reduced a ..........
        do  261  j = 1, l
          f = ar(i,j)
          g = e(j) - hh * f
          e(j) = g
          fi = -ai(i,j)
          gi = tau(2,j) - hh * fi
          tau(2,j) = -gi
c
          do  260  k = 1, j
            ar(j,k) = ar(j,k) - f * e(k) - g * ar(i,k)
     .      + fi * tau(2,k) + gi * ai(i,k)
            ai(j,k) = ai(j,k) - f * tau(2,k) - g * ai(i,k)
     .      - fi * e(k) - gi * ar(i,k)
  260   continue
  261   continue
c
 270    continue
        do  280  k = 1, l
          ar(i,k) = scale * ar(i,k)
          ai(i,k) = scale * ai(i,k)
  280   continue
c
        tau(2,l) = -si
 290    continue
        hh = d(i)
        d(i) = ar(i,i)
        ar(i,i) = hh
        ai(i,i) = scale * dsqrt(h)
  300 continue
c
      return
      end

