subroutine htridi(nm,n,ar,ai,d,e,e2,tau)

  !     implicit none
  integer :: i,j,k,l,n,nm
  double precision :: ar(nm,n),ai(nm,n),d(n),e(n),e2(n),tau(2,n)
  double precision :: f,g,h,fi,gi,hh,si,scale,pythag

  !     this subroutine is a translation of a complex analogue of
  !     the algol procedure tred1, num. math. 11, 181-195(1968)
  !     by martin, reinsch, and wilkinson.
  !     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).

  !     this subroutine reduces a complex hermitian matrix
  !     to a real symmetric tridiagonal matrix using
  !     unitary similarity transformations.

  !     on input

  !        nm must be set to the row dimension of two-dimensional
  !          array parameters as declared in the calling program
  !          dimension statement.

  !        n is the order of the matrix.

  !        ar and ai contain the real and imaginary parts,
  !          respectively, of the complex hermitian input matrix.
  !          only the lower triangle of the matrix need be supplied.

  !     on output

  !        ar and ai contain information about the unitary trans-
  !          formations used in the reduction in their full lower
  !          triangles.  their strict upper triangles and the
  !          diagonal of ar are unaltered.

  !        d contains the diagonal elements of the the tridiagonal matrix.

  !        e contains the subdiagonal elements of the tridiagonal
  !          matrix in its last n-1 positions.  e(1) is set to zero.

  !        e2 contains the squares of the corresponding elements of e.
  !          e2 may coincide with e if the squares are not needed.

  !        tau contains further information about the transformations.

  !     calls pythag for  dsqrt(a*a + b*b) .

  !     questions and comments should be directed to burton s. garbow,
  !     mathematics and computer science div, argonne national laboratory

  !     this version dated august 1983.

  !     ------------------------------------------------------------------

  tau(1,n) = 1.0d0
  tau(2,n) = 0.0d0

  do   i = 1, n
     d(i) = ar(i,i)
  enddo
  do  300  i = n, 1, -1
     l = i - 1
     h = 0.0d0
     scale = 0.0d0
     if (l < 1) go to 130
     !     .......... scale row (algol tol then not needed) ..........
     do   k = 1, l
        scale = scale + dabs(ar(i,k)) + dabs(ai(i,k))
     enddo

     if (scale /= 0.0d0) go to 140
     tau(1,l) = 1.0d0
     tau(2,l) = 0.0d0
130  continue
     e(i) = 0.0d0
     e2(i) = 0.0d0
     go to 290

140  continue
     do 150  k = 1, l
        ar(i,k) = ar(i,k) / scale
        ai(i,k) = ai(i,k) / scale
        h = h + ar(i,k) * ar(i,k) + ai(i,k) * ai(i,k)
150  enddo

     e2(i) = scale * scale * h
     g = dsqrt(h)
     e(i) = scale * g
     f = pythag(ar(i,l),ai(i,l))
     !     .......... form next diagonal element of matrix t ..........
     if (f == 0.0d0) go to 160
     tau(1,l) = (ai(i,l) * tau(2,i) - ar(i,l) * tau(1,i)) / f
     si = (ar(i,l) * tau(2,i) + ai(i,l) * tau(1,i)) / f
     h = h + f * g
     g = 1.0d0 + g / f
     ar(i,l) = g * ar(i,l)
     ai(i,l) = g * ai(i,l)
     if (l == 1) go to 270
     go to 170
160  continue
     tau(1,l) = -tau(1,i)
     si = tau(2,i)
     ar(i,l) = g
170  continue
     f = 0.0d0
     g = tau(1,l)
     ! ... use tau1,1:l) and tau(2,1:l) as work array ...
     !        do  180  j = 1, l
     !          tau(1,j) = 0
     !          tau(2,j) = 0
     !     180   continue
     tau = 0d0
     do  220  k = 1, l
        !          call yyxcpy(k-1,ar(i,k),-ai(i,k),ar(k,1),ai(k,1),nm,
        !     .          tau,tau(2,1),2,.true.)
        !          call yyaxpy(l-k+1,ar(i,k),-ai(i,k),ar(k,k),ai(k,k),1,
        !     .         tau(1,k),tau(2,k),2,.true.)
        tau(1,1:k-1) =  ar(i,k)*ar(k,1:k-1) + ai(i,k)*ai(k,1:k-1) + tau(1,1:k-1)
        tau(2,1:k-1) = -ar(i,k)*ai(k,1:k-1) - ai(i,k)*ar(k,1:k-1) + tau(2,1:k-1)
        tau(1,k:l)   =  ar(i,k)*ar(k,k:l)   + ai(i,k)*ai(k,k:l)   + tau(1,k:l)
        tau(2,k:l)   =  ar(i,k)*ai(k,k:l)   - ai(i,k)*ar(k,k:l)   + tau(2,k:l)
220  enddo
     !        call dcopy(l,tau,2,e,1)
     e(1:l) = tau(1,1:l)
     tau(1,l) = g
     do  240  j = 1, l
        e(j) = e(j)/h
        tau(2,j) = tau(2,j)/h
        f = f + e(j) * ar(i,j) - tau(2,j) * ai(i,j)
240  enddo
     hh = f / (h + h)
     !     .......... form reduced a ..........
     do  261  j = 1, l
        f = ar(i,j)
        g = e(j) - hh * f
        e(j) = g
        fi = -ai(i,j)
        gi = tau(2,j) - hh * fi
        tau(2,j) = -gi

        do  260  k = 1, j
           ar(j,k) = ar(j,k) - f * e(k) - g * ar(i,k) &
                + fi * tau(2,k) + gi * ai(i,k)
           ai(j,k) = ai(j,k) - f * tau(2,k) - g * ai(i,k) &
                - fi * e(k) - gi * ar(i,k)
260     enddo
261  enddo

270  continue
     do  280  k = 1, l
        ar(i,k) = scale * ar(i,k)
        ai(i,k) = scale * ai(i,k)
280  enddo

     tau(2,l) = -si
290  continue
     hh = d(i)
     d(i) = ar(i,i)
     ar(i,i) = hh
     ai(i,i) = scale * dsqrt(h)
300 enddo

  return
end subroutine htridi

