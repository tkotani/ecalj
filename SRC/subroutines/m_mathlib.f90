! Routines in m_handbook should be replaced by things in blas and lapack (2022-6-13,tk)
module m_mathlib
  public rs!,chkhss !, htridi,imtql2,htribk
  private
contains
  subroutine tql2(nm,n,d,e,z,ierr)
    integer :: i,j,k,l,m,n,ii,l1,l2,nm,mml,ierr
    double precision :: d(n),e(n),z(nm,n)
    double precision :: c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2!,pythag
    !     this subroutine is a translation of the algol procedure tql2,
    !     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
    !     wilkinson.
    !     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).

    !     this subroutine finds the eigenvalues and eigenvectors
    !     of a symmetric tridiagonal matrix by the ql method.
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

    ierr = 0
    if (n == 1) go to 1001

    do  i = 2, n
       e(i-1) = e(i)
    enddo

    f = 0.0d0
    tst1 = 0.0d0
    e(n) = 0.0d0

    do 240 l = 1, n
       j = 0
       h = dabs(d(l)) + dabs(e(l))
       if (tst1 < h) tst1 = h
       !     .......... look for small sub-diagonal element ..........
       do 110 m = l, n
          tst2 = tst1 + dabs(e(m))
          if (tst2 == tst1) go to 120
          !     .......... e(n) is always zero, so there is no exit
          !                through the bottom of the loop ..........
110    enddo

120    if (m == l) go to 220
130    if (j == 30) go to 1000
       j = j + 1
       !     .......... form shift ..........
       l1 = l + 1
       l2 = l1 + 1
       g = d(l)
       p = (d(l1) - g) / (2.0d0 * e(l))
       r = pythag(p,1.0d0)
       d(l) = e(l) / (p + dsign(r,p))
       d(l1) = e(l) * (p + dsign(r,p))
       dl1 = d(l1)
       h = g - d(l)
       if (l2 > n) go to 145

       do i = l2, n
          d(i) = d(i) - h
       enddo

145    f = f + h
       !     .......... ql transformation ..........
       p = d(m)
       c = 1.0d0
       c2 = c
       el1 = e(l1)
       s = 0.0d0
       mml = m - l
       !     .......... for i=m-1 step -1 until l do -- ..........
       do 200 ii = 1, mml
          c3 = c2
          c2 = c
          s2 = s
          i = m - ii
          g = c * e(i)
          h = c * p
          r = pythag(p,e(i))
          e(i+1) = s * r
          s = e(i) / r
          c = p / r
          p = c * d(i) - s * g
          d(i+1) = h + s * (c * g + s * d(i))
          !     .......... form vector ..........
          do 180 k = 1, n
             h = z(k,i+1)
             z(k,i+1) = s * z(k,i) + c * h
             z(k,i) = c * z(k,i) - s * h
180       enddo

200    enddo

       p = -s * s2 * c3 * el1 * e(l) / dl1
       e(l) = s * p
       d(l) = c * p
       tst2 = tst1 + dabs(e(l))
       if (tst2 > tst1) go to 130
220    d(l) = d(l) + f
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
260    enddo

       if (k == i) go to 300
       d(k) = d(i)
       d(i) = p

       do 280 j = 1, n
          p = z(j,i)
          z(j,i) = z(j,k)
          z(j,k) = p
280    enddo

300 enddo

    go to 1001
    !     .......... set error -- no convergence to an
    !                eigenvalue after 30 iterations ..........
1000 ierr = l
1001 return
  end subroutine tql2

  double precision function pythag(a,b)
    !     implicit none
    double precision :: a,b
    !     finds dsqrt(a**2+b**2) without overflow or destructive underflow
    double precision :: p,r,s,t,u
    p = dmax1(dabs(a),dabs(b))
    if (p == 0.0d0) go to 20
    r = (dmin1(dabs(a),dabs(b))/p)**2
10  continue
    t = 4.0d0 + r
    if (t == 4.0d0) go to 20
    s = r/t
    u = 1.0d0 + 2.0d0*s
    p = u*p
    r = (s/u)**2 * r
    go to 10
20  pythag = p
    return
  END function  pythag

!   subroutine htribk(nm,n,ar,ai,tau,m,zr,zi)

!     integer :: i,j,k,l,m,n,nm
!     double precision :: ar(nm,n),ai(nm,n),tau(2,n),zr(nm,m),zi(nm,m)
!     double precision :: h,s,si

!     !     this subroutine is a translation of a complex analogue of
!     !     the algol procedure trbak1, num. math. 11, 181-195(1968)
!     !     by martin, reinsch, and wilkinson.
!     !     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).

!     !     this subroutine forms the eigenvectors of a complex hermitian
!     !     matrix by back transforming those of the corresponding
!     !     real symmetric tridiagonal matrix determined by  htridi.

!     !     on input

!     !        nm must be set to the row dimension of two-dimensional
!     !          array parameters as declared in the calling program
!     !          dimension statement.

!     !        n is the order of the matrix.

!     !        ar and ai contain information about the unitary trans-
!     !          formations used in the reduction by  htridi  in their
!     !          full lower triangles except for the diagonal of ar.

!     !        tau contains further information about the transformations.

!     !        m is the number of eigenvectors to be back transformed.

!     !        zr contains the eigenvectors to be back transformed
!     !          in its first m columns.

!     !     on output

!     !        zr and zi contain the real and imaginary parts,
!     !          respectively, of the transformed eigenvectors
!     !          in their first m columns.

!     !     note that the last component of each returned vector
!     !     is real and that vector euclidean norms are preserved.

!     !     questions and comments should be directed to burton s. garbow,
!     !     mathematics and computer science div, argonne national laboratory

!     !     this version dated august 1983.
!     !     Aug 1990 MvS altered into daxpy-style loops
!     !     ------------------------------------------------------------------

!     if (m == 0) return
!     !     .......... transform the eigenvectors of the real symmetric
!     !                tridiagonal matrix to those of the hermitian
!     !                tridiagonal matrix. ..........
!     do  k = 1, n
!        do j = 1, m
!           zi(k,j) = -zr(k,j) * tau(2,k)
!           zr(k,j) = zr(k,j) * tau(1,k)
!        enddo
!     enddo

!     if (n == 1) return
!     !     .......... recover and apply the householder matrices ..........
!     do  140  i = 2, n
!        l = i - 1
!        h = ai(i,i)
!        if (h == 0.0d0) go to 140
!        !        do 100 j = 1, m
!        !          tau(1,j) = 0.0d0
!        !          tau(2,j) = 0.0d0
!        !     100   continue
!        tau=0d0
!        do  110  k = 1, l
!           !          call yyaxpy(m,ar(i,k),ai(i,k),zr(k,1),zi(k,1),nm,
!           !     .          tau,tau(2,1),2,.true.)
!           tau(1,:)= ar(i,k)*zr(k,:) - ai(i,k)*zi(k,:) + tau(1,:)
!           tau(2,:)= ar(i,k)*zi(k,:) + ai(i,k)*zr(k,:) + tau(2,:)
! 110    enddo
!        !     ... double divisions avoid possible underflow ...
!        tau=tau/h**2
!        !        call dscal(m,1/h,tau,2)
!        !        call dscal(m,1/h,tau,2)
!        !        call dscal(m,1/h,tau(2,1),2)
!        !        call dscal(m,1/h,tau(2,1),2)
!        do  120  j = 1, m
!           !          call yyxcpy(l,-tau(1,j),-tau(2,j),ar(i,1),ai(i,1),nm,
!           !     .    zr(1,j),zi(1,j),1,.true.)
!           zr(:,j)= -tau(1,j)*ar(i,:) + tau(2,j)*ai(i,:) + zr(:,j)
!           zi(:,j)= +tau(1,j)*ai(i,:) - tau(2,j)*ar(i,:) + zi(:,j)
! 120    enddo
! 140 enddo
!   end subroutine htribk

!   subroutine htridi(nm,n,ar,ai,d,e,e2,tau)

!     !     implicit none
!     integer :: i,j,k,l,n,nm
!     double precision :: ar(nm,n),ai(nm,n),d(n),e(n),e2(n),tau(2,n)
!     double precision :: f,g,h,fi,gi,hh,si,scale!,pythag

!     !     this subroutine is a translation of a complex analogue of
!     !     the algol procedure tred1, num. math. 11, 181-195(1968)
!     !     by martin, reinsch, and wilkinson.
!     !     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).

!     !     this subroutine reduces a complex hermitian matrix
!     !     to a real symmetric tridiagonal matrix using
!     !     unitary similarity transformations.

!     !     on input

!     !        nm must be set to the row dimension of two-dimensional
!     !          array parameters as declared in the calling program
!     !          dimension statement.

!     !        n is the order of the matrix.

!     !        ar and ai contain the real and imaginary parts,
!     !          respectively, of the complex hermitian input matrix.
!     !          only the lower triangle of the matrix need be supplied.

!     !     on output

!     !        ar and ai contain information about the unitary trans-
!     !          formations used in the reduction in their full lower
!     !          triangles.  their strict upper triangles and the
!     !          diagonal of ar are unaltered.

!     !        d contains the diagonal elements of the the tridiagonal matrix.

!     !        e contains the subdiagonal elements of the tridiagonal
!     !          matrix in its last n-1 positions.  e(1) is set to zero.

!     !        e2 contains the squares of the corresponding elements of e.
!     !          e2 may coincide with e if the squares are not needed.

!     !        tau contains further information about the transformations.

!     !     calls pythag for  dsqrt(a*a + b*b) .

!     !     questions and comments should be directed to burton s. garbow,
!     !     mathematics and computer science div, argonne national laboratory

!     !     this version dated august 1983.

!     !     ------------------------------------------------------------------

!     tau(1,n) = 1.0d0
!     tau(2,n) = 0.0d0

!     do   i = 1, n
!        d(i) = ar(i,i)
!     enddo
!     do  300  i = n, 1, -1
!        l = i - 1
!        h = 0.0d0
!        scale = 0.0d0
!        if (l < 1) go to 130
!        !     .......... scale row (algol tol then not needed) ..........
!        do   k = 1, l
!           scale = scale + dabs(ar(i,k)) + dabs(ai(i,k))
!        enddo

!        if (scale /= 0.0d0) go to 140
!        tau(1,l) = 1.0d0
!        tau(2,l) = 0.0d0
! 130    continue
!        e(i) = 0.0d0
!        e2(i) = 0.0d0
!        go to 290

! 140    continue
!        do 150  k = 1, l
!           ar(i,k) = ar(i,k) / scale
!           ai(i,k) = ai(i,k) / scale
!           h = h + ar(i,k) * ar(i,k) + ai(i,k) * ai(i,k)
! 150    enddo

!        e2(i) = scale * scale * h
!        g = dsqrt(h)
!        e(i) = scale * g
!        f = pythag(ar(i,l),ai(i,l))
!        !     .......... form next diagonal element of matrix t ..........
!        if (f == 0.0d0) go to 160
!        tau(1,l) = (ai(i,l) * tau(2,i) - ar(i,l) * tau(1,i)) / f
!        si = (ar(i,l) * tau(2,i) + ai(i,l) * tau(1,i)) / f
!        h = h + f * g
!        g = 1.0d0 + g / f
!        ar(i,l) = g * ar(i,l)
!        ai(i,l) = g * ai(i,l)
!        if (l == 1) go to 270
!        go to 170
! 160    continue
!        tau(1,l) = -tau(1,i)
!        si = tau(2,i)
!        ar(i,l) = g
! 170    continue
!        f = 0.0d0
!        g = tau(1,l)
!        ! ... use tau1,1:l) and tau(2,1:l) as work array ...
!        !        do  180  j = 1, l
!        !          tau(1,j) = 0
!        !          tau(2,j) = 0
!        !     180   continue
!        tau = 0d0
!        do  220  k = 1, l
!           !          call yyxcpy(k-1,ar(i,k),-ai(i,k),ar(k,1),ai(k,1),nm,
!           !     .          tau,tau(2,1),2,.true.)
!           !          call yyaxpy(l-k+1,ar(i,k),-ai(i,k),ar(k,k),ai(k,k),1,
!           !     .         tau(1,k),tau(2,k),2,.true.)
!           tau(1,1:k-1) =  ar(i,k)*ar(k,1:k-1) + ai(i,k)*ai(k,1:k-1) + tau(1,1:k-1)
!           tau(2,1:k-1) = -ar(i,k)*ai(k,1:k-1) - ai(i,k)*ar(k,1:k-1) + tau(2,1:k-1)
!           tau(1,k:l)   =  ar(i,k)*ar(k,k:l)   + ai(i,k)*ai(k,k:l)   + tau(1,k:l)
!           tau(2,k:l)   =  ar(i,k)*ai(k,k:l)   - ai(i,k)*ar(k,k:l)   + tau(2,k:l)
! 220    enddo
!        !        call dcopy(l,tau,2,e,1)
!        e(1:l) = tau(1,1:l)
!        tau(1,l) = g
!        do  240  j = 1, l
!           e(j) = e(j)/h
!           tau(2,j) = tau(2,j)/h
!           f = f + e(j) * ar(i,j) - tau(2,j) * ai(i,j)
! 240    enddo
!        hh = f / (h + h)
!        !     .......... form reduced a ..........
!        do  261  j = 1, l
!           f = ar(i,j)
!           g = e(j) - hh * f
!           e(j) = g
!           fi = -ai(i,j)
!           gi = tau(2,j) - hh * fi
!           tau(2,j) = -gi

!           do  260  k = 1, j
!              ar(j,k) = ar(j,k) - f * e(k) - g * ar(i,k) &
!                   + fi * tau(2,k) + gi * ai(i,k)
!              ai(j,k) = ai(j,k) - f * tau(2,k) - g * ai(i,k) &
!                   - fi * e(k) - gi * ar(i,k)
! 260       enddo
! 261    enddo

! 270    continue
!        do  280  k = 1, l
!           ar(i,k) = scale * ar(i,k)
!           ai(i,k) = scale * ai(i,k)
! 280    enddo

!        tau(2,l) = -si
! 290    continue
!        hh = d(i)
!        d(i) = ar(i,i)
!        ar(i,i) = hh
!        ai(i,i) = scale * dsqrt(h)
! 300 enddo

!     return
!   end subroutine htridi
  subroutine imtql2(nm,n,d,e,z,ierr)
    integer :: i,j,k,l,m,n,ii,nm,mml,ierr
    double precision :: d(n),e(n),z(nm,n)
    double precision :: b,c,f,g,p,r,s,tst1,tst2,d1mach,tol !,pythag
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
105    do 110 m = l, n
          if (m == n) go to 120
          tst1 = dabs(d(m)) + dabs(d(m+1))
          tst2 = tst1 + dabs(e(m))
          if (dabs(tst2-tst1) < tol*dabs(tst2)) go to 120
110    enddo

120    p = d(l)
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
180       enddo

200    enddo

       d(l) = d(l) - p
       e(l) = g
       e(m) = 0.0d0
       go to 105
       !     .......... recover from underflow ..........
210    d(i+1) = d(i+1) - p
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
260    enddo

       if (k == i) go to 300
       d(k) = d(i)
       d(i) = p

       do 280 j = 1, n
          p = z(j,i)
          z(j,i) = z(j,k)
          z(j,k) = p
280    enddo

300 enddo

    go to 1001
    !     .......... set error -- no convergence to an
    !                eigenvalue after 30 iterations ..........
1000 ierr = l
1001 continue
    call tcx('imtql2')
  end subroutine imtql2


  integer function chkhss(hess,n,w5,evtol,iopt,z,e)
    !- Returns the number of eigenvalues greater than a specified value
    !o z: eigenvectors if iopt > 1; otherwise, z unchanged
    !o e: eigenvalues
    !o hess: DESTROYED on output
    !     implicit none
    integer :: n,iopt
    double precision :: hess(n,n),w5(n,5),evtol,z(n,n),e(n)
    double precision :: dum(1,1)
    integer :: i,j

    j = 0
    chkhss = j
    if (iopt == 0) return
    if (iopt > 1) j = n
    call pshpr(0)
    !     call prm(.false.,1,6,' ',w5,n,5)

    call dsev1(n,hess,dum,w5,0,.false.,.false.,0, &
         j,9d9,i,z,e)
    call poppr
    j = 0
    do  i = 1, n
       if (e(i) < evtol) j = j+1
    enddo
    if (j > 0 .AND. iopt == 1) then
       print 333, j,e(1)
333    format(/' chkhss: Hessian has',i3,' negative evals; e(1)=', &
            f8.4,'.  Not updated.')
    elseif (j > 0) then
       print 332, j,e(1)
332    format(/' chkhss: Hessian has',i3,' negative evals; e(1)=', &
            f8.4,'.  Project unwanted part.')
       !        call prmx('e',e(1),n,n,1)
       !        call prmx('z',z,n,n,n)
       call rx('not implemented')
    endif
    chkhss = j
  end function chkhss

  subroutine dsev1(n,h,o,wk,ipr,lx,lov,linv,nmx,emx,nev,z,e)
    !- Diagonalize secular equation with overlap
    !----------------------------------------------------------------------
    !i Inputs
    !i    n:    dimension of hamiltonian
    !i    h,o:  hamiltonian, overlap matrices
    !i    wk:   work array length at least (5,11)*ndim (linv = (0,>0))
    !i    ipr:  verbosity
    !i    nmx:  maximum number of eigenvectors to be found
    !i    emx:  eigenvalue limit for eigenvectors to be found
    !i    lov:  true: overlap matrix, false, no overlap matrix
    !i   linv:  1 for inverse iteration
    !i    lx:   true, call x version for overlap handling
    !o Outputs
    !o    z:    eigenvectors; e, eigenvalues
    !o    nev:  number of eigenvectors found
    !u Updates
    !u   07 Apr 07 Bug fix: returns true nev, not nmx in tql2 branch
    !r Remarks
    !r    h,o,z are dimensioned (n,n)
    !r    h,o are OVERWRITTEN in this routine
    !----------------------------------------------------------------------
    !     implicit none
    ! Passed parameters
    logical :: lov,lx
    integer :: n,ipr,nmx,nev,linv
    double precision :: h(n,n),o(n,n),z(n,n),e(n),wk(n,11),emx
    ! Local variables
    integer :: ierr,j,iprint,iwk(n)

    call tcn('dsev1')
    nev = 0

    ! --- Eigenvalues of O^-1/2  H  O^-1/2 ---
    if (lov) then
       call dschd(n,n,o,wk,lx,ierr)
       call rxx(ierr.ne.0,'DSEV1: error in dschd')
       if (lx) then
          call dsredx(n,n,h,o,z)
       else
          call dsred(n,n,h,o)
       endif
    endif

    if (linv == 1 .AND. nmx > 0) then
       call dtridx(n,n,h,wk,wk(1,4),wk(1,5),wk(1,2))
    endif

    if (nmx <= 0) then
       call dtridx(n,n,h,e,wk,wk(1,4),wk(1,2))
       do   j = 1, n
          wk(j,1) = wk(j,1)**2
       enddo
       call tqlrat(n,e,wk,ierr)
       call rxx(ierr.ne.0,'DSEV1: tqlrat cannot find all evals')
       goto 100
    else if (linv == 1) then
       call imtqlv(n,wk,wk(1,4),wk(1,5),e,iwk,ierr,wk(1,6))
       call rxx(ierr.ne.0,'DSEV1: imtqlv cannot find all evals')
       !   ... Determine number of eigenvectors to be calculated
       nev = 1
       do   j = 2, n
          if (j <= nmx .AND. e(j-1) <= emx) nev = j
       enddo
       call tinvit(n,n,wk(1,1),wk(1,4),wk(1,5),nev,e,iwk,z, & !wk(1,11),z, &
            ierr,wk(1,6),wk(1,7),wk(1,8),wk(1,9),wk(1,10))
       call rxx(ierr.ne.0,'DSEV1: tinvit cannot find all evecs')
       call dtribx(n,n,h,wk(1,2),nev,z)
    else
       call tred2(n,n,h,e,wk(1,2),z)
       call tql2(n,n,e,wk(1,2),z,ierr)
       call rxx(ierr.ne.0,'DSEV1: tql2 cannot find all evecs')
       nev = n
    endif

    ! --- Get the eigenvectors of H - E O ---
    if ( .NOT. lov .OR. nmx <= 0) goto 100
    if (lx) then
       call dcopy(n*n,z,1,h,1)
       call dmpy(o,n,1,h,n,1,z,n,1,n,nev,n)
    else
       call dsbak(n,n,o,nev,z)
    endif

    ! --- Exit ---
100 continue
    if(iprint() >= 60) print 600, e
600 format(' evl='/(1x,8f10.5))
    call tcx('dsev1')
  end subroutine dsev1



  subroutine dtribx(nm,n,ar,tau,m,zr)

    implicit none
    integer :: i,j,k,l,m,n,nm
    double precision :: ar(nm,n),tau(n,2),zr(nm,m)
    double precision :: h,s

    !     this subroutine is a translation of a complex analogue of
    !     the algol procedure trbak1, num. math. 11, 181-195(1968)
    !     by martin, reinsch, and wilkinson.
    !     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).

    !     this subroutine forms the eigenvectors of a complex hermitian
    !     matrix by back transforming those of the corresponding
    !     real symmetric tridiagonal matrix determined by  htridi.

    !     on input

    !        nm must be set to the row dimension of two-dimensional
    !          array parameters as declared in the calling program
    !          dimension statement.

    !        n is the order of the matrix.

    !        ar and ai contain information about the unitary trans-
    !          formations used in the reduction by  htridi  in their
    !          full lower triangles except for the diagonal of ar.

    !        tau contains further information about the transformations.

    !        m is the number of eigenvectors to be back transformed.

    !        zr contains the eigenvectors to be back transformed
    !          in its first m columns.

    !     on output

    !        zr and zi contain the real and imaginary parts,
    !          respectively, of the transformed eigenvectors
    !          in their first m columns.

    !     note that the last component of each returned vector
    !     is real and that vector euclidean norms are preserved.

    !     questions and comments should be directed to burton s. garbow,
    !     mathematics and computer science div, argonne national laboratory

    !     this version dated august 1983.
    !     Aug 1990 MvS altered into daxpy-style loops.  Use this
    !     version with htridx, for unit strides
    !     ------------------------------------------------------------------

    if (m == 0) return
    !     .......... transform the eigenvectors of the real symmetric
    !                tridiagonal matrix to those of the full symmetric
    !                tridiagonal matrix. ..........
    do  k = 1, n
       do j = 1, m
          zr(k,j) = zr(k,j) * tau(k,1)
       enddo
    enddo

    if (n == 1) return
    do 140 i = 2, n
       l = i - 1
       h = -tau(i,2)
       if (h == 0.0d0) go to 140
       do 130 j = 1, m
          s = 0.0d0
          do 110 k = 1, l
             s = s + ar(k,i) * zr(k,j)
110       enddo
          !     .......... double divisions avoid possible underflow ..........
          s = (s / h) / h
          do k = 1, l
             zr(k,j) = zr(k,j) - s * ar(k,i)
          enddo
130    enddo
140 enddo
  end subroutine dtribx


  subroutine dtridx(nm,n,ar,d,e,e2,tau)

    implicit none
    !     implicit real*8 (v)
    integer :: i,j,k,l,n,nm
    double precision :: ar(nm,n),d(n),e(n),e2(n),tau(n,2)
    double precision :: f,g,h,hh,scale

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

    !     this version adapted from august 1983 htridi.  Differences
    !     with htridi are that indices of tau and a are permuted (uses
    !     upper triangle of a)
    !     ------------------------------------------------------------------

    call tcn('dtridx')
    tau(n,1) = 1.0d0

    do  i = 1, n
       d(i) = ar(i,i)
    enddo
    do  300  i = n, 1, -1
       l = i - 1
       h = 0.0d0
       scale = 0.0d0
       if (l < 1) go to 130
       !     .......... scale row (algol tol then not needed) ..........
       do   k = 1, l
          scale = scale + dabs(ar(k,i))
       enddo
       if (scale /= 0.0d0) go to 140
       tau(l,1) = 1.0d0
130    e(i) = 0.0d0
       e2(i) = 0.0d0
       go to 290

140    do  150  k = 1, l
          ar(k,i) = ar(k,i) / scale
          h = h + ar(k,i)**2
150    enddo

       e2(i) = scale * scale * h
       g = dsqrt(h)
       e(i) = scale * g
       f = dabs(ar(l,i))
       !     .......... form next diagonal element of matrix t ..........
       if (f == 0.0d0) go to 160
       tau(l,1) = ( - ar(l,i) * tau(i,1)) / f
       h = h + f * g
       g = 1.0d0 + g / f
       ar(l,i) = g * ar(l,i)
       if (l == 1) go to 270
       go to 170
160    tau(l,1) = -tau(i,1)
       ar(l,i) = g
170    f = 0.0d0
       do  240  j = 1, l
          !     .......... form element of a*u ..........
          g = 0.d0
          do  180  k = 1, j
             g = g + ar(k,j) * ar(k,i)
180       enddo
          do  200  k = j+1, l
             g = g + ar(j,k) * ar(k,i)
200       enddo
          !     .......... form element of p ..........
220       e(j) = g / h
          f = f + e(j) * ar(j,i)
240    enddo

       hh = f / (h + h)
       !     .......... form reduced a ..........
       do  j = 1, l
          f = ar(j,i)
          g = e(j) - hh * f
          e(j) = g
          do   k = 1, j
             ar(k,j) = ar(k,j) - f * e(k) - g * ar(k,i)
          enddo
       enddo
270    do  280  k = 1, l
          ar(k,i) = scale * ar(k,i)
280    enddo
290    hh = d(i)
       d(i) = ar(i,i)
       tau(i,2) = -scale * dsqrt(h)
300 enddo
    call tcx('dtridx')
  end subroutine dtridx

  subroutine dsredx(nm,n,hr,ar,wkr)
    !- Reduction of nonorthogonal hermitian matrix to orthogonal form
    ! ----------------------------------------------------------------
    !i Inputs
    !i   h,nm: hermitian matrix, declared as h(nm,*).  (Lower triangle only)
    !i   a: nonorthogonality matrix, Cholesky-decomposed by dschd into L(L+)
    !i      and L inverted.
    !i   wk: work array of same dimensions as h,a
    !i   n:  order of h and a
    !i   sw: false if h and a matrix are real.
    !o Outputs
    !o   H replaced by H'' = L^-1 H (L+)^-1
    !r Remarks
    !r   This version has more floating point operations and uses more
    !r   memory than yyhred, but calls zmpy for the n^3 operations.
    ! ----------------------------------------------------------------
    !     implicit none
    ! Passed parameters
    integer :: n,nm
    double precision :: hr(nm,n),ar(nm,n),wkr(nm,n)
    ! Local parameters
    integer :: i,j,mssiz
    parameter (mssiz = 48)
    call tcn('dsredx')
    ! --- Make a strictly L^-1 ---
    do  i = 1, n
       do j = 1, i-1
          ar(j,i) = 0
       enddo
    enddo
    ! --- wk <-  L^-1 H ----
    call dmpyt(2,mssiz,n,n,ar,nm,hr,nm,wkr,nm)
    ! --- Copy L^-1 to (L+)^-1 ---
    do    i = 1, n
       do j = i+1, n
          ar(i,j) =  ar(j,i)
          ar(j,i) = 0
       enddo
    enddo
    ! --- New H <-  L^-1 H L+^-1 ----
    call dmpyt(11,mssiz,n,n,wkr,nm,ar,nm,hr,nm)
    call tcx('dsredx')
  end subroutine dsredx

  subroutine dmpyt(case,nrow,n,m,ar,lda,br,ldb,cr,ldc)
    !- Multiplication of a triangular matrix by regular matrix
    ! ----------------------------------------------------------------
    !i Inputs:
    !i   n: number of rows of destination matrix to calculate
    !i   m: number of columns of destination matrix to calculate
    !i   case:  1(ones digit), rhs upper triangular
    !i          2(ones digit), lhs lower triangular
    !i          1(tens digit), result assumed symmetric
    !o Outputs:
    !o   product matrix stored in c
    !r Remarks:
    !r   This version uses BLAS3 dgemm.
    ! ----------------------------------------------------------------
    !     implicit none
    integer :: n,m,lda,ldb,ldc,nrow,case
    double precision :: ar(lda,1),br(ldb,1),cr(ldc,1)
    integer :: ic,ncol,nc,lc,ir,lr,nr,k
    logical :: lherm

    lherm = case/10 .gt. 0
    goto (1,2), mod(case,10)
    call rx('bad case in dmpyt')

    ! --- Case rhs is upper triangular ---
1   continue
    ncol = nrow
    k = n
    do  10  ic = 1, m, ncol
       nc = min(m-ic+1,ncol)
       lc = nc+ic-1
       if (lherm) k = lc
       call dmpy(ar,lda,1,br(1,ic),ldb,1,cr(1,ic),ldc,1, k,nc,lc)
10  enddo
    ! --- Copy to upper triangle ---
    if ( .NOT. lherm) return
    do    ir = 1, n
       do    ic = ir+1, n
          cr(ic,ir) =  cr(ir,ic)
       enddo
    enddo
    return
    ! --- Case lhs is lower triangular ---
2   continue
    k = m
    do  20  ir = 1, n, nrow
       nr = min(n-ir+1,nrow)
       lr = nr+ir-1
       if (lherm) k = lr
       call dmpy(ar(ir,1),lda,1,br,ldb,1,cr(ir,1),ldc,1,nr,k,lr)
20  enddo
    ! --- Copy to lower triangle ---
    if ( .NOT. lherm) return
    do    ir = 1, n
       do  ic = ir+1, n
          cr(ir,ic) =  cr(ic,ir)
       enddo
    enddo
  end subroutine dmpyt

  subroutine dsbak(nm,n,ar,m,zr)
    !- Back-transforms eigenvectors to nonorthogonal representation
    ! ----------------------------------------------------------------
    !i Inputs
    !i   z,nm: eigenvectors, declared as z(nm,*)
    !i   n: order of a and z
    !i   a: nonorthogonality matrix, Cholesky-decomposed by dschd into L(L+)
    !i   m: number of eigenvectors to be back transformed.
    !o Outputs
    !o   z transformed eigenvectors
    !r Remarks
    !r   Nonorthogonal eigenvectors are given by z <- (L+)^-1 z
    !r   This version uses vectorizable BLAS-style daxpy loops.
    ! ----------------------------------------------------------------
    integer :: m,n,nm
    double precision :: ar(nm,n),zr(nm,m)
    integer :: nmi,k
    do  10  nmi = n, 1, -1
       do    k = n, nmi+1, -1
          call daxpy(m,-ar(k,nmi),zr(k,1),nm,zr(nmi,1),nm)
       enddo
       call dscal(m,1/ar(nmi,nmi),zr(nmi,1),nm)
10  enddo
  end subroutine dsbak

  subroutine dschd(nm,n,ar,wk,swI,ierr)
    !- Cholesky decomposition of hermitian matrix
    ! ----------------------------------------------------------------
    !i Inputs
    !i   a,nm: hermitian matrix, declared as a(nm,*).  (Lower triangle only)
    !i   n:  order of a.
    !i   swI:true if to return L^-1
    !i   wk: real work array of dimension n (uneeded if swI is false)
    !o Outputs
    !o   A replaced by L or L^-1 if swI true
    !o   ierr:nonzero if matrix not positive definite.
    !r Remarks
    !r   Makes ljj = (ajj - sum_k<j ljk (l+)jk)^1/2
    !r         lij = (aij - sum_k<j lik (l+)jk)/ljj for i>j
    !r   The strict upper triangle is unused.
    !r   This version uses BLAS-style daxpy loops with unit stride.
    ! ----------------------------------------------------------------
    !     implicit none
    ! Passed parameters
    logical :: swI
    integer :: ierr,n,nm
    double precision :: ar(nm,n),wk(n)
    ! Local parameters
    integer :: i,j,k
    double precision :: xr
    call tcn('dschd')
    ! --- Cholesky decomposition of a into L(L+) (lower triangle only) ---
    do  10  j = 1, n
       do  20  k = 1, j-1
          xr = ar(j,k)
          do    i = j, n
             ar(i,j) = ar(i,j) - xr*ar(i,k)
          enddo
20     enddo
       ierr = j
       if (ar(j,j) <= 0) return
       ar(j,j) = dsqrt(ar(j,j))
       call dscal(n-j,1/ar(j,j),ar(j+1,j),1)
10  enddo
    ierr = 0
    if ( .NOT. swI) return
    ! --- Inversion of L (lower triangle only) ---
    do    j = 1, n
       ar(j,j) = 1/ar(j,j)
    enddo
    do  40  j = 2, n
       call dcopy(n-j+1,ar(j,j-1),1,wk,1)
       call dpzero(ar(j,j-1),n-j+1)
       ! This loop parallelizable ...
       do  50  k = 1, j-1
          xr = ar(j-1,k)
          do  i = j, n
             ar(i,k) = ar(i,k) - xr*wk(i-j+1)
          enddo
          ar(j,k) = ar(j,k)*ar(j,j)
50     enddo
40  enddo
    call tcx('dschd')
  end subroutine dschd


  subroutine dmpy(a,nca,nra,b,ncb,nrb,c,ncc,nrc,n,m,l)
    !- matrix multiplication
    ! ----------------------------------------------------------------
    !i Inputs:
    !i   a,nca,nra is the left matrix and respectively the spacing
    !i      between elements in adjacent columns and rows.
    !i   b,ncb,nrb is the right matrix and respectively the spacing
    !i      between elements in adjacent columns and rows.
    !i   c,ncc,nrc is the product matrix and respectively the spacing
    !i      between elements in adjacent columns and rows.
    !i   n,m: the number of rows and columns, respectively, to calculate
    !i   l:   length of vector for matrix multiply
    !o Outputs:
    !o   product matrix stored in c
    !r Remarks:
    !r   This is a general-purpose matrix multiplication routine,
    !r   multiplying a subblock of matrix a by a subblock of matrix b.
    !r   Normally matrix nc{a,b,c} is the row dimension of matrix {a,b,c}
    !r   and nr{a,b,c} is 1.  Reverse nr and nc for a transposed matrix.
    !r   Arrays are locally one-dimensional so as to optimize inner loop,
    !r   which is executed n*m*l times.  No attempt is made to optimize
    !r   the outer loops, executed n*m times.
    !r     Examples: product of (n,l) subblock of a into (l,m) subblock of b
    !r   call dmpy(a,nrowa,1,b,nrowb,1,c,nrowc,1,n,m,l)
    !r     nrowa, nrowb, and nrowc are the leading dimensions of a, b and c.
    !r     To generate the tranpose of that product, use:
    !r   call dmpy(a,nrowa,1,b,nrowb,1,c,1,nrowc,n,m,l)
    ! ----------------------------------------------------------------
    !     implicit none
    ! Passed parameters
    integer :: nca,nra,ncb,nrb,ncc,nrc,n,m,l
    double precision :: a(0:*), b(0:*), c(0:*)
    ! Local parameters
    double precision :: ar
    integer :: i,j,k,nccj,nrci,nrcicj,ncbj
    integer :: lda,ldb
    character(1) :: transa,transb

    if (nra == 1) then
       lda = nca
       transa = 'n'
    elseif (nca == 1) then
       lda = nra
       transa = 't'
    else
       lda = -1
    endif
    if (nrb == 1) then
       ldb = ncb
       transb = 'n'
    elseif (ncb == 1) then
       ldb = nrb
       transb = 't'
    else
       ldb = -1
    endif
    if (min(lda,ldb) < 0 .OR. nrc /= 1) goto 11
    ! if PARALLEL
    !      call pp_$dgemm(transa,transb,n,m,l,1d0,a,lda,b,ldb,0d0,c,ncc)
    ! else
    call dgemm(transa,transb,n,m,l,1d0,a,lda,b,ldb,0d0,c,ncc)
    ! endif
    return
11  continue
    ! endif
    ! --- Initialize array to zero ---
    do  10  i = n-1, 0, -1
       nrci = nrc*i
       nccj = -ncc
       do   j = m-1, 0, -1
          nccj = nccj + ncc
          nrcicj = nrci + nccj
          c(nrcicj) = 0
       enddo
10  enddo

    ! --- Do multiplication ---
    do  201  k = l-1, 0, -1
       do  20  i = n-1, 0, -1
          ar = a(nra*i + nca*k)
          if (ar == 0) goto 20
          call daxpy(m,ar,b(nrb*k),ncb,c(nrc*i),ncc)
20     enddo
201 enddo
  end subroutine dmpy


  subroutine tqlrat(n,d,e2,ierr)
    integer :: i,j,l,m,n,ii,l1,mml,ierr
    double precision :: d(n),e2(n)
    double precision :: b=1d99,c=1d99,f,g=1d99,h,p,r,s,t!,pythag!,epslon
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
105    do 110 m = l, n
          if (e2(m) <= c) go to 120
          !     .......... e2(n) is always zero, so there is no exit
          !                through the bottom of the loop ..........
110    enddo

120    if (m == l) go to 210
130    if (j == 30) go to 1000
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
200    enddo

       e2(l) = s * g
       d(l) = h
       !     .......... guard against underflow in convergence test ..........
       if (h == 0.0d0) go to 210
       if (dabs(e2(l)) <= dabs(c/h)) go to 210
       e2(l) = h * e2(l)
       if (e2(l) /= 0.0d0) go to 130
210    p = d(l) + f
       !     .......... order eigenvalues ..........
       if (l == 1) go to 250
       !     .......... for i=l step -1 until 2 do -- ..........
       do 230 ii = 2, l
          i = l + 2 - ii
          if (p >= d(i-1)) go to 270
          d(i) = d(i-1)
230    enddo

250    i = 1
270    d(i) = p
290 enddo
    go to 1001
    !     .......... set error -- no convergence to an
    !                eigenvalue after 30 iterations ..........
1000 ierr = l
1001 return
  end subroutine tqlrat

  subroutine tinvit(nm,n,d,e,e2,m,w,ind,z, &
       ierr,rv1,rv2,rv3,rv4,rv6)
    integer :: i,j,m,n,p,q,r,s,ii,ip,jj,nm,its,tag,ierr,group
    double precision :: d(n),e(n),e2(n),w(m),z(nm,m), &
         rv1(n),rv2(n),rv3(n),rv4(n),rv6(n)
    double precision :: u,v,uk,xu,x0,x1,eps2,eps3,eps4,norm,order!,pythag !,epslon, &       pythag
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
490    continue
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
505    group = 0
       go to 520
       !     .......... look for close or coincident roots ..........
510    if (dabs(x1-x0) >= eps2) go to 505
       group = group + 1
       if (order * (x1 - x0) <= 0.0d0) x1 = x0 + order * eps3
       !     .......... elimination with interchanges and
       !                initialization of vector ..........
520    v = 0.0d0

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
540       xu = e(i) / u
          rv4(i) = xu
          rv1(i-1) = u
          rv2(i-1) = v
          rv3(i-1) = 0.0d0
560       u = d(i) - x1 - xu * v
          if (i /= q) v = e(i+1)
580    enddo

       if (u == 0.0d0) u = eps3
       rv1(q) = u
       rv2(q) = 0.0d0
       rv3(q) = 0.0d0
       !     .......... back substitution
       !                for i=q step -1 until p do -- ..........
600    do 620 ii = p, q
          i = p + q - ii
          rv6(i) = (rv6(i) - u * rv2(i) - v * rv3(i)) / rv1(i)
          v = u
          u = rv6(i)
620    enddo
       !     .......... orthogonalize with respect to previous
       !                members of group ..........
       if (group == 0) go to 700
       j = r

       do 680 jj = 1, group
630       j = j - 1
          if (ind(j) /= tag) go to 630
          xu = 0.0d0

          do i = p, q
             xu = xu + rv6(i) * z(i,j)
          enddo

          do i = p, q
             rv6(i) = rv6(i) - xu * z(i,j)
          enddo
680    enddo

700    norm = 0.0d0

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
740    continue
       xu = eps4 / norm

       do i = p, q
          rv6(i) = rv6(i) * xu
       enddo
       !     .......... elimination operations on next vector
       !                iterate ..........
780    continue
       do 820 i = ip, q
          u = rv6(i)
          !     .......... if rv1(i-1) .eq. e(i), a row interchange
          !                was performed earlier in the
          !                triangularization process ..........
          if (rv1(i-1) /= e(i)) go to 800
          u = rv6(i-1)
          rv6(i-1) = rv6(i)
800       continue
          rv6(i) = u - rv4(i) * rv6(i-1)
820    enddo

       its = its + 1
       go to 600
       !     .......... set error -- non-converged eigenvector ..........
830    continue
       ierr = -r
       xu = 0.0d0
       go to 870
       !     .......... normalize so that sum of squares is
       !                1 and expand to full order ..........
840    u = 0.0d0

       do i = p, q
          u = pythag(u,rv6(i))
       enddo

       xu = 1.0d0 / u

870    continue
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

  real(8) function epslon (x)
    double precision :: x
    double precision :: a,b,c,eps

    !     this program should function properly on all systems
    !     satisfying the following two assumptions,
    !        1.  the base used in representing floating point
    !            numbers is not a power of three.
    !        2.  the quantity  a  in statement 10 is represented to
    !            the accuracy used in floating point variables
    !            that are stored in memory.
    !     the statement number 10 and the go to 10 are intended to
    !     force optimizing compilers to generate code satisfying
    !     assumption 2.
    !     under these assumptions, it should be true that,
    !            a  is not exactly equal to four-thirds,
    !            b  has a zero for its last bit or digit,
    !            c  is not exactly equal to one,
    !            eps  measures the separation of 1.0 from
    !                 the next larger floating point number.
    !     the developers of eispack would appreciate being informed
    !     about any systems where these assumptions do not hold.

    !     this version dated 4/6/83.

    a = 4.0d0/3.0d0
10  b = a - 1.0d0
    c = b + b + b
    eps = dabs(c-1.0d0)
    if (eps == 0.0d0) go to 10
    epslon = eps*dabs(x)
    return
  END function epslon

  subroutine imtqlv(n,d,e,e2,w,ind,ierr,rv1)
    integer :: i,j,k,l,m,n,ii,mml,tag,ierr
    double precision :: d(n),e(n),e2(n),w(n),rv1(n)
    double precision :: b,c,f,g,p,r,s,tst1,tst2,d1mach,tol !,pythag
    integer :: ind(n)

    !     this subroutine is a variant of  imtql1  which is a translation of
    !     algol procedure imtql1, num. math. 12, 377-383(1968) by martin and
    !     wilkinson, as modified in num. math. 15, 450(1970) by dubrulle.
    !     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).

    !     this subroutine finds the eigenvalues of a symmetric tridiagonal
    !     matrix by the implicit ql method and associates with them
    !     their corresponding submatrix indices.

    !     on input

    !        n is the order of the matrix.

    !        d contains the diagonal elements of the input matrix.

    !        e contains the subdiagonal elements of the input matrix
    !          in its last n-1 positions.  e(1) is arbitrary.

    !        e2 contains the squares of the corresponding elements of e.
    !          e2(1) is arbitrary.

    !     on output

    !        d and e are unaltered.

    !        elements of e2, corresponding to elements of e regarded
    !          as negligible, have been replaced by zero causing the
    !          matrix to split into a direct sum of submatrices.
    !          e2(1) is also set to zero.

    !        w contains the eigenvalues in ascending order.  if an
    !          error exit is made, the eigenvalues are correct and
    !          ordered for indices 1,2,...ierr-1, but may not be
    !          the smallest eigenvalues.

    !        ind contains the submatrix indices associated with the
    !          corresponding eigenvalues in w -- 1 for eigenvalues
    !          belonging to the first submatrix from the top,
    !          2 for those belonging to the second submatrix, etc..

    !        ierr is set to
    !          zero       for normal return,
    !          j          if the j-th eigenvalue has not been
    !                     determined after 30 iterations.

    !        rv1 is a temporary storage array.

    !     calls pythag for  dsqrt(a*a + b*b) .

    !     questions and comments should be directed to burton s. garbow,
    !     mathematics and computer science div, argonne national laboratory

    !     this version dated august 1983.

    !     ------------------------------------------------------------------

    call tcn('imtqlv')
    ierr = 0
    k = 0
    tag = 0
    tol = 2*d1mach(3)

    do 100 i = 1, n
       w(i) = d(i)
       if (i /= 1) rv1(i-1) = e(i)
100 enddo

    e2(1) = 0.0d0
    rv1(n) = 0.0d0

    do 290 l = 1, n
       j = 0
       !     .......... look for small sub-diagonal element ..........
105    do 110 m = l, n
          if (m == n) go to 120
          tst1 = dabs(w(m)) + dabs(w(m+1))
          tst2 = tst1 + dabs(rv1(m))
          if (dabs(tst2-tst1) < tol*dabs(tst2)) go to 120
          !     .......... guard against underflowed element of e2 ..........
          if (e2(m+1) == 0.0d0) go to 125
110    enddo

120    if (m <= k) go to 130
       if (m /= n) e2(m+1) = 0.0d0
125    k = m
       tag = tag + 1
130    p = w(l)
       if (m == l) go to 215
       if (j == 30) go to 1000
       j = j + 1
       !     .......... form shift ..........
       g = (w(l+1) - p) / (2.0d0 * rv1(l))
       r = pythag(g,1.0d0)
       g = w(m) - p + rv1(l) / (g + dsign(r,g))
       s = 1.0d0
       c = 1.0d0
       p = 0.0d0
       mml = m - l
       !     .......... for i=m-1 step -1 until l do -- ..........
       do 200 ii = 1, mml
          i = m - ii
          f = s * rv1(i)
          b = c * rv1(i)
          r = pythag(f,g)
          rv1(i+1) = r
          if (r == 0.0d0) go to 210
          s = f / r
          c = g / r
          g = w(i+1) - p
          r = (w(i) - g) * s + 2.0d0 * c * b
          p = s * r
          w(i+1) = g + p
          g = c * r - b
200    enddo

       w(l) = w(l) - p
       rv1(l) = g
       rv1(m) = 0.0d0
       go to 105
       !     .......... recover from underflow ..........
210    w(i+1) = w(i+1) - p
       rv1(m) = 0.0d0
       go to 105
       !     .......... order eigenvalues ..........
215    if (l == 1) go to 250
       !     .......... for i=l step -1 until 2 do -- ..........
       do 230 ii = 2, l
          i = l + 2 - ii
          if (p >= w(i-1)) go to 270
          w(i) = w(i-1)
          ind(i) = ind(i-1)
230    enddo

250    i = 1
270    w(i) = p
       ind(i) = tag
290 enddo

    go to 1001
    !     .......... set error -- no convergence to an
    !                eigenvalue after 30 iterations ..........
1000 ierr = l
1001 continue
    call tcx('imtqlv')
  end subroutine imtqlv

  !      subroutine rs(nm,n,a,w,matz,z,fv1,fv2,ierr)
  subroutine rs(n,a,w,z,ierr)

    integer :: n,nm,ierr,matz
    !      double precision a(nm,n),w(n),z(nm,n),fv1(n),fv2(n)
    double precision :: a(n,n),w(n),z(n,n),fv1(n),fv2(n)

    !     this subroutine calls the recommended sequence of
    !     subroutines from the eigensystem subroutine package (eispack)
    !     to find the eigenvalues and eigenvectors (if desired)
    !     of a real symmetric matrix.

    !     on input

    !        nm  must be set to the row dimension of the two-dimensional
    !        array parameters as declared in the calling program
    !        dimension statement.

    !        n  is the order of the matrix  a.

    !        a  contains the real symmetric matrix.

    !        matz  is an integer variable set equal to zero if
    !        only eigenvalues are desired.  otherwise it is set to
    !        any non-zero integer for both eigenvalues and eigenvectors.

    !     on output

    !        w  contains the eigenvalues in ascending order.

    !        z  contains the eigenvectors if matz is not zero.

    !        ierr  is an integer output variable set equal to an error
    !           completion code described in the documentation for tqlrat
    !           and tql2.  the normal completion code is zero.

    !        fv1  and  fv2  are temporary storage arrays.

    !     questions and comments should be directed to burton s. garbow,
    !     mathematics and computer science div, argonne national laboratory

    !     this version dated august 1983.

    !     ------------------------------------------------------------------

    nm=n
    matz=1

    if (n <= nm) go to 10
    ierr = 10 * n
    go to 50

10  if (matz /= 0) go to 20
    !     .......... find eigenvalues only ..........
    call  tred1(nm,n,a,w,fv1,fv2)
    !  tqlrat encounters catastrophic underflow on the Vax
    !     call  tqlrat(n,w,fv2,ierr)
    call  tql1(n,w,fv1,ierr)
    go to 50
    !     .......... find both eigenvalues and eigenvectors ..........
20  call tred2(nm,n,a,w,fv1,z)
    call tql2(nm,n,w,fv1,z,ierr)
50  return
  end subroutine rs
  subroutine tql1(n,d,e,ierr)

    integer :: i,j,l,m,n,ii,l1,l2,mml,ierr
    double precision :: d(n),e(n)
    double precision :: c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2!,pythag

    !     this subroutine is a translation of the algol procedure tql1,
    !     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
    !     wilkinson.
    !     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).

    !     this subroutine finds the eigenvalues of a symmetric
    !     tridiagonal matrix by the ql method.

    !     on input

    !        n is the order of the matrix.

    !        d contains the diagonal elements of the input matrix.

    !        e contains the subdiagonal elements of the input matrix
    !          in its last n-1 positions.  e(1) is arbitrary.

    !      on output

    !        d contains the eigenvalues in ascending order.  if an
    !          error exit is made, the eigenvalues are correct and
    !          ordered for indices 1,2,...ierr-1, but may not be
    !          the smallest eigenvalues.

    !        e has been destroyed.

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

    do 100 i = 2, n
       e(i-1) = e(i)
100 enddo

    f = 0.0d0
    tst1 = 0.0d0
    e(n) = 0.0d0

    do 290 l = 1, n
       j = 0
       h = dabs(d(l)) + dabs(e(l))
       if (tst1 < h) tst1 = h
       !     .......... look for small sub-diagonal element ..........
       do 110 m = l, n
          tst2 = tst1 + dabs(e(m))
          if (tst2 == tst1) go to 120
          !     .......... e(n) is always zero, so there is no exit
          !                through the bottom of the loop ..........
110    enddo

120    if (m == l) go to 210
130    if (j == 30) go to 1000
       j = j + 1
       !     .......... form shift ..........
       l1 = l + 1
       l2 = l1 + 1
       g = d(l)
       p = (d(l1) - g) / (2.0d0 * e(l))
       r = pythag(p,1.0d0)
       d(l) = e(l) / (p + dsign(r,p))
       d(l1) = e(l) * (p + dsign(r,p))
       dl1 = d(l1)
       h = g - d(l)
       if (l2 > n) go to 145

       do 140 i = l2, n
          d(i) = d(i) - h
140    enddo

145    f = f + h
       !     .......... ql transformation ..........
       p = d(m)
       c = 1.0d0
       c2 = c
       el1 = e(l1)
       s = 0.0d0
       mml = m - l
       !     .......... for i=m-1 step -1 until l do -- ..........
       do 200 ii = 1, mml
          c3 = c2
          c2 = c
          s2 = s
          i = m - ii
          g = c * e(i)
          h = c * p
          r = pythag(p,e(i))
          e(i+1) = s * r
          s = e(i) / r
          c = p / r
          p = c * d(i) - s * g
          d(i+1) = h + s * (c * g + s * d(i))
200    enddo

       p = -s * s2 * c3 * el1 * e(l) / dl1
       e(l) = s * p
       d(l) = c * p
       tst2 = tst1 + dabs(e(l))
       if (tst2 > tst1) go to 130
210    p = d(l) + f
       !     .......... order eigenvalues ..........
       if (l == 1) go to 250
       !     .......... for i=l step -1 until 2 do -- ..........
       do 230 ii = 2, l
          i = l + 2 - ii
          if (p >= d(i-1)) go to 270
          d(i) = d(i-1)
230    enddo

250    i = 1
270    d(i) = p
290 enddo

    go to 1001
    !     .......... set error -- no convergence to an
    !                eigenvalue after 30 iterations ..........
1000 ierr = l
1001 return
  end subroutine tql1
  subroutine tred1(nm,n,a,d,e,e2)

    integer :: i,j,k,l,n,ii,nm,jp1
    double precision :: a(nm,n),d(n),e(n),e2(n)
    double precision :: f,g,h,scale

    !     this subroutine is a translation of the algol procedure tred1,
    !     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
    !     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).

    !     this subroutine reduces a real symmetric matrix
    !     to a symmetric tridiagonal matrix using
    !     orthogonal similarity transformations.

    !     on input

    !        nm must be set to the row dimension of two-dimensional
    !          array parameters as declared in the calling program
    !          dimension statement.

    !        n is the order of the matrix.

    !        a contains the real symmetric input matrix.  only the
    !          lower triangle of the matrix need be supplied.

    !     on output

    !        a contains information about the orthogonal trans-
    !          formations used in the reduction in its strict lower
    !          triangle.  the full upper triangle of a is unaltered.

    !        d contains the diagonal elements of the tridiagonal matrix.

    !        e contains the subdiagonal elements of the tridiagonal
    !          matrix in its last n-1 positions.  e(1) is set to zero.

    !        e2 contains the squares of the corresponding elements of e.
    !          e2 may coincide with e if the squares are not needed.

    !     questions and comments should be directed to burton s. garbow,
    !     mathematics and computer science div, argonne national laboratory

    !     this version dated august 1983.

    !     ------------------------------------------------------------------

    do 100 i = 1, n
       d(i) = a(n,i)
       a(n,i) = a(i,i)
100 enddo
    !     .......... for i=n step -1 until 1 do -- ..........
    do 300 ii = 1, n
       i = n + 1 - ii
       l = i - 1
       h = 0.0d0
       scale = 0.0d0
       if (l < 1) go to 130
       !     .......... scale row (algol tol then not needed) ..........
       do 120 k = 1, l
          scale = scale + dabs(d(k))
120    enddo

       if (scale /= 0.0d0) go to 140

       do 125 j = 1, l
          d(j) = a(l,j)
          a(l,j) = a(i,j)
          a(i,j) = 0.0d0
125    enddo

130    e(i) = 0.0d0
       e2(i) = 0.0d0
       go to 300

140    do 150 k = 1, l
          d(k) = d(k) / scale
          h = h + d(k) * d(k)
150    enddo

       e2(i) = scale * scale * h
       f = d(l)
       g = -dsign(dsqrt(h),f)
       e(i) = scale * g
       h = h - f * g
       d(l) = f - g
       if (l == 1) go to 285
       !     .......... form a*u ..........
       do 170 j = 1, l
          e(j) = 0.0d0
170    enddo

       do 240 j = 1, l
          f = d(j)
          g = e(j) + a(j,j) * f
          jp1 = j + 1
          if (l < jp1) go to 220

          do 200 k = jp1, l
             g = g + a(k,j) * d(k)
             e(k) = e(k) + a(k,j) * f
200       enddo

220       e(j) = g
240    enddo
       !     .......... form p ..........
       f = 0.0d0

       do 245 j = 1, l
          e(j) = e(j) / h
          f = f + e(j) * d(j)
245    enddo

       h = f / (h + h)
       !     .......... form q ..........
       do 250 j = 1, l
          e(j) = e(j) - h * d(j)
250    enddo
       !     .......... form reduced a ..........
       do 280 j = 1, l
          f = d(j)
          g = e(j)

          do 260 k = j, l
             a(k,j) = a(k,j) - f * e(k) - g * d(k)
260       enddo

280    enddo

285    do 290 j = 1, l
          f = d(j)
          d(j) = a(l,j)
          a(l,j) = a(i,j)
          a(i,j) = f * scale
290    enddo

300 enddo

    return
  end subroutine tred1

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
130    e(i) = d(l)

       do 135 j = 1, l
          d(j) = z(l,j)
          z(i,j) = 0.0d0
          z(j,i) = 0.0d0
135    enddo

       go to 290

140    do 150 k = 1, l
          d(k) = d(k) / scale
          h = h + d(k) * d(k)
150    enddo

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
200       enddo

220       e(j) = g
240    enddo
       !     .......... form p ..........
       f = 0.0d0

       do 245 j = 1, l
          e(j) = e(j) / h
          f = f + e(j) * d(j)
245    enddo

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
280    enddo

290    d(i) = h
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
360       enddo
361    enddo

380    continue
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

end module m_mathlib
