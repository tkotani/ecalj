      subroutine eigch(a,n,jobn,d,z,iz,wk,ier)
C- Eigenvalues of Hermitian matrix -- partial simulation of IMSL
C ----------------------------------------------------------------
Ci Inputs:
Ci   a,n:  hermitian matrix, declared as a(n,*)
Ci   jobn: 10:  eigenvalues only
Ci         11:  eigenvalues and eigenvectors
Ci              Note: other IMSL switches aren't supported here
Ci   wk:   work array of length at least 3n (4n for APOLLO_BUG)
Ci   iz:   eigenvector matrix declared as z(iz,*)
Co Outputs:
Co   d:    eigenvalues
Co   z:    eigenvector matrix
Co   ier:  0, if all eigenvectors found
Co         k  if the eigenvalue iteration fails to converge; but
Co         eigenvalues (and eigenvectors) 1 through k-1 are correct.
Co         Note: this convention differs from IMSL
Cr Remarks:
Cr   Adapted from licepack:
Cr     david kahaner, cleve moler, g. w. stewart,
Cr       n.b.s.         u.n.m.      n.b.s./u.md.
Cp Procedures used:
Cp   (eispack) htribk, htridi, imtql2, tqlrat
Cp    (blas)   dcopy, dcopym
C ----------------------------------------------------------------
      integer i,ier,j,jobn,k,l,n,iz,m,mdim,min0,ijob
      double precision a(1),d(1),wk(1),z(2)
      
C --- job numbers to match IMSL ---
      ijob = mod(jobn,10)
      if (jobn-ijob .ne. 10 .or. ijob .gt. 1) stop 'eigch: bad jobn'

      if (n .eq. 1 .and. ijob .eq. 0) goto 35
      mdim = 2 * n
      if (ijob .eq. 0) goto 5
      if (n .gt. iz) stop 'eigch: ijob .ne. 0  and  n .gt. iz'
      if (n .eq. 1) goto 35

C --- Rearrange a if necessary when n .gt. iz and ijob .ne. 0 ---
      mdim = min0(mdim,2*iz)
      if (n .le. iz) goto 5
      l = n - 1
      do  4  j = 1, l
        m = 1+j*2*iz
        k = 1+j*2*n
        call dcopy(2*n,a(k),1,a(m),1)
    4 continue

C --- Separate real and imaginary parts ---
    5 continue
      do  10  j = 1, n
        k = (j-1) * mdim + 1
        l = k + n
        call dcopy(n,a(k+1),2,wk,1)
        call dcopy(n,a(k),2,a(k),1)
        call dcopy(n,wk,1,a(l),1)
   10 continue

C --- Eigenvalues and/or eigenvectors ---
#ifndef APOLLO_BUG
      call htridi(mdim,n,a(1),a(n+1),d,wk,wk(3*n+1),wk(n+1))
#else
      call htridi(mdim,n,a(1),a(n+1),d,wk,wk,wk(n+1))
#endif
      if (ijob .eq. 0) then
        do  12  j = 1, n
   12   wk(j) = wk(j)**2
        call tqlrat(n,d,wk,ier)
        return
      endif
      do  17  j = 1, n
        k = (j-1) * mdim + 1
        m = k + n - 1
        do  16  i = k, m
          z(i) = 0d0
   16   continue
        i = k + j - 1
        z(i) = 1d0
   17 continue
      call imtql2(mdim,n,d,wk,z,ier)
      if (ier .ne. 0) return
      call htribk(mdim,n,a(1),a(n+1),wk(n+1),n,z(1),z(n+1))

C --- Convert eigenvectors to complex storage ---
      do  20  j = 1, n
        k = (j-1) * mdim + 1
        i = (j-1) * 2 * iz + 1
        l = k + n
        call dcopy(n,z(k),1,wk,1)
        call dcopy(n,z(l),1,z(i+1),2)
        call dcopy(n,wk,1,z(i),2)
   20 continue
      return

C --- Take care of n=1 case ---
   35 continue
      d(1) = a(1)
      ier = 0
      if (ijob .eq. 0) return
      z(1) = a(1)
      z(2) = 0
      end
      subroutine dcopym(n,dx,incx,dy,incy)
c
c     copies negative of vector, x, to a vector, y.  adapted from:
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dy(1)
      integer i,incx,incy,ix,iy,n
c
      ix = 1
      iy = 1
      if (incx .lt. 0) ix = (1-n)*incx + 1
      if (incy .lt. 0) iy = (1-n)*incy + 1
      do  10  i = 1, n
        dy(iy) = -dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      end

