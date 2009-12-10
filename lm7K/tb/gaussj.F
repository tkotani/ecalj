      subroutine gaussj(n,nmax,m,mmax,ipiv,indxr,indxc,a,b)
C- Linear equation solution by Gauss-Jordan elimination
C     implicit none
C Passed parameters
      integer n,nmax,m,mmax
      integer ipiv(nmax),indxr(nmax),indxc(nmax)
      double precision a(nmax,nmax),b(nmax,mmax)
C Local parameters
      integer i,j,k,l,ll,irow,icol
      double precision big,dum,pivinv

      do  10  j = 1, n
        ipiv(j) = 0
   10 continue

      do  110  i = 1, n
        big = 0d0
        do  30  j = 1, n
          if (ipiv(j) .ne. 1) then
            do  20  k = 1, n
              if (ipiv(k) .eq. 0) then
                if (dabs(a(j,k)) .ge. big) then
                  big = dabs(a(j,k))
                  irow = j
                  icol = k
                endif
              else if (ipiv(k) .gt. 1) then
                call rx('GAUSSJ: singular matrix')
              endif
   20       continue
          endif
   30   continue
        ipiv(icol) = ipiv(icol) + 1

        if (irow .ne. icol) then
          do  40  l = 1, n
            dum = a(irow,l)
            a(irow,l) = a(icol,l)
            a(icol,l) = dum
   40     continue
          do  50  l = 1, m
            dum = b(irow,l)
            b(irow,l) = b(icol,l)
            b(icol,l) = dum
   50     continue
        endif

        indxr(i) = irow
        indxc(i) = icol
        call rxx(a(icol,icol) .eq. 0d0,'GAUSSJ: singular matrix')
        pivinv = 1d0 / a(icol,icol)
        a(icol,icol) = 1d0
        do  60  l = 1, n
          a(icol,l) = a(icol,l)*pivinv
   60   continue
        do  70  l = 1, m
          b(icol,l) = b(icol,l)*pivinv
   70   continue

        do  100  ll = 1, n
          if (ll .ne. icol) then
            dum = a(ll,icol)
            a(ll,icol) = 0d0
            do  80  l = 1, n
              a(ll,l) = a(ll,l) - a(icol,l)*dum
   80       continue
            do  90  l = 1, m
              b(ll,l) = b(ll,l) - b(icol,l)*dum
   90       continue
          endif
  100   continue
  110 continue

      do  130  l = n, 1, -1
        if (indxr(l) .ne. indxc(l)) then
          do  120  k = 1, n
            dum = a(k,indxr(l))
            a(k,indxr(l)) = a(k,indxc(l))
            a(k,indxc(l)) = dum
  120     continue
        endif
  130 continue

      end
