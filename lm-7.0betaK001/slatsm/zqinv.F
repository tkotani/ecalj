      subroutine zqinv(cs,a,lda,ncut,n,w,ldw,ierr)
C- Inversion of a complex matrix using Strassen's algorithm
C ----------------------------------------------------------------
Ci Inputs:
Ci   cs:   :if 'h', a is assumed hermitian
Ci   a     :matrix to be inverted
Ci   lda   :leading dimension of a
Ci   n     :rank of the matrix to be inverted
Ci   ncut  :If positive, call yyqinv for n>=ncut
Ci         :If zero, call yyqinv for n>=lncut, lncut made internally
Ci         :If negative, always call lapack routines.  In which case,
Ci         :-ncut = size of work array
Ci   w     :If ncut >=0, d.p. work array of dimension ldw*(n+1)
Ci         :If ncut <0,  d.p. work array of dimension at least 3*n,
Ci         :but preferably at least 66*n
Ci   ldw   :leading dimension of w, at least as large as n
Co Outputs:
Co   a     :is overwritten by inverse of input a
Co   ierr  :returned nonzero if matrix was not fully inverted.
Cr Remarks
Cr   zqinv is a front end for either LAPACK matrix inverter, or
Cr   yyqinv, which see.
Cb Bugs
Cb   The algorithm fails if a22 is singular, even if a is not.
Cb   Similarly, if smaller subblocks are singular, yqinv may fail
Cb   when called recursively.
Cr Remarks:
Cr   This is a front end for yyqinv, which which is passed the
Cr   real and imaginary parts of arrays.
Cr
Cr   The Lapack branch is recommended: it appears to be faster
Cr   and it is more accurate.
C ----------------------------------------------------------------
C     implicit none
C ... Passed parameters
      character*(*) cs
      integer n,ncut,lda,ldw,ierr
      double precision a(2,lda,n),w(ldw,1)
C ... Local parameters
      integer lncut,lwrk,nlev,i,j

C     call zprm('zqinv',2,a,lda,n,n)

      lncut = ncut
      if (ncut .eq. 0) lncut = 64

      if (n .gt. lncut .and. lncut .gt. 0) then

        nlev = 2
        call ztoy(a,lda,n,n,0)
        call yyqinv(cs,a,a(1+lda,1,1),2*lda,nlev,n,w,ldw,ierr)
        if (ierr .ne. 0) return
        call ztoy(a,lda,n,n,1)

      else

        if (ncut .lt. 0) then
          lwrk = (-ncut - n)/2
        else
          lwrk = ldw*n/2
        endif
        if (n .gt. 1 .and. lwrk .lt. n) call rx('zqinv: work too small')
        if (cs .eq. 'h') then
          call zhetrf('L',n,a,lda,w,w(1+n,1),lwrk,ierr)
          if (ierr .ne. 0) call rx('zqinv: singular matrix')
          call zhetri('L',n,a,lda,w,w(1+n,1),ierr)
          do  i = 1, n
          do  j = i+1, n
            a(1,i,j) = a(1,j,i)
            a(2,i,j) =-a(2,j,i)
          enddo
          enddo
        else
          call zgetrf(n,n,a,lda,w,ierr)
          if (ierr .ne. 0) call rx('zqinv: singular matrix')
          call zgetri(n,a,lda,w,w(1+n,1),lwrk,ierr)

        endif

      endif

C     call zprm('zqinv after inverse',2,a,lda,n,n)

      end

