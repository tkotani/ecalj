      subroutine dsredx(nm,n,hr,ar,wkr)
C- Reduction of nonorthogonal hermitian matrix to orthogonal form
C ----------------------------------------------------------------
Ci Inputs
Ci   h,nm: hermitian matrix, declared as h(nm,*).  (Lower triangle only)
Ci   a: nonorthogonality matrix, Cholesky-decomposed by dschd into L(L+)
Ci      and L inverted.
Ci   wk: work array of same dimensions as h,a
Ci   n:  order of h and a
Ci   sw: false if h and a matrix are real.
Co Outputs
Co   H replaced by H'' = L^-1 H (L+)^-1
Cr Remarks
Cr   This version has more floating point operations and uses more
Cr   memory than yyhred, but calls zmpy for the n^3 operations.
C ----------------------------------------------------------------
C     implicit none
C Passed parameters 
      integer n,nm
      double precision hr(nm,n),ar(nm,n),wkr(nm,n)
C Local parameters 
      integer i,j,mssiz
      parameter (mssiz = 48)

      call tcn('dsredx')

C --- Make a strictly L^-1 ---
      do  i = 1, n
        do j = 1, i-1
           ar(j,i) = 0
        enddo
      enddo

C --- wk <-  L^-1 H ---- 
      call dmpyt(2,mssiz,n,n,ar,nm,hr,nm,wkr,nm)

C --- Copy L^-1 to (L+)^-1 ---
      do    i = 1, n
        do j = i+1, n
          ar(i,j) =  ar(j,i)
          ar(j,i) = 0
       enddo
      enddo

C --- New H <-  L^-1 H L+^-1 ---- 
      call dmpyt(11,mssiz,n,n,wkr,nm,ar,nm,hr,nm)

c      print 337, ((hr(i,j), j=1,n), i=1,n)
c      pause
C  337 format(4f10.6)

      call tcx('dsredx')
      end

