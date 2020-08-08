      subroutine dmpyt(case,nrow,n,m,ar,lda,br,ldb,cr,ldc)
C- Multiplication of a triangular matrix by regular matrix
C ----------------------------------------------------------------
Ci Inputs:
Ci   n: number of rows of destination matrix to calculate
Ci   m: number of columns of destination matrix to calculate
Ci   case:  1(ones digit), rhs upper triangular
Ci          2(ones digit), lhs lower triangular
Ci          1(tens digit), result assumed symmetric
Co Outputs:
Co   product matrix stored in c
Cr Remarks:
Cr   This version uses BLAS3 dgemm.
C ----------------------------------------------------------------
C     implicit none
      integer n,m,lda,ldb,ldc,nrow,case
      double precision ar(lda,1),br(ldb,1),cr(ldc,1)
      integer ic,ncol,nc,lc,ir,lr,nr,k
      logical lherm

      lherm = case/10 .gt. 0
      goto (1,2), mod(case,10)
      call rx('bad case in dmpyt')

C --- Case rhs is upper triangular ---
    1 continue
      ncol = nrow
      k = n
      do  10  ic = 1, m, ncol
        nc = min(m-ic+1,ncol)
        lc = nc+ic-1
        if (lherm) k = lc
        call dmpy(ar,lda,1,br(1,ic),ldb,1,cr(1,ic),ldc,1,k,nc,lc)

   10 continue
C --- Copy to upper triangle ---
      if (.not. lherm) return
      do  12  ir = 1, n
        do  12  ic = ir+1, n
   12 cr(ic,ir) =  cr(ir,ic)
      return

C --- Case lhs is lower triangular ---
    2 continue
      k = m
      do  20  ir = 1, n, nrow
        nr = min(n-ir+1,nrow)
        lr = nr+ir-1
        if (lherm) k = lr
        call dmpy(ar(ir,1),lda,1,br,ldb,1,cr(ir,1),ldc,1,nr,k,lr)
   20 continue

C --- Copy to lower triangle ---
      if (.not. lherm) return
      do  22  ir = 1, n
        do  22  ic = ir+1, n
   22 cr(ir,ic) =  cr(ic,ir)
      end

