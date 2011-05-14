      subroutine yysbmm(mode,i1,i2,j1,j2,nc,ar,ai,lda,ija,offs,
     .xr,xi,ldx,br,bi,ldb)
C- Complex sparse-block-matrix dense-matrix multiply
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci          0  multiply contents into b
Ci          1  add product  a x to existing contents of b
Ci          2  add product -a x to existing contents of b
Ci          10s digit
Ci          0  a is real
Ci          1  only the diagonal parts of a are complex
Ci          2  a is complex
Ci   i1,i2 :calculate a(i,j)*x(j) for row subblocks i = i1..i2
Ci   j1,j2 :calculate a(i,j)*x(j) for col subblocks j = j1..j2
Ci   nc    :number of columns in x and the result matrix b
Ci   a     :sparse matrix, stored in block form by rows.
Ci          a consists of a vector of matrix subblocks:
Ci          a(*,*,i) = matrix subblock i
Ci   lda   :leading dimension of a
Ci   ija   :column index packed array pointer data to array a
Ci         ija(1,*) follows essentially the same conventions
Ci         as for scalar packed arrays (see da2spr)
Ci         except that indices now refer to matrix subblocks.
Ci         ija(1,1)= n+2, where n = max(number of rows, number of cols)
Ci         ija(1,i), i = 1,..n+1 = points to first entry in a for row i
Ci         ija(1,i), i = n+2... column index element a(i).  Thus
Co                   for row i, k ranges from ija(i) <= k < ija(i+1) and
Co                   sum_j a_ij x_j -> sum_k a_(ija(2,k)) x_(ija(1,k))
Ci         ija(2,*)  pointers to the matrix subblocks blocks in a:
Ci         ija(2,i), i=1..n  pointers to blocks on the diagonal of a
Ci         ija(2,i), i=n+2.. pointers to elements of a, grouped by rows
Ci   offs  :offsets to first entries in matrix subblocks
Ci          Thus the dimension of row i = offs(i+1) - offs(i)
Ci          If a consists of scalar subblocks, offs(i) = i-1.
Ci   xr,xi :dense matrix, and second operand
Ci   ldx   :leading dimension of x
Co Outputs
Co   br,bi :result matrix
Co   ldb   :leading dimension of b
Cr Remarks
Cr   This routine multiplies a sparse matrix whose elements
Cr   are matrix subblocks, by a dense matrix.
Cu Updates
Cb Bugs
Cb   Never checked for mode=2
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer i1,i2,j1,j2,nc,lda,ldb,ldx,ija(2,*),offs(i2),mode
      double precision ar(lda,lda,1),ai(lda,lda,1),
     .xr(ldx,nc),xi(ldx,nc),br(ldb,nc),bi(ldb,nc)
C ... Local parameters
      integer ir,k,pa,ofxb,ofbb,nra,nca,ic,i,j,mode0,mode1
      double precision alp
C     integer ofx,ofb
      logical:: isanrg, l_dummy_isanrg

C --- Setup ---
      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
Ckino isanrg is logical function,       call isanrg(mode0,0,2,'yysbmm:','1s digit mode',.true.)
      l_dummy_isanrg=isanrg(mode0,0,2,'yysbmm:','1s digit mode',.true.)
Ckino isanrg is logical function,       call isanrg(mode1,0,2,'yysbmm:','10s digit mode',.true.)
      l_dummy_isanrg=isanrg(mode1,0,2,'yysbmm:','10s digit mode',.true.)
      if (mode1 .eq. 0) then
        call dspbmm(mode0,i1,i2,j1,j2,nc,ar,lda,ija,offs,xr,ldx,br,ldb)
        return
      endif
      alp = 1
      if (mode0 .eq. 2) alp = -1
C     offsets shifting origin of x and b.
C     ofx = 0
C     ofb = 0

      call tcn('yysbmm')

C --- Initialize contents of b ---
      if (mode0 .eq. 0) then
C       do  6  i = offs(i1)-ofb+1, offs(i2+1)-ofb
        do  6  j = 1, nc
          do  6  i = offs(i1)+1, offs(i2+1)
            br(i,j) = 0
            bi(i,j) = 0
    6   continue
      endif

C --- For each row ir, multiply a(ir,ic) x(ic) ---
      do  10  ir = i1, i2
C       offset to b for this subblock
C       ofbb = offs(ir) - ofb
        ofbb = offs(ir)
C       row dimension of this subblock
        nra = offs(ir+1) - offs(ir)
C       pointer to diagonal subblock in a
        pa  = ija(2,ir)

C   ... b(ir) += a(ir)*x(ir).  Skip if missing diagonal element
        if (pa .ne. 0 .and. ir .ge. j1 .and. ir .le. j2) then
C         offset to x for this subblock
C         ofxb = offs(ir) - ofx
          ofxb = offs(ir)
          call yygemm('N4','N4',nra,nc,nra,alp,ar(1,1,pa),ai(1,1,pa),
     .    lda,xr(1+ofxb,1),xi(1+ofxb,1),ldx,1d0,
     .    br(ofbb+1,1),bi(ofbb+1,1),ldb)
        endif

C  ...  b(ir) = b(ir) + a(ir,ija(k))*x(ija(k)) for nonzero blocks in a
        do  11  k = ija(1,ir), ija(1,ir+1)-1
C         column index to a and row index to x
          ic  = ija(1,k)
          if (ic .lt. j1 .or. ic .gt. j2) goto 11
C         offset to row x for this subblock
C         ofxb = offs(ic) - ofx
          ofxb = offs(ic)
C         col dimension of subblock a and row dimension of x
          nca = offs(ic+1) - offs(ic)
C         pointer to subblock in a
          pa = ija(2,k)
C         b(ir) = b(ir) + a(ir,ija(k))*x(ija(k))
          if (mode1 .eq. 1) then
            call dgemm('N','N',nra,nc,nca,alp,ar(1,1,pa),lda,
     .      xr(1+ofxb,1),ldx,1d0,br(1+ofbb,1),ldb)
            call dgemm('N','N',nra,nc,nca,alp,ar(1,1,pa),lda,
     .      xi(1+ofxb,1),ldx,1d0,bi(1+ofbb,1),ldb)
          else
            call yygemm('N4','N4',nra,nc,nca,alp,ar(1,1,pa),ai(1,1,pa),
     .      lda,xr(1+ofxb,1),xi(1+ofxb,1),ldx,1d0,
     .      br(1+ofbb,1),bi(1+ofbb,1),ldb)
          endif
   11   continue
   10 continue
      call tcx('yysbmm')

      end
      subroutine yymsbm(mode,i1,i2,j1,j2,nr,ar,ai,lda,ija,offs,
     .xr,xi,ldx,br,bi,ldb)
C- Complex dense-matrix sparse-block-matrix multiply
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit
Ci          0  multiply contents into b
Ci          1  add product x a to existing contents of b
Ci          2  add product -x a to existing contents of b
Ci          10s digit
Ci          0  a is real
Ci          1  only the diagonal parts of a are complex
Ci          2  a is complex
Ci   i1,i2 :calculate x(i)*a(i,j) for row subblocks i = i1..i2
Ci   j1,j2 :calculate x(i)*a(i,j) for col subblocks j = j1..j2
Ci   nr    :number of rows in x and the result matrix b
Ci   a     :sparse matrix, stored in block form by rows.
Ci          a consists of a vector of matrix subblocks:
Ci          a(*,*,i) = matrix subblock i
Ci   lda   :leading dimension of a
Ci   ija   :column index packed array pointer data to array a
Ci         ija(1,*) follows essentially the same conventions
Ci         as for scalar packed arrays (see da2spr)
Ci         except that indices now refer to matrix subblocks.
Ci         ija(1,1)= n+2, where n = max(number of rows, number of cols)
Ci         ija(1,i), i = 1,..n+1 = points to first entry in a for row i
Ci         ija(1,i), i = n+2... column index element a(i).  Thus
Co                   for row i, k ranges from ija(i) <= k < ija(i+1) and
Co                   sum_j a_ij x_j -> sum_k a_(ija(2,k)) x_(ija(1,k))
Ci         ija(2,*)  pointers to the matrix subblocks blocks in a:
Ci         ija(2,i), i=1..n  pointers to blocks on the diagonal of a
Ci         ija(2,i), i=n+2.. pointers to elements of a, grouped by rows
Ci   offs  :offsets to first entries in matrix subblocks
Ci          offs(i,i=1..n) offset to first row in x and b for subblock i
Ci          Thus the dimension of row i = offs(i+1) - offs(i)
Ci          If a consists of scalar subblocks, offs(i) = i-1.
Ci   x     :dense matrix, and first operand
Ci   ldx   :leading dimension of x
Co Outputs
Co   b     :result matrix
Co   ldb   :leading dimension of b
Cr Remarks
Cr   This routine multiplies x a, with a=sparse matrix whose elements
Cr   are matrix subblocks
Cb Bugs
Cb   Never checked for mode=2
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer i1,i2,j1,j2,nr,lda,ldb,ldx,ija(2,*),offs(i2),mode
      double precision xr(ldx,nr),xi(ldx,nr),br(ldb,nr),bi(ldb,nr)
      double precision ar(lda,lda,1),ai(lda,lda,1)
C ... Local parameters
      integer ir,k,pa,ofxb,ofbb,nra,nca,ic,i,j,mode0,mode1,k1,k2
      double precision alp
C     integer ofx,ofb
      logical:: isanrg, l_dummy_isanrg

C --- Setup ---
      mode0 = mod(mode,10)
      mode1 = mod(mode/10,10)
Ckino isanrg is logical function,       call isanrg(mode0,0,2,'yymsbm:','1s digit mode',.true.)
      l_dummy_isanrg=isanrg(mode0,0,2,'yymsbm:','1s digit mode',.true.)
Ckino isanrg is logical function,       call isanrg(mode1,0,2,'yymsbm:','10s digit mode',.true.)
      l_dummy_isanrg=isanrg(mode1,0,2,'yymsbm:','10s digit mode',.true.)
      if (mode1 .eq. 0) then
        call dmspbm(mode0,i1,i2,j1,j2,nr,ar,lda,ija,offs,xr,ldx,br,ldb)
        return
      endif
      alp = 1
      if (mode0 .eq. 2) alp = -1
C     offsets shifting origin of x and b.
C     ofx = 0
C     ofb = 0


C ... Initialize contents of b
      if (mode0 .eq. 0) then
C       do  6  j = offs(j1)-ofb+1, offs(j2+1)-ofb
        do  6  j = offs(j1)+1, offs(j2+1)
          do  6  i = 1, nr
            br(i,j) = 0
            bi(i,j) = 0
    6   continue
      endif

C ... Product of diagonal elements x(ir)*a(ir). Skip if missing.
      k1 = max(i1,j1)
      k2 = min(i2,j2)
      do  8  ir = k1, k2
C       Pointer to diagonal subblock in a
        pa  = ija(2,ir)
C       offset to x and b for this subblock
C       ofxb = offs(ir) - ofx
C       ofbb = offs(ir) - ofb
        ofxb = offs(ir)
        ofbb = offs(ir)
C       row and column dimension of this subblock
        nra = offs(ir+1) - offs(ir)
C   ... b(ir) <- x(ir)*a(ir).
        if (pa .ne. 0) then
          call yygemm('N4','N4',nr,nra,nra,alp,xr(1,ofxb+1),xi(1,ofxb+1)
     .    ,ldx,ar(1,1,pa),ai(1,1,pa),lda,1d0,br(1,ofbb+1),bi(1,ofbb+1)
     .    ,ldb)
        endif
    8 continue

C ... b(ija(k)) += x(ir) * a(ir,ija(k)) for nonzero blocks ija
      do  10  ir = i1, i2
C       Pointer to diagonal subblock in a
        pa  = ija(2,ir)
C       offset to x for this subblock
C       ofxb = offs(ir) - ofx
        ofxb = offs(ir)
C       col dimension of x and row dimension of a in this subblock
        nra = offs(ir+1) - offs(ir)
        do  11  k = ija(1,ir), ija(1,ir+1)-1
C         column index to a and row index to x
          ic  = ija(1,k)
          if (ic .lt. j1 .or. ic .gt. j2) goto 11
C         offset to row x for this subblock
C         ofbb = offs(ic) - ofb
          ofbb = offs(ic)
C         row dimension of subblocks a and col dimension of b
          nca = offs(ic+1) - offs(ic)
C         pointer to subblock in a
          pa = ija(2,k)
C         b(ija(k)) += x(ir) * a(ir,ija(k))
          if (mode1 .eq. 1) then
            call dgemm('N','N',nr,nca,nra,alp,xr(1,1+ofxb),ldx,
     .      ar(1,1,pa),lda,1d0,br(1,1+ofbb),ldb)
            call dgemm('N','N',nr,nca,nra,alp,xi(1,1+ofxb),ldx,
     .      ar(1,1,pa),lda,1d0,bi(1,1+ofbb),ldb)
          else
            call yygemm('N','N',nr,nca,nra,alp,xr(1,1+ofxb),xi(1,1+ofxb)
     .      ,ldx,ar(1,1,pa),ai(1,1,pa),lda,1d0,br(1,1+ofbb),bi(1,1+
     .      ofbb),ldb)
          endif
   11   continue
   10 continue
      end

