      subroutine yygemm(transa,transb,m,n,k,alpha,ar,ai,lda,br,bi,ldb,
     .  beta,cr,ci,ldc)
C- Analog of zgemm, using real arithmetic
C ----------------------------------------------------------------
Ci Inputs:
Ci   transa: specifies the form of op( A ) to be used in multiplication
Ci           'N' or 'n',  op( A ) = A.
Ci           'T' or 't',  op( A ) = A', with A' = transpose A
Ci           'C' or 'c',  op( A ) = conjg( A' ).
Ci   transb: specifies the form of op( B ) to be used in multiplication
Ci           'N' or 'n',  op( B ) = B.
Ci           'T' or 't',  op( B ) = B', with B' = transpose B
Ci           'C' or 'c',  op( B ) = conjg( B' ).
Ci   m     :number of rows of the matrices  op( A )  and  C
Ci   n     :number of columns of the matrices  op( B )  and  C
Ci   k     :The length of the inner product
Ci   alpha :scaling factor scaling multiplication; see Outputs
Ci   ar,ai :first matrix in multiplication, real and imaginary parts
Ci   lda   :leading dimension of a
Ci   br,bi :second matrix in multiplication, real and imaginary parts
Ci   ofbi  :offset to imaginary part of b
Ci   ldb   :leading dimension of b
Ci   beta  :amount of c to add into result; see Outputs
Ci   cr,ci :to be added to matrix product, real and imaginary parts
Ci          When beta is zero, cr,ci need not be set on input
Ci   ldc   :leading dimension of c
Co Outputs:
Co   c     :overwritten by alpha*op( A )*op( B ) + beta*C
C Remarks:
Cr   Unless to tell it not to, yygemm may perform the multiplication
Cr   in three real multiplications.  In doing so, it temporarily
Cr   overwrites ar,ai,br,bi with linear combinations eg ar+ai.
Cr   This means that the input arrays are NOT restored exactly intact
Cr   and some precision may be lost.
Cr
Cr   You can make yygemm do the multiplication in the normal 4 real
Cr   multiplications; this leaves untouched ar,ai,br,bi.  To accomplish
Cr   this, make transa two characters long, and put '4' in the second
Cr   character.  NB  if a and b point to the same address,
Cr   you MUST put '4' the second character of transa.
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      character*(*)      transa, transb
      integer            m, n, k, lda, ldb, ldc
      double precision   alpha, beta
      double precision   ar(lda,1), br(ldb,1), cr(ldc,1),
     .                   ai(lda,1), bi(ldb,1), ci(ldc,1)
C Local variables
      logical ls,lno3
      integer im,in,kk,mm,nn,ncut
      parameter (ncut=10)
      double precision xx,s
      character*1 trnsa,trnsb

      if (m .eq. 0 .or. n .eq. 0 .or. k .eq. 0) return
C     call tcn('yygemm')
      trnsa = transa
      trnsb = transb
C ... Check cases where multiplication cannot be done fast way
      lno3 = beta .ne. 0 .or.
     .       m .lt. ncut .or. n .lt. ncut .or. k .lt. ncut
      if (len(transa) .gt. 1) then
        lno3 = lno3 .or. transa(2:2) .eq. '4'
      endif
      s = alpha
C ... If beta nonzero, no work space available; do normal way
      if (lno3) then

C ---   Real part ---
        ls = trnsa .eq. 'C' .eqv. trnsb .eq. 'C'
        if (ls .and. beta .ne. 0) then
          do  10  in = 1, n
          do  10  im = 1, m
   10     cr(im,in) = -cr(im,in)
        endif

        call dgemm(trnsa,trnsb,m,n,k,s,ai,lda,bi,ldb,beta,cr,ldc)

        if (ls) then
          do  12  in = 1, n
          do  12  im = 1, m
   12     cr(im,in) = -cr(im,in)
        endif
        call dgemm(trnsa,trnsb,m,n,k,s,ar,lda,br,ldb,1d0,cr,ldc)

C ---   Imaginary part ---
        ls = trnsb .eq. 'C'
        if (ls .and. beta .ne. 0d0) then
          do  14  in = 1, n
          do  14  im = 1, m
   14     ci(im,in) = -ci(im,in)
        endif

        call dgemm(trnsa,trnsb,m,n,k,s,ar,lda,bi,ldb,beta,ci,ldc)

        if (ls .neqv. trnsa .eq. 'C') then
          do  16  in = 1, n
          do  16  im = 1, m
   16     ci(im,in) = -ci(im,in)
        endif

        call dgemm(trnsa,trnsb,m,n,k,s,ai,lda,br,ldb,1d0,ci,ldc)

        if (trnsa .eq. 'C') then
          do  18  in = 1, n
          do  18  im = 1, m
   18     ci(im,in) = -ci(im,in)
        endif

C       call tcx('yygemm')
        return
      endif

C --- Matrix multiply in three real multiplications ---
C ... ci = ai *  bi ...
      call dgemm(trnsa,trnsb,m,n,k,s,ai,lda,bi,ldb,0d0,ci,ldc)

C ... cr = ar *  br ...
      call dgemm(trnsa,trnsb,m,n,k,s,ar,lda,br,ldb,0d0,cr,ldc)

C ... cr <- ar*br - ai*bi,  ci <- -ar*br - ai*bi  ...
      if (trnsa .eq. 'C' .eqv. trnsb .eq. 'C') then
        do  20  in = 1, n
        do  20  im = 1, m
          xx = cr(im,in)
          cr(im,in) =  xx - ci(im,in)
          ci(im,in) = -xx - ci(im,in)
   20   continue
      else
        do  120  in = 1, n
        do  120  im = 1, m
          xx = cr(im,in)
          cr(im,in) =  xx + ci(im,in)
          ci(im,in) = -xx + ci(im,in)
  120   continue
      endif

C ... ci += (ar+ai)(br+bi) ...
      kk = k
      mm = m
      if (trnsa .eq. 'C' .or. trnsa .eq. 'T') then
        kk = m
        mm = k
      endif
      if (trnsa .ne. 'C') then
        do  22  in = 1, kk
        do  22  im = 1, mm
   22   ar(im,in) = ar(im,in) + ai(im,in)
      else
        do  122  in = 1, kk
        do  122  im = 1, mm
  122   ar(im,in) = ar(im,in) - ai(im,in)
      endif

      kk = k
      nn = n
      if (trnsb .eq. 'C' .or. trnsb .eq. 'T') then
        kk = n
        nn = k
      endif
      if (trnsb .ne. 'C') then
        do  24  in = 1, nn
        do  24  im = 1, kk
   24   br(im,in) = br(im,in) + bi(im,in)
      else
        do  124  in = 1, nn
        do  124  im = 1, kk
  124   br(im,in) = br(im,in) - bi(im,in)
      endif
      call dgemm(trnsa,trnsb,m,n,k,s,ar,lda,br,ldb,1d0,ci,ldc)
      kk = k
      mm = m
      if (trnsa .eq. 'C' .or. trnsa .eq. 'T') then
        kk = m
        mm = k
      endif
      if (trnsa .ne. 'C') then
        do  26  in = 1, kk
        do  26  im = 1, mm
   26   ar(im,in) = ar(im,in) - ai(im,in)
      else
        do  126  in = 1, kk
        do  126  im = 1, mm
  126   ar(im,in) = ar(im,in) + ai(im,in)
      endif
      kk = k
      nn = n
      if (trnsb .eq. 'C' .or. trnsb .eq. 'T') then
        kk = n
        nn = k
      endif
      if (trnsb .ne. 'C') then
        do  28  in = 1, nn
        do  28  im = 1, kk
   28   br(im,in) = br(im,in) - bi(im,in)
      else
        do  128  in = 1, nn
        do  128  im = 1, kk
  128   br(im,in) = br(im,in) + bi(im,in)
      endif

C     call tcx('yygemm')
      end

