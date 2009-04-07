      subroutine decmp2(a11,a12,a21,a22,n,kpvt,wk)
C- Decomposes symmetric matrix in the special case a12, a21 are diagonal
C ----------------------------------------------------------------
Ci Inputs
Ci   a11   :11 block of matrix to be decomposed
Ci   a12   :12 block of matrix to be decomposed (diagonal)
Ci   a21   :21 block of matrix to be decomposed (diagonal)
Ci   a22   :22 block of matrix to be decomposed
Ci   n     :dimension of a11,a12,a21,a22
Ci   wk:   work array of length n*64
Co Outputs
Co   kpvt  :an integer vector of pivot indices
Co   a11   :the inverse of a11 is returned
Co   a22   :returns factorized form of (a22 - a21 a11^-1 a12)
Cr Remarks
Cr   a11, a22 are symmetric; only the upper triangle is used.
Cr   a12 and a21 are diagonal and stored as vectors
Cr   It is assumed that the a21(i)/a12(i) is constant
Cr   Use with dsolv2 to solve simultaneous equations.
Cu Updates
Cu   06 Aug 06 Redesigned to work with lapack
C ----------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer n,kpvt(n)
      double precision a11,a12,a21,a22,wk
      dimension a11(n,n),a12(1),a21(1),a22(n,n),wk(n)
C ... Local parameters
      integer ndef,i,j

      call tcn('decmp2')

C      call dpack(a11,n)
C      call dpack(a22,n)

C --- Make a11^-1 ---
C ... Packed Linpack version
C     call dspfa(a11,n,kpvt,ndef)
C     if (ndef .ne. 0) call rx('decmp2: matrix a11 singular')
C     call dspdi(a11,n,kpvt,0d0,0d0,wk,1)
C ... Linpack version
C     call dsifa(a11,n,n,kpvt,ndef)
C     if (ndef .ne. 0) call rx('decmp2: matrix a11 singular')
C     call dsidi(a11,n,n,kpvt,wk,wk,wk,001)
C     do  10  j = 1, n
C     do  10  i = 1, j
C       a11(j,i) = a11(i,j)
C  10 continue
C ... Lapack version
      call dsytrf('U',n,a11,n,kpvt,wk,n*64,ndef)
      if (ndef .ne. 0) call rx('decmp2: matrix a11 singular')
      call dsytri('U',n,a11,n,kpvt,wk,ndef)
      if (ndef .ne. 0) call rx('decmp2: matrix a11 singular')
      do  10  j = 1, n
      do  10  i = 1, j
        a11(j,i) = a11(i,j)
   10 continue
C     call prmx('a11^-1',a11,n,n,n)

C --- Make a22 - a21 a11^-1 a12 ---
C     k = 0
C     do  30  j = 1, n
C     do  30  i = 1, j
C       k = k+1
C       a22(k) = a22(k) - a21(i)*a11(k)*a12(j)
C  30 continue
      do  30  j = 1, n
      do  30  i = 1, j
        a22(i,j) = a22(i,j) - a21(i)*a11(i,j)*a12(j)
   30 continue

C     call dpack(a11,n)
C     call dpack(a22,n)

C --- Decompose (a22 - a21 a11^-1 a12) ---
C ... Packed Linpack version
C     call dspfa(a22,n,kpvt,ndef)
C     if (ndef .ne. 0) call rx('decmp2: matrix a22 singular')
C ... Linpack version
C     call dsifa(a22,n,n,kpvt,ndef)
C     if (ndef .ne. 0) call rx('decmp2: matrix a22 singular')
C ... Lapack version
      call dsytrf('U',n,a22,n,kpvt,wk,n*64,ndef)
      if (ndef .ne. 0) call rx('decmp2: matrix a22 singular')

      call tcx('decmp2')
      end


