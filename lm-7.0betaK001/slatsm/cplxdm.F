      subroutine cplxdm(kcplx,lda,ldb,ld1,ld2,nca,oi)
C- True first and second dimensions a complex matrix
C ----------------------------------------------------------------------
Ci Inputs
Ci   kcplx : specifies how complex matrix is matrix stored
Ci           0: a has real, imaginary separated with storage:
Ci              a = a(lda,ldb,2), with a(*,*,1..2) = real..imag
Ci           1: a is returned complex*16 format:
Ci              a = a(2,lda,ldb), with a(1..2,*,*) = real..imag
Ci           2: a has real, imaginary separated by columns
Ci              a = a(lda,2,ldb), with a(*,1..2,*) = real..imag
Ci   lda:    first formal dimension; see kcplx
Ci   ldb:    second formal dimension; see kcplx
Ci           ldb is only used if kcplx=0
Co Outputs
Co   ld1:    first true dimension:
Co           lda if kcplx=0
Co             2 if kcplx=1
Co           lda if kcplx=2
Co   ld2:    second true dimension
Co           ldb if kcplx=0
Co           lda if kcplx=1
Co             2 if kcplx=2
Co   nca:    Spacing between column elements
Co           lda    kcplx=0
Co           lda*2  kcplx=1
Co           lda*2  kcplx=2
Co    oi:    offset to first imaginary element:
Co           lda*ldb kcplx=0
Co           1       kcplx=1
Co           lda     kcplx=2
Cr Remarks
Cr   cplxdm attempts to isolate, and thus facilitate reconciling
Cr   incompatibilities in the differing complex representations.
Cr
Cr   ld1 and ld2 are the true first and second dimensions of a complex
Cr   matrix a(lda,ldb), whose formal leading dimensions are lda,ldb.
Cr   nca is the spacking between column elements.
Cu Updates
Cu   29 Jan 01 Changed meaning of nca (2nd to last argument)
C ----------------------------------------------------------------------
C     implicit none
      integer kcplx,lda,ldb,ld1,ld2,nca,oi

      if (kcplx .eq. 0) then
        ld1 = lda
        ld2 = ldb
        nca = lda
        oi  = lda*ldb
      elseif (kcplx .eq. 1) then
        ld1 = 2
        ld2 = lda
        nca = lda*2
        oi  = 1
      elseif (kcplx .eq. 2) then
        ld1 = lda
        ld2 = 2
        nca = lda*2
        oi  = lda
      else
        call isanrg(kcplx,0,2,'cplxdm:','kcplx',.true.)
      endif

      end

