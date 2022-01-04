c$$$      subroutine cplxdm(kcplx,lda,ldb,ld1,ld2,nca,oi)
c$$$C- True first and second dimensions a complex matrix
c$$$C ----------------------------------------------------------------------
c$$$Ci Inputs
c$$$Ci   kcplx : specifies how complex matrix is matrix stored
c$$$Ci           0: a has real, imaginary separated with storage:
c$$$Ci              a = a(lda,ldb,2), with a(*,*,1..2) = real..imag
c$$$Ci           1: a is returned complex*16 format:
c$$$Ci              a = a(2,lda,ldb), with a(1..2,*,*) = real..imag
c$$$Ci           2: a has real, imaginary separated by columns
c$$$Ci              a = a(lda,2,ldb), with a(*,1..2,*) = real..imag
c$$$Ci   lda:    first formal dimension; see kcplx
c$$$Ci   ldb:    second formal dimension; see kcplx
c$$$Ci           ldb is only used if kcplx=0
c$$$Co Outputs
c$$$Co   ld1:    first true dimension:
c$$$Co           lda if kcplx=0
c$$$Co             2 if kcplx=1
c$$$Co           lda if kcplx=2
c$$$Co   ld2:    second true dimension
c$$$Co           ldb if kcplx=0
c$$$Co           lda if kcplx=1
c$$$Co             2 if kcplx=2
c$$$Co   nca:    Spacing between column elements
c$$$Co           lda    kcplx=0
c$$$Co           lda*2  kcplx=1
c$$$Co           lda*2  kcplx=2
c$$$Co    oi:    offset to first imaginary element:
c$$$Co           lda*ldb kcplx=0
c$$$Co           1       kcplx=1
c$$$Co           lda     kcplx=2
c$$$Cr Remarks
c$$$Cr   cplxdm attempts to isolate, and thus facilitate reconciling
c$$$Cr   incompatibilities in the differing complex representations.
c$$$Cr
c$$$Cr   ld1 and ld2 are the true first and second dimensions of a complex
c$$$Cr   matrix a(lda,ldb), whose formal leading dimensions are lda,ldb.
c$$$Cr   nca is the spacking between column elements.
c$$$Cu Updates
c$$$Cu   29 Jan 01 Changed meaning of nca (2nd to last argument)
c$$$C ----------------------------------------------------------------------
c$$$C     implicit none
c$$$      integer kcplx,lda,ldb,ld1,ld2,nca,oi
c$$$      logical:: l_dummy_isanrg,isanrg
c$$$
c$$$      if (kcplx .eq. 0) then
c$$$        ld1 = lda
c$$$        ld2 = ldb
c$$$        nca = lda
c$$$        oi  = lda*ldb
c$$$      elseif (kcplx .eq. 1) then
c$$$        ld1 = 2
c$$$        ld2 = lda
c$$$        nca = lda*2
c$$$        oi  = 1
c$$$      elseif (kcplx .eq. 2) then
c$$$        ld1 = lda
c$$$        ld2 = 2
c$$$        nca = lda*2
c$$$        oi  = lda
c$$$      else
c$$$Ckino isanrg is logical function,         call isanrg(kcplx,0,2,'cplxdm:','kcplx',.true.)
c$$$        l_dummy_isanrg=isanrg(kcplx,0,2,'cplxdm:','kcplx',.true.)
c$$$      endif
c$$$
c$$$      end
c$$$
