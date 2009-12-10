      subroutine ymtrns(mode,a,nca,nra,ofai,b,ncb,nrb,ofbi,n1,n2,m1,m2)
C- Transpose of a complex matrix
C ----------------------------------------------------------------
Ci Inputs:
Ci   mode  :0 in-place tranpose (overwrite a by its transpose)
Ci         :  a(m1:m2; n1:n2) is overwritten by a(n1:n2; m1:m2)
Ci            a(n1:n2; m1:m2) is overwritten by a(m1:m2; n1:n2)
Ci         :1 out-of-place transpose: write transpose of a into b
Ci         :  b(m1:m2; n1:n2) is written
Ci   nca   :number of elements separating columns of a
Ci         : = leading dimension of a if real, imag fully separated
Ci         : = 2* leading dimension of a if real, imag sep by cols
Ci   nra   :number of elements separating rows of a (usu. 1)
Ci   ofai  :offset to imaginary part of a
Ci   n1,n2 :transpose rows n1..n2
Ci   m2,m2 :transpose columns m1..m2
Ci     ... The following are used for mode=1
Ci   ncb   :number of elements separating columns of b
Ci         : = leading dimension of b if real, imag fully separated
Ci         : = 2* leading dimension of b if real, imag sep by cols
Ci   nrb   :number of elements separating rows of b (usu. 1)
Ci   ofbi  :offset to imaginary part of b
Cio Inputs/Outputs
Cio  a     :on input, matrix to be transposed
Cio        :on output, subblock of a is transposed (mode=0)
Co   b     :tranpose of (n1:n2,m1:m2) subblock of matrix a (mode=1)
Cr Remarks
Cr   Transpose is complicated by the avoidance of double-counting
Cr   in the 'diagonal' block of tranpose (region 'x' in figure below)
Cr
Cr                     m1  n1     n2        m2
Cr               _________________________________
Cr               |\    .    .      .           .  |
Cr               |  \  .    .      .           .  |
Cr               |    \.    .      .           .  |
Cr               |     .\   .      .           .  |
Cr            n1 |     .__\_.______._____ _____.  |
Cr               |     |    \      .           |  |
Cr               |     | 2  .x\  1 .  1        |  |
Cr               |     |    .xxx\  .           |  |
Cr            n2 |     |____.xxxxx\.___________|  |
Cr               |                  \             |
Cr               |                    \           |
Cr               |                      \         |
Cr               |                        \       |
Cr               |                          \     |
Cr               |                            \   |
Cr               |______________________________\_|
Cr
C ----------------------------------------------------------------
C     implicit none
      integer mode,nra,nca,nrb,ncb,n1,n2,m1,m2,ofai,ofbi
      double precision a(0:*),b(0:*)
C ... Local parameters
      integer ofij,ofji,n,m,mmin1
      double precision wkr,wki

      call isanrg(mode,0,1,'ymtrns','mode',.true.)
      do  n = n1, n2

C   --- Out-of-place transpose ---
        if (mode .eq. 1) then
          ofij = nra*(n-1) + nca*(m1-1)
          ofji = ncb*(n-1) + nrb*(m1-1)
          do  m = m1, m2
            b(ofji)      = a(ofij)
            b(ofji+ofbi) = a(ofij+ofai)
            ofij = ofij + nca
            ofji = ofji + nrb
          enddo

C   --- In-place transpose ---
        else

C     ... Region 1 in diagram: loop over block m = m1 ... m2
C         unless m1<=n<=m2.  Then m = n ... m2 to exclude region x
          if (n .ge. m1 .and. n .le. m2) then
            mmin1 = n
          else
            mmin1 = m1
          endif
          ofij = nra*(n-1) + nca*(mmin1-1)
          ofji = nca*(n-1) + nra*(mmin1-1)
          do  m = mmin1, m2
            wkr = a(ofij)
            wki = a(ofij+ofai)
            a(ofij)      = a(ofji)
            a(ofij+ofai) = a(ofji+ofai)
            a(ofji)      = wkr
            a(ofji+ofai) = wki
            ofij = ofij + nca
            ofji = ofji + nra
          enddo

C     ... Region 2 in diagram: m = m1...n1 when m1<=n<=m2
          if (n .ge. m1 .and. n .le. m2) then
            ofij = nra*(n-1) + nca*(m1-1)
            ofji = nca*(n-1) + nra*(m1-1)
            do  m = m1, n1-1
              wkr = a(ofij)
              wki = a(ofij+ofai)
              a(ofij)      = a(ofji)
              a(ofij+ofai) = a(ofji+ofai)
              a(ofji)      = wkr
              a(ofji+ofai) = wki
              ofij = ofij + nca
              ofji = ofji + nra
            enddo
          endif
        endif
      enddo

      end

#if TEST
      subroutine fmain
      implicit none
      integer lda,ldb,nr1,nc1,ir,ic,i1mach,ld1,ld2,nca,oi,
     .  kcplxi,so,nr2,nc2,nra
      parameter (lda=7,ldb=8)
      double precision yy(lda,ldb,2),yy2(lda,ldb,2)
      double precision yy0(lda,ldb,2)
      character*8 fm

      so = i1mach(2)
      fm = '(8f8.1)'
      nr1 = 2
      nr2 = 5
      nc1 = 4
      nc2 = 7
      print 333,  nr1,nr2,nc1,nc2
  333 format('#nr1,nr2,nc1,nc2=? (default=',4i3,')')
      read(*,*) nr1,nr2,nc1,nc2

      do  10  ir = 1, lda
      do  10  ic = 1, ldb
        yy0(ir,ic,1) = 100*ir + ic
        yy0(ir,ic,2) = -(100*ir + ic)

   10 continue

      do  30  kcplxi = 0, 0, 1
      call cplxdm(kcplxi,lda,ldb,ld1,ld2,nca,oi)

      call dcopy(lda*ldb*2,yy0,1,yy,1)
      call ztoyy(yy,lda,ldb,lda,ldb,0,kcplxi)
      call dcopy(lda*ldb*2,yy,1,yy2,1)

C ... straight transpose, mode=0
      nra = 1
      if (kcplxi .eq. 1) nra = 2

      call ymtrns(0,yy,nca,nra,oi,yy2,nca,nra,oi,nr1,nr2,nc1,nc2)
      call ywrm(0,'complex copy yy',kcplxi+2,so,fm,yy,oi,lda,lda,ldb)
      call ymtrns(0,yy,nca,nra,oi,yy2,nca,nra,oi,nr1,nr2,nc1,nc2)
      call ztoyy(yy,lda,ldb,lda,ldb,kcplxi,0)
      do  ir = 1, lda
        do  ic = 1, ldb
          if (yy(ir,ic,1) .ne. yy0(ir,ic,1)) stop 'oops'
          if (yy(ir,ic,2) .ne. yy0(ir,ic,2)) stop 'oops'
        enddo
      enddo

C      call ymtrns(1,yy,nca,nra,oi,yy2,nca,nra,oi,nr1,nr2,nc1,nc2)
C      call ywrm(0,'transpose yy2',kcplxi+2,so,fm,yy2,oi,lda,lda,ldb)
C      call ymtrns(0,yy2,nca,nra,oi,yy2,nca,nra,oi,nr1,nr2,nc1,nc2)
C      call ztoyy(yy2,lda,ldb,lda,ldb,kcplxi,0)
Cc     this does not work (nor should it)
C      do  ir = 1, lda
C        do  ic = 1, ldb
C          if (yy2(ir,ic,1) .ne. yy0(ir,ic,1)) print *, ir,ic
C          if (yy2(ir,ic,1) .ne. yy0(ir,ic,1)) stop 'oops'
C          if (yy2(ir,ic,2) .ne. yy0(ir,ic,2)) print *, ir,ic
C          if (yy2(ir,ic,2) .ne. yy0(ir,ic,2)) stop 'oops'
C        enddo
C      enddo


   30 continue
      end
#endif
CTEST   echo 2  5  4  7 | a.out >out
C       diff out out.ymtrns
C       echo 4  7  2  5 | a.out >out
C       diff out out.ymtrns

