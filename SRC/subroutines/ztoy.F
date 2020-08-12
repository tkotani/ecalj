      subroutine ztoy(a,lda,nr,nc,lbak)
C- Separates real and imaginary parts of subblock of a complex*16 matrix
C ----------------------------------------------------------------------
Ci Inputs
Ci   a     :the matrix to separate or recombine
Ci   lda   :leading dimension of a
Ci   nr    :number of rows in a
Ci   nc    :number of columns in a
Ci   wk    :real work array of length nr
Ci   lbak  :0 separates real from imaginary, as described in Remarks
Ci         :1 recombines into double complex format
Co Outputs
Co   a     :transformed between complex*16 and real storage; see Remarks
Cr Remarks
Cr   Case lbak=0: 
Cr     Input a is stored in complex*16 format, i.e. 
Cr     a is effectively dimensioned a(2,lda,nc)
Cr       with a(1,*,*) = real part and a(2,*,*) = imaginary part
Cr     Output a is stored with real and imaginary separated by columns
Cr     a is effectively dimensioned a(lda,2,nc)
Cr       with a(*,1,*) = real part and a(*,2,*) = imaginary part
Cr
Cr   Case lbak=1: the process is reversed.
Cu Updates
Cu   28 Aug 04 (Kino) Added F90 code to avoid problems with dcopy
Cu   20 Jul 99 changed convention in cases lda>nr
C ----------------------------------------------------------------------
C     implicit none
      integer lda,lbak,nr,nc
      double precision a(1)
      real(8),allocatable:: wk(:)
      allocate(wk(nr))
      call ztoy1(a,lda,nr,nc,wk,lbak)
      deallocate(wk)
      end
      subroutine ztoy1(a,lda,nr,nc,wk,lbak)
C- Separates real and imaginary parts of subblock of a complex*16 matrix
C ----------------------------------------------------------------------
Ci Inputs
Ci   a     :the matrix to separate or recombine
Ci   lda   :leading dimension of a
Ci   nr    :number of rows in a
Ci   nc    :number of columns in a
Ci   wk    :real work array of length nr
Ci   lbak  :0 separates real from imaginary, as described in Remarks
Ci         :1 recombines into double complex format
Co Outputs
Co   a     :transformed between complex*16 and real storage; see Remarks
Cr Remarks
Cr   Case lbak=0: 
Cr     Input a is stored in complex*16 format, i.e. 
Cr     a is effectively dimensioned a(2,lda,nc)
Cr       with a(1,*,*) = real part and a(2,*,*) = imaginary part
Cr     Output a is stored with real and imaginary separated by columns
Cr     a is effectively dimensioned a(lda,2,nc)
Cr       with a(*,1,*) = real part and a(*,2,*) = imaginary part
Cr
Cr   Case lbak=1: the process is reversed.
Cu Updates
Cu  20 Jul 1999 changed convention in cases lda>nr
C ----------------------------------------------------------------------
C     implicit none
      integer lda,lbak,nr,nc
      double precision a(1),wk(nr)
      double precision,allocatable:: wk2(:)
      integer j,k,l,mdim

      mdim = 2*lda
      if (nr .gt. lda) call rx('ztoy: subblock exceeds array dimension')
      if (lbak .ne. 0) goto 20

C --- Separate real and imaginary parts ---
C ..  k is offset to col j, l to imaginary part
      allocate(wk2(nr))
      do  10  j = 1, nc
        k = (j-1)*mdim + 1
        l = k + lda
        call dcopy(nr,a(k+1),2,wk,1)
        call dcopy(nr,a(k),2,wk2,1)
        call dcopy(nr,wk2,1,a(k),1)
        call dcopy(nr,wk,1,a(l),1)
   10 continue
      deallocate(wk2)
      return

C --- Restore to complex storage ---
   20 continue
      allocate(wk2(nr))
      do  30  j = 1, nc
        k = (j-1)*mdim + 1
        l = k + lda
        call dcopy(nr,a(k),1,wk,1)
        call dcopy(nr,a(l),1,wk2,1)
        call dcopy(nr,wk2,1,a(k+1),2)
        call dcopy(nr,wk,1,a(k),2)
   30 continue
      deallocate(wk2)
      end

C$$$#if TEST
C$$$      subroutine fmain
C$$$      double precision a(2,4,2),z(8),aa(16)
C$$$      data a /1d0,11d0,2d0,12d0,3d0,13d0,4d0,14d0,
C$$$     .        5d0,15d0,6d0,16d0,7d0,17d0,8d0,18d0/
C$$$      equivalence(a,aa)
C$$$      integer w(1000)
C$$$      call wkinit(1000)
C$$$      lda = 4
C$$$      nr  = 3
C$$$      nc  = 2
C$$$      call zwrm('(2f12.6)',a,lda,nr,nc)
C$$$      print *, 'call ztoy with lbak=',0
C$$$      call ztoy(a,lda,nr,nc,0)
C$$$      call ywrm('(2f12.6)',a,lda,nr,nc)
C$$$      print *, 'call ztoy with lbak=',1
C$$$      call ztoy(a,lda,nr,nc,1)
C$$$      call zwrm('(2f12.6)',a,lda,nr,nc)
C$$$      end

C$$$      subroutine zwrm(fmt,s,ns,nr,nc)
C$$$      double precision s(0:1,ns,nc)
C$$$      character*(8) fmt
C$$$      print *, 'complex*16 s'
C$$$      print *, nr, nc
C$$$      do  10  i = 1, nr
C$$$   10 print fmt, (s(0,i,j), j=1,nc)
C$$$      do  20  i = 1, nr
C$$$   20 print fmt, (s(1,i,j), j=1,nc)
C$$$      end

C$$$      subroutine ywrm(fmt,s,lds,nr,nc)
C$$$      double precision s(lds,2,nc)
C$$$      character*(8) fmt
C$$$      print *, 'real,imaginary parts separated'
C$$$      print *, nr, nc
C$$$      do  10  i = 1, nr
C$$$   10 print fmt, (s(i,1,j), j=1,nc)
C$$$      do  20  i = 1, nr
C$$$   20 print fmt, (s(i,2,j), j=1,nc)
C$$$      end
C$$$#endif


