      subroutine ztoyy(a,lda,ldb,nr,nc,kcplxi,kcplxf)
C- Separates real and imaginary parts of subblock of a complex*16 matrix
C ----------------------------------------------------------------------
Ci Inputs
Ci   a     :the matrix to separate or recombine
Ci   lda   :formal leading dimension of a
Ci   ldb   :formal second dimension of a
Ci   nr    :number of rows in a
Ci   nc    :number of columns in a
Ci   kcplxi:initial complex storage mode.  Matrix a is formally
Ci          dimensioned a(lda,nc).  The true dimensions of a are:
Ci          0: a has real, imaginary separated and
Ci             a = a(lda,nc,1..2) with Re a = a(lda,nc)
Ci          1: a is in complex*16 format, i.e. a is dimensioned
Ci             a = a(1..2,lda,nc)
Ci          2: Real, imaginary parts separated by columns; a dimensioned
Ci             a = a(lda,2,nc), with a(*,1..2,*) = real..imag
Ci   kcplxf:final complex storage mode, with conventions as in kcplxi
Co Outputs
Co   a     :transformed from complex storage mode kcplxi to kcplxf.
Cr Remarks
Cr   subblock a(nr+1..lda,*) is overwritten when converting kcplx=0->1
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
      integer lda,ldb,nr,nc,kcplxi,kcplxf
      double precision a(1)
      integer lbak,kcplxl
      double precision wkr(nr),wkc(nc)
      if ( kcplxi .eq. kcplxf ) return

      call isanrg(kcplxi,0,2,'ztoyy','kcplxi',.true.)
      call isanrg(kcplxf,0,2,'ztoyy','kcplxf',.true.)

      kcplxl = kcplxi

C ... Convert kcplx=0 -> kcplx=2
      if (kcplxl .eq. 0) then
        call ztoy2(a,lda,ldb,nr,nc,wkc,1)
        kcplxl = 2
      endif
      if (kcplxl .eq. kcplxf) goto 100

C ... Convert kcplx=1 -> kcplx=2 or kcplx=2 -> kcplx=1
      if (kcplxl .ne. 0 .and. (kcplxl .eq. 1 .or. kcplxf .eq. 1)) then
        lbak = kcplxl-1
        call ztoy1(a,lda,nr,nc,wkr,lbak)
        kcplxl = 3-kcplxl
      endif
      if (kcplxl .eq. kcplxf) goto 100

C ... Convert kcplx=2 -> kcplx=0
      if (kcplxl .eq. 2 .and. kcplxf .eq. 0) then
        call ztoy2(a,lda,ldb,nr,nc,wkc,0)
        kcplxl = 0
      endif
C ... cleanup
  100 continue
      end

      subroutine ztoy2(a,lda,ldb,nr,nc,wk,lbak)
C- Swaps imaginary part of a into second index
C ----------------------------------------------------------------------
Ci Inputs
Ci   a     :the matrix to separate or recombine
Ci   lda   :dimension of a
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
Cr     a is effectively dimensioned a(2,lda,ldb)
Cr       with a(1,*,*) = real part and a(2,*,*) = imaginary part
Cr     Output a is stored with real and imaginary separated by columns
Cr     a is effectively dimensioned a(lda,2,ldb)
Cr       with a(*,1,*) = real part and a(*,2,*) = imaginary part
Cr
Cr   Case lbak=1: the process is reversed.
Cu Updates
Cu  20 Jul 1999 changed convention in cases lda>nr
C ----------------------------------------------------------------------
C     implicit none
      integer lda,ldb,nr,nc,lbak
      double precision a(lda,ldb),wk(nc)
      integer i

      if (nr .gt. lda) call rx('ztoy2: subblock exceeds lda')
      if (lbak .ne. 0) goto 20

C --- Change from a(lda,2,nc) to a(lda,nc,2) ---
      do  10  i = 1, nr
C       Copy imaginary part into wk
        call dcopy(nc,a(i,2),2*lda,wk,1)
C       Copy real part into a
        call dcopy(nc,a(i,1),2*lda,a(i,1),lda)
C       Copy imaginary part into a
        call dcopy(nc,wk,1,a(i,ldb+1),lda)
   10 continue
      return

C --- Change from a(lda,ldb,2) to a(lda,2,ldb) ---
   20 continue
      do  30  i = 1, nr
C       Copy real part a(*,1..nc) into wk
        call dcopy(nc,a(i,1),lda,wk,1)
C       Copy imaginary part a(*,1+ldb..nc+ldb) into a(*,2,1..nc)
        call dcopy(nc,a(i,ldb+1),lda,a(i,2),2*lda)
C       Copy real part into a(*,1,1..nc)
        call dcopy(nc,wk,1,a(i,1),2*lda)
   30 continue
      end

C$$$#if TEST
C$$$      subroutine fmain
C$$$C      implicit none
C$$$      integer lda,ldb,nr,nc,ir,ic,i1mach,ld1,ld2,ldr,oi,kcplxi,kcplxf,
C$$$     .  lerr
C$$$      parameter (lda=5,ldb=6,nr=4,nc=3)
C$$$      double precision yy(lda,ldb,2),y(lda,2,ldb),z(2,lda,ldb)
C$$$      double precision yy0(lda,ldb,2),diff
C$$$      equivalence (yy,y)
C$$$      equivalence (yy,z)

C$$$      integer w(1000)
C$$$      common /w/ w
C$$$      call wkinit(1000)

C$$$      lerr = 0

C$$$      do  10  ir = 1, lda
C$$$      do  10  ic = 1, ldb
C$$$        yy(ir,ic,1) = 100*ir + ic
C$$$        yy(ir,ic,2) = -(100*ir + ic)
C$$$        yy0(ir,ic,1) = 100*ir + ic
C$$$        yy0(ir,ic,2) = -(100*ir + ic)
C$$$   10 continue
C$$$      kcplxi = 0

C$$$      call cplxdm(kcplxi,lda,ldb,ld1,ld2,ldr,oi)
C$$$      call ywrm(0,'yy',kcplxi+2,i1mach(2),'(6f8.2)',yy,oi,lda,nr,ldb)

C$$$C ... 0 -> 2
C$$$      kcplxf = 2
C$$$      call ztoyy(yy,lda,ldb,nr,nc,kcplxi,kcplxf)
C$$$      diff = 0
C$$$      do  20  ir = 1, nr
C$$$      do  20  ic = 1, nc
C$$$        diff = diff + abs(y(ir,1,ic)-yy0(ir,ic,1)) +
C$$$     .                abs(y(ir,2,ic)-yy0(ir,ic,2))
C$$$   20 continue
C$$$      print 333, kcplxi,kcplxf,diff,yy(5,1,1),yy(5,1,2),yy(5,2,1)
C$$$      if (abs(diff) .gt. 1d-10) lerr = lerr+1
C$$$  333 format(' diff',i2,' ->',i2,' =',f10.2,2x,3f10.2)
C$$$      kcplxi = kcplxf
C$$$C     call cplxdm(kcplxi,lda,ldb,ld1,ld2,ldr,oi)
C$$$C     call ywrm(0,'y',kcplxi+2,i1mach(2),'(6f8.2)',yy,oi,lda,nr,nc)

C$$$C ... 2 -> 0
C$$$      kcplxf = 0
C$$$      call ztoyy(yy,lda,ldb,nr,nc,kcplxi,kcplxf)
C$$$      diff = 0
C$$$      do  50  ir = 1, nr
C$$$      do  50  ic = 1, nc
C$$$        diff = diff + abs(yy(ir,ic,1)-yy0(ir,ic,1)) +
C$$$     .                abs(yy(ir,ic,2)-yy0(ir,ic,2))
C$$$   50 continue
C$$$      print 333, kcplxi,kcplxf,diff,yy(5,1,1),yy(5,1,2),yy(5,2,1)
C$$$      if (abs(diff) .gt. 1d-10) lerr = lerr+1
C$$$      kcplxi = kcplxf
C$$$C      call cplxdm(kcplxi,lda,ldb,ld1,ld2,ldr,oi)
C$$$C      call ywrm(0,'y',kcplxi+2,i1mach(2),'(6f8.2)',yy,oi,lda,nr,nc)

C$$$C ... 0 -> 1
C$$$      kcplxf = 1
C$$$      call ztoyy(yy,lda,ldb,nr,nc,kcplxi,kcplxf)
C$$$      diff = 0
C$$$      do  60  ir = 1, nr
C$$$      do  60  ic = 1, nc
C$$$        diff = diff + abs(z(1,ir,ic)-yy0(ir,ic,1)) +
C$$$     .                abs(z(2,ir,ic)-yy0(ir,ic,2))
C$$$   60 continue
C$$$      print 333, kcplxi,kcplxf,diff,yy(5,1,1),yy(5,1,2),yy(5,2,1)
C$$$      if (abs(diff) .gt. 1d-10) lerr = lerr+1
C$$$      kcplxi = kcplxf
C$$$C      call cplxdm(kcplxi,lda,ldb,ld1,ld2,ldr,oi)
C$$$C      call ywrm(0,'y',kcplxi+2,i1mach(2),'(6f8.2)',yy,oi,lda,nr,nc)

C$$$C ... 1 -> 2
C$$$      kcplxf = 2
C$$$      call ztoyy(yy,lda,ldb,nr,nc,kcplxi,kcplxf)
C$$$      diff = 0
C$$$      do  40  ir = 1, nr
C$$$      do  40  ic = 1, nc
C$$$        diff = diff + abs(y(ir,1,ic)-yy0(ir,ic,1)) +
C$$$     .                abs(y(ir,2,ic)-yy0(ir,ic,2))
C$$$   40 continue
C$$$      print 333, kcplxi,kcplxf,diff,yy(5,1,1),yy(5,1,2),yy(5,2,1)
C$$$      if (abs(diff) .gt. 1d-10) lerr = lerr+1
C$$$      kcplxi = kcplxf
C$$$C     call cplxdm(kcplxi,lda,ldb,ld1,ld2,ldr,oi)
C$$$C     call ywrm(0,'y',kcplxi+2,i1mach(2),'(6f8.2)',yy,oi,lda,nr,nc)

C$$$C ... 2 -> 1
C$$$      kcplxf = 1
C$$$      call ztoyy(yy,lda,ldb,nr,nc,kcplxi,kcplxf)
C$$$      diff = 0
C$$$      do  30  ir = 1, nr
C$$$      do  30  ic = 1, nc
C$$$        diff = diff + abs(z(1,ir,ic)-yy0(ir,ic,1)) +
C$$$     .                abs(z(2,ir,ic)-yy0(ir,ic,2))
C$$$   30 continue
C$$$      print 333, kcplxi,kcplxf,diff,yy(5,1,1),yy(5,1,2),yy(5,2,1)
C$$$      if (abs(diff) .gt. 1d-10) lerr = lerr+1
C$$$      kcplxi = kcplxf
C$$$C     call cplxdm(kcplxi,lda,ldb,ld1,ld2,ldr,oi)
C$$$C     call ywrm(0,'y',kcplxi+2,i1mach(2),'(6f8.2)',yy,oi,lda,nr,nc)

C$$$C ... 1 -> 0
C$$$      kcplxf = 0
C$$$      call ztoyy(yy,lda,ldb,nr,nc,kcplxi,kcplxf)
C$$$      diff = 0
C$$$      do  70  ir = 1, nr
C$$$      do  70  ic = 1, nc
C$$$        diff = diff + abs(yy(ir,ic,1)-yy0(ir,ic,1)) +
C$$$     .                abs(yy(ir,ic,2)-yy0(ir,ic,2))
C$$$   70 continue
C$$$      print 333, kcplxi,kcplxf,diff,yy(5,1,1),yy(5,1,2),yy(5,2,1)
C$$$      if (abs(diff) .gt. 1d-10) lerr = lerr+1
C$$$      kcplxi = kcplxf
C$$$      call cplxdm(kcplxi,lda,ldb,ld1,ld2,ldr,oi)
C$$$      call ywrm(0,'y',kcplxi+2,i1mach(2),'(6f8.2)',yy,oi,lda,nr,nc)

C$$$      print *, ' '
C$$$      if (lerr .eq. 0) then
C$$$        print *, 'ztoyy : all conversions successful'
C$$$      else
C$$$        print *, 'ztoyy : ',lerr, 'erroneous conversions'
C$$$      endif

C$$$      call cexit(lerr,1)
C$$$      end

C$$$#endif

