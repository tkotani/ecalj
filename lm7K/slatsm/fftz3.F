c takao modified this based on delw_2.2.
      subroutine fftz3(c,n1,n2,n3,k1,k2,k3,nfft,iset,isig)
C     implicit none
      integer n1,n2,n3,k1,k2,k3,nfft,iset,isig
      double complex c(k1,k2,k3,nfft)
C Local variables
      integer i1,i2,i3,id,iord,iopt,ow1,ow2,ow3,oiwk,ierr,ifft
      double precision scale
      save ow1,ow2,ow3,oiwk

C --- A public-domain fft package.  See http://www.fftw.org/ ---
C ... Start of include file fftw_f77.i that comes with the fftw package
c     This file contains PARAMETER statements for various constants
c     that can be passed to FFTW routines.  You should include
c     this file in any FORTRAN program that calls the fftw_f77
c     routines (either directly or with an #include statement
c     if you use the C preprocessor).
      integer FFTW_FORWARD,FFTW_BACKWARD
      parameter (FFTW_FORWARD=-1,FFTW_BACKWARD=1)

      integer FFTW_REAL_TO_COMPLEX,FFTW_COMPLEX_TO_REAL
      parameter (FFTW_REAL_TO_COMPLEX=-1,FFTW_COMPLEX_TO_REAL=1)

      integer FFTW_THREADSAFE
      parameter (FFTW_THREADSAFE=128)
C ... End of include file fftw_f77.i that comes with the fftw package
      integer plan(2),jsig

      integer FFTW_ESTIMATE,FFTW_MEASURE
      INTEGER FFTW_PRESERVE_INPUT
      integer FFTW_OUT_OF_PLACE,FFTW_IN_PLACE,FFTW_USE_WISDOM

      parameter (FFTW_ESTIMATE=64,FFTW_MEASURE=0)
      parameter (FFTW_PRESERVE_INPUT=16)

      call tcn('fftz3')
      if (n1 .eq. 1 .and. n2 .eq. 1 .and. n3 .eq. 1) goto 99

      jsig = 0
      if (isig .eq. -1) jsig = FFTW_FORWARD
      if (isig .eq.  1) jsig = FFTW_BACKWARD
      do  10  ifft = 1, nfft
      if (n2 .eq. 1 .and. n3 .eq. 1) then
        call dfftw_plan_dft_1d(plan,n1,c(1,1,1,ifft),c(1,1,1,ifft),
     .    jsig,FFTW_ESTIMATE)
      elseif (n3 .eq. 1) then
        call dfftw_plan_dft_2d(plan,n1,n2,c(1,1,1,ifft),c(1,1,1,ifft),
     .    jsig,FFTW_ESTIMATE)
      else
        call dfftw_plan_dft_3d(plan,n1,n2,n3,c(1,1,1,ifft),c(1,1,1,ifft)
     .    ,jsig,FFTW_ESTIMATE)
      endif
      call dfftw_execute(plan)
      call dfftw_destroy_plan(plan)
   10 continue

C ... Renormalize forward transform
      if (isig .gt. 0 .or. isig .lt. -10) goto 99
      scale = 1/dble(n1*n2*n3)
      do  20  ifft = 1, nfft
      do  20  i3 = 1, n3
      do  20  i2 = 1, n2
      do  20  i1 = 1, n1
   20 c(i1,i2,i3,ifft) = scale*c(i1,i2,i3,ifft)
C     call zprm3('c',0,c,n1,n2,n3)
   99 call tcx('fftz3')
      return

c      entry fftz3s(n1,n2,n3,k1,k2,k3,iset)
      end

C$$$      subroutine fftz3c(c,f,n1,n2,n3,k1,k2,k3,isw,isig)
C$$$C- Convolution of double complex function with a filter by FFT
C$$$C ----------------------------------------------------------------------
C$$$Ci Inputs
C$$$Ci   isw   1s digit
C$$$Ci         1  do not FT input f
C$$$Ci         2  do not FT input c
C$$$Ci         3  combination of 1+2
C$$$Ci        10s digit
C$$$Ci         1  on output, do not inverse FT f, but return f~
C$$$Ci         2  on output, do not inverse FT c, but return c~
C$$$Ci         3  combination of 1+2
C$$$Ci isig    sense of 'forward' and 'inverse' transforms
C$$$Cr         Use isig=-1 for usual definition of 'forward'
C$$$Cr            f~ = sum_x f(x) exp(-i q x) / (n1*n2*n3)
C$$$Cr         Use isig=1 for reverse definition:
C$$$Cr            f~ = sum_x f(x) exp( i q x)
C$$$Co Outputs
C$$$Co   c is overwritten by (c~ * f~)~, where ~ is FT
C$$$C ----------------------------------------------------------------------
C$$$C     implicit none
C$$$      integer n1,n2,n3,k1,k2,k3,isw,isig
C$$$      double complex c(k1,k2,k3),f(k1,k2,k3)
C$$$C Local variables
C$$$      integer isw0,isw1,i1,i1m,i2,i2m,i3,i3m
C$$$      isw0 = mod(isw,10)
C$$$      isw1 = mod(isw/10,10)

C$$$      if (isw0/2 .eq. 0)      call fftz3(c,n1,n2,n3,k1,k2,k3,1,0,isig)
C$$$      if (mod(isw0,2) .eq. 0) call fftz3(f,n1,n2,n3,k1,k2,k3,1,0,isig)

C$$$      do  10  i3 = 1, n3
C$$$      i3m = n3+2-i3
C$$$      if (i3 .eq. 1) i3m = 1

C$$$        do  20  i2 = 1, n2
C$$$        i2m = n2+2-i2
C$$$        if (i2 .eq. 1) i2m = 1

C$$$          do  30  i1 = 1, n1
C$$$          i1m = n1+2-i1
C$$$          if (i1 .eq. 1) i1m = 1

C$$$            c(i1,i2,i3) = c(i1,i2,i3)*f(i1m,i2m,i3m)
C$$$   30     continue
C$$$   20   continue
C$$$   10 continue

C$$$      if (isw1/2 .eq. 0)      call fftz3(c,n1,n2,n3,k1,k2,k3,1,0,-isig)
C$$$      if (mod(isw1,2) .eq. 0) call fftz3(f,n1,n2,n3,k1,k2,k3,1,0,-isig)
C$$$      end


