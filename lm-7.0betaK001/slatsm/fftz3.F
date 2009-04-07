#define FFTW3 1
#define FFTW 1
#define IN_PLACE 1
      subroutine fftz30(n1,n2,n3,k1,k2,k3)
C- Three-dimensional double complex-to-complex FFT
C ----------------------------------------------------------------------
Ci Inputs and parameters for fft3, fft30, fft3s
Ci    c(k1,k2,k3)   array to transform; overwritten on output (fft3)
Ci    n1,n2,n3      number of elements of the FFT
Ci    k1,k2,k3      dimensions of array c
Ci    isig         -1 for forward transform (real to recip)
Ci                  1 for inverse (assemble an real-space mesh)
Ci                 10s digit add 1 to suppress  / (n1*n2*n3) scaling
Ci    iset          a switch specifying information about setup;
Ci                  see Remarks.
Co Outputs
Cr   See Remarks.
Co   Forward transform: f(q) = sum_x f(x) exp(-i q x) / (n1*n2*n3)
Co   Reverse transform: f(q) = sum_x f(x) exp( i q x)
Cr Remarks
Cr  This code is a generic interface for several FFT routines.
Cr    fftz30(n1,n2,n3, k1,k2,k3)
Cr      Select the minimum dimensions k1,k2,k3 needed for fftz3.
Cr      This was included as an insurance in case some future
Cr      implementation requires k1,k2,k3 different from n1,n2,n3.
Cr      None of the existing ones here do, though some require
Cr      k1,k2,k3 are equal to n1,n2,n3.
Cr
Cr    fftz3s(n1,n2,n3,k1,k2,k3, iset)
Cr      This is an optional, implementation-specific setup routine.
Cr      Some implementations an initialization of internal coefficients
Cr      which need not repeated for subsequent calls.  For repeated
Cr      FFTs, the setup need not be repeated.
Cr      Calling fftz3s with iset=0 allocates memory for, and
Cr      initializes these coeficients.  To free this memory,
Cr      call fftz3s with the value of iset returned by the initial call.
Cr
Cr    fftz3(c,n1,n2,n3,k1,k2,k3,nfft,iset,isig)
Cr      This is the FFT routine.
Cr      Call with isig=-1 for forward, isig=1 for reverse transform:
Cr      Forward transform: f(q) = sum_x f(x) exp(-i q x) / (n1*n2*n3)
Cr      Reverse transform: f(q) = sum_x f(x) exp( i q x)
Cr      Call with iset=0 if fftz3s has not been invoked.
Cr
Cr    fftz3c(c,f,n1,n2,n3,k1,k2,k3,isw,isig)
Cr      Convolution of function c with f using the FFT routine.
C ----------------------------------------------------------------------
C     implicit none
      integer n1,n2,n3,k1,k2,k3
#if GOEDECKER
#if SGI | SGI8
      k1 = n1+1
      k2 = n2+1
      k3 = n3*2
#elif DECA
      k1 = n1
      k2 = n2
      k3 = n3*2
#else
      k1 = n1
      k2 = n2
      k3 = n3*2
#endif
#elif ESSL
c ... For info, here are the strides suggested by essl
C      call stride(n2, n1, inc2x, 'c', 0)
C      call stride(n3, n2*inc2x, inc3x, 'c', 0)
C      write(6,300) n1,n2,n3,inc2x,inc3x
C  300 format('n1,n2,n3=',3i5,'    stride  suggests inc2,inc3=',2i5)
c ... but we just take k1=n1+1
      k1 = n1+1
      k2 = n2
      k3 = n3
#else
      k1 = n1
      k2 = n2
      k3 = n3
#endif
      end

      subroutine fftz3(c,n1,n2,n3,k1,k2,k3,nfft,iset,isig)
C     implicit none
      integer n1,n2,n3,k1,k2,k3,nfft,iset,isig
      double complex c(k1,k2,k3,nfft)
C Local variables
      integer i1,i2,i3,id,iord,iopt,ow1,ow2,ow3,oiwk,ierr,ifft
#if GOEDECKER
      integer ncache,iz
      integer w(1)
      common /w/ w
#endif
      double precision scale
      save ow1,ow2,ow3,oiwk

#if GOEDECKER
      call tcn('fftz3')
      if (n1 .eq. 1 .and. n2 .eq. 1 .and. n3 .eq. 1) goto 99
#if SGI | SGI8
      ncache=5*1024
#elif DECA
      ncache=5*1024
#else
      ncache=5*1024
#endif
      call defrr(ow1, ncache+1)
      iz = 1
      do  10  ifft = 1, nfft
   10 call gdfft(n1,n2,n3,k1,k2,k3/2,c(1,1,1,ifft),mod(isig,10),iz,
     .    w(ow1),ncache)
C ... Move and renormalize fft
      if (iz .eq. 2 .or. isig .gt. 0) then
        id = 0
        if (iz .eq. 2) id = k3/2
        if (isig .gt. 0) then
          do  20  ifft = 1, nfft
          do  20  i3 = 1, n3
          do  20  i2 = 1, n2
          do  20  i1 = 1, n1
   20     c(i1,i2,i3,ifft) = c(i1,i2,i3+id,ifft)
        else
          scale = 1/dble(n1*n2*n3)
          do  22  ifft = 1, nfft
          do  22  i3 = 1, n3
          do  22  i2 = 1, n2
          do  22  i1 = 1, n1
   22     c(i1,i2,i3,ifft) = scale*c(i1,i2,i3+id,ifft)
        endif
        call rlse(ow1)
      endif
#elif SGI | SGI8
      integer w(1)
      common /w/ w

      call tcn('fftz3')
      if (n1 .eq. 1 .and. n2 .eq. 1 .and. n3 .eq. 1) goto 99

C --- The SGI library routine ---
      if (iset .eq. 0) then
        call defcc(ow1, 45+n1+n2+n3)
        call zfft3di(n1,n2,n3,w(ow1))
      endif
      do  10  ifft = 1, nfft
   10 call zfft3d(mod(isig,10),n1,n2,n3,c(1,1,1,ifft),n1,n2,w(ow1))
      if (iset .eq. 0) call rlse(ow1)
C ... Renormalize forward transform
      if (isig .gt. 0 .or. isig .lt. -10) goto 99
      scale = 1/dble(n1*n2*n3)
      do  20  ifft = 1, nfft
      do  20  i3 = 1, n3
      do  20  i2 = 1, n2
      do  20  i1 = 1, n1
   20 c(i1,i2,i3,ifft) = scale*c(i1,i2,i3,ifft)

#elif ESSL
      integer w(1)
      common /w/ w
      integer naux,isign,inc2,inc3
      call tcn('fftz3')

      call rx('FFTZ3: this version does not seem to work')
c ... Work array
      naux = 60000
      if (max(n2,n3).ge.252) call rx('fftz3(essl): bad naux')
      call defrr(ow1, naux)

c ... Do the FFT
      isign = -isig
      inc2 = k1
      inc3 = k1*k2
      scale = 1d0
      if (isig .lt. 0) scale = 1d0/dble(n1*n2*n3)
      if (isig .lt. -10) scale = 1
      do  10  ifft = 1, nfft
      call dcft3(c(1,1,1,ifft),inc2,inc3,c,inc2,inc3,n1,n2,n3,isign,
     .           scale,w(ow1),naux)
 10   continue
      call rlse(ow1)

C ... Renormalize forward transform ... already done by scale
C      if (isig .gt. 0 .or. isig .lt. -10) goto 99
C      scale = 1/dble(n1*n2*n3)
C      do  20  ifft = 1, nfft
C      do  20  i3 = 1, n3
C      do  20  i2 = 1, n2
C      do  20  i1 = 1, n1
C   20 c(i1,i2,i3,ifft) = scale*c(i1,i2,i3,ifft)

#elif FFTW | FFTW3
C
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

#if FFTW3
      integer FFTW_ESTIMATE,FFTW_MEASURE
      parameter (FFTW_ESTIMATE=64,FFTW_MEASURE=0)
      INTEGER FFTW_PRESERVE_INPUT
      parameter (FFTW_PRESERVE_INPUT=16)
#else
      integer FFTW_ESTIMATE,FFTW_MEASURE
      parameter (FFTW_ESTIMATE=0,FFTW_MEASURE=1)
      integer FFTW_OUT_OF_PLACE,FFTW_IN_PLACE,FFTW_USE_WISDOM
      parameter (FFTW_OUT_OF_PLACE=0)
      parameter (FFTW_IN_PLACE=8,FFTW_USE_WISDOM=16)
#endif

      integer FFTW_THREADSAFE
      parameter (FFTW_THREADSAFE=128)
C ... End of include file fftw_f77.i that comes with the fftw package
      integer plan(2),jsig

      call tcn('fftz3')
      if (n1 .eq. 1 .and. n2 .eq. 1 .and. n3 .eq. 1) goto 99

#if FFTW3
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
#else
#if IN_PLACE
      call fftw3d_f77_create_plan(plan,n1,n2,n3,isig,
     .                            FFTW_ESTIMATE+FFTW_IN_PLACE)
      do  10  ifft = 1, nfft
   10 call fftwnd_f77_one(plan,c(1,1,1,ifft),0)
#else
      call fftw3d_f77_create_plan(plan,n1,n2,n3,isig,
     .                            FFTW_ESTIMATE+FFTW_OUT_OF_PLACE)
      call defcc(ow1, n1*n2*n3)
      call fftwnd_f77_one(plan,c(1,1,1,ifft),w(ow1))
      call dcopy(2*n1*n2*n3,w(ow1),1,c(1,1,1,ifft),1)
      call rlse(ow1)
#endif
      call fftwnd_f77_destroy_plan(plan)
#endif

C ... Renormalize forward transform
      if (isig .gt. 0 .or. isig .lt. -10) goto 99
      scale = 1/dble(n1*n2*n3)
      do  20  ifft = 1, nfft
      do  20  i3 = 1, n3
      do  20  i2 = 1, n2
      do  20  i1 = 1, n1
   20 c(i1,i2,i3,ifft) = scale*c(i1,i2,i3,ifft)

#else
C --- A package taken from Prof. Ho, Davis. ---
      integer w(1)
      common /w/ w

      call tcn('fftz3')
      if (n1 .eq. 1 .and. n2 .eq. 1 .and. n3 .eq. 1) goto 99
      id = k1
      iord = 1

C ... Work arrays for transform
      iopt = 0
      if (n3 .lt. 32) iopt = 1
      if (iset .eq. 0) then
        call defrr(ow1, 4*n1+14)
        call defrr(ow2, 4*n2*((1-iopt)+iopt*(id+1))+14)
        call defrr(ow3, 4*n3+14)
        call defrr(oiwk,max(n1,n2,n3))
C   ... Set up for FFT
        call c3fft(c,id,n1,n2,n3,w(ow1),w(ow2),w(ow3),iopt,0,
     .    iord,w(oiwk),ierr)
        if (ierr .eq. 2) call rx('c3fft needs prime factors 2,3,5')
        if (ierr .ne. 0) call rx('c3fft called improperly for isig=0')
      endif
      do  10  ifft = 1, nfft
   10 call c3fft(c(1,1,1,ifft),id,n1,n2,n3,w(ow1),w(ow2),w(ow3),iopt,
     .  mod(isig,10),iord,w(oiwk),ierr)
      if (ierr .ne. 0) call rx('c3fft called improperly')
      if (iset .eq. 0) call rlse(ow1)

C ... Renormalize forward transform
      if (isig .gt. 0 .or. isig .lt. -10) goto 99
      scale = 1/dble(n1*n2*n3)
      do  20  ifft = 1, nfft
      do  20  i3 = 1, n3
      do  20  i2 = 1, n2
      do  20  i1 = 1, n1
   20 c(i1,i2,i3,ifft) = scale*c(i1,i2,i3,ifft)
#endif

C     call zprm3('c',0,c,n1,n2,n3)
   99 call tcx('fftz3')
      return

      entry fftz3s(n1,n2,n3,k1,k2,k3,iset)

#if SGI | SGI8
      if (iset .eq. 0) then
        call defcc(ow1, 45+n1+n2+n3)
        call zfft3di(n1,n2,n3,w(ow1))
        iset = 1
      else
        call rlse(ow1)
        iset = 0
      endif
#elif FFTW | FFTW3
C     Nothing needed using FFTW
#else
      if (iset .eq. 0) then
        id = k1
        iopt = 0
        if (n3 .lt. 32) iopt = 1
        iord = 1
        call defrr(ow1, 4*n1+14)
        call defrr(ow2, 4*n2*((1-iopt)+iopt*(id+1))+14)
        call defrr(ow3, 4*n3+14)
        call defrr(oiwk,max(n1,n2,n3))
        call c3fft(c,id,n1,n2,n3,w(ow1),w(ow2),w(ow3),iopt,0,
     .    iord,w(oiwk),ierr)
        if (ierr .eq. 2) call rx('c3fft needs prime factors 2,3,5')
        if (ierr .ne. 0) call rxi('c3fft called improperly, err=',ierr)
        iset = 1
      else
        call rlse(ow1)
        iset = 0
      endif
#endif

      end

      subroutine gtpfac(npfac,pfac)
C- Returns allowed prime factors in integers for uniform mesh
C     implicit none
      integer npfac,pfac(5)
#if SGI | SGI8 | FFTW | FFTW3
      npfac = 5
      pfac(1) = 2
      pfac(2) = 3
      pfac(3) = 5
      pfac(4) = 7
      pfac(5) = 11
#else
      npfac = 3
      pfac(1) = 2
      pfac(2) = 3
      pfac(3) = 5
#endif
      end

      subroutine fftz3c(c,f,n1,n2,n3,k1,k2,k3,isw,isig)
C- Convolution of double complex function with a filter by FFT
C ----------------------------------------------------------------------
Ci Inputs
Ci   isw   1s digit
Ci         1  do not FT input f
Ci         2  do not FT input c
Ci         3  combination of 1+2
Ci        10s digit
Ci         1  on output, do not inverse FT f, but return f~
Ci         2  on output, do not inverse FT c, but return c~
Ci         3  combination of 1+2
Ci isig    sense of 'forward' and 'inverse' transforms
Cr         Use isig=-1 for usual definition of 'forward'
Cr            f~ = sum_x f(x) exp(-i q x) / (n1*n2*n3)
Cr         Use isig=1 for reverse definition:
Cr            f~ = sum_x f(x) exp( i q x)
Co Outputs
Co   c is overwritten by (c~ * f~)~, where ~ is FT
C ----------------------------------------------------------------------
C     implicit none
      integer n1,n2,n3,k1,k2,k3,isw,isig
      double complex c(k1,k2,k3),f(k1,k2,k3)
C Local variables
      integer isw0,isw1,i1,i1m,i2,i2m,i3,i3m
      isw0 = mod(isw,10)
      isw1 = mod(isw/10,10)

      if (isw0/2 .eq. 0)      call fftz3(c,n1,n2,n3,k1,k2,k3,1,0,isig)
      if (mod(isw0,2) .eq. 0) call fftz3(f,n1,n2,n3,k1,k2,k3,1,0,isig)

      do  10  i3 = 1, n3
      i3m = n3+2-i3
      if (i3 .eq. 1) i3m = 1

        do  20  i2 = 1, n2
        i2m = n2+2-i2
        if (i2 .eq. 1) i2m = 1

          do  30  i1 = 1, n1
          i1m = n1+2-i1
          if (i1 .eq. 1) i1m = 1

            c(i1,i2,i3) = c(i1,i2,i3)*f(i1m,i2m,i3m)
   30     continue
   20   continue
   10 continue

      if (isw1/2 .eq. 0)      call fftz3(c,n1,n2,n3,k1,k2,k3,1,0,-isig)
      if (mod(isw1,2) .eq. 0) call fftz3(f,n1,n2,n3,k1,k2,k3,1,0,-isig)
      end

