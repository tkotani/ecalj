c takao modified this based on delw_2.2.
      subroutine fftz30(n1,n2,n3,k1,k2,k3)
      integer:: n1,n2,n3,k1,k2,k3
      k1=n1
      k2=n2
      k3=n3
      end

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

