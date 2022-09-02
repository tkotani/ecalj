subroutine fftz30(n1,n2,n3,k1,k2,k3)
  integer:: n1,n2,n3,k1,k2,k3
  k1=n1
  k2=n2
  k3=n3
end subroutine fftz30

subroutine fftz3(c,n1,n2,n3,k1,k2,k3,nfft,iset,isig)
  implicit none
  integer :: n1,n2,n3,k1,k2,k3,nfft,iset,isig
  complex(8):: c(k1,k2,k3,nfft)
  integer :: i1,i2,i3,id,iord,iopt,ow1,ow2,ow3,oiwk,ierr,ifft
  real(8):: scale
  save ow1,ow2,ow3,oiwk
  ! --- A public-domain fft package.  See http://www.fftw.org/ ---
  ! ... Start of include file fftw_f77.i that comes with the fftw package
  !     This file contains PARAMETER statements for various constants
  !     that can be passed to FFTW routines.  You should include
  !     this file in any FORTRAN program that calls the fftw_f77
  !     routines (either directly or with an #include statement
  !     if you use the C preprocessor).
  integer :: FFTW_FORWARD,FFTW_BACKWARD
  parameter (FFTW_FORWARD=-1,FFTW_BACKWARD=1)
  integer :: FFTW_REAL_TO_COMPLEX,FFTW_COMPLEX_TO_REAL
  parameter (FFTW_REAL_TO_COMPLEX=-1,FFTW_COMPLEX_TO_REAL=1)
  integer :: FFTW_THREADSAFE
  parameter (FFTW_THREADSAFE=128)
  ! ... End of include file fftw_f77.i that comes with the fftw package
  integer :: jsig
  integer(8):: plan
  integer :: FFTW_ESTIMATE,FFTW_MEASURE
  INTEGER :: FFTW_PRESERVE_INPUT
  integer :: FFTW_OUT_OF_PLACE,FFTW_IN_PLACE,FFTW_USE_WISDOM
  parameter (FFTW_ESTIMATE=64,FFTW_MEASURE=0)
  parameter (FFTW_PRESERVE_INPUT=16)
  call tcn('fftz3')
  if (n1 == 1 .AND. n2 == 1 .AND. n3 == 1) goto 99
  jsig = 0
  if (isig == -1) jsig = FFTW_FORWARD
  if (isig ==  1) jsig = FFTW_BACKWARD
  do  10  ifft = 1, nfft
     if (n2 == 1 .AND. n3 == 1) then
        call dfftw_plan_dft_1d(plan,n1,c(1,1,1,ifft),c(1,1,1,ifft),      jsig,FFTW_ESTIMATE)
     elseif (n3 == 1) then
        call dfftw_plan_dft_2d(plan,n1,n2,c(1,1,1,ifft),c(1,1,1,ifft),   jsig,FFTW_ESTIMATE)
     else
        call dfftw_plan_dft_3d(plan,n1,n2,n3,c(1,1,1,ifft),c(1,1,1,ifft),jsig,FFTW_ESTIMATE)
     endif
     call dfftw_execute(plan,c(1,1,1,ifft),c(1,1,1,ifft))
     call dfftw_destroy_plan(plan)
10 enddo
  ! ... Renormalize forward transform
  if (isig > 0 .OR. isig < -10) goto 99
  scale = 1/dble(n1*n2*n3)
  do    ifft = 1, nfft
     do   i3 = 1, n3
        do  i2 = 1, n2
           do  i1 = 1, n1
              c(i1,i2,i3,ifft) = scale*c(i1,i2,i3,ifft)
           enddo
        enddo
     enddo
  enddo
99 continue
  call tcx('fftz3')
  return
end subroutine fftz3

