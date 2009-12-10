      subroutine prterr(ndim,nlm,nlm1,xi1,xi2,
     .  err,irerr,ilerr,eps,ireps,ileps)
C- Check the abs. and relative difference between arrays xi1 and xi2
C  xi1 is xi1-xi2 on the output
C ----------------------------------------------------------------
Ci Inputs
Ci   ndim  :leading dimension of xi1 and xi2
Ci   nlm   :first actual dimension of xi1 and xi2 
Ci   nlm1  :second actual dimension of xi1 and xi2 
Ci Inputs/Outputs
Cio  xi1,xi2 :arrays to compare. On output xi1 is xi1-xi2
Co Outputs
Co   err   :max abs error, |xi1(L,L')-xi2(L,L')|
Co   irerr,ilerr :max abs error occurs at L = irerr and L'=ilerr
Co   eps   :max rel error, |xi1(L,L')-xi2(L,L')|/|xi2(L,L')|
Co   ireps,ileps :max rel error occurs at L = ireps and L'=ileps
Cr Remarks
Cr   Whenever |xi2(L,L')|=0, eps is assumed to be zero  
Cu Updates
Cu   27 Feb 08 First written
C ----------------------------------------------------------------
C     implicit none
      integer ndim,nlm,nlm1
      double precision xi1(ndim,nlm1),xi2(ndim,nlm1)
      double precision err,eps,xx
      integer ilerr,irerr,ileps,ireps,ilm1,ilm

      call daxpy(ndim*nlm1,-1d0,xi2,1,xi1,1)
      err = -1d0
      eps = -1d0
      do ilm1 = 1, nlm1 
        do ilm = 1, nlm 
          if(dabs(xi1(ilm,ilm1)) .gt. err) then
             irerr = ilm
             ilerr = ilm1
             err = dabs(xi1(ilm,ilm1))
          endif
          xx = dabs(xi2(ilm,ilm1))
          if(xx .ne. 0d0 .and. dabs(xi1(ilm,ilm1))/xx .gt. eps) then
             ireps = ilm
             ileps = ilm1
             eps = dabs(xi1(ilm,ilm1)/xi2(ilm,ilm1))
          endif
        enddo
      enddo
c300  format(11x,a,g8.2,a,i2,a,i4,a,f9.6,a) 
c     write(6,300) 'Maximal difference = ',err,' at L'' =',ilerr,
c    .   ' and L =',irerr   
c     write(6,300) 'Maximal rel. error = ',eps,' at L'' =',ileps,
c    .   ' and L =',ireps   
c     if(err .eq. 0d0 .and. eps .eq. 0d0) then 
c        write(6,300) 'prterr: Arrays are identical!' 
c     endif

      end 

