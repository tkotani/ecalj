      subroutine ppfac(npfac,pfac,fmax,job,fac,nfac)
C- Find all products of prime factors within some maximum
C ----------------------------------------------------------------------
Ci Inputs
Ci   npfac :number of prime factors
Ci   pfac  :list of prime factors
Ci   fmax  :maximum value of product
Ci   job   :0 list is returned unsorted
Ci         :1 list is returned sorted
Co Outputs
Co   fac   :all products of pfac <= fmax
Co         :List is returned sorted.
Co         :fac must be dimensioned at least nfac
Co         :or 3*nfac if fac is returned sorted.
Co   nfac  :length of fac
Cr Remarks
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer npfac,pfac(npfac),fac(*),fmax,job,nfac
C ... Local parameters
      integer i,m,k,fack,nfaci

      nfac = 0
      nfaci = 0
      do  i = 1, npfac
        nfac = nfac+1
        nfaci = nfac
        fac(nfac) = pfac(i)
        do  k  = 1, nfaci
        fack = fac(k)
        do  m = 1, fmax
          fack = fack*pfac(i)
          if (fack .le. fmax) then
            nfac = nfac+1
            fac(nfac) = fack
          else
            goto 10
          endif
        enddo
   10   continue
        enddo
      enddo

      if (job .ne. 0) call ivheap(1,nfac,fac,fac(1+2*nfac),0)
      end
C     Testing
C      subroutine fmain
C      integer fmax,nfac,nprime
C      parameter (nprime=4,fmax=100)
C      integer fac(3*fmax),pfac(nprime)
C      data pfac /2,3,5,7/
C      call ppfac(nprime,pfac,fmax,1,fac,nfac)
C      print 333, (fac(i), i=1,nfac)
C  333 format(i4)
C      end

