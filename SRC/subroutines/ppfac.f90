subroutine ppfac(npfac,pfac,fmax,job,fac,nfac)
  !- Find all products of prime factors within some maximum
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   npfac :number of prime factors
  !i   pfac  :list of prime factors
  !i   fmax  :maximum value of product
  !i   job   :0 list is returned unsorted
  !i         :1 list is returned sorted
  !o Outputs
  !o   fac   :all products of pfac <= fmax
  !o         :List is returned sorted.
  !o         :fac must be dimensioned at least nfac
  !o         :or 3*nfac if fac is returned sorted.
  !o   nfac  :length of fac
  !r Remarks
  !u Updates
  ! ----------------------------------------------------------------------
  implicit none
  integer :: npfac,pfac(npfac),fac(*),fmax,job,nfac
  integer :: i,m,k,fack,nfaci
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
           if (fack <= fmax) then
              nfac = nfac+1
              fac(nfac) = fack
           else
              goto 10
           endif
        enddo
10      continue
     enddo
  enddo
  if (job /= 0) call ivheap(1,nfac,fac,fac(1+2*nfac),0)
end subroutine ppfac
!     Testing
!      subroutine fmain
!      integer fmax,nfac,nprime
!      parameter (nprime=4,fmax=100)
!      integer fac(3*fmax),pfac(nprime)
!      data pfac /2,3,5,7/
!      call ppfac(nprime,pfac,fmax,1,fac,nfac)
!      print 333, (fac(i), i=1,nfac)
!  333 format(i4)
!      end

