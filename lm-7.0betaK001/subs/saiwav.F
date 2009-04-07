      subroutine saiwav(ldim,lidim,nd1,nd2,a)
C- Average of up- and down- channels of an array for iwaves
C ----------------------------------------------------------------
Ci Inputs
Ci   ldim  :dimension of lower block
Ci   lidim :dimension of lower+intermediate block
Ci   nd1   :leading dimension of array a
Ci   nd2   :second dimension of array a
Co Outputs
Co   a     :spin-up and down-down channels of (ldim+1..) are averaged
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer ldim,lidim,nd1,nd2
      double precision a(nd1,nd2,2)
C Local parameters
      integer i,k
      double precision xx

      do  10  i = ldim+1, lidim
        do  12  k = 1, nd1
          xx = (a(k,i,1) + a(k,i,2))/2
          a(k,i,1) = xx
          a(k,i,2) = xx
   12   continue
   10 continue
      end

