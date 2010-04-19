      subroutine meptr(i,j,iam,npm,k)
C- Returns index to matrix element table, if one exists
C ----------------------------------------------------------------
Ci Inputs
Ci   i,j: two classes for which matrix element index is sought
Ci   iam,npm: see remarks
Co Outputs
Co   k:  index, see remarks
Cr Remarks
Cr   iam(1,kk) and iam(2,kk) are the two classes for the kkth ME pair
Cr   k = iam(3,kk) is a pointer to tabme(*,k) which holds the
Cr   matrix elements <iam(1,kk) | H | iam(2,kk)>
Cr   npm(0,i) is number matrix elts tabulated for ith class
Cr   npm(1,i) is number matrix elts in all sites preceding ith class
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer i,j,iam(3,1),npm(0:1,1),k,iprint
C Local variables
      integer l, lmax

      l = npm(1,i)
      k = 0
      lmax = l+npm(0,i)
   10 l = l+1
      if (l .gt. lmax) then
        if (iprint() .ge. 20) print 333, i,j
  333   format(' MEPTR: ***WARNING*** no melts recorded for i,j=',2i3)
        return
      endif
      if (iam(2,l) .ne. j) goto 10
      k = iam(3,l)
      end
