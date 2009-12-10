      subroutine csmaln(v,n,machep,lunit,iopt,nzero)
C- Sets to zero numerically small numbers in a vector
C----------------------------------------------------------------------
Ci Inputs
Ci    v,n:  array of length n
Ci    machep: machine precision
Ci    iopt: not used now
Co Outputs
Co    elements 1..n for which abs(v)<machep are set to zero.
Co    nzero: number of elements changed
Cr Remarks
C----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer n,iopt,nzero(2),lunit
      double precision v(n),machep
C Local variables:
      integer i

      nzero(1) = 0
      nzero(2) = 0
      do  10  i = 1, n
        if (dabs(v(i)) .gt. machep) goto 10
        if (v(i) .eq. 0d0) then
          nzero(1) = nzero(1)+1
        else
          nzero(2) = nzero(2)+1
          v(i) = 0
        endif
   10 continue

      if (lunit .ge. 0)
     .  write(lunit,333) nzero(1),nzero(2),nzero(1)+nzero(2),n
  333 format(' csmaln:',i6,' +',i6,' =',i6,'  small numbers of',i6)
      end

