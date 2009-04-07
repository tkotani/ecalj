      subroutine covsrt(nvar,npar,ivar,cov)
C- Repack the covariance matrix to the true order
C ----------------------------------------------------------------------
Ci Inputs
Ci   nvar: number of parameters to vary out of npar total parameters
Ci   npar: total number of parameters (including fixed parameters)
Ci   ivar(1,i): points to the position in the full list of ith variable
Co Outputs
Co   cov: covariance matrix, repacked to true order on output
Cr Remarks
Cr   Elements of cov associated with fixed parameters set to zero.
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nvar,npar
      integer ivar(2,nvar)
      double precision cov(npar,npar)
C Local parameters
      integer i,j
      double precision swap

C --- Zero all elements below diagonal ---
      do  20  j = 1, npar-1
        do  10  i = j+1, npar
          cov(i,j) = 0d0
   10   continue
   20 continue

C --- Repack off-diagonal elements of fit into correct locations ---
      do  40  i = 1, nvar-1
        do  30  j = i+1, nvar
          if (ivar(1,j) .gt. ivar(1,i)) then
            cov(ivar(1,j),ivar(1,i)) = cov(i,j)
          else
            cov(ivar(1,i),ivar(1,j)) = cov(i,j)
          endif
   30   continue
   40 continue

C --- Temporarily store original diagonal elements in top row ---
      swap = cov(1,1)
      do  50  j = 1, npar
        cov(1,j) = cov(j,j)
        cov(j,j) = 0d0
   50 continue
      cov(ivar(1,1),ivar(1,1)) = swap

C --- Sort elements into proper order on diagonal ---
      do  60  j = 2, nvar
        cov(ivar(1,j),ivar(1,j)) = cov(1,j)
   60 continue

C --- Fill in above diagonal by symmetry ---
      do  80  j = 2, npar
        do  70  i = 1, j-1
          cov(i,j) = cov(j,i)
   70   continue
   80 continue

      end
