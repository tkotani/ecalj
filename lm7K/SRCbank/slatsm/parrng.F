      subroutine parrng(n1,n2,nprocs,irank,ista,iend)
C- returns the working range for parallelization
C ----------------------------------------------------------------------
Ci Inputs
Ci   n1      :starting index of the entire range
Ci   n2      :ending index of the entire range
Ci   nprocs  :number of parallel processes
Ci   irank   :rank of the calling process
Co Outputs
Co   ista    :starting index for the calling process
Co   iend    :ending index for the calling process
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer n1,n2,nprocs,irank,ista,iend
C ... Local parameters
      integer iwork1,iwork2

C --- Calculating the working range ---
      iwork1 = (n2 - n1 + 1) / nprocs
      iwork2 = mod(n2 - n1 + 1, nprocs)
      ista = irank * iwork1 + n1 + min(irank, iwork2)
      iend = ista + iwork1 - 1
      if (iwork2 .gt. irank) iend = iend + 1
      end


