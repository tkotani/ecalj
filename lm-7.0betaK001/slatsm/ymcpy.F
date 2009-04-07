      subroutine ymcpy(a,nca,nra,offa,b,ncb,nrb,offb,n,m)
C- General matrix copy
C ----------------------------------------------------------------
Ci Inputs:
Ci   a,nca,nra is the left matrix and respectively the number of
Ci      elements separating columns and rows.
Ci   offa      is the offset to the imaginary part of a
Ci   b,ncb,nrb is the right matrix and respectively the number of
Ci      elements separating columns and rows.
Ci   offb      is the offset to the imaginary part of b
Ci   n,m: the number of columns and rows, respectively, to calculate
Co Outputs:
Co   result matrix stored in c
Cr Remarks:
Cr   This is a complex analog of dmcpy
C ----------------------------------------------------------------
C     implicit none
      integer nca,nra,ncb,nrb,n,m,offa,offb
      double precision a(0:*), b(0:*)

      call dmcpy(a,nca,nra,b,ncb,nrb,n,m)
      call dmcpy(a(offa),nca,nra,b(offb),ncb,nrb,n,m)

      end

