      subroutine dmx22(a,nca,nra,b,ncb,nrb,c,ncc,nrc)
C- General multiplication of 2x2 matrices a*b
C ----------------------------------------------------------------------
Ci Inputs
Ci   a,nca,nra is the left matrix and respectively the spacing
Ci      between elements in adjacent columns and rows.
Ci      Typically nca=2 and nra=1
Ci   b,ncb,nrb is the right matrix and respectively the spacing
Ci      between elements in adjacent columns and rows.
Ci      Typically ncb=2 and nrb=1
Co Outputs
Co   c,ncc,nrc is the product matrix and respectively the spacing
Co      between elements in adjacent columns and rows.
Cr Remarks
Cr   It is permissible for any of a,b,c to use the same address space
Cu Updates
Cu   17 Mar 03 First created (from A. Chantis)
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nca,nra,ncb,nrb,ncc,nrc
C      double precision a(2,2), b(2,2), c(2,2)
      double precision a(0:*), b(0:*), c(0:*)
C ... Local parameters
      double precision cloc(2,2)

C      cloc(1,1) = a(1,1)*b(1,1) + a(1,2)*b(2,1)
C      cloc(1,2) = a(1,1)*b(1,2) + a(1,2)*b(2,2)
C      cloc(2,1) = a(2,1)*b(1,1) + a(2,2)*b(2,1)
C      cloc(2,2) = a(2,1)*b(1,2) + a(2,2)*b(2,2)

      cloc(1,1) = a(nra*0+nca*0)*b(nrb*0+ncb*0) +
     .a(nra*0+nca*1)*b(nrb*1+ncb*0)
      cloc(1,2) = a(nra*0+nca*0)*b(nrb*0+ncb*1) +
     .a(nra*0+nca*1)*b(nrb*1+ncb*1)
      cloc(2,1) = a(nra*1+nca*0)*b(nrb*0+ncb*0) +
     .a(nra*1+nca*1)*b(nrb*1+ncb*0)
      cloc(2,2) = a(nra*1+nca*0)*b(nrb*0+ncb*1) +
     .a(nra*1+nca*1)*b(nrb*1+ncb*1)

      c(nrc*0+ncc*0) = cloc(1,1)
      c(nrc*0+ncc*1) = cloc(1,2)
      c(nrc*1+ncc*0) = cloc(2,1)
      c(nrc*1+ncc*1) = cloc(2,2)

      end

C     Test
C      subroutine fmain
C
C      double precision a(2,2), b(2,2), c(2,2)
C      double precision aa(3,2,4,2),bb(2,2,3,2),cc(5,2,7,2)
C
C      a(1,1) = 1d0
C      a(1,2) = 2d0
C      a(2,1) = 4d0
C      a(2,2) = 3d0
C
C      b(1,1) = -1d0
C      b(1,2) = 2d0
C      b(2,1) = -4d0
C      b(2,2) = 1d0
C
CC     call prmx('a',a,2,2,2)
CC     call prmx('b',b,2,2,2)
C
CC     call dmpy22(a,b,c)
CC     call prmx('c',c,2,2,2)
C
C      aa(1,1,1,1) = a(1,1)
C      aa(1,1,1,2) = a(1,2)
C      aa(1,2,1,1) = a(2,1)
C      aa(1,2,1,2) = a(2,2)
C
C      bb(1,1,1,1) = b(1,1)
C      bb(1,1,1,2) = b(1,2)
C      bb(1,2,1,1) = b(2,1)
C      bb(1,2,1,2) = b(2,2)
C
C      call dmx22(aa,3*2*4,3,bb,2*2*3,2,cc,5*2*7,5)
C
C      print *, 'the following should be zero'
C      print *, cc(1,1,1,1) - a(1,1)*b(1,1) - a(1,2)*b(2,1)
C      print *, cc(1,1,1,2) - a(1,1)*b(1,2) - a(1,2)*b(2,2)
C      print *, cc(1,2,1,1) - a(2,1)*b(1,1) - a(2,2)*b(2,1)
C      print *, cc(1,2,1,2) - a(2,1)*b(1,2) - a(2,2)*b(2,2)
C
CC     call prmx('c',c,2,2,2)
C
C      end

