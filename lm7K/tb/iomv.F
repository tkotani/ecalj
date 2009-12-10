      subroutine iomv(nbas,pos,alat,ifi)
C- Append positions to mv file
C ----------------------------------------------------------------------
Ci Inputs:
Ci   nbas,pos,ifi
Co Outputs:
Co   
Cr Remarks
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      integer nbas,ifi
      double precision pos(3,nbas),alat
C Local Variables
      integer ib,m
      write (ifi,1) ((alat*pos(m,ib),m=1,3),ib=1,nbas)
    1 format(8f10.5)
          
      end
