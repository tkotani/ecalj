      subroutine shftpp(nc,nlsp,pp,vold,vnew,oshft,nshft)
C- Shift or undo shift of pp's by constant potential
C ----------------------------------------------------------------------
Ci Inputs:
Ci   nc,pp,ves:
Cr Remarks:   shift pp's as follows:
Cr     oshft\nshft   F                 T
Cr           ---------------------------------------
Cr       F  |   do nothing       shift by vnew
Cr       T  | shift by -vold   shift by vnew-vold
Cr           ---------------------------------------
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      logical oshft,nshft
      integer nc,nlsp
      double precision pp(6,nlsp,nc),vold(nc),vnew(nc)
C Local variables
      integer ic

C --- Shift enu and c for each class ---
      do  10  ic = 1, nc
        if (oshft) then
          call daxpy(nlsp,-1d0,vold(ic),0,pp(1,1,ic),6)
          call daxpy(nlsp,-1d0,vold(ic),0,pp(2,1,ic),6)
        endif

        if (nshft) then
          call daxpy(nlsp,1d0,vnew(ic),0,pp(1,1,ic),6)
          call daxpy(nlsp,1d0,vnew(ic),0,pp(2,1,ic),6)
        endif
   10 continue

      end





