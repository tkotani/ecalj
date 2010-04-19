      subroutine ioqm(nsp,nelts,neltsm,qmpol,mmom,ierr)
C- I/O multipole and magnetic moments
C ----------------------------------------------------------------------
Ci Inputs:
Ci   nsp,nelts,neltsm,qmpol,mmom,ierr
Co Outputs:
Co   ierr
Cr Remarks
Cr   for ierr=1, write moments to disc
Cr   for ierr=-1, attempt to read. If successful return ierr=0
C ----------------------------------------------------------------------
C     implicit none
C Passed Parameters
      integer nsp,nelts,neltsm,ierr
      double precision qmpol(nelts),mmom(neltsm)
C Local Variables
      integer ifi,fopn,iprint

      ifi = fopn('QMOM')
      rewind ifi
      if (ierr .eq. 1) then
        write (ifi) qmpol
        if (nsp .eq. 2) then
          write (ifi) mmom
        endif
      else
        read (ifi, end=1, err=1) qmpol
        if (nsp .eq. 2) then
          read (ifi, end=1, err=1) mmom
        endif
        ierr = 0
      endif
      return
    1 continue 
      ierr = -1
      if (iprint() .gt. 30) then
        print *, ' '
        print *, 'IOQM: error reading moments, using defaults ..'
      endif
      end
