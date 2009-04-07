      subroutine fermi(qval,dos,ndos,emin,emax,nsp,eferm,e1,e2,dosef)
C- Makes fermi energy from integrated density
C ----------------------------------------------------------------------
Ci Inputs
Ci   qval, number of electrons to fermi level; dos(i) integrated
Ci   density at bin i; ndos, number of bins + 1; emin, emax, energy
Ci   window.
Co Outputs
Co   Eferm, Fermi energy; e1, e2, confidence limits on Fermi energy
Co   i.e., Fermi energy lies between e1 and e2.
Co   dosef:  density of states at fermi level
Cr Remarks
Cr   emin and e1 (and emax and e2) may point to the same address.
Cr   This version uses idos decomposed into spin up, down for nsp=2
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer ndos,nsp
      double precision qval,dos(ndos,1),emin,emax,eferm,e1,e2,dosef
C Local parameters
      integer i1,ie,i2
      double precision de,q1,q2,d1mach,wt
C External procedures

C --- Check bounds of DOS ---
      wt = 1d0/2
      i2 = 1
      if (nsp .eq. 2) then
        wt = 1
        i2 = 2
      endif
      q1 = (dos(1,1)+dos(1,i2))*wt
      q2 = (dos(ndos,1)+dos(ndos,i2))*wt
      if (q1.gt.qval .or. q2.lt.qval) call fexit3(-1,111,
     .  ' Exit -1 FERMI: NOS ( %1,6;6g %1,6;6g ) does not '//
     .  'encompass Q = %1;6g',q1,q2,qval)

C --- Find bin that boxes E_f ---
      de = (emax-emin)/(ndos-1)
      i1 = 1
      q1 = (qval + d1mach(3))/wt
      do  1  ie = 2, ndos
        if (dos(ie,1)+dos(ie,i2) .gt. q1) goto 2
        i1 = ie
    1 continue
      call rx('bug in FERMI')
    2 continue

C --- Linear interpolation for the Fermi level ---
      q1 = (dos(i1,1)+dos(i1,i2))*wt
      q2 = (dos(i1+1,1)+dos(i1+1,i2))*wt
      e1 = emin + de*(i1-1)
      e2 = e1 + de
      eferm = e1 + (qval-q1)/(q2-q1)*de
      dosef = (q2-q1)/de

      end

