      subroutine intnos(ndos,dos,emin,emax,qval,efermi,dosef,eband)
C- Finds E_F from a tabulated number of states function
C-----------------------------------------------------------------------
Ci  Input
Ci    ndos : number of tabulated points; dos : integrated DOS
Ci    emin, emax : energy range of tabulation;
Ci    qval : number of valence electrons
Co  Output
Co    efermi : Fermi energy, dosef : DOS at E_f
C-----------------------------------------------------------------------
C     implicit none
      integer ndos
      double precision dos(0:ndos-1),emin,emax,qval,efermi,eband,dosef
      integer i,meshpt,iprint,lgunit
      double precision step,sum,q,q1,q2,e1,eps,d1mach
      eps = d1mach(3)
C --- make Fermi energy ---
      step = (emax - emin) / (ndos - 1)
      meshpt = 0
      q = qval + eps
      do  1  i = 1, ndos-1
        if ( dos(i) .ge. q ) goto 2
        meshpt = i
    1 continue
    2 continue
      if (meshpt .eq. ndos-1)
     .call rx('INTNOS : Fermi energy lies above emax')
C E_F lies between mesh points meshpt and meshpt+1 -- interpolate :
      q1 = dos(meshpt)
      q2 = dos(meshpt+1)
      e1 = emin + step * meshpt
      efermi = e1 + ( qval-q1 ) / ( q2-q1 ) * step
      dosef = (q2 - dos(meshpt-1)) / (2*step)
C --- make band energy by partial integration ---
      sum = .5d0 * q1
      do  3  i = 1, meshpt-1
        sum = sum + dos(i)
    3 continue
      sum = sum * step
      sum = sum + .5d0 * (efermi - e1) * (qval + q1)
      eband = efermi * qval - sum
      if (iprint() .ge. 30) then
c        do  12  i = 1, 1
         write(lgunit(1),"(' INTNOS: Fermi energy=,f15.8',
     .    '  band energy=',f15.8,' DOS(E_f)=',f15.8)") efermi, eband, dosef
c   12   continue
C        write(*,10) efermi, eband, dosef
C        write(fopn('LOG'),10) efermi, eband, dosef
C   10 format(' INTNOS: Fermi energy =',f10.6,'; band energy =',f11.6/
C     .       '          DOS(E_f) =',f11.6)
      endif
      end

