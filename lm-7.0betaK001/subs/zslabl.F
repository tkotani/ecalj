      subroutine zslabl(iopt,slabl,iz)
C- Matches element label with atomic number
C ----------------------------------------------------------------------
Ci Inputs:
Ci   iopt  :-1  iz   =input  ; slabl=output
Ci           1  slabl=input  ; iz   =output
Ci Inputs/Outputs:
Cio  slabl :name of chemical formula.
Cio  iz    :nuclear charge.
Cio        :If slabl is input but does not match any known unit,
Cio        :iz is returned -1
Cu Updates
Cu   02 Nov 01 Adapted from Stuttgart lmto56 zclabl
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer iopt,iz
      character*(*) slabl
C Local variables:
      integer kz
      character*2 aslabl(0:100), ccopy*2
C External calls:
      external  chcase
C Intrinsic functions:
      intrinsic  ichar
C Data statements:
      data aslabl/'E ',
     .            'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne',
     .            'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca',
     .            'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn',
     .            'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr',
     .            'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     .            'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd',
     .            'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     .            'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg',
     .            'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
     .            'Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm'/


C --- Return label corresponding to supplied atomic number ---
      if (iopt .eq. -1) then
        if (iz.lt.0 .or. iz.gt.100)
     .    call rxi('ZSLABL: bad atomic number',iz)
        slabl = aslabl(iz)

C --- Find atomic number from label ---
      elseif (iopt .eq. 1) then
        ccopy = slabl
        call  chcase( 1,1,ccopy(1:1))
        call  chcase(-1,1,ccopy(2:2))
        if (ccopy(2:2).ge.'0' .and. ccopy(2:2).le.'9') ccopy(2:2) = ' '
        do  kz = 0, 100
          iz = kz
          if (aslabl(iz) .eq. ccopy(1:2)) return
        enddo
C       write(lgunit(2),401) slabl
        iz = -1

      else
        call rxi(' ZSLABL: bad iopt',iopt)
      endif

C  400 format(/' ZSLABL: bad Z=',i5,' set to -1')
C  401 format(/' ZSLABL: could not find atom ',a4,' in my list,',
C     .       ' Z set to -1')
      end

