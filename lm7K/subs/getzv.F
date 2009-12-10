      subroutine getzv(nclass,nrc,z,qc,sbz,sctrl,zval)
C- Calculate and save number of valence electrons
C     implicit none
      integer nclass
      double precision qc(1),sbz(1),sctrl(1),z(nclass)
      integer nrc(1)
C Local variables
      double precision zval,dgets,zbak(2)
      integer ic,lgunit,iprint,isw
      logical lgors

      zval = dgets('bz zval',sbz)
      if (zval .eq. 0) then
        do  10  ic = 1, nclass
   10   zval = zval + nrc(ic) * (z(ic)-qc(ic))
      endif
      call pack1('bz zval',sbz,zval)

C --- Printout ---
      call upack1('ctrl zbak',sctrl,zbak)
      if (iprint() .ge. 30) then
        call awrit4('%N GETZV:  %d valence electrons'//
     .    '%?#n#  zbak=%d  qbak=%d##',' ',80,lgunit(1),zval,
     .    isw(lgors('ctrl lasa,64 ',sctrl)),zbak,zbak(2))
      endif
      if (iprint() .ge. 20) then
        call awrit4('     zval %d'//'%?#n#  zbak %;2d  qbak %;2d##',
     .    ' ',80,lgunit(2),zval,isw(lgors('ctrl lasa,64 ',sctrl)),
     .    zbak,zbak(2))
      endif

      end

