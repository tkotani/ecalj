Cgetarg...info...           structure ['getzv', 'bz', 'sbz']
Cgetarg...info...           structure ['getzv', 'ctrl', 'sctrl']
Cgetarg...info...           use_to_add ['getzv', 'm_struc_def']
      subroutine getzv(nclass,nrc,z,qc,sbz,sctrl,zval)
       
       use m_struc_def  !Cgetarg

C- Calculate and save number of valence electrons
C     implicit none
      integer nclass
Cgetarg       double precision qc(1),sbz(1),sctrl(1),z(nclass)
       real(8):: qc(1) , z(nclass) 
       type(s_bz)::sbz
       type(s_ctrl)::sctrl

      integer nrc(1)
C Local variables
      double precision zval,dgets,zbak(2)
      integer ic,lgunit,iprint,isw
      logical lgors

Cgetarg       zval = dgets('bz zval',sbz)
       zval = sbz%zval 

      if (zval .eq. 0) then
        do  10  ic = 1, nclass
   10   zval = zval + nrc(ic) * (z(ic)-qc(ic))
      endif
Cgetarg...info...           real(8) :: zval , real(8) :: sbz%zval 1
Cgetarg       call pack1('bz zval',sbz,zval)
       
       sbz%zval=zval 


C --- Printout ---
Cgetarg...info...           real(8) :: zbak(2) , real(8) :: sctrl%zbak 2
Cgetarg       call upack1('ctrl zbak',sctrl,zbak)
       
       i_copy_size=size(sctrl%zbak) 
       call dcopy(i_copy_size,sctrl%zbak,1,zbak,1) 

      if (iprint() .ge. 30) then
Cgetarg         call awrit4('%N GETZV:  %d valence electrons'//
Cgetarg      .    '%?#n#  zbak=%d  qbak=%d##',' ',80,lgunit(1),zval,
Cgetarg      .    isw(lgors('ctrl lasa,64 ',sctrl)),zbak,zbak(2))
         call awrit4 ( '%N GETZV:  %d valence electrons' // '%?#n#  zbak=%d  qbak=%d##' 
     .   , ' ' , 80 , lgunit ( 1 ) , zval , isw ( iand(64,int(sctrl%lasa)) .ne.0 
     .   ) , zbak , zbak ( 2 ) ) 

      endif
      if (iprint() .ge. 20) then
Cgetarg         call awrit4('     zval %d'//'%?#n#  zbak %;2d  qbak %;2d##',
Cgetarg      .    ' ',80,lgunit(2),zval,isw(lgors('ctrl lasa,64 ',sctrl)),
Cgetarg      .    zbak,zbak(2))
         call awrit4 ( '     zval %d' // '%?#n#  zbak %;2d  qbak %;2d##' 
     .   , ' ' , 80 , lgunit ( 2 ) , zval , isw ( iand(64,int(sctrl%lasa)) .ne.0 
     .   ) , zbak , zbak ( 2 ) ) 

      endif

Cgetarg       end
       end subroutine getzv 


