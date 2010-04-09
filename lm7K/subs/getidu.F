Cgetarg...info...           structure ['getidu', 'site', 'ssite']
Cgetarg...info...           structure ['getidu', 'spec', 'sspec']
Cgetarg...info...           use_to_add ['getidu', 'm_struc_def']
      subroutine getidu(nbas,sspec,ssite,idu)
       
       use m_struc_def  !Cgetarg

C- Get idu for all sites
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: lmxa idu
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: spec
Co Outputs
Co   idu   :idu(l+1,ib)=1 => this l and site has a nonlocal U matrix
Cr Remarks
Cr
Cu Updates
Cu   Mar 29, 07 (Jialei) first created
C------------------------------------
C     implicit none
C ... Passed parameters
      integer nbas,idu(4,nbas)
Cgetarg       double precision sspec(1), ssite(1)
       type(s_spec)::sspec(*)
       type(s_site)::ssite(*)

C ... Local parameters
      integer igetss,is,ib

      call iinit(idu,4*nbas)
      do  ib = 1, nbas
Cgetarg         is = igetss('site spec',ib,ssite)
         is = int(ssite(ib)%spec) 

Cgetarg...info...           integer :: idu(4,nbas) , integer(8) :: sspec%idu 4
Cgetarg         call upack('spec idu',sspec,is,idu(1,ib),0,0,0)
         
         i_copy_size=size(sspec(is)%idu) 
         call i8icopy(i_copy_size,sspec(is)%idu,1,idu(1,ib),1) 

      enddo

Cgetarg       end
       end subroutine getidu 





