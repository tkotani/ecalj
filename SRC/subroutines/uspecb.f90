      module m_uspecb
      contains
      subroutine uspecb(is,rsmh,eh)
      use m_lmfinit,only: nkaph,nspec,nkapii,n0,nkap0,
     &     rsmh1ss=>rsmh1,rsmh2ss=>rsmh2,eh1ss=>eh1,eh2ss=>eh2,lhh,lpzex,sspec=>v_sspec 
      use m_elocp,only: rsml,ehl
      implicit none
      intent(in)::     is
      intent(out)::        rsmh,eh
C output only -----------------------------------------------------------------
Ci         :1,2,3: parameters consist of rsmh,eh, where:
Ci         rsmh,eh         first orbital in basis
Ci         rsmh2,eh2       if nkap=2
Ci         rsmhl,ehl       if extended local orbitals
!NOTE: nkaph is index for local orbital globally. 
! Followings are old data handled in this routine (obsolate) 
!  lh    :lh(1..nkape) = max l for which envelope basis fn defined
!        :lh(1+nkape..nkap0) are initialized to -1 (see m_lmfinit
!  rsmh  :rsmh(1..l+1,1..nkape,is1..is2), rsmh(1..l+1,nkaph,is1..is2)
!        :smoothing radii for envelope functions, extended loc. orbitals
!        :(rsmh is input for lpack=1)
!  eh    :eh(1..l+1,1..nkape,is1..is2), eh(1..l+1,nkaph,is1..is2)
!        :energies of envelope functions and extended local orbitals
C ----------------------------------------------------------------------
      integer is 
      real(8):: rsmh(n0,nkap0,is:is) , eh(n0,nkap0,is:is)
      eh  = 0d0
      rsmh= 0d0
      rsmh(:,1,is) = rsmh1ss(:,is)
      eh(:,1,is)   = eh1ss(:,is)
      if(nkapii(is)==2) rsmh(:,2,is)=rsmh2ss(:,is)
      if(nkapii(is)==2) eh(:,2,is)=eh2ss(:,is)
      if(lpzex(is)==1.and.allocated(rsml)) rsmh(:,nkaph,is)=rsml(:,is)
      if(lpzex(is)==1.and.allocated(ehl))  eh(:,nkaph,is)=ehl(:,is)
      end subroutine uspecb
      end module
