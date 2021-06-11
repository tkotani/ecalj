      module m_uspecb
      integer,allocatable,protected:: nkapii(:)
      
      contains
      subroutine uspecb(is,lh,rsmh,eh,nkape) !lpack,mode,sspec,is1,is2,lh,rsmh,eh,nkape)
      use m_struc_def  
      use m_lmfinit,only: nkaph,nspec,nkapii_=>nkapii,nkaphh_=>nkaphh,
     &     rsmhss=>rsmh,rsmh2ss=>rsmh2,ehss=>eh,eh2ss=>eh2,lhh,lpzex,sspec=>v_sspec !,rsml,ehl
       use m_elocp,only: rsml,ehl
C output only -----------------------------------------------------------------
Ci         :1,2,3: parameters consist of rsmh,eh, where:
Ci              rsmh,eh         first orbital in basis
Ci              rsmh2,eh2       if nkap=2
Ci              rsmhl,ehl       if extended local orbitals
!! Outputs. We generate all quanitites in m_lmfinit and m_elocp
!   lh    :lh(1..nkape) = max l for which envelope basis fn defined
!         :lh(1+nkape..nkap0) are initialized to -1 (see m_lmfinit
Ci  rsmh  :rsmh(1..l+1,1..nkape,is1..is2), rsmh(1..l+1,nkaph,is1..is2)
Cio        :smoothing radii for envelope functions, extended loc. orbitals
Cio        :(rsmh is input for lpack=1)
Ci  eh    :eh(1..l+1,1..nkape,is1..is2), eh(1..l+1,nkaph,is1..is2)
Cio        :energies of envelope functions and extended local orbitals
Cio        :(eh is input for lpack=1)
!  nkape :=nkaphh(i)
!       number of envelope function types per l q.n. for is1..is2
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nkap0,n0,i_copy_size,nkape
      parameter (nkap0=3,n0=10)
      integer is,lh(nkap0,is:is)
      real(8):: rsmh(n0,nkap0,is:is) , eh(n0,nkap0,is:is)
c      type(s_spec)::sspec(*)
      character spid*8
      integer lmxb,ikap,kkap,nkapi,i,j,k,l,npzi,lpzi,lpzl,npze,lmxa
      double precision orbp(n0,2,nkap0),dasum,pnz(n0,2)
      logical:: init=.true.
      if(init) then
         allocate(nkapii(nspec))
         init=.false.
      endif
c      if(lpack==0) then
          npze = 0
          eh  = 0d0
          rsmh= 0d0
          nkapii(is)= nkapii_(is)
          nkape     = nkaphh_(is)
          lh(:,is)  = lhh(:,is)
          rsmh(:,1,is) = rsmhss(:,is)
          eh(:,1,is) = ehss(:,is)
          if(nkapii(is)==2) rsmh(:,2,is)=rsmh2ss(:,is)
          if(nkapii(is)==2) eh(:,2,is)=eh2ss(:,is)
          if(lpzex(is)==1.and.allocated(rsml)) rsmh(:,nkaph,is)=rsml(:,is)
          if(lpzex(is)==1.and.allocated(ehl))  eh(:,nkaph,is)=ehl(:,is)
c          print *,"qqqqqqq",lh(:,is)
c          print *,"qqqqqqq",rsmh(:,1,is)
c          print *,"qqqqqqq",eh(:,1,is)
c          print *,"qqqqqqq",rsmh(:,2,is)
c          print *,"qqqqqqq",eh(:,2,is)
c          print *,"qqqqqqq",mode,is,nkapii(is),nkape
c          return
c      endif   
      end subroutine uspecb
      end module
