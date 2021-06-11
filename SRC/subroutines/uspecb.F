      module m_uspecb
      integer,allocatable,protected:: nkapii(:)
      
      contains
      subroutine uspecb(lpack,mode,sspec,is1,is2,lh,rsmh,eh,nkape)
      use m_struc_def  !Cgetarg
      use m_lmfinit,only: nkaph,nspec,nkapii_=>nkapii,nkaphh,
     &     rsmhss=>rsmh,rsmh2ss=>rsmh2,ehss=>eh,eh2ss=>eh2,lhh,lpzex !,rsml,ehl
       use m_elocp,only: rsml,ehl
C- Pack/unpack parameters related to basis
C ----------------------------------------------------------------------
Ci Inputs
Ci   lpack :0 for unpack, 1 for pack
Ci         :NB packing case, it is the caller's responsibility
Ci         :to ensure rsmh and eh are nonzero only for l and energies
Ci         :1..nkape which have envelope functions.
Ci   mode  :At present, this is the only implemented mode
Ci         :-1 Only forces uspecb to read nkaph from global variable
Ci         :0 No reference to envelope function parameters.
Ci         :  Only nkape is returned; see below
Ci         :1,2,3: parameters consist of rsmh,eh, where:
Ci              rsmh,eh         first orbital in basis
Ci              rsmh2,eh2       if nkap=2
Ci              rsmhl,ehl       if extended local orbitals
Ci            Data taken from spec->orbp (see Remarks)
Ci            modes 1..4 are similar except that:
Ci            mode=1 does not return envelope information about
Ci                   local orbitals
Ci            mode=2 returns rsmh,eh in last channel corresponding to
Ci                   loc. orb. (rsmh>0 for extended loc. orbs; else 0)
Ci                   and sets nkape to nkaph
Ci            xxxmode=3 intended for local orbitals of the second type
Ci                   NOT IMPLEMENTED
Ci            mode=4 only returns envelope information about
Ci                   local orbitals in rsmh(*,1),eh(*,1)
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: lmxb pz orbp
Ci     Stored:    name orbp
Ci     Passed to: spacks
Ci  is1,is2:range of species for which to extract parameters
Cio Inputs/Outputs
Co   ... All of the following are generated for lpack=0
Co   lh    :lh(1..nkape) = max l for which envelope basis fn defined
Co         :Also lh(nkaph) = max l for which a local orbital defined
Co         :lh(1+nkape..nkap0) are initialized to -1
Cio  rsmh  :rsmh(1..l+1,1..nkape,is1..is2), rsmh(1..l+1,nkaph,is1..is2)
Cio        :smoothing radii for envelope functions, extended loc. orbitals
Cio        :(rsmh is input for lpack=1)
Cio  eh    :eh(1..l+1,1..nkape,is1..is2), eh(1..l+1,nkaph,is1..is2)
Cio        :energies of envelope functions and extended local orbitals
Cio        :(eh is input for lpack=1)
Cio  nkape :number of envelope function types per l q.n. for is1..is2
Cio        :1s digit contains number of energy channels in rsmh,eh
Cio        :which contain valence orbitals with envelope functions.
Cio        :nkape may be modified depending on mode:
Cio        :mode=1: nkape not affect by local orbitals
Cio        :mode=2: 1 added to nkape if local orbitals are included.
Co         :When called with mode=0, nkape returns the following:
Co         :1s    digit always zero
Co         :10s   digit 1 if local orbital of first type exists
Co         :100s  digit 1 if local orbital of 2nd type exists +10 mode
Cl Local variables
Cl   npzi  :1 if there are extended local orbitals AND mode>1
Cl         :0 otherwise
Cl   nkapi :number of envelope functions EXCLUDING
Cl         :extended local orbitals
Cr Remarks
Cr   Case parameters contained in sspec->orbp
Cr     orbp is dimensioned orbp(n0,2,nkap0)
Cr     1st index for orbital types ?
Cr     2nd index :1 for rsmh, 2 for eh
Cr     3rd index :1 for (rsmh,eh) 2 for (rsmh2,eh2), 3 loc orb.
Cr     Or, if no rsmh2,eh2, loc. orb occupy 2nd spot
Cu Updates
Cu    8 Jul 05 New mode=-1
Cu   25 Jul 04 Redesigned for extended local orbitals
Cu   10 Apr 02 Redimensioned eh,rsmh to accomodate larger lmax
Cu   19 Feb 02 mode 0 implemented.  Convention for npqni changed.
Cu   12 Feb 02 npqn,nkape are returned for range is1..is2
Cu   25 Aug 01 Extended to return npqn,lh for local orbitals
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer nkap0,n0,nkape,i_copy_size
      parameter (nkap0=3,n0=10)
      integer lpack,mode,is1,is2,lh(nkap0,is1:is2)
      real(8):: rsmh(n0,nkap0,is1:is2) , eh(n0,nkap0,is1:is2)
      type(s_spec)::sspec(*)
      character spid*8
      integer is,lmxb,ikap,kkap,nkapi,i,j,k,l,npzi,lpzi,lpzl,npze,lmxa
      double precision orbp(n0,2,nkap0),dasum,pnz(n0,2)
      logical:: init=.true.
      if(init) then
         allocate(nkapii(nspec))
         init=.false.
      endif
       
c$$$      if(lpack/=0) then
c$$$          do  is = is1, is2
c$$$            nkapi = mod(nkape,10)
c$$$            npzi = 0
c$$$            if (nkape .gt. 10) npzi = 1
c$$$            kkap = 0
c$$$            do  ikap = 1, nkaph
c$$$              if (ikap.le.nkapi .or. ikap.eq.nkaph .and. npzi.ne.0) then
c$$$                kkap = kkap+1
c$$$                call dcopy(n0,rsmh(1,kkap,is),1,orbp(1,1,ikap),1)
c$$$                call dcopy(n0,eh(1,kkap,is),  1,orbp(1,2,ikap),1)
c$$$              endif
c$$$            enddo
c$$$            i_copy_size=size(sspec(is)%orbp)
c$$$            call dcopy(i_copy_size,orbp,1,sspec(is)%orbp,1)
c$$$        enddo
c$$$      endif
      
      if(lpack==0) then
          nkape = 0
          npze = 0
          eh  = 0d0
          rsmh= 0d0
          is=is1
          nkape     = nkaphh(is)
          nkapii(is)= nkapii_(is)
          lh(:,is)  = lhh(:,is)
          rsmh(:,1,is) = rsmhss(:,is)
          eh(:,1,is) = ehss(:,is)
          if(nkapii(is)==2) rsmh(:,2,is)=rsmh2ss(:,is)
          if(nkapii(is)==2) eh(:,2,is)=eh2ss(:,is)
          if(lpzex(is)==1.and.allocated(rsml))  rsmh(:,nkaphh(is),is)=rsml(:,is)
          if(lpzex(is)==1.and.allocated(ehl))  eh(:,nkaphh(is),is)=ehl(:,is)
c          print *,"qqqqqqq",lh(:,is)
c          print *,"qqqqqqq",rsmh(:,1,is)
c          print *,"qqqqqqq",eh(:,1,is)
c          print *,"qqqqqqq",rsmh(:,2,is)
c          print *,"qqqqqqq",eh(:,2,is)
c          print *,"qqqqqqq",mode,is,nkapii(is),nkape
c          return
      endif   
ccccccccccccccccccccccccccccccccccccccccccccccc
      end subroutine uspecb
      end module
