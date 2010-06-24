Cchp1 cmdvar= real(8),rv_p_,ogv
      subroutine suham2(sctrl,slat,sspec,ssite,sham,spot,sstrn)

      use m_struc_def  !Cgetarg

C- Further hamiltonian setup after potential generated
C ----------------------------------------------------------------------
Ci Inputs
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read: lfp
Ci     Stored:
Ci     Passed to:
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: ng ogv alat tolft
Ci     Stored:
Ci     Passed to:
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read:
Ci     Stored:
Ci     Passed to: sugcut
Ci   ssite :struct for site-specific information; see routine usite
Ci   sham  :struct for parameters defining hamiltonian; see routine uham
Ci   spot  :struct for information about the potential; see routine upot
Ci   sstrn :struct for global strings
Co Outputs
Cl Local variables
Cr Remarks
Cr   This routine completes energy-independent hamiltonian setup,
Cr   which requires information not when suham is called.
Cr   The basis, e.g. local orbitals, depend on the potential
Cu Updates
Cu   16 Aug 04 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      type(s_ctrl)::sctrl
      type(s_lat)::slat
      type(s_spec)::sspec(*)
      type(s_site)::ssite(*)
      type(s_ham)::sham
      type(s_pot)::spot

      character*(*) sstrn
C ... Local parameters
      integer nspec,lfp,nglob
C ... Needed for
Cchp1       integer ng,ogv,i
       integer:: ng , i 
       real(8),pointer :: rv_p_ogv(:)

      double precision alat,tolgv,dum

      integer w(1)
      common /w/ w

C --- Setup ---

      lfp=sctrl%lfp

      nspec = nglob('nspec')

C --- FP setup ---
      if (lfp .ne. 0) then
        call uspecb(0,0,sspec,1,nspec,dum,dum,dum,i)
        if (i .ge. 100) then

          ng=slat%ng
Cchp1           ogv=slat%ogv
           rv_p_ogv => slat%rv_p_ogv 

          alat=slat%alat
          tolgv=slat%tolft

Cchp1           call sugcut(2,nspec,sspec,alat,ng,w(ogv),tolgv)
           call sugcut ( 2 , nspec , sspec , alat , ng , rv_p_ogv , tolgv 
     .     ) 

        endif
      endif

      end subroutine suham2


