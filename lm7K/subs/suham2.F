Cgetarg...info...           structure ['suham2', 'ctrl', 'sctrl']
Cgetarg...info...           structure ['suham2', 'lat', 'slat']
Cgetarg...info...           use_to_add ['suham2', 'm_struc_def']
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
Cgetarg       double precision sctrl(1),slat(1),sspec(1),ssite(1),
Cgetarg      .  sham(1),spot(1)
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
      integer ng,ogv,i
      double precision alat,tolgv,dum

      integer w(1)
      common /w/ w

C --- Setup ---
Cgetarg...info...           integer :: lfp , integer(8) :: sctrl%lfp 1
Cgetarg       call upack('ctrl lfp',sctrl,lfp,0,0,0,0)
       
       lfp=sctrl%lfp

      nspec = nglob('nspec')

C --- FP setup ---
      if (lfp .ne. 0) then
        call uspecb(0,0,sspec,1,nspec,dum,dum,dum,i)
        if (i .ge. 100) then
Cgetarg...info...           integer :: ng , integer(8) :: slat%ng 1
Cgetarg...info...           integer :: ogv , integer(8) :: slat%ogv 1
Cgetarg...info...           real(8) :: alat , real(8) :: slat%alat 1
Cgetarg...info...           real(8) :: tolgv , real(8) :: slat%tolft 1
Cgetarg           call upack('lat ng ogv alat tolft',slat,ng,ogv,alat,tolgv,0)
           
           ng=slat%ng
           ogv=slat%ogv
           alat=slat%alat
           tolgv=slat%tolft

          call sugcut(2,nspec,sspec,alat,ng,w(ogv),tolgv)
        endif
      endif

Cgetarg       end
       end subroutine suham2 


