      module m_rdctrlchk
C     implicit none

C ... Fixed parameters
      integer,parameter:: NULLI =-99999
      real(8),parameter::    NULLR =-99999
      logical,parameter:: T=.true.,F=.false.
      integer,parameter :: n0=10
      integer,parameter :: nkap0=3
      integer(2):: nono

C ... IO
      integer:: io_show=0,io_help=0,nvario=0
c      character,allocatable:: jobid
c      character(128*20) :: cmd,const
c      character*(10000) :: jobdat,master
c      integer:: jobc,ix,nout

C ... HEADER, SYMGRP
      character(256):: header,symg=' '

C ... HAM
      logical :: frzwf,ham_ewald
      integer:: ctrl_lfrce,ham_lxcf,gga,ftmesh(3),nmto=0,lrsig=0,
     .  nsp=1,ham_qasa,lrel=1,lso=0
      logical:: ltbe  ! set to T if TBE program (used for defaults)
      integer lfp  ! set to T if FP program (used for defaults)
      real(8):: lat_gmax,tolft,elind,dqval,kmto(10),rsrnge,vmtz
      real(8):: alfsi=nullr,dabc(3)=nullr,rsstol
      real(8):: pmin(n0),pmax(n0)
C     sigp holds parameters for approximating self-energy sigma
C       arg 1: mode : specifies how to set its diagonal part
C              for states above the high-energy cutoff
C              0 constrain sigii to be > asig+bsig*e
C              1 constrain sigii to be = asig+bsig*e
C              2 constrain sigii to be > asig and < bsig
C              3 constraint same as mode 1.
C                Mode 3 differs in that the least-squares fit to
C                sigii (for informational purposes only, to help
C                estimate asig and bsig) is done for states between
C                efit and nmax or emax
C       arg 2: nmin : sigma for states 1..nmin are approximated by sigii
C       arg 3: emin : (used only if nmin<0)
C                   : sigma for levels e<emin are approximated by sigii
C       arg 4: nmax : sigma for levels i>nmax are approximated by
C                     sigii AND constrained according to mode
C       arg 5: emax : (used only if nmax<=0)
C                   : sigma for levels e<emax are approximated by
C                     sigii AND constrained according to mode
C       arg 6: asig : constraint used to approximate
C                     sigii = asig + E * bsig  or
C                     asig < sigii < bsig
C       arg 7: bsig : constraint used to approximate
C                     sigii = asig + E * bsig  or
C                     asig < sigii < bsig
      real(8):: sigp(10)=0,sigp_emin,sigp_emax,sigp_a,sigp_b,sigp_efit
      integer:: sigp_mode=3,sigp_nmin=0,sigp_nmax=0
      equivalence (sigp_emin,sigp(3)),(sigp_emax,sigp(5)),
     .             (sigp_a,sigp(6)),(sigp_b,sigp(7)),(sigp_efit,sigp(8))
C      real(8):: sigp_mode,sigp_nmin,sigp_nmax
C      equivalence (sigp_mode,sigp(1)),(sigp_nmin,sigp(2)),
C     .  (sigp_nmax,sigp(4))

C ... OPTIONS
      logical ::
     .  lcd1,lcd2,lcd4=F,lcd8=F,lcd64,lasa4=T,lasa8=F,lasa32=F,lasa64=F,
     .  lham1,lham4=F,lham8,lham16,lham32=F,lham64,lham128,lham256,lves
C     Initialize, since read is only executed when NSPIN=2
      logical:: lncol1=F,lncol2=F,lncol4=F,lncol8=F,lncol16=F,
     .  lncol32=F,lncol64=F
      character(256):: sxopt=' '
      integer:: lasa3=1,lsx,lscr,smalit(2)=(/80,2/),nesabc(3),
     .  lstonr(3)=0,quit,nl,lham3,lpfloat
      real(8):: asa_elin,rmines,rmaxes,ham_qss(4)

C ... STRUC
C     Initializing nbas and alat flags no nbas or alat has been input
      real(8):: lat_slat(9),dlat,alat=NULLR,plat(9),lat_gam(4),dalat
      real(8):: vol,avw,lat_dist(9)
      integer:: nbas=NULLI,nbasp=NULLI,nsite=NULLI,nspec,nclass
      integer:: lat_ldist=0

C ... SPEC

C ... SITE

C ... Iterations

C ... BZ

C ... TB

C ... Ewald

C ... STR

C ... TCF

C ... CGF

C ... PGF

C ... DYN

C ... GW

C ... PLANE

C ... ORDERN

C ... OPTICS

      contains

      subroutine readctrl67(rcd_in,prgn,vstrn,vn)
C- Reads data parsed from ctrl file into m_rdctrlchk
C ----------------------------------------------------------------------
Ci Inputs
Ci   rcd_in:string with contents ctrl file
Ci   prgn:name of main program
Ci   vstrn :version
Ci   vn    :
Co Outputs
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   09 Aug 07
C ----------------------------------------------------------------------
      use m_toksw
      use m_gtv
C     implicit none
C ... Passed parameters
      character*(*):: rcd_in, prgn
      character(6):: vstrn(2)
      double precision vn(2)
C ... Local variables
      integer::j,lgunit
      integer:: sw
      integer:: jj(2)
      character(256):: a,outs
      logical :: cmdopt
      double precision vsn,xv(2*n0)
      logical :: ltmp
      logical,parameter:: T=.true.,F=.false.
      logical::  debug=.false.

      character(128) :: nm
      integer :: stdo,stdl,stde

C     call shosyv(0,0,0,6)

C     Used to set tbe-dependent or f-dependent defaults
      ltbe = prgn .eq. 'TBE'
      lfp = 0
      if (prgn=='LMF' .or. prgn=='LMFGWD' .or. prgn=='LMMC') lfp = 1

C --- Initialize ---
      debug = cmdopt('--debug',6,0,a)
      call toksw_init_chk(debug)
      stdo = lgunit(1)
      stdl = lgunit(2)
      stde = lgunit(1)
C     call gtv_setst(,stdo,stdl,stde)

      write(stdo,331) trim(prgn)
  331 format(/' --- Convert input v6 to v7 for program ',a,' ---'/)

C     In this program, never show, but there is a help mode
      io_show = 0
      if (cmdopt('--input',7,0,a)) io_help = 1
      call gtv_setio(debug,io_show,io_help)
      if (io_help .eq. 1) then
        write(stdo,332)
  332   format(/' Token',t19,'Input   cast  (size,min)'/
     .    ' ------------------------------------------')
      elseif (io_show .ne. 0) then
        call gtv_setio(debug,io_show,io_help)
        write(stdo,333)
  333   format(/' Token',t19,
     .    'Input   cast  (size,min,read,def)     result')
      endif

C     goto 1234

C --- Struc ---
      nm='STRUC_NSPEC'; sw = tksw(prgn,nm)
      if (sw .eq. 1) sw = 0
      call gtv(trim(nm),sw,nspec,def_i4=1,
     .  note='Number of species to read from SPEC category.'//
     .  '%N%3fIf not present, NSPEC assumes the value 1')
      if (nspec == NULLI) nspec = 1

      call gtvchk(2,prgn,'STRUC_NCLASS',
     .  'Number of classes',
     .  nwtok='STRUC_NSPEC')


C --- IO ---
C      nm='IO_HELP'; call gtv(trim(nm),tksw(prgn,nm),io_help,
C     .  def_i4=0, note='Show what input would be sought, '//
C     .  'without attempting to read data')
C      if (cmdopt('--input',7,0,a)) io_help = 1  !optio=0 in old code
C      if (io_help == 1) io_show = 1
C      call gtv_setio(debug,io_show,io_help)

C --- Version ---
      do  j = 1, 2
      vsn = dble(int(vn(j)*100))/100
      if (vsn .ne. 0d0) then
        nm = 'VERS_'//vstrn(j)
        outs = 'Input style'
        if (j .eq. 2) outs='Program'
        call gtvchk(0,prgn,nm,trim(outs)//' version check')
      endif
      enddo

C --- Options ---
C --- Hamiltonian parameters ---
      call gtvchk(2,prgn,'OPTIONS_NSPIN',
     .  'Set to 2 for spin polarized calculations',nwtok='HAM_NSPIN')
C     .  note2='... Replace with `HAM_NSPIN''')

      call gtvchk(2,prgn,'OPTIONS_REL',
     .  'for relativistic Schrodinger equation',nwtok='HAM_REL')

      call gtvchk(2,prgn,'OPTIONS_SX',
     .  'Screened exchange:',nwtok='HAM_SX')

      call gtvchk(2,prgn,'OPTIONS_SXOPTS',
     .  'Screened exchange:',nwtok='HAM_SXOPTS')

      call gtvchk(2,prgn,'OPTIONS_NRMIX',
     .  'Small iterations, sphere program',nwtok='ITER_NRMIX')

      call gtvchk(2,prgn,'OPTIONS_XCFUN','LD EXC functional',
     .  nwtok='HAM_XCFUN')

      call gtvchk(2,prgn,'OPTIONS_SO','for spin-orbit coupling',
     .  nwtok='HAM_SO')

      call gtvchk(2,prgn,'OPTIONS_NONCOL','Noncollinear magnetism',
     .  nwtok='HAM_NONCOL')

      call gtvchk(2,prgn,'OPTIONS_BFIELD','Applied magnetic field',
     .  nwtok='HAM_BFIELD')

      call gtvchk(2,prgn,'OPTIONS_SS','Magnetic spin spiral',
     .  nwtok='HAM_SS')

      call gtvchk(2,prgn,'OPTIONS_ZBAK',
     .  'turns on Ewald MT correction;  adds hom. back. charge',
     .  nwtok='OPTIONS_MTCOR',note2='... use BZ_ZBAK for bg. chg,'//
     .  ' OPTIONS_MTCOR for MT corr')

      call gtvchk(0,prgn,'OPTIONS_ASA','Contains ASA-specific tokens')

      call gtvchk(2,prgn,'OPTIONS_ADNF','Turn on automatic downfolding',
     .  nwtok='OPTIONS_ASA_ADNF')

      call gtvchk(2,prgn,'OPTIONS_NSPH',
     .  'generate multipole moments in output density',nwtok=
     .  'OPTIONS_ASA_NSPH')

      call gtvchk(2,prgn,'OPTIONS_TWOC','Two-center ASA hamiltonian',
     .  nwtok='OPTIONS_ASA_TWOC')

      call gtvchk(2,prgn,'OPTIONS_GAMMA',
     .  'Gamma representation (2nd generation ASA)',nwtok=
     .  'OPTIONS_ASA_GAMMA')

      call gtvchk(2,prgn,'OPTIONS_CCOR','Turn on combined correction',
     .  nwtok='OPTIONS_ASA_CCOR')

      call gtvchk(2,prgn,'OPTIONS_ELIN',
     .  'Energy to linearize for CCOR (2nd gen ASA 2C hamiltonian)',
     .  nwtok='OPTIONS_ASA_ELIN')

      call gtvchk(2,prgn,'OPTIONS_NEWREP',
     .  'Interactively transform representation, (2nd gen ASA)
     .  ',nwtok='OPTIONS_ASA_NEWREP')

      call gtvchk(2,prgn,'OPTIONS_NOHYB','Turns off hybridisation ',
     .  nwtok='OPTIONS_ASA_NOHYB')

      call gtvchk(2,prgn,'OPTIONS_MTCOR','Turns on Ewald MT correction',
     .  nwtok='OPTIONS_ASA_MTCOR')

      call gtvchk(2,prgn,'OPTIONS_MTCOR_QMT',
     .  'Override standard background charge',nwtok=
     .  'OPTIONS_ASA_QMT')

      call gtvchk(3,prgn,'HAM_QASA','conventions for 2nd gen '//
     .  'ASA moments',note2='default value has changed to 3')

      nm = 'HAM_SIGP'
      call gtv(nm,tksw(prgn,nm),nono,Texist=ltmp) ! see if new style is present
      call gtvchk(4,prgn,
     .  'HAM_SIGP','Parameters for spin statics',
     .  note2='Replace `HAM_SIGP: ...'' with `HAM_SIGP[tokens]''',
     .  news=ltmp)

C --- Symmetry group ---
C --- Crystal Green's function ---
C --- Planar Green's function ---
C --- Species (old CLASS) ---
      if (io_help == 0) then
        print *,' '
      else
        call info0(0,1,0,
     .    ' ... the following checks apply to each species')
      endif
      do  j = 1, nspec
        if (io_help == 0) then
        call info2(0,0,0,' ... checking species %i of %i',j,nspec)
        endif
        jj= (/1,j/)
        call gtvchk(3,prgn,'SPEC_ATOM_A',
     .    'Radial mesh point spacing parameter',cindx=jj,
     .    note2='default value has changed from 0.03 to 0.025')

        call gtvchk(3,prgn,'SPEC_ATOM_LMXA',
     .    'l-cutoff for augmentation',cindx=jj,
     .    note2='default value may be different (depends on NL)')

      enddo
      call info0(0,0,1,' ... done with species')
        
C --- Site ---
      call gtvchk(1,prgn,'SITE_MODE',
     .  'Set to T if input positions are fractions of plat',note2=
     .  '... SITE_ATOM_XPOS reads site positions as fractions of plat')

C --- TCF ---
      call gtvchk(2,prgn,'STR','Parameters for two-center fit',
     .  nwtok='TCF')

C --- Structure constants ---
      nm = 'STR_ENV_MODE'
      call gtv(nm,tksw(prgn,nm),nono,Texist=ltmp) ! see if new style is present
      call gtvchk(4,prgn,
     .  'STR_MODE','Type of envelope functions (structure constants)',
     .  note2='Replace `STR_MODE ...'' with `STR_ENV_MODE''',
     .  news=ltmp)

C --- Two-center fit ---
C --- Brillouin Zone ---
      call gtvchk(1,prgn,'BZ_N.W',
     .  'Polynomial order and broadening for M-P sampling integration',
     .  note2='... Use BZ_N and BZ_W separately')

      call gtvchk(3,prgn,'BZ_EFMAX','Find evecs up to efmax',
     .  note2='default value has changed to 2')

C --- Ewald sums ---
C --- Iterations (formerly MIX) ---
      if (tksw(prgn,'MIX')/=2) then
        call gtvchk(2,prgn,'MIX','Parameters for mixing, iterations',
     .  nwtok='ITER')

        call gtvchk(2,prgn,'MIX_MODE','Mixing rules for charge mixing',
     .    nwtok='ITER_MIX')

        call gtvchk(2,prgn,'MIX_CONV','Tolerance in energy change '//
     .    'from prior iteration for self-consistency',nwtok='ITER_CONV')

        call gtvchk(2,prgn,'MIX_CONVC','Tolerance in output-input'//
     .    ' charge for self-consistency',nwtok='ITER_CONVC')

        call gtvchk(2,prgn,'MIX_AMODE','Mixing rules for Euler '//
     .    'angle mixing',nwtok='ITER_AMIX')

C       call snit
        call gtvchk(3,prgn,'ITER_CONV','Tolerance in energy change '//
     .    'from prior iteration for self-consistency',
     .    note2='default value has changed to 0 to 1E-4')

      endif
C --- Master ---
      call gtvchk(1,prgn,'MASTER','Sets job-dependent variables',
     .  note2='... Select one job; variables in CONST')

C --- TB ---
C --- Map (do not maintain for now) ---
C --- Optics ---
C --- Dynamics ---
      nm = 'DYN_MSTAT_MODE'
      call gtv(nm,tksw(prgn,nm),nono,Texist=ltmp) ! see if new style is present
      call gtvchk(4,prgn,'DYN_MSTAT','Parameters for molecular statics',
     .  note2='Replace `MSTAT: ...'' with `MSTAT[tokens]''',news=ltmp)

      nm = 'DYN_SSTAT_MODE'
      call gtv(nm,tksw(prgn,nm),nono,Texist=ltmp) ! see if new style is present
      call gtvchk(4,prgn,
     .  'OPTIONS_SDYN','Parameters for spin statics',
     .  note2='Replace `OPTIONS_SDYN: ...'' with `DYN_SSTAT[tokens]''',
     .  news=ltmp)

CC --- GW ---
CC --- Plane ---
CC --- Start ---
      call gtvchk(2,prgn,'START_NIT',
     .  'maximum number of iterations in self-consistency cycle',
     .  nwtok='ITER_NIT')

      call gtvchk(2,prgn,'START_CNVG','Tolerance in output-input'//
     .  ' charge for self-consistency',nwtok='ITER_CONVC')

      stop

      end subroutine

      subroutine gtvchk(mode,prgn,token,note1,nwtok,note2,news,cindx)
C- Find token and print message about how to update input
C ----------------------------------------------------------------------
Ci Inputs
Ci  mode   :0 token is now input (not previously)
Ci         :1 token is no longer used
Ci         :2 token is no longer used; replace by nwtok
Ci         :3 token's default value has been changed.
Ci         :4 input style has been changed
Ci  prgn   :name of program
Ci  token  :name of token, including parents.
Ci         :mode=0: token is now sought
Ci         :mode=1,2: token is no longer sought
Ci         :mode=2: token is to be replaced by nwtok
Co  Texist :.true. if the token could be matched; otherwise .false.
Ci Inputs (optional)
Ci  note1  :printout string for token
Ci  nwtok  :name of new token (mode 2)
Ci  note2  :printout string for nwtok
Ci  cindx  :used to indicate multiple occurences of a token
Ci  news   :(mode 4 only): T if new syntax is being used.
Co Outputs (optional)
Cl Local variables
Cr Remarks
Cu Updates
Cu   10 Oct 07
C ----------------------------------------------------------------------
      use m_toksw
      use m_gtv
C     implicit none
C ... Passed parameters
      character(*),intent(in):: prgn,token,note1
      character(*),optional,intent(in):: note2,nwtok
      integer,intent(in):: mode
      integer,optional,intent(in):: cindx(1:)
      logical,optional,intent(in):: news
C ... Local parameters
      logical:: ltmp,ltmp2,lnote2,lnews
      integer:: sw,sw1
      character(32):: strn1,strn2

      sw = tksw(prgn,token)
      if (sw .eq. 2) return
      sw1 = sw
      if (io_help .eq. 0 .and. sw .eq. 1) sw1 = 0
      lnote2 = present(note2)
      if (present(news)) then
        lnews = news
      else
        lnews = .false.
      endif

      if (mode .eq. 0) then
        call gtv_none(token,sw1,nono,Texist=ltmp,cindx=cindx,
     .    note='(new): '//trim(note1))
        if (io_help == 0 .and. .not. ltmp) then
          if (sw == 0) then
            write(*,333) trim(token)
          else
            write(*,333) trim(token),'(required)'
          endif
  333     format(' Missing new token  ',a:,2x,a)
        endif
      elseif (mode .eq. 1 .or. mode .eq. 4) then
        strn1 = '(no longer used)'
        strn2 = 'is no longer used'
        if (mode .eq. 4) then
          strn1 = '(syntax changed)'
          strn2 = 'syntax changed'
        endif
        call gtv_none(token,sw1,nono,Texist=ltmp,cindx=cindx,
     .    note=trim(strn1)//':  '//trim(note1))
        if (io_help == 1 .and. lnote2) then ! Help mode
          write(*,'(3x,a)') note2
C       New syntax present, old syntax absent
        elseif (io_help == 0 .and. lnews .and. .not. ltmp) then 
          write(*,338) trim(token),'ok (new syntax)'
        elseif (io_help == 0 .and. lnews .and. lnote2) then ! Using new syntax
          write(*,338) trim(token),'ok (new syntax)',trim(note2)
  338     format('  Token `',a,''' ',a:' ...',t45,a:1x,a)
        elseif (io_help == 0 .and. lnews) then ! Using new syntax
          write(*,338) trim(token),'ok (new syntax)'
        elseif (io_help == 0 .and. ltmp .and. lnote2) then ! old token or syntax
          write(*,334) trim(token),trim(strn2),note2
        elseif (io_help == 0 .and. ltmp) then ! old token or syntax
          write(*,334) trim(token),trim(strn2)
  334     format(' *Token `',a,''' ',a:' ...',t45,a:1x,a)
        endif
      elseif (mode .eq. 2) then
        if (.not. present(nwtok)) call rx(
     .    'gtvchk requires nwtok for mode=2')
          call gtv_none(token,sw1,nono,Texist=ltmp,cindx=cindx,
     .      note='(replace with `'//trim(nwtok)//'''): '//trim(note1))
        if (io_help == 0) then
          strn1 = '(no longer used): '
          call gtv_none(nwtok,sw1,nono,Texist=ltmp2,cindx=cindx)
          if (ltmp2 .and. ltmp) then
            write(*,337) trim(nwtok),trim(token)
  337       format('  Token `',a,''' found ... ',t45,
     .        'remove obsolete ','`',a,'''')
          elseif (ltmp .and. lnote2) then
            write(*,334) trim(token),trim(strn1),note2
          elseif (ltmp) then
            write(*,335) trim(token),trim(nwtok)
  335       format(' *Token `',a,
     .        ''' no longer used ...',t45,'move to `',a,'''')
          endif
        endif
      elseif (mode .eq. 3) then
        call gtv_none(token,sw1,nono,Texist=ltmp,cindx=cindx,
     .    note='(new default value): '//trim(note1))
        if (io_help == 0 .and. .not. ltmp) then
        if (lnote2) then
          write(*,336) trim(token),note2
        else
          write(*,336) trim(token)
        endif
        endif
  336   format(' *Token `',a,''' is missing':' ...',t45,a)
      else
        call rx('gtvchk: bad mode')
      endif

      end subroutine gtvchk

      end module

      subroutine toksw_init_chk(debug)
      use m_toksw
      character*(200):: io_sw, asa_sw
      logical:: debug

C --- Store data into arrays above ---
      call clear_swtok(debug)
      io_sw = ' HEADER~ CONST~ CMD~ IO_SHOW~ IO_HELP~ IO_VERBOS~ '//
     .  'IO_WKP~ IO_IACTIV~ IO_TIM~' // ' STRUC_NSPEC~ STRUC_NCLASS~'
      asa_sw = ' OPTIONS~ OPTIONS_ASA OPTIONS_ADNF~ '
     .  //'OPTIONS_NSPH~ OPTIONS_TWOC~ OPTIONS_CCOR~ '
     .  //'OPTIONS_GAMMA~ OPTIONS_ELIN~ OPTIONS_NEWREP~ '
     .  //'OPTIONS_NOHYB~ OPTIONS_MTCOR~ '
     .  //'OPTIONS_MTCOR_QMT~'
C      asa_sw = ' OPTIONS~ OPTIONS_ASA OPTIONS_ASA_ADNF~ '
C     .  //'OPTIONS_ASA_NSPH~ OPTIONS_ASA_TWOC~ OPTIONS_ASA_CCOR~ '
C     .  //'OPTIONS_ASA_GAMMA~ OPTIONS_ASA_ELIN~ OPTIONS_ASA_NEWREP~ '
C     .  //'OPTIONS_ASA_NOHYB~ OPTIONS_ASA_MTCOR~ '
C     .  //'OPTIONS_ASA_MTCOR_QMT~ OPTIONS_ASA_SCR~'
C      struc_sw =
C     .  ' STRUC STRUC_NSPEC~ STRUC_NCLASS STRUC_FILE~ STRUC_NBAS '//
C     .  'STRUC_PLAT STRUC_ALAT STRUC_DALAT~ STRUC_NL~ STRUC_SHEAR~ '//
C     .  'STRUC_ROT~ STRUC_DEFGRD~ STRUC_STRAIN~ STRUC_ALPHA~ '//
C     .  'STRUC_LOADV~ '

C --- LM switches ---
      call nswadd()
      call tkadd(" LM::" )
      call tkadd(" VERS_LM")
      call tkadd(" VERS_ASA")
      call tkadd(trim(io_sw))
      call tkadd(trim(asa_sw))
      call tkadd(" OPTIONS_NSPIN~ OPTIONS_REL~ OPTIONS_SO~")
      call tkadd(" OPTIONS_XCFUN~")
      call tkadd(" OPTIONS_NONCOL~ OPTIONS_SS~ OPTIONS_BFIELD~")
      call tkadd(" OPTIONS_SX~  OPTIONS_SXOPTS~")
C     call tkadd(" OPTIONS_QASA~")
      call tkadd(" OPTIONS_ZBAK~")
      call tkadd(" OPTIONS_SDYN~")
      call tkadd(" OPTIONS_NRMIX~")
      call tkadd(" HAM_QASA~")
      call tkadd(" SITE_MODE~")
      call tkadd(" SPEC_ATOM_A~")
      call tkadd(" BZ_N.W~ BZ_EFMAX~")
      call tkadd(" MIX~")
      call tkadd(" MIX_MODE~")
      call tkadd(" MIX_CONV~")
      call tkadd(" MIX_CONVC~")
      call tkadd(" MIX_AMODE~")
      call tkadd(" DYN_SSTAT~")
      call tkadd(" DYN_SSTAT_MODE~")
      call tkadd(" ITER_CONV~ ITER_NIT~")
      call tkadd(" MASTER~")
      call tkadd(" START_NIT~")
      call tkadd(" START_CNVG~")
C --- LMFA switches ---
      call nswadd()
      call tkadd(" LMFA::" )
      call tkadd(" VERS_LM")
C     call tkadd(" VERS_FP")
      call tkadd(trim(io_sw))
      call tkadd(" OPTIONS_NSPIN~ OPTIONS_REL~ OPTIONS_SO~")
      call tkadd(" OPTIONS_XCFUN~")
      call tkadd(" SPEC_ATOM_A~")
      call tkadd(" SITE_MODE~")
      call tkadd(" MIX~")
C --- LMFGWD switches ---
      call nswadd()
      call tkadd(" LMFGWD::" )
      call tkadd(" VERS_LM")
      call tkadd(" VERS_FP")
      call tkadd(trim(io_sw))
      call tkadd(" OPTIONS_NSPIN~ OPTIONS_REL~ OPTIONS_SO~")
      call tkadd(" OPTIONS_XCFUN~")
      call tkadd(" HAM_SIGP~")
      call tkadd(" SPEC_ATOM_A~ SPEC_ATOM_LMXA~")
      call tkadd(" SITE_MODE~")
C --- LMF switches ---
      call nswadd()
      call tkadd(" LMF::" )
      call tkadd(" VERS_LM")
      call tkadd(" VERS_FP")
      call tkadd(trim(io_sw))
      call tkadd(" OPTIONS_NSPIN~ OPTIONS_REL~ OPTIONS_SO~")
      call tkadd(" OPTIONS_XCFUN~")
      call tkadd(" SPEC_ATOM_A~ SPEC_ATOM_LMXA~")
      call tkadd(" SITE_MODE~")
      call tkadd(" BZ_N.W~")
      call tkadd(" MIX~")
      call tkadd(" MIX_MODE~")
      call tkadd(" MIX_CONV~")
      call tkadd(" MIX_CONVC~")
C     call tkadd(" MIX_AMODE~")
      call tkadd(" HAM_SIGP~")
      call tkadd(" DYN_MSTAT~")
      call tkadd(" DYN_NIT~")
      call tkadd(" DYN_MSTAT_MODE~")
      call tkadd(" START_NIT~")
C --- LMMC switches ---
      call nswadd()
      call tkadd(" LMMC::" )
      call tkadd(" VERS_LM")
      call tkadd(" VERS_MOL")
      call tkadd(trim(io_sw))
      call tkadd(" OPTIONS_NSPIN~ OPTIONS_REL~ OPTIONS_SO~")
      call tkadd(" OPTIONS_XCFUN~")
      call tkadd(" SPEC_ATOM_A~")
      call tkadd(" STR")
      call tkadd(" MIX~")
      call tkadd(" MIX_MODE~")
      call tkadd(" MIX_CONV~")
      call tkadd(" MIX_CONVC~")
C     call tkadd(" MIX_AMODE~")
      call tkadd(" START_NIT~")
C --- LMGF switches ---
      call nswadd()
      call tkadd(" LMGF::" )
      call tkadd(" VERS_LM")
      call tkadd(" VERS_ASA")
      call tkadd(trim(io_sw))
      call tkadd(" OPTIONS_NSPIN~ OPTIONS_REL~ OPTIONS_SO~")
      call tkadd(" OPTIONS_XCFUN~")
      call tkadd(" OPTIONS_NONCOL~ OPTIONS_SS~ OPTIONS_BFIELD~")
      call tkadd(" OPTIONS_NRMIX~")
      call tkadd(" SPEC_ATOM_A~")
      call tkadd(" SITE_MODE~")
C     call tkadd(" BZ_N.W~ BZ_EFMAX~")
      call tkadd(" MIX~")
      call tkadd(" MIX_MODE~")
      call tkadd(" MIX_CONV~")
      call tkadd(" MIX_CONVC~")
      call tkadd(" MIX_AMODE~")
      call tkadd(" ITER_CONV~")
      call tkadd(" START_NIT~")
      call tkadd(" START_CNVG~")
C --- LMDOS switches ---
      call nswadd()
      call tkadd(" LMDOS::" )
      call tkadd(" VERS_LM")
      call tkadd(trim(io_sw))
      call tkadd(" OPTIONS_NSPIN~")
      call tkadd(" SITE_MODE~")
      call tkadd(" BZ_N.W~ BZ_EFMAX~")
C --- LMCHK switches ---
      call nswadd()
      call tkadd(" LMCHK::" )
      call tkadd(" VERS_LM")
      call tkadd(trim(io_sw))
      call tkadd(" OPTIONS_NSPIN~ OPTIONS_NONCOL~ OPTIONS_SS~")
      call tkadd(" OPTIONS_NONCOL~ OPTIONS_SS~")
      call tkadd(" SPEC_ATOM_A~")
      call tkadd(" SITE_MODE~")
C --- LMSCELL switches ---
      call nswadd()
      call tkadd(" LMSCELL::" )
      call tkadd(" VERS_LM")
      call tkadd(trim(io_sw))
      call tkadd(" OPTIONS_NSPIN~ OPTIONS_NONCOL~")
      call tkadd(" SITE_MODE~")
C --- LMSTR switches ---
      call nswadd()
      call tkadd(" LMSTR::" )
      call tkadd(" VERS_LM")
      call tkadd(" VERS_ASA")
      call tkadd(trim(io_sw))
      call tkadd(" OPTIONS_NSPIN~")
      call tkadd(" SITE_MODE~")
      call tkadd(" STR_MODE~")
      call tkadd(" STR_ENV_MODE~")
C --- LMXBS, switches ---
      call nswadd()
      call tkadd(" LMXBS,::" )
      call tkadd(" VERS_LM")
      call tkadd(trim(io_sw))
      call tkadd(" SITE_MODE~")
C --- LMCTL switches ---
      call nswadd()
      call tkadd(" LMCTL::" )
      call tkadd(" VERS_LM")
      call tkadd(" VERS_ASA")
      call tkadd(trim(io_sw))
      call tkadd(" OPTIONS_NSPIN~")
      call tkadd(" OPTIONS_NONCOL~ OPTIONS_SS~")
      call tkadd(" SPEC_ATOM_A~")
      call tkadd(" SITE_MODE~")
C --- LMPG switches ---
      call nswadd()
      call tkadd(" LMPG::" )
      call tkadd(" VERS_LM")
      call tkadd(" VERS_ASA")
      call tkadd(trim(io_sw))
      call tkadd(" OPTIONS_NSPIN~ OPTIONS_REL~ OPTIONS_SO~")
      call tkadd(" OPTIONS_XCFUN~")
      call tkadd(" OPTIONS_NONCOL~ OPTIONS_SS~ OPTIONS_BFIELD~")
      call tkadd(" OPTIONS_NRMIX~")
      call tkadd(" SPEC_ATOM_A~")
      call tkadd(" SITE_MODE~")
C     call tkadd(" BZ_N.W~ BZ_EFMAX~")
      call tkadd(" MIX~")
      call tkadd(" MIX_MODE~")
      call tkadd(" MIX_CONV~")
      call tkadd(" MIX_CONVC~")
      call tkadd(" MIX_AMODE~")
      call tkadd(" ITER_CONV~")
      call tkadd(" START_NIT~")
      call tkadd(" START_CNVG~")
C --- LMPLAN switches ---
      call nswadd()
      call tkadd(" LMPLAN::" )
      call tkadd(" VERS_LM")
      call tkadd(trim(io_sw))
      call tkadd(" OPTIONS_NSPIN~")
      call tkadd(" OPTIONS_NONCOL~ OPTIONS_SS~")
C     call tkadd(" SPEC_ATOM_A~")
      call tkadd(" SITE_MODE~")
C --- LMIMP switches ---
      call nswadd()
      call tkadd(" LMIMP::" )
      call tkadd(" VERS_LM")
      call tkadd(" VERS_ASA")
      call tkadd(trim(io_sw))
      call tkadd(" OPTIONS_NSPIN~")
      call tkadd(" SPEC_ATOM_A~")
      call tkadd(" SITE_MODE~")
C --- LMAVGM switches ---
      call nswadd()
      call tkadd(" LMAVGM::" )
      call tkadd(" VERS_LM")
      call tkadd(" VERS_ASA")
      call tkadd(trim(io_sw))
      call tkadd(" OPTIONS_NSPIN~")
      call tkadd(" SITE_MODE~")
C --- TBE switches ---
      call nswadd()
      call tkadd(" TBE::" )
      call tkadd(" VERS_LM")
      call tkadd(" VERS_TB")
      call tkadd(trim(io_sw))
      call tkadd(" OPTIONS_NSPIN~ OPTIONS_SO~ OPTIONS_XCFUN~")
      call tkadd(" SITE_MODE~")
      call tkadd(" BZ_N.W~")
      call tkadd(" MIX~")
      call tkadd(" MIX_MODE~")
      call tkadd(" MIX_CONV~")
      call tkadd(" MIX_CONVC~")
C     call tkadd(" MIX_AMODE~")
      call tkadd(" START_NIT~")
      call tkadd(" START_CNVG~")
C --- LMMAG switches ---
      call nswadd()
      call tkadd(" LMMAG::" )
      call tkadd(" VERS_LM")
      call tkadd(" VERS_MM")
      call tkadd(trim(io_sw))
      call tkadd(" OPTIONS_SO~ OPTIONS_SS~ OPTIONS_BFIELD~")
      call tkadd(" OPTIONS_XCFUN~")
      call tkadd(" SITE_MODE~")
      end

C      subroutine snit
C      return
C      end

