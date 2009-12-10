      subroutine rdctrlchk(recrd,recln,nrecs,
     . prgnam,vrsion,vn,vn2,
     . pass2,slabl_,osbz,osctrl,osham,ospot,oslat,osmix,osspec,ossite,
     . osstr,osarry,osmove,ostb,sstrn)
C- Main input for LMTO programs
C ----------------------------------------------------------------------
Ci Inputs
Ci   recrd (recln*nrecs) : preprocessed input
Ci   prgnam:name of main program
Ci   vrsion:string specifying expected program version
Ci   vn,vn2:major and minor versions
Ci   pass2 :flags whether call is 1st or 2nd pass. (2nd pass is sometimes
Ci         :used to read in class-specific info, e.g. ASA moments)
Co Outputs
Co   Input file is read and data is packed into these structures:
Co   slabl :vector of species labels
Co   sbz   :struct for the Brillouin Zone; see routine ubz
Co     Elts read: lmet lio,18 lmull fsmom
Co     Stored:    n w efmax lmet semsh zval ndos ef def range lio dosw
Co     Passed to: ubz dval rdccat
Co   sctrl :struct for program flow parameters; see routine uctrl
Co     Elts read: nbas nclass nspec nspin nl lncol lsx lscr lmet lrel
Co                lordn loptc lpgf mdprm lham,4 lxcf lfrce sdmod
Co                lasa lcd ltb lqp,2
Co     Stored:    lasa lfp lbas lcd lmet lqp lrel nspin nitmv lrs lxcf
Co                nl lpgf maxit smalit tol ltb zbak lncol sclwsr omax1
Co                omax2 nvario nsite nbas nspec modep
Co     Passed to: uctrl dval rdccat lgors lsets susite
Co   sham  :struct for parameters defining hamiltonian; see routine uham
Co     Elts read: lsig
Co     Stored:    rsrnge sigp rsstol lncol lxcf lham
Co     Passed to: uham dval susite
Co   spot  :struct for information about the potential; see routine upot
Co     Elts read: opnu oqnu oves opp osoptc
Co     Stored:    osoptc osgw
Co     Passed to: upot dval rdccat susite
Co   slat  :struct for lattice information; see routine ulat
Co     Elts read: alat avw
Co     Stored:    as nkdmx nkqmx tol gam tolft
Co     Passed to: ulat dval rdccat susite
Co   smix  :struct for charge mixing parameters; see routine umix
Co     Elts read: lxpot,3
Co     Stored:    fn r b bv wc w mmix nsave
Co     Passed to: umix dval spacks rdccat
Co   sspec :struct for species-specific information; see routine uspec
Co     Elts read: rmt
Co     Stored:    norp lmxa lmxpb hcr lmxf coreq pb1 pb2 coreh etf idxdn
Co     Passed to: uspec dval spackv spacks ioorbp scalss suidx
Co   ssite :struct for site-specific information; see routine usite
Co     Elts read:
Co     Stored:    relax
Co     Passed to: rdccat usite dval spackv
Co   sstr  :struct for parameters for screened strux; see routine ustr
Co     Elts read: skmsh n symg rmax
Co     Stored:    nkaps rmax rfit kaps lmaxw loka drwats
Co     Passed to: ustr dval rdccat
Co   sarry
Co     Elts read:
Co     Stored:
Co     Passed to: uarray dval susite
Co   smove :struct for dynamics information; see routine umove
Co     Elts read:
Co     Stored:    gyro prmint
Co     Passed to: umove dval rdccat
Co   sstrn :struct for global strings
Co     Elts read: symg
Co     Stored:
Co     Passed to: len rdccat parstr
Cg Global variables
Cg   The following global variables are set by rdctrl and may be accessed by
Cg   any routine via function call 'dglob' (for double) or 'nglob' (for int)
Cg   avw   :global length scale, usu. the average Wigner-Seitz radius,
Cg         :used in various places to set a length scale for a range,
Cg         :sometimes in generating structure constants, etc.
Cg   lrel  :specifies type of Schrodinger equation
Cg         :0 nonrelativistic Schrodinger equation
Cg         :1 scalar relativistic Schrodinger equation
Cg         :2 Dirac equation
Cg   lxcf  :specifies type of XC potential.  1s digit specifies local XC:
Cg         :1 for Ceperly-Alder
Cg         :2 for Barth-Hedin (ASW fit)
Cg         :3 for PW91
Cg         :4 for PBE
Cg         :10s digit specifies type of gradient correction
Cg         :0 no gradient correction
Cg         :1 Langreth-Mehl
Cg         :2 PW91
Cg         :3 PBE
Cg         :4 PBE with Becke exchange
Cg   mxorb :nkaph * (maximum number of lm channels in any sphere)
Cg         :Used for dimensioning the indexing arrays involved in
Cg         :assembling the hamiltonian;
Cg   nbas  :number of atoms in the basis
Cg   nbasp :number of atoms in the padded basis
Cg         :(when extensions are needed, e.g. in layer GF code)
Cg   nkape :NOT USED The maximum number of envelope functions centered at
Cg         :particular R and l channel
Cg         :NB: nkape is not used now.
Cg   nkaph :The maximum number of radial functions centered at
Cg         :particular R and l channel used in the lmto basis.
Cg   nl    :1+Maximum l-cutoff for augmentation
Cg   npl   :(not set by rdctrl) number of principal layers (layer geometries)
Cg   nkaph :The maximum number of "principal quantum" numbers centered
Cg         :at a particular R and l channel --- energies for one Rl
Cg         :at which augmentation (phi-phidot) functions are made.
Cg   nsp   :1 if not spin-polarized; otherwise 2
Cg   nspec :number of species
Cg   stde  :standard error file
Cg   stdl  :standard log file
Cg   stdo  :standard output file
Cr Remarks
Cr rdctrl does:
Cr  1. allocate the following structure arrays
Cr     osbz,osctrl,osham,ospot,oslat,osmix,osspec,ossite,osstr,osarry
Cr  2. read input data specified by tokens
Cr  3. If pass2, read class parameters from START
Cu Updates
Cu   19 Sep 07 (TK+MvS) Adapted from rdctrl, 1st cut at new input
Cu   20 Oct 06 Broadcast species so floating sites work properly in MPI
Cu   06 Aug 06 Remove defaults for STR RMAX and HCR
Cu   24 Nov 05 Remove mpi-specific calls
Cu   08 Jul 05 Assign nat as global variable
Cu             fix bug so --rdbasp works again
Cu   27 Mar 05 Add read option --rs=.,.,2,.. -> add 512 to lrs
Cu   21 Dec 04 Add switch to rotate FP local density on file read
Cu   16 Aug 04 Changes for extended local orbitals
Cu   18 Jun 04 printout of correct LDA+GGA functional
Cu   20 Dec 03 --rs rechecked in case made part of CMD in ctrl file
Cu   07 Sep 03 (lmf) rdctrl can read basis info from basis file
Cu   21 May 03 Added setup for sham->sigp
Cu   20 Mar 03 Change default for ctrl->tol:3 = etol
Cu   18 Mar 03 Added handling for fully relativistic case
Cu   15 Feb 02 (ATP) Added MPI parallelization
Cu   24 Aug 01 Extended to handle local orbitals.
Cu   28 Apr 98 code for new category 'OPTICS'
C ----------------------------------------------------------------------
      use m_rdctrlchk
      use m_gtv
C     implicit none
C ... Passed parameters
      integer:: recln,nrecs,ioff
      character strn*(recln)
      character*(1000) recrd
      logical pass2
      integer osbz,osctrl,osham,ospot,oslat,osmix,osspec,
     .  ossite,osstr,osarry,osmove,ostb
      character  prgnam*(*), sstrn*(*)
c      character toksw(0:30)*(*), vrsion*6
C     character  vrsion*6
      character(6):: vrsion(2)
      double precision vn(2),vn2(2)
C ... Heap
      integer w(1)
      common /w/ w
C ... Local parameters
      integer procid,master
      logical lgors,cmdopt,bittst,ltmp,parstr,ioorbp
      integer a2vec,bitand,fopna,getdig,i,igets,iprint,
     .  iprt,irs(5),ishow,isw,ifi,ix(n0*nkap0),j,k,l,lasa,lbas,lcd,
     .  lfrce,lfrzw,lgunit,lham,lmet,lncol,lnsph,loptc,llordn,lpgf,
     .  llrel,lrs,lstsym,lsx1,ltb,lxcf,nat,nbasl,nclasp,nlibu,
C     .  nclass,nl,nsite,nsp,nspec,
     .  nglob,nkap,nspc,nlmax,
     .  scrwid,stdo,stdl,stde,k1,k2,mpipid,lsx_
C ... BZ
      integer bz_lio
C ... basis
      double precision orbp(n0,2,nkap0)


c lscr,lsx,
      integer o,oclabl,ohave,oics,opnu,opp,oqnu,osgw,osoptc,
     .  osordn,oves,owk

      character slabl_(1)*8, fileid*64, prgn*12
      double precision dval,dglob,xx(n0),dgets
      double precision ekap(6)
      save ishow

      integer:: ioff2,irec

      character(600000):: rcd  !ifc henry is strange
c      character(60000):: rcd  !ifc henry is strange

      scrwid = 80
      stdo = lgunit(1)
      stdl = lgunit(2)
      stde = stdo

C --- Initialize gtv; copy recrd to rcd ---
      call gtv_setst(stdo,stdl,stde)
      call gtv_setrcd(recrd,nrecs,recln)

C --- Copy recrd to rcd
C      rcd = ' '
C      ioff2 = 1
C      do  irec = 1, nrecs
C        ioff = (irec-1)*recln + 1
CC       Blank line; skip
C        if (recrd(ioff:ioff)=='#' .or. recrd(ioff:ioff)=='!') then
CC       No character in 1st column
C        elseif (recrd(ioff:ioff)==' ') then
CC         rcd(ioff2:ioff2) = ' '
C          rcd(ioff2+1:) = adjustl(recrd(ioff:ioff+recln-1))
C          ioff2 = len_trim(rcd)+1
CC       Character in 1st col: add @CAT@ to mark category delimiter
C        else
CC         Replace CLASS with SPEC
C          if (recrd(ioff:ioff+4) .eq. 'CLASS') then
C            recrd(ioff:ioff+4) = 'SPEC'
C          endif
C          write(rcd(ioff2:ioff2+recln-1+6),"(1x,a,a)")
C     .    '@CAT@',recrd(ioff:ioff+recln-1)
C          ioff2 = len_trim(rcd)+1
C        endif
CC       write(*,"(i4,a)") ioff, rcd(1:ioff2)
C      enddo
C      rcd(ioff2:) = ' @CAT@EOF'

C --- Read input parameters from contents of rcd ---
      if (cmdopt('--prog=',7,0,fileid)) then
        prgn = fileid(8:)
        vrsion(2) = ' '
      else
        prgn = 'LM'
      endif
      if (prgn == 'LM' .or. prgn == 'LMGF' .or.  prgn == 'LMPG' .or.
     .  prgn == 'LMSTR' .or. prgn == 'LMCTL' .or. prgn == 'LMAVGM'
     .  .or. prgn == 'LMPLAN')
     .  vrsion(2) = 'ASA'          
      if (prgn=='LMF' .or. prgn=='LMFGWD')
     .  vrsion(2) = 'FP'
      if (prgn=='LMMC') vrsion(2) = 'MOL'
      if (prgn=='TBE') vrsion(2) = 'TB'
      if (prgn=='LMMAG') vrsion(2) = 'MM'
      call readctrl67(rcd,prgn,vrsion(1),vn(1))
      call rx0('done')

      end

