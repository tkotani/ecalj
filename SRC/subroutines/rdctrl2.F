!! note. all initial data is pushed into m_lmfinit
      subroutine rdctrl2( recrd,recln,nrecs, prgnam)
      use m_work,only: iv_a_oips,rv_a_opos
      use m_work,only: rv_a_ocy,rv_a_ocg,iv_a_ojcg,iv_a_oidxcg !CG coefficients
      use m_globalvariables !global variables
      use m_lmfinit !contains main initial data
      use m_struc_func
      use m_toksw,only:tksw
      use m_gtv 
!! followings are old document ------
!!
c      use m_susite_func, only: susite
C- Main input for LMTO programs
C ----------------------------------------------------------------------
Ci Inputs
Ci   recrd (recln*nrecs) : preprocessed input
Ci   prgnam:name of main program
Ci   vrsion:string specifying expected program version
Ci   vn,vn2:major and minor versions
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
Cg         :103 for PBE
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
Cr     v_sbz,v_sctrl,v_sham,v_spot,v_slat,v_smix,v_sspec,v_ssite,v_sstr,v_sarry
Cr  2. read input data specified by tokens
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
      implicit none
#if MPI | MPIK
      include "mpif.h"
#endif
      integer,intent(in):: recln,nrecs
      character*(recln),intent(in):: recrd(nrecs)
      character,intent(in)::  prgnam*(*)
      character strn*(recln)
      integer:: i_spec
      character fileid*64
      integer:: procid,nproc,master=0
      logical lgors,cmdopt,bittst,ltmp,ioorbp!,asa
      double precision dval,dglob,xx(n0*2),dgets !,ekap(6)
      integer a2vec,bitand,fopna,getdig,i,is,igets,iprint,
     .iprt,irs(5),isw,ifi,ix(n0*nkap0),j,k,l,lbas,noinv, !noinv was iand(lqp,1)
     .lfrzw,lgunit,lmet,
c lncol,lham
     .lrs,lstsym,ltb,lxcf,nat,nlibu,nclasp,
     .nglob,nspc,nlmax,scrwid,stdo,stdl,stde,k1,k2,mpipid !nkap,lsx1,
      character*(8),allocatable::clabl(:)
      integer,allocatable:: ipc(:),initc(:),ics(:)
      real(8),allocatable:: pnuc(:,:,:),qnuc(:,:,:,:),
     .pp(:,:,:,:),ves(:),zc(:)
      integer:: dvec1(3)=1, dvec2(3)=0
C ... basis
      double precision orbp(n0,2,nkap0)
      integer o,oclabl,ohave,oics,opnu,opp,oqnu,osgw,osoptc,oves,owk !osordn,
      real(8):: pnux(20)
      integer:: nnn
      integer:: i_copy_size,i_spacks,iendx,inix,i_spackv
      real(8):: seref
      integer:: ib , ispec
      integer,allocatable:: wowk(:)
      logical:: isanrg,l_dummy_isanrg
      integer:: lmxcg,lmxcy,lnjcg,lnxcg,nlm
      
ccccccccccccccccccccccccccccccccccccccccccccccccc
      integer::nout,nn,i0,ivec(10),iosite
      integer:: io_iactive,io_tim(2),verbos(5)
      character(256)::  a,outs
      logical::  mlog=.false.
      integer:: lmxbj,lmxaj,nlbj,nlbji,nlaj,nlaji
      double precision vsn,vers,xv(2*n0)
      character(256*16) :: bigstr=' '
      integer:: it
      logical :: parmxp
c      logical,parameter:: T=.true.,F=.false.
      logical::  debug=.false.

      integer:: lp1,lpzi,nrmx=1501
c      integer,parameter:: mxspec=256
      real(8):: xxx
c      type (s_mix) :: vamix
      real(8)::  avwsr
      integer::izerv(n0)=(/(0,i=1,n0)/)
      real(8)::   zerov(n0)=(/(0d0,i=1,n0)/)

c      real(8):: tbalph(5,5)
      integer:: ii,sw
      real(8):: dasum!,dglob
      character(128) :: nm
      real(8):: nullrv(256)
      integer:: nulliv(256),jj(2) !,nkapsi
      logical:: noelind

c$$$      data tbalph/.214286D0,.000001D0,.000001D0,.000001D0,.001D-3,
c$$$     ..287234D0,.025822D0,.000001D0,.000001D0,.001D-3,
c$$$     ..348485D0,.053030D0,.010714D0,.000001D0,.001D-3,
c$$$     ..385057D0,.073209D0,.022481D0,.006069D0,.001D-3,
c$$$     ..404761D0,.087927D0,.032886D0,.012257D0,4.2735D-3/

      integer:: levelinit=0
      integer:: lx,lxx
      character*256:: sss
      logical:: sexist
      integer:: ibas,ierr,lc
      logical:: ltet
      
      procid = mpipid(1)
      nproc  = mpipid(0)
      scrwid = 80
      stdo = lgunit(1)
      stdl = lgunit(2)
      stde = stdo
      lfp = 1
      nullrv = nullr
      nulliv  =nulli
      debug = cmdopt('--debug',6,0,a)
!! --- Copy recrd to rcd in m_gtv 
      call gtv_setrcd(recrd,nrecs,recln,stdo,stdl,stde)
C --- Initialize ---
      call toksw_init(debug)
      if (cmdopt('--show',6,0,a))  io_show = 1 !iosw=3 in old code
      if (procid .ne. master) io_show = 0
      if (cmdopt('--input',7,0,a)) io_help = 1 !optio=0 in old code
      call gtv_setio(debug,io_show,io_help)
      if (io_help .eq. 1) then
        write(stdo,332)
  332   format(/' Token',t19,'Input   cast  (size,min)'/
     .  ' ------------------------------------------')
      elseif (io_show .ne. 0) then
        write(stdo,333)
  333   format(/' Token',t19,
     .  'Input   cast  (size,min,read,def)     result')
      endif
C     --- IO ---
      nm='IO_SHOW'; call gtv(trim(nm),tksw(prgnam,nm),io_show,def_i4=
     .     0,note='Echo data as it is read from input file')
      if (cmdopt('--show',6,0,a)) io_show = 1
      if (procid .ne. master) io_show = 0
      if (cmdopt('--show=',7,0,a)) then
         i = 7
         i = a2vec(a,len(a),i,2,' ',1,1,1,it,io_show)
      endif
      call gtv_setio(debug,io_show,io_help) ! In case io_show changed
      nm='IO_HELP'; call gtv(trim(nm),tksw(prgnam,nm),io_help,
     .     def_i4=0, note='Show what input would be sought, '//
     .     'without attempting to read data')
      if (cmdopt('--input',7,0,a)) io_help = 1 !optio=0 in old code
      if (io_help == 1) io_show = 1
      call gtv_setio(debug,io_show,io_help) ! In case io_help changed

!     !--- CONST --- for backword compatibility. We will remove this.
      call numsyv(nvario)
      nm='CONST'; call gtv(trim(nm),tksw(prgnam,nm), bigstr, note=
     .     'Constants may declared for use in expressions'//
     .     '%N   Variables may also be set from the command-line:  -vnam=#',
     .     nout=nout)
      if (nout .eq. 1) then
         i = 0
         call parsyv(bigstr,len_trim(bigstr),1999,0,i)
         if (io_show .ne. 0) call shosyv(0,0,0,stdo)
      endif

!     ! begining of log file.
      if (stdl>0 .and. iprint()>0) call headl2(prgnam,0,stdl)

!     ! --- Header ---
      nm='HEADER'; call gtv(trim(nm),tksw(prgnam,nm),header,note=
     .     'Contents displayed at beginning of program execution',
     .     nout=nout)
      if (nout .gt. 0) then
         outs  = 'HEADER '//trim(header)
         if (procid .eq. master) then
            if (io_show .ne. 0) write(stdo,'(1x)')
            write(stdo,'(1x,a)') trim(outs)
            if (io_show .ne. 0) write(stdo,'(1x)')
            if (stdl>0 .and. iprint()>0) write(stdl,'(1x,a)') trim(outs)
         endif
      endif

C     --- IO ---
      nm='IO_VERBOS'; call gtv(trim(nm),tksw(prgnam,nm),verbos,
     .     note='Verbosity stack for printout'//
     .     '%N   May also be set from the command-line: --pr#1[,#2]',
     .     def_i4v=(/30/),nout=i0)
C     ... Override verbos w/ -pr commmand-line arg
      if ((cmdopt('--pr',4,0,a) .or. cmdopt('-pr',3,0,a)) ) then
         i = 4
         if (cmdopt('-pr',3,0,a)) i = 3
         i = a2vec(a,len(a),i,2,', ',2,2,5,it,verbos)
         if (i .lt. 0) call rxs('error parsing switch ',a)
         i0 = max(i0,i)
      endif
C     Copy verbosities to print stack
      if (i0 >= 1) then
         do  i = 0, 5
            call sprt(i,verbos(min(i+1,i0)))
         enddo
      endif
      if (procid .ne. master) then
         call pshpr(0)
         call pshpr(0)
         do  i = 1, 4
            call sprt(i,0)
         enddo
      endif
      nm='IO_IACTIV'; call gtv(trim(nm),tksw(prgnam,nm),
     .     io_iactive, note='Turn on interactive mode'//
     .     '%N   May also be controlled from the command-line:'//
     .     '  --iactiv  or  --iactiv=no',
     .     def_i4=0)
      if( cmdopt('--no-iactiv ',7,0,a)) io_iactive = 0
      if( cmdopt('--iactiv',7,0,a))     io_iactive = 1
      if( cmdopt('--iactiv=no',11,0,a)) io_iactive = 0
      if (io_iactive .gt. 1) io_iactive = 1
      call initqu(io_iactive)
      nm='IO_TIM'; call gtv(trim(nm),tksw(prgnam,nm),io_tim,
     .     note='Turns CPU timing log.  Value sets tree depth.'//
     .     '%N   Optional 2nd arg prints CPU times as routines execute.'//
     .     '%N   Args may be set through command-line: --time=#1[,#2]',
     .     def_i4v=(/1,1/),nmin=1,nout=i0)
      if(i0==1) io_tim(2)=io_tim(1)
C     Override with '--time=' commmand-line arg
      if ( cmdopt('--time=',7,0,a) ) then
         i = 7
!     !  '--time did not work because a2vec did not work. Set IO_TIM=3,5 for example.'
         i = a2vec(a,len(a),i,2,', ',2,2,2,it,io_tim) !this does not work...
         if (i .lt. 0) call rxs('error parsing switch ',a)
         i0 = max(i0,i)
      endif
      if ( i0 >=1 ) call tcinit(io_tim(2),io_tim(1),levelinit)

!     ! --- Struc ---
      if (tksw(prgnam,'STRUC') == 2) goto 59
      if (io_show+io_help/=0) call info0(1,1,
     .     0,' --- Parameters for crystal structure ---')
      nm='STRUC_FILE'; call gtv(trim(nm),tksw(prgnam,nm),outs,nmin=10,
     .     nout=nout, note=
     .     'Name of site file containing basis and lattice information.'//
     .     '%N%3fRead NBAS, PLAT, and optionally ALAT from site file, '//
     .     'if specified.'//
     .     '%N%3fOtherwise, they are read from the ctrl file.')
      if (nout == 1) then
         xxx = nullr
         if (procid .eq. master) then
            j = iosite(7000,3d0,0,trim(outs),i,slabl,alat,plat,nbas,
     .           nspec,xxx,xxx,xxx,xxx,xxx,xxx,xxx)
         endif
         call mpibc1(nbas,1,2,mlog,'readctrl','nbas')
         call mpibc1(alat,1,4,mlog,'readctrl','alat')
         call mpibc1(plat,9,4,mlog,'readctrl','plat')
      endif
      if (alat .eq. NULLR) then
         nm='STRUC_ALAT'; call gtv(trim(nm),tksw(prgnam,nm),
     .        alat, note= 'Scaling of lattice vectors, in a.u.')
      endif
      if (nbas .eq. NULLI) then
         nm='STRUC_NBAS'; call gtv(trim(nm),tksw(prgnam,nm),
     .        nbas,note='Size of basis')
         nm='STRUC_PLAT'; call gtv(trim(nm),tksw(prgnam,nm),
     .        temp33, nmin=9, nout=nout, note=
     .        'Primitive lattice vectors, in units of alat')
         plat(1:3,1)=temp33(1:3)
         plat(1:3,2)=temp33(4:6)
         plat(1:3,3)=temp33(7:9)
      endif
      avw = avwsr(plat,alat,vol,nbas)
      nm='STRUC_NSPEC'; call gtv(trim(nm),tksw(prgnam,nm),nspec,
     .     note='Number of species to read from SPEC category.'//
     .     '%N%3fIf not present, NSPEC will be obtained '//
     .     'by counting entries in the SPEC category' )
      if (io_help == 0 .and. nspec .eq. NULLI) then
!!    ... Count number of species in SPEC category
         nm='SPEC_ATOM'; sw = tksw(prgnam,nm)
         if (sw .ne. 2) then
            j = 0; nspec = 0
            do  while (nspec .le. 0)
               j = j+1; jj= (/1,j/)
               if (.not. debug) call pshpr(0)
               call gtv(trim(nm),0,nono,Texist=ltmp,cindx=jj)
               if (.not. debug) call poppr
               if (.not. ltmp) nspec = j-1
            enddo
            if (io_show>0) call info2(1,0,0,
     .           ' ... found %i species in SPEC category',nspec,0)
         endif
      endif
C     Extra site positions for e.g. point multipoles
      nm='STRUC_NBASP'; call gtv(trim(nm),tksw(prgnam,nm),nbasp,
     .     def_i4=nbas,note='nbas + no. of point multipoles')
      if (nbasp .eq. NULLI) nbasp = nbas
      nsite = nbasp
      nm='STRUC_SLAT'; call gtv(trim(nm),tksw(prgnam,nm),temp33,
     .     nmin=9,note='Supercell lattice vectors')
      lat_slat(1:3,1)=temp33(1:3)
      lat_slat(1:3,2)=temp33(4:6)
      lat_slat(1:3,3)=temp33(7:9)
      nm='STRUC_DALAT'; call gtv(trim(nm),tksw(prgnam,nm),
     .     dalat, def_r8=0d0, note='added to alat after input is read')
      nm='STRUC_NL'; call gtv(trim(nm),tksw(prgnam,nm),nl,def_i4=3,
     .     note='global default LMXA+1 for basis and augmentation.')
C     ... Lattice distortion or rotation
      sw = tksw(prgnam,'STRUC_SHEAR')
      lat_gam(1) = NULLR
      if (sw/=2) then
         j = 0
         nm='STRUC_SHEAR'; call gtv(trim(nm),0,lat_gam,nmin=4,or=T,
     .        note='Volume-conserving shear of PLAT (1=ideal)')
         if (lat_gam(1) /= NULLR .and. io_help == 0) then
            j = 0
            goto 880
         endif
         nm='STRUC_ROT'; call gtv(trim(nm),0,outs,nmin=10,nout=nout,or=T)
         if (nout /= 0 .and. io_help == 0) then
            call a2rotm(outs,.false.,iprint()-10,lat_dist)
            j = 2
            goto 880
         endif
         nm='STRUC_DEFGRD'; call gtv(trim(nm),0,lat_dist,nmin=9, or=T)
         if (lat_dist(1) /= NULLR .and. io_help == 0) then
            j = 2
            goto 880
         endif
         nm='STRUC_STRAIN'; call gtv(trim(nm),sw,lat_dist(1:6),nmin=6,
     .        note='6 Voigt strain matrix elements')
         if (lat_dist(1) /= NULLR .or. io_help /= 0) then
            nm='STRUC_ALPHA'; call gtv(trim(nm),1,lat_dist(7),nout=nout,
     .           note=
     .           'Amplitude of (Voigt) strain.  Only read if STRAIN input')
            if (nout /= 0) j = 3
         endif
 880     continue
         lat_ldist = j
      endif
      if (lat_gam(1) == NULLR) lat_gam=(/0d0,0d0,1d0,1d0/)
C     ... End of struc tokens
   59 continue

      
C     --- Options ---
      if (tksw(prgnam,'OPTIONS') == 2) goto 69
      if (io_show+io_help/=0) call info0(1,1,
     .     0,' --- Program Options ---')
      nm='OPTIONS_HF'; call gtv(trim(nm),tksw(prgnam,nm),lhf, !lcd2,
     .     def_lg=F,note='T for non-self-consistent Harris')
      nm='OPTIONS_RMINES'; call gtv(trim(nm),tksw(prgnam,nm),rmines,
     .     def_r8=1d0,note='Minimum MT radius when finding new ES')
      nm='OPTIONS_RMAXES'; call gtv(trim(nm),tksw(prgnam,nm),rmaxes,
     .     def_r8=2d0,note='Maximum MT radius when finding new ES')
      nm='OPTIONS_NESABC'; call gtv(trim(nm),tksw(prgnam,nm),nesabc,
     .     nmin=3,def_i4v=(/100,100,100/),
     .     note='No. divisions when searching for empty spheres')
      outs = ' '
      nm='OPTIONS_Q'; call gtv(trim(nm),tksw(prgnam,nm),outs,nmin=10,
     .     note='Use Q=show, Q=atom, or Q=band to quit after '//
     .     'input, sphere calc or band pass.'//
     .     '%N%3fCommand-line `--quit=string'' overrides input file')
      if (outs=='SHOW' .or. outs=='show' .or.
     .     cmdopt('--quit=show',11,0,a)) then; quit=1
      elseif (outs=='ATOM' .or. outs=='atom' .or.
     .        cmdopt('--quit=atom',11,0,a)) then; quit=2
      elseif (outs=='BAND' .or. outs=='band' .or.
     .        cmdopt('--quit=band',11,0,a)) then; quit=4
      elseif (outs=='HAM' .or. outs=='ham' .or.
     .        cmdopt('--quit=ham',10,0,a)) then; quit=8
      elseif (outs=='dmat' .or.
     .        cmdopt('--quit=dmat',11,0,a)) then; quit=16
      elseif (outs==' ') then
      else
         call rx('OPTIONS_Q= must contain one of SHOW, ATOM, or BAND')
      endif
      nm='OPTIONS_PFLOAT';call gtv(trim(nm),tksw(prgnam,nm),lpfloat,
     .     def_i4=1,note=       !feb2012takao
     .     'Controls how band CG is determined in floating Pnu.'//
     .     '%N%5f0 Band CG found by traditional method'//
     .     '%N%5f1 Band CG found from true energy moment of density')
      nm='OPTIONS_SCR'; call gtv(trim(nm),tksw(prgnam,nm),lscr,
     .     def_i4=0,note='Use scr to accelerate convergence:'//
     .     '%N%3f0 do nothing'//
     .     '%N%3f1 Make ASA static response function (see documentation)'//
     .     '%N%3f2 Use response to screen output q and ves'//
     .     '%N%3f4 Use model response to screen output q'//
     .     '%N%3f6 Use response to screen output ves only'//
     .     '%N%5fAdd 1 to combine mode 1 with another mode'//
     .     '%N%5fAdd 10*k to compute intra-site contribution to vbare'//
     .     ' each kth iteration'//
     .     '%N%5fAdd 100*k to compute response function'//
     .     ' on every kth iteration')
   69 continue                  ! end of OPTIONS

C     --- Hamiltonian parameters ---
      if (io_show+io_help/=0 .and. tksw(prgnam,'HAM')/=2) call info0(1,1,
     .     0,' --- Parameters for hamiltonian ---')
      nsp = 1
      nm='HAM_NSPIN'; call gtv(trim(nm),tksw(prgnam,nm),nsp,
     .     def_i4=1,note='Set to 2 for spin polarized calculations')
      if (io_help .eq. 0) then
         l_dummy_isanrg=isanrg(nsp,1,2,'rdctrl','nsp',T)
      else
         nsp = 1
      endif
c     takao !need inversion check...
c     if(prgnam=='LMCHK') then
c     write(6,*) '### ---> ENFORCE nsp=1 for LMCHK ###'
c     nsp=1
c     endif
C     Flags representing how cd is represented: default values
      lcd4=F 
      if (prgnam == 'LMF' .or. prgnam == 'LMFGWD') lcd4=T
      nm='HAM_REL'; call gtv(trim(nm),tksw(prgnam,nm),lrel,def_i4=1,
     .     note='0 for nonrelativistic Schrodinger equation'//
     .     '%N%3f1 for scalar relativistic Schrodinger equation'//
     .     '%N%3f2 for Dirac equation')
C     Fully relativistic => spin-orbit coupling
c     if (lrel==2) lncol4 = T
      if (lrel==2) lso=1
      if (nsp==2 .or. io_help/=0) then
         if (io_help /= 0) call info0(1,0,0,' * To read the magnetic '//
     .        'parameters below, HAM_NSPIN must be 2')
         nm='HAM_SO'; call gtv(trim(nm),tksw(prgnam,nm),
     .        lso,def_i4=0,note=
     .        'Spin-orbit coupling (for REL=1)'//
     .        '%N%3f0 : no SO coupling'//
     .        '%N%3f1 : Add L.S to hamiltonian'//
     .        '%N%3f2 : Add Lz.Sz only to hamiltonian'//
     .        '%N%3f3 : Like 2, but also compute <L.S-LzSz> by perturbation')
         if (io_help==0) l_dummy_isanrg=isanrg(lso,0,3,' rdctrl:','SO',T)
      endif
      sw = tksw(prgnam,'HAM_GMAX')
      if (sw/=2) then
         nm='HAM_GMAX'; call gtv(trim(nm),sw,lat_gmax,nmin=1,nout=nout,
     .        note='Energy cutoff for plane-wave mesh',or=T)
         if (nout .ne. 0) then
            sw = 2
         else
            lat_gmax = 0
         endif
         nm='HAM_FTMESH'; call gtv(trim(nm),sw,ftmesh,nout=nout,
     .        note='No. divisions for plane-wave mesh '//
     .        'along each of 3 lattice vectors.'//
     .        '%N%3fSupply one number for all vectors or a separate '//
     .        'number for each vector.')
         call fill3in(nout,ftmesh)
      endif
      nm='HAM_TOL'; call gtv(trim(nm),tksw(prgnam,nm),tolft,
     .     def_r8=1d-6, note='w.f. tolerance for FT mesh')
      nm='HAM_FRZWF'; call gtv(trim(nm),tksw(prgnam,nm),frzwf,def_lg=F,
     .     note='Set to freeze augmentation wave functions for all species')
      nm='HAM_FORCES'; call gtv(trim(nm),tksw(prgnam,nm),ctrl_lfrce,
     .     def_i4=0,note=
     .     'Controls the ansatz for density shift in force calculation.'//
     .     '%N%3f-1 no force%3f0 no shift'//
     .     '%N%3f 1 free-atom shift  12 screened core+nucleus')
      if(.not.noelind()) then
         nm='HAM_ELIND'; call gtv(trim(nm),tksw(prgnam,nm),elind,
     .        def_r8=0d0,note='Lindhard energy for model screening')
      else
         elind=0d0
      endif
      nm='HAM_XCFUN'; call gtv(trim(nm),tksw(prgnam,nm),ham_lxcf,def_i4=2,
     .     note='Specifies local exchange correlation functional:'//
     .     '%N%3f1 for Ceperly-Alder (VWN)'//
     .     '%N%3f2 for Barth-Hedin (ASW fit)'//
     .     '%N%3f103 for PBE-GGA (use xcpbe.F in ABINIT')
      nm='HAM_GGA'; call gtv(trim(nm),tksw(prgnam,nm),ggax,
     .     nout=nout, def_i4=0, note='UnUsed now. See XCFUN...')
      if(ggax>0) call rx('Use HAM_XCFUN instead of HAM_GGA')
      nm='HAM_ALFSI'; call gtv(trim(nm),tksw(prgnam,nm),
     .     alfsi, def_r8=0d0,note=
     .     'Coefficient to artificial addition of overlap to hamiltonian')
      nm='HAM_DABC' ; call gtv(trim(nm),tksw(prgnam,nm),
     .     dabc,nout=nout,
     .     note='Spacings for real-space interstital mesh')
      if (nout == 1) dabc(2) = dabc(1)
      if (nout <  3) dabc(3) = dabc(2)
!     !    Merge with ?
      nm='HAM_DQVAL'; call gtv(trim(nm),tksw(prgnam,nm),dqval,
     .     note='Total charge')
      nm='HAM_RDSIG'; call gtv(trim(nm),tksw(prgnam,nm),lrsig,def_i4=0,
     .     note=
     .     'Controls how self-energy is added to '//
     .     'local exchange correlation functional:'//
     .     '%N%3f   0 do not read Sigma'//
     .     '%N%3f  12 read ESEAVR as Sigma-Vxc')
      nm='HAM_RSRNGE'; call gtv(trim(nm),tksw(prgnam,nm),rsrnge,
     .     def_r8=5d0, note='Maximum range in connecting vectors for '//
     .     'r.s. sigma (units of alat)')
      nm='HAM_RSSTOL'; call gtv(trim(nm),tksw(prgnam,nm),
     .     rsstol, def_r8=5d-6, note=
     .     'Max tolerance in Bloch sum error for r.s. sigma ')
      nm='HAM_ScaledSigma'; call gtv(trim(nm),tksw(prgnam,nm),scaledsigma,
     .     def_r8=1d0, note='=\alpha_Q for QSGW-LDA hybrid. \alpha \times (\Sigma-Vxc^LDA) is added to LDA/GGA Hamiltonian.')
      nm='HAM_EWALD'; call gtv(trim(nm),tksw(prgnam,nm),ham_ewald,
     .     def_lg=.false.,note='Make strux by Ewald summation')
      nm='HAM_VMTZ'; call gtv(trim(nm),tksw(prgnam,nm),vmtz,def_r8=0d0,
     .     note='Muffin-tin zero defining wave functions')
      nm='HAM_PMIN'; call gtv(trim(nm),tksw(prgnam,nm),pmin,
     .     def_r8v=zerov,nout=nout,note=
     .     'Global minimum in fractional part of P-functions.'//
     .     '%N%3fEnter values for l=0..nl:'//
     .     '%N%3f0: no minimum constraint'//
     .     '%N%3f#: with #<1, floor of fractional P is #'//
     .     '%N%3f1: use free-electron value as minimum')
      nm='HAM_PMAX'; call gtv(trim(nm),tksw(prgnam,nm),
     .     pmax, def_r8v=zerov, nout=nout, note=
     .     'Global maximum in fractional part of P-functions.'//
     .     '%N%3fEnter values for l=0..nl:'//
     .     '%N%3f0: no maximum constraint'//
     .     '%N%3f#: with #<1, ceiling of fractional P is #')
!     ! We set default oveps=1d-7 16Nov2015. This was zero before the data.
      nm='HAM_OVEPS'; call gtv(trim(nm),tksw(prgnam,nm),
     .     oveps, def_r8=1d-7, nout=nout, note=
     .     'Diagonalize hamiltonian in reduced hilbert space,'//
     .     '%N%3fdiscarding part with evals of overlap < OVEPS')
      nm='HAM_STABILIZE'; call gtv(trim(nm),tksw(prgnam,nm),
     .     delta_stabilize, def_r8=-1d0, nout=nout,note=
     .     'Experimental. Stabilizer for Diagonalize hamiltonian (negative means unused),'//
     .     '%N%3f "H --> H + HAM_STABILIZE*O^-1" in zhev_tk(diagonalization)')
C     ... APW basis
      nm='HAM_PWMODE'; call gtv(trim(nm),tksw(prgnam,nm),pwmode,
     .     def_i4=0,note=
     .     'Controls APW addition to LMTO basis'//
     .     '%N%3f1s digit:'//
     .     '%N%6f0: LMTO basis only'//
     .     '%N%6f1: Mixed LMTO+PW'//
     .     '%N%6f2: PW basis only'//
     .     '%N%3f10s digit:'//
     .     '%N%6f0: PW basis fixed'//
     .     '%N%6f1: PW basis q-dependent')
      if(pwmode==10) pwmode=0   !takao added. corrected Sep2011
      nm='HAM_PWEMIN'; call gtv(trim(nm),tksw(prgnam,nm),
     .     pwemin, def_r8=0d0, nout=nout, note=
     .     'Include APWs with energy E > PWEMIN (Ry)')
      nm='HAM_PWEMAX'; call gtv(trim(nm),tksw(prgnam,nm),
     .     pwemax, def_r8=0d0, nout=nout, note=
     .     'Include APWs with energy E < PWEMAX (Ry)')
      nm='HAM_NPWPAD'; call gtv(trim(nm),tksw(prgnam,nm),
     .     npwpad, def_i4=-1, note=
     .     'Overrides default padding of variable basis dimension')
      nm='HAM_UDIAG'; call gtv(trim(nm),tksw(prgnam,nm),ham_udiag,
     .     def_i4=0,note='nonzero => diagonal-only LDA+U')
     
C     --- Symmetry group ---
      if (io_show+io_help/=0 .and. tksw(prgnam,'SYMGRP')/=2)
     .     call info0(1,1,0,' --- Symmetry group operations ---')
      nm='SYMGRP'; call gtv(trim(nm),tksw(prgnam,nm),symg,
     .     note='Generators for symmetry group')
C     --- Symmetry group for AF --- !june2015
      if (io_show+io_help/=0 .and. tksw(prgnam,'SYMGRPAF')/=2)
     .     call info0(1,1,0,' --- Symmetry group operations ---')
      nm='SYMGRPAF'; call gtv(trim(nm),tksw(prgnam,nm),symgaf,
     .     note='An (or multiple) Extra Generator for adding anti ferro symmetry')
C     --- Species (old CLASS) ---
      if (io_help == 1) nspec = 1
      if (tksw(prgnam,'SPEC') == 2) goto 79
      if (io_show+io_help/=0) call info0(1,1,
     .     0,' --- Parameters for species data ---')
      if (io_help /= 0) call info0(1,0,0,' * The next four tokens '//
     .     'apply to the automatic sphere resizer')
      nm='SPEC_SCLWSR'; call gtv(trim(nm),tksw(prgnam,nm),
     .     sclwsr, def_r8=0d0, note=
     .     'Scales sphere radii, trying to reach volume = '//
     .     'SCLWSR * cell volume'//
     .     '%N%3fSCLWSR=0 turns off this option.'//
     .     '%N%3fAdd  10  to initially scale non-ES first;'//
     .     '%N%3f or  20  to scale ES independently.')
      nm='SPEC_OMAX1'; call gtv(trim(nm),tksw(prgnam,nm),omax1,
     .     def_r8v=(/0d0,0d0,0d0/),note=
     .     'Limits max sphere overlaps when adjusting MT radii')
      nm='SPEC_OMAX2'; call gtv(trim(nm),tksw(prgnam,nm),
     .     omax2, def_r8v=(/0d0,0d0,0d0/),note=
     .     'Sphere overlap constraints of second type',nout=nout)
      nm='SPEC_WSRMAX'; call gtv(trim(nm),tksw(prgnam,nm),wsrmax,
     .     def_r8=0d0,note=
     .     'If WSRMAX is nonzero, no sphere radius may exceed its value')
      if (io_help .eq. 1) then
         write(*,382)
 382     format(/' * ',
     .        'The following tokens are input for each species. ',
     .        'Data sandwiched'/3x,'between successive occurences of ',
     .        'token ATOM apply to one species.')
         nspec = 1
      endif

      if (nspec .eq. 0) goto 79

      allocate(pnu(n0,nsp,nspec),qnu(n0,nsp,nspec),
     .     pz(n0,nsp,nspec),amom(n0,nspec),idmod(n0,nspec),
     .     rsmh(n0,nspec),eh(n0,nspec),rsmh2(n0,nspec),eh2(n0,nspec),
     .     pb1(nspec),pb2(nspec),lmxpb(nspec),
     .     ehvl(n0,nspec),
     .     qpol(n0,nspec),stni(nspec),tbvso(4,nspec),
     .     iq1(n0,nspec),ivso(n0,nspec),
     .     rg(nspec),rsma(nspec),rfoca(nspec),rsmfa(nspec),rcfa(2,nspec),
     .     rs3(nspec),rham(nspec),rmt(nspec),rsmv(nspec),
     .     nxi(nspec),exi(n0,nspec),rint(nspec),rcut(nspec),
     .     spec_a(nspec),z(nspec),nr(nspec),mass(nspec),eref(nspec),
     .     coreh(nspec),coreq(2,nspec),
     .     colxbs(3,nspec),radxbs(nspec),
     .     idxdn(n0,nspec),
     .     hcr(n0,nspec),       !alpha(n0,nspec),
     .     idu(4,nspec),uh(4,nspec),jh(4,nspec),
     .     dv(nspec),grp(nspec),grp2(nspec),
     .     mxcst1(nspec),mxcst2(nspec),mxcst4(nspec),
     .     kmxt(nspec),kmxv(nspec),
     .     lfoca(nspec),lmxl(nspec),lxi(nspec),lmxa(nspec),lmxb(nspec),nmcore(nspec))

      mxcst1=F
      mxcst2=F
      mxcst4=F

      nkaph = 1
      lpzi = 1
      qpol = NULLR
      rsmh = 0d0
      rsmh2 = 0d0
      eh = NULLR
      eh2 = NULLR
      hcr = NULLR
      idmod = NULLI
      ehvl = NULLR
      do  j = 1, nspec
        if(debug) print *,'nspec mxcst j-loop j nspec',j,nspec
        colxbs(:,j) = NULLR; radxbs(j) = NULLR
        rcfa(:,j) = NULLR; rfoca(j) = 0d0; rg(j) = 0d0
        rham(j) = NULLR; rsma(j) = 0d0; rsmfa(j) = 0d0
        spec_a(j) = NULLR; nr(j) = NULLI
        exi(:,j) = NULLR
        coreh(j) = ' '; coreq(:,j) = NULLR
        eref(j) = NULLR
        mass(j) = NULLR

        if (io_help /= 0) then
          write(stdo,'(1x)')
        elseif (io_help == 0 .and. io_show>0) then
          call info(1,0,0,' ... Species %i',j,0)
        endif

        jj= (/1,j/)
        nm='SPEC_ATOM'; call gtv(trim(nm),tksw(prgnam,nm),slabl(j),
     .       nmin=10,cindx=jj,note='Species label')

        nm='SPEC_ATOM_Z'; call gtv(trim(nm),tksw(prgnam,nm),z(j),
     .       cindx=jj,note='Atomic number')

        if(debug) print *,'nspec aaa mxcst j-loop j nspec',j,nspec

        sw = tksw(prgnam,'SPEC_ATOM_R')
        if (sw /= 2) then
           nout = 0
           nm='SPEC_ATOM_R'; call gtv(trim(nm),sw,rmt(j),cindx=jj,
     .          nout=nout,note= 'Augmentation sphere radius rmax',or=T)
           if (nout .ne. 1) then
              nm='SPEC_ATOM_R/W';call gtv(trim(nm),sw,rmt(j),cindx=jj,
     .             nout=nout,note='rmax relative to average WS radius',or=T)
              if (nout == 1) then
                 rmt(j) =rmt(j)*avw
              else
                 nm='SPEC_ATOM_R/A';call gtv(trim(nm),sw,rmt(j),cindx=jj,
     .                nout=nout,note='rmax relative to lattice constant')
                 if (nout == 1) then
                    rmt(j) =rmt(j)*alat
                 else
c     takao for lmchk even when R is not given. See default of LMCHK
                    rmt(j) = 0d0
                 endif
              endif
           endif
        endif
c     if (ltbe) rmt(j) = 1d0

        if(debug) print *,'nspec bbb mxcst j-loop j nspec',j,nspec

C     ... Radial mesh parameters: determine default value of a
        i0 = NULLI
        xxx = NULLR
        if (io_help .ne. 1) then
           call pshpr(0)
           call rmesh(z(j),rmt(j),lrel,.false.,nrmx,xxx,i0)
           call poppr
           if (xxx .eq. .03d0) xxx = .015d0 !.025d0 jun2012 .025 to .015 as default.
        endif
        nm='SPEC_ATOM_A'; call gtv(trim(nm),tksw(prgnam,nm),spec_a(j),
     .       def_r8=xxx,cindx=jj,nout=nout,
     .       note='Radial mesh point spacing parameter')
C     Determine default NR
        if (tksw(prgnam,'SPEC_ATOM_NR') .ne. 2) then
           i0 = 0
           call pshpr(0)
           call rmesh(z(j),rmt(j),lrel,.false.,nrmx,spec_a(j),i0)
           call poppr
        endif
        nm='SPEC_ATOM_NR'; call gtv(trim(nm),tksw(prgnam,nm),nr(j),
     .       def_i4=i0,cindx=jj, note='Number of radial mesh points')
        if (nr(j) .eq. 0) nr(j) = i0
        nm='SPEC_ATOM_LMX'; call gtv(trim(nm),tksw(prgnam,nm),lmxb(j),
     .       def_i4=max(nl-1,NULLI),cindx=jj,note='l-cutoff for basis')
C     Running account of maximum lmxb
        lmxbj = NULLI
        if (io_help .ne. 1) then
           lmxbx = max(lmxbx,lmxb(j))
           lmxbj = lmxb(j)
        endif
C     nlbj = number of elements associated with lmxb
C     0 => no elements
C     nlbji: ditto, but used to specify number of default values
        nlbj = 1+lmxbj
        nlbji = nlbj
        if (io_help .eq. 1) then
           nlbji = NULLI
           nlbj = 1
        elseif (lmxbj == NULLI) then
           nlbji = NULLI
           nlbj = 0
        endif
        if(debug) print *,'nspec ccc mxcst j-loop j nspec',j,nspec

C     ... Basis set for lmf
        nm='SPEC_ATOM_RSMH'; call gtv(trim(nm),tksw(prgnam,nm),
     .       rsmh(1:nlbj,j),cindx=jj,nout=nout,nmin=nlbj,def_r8v=zerov,
     .       note='Smoothing radii for basis')
        nn = NULLI
        if (nout .gt. 0) nn = nout
C     Find nn=lmax+1 for which rsmh ne 0 => reduce # EH required
        do  i = 1, nlbj
           if (rsmh(i,j) .eq. 0) cycle
           nn = i
        enddo
        nm='SPEC_ATOM_EH'; call gtv(trim(nm),tksw(prgnam,nm),
     .       eh(1:nlbj,j),nmin=nn,cindx=jj,
     .       note='Kinetic energies for basis')
        nn = NULLI
        if (nout>0) nn = nout
        nm='SPEC_ATOM_RSMH2'; call gtv(trim(nm),tksw(prgnam,nm),
     .       rsmh2(1:nlbj,j),nmin=nlbj,def_r8v=zerov,cindx=jj,nout=nout,
     .       Texist=ltmp,note='Basis smoothing radii, second group')
        if (ltmp) then
           nn = NULLI
           sw = tksw(prgnam,nm)
           if (nout>0) then
              nkaph=2
              nn = nout
              sw = 1
C     Find nn=lmax+1 for which rsmh ne 0 => reduce # EH required
              do  i = 1, nlbj
                 if (rsmh2(i,j) .eq. 0) cycle
                 nn = i
              enddo
           endif
           nm='SPEC_ATOM_EH2'; call gtv(trim(nm),sw,
     .          eh2(1:nlbj,j),nmin=nn,cindx=jj,
     .          note='Basis kinetic energies, second group')
        endif
        if(debug) print *,'nspec ddd mxcst j-loop j nspec',j,nspec
        nm='SPEC_ATOM_EHVL'; call gtv(trim(nm),tksw(prgnam,nm),
     .       ehvl(1:nlbj,j),cindx=jj,def_r8v=(/(-0.5d0,i=1,n0)/),
     .       nout=nout,nmin=nlbj, note='val-lap fit energies')

C     ... Determine lmxa: floating orbitals have no lmxa
        nm='SPEC_ATOM_LMXA'; sw = tksw(prgnam,nm)
        if(debug) print *,'nspec ddd000 mxcst j rmt io_help lmxa=',j,rmt(j),io_help,lmxa(j)
        if (rmt(j) == 0.and.io_help/=1) then !takao iohelp/=1 added.
           lmxa(j) = -1
        elseif (sw == 2) then   !lmxa not read: look for subsitute
           lmxa(j) = 4
        else
C     Possibly replace with a formula that depends on rmt
           call gtv(trim(nm),sw,lmxa(j),
     .          def_i4=max(nl-1,lmxbj,NULLI),cindx=jj,Texist=ltmp,note=
     .          'l-cutoff for augmentation')
C     lmxb may not exceed lmxa
           if (io_help == 0 .and. lmxbj .gt. lmxa(j))
     .          call rx2('species '//trim(slabl(j))//' : LMX=%i '//
     .          'exceeds LMXA=%i.  Revise input so that LMX<=LMXA.',
     .          lmxbj,lmxa(j))
        endif
        lmxaj = lmxa(j)
C     nlaj = number of elements associated with lmxa
C     0 => no elements
C     nlaji: ditto, but used to specify number of default values
        nlaj = 1+lmxaj
        if(debug) print *,'nspec ddd111 mxcst j-loop j nspec nlaj lmxaj=',j,nspec,nlaj,lmxaj
        nlaji = nlaj
        if (io_help .eq. 1) then
           nlaji = NULLI
           nlaj = 1
        elseif (lmxaj == NULLI) then
           nlaji = NULLI
           nlaj = 0
        endif
C     if (lmxaj == NULLI) nlaj = NULLI
C     if (lmxaj == NULLI) nlaji = NULLI


C     ... Parameters that depend on the existence of an augmentation sphere
C     lmxl = l-cutoff for numerical rep'sn of density in sphere
C     lfoca = mode for treating core
C     kmxt = kmax for expansion of envelope wf tails
C     kmxv = cutoff to expand smoothed potential
C     hardwired for now
        kmxv(j) = 15            !Not input
C     Cannot set default here: need set after rescaling of rmt
C     rsmv(j) = rmt(j)*.5d0   !Not input
        rsmv(j) = 0d0           !Not input
        kmxt(j) = -1            !If sought, default will be set below
        lfoca(j) = 0            !If sought, default will be reset below
        nmcore(j)=0
        lmxl(j) =  lmxaj        !Use lmxaj in case not sought (ASA:mpol)
        lmxpb(j) = NULLI        !If sought, default will be set below
        lxi(j) = NULLI          !If sought, default will be set below
        nxi(j) = NULLI          !If sought, default will be set below
        pb1(j) = ' '
        pb2(j) = ' '
        call dpzero(pnu(1,1,j),n0*nsp)
        call dpzero(pz(1,1,j),n0*nsp)
        rint(j) = NULLI         !No default for now
        rcut(j) = NULLI         !No default for now
        rs3(j) = NULLI          !If sought, default will be set below
        idxdn(:,j) = 1
        pnu(1,1,j) = NULLI      !If sought, default will be set below
        qnu(1,1,j) = NULLI      !If sought, default will be set below

        if(debug) print *,'nspec ddd222 mxcst j-loop j nspec nlaj=',j,nspec,nlaj

        if (nlaj .ne. 0) then
           nm='SPEC_ATOM_LMXL'; call gtv(trim(nm),tksw(prgnam,nm),lmxl(j),
     .          cindx=jj,def_i4=lmxaj,note=
     .          'lmax for which to accumulate rho,V in sphere')

C     ... Set up default P,Q in absence of explicit specification
           call dpzero(pnu(1,1,j),nsp*n0)
           call dpzero(qnu(1,1,j),nsp*n0)

!     !       This branch sets defaults beforehand ... too complicated
C     if (io_help == 0) then
C     C         Default P,Q, nsp=1
C     if (debug) call pshpr(100)
C     call defpq(z(j),lmxaj,1,pnu(1,1,j),qnu(1,1,j))
C     if (debug) call poppr
C     endif
C     nm='SPEC_ATOM_P'; call gtv(trim(nm),tksw(prgnam,nm),
C     .    pnu(1:nlaj,1,j),def_r8v=pnu(:,1,j),nmin=nlaji,cindx=jj,note=
C     .    'Starting log der. parameters for each l')
C     nm='SPEC_ATOM_Q'; call gtv(trim(nm),tksw(prgnam,nm),
C     .    qnu(1:nlaj,1,j),def_r8v=qnu(:,1,j),nmin=nlaji,cindx=jj,note=
C     .    'Starting sphere charges for each l channel')
C     call snit

!     ! -- takao move back default value of dev_r8v to zero june2012 --
           nm='SPEC_ATOM_P'; call gtv(trim(nm),tksw(prgnam,nm),
     .          pnu(1:nlaj,1,j),def_r8v=zerov,cindx=jj,note=
     .          'Starting log der. parameters for each l')
           nm='SPEC_ATOM_Q'; call gtv(trim(nm),tksw(prgnam,nm),
     .          qnu(1:nlaj,1,j),def_r8v=zerov,cindx=jj,note=
     .          'Starting valence charges for each l channel.'//
     .          '%N%2f Q do not include semicore(PZ) electrons.'//
     &          '%N%2f Charge configuration is shown by lmfa %N'//
     &          '%N%2f WARN: This version cannot treat two valence channels'//
     &          '%N%2f per l (Q for a l-channl is zero if the l is with PZ).'//
     &          '%N%2f This causes a problem typically in Li; then we '//
     &          '%N%2f can not treat both of PZ=1.9 and P=2.2 as valence.'//
     &          '%N%2f To avoid this, use Q=0,1 together.'//
     &          ' This trick supply an '//
     &          '%N%2f electron to 2p channel; this trick works fine.')

!     ! ==== Reset default P,Q in absence of explicit specification ====
           if (io_help == 0) then
              if (io_show .ne. 0) call pshpr(50)

!     ! -- takao jun2012. qnu is set by default p. --
!     ! This looks too complicated. Fix this in future.
!     ! In anyway, we expect pnu and qnu are correctly returned (qnu does not care value of given P).
c     print *,'qnuin ',sum(abs(qnu(:,:,j))),qnu(:,:,j)
!     ! set default pnu. See the following section 'correct qnu'
!     ! isp=1 means charge. isp=2 means mmom
              if(allocated(pnudefault)) deallocate(pnudefault,qnudefault,qnudummy)
              allocate(pnudefault(n0,nsp),qnudefault(n0,nsp),qnudummy(n0,nsp))
              pnudefault=0d0
              qnudefault=0d0
              qnudummy=0d0
              iqnu=1
              if(sum(abs(qnu(:,1,j)))<1d-8) iqnu=0 !check initial Q is given or not.
              call defpq(z(j),lmxaj,1,pnudefault,qnudefault) ! qnu is given here for default pnu.
              call defpq(z(j),lmxaj,1,pnu(1,1,j),qnudummy) ! set pnu. qnu is kept (but not used here).
              if(iqnu==0) qnu(:,1,j)=qnudefault(:,1)
!     !

C     call snit
              if (io_show .ne. 0) call poppr
           endif
           if (nsp .eq. 2) call dcopy(n0,pnu(1,1,j),1,pnu(1,2,j),1)
           if (nsp .eq. 2 .or. io_help == 1) then
              nm='SPEC_ATOM_MMOM'; call gtv(trim(nm),tksw(prgnam,nm),
     .             qnu(1:nlaj,2,j),def_r8v=zerov,nmin=nlaji,cindx=jj,note=
     .             'Starting mag. moms for each l channel.'//
     .             '%N%2f For a chanel with PZ, this is enforced to be zero.'//
     .             '%N%2f See explanation for SPEC_ATOM_Q.')
           endif

           nm='SPEC_ATOM_NMCORE'; call gtv(trim(nm),tksw(prgnam,nm),nmcore(j), !takao jun2012
     .          def_i4=0,cindx=jj,note='spin-averaged core: jun2012takao'//
     .          '%N%3f0(default): spin-polarized core'//
     .          '%N%3f1         : spin-averaged core density is from spin-averaged potential')

           nm='SPEC_ATOM_PZ'; call gtv(trim(nm),tksw(prgnam,nm),
     .          pz(1:nlaj,1,j),def_r8v=zerov,nmin=nlaji,cindx=jj,note=
     .          'Starting semicore log der. parameters'//
     .          '%N%10fAdd 10 to attach Hankel tail'//
     .          '%N%10fAdd 20 to include perturbatively',nout=nout)
           if (nout>0) then
              if (dasum(nlaj,pz(1,1,j),1) /= 0) lpzi = max(lpzi,2)
           endif

!     ! correct qnu jun2012takao  2012july->mod(int(pz...,10)
!     ! our four cases are
!     !  P=Pdefault      ! qnu
!     !  Pdefault < P    ! Pdefault is filled as core
!     !  Pz < P=Pdefault ! qnu + 2*(2l+1)
!     !  Pz=Pdefault < P ! qnu
           if(iqnu==0) then
              do lx=0,lmxaj     !correct valence number of electrons.
                 if(pz(lx+1,1,j)<1d-8) then ! PZ(local orbital) not exist
                    if( int(pnudefault(lx+1,1)) < int(pnu(lx+1,1,j)) ) then
                       qnu(lx+1,1,j)= 0d0 ! pnudefault is filled and no q for pnu. (core hole case or so)
                    endif
                 else           !PZ exist
c     print *,'qnu=',lx,qnu(lx+1,1,j)
                    if( mod(int(pz(lx+1,1,j)),10)<int(pnudefault(lx+1,1)) ) then
                       qnu(lx+1,1,j)= qnu(lx+1,1,j)+ 2d0*(2d0*lx+1d0)
                    endif
                 endif
              enddo
           endif

           i0 = 1
           if (z(j) .le. 8) i0 = 0
           if (io_help .eq. 1) i0 = NULLI
           nm='SPEC_ATOM_LFOCA';call gtv(trim(nm),tksw(prgnam,nm),lfoca(j),
c     .    def_i4=i0,cindx=jj,note='FOCA switch 0, 1 or 2 (see docs)')
     .          def_i4=i0,cindx=jj,note='FOCA switch 0(within MT):'
     .          //'=1(frozenCore). Default: 1 for z>8;0 for z<=8') !takao Aug2010
           if(io_help==0) then
              if(lfoca(j)/=0.and.lfoca(j)/=1) then
                 call rx('LFOCA should be 0 or 1 (Aug2010): 2 is not allowed')
              endif
           endif

C     call snot
           nm='SPEC_ATOM_KMXA'; call gtv(trim(nm),tksw(prgnam,nm),
     .          kmxt(j),def_i4=3,cindx=jj,note=
     .          'k-cutoff for projection of wave functions in sphere.')
C     Cannot set default here: need set after rescaling of rmt
C     xxx = 0.40d0*rmt(j); if (io_help .eq. 1) xxx = NULLI
           xxx = 0d0; if (io_help .eq. 1) xxx = NULLI
           nm='SPEC_ATOM_RSMA'; call gtv(trim(nm),tksw(prgnam,nm),rsma(j),
     .          def_r8=xxx,cindx=jj,note=
     .          'Smoothing for projection of wave functions in sphere.'//
     .          '%N%3finput<0 => choose default * -input')

C     Cannot set default here: need set after rescaling of rmt
C     xxx = 0.25d0*rmt(j); if (io_help .eq. 1) xxx = NULLI
           xxx = 0d0; if (io_help .eq. 1) xxx = NULLI
           nm='SPEC_ATOM_RSMG'; call gtv(trim(nm),tksw(prgnam,nm),rg(j),
     .          def_r8=xxx,cindx=jj,note=
     .          'Smoothing for projection of charge in sphere.'//
     .          '%N%3finput<0 => choose default * -input')
C     Cannot set default here: need set after rescaling of rmt
C     xxx = 0.40d0*rmt(j); if (io_help .eq. 1) xxx = NULLI
           xxx = 0d0; if (io_help .eq. 1) xxx = NULLI
           nm='SPEC_ATOM_RFOCA'; call gtv(trim(nm),tksw(prgnam,nm),rfoca(j),
     .          def_r8=xxx,cindx=jj,note=
     .          'Smoothing for core tail.  input<0 => choose default * -input')
C     Cannot set default here: need set after rescaling of rmt
C     xxx = 0.50d0*rmt(j); if (io_help .eq. 1) xxx = NULLI
           xxx = 0d0; if (io_help .eq. 1) xxx = NULLI
           nm='SPEC_ATOM_RSMFA'; call gtv(trim(nm),tksw(prgnam,nm),rsmfa(j),
     .          def_r8=xxx,cindx=jj,note=
     .          'Smoothing for free atom.  input<0 => choose default * -input')

           nm='SPEC_ATOM_RCFA'; call gtv(trim(nm),tksw(prgnam,nm),rcfa(1:2,j)
     .          ,def_r8v=zerov,nmin=2,cindx=jj,note=
     .          'Cutoff radius for renormalization of free atom density'//
     .          '(WARN:takao rnatm.F is not tested).'//
     .          '%N%3fOptional 2nd argument = width'//
     .          '%N%3fRCFA<0 => renormalize potential instead of density')

C     Negative radii: convert to actual numbers
           if (rg(j)    .lt. 0) rg(j)    = -rg(j)*0.25d0*rmt(j)
           if (rsma(j)  .lt. 0) rsma(j)  = -rsma(j)*0.4d0*rmt(j)
           if (rfoca(j) .lt. 0) rfoca(j) = -rfoca(j)*0.4d0*rmt(j)
           if (rsmfa(j) .lt. 0) rsmfa(j) = -rsmfa(j)*0.5d0*rmt(j)

           nm='SPEC_ATOM_RS3'; call gtv(trim(nm),tksw(prgnam,nm),rs3(j),
c     .    def_r8=1d0,cindx=jj,
     .          def_r8=0.5d0,cindx=jj, !takao Nov26,2010
     .          note='Minimum smoothing radius for local orbital')

C     Downfolding switches: if auto DNF turned on, default is zero.
c     if (lham4) then
c     nn = 0
c     C       Help mode: If ADNF COULD have been turned on, default unknown
c     elseif (tksw(prgnam,'OPTIONS_ADNF')/= 2.and. io_help==1) then
c     nn = NULLI
c     C       All other cases: default is 1
c     else
           nn = 1
c     endif
           nm='SPEC_ATOM_IDXDN'; call gtv(trim(nm),tksw(prgnam,nm),
     .          idxdn(1:nlbj,j),nmin=nlbj,
     .          def_i4v=(/(nn,i=1,nlbj)/),cindx=jj,
     .          note='downfolding index: 0, auto; 1, no dnf; 2, fold down;'//
     .          ' 3, neglect')

C     'Hard Core' radius entering into strux
           nm='SPEC_ATOM_HCR'; call gtv(trim(nm),tksw(prgnam,nm),
     .          hcr(1:nlbj,j),nmin=nlbj,cindx=jj,nout=nout,or=T,
     .          note='Hard sphere radii for structure constants')
           if (nout == 0) then  !nout=-1 if sw=2; otherwise nout=0 unless data was read
              nm='SPEC_ATOM_HCR/R'; call gtv(trim(nm),tksw(prgnam,nm),
     .             hcr(1:nlbj,j),nmin=nlbj,cindx=jj,nout=nout,
     .             def_r8v=(/(0.7d0,I=1,n0)/),note=
     .             'Hard sphere radii for structure constants, units of R')
              call dscal(nout,rmt(j),hcr(1,j),1)
           endif
           if(debug) print *,'nspec eee mxcst j-loop j nspec',j,nspec

          nm='SPEC_ATOM_IDMOD'; call gtv(trim(nm),tksw(prgnam,nm),
     .         idmod(1:nlaj,j),nmin=nlaji,def_i4v=(/(0,i=1,n0)/),
     .         cindx=jj,note=
     .         'idmod=0 floats P to band CG, 1 freezes P, 2 freezes enu')

          nm='SPEC_ATOM_DV'; call gtv(trim(nm),tksw(prgnam,nm),dv(j),
     .         def_r8=0d0,cindx=jj,note='Artificial constant potential '//
     .         'shift added to spheres belonging to this species')
          nm='SPEC_ATOM_MIX'; call gtv(trim(nm),tksw(prgnam,nm),mxcst1(j),
     .         def_lg=F,cindx=jj,note='Set to suppress '//
     .         'self-consistency of classes in this spec')
          if(debug) print *,' mxcst nspec 111 =',nspec,j,jj,mxcst1(j)
          nm='SPEC_ATOM_CSTRMX'; call gtv(trim(nm),tksw(prgnam,nm),
     .         mxcst2(j),cindx=jj,def_lg=F,note='Set to exclude this'//
     .         ' species when automatically resizing sphere radii'//
     .         ' (SCLWSR>0)')
          if (sclwsr .eq. 0) mxcst2(j) = F
          nm='SPEC_ATOM_FRZWF'; call gtv(trim(nm),tksw(prgnam,nm),
     .         mxcst4(j),cindx=jj,def_lg=F,note=
     .         'Set to freeze augmentation wave functions for this species')
        endif    ! end of input dependent on presence of aug sphere.
        if(debug) print *,'nspec eee222 mxcst j-loop j nspec',j,nspec
!!
        sw = tksw(prgnam,'SPEC_ATOM_IDU')
        if (io_help .gt. 0 .and. sw .lt. 2) then
           call info0(1,0,0,' * ... The next three tokens are for LDA+U')
        endif
        nm='SPEC_ATOM_IDU'; call gtv(trim(nm),sw,idu(:,j),cindx=jj,
     .       def_i4v=(/(0,i=1,n0)/),note=
     .       'LDA+U mode:  0 nothing, 1 AMF, 2 FLL, 3 mixed; +10: nothin if sigm.* exist')
        nm='SPEC_ATOM_UH'; call gtv(trim(nm),sw,uh(:,j),cindx=jj,
     .       def_r8v=zerov,note='Hubbard U for LDA+U')
        nm='SPEC_ATOM_JH'; call gtv(trim(nm),sw,jh(:,j),cindx=jj,
     .       def_r8v=zerov,note='Exchange parameter J for LDA+U')

!     ! 2019 auto off lda+u uh jh if sigm exist
        if(procid==master .and. sum(abs(idu(:,j)))/=0) then 
!     ! get sexist if sigm.* exist.
c     ifi=fopn('tempyyy') 
c     sss=''
c     INQUIRE (ifi, NAME=sss) !second run, 'sss' must contain extention, eg, '.copt'
c     call fclose(ifi)
c     ipos = scan(sss,'.',back=.true.) !backward search of '.' from the last.
c     sname= sss(ipos+1:256) !sname contains extention such as 'copt'.
c     open(newunit=ifi, file='tempyyy.'//trim(sname)) 
c     close(ifi,status='delete')
           inquire(file='sigm.'//trim(sname),exist=sexist)
           if(sexist) then
              write(6,*)'We find sigm file --> noLDA+U-> set UH=JH=0'
              do lxx=0+1,3+1
                 if(idu(lxx,j)>10) then
                    uh(lxx,j) = 0d0
                    jh(lxx,j) = 0d0
                 endif
              enddo  
           endif
           do lxx=0+1,3+1
              idu(lxx,j) = mod(idu(lxx,j),10)
           enddo  
        endif
        call mpibc1(idu(:,j),4,2,mlog,'readctrl','idu')
        call mpibc1(uh(:,j),4,4,mlog,'readctrl','uh')
        call mpibc1(jh(:,j),4,4,mlog,'readctrl','jh')
c     write(6,"('j idu=',i5,256i5)")j, idu(:,j)
c     write(6,"('j uh =',i5,256f10.5)")j, uh(:,j)
c     write(6,"('j jh =',i5,256f10.5)")j, jh(:,j)

        nm='SPEC_ATOM_PBAS'; call gtv(trim(nm),tksw(prgnam,nm),pb1(j),
     .       nmin=10,cindx=jj,note='product basis for GW')
        nm='SPEC_ATOM_PBAS2'; call gtv(trim(nm),tksw(prgnam,nm),
     .       pb2(j),cindx=jj,note='second product basis for GW')
        nm='SPEC_ATOM_LMXPB';call gtv(trim(nm),tksw(prgnam,nm),lmxpb(j),
     .       def_i4=4,cindx=jj,note='l-cutoff for product basis')
c     if (nout==0) spec_lmxpb(j)='        '

C     Some sanity checks
        if (io_help == 0 .and. lmxaj .ge. 0) then
           if (tksw(prgnam,'SPEC_ATOM_LFOCA') /= 2)
     .          l_dummy_isanrg=isanrg(lfoca(j),0,2,'rdctrl','lfoca',T)
           if (tksw(prgnam,'SPEC_ATOM_LMXL') /= 2) l_dummy_isanrg=
     .          isanrg(lmxl(j),min(0,lmxaj),max(0,lmxaj),'rdctrl','lmxl',T)
           if (tksw(prgnam,'SPEC_ATOM_KMXA') /= 2)
     .          l_dummy_isanrg=isanrg(kmxt(j),2,25,' rdctrl (warning):','kmxa',F)
        endif
        coreh(j) = ' '
        nm='SPEC_ATOM_C-HOLE'; call gtv(trim(nm),tksw(prgnam,nm),coreh(j),
     .       nmin=10,cindx=jj,note='Channel for core hole')
        nm='SPEC_ATOM_C-HQ'; call gtv(trim(nm),tksw(prgnam,nm),coreq(:,j),
     .       def_r8v=(/-1d0,0d0/),cindx=jj,nmin=2,note=
     .       'Charge in core hole.  '//
     .       'Optional 2nd entry is moment of core hole:'//
     .       '%N%5fQ(spin1) = full + C-HQ(1)/2 + C-HQ(2)/2'//
     .       '%N%5fQ(spin2) = full + C-HQ(1)/2 - C-HQ(2)/2')
        nm='SPEC_ATOM_EREF'; call gtv(trim(nm),tksw(prgnam,nm),eref(j),
     .       def_r8=0d0,cindx=jj,note=
     .       'Reference energy subtracted from total energy')
        nm='SPEC_ATOM_AMASS'; call gtv(trim(nm),tksw(prgnam,nm),mass(j),
     .       cindx=jj,note='Nuclear mass in a.u. (for dynamics)')
c     nm='SPEC_ATOM_COLOUR'; call gtv(trim(nm),tksw(prgnam,nm),
c     .  colxbs(:,j),def_r8v=zerov,cindx=jj,note=
c     .  'Colour for xbs') !,-3,4,i1,TF(37))
c     xxx = rmt(j); if (io_help .eq. 1) xxx = NULLI
c     nm='SPEC_ATOM_RADIUS'; call gtv(trim(nm),tksw(prgnam,nm),
c     .  radxbs(j),def_r8=xxx,cindx=jj,note='Radius for xbs')
      enddo                     ! Loop over species
C     ... Cleanup after looping over species data
      if (io_help==0) then
C     Maximum L-cutoff
         call imxmn(nspec,lmxa,1,i0,lmxax)
         nlmax = (lmxax+1)**2
         nkaph = nkaph+lpzi-1
         globalvariables%nkaph = dble(nkaph); globalvariables%l_nkaph = globalvariables%l_nkaph +1
         globalvariables%nlmax = dble(nlmax); globalvariables%l_nlmax = globalvariables%l_nlmax +1
         globalvariables%mxorb = dble(nkaph)*nlmax; globalvariables%l_mxorb = globalvariables%l_mxorb +1
      endif
 79   continue

!     ! --- Site ---
      if (io_show+io_help/=0 .and. tksw(prgnam,'SITE')/=2) call info0(1,1,0,' --- Parameters for site data ---')
      if (io_help .eq. 1) then
         nbas = 1
         nsite =2
      endif
      allocate(pos(3,nsite),vel(3,nsite),eula(3,nsite),vshft(nsite),
     .     ips(nsite),ipl(nsite),plv(nsite),irlx(3,nsite),
     .     iantiferro(nsite))   !mpole(nsite),dpole(3,nsite),
c     if (ltbe) allocate(delta(n0,nsite),ndelta(nsite))

C     Default values
      vel  = 0d0
      eula = 0d0
      vshft = 0d0
      ipl  = 0
      plv  = 0
      irlx = 0
C     Should always be set, if used
      ips  = NULLI
      pos  = NULLR
c     mpole = 0d0
c     dpole = 0d0

      if (io_help .eq. 1) then
         if (iprint() .gt. 0) write(*,383)
 383     format(/' * ',
     .        'The following tokens are input for each site. ',
     .        'Data sandwiched'/3x,'between successive occurences of ',
     .        'token ATOM apply to one site.'/3x,'Alternatively, all ',
     .        'site data can be read in via the SITE file.')
         nbas = 1
         nsite = 2
      endif

!     ! ... Site data, one pass for each atom
      do  j = 1, nbas
         if (io_help /= 0) then
            write(stdo,'(1x)')
         elseif (io_help == 0 .and. io_show>0) then
            call info(1,0,0,' ... Site %i',j,0)
         endif
         jj=(/1,j/)
         nm='SITE_ATOM'; call gtv(trim(nm),tksw(prgnam,nm),alabl,nmin=10,
     .        cindx=jj,note='Species label')
         if (io_help .ne. 1) then
            do  i = 1, nspec
c     print *, alabl,' xxx ', slabl(i)
               if (trim(alabl) == trim(slabl(i)) ) then
                  ips(j) = i
c     print *, alabl,' xxx ', i,ips(j)
                  goto 881
               endif
            enddo
            call rxs('Category SITE referred to'//
     .           ' nonexistent species: ',alabl)
         endif
 881     continue
!     !  ... Site positions
         sw = tksw(prgnam,'SITE_ATOM_XPOS')
         nm='SITE_ATOM_POS'; call gtv(trim(nm),tksw(prgnam,nm),pos(:,j),
     .        nout=nout,cindx=jj,note='Atom coordinates, in units of alat',
     .        or=(sw.ne.2))
         if (nout == 0 .or. tksw(prgnam,'SITE_ATOM_POS') == 2) then !nout=-1 if sw=2; otherwise nout=0 unless data was read
            nm='SITE_ATOM_XPOS'; call gtv(trim(nm),tksw(prgnam,nm),pos(:,j),
     .           cindx=jj,note='Atom coordinates, as (fractional) '//
     .           'multiples of the lattice vectors')
            call dcopy(3,pos(1,j),1,xv,1)
            call dmpy(plat,3,1,xv,3,1,pos(1,j),3,1,3,1,3)
         endif
         nm='SITE_ATOM_DPOS'; call gtv(trim(nm),tksw(prgnam,nm),
     .        xv(1:3),def_r8v=zerov,nout=nout,cindx=jj,note=
     .        'Shift in atom coordinates added to pos')
         if (nout .eq. 3) call daxpy(3,1d0,xv,1,pos(1,j),1)

         nm='SITE_ATOM_V0'; call gtv(trim(nm),tksw(prgnam,nm),vel(:,j),
     .        def_r8v=zerov,cindx=jj,
     .        note='Initial velocity for molecular dynamics')

         nm='SITE_ATOM_RELAX'; call gtv(trim(nm),tksw(prgnam,nm),irlx(:,j),
     .        def_i4v=(/(1,i=1,n0)/),cindx=jj,note=
     .        'relax site positions (lattice dynamics) or '//
     .        'Euler angles (spin dynamics)')

         nm='SITE_ATOM_VSHFT'; call gtv(trim(nm),tksw(prgnam,nm),vshft(j),
     .        def_r8=0d0,cindx=jj,
     .        note='Constant potential shift for this site')

         nm='SITE_ATOM_ROT'; call gtv(trim(nm),tksw(prgnam,nm),outs,nmin=
     .        10,cindx=jj,nout=nout,note=
     .        'Rotation of spin quantization axis at this site')
         if (io_help .ne. 1 .and. nout .eq. 1) then
            call numsyv(i0)
            call lodsyv('ib',0,dble(j),ii)
            call lodsyv('x', 0,pos(1,j),ii)
            call lodsyv('y', 0,pos(2,j),ii)
            call lodsyv('z', 0,pos(3,j),ii)
            call a2rotm(outs,F,iprint()-10,xv)
            call clrsyv(i0)
            call rm2eua(xv,eula(1,j),eula(2,j),eula(3,j))
            if (io_show>0) call info5(0,0,0,
     .           '%15pROT= : Euler alpha = %1;6d  beta = %1;6d'//
     .           '  gamma = %1;6d',eula(1,j),eula(2,j),eula(3,j),0,0)
         endif

         nm='SITE_ATOM_PL'; call gtv(trim(nm),tksw(prgnam,nm),ipl(j),
     .        def_i4=0,cindx=jj,
     .        note='Assign principal layer number to this site')

         nm='SITE_ATOM_PLV'; call gtv(trim(nm),tksw(prgnam,nm),plv(j),
     .        def_i4=0,cindx=jj,
     .        note='Assign PL potential index to this site')
!     ! takao
         nm='SITE_ATOM_AF'; call gtv(trim(nm),tksw(prgnam,nm),iantiferro(j),
     .        cindx=jj,def_i4=0,
     .        note='antiferro ID:=i and -i should be af-pair, we look for space-group operation with spin-flip')
      enddo
 89   continue

!! --- Structure constants ---
      if (io_show+io_help/=0 .and. tksw(prgnam,'STR')/=2) call info0(1,1,
     .     0,' --- Parameters for structure constants ---')
      nm='STR_RMAXS'; call gtv(trim(nm),tksw(prgnam,nm),str_rmax,
     .     nout=nout,note='Radial cutoff for strux, in a.u.',or=T)
      if (nout == 0) then       !nout=-1 if sw=2; otherwise nout=0 unless data was read
         nm='STR_RMAX'; call gtv(trim(nm),tksw(prgnam,nm),str_rmax,
     .        def_r8=0d0,note='Radial cutoff for strux, in units of avw')
         str_rmax = str_rmax*avw
      endif
      nm='STR_MXNBR'; call gtv(trim(nm),tksw(prgnam,nm),str_mxnbr,
     .     def_i4=0,note='Max number of nbrs (for dimensioning arrays)')

!     ! --- Brillouin Zone ---
      if (io_show+io_help/=0 .and. tksw(prgnam,'BZ')/=2) call info0(1,1,
     .     0,' --- Parameters for Brillouin zone integration ---')
      nm='BZ_GETQP'; call gtv(trim(nm),tksw(prgnam,nm),bz_lio1,
     .     def_lg=F,note='Read qp from disk',or=T)
      nm='BZ_NKABC'; sw=tksw(prgnam,nm); if (bz_lio1) sw = 2
      call gtv(trim(nm),sw,bz_nabc,nout=nout,
     .     note='No. qp along each of 3 lattice vectors.'//
     .     '%N%3fSupply one number for all vectors or a separate '//
     .     'number for each vector.')
      call fill3in(nout,bz_nabc)
      nm='BZ_PUTQP'; call gtv(trim(nm),tksw(prgnam,nm),bz_lio2,
     .     def_lg=F,note='Write qp to disk')
      nm='BZ_BZJOB';call gtv(trim(nm),tksw(prgnam,nm),bz_lshft,nout=nout,
     .     def_i4v=izerv(1:1),note=
     .     '0 centers BZ mesh at origin, 1 centers off origin'//
     .     '%N%3fSupply one number for all vectors or a separate '//
     .     'number for each vector.')
      call fill3in(nout,bz_lshft)
      nm='BZ_METAL'; call gtv(trim(nm),tksw(prgnam,nm),bz_lmet,
     .     def_i4=1,note=
     .     '0 assume insulator; 1 save evecs on disk; '//
     .     '2 use wgt from prior iter%N%3f3 always make two band passes;'//
     .     ' 4 BZ integration with 3-point scheme')
      nm='BZ_TETRA'; call gtv(trim(nm),tksw(prgnam,nm),ctrl_lmet2,
     .     def_lg=T,note='Tetrahedron integration')
      if( cmdopt('--tdos',6,0,a).or.cmdopt('--pdos',6,0,a)) then
         print *,' --tdos or --pdos enforces BZ_METAL=3 and BZ_TETRA=1'
         bz_lmet=3
         ctrl_lmet2=.true.
      endif
      nm='BZ_N'; call gtv(trim(nm),tksw(prgnam,nm),bz_n,
     .     def_i4=0,note=
     .     'N>0: Polynomial order for Methfessel-Paxton sampling%N%3f'//
     .     'N=0: Conventional Gaussian sampling%N%3f'//
     .     'N<0: Broadening by Fermi-Dirac distribution%N%3f'//
     .     'To be used in conjunction with W= ; see next')
      nm='BZ_W'; call gtv(trim(nm),tksw(prgnam,nm),bz_w,
     .     def_r8=5d-3,note=
     .     'N>=0: Line broadening for sampling integration%N%3f'//
     .     'N<0 : Temperature for Fermi distribution (Ry)')
      nm='BZ_EF0'; call gtv(trim(nm),tksw(prgnam,nm),bz_ef,
     .     def_r8=0d0,note='Initial guess at Fermi energy')
      nm='BZ_DELEF'; call gtv(trim(nm),tksw(prgnam,nm),bz_def,
     .     def_r8=0.05d0,note='Initial uncertainty in Fermi energy')
      nm='BZ_ZBAK'; call gtv(trim(nm),tksw(prgnam,nm),zbak(1),
     .     def_r8=0d0,note='Homogeneous background charge')
      nm='BZ_SAVDOS'; call gtv(trim(nm),tksw(prgnam,nm),ctrl_ldos,
     .     def_i4=0,note='Choose some combination of the following:'
     .     //'%N%3f1 Write DOS to directly disk (NPTS and DOS also needed)'
     .     //'%N%3f2 Write weights for partial DOS'
     .     //'%N%3f4 Same as (2), but weights m-resolved')
!     ! -- Sanity checks and other initialization ---
      if( iand(8,ctrl_ldos) /=0 .and. iand(4+2,ctrl_ldos) /=0 ) 
     &     call rx ( 'inconsistent BZ_SAVDOS options' )
      nm='BZ_DOS'; call gtv(trim(nm),tksw(prgnam,nm),bz_dosw,def_r8v=
     .     (/-1d0,0d0/),note='Energy window over which DOS accumulated')
      nm='BZ_NPTS'; call gtv(trim(nm),tksw(prgnam,nm),bz_ndos,def_i4=1001,
     .     note='No. DOS points (sampling integration, and lmdos)')
!     ! mar2013
      nm='BZ_DOSMAX'; call gtv(trim(nm),tksw(prgnam,nm),bz_dosmax,def_r8=1.5d0,
     &     note='Maximum energy to which DOS accumulated, relative to Efermi')
      xxx = 5d0
      nm='BZ_EFMAX'; call gtv(trim(nm),tksw(prgnam,nm),bz_efmax,
     .     def_r8=xxx,note='Find evecs up to efmax')
      nm='BZ_NEVMX'; call gtv(trim(nm),tksw(prgnam,nm),bz_nevmx,
     .     def_i4=0,note='Find at most nevmx eigenvectors'//
     .     '%N%3fIf NEVMX=0, program uses internal default'//
     .     '%N%3fIf NEVMX<0, no eigenvectors are generated')
      nm='BZ_NOINV'; call gtv(trim(nm),tksw(prgnam,nm),ctrl_noinv,
     .     def_lg=F,note=
     .     'Suppress automatic inclusion of inversion symmetry for BZ')
      nm='BZ_FSMOM'; call gtv(trim(nm),tksw(prgnam,nm),bz_fsmom,
     .     def_r8=NULLR,note='Fixed-spin moment (fixed-spin moment method)')
      nm='BZ_FSMOMMETHOD';call gtv(trim(nm),tksw(prgnam,nm),bz_fsmommethod,
     .     def_i4=0,note='Method of Fixed-spin moment 0:original 1:discrete')

!! --- Ewald sums ---
      if (io_show+io_help/=0 .and. tksw(prgnam,'EWALD')/=2)
     .     call info0(1,1,0,' --- Parameters for Ewald sums ---')
      nm='EWALD_AS'; call gtv(trim(nm),tksw(prgnam,nm),lat_as,
     .     def_r8=2d0,note='Ewald smoothing parameter')
      nm='EWALD_TOL'; call gtv(trim(nm),tksw(prgnam,nm),lat_tol,
     .     def_r8=1d-8,note='Ewald tolerance')
      nm='EWALD_NKDMX'; call gtv(trim(nm),tksw(prgnam,nm),lat_nkdmx,
     .     def_i4=800,note='Ewald tolerance')
      nm='EWALD_RPAD'; call gtv(trim(nm),tksw(prgnam,nm),lat_rpad,
     .     def_r8=0d0,note='Scale rcutoff by rpad when lattice vectors'//
     .     ' padded in oblong geometries')
      
!! --- Iterations (formerly MIX) ---
      vmix%b = NULLI            ! Not set
      if (tksw(prgnam,'ITER')/=2) then
         if (io_show+io_help/=0)
     .        call info0(1,1,0,' --- Parameters for iterations ---')
C     Default values for smix (array has same same structure as lstra smix)
         call umix_init(vmix)
         vmix%b = 1             ! beta
         vmix%bv = 1            ! bv
         vmix%fn='mixm'
         vmix%kill =  0         ! nkill
         vmix%lxpot =  0        ! lxpot
         vmix%mmix = -1         ! mmix
         vmix%mode = 0          ! mode (0=Anderson)
         vmix%nsave = 8         ! nsave = # iter to save on disk
         vmix%r=' '
         vmix%tolu = 0          ! tolu
         vmix%umix = 1          ! umix (mixing parm for LDA+U)
         vmix%w(1) = 1          ! w(1)
         vmix%w(2) = 1          ! w(2)
         vmix%wc = -1           ! wc
         smalit = NULLI
         nm='ITER_NIT';call gtv(trim(nm),tksw(prgnam,nm),iter_maxit,
     .        def_i4=30,
     .        note='maximum number of iterations in self-consistency cycle')
         nm='ITER_NRMIX'; call gtv(trim(nm),tksw(prgnam,nm),smalit,
     .        def_i4v=(/80,2/),note='Sphere program, (1) max iter; '//
     .        '(2) no. prior iter for Anderson mixing ')
         nm='ITER_MIX'; sw=tksw(prgnam,nm); call gtv(trim(nm),sw,iter_mix,
     .        nmin=10,nout=nout,note='Mixing rules for charge mixing.  Syntax:')
         if (io_help/=0 .and. tksw(prgnam,nm)/=2) print 345
 345     format(3x,
     .        'A[nmix][,b=beta][,bv=betv][,n=nit][,w=w1,w2][,nam=fn][,k=nkill]',
     .        '[;...] or'/
     .        3x,'B[nmix][,b=beta][,bv=betv][,wc=wc][,n=#][,w=w1,w2][,nam=fn]',
     .        '[,k=nkill]')
c     c takao
c     mmm=adjustl(iter_mix)
c     if(mmm(1:1)=='B')
c     &      call rx(" Broyden mixing is problematic."//
c     &              " If necessary, you have to fix it (mixrho.F)")
         if (nout==1 .and. io_help==0 .and. sw/=2) then
            alabl=vmix%fn
            if (io_show==0) call pshpr(0)
            if (.not. parmxp(-1,iter_mix,len_trim(iter_mix),
     .           int(vmix%mode),int(vmix%nmix),vmix%w(1),vmix%b,vmix%elind,
     .           alabl,vmix%wc,vmix%kill,vmix%bv,xxx)) then
               if (io_show==0) call poppr
               call rx( 'RDCCAT: parse in parmxp failed 1')
            endif
            if (io_show==0) call poppr
         endif

         nm='ITER_CONV'; call gtv(trim(nm),tksw(prgnam,nm),ctrl_tol(3),
     .        def_r8=1d-4,note='Tolerance in energy change from prior'//
     .        ' iteration for self-consistency')
         nm='ITER_CONVC'; call gtv(trim(nm),tksw(prgnam,nm),ctrl_tol(1),
     .        def_r8=1d-4,note='Tolerance in output-input charge for'//
     .        ' self-consistency')
!     !
         i0 = 0
         nm='ITER_XIPMX'; call gtv(trim(nm),tksw(prgnam,nm),i0,def_i4=0,
     .        nout=nout,note='Mix potential independently of charge:'//
     .        '%N%3fXIPMX=1: mix vin and v(qmix)'//
     .        '%N%3fXIPMX=2: mix vin and v(qout)')
         if (nout .gt. 0) then
            l_dummy_isanrg = isanrg(i0,0,3,'rdctrl','XIPMX',T)
            if (i0 .gt. 0) lves = 1 
         endif
         vmix%lxpot = i0

         nm='ITER_UMIX'; call gtv(trim(nm),tksw(prgnam,nm),vmix%umix,
     .        def_r8=1d0,note='Mixing parameter for densmat in LDA+U')
         nm='ITER_TOLU'; call gtv(trim(nm),tksw(prgnam,nm),vmix%tolu,
     .        def_r8=0d0,note='Tolerance for densmat in LDA+U')
         nm='ITER_NITU'; call gtv(trim(nm),tksw(prgnam,nm),i0,
     .        def_i4=0,note='Max number of LDA+U iterations of densmat')
         vmix%nitu = i0
      endif                     ! iterations category
!! --- Dynamics ---
      if (io_show+io_help/=0 .and. tksw(prgnam,'DYN')/=2) call info0(1,1,
     .     0,' --- Parameters for dynamics and statics ---')
      nm='DYN_NIT'; call gtv(trim(nm),tksw(prgnam,nm),nitmv,def_i4=1,
     .     note='maximum number of relaxation steps (statics)'//
     .     ' or time steps (dynamics)')
      nm='DYN_MSTAT'; call gtv(trim(nm),tksw(prgnam,nm),nono,Texist=ltmp,
     .     note='Parameters for molecular statics')
      mdprm(1) = 0
      if (io_help .ne. 0 .or. ltmp) then
         if(ctrl_lfrce==0) ctrl_lfrce=1 ! takao july2012
         nm='DYN_MSTAT_MODE'; call gtv(trim(nm),tksw(prgnam,nm),i0,
     .        def_i4=0,note=
     .        '0: no relaxation  '//
     .        '4: conjugate gradients  '//
     .        '5: Fletcher-Powell  '//
     .        '6: Broyden')
         mdprm(1) = i0
         nm='DYN_MSTAT_HESS'; call gtv(trim(nm),tksw(prgnam,nm),ltmp,
     .        def_lg=T,note='Read hessian matrix')
         mdprm(2) = isw(ltmp)
         nm='DYN_MSTAT_XTOL'; call gtv(trim(nm),tksw(prgnam,nm),mdprm(3),
     .        def_r8=1d-3,note=
     .        'Convergence criterion in displacements'//
     .        '%N%3fXTOL>0: use length;  <0: use max val;  =0: do not use')
         nm='DYN_MSTAT_GTOL'; call gtv(trim(nm),tksw(prgnam,nm),mdprm(4),
     .        def_r8=0d0,note=
     .        'Convergence criterion in gradients'//
     .        '%N%3fGTOL>0: use length;  <0: use max val;  =0: do not use')
         nm='DYN_MSTAT_STEP'; call gtv(trim(nm),tksw(prgnam,nm),mdprm(5),
     .        def_r8=0.015d0,note='Initial (and maximum) step length')
         nm='DYN_MSTAT_NKILL'; call gtv(trim(nm),tksw(prgnam,nm),i0,
     .        def_i4=0,note='Remove hessian after NKILL iter')
         mdprm(6) = i0
         nm='DYN_MSTAT_PDEF'; call gtv(trim(nm),tksw(prgnam,nm),lat_defm,
     .        def_r8v=zerov,note='Lattice deformation modes')
      endif
      if (io_help > 0) call cexit(0,1)
cccc  end of read input parameter ccccccccccccccccccccccccccxxxxxxxxxxxxxxxx
      
c!! Optionally read positions from pos file
c      if (cmdopt('--rpos=',7,0,fileid)) then
c        call iopos(.false.,-1,fileid(8:),nbasp,pos)
c      endif

c$$$!! --- Allocate and copy input to sbz ---
c$$$      call ubz_init(v_sbz)
c$$$      v_sbz%size=ubz_size()
c$$$      v_sbz%def=bz_def
c$$$      v_sbz%dosw=bz_dosw(:)
c$$$      v_sbz%dosmax=bz_dosmax
c$$$      v_sbz%ef=bz_ef
c$$$      v_sbz%efmax=bz_efmax
c$$$      v_sbz%fsmom=bz_fsmom
c$$$      v_sbz%fsmommethod=bz_fsmommethod
c$$$      v_sbz%lcond=bz_lcond(:)
c$$$      v_sbz%lmet=bz_lmet
c$$$      v_sbz%lmull=bz_lmull
c$$$      v_sbz%lshft = bz_lshft
c$$$      v_sbz%n     = bz_n
c$$$      v_sbz%ndos  = bz_ndos
c$$$      v_sbz%nevmx = bz_nevmx
c$$$      v_sbz%nkabc = bz_nabc
c$$$      v_sbz%range=bz_range
c$$$      v_sbz%w=bz_w
!! Allocate and copy input to sctrl ---
      v_sctrl%size=uctrl_size()
!! lbas: 1 Hamiltonian has no screening transformation
!!           2 Hamiltonian is nfp style
!!          16 freeze phi,phidot for all species
!! lbas = 3 = 1+2 in lmf
!! frzwf: freeze basis or not. j= 1+2+16 if frzwf=T
      lbas = 3 !!   For now, LMF => fp; no screening; nfp-style hamiltonian
      j = lbas + 16*isw(frzwf)
      v_sctrl%lbas=j
      v_sctrl%defm=lat_defm(:)
      v_sctrl%ldos=ctrl_ldos
      v_sctrl%lfp=lfp
      v_sctrl%lfrce=ctrl_lfrce
!! takao. we may need to remove lmet4 or lmet8... but a little complicated...
!! Now ctrl_lmet4 and ctrl_lmet8 are fixed as .true.in m_rdctl.F
!!     ctrl_lmet    1 metal  2 tetrahedron
!!                 4 (PGF) V-shift1 is zero
!!                 8 (PGF) V-shift2 is zero
      lmet = isw(bz_lmet.ne.0) + 2*isw(ctrl_lmet2) != 1*lmetal + 2*ltetrahedron
     .+4*isw(ctrl_lmet4) + 8*isw(ctrl_lmet8)

C takao in lmf only lncol4, lncol32, lncol64 are used.
C     lncol 1 noncollinear magnetism                      !=F in lmf
C           2 spin spirals                                !=F in lmf
C           4 spin-orbit coupling                         !Used in lmf.
C           8 External magnetic field                     !=F in lmf
C          16 mag. forces                                 !=F in lmf
C          32 spin-orbit coupling, LzSz only              !Used in lmf
C          64 spin-orbit coupling, LzSz + (L.S-LzSz) pert !Used in lmf
C     If spin-orbit or SS, also turn on noncollinear

c      k = 2*isw(lncol2)+4*isw(lncol4)+8*isw(lncol8)+16*isw(lncol16)
c      if (k .ne. 0) lncol1=T
c      lncol = 1*isw(lncol1)+2*isw(lncol2)+4*isw(lncol4)+8*isw(lncol8)+
c     .  16*isw(lncol16)+32*isw(lncol32)+64*isw(lncol64)

c takao lncol is used only in limited routines.
c One of lncol4,lncol32, or lncor64 is true, others are false.
c      lncol= 4*isw(lncol4)+32*isw(lncol32)+64*isw(lncol64)
      v_sctrl%lmet=lmet
c      v_sctrl%lncol=lncol
C     lrs  switches concerning restart mode.
      irs=0
      irs(1) = 1
      irs(2) = 1
      if (cmdopt('--rs=',5,0,strn)) then
        i = 5
        j = a2vec(strn,len(strn),i,2,', ',2,2,5,ix,irs)
      endif
      irs(1) = mod(mod(irs(1),10),4) + 4*getdig(irs(1),1,10)   + 8*getdig(irs(1),1,100)
      lrs = 1*mod(irs(1),8) + 8*irs(2)+32*irs(3) + 64*irs(4) + 128*irs(5) + 256*mod(irs(1)/8,2)
!!
!!     irs(1) tells what to read and whether to invoke smshft.
!!     4s' bit of irs(1) -> invoke smshft after file read.
!!     8s' bit of irs(1) -> rotate local density after file read
!!     0+1's bits irs(1)     action
!!    0              read from atom file
!!    1              read from binary rst file
!!    2              read from ascii rsta file
!!    3              read nothing (data already input)

!!         1 Read from restart file
!!         2 Read from restart file, ascii mode
!!         4 Read from restart file, invoke smshft
!!       256 rotate local density after reading
      irs1 = 1*mod(irs(1),8) + 256*mod(irs(1)/8,2) !  iand(1+2+4,int(lrs)) + 8*isw( iand(256,int(lrs)) .ne.0 )
!!         8 Write new density to restart file
!!        16 Write new density to restart file, ascii format
      irs2 = irs(2) !iand( 8+16,int(lrs)) /8
!!        32 read site positions from input file
      irs3 = irs(3) ! isw ( iand(32,int(lrs)) .ne.0 )
!!        64 read starting fermi level from input file
      irs4 = irs(4) !isw ( iand(64,int(lrs)) .ne.0 )
!!       128 read starting pnu level from input file
      irs5 = irs(5) !isw ( iand(128,int(lrs)) .ne.0 )
      irs11 = 1*mod(irs(1),8) !iand(1+2+4,int(lrs))

C     lqp 1 do not add inversion 2 inverse iteration
c      lqp = 1*isw(ctrl_noinv)+2*isw(ctrl_lqp2)  !lqp= ctrl_noinv + 2*ctrl_lqp2
      noinv = isw (ctrl_noinv)  ! T->1 F->0
c      lpq2 = isw (ctrl_lqp2)  ! T->1 F->0
C     lscr 0 do nothing
C          1 Make P0(0)
C          2 Screen output q and ves
C          3 Screen output ves only
C            Add 10*k to compute intra-site contribution to
C            vbare each kth iteration
C            Add 100*k to compute response function only
C            each kth iteration
C          4 Use model response to screen output q
C            Add 1 to combine mode 1 with another mode
C            Add 10*k to compute intra-site contribution to
C            vbare each kth iteration
C            Add 100*k to compute response function only each kth iteration
      v_sctrl%noinv=noinv
      v_sctrl%lrel=lrel
      v_sctrl%lrs=lrs
      v_sctrl%lves=lves !!!  lves    1 take ves as input
!!  lxcf   parameter defining XC functional
!!      1s digit:
!!      1 for Ceperly-Alder (VWN)
!!      2 for Barth-Hedin (ASW fit)
!!      103 for PBE-GGA
      v_sctrl%lxcf= ham_lxcf
      v_sctrl%maxit=iter_maxit
      v_sctrl%mdprm=mdprm
      v_sctrl%modep=2
      v_sctrl%nbas=nbas

      nl = max(lmxbx,lmxax)+1

      v_sctrl%nbasp=nbasp
      v_sctrl%nesabc=nesabc
      v_sctrl%nitmv=nitmv
      v_sctrl%nl=nl
      v_sctrl%nsite=nsite
      v_sctrl%nspec=nspec

      v_sctrl%nspin=nsp
      v_sctrl%nvario=nvario
      v_sctrl%omax1 = omax1
      v_sctrl%omax2 = omax2
      v_sctrl%quit  = quit
      v_sctrl%rmaxes= rmaxes
      v_sctrl%rmines= rmines
      v_sctrl%sclwsr= sclwsr
      v_sctrl%sdmod = sdmod
      v_sctrl%sdprm = sdprm
      v_sctrl%smalit= smalit
      v_sctrl%tol   = ctrl_tol
      v_sctrl%wsrmax= wsrmax
      v_sctrl%zbak  = zbak
      v_sctrl%pfloat= lpfloat

      v_slat%size=ulat_size()
      if (dalat == NULLR) dalat=0
      v_slat%alat=alat+dalat
      v_slat%as=lat_as
      v_slat%avw=avw
      v_slat%nkdmx=lat_nkdmx
      v_slat%nkqmx=lat_nkdmx
      v_slat%gam =lat_gam
      v_slat%gmax=lat_gmax
      v_slat%nabc=ftmesh
      v_slat%ldist= lat_ldist
      v_slat%dist = lat_dist
      v_slat%plat=plat
      v_slat%rpad=lat_rpad
      v_slat%slat=lat_slat
      v_slat%tol=lat_tol
      v_slat%tolft=tolft
      v_slat%vol=vol
!! setcg
!!... Choose dimensions for arrays
      lmxcg=8
      lmxcy=12
      if (lmxcg .le. 6) then
        lnjcg = 6500
        lnxcg = 1300
      else if (lmxcg .le. 8) then
        lnjcg = 22700
        lnxcg = 3400
      else if (lmxcg .le. 10) then
        lnjcg = 62200
        lnxcg = 7400
      else
        call rxi('setcg: cannot handle lmxcg=',lmxcg)
      endif
C ... Allocate and occupy the arrays
      nlm=(lmxcy+1)**2
      allocate(rv_a_ocy(abs(nlm)))
      allocate(rv_a_ocg(abs(lnjcg)))
      allocate(iv_a_ojcg(abs(lnjcg)))
      allocate(iv_a_oidxcg(abs(lnxcg)))
      call sylmnc ( rv_a_ocy , lmxcy )
      call scg ( lmxcg , rv_a_ocg , iv_a_oidxcg , iv_a_ojcg )
!!

C --- Allocate and copy input to sham ---
      v_sham%size=uham_size()
      v_sham%alfsi=alfsi
      i_copy_size=size(v_sham%dabc)
      call dcopy(i_copy_size,dabc,1,v_sham%dabc,1)
      v_sham%elind=elind
c      v_sham%nmto=nmto
c      i_copy_size=size(v_sham%kmto)
c      call dcopy(i_copy_size,kmto,1,v_sham%kmto,1)

C     Replicate ctrl->lncol in ham->lncol, ditto for lham,lgen3
c      v_sham%lncol=(int(v_sctrl%lncol))

C     Mop up integer parts of sigp
c      sigp(1) = sigp_mode
c      sigp(2) = sigp_nmin
c      sigp(4) = sigp_nmax

      v_sham%nkaph=nkaph
      i_copy_size=size(v_sham%pmax)
      call dcopy(i_copy_size,pmax,1,v_sham%pmax,1)
      i_copy_size=size(v_sham%pmin)
      call dcopy(i_copy_size,pmin,1,v_sham%pmin,1)
      v_sham%lsig=lrsig
c      i_copy_size=size(v_sham%sigp)
c      call dcopy(i_copy_size,sigp,1,v_sham%sigp,1)
      i_copy_size=size(v_sham%qss)
      call dcopy(i_copy_size,ham_qss,1,v_sham%qss,1)

      v_sham%rsrnge=rsrnge
      v_sham%rsstol=rsstol
      v_sham%udiag=ham_udiag
      v_sham%scaledsigma=scaledsigma

!!  Parameters for APW
      v_sham%pwmode=pwmode
      v_sham%npwpad=npwpad
      v_sham%pwemin=pwemin
      v_sham%pwemax=pwemax
      v_sham%oveps=oveps
      v_sham%delta_stabilize=delta_stabilize !takao sep2010

!!... Allocate and initialize sarray
      v_sarry%size=uarray_size()

!! ... Allocate and initialize spot
      v_spot%size=upot_size()
      v_spot%vmtz0=vmtz

C --- Allocate and copy input to smix ---
      v_smix%size=umix_size()
      v_smix=vmix

cC --- Allocate and copy input to sstr ---
c      v_sstr%size=ustr_size()
c      v_sstr%mxnbr=str_mxnbr
c      v_sstr%rmax=str_rmax
C --- Allocate and copy input to sspec ---
      allocate(v_sspec(nspec))
        do i_spackv=1,nspec
          v_sspec(i_spackv)%size=uspec_size()
          v_sspec(i_spackv)%eh3=-0.5d0
          v_sspec(i_spackv)%etf=-1d0
          v_sspec(i_spackv)%lmxf= 2*nl-2
          v_sspec(i_spackv)%norp=2
          v_sspec(i_spackv)%vmtz=-0.5d0
        enddo

        do  j = 1, nspec
          slabl_(j) = slabl(j)
          v_sspec(j)%pb1=pb1(j)
          v_sspec(j)%pb2=pb2(j)
          v_sspec(j)%coreh=coreh(j)
          v_sspec(j)%name=slabl(j)
          v_sspec(j)%nmcore=nmcore(j)
          v_sspec(j)%a=spec_a(j)
          v_sspec(j)%nr=nr(j)
c          v_sspec(j)%alpha=alpha(:,j)
          v_sspec(j)%coreq=coreq(:,j)
          v_sspec(j)%lxi=lxi(j)
          v_sspec(j)%nxi=nxi(j)
          v_sspec(j)%exi=exi(:,j)
          v_sspec(j)%group=grp(j)
          v_sspec(j)%grp2=grp2(j)
          v_sspec(j)%hcr=hcr(:,j)
          v_sspec(j)%idmod=idmod(:,j)
          v_sspec(j)%ehvl=ehvl(:,j)
          ix = 1
          ix(1:1+lmxb(j)) = idxdn(1:1+lmxb(j),j)
          v_sspec(j)%idxdn = ix
          v_sspec(j)%idu = idu(:,j)
          v_sspec(j)%jh=jh(:,j)
          v_sspec(j)%uh=uh(:,j)
          v_sspec(j)%kmxt=kmxt(j)
          v_sspec(j)%kmxv=kmxv(j)
          v_sspec(j)%lfoca=lfoca(j)
          v_sspec(j)%rsmv=rsmv(j)
          v_sspec(j)%lmxa=lmxa(j)
          v_sspec(j)%lmxb=lmxb(j)
          v_sspec(j)%lmxl=lmxl(j)
          v_sspec(j)%lmxpb=lmxpb(j)

          i = 1*isw(mxcst1(j))+2*isw(mxcst2(j))+4*isw(mxcst4(j))
          orbp=0d0 !call dpzero(orbp,n0*2*nkap0)
          call dcopy(n0,rsmh(1,j),1,orbp(1,1,1),1)
          call dcopy(n0,eh(1,j),1,orbp(1,2,1),1)
          call dcopy(n0,rsmh2(1,j),1,orbp(1,1,2),1)
          call dcopy(n0,eh2(1,j),1,orbp(1,2,2),1)
          i_copy_size=size(v_sspec(j)%orbp)
          call dcopy(i_copy_size,orbp,1,v_sspec(j)%orbp,1)

          v_sspec(j)%mass=mass(j)
          v_sspec(j)%mxcst=i
!!     Pack P,Q,PZ for both spins
          v_sspec(j)%p = pnu(1:n0*nsp,1,j) 
          v_sspec(j)%pz= pz(1:n0*nsp,1,j)  
          v_sspec(j)%q = qnu(1:n0*nsp,1,j) 
          v_sspec(j)%z=z(j)
          v_sspec(j)%colxbs=colxbs(:,j)
          v_sspec(j)%radxbs=radxbs(j)
          v_sspec(j)%rcfa=rcfa(:,j)
          v_sspec(j)%rcut=rcut(j)
          v_sspec(j)%rfoca=rfoca(j)
          v_sspec(j)%rg=rg(j)
          v_sspec(j)%rmt=rmt(j)
          v_sspec(j)%rs3=rs3(j)
          v_sspec(j)%rsma=rsma(j)
          v_sspec(j)%rsmfa=rsmfa(j)
          v_sspec(j)%dv=dv(j)
          v_sspec(j)%eref=eref(j)
          v_sspec(j)%rham=rham(j)
          v_sspec(j)%rint=rint(j)
        enddo

        allocate(v_ssite(nsite))
        do i=1,nsite
          v_ssite(i)%size=usite_size()
        enddo
        do  j = 1, nsite
          v_ssite(j)%spec=ips(j)
          v_ssite(j)%class=ips(j)
          v_ssite(j)%clabel=slabl(ips(j))
          v_ssite(j)%pl=ipl(j)
          v_ssite(j)%plv=plv(j)
          v_ssite(j)%pos=pos(1:3,j)
          v_ssite(j)%vel=vel(1:3,j)
          v_ssite(j)%vshft=vshft(j)
          v_ssite(j)%relax=irlx(:,j)
          v_ssite(j)%eula=eula(:,j)
          v_ssite(j)%iantiferro=iantiferro(j)
        enddo
        sstrnmix=trim(iter_mix)

!! ... Suppress symmetry operations for special circumstances
!     !     Switches that automatically turn of all symops
!     ! --pdos mar2003 added. Also in lmv7.F
          lstsym = 0
          if ((mdprm(1) .ge. 1 .and. mdprm(1) .le. 3) .or.   
     .         cmdopt('--cls',5,0,strn) .or. cmdopt('--nosym',7,0,strn).or.cmdopt('--pdos',6,0,strn)) then
             symg = 'e'
             lstsym = 2         !lstsym=2: turn off symops
             noinv=1
          endif
          sstrnsymg=trim(symg)
!!
        if (cmdopt('--rdbasp',8,0,strn)) call rx('not support --rdbasp')

!!     Add dalat to alat
!!  v_slat%alat=(v_slat%alat)+dalat !this is a bug; this should be commented out. This
!!  is a bug for lm7K when I started lm7K. I had included this bug here. Fixed at 28May2010.

!! Dirac equation requires spin polarization
        if( nsp==1 .and. int(v_sctrl%lrel)==2 ) call rx('rdccat: Dirac equation requires NSPIN=2')
!! ... Suppress inversion when noncollinear magnetism, SX, NLO
        if(lso.ne. 0) noinv=1   !lqp = lqp-bitand(lqp,1)+1
        v_sctrl%noinv=noinv
!     ! ... Setup for idxdn ... ctrl->lham,4 is automatic downfolding switch
C     No screening => no downfolding; j=3
        j=3
C     nfp-style basis:
        j=j+10
        call suidx(globalvariables%nkaph,j,nspec,v_sspec)

!! ... Set some global variables
          globalvariables%nspec = dble(nspec); globalvariables%l_nspec = globalvariables%l_nspec +1
          globalvariables%nbas = dble(nbas);   globalvariables%l_nbas = globalvariables%l_nbas +1
          globalvariables%nbasp= dble(nbasp);  globalvariables%l_nbasp = globalvariables%l_nbasp +1
          globalvariables%nsp =  dble(nsp);    globalvariables%l_nsp = globalvariables%l_nsp +1
          globalvariables%nl =   dble(nl);     globalvariables%l_nl = globalvariables%l_nl +1
          globalvariables%lrel = dble(avw);    globalvariables%l_avw = globalvariables%l_avw +1
          globalvariables%lrel = dble(lrel);   globalvariables%l_lrel = globalvariables%l_lrel +1
          lxcf = int(v_sctrl%lxcf)
          globalvariables%lxcf = dble(lxcf); globalvariables%l_lxcf = globalvariables%l_lxcf +1
          globalvariables%stdo = dble(stdo); globalvariables%l_stdo = globalvariables%l_stdo +1
          globalvariables%stdl = dble(stdl); globalvariables%l_stdl = globalvariables%l_stdl +1
          globalvariables%stde = dble(stde); globalvariables%l_stde = globalvariables%l_stde +1
          nspc = 1
          if( lso==1 ) nspc = 2
c          if (iand(lncol,4) .ne. 0) nspc = 2
          globalvariables%nspc = dble(nspc); globalvariables%l_nspc = globalvariables%l_nspc +1; xx(1) = dble(nspc)
C     Make nat = number of real atoms as nbas - # sites w/ floating orbitals
          if (procid .eq. master) then
            nat = nbas
            do  i = 1, nbas
              j=v_ssite(i)%spec
              l=v_sspec(j)%lmxa
              if (l .eq. -1) nat = nat-1
            enddo
          endif
          call mpibc1(nat,1,2,0,'rdctrl','nat')
          globalvariables%nat = dble(nat); globalvariables%l_nat = globalvariables%l_nat +1; xx(1) = dble(nat)
C ... Set modep
          v_sctrl%modep=2
C ... Count LDA+U blocks (printout only)
          allocate(wowk(nbas))
          wowk=0
          call pshpr(0)
          call suldau(nbas,v_sspec,v_ssite,nlibu,k,wowk)
          v_sham%nlibu=nlibu
          call poppr
          deallocate(wowk,pnu,qnu,pz,amom,idmod,rsmh,eh,rsmh2,eh2,pb1,pb2,
     .     lmxpb,qpol,stni,tbvso,iq1,ivso,rg,rsma,rfoca,rsmfa,rcfa,nxi,
     .     exi,rint,rcut,coreq,mass,colxbs,radxbs,rs3,rham,idxdn,hcr,
     .     rmt,idu,uh,jh,dv,grp,grp2,mxcst1,mxcst2,mxcst4,kmxt,kmxv, !alpha,
     .     lfoca,eref,lmxl,lxi,coreh,lmxa,lmxb,spec_a,z,nr,rsmv,
     .     pos,vel,eula,vshft,ips,ipl,plv,irlx)

!! --- Printout ---
          if (iprint() .ge. 20 .and. procid .eq. master) then
            do  90  k = 1, 2
              strn = '  '//prgnam
              lxcf=v_sctrl%lxcf
c              lfrzw = isw ( iand(16,int(v_sctrl%lbas)) .ne.0 )
              lfrzw = isw(frzwf)
              call awrit8(' special:%10p'//
     .         '%?;n; forces,;;'//
     .         '%?;n==2; Dirac equation,;;'//
     .         '%?;n; spherical-harmonics,;;'//
     .         '%?;n>0; Im(eps(w)),;;%-1j%?;n<0; JDOS,;;'//
     .         '%?;n==1; eps^-1,;;%-1j%?;(n>=2); scr-rho-out,;;'//
     .         '%-1j%?;(n>=4);%b(model eps),;;'//
     .         '%?;n; Order-N:?,;;%-1j'//
     .         '%?;n==1;%2bEmbedded-Cluster,;;%-1j'//
     .         '%?;n==2;%2bVanderbuilt,;;'//
     .         '%?;n; APW basis,;;'//
     .         '%?;n; oveps,;;'//
     .         '%b %b',strn,scrwid,0,ctrl_lfrce,lrel,0,0, !ctrl_loptc,
     .         mod(lscr,10),0,pwmode,isw(oveps.ne.0))
              if (strn .ne. ' special:') call awrit0(strn,' ',-80,lgunit(k))
              if (prgnam .eq. 'LMF' .or. prgnam .eq. 'LMFA') then
                if(lxcf/=1.and.lxcf/=2.and.lxcf/=103) then
                  call rx('only XCFUN=1,2,103 are allowed now')
                endif
              endif
              if (prgnam(1:3) .eq. 'LMF') then
                call awrit1 ( ' special:%10p' // '%?;n==1; core-level-optics,;;' 
     .      // '%b %b' , strn , scrwid , 0 , 0 )
                if (strn .ne. ' special:') call awrit0(strn,' ',-80,lgunit(k))
              endif
   90       continue
          endif

!     ! --- takao embed contents in susite here. This is only for lmf and lmfgw.
          seref = 0d0
          do ib = 1, nsite
             ispec = v_ssite(ib)%spec
             if (ispec <= 0) cycle
             v_ssite(ispec)%clabel = v_sspec(ispec)%name
             seref = seref +     v_sspec(ispec)%eref
          enddo
          v_sham%seref=seref
          allocate(rv_a_opos(abs(3*nsite)))
          if (allocated(iv_a_oips)) deallocate(iv_a_oips) !v_sarry%
          allocate(iv_a_oips(abs(nsite)))
          do i_spackv=1,nsite
             rv_a_opos(3*(i_spackv-1)+1:3*(i_spackv-1)+3)= v_ssite( i_spackv )%pos
             iv_a_oips(i_spackv) =                        v_ssite( i_spackv )%spec
          enddo
          if (procid .eq. master) then
             if (iprint() .ge. 20) then
                if (lstsym .eq. 1) then
                   write(stdo,357) trim(sstrnsymg) !sstrn(i:j)
 357               format(/' Automatic symmetry finder turned off.  Use: ',a)
                elseif (lstsym .eq. 2) then
                   write(stdo,358)
 358               format(/' Symmetry operations suppressed')
                endif
             endif
          endif
          do i_spec=1,nspec
             call mpibc1_s_spec(v_sspec(i_spec),.false.,'rdctrl','sspec')
          enddo
!     ! Debugging printout ---
          if (io_help == 0 .and. io_show > 1) then
             print *, '---------- contents of sstrn ------------'
             print *, 'mix: ', trim(sstrnmix)
             print *, 'symg:', trim(sstrnsymg)
             call rx0('done show')
          endif
!     !     . slabl_,v_sbz,v_sctrl,v_sham,v_spot,v_slat,v_smix,v_sspec,v_ssite, !v_sstr,
!     !     . v_sarry,sstrnmix,sstrnsymg
          write(6,*) '--- End of readin ctrl file. ---'
          print *,'end of rdctrl2'
!!  -- End of readin ctrl file. -----------




           
#if MPI | MPIK
      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
      print *,' mpipid=',procid
      call MPI_BARRIER( MPI_COMM_WORLD, ierr )
#endif

      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc Hereafter is originally in lmv7.F
!     !--- Lattice setup:deformation ---
      !!v_* is going to be overwritten 
      call lattic() !this gives v_slat
      do ibas=1,v_sctrl%nbas
        v_ssite(ibas)%pos = rv_a_opos(3*(ibas-1)+1:3*(ibas-1)+3) !atom positions in the primitive cell
      enddo
!! --- We store all settings in module variables of m_rdctrl.
      print *,'end of lattic'
!!--- Generate symmetry operations; split species into classes  ---
      strn = 'find'
      if(len_trim(sstrnsymg)>0) strn=trim(sstrnsymg)
      if (cmdopt('--nosym',7,0,outs).or.cmdopt('--pdos',6,0,outs)) strn = ' '
!! when lmfgw, 1st digit of lc is zero--> no inversion added in mksym.
      lc = 10
      if (.not. prgnam == 'LMFGWD') lc = lc+2 ! inversion for lmf
      if (prgnam=='LMFA'.or.prgnam=='LMCHK') then
         lc=20
         if ( .not. iand(1,int(v_sctrl%noinv)) .ne.0 ) lc = lc+2
      endif   
c#ifndef LMFGWD 
!! iand(2**n,i) means 2**n ditit of i.
c      if( v_sctrl%noinv==0 ) lc = lc+2 ! inversion for lmf
c#endif
c      
c$$$!! Read AF case 2015, SPGGRAF category and SITE_ATOM_AF
c$$$      ipr10= iprint()>10 !this is only for master
c$$$      if(len_trim(symgaf)>0) then
c$$$        nbas = v_sctrl%nbas
c$$$        if(ipr10) then
c$$$          write(6,*)
c$$$          write(6,"(a)")       ' AF: ======================================== '
c$$$          write(6,"(a)")       ' AF: Antiferro mode: SPGGRAF='//trim(symgaf)
c$$$          write(6,"(a)")       ' AF:  (neglct waring in GENSYM) '
c$$$          write(6,"(a)")       ' AF:  Caution: I use v_sarry2=v_arry. Is this allowed in fortran?'
c$$$          do j=1,nbas
c$$$            write(6,"(a,2i3)") ' AF:  ibas,AF=',j,v_ssite(j)%iantiferro
c$$$          enddo
c$$$        endif
c$$$c        v_sarry2=v_sarry
c$$$        allocate(v_ssite2(nbas))
c$$$        v_ssite2=v_ssite
c$$$        inumaf=0
c$$$        do j=1,nbas
c$$$          do k=j,nbas
c$$$            if( v_ssite2(j)%iantiferro+v_ssite2(k)%iantiferro==0) then
c$$$              v_sarry2%iv_a_oips(k) =v_sarry2%iv_a_oips(j)
c$$$              inumaf=inumaf+1
c$$$              exit
c$$$            endif
c$$$          enddo
c$$$        enddo
c$$$        v_sctrl2 = v_sctrl
c$$$        v_sctrl2%nspec= v_sctrl%nspec - inumaf
c$$$        v_slat2=v_slat
c$$$        strn2=trim(strn)//' '//trim(symgaf)
c$$$        if(procid==master) call pshpr(60)
c$$$        call mksym(lc,slabl_,strn2 ,v_sctrl2,v_slat2,v_ssite2,v_sarry2)
c$$$        if(procid==master) call poppr()
c$$$        if(ipr10) write(6,"(a)")' AF: mksym, generator= SYMGRP+SYMGRPAF= '//trim(strn2)
c$$$        ngrpaf=v_slat2%nsgrp
c$$$        allocate(iclasstaf(nbas),symops_af(3,3,ngrpaf),ag_af(3,ngrpaf))
c$$$        call dcopy ( ngrpaf * 9 , v_slat2%rv_a_osymgr , 1 , symops_af , 1 )
c$$$        call dcopy ( ngrpaf * 3 , v_slat2%rv_a_oag ,    1 , ag_af , 1 )
c$$$        if(ipr10) write(6,"(a,i3)") ' AF: ngrpaf=',ngrpaf
c$$$        do ib=1,nbas
c$$$          iclasstaf(ib)=v_ssite2(ib)%class
c$$$        enddo
c$$$!! We use  symops_af, ngrpaf and so on in suham.F to make a table in sugw.F-gen_hamindex.F.
c$$$!!   It is passed through m_hamindex
c$$$!!   The table is for rotation of the self-energy file sigm.*, used in rdsigm2.F-hamfb3k-rotsig.
c$$$!!
c$$$c$$$        call mptauof ( symops_af , ngrpaf , v_slat2%plat , nbas , v_slat2%rv_a_opos , iclasstaf
c$$$c$$$     .   , miat , tiat , invgx , shtvg )
c$$$c$$$        do ig = 1,ngrpaf
c$$$c$$$          write(*,"(2i4,3e24.16)") ig, invgx(ig), shtvg(1:3,ig)
c$$$c$$$          do i=1,3
c$$$c$$$            write(*,"(3e24.16)") symops_af(i,1:3,ig)
c$$$c$$$          enddo
c$$$c$$$        enddo
c$$$          if(ipr10) write(6,"(a)")       ' AF: ===== end of AF section================= '
c$$$          if(ipr10) write(6,"(a)")       
c$$$       endif
       
!! ---------------------------
      if(procid==master) call pshpr(60)
      call mksym(lc,slabl_,strn,v_sctrl,v_slat,v_ssite) !,v_sarry)
      if(procid==master) call poppr
c      print *,'eeeeeeeeeeeeeeevvvvvvvvvvvvvvvvv1111qqqqqqqqqqqq'
c      call rx0(prgnam)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c$$$!! test1 ------------------
c$$$      ngrp=v_slat%nsgrp
c$$$      nbas=v_sctrl%nbas
c$$$      allocate(iclasst(nbas),invgx(ngrp),miat(nbas,ngrp),tiat(3,nbas,ngrp),shtvg(3,ngrp),symops(3,3,ngrp),ag(3,ngrp))
c$$$      call dcopy ( ngrp * 9 , v_slat%rv_a_osymgr , 1 , symops , 1 )
c$$$      call dcopy ( ngrp * 3 , v_slat%rv_a_oag , 1 , ag , 1 )
c$$$      write(*,*) 'ngrp=',ngrp
c$$$      do ib=1,nbas
c$$$        iclasst(ib)=v_ssite(ib)%class
c$$$      enddo
c$$$      call mptauof ( symops , ngrp , v_slat%plat , nbas , v_slat%rv_a_opos , iclasst 
c$$$     . , miat , tiat , invgx , shtvg )
c$$$      do ig = 1,ngrp
c$$$        write(*,"(2i4,3e24.16)") ig, invgx(ig), shtvg(1:3,ig)
c$$$        do i=1,3
c$$$          write(*,"(3e24.16)") symops(i,1:3,ig)
c$$$        enddo
c$$$      enddo
c
!!  Set default values for species data
      call defspc(v_sspec)
c#ifndef LMFGWD
c      ltet = iand(1,int(v_sctrl%lmet)) .ne.0 .or. iand(4+2+1,int(v_sctrl%ldos)) .ne.0
c!! iand(2**n,i) means 2**n ditit of i.
c!! iand(4+1,i) = 1,0,0,4,5,4,5,0, 1,0,0,4,5,4,5,0, ... for i=1,2,3,4,5,6 (I think).
c      call mkqp(v_sctrl,v_sbz,v_slat,ltet,F,1,-2)
c#endif
c      if(.not.prgnam == 'LMFGWD') then
      if(prgnam == 'LMF') then
         ltet = iand(1,int(v_sctrl%lmet)) .ne.0 .or. iand(4+2+1,int(v_sctrl%ldos)) .ne.0
!! iand(2**n,i) means 2**n ditit of i.
!! iand(4+1,i) = 1,0,0,4,5,4,5,0, 1,0,0,4,5,4,5,0, ... for i=1,2,3,4,5,6 (I think).
         call mkqp(v_sctrl,v_slat,ltet,F,1,-2)
      endif
!!    show and quit
      if( int(v_sctrl%quit) == 1 ) then
        write(6,"(a)") prgnam//"Q=SHOW encountered"
        call rx0(prgnam)
      endif
      end subroutine rdctrl2
