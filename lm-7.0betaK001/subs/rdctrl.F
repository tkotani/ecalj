      subroutine rdctrl(prgnam,wksize,toksw,ncount,mcount,vrsion,vn,vn2,
     .  pass2,slabl,osbz,osctrl,osham,ospot,oslat,osmix,osspec,ossite,
     .  osstr,osarry,osmove,ostb,sstrn)
C- Main input for ASA  and related programs
C ----------------------------------------------------------------------
Ci Inputs
Ci   prgnam:name of main program
Ci   wksize:size of global work array
Ci   toksw :set of strings specifying whether each token to be read are
Ci         :required, optional, to be ignored, or whether internal default
Ci         :is taken (see rdccat.f).
Ci   ncount:table of categories to read (see rdccat for the association
Ci         :of categories to table entries)
Ci   mcount:number of categories to read
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
Cu   29 Jul 08 Print zbak in header information
Cu   05 Jul 08 new default ham->pwpad
Cu   27 Jun 08 New defaults str->rmaxg,ivl spec->ehvl
Cu   05 Nov 07 some cleanup in preparation for new v7 input
Cu   08 Jun 07 New structure stb
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
C     implicit none
C ... Passed parameters
      logical pass2
      integer wksize,osbz,osctrl,osham,ospot,oslat,osmix,osspec,
     .  ossite,osstr,osarry,osmove,ostb,mcount,ncount(0:mcount)
      character  prgnam*(*), sstrn*(*)
      character toksw(0:30)*(*), vrsion*6
      double precision vn,vn2
C ... Heap
      integer w(1)
      common /w/ w
C ... Local parameters
      integer procid,master
      logical lgors,cmdopt,bittst,ltmp,parstr,ioorbp
      integer n0,nkap0
      parameter(n0=10,nkap0=3)
      integer TF(0:40),a2vec,bitand,fopna,getdig,i,i1mach,igets,iprint,
     .  iprt,irs(5),ishow,isw,ifi,ix(n0*nkap0),j,k,l,lasa,lbas,lcd,lfp,
     .  lfrce,lfrzw,lgunit,lham,lmet,lncol,lnsph,loptc,lordn,lpgf,lrel,
     .  lrs,lscr,lstsym,lsx,lsx1,lrsig,ltb,lxcf,mxrecs,nat,nbas,nclasp,
     .  nlibu,nclass,nfilin,nglob,nitmv,nkap,nl,nsite,nsp,nspc,nspec,
     .  nlmax,recln0,scrwid,stdo,stdl,stde,k1,k2,mpipid,nbasp,pwmode
      integer o,oclabl,offset,ohave,oics,opnu,opp,oqnu,osgw,osoptc,
     .  osordn,oves,owk
      parameter (mxrecs=10000,recln0=120)
      character strn*(recln0)
      character slabl(1)*8, fileid*64
      double precision dval,dglob,xx(n0),dalat,dgets,zbak(2)
      double precision ekap(6),prmint(20),oveps
      integer NULLI
      parameter (NULLI=-99999)
      save ishow
C     Default parameters for Bulirsch-Stoer integration
      data prmint /1d0,1d0,0d0,2d-4,7d0,11d0,2d0,4d0,6d0,8d0,12d0,16d0,
     .  24d0,32d0,48d0,64d0,96d0,3*0d0/

C     call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
      procid = mpipid(1)
      master = 0

      scrwid = 80
      stdo = lgunit(1)
      stdl = lgunit(2)
      stde = stdo

      if (pass2) then
        call upack('ctrl nbas nclass nspec nspin nl',w(osctrl),
     .    nsite,nclass,nspec,nsp,nl)
        nbas = nsite
        call upack('pot opnu oqnu oves',w(ospot),opnu,oqnu,oves,0,0)
        call upack2('pot opp osoptc',w(ospot),opp,osoptc)
        call upack('array oics oclabl nclasp ohave',w(osarry),
     .    oics,oclabl,nclasp,ohave,0)
        goto 5
      endif

C ... initialize sstrn
      call ustrn(w,0,len(sstrn),0,0,0)
      sstrn = ' '

C ... Allocate and default values for sbz
      call defdr(osbz,1000)
      call ubz(w(osbz),-1,-000,0,0,0)
      call redfrr(osbz,nint(dval(w(osbz),1)))
      call dpzero(xx,10)
      xx(1) = 10
      xx(3) = -1
      xx(5) = .01d0
      xx(9) = 1d-10
C     xx(10) = .005d0
      call pack5('bz n w efmax lmet semsh',w(osbz),0,5d-3,.5d0,1,xx)
      if (prgnam(1:2) .eq. 'TB')
     .  call pack3('bz efmax zval lmull',w(osbz),5d0,0d0,-1)
      if (prgnam(1:3) .eq. 'LMF') call pack1('bz efmax',w(osbz),2d0)
      call pack5('bz ndos ef def range',w(osbz),
     .  1001,0d0,.05d0,5d0,0d0)
C     call pack5('bz ndos ef def range',w(osbz),
C    .  1001,0d0,.05d0,9d0,0d0)
      xx(1) = -1
      xx(2) =  0
      call pack2('bz lio dosw',w(osbz),0,xx)

C ... Allocate and default values for sarray
      call defrr(osarry,1000)
      call uarray(w(osarry),-1,NULLI,0,0,0)
      call redfrr(osarry,nint(dval(w(osarry),1)))

C ... Allocate sctrl, specify some default values
      call defdr(osctrl,1000)
      call uctrl(w(osctrl),-1,0,0,0,0)
      call redfrr(osctrl,nint(dval(w(osctrl),1)))
      call uctrl(w(osctrl),-1,0,0,0,0)
      lasa = 1+4
      lfp = 0
      lbas = 0
C     For now, LMF => fp; no screening; nfp-style hamiltonian
      if (vrsion .eq. 'LMF-') then
        lfp = 1
        lbas = 3
      endif
C     Default for ES finder
      call ivset(irs,1,3,100)
      call pack5('ctrl rmines rmaxes nesabc',w(osctrl),1d0,2d0,irs,0,0)
C     Default: cd represented in plane waves
      lcd = 0
      if (lfp .eq. 1) lcd = 4
      call pack5('ctrl lasa lfp lbas lcd',w(osctrl),lasa,lfp,lbas,lcd,0)
      lrel = 1
      nsp = 1
      lmet = 1+2+4+8
      call pack5('ctrl lmet lqp lrel nspin',w(osctrl),lmet,2,lrel,nsp,0)
      nitmv = 1
      call ivset(irs,1,5,0)
      irs(1) = 1
      irs(2) = 1
      if (cmdopt('--rs=',5,0,strn)) then
        i = 5
        j = a2vec(strn,len(strn),i,2,', ',2,2,5,ix,irs)
      endif
      irs(1) = mod(mod(irs(1),10),4) + 4*getdig(irs(1),1,10)
     .       + 8*getdig(irs(1),1,100)
      lrs = 1*mod(irs(1),8)+8*irs(2)+32*irs(3)+64*irs(4)+128*irs(5)
     .    + 256*mod(irs(1)/8,2)
      call pack2('ctrl nitmv lrs',w(osctrl),nitmv,lrs)
      dalat = 0
      ix(1) = 80
      ix(2) = 2
      ix(3) = 0
      ix(4) = 0
      if (prgnam .eq. 'LMPG' .or. prgnam .eq. 'TBPG') ix(3) = 1
      nl = 3
      call pack5('ctrl lxcf nl lpgf maxit smalit',w(osctrl),2,nl,ix(3),
     .  1,ix)
      call dvset(xx,1,3,1d-4)
      xx(3) = 0d0
      if (prgnam .eq. 'LMF') xx(3) = 1d-4
      call pack5('ctrl tol ltb',w(osctrl),xx,256,0,0,0)
      zbak(1) = 0d0
      zbak(2) = 0d0
      call pack1('ctrl zbak',w(osctrl),zbak)
      if (vrsion .eq. 'MM-') call pack1('ctrl lncol',w(osctrl),17)
C     Default sphere overlap constraints
      xx(1) = 0.16d0
      xx(2) = 0.18d0
      xx(3) = 0.20d0
      xx(4) = 0.40d0
      xx(5) = 0.45d0
      xx(6) = 0.50d0
      call pack5('ctrl sclwsr omax1 omax2',w(osctrl),0d0,xx,xx(4),0,0)
C     call shstru('ctrl',w(osctrl),0,0)
C     Default lattice deformation mode
      call dvset(xx,1,6,1d0)
      xx(1) = 0
      call pack5('ctrl defm',w(osctrl),xx,0,0,0,0)

C ... Allocate and default values for slat
      call defdr(oslat,1000)
      call ulat(w(oslat),-1,0,0,0,0)
      call redfrr(oslat,nint(dval(w(oslat),1)))
      xx(1) = 0
      xx(2) = 0
      xx(3) = 1
      xx(4) = 1
      call pack5('lat as nkdmx nkqmx tol gam',w(oslat),
     .  2d0,800,800,1d-8,xx)
      call pack1('lat tolft',w(oslat),1d-6)
C ... Allocate and initialize sham
      call defdr(osham,1000)
      call uham(w(osham),-1,-000,0,0,0)
      call redfrr(osham,nint(dval(w(osham),1)))
      call dpzero(xx,n0)
C     xx(5) = 2
C     xx(7) = 0.08d0
      xx(1) = 3
      xx(5) = 2
      xx(6) = 0.02d0
      xx(7) = 0.06d0
C     xx(8) = .001d0
C     xx(9) = .0013d0
C     xx(10)= .0012d0
      call pack5('ham rsrnge sigp rsstol',w(osham),5d0,xx,5d-6,0,0)
      call pack2('ham npwpad pwemin',w(osham),-1,0d0)
C     call shstru('ham',w(osham),0,0)
C     stop

C ... Allocate and initialize spot
      call defdr(ospot,1000)
      call upot(w(ospot),-1,-000,0,0,0)
      call redfrr(ospot,nint(dval(w(ospot),1)))
C ... Allocate and initialize smix
      call defrr(osmix,1000)
      call umix(w(osmix),-1,-000,0,0,0)
      call redfrr(osmix,nint(dval(w(osmix),1)))
      call spacks(1,'mix fn',w(osmix),'mix',0,0)
      call spacks(1,'mix r',w(osmix),' ',0,0)
      xx(1) = 1
      xx(2) = 1
      xx(3) = 0
      call pack5('mix b bv wc w mmix',w(osmix),1d0,1d0,-1d0,xx,-1)
      call pack1('mix nsave',w(osmix),8)
C     call shstru('mix',w(osmix),0,0)
C     stop
C ... Allocate and initialize soptc
      call defrr(osoptc,1000)
      call uoptic(w(osoptc),-1,0,0,0,0)
      call redfrr(osoptc,nint(dval(w(osoptc),1)))
      call pack1('pot osoptc',w(ospot),osoptc)
      call pack1('optic ne',w(osoptc),501)
      call pack1('optic window:2',w(osoptc),1d0)
C ... Allocate and initialize sordn
      call defrr(osordn,1000)
      call uordn(w(osordn),-1,0,0,0,0)
      call redfrr(osordn,nint(dval(w(osordn),1)))
      call pack1('array osordn',w(osarry),osordn)
C ... Allocate and initialize smove
      call defrr(osmove,1000)
      call umove(w(osmove),-1,0,0,0,0)
      call redfrr(osmove,nint(dval(w(osmove),1)))
C     set defaults to gyro, Bulirsch-Stoer integration
      call pack1('move gyro',w(osmove),2d0)
      call pack1('move prmint',w(osmove),prmint)
C     call shstru('move',w(osmove),0,0)
C     call pack1('array osmove',w(osarry),osmove)
C ... Allocate and initialize sgw
      call defrr(osgw,1000)
      call ugw(w(osgw),-1,0,0,0,0)
      call redfrr(osgw,nint(dval(w(osgw),1)))
      call pack5('gw gcutb gcutx qoffp nband gsmear',w(osgw),2.7d0,
     .  2.2d0,1d0,9999,.003d0)
C     Default : Faleev's fast energy integration
      xx(1) = .01d0
      xx(2) = .04d0
      call pack5('gw lgw delre ecuts',w(osgw),1,xx,2.5d0,0,0)
C     xx(1) = .04d0
C     xx(2) = .01d0
C     call pack2('gw lgw delre',w(osgw),0,xx)
      call pack5('gw nime delre deltax deltaw pbtol',w(osgw),6,xx,
     .  -1d-4,.02d0,1d-3)
      call spacks(1,'gw pb1',w(osgw),'111',0,0)
      call spacks(1,'gw pb2',w(osgw),'1111',0,0)
      call pack1('pot osgw',w(ospot),osgw)

C ... Allocate sstr, initialize with some default values
      call defrr(osstr,1000)
      call ustr(w(osstr),-1,0,0,0,0)
      call redfrr(osstr,nint(dval(w(osstr),1)))
      xx(1) = 0
      xx(2) = -1
      xx(3) = -2.3d0
      if (prgnam .eq. 'TBE') then
        call pack5('str nkaps rmax rfit kaps lmaxw',w(osstr),1,
     .    0d0,.8d0,xx,-1)
      else
        call pack5('str nkaps rmax rfit kaps lmaxw',w(osstr),1,
     .    0.0d0,.8d0,xx,-1)
        call pack2('str rmaxg ivl',w(osstr),0.7d0,0)
      endif
      call pack2('str loka drwats',w(osstr),1,.1d0)
C     call shstru('str',w(osstr),0,0)
C     stop

C ... Allocate and default values for stb
      call defrr(ostb,1000)
      call utb(w(ostb),-1,NULLI,0,0,0)
      call redfrr(ostb,nint(dval(w(ostb),1)))
      call pack1('tb alam',w(ostb),.001d0)
      call pack1('tb alsc',w(ostb),10d0)
C     call shstru('tb',w(ostb),0,0)
C     stop

C ... Other initialization
      ossite = 1
      osspec = 1
      nbas = 1
      nspec = 1
      nsite = 1
      ishow = 1

    5 nfilin = 0
      if (procid .eq. master) then
      if (cmdopt('--input',7,0,strn)) then
        ishow = 0
      else
        nfilin = fopna('CTRL',-1,1)
      endif
      endif

      offset = wksize*i1mach(6) - mxrecs*recln0 - 100

C --- Read from mcount input categories ---
      do  10  i = 0, mcount
        if (.not. pass2 .and. ncount(i) .eq. 4) then
          call numsyv(j)
          call pack1('ctrl nvario',w(osctrl),j)
        endif
        if (pass2 .and. ncount(i).ne.15) goto 10
        k = mpipid(2)
C       call MPI_BARRIER( MPI_COMM_WORLD, ierr )

        call rdccat(mxrecs,recln0,nfilin,offset,ishow,ncount(i),pass2,
     .    vrsion,vn,prgnam,TF,toksw(max(ncount(i),0)),strn,nspec,nsite,
     .    nbas,nl,nsp,dalat,slabl,nclasp,w(opnu),w(oqnu),w(opp),w(oics),
     .    w(oclabl),w(ohave),w(oves),w(osbz),w(osctrl),
     .    w(osham),w(ospot),w(oslat),w(osmix),w(osoptc),w(osordn),
     .    w(osspec),w(ossite),w(osstr),w(osmove),w(osgw),w(ostb),sstrn)
C       really belongs in the call to rdccat
        call upack1('ctrl nbasp',w(osctrl),nbasp)

C       Repeat, in case --rs specified in CMD
        if (ncount(i) .eq. 25 .and. cmdopt('--rs=',5,0,strn)) then
          k = 5
          j = a2vec(strn,len(strn),k,2,', ',2,2,5,ix,irs)
          irs(1) = mod(mod(irs(1),10),4) + 4*getdig(irs(1),1,10)
     .           + 8*getdig(irs(1),1,100)
          lrs = 1*mod(irs(1),8)+8*irs(2)+32*irs(3)+64*irs(4)+128*irs(5)
     .        + 256*mod(irs(1)/8,2)
          call pack1('ctrl lrs',w(osctrl),lrs)
        endif

        call defi(o,1)
        if (o .gt. offset) call rx('increase wksize')
        call rlse(o)
C   ... Some structure-dependent initialization
        if (ncount(i) .eq. 6) then
          call pack5('ctrl nsite nbas nbasp nspec',w(osctrl),
     .      nsite,nbas,nbasp,nspec,0)
C     ... Allocate and initialize sspec and ssite
          call defrr(osspec,1000)
          call uspec(w(osspec),-1,0,1,0,0,0)
          call redfrr(osspec,nspec*nint(dval(w(osspec),1)))
          call uspec(w(osspec),-1,-000,nspec,0,0,0)
          call spackv(1,'spec norp',w(osspec),1,nspec,2)
C     ... Init default parameters for local orbitals
          call spackv(1,'spec rs3',w(osspec),1,nspec,1d0)
          call spackv(1,'spec eh3',w(osspec),1,nspec,-0.5d0)
          call dvset(xx,1,10,-.5d0)
          call spackv(1,'spec ehvl',w(osspec),1,nspec,xx)
          call spackv(1,'spec vmtz',w(osspec),1,nspec,-0.5d0)
C     ... Init default for sigma,lmxa,lmxf
          call dvset(xx,1,10,-1d0)
C     ... First set lmxa to 4 to set def. hcr to l=4
          call spackv(1,'spec lmxa',w(osspec),1,nspec,4)
          call spackv(1,'spec lmxpb',w(osspec),1,nspec,3)
          call spackv(1,'spec hcr',w(osspec),1,nspec,xx)
          call spackv(1,'spec lmxa',w(osspec),1,nspec,nl-1)
          call spackv(1,'spec lmxf',w(osspec),1,nspec,2*nl-2)
          xx(1) = -1
          xx(2) = 0
          call spackv(1,'spec coreq',w(osspec),1,nspec,xx)
          if (prgnam(1:2) .eq. 'TB')
     .    call spackv(1,'spec rmt',w(osspec),1,nspec,1d0)
          do  12  j = 1, nspec
C           call spacks(1,'spec pb1',w(osspec),'111     ',j,j)
C           call spacks(1,'spec pb2',w(osspec),'1111    ',j,j)
            call spacks(1,'spec pb1',w(osspec),'        ',j,j)
            call spacks(1,'spec pb2',w(osspec),'        ',j,j)
            call spacks(1,'spec coreh',w(osspec),'        ',j,j)
   12     continue
C         call shstru('spec',w(osspec),1,nspec)
C         stop

C     ... NB: essential that ssite is LAST array allocated here
          call defrr(ossite,1000)
          call usite(w(ossite),-1,0,1,0,0,0)
          call redfrr(ossite,nsite*nint(dval(w(ossite),1)))
          call usite(w(ossite),-1,0,nsite,0,0,0)
          call ivset(ix,1,n0*nkap0,1)
          call spackv(1,'site relax',w(ossite),1,nsite,ix)
C     ... Default lmxa
          call spackv(1,'spec lmxa',w(osspec),1,nspec,nl-1)
          call spackv(1,'spec etf',w(osspec),1,nspec,-1d0)
C     ... When no automatic downfolding, set default idxdn to 1
          if (.not. lgors('ctrl lham,4',w(osctrl)))
     .      call spackv(1,'spec idxdn',w(osspec),1,nspec,ix)
C        elseif (ncount(i) .eq. 8) then
C.          if (cmdopt('--rpos=',7,0,strn)) then
C            call defrr(opos,3*nsite)
C            call spackv(10,'site pos',w(ossite),1,nsite,w(opos))
C            call iopos(.false.,-1,strn(8:),nbas,w(opos))
C            call spackv(11,'site pos',w(ossite),1,nsite,w(opos))
C            call rlse(opos)
C          endif
        endif

   10 continue
      if (ishow .eq. 0) call cexit(0,1)
      if (pass2) return

C --- Miscellaneous cleanup and initialization ---
      if (lgors('ctrl lgen3',w(osctrl))) then
        nkap = dgets('str skmsh:1.1',w(osstr))
        xx(1) = dgets('str skmsh:2.2',w(osstr))
        xx(2) = dgets('str skmsh:3.3',w(osstr))
        call gausq(nkap,xx(1),xx(2),ekap,xx(3),0,0)
      endif

C ... Reset default kmxa if PW turned on (should be done elsewhere)
C      call upack2('ham pwmode pwemax',w(osham),pwmode,pwemax)
C      if (pwmode .gt. 0 .and. pwemax .gt. 3) then
C      do  i = 1, nspec
C        print *,i
C      enddo
C      endif
C      stop

C ... Read the basis from the atm file
      if (cmdopt('--rdbasp',8,0,strn)) then
        fileid = 'basp'
        if (strn(9:12) .eq. ':fn=') then
          fileid = strn(13:)
        else
        endif

C       Number of envelope functions before modifications
C       Takao says: I think xx(n0*nkap0*nspec) is required.
        call rx('check dimension of xx for this call.')
        call uspecb(0,1,w(osspec),1,nspec,xx,xx,xx,k1)

C       Read new parms; return in k2 max no. env. fns/site read 
        call strip(fileid,i,j)
        ifi = fopna(fileid(1:j),-1,0)
        rewind ifi
        call upack1('ctrl nspec',w(osctrl),nspec)
        if (.not. ioorbp(111,2,1,nspec,w(osspec),k2,ifi))
     .    call rxs2('lmfp: failed to find BASIS: token in file "',
     .    fileid(1:j),'"')
        k2 = mod(k2,10)
        call fclr(' ',ifi)

C       File read cause number of envelope functions to increase?
C       If so, increment nkaph by 1; keep old nkaph in k1 for later
        if (k2 .gt. k1) then
          k1 = nglob('nkaph')
          xx(1) = dglob('nkaph',dble(k1+1),1)
          k = nglob('nkaph')
          call uspecb(0,-1,w(osspec),1,nspec,xx,xx,xx,xx)
          nlmax = igets('ctrl nl',w(osctrl))**2
          xx(1) = dglob('mxorb',dble(k)*nlmax,1)
        endif
      endif

      call pack1('lat alat',w(oslat),dgets('lat alat',w(oslat))+dalat)
      call upack('ctrl lncol lsx lscr lmet lrel',
     .  w(osctrl),lncol,lsx,lscr,lmet,lrel)
      lsx1 = mod(lsx,2)
      lrsig = igets('ham lsig',w(osham))
      call upack('ctrl lordn loptc',w(osctrl),lordn,loptc,0,0,0)
      lpgf = igets('ctrl lpgf',w(osctrl))
C     lcgf = igets('ctrl lcgf',w(osctrl))
      lnsph = isw(lgors('ctrl lasa,32',w(osctrl)))
      if (lnsph .ne. 0) call pack1('mix nsave',w(osmix),3)
C     replicate ctrl->lncol in ham->lncol
      call pack1('ham lncol',w(osham),lncol)
C ... Use true spherical harmonics
      if (bittst(lncol,4) .or. loptc .ne. 0 .or. lrel .eq. 2)
     .  call lsets('ctrl lham',w(osctrl),.true.,256)

C ... Dirac equation requires spin polarization
      if (nsp.eq.1 .and. igets('ctrl lrel',w(osctrl)).eq.2) then
        call rx('rdccat: Dirac equation requires NSPIN=2')
      endif

C ... Suppress symmetry operations for:
      lstsym = 0
C     noncollinear magnetism
      if (lncol .ne. 0) lstsym = 1
C     molecular dynamics
      if (dgets('ctrl mdprm',w(osctrl)) .eq. 1 .or.
     .    dgets('ctrl mdprm',w(osctrl)) .eq. 2 .or.
     .    dgets('ctrl mdprm',w(osctrl)) .eq. 3) lstsym = 2
C     core-level optics
      if (cmdopt('--cls',5,0,strn)) lstsym = 2
      if (cmdopt('--nosym',7,0,strn)) lstsym = 2
      if (lstsym .ne. 0) then
        call upacks('strn symg',i,j)
        if (j .ge. i) then
C         Discard all symops
          if (lstsym .eq. 2) then
            sstrn(i:j) = 'e'
C         Suppress automatic finder
          else
            k = i-1
            ltmp = parstr(sstrn,'find ',j-3,5,' ',k,l)
            if (ltmp) then
              sstrn(k+1:l) = 'e'
            endif
          endif
C       No symops supplied: suppress automatic finder
        else
          call lstra('strn symg',i,o,j)
          call ustrn(w,-o,1,i,ix,3)
          sstrn(i:i) = 'e'
        endif
        call lsets('ctrl lqp',w(osctrl),1,1)
      endif
C ... Suppress inversion when noncollinear magnetism, SX, NLO
      if (lncol+lsx1 .ne. 0 .or. loptc .ge. 10)
     .  call lsets('ctrl lqp',w(osctrl),1,1)
      if (mod(lscr,10) .eq. 1 .and. prgnam .eq. 'LM')
     .  call lsets('ctrl lqp',w(osctrl),1,1)
C ... No Euler angle motion unless forces are set
C     if (.not. bittst(lncol,16)) call pack1('ctrl sdmod',w(osctrl),-1)

C ... Special pgf initialization
      if (lpgf .ne. 0) then
        call lsets('ctrl lmet',w(osctrl),.false.,2)
        call redfrr(ossite,3*nsite*nint(dval(w(ossite),1)))
      endif

C ... Setup for idxdn ... ctrl->lham,4 is automatic downfolding switch
      j = igets('ctrl lham,4',w(osctrl))/4
      j = 2*(1-j)
C     No screening => no downfolding; also lmxb<l<=lmxa => 'high'
C     Probably ought to have lmxb<l<=lmxa => 'high' always
      if (lgors('ctrl lbas,1',w(osctrl))) j = 3
C     nfp-style basis:
      if (lgors('ctrl lbas,2',w(osctrl))) j = j+10
      call suidx(nglob('nkaph'),j,nspec,w(osspec))

C --- ASA-specific initialization ---
C ... nothing ASA-specific

C ... Set some global variables
      xx(1) = dglob('nspec',dble(nspec),1)
      xx(1) = dglob('nbas',dble(nbas),1)
      xx(1) = dglob('nbasp',dble(nbasp),1)
      xx(1) = dglob('nsp',dble(nsp),1)
      xx(1) = dglob('nl',dble(nl),1)
      xx(1) = dglob('avw',dgets('lat avw',w(oslat)),1)
      xx(1) = dglob('lrel',dble(lrel),1)
      lxcf  = igets('ctrl lxcf',w(osctrl))
      xx(1) = dglob('lxcf',dble(lxcf),1)
      xx(1) = dglob('stdo',dble(stdo),1)
      xx(1) = dglob('stdl',dble(stdl),1)
      xx(1) = dglob('stde',dble(stde),1)
      nspc = 1
      if (bitand(lncol,1+2+4) .ne. 0) nspc = 2
      xx(1) = dglob('nspc',dble(nspc),1)
C     nat = number of real atoms as nbas - # sites w/ floating orbitals
      if (procid .eq. master) then
      nat = nbas
      do  i = 1, nbas
        call upack('site spec',w(ossite),i,j,0,0,0)
        call upack('spec lmxa',w(osspec),j,l,0,0,0)
        if (l .eq. -1) nat = nat-1
        call spacks(1,'site clabel',w(ossite),' ',i,i)
      enddo
      endif

      call mpibc1(nat,1,2,0,'rdctrl','nat')
      xx(1) = dglob('nat',dble(nat),1)

C ... Set modep
      ix(1) = 2
      ix(2) = 2
      ix(3) = 2
      if (lgors('ctrl lpgf,-1',w(osctrl))) ix(3) = 0
      call pack1('ctrl modep',w(osctrl),ix)
      call pack1('ham lxcf',w(osham),lxcf)

C ... End of copy

C ... Count LDA+U blocks (NB: doesn't count padding layers in PGF case)
      call defi(owk,-nbas)
      call pshpr(0)
      call suldau(nbas,w(osspec),w(ossite),nlibu,k,w(owk))
      call pack5('ham nlibu lmaxu',w(osham),nlibu,k,0,0,0)
      call poppr
      call rlse(owk)

C --- Printout ---
      if (iprint() .ge. 20 .and. procid .eq. master) then
        do  90  k = 1, 2

        strn = '  '//prgnam
        call awrit5(
     .    '%N %a:%12palat = %;5d  '//
     .    'nbas = %i%?#n#%-1j(+%i)##%?#n#%-1j+%i##'//
     .    '  nspec = %i',strn,scrwid,0,dgets('lat alat',w(oslat)),
     .    nat,nbas-nat,nsite-nbas,nspec)
        call awrit6('%a  vn %,2d  verb %i%?#n#,%i#%j#%?#n#,%i#%j#',
     .    strn,scrwid,-lgunit(k),vn-vn2/1d5,
     .    iprint(),iprint()-iprt(1),iprt(1),iprt(2)-iprt(1),iprt(2))
        call upack('ctrl lxcf lfrce',w(osctrl),lxcf,lfrce,0,0,0)
        i = igets('ctrl lham,256',w(osctrl))
        lfrzw = isw(lgors('ctrl lbas,16',w(osctrl)))
        call upack1('ctrl zbak',w(osctrl),zbak)


        call upack2('ham pwmode oveps',w(osham),pwmode,oveps)

        call awrit8(' special:%10p'//
     .    '%?;n; forces,;;'//
     .    '%?;n==2; Dirac equation,;;'//
     .    '%?;n; spherical-harmonics,;;'//
     .    '%?;n>0; Im(eps(w)),;;%-1j%?;n<0; JDOS,;;'//
     .    '%?;n==1; eps^-1,;;%-1j%?;(n>=2); scr-rho-out,;;'//
     .    '%-1j%?;(n>=4);%b(model eps),;;'//
     .    '%?;n; Order-N:?,;;%-1j'//
     .    '%?;n==1;%2bEmbedded-Cluster,;;%-1j'//
     .    '%?;n==2;%2bVanderbuilt,;;'//
     .    '%?;n; APW basis,;;'//
     .    '%?;n; oveps,;;',
     .    strn,scrwid,0,
     .    lfrce,lrel,i,loptc,mod(lscr,10),lordn,pwmode,isw(oveps.ne.0))
        if (prgnam(1:3) .eq. 'LMF') then
          call awrit1('%a'//
     .      '%?;n==1; core-level-optics,;;',strn,scrwid,0,
     .      igets('optic cls',w(osoptc)))
        endif
        call awrit8('%a'//
     .    '%?;n==1; core-level-optics,;;'//
     .    '%?;n; zbak=%d,;%j;'//
     .    '%b %b',
     .    strn,scrwid,0,
     .    isw(prgnam(1:3).eq.'LMF')*igets('optic cls',w(osoptc)),
     .    isw(zbak.ne.0),zbak,0,0,0,0,0)
        if (strn .ne. ' special') call awrit0(strn,' ',-80,lgunit(k))

        if (prgnam(1:2) .ne. 'TB') then
          call awrit1('%x pot:%10p'//'%?;n==0; non-rel,;.;',strn,
     .      scrwid,0,lrel)
          call awrit8('%a%?;%c==,;;%b%a%b;%-1j'//
     .    '%?;n>1; spin-pol,;;'//
     .    '%?;n; LDA+U,;;'//
     .    '%?;n==1; XC:CA,;;%-1j'//
     .    '%?;n==2; XC:BH,;;%-1j'//
     .    '%?;n==3; XC:PW91,;;%-1j'//
     .    '%?;n==4; XC:PBE,;;'//
     .    '%?;n==1;%b+LMH(gga),;;%-1j'//
     .    '%?;n==2;%b+PW91(gga),;;%-1j'//
     .    '%?;n==3;%b+PBE(gga),;;%-1j'//
     .    '%?;n==4;%b+Becke(gga),;;'//
     .    '%?;n; frozen-wave-functions,;;'//
     .    '%?;n==1; nsph-mpol,;;'//
     .    '%?;n; read Sigma,;;'//
     .    '%?;n; make SX,;;'//
     .    '%b %b',
     .      strn,scrwid,0,nsp,nlibu,mod(lxcf,10),lxcf/10,lfrzw,lnsph,
     .      lrsig,lsx)
          call awrit0(strn,' ',-80,lgunit(k))
          call isanrg(mod(lxcf,10),0,4,prgnam,'XC functional',.true.)
          call isanrg(lxcf/10,0,4,prgnam,'GGA functional',.true.)
          if (mod(lxcf,10) .eq. 2 .and. lxcf/10 .gt. 1 .or.
     .      mod(lxcf,10) .ne. 2 .and. lxcf/10 .eq. 1) call info0
     .      (10,0,0,'%10f(warning) mixing incompatible functionals')
        endif

        if (bitand(lncol,1+2+4+8+32+64) .ne. 0) then
          i = igets('ctrl sdmod',w(osctrl))
          if (.not. bittst(lncol,16)) i = -1
          call awrit8(' noncoll: '//
     .    '%?;n; Non-coll,;;'//
     .    '%?;n; Spin-orbit,;;'//
     .    '%?;n; LzSz,;;%-1j%?;(n>32);%b+L.S(pert),;;'//
     .    '%?;n; B-field,;;'//
     .    '%?;n; spin-spiral,;;'//
     .    '%?;n; mag-forces:;;'//
     .    '%?;(n%10)<2&(n%10)>=0; relax,;;'//
     .    '%?;(n%10)>1; spin-dynamics,;;'//
     .    '%a%b %a',strn,scrwid,lgunit(k),
     .      bitand(lncol,1),bitand(lncol,4),bitand(lncol,32+64),
     .      bitand(lncol,8),bitand(lncol,2),bitand(lncol,16),i,i)
        endif

        ltmp = prgnam .eq. 'LM' .or. prgnam .eq. 'LMPG'
     .    .or. prgnam .eq. 'LMGF'
        if (ltmp) then
          call upack('ctrl lasa lham lcd',w(osctrl),lasa,lham,lcd,0,0)
          call awrit8('%x asa:%10p'//
     .      '%?;n; no-ccor,;;'//
     .      '%?;n; gam-rep,;;'//
     .      '%?;n; two-c-H,;;'//
     .      '%?;n;%b + pert-ev,;;'//
     .      '%?;n; map,;;'//
     .      '%?;n; indep-vmix,;;'//
     .      '%?;n; mt-corr,;;'//
     .      '%?;n; frozen core,;;'//
     .      '%a%b ',strn,scrwid,0,
     .      isw(bitand(lasa,4).eq.0),
     .      bitand(lham,128),
     .      bitand(lham,3),bitand(lham,2),
     .      bitand(lasa,16),igets('mix lxpot,3',w(osmix)),
     .      bitand(lasa,64),bitand(lcd,1))
          if (strn .ne. ' asa') call awrit0(strn,' ',-80,lgunit(k))
        endif

        if (prgnam(1:2) .eq. 'TB') then
          call upack2('ctrl mdprm ltb',w(osctrl),xx,ltb)
          ix(1) = nint(xx(1))
          ix(2) = nint(xx(2))
          call upack1('str rmax',w(osstr),xx)
          call awrit1('%x TB: %11prmaxh = %d,',strn,scrwid,0,xx)
          call awrit6('%a'//
     .      '%?;n; s-c multipoles,;;'//
     .      '%?;n; s-c multipoles: MRS theory,;;'//
     .      '%?;n; read del,;;'//
     .      '%?;n; n-orthog TB,;;'//
     .      '%?;n==4; m-stat: Conj. grad.;;%-1j'//
     .      '%?;n==5; m-stat: F-P;;%-1j'//
     .      '%?;n==6; m-stat: Broy;;%-1j'//
     .      '%?;n==1; MD (NVE),;;%-1j'//
     .      '%?;n==2; MD (NVT),;;%-1j'//
     .      '%?;n==3; MD (NPT),;;'//
     .      '%?;n; pair-only,;;',
     .      strn,scrwid,0,bitand(ltb,2**15),bitand(ltb,2**12),
     .      bitand(ltb,2**16),bitand(ltb,1),ix(1),bitand(ltb,512))
          if (ix(1) .ge. 4) then
            call awrit2('%a'//
     .      '%?;n; rlx-vol,;;'//
     .      '%?;n; i/o-hess,;;',
     .      strn,scrwid,0,bitand(ix(2),1),bitand(ix(2),2))
          endif
          call awrit5('%a'//
     .      '%?;n; trh,;;'//
     .      '%?;n; rho,;;'//
     .      '%?;n; spin-orb,;;'//
     .      '%?;n; crysf,;;'//
     .      '%?;n;%b+ovlp-cf,;;'//
     .      '%a%b %a',strn,-scrwid,-lgunit(k),
     .      bitand(ltb,1024),bitand(ltb,2048),
     .      bitand(lncol,4),bitand(ltb,2),bitand(ltb,4))
        endif

        call awrit7(' bz:%10p'//
     .    '%?;n; metal,; nonmetal,;'//
     .    '%-1j%?;n>1;%b(%-1j%i),;;'//
     .    '%?;n; tetra,;;'//
     .    '%?;n; get-qp,;;'//
     .    '%?;n; invit,;;'//
     .    '%?;n; dens-mat,;;'//
     .    '%?;(n>0); %-1jmull=%i,;;'//
     .    '%?;n; fixed-spin-mom,;;%b ',
     .    strn,scrwid,lgunit(k),
     .    igets('bz lmet',w(osbz)),
     .    bitand(lmet,2),igets('bz lio,1',w(osbz)),
     .    igets('ctrl lqp,2',w(osctrl)),
     .    igets('bz lio,8',w(osbz)),
     .    igets('bz lmull',w(osbz)),
     .    isw(dgets('bz fsmom',w(osbz)).ne.0))

   90   continue
      endif

C --- Sanity checks and other initialization ---
      if (lgors('ctrl ldos,8',w(osctrl)) .and.
     .    lgors('ctrl ldos,4+2',w(osctrl)))
     .  call rx('inconsistent DOS options')
      call rxx(lgors('ctrl lstonr,-1',w(osctrl)) .and. nsp .eq. 2,
     .  'Stoner model not compatible with nsp=2')
      call rxx(lgors('ctrl lstonr,-1',w(osctrl)) .and.
     . lgors('ctrl lham,4',w(osctrl)),'Stoner not compatible with ADNF')

C --- Check and order principal layers, and sites by PL ---
      call susite(w(osctrl),w(osham),w(ospot),w(osarry),w(oslat),
     .  w(osspec),ossite)

      call pack1('ham lham',w(osham),igets('ctrl lham',w(osctrl)))

      if (procid .eq. master) then
      if (iprint() .ge. 20) then
        if (lstsym .eq. 1) then
          call upacks('strn symg',i,j)
          write(stdo,357) sstrn(i:j)
  357     format(/' Automatic symmetry finder turned off.  Use: ',a)
        elseif (lstsym .eq. 2) then
          write(stdo,358)
  358     format(/' Symmetry operations suppressed')
        endif
      endif
      endif

C     Broadcast spec structure to fix floating orbitals case
      call mpibc1(w(osspec),nspec*nint(dval(w(osspec),1)),4,0,
     .  'rdctrl','sspec')

C      ifi = fopna('out',-1,0)
C      k = nspec*nint(dval(w(osspec),1))
C      print *, procid,k
C      call ywrm(0,'sspec',1,ifi,'(9f15.6)',w(osspec),0,k,k,1)
C      k = nsite*nint(dval(w(ossite),1))
C      print *, procid,k
C      call ywrm(0,'ssite',1,ifi,'(9f15.6)',w(ossite),0,k,k,1)
C      k = nint(dval(w(osbz),1))
C      print *, procid,k
C      call ywrm(0,'sbz',1,ifi,'(9f15.6)',w(osbz),0,k,k,1)
C      k = nint(dval(w(ospot),1))
C      print *, procid,k
C      call ywrm(0,'spot',1,ifi,'(9f15.6)',w(ospot),0,k,k,1)
C      k = nint(dval(w(oslat),1))
C      print *, procid,k
C      call ywrm(0,'slat',1,ifi,'(9f15.6)',w(oslat),0,k,k,1)
C      k = nint(dval(w(osham),1))
C      print *, procid,k
C      call ywrm(0,'sham',1,ifi,'(9f15.6)',w(osham),0,k,k,1)
C      print *, procid
C      call shstru('spec',w(osspec),nspec,nspec)
C      call rx0('done')

      if (cmdopt('--show=',7,0,fileid)) then
      i = 7
      i = a2vec(fileid,len(fileid),i,2,' ',1,1,1,ix,j)
      if (j .gt. 1) then
        print *, '---------- contents of sstrn ------------'
        call upacks('strn amix',i,j)
        print *, 'amix:', sstrn(i:j)
        call upacks('strn mix',i,j)
        print *, 'mix:', sstrn(i:j)
        call upacks('strn gfopt',i,j)
        print *, 'gfopt:', sstrn(i:j)
        call upacks('strn mmham',i,j)
        print *, 'mmham:', sstrn(i:j)
        call upacks('strn symg',i,j)
        print *, 'symg:', sstrn(i:j)
        call upacks('strn sxopt',i,j)
        print *, 'sxopt:', sstrn(i:j)

        call shstru('ctrl',w(osctrl),0,0)
        call shstru('bz',w(osbz),0,0)
        call shstru('lat',w(oslat),0,0)
        call shstru('array',w(osarry),0,0)
        call shstru('ham',w(osham),0,0)
C       call shstru('pot',w(ospot),0,0)
        call shstru('mix',w(osmix),0,0)
        call shstru('move',w(osmove),0,0)
        call shstru('optic',w(osoptc),0,0)
        call shstru('gw',w(osgw),0,0)
        call shstru('str',w(osstr),0,0)
        call shstru('tb',w(ostb),0,0)
        call shstru('spec',w(osspec),1,nspec)
        call shstru('site',w(ossite),1,nsite)
        call rx0('done show')
      endif
      endif

      end

