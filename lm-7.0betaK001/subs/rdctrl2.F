      subroutine rdctrl2(recrd,recln,nrecs,
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
      use m_rdctrl
      use m_gtv
C     implicit none
C     include "mpif.h"
C ... Passed parameters
      integer:: recln,nrecs
      character slabl_(1)*8
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
C ... Local parameters
      character fileid*64
      integer procid,nproc,master
      logical lgors,cmdopt,bittst,ltmp,ioorbp,asa
      double precision dval,dglob,xx(n0*2),dgets,ekap(6)
      integer a2vec,bitand,fopna,getdig,i,is,igets,iprint,
     .  iprt,irs(5),isw,ifi,ix(n0*nkap0),j,k,l,lasa,lbas,lcd,lqp,
     .  lfrzw,lgunit,lham,lmet,lncol,lordn,
     .  lrs,lstsym,lsx1,ltb,lxcf,nat,nlibu,nclasp,
     .  nglob,nkap,nspc,nlmax,scrwid,stdo,stdl,stde,k1,k2,mpipid
        character*(8),allocatable::clabl(:)
        integer,allocatable:: ipc(:),initc(:),ics(:)
        real(8),allocatable:: pnuc(:,:,:),qnuc(:,:,:,:),
     .    pp(:,:,:,:),ves(:),zc(:)
      integer:: dvec1(3)=1, dvec2(3)=0

C ... basis
      double precision orbp(n0,2,nkap0)
      integer o,oclabl,ohave,oics,opnu,opp,oqnu,osgw,osoptc,
     .  osordn,oves,owk
C ... Heap
      integer w(1)
      common /w/ w

c takao 
      real(8):: pnux(20)
      integer:: nnn

      procid = mpipid(1)
      nproc  = mpipid(0)
      master = 0

      scrwid = 80
      stdo = lgunit(1)
      stdl = lgunit(2)
      stde = stdo

C --- Initialize gtv; copy recrd to rcd ---
      call gtv_setst(stdo,stdl,stde)
      call gtv_setrcd(recrd,nrecs,recln)

c$$$      if (pass2) then
c$$$        call upack('ctrl nbas nclass nspec nspin nl',w(osctrl),
c$$$     .    nsite,nclass,nspec,nsp,nl)
c$$$        nbas = nsite
c$$$        call upack('pot opnu oqnu oves',w(ospot),opnu,oqnu,oves,0,0)
c$$$        call upack2('pot opp osoptc',w(ospot),opp,osoptc)
c$$$        call upack('array oics oclabl nclasp ohave',w(osarry),
c$$$     .    oics,oclabl,nclasp,ohave,0)
c$$$        allocate(clabl(nclasp))
c$$$        allocate(lmxa(nclasp),pnuc(nl,nsp,nclasp),qnuc(3,nl,nsp,nclasp),
c$$$     .    pp(6,nl,nsp,nclasp),ves(nclasp),zc(nclasp),
c$$$     .    initc(nclasp),ics(nclasp))
c$$$        call dpcopy(w(opnu),pnuc,1,nl*nsp*nclasp,1d0)
c$$$        call dpcopy(w(oqnu),qnuc,1,3*nl*nsp*nclasp,1d0)
c$$$        call dpcopy(w(opp),pp,1,6*nl*nsp*nclasp,1d0)
c$$$        call dpcopy(w(oves),ves,1,nclasp,1d0)
c$$$        call icopy(nclasp,w(ohave),1,initc,1)
c$$$        call icopy(nclasp,w(oics),1,ics,1)
c$$$        do  j = 1, nclasp
c$$$          xx = dval(w(oclabl),j)
c$$$          call r8tos8(xx,clabl(j))
c$$$          is = w(oics+j-1)
c$$$          call upack('spec lmxa z',w(osspec),is,
c$$$     .      lmxa(j),zc(j),0,0)
c$$$        enddo
c$$$
c$$$
c$$$        call readctrlpq(prgnam,nclasp,nl,nsp,pnuc,qnuc,pp,zc,
c$$$     .    ves,initc,ics,clabl)
c$$$
c$$$        call dpcopy(pnuc,w(opnu),1,nl*nsp*nclasp,1d0)
c$$$        call dpcopy(qnuc,w(oqnu),1,3*nl*nsp*nclasp,1d0)
c$$$        call dpcopy(pp,w(opp),1,6*nl*nsp*nclasp,1d0)
c$$$        call dpcopy(ves,w(oves),1,nclasp,1d0)
c$$$        call icopy(nclasp,initc,1,w(ohave),1)
c$$$        deallocate(lmxa,pnuc,qnuc,pp,ves,zc,initc,ics,clabl)
c$$$        call upack('ctrl lasa',w(osctrl),lasa,0,0,0,0)
c$$$        lasa = lasa + lasa3 + 8*isw(lasa8)
c$$$        call pack5('ctrl lasa',w(osctrl),lasa,0,0,0,0)
c$$$        return
c$$$      endif

C --- Read input parameters from contents of rcd ---
      call readctrl(prgnam,vrsion(1),vn(1))
      if (io_help > 0) then
        call readctrlpq(prgnam,nclasp,nl,nsp,pnuc,qnuc,pp,zc,
     .    ves,initc,ics,clabl)
        call cexit(0,1)
      endif

      call defrr(owk,1)
C     Set switches depending type of program
      lbas = 0
      lcd = 0
C     For now, LMF => fp; no screening; nfp-style hamiltonian
C     cd represented in plane waves
      if (lfp .ne. 0) then
        lbas = 3
      endif
      if (trim(prgnam) .eq. 'LMMC') then
        lbas = 1
      endif
      asa = .false.
      if (prgnam .eq. 'LM' .or. prgnam .eq. 'LMGF' .or.
     .    prgnam .eq. 'LMPG' .or. prgnam .eq. 'LMCTL') then
        asa = .true.
      endif

C ... Optionally read positions from pos file
      if (cmdopt('--rpos=',7,0,fileid)) then
        call iopos(.false.,-1,fileid(8:),nbasp,pos)
      endif

C ------------------- Copy to structures ----------------------
C     (query-replace-regexp "\\([a-z_]+\\)\\([0-9]+\\)" "\\2*isw(\\1\\2)" nil)

C ... Initialize sstrn
      call ustrn(w,0,len(sstrn),0,0,0)
      sstrn = ' '

C --- Allocate and copy input to sbz ---
      call defdr(osbz,1000)
      call ubz(w(osbz),-1,-000,0,0,0)
      call redfrr(osbz,nint(dval(w(osbz),1)))
      call pack5('bz def dosw ef efmax fsmom',w(osbz),
     .  bz_def,bz_dosw,bz_ef,bz_efmax,bz_fsmom)
      i = 1*isw(bz_lio1)+2*isw(bz_lio2)+8*isw(bz_lio8)
      call pack5('bz lcond lio lmet lmull lopt',w(osbz),
     .  bz_lcond,i,bz_lmet,bz_lmull,0)
      call pack5('bz lshft n ndos nevmx nkabc',w(osbz),
     .  bz_lshft,bz_n,bz_ndos,bz_nevmx,bz_nabc)
      call pack5('bz range semsh w zval',w(osbz),
     .  bz_range,bz_semsh,bz_w,bz_zval,0)

C --- Allocate and copy input to sctrl ---
      call defdr(osctrl,1000)
      call uctrl(w(osctrl),-1,0,0,0,0)
      call redfrr(osctrl,nint(dval(w(osctrl),1)))
      call uctrl(w(osctrl),-1,0,0,0,0)

C     lasa: 1 Make V from P,Q  2 Make pp  4 ccor  8 free atm
C          16 map  32 nonspherical mpol moms 64 MT corr
C         128 interpretation of sphere Q2; see newrho.f
C         256 how ASA Q1,Q2 are accumulated; see makwts.f
C         512 (spin pol) alpha repsn = (gamma(1) + gamma(nsp))/2
      lasa=4*isw(lasa4)+32*isw(lasa32)+64*isw(lasa64)+128*ham_qasa
     .    +512*isw(lasa512)
      if (.not. asa) lasa=0
C     lbas: 1 Hamiltonian has no screening transformation
C           2 Hamiltonian is nfp style
C          16 freeze phi,phidot for all species
      j = lbas + 16*isw(frzwf)
C     lcd: 1 freeze core
C          2 non-self-consistent Harris
C          4 represent full potential density on a uniform mesh
C          8 represent full potential density via TCF
C         16 unused
C         32 unused
C         64 (molecules) XC potential by FFT
      k = 1*isw(lcd1)+2*isw(lcd2)+4*isw(lcd4)+8*isw(lcd8)+64*isw(lcd64)
      call pack5('ctrl defm elin lasa lbas lcd',w(osctrl),
     .  lat_defm,asa_elin,lasa,j,k)
      k = isw(nmto.gt.1)
      if (nmto.gt.1) k = k+2*isw(ham_ewald)
      call pack5('ctrl lcgf ldos lfp lfrce lgen3',w(osctrl),
     . ctrl_lcgf,ctrl_ldos,lfp,ctrl_lfrce,k)
C     lham  1 (ASA) 2-center
C           1 (molecules) two-panel
C           2 (ASA) 2-c + pert. corr
C           4 (ASA) auto-down-fold
C           8 (ASA) change rep interactively
C          16 (ASA) suppress d hybridization
C          32 (ASA) preserve ortho. evecs
C          64 (ASA) save evecs to disk
C         128 (ASA) gamma-rep
C         256       use true spherical harmonics
      lham = 1*isw(lham1)+4*isw(lham4)+8*isw(lham8)+16*isw(lham16)+
     .      32*isw(lham32)+64*isw(lham64)+128*isw(lham128)+
     .     256*isw(lham256) + lham3
C     ctrl_lmet    1 metal  2 tetrahedron
Cr                 4 (PGF) V-shift1 is zero
Cr                 8 (PGF) V-shift2 is zero
      lmet = isw(bz_lmet.ne.0) + 2*isw(ctrl_lmet2) +
     .     4*isw(ctrl_lmet4) + 8*isw(ctrl_lmet8)
C     lncol 1 noncollinear magnetism
C           2 spin spirals
C           4 spin-orbit coupling
C           8 External magnetic field
C          16 mag. forces
C          32 spin-orbit coupling, LzSz only
C          64 spin-orbit coupling, LzSz + (L.S-LzSz) pert
C     If spin-orbit or SS, also turn on noncollinear
      k = 2*isw(lncol2)+4*isw(lncol4)+8*isw(lncol8)+16*isw(lncol16)
      if (k .ne. 0) lncol1=T
      lncol = 1*isw(lncol1)+2*isw(lncol2)+4*isw(lncol4)+8*isw(lncol8)+
     .  16*isw(lncol16)+32*isw(lncol32)+64*isw(lncol64)
      lordn = 0
      call pack5('ctrl lham lmet lncol loptc lordn',w(osctrl),
     .  lham,lmet,lncol,ctrl_loptc,lordn)
C     lrs  switches concerning restart mode.
C         1 Read from restart file
C         2 Read from restart file, ascii mode
C         4 Read from restart file, invoke smshft
C         8 Write new density to restart file
C        16 Write new density to restart file, ascii format
C        32 read site positions from input file
C        64 read starting fermi level from input file
C       128 read starting pnu level from input file
C       256 rotate local density after reading
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
C     lqp 1 do not add inversion 2 inverse iteration
      lqp = 1*isw(ctrl_lqp1)+2*isw(ctrl_lqp2)
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
       call pack5('ctrl lpgf lqp lrel lrs lscr',w(osctrl),
     .  ctrl_lpgf,lqp,lrel,lrs,lscr)
       
C      lstr  no longer used.
C      ltb switches for empirical tight-binding
C         1 overlap        2 crystal-field     4 ovlp+CF
C         8 add ebarLL    16 forces           32 fij
C        64 not used     128 pressure        256 evdisc
C       512 pair pot    1024 TrH & local E  2048 local rho
C      2^12 Hubbard U   2^13 No Madelung    2^14 wgt avg U
C      2^15 L>0 estat   2^16 disc read incr 2^17 gamma-pt
       i = 1*isw(ltb1)+2*isw(ltb2)+4*isw(ltb4)+8*isw(ltb8)+16*isw(ltb16)
     .   +32*isw(ltb32)+64*isw(ltb64)+128*isw(ltb128)+256*isw(ltb256)
     .   +512*isw(ltb512)+1024*isw(ltb1024)+2048*isw(ltb2048)
     .   +4096*isw(ltb4096)+2**13*isw(ltb213)+2**14*isw(ltb214)
     .   +2**15*isw(ltb215)+2**16*isw(ltb216)+2**17*isw(ltb217)
     .   +2**18*isw(ltb218)
C      lves    1 take ves as input
       call pack5('ctrl lstonr lstr lsx ltb lves',w(osctrl),
     .   lstonr,0,lsx,i,lves)
C      lxcf   parameter defining XC functional
C      1s digit:
C      1 for Ceperly-Alder
C      2 for Barth-Hedin (ASW fit)
C      3 for PW91
C      4 for PBE
C      10s digit for GGAs
C      0 for LSDA
C      1 for LMH
C      2 for PW91
C      3 for PBE
C      4 for PBE with Becke exchange
       i = ham_lxcf + 10*gga
C      Set modep
       ix(1:3) = 2
       if (prgnam == 'LMPG') ix(3) = 0
       if (prgnam == 'LMMC') ix(1:3) = 0
       call pack5('ctrl lxcf maxit mdprm modep nbas',w(osctrl),
     .   i,iter_maxit,mdprm,ix(1:3),nbas)
C      Reset nl
       if (nl .ne. max(lmxbx,lmxax)+1 .and. io_help == 0) then
         call info2(20,1,0,' rdctrl: reset global max nl from %i to %i',
     .     nl,max(lmxbx,lmxax)+1)
         nl = max(lmxbx,lmxax)+1
       endif
C      call pack5('ctrl nbasp nclass nesabc nitmv nl',w(osctrl),
       call pack5('ctrl nbasp nesabc nitmv nl',w(osctrl),
     .   nbasp,nesabc,nitmv,nl,0)
C      call pack5('ctrl nlibu nmap npl nsite nspec',w(osctrl),
       call pack5('ctrl nsite nspec',w(osctrl),
     .   nsite,nspec,0,0,0)
C      call shosyv(0,0,0,6)
       call pack5('ctrl nspin nvario omax1 omax2 quit',w(osctrl),
     .   nsp,nvario,omax1,omax2,quit)
       call pack5('ctrl rmaxes rmines sclwsr sdmod sdprm',w(osctrl),
     .   rmaxes,rmines,sclwsr,sdmod,sdprm)
       call pack5('ctrl smalit tol wsrmax zbak pfloat',w(osctrl),
     .   smalit,ctrl_tol,wsrmax,zbak,lpfloat)

C --- Allocate and copy input to slat ---
      call defdr(oslat,1000)
      call ulat(w(oslat),-1,0,0,0,0)
      call redfrr(oslat,nint(dval(w(oslat),1)))
      if (dalat == NULLR) dalat=0
      call pack5('lat alat as avw nkdmx nkqmx',w(oslat),
     .  alat+dalat,lat_as,avw,lat_nkdmx,lat_nkdmx)
      call pack5('lat gam gmax nabc',w(oslat),
     .  lat_gam,lat_gmax,ftmesh,0,0)
      call pack5('lat gam ldist dist',w(oslat),
     .  lat_gam,lat_ldist,lat_dist,0,0)
      call pack5('lat plat plat2',w(oslat),plat,slat_plat2,0,0,0)
      call pack5('lat rpad slat tol tolft vol',w(oslat),
     .  lat_rpad,lat_slat,lat_tol,tolft,vol)
      call dcopy(6,plat,1,xx,1)
      call dcopy(3,platl,1,xx(7),1)
      call pack5('lat platl',w(oslat),xx,0,0,0,0)
      call dcopy(3,platr,1,xx(7),1)
      call pack5('lat platr',w(oslat),xx,0,0,0,0)

C --- Allocate and copy input to sham ---
      call defdr(osham,1000)
      call uham(w(osham),-1,-000,0,0,0)
      call redfrr(osham,nint(dval(w(osham),1)))
      call pack5('ham alfsi dabc elind nmto kmto',w(osham),
     .  alfsi,dabc,elind,nmto,kmto)
C     Replicate ctrl->lncol in ham->lncol, ditto for lham,lgen3
      call pack1('ham lncol',w(osham),igets('ctrl lncol',w(osctrl)))
      call pack1('ham lham',w(osham),igets('ctrl lham',w(osctrl)))
      call pack1('ham lgen3',w(osham),igets('ctrl lgen3',w(osctrl)))
C     Mop up integer parts of sigp
      sigp(1) = sigp_mode
      sigp(2) = sigp_nmin
      sigp(4) = sigp_nmax
      call pack5('ham nkaph pmax pmin lsig sigp',w(osham),
     .  nkaph,pmax,pmin,lrsig,sigp)
      call pack5('ham qss rsrnge rsstol udiag',w(osham),
     .  ham_qss,rsrnge,rsstol,ham_udiag,0)
C     Parameters for APW
      call pack5('ham pwmode npwpad pwemin pwemax',w(osham),pwmode,
     .  npwpad,pwemin,pwemax,0)
      call pack1('ham oveps',w(osham),oveps)

C ... Allocate and initialize sarray
      call defrr(osarry,1000)
      call uarray(w(osarry),-1,NULLI,0,0,0)
      call redfrr(osarry,nint(dval(w(osarry),1)))

C ... Allocate and initialize spot
      call defdr(ospot,1000)
      call upot(w(ospot),-1,-000,0,0,0)
      call redfrr(ospot,nint(dval(w(ospot),1)))
      call pack5('pot vmtz0',w(ospot),vmtz,0,0,0,0)

C --- Allocate and copy input to smix ---
      call defrr(osmix,1000)
      call umix(w(osmix),-1,-000,0,0,0)
      call redfrr(osmix,nint(dval(w(osmix),1)))
ckino's fix for SR11000 size() return's integer(8), right?
      nnn = size(smix)
      call dcopy(nnn,smix,1,w(osmix),1)

C ... Allocate and initialize sordn
      call defrr(osordn,1000)
      call uordn(w(osordn),-1,0,0,0,0)
      call redfrr(osordn,nint(dval(w(osordn),1)))
      call pack1('array osordn',w(osarry),osordn)

C --- Allocate and copy input to smove ---
      call defrr(osmove,1000)
      call umove(w(osmove),-1,0,0,0,0)
      call redfrr(osmove,nint(dval(w(osmove),1)))
      if (lbsprm) then          !Load Bulirsch-Stoer parameters into structure
        prmint(2) = isw(prmint_new)
        prmint(3) = prmint_ts0
        prmint(4) = prmint_tol
        prmint(5) = prmint_mx
        prmint(6) = prmint_mi
        prmint(7:6+prmint_mi) = prmint_nseq(1:prmint_mi)
C        print 331, prmint(1:17)
C  331   format(17f8.4)
      endif
      call pack5('move gyro nmodt modt ct prmint',w(osmove),
     .  2d0,gd_nmodt,gd_modt,gd_ct,prmint)
      call pack5('move kt ts tsequ tstot',w(osmove),
     .  move_kt,move_ts,move_tsequ,move_tstot,0)

C --- Allocate and copy input to soptic ---
      call defrr(osoptc,1000)
      call uoptic(w(osoptc),-1,0,0,0,0)
      call redfrr(osoptc,nint(dval(w(osoptc),1)))
      call pack1('pot osoptc',w(ospot),osoptc)
      call pack5('optic ne window ocrng unrng esciss',w(osoptc),
     .  optic_ne,optic_window,optic_ocrng,optic_unrng,optic_esciss)
      call pack2('optic mode nchi2',w(osoptc),isw(optic_mode1),
     .  optic_nchi2)
      call pack1('optic axes',w(osoptc),optic_axes) ! nchi2 must be
                                                    ! packed beforehand
C --- Allocate and copy input to sgw ---
      call defrr(osgw,1000)
      call ugw(w(osgw),-1,0,0,0,0)
      call redfrr(osgw,nint(dval(w(osgw),1)))
C     Old defaults
C     call pack5('gw gcutb gcutx qoffp nband gsmear',w(osgw),2.7d0,
C    .  2.2d0,1d0,9999,.003d0)
CC    Default : Faleev's fast energy integration
C     xx(1) = .01d0
C     xx(2) = .04d0
C     call pack5('gw lgw delre ecuts',w(osgw),1,xx,2.5d0,0,0)
C     call pack5('gw nime delre deltax deltaw pbtol',w(osgw),6,xx,
C    .  -1d-4,.02d0,1d-3)
C     call spacks(1,'gw pb1',w(osgw),'111',0,0)
C     call spacks(1,'gw pb2',w(osgw),'1111',0,0)
      call pack5('gw gcutb gcutx qoffp nband gsmear',w(osgw),
     .  gw_gcutb,gw_gcutx,gw_qoffp,gw_nband,gw_gsmear)
      call pack5('gw lgw nkabc delre ecuts',w(osgw),
     .  1,gw_nabc,gw_delre,gw_ecuts,0)
      call pack5('gw nime deltax deltaw pbtol',w(osgw),
     .  gw_nime,gw_deltax,gw_deltaw,gw_pbtol,0)
      call spacks(1,'gw pb1',w(osgw),'111',0,0)
      call spacks(1,'gw pb2',w(osgw),'1111',0,0)
      call pack1('pot osgw',w(ospot),osgw)

C --- Allocate and copy input to sstr ---
      call defrr(osstr,1000)
      call ustr(w(osstr),-1,0,0,0,0)
      call redfrr(osstr,nint(dval(w(osstr),1)))
C     Old defaults
C     xx(1) = 0
C     xx(2) = -1
C     xx(3) = -2.3d0
C     if (prgnam .eq. 'TBE') then
C       call pack5('str nkaps rmax rfit kaps lmaxw',w(osstr),1,
C    .    0d0,.8d0,xx,-1)
C     else
C       call pack5('str nkaps rmax rfit kaps lmaxw',w(osstr),1,
C    .    0d0,.8d0,xx,-1)
C     endif
C     call pack2('str loka drwats',w(osstr),1,.1d0)
C    loka switch
      i = 1
      if (str_mode .eq. 1 .or. str_mode .eq. 3) i = 0
C     Pack iinv parameters
      call dpzero(xx,5)
      xx(1) = iinv_nit 
      xx(2) = iinv_ncut
      xx(3) = iinv_tol
      call pack5('str adec amode drwats iinv',w(osstr),
     .  tcf_adec,str_mode,str_drwats,xx,0)
      call pack5('str lequiv loka lmaxw lshow mxnbr',w(osstr),
     .  isw(str_lequiv1),i,str_lmaxw,isw(str_lshow1),str_mxnbr)
      call pack5('str nalf nbisi ncupl ndust wztcf',w(osstr),
     .  tcf_nalf,tcf_nbisi,tcf_ncupl,tcf_ndust,tcf_wztcf)
      call pack5('str nkaps rmax rfit kaps',w(osstr),
     .  str_nkaps,str_rmax,0.8d0,str_kaps,0)
      call pack5('str rmaxg ivl',w(osstr),str_rmaxg,str_ivl,0,0,0)

C --- Allocate and copy input to stb ---
      call defrr(ostb,1000)
      call utb(w(ostb),-1,NULLI,0,0,0)
      call redfrr(ostb,nint(dval(w(ostb),1)))
C     old defaults
C      call pack1('tb alam',w(ostb),.001d0)
C      call pack1('tb alsc',w(ostb),10d0)
      call pack5('tb alam alsc',w(ostb),
     .  tb_alam,tb_alsc,0,0,0)
      call pack5('tb ebfit fmode nbfit rmfit',w(ostb),
     .  tb_ebfit,tb_fmode,tb_nbfit,tb_rmfit,0)


C --- Allocate and copy input to sspec ---
      call defrr(osspec,1000)
      call uspec(w(osspec),-1,0,1,0,0,0)
      call redfrr(osspec,nspec*nint(dval(w(osspec),1)))
      call uspec(w(osspec),-1,-000,nspec,0,0,0)
C ... old defaults
C      call spackv(1,'spec rs3',w(osspec),1,nspec,1d0)
C      call spackv(1,'spec vmtz',w(osspec),1,nspec,-0.5d0)
CC ... Init default for sigma,lmxa,lmxf
C      call dvset(xx,1,10,-1d0)
CC ... First set lmxa to 4 to set def. hcr to l=4
C      call spackv(1,'spec lmxa',w(osspec),1,nspec,4)
C      call spackv(1,'spec lmxpb',w(osspec),1,nspec,3)
C      call spackv(1,'spec hcr',w(osspec),1,nspec,xx)
C      call spackv(1,'spec lmxa',w(osspec),1,nspec,nl-1)
C      call spackv(1,'spec lmxf',w(osspec),1,nspec,2*nl-2)
C      xx(1) = -1
C      xx(2) = 0
C      call spackv(1,'spec coreq',w(osspec),1,nspec,xx)
C      if (prgnam(1:2) .eq. 'TB')
C     .  call spackv(1,'spec rmt',w(osspec),1,nspec,1d0)
C      do  j = 1, nspec
CC       call spacks(1,'spec pb1',w(osspec),'111     ',j,j)
CC       call spacks(1,'spec pb2',w(osspec),'1111    ',j,j)
C        call spacks(1,'spec pb1',w(osspec),'        ',j,j)
C        call spacks(1,'spec pb2',w(osspec),'        ',j,j)
C        call spacks(1,'spec coreh',w(osspec),'        ',j,j)
C      enddo
      call spackv(1,'spec eh3',w(osspec),1,nspec,-0.5d0)
      call spackv(1,'spec etf',w(osspec),1,nspec,-1d0)
      call spackv(1,'spec lmxf',w(osspec),1,nspec,2*nl-2)
      call spackv(1,'spec norp',w(osspec),1,nspec,2)
      call spackv(1,'spec vmtz',w(osspec),1,nspec,-0.5d0)

      do  j = 1, nspec
        slabl_(j) = slabl(j)
        call spacks(1,'spec pb1',w(osspec),pb1(j),j,j)
        call spacks(1,'spec pb2',w(osspec),pb2(j),j,j)
        call spacks(1,'spec coreh',w(osspec),coreh(j),j,j)
        call spacks(1,'spec name',w(osspec),slabl(j),j,j)
        call pack5('spec a nr alpha coreq ',w(osspec),j,
     .    spec_a(j),nr(j),alpha(1,j),coreq(1,j))
        call pack5('spec lxi nxi exi',w(osspec),j,
     .    lxi(j),nxi(j),exi(1,j),0)
        call pack5('spec group grp2 hcr idmod',w(osspec),j,
     .    grp(j),grp2(j),hcr(1,j),idmod(1,j))
        call pack5('spec ehvl',w(osspec),j,ehvl(1,j),0,0,0)
C       Set idxdn
        call ivset(ix,1,n0*nkap0,1)
        call icopy(1+lmxb(j),idxdn(1,j),1,ix,1)
        call pack5('spec idxdn idu jh uh',w(osspec),j,
     .    ix,idu(1,j),jh(1,j),uh(1,j))
        call pack5('spec kmxt kmxv lfoca rsmv',w(osspec),j,
     .    kmxt(j),kmxv(j),lfoca(j),rsmv(j))
        call pack5('spec lmxa lmxb lmxl lmxpb',w(osspec),j,
     .    lmxa(j),lmxb(j),lmxl(j),lmxpb(j))
C       pack mxcst(j)        
        i = 1*isw(mxcst1(j))+2*isw(mxcst2(j))+4*isw(mxcst4(j))
C       pack orbp(j)
        call dpzero(orbp,n0*2*nkap0)
        call dcopy(n0,rsmh(1,j),1,orbp(1,1,1),1)
        call dcopy(n0,eh(1,j),1,orbp(1,2,1),1)
        call dcopy(n0,rsmh2(1,j),1,orbp(1,1,2),1)
        call dcopy(n0,eh2(1,j),1,orbp(1,2,2),1)
        call pack5('spec mass mxcst orbp',w(osspec),j,mass(j),i,orbp,0)
C       Pack P,Q,PZ for both spins        
        call dpzero(orbp,n0*6)
        call dcopy(n0*nsp,pnu(1,1,j),1,orbp(1,1,1),1)
        call dcopy(n0*nsp,pz(1,1,j),1,orbp(1,1,2),1)
        call dcopy(n0*nsp,qnu(1,1,j),1,orbp(1,1,3),1)
        call pack5('spec p pz q z',w(osspec),j,orbp,
     .    orbp(1,1,2),orbp(1,1,3),z(j))
        call pack5('spec colxbs radxbs',w(osspec),j,
     .    colxbs(1,j),radxbs(j),0,0)
        call pack5('spec rcfa rcut rfoca rg rham',w(osspec),j,
     .    rcfa(1,j),rcut(j),rfoca(j),rg(j))
        call pack5('spec rmt rs3 rsma rsmfa',w(osspec),j,
     .    rmt(j),rs3(j),rsma(j),rsmfa(j))
        if (ltbe) then
          call pack5('spec iq1 ivso stni vso',w(osspec),j,
     .      iq1(1,j),ivso(1,j),stni(j),tbvso(1,j))
          call pack5('spec qpol',w(osspec),j,qpol(1,j),0,0,0)
          if (ltb217 .or. ltb218) then
           call pack1('bz nkabc',w(osbz),dvec1)
           call pack1('bz lshft',w(osbz),dvec2)
          endif
        endif
        call pack5('spec dv eref rham rint',w(osspec),j,
     .    dv(j),eref(j),rham(j),rint(j))
      enddo
C     call shstru('spec',w(osspec),1,nspec);stop
cccccccccccccc
c       call upack('spec p ',w(osspec),1,pnux,0,0,0)
c       print *,'xxx:111 222 ini pnu=',pnu
cccccccccccccc


C --- Allocate and copy input to ssite ---
C     NB: essential that ssite is LAST array allocated here
      call defrr(ossite,1000)
      call usite(w(ossite),-1,0,1,0,0,0)
      call redfrr(ossite,nsite*nint(dval(w(ossite),1)))
      call usite(w(ossite),-1,0,nsite,0,0,0)
C     Old defaults
C     call ivset(ix,1,n0*nkap0,1)
C     call spackv(1,'site relax',w(ossite),1,nsite,ix)
      do  j = 1, nsite
C       Pack spec into class in case spec info is sufficient
C       and spid into class label
        call pack5('site spec class dpole mpole',w(ossite),j,
     .    ips(j),ips(j),dpole(1,j),mpole(j))
        call spacks(1,'site clabel',w(ossite),slabl(ips(j)),j,j)
        call pack5('site pl plv pos vel',w(ossite),j,
     .    ipl(j),plv(j),pos(1,j),vel(1,j))
        call pack5('site vshft relax eula',w(ossite),j,
     .    vshft(j),irlx(1,j),eula(1,j),0)
        if (ltbe) then
          call pack2('site ndelta',w(ossite),j,ndelta(j))
          call pack2('site delta',w(ossite),j,delta(1,j))
        endif
      enddo

C ... Copy string outputs to sstrn : amix, gfopt, jobid, mix, mmham, sxopt, symg
      j = len_trim(iter_amix)          ! Euler angle mixing amix
      if (j .gt. 0) then
        call lstra('strn amix',i,o,k)  ! o = index amix has (see ustrn),
        call ustrn(w,-o,1,i,k,j)       ! Copy to sstrn(i:i+j-1)
        sstrn(i:i+j-1) = iter_amix
      endif
      j = len_trim(iter_mix)           ! density mixing mix
      if (j .gt. 0) then
        call lstra('strn mix',i,o,k)   ! o = index mix has (see ustrn),
        call ustrn(w,-o,1,i,k,j)       ! Copy to sstrn(i:i+j-1)
        sstrn(i:i+j-1) = iter_mix
      endif
      j = len_trim(mmham)              ! micromagnetics hamiltonian
      if (j .gt. 0) then
        call lstra('strn mmham',i,o,k) ! o = index mmham has (see ustrn),
        call ustrn(w,-o,1,i,k,j)       ! string should be copied to sstrn(i:j)
        sstrn(i:i+j-1) = mmham
      endif
      j = len_trim(header)             ! jobid
      if (j .gt. 0) then
        call lstra('strn jobid',i,o,k) ! o = index jobid has (see ustrn),
        call ustrn(w,-o,1,i,k,j)       ! Copy to sstrn(i:i+j-1)
        sstrn(i:i+j-1) = header
      endif
      j = len_trim(gfopt)              ! Green's function gfopt
      if (j .gt. 0) then
        call lstra('strn gfopt',i,o,k) ! o = index mix has (see ustrn),
        call ustrn(w,-o,1,i,k,j)       ! Copy to sstrn(i:i+j-1)
        sstrn(i:i+j-1) = gfopt
      endif
      j = len_trim(sxopt)              ! Screened exchange sxopt
      if (j .gt. 0) then
        call lstra('strn sxopt',i,o,k) ! o = index mix has (see ustrn),
        call ustrn(w,-o,1,i,k,j)       ! Copy to sstrn(i:i+j-1)
        sstrn(i:i+j-1) = sxopt
      endif

C ... Suppress symmetry operations for special circumstances
      lstsym = 0                
      if (lncol .ne. 0) lstsym=1 !lstsym=1: noncollinear case
                                 !      =2: turn off symops
C     Switches that automatically turn of all symops 
      if ((mdprm(1) .ge. 1 .and. mdprm(1) .le. 3) .or.   
     .  cmdopt('--cls',5,0,strn) .or. cmdopt('--nosym',7,0,strn)) then
        symg = 'e'
        lstsym = 2               !lstsym=2: turn off symops
      endif
C     Switches that turn off automatic finder, incl. inversion
      if (lstsym .ne. 0) then
        i = 1
        do while (i .ne. 0)
          i = index(symg,'find')
          if (i .ne. 0) then
            symg(i:i+3) = ' '
          endif
        enddo
        if (symg .eq. ' ') symg = 'e'  ! suppress ops if none explicit
        lqp = lqp-bitand(lqp,1)+1
      endif
      j = len_trim(symg)               ! Symmetry group symg
      if (j .gt. 0) then
        call lstra('strn symg',i,o,k)  ! o = index symg has (see ustrn),
        call ustrn(w,-o,1,i,k,j)       ! Copy to sstrn(i:i+j-1)
        sstrn(i:i+j-1) = symg
      endif
C ... End of copy

C --- Miscellaneous cleanup and initialization ---
      if (lgors('ctrl lgen3',w(osctrl))) then
        nkap = dgets('str skmsh:1.1',w(osstr))
        xx(1) = dgets('str skmsh:2.2',w(osstr))
        xx(2) = dgets('str skmsh:3.3',w(osstr))
        call gausq(nkap,xx(1),xx(2),ekap,xx(3),0,0)
      endif

C ... Read the basis from the atm file
      if (cmdopt('--rdbasp',8,0,strn)) then
        fileid = 'basp'
        if (strn(9:12) .eq. ':fn=') then
          fileid = strn(13:13+63) !tk to avoid complaint; replace fileid = strn(13:)
        else
        endif

C       Number of envelope functions before modifications
        call uspecb(0,1,w(osspec),1,nspec,xx,xx,xx,k1)

C       Read new parms; return in k2 max no. env. fns/site read
        call strip(fileid,i,j)
        ifi = fopna(fileid(1:j),-1,0)
        rewind ifi
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

C     Add dalat to alat
      call pack1('lat alat',w(oslat),dgets('lat alat',w(oslat))+dalat)
C     Replicate ctrl->lncol in ham->lncol
      call pack1('ham lncol',w(osham),lncol)
      lsx1 = mod(lsx,2)
      if (lasa32) call pack1('mix nsave',w(osmix),3)
C ... Use true spherical harmonics
      if (bittst(lncol,4) .or. ctrl_loptc .ne. 0 .or. lrel .eq. 2)
     .  call lsets('ctrl lham',w(osctrl),.true.,256)

C ... Dirac equation requires spin polarization
      if (nsp.eq.1 .and. igets('ctrl lrel',w(osctrl)).eq.2) then
        call rx('rdccat: Dirac equation requires NSPIN=2')
      endif

C ... Suppress inversion when noncollinear magnetism, SX, NLO
      if (lncol+lsx1 .ne. 0 .or. ctrl_loptc .ge. 10 .or.
     .   (mod(lscr,10) .eq. 1 .and. prgnam .eq. 'LM'))
     .  lqp = lqp-bitand(lqp,1)+1
      call pack5('ctrl lqp',w(osctrl),lqp,0,0,0,0) !repack lqp

C ... Special pgf initialization
      if (ctrl_lpgf(1) .ne. 0) then
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
C     Make nat = number of real atoms as nbas - # sites w/ floating orbitals
      if (procid .eq. master) then
      nat = nbas
      do  i = 1, nbas
        call upack('site spec',w(ossite),i,j,0,0,0)
        call upack('spec lmxa',w(osspec),j,l,0,0,0)
        if (l .eq. -1) nat = nat-1
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

C ... Count LDA+U blocks (printout only)
      call defi(owk,-nbas)
      call pshpr(0)
      call suldau(nbas,w(osspec),w(ossite),nlibu,k,w(owk))
      call pack5('ham nlibu lmaxu',w(osham),nlibu,k,0,0,0)
      call poppr
      call rlse(owk)

C     Free arrays used to read input
      deallocate(pnu,qnu,pz,amom,idmod,rsmh,eh,rsmh2,eh2,pb1,pb2,
     .  lmxpb,qpol,stni,tbvso,iq1,ivso,rg,rsma,rfoca,rsmfa,rcfa,nxi,
     .  exi,rint,rcut,coreq,mass,colxbs,radxbs,rs3,rham,idxdn,hcr,
     .  rmt,alpha,idu,uh,jh,dv,grp,grp2,mxcst1,mxcst2,mxcst4,kmxt,kmxv,
     .  lfoca,eref,lmxl,lxi,coreh,lmxa,lmxb,spec_a,z,nr,rsmv)
      deallocate(pos,vel,eula,vshft,ips,ipl,plv,irlx,mpole,dpole)
      if (ltbe) deallocate(delta,ndelta)

C      call rx0('done copying input to structures')

C --- Printout ---
      if (iprint() .ge. 20 .and. procid .eq. master) then
        do  90  k = 1, 2

        strn = '  '//prgnam
        if (prgnam .eq. 'LMMC') then
          call awrit4(
     .    '%N %a:%12palat = %;5d  '//
     .    'nbas = %i%?#n#%-1j(+%i)##'//
     .    '  nspec = %i',strn,scrwid,0,dgets('lat alat',w(oslat)),
     .      nbas,nsite-nbas,nspec)
        else
          call awrit6(
     .    '%N %a:%12palat = %;5d  '//
     .    'nbas = %i%?#n#%-1j(+%i)##%?#n#+%i#%j#'//
     .    '  nspec = %i',strn,scrwid,0,dgets('lat alat',w(oslat)),
     .    nat,nbas-nat,nbas-nsite,nsite,nspec)
        endif
        if (vn(2) .eq. 0) then
          call awrit6('%a  vn %,2d  verb %i%?#n#,%i#%j#%?#n#,%i#%j#',
     .    strn,scrwid,-lgunit(k),vn(1)-vn2(1)/1d5,
     .    iprint(),iprint()-iprt(1),iprt(1),iprt(2)-iprt(1),iprt(2))
        else
          call awrit7('%a  vn %,2d('//trim(prgnam)//' %,1d)  '//
     .      'verb %i%?#n#,%i#%j#%?#n#,%i#%j#',
     .    strn,scrwid,-lgunit(k),vn(1)-vn2(1)/1d5,vn(2),
     .    iprint(),iprint()-iprt(1),iprt(1),iprt(2)-iprt(1),iprt(2))
        endif

        call upack('ctrl lxcf',w(osctrl),lxcf,0,0,0,0)
        i = igets('ctrl lham,256',w(osctrl))
        lfrzw = isw(lgors('ctrl lbas,16',w(osctrl)))

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
     .    '%?;n; oveps,;;'//
     .    '%b %b',strn,scrwid,0,ctrl_lfrce,lrel,i,ctrl_loptc,
     .    mod(lscr,10),lordn,pwmode,isw(oveps.ne.0))
        if (strn .ne. ' special:') call awrit0(strn,' ',-80,lgunit(k))

        ltmp = prgnam.eq.'LM'  .or. prgnam.eq.'LMMC'  .or.
     .         prgnam.eq.'LMF' .or. prgnam.eq.'LMFGWD' .or.
     .         prgnam.eq.'LMGF' .or. prgnam.eq.'LMPG' .or.
     .         prgnam.eq.'LMFA' .or. prgnam.eq.'LMCHK'
        if (ltmp) then
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
     .      strn,scrwid,0,nsp,nlibu,mod(lxcf,10),lxcf/10,lfrzw,
     .      isw(lasa32),lrsig,lsx)
          call awrit0(strn,' ',-80,lgunit(k))
        endif
        call isanrg(mod(lxcf,10),0,4,prgnam,'XC functional',.true.)
        call isanrg(lxcf/10,0,4,prgnam,'GGA functional',.true.)
        if (mod(lxcf,10) .eq. 2 .and. lxcf/10 .gt. 1 .or.
     .      mod(lxcf,10) .ne. 2 .and. lxcf/10 .eq. 1) call info0
     .    (10,0,0,'%10f(warning) mixing incompatible functionals')

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

        if (asa) then
          call upack('ctrl lham lcd',w(osctrl),lham,lcd,0,0,0)
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
     .      bitand(lasa,4).eq.0,
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

        if (prgnam(1:3) .eq. 'LMF') then
          call awrit1(' special:%10p'//
     .      '%?;n==1; core-level-optics,;;'//
     .      '%b %b',strn,scrwid,0,igets('optic cls',w(osoptc)))
          if (strn .ne. ' special:') call awrit0(strn,' ',-80,lgunit(k))
        endif

c       ltmp = prgnam.eq.'LM'   .or. prgnam.eq.'LMMC'  .or.
        ltmp = prgnam.eq.'LM'   .or. 
     .         prgnam.eq.'LMF'  .or. prgnam.eq.'LMFGWD' .or.
     .         prgnam.eq.'LMGF' .or. prgnam.eq.'LMPG' .or. 
     .         prgnam.eq.'LMDOS'.or. prgnam.eq.'TBE'

        if (ltmp) call awrit7(' bz:%10p'//
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

C --- Debugging printout ---
      if (io_help == 0 .and. io_show > 1) then
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
C       call shstru('pot',w(ospot),0,0)
        call shstru('ham',w(osham),0,0)
        if (smix(2) /= NULLI)
     .  call shstru('mix',w(osmix),0,0)
        call shstru('move',w(osmove),0,0)
        call shstru('optic',w(osoptc),0,0)
        call shstru('gw',w(osgw),0,0)
        call shstru('str',w(osstr),0,0)
        call shstru('tb',w(ostb),0,0)
        call shstru('spec',w(osspec),1,nspec)
        call shstru('site',w(ossite),1,nsite)
        call rx0('done show')
      endif
C     call rx0('done')
      end

C      subroutine string_replace(rcd,aaa,bbb)
C      character*(*):: rcd,aaa,bbb
C      integer ii,lenaaa
C      ii = index(rcd, aaa)
C      lenaaa= len(aaa)
Cc      print *,'aaa',aaa,'!', aaalen
Cc      print *,'nrecsnrecs',ii, '!'//rcd(ii:ii+aaalen-1)//'!'
C      rcd(ii:ii+lenaaa-1)=bbb
C      end
C

