      module m_rdctrl2_func
      contains
      subroutine rdctrl2(recrd,recln,nrecs,
     .prgnam,vrsion,vn,vn2,
     .pass2,slabl_,v_sbz,v_sctrl,v_sham,v_spot,v_slat,v_smix,v_sspec,v_ssite,
     .v_sstr,v_sarry,v_smove,v_stb,sstrn)
      use m_rdctrl
      use m_gtv
      use m_struc_func
      use m_globalvariables
c      use m_susite_func, only: susite
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
      implicit none
C     include "mpif.h"
C ... Passed parameters
      integer:: recln,nrecs
      character slabl_(1)*8
      character strn*(recln)
      character*(*) recrd
      logical pass2
Cki      integer v_sbz,v_sctrl,v_sham,v_spot,v_slat,v_smix,v_sspec,
Cki     .  v_ssite,v_sstr,v_sarry,v_smove,v_stb
      type(s_bz):: v_sbz
      type(s_ctrl):: v_sctrl
      type(s_ham):: v_sham
      type(s_pot):: v_spot
      type(s_lat):: v_slat
      type(s_mix):: v_smix
      type(s_str):: v_sstr
      type(s_array):: v_sarry
      type(s_move):: v_smove
      type(s_tb):: v_stb
      type(s_spec),pointer:: v_sspec(:)
      type(s_site),pointer:: v_ssite(:)

      integer:: i_spec

      character  prgnam*(*), sstrn*(*)
c      character toksw(0:30)*(*), vrsion*6
C     character  vrsion*6
      character(6):: vrsion(2)
      double precision vn(2),vn2(2)
C ... Local parameters
      character fileid*64
      integer procid,nproc,master
      logical lgors,cmdopt,bittst,ltmp,ioorbp,asa
      double precision dval,dglob,xx(n0*2),dgets !,ekap(6)
      integer a2vec,bitand,fopna,getdig,i,is,igets,iprint,
     .iprt,irs(5),isw,ifi,ix(n0*nkap0),j,k,l,lasa,lbas,lcd,lqp,
     .lfrzw,lgunit,lmet,lncol, !,lham
     .lrs,lstsym,ltb,lxcf,nat,nlibu,nclasp,
     .nglob,nspc,nlmax,scrwid,stdo,stdl,stde,k1,k2,mpipid !nkap,lsx1,
      character*(8),allocatable::clabl(:)
      integer,allocatable:: ipc(:),initc(:),ics(:)
      real(8),allocatable:: pnuc(:,:,:),qnuc(:,:,:,:),
     .pp(:,:,:,:),ves(:),zc(:)
      integer:: dvec1(3)=1, dvec2(3)=0

C ... basis
      double precision orbp(n0,2,nkap0)
      integer o,oclabl,ohave,oics,opnu,opp,oqnu,osgw,osoptc,
     .oves,owk !osordn,
C ... Heap

c takao
      real(8):: pnux(20)
      integer:: nnn
c      integer:: umix_size,umove_size,uspec_size,usite_size,ubz_size,
c     &  uctrl_size,ulat_size,uham_size,uarray_size,upot_size,ustr_size
c
      integer:: i_copy_size,i_spacks,iendx,inix,i_spackv
      real(8):: seref
      integer:: ib , ispec 
      real(8),pointer :: rv_p_opos(:)
      integer,pointer :: iv_p_oips(:)
      integer,allocatable:: wowk(:)

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
c$$$        call upack('ctrl nbas nclass nspec nspin nl',v_sctrl,
c$$$     .    nsite,nclass,nspec,nsp,nl)
c$$$        nbas = nsite
c$$$        call upack('pot opnu oqnu oves',v_spot,opnu,oqnu,oves,0,0)
c$$$        call upack2('pot opp osoptc',v_spot,opp,osoptc)
c$$$        call upack('array oics oclabl nclasp ohave',v_sarry,
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
c$$$          call upack('spec lmxa z',v_sspec,is,
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
c$$$        call upack('ctrl lasa',v_sctrl,lasa,0,0,0,0)
c$$$        lasa = lasa + lasa3 + 8*isw(lasa8)
c$$$        call pack5('ctrl lasa',v_sctrl,lasa,0,0,0,0)
c$$$        return
c$$$      endif

C --- Read input parameters from contents of rcd ---
      call readctrl(prgnam,vrsion(1),vn(1))
      if (io_help > 0) then
        call readctrlpq(prgnam,nclasp,nl,nsp,pnuc,qnuc,pp,zc,
     .  ves,initc,ics,clabl)
        call cexit(0,1)
      endif

c      call defrr(owk,1)
C     Set switches depending type of program
      lbas = 0
      lcd = 0
C     For now, LMF => fp; no screening; nfp-style hamiltonian
C     cd represented in plane waves
      if (lfp .ne. 0) then
        lbas = 3
      endif
c      if (trim(prgnam) .eq. 'LMMC') then
c        lbas = 1
c      endif
      asa = .false.
c      if (prgnam .eq. 'LM' .or. prgnam .eq. 'LMGF' .or.
c     .    prgnam .eq. 'LMPG' .or. prgnam .eq. 'LMCTL') then
c        asa = .true.
c      endif

C ... Optionally read positions from pos file
      if (cmdopt('--rpos=',7,0,fileid)) then
        call iopos(.false.,-1,fileid(8:),nbasp,pos)
      endif

C ------------------- Copy to structures ----------------------
C     (query-replace-regexp "\\([a-z_]+\\)\\([0-9]+\\)" "\\2*isw(\\1\\2)" nil)

C ... Initialize sstrn
c      call ustrn(w,0,len(sstrn),0,0,0) !takao this does nothing.
      sstrn = ' '

C --- Allocate and copy input to sbz ---
Cki      call defdr(v_sbz,1000)
Cki      call ubz(v_sbz,-1,-000,0,0,0)
Cki      call redfrr(v_sbz,nint(dval(v_sbz,1)))

ctakao Kino's Cki error
      call ubz_init(v_sbz)
      v_sbz%size=ubz_size()


      v_sbz%def=bz_def
      i_copy_size=size(v_sbz%dosw)
      call dcopy(i_copy_size,bz_dosw,1,v_sbz%dosw,1)
      v_sbz%ef=bz_ef
      v_sbz%efmax=bz_efmax
      v_sbz%fsmom=bz_fsmom
      v_sbz%fsmommethod=bz_fsmommethod

      i = 1*isw(bz_lio1)+2*isw(bz_lio2) !+8*isw(bz_lio8)

      i_copy_size=size(v_sbz%lcond)
      call dcopy(i_copy_size,bz_lcond,1,v_sbz%lcond,1)
      v_sbz%lio=i
      v_sbz%lmet=bz_lmet
      v_sbz%lmull=bz_lmull
      v_sbz%lopt=0


      i_copy_size=size(v_sbz%lshft)
      call ii8copy(i_copy_size,bz_lshft,1,v_sbz%lshft,1)
      v_sbz%n=bz_n
      v_sbz%ndos=bz_ndos
      v_sbz%nevmx=bz_nevmx
      i_copy_size=size(v_sbz%nkabc)
      call ii8copy(i_copy_size,bz_nabc,1,v_sbz%nkabc,1)


      v_sbz%range=bz_range
      i_copy_size=size(v_sbz%semsh)
      call dcopy(i_copy_size,bz_semsh,1,v_sbz%semsh,1)
      v_sbz%w=bz_w
      v_sbz%zval=bz_zval


C --- Allocate and copy input to sctrl ---
Cki      call defdr(v_sctrl,1000)
Cki      call uctrl(v_sctrl,-1,0,0,0,0)
Cki      call redfrr(v_sctrl,nint(dval(v_sctrl,1)))
Cki      call uctrl(v_sctrl,-1,0,0,0,0)
      v_sctrl%size=uctrl_size()

C     lasa: 1 Make V from P,Q  2 Make pp  4 ccor  8 free atm
C          16 map  32 nonspherical mpol moms 64 MT corr
C         128 interpretation of sphere Q2; see newrho.f
C         256 how ASA Q1,Q2 are accumulated; see makwts.f
C         512 (spin pol) alpha repsn = (gamma(1) + gamma(nsp))/2
C$$$      lasa=4*isw(lasa4)+32*isw(lasa32)+64*isw(lasa64)+128*ham_qasa
C$$$     .    +512*isw(lasa512)
C$$$      if (.not. asa) lasa=0
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
      k = 2*isw(lcd2)+4*isw(lcd4) !+ 1*isw(lcd1)   +8*isw(lcd8)+64*isw(lcd64)

      i_copy_size=size(v_sctrl%defm)
      call dcopy(i_copy_size,lat_defm,1,v_sctrl%defm,1)
C$$$       v_sctrl%elin=asa_elin
C$$$       v_sctrl%lasa=lasa
      v_sctrl%lbas=j
      v_sctrl%lcd=k

c      k = isw(nmto.gt.1)
c      if (nmto.gt.1) k = k+2*isw(ham_ewald)

c       v_sctrl%lcgf=ctrl_lcgf
      v_sctrl%ldos=ctrl_ldos
      v_sctrl%lfp=lfp
      v_sctrl%lfrce=ctrl_lfrce
c       v_sctrl%lgen3= lgen3

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
c      lham = 1*isw(lham1)+4*isw(lham4)+8*isw(lham8)+16*isw(lham16)+
c     .      32*isw(lham32)+64*isw(lham64)+128*isw(lham128)+
c     .     256*isw(lham256) + lham3
c     lham = 256*isw(lham256)

c takao. we may need to remove lmet4 or lmet8... but a little complicated...
c Now ctrl_lmet4 and ctrl_lmet8 are fixed as .true.in m_rdctl.F
C     ctrl_lmet    1 metal  2 tetrahedron
Cr                 4 (PGF) V-shift1 is zero
Cr                 8 (PGF) V-shift2 is zero
      lmet = isw(bz_lmet.ne.0) + 2*isw(ctrl_lmet2) 
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
      lncol= 4*isw(lncol4)+32*isw(lncol32)+64*isw(lncol64)


c      lordn = 0

c       v_sctrl%lham=lham
      v_sctrl%lmet=lmet
      v_sctrl%lncol=lncol
c       v_sctrl%loptc=ctrl_loptc
c       v_sctrl%lordn=lordn

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
     .+ 8*getdig(irs(1),1,100)

ccccccccccccccccccccccccccccccccccccccccccc
c      print *,'xxxxxxxxxx irs(1)=',irs(1)
ccccccccccccccccccccccccccccccccccccccccccc

      lrs = 1*mod(irs(1),8)+8*irs(2)+32*irs(3)+64*irs(4)+128*irs(5)
     .+ 256*mod(irs(1)/8,2)
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

c        i_copy_size=size(v_sctrl%lpgf)
c        call ii8copy(i_copy_size,ctrl_lpgf,1,v_sctrl%lpgf,1)
      v_sctrl%lqp=lqp
      v_sctrl%lrel=lrel
      v_sctrl%lrs=lrs
      v_sctrl%lscr=lscr


C      lstr  no longer used.
C      ltb switches for empirical tight-binding
C         1 overlap        2 crystal-field     4 ovlp+CF
C         8 add ebarLL    16 forces           32 fij
C        64 not used     128 pressure        256 evdisc
C       512 pair pot    1024 TrH & local E  2048 local rho
C      2^12 Hubbard U   2^13 No Madelung    2^14 wgt avg U
C      2^15 L>0 estat   2^16 disc read incr 2^17 gamma-pt
C$$$       i = 1*isw(ltb1)+2*isw(ltb2)+4*isw(ltb4)+8*isw(ltb8)+16*isw(ltb16)
C$$$     .   +32*isw(ltb32)+64*isw(ltb64)+128*isw(ltb128)+256*isw(ltb256)
C$$$     .   +512*isw(ltb512)+1024*isw(ltb1024)+2048*isw(ltb2048)
C$$$     .   +4096*isw(ltb4096)+2**13*isw(ltb213)+2**14*isw(ltb214)
C$$$     .   +2**15*isw(ltb215)+2**16*isw(ltb216)+2**17*isw(ltb217)
C$$$     .   +2**18*isw(ltb218)
C      lves    1 take ves as input

C$$$        i_copy_size=size(v_sctrl%lstonr)
C$$$        call ii8copy(i_copy_size,lstonr,1,v_sctrl%lstonr,1)
      v_sctrl%lstr=0
c        v_sctrl%lsx=lsx
c        v_sctrl%ltb=i
      v_sctrl%lves=lves

C      lxcf   parameter defining XC functional
C      1s digit:
C      1 for Ceperly-Alder (VWN)
C      2 for Barth-Hedin (ASW fit)
C      103 for PBE-GGA
      v_sctrl%lxcf= ham_lxcf


!! ==  Set modep ===
      ix(1:3) = 2
c      if (prgnam == 'LMPG') ix(3) = 0
c      if (prgnam == 'LMMC') ix(1:3) = 0

      v_sctrl%maxit=iter_maxit
      i_copy_size=size(v_sctrl%mdprm)
      call dcopy(i_copy_size,mdprm,1,v_sctrl%mdprm,1)
      i_copy_size=size(v_sctrl%modep)
      call ii8copy(i_copy_size,ix(1:3),1,v_sctrl%modep,1)
      v_sctrl%nbas=nbas

C      Reset nl
      if (nl .ne. max(lmxbx,lmxax)+1 .and. io_help == 0) then
        call info2(20,1,0,' rdctrl: reset global max nl from %i to %i',
     .  nl,max(lmxbx,lmxax)+1)
        nl = max(lmxbx,lmxax)+1
      endif
C      call pack5('ctrl nbasp nclass nesabc nitmv nl',v_sctrl,

      v_sctrl%nbasp=nbasp
      i_copy_size=size(v_sctrl%nesabc)
      call ii8copy(i_copy_size,nesabc,1,v_sctrl%nesabc,1)
      v_sctrl%nitmv=nitmv
      v_sctrl%nl=nl

C      call pack5('ctrl nlibu nmap npl nsite nspec',v_sctrl,

      v_sctrl%nsite=nsite
      v_sctrl%nspec=nspec

C      call shosyv(0,0,0,6)

      v_sctrl%nspin=nsp
      v_sctrl%nvario=nvario
      i_copy_size=size(v_sctrl%omax1)
      call dcopy(i_copy_size,omax1,1,v_sctrl%omax1,1)
      i_copy_size=size(v_sctrl%omax2)
      call dcopy(i_copy_size,omax2,1,v_sctrl%omax2,1)
      v_sctrl%quit=quit


      v_sctrl%rmaxes=rmaxes
      v_sctrl%rmines=rmines
      v_sctrl%sclwsr=sclwsr
      v_sctrl%sdmod=sdmod
      i_copy_size=size(v_sctrl%sdprm)
      call dcopy(i_copy_size,sdprm,1,v_sctrl%sdprm,1)


      i_copy_size=size(v_sctrl%smalit)
      call ii8copy(i_copy_size,smalit,1,v_sctrl%smalit,1)
      i_copy_size=size(v_sctrl%tol)
      call dcopy(i_copy_size,ctrl_tol,1,v_sctrl%tol,1)
      v_sctrl%wsrmax=wsrmax
      i_copy_size=size(v_sctrl%zbak)
      call dcopy(i_copy_size,zbak,1,v_sctrl%zbak,1)
      v_sctrl%pfloat=lpfloat


C --- Allocate and copy input to slat ---
CKi      call defdr(v_slat,1000)
CKi      call ulat(v_slat,-1,0,0,0,0)
CKi      call redfrr(v_slat,nint(dval(v_slat,1)))
      v_slat%size=ulat_size()

      if (dalat == NULLR) dalat=0

      v_slat%alat=alat+dalat
      v_slat%as=lat_as
      v_slat%avw=avw
      v_slat%nkdmx=lat_nkdmx
      v_slat%nkqmx=lat_nkdmx


      i_copy_size=size(v_slat%gam)
      call dcopy(i_copy_size,lat_gam,1,v_slat%gam,1)
      v_slat%gmax=lat_gmax
      i_copy_size=size(v_slat%nabc)
      call ii8copy(i_copy_size,ftmesh,1,v_slat%nabc,1)


c takao I think gam is alreay setted above.
c       i_copy_size=size(v_slat%gam)
c       call dcopy(i_copy_size,lat_gam,1,v_slat%gam,1)

      v_slat%ldist=lat_ldist
      i_copy_size=size(v_slat%dist)
      call dcopy(i_copy_size,lat_dist,1,v_slat%dist,1)


      i_copy_size=size(v_slat%plat)
      call dcopy(i_copy_size,plat,1,v_slat%plat,1)
c       i_copy_size=size(v_slat%plat2)
c       call dcopy(i_copy_size,slat_plat2,1,v_slat%plat2,1)


      v_slat%rpad=lat_rpad
      i_copy_size=size(v_slat%slat)
      call dcopy(i_copy_size,lat_slat,1,v_slat%slat,1)
      v_slat%tol=lat_tol
      v_slat%tolft=tolft
      v_slat%vol=vol

c      call dcopy(6,plat,1,xx,1)
c      call dcopy(3,platl,1,xx(7),1)

c       i_copy_size=size(v_slat%platl)
c       call dcopy(i_copy_size,xx,1,v_slat%platl,1)

c      call dcopy(3,platr,1,xx(7),1)

c       i_copy_size=size(v_slat%platr)
c       call dcopy(i_copy_size,xx,1,v_slat%platr,1)


C --- Allocate and copy input to sham ---
CKi      call defdr(v_sham,1000)
CKi      call uham(v_sham,-1,-000,0,0,0)
CKi      call redfrr(v_sham,nint(dval(v_sham,1)))
      v_sham%size=uham_size()


      v_sham%alfsi=alfsi
      i_copy_size=size(v_sham%dabc)
      call dcopy(i_copy_size,dabc,1,v_sham%dabc,1)
      v_sham%elind=elind
      v_sham%nmto=nmto
      i_copy_size=size(v_sham%kmto)
      call dcopy(i_copy_size,kmto,1,v_sham%kmto,1)

C     Replicate ctrl->lncol in ham->lncol, ditto for lham,lgen3

      v_sham%lncol=(int(v_sctrl%lncol))


c       v_sham%lham=(int(v_sctrl%lham))

c       v_sham%lgen3=(int(v_sctrl%lgen3))

C     Mop up integer parts of sigp
      sigp(1) = sigp_mode
      sigp(2) = sigp_nmin
      sigp(4) = sigp_nmax

      v_sham%nkaph=nkaph
      i_copy_size=size(v_sham%pmax)
      call dcopy(i_copy_size,pmax,1,v_sham%pmax,1)
      i_copy_size=size(v_sham%pmin)
      call dcopy(i_copy_size,pmin,1,v_sham%pmin,1)
      v_sham%lsig=lrsig
      i_copy_size=size(v_sham%sigp)
      call dcopy(i_copy_size,sigp,1,v_sham%sigp,1)


      i_copy_size=size(v_sham%qss)
      call dcopy(i_copy_size,ham_qss,1,v_sham%qss,1)
      v_sham%rsrnge=rsrnge
      v_sham%rsstol=rsstol
      v_sham%udiag=ham_udiag

C     Parameters for APW

      v_sham%pwmode=pwmode
      v_sham%npwpad=npwpad
      v_sham%pwemin=pwemin
      v_sham%pwemax=pwemax


      v_sham%oveps=oveps
      v_sham%delta_stabilize=delta_stabilize !takao sep2010


C ... Allocate and initialize sarray
CKi      call defrr(v_sarry,1000)
CKi      call uarray(v_sarry,-1,NULLI,0,0,0)
CKi      call redfrr(v_sarry,nint(dval(v_sarry,1)))
      v_sarry%size=uarray_size()

C ... Allocate and initialize spot
CKi      call defdr(v_spot,1000)
CKi      call upot(v_spot,-1,-000,0,0,0)
CKi      call redfrr(v_spot,nint(dval(v_spot,1)))
      v_spot%size=upot_size()


      v_spot%vmtz0=vmtz


C --- Allocate and copy input to smix ---
CKi      call defrr(v_smix,1000)
CKi      call umix(v_smix,-1,-000,0,0,0)
CKi      call redfrr(v_smix,nint(dval(v_smix,1)))
      v_smix%size=umix_size()

ckino's fix for SR11000 size() return's integer(8), right?
c      nnn = size(smix)
c      nnn = umix_size()
c      call dcopy(nnn,vmix,1,v_smix,1)
      v_smix=vmix
      

C ... Allocate and initialize sordn
Cki delete sordn
c#if 1
c       v_sarry%osordn=osordn
c#else
c      call defrr(osordn,1000)
c      call uordn(w(osordn),-1,0,0,0,0)
c      call redfrr(osordn,nint(dval(w(osordn),1)))
cCgetarg...info...           integer :: osordn , integer(8) :: v_sarry%osordn 1
cCgetarg       call pack1('array osordn',v_sarry,osordn)
c
c       v_sarry%osordn=osordn
c
c#endif

C --- Allocate and copy input to smove ---
Cki      call defrr(v_smove,1000)
Cki      call umove(v_smove,-1,0,0,0,0)
Cki      call redfrr(v_smove,nint(dval(v_smove,1)))
      v_smove%size=umove_size()

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

      v_smove%gyro=2d0
      v_smove%nmodt=gd_nmodt
      i_copy_size=size(v_smove%modt)
      call ii8copy(i_copy_size,gd_modt,1,v_smove%modt,1)
      i_copy_size=size(v_smove%ct)
      call dcopy(i_copy_size,gd_ct,1,v_smove%ct,1)
      i_copy_size=size(v_smove%prmint)
      call dcopy(i_copy_size,prmint,1,v_smove%prmint,1)


      v_smove%kt=move_kt
      v_smove%ts=move_ts
      v_smove%tsequ=move_tsequ
      v_smove%tstot=move_tstot


C --- Allocate and copy input to soptic ---
Cki  delete soptc
C$$$#if 1
C$$$       osoptc=0
C$$$       v_spot%osoptc=osoptc
C$$$#else
C$$$      call defrr(osoptc,1000)
C$$$      call uoptic(w(osoptc),-1,0,0,0,0)
C$$$      call redfrr(osoptc,nint(dval(w(osoptc),1)))
C$$$Cgetarg...info...           integer :: osoptc , integer(8) :: v_spot%osoptc 1
C$$$Cgetarg       call pack1('pot osoptc',v_spot,osoptc)

C$$$       v_spot%osoptc=osoptc

C$$$Cgetarg...info...           NG:: optic_ne , integer(8) :: w(osoptc)%ne 1
C$$$Cgetarg...info...           NG:: optic_window , real(8) :: w(osoptc)%window 2
C$$$Cgetarg...info...           NG:: optic_ocrng , integer(8) :: w(osoptc)%ocrng 2
C$$$Cgetarg...info...           NG:: optic_unrng , integer(8) :: w(osoptc)%unrng 2
C$$$Cgetarg...info...           NG:: optic_esciss , real(8) :: w(osoptc)%esciss 1
C$$$Cgetarg       call pack5('optic ne window ocrng unrng esciss',w(osoptc),
C$$$Cgetarg      .  optic_ne,optic_window,optic_ocrng,optic_unrng,optic_esciss)

C$$$       w(osoptc)%ne=optic_ne
C$$$       i_copy_size=size(w(osoptc)%window)
C$$$       call dcopy(i_copy_size,optic_window,1,w(osoptc)%window,1)
C$$$       i_copy_size=size(w(osoptc)%ocrng)
C$$$       call ii8copy(i_copy_size,optic_ocrng,1,w(osoptc)%ocrng,1)
C$$$       i_copy_size=size(w(osoptc)%unrng)
C$$$       call ii8copy(i_copy_size,optic_unrng,1,w(osoptc)%unrng,1)
C$$$       w(osoptc)%esciss=optic_esciss

C$$$Cgetarg...info...           integer :: isw , integer(8) :: w(osoptc)%mode 1
C$$$Cgetarg...info...           NG:: optic_nchi2 , integer(8) :: w(osoptc)%nchi2 1
C$$$Cgetarg       call pack2('optic mode nchi2',w(osoptc),isw(optic_mode1),
C$$$Cgetarg      .  optic_nchi2)

C$$$       w(osoptc)%mode=isw(optic_mode1)
C$$$       w(osoptc)%nchi2=optic_nchi2

C$$$Cgetarg...info...           NG:: optic_axes , integer(8) :: w(osoptc)%axes 18
C$$$Cgetarg       call pack1('optic axes',w(osoptc),optic_axes) ! nchi2 must be

C$$$       i_copy_size=size(w(osoptc)%axes)
C$$$       call ii8copy(i_copy_size,optic_axes,1,w(osoptc)%axes,1)

C$$$#endif
      ! packed beforehand
C --- Allocate and copy input to sgw ---
Cki delete sgw
C$$$#if 1
C$$$       osgw=0
C$$$       v_spot%osgw=osgw
C$$$#else
C$$$      call defrr(osgw,1000)
C$$$      call ugw(w(osgw),-1,0,0,0,0)
C$$$      call redfrr(osgw,nint(dval(w(osgw),1)))
C$$$C     Old defaults
C$$$C     call pack5('gw gcutb gcutx qoffp nband gsmear',w(osgw),2.7d0,
C$$$C    .  2.2d0,1d0,9999,.003d0)
C$$$CC    Default : Faleev's fast energy integration
C$$$C     xx(1) = .01d0
C$$$C     xx(2) = .04d0
C$$$C     call pack5('gw lgw delre ecuts',w(osgw),1,xx,2.5d0,0,0)
C$$$C     call pack5('gw nime delre deltax deltaw pbtol',w(osgw),6,xx,
C$$$C    .  -1d-4,.02d0,1d-3)
C$$$C     call spacks(1,'gw pb1',w(osgw),'111',0,0)
C$$$C     call spacks(1,'gw pb2',w(osgw),'1111',0,0)
C$$$Cgetarg...info...           NG:: gw_gcutb , real(8) :: w(osgw)%gcutb 1
C$$$Cgetarg...info...           NG:: gw_gcutx , real(8) :: w(osgw)%gcutx 1
C$$$Cgetarg...info...           NG:: gw_qoffp , real(8) :: w(osgw)%qoffp 1
C$$$Cgetarg...info...           NG:: gw_nband , integer(8) :: w(osgw)%nband 1
C$$$Cgetarg...info...           NG:: gw_gsmear , real(8) :: w(osgw)%gsmear 1
C$$$Cgetarg       call pack5('gw gcutb gcutx qoffp nband gsmear',w(osgw),
C$$$Cgetarg      .  gw_gcutb,gw_gcutx,gw_qoffp,gw_nband,gw_gsmear)

C$$$       w(osgw)%gcutb=gw_gcutb
C$$$       w(osgw)%gcutx=gw_gcutx
C$$$       w(osgw)%qoffp=gw_qoffp
C$$$       w(osgw)%nband=gw_nband
C$$$       w(osgw)%gsmear=gw_gsmear

C$$$Cgetarg...info...           NG:: gw_nabc , integer(8) :: w(osgw)%nkabc 3
C$$$Cgetarg...info...           NG:: gw_delre , real(8) :: w(osgw)%delre 2
C$$$Cgetarg...info...           NG:: gw_ecuts , real(8) :: w(osgw)%ecuts 1
C$$$Cgetarg       call pack5('gw lgw nkabc delre ecuts',w(osgw),
C$$$Cgetarg      .  1,gw_nabc,gw_delre,gw_ecuts,0)

C$$$       w(osgw)%lgw=1
C$$$       i_copy_size=size(w(osgw)%nkabc)
C$$$       call ii8copy(i_copy_size,gw_nabc,1,w(osgw)%nkabc,1)
C$$$       i_copy_size=size(w(osgw)%delre)
C$$$       call dcopy(i_copy_size,gw_delre,1,w(osgw)%delre,1)
C$$$       w(osgw)%ecuts=gw_ecuts

C$$$Cgetarg...info...           NG:: gw_nime , integer(8) :: w(osgw)%nime 1
C$$$Cgetarg...info...           NG:: gw_deltax , real(8) :: w(osgw)%deltax 1
C$$$Cgetarg...info...           NG:: gw_deltaw , real(8) :: w(osgw)%deltaw 1
C$$$Cgetarg...info...           NG:: gw_pbtol , real(8) :: w(osgw)%pbtol 1
C$$$Cgetarg       call pack5('gw nime deltax deltaw pbtol',w(osgw),
C$$$Cgetarg      .  gw_nime,gw_deltax,gw_deltaw,gw_pbtol,0)

C$$$       w(osgw)%nime=gw_nime
C$$$       w(osgw)%deltax=gw_deltax
C$$$       w(osgw)%deltaw=gw_deltaw
C$$$       w(osgw)%pbtol=gw_pbtol

C$$$Cgetarg...info...          check_variable '111' not found
C$$$Cgetarg       call spacks(1,'gw pb1',w(osgw),'111',0,0)
C$$$       do i_spacks=0,0
C$$$       call spacks_copy('p',w(osgw)(i_spacks)%pb1,0,0,'111',i_spacks)
C$$$       enddo

C$$$Cgetarg...info...          check_variable '1111' not found
C$$$Cgetarg       call spacks(1,'gw pb2',w(osgw),'1111',0,0)
C$$$       do i_spacks=0,0
C$$$       call spacks_copy('p',w(osgw)(i_spacks)%pb2,0,0,'1111',i_spacks)
C$$$       enddo

C$$$Cgetarg...info...           integer :: osgw , integer(8) :: v_spot%osgw 1
C$$$Cgetarg       call pack1('pot osgw',v_spot,osgw)

C$$$       v_spot%osgw=osgw

C$$$#endif
C --- Allocate and copy input to sstr ---
CKi      call defrr(v_sstr,1000)
CKi      call ustr(v_sstr,-1,0,0,0,0)
CKi      call redfrr(v_sstr,nint(dval(v_sstr,1)))
      v_sstr%size=ustr_size()

C     Old defaults
C     xx(1) = 0
C     xx(2) = -1
C     xx(3) = -2.3d0
C     if (prgnam .eq. 'TBE') then
C       call pack5('str nkaps rmax rfit kaps lmaxw',v_sstr,1,
C    .    0d0,.8d0,xx,-1)
C     else
C       call pack5('str nkaps rmax rfit kaps lmaxw',v_sstr,1,
C    .    0d0,.8d0,xx,-1)
C     endif
C     call pack2('str loka drwats',v_sstr,1,.1d0)
C    loka switch
c      i = 1
c      if (str_mode .eq. 1 .or. str_mode .eq. 3) i = 0
C     Pack iinv parameters
c      call dpzero(xx,5)
c      xx(1) = iinv_nit
c      xx(2) = iinv_ncut
c      xx(3) = iinv_tol

c       v_sstr%adec=tcf_adec
c       v_sstr%amode=str_mode ! amode is dummy now
c       v_sstr%drwats=str_drwats
c       i_copy_size=size(v_sstr%iinv)
c       call dcopy(i_copy_size,xx,1,v_sstr%iinv,1)


c       v_sstr%lequiv=isw(str_lequiv1)
c       v_sstr%loka=i
c       v_sstr%lmaxw=str_lmaxw
c       v_sstr%lshow=isw(str_lshow1)
      v_sstr%mxnbr=str_mxnbr


c       v_sstr%nalf=tcf_nalf
c       i_copy_size=size(v_sstr%nbisi)
c       call ii8copy(i_copy_size,tcf_nbisi,1,v_sstr%nbisi,1)
c       v_sstr%ncupl=tcf_ncupl
c       v_sstr%ndust=tcf_ndust
c       v_sstr%wztcf=tcf_wztcf


c       v_sstr%nkaps=str_nkaps
      v_sstr%rmax=str_rmax
c       v_sstr%rfit=0.8d0
c       i_copy_size=size(v_sstr%kaps)
c       call dcopy(i_copy_size,str_kaps,1,v_sstr%kaps,1)


c       v_sstr%rmaxg=str_rmaxg
c       v_sstr%ivl=str_ivl


C$$$C --- Allocate and copy input to stb ---
C$$$Cki      call defrr(v_stb,1000)
C$$$Cki      call utb(v_stb,-1,NULLI,0,0,0)
C$$$Cki      call redfrr(v_stb,nint(dval(v_stb,1)))
C$$$        v_stb%size=utb_size()

C$$$C     old defaults
C$$$C      call pack1('tb alam',v_stb,.001d0)
C$$$C      call pack1('tb alsc',v_stb,10d0)
C$$$Cgetarg...info...           NG:: tb_alam , real(8) :: v_stb%alam 1
C$$$Cgetarg...info...           NG:: tb_alsc , real(8) :: v_stb%alsc 1
C$$$Cgetarg       call pack5('tb alam alsc',v_stb,
C$$$Cgetarg      .  tb_alam,tb_alsc,0,0,0)

C$$$       v_stb%alam=tb_alam
C$$$       v_stb%alsc=tb_alsc

C$$$Cgetarg...info...           NG:: tb_ebfit , real(8) :: v_stb%ebfit 2
C$$$Cgetarg...info...           NG:: tb_fmode , integer(8) :: v_stb%fmode 1
C$$$Cgetarg...info...           NG:: tb_nbfit , integer(8) :: v_stb%nbfit 2
C$$$Cgetarg...info...           NG:: tb_rmfit , real(8) :: v_stb%rmfit 2
C$$$Cgetarg       call pack5('tb ebfit fmode nbfit rmfit',v_stb,
C$$$Cgetarg      .  tb_ebfit,tb_fmode,tb_nbfit,tb_rmfit,0)

C$$$       i_copy_size=size(v_stb%ebfit)
C$$$       call dcopy(i_copy_size,tb_ebfit,1,v_stb%ebfit,1)
C$$$       v_stb%fmode=tb_fmode
C$$$       i_copy_size=size(v_stb%nbfit)
C$$$       call ii8copy(i_copy_size,tb_nbfit,1,v_stb%nbfit,1)
C$$$       i_copy_size=size(v_stb%rmfit)
C$$$       call dcopy(i_copy_size,tb_rmfit,1,v_stb%rmfit,1)



C --- Allocate and copy input to sspec ---
Cki      call defrr(v_sspec,1000)
Cki      call uspec(v_sspec,-1,0,1,0,0,0)
Cki      call redfrr(v_sspec,nspec*nint(dval(v_sspec,1)))
Cki      call uspec(v_sspec,-1,-000,nspec,0,0,0)
      allocate(v_sspec(nspec))
      do i=1,nspec; v_sspec(i)%size=uspec_size(); enddo

C ... old defaults
C      call spackv(1,'spec rs3',v_sspec,1,nspec,1d0)
C      call spackv(1,'spec vmtz',v_sspec,1,nspec,-0.5d0)
CC ... Init default for sigma,lmxa,lmxf
C      call dvset(xx,1,10,-1d0)
CC ... First set lmxa to 4 to set def. hcr to l=4
C      call spackv(1,'spec lmxa',v_sspec,1,nspec,4)
C      call spackv(1,'spec lmxpb',v_sspec,1,nspec,3)
C      call spackv(1,'spec hcr',v_sspec,1,nspec,xx)
C      call spackv(1,'spec lmxa',v_sspec,1,nspec,nl-1)
C      call spackv(1,'spec lmxf',v_sspec,1,nspec,2*nl-2)
C      xx(1) = -1
C      xx(2) = 0
C      call spackv(1,'spec coreq',v_sspec,1,nspec,xx)
C      if (prgnam(1:2) .eq. 'TB')
C     .  call spackv(1,'spec rmt',v_sspec,1,nspec,1d0)
C      do  j = 1, nspec
CC       call spacks(1,'spec pb1',v_sspec,'111     ',j,j)
CC       call spacks(1,'spec pb2',v_sspec,'1111    ',j,j)
C        call spacks(1,'spec pb1',v_sspec,'        ',j,j)
C        call spacks(1,'spec pb2',v_sspec,'        ',j,j)
C        call spacks(1,'spec coreh',v_sspec,'        ',j,j)
C      enddo
      i_copy_size=1;
      do i_spackv=1,nspec
        call spackv_array_copy_r8_r8('p',v_sspec(i_spackv)%eh3,i_copy_size,1,-0.5d0)
      enddo

      i_copy_size=1;
      do i_spackv=1,nspec
        call spackv_array_copy_r8_r8('p',v_sspec(i_spackv)%etf,i_copy_size,1,-1d0)
      enddo

      i_copy_size=1;
      do i_spackv=1,nspec
        call spackv_array_copy_i8_i('p',v_sspec(i_spackv)%lmxf,i_copy_size,1,2*nl-2)
      enddo

      i_copy_size=1;
      do i_spackv=1,nspec
        call spackv_array_copy_i8_i('p',v_sspec(i_spackv)%norp,i_copy_size,1,2)
      enddo

      i_copy_size=1;
      do i_spackv=1,nspec
        call spackv_array_copy_r8_r8('p',v_sspec(i_spackv)%vmtz,i_copy_size,1,-0.5d0)
      enddo


      do  j = 1, nspec
        slabl_(j) = slabl(j)
        do i_spacks=j,j
          call spacks_copy('p',v_sspec(i_spacks)%pb1,j,j,pb1(j),i_spacks)
        enddo

        do i_spacks=j,j
          call spacks_copy('p',v_sspec(i_spacks)%pb2,j,j,pb2(j),i_spacks)
        enddo

        do i_spacks=j,j
          call spacks_copy('p',v_sspec(i_spacks)%coreh,j,j,coreh(j),i_spacks)
        enddo

        do i_spacks=j,j
          call spacks_copy('p',v_sspec(i_spacks)%name,j,j,slabl(j),i_spacks)
        enddo


        v_sspec(j)%a=spec_a(j)
        v_sspec(j)%nr=nr(j)
        i_copy_size=size(v_sspec(j)%alpha)
        call dcopy(i_copy_size,alpha(1,j),1,v_sspec(j)%alpha,1)
        i_copy_size=size(v_sspec(j)%coreq)
        call dcopy(i_copy_size,coreq(1,j),1,v_sspec(j)%coreq,1)


        v_sspec(j)%lxi=lxi(j)
        v_sspec(j)%nxi=nxi(j)
        i_copy_size=size(v_sspec(j)%exi)
        call dcopy(i_copy_size,exi(1,j),1,v_sspec(j)%exi,1)


        v_sspec(j)%group=grp(j)
        v_sspec(j)%grp2=grp2(j)
        i_copy_size=size(v_sspec(j)%hcr)
        call dcopy(i_copy_size,hcr(1,j),1,v_sspec(j)%hcr,1)
        i_copy_size=size(v_sspec(j)%idmod)
        call ii8copy(i_copy_size,idmod(1,j),1,v_sspec(j)%idmod,1)


        i_copy_size=size(v_sspec(j)%ehvl)
        call dcopy(i_copy_size,ehvl(1,j),1,v_sspec(j)%ehvl,1)

C       Set idxdn
        call ivset(ix,1,n0*nkap0,1)
        call icopy(1+lmxb(j),idxdn(1,j),1,ix,1)

        i_copy_size=size(v_sspec(j)%idxdn)
        call ii8copy(i_copy_size,ix,1,v_sspec(j)%idxdn,1)
        i_copy_size=size(v_sspec(j)%idu)
        call ii8copy(i_copy_size,idu(1,j),1,v_sspec(j)%idu,1)
        i_copy_size=size(v_sspec(j)%jh)
        call dcopy(i_copy_size,jh(1,j),1,v_sspec(j)%jh,1)
        i_copy_size=size(v_sspec(j)%uh)
        call dcopy(i_copy_size,uh(1,j),1,v_sspec(j)%uh,1)


        v_sspec(j)%kmxt=kmxt(j)
        v_sspec(j)%kmxv=kmxv(j)
        v_sspec(j)%lfoca=lfoca(j)
        v_sspec(j)%rsmv=rsmv(j)


        v_sspec(j)%lmxa=lmxa(j)
        v_sspec(j)%lmxb=lmxb(j)
        v_sspec(j)%lmxl=lmxl(j)
        v_sspec(j)%lmxpb=lmxpb(j)

C       pack mxcst(j)
        i = 1*isw(mxcst1(j))+2*isw(mxcst2(j))+4*isw(mxcst4(j))
ccccccccccccccccccccccc
        print *,' mxcst switch =',j,i,mxcst1(j),mxcst2(j),mxcst4(j)
c        i=0
ccccccccccccccccccccccc
C       pack orbp(j)
        call dpzero(orbp,n0*2*nkap0)
        call dcopy(n0,rsmh(1,j),1,orbp(1,1,1),1)
        call dcopy(n0,eh(1,j),1,orbp(1,2,1),1)
        call dcopy(n0,rsmh2(1,j),1,orbp(1,1,2),1)
        call dcopy(n0,eh2(1,j),1,orbp(1,2,2),1)

        v_sspec(j)%mass=mass(j)
        v_sspec(j)%mxcst=i
        i_copy_size=size(v_sspec(j)%orbp)
        call dcopy(i_copy_size,orbp,1,v_sspec(j)%orbp,1)

C       Pack P,Q,PZ for both spins
        call dpzero(orbp,n0*6)
        call dcopy(n0*nsp,pnu(1,1,j),1,orbp(1,1,1),1)
        call dcopy(n0*nsp,pz(1,1,j),1,orbp(1,1,2),1)
        call dcopy(n0*nsp,qnu(1,1,j),1,orbp(1,1,3),1)

        i_copy_size=size(v_sspec(j)%p)
        call dcopy(i_copy_size,orbp,1,v_sspec(j)%p,1)
        i_copy_size=size(v_sspec(j)%pz)
        call dcopy(i_copy_size,orbp(1,1,2),1,v_sspec(j)%pz,1)
        i_copy_size=size(v_sspec(j)%q)
        call dcopy(i_copy_size,orbp(1,1,3),1,v_sspec(j)%q,1)
        v_sspec(j)%z=z(j)


        i_copy_size=size(v_sspec(j)%colxbs)
        call dcopy(i_copy_size,colxbs(1,j),1,v_sspec(j)%colxbs,1)
        v_sspec(j)%radxbs=radxbs(j)

c        call pack5('spec rcfa rcut rfoca rg rham',v_sspec,j,

        i_copy_size=size(v_sspec(j)%rcfa)
        call dcopy(i_copy_size,rcfa(1,j),1,v_sspec(j)%rcfa,1)
        v_sspec(j)%rcut=rcut(j)
        v_sspec(j)%rfoca=rfoca(j)
        v_sspec(j)%rg=rg(j)


        v_sspec(j)%rmt=rmt(j)
        v_sspec(j)%rs3=rs3(j)
        v_sspec(j)%rsma=rsma(j)
        v_sspec(j)%rsmfa=rsmfa(j)

C$$$        if (ltbe) then
C$$$Cgetarg...info...           NG:: iq1(1,j) , integer(8) :: v_sspec%iq1 4
C$$$Cgetarg...info...           NG:: ivso(1,j) , integer(8) :: v_sspec%ivso 4
C$$$Cgetarg...info...           NG:: stni(j) , real(8) :: v_sspec%stni 1
C$$$Cgetarg...info...           NG:: tbvso(1,j) , real(8) :: v_sspec%vso 4
C$$$Cgetarg           call pack5('spec iq1 ivso stni vso',v_sspec,j,
C$$$Cgetarg      .      iq1(1,j),ivso(1,j),stni(j),tbvso(1,j))

C$$$           i_copy_size=size(v_sspec(j)%iq1)
C$$$           call ii8copy(i_copy_size,iq1(1,j),1,v_sspec(j)%iq1,1)
C$$$           i_copy_size=size(v_sspec(j)%ivso)
C$$$           call ii8copy(i_copy_size,ivso(1,j),1,v_sspec(j)%ivso,1)
C$$$           v_sspec(j)%stni=stni(j)
C$$$           i_copy_size=size(v_sspec(j)%vso)
C$$$           call dcopy(i_copy_size,tbvso(1,j),1,v_sspec(j)%vso,1)

C$$$Cgetarg...info...           NG:: qpol(1,j) , real(8) :: v_sspec%qpol 10
C$$$Cgetarg           call pack5('spec qpol',v_sspec,j,qpol(1,j),0,0,0)

C$$$           i_copy_size=size(v_sspec(j)%qpol)
C$$$           call dcopy(i_copy_size,qpol(1,j),1,v_sspec(j)%qpol,1)

C$$$          if (ltb217 .or. ltb218) then
C$$$Cgetarg...info...           integer :: dvec1(3) , integer(8) :: v_sbz%nkabc 3
C$$$Cgetarg            call pack1('bz nkabc',v_sbz,dvec1)

C$$$            i_copy_size=size(v_sbz%nkabc)
C$$$            call ii8copy(i_copy_size,dvec1,1,v_sbz%nkabc,1)

C$$$Cgetarg...info...           integer :: dvec2(3) , integer(8) :: v_sbz%lshft 3
C$$$Cgetarg            call pack1('bz lshft',v_sbz,dvec2)

C$$$            i_copy_size=size(v_sbz%lshft)
C$$$            call ii8copy(i_copy_size,dvec2,1,v_sbz%lshft,1)

C$$$          endif
C$$$        endif

        v_sspec(j)%dv=dv(j)
        v_sspec(j)%eref=eref(j)
        v_sspec(j)%rham=rham(j)
        v_sspec(j)%rint=rint(j)

      enddo
C     call shstru('spec',v_sspec,1,nspec);stop
cccccccccccccc
c       call upack('spec p ',v_sspec,1,pnux,0,0,0)
c       print *,'xxx:111 222 ini pnu=',pnu
cccccccccccccc


C --- Allocate and copy input to ssite ---
C     NB: essential that ssite is LAST array allocated here
CKi      call defrr(v_ssite,1000)
CKi      call usite(v_ssite,-1,0,1,0,0,0)
CKi      call redfrr(v_ssite,nsite*nint(dval(v_ssite,1)))
CKi      call usite(v_ssite,-1,0,nsite,0,0,0)
      allocate(v_ssite(nsite))
      do i=1,nsite; v_ssite(i)%size=usite_size(); enddo

C     Old defaults
C     call ivset(ix,1,n0*nkap0,1)
C     call spackv(1,'site relax',v_ssite,1,nsite,ix)
      do  j = 1, nsite
C       Pack spec into class in case spec info is sufficient
C       and spid into class label

        v_ssite(j)%spec=ips(j)
        v_ssite(j)%class=ips(j)
        i_copy_size=size(v_ssite(j)%dpole)
        call dcopy(i_copy_size,dpole(1,j),1,v_ssite(j)%dpole,1)
        v_ssite(j)%mpole=mpole(j)

        do i_spacks=j,j
          call spacks_copy('p',v_ssite(i_spacks)%clabel,j,j,slabl(ips(j)),i_spacks)
        enddo


        v_ssite(j)%pl=ipl(j)
        v_ssite(j)%plv=plv(j)
        i_copy_size=size(v_ssite(j)%pos)
        call dcopy(i_copy_size,pos(1,j),1,v_ssite(j)%pos,1)
        i_copy_size=size(v_ssite(j)%vel)
        call dcopy(i_copy_size,vel(1,j),1,v_ssite(j)%vel,1)


        v_ssite(j)%vshft=vshft(j)
        i_copy_size=size(v_ssite(j)%relax)
        call ii8copy(i_copy_size,irlx(1,j),1,v_ssite(j)%relax,1)
        i_copy_size=size(v_ssite(j)%eula)
        call dcopy(i_copy_size,eula(1,j),1,v_ssite(j)%eula,1)

C$$$        if (ltbe) then
C$$$Cgetarg...info...           NG:: ndelta(j) , integer(8) :: v_ssite%ndelta 1
C$$$Cgetarg           call pack2('site ndelta',v_ssite,j,ndelta(j))

C$$$           v_ssite(j)%ndelta=ndelta(j)

C$$$Cgetarg...info...           NG:: delta(1,j) , real(8) :: v_ssite%delta 6
C$$$Cgetarg           call pack2('site delta',v_ssite,j,delta(1,j))

C$$$           i_copy_size=size(v_ssite(j)%delta)
C$$$           call dcopy(i_copy_size,delta(1,j),1,v_ssite(j)%delta,1)

C$$$        endif
      enddo

C ... Copy string outputs to sstrn : amix, gfopt, jobid, mix, mmham, sxopt, symg
      j = len_trim(iter_amix)          ! Euler angle mixing amix
      if (j .gt. 0) then
c        call lstra('strn amix',i,o,k)  ! o = index amix has (see ustrn),
c        call ustrn(w,-o,1,i,k,j)       ! Copy to sstrn(i:i+j-1)
c        sstrn(i:i+j-1) = iter_amix
        call packs('strn amix',sstrn,iter_amix)
        call upacks('strn amix',inix,iendx)
c         print *,'vvvvv ',sstrn(inix:iendx)
c         stop 'xxxxxxxxxxxxxxx'
      endif
      j = len_trim(iter_mix)           ! density mixing mix
      if (j .gt. 0) then
c        call lstra('strn mix',i,o,k)   ! o = index mix has (see ustrn),
c        call ustrn(w,-o,1,i,k,j)       ! Copy to sstrn(i:i+j-1)
c        sstrn(i:i+j-1) = iter_mix
        call packs('strn mix',sstrn,iter_mix)
      endif
c      j = len_trim(mmham)              ! micromagnetics hamiltonian
c      if (j .gt. 0) then
c        call lstra('strn mmham',i,o,k) ! o = index mmham has (see ustrn),
c        call ustrn(w,-o,1,i,k,j)       ! string should be copied to sstrn(i:j)
c        sstrn(i:i+j-1) = mmham
c         call packs('strn mmham',sstrn,mmham)
c      endif
      j = len_trim(header)             ! jobid
      if (j .gt. 0) then
c        call lstra('strn jobid',i,o,k) ! o = index jobid has (see ustrn),
c        call ustrn(w,-o,1,i,k,j)       ! Copy to sstrn(i:i+j-1)
c        sstrn(i:i+j-1) = header
        call packs('strn jobid',sstrn,header)
      endif

c      j = len_trim(gfopt)              ! Green's function gfopt
c      if (j .gt. 0) then
c        call lstra('strn gfopt',i,o,k) ! o = index mix has (see ustrn),
c        call ustrn(w,-o,1,i,k,j)       ! Copy to sstrn(i:i+j-1)
c        sstrn(i:i+j-1) = gfopt
c        call packs('strn gfopt',sstrn,gfopt)
c      endif

c      j = len_trim(sxopt)              ! Screened exchange sxopt
c      if (j .gt. 0) then
c        call lstra('strn sxopt',i,o,k) ! o = index mix has (see ustrn),
c        call ustrn(w,-o,1,i,k,j)       ! Copy to sstrn(i:i+j-1)
c        sstrn(i:i+j-1) = sxopt
c         call packs('strn sxopt',sstrn,sxopt)
c      endif

C ... Suppress symmetry operations for special circumstances
      lstsym = 0
      if (lncol .ne. 0) lstsym=1 !lstsym=1: noncollinear case
      !      =2: turn off symops
C     Switches that automatically turn of all symops
      if ((mdprm(1) .ge. 1 .and. mdprm(1) .le. 3) .or.   
     .cmdopt('--cls',5,0,strn) .or. cmdopt('--nosym',7,0,strn)) then
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
c        call lstra('strn symg',i,o,k)  ! o = index symg has (see ustrn),
c        call ustrn(w,-o,1,i,k,j)       ! Copy to sstrn(i:i+j-1)
c        sstrn(i:i+j-1) = symg
        call packs('strn symg',sstrn,symg)
      endif
C ... End of copy

C$$$C --- Miscellaneous cleanup and initialization ---
C$$$Cgetarg       if (lgors('ctrl lgen3',v_sctrl)) then
C$$$       if ( iand(0,int(v_sctrl%lgen3)) .ne.0 ) then
C$$$Cgetarg         nkap = dgets('str skmsh:1.1',v_sstr)
C$$$         nkap = v_sstr%skmsh(1)
C$$$Cgetarg         xx(1) = dgets('str skmsh:2.2',v_sstr)
C$$$         xx ( 1 ) = v_sstr%skmsh(2)
C$$$Cgetarg         xx(2) = dgets('str skmsh:3.3',v_sstr)
C$$$         xx ( 2 ) = v_sstr%skmsh(3)
C$$$        call gausq(nkap,xx(1),xx(2),ekap,xx(3),0,0)
C$$$      endif

C ... Read the basis from the atm file
      if (cmdopt('--rdbasp',8,0,strn)) then
         call rx('not support --rdbasp')
C$$$        fileid = 'basp'
C$$$        if (strn(9:12) .eq. ':fn=') then
C$$$          fileid = strn(13:13+63) !tk to avoid complaint; replace fileid = strn(13:)
C$$$        else
C$$$        endif

C$$$C       Number of envelope functions before modifications
C$$$        call uspecb(0,1,v_sspec,1,nspec,xx,xx,xx,k1)

C$$$C       Read new parms; return in k2 max no. env. fns/site read
C$$$        call strip(fileid,i,j)
C$$$        ifi = fopna(fileid(1:j),-1,0)
C$$$        rewind ifi
C$$$        if (.not. ioorbp(111,2,1,nspec,v_sspec,k2,ifi))
C$$$     .  call rxs2('lmfp: failed to find BASIS: token in file "',
C$$$     .  fileid(1:j),'"')
C$$$        k2 = mod(k2,10)
C$$$        call fclr(' ',ifi)

C$$$C       File read cause number of envelope functions to increase?
C$$$C       If so, increment nkaph by 1; keep old nkaph in k1 for later
C$$$        if (k2 .gt. k1) then
C$$$          k1 = nglob('nkaph')
C$$$          xx(1) = dglob('nkaph',dble(k1+1),1)
C$$$          k = nglob('nkaph')
C$$$          call uspecb(0,-1,v_sspec,1,nspec,xx,xx,xx,xx)
C$$$          nlmax = int(v_sctrl%nl) * * 2

C$$$          xx(1) = dglob('mxorb',dble(k)*nlmax,1)
C$$$        endif
      endif

C     Add dalat to alat

CXX  v_slat%alat=(v_slat%alat)+dalat !this is a bug; this should be commented out. This 
CXX  is a bug for lm7K when I started lm7K. I had included this bug here. Fixed at 28May2010.

C     Replicate ctrl->lncol in ham->lncol

      v_sham%lncol=lncol

c      lsx1 = mod(lsx,2)
C$$$       if ( lasa32 ) then

C$$$       v_smix%nsave=3
C$$$       endif

C ... Use true spherical harmonics

c       if ( bittst ( lncol , 4 ) .or. ctrl_loptc .ne. 0 .or. lrel .eq.
c     . 2 ) then
c       if ( bittst ( lncol , 4 ) .or. lrel .eq. 2 ) then
c       call lsets_bitop_i8(v_sctrl%lham,1,.true.,256)
c       endif


C ... Dirac equation requires spin polarization
      if ( nsp.eq.1 .and. int(v_sctrl%lrel) .eq.2 ) then

        call rx('rdccat: Dirac equation requires NSPIN=2')
      endif

C ... Suppress inversion when noncollinear magnetism, SX, NLO
c      if (lncol+lsx1 .ne. 0 .or. ctrl_loptc .ge. 10 .or.
c     .   (mod(lscr,10) .eq. 1 .and. prgnam .eq. 'LM'))
      if (lncol .ne. 0)
     .lqp = lqp-bitand(lqp,1)+1

      v_sctrl%lqp=lqp


C ... Special pgf initialization
C$$$      if (ctrl_lpgf(1) .ne. 0) then
C$$$Cgetarg         call lsets('ctrl lmet',v_sctrl,.false.,2)

C$$$         call lsets_bitop_i8(v_sctrl%lmet,1,.false.,2)

C$$$        call redfrr(v_ssite,3*nsite*nint(dval(v_ssite,1)))
C$$$      endif

C ... Setup for idxdn ... ctrl->lham,4 is automatic downfolding switch
c     j = iand(4,int(v_sctrl%lham)) /4
c     j = 2*(1-j)
      j= 2

C     No screening => no downfolding; also lmxb<l<=lmxa => 'high'
C     Probably ought to have lmxb<l<=lmxa => 'high' always
      if ( iand(1,int(v_sctrl%lbas)) .ne.0 ) j = 3

C     nfp-style basis:
      if ( iand(2,int(v_sctrl%lbas)) .ne.0 ) j = j+10

Changenglob      call suidx(nglob('nkaph'),j,nspec,v_sspec)
      call suidx(globalvariables%nkaph,j,nspec,v_sspec)

C ... Set some global variables
Changenglob      xx(1) = dglob('nspec',dble(nspec),1)
      globalvariables%nspec = dble(nspec); globalvariables%l_nspec = globalvariables%l_nspec +1; xx(1) = dble(nspec)
Changenglob      xx(1) = dglob('nbas',dble(nbas),1)
      globalvariables%nbas = dble(nbas); globalvariables%l_nbas = globalvariables%l_nbas +1; xx(1) = dble(nbas)
Changenglob      xx(1) = dglob('nbasp',dble(nbasp),1)
      globalvariables%nbasp = dble(nbasp); globalvariables%l_nbasp = globalvariables%l_nbasp +1; xx(1) = dble(nbasp)
Changenglob      xx(1) = dglob('nsp',dble(nsp),1)
      globalvariables%nsp = dble(nsp); globalvariables%l_nsp = globalvariables%l_nsp +1; xx(1) = dble(nsp)
Changenglob      xx(1) = dglob('nl',dble(nl),1)
      globalvariables%nl = dble(nl); globalvariables%l_nl = globalvariables%l_nl +1; xx(1) = dble(nl)
Changenglob      xx ( 1 ) = dglob ( 'avw' , v_slat%avw , 1 )
      globalvariables%lrel = dble(avw); globalvariables%l_avw = globalvariables%l_avw +1; xx(1) = dble(avw)

Changenglob      xx(1) = dglob('lrel',dble(lrel),1)
      globalvariables%lrel = dble(lrel); globalvariables%l_lrel = globalvariables%l_lrel +1; xx(1) = dble(lrel)
      lxcf = int(v_sctrl%lxcf)

Changenglob      xx(1) = dglob('lxcf',dble(lxcf),1)
      globalvariables%lxcf = dble(lxcf); globalvariables%l_lxcf = globalvariables%l_lxcf +1; xx(1) = dble(lxcf)
Changenglob      xx(1) = dglob('stdo',dble(stdo),1)
      globalvariables%stdo = dble(stdo); globalvariables%l_stdo = globalvariables%l_stdo +1; xx(1) = dble(stdo)
Changenglob      xx(1) = dglob('stdl',dble(stdl),1)
      globalvariables%stdl = dble(stdl); globalvariables%l_stdl = globalvariables%l_stdl +1; xx(1) = dble(stdl)
Changenglob      xx(1) = dglob('stde',dble(stde),1)
      globalvariables%stde = dble(stde); globalvariables%l_stde = globalvariables%l_stde +1; xx(1) = dble(stde)
      nspc = 1
c      if (bitand(lncol,1+2+4) .ne. 0) nspc = 2
      if (iand(lncol,4) .ne. 0) nspc = 2
Changenglob      xx(1) = dglob('nspc',dble(nspc),1)
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
Changenglob      xx(1) = dglob('nat',dble(nat),1)
      globalvariables%nat = dble(nat); globalvariables%l_nat = globalvariables%l_nat +1; xx(1) = dble(nat)

C ... Set modep
      ix(1) = 2
      ix(2) = 2
      ix(3) = 2
C$$$       if ( iand(-1,int(v_sctrl%lpgf(1))) .ne.0 ) ix ( 3 ) = 0


      i_copy_size=size(v_sctrl%modep)
      call ii8copy(i_copy_size,ix,1,v_sctrl%modep,1)


c      v_sham%lxcf=lxcf


C ... Count LDA+U blocks (printout only)
c      call defi(owk,-nbas)
      allocate(wowk(nbas))
      wowk=0
      call pshpr(0)
      call suldau(nbas,v_sspec,v_ssite,nlibu,k,wowk)
      v_sham%nlibu=nlibu
      v_sham%lmaxu=k
      call poppr
c      call rlse(owk)
      deallocate(wowk)

C     Free arrays used to read input
      deallocate(pnu,qnu,pz,amom,idmod,rsmh,eh,rsmh2,eh2,pb1,pb2,
     .lmxpb,qpol,stni,tbvso,iq1,ivso,rg,rsma,rfoca,rsmfa,rcfa,nxi,
     .exi,rint,rcut,coreq,mass,colxbs,radxbs,rs3,rham,idxdn,hcr,
     .rmt,alpha,idu,uh,jh,dv,grp,grp2,mxcst1,mxcst2,mxcst4,kmxt,kmxv,
     .lfoca,eref,lmxl,lxi,coreh,lmxa,lmxb,spec_a,z,nr,rsmv)
      deallocate(pos,vel,eula,vshft,ips,ipl,plv,irlx,mpole,dpole)
c      if (ltbe) deallocate(delta,ndelta)

C      call rx0('done copying input to structures')

C --- Printout ---
      if (iprint() .ge. 20 .and. procid .eq. master) then
        do  90  k = 1, 2

          strn = '  '//prgnam
C$$$        if (prgnam .eq. 'LMMC') then
C$$$Cgetarg           call awrit4(
C$$$Cgetarg      .    '%N %a:%12palat = %;5d  '//
C$$$Cgetarg      .    'nbas = %i%?#n#%-1j(+%i)##'//
C$$$Cgetarg      .    '  nspec = %i',strn,scrwid,0,dgets('lat alat',v_slat),
C$$$Cgetarg      .      nbas,nsite-nbas,nspec)
C$$$           call awrit4 ( '%N %a:%12palat = %;5d  ' // 'nbas = %i%?#n#%-1j(+%i)##'
C$$$     .     // '  nspec = %i' , strn , scrwid , 0 , v_slat%alat , nbas ,
C$$$     .     nsite-nbas , nspec )

C$$$        else
C$$$Cgetarg           call awrit6(
C$$$Cgetarg      .    '%N %a:%12palat = %;5d  '//
C$$$Cgetarg      .    'nbas = %i%?#n#%-1j(+%i)##%?#n#+%i#%j#'//
C$$$Cgetarg      .    '  nspec = %i',strn,scrwid,0,dgets('lat alat',v_slat),
C$$$Cgetarg      .    nat,nbas-nat,nbas-nsite,nsite,nspec)
C$$$           call awrit6 ( '%N %a:%12palat = %;5d  ' // 'nbas = %i%?#n#%-1j(+%i)##%?#n#+%i#%j#'
C$$$     .     // '  nspec = %i' , strn , scrwid , 0 , v_slat%alat , nat , nbas-nat
C$$$     .     , nbas-nsite , nsite , nspec )

C$$$        endif
          if (vn(2) .eq. 0) then
            call awrit6('%a  vn %,2d  verb %i%?#n#,%i#%j#%?#n#,%i#%j#',
     .      strn,scrwid,-lgunit(k),vn(1)-vn2(1)/1d5,
     .      iprint(),iprint()-iprt(1),iprt(1),iprt(2)-iprt(1),iprt(2))
          else
            call awrit7('%a  vn %,2d('//trim(prgnam)//' %,1d)  '//
     .      'verb %i%?#n#,%i#%j#%?#n#,%i#%j#',
     .      strn,scrwid,-lgunit(k),vn(1)-vn2(1)/1d5,vn(2),
     .      iprint(),iprint()-iprt(1),iprt(1),iprt(2)-iprt(1),iprt(2))
          endif


          lxcf=v_sctrl%lxcf

c         i = iand(256,int(v_sctrl%lham))

          lfrzw = isw ( iand(16,int(v_sctrl%lbas)) .ne.0 )


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
c     .    '%b %b',strn,scrwid,0,ctrl_lfrce,lrel,i,0, !ctrl_loptc,
c     .    mod(lscr,10),lordn,pwmode,isw(oveps.ne.0))
     .    '%b %b',strn,scrwid,0,ctrl_lfrce,lrel,0,0, !ctrl_loptc,
     .    mod(lscr,10),0,pwmode,isw(oveps.ne.0))
          if (strn .ne. ' special:') call awrit0(strn,' ',-80,lgunit(k))

C$$$          ltmp = prgnam.eq.'LM'  .or. prgnam.eq.'LMMC'  .or.
C$$$     .    prgnam.eq.'LMF' .or. prgnam.eq.'LMFGWD' .or.
C$$$     .    prgnam.eq.'LMGF' .or. prgnam.eq.'LMPG' .or.
C$$$     .    prgnam.eq.'LMFA' .or. prgnam.eq.'LMCHK'
C$$$          if (ltmp) then
C$$$            call awrit1('%x pot:%10p'//'%?;n==0; non-rel,;.;',strn,
C$$$     .      scrwid,0,lrel)
C$$$            call awrit8('%a%?;%c==,;;%b%a%b;%-1j'//
C$$$     .      '%?;n>1; spin-pol,;;'//
C$$$     .      '%?;n; LDA+U,;;'//
C$$$     .      '%?;n==1; XC:CA,;;%-1j'//
C$$$     .      '%?;n==2; XC:BH,;;%-1j'//
C$$$     .      '%?;n==3; XC:PW91,;;%-1j'//
C$$$     .      '%?;n==4; XC:PBE,;;'//
C$$$     .      '%?;n==1;%b+LMH(gga),;;%-1j'//
C$$$     .      '%?;n==2;%b+PW91(gga),;;%-1j'//
C$$$     .      '%?;n==3;%b+PBE(gga),;;%-1j'//
C$$$     .      '%?;n==4;%b+Becke(gga),;;'//
C$$$     .      '%?;n; frozen-wave-functions,;;'//
C$$$     .      '%?;n==1; nsph-mpol,;;'//
C$$$     .      '%?;n; read Sigma,;;'//
C$$$     .      '%?;n; make SX,;;'//
C$$$     .      '%b %b',
C$$$clm72
C$$$     .      strn,scrwid,0,nsp,nlibu,mod(lxcf,100),lxcf/100,lfrzw, !lm72
C$$$c     .      strn,scrwid,0,nsp,nlibu,mod(lxcf,10),lxcf/10,lfrzw,
C$$$c     .      isw(lasa32),lrsig,lsx)
C$$$     .      0,lrsig,0)
C$$$            call awrit0(strn,' ',-80,lgunit(k))
C$$$          endif

c          call isanrg(mod(lxcf,10),0,4,prgnam,'XC functional',.true.)
c          call isanrg(lxcf/10,0,4,prgnam,'GGA functional',.true.)
c          if (mod(lxcf,10) .eq. 2 .and. lxcf/10 .gt. 1 .or.
c     .    mod(lxcf,10) .ne. 2 .and. lxcf/10 .eq. 1) call info0
c     .    (10,0,0,'%10f(warning) mixing incompatible functionals')
clm72
      if (prgnam .eq. 'LMF' .or. prgnam .eq. 'LMFA') then
        if(lxcf/=1.and.lxcf/=2.and.lxcf/=103) then
          call rx('only XCFUN=1,2,103 are allowed now')
        endif
      endif  
C$$$        call isanrg(mod(lxcf,100),0,3,prgnam,'XC functional',.true.)
C$$$        call isanrg(lxcf/100,0,1,prgnam,'GGA functional',.true.)
C$$$        if (mod(lxcf,100) .eq. 2 .and. lxcf/100 .gt. 1 .or.
C$$$     .      mod(lxcf,100) .ne. 2 .and. lxcf/100 .eq. 1) call info0
C$$$     .    (10,0,0,'%10f(warning) mixing incompatible functionals')

C$$$        if (bitand(lncol,1+2+4+8+32+64) .ne. 0) then
C$$$Cgetarg           i = igets('ctrl sdmod',v_sctrl)
C$$$           i = int(v_sctrl%sdmod)

C$$$          if (.not. bittst(lncol,16)) i = -1
C$$$          call awrit8(' noncoll: '//
C$$$     .    '%?;n; Non-coll,;;'//
C$$$     .    '%?;n; Spin-orbit,;;'//
C$$$     .    '%?;n; LzSz,;;%-1j%?;(n>32);%b+L.S(pert),;;'//
C$$$     .    '%?;n; B-field,;;'//
C$$$     .    '%?;n; spin-spiral,;;'//
C$$$     .    '%?;n; mag-forces:;;'//
C$$$     .    '%?;(n%10)<2&(n%10)>=0; relax,;;'//
C$$$     .    '%?;(n%10)>1; spin-dynamics,;;'//
C$$$     .    '%a%b %a',strn,scrwid,lgunit(k),
C$$$     .      bitand(lncol,1),bitand(lncol,4),bitand(lncol,32+64),
C$$$     .      bitand(lncol,8),bitand(lncol,2),bitand(lncol,16),i,i)
C$$$        endif

C$$$        if (asa) then
C$$$Cgetarg...info...           integer :: lham , integer(8) :: v_sctrl%lham 1
C$$$Cgetarg...info...           integer :: lcd , integer(8) :: v_sctrl%lcd 1
C$$$Cgetarg           call upack('ctrl lham lcd',v_sctrl,lham,lcd,0,0,0)

C$$$           lham=v_sctrl%lham
C$$$           lcd=v_sctrl%lcd

C$$$Cgetarg           call awrit8('%x asa:%10p'//
C$$$Cgetarg      .      '%?;n; no-ccor,;;'//
C$$$Cgetarg      .      '%?;n; gam-rep,;;'//
C$$$Cgetarg      .      '%?;n; two-c-H,;;'//
C$$$Cgetarg      .      '%?;n;%b + pert-ev,;;'//
C$$$Cgetarg      .      '%?;n; map,;;'//
C$$$Cgetarg      .      '%?;n; indep-vmix,;;'//
C$$$Cgetarg      .      '%?;n; mt-corr,;;'//
C$$$Cgetarg      .      '%?;n; frozen core,;;'//
C$$$Cgetarg      .      '%a%b ',strn,scrwid,0,
C$$$Cgetarg      .      bitand(lasa,4).eq.0,
C$$$Cgetarg      .      bitand(lham,128),
C$$$Cgetarg      .      bitand(lham,3),bitand(lham,2),
C$$$Cgetarg      .      bitand(lasa,16),igets('mix lxpot,3',v_smix),
C$$$Cgetarg      .      bitand(lasa,64),bitand(lcd,1))
C$$$           call awrit8 ( '%x asa:%10p' // '%?;n; no-ccor,;;' // '%?;n; gam-rep,;;'
C$$$     .     // '%?;n; two-c-H,;;' // '%?;n;%b + pert-ev,;;' // '%?;n; map,;;'
C$$$     .     // '%?;n; indep-vmix,;;' // '%?;n; mt-corr,;;' // '%?;n; frozen core,;;'
C$$$     .     // '%a%b ' , strn , scrwid , 0 , bitand ( lasa , 4 ) .eq.0 ,
C$$$     .     bitand ( lham , 128 ) , bitand ( lham , 3 ) , bitand ( lham ,
C$$$     .     2 ) , bitand ( lasa , 16 ) , iand(3,int(v_smix%lxpot)) , bitand
C$$$     .     ( lasa , 64 ) , bitand ( lcd , 1 ) )

C$$$          if (strn .ne. ' asa') call awrit0(strn,' ',-80,lgunit(k))
C$$$        endif

C$$$        if (prgnam(1:2) .eq. 'TB') then
C$$$Cgetarg...info...           real(8) :: xx(n0*2) , real(8) :: v_sctrl%mdprm 6
C$$$Cgetarg...info...           integer :: ltb , integer(8) :: v_sctrl%ltb 1
C$$$Cgetarg           call upack2('ctrl mdprm ltb',v_sctrl,xx,ltb)

C$$$           i_copy_size=size(v_sctrl%mdprm)
C$$$           call dcopy(i_copy_size,v_sctrl%mdprm,1,xx,1)
C$$$c           ltb=v_sctrl%ltb

C$$$          ix(1) = nint(xx(1))
C$$$          ix(2) = nint(xx(2))
C$$$Cgetarg...info...           real(8) :: xx(n0*2) , real(8) :: v_sstr%rmax 1
C$$$Cgetarg           call upack1('str rmax',v_sstr,xx)

C$$$           xx=v_sstr%rmax

C$$$          call awrit1('%x TB: %11prmaxh = %d,',strn,scrwid,0,xx)
C$$$          call awrit6('%a'//
C$$$     .      '%?;n; s-c multipoles,;;'//
C$$$     .      '%?;n; s-c multipoles: MRS theory,;;'//
C$$$     .      '%?;n; read del,;;'//
C$$$     .      '%?;n; n-orthog TB,;;'//
C$$$     .      '%?;n==4; m-stat: Conj. grad.;;%-1j'//
C$$$     .      '%?;n==5; m-stat: F-P;;%-1j'//
C$$$     .      '%?;n==6; m-stat: Broy;;%-1j'//
C$$$     .      '%?;n==1; MD (NVE),;;%-1j'//
C$$$     .      '%?;n==2; MD (NVT),;;%-1j'//
C$$$     .      '%?;n==3; MD (NPT),;;'//
C$$$     .      '%?;n; pair-only,;;',
C$$$     .      strn,scrwid,0,bitand(ltb,2**15),bitand(ltb,2**12),
C$$$     .      bitand(ltb,2**16),bitand(ltb,1),ix(1),bitand(ltb,512))
C$$$          if (ix(1) .ge. 4) then
C$$$            call awrit2('%a'//
C$$$     .      '%?;n; rlx-vol,;;'//
C$$$     .      '%?;n; i/o-hess,;;',
C$$$     .      strn,scrwid,0,bitand(ix(2),1),bitand(ix(2),2))
C$$$          endif
C$$$          call awrit5('%a'//
C$$$     .      '%?;n; trh,;;'//
C$$$     .      '%?;n; rho,;;'//
C$$$     .      '%?;n; spin-orb,;;'//
C$$$     .      '%?;n; crysf,;;'//
C$$$     .      '%?;n;%b+ovlp-cf,;;'//
C$$$     .      '%a%b %a',strn,-scrwid,-lgunit(k),
C$$$     .      bitand(ltb,1024),bitand(ltb,2048),
C$$$     .      bitand(lncol,4),bitand(ltb,2),bitand(ltb,4))
C$$$        endif

          if (prgnam(1:3) .eq. 'LMF') then
            call awrit1 ( ' special:%10p' // '%?;n==1; core-level-optics,;;' 
     .      // '%b %b' , strn , scrwid , 0 , 0 )

            if (strn .ne. ' special:') call awrit0(strn,' ',-80,lgunit(k))
          endif

c       ltmp = prgnam.eq.'LM'   .or. prgnam.eq.'LMMC'  .or.
          ltmp = prgnam.eq.'LM'   .or. 
     .    prgnam.eq.'LMF'  .or. prgnam.eq.'LMFGWD' .or.
     .    prgnam.eq.'LMGF' .or. prgnam.eq.'LMPG' .or. 
     .    prgnam.eq.'LMDOS'.or. prgnam.eq.'TBE'

          if ( ltmp ) call awrit7 ( ' bz:%10p' // '%?;n; metal,; nonmetal,;' 
     .    // '%-1j%?;n>1;%b(%-1j%i),;;' // '%?;n; tetra,;;' // '%?;n; get-qp,;;' 
     .    // '%?;n; invit,;;' // '%?;n; dens-mat,;;' // '%?;(n>0); %-1jmull=%i,;;' 
     .    // '%?;n; fixed-spin-mom,;;%b ' , strn , scrwid , lgunit ( k 
     .    ) , int(v_sbz%lmet) , bitand ( lmet , 2 ) , iand(1,int(v_sbz%lio))
c    .   , iand(2,int(v_sctrl%lqp)) , iand(8,int(v_sbz%lio)) , int(v_sbz%lmull) 
     .    , iand(2,int(v_sctrl%lqp)) , .false.                , int(v_sbz%lmull) 
     .    , isw ( v_sbz%fsmom .ne.0 ) )


   90   continue
      endif

C --- Sanity checks and other initialization ---
      if ( iand(8,int(v_sctrl%ldos)) .ne.0 .and. iand(4+2,int(v_sctrl%ldos)) .ne.0 
     .) call rx ( 'inconsistent DOS options' )

C$$$Cgetarg       call rxx(lgors('ctrl lstonr,-1',v_sctrl) .and. nsp .eq. 2,
C$$$Cgetarg      .  'Stoner model not compatible with nsp=2')
C$$$       call rxx ( iand(-1,int(v_sctrl%lstonr(1))) .ne.0 .and. nsp .eq.
C$$$     . 2 , 'Stoner model not compatible with nsp=2' )

C$$$Cgetarg       call rxx(lgors('ctrl lstonr,-1',v_sctrl) .and.
C$$$Cgetarg      . lgors('ctrl lham,4',v_sctrl),'Stoner not compatible with ADNF')
C$$$       call rxx ( iand(-1,int(v_sctrl%lstonr(1))) .ne.0 .and. iand(4,int(v_sctrl%lham)) .ne.0
C$$$     . , 'Stoner not compatible with ADNF' )


C --- takao embed contents in susite here. This is only for lmf and lmfgw.
c      call susite(v_sctrl,v_sham,v_spot,v_sarry,v_slat,v_sspec,v_ssite)
      seref = 0d0
      do ib = 1, nsite
        ispec = v_ssite(ib)%spec
        if (ispec <= 0) cycle
c         call spacks_copy('u',sspec(i)%name,i,i,spid,ispec)
c         call spacks_copy('p',v_ssite(i)%clabel,ib,ib,spid,ispec)
        v_ssite(ispec)%clabel = v_sspec(ispec)%name
        seref = seref +     v_sspec(ispec)%eref
      enddo
      v_sham%seref=seref
c
       allocate(rv_p_opos(abs(3*nsite)))
       if (3*nsite<0) rv_p_opos(:)=0.0d0

       allocate(iv_p_oips(abs(nsite)))
       if (nsite<0) iv_p_oips(:)=0

      do i_spackv=1,nsite
        i_copy_size=3
         call spackv_array_copy_r8_r8 ( 'u' , v_ssite ( i_spackv ) %pos 
     .   , i_copy_size , i_spackv + 1 - 1 , rv_p_opos ) 

        i_copy_size=1
         call spackv_array_copy_i8_i ( 'u' , v_ssite ( i_spackv ) %spec 
     .   , i_copy_size , i_spackv + 1 - 1 , iv_p_oips ) 

      enddo
       v_slat%rv_p_opos => rv_p_opos 

       v_sarry%rv_p_opos => rv_p_opos 

       v_sarry%iv_p_oips => iv_p_oips 



c       v_sham%lham=(int(v_sctrl%lham))

      if (procid .eq. master) then
        if (iprint() .ge. 20) then
          if (lstsym .eq. 1) then
            call upacks('strn symg',i,j)
            write(stdo,357) sstrn(i:j)
  357       format(/' Automatic symmetry finder turned off.  Use: ',a)
          elseif (lstsym .eq. 2) then
            write(stdo,358)
  358       format(/' Symmetry operations suppressed')
          endif
        endif
      endif
C     Broadcast spec structure to fix floating orbitals case
cki      call mpibc1(v_sspec,nspec*nint(dval(v_sspec,1)),4,0,
       do i_spec=1,nspec
       call mpibc1_s_spec(v_sspec(i_spec),.false.,
     .'rdctrl','sspec')
       enddo
       

C --- Debugging printout ---
      if (io_help == 0 .and. io_show > 1) then
        print *, '---------- contents of sstrn ------------'
        call upacks('strn amix',i,j)
        print *, 'amix:', sstrn(i:j)
        call upacks('strn mix',i,j)
        print *, 'mix:', sstrn(i:j)
c        call upacks('strn gfopt',i,j)
c        print *, 'gfopt:', sstrn(i:j)
c        call upacks('strn mmham',i,j)
c        print *, 'mmham:', sstrn(i:j)
        call upacks('strn symg',i,j)
        print *, 'symg:', sstrn(i:j)
c        call upacks('strn sxopt',i,j)
c        print *, 'sxopt:', sstrn(i:j)

C$$$        call shstru('ctrl',v_sctrl,0,0)
C$$$        call shstru('bz',v_sbz,0,0)
C$$$        call shstru('lat',v_slat,0,0)
C$$$        call shstru('array',v_sarry,0,0)
C$$$C       call shstru('pot',v_spot,0,0)
C$$$        call shstru('ham',v_sham,0,0)
C$$$c        if (smix(2) /= NULLI)
C$$$        if (vmix%b /= NULLI)
C$$$     .  call shstru('mix',v_smix,0,0)
C$$$        call shstru('move',v_smove,0,0)
C$$$Cki        call shstru('optic',w(osoptc),0,0)
C$$$Cki        call shstru('gw',w(osgw),0,0)
C$$$        call shstru('str',v_sstr,0,0)
C$$$        call shstru('tb',v_stb,0,0)
C$$$        call shstru('spec',v_sspec,1,nspec)
C$$$        call shstru('site',v_ssite,1,nsite)
        call rx0('done show')

      endif
c      call rx0('done')
      end subroutine rdctrl2

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
      end module m_rdctrl2_func
