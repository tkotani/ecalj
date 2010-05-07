Cgetarg...info...           structure ['rdctrl2', 'bz', 'v_sbz']
Cgetarg...info...           structure ['rdctrl2', 'ctrl', 'v_sctrl']
Cgetarg...info...           structure ['rdctrl2', 'lat', 'v_slat']
Cgetarg...info...           structure ['rdctrl2', 'ham', 'v_sham']
Cgetarg...info...           structure ['rdctrl2', 'pot', 'v_spot']
Cgetarg...info...           structure ['rdctrl2', 'array', 'v_sarry']
Cgetarg...info...           structure ['rdctrl2', 'move', 'v_smove']
Cgetarg...info...           structure ['rdctrl2', 'optic', 'w(osoptc)']
Cgetarg...info...           structure ['rdctrl2', 'gw', 'w(osgw)']
Cgetarg...info...           structure ['rdctrl2', 'str', 'v_sstr']
Cgetarg...info...           structure ['rdctrl2', 'tb', 'v_stb']
Cgetarg...info...           structure ['rdctrl2', 'spec', 'v_sspec']
Cgetarg...info...           structure ['rdctrl2', 'site', 'v_ssite']
Cgetarg...info...           structure ['rdctrl2', 'mix', 'v_smix']
      module m_rdctrl2_func
      contains
      subroutine rdctrl2(recrd,recln,nrecs,
     . prgnam,vrsion,vn,vn2,
     . pass2,slabl_,v_sbz,v_sctrl,v_sham,v_spot,v_slat,v_smix,v_sspec,v_ssite,
     . v_sstr,v_sarry,v_smove,v_stb,sstrn)

      use m_rdctrl
      use m_gtv
c      use m_struc_size !   use m_struc_func
      use m_susite_func, only: susite
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
      character*(1000) recrd
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
     .  iprt,irs(5),isw,ifi,ix(n0*nkap0),j,k,l,lasa,lbas,lcd,lqp,
     .  lfrzw,lgunit,lmet,lncol, !,lham
     .  lrs,lstsym,ltb,lxcf,nat,nlibu,nclasp,
     .  nglob,nspc,nlmax,scrwid,stdo,stdl,stde,k1,k2,mpipid !nkap,lsx1,
        character*(8),allocatable::clabl(:)
        integer,allocatable:: ipc(:),initc(:),ics(:)
        real(8),allocatable:: pnuc(:,:,:),qnuc(:,:,:,:),
     .    pp(:,:,:,:),ves(:),zc(:)
      integer:: dvec1(3)=1, dvec2(3)=0

C ... basis
      double precision orbp(n0,2,nkap0)
      integer o,oclabl,ohave,oics,opnu,opp,oqnu,osgw,osoptc,
     .  oves,owk !osordn,
C ... Heap
      integer w(1)
      common /w/ w

c takao 
      real(8):: pnux(20)
      integer:: nnn
      integer:: umix_size,umove_size,uspec_size,usite_size,ubz_size,
     &  uctrl_size,ulat_size,uham_size,uarray_size,upot_size,ustr_size
c
      integer:: i_copy_size,i_spacks,iendx,inix,i_spackv

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

Cgetarg...info...           NG:: bz_def , real(8) :: v_sbz%def 1
Cgetarg...info...           NG:: bz_dosw , real(8) :: v_sbz%dosw 2
Cgetarg...info...           NG:: bz_ef , real(8) :: v_sbz%ef 1
Cgetarg...info...           NG:: bz_efmax , real(8) :: v_sbz%efmax 1
Cgetarg...info...           NG:: bz_fsmom , real(8) :: v_sbz%fsmom 1
Cgetarg       call pack5('bz def dosw ef efmax fsmom',v_sbz,
Cgetarg      .  bz_def,bz_dosw,bz_ef,bz_efmax,bz_fsmom)
       
       v_sbz%def=bz_def 
       i_copy_size=size(v_sbz%dosw) 
       call dcopy(i_copy_size,bz_dosw,1,v_sbz%dosw,1) 
       v_sbz%ef=bz_ef 
       v_sbz%efmax=bz_efmax 
       v_sbz%fsmom=bz_fsmom 

      i = 1*isw(bz_lio1)+2*isw(bz_lio2) !+8*isw(bz_lio8)
Cgetarg...info...           NG:: bz_lcond , real(8) :: v_sbz%lcond 4
Cgetarg...info...           integer :: i , integer(8) :: v_sbz%lio 1
Cgetarg...info...           NG:: bz_lmet , integer(8) :: v_sbz%lmet 1
Cgetarg...info...           NG:: bz_lmull , integer(8) :: v_sbz%lmull 1
Cgetarg       call pack5('bz lcond lio lmet lmull lopt',v_sbz,
Cgetarg      .  bz_lcond,i,bz_lmet,bz_lmull,0)
       
       i_copy_size=size(v_sbz%lcond) 
       call dcopy(i_copy_size,bz_lcond,1,v_sbz%lcond,1) 
       v_sbz%lio=i 
       v_sbz%lmet=bz_lmet 
       v_sbz%lmull=bz_lmull 
       v_sbz%lopt=0 

Cgetarg...info...           NG:: bz_lshft , integer(8) :: v_sbz%lshft 3
Cgetarg...info...           NG:: bz_n , integer(8) :: v_sbz%n 1
Cgetarg...info...           NG:: bz_ndos , integer(8) :: v_sbz%ndos 1
Cgetarg...info...           NG:: bz_nevmx , integer(8) :: v_sbz%nevmx 1
Cgetarg...info...           NG:: bz_nabc , integer(8) :: v_sbz%nkabc 3
Cgetarg       call pack5('bz lshft n ndos nevmx nkabc',v_sbz,
Cgetarg      .  bz_lshft,bz_n,bz_ndos,bz_nevmx,bz_nabc)
       
       i_copy_size=size(v_sbz%lshft) 
       call ii8copy(i_copy_size,bz_lshft,1,v_sbz%lshft,1) 
       v_sbz%n=bz_n 
       v_sbz%ndos=bz_ndos 
       v_sbz%nevmx=bz_nevmx 
       i_copy_size=size(v_sbz%nkabc) 
       call ii8copy(i_copy_size,bz_nabc,1,v_sbz%nkabc,1) 

Cgetarg...info...           NG:: bz_range , real(8) :: v_sbz%range 1
Cgetarg...info...           NG:: bz_semsh , real(8) :: v_sbz%semsh 10
Cgetarg...info...           NG:: bz_w , real(8) :: v_sbz%w 1
Cgetarg...info...           NG:: bz_zval , real(8) :: v_sbz%zval 1
Cgetarg       call pack5('bz range semsh w zval',v_sbz,
Cgetarg      .  bz_range,bz_semsh,bz_w,bz_zval,0)
       
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
Cgetarg...info...           NG:: lat_defm , real(8) :: v_sctrl%defm 6
Cgetarg...info...           NG:: asa_elin , real(8) :: v_sctrl%elin 1
Cgetarg...info...           integer :: lasa , integer(8) :: v_sctrl%lasa 1
Cgetarg...info...           integer :: j , integer(8) :: v_sctrl%lbas 1
Cgetarg...info...           integer :: k , integer(8) :: v_sctrl%lcd 1
Cgetarg       call pack5('ctrl defm elin lasa lbas lcd',v_sctrl,
Cgetarg      .  lat_defm,asa_elin,lasa,j,k)
       
       i_copy_size=size(v_sctrl%defm) 
       call dcopy(i_copy_size,lat_defm,1,v_sctrl%defm,1) 
C$$$       v_sctrl%elin=asa_elin 
C$$$       v_sctrl%lasa=lasa 
       v_sctrl%lbas=j 
       v_sctrl%lcd=k 

c      k = isw(nmto.gt.1)
c      if (nmto.gt.1) k = k+2*isw(ham_ewald)
Cgetarg...info...           NG:: ctrl_lcgf , integer(8) :: v_sctrl%lcgf 1
Cgetarg...info...           NG:: ctrl_ldos , integer(8) :: v_sctrl%ldos 1
Cgetarg...info...           NG:: lfp , integer(8) :: v_sctrl%lfp 1
Cgetarg...info...           NG:: ctrl_lfrce , integer(8) :: v_sctrl%lfrce 1
Cgetarg...info...           integer :: k , integer(8) :: v_sctrl%lgen3 1
Cgetarg       call pack5('ctrl lcgf ldos lfp lfrce lgen3',v_sctrl,
Cgetarg      . ctrl_lcgf,ctrl_ldos,lfp,ctrl_lfrce,k)
       
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
     .      +4*isw(ctrl_lmet4) + 8*isw(ctrl_lmet8)


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
Cgetarg...info...           integer :: lham , integer(8) :: v_sctrl%lham 1
Cgetarg...info...           integer :: lmet , integer(8) :: v_sctrl%lmet 1
Cgetarg...info...           integer :: lncol , integer(8) :: v_sctrl%lncol 1
Cgetarg...info...           NG:: ctrl_loptc , integer(8) :: v_sctrl%loptc 1
Cgetarg...info...           integer :: lordn , integer(8) :: v_sctrl%lordn 1
Cgetarg       call pack5('ctrl lham lmet lncol loptc lordn',v_sctrl,
Cgetarg      .  lham,lmet,lncol,ctrl_loptc,lordn)
       
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
     .       + 8*getdig(irs(1),1,100)

ccccccccccccccccccccccccccccccccccccccccccc
c      print *,'xxxxxxxxxx irs(1)=',irs(1)
ccccccccccccccccccccccccccccccccccccccccccc

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
Cgetarg...info...           NG:: ctrl_lpgf , integer(8) :: v_sctrl%lpgf 2
Cgetarg...info...           integer :: lqp , integer(8) :: v_sctrl%lqp 1
Cgetarg...info...           NG:: lrel , integer(8) :: v_sctrl%lrel 1
Cgetarg...info...           integer :: lrs , integer(8) :: v_sctrl%lrs 1
Cgetarg...info...           NG:: lscr , integer(8) :: v_sctrl%lscr 1
Cgetarg        call pack5('ctrl lpgf lqp lrel lrs lscr',v_sctrl,
Cgetarg      .  ctrl_lpgf,lqp,lrel,lrs,lscr)
        
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
Cgetarg...info...           NG:: lstonr , integer(8) :: v_sctrl%lstonr 3
Cgetarg...info...           NG:: lsx , integer(8) :: v_sctrl%lsx 1
Cgetarg...info...           integer :: i , integer(8) :: v_sctrl%ltb 1
Cgetarg...info...           NG:: lves , integer(8) :: v_sctrl%lves 1
Cgetarg        call pack5('ctrl lstonr lstr lsx ltb lves',v_sctrl,
Cgetarg      .   lstonr,0,lsx,i,lves)
        
C$$$        i_copy_size=size(v_sctrl%lstonr) 
C$$$        call ii8copy(i_copy_size,lstonr,1,v_sctrl%lstonr,1) 
        v_sctrl%lstr=0 
c        v_sctrl%lsx=lsx 
c        v_sctrl%ltb=i 
        v_sctrl%lves=lves 

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
Cgetarg...info...           integer :: i , integer(8) :: v_sctrl%lxcf 1
Cgetarg...info...           NG:: iter_maxit , integer(8) :: v_sctrl%maxit 1
Cgetarg...info...           NG:: mdprm , real(8) :: v_sctrl%mdprm 6
Cgetarg...info...           integer :: ix(n0*nkap0) , integer(8) :: v_sctrl%modep 3
Cgetarg...info...           NG:: nbas , integer(8) :: v_sctrl%nbas 1
Cgetarg        call pack5('ctrl lxcf maxit mdprm modep nbas',v_sctrl,
Cgetarg      .   i,iter_maxit,mdprm,ix(1:3),nbas)
        
        v_sctrl%lxcf=i 
        v_sctrl%maxit=iter_maxit 
        i_copy_size=size(v_sctrl%mdprm) 
        call dcopy(i_copy_size,mdprm,1,v_sctrl%mdprm,1) 
        i_copy_size=size(v_sctrl%modep) 
        call ii8copy(i_copy_size,ix(1:3),1,v_sctrl%modep,1) 
        v_sctrl%nbas=nbas 

C      Reset nl
       if (nl .ne. max(lmxbx,lmxax)+1 .and. io_help == 0) then
         call info2(20,1,0,' rdctrl: reset global max nl from %i to %i',
     .     nl,max(lmxbx,lmxax)+1)
         nl = max(lmxbx,lmxax)+1
       endif
C      call pack5('ctrl nbasp nclass nesabc nitmv nl',v_sctrl,
Cgetarg...info...           NG:: nbasp , integer(8) :: v_sctrl%nbasp 1
Cgetarg...info...           NG:: nesabc , integer(8) :: v_sctrl%nesabc 3
Cgetarg...info...           NG:: nitmv , integer(8) :: v_sctrl%nitmv 1
Cgetarg...info...           NG:: nl , integer(8) :: v_sctrl%nl 1
Cgetarg        call pack5('ctrl nbasp nesabc nitmv nl',v_sctrl,
Cgetarg      .   nbasp,nesabc,nitmv,nl,0)
        
        v_sctrl%nbasp=nbasp 
        i_copy_size=size(v_sctrl%nesabc) 
        call ii8copy(i_copy_size,nesabc,1,v_sctrl%nesabc,1) 
        v_sctrl%nitmv=nitmv 
        v_sctrl%nl=nl 

C      call pack5('ctrl nlibu nmap npl nsite nspec',v_sctrl,
Cgetarg...info...           NG:: nsite , integer(8) :: v_sctrl%nsite 1
Cgetarg...info...           NG:: nspec , integer(8) :: v_sctrl%nspec 1
Cgetarg        call pack5('ctrl nsite nspec',v_sctrl,
Cgetarg      .   nsite,nspec,0,0,0)
        
        v_sctrl%nsite=nsite 
        v_sctrl%nspec=nspec 

C      call shosyv(0,0,0,6)
Cgetarg...info...           NG:: nsp , integer(8) :: v_sctrl%nspin 1
Cgetarg...info...           NG:: nvario , integer(8) :: v_sctrl%nvario 1
Cgetarg...info...           NG:: omax1 , real(8) :: v_sctrl%omax1 3
Cgetarg...info...           NG:: omax2 , real(8) :: v_sctrl%omax2 3
Cgetarg...info...           NG:: quit , integer(8) :: v_sctrl%quit 1
Cgetarg        call pack5('ctrl nspin nvario omax1 omax2 quit',v_sctrl,
Cgetarg      .   nsp,nvario,omax1,omax2,quit)
        
        v_sctrl%nspin=nsp 
        v_sctrl%nvario=nvario 
        i_copy_size=size(v_sctrl%omax1) 
        call dcopy(i_copy_size,omax1,1,v_sctrl%omax1,1) 
        i_copy_size=size(v_sctrl%omax2) 
        call dcopy(i_copy_size,omax2,1,v_sctrl%omax2,1) 
        v_sctrl%quit=quit 

Cgetarg...info...           NG:: rmaxes , real(8) :: v_sctrl%rmaxes 1
Cgetarg...info...           NG:: rmines , real(8) :: v_sctrl%rmines 1
Cgetarg...info...           NG:: sclwsr , real(8) :: v_sctrl%sclwsr 1
Cgetarg...info...           NG:: sdmod , integer(8) :: v_sctrl%sdmod 1
Cgetarg...info...           NG:: sdprm , real(8) :: v_sctrl%sdprm 5
Cgetarg        call pack5('ctrl rmaxes rmines sclwsr sdmod sdprm',v_sctrl,
Cgetarg      .   rmaxes,rmines,sclwsr,sdmod,sdprm)
        
        v_sctrl%rmaxes=rmaxes 
        v_sctrl%rmines=rmines 
        v_sctrl%sclwsr=sclwsr 
        v_sctrl%sdmod=sdmod 
        i_copy_size=size(v_sctrl%sdprm) 
        call dcopy(i_copy_size,sdprm,1,v_sctrl%sdprm,1) 

Cgetarg...info...           NG:: smalit , integer(8) :: v_sctrl%smalit 2
Cgetarg...info...           NG:: ctrl_tol , real(8) :: v_sctrl%tol 3
Cgetarg...info...           NG:: wsrmax , real(8) :: v_sctrl%wsrmax 1
Cgetarg...info...           NG:: zbak , real(8) :: v_sctrl%zbak 2
Cgetarg...info...           NG:: lpfloat , integer(8) :: v_sctrl%pfloat 1
Cgetarg        call pack5('ctrl smalit tol wsrmax zbak pfloat',v_sctrl,
Cgetarg      .   smalit,ctrl_tol,wsrmax,zbak,lpfloat)
        
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
Cgetarg...info...           NG:: alat+dalat , real(8) :: v_slat%alat 1
Cgetarg...info...           NG:: lat_as , real(8) :: v_slat%as 1
Cgetarg...info...           NG:: avw , real(8) :: v_slat%avw 1
Cgetarg...info...           NG:: lat_nkdmx , integer(8) :: v_slat%nkdmx 1
Cgetarg...info...           NG:: lat_nkdmx , integer(8) :: v_slat%nkqmx 1
Cgetarg       call pack5('lat alat as avw nkdmx nkqmx',v_slat,
Cgetarg      .  alat+dalat,lat_as,avw,lat_nkdmx,lat_nkdmx)
       
       v_slat%alat=alat+dalat 
       v_slat%as=lat_as 
       v_slat%avw=avw 
       v_slat%nkdmx=lat_nkdmx 
       v_slat%nkqmx=lat_nkdmx 

Cgetarg...info...           NG:: lat_gam , real(8) :: v_slat%gam 4
Cgetarg...info...           NG:: lat_gmax , real(8) :: v_slat%gmax 1
Cgetarg...info...           NG:: ftmesh , integer(8) :: v_slat%nabc 3
Cgetarg       call pack5('lat gam gmax nabc',v_slat,
Cgetarg      .  lat_gam,lat_gmax,ftmesh,0,0)
       
       i_copy_size=size(v_slat%gam) 
       call dcopy(i_copy_size,lat_gam,1,v_slat%gam,1) 
       v_slat%gmax=lat_gmax 
       i_copy_size=size(v_slat%nabc) 
       call ii8copy(i_copy_size,ftmesh,1,v_slat%nabc,1) 

Cgetarg...info...           NG:: lat_gam , real(8) :: v_slat%gam 4
Cgetarg...info...           NG:: lat_ldist , integer(8) :: v_slat%ldist 1
Cgetarg...info...           NG:: lat_dist , real(8) :: v_slat%dist 9
Cgetarg       call pack5('lat gam ldist dist',v_slat,
Cgetarg      .  lat_gam,lat_ldist,lat_dist,0,0)
       
c takao I think gam is alreay setted above.
c       i_copy_size=size(v_slat%gam) 
c       call dcopy(i_copy_size,lat_gam,1,v_slat%gam,1) 

       v_slat%ldist=lat_ldist 
       i_copy_size=size(v_slat%dist) 
       call dcopy(i_copy_size,lat_dist,1,v_slat%dist,1) 

Cgetarg...info...           NG:: plat , real(8) :: v_slat%plat 9
Cgetarg...info...           NG:: slat_plat2 , real(8) :: v_slat%plat2 9
Cgetarg       call pack5('lat plat plat2',v_slat,plat,slat_plat2,0,0,0)
       
       i_copy_size=size(v_slat%plat) 
       call dcopy(i_copy_size,plat,1,v_slat%plat,1) 
c       i_copy_size=size(v_slat%plat2) 
c       call dcopy(i_copy_size,slat_plat2,1,v_slat%plat2,1) 

Cgetarg...info...           NG:: lat_rpad , real(8) :: v_slat%rpad 1
Cgetarg...info...           NG:: lat_slat , real(8) :: v_slat%slat 9
Cgetarg...info...           NG:: lat_tol , real(8) :: v_slat%tol 1
Cgetarg...info...           NG:: tolft , real(8) :: v_slat%tolft 1
Cgetarg...info...           NG:: vol , real(8) :: v_slat%vol 1
Cgetarg       call pack5('lat rpad slat tol tolft vol',v_slat,
Cgetarg      .  lat_rpad,lat_slat,lat_tol,tolft,vol)
       
       v_slat%rpad=lat_rpad 
       i_copy_size=size(v_slat%slat) 
       call dcopy(i_copy_size,lat_slat,1,v_slat%slat,1) 
       v_slat%tol=lat_tol 
       v_slat%tolft=tolft 
       v_slat%vol=vol 

c      call dcopy(6,plat,1,xx,1)
c      call dcopy(3,platl,1,xx(7),1)
Cgetarg...info...           real(8) :: xx(n0*2) , real(8) :: v_slat%platl 9
Cgetarg       call pack5('lat platl',v_slat,xx,0,0,0,0)
       
c       i_copy_size=size(v_slat%platl) 
c       call dcopy(i_copy_size,xx,1,v_slat%platl,1) 

c      call dcopy(3,platr,1,xx(7),1)
Cgetarg...info...           real(8) :: xx(n0*2) , real(8) :: v_slat%platr 9
Cgetarg       call pack5('lat platr',v_slat,xx,0,0,0,0)
       
c       i_copy_size=size(v_slat%platr) 
c       call dcopy(i_copy_size,xx,1,v_slat%platr,1) 


C --- Allocate and copy input to sham ---
CKi      call defdr(v_sham,1000)
CKi      call uham(v_sham,-1,-000,0,0,0)
CKi      call redfrr(v_sham,nint(dval(v_sham,1)))
      v_sham%size=uham_size()

Cgetarg...info...           NG:: alfsi , real(8) :: v_sham%alfsi 1
Cgetarg...info...           NG:: dabc , real(8) :: v_sham%dabc 3
Cgetarg...info...           NG:: elind , real(8) :: v_sham%elind 1
Cgetarg...info...           NG:: nmto , integer(8) :: v_sham%nmto 1
Cgetarg...info...           NG:: kmto , real(8) :: v_sham%kmto 6
Cgetarg       call pack5('ham alfsi dabc elind nmto kmto',v_sham,
Cgetarg      .  alfsi,dabc,elind,nmto,kmto)
       
       v_sham%alfsi=alfsi 
       i_copy_size=size(v_sham%dabc) 
       call dcopy(i_copy_size,dabc,1,v_sham%dabc,1) 
       v_sham%elind=elind 
       v_sham%nmto=nmto 
       i_copy_size=size(v_sham%kmto) 
       call dcopy(i_copy_size,kmto,1,v_sham%kmto,1) 

C     Replicate ctrl->lncol in ham->lncol, ditto for lham,lgen3
Cgetarg       call pack1('ham lncol',v_sham,igets('ctrl lncol',v_sctrl))
       
       v_sham%lncol=(int(v_sctrl%lncol)) 

Cgetarg       call pack1('ham lham',v_sham,igets('ctrl lham',v_sctrl))
       
c       v_sham%lham=(int(v_sctrl%lham)) 

Cgetarg       call pack1('ham lgen3',v_sham,igets('ctrl lgen3',v_sctrl))
c       v_sham%lgen3=(int(v_sctrl%lgen3)) 

C     Mop up integer parts of sigp
      sigp(1) = sigp_mode
      sigp(2) = sigp_nmin
      sigp(4) = sigp_nmax
Cgetarg...info...           NG:: nkaph , integer(8) :: v_sham%nkaph 1
Cgetarg...info...           NG:: pmax , real(8) :: v_sham%pmax 10
Cgetarg...info...           NG:: pmin , real(8) :: v_sham%pmin 10
Cgetarg...info...           NG:: lrsig , integer(8) :: v_sham%lsig 1
Cgetarg...info...           NG:: sigp , real(8) :: v_sham%sigp 10
Cgetarg       call pack5('ham nkaph pmax pmin lsig sigp',v_sham,
Cgetarg      .  nkaph,pmax,pmin,lrsig,sigp)
       
       v_sham%nkaph=nkaph 
       i_copy_size=size(v_sham%pmax) 
       call dcopy(i_copy_size,pmax,1,v_sham%pmax,1) 
       i_copy_size=size(v_sham%pmin) 
       call dcopy(i_copy_size,pmin,1,v_sham%pmin,1) 
       v_sham%lsig=lrsig 
       i_copy_size=size(v_sham%sigp) 
       call dcopy(i_copy_size,sigp,1,v_sham%sigp,1) 

Cgetarg...info...           NG:: ham_qss , real(8) :: v_sham%qss 4
Cgetarg...info...           NG:: rsrnge , real(8) :: v_sham%rsrnge 1
Cgetarg...info...           NG:: rsstol , real(8) :: v_sham%rsstol 1
Cgetarg...info...           NG:: ham_udiag , integer(8) :: v_sham%udiag 1
Cgetarg       call pack5('ham qss rsrnge rsstol udiag',v_sham,
Cgetarg      .  ham_qss,rsrnge,rsstol,ham_udiag,0)
       
       i_copy_size=size(v_sham%qss) 
       call dcopy(i_copy_size,ham_qss,1,v_sham%qss,1) 
       v_sham%rsrnge=rsrnge 
       v_sham%rsstol=rsstol 
       v_sham%udiag=ham_udiag 

C     Parameters for APW
Cgetarg...info...           NG:: pwmode , integer(8) :: v_sham%pwmode 1
Cgetarg...info...           NG:: npwpad , integer(8) :: v_sham%npwpad 1
Cgetarg...info...           NG:: pwemin , real(8) :: v_sham%pwemin 1
Cgetarg...info...           NG:: pwemax , real(8) :: v_sham%pwemax 1
Cgetarg       call pack5('ham pwmode npwpad pwemin pwemax',v_sham,pwmode,
Cgetarg      .  npwpad,pwemin,pwemax,0)
       
       v_sham%pwmode=pwmode 
       v_sham%npwpad=npwpad 
       v_sham%pwemin=pwemin 
       v_sham%pwemax=pwemax 

Cgetarg...info...           NG:: oveps , real(8) :: v_sham%oveps 1
Cgetarg       call pack1('ham oveps',v_sham,oveps)
       
       v_sham%oveps=oveps 


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

Cgetarg...info...           NG:: vmtz , real(8) :: v_spot%vmtz0 1
Cgetarg       call pack5('pot vmtz0',v_spot,vmtz,0,0,0,0)
       
       v_spot%vmtz0=vmtz 


C --- Allocate and copy input to smix ---
CKi      call defrr(v_smix,1000)
CKi      call umix(v_smix,-1,-000,0,0,0)
CKi      call redfrr(v_smix,nint(dval(v_smix,1)))
       v_smix%size=umix_size()

ckino's fix for SR11000 size() return's integer(8), right?
c      nnn = size(smix)
      nnn = umix_size()
      call dcopy(nnn,vmix,1,v_smix,1)

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
Cgetarg...info...           NG:: gd_nmodt , integer(8) :: v_smove%nmodt 1
Cgetarg...info...           NG:: gd_modt , integer(8) :: v_smove%modt 3
Cgetarg...info...           NG:: gd_ct , real(8) :: v_smove%ct 3
Cgetarg...info...           NG:: prmint , real(8) :: v_smove%prmint 20
Cgetarg       call pack5('move gyro nmodt modt ct prmint',v_smove,
Cgetarg      .  2d0,gd_nmodt,gd_modt,gd_ct,prmint)
       
       v_smove%gyro=2d0 
       v_smove%nmodt=gd_nmodt 
       i_copy_size=size(v_smove%modt) 
       call ii8copy(i_copy_size,gd_modt,1,v_smove%modt,1) 
       i_copy_size=size(v_smove%ct) 
       call dcopy(i_copy_size,gd_ct,1,v_smove%ct,1) 
       i_copy_size=size(v_smove%prmint) 
       call dcopy(i_copy_size,prmint,1,v_smove%prmint,1) 

Cgetarg...info...           NG:: move_kt , real(8) :: v_smove%kt 1
Cgetarg...info...           NG:: move_ts , real(8) :: v_smove%ts 1
Cgetarg...info...           NG:: move_tsequ , real(8) :: v_smove%tsequ 1
Cgetarg...info...           NG:: move_tstot , real(8) :: v_smove%tstot 1
Cgetarg       call pack5('move kt ts tsequ tstot',v_smove,
Cgetarg      .  move_kt,move_ts,move_tsequ,move_tstot,0)
       
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
Cgetarg...info...           NG:: tcf_adec , real(8) :: v_sstr%adec 1
Cgetarg...info...           NG:: str_mode , integer(8) :: v_sstr%amode 1
Cgetarg...info...           NG:: str_drwats , real(8) :: v_sstr%drwats 1
Cgetarg...info...           real(8) :: xx(n0*2) , real(8) :: v_sstr%iinv 5
Cgetarg       call pack5('str adec amode drwats iinv',v_sstr,
Cgetarg      .  tcf_adec,str_mode,str_drwats,xx,0)
       
c       v_sstr%adec=tcf_adec 
c       v_sstr%amode=str_mode ! amode is dummy now
c       v_sstr%drwats=str_drwats 
c       i_copy_size=size(v_sstr%iinv) 
c       call dcopy(i_copy_size,xx,1,v_sstr%iinv,1) 

Cgetarg...info...           integer :: isw , integer(8) :: v_sstr%lequiv 1
Cgetarg...info...           integer :: i , integer(8) :: v_sstr%loka 1
Cgetarg...info...           NG:: str_lmaxw , integer(8) :: v_sstr%lmaxw 1
Cgetarg...info...           integer :: isw , integer(8) :: v_sstr%lshow 1
Cgetarg...info...           NG:: str_mxnbr , integer(8) :: v_sstr%mxnbr 1
Cgetarg       call pack5('str lequiv loka lmaxw lshow mxnbr',v_sstr,
Cgetarg      .  isw(str_lequiv1),i,str_lmaxw,isw(str_lshow1),str_mxnbr)
       
c       v_sstr%lequiv=isw(str_lequiv1) 
c       v_sstr%loka=i 
c       v_sstr%lmaxw=str_lmaxw 
c       v_sstr%lshow=isw(str_lshow1) 
       v_sstr%mxnbr=str_mxnbr 

Cgetarg...info...           NG:: tcf_nalf , integer(8) :: v_sstr%nalf 1
Cgetarg...info...           NG:: tcf_nbisi , integer(8) :: v_sstr%nbisi 3
Cgetarg...info...           NG:: tcf_ncupl , integer(8) :: v_sstr%ncupl 1
Cgetarg...info...           NG:: tcf_ndust , integer(8) :: v_sstr%ndust 1
Cgetarg...info...           NG:: tcf_wztcf , real(8) :: v_sstr%wztcf 1
Cgetarg       call pack5('str nalf nbisi ncupl ndust wztcf',v_sstr,
Cgetarg      .  tcf_nalf,tcf_nbisi,tcf_ncupl,tcf_ndust,tcf_wztcf)
       
c       v_sstr%nalf=tcf_nalf 
c       i_copy_size=size(v_sstr%nbisi) 
c       call ii8copy(i_copy_size,tcf_nbisi,1,v_sstr%nbisi,1) 
c       v_sstr%ncupl=tcf_ncupl 
c       v_sstr%ndust=tcf_ndust 
c       v_sstr%wztcf=tcf_wztcf 

Cgetarg...info...           NG:: str_nkaps , integer(8) :: v_sstr%nkaps 1
Cgetarg...info...           NG:: str_rmax , real(8) :: v_sstr%rmax 1
Cgetarg...info...           NG:: str_kaps , real(8) :: v_sstr%kaps 6
Cgetarg       call pack5('str nkaps rmax rfit kaps',v_sstr,
Cgetarg      .  str_nkaps,str_rmax,0.8d0,str_kaps,0)
       
c       v_sstr%nkaps=str_nkaps 
       v_sstr%rmax=str_rmax 
c       v_sstr%rfit=0.8d0 
c       i_copy_size=size(v_sstr%kaps) 
c       call dcopy(i_copy_size,str_kaps,1,v_sstr%kaps,1) 

Cgetarg...info...           NG:: str_rmaxg , real(8) :: v_sstr%rmaxg 1
Cgetarg...info...           NG:: str_ivl , integer(8) :: v_sstr%ivl 1
Cgetarg       call pack5('str rmaxg ivl',v_sstr,str_rmaxg,str_ivl,0,0,0)
       
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
Cgetarg       call spackv(1,'spec eh3',v_sspec,1,nspec,-0.5d0)
       i_copy_size=1; 
       do i_spackv=1,nspec 
       call spackv_array_copy_r8_r8('p',v_sspec(i_spackv)%eh3,i_copy_size,1,-0.5d0)
       enddo

Cgetarg       call spackv(1,'spec etf',v_sspec,1,nspec,-1d0)
       i_copy_size=1; 
       do i_spackv=1,nspec 
       call spackv_array_copy_r8_r8('p',v_sspec(i_spackv)%etf,i_copy_size,1,-1d0)
       enddo

Cgetarg       call spackv(1,'spec lmxf',v_sspec,1,nspec,2*nl-2)
       i_copy_size=1; 
       do i_spackv=1,nspec 
       call spackv_array_copy_i8_i('p',v_sspec(i_spackv)%lmxf,i_copy_size,1,2*nl-2)
       enddo

Cgetarg       call spackv(1,'spec norp',v_sspec,1,nspec,2)
       i_copy_size=1; 
       do i_spackv=1,nspec 
       call spackv_array_copy_i8_i('p',v_sspec(i_spackv)%norp,i_copy_size,1,2)
       enddo

Cgetarg       call spackv(1,'spec vmtz',v_sspec,1,nspec,-0.5d0)
       i_copy_size=1; 
       do i_spackv=1,nspec 
       call spackv_array_copy_r8_r8('p',v_sspec(i_spackv)%vmtz,i_copy_size,1,-0.5d0)
       enddo


      do  j = 1, nspec
        slabl_(j) = slabl(j)
Cgetarg...info...          check_variable pb1(j) not found
Cgetarg         call spacks(1,'spec pb1',v_sspec,pb1(j),j,j)
         do i_spacks=j,j 
         call spacks_copy('p',v_sspec(i_spacks)%pb1,j,j,pb1(j),i_spacks)
         enddo

Cgetarg...info...          check_variable pb2(j) not found
Cgetarg         call spacks(1,'spec pb2',v_sspec,pb2(j),j,j)
         do i_spacks=j,j 
         call spacks_copy('p',v_sspec(i_spacks)%pb2,j,j,pb2(j),i_spacks)
         enddo

Cgetarg...info...          check_variable coreh(j) not found
Cgetarg         call spacks(1,'spec coreh',v_sspec,coreh(j),j,j)
         do i_spacks=j,j 
         call spacks_copy('p',v_sspec(i_spacks)%coreh,j,j,coreh(j),i_spacks)
         enddo

Cgetarg...info...          check_variable slabl(j) not found
Cgetarg         call spacks(1,'spec name',v_sspec,slabl(j),j,j)
         do i_spacks=j,j 
         call spacks_copy('p',v_sspec(i_spacks)%name,j,j,slabl(j),i_spacks)
         enddo

Cgetarg...info...           NG:: spec_a(j) , real(8) :: v_sspec%a 1
Cgetarg...info...           NG:: nr(j) , integer(8) :: v_sspec%nr 1
Cgetarg...info...           NG:: alpha(1,j) , real(8) :: v_sspec%alpha 10
Cgetarg...info...           NG:: coreq(1,j) , real(8) :: v_sspec%coreq 2
Cgetarg         call pack5('spec a nr alpha coreq ',v_sspec,j,
Cgetarg      .    spec_a(j),nr(j),alpha(1,j),coreq(1,j))
         
         v_sspec(j)%a=spec_a(j) 
         v_sspec(j)%nr=nr(j) 
         i_copy_size=size(v_sspec(j)%alpha) 
         call dcopy(i_copy_size,alpha(1,j),1,v_sspec(j)%alpha,1) 
         i_copy_size=size(v_sspec(j)%coreq) 
         call dcopy(i_copy_size,coreq(1,j),1,v_sspec(j)%coreq,1) 

Cgetarg...info...           NG:: lxi(j) , integer(8) :: v_sspec%lxi 1
Cgetarg...info...           NG:: nxi(j) , integer(8) :: v_sspec%nxi 1
Cgetarg...info...           NG:: exi(1,j) , real(8) :: v_sspec%exi 10
Cgetarg         call pack5('spec lxi nxi exi',v_sspec,j,
Cgetarg      .    lxi(j),nxi(j),exi(1,j),0)
         
         v_sspec(j)%lxi=lxi(j) 
         v_sspec(j)%nxi=nxi(j) 
         i_copy_size=size(v_sspec(j)%exi) 
         call dcopy(i_copy_size,exi(1,j),1,v_sspec(j)%exi,1) 

Cgetarg...info...           NG:: grp(j) , integer(8) :: v_sspec%group 1
Cgetarg...info...           NG:: grp2(j) , integer(8) :: v_sspec%grp2 1
Cgetarg...info...           NG:: hcr(1,j) , real(8) :: v_sspec%hcr 10
Cgetarg...info...           NG:: idmod(1,j) , integer(8) :: v_sspec%idmod 10
Cgetarg         call pack5('spec group grp2 hcr idmod',v_sspec,j,
Cgetarg      .    grp(j),grp2(j),hcr(1,j),idmod(1,j))
         
         v_sspec(j)%group=grp(j) 
         v_sspec(j)%grp2=grp2(j) 
         i_copy_size=size(v_sspec(j)%hcr) 
         call dcopy(i_copy_size,hcr(1,j),1,v_sspec(j)%hcr,1) 
         i_copy_size=size(v_sspec(j)%idmod) 
         call ii8copy(i_copy_size,idmod(1,j),1,v_sspec(j)%idmod,1) 

Cgetarg...info...           NG:: ehvl(1,j) , real(8) :: v_sspec%ehvl 10
Cgetarg         call pack5('spec ehvl',v_sspec,j,ehvl(1,j),0,0,0)
         
         i_copy_size=size(v_sspec(j)%ehvl) 
         call dcopy(i_copy_size,ehvl(1,j),1,v_sspec(j)%ehvl,1) 

C       Set idxdn
        call ivset(ix,1,n0*nkap0,1)
        call icopy(1+lmxb(j),idxdn(1,j),1,ix,1)
Cgetarg...info...           integer :: ix(n0*nkap0) , integer(8) :: v_sspec%idxdn 30
Cgetarg...info...           NG:: idu(1,j) , integer(8) :: v_sspec%idu 4
Cgetarg...info...           NG:: jh(1,j) , real(8) :: v_sspec%jh 4
Cgetarg...info...           NG:: uh(1,j) , real(8) :: v_sspec%uh 4
Cgetarg         call pack5('spec idxdn idu jh uh',v_sspec,j,
Cgetarg      .    ix,idu(1,j),jh(1,j),uh(1,j))
         
         i_copy_size=size(v_sspec(j)%idxdn) 
         call ii8copy(i_copy_size,ix,1,v_sspec(j)%idxdn,1) 
         i_copy_size=size(v_sspec(j)%idu) 
         call ii8copy(i_copy_size,idu(1,j),1,v_sspec(j)%idu,1) 
         i_copy_size=size(v_sspec(j)%jh) 
         call dcopy(i_copy_size,jh(1,j),1,v_sspec(j)%jh,1) 
         i_copy_size=size(v_sspec(j)%uh) 
         call dcopy(i_copy_size,uh(1,j),1,v_sspec(j)%uh,1) 

Cgetarg...info...           NG:: kmxt(j) , integer(8) :: v_sspec%kmxt 1
Cgetarg...info...           NG:: kmxv(j) , integer(8) :: v_sspec%kmxv 1
Cgetarg...info...           NG:: lfoca(j) , integer(8) :: v_sspec%lfoca 1
Cgetarg...info...           NG:: rsmv(j) , real(8) :: v_sspec%rsmv 1
Cgetarg         call pack5('spec kmxt kmxv lfoca rsmv',v_sspec,j,
Cgetarg      .    kmxt(j),kmxv(j),lfoca(j),rsmv(j))
         
         v_sspec(j)%kmxt=kmxt(j) 
         v_sspec(j)%kmxv=kmxv(j) 
         v_sspec(j)%lfoca=lfoca(j) 
         v_sspec(j)%rsmv=rsmv(j) 

Cgetarg...info...           NG:: lmxa(j) , integer(8) :: v_sspec%lmxa 1
Cgetarg...info...           NG:: lmxb(j) , integer(8) :: v_sspec%lmxb 1
Cgetarg...info...           NG:: lmxl(j) , integer(8) :: v_sspec%lmxl 1
Cgetarg...info...           NG:: lmxpb(j) , integer(8) :: v_sspec%lmxpb 1
Cgetarg         call pack5('spec lmxa lmxb lmxl lmxpb',v_sspec,j,
Cgetarg      .    lmxa(j),lmxb(j),lmxl(j),lmxpb(j))
         
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
Cgetarg...info...           NG:: mass(j) , real(8) :: v_sspec%mass 1
Cgetarg...info...           integer :: i , integer(8) :: v_sspec%mxcst 1
Cgetarg...info...           real(8) :: orbp(n0,2,nkap0) , real(8) :: v_sspec%orbp 60
Cgetarg         call pack5('spec mass mxcst orbp',v_sspec,j,mass(j),i,orbp,0)
         
         v_sspec(j)%mass=mass(j) 
         v_sspec(j)%mxcst=i 
         i_copy_size=size(v_sspec(j)%orbp) 
         call dcopy(i_copy_size,orbp,1,v_sspec(j)%orbp,1) 

C       Pack P,Q,PZ for both spins        
        call dpzero(orbp,n0*6)
        call dcopy(n0*nsp,pnu(1,1,j),1,orbp(1,1,1),1)
        call dcopy(n0*nsp,pz(1,1,j),1,orbp(1,1,2),1)
        call dcopy(n0*nsp,qnu(1,1,j),1,orbp(1,1,3),1)
Cgetarg...info...           real(8) :: orbp(n0,2,nkap0) , real(8) :: v_sspec%p 20
Cgetarg...info...           real(8) :: orbp(n0,2,nkap0) , real(8) :: v_sspec%pz 20
Cgetarg...info...           real(8) :: orbp(n0,2,nkap0) , real(8) :: v_sspec%q 20
Cgetarg...info...           NG:: z(j) , real(8) :: v_sspec%z 1
Cgetarg         call pack5('spec p pz q z',v_sspec,j,orbp,
Cgetarg      .    orbp(1,1,2),orbp(1,1,3),z(j))
         
         i_copy_size=size(v_sspec(j)%p) 
         call dcopy(i_copy_size,orbp,1,v_sspec(j)%p,1) 
         i_copy_size=size(v_sspec(j)%pz) 
         call dcopy(i_copy_size,orbp(1,1,2),1,v_sspec(j)%pz,1) 
         i_copy_size=size(v_sspec(j)%q) 
         call dcopy(i_copy_size,orbp(1,1,3),1,v_sspec(j)%q,1) 
         v_sspec(j)%z=z(j) 

Cgetarg...info...           NG:: colxbs(1,j) , real(8) :: v_sspec%colxbs 3
Cgetarg...info...           NG:: radxbs(j) , real(8) :: v_sspec%radxbs 1
Cgetarg         call pack5('spec colxbs radxbs',v_sspec,j,
Cgetarg      .    colxbs(1,j),radxbs(j),0,0)
         
         i_copy_size=size(v_sspec(j)%colxbs) 
         call dcopy(i_copy_size,colxbs(1,j),1,v_sspec(j)%colxbs,1) 
         v_sspec(j)%radxbs=radxbs(j) 

c        call pack5('spec rcfa rcut rfoca rg rham',v_sspec,j,
Cgetarg...info...           NG:: rcfa(1,j) , real(8) :: v_sspec%rcfa 2
Cgetarg...info...           NG:: rcut(j) , real(8) :: v_sspec%rcut 1
Cgetarg...info...           NG:: rfoca(j) , real(8) :: v_sspec%rfoca 1
Cgetarg...info...           NG:: rg(j) , real(8) :: v_sspec%rg 1
Cgetarg         call pack5('spec rcfa rcut rfoca rg',v_sspec,j,
Cgetarg      .    rcfa(1,j),rcut(j),rfoca(j),rg(j))
         
         i_copy_size=size(v_sspec(j)%rcfa) 
         call dcopy(i_copy_size,rcfa(1,j),1,v_sspec(j)%rcfa,1) 
         v_sspec(j)%rcut=rcut(j) 
         v_sspec(j)%rfoca=rfoca(j) 
         v_sspec(j)%rg=rg(j) 

Cgetarg...info...           NG:: rmt(j) , real(8) :: v_sspec%rmt 1
Cgetarg...info...           NG:: rs3(j) , real(8) :: v_sspec%rs3 1
Cgetarg...info...           NG:: rsma(j) , real(8) :: v_sspec%rsma 1
Cgetarg...info...           NG:: rsmfa(j) , real(8) :: v_sspec%rsmfa 1
Cgetarg         call pack5('spec rmt rs3 rsma rsmfa',v_sspec,j,
Cgetarg      .    rmt(j),rs3(j),rsma(j),rsmfa(j))
         
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
Cgetarg...info...           NG:: dv(j) , real(8) :: v_sspec%dv 1
Cgetarg...info...           NG:: eref(j) , real(8) :: v_sspec%eref 1
Cgetarg...info...           NG:: rham(j) , real(8) :: v_sspec%rham 1
Cgetarg...info...           NG:: rint(j) , real(8) :: v_sspec%rint 1
Cgetarg         call pack5('spec dv eref rham rint',v_sspec,j,
Cgetarg      .    dv(j),eref(j),rham(j),rint(j))
         
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
Cgetarg...info...           NG:: ips(j) , integer(8) :: v_ssite%spec 1
Cgetarg...info...           NG:: ips(j) , integer(8) :: v_ssite%class 1
Cgetarg...info...           NG:: dpole(1,j) , real(8) :: v_ssite%dpole 3
Cgetarg...info...           NG:: mpole(j) , real(8) :: v_ssite%mpole 1
Cgetarg         call pack5('site spec class dpole mpole',v_ssite,j,
Cgetarg      .    ips(j),ips(j),dpole(1,j),mpole(j))
         
         v_ssite(j)%spec=ips(j) 
         v_ssite(j)%class=ips(j) 
         i_copy_size=size(v_ssite(j)%dpole) 
         call dcopy(i_copy_size,dpole(1,j),1,v_ssite(j)%dpole,1) 
         v_ssite(j)%mpole=mpole(j) 

Cgetarg...info...          check_variable slabl(ips(j)) not found
Cgetarg         call spacks(1,'site clabel',v_ssite,slabl(ips(j)),j,j)
         do i_spacks=j,j 
         call spacks_copy('p',v_ssite(i_spacks)%clabel,j,j,slabl(ips(j)),i_spacks)
         enddo

Cgetarg...info...           NG:: ipl(j) , integer(8) :: v_ssite%pl 1
Cgetarg...info...           NG:: plv(j) , integer(8) :: v_ssite%plv 1
Cgetarg...info...           NG:: pos(1,j) , real(8) :: v_ssite%pos 3
Cgetarg...info...           NG:: vel(1,j) , real(8) :: v_ssite%vel 3
Cgetarg         call pack5('site pl plv pos vel',v_ssite,j,
Cgetarg      .    ipl(j),plv(j),pos(1,j),vel(1,j))
         
         v_ssite(j)%pl=ipl(j) 
         v_ssite(j)%plv=plv(j) 
         i_copy_size=size(v_ssite(j)%pos) 
         call dcopy(i_copy_size,pos(1,j),1,v_ssite(j)%pos,1) 
         i_copy_size=size(v_ssite(j)%vel) 
         call dcopy(i_copy_size,vel(1,j),1,v_ssite(j)%vel,1) 

Cgetarg...info...           NG:: vshft(j) , real(8) :: v_ssite%vshft 1
Cgetarg...info...           NG:: irlx(1,j) , integer(8) :: v_ssite%relax 3
Cgetarg...info...           NG:: eula(1,j) , real(8) :: v_ssite%eula 3
Cgetarg         call pack5('site vshft relax eula',v_ssite,j,
Cgetarg      .    vshft(j),irlx(1,j),eula(1,j),0)
         
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
        fileid = 'basp'
        if (strn(9:12) .eq. ':fn=') then
          fileid = strn(13:13+63) !tk to avoid complaint; replace fileid = strn(13:)
        else
        endif

C       Number of envelope functions before modifications
        call uspecb(0,1,v_sspec,1,nspec,xx,xx,xx,k1)

C       Read new parms; return in k2 max no. env. fns/site read
        call strip(fileid,i,j)
        ifi = fopna(fileid(1:j),-1,0)
        rewind ifi
        if (.not. ioorbp(111,2,1,nspec,v_sspec,k2,ifi))
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
          call uspecb(0,-1,v_sspec,1,nspec,xx,xx,xx,xx)
Cgetarg           nlmax = igets('ctrl nl',v_sctrl)**2
           nlmax = int(v_sctrl%nl) * * 2 

          xx(1) = dglob('mxorb',dble(k)*nlmax,1)
        endif
      endif

C     Add dalat to alat
Cgetarg       call pack1('lat alat',v_slat,dgets('lat alat',v_slat)+dalat)
       
       v_slat%alat=(v_slat%alat)+dalat 

C     Replicate ctrl->lncol in ham->lncol
Cgetarg...info...           integer :: lncol , integer(8) :: v_sham%lncol 1
Cgetarg       call pack1('ham lncol',v_sham,lncol)
       
       v_sham%lncol=lncol 

c      lsx1 = mod(lsx,2)
Cgetarg       if (lasa32) call pack1('mix nsave',v_smix,3)
C$$$       if ( lasa32 ) then 
       
C$$$       v_smix%nsave=3 
C$$$       endif 

C ... Use true spherical harmonics
Cgetarg       if (bittst(lncol,4) .or. ctrl_loptc .ne. 0 .or. lrel .eq. 2)
Cgetarg      .  call lsets('ctrl lham',v_sctrl,.true.,256)

c       if ( bittst ( lncol , 4 ) .or. ctrl_loptc .ne. 0 .or. lrel .eq. 
c     . 2 ) then 
c       if ( bittst ( lncol , 4 ) .or. lrel .eq. 2 ) then 
c       call lsets_bitop_i8(v_sctrl%lham,1,.true.,256)
c       endif 


C ... Dirac equation requires spin polarization
Cgetarg       if (nsp.eq.1 .and. igets('ctrl lrel',v_sctrl).eq.2) then
       if ( nsp.eq.1 .and. int(v_sctrl%lrel) .eq.2 ) then 

        call rx('rdccat: Dirac equation requires NSPIN=2')
      endif

C ... Suppress inversion when noncollinear magnetism, SX, NLO
c      if (lncol+lsx1 .ne. 0 .or. ctrl_loptc .ge. 10 .or.
c     .   (mod(lscr,10) .eq. 1 .and. prgnam .eq. 'LM'))
      if (lncol .ne. 0)
     .  lqp = lqp-bitand(lqp,1)+1

Cgetarg...info...           integer :: lqp , integer(8) :: v_sctrl%lqp 1
Cgetarg       call pack5('ctrl lqp',v_sctrl,lqp,0,0,0,0) !repack lqp
       v_sctrl%lqp=lqp 


C ... Special pgf initialization
C$$$      if (ctrl_lpgf(1) .ne. 0) then
C$$$Cgetarg         call lsets('ctrl lmet',v_sctrl,.false.,2)
         
C$$$         call lsets_bitop_i8(v_sctrl%lmet,1,.false.,2)

C$$$        call redfrr(v_ssite,3*nsite*nint(dval(v_ssite,1)))
C$$$      endif

C ... Setup for idxdn ... ctrl->lham,4 is automatic downfolding switch
Cgetarg       j = igets('ctrl lham,4',v_sctrl)/4
c     j = iand(4,int(v_sctrl%lham)) /4 
c     j = 2*(1-j)
      j= 2

C     No screening => no downfolding; also lmxb<l<=lmxa => 'high'
C     Probably ought to have lmxb<l<=lmxa => 'high' always
Cgetarg       if (lgors('ctrl lbas,1',v_sctrl)) j = 3
       if ( iand(1,int(v_sctrl%lbas)) .ne.0 ) j = 3 

C     nfp-style basis:
Cgetarg       if (lgors('ctrl lbas,2',v_sctrl)) j = j+10
       if ( iand(2,int(v_sctrl%lbas)) .ne.0 ) j = j+10 

      call suidx(nglob('nkaph'),j,nspec,v_sspec)

C ... Set some global variables
      xx(1) = dglob('nspec',dble(nspec),1)
      xx(1) = dglob('nbas',dble(nbas),1)
      xx(1) = dglob('nbasp',dble(nbasp),1)
      xx(1) = dglob('nsp',dble(nsp),1)
      xx(1) = dglob('nl',dble(nl),1)
Cgetarg       xx(1) = dglob('avw',dgets('lat avw',v_slat),1)
       xx ( 1 ) = dglob ( 'avw' , v_slat%avw , 1 ) 

      xx(1) = dglob('lrel',dble(lrel),1)
Cgetarg       lxcf  = igets('ctrl lxcf',v_sctrl)
       lxcf = int(v_sctrl%lxcf) 

      xx(1) = dglob('lxcf',dble(lxcf),1)
      xx(1) = dglob('stdo',dble(stdo),1)
      xx(1) = dglob('stdl',dble(stdl),1)
      xx(1) = dglob('stde',dble(stde),1)
      nspc = 1
c      if (bitand(lncol,1+2+4) .ne. 0) nspc = 2
      if (iand(lncol,4) .ne. 0) nspc = 2
      xx(1) = dglob('nspc',dble(nspc),1)
C     Make nat = number of real atoms as nbas - # sites w/ floating orbitals
      if (procid .eq. master) then
      nat = nbas
      do  i = 1, nbas
Cgetarg...info...           integer :: j , integer(8) :: v_ssite%spec 1
Cgetarg         call upack('site spec',v_ssite,i,j,0,0,0)
         
         j=v_ssite(i)%spec

Cgetarg...info...           integer :: l , integer(8) :: v_sspec%lmxa 1
Cgetarg         call upack('spec lmxa',v_sspec,j,l,0,0,0)
         
         l=v_sspec(j)%lmxa

        if (l .eq. -1) nat = nat-1
      enddo
      endif
      call mpibc1(nat,1,2,0,'rdctrl','nat')
      xx(1) = dglob('nat',dble(nat),1)

C ... Set modep
      ix(1) = 2
      ix(2) = 2
      ix(3) = 2
Cgetarg       if (lgors('ctrl lpgf,-1',v_sctrl)) ix(3) = 0
C$$$       if ( iand(-1,int(v_sctrl%lpgf(1))) .ne.0 ) ix ( 3 ) = 0 

Cgetarg...info...           integer :: ix(n0*nkap0) , integer(8) :: v_sctrl%modep 3
Cgetarg       call pack1('ctrl modep',v_sctrl,ix)
       
       i_copy_size=size(v_sctrl%modep) 
       call ii8copy(i_copy_size,ix,1,v_sctrl%modep,1) 

Cgetarg...info...           integer :: lxcf , integer(8) :: v_sham%lxcf 1
Cgetarg       call pack1('ham lxcf',v_sham,lxcf)
       
       v_sham%lxcf=lxcf 


C ... Count LDA+U blocks (printout only)
      call defi(owk,-nbas)
      call pshpr(0)
      call suldau(nbas,v_sspec,v_ssite,nlibu,k,w(owk))
Cgetarg...info...           integer :: nlibu , integer(8) :: v_sham%nlibu 1
Cgetarg...info...           integer :: k , integer(8) :: v_sham%lmaxu 1
Cgetarg       call pack5('ham nlibu lmaxu',v_sham,nlibu,k,0,0,0)
       
       v_sham%nlibu=nlibu 
       v_sham%lmaxu=k 

      call poppr
      call rlse(owk)

C     Free arrays used to read input
      deallocate(pnu,qnu,pz,amom,idmod,rsmh,eh,rsmh2,eh2,pb1,pb2,
     .  lmxpb,qpol,stni,tbvso,iq1,ivso,rg,rsma,rfoca,rsmfa,rcfa,nxi,
     .  exi,rint,rcut,coreq,mass,colxbs,radxbs,rs3,rham,idxdn,hcr,
     .  rmt,alpha,idu,uh,jh,dv,grp,grp2,mxcst1,mxcst2,mxcst4,kmxt,kmxv,
     .  lfoca,eref,lmxl,lxi,coreh,lmxa,lmxb,spec_a,z,nr,rsmv)
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
     .    strn,scrwid,-lgunit(k),vn(1)-vn2(1)/1d5,
     .    iprint(),iprint()-iprt(1),iprt(1),iprt(2)-iprt(1),iprt(2))
        else
          call awrit7('%a  vn %,2d('//trim(prgnam)//' %,1d)  '//
     .      'verb %i%?#n#,%i#%j#%?#n#,%i#%j#',
     .    strn,scrwid,-lgunit(k),vn(1)-vn2(1)/1d5,vn(2),
     .    iprint(),iprint()-iprt(1),iprt(1),iprt(2)-iprt(1),iprt(2))
        endif

Cgetarg...info...           integer :: lxcf , integer(8) :: v_sctrl%lxcf 1
Cgetarg         call upack('ctrl lxcf',v_sctrl,lxcf,0,0,0,0)
         
         lxcf=v_sctrl%lxcf

Cgetarg         i = igets('ctrl lham,256',v_sctrl)
c         i = iand(256,int(v_sctrl%lham)) 

Cgetarg         lfrzw = isw(lgors('ctrl lbas,16',v_sctrl))
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
c     .      isw(lasa32),lrsig,lsx)
     .      0,lrsig,0)
          call awrit0(strn,' ',-80,lgunit(k))
        endif
        call isanrg(mod(lxcf,10),0,4,prgnam,'XC functional',.true.)
        call isanrg(lxcf/10,0,4,prgnam,'GGA functional',.true.)
        if (mod(lxcf,10) .eq. 2 .and. lxcf/10 .gt. 1 .or.
     .      mod(lxcf,10) .ne. 2 .and. lxcf/10 .eq. 1) call info0
     .    (10,0,0,'%10f(warning) mixing incompatible functionals')

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
Cgetarg           call awrit1(' special:%10p'//
Cgetarg      .      '%?;n==1; core-level-optics,;;'//
Cgetarg      .      '%b %b',strn,scrwid,0,igets('optic cls',w(osoptc)))
           call awrit1 ( ' special:%10p' // '%?;n==1; core-level-optics,;;' 
     .     // '%b %b' , strn , scrwid , 0 , 0 )

          if (strn .ne. ' special:') call awrit0(strn,' ',-80,lgunit(k))
        endif

c       ltmp = prgnam.eq.'LM'   .or. prgnam.eq.'LMMC'  .or.
        ltmp = prgnam.eq.'LM'   .or. 
     .         prgnam.eq.'LMF'  .or. prgnam.eq.'LMFGWD' .or.
     .         prgnam.eq.'LMGF' .or. prgnam.eq.'LMPG' .or. 
     .         prgnam.eq.'LMDOS'.or. prgnam.eq.'TBE'

Cgetarg         if (ltmp) call awrit7(' bz:%10p'//
Cgetarg      .    '%?;n; metal,; nonmetal,;'//
Cgetarg      .    '%-1j%?;n>1;%b(%-1j%i),;;'//
Cgetarg      .    '%?;n; tetra,;;'//
Cgetarg      .    '%?;n; get-qp,;;'//
Cgetarg      .    '%?;n; invit,;;'//
Cgetarg      .    '%?;n; dens-mat,;;'//
Cgetarg      .    '%?;(n>0); %-1jmull=%i,;;'//
Cgetarg      .    '%?;n; fixed-spin-mom,;;%b ',
Cgetarg      .    strn,scrwid,lgunit(k),
Cgetarg      .    igets('bz lmet',v_sbz),
Cgetarg      .    bitand(lmet,2),igets('bz lio,1',v_sbz),
Cgetarg      .    igets('ctrl lqp,2',v_sctrl),
Cgetarg      .    igets('bz lio,8',v_sbz),
Cgetarg      .    igets('bz lmull',v_sbz),
Cgetarg      .    isw(dgets('bz fsmom',v_sbz).ne.0))
         if ( ltmp ) call awrit7 ( ' bz:%10p' // '%?;n; metal,; nonmetal,;' 
     .   // '%-1j%?;n>1;%b(%-1j%i),;;' // '%?;n; tetra,;;' // '%?;n; get-qp,;;' 
     .   // '%?;n; invit,;;' // '%?;n; dens-mat,;;' // '%?;(n>0); %-1jmull=%i,;;' 
     .   // '%?;n; fixed-spin-mom,;;%b ' , strn , scrwid , lgunit ( k 
     .   ) , int(v_sbz%lmet) , bitand ( lmet , 2 ) , iand(1,int(v_sbz%lio)) 
c    .   , iand(2,int(v_sctrl%lqp)) , iand(8,int(v_sbz%lio)) , int(v_sbz%lmull) 
     .   , iand(2,int(v_sctrl%lqp)) , .false.                , int(v_sbz%lmull) 
     .   , isw ( v_sbz%fsmom .ne.0 ) ) 


   90   continue
      endif

C --- Sanity checks and other initialization ---
Cgetarg       if (lgors('ctrl ldos,8',v_sctrl) .and.
Cgetarg      .    lgors('ctrl ldos,4+2',v_sctrl))
Cgetarg      .  call rx('inconsistent DOS options')
       if ( iand(8,int(v_sctrl%ldos)) .ne.0 .and. iand(4+2,int(v_sctrl%ldos)) .ne.0 
     . ) call rx ( 'inconsistent DOS options' ) 

C$$$Cgetarg       call rxx(lgors('ctrl lstonr,-1',v_sctrl) .and. nsp .eq. 2,
C$$$Cgetarg      .  'Stoner model not compatible with nsp=2')
C$$$       call rxx ( iand(-1,int(v_sctrl%lstonr(1))) .ne.0 .and. nsp .eq. 
C$$$     . 2 , 'Stoner model not compatible with nsp=2' ) 

C$$$Cgetarg       call rxx(lgors('ctrl lstonr,-1',v_sctrl) .and.
C$$$Cgetarg      . lgors('ctrl lham,4',v_sctrl),'Stoner not compatible with ADNF')
C$$$       call rxx ( iand(-1,int(v_sctrl%lstonr(1))) .ne.0 .and. iand(4,int(v_sctrl%lham)) .ne.0 
C$$$     . , 'Stoner not compatible with ADNF' ) 

C --- takao's simplified susite only for lmf, lmfgw.
c      print *,'www000 spec=',v_ssite(1)%spec
      call susite(v_sctrl,v_sham,v_spot,v_sarry,v_slat,
     .  v_sspec,v_ssite)


Cgetarg       call pack1('ham lham',v_sham,igets('ctrl lham',v_sctrl))
c       v_sham%lham=(int(v_sctrl%lham)) 


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
      call mpibc1(v_sspec,nspec*nint(dval(v_sspec,1)),4,0,
     .  'rdctrl','sspec')

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
