      module m_rdctrl2_func
      contains
      subroutine rdctrl2(recrd,recln,nrecs,
     . prgnam, !vrsion,vn,vn2,
     . pass2,slabl_,v_sbz,v_sctrl,v_sham,v_spot,v_slat,v_smix,v_sspec,v_ssite,
     . v_sstr,v_sarry,sstrnmix,sstrnsymg)  !,v_stb,v_smove
      use m_gtv,only: gtv_setst,gtv_setrcd

      use m_struc_func
      use m_globalvariables

      use m_rdctrl,only:
     &  readctrl,readctrlpq,
     &  NULLR,n0,nkap0,s_bz,s_ctrl,s_ham,s_pot,s_lat,s_mix,s_str,s_array,s_spec,s_site, !s_tb,,s_move
!!
     &  io_show,io_help,nvario,header,symg,
     &  frzwf,ham_ewald,ctrl_lfrce,ham_lxcf,gga,ftmesh,lrsig, !nmto,
     &  nsp,lrel,lso,ham_udiag, lfp, lat_gmax,tolft,elind,dqval,kmto,rsrnge,vmtz,scaledsigma,
     &  alfsi,dabc,rsstol,pmin,pmax,
     &  sigp,sigp_emin,sigp_emax,sigp_a,sigp_b,sigp_efit,
     &    sigp_mode,sigp_nmin,sigp_nmax,
!! OPTIONAL
     &    lves,                        
c lcd2,lcd4,
c     &    lncol1,lncol2,lncol4,lncol8,lncol16,lncol32,lncol64,
     &    lscr,smalit,nesabc, 
     &    lstonr,quit,nl,lpfloat, rmines,rmaxes,ham_qss,
!! STRUC
     &  lat_slat,dlat,alat,plat,lat_gam,dalat, 
     &     vol,avw,lat_dist,
     &     nbas,nbasp,nsite,nspec,nclass,lat_ldist,
!! SPEC
     &  omax1,omax2,wsrmax,sclwsr,   
     &    mxspec,slabl,lmxbx,lmxax,nkaph,
     &    mxcst1,mxcst2,mxcst4,idxdn, grp,grp2,idu,lmxb,lmxa,idmod,iq1,ivso,
     .    kmxt,kmxv,lfoca,lmxl,lxi,nxi,nr,lmxpb,
     .    rsmh,rsmh2,eh,eh2,
     .    hcr,rs3,rham,alpha,ehvl,
     .    dv, uh,jh,
     .    qpol,stni,tbvso,
     .    pnu,qnu,
     .    coreq,mass,colxbs,radxbs,
     .    rg,rsma,rfoca,
     .    rsmfa,rcfa,
     .    exi,rint,rcut,rmt,pz,
     .    amom,spec_a,z,eref,rsmv,pb1,pb2, coreh,nmcore,
!! ... SITE
     .  alabl, pos,vel,eula,vshft, ips,ipl,plv,irlx,ndelta,delta,mpole,dpole,   iantiferro,
!! ... Iterations
     . ctrl_tol,vmix,iter_maxit,iter_mix, !iter_amix,
!! ... BZ
     . bz_nabc, bz_lshft, bz_lmet,bz_n,bz_nevmx,bz_lmull,ctrl_ldos,bz_ndos,bz_fsmommethod,
     . bz_w,bz_ef,bz_def,bz_efmax,bz_zval,bz_fsmom,
     . bz_semsh,zbak,bz_dosw,bz_lcond,bz_range, bz_lio1,bz_lio2,ctrl_lmet2, 
     . ctrl_lmet4,ctrl_lmet8,ctrl_noinv,bz_dosmax, !ctrl_lqp2,
!! ... Ewald
     . lat_as,lat_tol,lat_rpad, lat_nkdmx,
!! ... STR
     . str_rmax, str_mxnbr,
!! ... DYN
     . sdmod, mdprm,lat_defm,sdprm,
c     . prmint,
c     . prmint_ts0,prmint_tol,gd_ct,
c     . move_kt,move_ts,move_tsequ,move_tstot,
c     . prmint_new,lbsprm,lsdyn,
c     . prmint_mi,prmint_mx,prmint_nseq,
c     . gd_modt,gd_nmodt,
c      real(8), parameter:: fs = 20.67098d0, degK = 6.3333d-6   ! defaults for MD
     . nitmv,
!! PW
     . pwmode,npwpad,ncutovl,pwemax,pwemin,oveps,delta_stabilize !PW

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
      type(s_bz):: v_sbz
      type(s_ctrl):: v_sctrl
      type(s_ham):: v_sham
      type(s_pot):: v_spot
      type(s_lat):: v_slat
      type(s_mix):: v_smix
      type(s_str):: v_sstr
      type(s_array):: v_sarry
c      type(s_move):: v_smove
c      type(s_tb):: v_stb
      type(s_spec),allocatable:: v_sspec(:)
      type(s_site),allocatable:: v_ssite(:)

      integer:: i_spec

      character  prgnam*(*) !, sstrn*(*)
      character*(*)::sstrnmix,sstrnsymg
c      character(6):: vrsion(2)
c      double precision vn(2),vn2(2)
C ... Local parameters
      character fileid*64
      integer procid,nproc,master
      logical lgors,cmdopt,bittst,ltmp,ioorbp!,asa
      double precision dval,dglob,xx(n0*2),dgets !,ekap(6)
      integer a2vec,bitand,fopna,getdig,i,is,igets,iprint,
     .iprt,irs(5),isw,ifi,ix(n0*nkap0),j,k,l,lbas,lcd,noinv, !noinv was iand(lqp,1)
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
      integer o,oclabl,ohave,oics,opnu,opp,oqnu,osgw,osoptc,
     .oves,owk !osordn,

      real(8):: pnux(20)
      integer:: nnn
      integer:: i_copy_size,i_spacks,iendx,inix,i_spackv
      real(8):: seref
      integer:: ib , ispec
      integer,allocatable:: wowk(:)
      logical:: isanrg,l_dummy_isanrg
      integer:: lmxcg,lmxcy,lnjcg,lnxcg,nlm

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

!! --- Read input parameters from contents of rcd ---
      call readctrl(prgnam) !,vrsion(1),vn(1))
      if (io_help > 0) then
        call readctrlpq(prgnam,nclasp,nl,nsp,pnuc,qnuc,pp,zc,
     .  ves,initc,ics,clabl)
        call cexit(0,1)
      endif

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
c      asa = .false.
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
c      sstrn = ' '

C --- Allocate and copy input to sbz ---
      call ubz_init(v_sbz)
      v_sbz%size=ubz_size()
      v_sbz%def=bz_def
      i_copy_size=size(v_sbz%dosw)
      call dcopy(i_copy_size,bz_dosw,1,v_sbz%dosw,1)
      v_sbz%dosmax=bz_dosmax
      v_sbz%ef=bz_ef
      v_sbz%efmax=bz_efmax
      v_sbz%fsmom=bz_fsmom
      v_sbz%fsmommethod=bz_fsmommethod

      i = 1*isw(bz_lio1)+2*isw(bz_lio2) !+8*isw(bz_lio8)
      i_copy_size=size(v_sbz%lcond)
      call dcopy(i_copy_size,bz_lcond,1,v_sbz%lcond,1)
c      v_sbz%lio=i
      v_sbz%lmet=bz_lmet
      v_sbz%lmull=bz_lmull
c      v_sbz%lopt=0


c      i_copy_size=size(v_sbz%lshft)
c      call ii8copy(i_copy_size,bz_lshft,1,v_sbz%lshft,1)
      v_sbz%lshft = bz_lshft
      v_sbz%n     = bz_n
      v_sbz%ndos  = bz_ndos
      v_sbz%nevmx = bz_nevmx
c      i_copy_size = size(v_sbz%nkabc)
c      call ii8copy(i_copy_size,bz_nabc,1,v_sbz%nkabc,1)
      v_sbz%nkabc = bz_nabc

      v_sbz%range=bz_range
c      i_copy_size=size(v_sbz%semsh)
c      call dcopy(i_copy_size,bz_semsh,1,v_sbz%semsh,1)
      v_sbz%w=bz_w
c      v_sbz%zval=bz_zval


C --- Allocate and copy input to sctrl ---
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
c      k = 2*isw(lcd2)+4*isw(lcd4) !+ 1*isw(lcd1)   +8*isw(lcd8)+64*isw(lcd64)

      i_copy_size=size(v_sctrl%defm)
      call dcopy(i_copy_size,lat_defm,1,v_sctrl%defm,1)
      v_sctrl%lbas=j
c      v_sctrl%lcd=k
      v_sctrl%ldos=ctrl_ldos
      v_sctrl%lfp=lfp
      v_sctrl%lfrce=ctrl_lfrce

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
c      lncol= 4*isw(lncol4)+32*isw(lncol32)+64*isw(lncol64)
      v_sctrl%lmet=lmet
c      v_sctrl%lncol=lncol

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
     . + 8*getdig(irs(1),1,100)
      lrs = 1*mod(irs(1),8)+8*irs(2)+32*irs(3)+64*irs(4)+128*irs(5)
     . + 256*mod(irs(1)/8,2)
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

c        i_copy_size=size(v_sctrl%lpgf)
c        call ii8copy(i_copy_size,ctrl_lpgf,1,v_sctrl%lpgf,1)
c      v_sctrl%lqp=lqp
      v_sctrl%noinv=noinv
c      v_sctrl%lqp2=lqp2
      v_sctrl%lrel=lrel
      v_sctrl%lrs=lrs
c      v_sctrl%lscr=lscr

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
c      v_sctrl%lstr=0
      v_sctrl%lves=lves

C      lxcf   parameter defining XC functional
C      1s digit:
C      1 for Ceperly-Alder (VWN)
C      2 for Barth-Hedin (ASW fit)
C      103 for PBE-GGA
      v_sctrl%lxcf= ham_lxcf

!! ==  Set modep ===
c      ix(1:3) = 2
c      if (prgnam == 'LMPG') ix(3) = 0
c      if (prgnam == 'LMMC') ix(1:3) = 0

      v_sctrl%maxit=iter_maxit
      i_copy_size=size(v_sctrl%mdprm)
      call dcopy(i_copy_size,mdprm,1,v_sctrl%mdprm,1)
      i_copy_size=size(v_sctrl%modep)
c      call ii8copy(i_copy_size,ix(1:3),1,v_sctrl%modep,1)
      v_sctrl%modep=2
      v_sctrl%nbas=nbas

C      Reset nl
      if (nl .ne. max(lmxbx,lmxax)+1 .and. io_help == 0) then
        call info2(20,1,0,' rdctrl: reset global max nl from %i to %i',
     .  nl,max(lmxbx,lmxax)+1)
        nl = max(lmxbx,lmxax)+1
      endif

      v_sctrl%nbasp=nbasp
c      i_copy_size=size(v_sctrl%nesabc)
c      call ii8copy(i_copy_size,nesabc,1,v_sctrl%nesabc,1)
      v_sctrl%nesabc=nesabc
      v_sctrl%nitmv=nitmv
      v_sctrl%nl=nl
      v_sctrl%nsite=nsite
      v_sctrl%nspec=nspec

      v_sctrl%nspin=nsp
      v_sctrl%nvario=nvario
c      i_copy_size=size(v_sctrl%omax1)
c      call dcopy(i_copy_size,omax1,1,v_sctrl%omax1,1)
c      i_copy_size=size(v_sctrl%omax2)
c      call dcopy(i_copy_size,omax2,1,v_sctrl%omax2,1)
      v_sctrl%omax1 = omax1
      v_sctrl%omax2 = omax2
      v_sctrl%quit  = quit
      v_sctrl%rmaxes= rmaxes
      v_sctrl%rmines= rmines
      v_sctrl%sclwsr= sclwsr
      v_sctrl%sdmod = sdmod
c      i_copy_size=size(v_sctrl%sdprm)
c      call dcopy(i_copy_size,sdprm,1,v_sctrl%sdprm,1)
c      i_copy_size=size(v_sctrl%smalit)
c      call ii8copy(i_copy_size,smalit,1,v_sctrl%smalit,1)
c      i_copy_size=size(v_sctrl%tol)
c      call dcopy(i_copy_size,ctrl_tol,1,v_sctrl%tol,1)
c      i_copy_size=size(v_sctrl%zbak)
c      call dcopy(i_copy_size,zbak,1,v_sctrl%zbak,1)
      v_sctrl%sdprm = sdprm
      v_sctrl%smalit= smalit
      v_sctrl%tol   = ctrl_tol
      v_sctrl%wsrmax= wsrmax
      v_sctrl%zbak  = zbak
      v_sctrl%pfloat= lpfloat

C --- Allocate and copy input to slat ---
      v_slat%size=ulat_size()
      if (dalat == NULLR) dalat=0
      v_slat%alat=alat+dalat
      v_slat%as=lat_as
      v_slat%avw=avw
      v_slat%nkdmx=lat_nkdmx
      v_slat%nkqmx=lat_nkdmx
      v_slat%gam =lat_gam
      v_slat%gmax=lat_gmax
c      i_copy_size=size(v_slat%gam)
c      call dcopy(i_copy_size,lat_gam,1,v_slat%gam,1)
c      i_copy_size=size(v_slat%nabc)
c      call ii8copy(i_copy_size,ftmesh,1,v_slat%nabc,1)
c takao I think gam is alreay setted above.
c       i_copy_size=size(v_slat%gam)
c       call dcopy(i_copy_size,lat_gam,1,v_slat%gam,1)
c      i_copy_size=size(v_slat%dist)
c      call dcopy(i_copy_size,lat_dist,1,v_slat%dist,1)
c      i_copy_size=size(v_slat%plat)
c      call dcopy(i_copy_size,plat,1,v_slat%plat,1)
c       i_copy_size=size(v_slat%plat2)
c       call dcopy(i_copy_size,slat_plat2,1,v_slat%plat2,1)
c      i_copy_size=size(v_slat%slat)
c      call dcopy(i_copy_size,lat_slat,1,v_slat%slat,1)
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
      allocate(v_slat%rv_a_ocy(abs(nlm)))
      allocate(v_slat%rv_a_ocg(abs(lnjcg)))
      allocate(v_slat%iv_a_ojcg(abs(lnjcg)))
      allocate(v_slat%iv_a_oidxcg(abs(lnxcg)))
      call sylmnc ( v_slat%rv_a_ocy , lmxcy )
      call scg ( lmxcg , v_slat%rv_a_ocg , v_slat%iv_a_oidxcg , v_slat%iv_a_ojcg )
!!

C --- Allocate and copy input to sham ---
      v_sham%size=uham_size()
      v_sham%alfsi=alfsi
      i_copy_size=size(v_sham%dabc)
      call dcopy(i_copy_size,dabc,1,v_sham%dabc,1)
      v_sham%elind=elind
c      v_sham%nmto=nmto
      i_copy_size=size(v_sham%kmto)
      call dcopy(i_copy_size,kmto,1,v_sham%kmto,1)

C     Replicate ctrl->lncol in ham->lncol, ditto for lham,lgen3
c      v_sham%lncol=(int(v_sctrl%lncol))

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
      v_sham%scaledsigma=scaledsigma

C     Parameters for APW

      v_sham%pwmode=pwmode
      v_sham%npwpad=npwpad
      v_sham%pwemin=pwemin
      v_sham%pwemax=pwemax


      v_sham%oveps=oveps
      v_sham%delta_stabilize=delta_stabilize !takao sep2010

C ... Allocate and initialize sarray
      v_sarry%size=uarray_size()

C ... Allocate and initialize spot
      v_spot%size=upot_size()
      v_spot%vmtz0=vmtz

C --- Allocate and copy input to smix ---
      v_smix%size=umix_size()

ckino's fix for SR11000 size() return's integer(8), right?
c      nnn = size(smix)
c      nnn = umix_size()
c      call dcopy(nnn,vmix,1,v_smix,1)
      v_smix=vmix

c$$$C --- Allocate and copy input to smove ---
c$$$      v_smove%size=umove_size()
c$$$
c$$$      if (lbsprm) then          !Load Bulirsch-Stoer parameters into structure
c$$$        prmint(2) = isw(prmint_new)
c$$$        prmint(3) = prmint_ts0
c$$$        prmint(4) = prmint_tol
c$$$        prmint(5) = prmint_mx
c$$$        prmint(6) = prmint_mi
c$$$        prmint(7:6+prmint_mi) = prmint_nseq(1:prmint_mi)
c$$$      endif
c$$$
c$$$      v_smove%gyro=2d0
c$$$      v_smove%nmodt=gd_nmodt
c$$$      i_copy_size=size(v_smove%modt)
c$$$      call ii8copy(i_copy_size,gd_modt,1,v_smove%modt,1)
c$$$      i_copy_size=size(v_smove%ct)
c$$$      call dcopy(i_copy_size,gd_ct,1,v_smove%ct,1)
c$$$      i_copy_size=size(v_smove%prmint)
c$$$      call dcopy(i_copy_size,prmint,1,v_smove%prmint,1)
c$$$
c$$$      v_smove%kt=move_kt
c$$$      v_smove%ts=move_ts
c$$$      v_smove%tsequ=move_tsequ
c$$$      v_smove%tstot=move_tstot

C --- Allocate and copy input to sstr ---
      v_sstr%size=ustr_size()
      v_sstr%mxnbr=str_mxnbr
      v_sstr%rmax=str_rmax
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
c          call spackv_array_copy_r8_r8('p',v_sspec(i_spackv)%eh3,i_copy_size,1,-0.5d0)

c        i_copy_size=1;
c        do i_spackv=1,nspec
c          call spackv_array_copy_r8_r8('p',v_sspec(i_spackv)%etf,i_copy_size,1,-1d0)
c        enddo
c        i_copy_size=1;
c        do i_spackv=1,nspec
c          call spackv_array_copy_i8_i('p',v_sspec(i_spackv)%lmxf,i_copy_size,1,2*nl-2)
c        enddo
c        i_copy_size=1;
c        do i_spackv=1,nspec
c          call spackv_array_copy_i8_i('p',v_sspec(i_spackv)%norp,i_copy_size,1,2)
c        enddo
c        i_copy_size=1;
c        do i_spackv=1,nspec
c          call spackv_array_copy_r8_r8('p',v_sspec(i_spackv)%vmtz,i_copy_size,1,-0.5d0)
c        enddo

        do  j = 1, nspec
          slabl_(j) = slabl(j)
c          do i_spacks=j,j
c            call spacks_copy('p',v_sspec(i_spacks)%pb1,j,j,pb1(j),i_spacks)
c            call spacks_copy('p',v_sspec(i_spacks)%pb2,j,j,pb2(j),i_spacks)
c            call spacks_copy('p',v_sspec(i_spacks)%coreh,j,j,coreh(j),i_spacks)
c            call spacks_copy('p',v_sspec(i_spacks)%name,j,j,slabl(j),i_spacks)
c          enddo
          v_sspec(j)%pb1=pb1(j)
          v_sspec(j)%pb2=pb2(j)
          v_sspec(j)%coreh=coreh(j)
          v_sspec(j)%name=slabl(j)

          v_sspec(j)%nmcore=nmcore(j)
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

c          i_copy_size=size(v_sspec(j)%idmod)
c          call ii8copy(i_copy_size,idmod(1,j),1,v_sspec(j)%idmod,1)
          v_sspec(j)%idmod=idmod(:,j)
          i_copy_size=size(v_sspec(j)%ehvl)
          call dcopy(i_copy_size,ehvl(1,j),1,v_sspec(j)%ehvl,1)

C       Set idxdn
c          call ivset(ix,1,n0*nkap0,1)
c          call icopy(1+lmxb(j),idxdn(1,j),1,ix,1)
          ix = 1
          ix(1:1+lmxb(j)) = idxdn(1:1+lmxb(j),j)
c          i_copy_size=size(v_sspec(j)%idxdn)
c          call ii8copy(i_copy_size,ix,1,v_sspec(j)%idxdn,1)
          v_sspec(j)%idxdn = ix
c          i_copy_size=size(v_sspec(j)%idu)
c          call ii8copy(i_copy_size,idu(1,j),1,v_sspec(j)%idu,1)
          v_sspec(j)%idu = idu(:,j)
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

          i_copy_size=size(v_sspec(j)%rcfa)
          call dcopy(i_copy_size,rcfa(1,j),1,v_sspec(j)%rcfa,1)
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

C --- Allocate and copy input to ssite ---
C     NB: essential that ssite is LAST array allocated here
CKi      call defrr(v_ssite,1000)
CKi      call usite(v_ssite,-1,0,1,0,0,0)
CKi      call redfrr(v_ssite,nsite*nint(dval(v_ssite,1)))
CKi      call usite(v_ssite,-1,0,nsite,0,0,0)
ckino Jan.04.2012:        allocate(v_ssite(nsite))
        allocate(v_ssite(nsite))
        do i=1,nsite
          v_ssite(i)%size=usite_size()
        enddo
        do  j = 1, nsite
          v_ssite(j)%spec=ips(j)
          v_ssite(j)%class=ips(j)
c          i_copy_size=size(v_ssite(j)%dpole)
c          call dcopy(i_copy_size,dpole(1,j),1,v_ssite(j)%dpole,1)
          v_ssite(j)%dpole=dpole(:,j)
          v_ssite(j)%mpole=mpole(j)
c          do i_spacks=j,j
c            call spacks_copy('p',v_ssite(i_spacks)%clabel,j,j,slabl(ips(j)),i_spacks)
c          enddo
          v_ssite(j)%clabel=slabl(ips(j))
          v_ssite(j)%pl=ipl(j)
          v_ssite(j)%plv=plv(j)
c            i_copy_size=size(v_ssite(j)%pos)
c            call dcopy(i_copy_size,pos(1,j),1,,1)
c            i_copy_size=size(v_ssite(j)%vel)
c            call dcopy(i_copy_size,vel(1,j),1,v_ssite(j)%vel,1)
          v_ssite(j)%pos=pos(1:3,j)
          v_ssite(j)%vel=vel(1:3,j)
          v_ssite(j)%vshft=vshft(j)
c          i_copy_size=size(v_ssite(j)%relax)
c          call ii8copy(i_copy_size,irlx(1,j),1,v_ssite(j)%relax,1)
          v_ssite(j)%relax=irlx(:,j)
c          i_copy_size=size(v_ssite(j)%eula)
c          call dcopy(i_copy_size,eula(1,j),1,v_ssite(j)%eula,1)
          v_ssite(j)%eula=eula(:,j)
          v_ssite(j)%iantiferro=iantiferro(j)
        enddo

C ... Copy string outputs to sstrn : amix, gfopt, jobid, mix, mmham, sxopt, symg
c          j = len_trim(iter_amix)          ! Euler angle mixing amix
c          if (j .gt. 0) then
c            v_sstrn%amix=trim(iter_amix)
c            call packs('strn amix',sstrn,iter_amix)
c            call upacks('strn amix',inix,iendx)
c          endif
          j = len_trim(iter_mix)           ! density mixing mix
          if (j .gt. 0) then
            sstrnmix=trim(iter_mix)
          else
            sstrnmix=''
c            call packs('strn mix',sstrn,iter_mix)
          endif
c          j = len_trim(header)             ! jobid
c          if (j .gt. 0) then
c            v_sstrn%header=trim(header)
c            call packs('strn jobid',sstrn,header)
c          endif

C ... Suppress symmetry operations for special circumstances
          lstsym = 0
c takao comment out next line. !! dec 2012
c          if (lncol .ne. 0) lstsym=1 !lstsym=1: noncollinear case
          !      =2: turn off symops
C     Switches that automatically turn of all symops
          if ((mdprm(1) .ge. 1 .and. mdprm(1) .le. 3) .or.   
!! --pdos mar2003 added. Also in lmv7.F
     .     cmdopt('--cls',5,0,strn) .or. cmdopt('--nosym',7,0,strn).or.cmdopt('--pdos',6,0,strn)) then
          symg = 'e'
          lstsym = 2            !lstsym=2: turn off symops
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
          if (symg .eq. ' ') symg = 'e' ! suppress ops if none explicit
c          lqp = lqp-bitand(lqp,1)+1 ! noinv=1 -->do nothing, noinv=0 --> noinv=1
          noinv=1
        endif
        j = len_trim(symg)               ! Symmetry group symg
        if (j .gt. 0) then
c            call packs('strn symg',sstrn,symg)
        sstrnsymg=trim(symg)
        else
        sstrnsymg=''
        endif
C ... End of copy


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
c          v_sham%lncol=lncol

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
c          if (lncol .ne. 0)
          if(lso.ne. 0) noinv=1 !lqp = lqp-bitand(lqp,1)+1
          v_sctrl%noinv=noinv

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
ckino Jan.04.2012:        call suidx(globalvariables%nkaph,j,nspec,v_sspec)
          call suidx(globalvariables%nkaph,j,nspec,v_sspec)

C ... Set some global variables
          globalvariables%nspec = dble(nspec); globalvariables%l_nspec = globalvariables%l_nspec +1; xx(1) = dble(nspec)
          globalvariables%nbas = dble(nbas); globalvariables%l_nbas = globalvariables%l_nbas +1; xx(1) = dble(nbas)
          globalvariables%nbasp = dble(nbasp); globalvariables%l_nbasp = globalvariables%l_nbasp +1; xx(1) = dble(nbasp)
          globalvariables%nsp = dble(nsp); globalvariables%l_nsp = globalvariables%l_nsp +1; xx(1) = dble(nsp)
          globalvariables%nl = dble(nl); globalvariables%l_nl = globalvariables%l_nl +1; xx(1) = dble(nl)
          globalvariables%lrel = dble(avw); globalvariables%l_avw = globalvariables%l_avw +1; xx(1) = dble(avw)
          globalvariables%lrel = dble(lrel); globalvariables%l_lrel = globalvariables%l_lrel +1; xx(1) = dble(lrel)
          lxcf = int(v_sctrl%lxcf)
          globalvariables%lxcf = dble(lxcf); globalvariables%l_lxcf = globalvariables%l_lxcf +1; xx(1) = dble(lxcf)
          globalvariables%stdo = dble(stdo); globalvariables%l_stdo = globalvariables%l_stdo +1; xx(1) = dble(stdo)
          globalvariables%stdl = dble(stdl); globalvariables%l_stdl = globalvariables%l_stdl +1; xx(1) = dble(stdl)
          globalvariables%stde = dble(stde); globalvariables%l_stde = globalvariables%l_stde +1; xx(1) = dble(stde)
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
c          ix(1) = 2
c          ix(2) = 2
c          ix(3) = 2
c          i_copy_size=size(v_sctrl%modep)
c          call ii8copy(i_copy_size,ix,1,v_sctrl%modep,1)
          v_sctrl%modep=2
C ... Count LDA+U blocks (printout only)
c      call defi(owk,-nbas)
          allocate(wowk(nbas))
          wowk=0
          call pshpr(0)
          call suldau(nbas,v_sspec,v_ssite,nlibu,k,wowk)
          v_sham%nlibu=nlibu
c          v_sham%lmaxu=k
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
c              if (vn(2) .eq. 0) then
c                call awrit6('%a  vn %,2d  verb %i%?#n#,%i#%j#%?#n#,%i#%j#',
c     .      strn,scrwid,-lgunit(k),vn(1)-vn2(1)/1d5,
c     .      iprint(),iprint()-iprt(1),iprt(1),iprt(2)-iprt(1),iprt(2))
c              else
c                call awrit7('%a  vn %,2d('//trim(prgnam)//' %,1d)  '//
c     .      'verb %i%?#n#,%i#%j#%?#n#,%i#%j#',
c     .      strn,scrwid,-lgunit(k),vn(1)-vn2(1)/1d5,vn(2),
c     .      iprint(),iprint()-iprt(1),iprt(1),iprt(2)-iprt(1),iprt(2))
c              endif
              lxcf=v_sctrl%lxcf
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

c$$$c       ltmp = prgnam.eq.'LM'   .or. prgnam.eq.'LMMC'  .or.
c$$$              ltmp = prgnam.eq.'LM'   .or. 
c$$$     .         prgnam.eq.'LMF'  .or. prgnam.eq.'LMFGWD' .or.
c$$$     .         prgnam.eq.'LMGF' .or. prgnam.eq.'LMPG' .or. 
c$$$     .         prgnam.eq.'LMDOS'.or. prgnam.eq.'TBE'
c$$$
c$$$              if ( ltmp ) call awrit7 ( ' bz:%10p' // '%?;n; metal,; nonmetal,;' 
c$$$     .         // '%-1j%?;n>1;%b(%-1j%i),;;' // '%?;n; tetra,;;' // '%?;n; get-qp,;;' 
c$$$     .         // '%?;n; invit,;;' // '%?;n; dens-mat,;;' // '%?;(n>0); %-1jmull=%i,;;' 
c$$$     .         // '%?;n; fixed-spin-mom,;;%b ' , strn , scrwid , lgunit ( k 
c$$$     .         ) , int(v_sbz%lmet) , bitand ( lmet , 2 ) , iand(1,int(v_sbz%lio))
c$$$c    .   , iand(2,int(v_sctrl%lqp)) , iand(8,int(v_sbz%lio)) , int(v_sbz%lmull) 
c$$$     .         , iand(2,int(v_sctrl%lqp)) , .false.                , int(v_sbz%lmull) 
c$$$c     .    , isw ( v_sbz%fsmom .ne.0 ) )
c$$$     .         , isw ( v_sbz%fsmom .ne.NULLR ) )


   90       continue
          endif

C --- Sanity checks and other initialization ---
          if ( iand(8,int(v_sctrl%ldos)) .ne.0 .and. iand(4+2,int(v_sctrl%ldos)) .ne.0 
     .) call rx ( 'inconsistent DOS options' )

C --- takao embed contents in susite here. This is only for lmf and lmfgw.
c      call susite(v_sctrl,v_sham,v_spot,v_sarry,v_slat,v_sspec,v_ssite)
          seref = 0d0
          do ib = 1, nsite
            ispec = v_ssite(ib)%spec
            if (ispec <= 0) cycle
            v_ssite(ispec)%clabel = v_sspec(ispec)%name
            seref = seref +     v_sspec(ispec)%eref
          enddo
          v_sham%seref=seref
          allocate(v_slat%rv_a_opos(abs(3*nsite)))
          if (3*nsite<0) v_slat%rv_a_opos(:)=0.0d0
          if (allocated(v_sarry%iv_a_oips)) deallocate(v_sarry%iv_a_oips)
          allocate(v_sarry%iv_a_oips(abs(nsite)))
          if (nsite<0) v_sarry%iv_a_oips(:)=0

          do i_spackv=1,nsite
c            i_copy_size=3
c            call spackv_array_copy_r8_r8 ( 'u' , v_ssite ( i_spackv ) %pos
c     .   , i_copy_size , i_spackv + 1 - 1 , v_slat%rv_a_opos )
c            i_copy_size=1
c            call spackv_array_copy_i8_i ( 'u' , v_ssite ( i_spackv ) %spec
c     .   , i_copy_size , i_spackv + 1 - 1 , v_sarry%iv_a_oips )
            v_slat%rv_a_opos(3*(i_spackv-1)+1:3*(i_spackv-1)+3)= v_ssite( i_spackv )%pos
            v_sarry%iv_a_oips(i_spackv) =                        v_ssite( i_spackv )%spec
          enddo
          if (procid .eq. master) then
            if (iprint() .ge. 20) then
              if (lstsym .eq. 1) then
c                call upacks('strn symg',i,j)
                write(stdo,357) trim(sstrnsymg) !sstrn(i:j)
  357           format(/' Automatic symmetry finder turned off.  Use: ',a)
              elseif (lstsym .eq. 2) then
                write(stdo,358)
  358           format(/' Symmetry operations suppressed')
              endif
            endif
          endif
          do i_spec=1,nspec
            call mpibc1_s_spec(v_sspec(i_spec),.false.,'rdctrl','sspec')
          enddo
C --- Debugging printout ---
          if (io_help == 0 .and. io_show > 1) then
            print *, '---------- contents of sstrn ------------'
            print *, 'mix: ', trim(sstrnmix)
            print *, 'symg:', trim(sstrnsymg)
            call rx0('done show')
c$$$            call upacks('strn amix',i,j)
c$$$            print *, 'amix:', sstrn(i:j)
c$$$            call upacks('strn mix',i,j)
c$$$            print *, 'mix:', sstrn(i:j)
c$$$            call upacks('strn symg',i,j)
c$$$            print *, 'symg:', sstrn(i:j)
c$$$            call rx0('done show')
          endif
          end subroutine rdctrl2
          end module m_rdctrl2_func
