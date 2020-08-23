      module m_rdctrl2_func
      integer,protected:: irs1,irs2,irs3,irs4,irs5,irs11
      real(8),protected::str_mxnbr,str_rmax

      contains
      subroutine rdctrl2(
     i  recrd,recln,nrecs, prgnam, 
     . slabl_,v_sbz,v_sctrl,v_sham,v_spot,v_slat,v_smix,v_sspec,v_ssite, 
     . v_sarry,sstrnmix,sstrnsymg)  !v_sstr,

      use m_gtv,only: gtv_setrcd
      use m_struc_func
      use m_globalvariables

      use m_rdctrl,only:
     &  readctrl,!readctrlpq,
     &  NULLR,n0,nkap0,s_bz,s_ctrl,s_ham,s_pot,s_lat,s_mix,s_array,s_spec,s_site, !s_tb,,s_move,s_str
!!
     &  io_show,io_help,nvario,header,symg,
     &  frzwf,ham_ewald,ctrl_lfrce,ham_lxcf,gga,ftmesh,lrsig, !nmto,
     &  nsp,lrel,lso,ham_udiag, lfp, lat_gmax,tolft,elind,dqval,rsrnge,vmtz,scaledsigma, !kmto,
     &  alfsi,dabc,rsstol,pmin,pmax,
     &  sigp,sigp_emin,sigp_emax,sigp_a,sigp_b,sigp_efit,
     &    sigp_mode,sigp_nmin,sigp_nmax,
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
!! ... Iterations,MIX
     . ctrl_tol,vmix,iter_maxit,iter_mix, !iter_amix,
!! ... BZ
     . bz_nabc, bz_lshft, bz_lmet,bz_n,bz_nevmx,bz_lmull,ctrl_ldos,bz_ndos,bz_fsmommethod,
     . bz_w,bz_ef,bz_def,bz_efmax,bz_zval,bz_fsmom,
     . bz_semsh,zbak,bz_dosw,bz_lcond,bz_range, bz_lio1,bz_lio2,ctrl_lmet2, 
     . ctrl_lmet4,ctrl_lmet8,ctrl_noinv,bz_dosmax, !ctrl_lqp2,
!! ... Ewald
     . lat_as,lat_tol,lat_rpad, lat_nkdmx,
c!! ... STR
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
      integer,intent(in):: recln,nrecs
      character*(recln),intent(in):: recrd(nrecs)
      character,intent(in)::  prgnam*(*) 

      character slabl_(1)*8
      character strn*(recln)
      type(s_bz):: v_sbz
      type(s_ctrl):: v_sctrl
      type(s_ham):: v_sham
      type(s_pot):: v_spot
      type(s_lat):: v_slat
      type(s_mix):: v_smix
c      type(s_str):: v_sstr
      type(s_array):: v_sarry
c      type(s_move):: v_smove
c      type(s_tb):: v_stb
      type(s_spec),allocatable:: v_sspec(:)
      type(s_site),allocatable:: v_ssite(:)

      integer:: i_spec
      character*(*)::sstrnmix,sstrnsymg
c      character(6):: vrsion(2)
c      double precision vn(2),vn2(2)
C ... Local parameters
      character fileid*64
      integer procid,nproc,master
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
!! --- Copy recrd to rcd in m_gtv 
      call gtv_setrcd(recrd,nrecs,recln,stdo,stdl,stde)
!! --- Read input parameters from rcd in m_gtv ---
      call readctrl(prgnam) 
      if (io_help > 0) call cexit(0,1)
!!   For now, LMF => fp; no screening; nfp-style hamiltonian
      lbas = 3 

c!! Optionally read positions from pos file
c      if (cmdopt('--rpos=',7,0,fileid)) then
c        call iopos(.false.,-1,fileid(8:),nbasp,pos)
c      endif

!! --- Allocate and copy input to sbz ---
      call ubz_init(v_sbz)
      v_sbz%size=ubz_size()
      v_sbz%def=bz_def
      v_sbz%dosw=bz_dosw(:)
      v_sbz%dosmax=bz_dosmax
      v_sbz%ef=bz_ef
      v_sbz%efmax=bz_efmax
      v_sbz%fsmom=bz_fsmom
      v_sbz%fsmommethod=bz_fsmommethod
      v_sbz%lcond=bz_lcond(:)
      v_sbz%lmet=bz_lmet
      v_sbz%lmull=bz_lmull
      v_sbz%lshft = bz_lshft
      v_sbz%n     = bz_n
      v_sbz%ndos  = bz_ndos
      v_sbz%nevmx = bz_nevmx
      v_sbz%nkabc = bz_nabc
      v_sbz%range=bz_range
      v_sbz%w=bz_w

!! Allocate and copy input to sctrl ---
      v_sctrl%size=uctrl_size()
!! lbas: 1 Hamiltonian has no screening transformation
!!           2 Hamiltonian is nfp style
!!          16 freeze phi,phidot for all species
!! lbas = 3 = 1+2 in lmf
!! frzwf: freeze basis or not. j= 1+2+16 if frzwf=T
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
      irs(1) = mod(mod(irs(1),10),4) + 4*getdig(irs(1),1,10)
     . + 8*getdig(irs(1),1,100)
      lrs = 1*mod(irs(1),8)+8*irs(2)+32*irs(3)+64*irs(4)+128*irs(5)
     . + 256*mod(irs(1)/8,2)
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

!!  lves    1 take ves as input
      v_sctrl%lves=lves

!!  lxcf   parameter defining XC functional
!!      1s digit:
!!      1 for Ceperly-Alder (VWN)
!!      2 for Barth-Hedin (ASW fit)
!!      103 for PBE-GGA
      v_sctrl%lxcf= ham_lxcf

!! ==  Set modep ===
c      ix(1:3) = 2
c      if (prgnam == 'LMPG') ix(3) = 0
c      if (prgnam == 'LMMC') ix(1:3) = 0

      v_sctrl%maxit=iter_maxit
      v_sctrl%mdprm=mdprm
      v_sctrl%modep=2
      v_sctrl%nbas=nbas

C      Reset nl
      if (nl .ne. max(lmxbx,lmxax)+1 .and. io_help == 0) then
        call info2(20,1,0,' rdctrl: reset global max nl from %i to %i',
     .  nl,max(lmxbx,lmxax)+1)
        nl = max(lmxbx,lmxax)+1
      endif

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
c      i_copy_size=size(v_sham%kmto)
c      call dcopy(i_copy_size,kmto,1,v_sham%kmto,1)

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
          v_sspec(j)%alpha=alpha(:,j)
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
          v_ssite(j)%dpole=dpole(:,j)
          v_ssite(j)%mpole=mpole(j)
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

        j = len_trim(iter_mix)  ! density mixing mix
        if (j .gt. 0) then
          sstrnmix=trim(iter_mix)
        else
          sstrnmix=''
        endif

!! ... Suppress symmetry operations for special circumstances
          lstsym = 0
!!     Switches that automatically turn of all symops
          if ((mdprm(1) .ge. 1 .and. mdprm(1) .le. 3) .or.   
!! --pdos mar2003 added. Also in lmv7.F
     .     cmdopt('--cls',5,0,strn) .or. cmdopt('--nosym',7,0,strn).or.cmdopt('--pdos',6,0,strn)) then
          symg = 'e'
          lstsym = 2            !lstsym=2: turn off symops
        endif
!!     Switches that turn off automatic finder, incl. inversion
        if (lstsym .ne. 0) then
          i = 1
          do while (i .ne. 0)
            i = index(symg,'find')
            if (i .ne. 0) then
              symg(i:i+3) = ' '
            endif
          enddo
          if (symg .eq. ' ') symg = 'e' ! suppress ops if none explicit
!!          lqp = lqp-bitand(lqp,1)+1 ! noinv=1 -->do nothing, noinv=0 --> noinv=1
          noinv=1
        endif
        j = len_trim(symg)               ! Symmetry group symg
        if (j .gt. 0) then
          sstrnsymg=trim(symg)
        else
          sstrnsymg=''
        endif
!!
        if (cmdopt('--rdbasp',8,0,strn)) call rx('not support --rdbasp')

!!     Add dalat to alat
!!  v_slat%alat=(v_slat%alat)+dalat !this is a bug; this should be commented out. This
!!  is a bug for lm7K when I started lm7K. I had included this bug here. Fixed at 28May2010.

!! Dirac equation requires spin polarization
        if( nsp==1 .and. int(v_sctrl%lrel)==2 ) call rx('rdccat: Dirac equation requires NSPIN=2')

!! ... Suppress inversion when noncollinear magnetism, SX, NLO
          if(lso.ne. 0) noinv=1 !lqp = lqp-bitand(lqp,1)+1
          v_sctrl%noinv=noinv
!! ... Setup for idxdn ... ctrl->lham,4 is automatic downfolding switch
          j= 2
C     No screening => no downfolding; also lmxb<l<=lmxa => 'high'
C     Probably ought to have lmxb<l<=lmxa => 'high' always

c          if ( iand(1,int(v_sctrl%lbas)) .ne.0 ) j = 3
          j=3
C     nfp-style basis:
c          if ( iand(2,int(v_sctrl%lbas)) .ne.0 ) j = j+10
          j=j+10
          call suidx(globalvariables%nkaph,j,nspec,v_sspec)

!! ... Set some global variables
          globalvariables%nspec = dble(nspec); globalvariables%l_nspec = globalvariables%l_nspec +1; xx(1) = dble(nspec)
          globalvariables%nbas = dble(nbas);   globalvariables%l_nbas = globalvariables%l_nbas +1; xx(1) = dble(nbas)
          globalvariables%nbasp= dble(nbasp);  globalvariables%l_nbasp = globalvariables%l_nbasp +1; xx(1) = dble(nbasp)
          globalvariables%nsp =  dble(nsp);    globalvariables%l_nsp = globalvariables%l_nsp +1; xx(1) = dble(nsp)
          globalvariables%nl =   dble(nl);     globalvariables%l_nl = globalvariables%l_nl +1; xx(1) = dble(nl)
          globalvariables%lrel = dble(avw);    globalvariables%l_avw = globalvariables%l_avw +1; xx(1) = dble(avw)
          globalvariables%lrel = dble(lrel);   globalvariables%l_lrel = globalvariables%l_lrel +1; xx(1) = dble(lrel)
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
          v_sctrl%modep=2
C ... Count LDA+U blocks (printout only)
          allocate(wowk(nbas))
          wowk=0
          call pshpr(0)
          call suldau(nbas,v_sspec,v_ssite,nlibu,k,wowk)
          v_sham%nlibu=nlibu
          call poppr
          deallocate(wowk)
C     Free arrays used to read input
          deallocate(pnu,qnu,pz,amom,idmod,rsmh,eh,rsmh2,eh2,pb1,pb2,
     .     lmxpb,qpol,stni,tbvso,iq1,ivso,rg,rsma,rfoca,rsmfa,rcfa,nxi,
     .     exi,rint,rcut,coreq,mass,colxbs,radxbs,rs3,rham,idxdn,hcr,
     .     rmt,alpha,idu,uh,jh,dv,grp,grp2,mxcst1,mxcst2,mxcst4,kmxt,kmxv,
     .     lfoca,eref,lmxl,lxi,coreh,lmxa,lmxb,spec_a,z,nr,rsmv,
     .     pos,vel,eula,vshft,ips,ipl,plv,irlx,mpole,dpole)

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

!! --- takao embed contents in susite here. This is only for lmf and lmfgw.
          seref = 0d0
          do ib = 1, nsite
            ispec = v_ssite(ib)%spec
            if (ispec <= 0) cycle
            v_ssite(ispec)%clabel = v_sspec(ispec)%name
            seref = seref +     v_sspec(ispec)%eref
          enddo
          v_sham%seref=seref
          allocate(v_slat%rv_a_opos(abs(3*nsite)))
          if (allocated(v_sarry%iv_a_oips)) deallocate(v_sarry%iv_a_oips)
          allocate(v_sarry%iv_a_oips(abs(nsite)))
          do i_spackv=1,nsite
            v_slat%rv_a_opos(3*(i_spackv-1)+1:3*(i_spackv-1)+3)= v_ssite( i_spackv )%pos
            v_sarry%iv_a_oips(i_spackv) =                        v_ssite( i_spackv )%spec
          enddo
          if (procid .eq. master) then
            if (iprint() .ge. 20) then
              if (lstsym .eq. 1) then
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
!! Debugging printout ---
          if (io_help == 0 .and. io_show > 1) then
            print *, '---------- contents of sstrn ------------'
            print *, 'mix: ', trim(sstrnmix)
            print *, 'symg:', trim(sstrnsymg)
            call rx0('done show')
          endif
          end subroutine rdctrl2
          end module m_rdctrl2_func
