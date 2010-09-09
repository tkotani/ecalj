Cdwdef1 have_common_w_in readctrlpq , #_of_w_access= 0
      module m_rdctrl
      use m_struc_def
      implicit none

C ... Fixed parameters
      integer,parameter:: NULLI =-99999
      real(8),parameter::    NULLR =-99999
      logical,parameter:: T=.true.,F=.false.
      integer,parameter :: n0=10
      integer,parameter :: nkap0=3
      integer(2):: nono

C ... IO
      integer:: io_show=0,io_help=0,nvario=0

C ... HEADER, SYMGRP
      character(256):: header,symg=' '

C ... HAM
      logical :: frzwf,ham_ewald
      integer:: ctrl_lfrce,ham_lxcf,gga,ftmesh(3),nmto=0,lrsig=0,
     .nsp=1,lrel=1,lso=0,ham_udiag=0,ggax !,ham_qasa
c      logical:: ltbe  ! set to T if TBE program (used for defaults)
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
     .(sigp_a,sigp(6)),(sigp_b,sigp(7)),(sigp_efit,sigp(8))
C      real(8):: sigp_mode,sigp_nmin,sigp_nmax
C      equivalence (sigp_mode,sigp(1)),(sigp_nmin,sigp(2)),
C     .  (sigp_nmax,sigp(4))

C ... OPTIONS
      logical ::
     .lcd2,lcd4=F !,lcd64, lasa32=F,lcd1,lcd8=Flasa32 is dummy not used but wor awrit8 in rdctrl2
C$$$    &  ,lasa4=T,lasa8=F,,lasa64=F,
C$$$     .  lasa512
c     &  ,lham1,lham4=F,lham8=F,lham16=F,lham32=F,lham64,lham128=F,
c     .  lham256 !,lves takao2010apr lves is integer
      integer ::lves
C     Initialize, since read is only executed when NSPIN=2
      logical:: lncol1=F,lncol2=F,lncol8=F,lncol16=F
c     .  lncol32=F,lncol64=F
      logical:: lncol4=F, lncol32=F,lncol64=F

c      character(256):: sxopt=' '
      integer:: lscr,smalit(2),nesabc(3), !lasa3=1,lsx,
     .lstonr(3)=0,quit,nl,lpfloat
      real(8):: rmines,rmaxes,ham_qss(4) !asa_elin,,lham3=0

C ... STRUC
C     Initializing nbas and alat flags no nbas or alat has been input
      real(8):: lat_slat(9),dlat,alat=NULLR,plat(9),lat_gam(4),dalat
      real(8):: vol,avw,lat_dist(9)
      integer:: nbas=NULLI,nbasp=NULLI,nsite=NULLI,nspec,nclass
      integer:: lat_ldist=0

C ... SPEC
      real(8):: omax1(3),omax2(3),wsrmax,sclwsr
      integer,parameter ::mxspec=256
      character*8 slabl(mxspec)
      integer:: lmxbx=-1,lmxax,nkaph
      logical,allocatable:: mxcst1(:),mxcst2(:),mxcst4(:)
      integer,allocatable:: idxdn(:,:), grp(:),grp2(:),
     .idu(:,:),lmxb(:),lmxa(:),idmod(:,:),iq1(:,:),ivso(:,:),
     .kmxt(:),kmxv(:),lfoca(:),lmxl(:),lxi(:),nxi(:),nr(:),lmxpb(:)
      real(8),allocatable:: rsmh(:,:),rsmh2(:,:),eh(:,:),eh2(:,:),
     .hcr(:,:),rs3(:),rham(:),alpha(:,:),ehvl(:,:),
     .dv(:), uh(:,:),jh(:,:),
     .qpol(:,:),stni(:),tbvso(:,:),
     .pnu(:,:,:),qnu(:,:,:),
     .coreq(:,:),mass(:),colxbs(:,:),radxbs(:),
     .rg(:),rsma(:),rfoca(:),
     .rsmfa(:),rcfa(:,:),
     .exi(:,:),rint(:),rcut(:),rmt(:),pz(:,:,:),
     .amom(:,:),spec_a(:),z(:),eref(:),rsmv(:)
      character*(8),allocatable::pb1(:),pb2(:), coreh(:)

C ... SITE
      character*8:: alabl
      real(8),allocatable   :: pos(:,:),vel(:,:),eula(:,:),vshft(:)
      integer,allocatable:: ips(:),ipl(:),plv(:),irlx(:,:),ndelta(:)
      real(8),allocatable   :: delta(:,:),mpole(:),dpole(:,:)

C ... Iterations
      character(128) :: iter_mix=' ',iter_amix=' '
      real(8):: ctrl_tol(3)=1d-4
      type (s_mix):: vmix
      integer:: iter_maxit=1

C ... BZ
      integer:: bz_nabc(3)=NULLI,bz_lshft(3)=0,
     .bz_lmet,bz_n,bz_nevmx,bz_lmull,ctrl_ldos,bz_ndos
      real(8):: bz_w,bz_ef,bz_def,bz_efmax,bz_zval,bz_fsmom,
     .bz_semsh(10),zbak(2),bz_dosw(2),bz_lcond(4),bz_range=5d0
      logical:: bz_lio1,bz_lio2,ctrl_lmet2, !,bz_lio8
     .ctrl_lmet4=.true.,ctrl_lmet8=.true.,ctrl_lqp1,ctrl_lqp2

C$$$C ... TB
C$$$      integer:: tb_fmode=NULLI,tb_nbfit(2)=0
C$$$      real(8):: tb_ebfit(2)=0d0,tb_rmfit(2)=0d0,tb_alam=NULLR,
C$$$     .  tb_alsc=NULLR
C$$$      logical:: ltb1,ltb2,ltb4,ltb8,ltb16,
C$$$     .  ltb32,ltb64,ltb128,ltb256,ltb512,
C$$$     .  ltb1024,ltb2048,ltb4096,ltb213,ltb214,
C$$$     .  ltb215,ltb216,ltb217,ltb218

C ... Ewald
      real(8):: lat_as,lat_tol,lat_rpad
      integer:: lat_nkdmx

C ... STR
      real(8):: str_rmax=nullr !, str_rmaxg=nullr,
c     .  str_drwats=NULLR !,iinv_tol                     !str_kaps(6)=NULLR,
      integer:: str_mxnbr !,str_lmaxw=NULLI !str_mode=nulli,str_nkaps=1,
c      integer:: iinv_nit,iinv_ncut,
c      integer:: str_ivl=0  !tcf_nalf,
c      logical:: str_lshow1 !,str_lequiv1

C ... TCF
c      integer tcf_nbisi(3),tcf_ncupl,tcf_ndust
c      real(8):: tcf_wztcf !tcf_adec,

C ... CGF
c      integer:: ctrl_lcgf

C ... PGF
c      integer:: ctrl_lpgf(2)=0
c      character(256):: gfopt=' '
c      real(8):: platl(3)=0d0,platr(3)=0d0

C ... DYN
C   structure of mdprm:
C   arg 1: 0 no relaxation or dynamics
C          1 new dynamics 2  restart dynamics
C          4 relax with conjugate gradients
C          5 relax with variable metric
C          6 relax with Broyden
C   arg 2: statics: switch
C          1 read hessian matrix
C          dynamics:
C            number of iterations between printouts.
C   arg 3: (stat) relaxation x-tolerance
C          (dyn)  temperature
C   arg 4: (stat) relaxation g-tolerance
C          (dyn)  time step
C   arg 5: (stat) step length
C          (dyn)  relaxation time
C   arg 6: (stat) Remove hessian after this many steps
C          (dyn)  --
C   Structure of sdprm :parameters governing dynamics (see magtrq.f)
C   arg 1: scale factor amplifying magnetic forces
C   arg 2: time step tau for Landau dynamics
C   arg 3: reference energy etot0 for Nose thermostat
C          (not used unless fscal is zero)
C   arg 4: maximum allowed change in angle (not used unless nonneg)
C   arg 5: etol: set tau=0 this iter if etot-ehf>etol

C  prmint  Parameters for numerical integration
C          For Bulirsch-Stoer integration:
C  arg  1:   mode: (1 for BS-integration)
C  arg  2:   rst:  1 for start new integration, 0 to continue
C  arg  3:   ts0:  minimum time step size
C  arg  4:   tol:  tolerance in integration errors
C  arg  5:   mx:   order of rational function extrapolation
C  arg  6:   mi:   number of midpoint rules
C  arg  7-17 nseq: sequence of no. midpoint divisions
C  arg  18:        offset to bs workspace
      integer:: sdmod=-1
      real(8):: mdprm(6),lat_defm(6),sdprm(5)
      real(8):: prmint(20)=(/1d0,1d0,0d0,2d-4,7d0,11d0,2d0,4d0,6d0,8d0,
     .12d0,16d0,24d0,32d0,48d0,64d0,96d0,0d0,0d0,0d0/),
     .prmint_ts0,prmint_tol,gd_ct(3),
     .move_kt,move_ts,move_tsequ,move_tstot
      logical:: prmint_new,lbsprm=F,lsdyn=F
      integer:: nitmv,prmint_mi,prmint_mx,prmint_nseq(11),
     .gd_modt(3),gd_nmodt=1
      character(1024):: mmham=' '
      real(8), parameter:: fs = 20.67098d0, degK = 6.3333d-6   ! defaults for MD

cC ... GW
c      integer:: gw_mksig=0,gw_nabc(3)=NULLI,gw_nband=-NULLI,
c     .  gw_nime=6
c      real(8):: gw_gcutb=NULLR,gw_gcutx=NULLR,gw_ecuts=NULLR,
c     .  gw_delre(2)=NULLR,gw_deltax=NULLR,gw_deltaw=NULLR,
c     .  gw_gsmear=NULLR,gw_pbtol=NULLR,gw_qoffp=NULLR

C ... PLANE
c      real(8):: slat_plat2(9)

C ... ORDERN
C      integer lordn
C      real(8):: plate(9)    !slat, 'lat plate'
C      real(8):: efre,ordn_rmaxg

C ... OPTICS
c      integer:: optic_ne,optic_ocrng(2),optic_unrng(2),ctrl_loptc,
c     .  optic_nchi2=0,optic_axes(18)
c      real(8):: optic_window(2),optic_esciss
c      logical :: optic_mode1

C ... PW takao
      integer:: pwmode,npwpad,ncutovl !ncutovl is by takao. not in lm-7.0beta.
      real(8)::    pwemax,pwemin,oveps
      character(1)::mmm
ccccccccccccccccccccccccccccccccccccccccc

Cr         1   amix    amix mode
cccCr         2   gemb    file names for embedding GF
cccCr         3   gfopt   GF options
Cr         4   jobid   job id
cccCr         5   map     mapping
Cr         6   mix     mixing mode
cccCr         7   mmham   micromagnetics hamiltonian
cccCr         8   sxopt   SX options
Cr         9   symg    symmetry group

ctest
c      real(8):: HAM_FORCES, MIX_MODE
c      real(8):: HAM_SIGP_b, MIX_b, MIX_k
c      real(8),allocatable:: SPEC_SITE(:,:)
c      real(8):: MIX_MODE_b,
c     &          SPEC_RSMH
c      real(8):: HAM_SIGP(8),z(100),r(100),rw(100)=-9999,ra(100)=-9999
c      integer :: VERBOS(5),nclasstest, n_verbos
c      character(16) :: atom(50)
!
      contains
      subroutine readctrl(prgn,vstrn,vn)
C- Reads data parsed from ctrl file into m_rdctrl
C ----------------------------------------------------------------------
Ci Inputs
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
      use m_globalvariables
C     implicit none
C ... Passed parameters
      character*(*):: prgn
      character(6):: vstrn(2)
      double precision vn(2)
C ... Local variables
      integer::nout,i,j,nn,i0,iprint,lgunit,ivec(10),iosite
      integer:: io_iactive,io_tim(2),verbos(5),isw
      character(256)::  a,outs
      logical :: cmdopt
      integer,parameter :: master=0
      integer:: procid,nproc
      logical::  mlog=.false.
      integer:: lmxbj,lmxaj,nlbj,nlbji,nlaj,nlaji
      integer:: mpipid,a2vec
      double precision vsn,vers,xv(2*n0)
      character(256*16) :: bigstr=' '
      integer:: it
      logical :: parmxp
      logical,parameter:: T=.true.,F=.false.
      logical::  debug=.false.

      integer:: lp1,lpzi,nlmax,nrmx=1501
      integer,parameter:: mxspec=256
      real(8):: xx
      type (s_mix) :: vamix
      real(8)::  avwsr
      integer::izerv(n0)=(/(0,i=1,n0)/)
      real(8)::   zerov(n0)=(/(0d0,i=1,n0)/)

      real(8):: tbalph(5,5)
      integer:: ii,sw
      real(8):: dasum,dglob
      character(128) :: nm
      integer :: stdo,stdl,stde

      real(8):: nullrv(256)
      integer:: nulliv(256),jj(2) !,nkapsi
      logical:: ltmp,noelind

      data tbalph/.214286D0,.000001D0,.000001D0,.000001D0,.001D-3,
     ..287234D0,.025822D0,.000001D0,.000001D0,.001D-3,
     ..348485D0,.053030D0,.010714D0,.000001D0,.001D-3,
     ..385057D0,.073209D0,.022481D0,.006069D0,.001D-3,
     ..404761D0,.087927D0,.032886D0,.012257D0,4.2735D-3/

      integer:: levelinit=0

C ... For MPI
      procid = mpipid(1)
      nproc  = mpipid(0)
C     mlog = cmdopt('--mlog',6,0,a)

C     Used to set tbe-dependent or f-dependent defaults
c      ltbe = prgn .eq. 'TBE'
      lfp = 0
c      if (prgn=='LMF' .or. prgn=='LMFGWD' .or. prgn=='LMMC') lfp = 1
      if (prgn=='LMF' .or. prgn=='LMFGWD') lfp = 1

C      print *, '!! choose default nl? ... NL now only sets defaults'
C      print *, '!! molecules stuff around LMXB not done yet '
C      print *, '!! fix EMESH'
C      print *, '!! change default xcf?'

C --- Initialize ---
      nullrv = nullr
      nulliv  =nulli
      debug = cmdopt('--debug',6,0,a)
      call toksw_init(debug)
      stdo = lgunit(1)
      stdl = lgunit(2)
      stde = lgunit(1)
C     call gtv_setst(debug,stdo,stdl,stde)
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

C      goto 1234

C --- IO ---
      nm='IO_SHOW'; call gtv(trim(nm),tksw(prgn,nm),io_show,def_i4=
     .0,note='Echo data as it is read from input file')
      if (cmdopt('--show',6,0,a)) io_show = 1
      if (procid .ne. master) io_show = 0
      if (cmdopt('--show=',7,0,a)) then
        i = 7
        i = a2vec(a,len(a),i,2,' ',1,1,1,it,io_show)
      endif
      call gtv_setio(debug,io_show,io_help) ! In case io_show changed

      nm='IO_HELP'; call gtv(trim(nm),tksw(prgn,nm),io_help,
     .def_i4=0, note='Show what input would be sought, '//
     .'without attempting to read data')
      if (cmdopt('--input',7,0,a)) io_help = 1  !optio=0 in old code
      if (io_help == 1) io_show = 1
      call gtv_setio(debug,io_show,io_help) ! In case io_help changed

C --- Version ---
      do  j = 1, 2
        vsn = dble(int(vn(j)*100))/100
        if (vsn .ne. 0d0) then
          nm = 'VERS_'//vstrn(j)
          outs = 'Input style'
          if (j .eq. 2) outs='Program'
          sw = tksw(prgn,nm)
          if (sw .ne. 2) then
            call gtv(trim(nm),sw,vers,note=trim(outs)//
     .      ' version check: integer part must be at least as large as:')
            if (io_help .eq. 1) then
              call info2(1,0,0,'   '//trim(outs)//' version: %d',vsn,0)
            else
              if (io_show .ne. 0) call info2(1,0,0,'  Expected '//
     .        trim(outs)//' version: %d'//' ... found in file: %d',vsn,
     .        vers)
              if ( vers .lt. int(vsn)) call fexit(-1,1,
     .        '  Exit -1 RDCTRL: control file earlier than v%d',vsn)
            endif
          endif
        endif
      enddo

C --- CMD ---
C Token ordering: CMD
c  250 continue
      nm='CMD'; call gtv(trim(nm),tksw(prgn,nm), bigstr, note=
     .'Contents are appended to command-line arguments',nout=nout)
      if (nout .eq. 1) then
        call acmdop(bigstr,len_trim(bigstr),0)
      endif

C --- CONST ---
      call numsyv(nvario)
      nm='CONST'; call gtv(trim(nm),tksw(prgn,nm), bigstr, note=
     .'Constants may declared for use in expressions'//
     .'%N   Variables may also be set from the command-line:  -vnam=#',
     .nout=nout)
      if (nout .eq. 1) then
        i = 0
        call parsyv(bigstr,len_trim(bigstr),1999,0,i)
        if (io_show .ne. 0) call shosyv(0,0,0,stdo)
      endif

C --- Header ---
C Token ordering: HEADER
C   10 continue
      nm='HEADER'; call gtv(trim(nm),tksw(prgn,nm),header,note=
     .'Contents displayed at beginning of program execution',
     .nout=nout)
      if (nout .gt. 0) then
        if (stdl>0 .and. iprint()>0) call headl2(prgn,0,stdl)
        outs  = 'HEADER '//trim(header)
        if (procid .eq. master) then
          if (io_show .ne. 0) write(stdo,'(1x)')
          write(stdo,'(1x,a)') trim(outs)
          if (io_show .ne. 0) write(stdo,'(1x)')
          if (stdl>0 .and. iprint()>0) write(stdl,'(1x,a)') trim(outs)
        endif
      endif

C --- IO ---
      nm='IO_VERBOS'; call gtv(trim(nm),tksw(prgn,nm),verbos,
     .note='Verbosity stack for printout'//
     .'%N   May also be set from the command-line: --pr#1[,#2]',
     .def_i4v=(/30/),nout=i0)
C ... Override verbos w/ -pr commmand-line arg
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

c      ltmp = F
c      nm='IO_WKP'; call gtv(trim(nm),tksw(prgn,nm),ltmp,def_lg=F,
c     .note='Print size of dynamic arrays when allocated')
c      if (ltmp) call wkprnt(1)
      nm='IO_IACTIV'; call gtv(trim(nm),tksw(prgn,nm),
     .io_iactive, note='Turn on interactive mode'//
     .'%N   May also be controlled from the command-line:'//
     .'  --iactiv  or  --iactiv=no',
     .def_i4=0)
      if( cmdopt('--no-iactiv ',7,0,a)) io_iactive = 0
      if( cmdopt('--iactiv',7,0,a))     io_iactive = 1
      if( cmdopt('--iactiv=no',11,0,a)) io_iactive = 0
      if (io_iactive .gt. 1) io_iactive = 1
      call initqu(io_iactive)

 1234 continue

      nm='IO_TIM'; call gtv(trim(nm),tksw(prgn,nm),io_tim,
     .note='Turns CPU timing log.  Value sets tree depth.'//
     .'%N   Optional 2nd arg prints CPU times as routines execute.'//
     .'%N   Args may be set through command-line: --time=#1[,#2]',
     .def_i4v=(/0,0/),nmin=2,nout=i0)
C     Override with '--time=' commmand-line arg
      if ( cmdopt('--time=',7,0,a) ) then
        i = 7
        i = a2vec(a,len(a),i,2,', ',2,2,2,it,io_tim)
        if (i .lt. 0) call rxs('error parsing switch ',a)
        i0 = max(i0,i)
      endif
      if ( i0 >=1 ) call tcinit(io_tim(2),io_tim(1),levelinit)

C --- Struc ---
      if (tksw(prgn,'STRUC') == 2) goto 59
      if (io_show+io_help/=0) call info0(1,1,
     .0,' --- Parameters for crystal structure ---')
      nm='STRUC_FILE'; call gtv(trim(nm),tksw(prgn,nm),outs,nmin=10,
     .nout=nout, note=
     .'Name of site file containing basis and lattice information.'//
     .'%N%3fRead NBAS, PLAT, and optionally ALAT from site file, '//
     .'if specified.'//
     .'%N%3fOtherwise, they are read from the ctrl file.')
      if (nout == 1) then
        xx = nullr
        if (procid .eq. master) then
          j = iosite(7000,3d0,0,trim(outs),i,slabl,alat,plat,nbas,
     .    nspec,xx,xx,xx,xx,xx,xx,xx)
        endif
        call mpibc1(nbas,1,2,mlog,'readctrl','nbas')
        call mpibc1(alat,1,4,mlog,'readctrl','alat')
        call mpibc1(plat,9,4,mlog,'readctrl','plat')
      endif

C     call shosyv(0,0,0,stdo)
      if (alat .eq. NULLR) then
        nm='STRUC_ALAT'; call gtv(trim(nm),tksw(prgn,nm),
     .  alat, note= 'Scaling of lattice vectors, in a.u.')
      endif
      if (nbas .eq. NULLI) then
        nm='STRUC_NBAS'; call gtv(trim(nm),tksw(prgn,nm),
     .  nbas,note='Size of basis')
        nm='STRUC_PLAT'; call gtv(trim(nm),tksw(prgn,nm),
     .  plat, nmin=9, nout=nout, note=
     .  'Primitive lattice vectors, in units of alat')
      endif
      avw = avwsr(plat,alat,vol,nbas)
      nm='STRUC_NSPEC'; call gtv(trim(nm),tksw(prgn,nm),nspec,
     .note='Number of species to read from SPEC category.'//
     .'%N%3fIf not present, NSPEC will be obtained '//
     .'by counting entries in the SPEC category' )
      if (io_help == 0 .and. nspec .eq. NULLI) then

C ... Count number of species in SPEC category
        nm='SPEC_ATOM'; sw = tksw(prgn,nm)
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
     .    ' ... found %i species in SPEC category',nspec,0)
        endif
      endif

C     Extra site positions for e.g. point multipoles
      nm='STRUC_NBASP'; call gtv(trim(nm),tksw(prgn,nm),nbasp,
     .def_i4=nbas,note='nbas + no. of point multipoles')
      if (nbasp .eq. NULLI) nbasp = nbas
      nsite = nbasp

      nm='STRUC_SLAT'; call gtv(trim(nm),tksw(prgn,nm),lat_slat,
     .nmin=9,note='Supercell lattice vectors')
      nm='STRUC_DALAT'; call gtv(trim(nm),tksw(prgn,nm),
     .dalat, def_r8=0d0, note='added to alat after input is read')
      nm='STRUC_NL'; call gtv(trim(nm),tksw(prgn,nm),nl,def_i4=3,
     .note='global default lmax+1 for basis and augmentation')

C ... Lattice distortion or rotation
      sw = tksw(prgn,'STRUC_SHEAR')
      lat_gam(1) = NULLR
      if (sw/=2) then
        j = 0
        nm='STRUC_SHEAR'; call gtv(trim(nm),0,lat_gam,nmin=4,or=T,
     .  note='Volume-conserving shear of PLAT (1=ideal)')
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
     .  note='6 Voigt strain matrix elements')
        if (lat_dist(1) /= NULLR .or. io_help /= 0) then
          nm='STRUC_ALPHA'; call gtv(trim(nm),1,lat_dist(7),nout=nout,
     .    note=
     .    'Amplitude of (Voigt) strain.  Only read if STRAIN input')
          if (nout /= 0) j = 3
        endif
  880   continue
        lat_ldist = j
      endif
      if (lat_gam(1) == NULLR) lat_gam=(/0d0,0d0,1d0,1d0/)

C ... End of struc tokens
   59 continue

C --- Options ---
      if (tksw(prgn,'OPTIONS') == 2) goto 69
      if (io_show+io_help/=0) call info0(1,1,
     .0,' --- Program Options ---')
C      nm='OPTIONS'; call gtv(trim(nm),tksw(prgn,nm),Texist=ltmp,
C     .  note='Various optional input')
      nm='OPTIONS_HF'; call gtv(trim(nm),tksw(prgn,nm),lcd2,
     .def_lg=F,note='T for non-self-consistent Harris')
c      nm='OPTIONS_SHARM'; call gtv(trim(nm),tksw(prgn,nm),lham256,
c     .  def_lg=F, note='Use true spherical harmonics')
c      nm='OPTIONS_FRZ'; call gtv(trim(nm),tksw(prgn,nm),lcd1,
c     .  def_lg=F, note= 'Freeze core')
c      nm='OPTIONS_XCQS'; call gtv(trim(nm),tksw(prgn,nm),lcd64,note=
c     .  'XC potential in reciprocal space')

c      nm='OPTIONS_TPAN'; call gtv(trim(nm),tksw(prgn,nm),lham1,
c     .  def_lg=F,note='two-panel')
c      nm='OPTIONS_SAVVEC'; call gtv(trim(nm),tksw(prgn,nm),lham64,
c     .  def_lg=F,note='Save eigenvectors on disk')

      nm='OPTIONS_RMINES'; call gtv(trim(nm),tksw(prgn,nm),rmines,
     .def_r8=1d0,note='Minimum MT radius when finding new ES')
      nm='OPTIONS_RMAXES'; call gtv(trim(nm),tksw(prgn,nm),rmaxes,
     .def_r8=2d0,note='Maximum MT radius when finding new ES')
      nm='OPTIONS_NESABC'; call gtv(trim(nm),tksw(prgn,nm),nesabc,
     .nmin=3,def_i4v=(/100,100,100/),
     .note='No. divisions when searching for empty spheres')
C      nm='OPTIONS_STONER'; call gtv(trim(nm),tksw(prgn,nm),
C     .  lstonr, note=
C     .  'Generalised Stoner rigid band calculation '//
C     .  '%N%3fSecond argument is number of points; '//
C     .  'third for graphical output')

      outs = ' '
      nm='OPTIONS_Q'; call gtv(trim(nm),tksw(prgn,nm),outs,nmin=10,
     .note='Use Q=show, Q=atom, or Q=band to quit after '//
     .'input, sphere calc or band pass.'//
     .'%N%3fCommand-line `--quit=string'' overrides input file')
      if (outs=='SHOW' .or. outs=='show' .or.
     .cmdopt('--quit=show',11,0,a)) then; quit=1
      elseif (outs=='ATOM' .or. outs=='atom' .or.
     .cmdopt('--quit=atom',11,0,a)) then; quit=2
      elseif (outs=='BAND' .or. outs=='band' .or.
     .cmdopt('--quit=band',11,0,a)) then; quit=4
      elseif (outs=='HAM' .or. outs=='ham' .or.
     .cmdopt('--quit=ham',10,0,a)) then; quit=8
      elseif (outs==' ') then
      else
        call rx('OPTIONS_Q= must contain one of SHOW, ATOM, or BAND')
      endif

      nm='OPTIONS_PFLOAT';call gtv(trim(nm),tksw(prgn,nm),lpfloat,
     .def_i4=0,note=
     .'Controls how band CG is determined in floating Pnu.'//
     .'%N%5f0 Band CG found by traditional method'//
     .'%N%5f1 Band CG found from true energy moment of density')

      nm='OPTIONS_SCR'; call gtv(trim(nm),tksw(prgn,nm),lscr,
     .def_i4=0,note='Use scr to accelerate convergence:'//
     .'%N%3f0 do nothing'//
     .'%N%3f1 Make ASA static response function (see documentation)'//
     .'%N%3f2 Use response to screen output q and ves'//
     .'%N%3f4 Use model response to screen output q'//
     .'%N%3f6 Use response to screen output ves only'//
     .'%N%5fAdd 1 to combine mode 1 with another mode'//
     .'%N%5fAdd 10*k to compute intra-site contribution to vbare'//
     .' each kth iteration'//
     .'%N%5fAdd 100*k to compute response function'//
     .' on every kth iteration')

C$$$C     ASA-specific
C$$$      sw = tksw(prgn,'OPTIONS_ASA')
C$$$      if (io_show+io_help/=0 .and. sw/=2)
C$$$     .  call info0(1,1,0,' * The following are ASA-specific')
C$$$      if (sw .ne. 2) then
C$$$C     This call does nothing except for printout
C$$$      nm='OPTIONS_ASA'; call gtv(trim(nm),sw,nono,Texist=ltmp,
C$$$     .  note= 'Contains the following ASA-specific tokens:')
c      nm='OPTIONS_ASA_ADNF'; call gtv(trim(nm),tksw(prgn,nm),lham4,
c     .  def_lg=F,note= 'Turn on automatic downfolding')
C$$$      nm='OPTIONS_ASA_NSPH'; call gtv(trim(nm),tksw(prgn,nm),lasa32,
C$$$     .  note='generate multipole moments in output density')
c      nm='OPTIONS_ASA_TWOC'; call gtv(trim(nm),tksw(prgn,nm),lham3,
c     .  def_i4=0,note='Two-center ASA hamiltonian.'//
c     .  '  Use TWOC=3 for 2C+pert corr')
C$$$      j = 0
C$$$      nm='OPTIONS_ASA_GAMMA'; call gtv(trim(nm),tksw(prgn,nm),j,
C$$$     .  def_i4=0,note=
C$$$     .  'gamma representation'//
C$$$     .  '%N%3f1: gamma'//
C$$$     . '%N%02f>1: (gamma(1)+gamma(nsp))/2'//
C$$$     .  '%N%3f : Some noncollinear algorithms automatically'//
C$$$     .  '%N%3f : average spins this latter representation')
C$$$      if (io_help==0) call isanrg(j,0,2,' rdccat:','GAMMA',.true.)
c      lham128 = j .gt. 0
C$$$      lasa512 = j .gt. 1
C$$$      nm='OPTIONS_ASA_CCOR'; call gtv(trim(nm),tksw(prgn,nm),lasa4,
C$$$     .  def_lg=T,note= 'Turn on combined correction')
C$$$      nm='OPTIONS_ASA_ELIN'; call gtv(trim(nm),tksw(prgn,nm),asa_elin,
C$$$     .  def_r8=0d0,note=
C$$$     .  'Energy to linearize for CCOR (2nd gen ASA 2C hamiltonian)')
c      nm='OPTIONS_ASA_NEWREP'; call gtv(trim(nm),tksw(prgn,nm),lham8,
c     .  def_lg=F,note=
c     .  'Interactively transform representation, (2nd gen ASA)')
c      nm='OPTIONS_ASA_NOHYB'; call gtv(trim(nm),tksw(prgn,nm),lham16,
c     .  def_lg=F,note='Turns off hybridisation ')
C$$$      nm='OPTIONS_ASA_MTCOR'; call gtv(trim(nm),tksw(prgn,nm),lasa64,
C$$$     .  def_lg=F,note='Turns on Ewald MT correction')
C$$$      nm='OPTIONS_ASA_QMT'; call gtv(trim(nm),tksw(prgn,nm),
C$$$     .  zbak(2),def_r8=0d0,note=
C$$$     .  'Override standard background charge for'//
C$$$     .  ' Ewald MT correction%N%3fInput only meaningful if MTCOR=T')
C$$$      endif

   69 continue                  ! end of OPTIONS

C --- Hamiltonian parameters ---
      if (io_show+io_help/=0 .and. tksw(prgn,'HAM')/=2) call info0(1,1,
     .0,' --- Parameters for hamiltonian ---')
      nsp = 1
c      if (prgn == 'MMAG') then
c        nsp = 2
c        lncol1 = T
c      endif
      nm='HAM_NSPIN'; call gtv(trim(nm),tksw(prgn,nm),nsp,
     .def_i4=1,note='Set to 2 for spin polarized calculations')
      if (io_help .eq. 0) then
        call isanrg(nsp,1,2,'rdctrl','nsp',T)
      else
        nsp = 1
      endif
ctakao
      if(prgn=='LMCHK') then
        write(6,*) '### ---> ENFORCE nsp=1 for LMCHK ###'
        nsp=1
      endif


C     Flags representing how cd is represented: default values
      if (prgn == 'LMF' .or. prgn == 'LMFGWD') lcd4=T
c      if (prgn == 'LMMC') lcd8=T

      nm='HAM_REL'; call gtv(trim(nm),tksw(prgn,nm),lrel,def_i4=1,
     .note='0 for nonrelativistic Schrodinger equation'//
     .'%N%3f1 for scalar relativistic Schrodinger equation'//
     .'%N%3f2 for Dirac equation')
C     Fully relativistic => spin-orbit coupling
      if (lrel==2) lncol4 = T

      if (nsp==2 .or. io_help/=0) then
        if (io_help /= 0) call info0(1,0,0,' * To read the magnetic '//
     .  'parameters below, HAM_NSPIN must be 2')
        nm='HAM_SO'; call gtv(trim(nm),tksw(prgn,nm),
     .  lso,def_i4=0,note=
     .  'Spin-orbit coupling (for REL=1)'//
     .  '%N%3f0 : no SO coupling'//
     .  '%N%3f1 : Add L.S to hamiltonian'//
     .  '%N%3f2 : Add Lz.Sz only to hamiltonian'//
     .  '%N%3f3 : Like 2, but also compute <L.S-LzSz> by perturbation')
        if (io_help==0) call isanrg(lso,0,3,' rdctrl:','SO',T)
        if (lso == 1) lncol4 = T
        if (lso == 2) lncol32= T
        if (lso == 3) lncol64= T
c      nm='HAM_NONCOL'; call gtv(trim(nm),tksw(prgn,nm),lncol1,
c     .  def_lg=F, note='Noncollinear magnetism')
c      nm='HAM_SS'; call gtv(trim(nm),tksw(prgn,nm),ham_qss,nmin=4,nout=
c     .  nout,note='Magnetic spin spiral, direction vector and angle')
c      lncol2 = (nout == 4)
C     If spin-orbit or SS, also turn on noncollinear
c      if (lncol2 .or. lncol4) lncol1 = T
c      if (lncol1 .or. io_help/=0) then
cC      nm='HAM_SPINTRQ'; call gtv(trim(nm),tksw(prgn,nm),lncol16,
cC     .  def_lg=F, note='Calculate spin torques (requires NONCOL=t)')
c      nm='HAM_BFIELD'; call gtv(trim(nm),tksw(prgn,nm),lncol8,
c     .  def_lg=F, note='Applied magnetic field (requires NONCOL=t)')
c      endif
      endif

      sw = tksw(prgn,'HAM_GMAX')
      if (sw/=2) then
        nm='HAM_GMAX'; call gtv(trim(nm),sw,lat_gmax,nmin=1,nout=nout,
     .  note='Energy cutoff for plane-wave mesh',or=T)
        if (nout .ne. 0) then
          sw = 2
        else
          lat_gmax = 0
        endif
        nm='HAM_FTMESH'; call gtv(trim(nm),sw,ftmesh,nout=nout,
     .  note='No. divisions for plane-wave mesh '//
     .  'along each of 3 lattice vectors.'//
     .  '%N%3fSupply one number for all vectors or a separate '//
     .  'number for each vector.')
        call fill3in(nout,ftmesh)
      endif
      nm='HAM_TOL'; call gtv(trim(nm),tksw(prgn,nm),tolft,
     .def_r8=1d-6, note='w.f. tolerance for FT mesh')
      nm='HAM_FRZWF'; call gtv(trim(nm),tksw(prgn,nm),frzwf,def_lg=F,
     .note='Set to freeze augmentation wave functions for all species')
      nm='HAM_FORCES'; call gtv(trim(nm),tksw(prgn,nm),ctrl_lfrce,
     .def_i4=0,note=
     .'Controls the ansatz for density shift in force calculation.'//
     .'%N%3f-1 no force%3f0 no shift'//
     .'%N%3f 1 free-atom shift  12 screened core+nucleus')
      if(.not.noelind()) then
        nm='HAM_ELIND'; call gtv(trim(nm),tksw(prgn,nm),elind,
     .  def_r8=0d0,note='Lindhard energy for model screening')
      else
        elind=0d0
      endif
      nm='HAM_XCFUN'; call gtv(trim(nm),tksw(prgn,nm),ham_lxcf,def_i4=2,
     .note='Specifies local exchange correlation functional:'//
     .'%N%3f1 for Ceperly-Alder (VWN)'//
     .'%N%3f2 for Barth-Hedin (ASW fit)'//
     .'%N%3f103 for PBE-GGA (use xcpbe.F in ABINIT')
      nm='HAM_GGA'; call gtv(trim(nm),tksw(prgn,nm),ggax,
     .nout=nout, def_i4=0, note='UnUsed now. See XCFUN...')
      if(ggax>0) call rx('Use HAM_XCFUN instead of HAM_GGA')

c      if (nout/=0) ham_lxcf= ham_lxcf
c      nm='HAM_EBAS'; call gtv(trim(nm),tksw(prgn,nm),
c     .  kmto,nout=nmto,note='LMTO envelope kinetic energies')
C      if (nmto>0) then
c        call packk('ham kmto', sham, kmto, n=nmto)
c        call packk('ham nmto', sham, nmto)
C      endif
      nm='HAM_ALFSI'; call gtv(trim(nm),tksw(prgn,nm),
     .alfsi, def_r8=0d0,note=
     .'Coefficient to artificial addition of overlap to hamiltonian')
      nm='HAM_DABC' ; call gtv(trim(nm),tksw(prgn,nm),
     .dabc,nout=nout,
     .note='Spacings for real-space interstital mesh')
      if (nout == 1) dabc(2) = dabc(1)
      if (nout <  3) dabc(3) = dabc(2)
!!    Merge with ?
      nm='HAM_DQVAL'; call gtv(trim(nm),tksw(prgn,nm),dqval,
     .note='Total charge')
c      call packk('ham alfsi',sham, alfsi)
c      call packk('ham dabc' ,sham, dabc,n=3)
c      call packk('ctrl zbak',sctrl,dqval)
      nm='HAM_RDSIG'; call gtv(trim(nm),tksw(prgn,nm),lrsig,def_i4=0,
     .note=
     .'Controls how self-energy is added to '//
     .'local exchange correlation functional:'//
     .'%N%3f 1s digit:'//
     .'%N%3f   0 do not read Sigma'//
     .'%N%3f   1 add sigma to potential'//
     .'%N%3f   Add 2 to symmetrize Sigma(T)'//
     .'%N%3f   Add 4 to include Re(Sigma(T)) only'//
     .'%N%3f 10s digit:'//
     .'%N%3f   0 simple interpolation'//
     .'%N%3f   1 approx high and low sigma by diagonal (see sigp)'//
     .'%N%3f   3 interpolate sigma by LDA evecs'//
     .'%N%3f     (In this mode use 100s digit for'//
     .' # interpolation points)'//
     .'%N%3f Add 10000 to indicate sigma has no symops'//
     .'%N%3f Add 20000 to use minimum neighbor table'//
     .'%N%3f Add 40000 to allow file qp mismatch')
      nm='HAM_RSRNGE'; call gtv(trim(nm),tksw(prgn,nm),rsrnge,
     .def_r8=5d0, note='Maximum range in connecting vectors for '//
     .'r.s. sigma (units of alat)')
      nm='HAM_RSSTOL'; call gtv(trim(nm),tksw(prgn,nm),
     .rsstol, def_r8=5d-6, note=
     .'Max tolerance in Bloch sum error for r.s. sigma ')
c      call packk('ham rsrnge',sham,rsrnge)
c      call packk('ham rsstol',sham,rsstol)

c      call packk('ham sigp',sham,sigp,n=8)
c      nm='HAM_NMTO'; sw = tksw(prgn,nm); call gtv(trim(nm),sw,
c     .  nmto,def_i4=0,Texist=ltmp, note=
c     .  'Order of polynomial approximation for NMTO hamiltonian')
c      if (nmto>0) call isanrg(nmto,2,5,'rdctrl','nmto',.true.)
c      if (nmto .gt. 0 .or. (io_help .ne. 0 .and. ltmp)) then
c!       Warning: kmto also used by lmmc; meaning is somewhat different
c        nm='HAM_KMTO'; call gtv(trim(nm),tksw(prgn,nm),
c     .    kmto(1:max(1,nmto)),nmin=nmto,
c     .    note='Corresponding NMTO kinetic energies'//
c     .    '%N%3fRead NMTO values, or skip if NMTO=0' )
c      endif
      nm='HAM_EWALD'; call gtv(trim(nm),tksw(prgn,nm),ham_ewald,
     .def_lg=.false.,note='Make strux by Ewald summation')
      nm='HAM_VMTZ'; call gtv(trim(nm),tksw(prgn,nm),vmtz,def_r8=0d0,
     .note='Muffin-tin zero defining wave functions')
c      i = igets('ctrl lasa,128+256',sctrl)/128
c      j = igets('ctrl lasa',sctrl) - 128*i
c      i0 = partok(recrd(catbeg),'QASA=','= ,',i,
c     .  '0 => Methfessel conventions for 2nd gen ASA moments '//
c     .             'Q0..2'//
c     .  '%N%3f1 => Q2 = coff to phidot**2 - p phi**2 in sphere'//
c     .  '%N%3f2 => Q1,Q2 accumulated as coffs to <phi*phidot> '//
c     .                 'and <phidot**2>'//
c     .  '%N%3f3 => 1+2 (Stuttgart conventions)',-1,2,0,TF(11))
c      call pack1('ctrl lasa',sctrl,j+128*mod(i,4))

C$$$      nm='HAM_QASA'; call gtv(trim(nm),tksw(prgn,nm),ham_qasa,
C$$$     .  def_i4=3,note='0 => Methfessel conventions for 2nd gen'//
C$$$     .  ' ASA moments Q0..2'//
C$$$     .  '%N%3f1 => Q2 = coff to phidot**2 - p phi**2 in sphere'//
C$$$     .  '%N%3f2 => Q1,Q2 accumulated as coffs to <phi*phidot> '//
C$$$     .                 'and <phidot**2>'//
C$$$     .  '%N%3f3 => 1+2 (Stuttgart conventions)')
      nm='HAM_PMIN'; call gtv(trim(nm),tksw(prgn,nm),pmin,
     .def_r8v=zerov,nout=nout,note=
     .'Global minimum in fractional part of P-functions.'//
     .'%N%3fEnter values for l=0..nl:'//
     .'%N%3f0: no minimum constraint'//
     .'%N%3f#: with #<1, floor of fractional P is #'//
     .'%N%3f1: use free-electron value as minimum')
C     print *, pmax
C     call snit
      nm='HAM_PMAX'; call gtv(trim(nm),tksw(prgn,nm),
     .pmax, def_r8v=zerov, nout=nout, note=
     .'Global maximum in fractional part of P-functions.'//
     .'%N%3fEnter values for l=0..nl:'//
     .'%N%3f0: no maximum constraint'//
     .'%N%3f#: with #<1, ceiling of fractional P is #')
C      print *, nout
C      print *,pmax
C      stop

      nm='HAM_OVEPS'; call gtv(trim(nm),tksw(prgn,nm),
     .oveps, def_r8=0d0, nout=nout, note=
     .'Diagonalize hamiltonian in reduced hilbert space,'//
     .'%N%3fdiscarding part with evals of overlap < OVEPS')

      nm='HAM_OVNCUT'; call gtv(trim(nm),tksw(prgn,nm),
     .ncutovl, def_i4=0, note=
     .'Diagonalize hamiltonian in reduced hilbert space,'//
     .'%N%3fdiscarding OVNCUT basis')

C ... APW basis
      nm='HAM_PWMODE'; call gtv(trim(nm),tksw(prgn,nm),pwmode,
     .def_i4=0,note=
     .'Controls APW addition to LMTO basis'//
     .'%N%3f1s digit:'//
     .'%N%6f0: LMTO basis only'//
     .'%N%6f1: Mixed LMTO+PW'//
     .'%N%6f2: PW basis only'//
     .'%N%3f10s digit:'//
     .'%N%6f0: PW basis fixed'//
     .'%N%6f1: PW basis q-dependent')
      nm='HAM_PWEMIN'; call gtv(trim(nm),tksw(prgn,nm),
     .pwemin, def_r8=0d0, nout=nout, note=
     .'Include APWs with energy E > PWEMIN (Ry)')
      nm='HAM_PWEMAX'; call gtv(trim(nm),tksw(prgn,nm),
     .pwemax, def_r8=0d0, nout=nout, note=
     .'Include APWs with energy E < PWEMAX (Ry)')
      nm='HAM_NPWPAD'; call gtv(trim(nm),tksw(prgn,nm),
     .npwpad, def_i4=-1, note=
     .'Overrides default padding of variable basis dimension')

C$$$      nm='HAM_SX'; call gtv(trim(nm),tksw(prgn,nm),lsx,
C$$$     .  def_i4=0,note='Screened exchange:'//
C$$$     .  '%N%3f 0 do nothing'//
C$$$     .  '%N%3f 1 Calculate SX Sigma'//
C$$$     .  '%N%3f11 Calculate SX Sigma w/ local W')
C$$$      nm='HAM_SXOPTS'; call gtv(trim(nm),tksw(prgn,nm),sxopt,nmin=10,
C$$$     .  note='Options for screened exchange,.e.g.  SXOPTS=rssig;nit=3')

      sw = tksw(prgn,'HAM_SIGP')
C     This call does nothing except for printout
      if (sw .ne. 2) then
        nm='HAM_SIGP'; call gtv(trim(nm),sw,nono,Texist=ltmp,
     .  note= 'Parameters for replacing high (and possibly low) '//
     .  'subblocks of sigma'//
     .  '%N%3fwith a linear function of LDA eigenvalues (diagonal)')
        nm='HAM_SIGP_MODE'; call gtv(trim(nm),sw,sigp_mode,def_i4=3,
     .  nout=nout,
     .  note= 'Specifies linear fitting function'//
     .  '%N%3f 0 sigii > a + b*e'//
     .  '%N%3f 1 sigii = a + b*e'//
     .  '%N%3f 2 sigii > a  and  sigii < b'//
     .  '%N%3f 3 sigii = a + b*e')
        nm='HAM_SIGP_A'; call gtv(trim(nm),sw,sigp_a,def_r8=0.02d0,
     .  note='a parameter in fit')
        nm='HAM_SIGP_B'; call gtv(trim(nm),sw,sigp_b,def_r8=0.06d0,
     .  note='b parameter in fit')
        nm='HAM_SIGP_NMIN'; call gtv(trim(nm),sw,sigp_nmin,def_i4=0,
     .  note='Replace subblock 1..nmin with fit')
        nm='HAM_SIGP_EMIN'; call gtv(trim(nm),sw,sigp_emin,def_r8=0d0,
     .  note='Replace subblock elda<emin with fit (not used if nmin>=0)')
        nm='HAM_SIGP_NMAX'; call gtv(trim(nm),sw,sigp_nmax,def_i4=0,
     .  note='Replace subblock nmax+1.. with fit')
        nm='HAM_SIGP_EMAX'; call gtv(trim(nm),sw,sigp_emax,def_r8=2d0,
     .  note='Replace subblock elda>emax with fit (not used if nmax>0)')
        nm='HAM_SIGP_EFIT'; call gtv(trim(nm),sw,sigp_efit,def_r8=0d0,
     .  note='Energy cutoff for fitting sigii '//
     .  '(not used in any calculation')
      endif

      nm='HAM_UDIAG'; call gtv(trim(nm),tksw(prgn,nm),ham_udiag,
     .def_i4=0,note='nonzero => diagonal-only LDA+U')

C --- Symmetry group ---
      if (io_show+io_help/=0 .and. tksw(prgn,'SYMGRP')/=2)
     .call info0(1,1,0,' --- Symmetry group operations ---')
      nm='SYMGRP'; call gtv(trim(nm),tksw(prgn,nm),symg,
     .note='Generators for symmetry group')

C$$$C --- Crystal Green's function ---
C$$$      if (tksw(prgn,'GF') == 2) goto 269
C$$$      if (io_show+io_help/=0) call info0(1,1,
C$$$     .  0,' --- Parameters for GF ---')
C$$$      nm='GF_MODE'; call gtv(trim(nm),tksw(prgn,nm),ctrl_lcgf,
C$$$     .  def_i4=0,note=' 0: do nothing'//
C$$$     . '%N%3f 1: self-consistent cycle'//
C$$$     .'%N%3f10: Transverse exchange interactions J(q), MST'//
C$$$     .'%N%3f11: Read J(q) from disk and print derivative properties'//
C$$$C    .'%N%3f12: generate G(R,R'')'//
C$$$     .'%N%3f14: Longitudinal exchange interactions J(q), MST'//
C$$$     .'%N%3f20: Transverse chi+- from ASA GF'//
C$$$     .'%N%3f21: Read chi from disk and print derivative properties'//
C$$$     .'%N%3f24: Transverse chi++,chi-- from ASA GF')

C$$$      nm='GF_GFOPTS'; call gtv(trim(nm),tksw(prgn,nm),gfopt,nmin=10,
C$$$     .  note='Switches governing execution of GF')
C$$$  269 continue

C$$$C --- Planar Green's function ---
C$$$      if (tksw(prgn,'PGF') == 2) goto 219
C$$$      if (io_show+io_help/=0) call info0(1,1,
C$$$     .  0,' --- Parameters for PGF ---')
C$$$C     nm='PGF'; call gtv(trim(nm),tksw(prgn,nm),Texist=ltmp,note=
C$$$C    .  'Parameters for planar Green''s function')
C$$$      ctrl_lpgf(1) = 0
C$$$      nm='PGF_MODE'; call gtv(trim(nm),tksw(prgn,nm),
C$$$     .  ctrl_lpgf(1),note='0: do nothing'//
C$$$     .'%N%3f1: diagonal layer GF'//
C$$$     .'%N%3f2: left- and right-bulk bulk GF'//
C$$$     .'%N%3f3: find k(E) for left bulk'//
C$$$     .'%N%3f4: find k(E) for right bulk'//
C$$$     .'%N%3f5: Calculate current')
C$$$      if (io_help == 0) then
C$$$        if (ctrl_lpgf(1) == NULLI) ctrl_lpgf(1) = 0
C$$$        if (ctrl_lpgf(1)==0) goto 219
C$$$      else
C$$$      if (io_show+io_help/=0) call info0(1,0,0,
C$$$     .    ' ... The following are read if MODE>0')
C$$$      endif
C$$$      nm='PGF_SPARSE'; call gtv(trim(nm),tksw(prgn,nm),ctrl_lpgf(2),
C$$$     .  def_i4=0,note='0: Calculate GF layer by layer'//
C$$$     .  '%N%3f1: Calculate GF with LU decomposition')
C$$$      nm='PGF_PLATL'; call gtv(trim(nm),tksw(prgn,nm),
C$$$     .  platl, nmin=3,note='third lattice vector of left bulk')
C$$$      nm='PGF_PLATR'; call gtv(trim(nm),tksw(prgn,nm),
C$$$     .  platr, nmin=3,note='third lattice vector of right bulk')
C$$$      gfopt = ' '
C$$$      nm='PGF_GFOPTS'; call gtv(trim(nm),tksw(prgn,nm),gfopt,nmin=10,
C$$$     .  note='Switches governing execution of PGF')
C$$$  219 continue

C --- Species (old CLASS) ---
      if (io_help == 1) nspec = 1
      if (tksw(prgn,'SPEC') == 2) goto 79
      if (io_show+io_help/=0) call info0(1,1,
     .0,' --- Parameters for species data ---')
C     SPEC call is only helpful for printout
C      if (io_help /= 0 .or. io_show>0) write(stdo,'(1x)')
C      nm='SPEC'; call gtv(trim(nm),tksw(prgn,nm),Texist=ltmp,note=
C     .  'Parameters for species data')
C      if (.not. ltmp) goto 79
      if (io_help /= 0) call info0(1,0,0,' * The next four tokens '//
     .'apply to the automatic sphere resizer')
      nm='SPEC_SCLWSR'; call gtv(trim(nm),tksw(prgn,nm),
     .sclwsr, def_r8=0d0, note=
     .'Scales sphere radii, trying to reach volume = '//
     .'SCLWSR * cell volume'//
     .'%N%3fSCLWSR=0 turns off this option.'//
     .'%N%3fAdd  10  to initially scale non-ES first;'//
     .'%N%3f or  20  to scale ES independently.')
C     if (sclwsr .ne. 0) then
      nm='SPEC_OMAX1'; call gtv(trim(nm),tksw(prgn,nm),omax1,
ctakao
c     .  def_r8v=(/.16d0,.18d0,.20d0/),note=
     .def_r8v=(/0d0,0d0,0d0/),note=
     .'Limits max sphere overlaps when adjusting MT radii')
      nm='SPEC_OMAX2'; call gtv(trim(nm),tksw(prgn,nm),
ctakao
c     .  omax2, def_r8v=(/.40d0,.45d0,.50d0/),note=
     .omax2, def_r8v=(/0d0,0d0,0d0/),note=
     .'Sphere overlap constraints of second type',nout=nout)
      nm='SPEC_WSRMAX'; call gtv(trim(nm),tksw(prgn,nm),wsrmax,
     .def_r8=0d0,note=
     .'If WSRMAX is nonzero, no sphere radius may exceed its value')
C     endif

      if (io_help .eq. 1) then
        write(*,382)
  382   format(/' * ',
     .  'The following tokens are input for each species. ',
     .  'Data sandwiched'/3x,'between successive occurences of ',
     .  'token ATOM apply to one species.')
        nspec = 1
      endif

      if (nspec .eq. 0) goto 79

      allocate(pnu(n0,nsp,nspec),qnu(n0,nsp,nspec),
     .pz(n0,nsp,nspec),amom(n0,nspec),idmod(n0,nspec),
     .rsmh(n0,nspec),eh(n0,nspec),rsmh2(n0,nspec),eh2(n0,nspec),
     .pb1(nspec),pb2(nspec),lmxpb(nspec),
     .ehvl(n0,nspec),
     .qpol(n0,nspec),stni(nspec),tbvso(4,nspec),
     .iq1(n0,nspec),ivso(n0,nspec),
     .rg(nspec),rsma(nspec),rfoca(nspec),rsmfa(nspec),rcfa(2,nspec),
     .rs3(nspec),rham(nspec),rmt(nspec),rsmv(nspec),
     .nxi(nspec),exi(n0,nspec),rint(nspec),rcut(nspec),
     .spec_a(nspec),z(nspec),nr(nspec),mass(nspec),eref(nspec),
     .coreh(nspec),coreq(2,nspec),
     .colxbs(3,nspec),radxbs(nspec),
     .idxdn(n0,nspec),
     .hcr(n0,nspec),alpha(n0,nspec),
     .idu(4,nspec),uh(4,nspec),jh(4,nspec),
     .dv(nspec),grp(nspec),grp2(nspec),
     .mxcst1(nspec),mxcst2(nspec),mxcst4(nspec),
     .kmxt(nspec),kmxv(nspec),
     .lfoca(nspec),lmxl(nspec),lxi(nspec),lmxa(nspec),lmxb(nspec))

ctakao Sep2009
      mxcst1=F
      mxcst2=F
      mxcst4=F

c      jc = 0
c      if (optio .eq. 2) nkaph = nglob('nkaph')
C     Loop over species
c
ctakao test
c      if(prgn=='LMF') then
c      do  j = 1, nspec
c        jj= (/1,j/)
c        nm='SPEC_ATOM_RSMH'; call gtv(trim(nm),0,
c     .    rsmh(1,j),cindx=jj,nout=nout,
c     .    note='Smoothing radii for basis')
c        print *,'qqq: nout=',nout
c      enddo
c      endif
c      stop 'xxxxxxxxxxxxxxxxxxxxxxxx'
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
cccccccccccccccccccccctakao
        if(debug) print *,'nspec mxcst j-loop j nspec',j,nspec
ccccccccccccccccccccccc
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
        nm='SPEC_ATOM'; call gtv(trim(nm),tksw(prgn,nm),slabl(j),
     .  nmin=10,cindx=jj,note='Species label')

        nm='SPEC_ATOM_Z'; call gtv(trim(nm),tksw(prgn,nm),z(j),
     .  cindx=jj,note='Atomic number')

cccccccccccccccccccccctakao
        if(debug) print *,'nspec aaa mxcst j-loop j nspec',j,nspec
        sw = tksw(prgn,'SPEC_ATOM_R')
        if (sw /= 2) then
          nout = 0
          nm='SPEC_ATOM_R'; call gtv(trim(nm),sw,rmt(j),cindx=jj,
     .    nout=nout,note= 'Augmentation sphere radius rmax',or=T)
          if (nout .ne. 1) then
            nm='SPEC_ATOM_R/W';call gtv(trim(nm),sw,rmt(j),cindx=jj,
     .      nout=nout,note='rmax relative to average WS radius',or=T)
            if (nout == 1) then
              rmt(j) =rmt(j)*avw
            else
              nm='SPEC_ATOM_R/A';call gtv(trim(nm),sw,rmt(j),cindx=jj,
     .        nout=nout,note='rmax relative to lattice constant')
              if (nout == 1) then
                rmt(j) =rmt(j)*alat
              else
c takao for lmchk even when R is not given. See default of LMCHK
                rmt(j) = 0d0
              endif
            endif
          endif
        endif
c        if (ltbe) rmt(j) = 1d0

cccccccccccccccccccccctakao
        if(debug) print *,'nspec bbb mxcst j-loop j nspec',j,nspec

C   ... Radial mesh parameters: determine default value of a
        i0 = NULLI
        xx = NULLR
        if (io_help .ne. 1) then
          call pshpr(0)
          call rmesh(z(j),rmt(j),lrel,.false.,nrmx,xx,i0)
          call poppr
          if (xx .eq. .03d0) xx = .025d0
        endif
        nm='SPEC_ATOM_A'; call gtv(trim(nm),tksw(prgn,nm),spec_a(j),
     .  def_r8=xx,cindx=jj,nout=nout,
     .  note='Radial mesh point spacing parameter')
C       Determine default NR
        if (tksw(prgn,'SPEC_ATOM_NR') .ne. 2) then
          i0 = 0
          call pshpr(0)
          call rmesh(z(j),rmt(j),lrel,.false.,nrmx,spec_a(j),i0)
          call poppr
        endif
        nm='SPEC_ATOM_NR'; call gtv(trim(nm),tksw(prgn,nm),nr(j),
     .  def_i4=i0,cindx=jj, note='Number of radial mesh points')
        if (nr(j) .eq. 0) nr(j) = i0
        nm='SPEC_ATOM_LMX'; call gtv(trim(nm),tksw(prgn,nm),lmxb(j),
     .  def_i4=max(nl-1,NULLI),cindx=jj,note='l-cutoff for basis')
C       Running account of maximum lmxb
        lmxbj = NULLI
        if (io_help .ne. 1) then
          if (prgn .eq. 'LMMC') then
            call lx2vec(lmxb(j),0,nn,ivec)
            do  i = 1, nn
              lmxbx = max(lmxbx,ivec(i))
              lmxbj = max(lmxbj,ivec(i))
            enddo
          else
            lmxbx = max(lmxbx,lmxb(j))
            lmxbj = lmxb(j)
          endif
        endif
C       nlbj = number of elements associated with lmxb
C              0 => no elements
C       nlbji: ditto, but used to specify number of default values
        nlbj = 1+lmxbj
        nlbji = nlbj
        if (io_help .eq. 1) then
          nlbji = NULLI
          nlbj = 1
        elseif (lmxbj == NULLI) then
          nlbji = NULLI
          nlbj = 0
        endif
cccccccccccccccccccccctakao
        if(debug) print *,'nspec ccc mxcst j-loop j nspec',j,nspec

C   ... Basis set for lmf
        nm='SPEC_ATOM_RSMH'; call gtv(trim(nm),tksw(prgn,nm),
     .  rsmh(1:nlbj,j),cindx=jj,nout=nout,nmin=nlbj,def_r8v=zerov,
     .  note='Smoothing radii for basis')
        nn = NULLI
        if (nout .gt. 0) nn = nout
C       Find nn=lmax+1 for which rsmh ne 0 => reduce # EH required
        do  i = 1, nlbj
          if (rsmh(i,j) .eq. 0) cycle
          nn = i
        enddo
        nm='SPEC_ATOM_EH'; call gtv(trim(nm),tksw(prgn,nm),
     .  eh(1:nlbj,j),nmin=nn,cindx=jj,
     .  note='Kinetic energies for basis')
        nn = NULLI
        if (nout>0) nn = nout
        nm='SPEC_ATOM_RSMH2'; call gtv(trim(nm),tksw(prgn,nm),
     .  rsmh2(1:nlbj,j),nmin=nlbj,def_r8v=zerov,cindx=jj,nout=nout,
     .  Texist=ltmp,note='Basis smoothing radii, second group')
        if (ltmp) then
          nn = NULLI
          sw = tksw(prgn,nm)
          if (nout>0) then
            nkaph=2
            nn = nout
            sw = 1
C         Find nn=lmax+1 for which rsmh ne 0 => reduce # EH required
            do  i = 1, nlbj
              if (rsmh2(i,j) .eq. 0) cycle
              nn = i
            enddo
          endif
          nm='SPEC_ATOM_EH2'; call gtv(trim(nm),sw,
     .    eh2(1:nlbj,j),nmin=nn,cindx=jj,
     .    note='Basis kinetic energies, second group')
        endif
cccccccccccccccccccccctakao
        if(debug) print *,'nspec ddd mxcst j-loop j nspec',j,nspec

        nm='SPEC_ATOM_EHVL'; call gtv(trim(nm),tksw(prgn,nm),
     .  ehvl(1:nlbj,j),cindx=jj,def_r8v=(/(-0.5d0,i=1,n0)/),
     .  nout=nout,nmin=nlbj, note='val-lap fit energies')

C   ... Determine lmxa: floating orbitals have no lmxa
        nm='SPEC_ATOM_LMXA'; sw = tksw(prgn,nm)
ccccccccccccccccccc takao ccccc
        if(debug) print *,'nspec ddd000 mxcst j rmt io_help lmxa=',j,rmt(j),io_help,lmxa(j)
cccccccccccccccccccccccccccccc
        if (rmt(j) == 0.and.io_help/=1) then !takao iohelp/=1 added.
          lmxa(j) = -1
        elseif (sw == 2) then   !lmxa not read: look for subsitute
          if (lfp > 0) then     !Possibly replace with a formula
            !that depends on rmt
            lmxa(j) = 4
          elseif (tksw(prgn,'SPEC_ATOM_LMX') /= 2) then
            lmxa(j) = lmxb(j)
          elseif (tksw(prgn,'STRUC_NL') /= 2) then
            lmxa(j) = nl-1
          else
            lmxa(j) = 0
          endif
        else
C         Possibly replace with a formula that depends on rmt
          call gtv(trim(nm),sw,lmxa(j),
     .    def_i4=max(nl-1,lmxbj,NULLI),cindx=jj,Texist=ltmp,note=
     .    'l-cutoff for augmentation')
C         lmxb may not exceed lmxa
          if (io_help == 0 .and. lmxbj .gt. lmxa(j))
     .    call rx2('species '//trim(slabl(j))//' : LMX=%i '//
     .    'exceeds LMXA=%i.  Revise input so that LMX<=LMXA.',
     .    lmxbj,lmxa(j))
        endif
        lmxaj = lmxa(j)
C       nlaj = number of elements associated with lmxa
C              0 => no elements
C       nlaji: ditto, but used to specify number of default values
        nlaj = 1+lmxaj
ccccccccccccccccccc takao ccccc
        if(debug) print *,'nspec ddd111 mxcst j-loop j nspec nlaj lmxaj=',j,nspec,nlaj,lmxaj
cccccccccccccccccccccccccccccc
        nlaji = nlaj
        if (io_help .eq. 1) then
          nlaji = NULLI
          nlaj = 1
        elseif (lmxaj == NULLI) then
          nlaji = NULLI
          nlaj = 0
        endif
C        if (lmxaj == NULLI) nlaj = NULLI
C        if (lmxaj == NULLI) nlaji = NULLI



C   ... Parameters that depend on the existence of an augmentation sphere
C       lmxl = l-cutoff for numerical rep'sn of density in sphere
C       lfoca = mode for treating core
C       kmxt = kmax for expansion of envelope wf tails
C       kmxv = cutoff to expand smoothed potential
C       hardwired for now
        kmxv(j) = 15            !Not input
C       Cannot set default here: need set after rescaling of rmt
C       rsmv(j) = rmt(j)*.5d0   !Not input
        rsmv(j) = 0d0           !Not input
        kmxt(j) = -1            !If sought, default will be set below
        lfoca(j) = 0            !If sought, default will be reset below
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

cccccccccccccccccccccctakao
        if(debug) print *,'nspec ddd222 mxcst j-loop j nspec nlaj=',j,nspec,nlaj
ccccccccccccccccccccccccccc

        if (nlaj .ne. 0) then

          nm='SPEC_ATOM_LMXL'; call gtv(trim(nm),tksw(prgn,nm),lmxl(j),
     .    cindx=jj,def_i4=lmxaj,note=
     .    'lmax for which to accumulate rho,V in sphere')

C   ... Set up default P,Q in absence of explicit specification
          call dpzero(pnu(1,1,j),nsp*n0)
          call dpzero(qnu(1,1,j),nsp*n0)
C       This branch sets defaults beforehand ... too complicated
C        if (io_help == 0) then
CC         Default P,Q, nsp=1
C          if (debug) call pshpr(100)
C          call defpq(z(j),lmxaj,1,pnu(1,1,j),qnu(1,1,j))
C          if (debug) call poppr
C        endif
C        nm='SPEC_ATOM_P'; call gtv(trim(nm),tksw(prgn,nm),
C     .    pnu(1:nlaj,1,j),def_r8v=pnu(:,1,j),nmin=nlaji,cindx=jj,note=
C     .    'Starting log der. parameters for each l')
C        nm='SPEC_ATOM_Q'; call gtv(trim(nm),tksw(prgn,nm),
C     .    qnu(1:nlaj,1,j),def_r8v=qnu(:,1,j),nmin=nlaji,cindx=jj,note=
C     .    'Starting sphere charges for each l channel')
C        call snit
          nm='SPEC_ATOM_P'; call gtv(trim(nm),tksw(prgn,nm),
     .    pnu(1:nlaj,1,j),cindx=jj,note=
     .    'Starting log der. parameters for each l')
          nm='SPEC_ATOM_Q'; call gtv(trim(nm),tksw(prgn,nm),
     .    qnu(1:nlaj,1,j),cindx=jj,note=
     .    'Starting valence charges for each l channel.'//
     .    '%N%2f Q do not include semicore(PZ) electrons.'//
     &    '%N%2f Charge configuration is shown by lmfa %N'//
     &    '%N%2f WARN: This version cannot treat two valence channels'//
     &    '%N%2f per l (Q for a l-channl is zero if the l is with PZ).'//
     &    '%N%2f This causes a problem typically in Li; then we '//
     &    '%N%2f can not treat both of PZ=1.9 and P=2.2 as valence.'//
     &    '%N%2f To avoid this, use Q=0,1 together.'//
     &    ' This trick supply an '//
     &    '%N%2f electron to 2p channel; this trick works fine.')

C       Reset default P,Q in absence of explicit specification
          if (io_help == 0) then
            if (io_show .ne. 0) call pshpr(50)
C         call snit
            call defpq(z(j),lmxaj,1,pnu(1,1,j),qnu(1,1,j))
C         call snit
            if (io_show .ne. 0) call poppr
          endif
          if (nsp .eq. 2) call dcopy(n0,pnu(1,1,j),1,pnu(1,2,j),1)
          if (nsp .eq. 2 .or. io_help == 1) then
            nm='SPEC_ATOM_MMOM'; call gtv(trim(nm),tksw(prgn,nm),
     .      qnu(1:nlaj,2,j),def_r8v=zerov,nmin=nlaji,cindx=jj,note=
     .      'Starting mag. moms for each l channel.'//
     .      '%N%2f For a chanel with PZ, this is enforced to be zero.'//
     .      '%N%2f See explanation for SPEC_ATOM_Q.')
          endif

          nm='SPEC_ATOM_PZ'; call gtv(trim(nm),tksw(prgn,nm),
     .    pz(1:nlaj,1,j),def_r8v=zerov,nmin=nlaji,cindx=jj,note=
     .    'Starting semicore log der. parameters'//
     .    '%N%10fAdd 10 to attach Hankel tail'//
     .    '%N%10fAdd 20 to include perturbatively',nout=nout)
          if (nout>0) then
            if (dasum(nlaj,pz(1,1,j),1) /= 0) lpzi = max(lpzi,2)
          endif

          i0 = 1
          if (z(j) .le. 8) i0 = 0
          if (io_help .eq. 1) i0 = NULLI
          nm='SPEC_ATOM_LFOCA';call gtv(trim(nm),tksw(prgn,nm),lfoca(j),
c     .    def_i4=i0,cindx=jj,note='FOCA switch 0, 1 or 2 (see docs)')
     .    def_i4=i0,cindx=jj,note='FOCA switch 0(within MT):'
     .    //'=1(frozenCore). Default: 1 for z>8;0 for z<=8')     !takao Aug2010
          if(io_help==0) then
          if(lfoca(j)/=0.and.lfoca(j)/=1) then
            call rx('LFOCA should be 0 or 1 (Aug2010): 2 is not allowed')
          endif
          endif

C       call snot
          nm='SPEC_ATOM_KMXA'; call gtv(trim(nm),tksw(prgn,nm),
     .    kmxt(j),def_i4=3,cindx=jj,note=
     .    'k-cutoff for projection of wave functions in sphere.')
C       Cannot set default here: need set after rescaling of rmt
C       xx = 0.40d0*rmt(j); if (io_help .eq. 1) xx = NULLI
          xx = 0d0; if (io_help .eq. 1) xx = NULLI
          nm='SPEC_ATOM_RSMA'; call gtv(trim(nm),tksw(prgn,nm),rsma(j),
     .    def_r8=xx,cindx=jj,note=
     .    'Smoothing for projection of wave functions in sphere.'//
     .    '%N%3finput<0 => choose default * -input')

C       Cannot set default here: need set after rescaling of rmt
C       xx = 0.25d0*rmt(j); if (io_help .eq. 1) xx = NULLI
          xx = 0d0; if (io_help .eq. 1) xx = NULLI
          nm='SPEC_ATOM_RSMG'; call gtv(trim(nm),tksw(prgn,nm),rg(j),
     .    def_r8=xx,cindx=jj,note=
     .    'Smoothing for projection of charge in sphere.'//
     .    '%N%3finput<0 => choose default * -input')
C       Cannot set default here: need set after rescaling of rmt
C       xx = 0.40d0*rmt(j); if (io_help .eq. 1) xx = NULLI
          xx = 0d0; if (io_help .eq. 1) xx = NULLI
          nm='SPEC_ATOM_RFOCA'; call gtv(trim(nm),tksw(prgn,nm),rfoca(j),
     .    def_r8=xx,cindx=jj,note=
     .    'Smoothing for core tail.  input<0 => choose default * -input')
C       Cannot set default here: need set after rescaling of rmt
C       xx = 0.50d0*rmt(j); if (io_help .eq. 1) xx = NULLI
          xx = 0d0; if (io_help .eq. 1) xx = NULLI
          nm='SPEC_ATOM_RSMFA'; call gtv(trim(nm),tksw(prgn,nm),rsmfa(j),
     .    def_r8=xx,cindx=jj,note=
     .    'Smoothing for free atom.  input<0 => choose default * -input')
          nm='SPEC_ATOM_RCFA'; call gtv(trim(nm),tksw(prgn,nm),rcfa(1:2,j)
     .    ,def_r8v=zerov,nmin=2,cindx=jj,note=
     .    'Cutoff radius for renormalization of free atom density.'//
     .    '%N%3fOptional 2nd argument = width'//
     .    '%N%3fRCFA<0 => renormalize potential instead of density')

C       Negative radii: convert to actual numbers
          if (rg(j)    .lt. 0) rg(j)    = -rg(j)*0.25d0*rmt(j)
          if (rsma(j)  .lt. 0) rsma(j)  = -rsma(j)*0.4d0*rmt(j)
          if (rfoca(j) .lt. 0) rfoca(j) = -rfoca(j)*0.4d0*rmt(j)
          if (rsmfa(j) .lt. 0) rsmfa(j) = -rsmfa(j)*0.5d0*rmt(j)

          nm='SPEC_ATOM_RS3'; call gtv(trim(nm),tksw(prgn,nm),rs3(j),
     .    def_r8=1d0,cindx=jj,
     .    note='Minimum smoothing radius for local orbital')

C       Downfolding switches: if auto DNF turned on, default is zero.
c        if (lham4) then
c          nn = 0
cC       Help mode: If ADNF COULD have been turned on, default unknown
c        elseif (tksw(prgn,'OPTIONS_ADNF')/= 2.and. io_help==1) then
c          nn = NULLI
cC       All other cases: default is 1
c        else
          nn = 1
c        endif
          nm='SPEC_ATOM_IDXDN'; call gtv(trim(nm),tksw(prgn,nm),
     .    idxdn(1:nlbj,j),nmin=nlbj,
     .    def_i4v=(/(nn,i=1,nlbj)/),cindx=jj,
     .    note='downfolding index: 0, auto; 1, no dnf; 2, fold down;'//
     .    ' 3, neglect')

C       'Hard Core' radius entering into strux
          nm='SPEC_ATOM_HCR'; call gtv(trim(nm),tksw(prgn,nm),
     .    hcr(1:nlbj,j),nmin=nlbj,cindx=jj,nout=nout,or=T,
     .    note='Hard sphere radii for structure constants')
          if (nout == 0) then !nout=-1 if sw=2; otherwise nout=0 unless data was read
            nm='SPEC_ATOM_HCR/R'; call gtv(trim(nm),tksw(prgn,nm),
     .      hcr(1:nlbj,j),nmin=nlbj,cindx=jj,nout=nout,
     .      def_r8v=(/(0.7d0,I=1,n0)/),note=
     .      'Hard sphere radii for structure constants, units of R')
            call dscal(nout,rmt(j),hcr(1,j),1)
          endif

cccccccccccccccccccccctakao
          if(debug) print *,'nspec eee mxcst j-loop j nspec',j,nspec

C       ASA tight-binding alpha parameters
          call dpzero(xv,n0)
          lp1 = nlbj
          nm='SPEC_ATOM_ALPHA'; sw=tksw(prgn,nm)
C       Set default values; poke into xv
          if (sw /= 2) then
C       Help mode, default -> unknown
            if (io_help == 1) then
              xv(1) = NULLI
            else
C         Determine effective 1+lmxb: neglected orbitals => reduce lmxb
              if (tksw(prgn,'SPEC_ATOM_IDXDN') /= 2) then
   75           if (lp1 .gt. 1 .and. idxdn(lp1,j) .ge. 4) then
                  lp1 = lp1-1
                  goto 75
                endif
              endif
              call dcopy(5,tbalph(1,lp1),1,xv,1)
            endif
          endif

          i0 = lp1
          if (lmxbj == NULLI) i0 = NULLI
C       Copy all 5 screening parameters, to avoid division zero
          call dcopy(5,tbalph(1,lp1),1,alpha(1,j),1)
          call gtv(trim(nm),sw,alpha(1:lp1,j),cindx=jj,nmin=i0,
     .    def_r8v=xv,note='Screening parameters for structure constants')

          nm='SPEC_ATOM_IDMOD'; call gtv(trim(nm),tksw(prgn,nm),
     .    idmod(1:nlaj,j),nmin=nlaji,def_i4v=(/(0,i=1,n0)/),
     .    cindx=jj,note=
     .    'idmod=0 floats P to band CG, 1 freezes P, 2 freezes enu')

          nm='SPEC_ATOM_DV'; call gtv(trim(nm),tksw(prgn,nm),dv(j),
     .    def_r8=0d0,cindx=jj,note='Artificial constant potential '//
     .    'shift added to spheres belonging to this species')
ccccccccccccccc
c        mxcst1(j)=F
ccccccccccccccc
          nm='SPEC_ATOM_MIX'; call gtv(trim(nm),tksw(prgn,nm),mxcst1(j),
     .    def_lg=F,cindx=jj,note='Set to suppress '//
     .    'self-consistency of classes in this spec')
ccccccccccccccccccccccccc
          if(debug) print *,' mxcst nspec 111 =',nspec,j,jj,mxcst1(j)
c        stop 'xxxxxxxxxxxx'
ccccccccccccccccccccccccc
          nm='SPEC_ATOM_CSTRMX'; call gtv(trim(nm),tksw(prgn,nm),
     .    mxcst2(j),cindx=jj,def_lg=F,note='Set to exclude this'//
     .    ' species when automatically resizing sphere radii'//
     .    ' (SCLWSR>0)')
          if (sclwsr .eq. 0) mxcst2(j) = F

          nm='SPEC_ATOM_FRZWF'; call gtv(trim(nm),tksw(prgn,nm),
     .    mxcst4(j),cindx=jj,def_lg=F,note=
     .    'Set to freeze augmentation wave functions for this species')

        endif    ! end of input dependent on presence of aug sphere.
cccccccccccccccccccccctakao
        if(debug) print *,'nspec eee222 mxcst j-loop j nspec',j,nspec

        nm='SPEC_ATOM_RHAM'; call gtv(trim(nm),tksw(prgn,nm),rham(j),
     .  cindx=jj,note='Radius delimiting range of r.s. hamiltonian')

        nm='SPEC_ATOM_GROUP'; call gtv(trim(nm),tksw(prgn,nm),
     .  grp(j),def_i4=0,cindx=jj)
        nm='SPEC_ATOM_GRP2'; call gtv(trim(nm),tksw(prgn,nm),
     .  grp2(j),def_i4=0,cindx=jj)

        sw = tksw(prgn,'SPEC_ATOM_IDU')
        if (io_help .gt. 0 .and. sw .lt. 2) then
          call info0(1,0,0,' * ... The next three tokens are for LDA+U')
        endif
        nm='SPEC_ATOM_IDU'; call gtv(trim(nm),sw,idu(:,j),cindx=jj,
     .  def_i4v=(/(0,i=1,n0)/),note=
     .  'LDA+U mode:  0 nothing, 1 AMF, 2 FLL, 3 mixed ')
        nm='SPEC_ATOM_UH'; call gtv(trim(nm),sw,uh(:,j),cindx=jj,
     .  def_r8v=zerov,note='Hubbard U for LDA+U')
        nm='SPEC_ATOM_JH'; call gtv(trim(nm),sw,jh(:,j),cindx=jj,
     .  def_r8v=zerov,note='Exchange parameter J for LDA+U')

        nm='SPEC_ATOM_PBAS'; call gtv(trim(nm),tksw(prgn,nm),pb1(j),
     .  nmin=10,cindx=jj,note='product basis for GW')
        nm='SPEC_ATOM_PBAS2'; call gtv(trim(nm),tksw(prgn,nm),
     .  pb2(j),cindx=jj,note='second product basis for GW')
        nm='SPEC_ATOM_LMXPB';call gtv(trim(nm),tksw(prgn,nm),lmxpb(j),
     .  def_i4=4,cindx=jj,note='l-cutoff for product basis')
c        if (nout==0) spec_lmxpb(j)='        '

C   ... tbe specific
        nm='SPEC_ATOM_QPOL'; call gtv(trim(nm),tksw(prgn,nm),
     .  qpol(:,j),def_r8v=zerov,cindx=jj,note=
     .  'up to ten polarisability parameters')
        nm='SPEC_ATOM_I'; call gtv(trim(nm),tksw(prgn,nm),stni(j),
     .  def_r8=0d0,cindx=jj,note='TB Stoner I (l=2)')
        nm='SPEC_ATOM_FRZQ1'; call gtv(trim(nm),tksw(prgn,nm),
     .  iq1(1:nlbj*nsp,j),nmin=nlbji*nsp,def_i4v=izerv,
     .  note='Switches to freeze site energies')
        nm='SPEC_ATOM_FRZVSO'; call gtv(trim(nm),tksw(prgn,nm),
     .  ivso(2:nlbj,j),nmin=max(nlbji-1,NULLI),def_i4v=izerv,
     .  note='Switches to freeze TB spin-orbit parameters')

        if (lncol4 .and. lmxaj .ge. 1 .or. io_help .ne. 0) then
c          call dpzero(tmp,n0)
c          if (optio .eq. 1) call pack2('spec alpha',sspec,j,tmp)
c          i0 = partks(recrd(catbeg),'VSO=','=','spec alpha:2',sspec,j,
c     .      'TB spin-orbit parameters',-lmxb,4,i1,TF(30))
c          tmp=0
          tbvso(:,j)=0
          nm='SPEC_ATOM_VSO'
          call gtv(trim(nm),tksw(prgn,nm),tbvso(2:nlbj,j),
     .    def_r8v=zerov,cindx=jj,note='TB spin-orbit parameters')
        endif

C   ... molecules-specific.  lmxb will change, so use TF(34) for now
c        i0 = partks(recrd(catbeg),'LMXB=','=','spec lmxpb',sspec,j,
c     .    'compound l-cutoffs for basis (molecules)',-1,2,i1,TF(34))
c        i0 = partks(recrd(catbeg),'LXI=','=','spec lxi',sspec,j,
c     .    'l-cutoffs for interstitial product basis',-1,2,i1,TF(34))
c        i = igetss('spec lxi',j,sspec)
c        call lx2vec(i,0,ii,ivec)
c        i0 = partks(recrd(catbeg),'EXI=','=','spec exi',sspec,j,
c     .    'Hankel energies for interstitial product basis',-ii,4,i1,
c     .    TF(34))
c        i0 = partks(recrd(catbeg),'RINT=','=','spec rint',sspec,j,
c     .    'range of interstitial product basis',-1,4,i1,TF(35))
c        i0 = partks(recrd(catbeg),'RCUT=','=','spec rcut',sspec,j,
c     .    'radius over which true,smooth Hankels differ',-1,4,i1,TF(36))

        nm='SPEC_ATOM_LMXB'; call gtv(trim(nm),tksw(prgn,nm),lmxpb(j),
     .  cindx=jj,note='compound l-cutoffs for basis (molecules)')
        nm='SPEC_ATOM_LXI'; call gtv(trim(nm),tksw(prgn,nm),lxi(j),
     .  cindx=jj,note='l-cutoffs for interstitial product basis')
        nm='SPEC_ATOM_EXI'; call gtv(trim(nm),tksw(prgn,nm),exi(:,j),
     .  cindx=jj,note='Hankel energies for interstitial product basis'
     .  ,nout=nout)
        if (nout .ge. 0) nxi(j) = nout
        nm='SPEC_ATOM_RINT'; call gtv(trim(nm),tksw(prgn,nm),rint(j),
     .  cindx=jj,note='range of interstitial product basis')
        nm='SPEC_ATOM_RCUT'; call gtv(trim(nm),tksw(prgn,nm),rcut(j),
     .  cindx=jj,note='radius over which true,smooth Hankels differ')

C       Some sanity checks
        if (io_help == 0 .and. lmxaj .ge. 0) then
          if (tksw(prgn,'SPEC_ATOM_LFOCA') /= 2)
     .    call isanrg(lfoca(j),0,2,'rdctrl','lfoca',T)
          if (tksw(prgn,'SPEC_ATOM_LMXL') /= 2) call
     .    isanrg(lmxl(j),min(0,lmxaj),max(0,lmxaj),'rdctrl','lmxl',T)
          if (tksw(prgn,'SPEC_ATOM_KMXA') /= 2)
     .    call isanrg(kmxt(j),2,25,' rdctrl (warning):','kmxa',F)
        endif

c        i0 = partks(recrd(catbeg),'C-HOLE ',' ','spec coreh',sspec,j,
c     .    'channel for core hole',1,1,i1,TF(19))
c        i0 = partks(recrd(catbeg),'C-HQ=','= ,','spec coreq',sspec,j,
c     .    'Charge in core hole.  '//
c     .    'Optional 2nd entry is moment of core hole:'//
c     .    '%N%12fQ(spin1) = full + C-HQ(1)/2 + C-HQ(2)/2'//
c     .    '%N%12fQ(spin2) = full + C-HQ(1)/2 - C-HQ(2)/2',2,4,i1,TF(19))
c        i0 = partks(recrd(catbeg),'EREF=','=','spec eref',sspec,j,
c     .    'Reference energy subtracted from total energy',1,4,i1,TF(4))

        coreh(j) = ' '
        nm='SPEC_ATOM_C-HOLE'; call gtv(trim(nm),tksw(prgn,nm),coreh(j),
     .  nmin=10,cindx=jj,note='Channel for core hole')
        nm='SPEC_ATOM_C-HQ'; call gtv(trim(nm),tksw(prgn,nm),coreq(:,j),
     .  def_r8v=(/-1d0,0d0/),cindx=jj,nmin=2,note=
     .  'Charge in core hole.  '//
     .  'Optional 2nd entry is moment of core hole:'//
     .  '%N%5fQ(spin1) = full + C-HQ(1)/2 + C-HQ(2)/2'//
     .  '%N%5fQ(spin2) = full + C-HQ(1)/2 - C-HQ(2)/2')

        nm='SPEC_ATOM_EREF'; call gtv(trim(nm),tksw(prgn,nm),eref(j),
     .  def_r8=0d0,cindx=jj,note=
     .  'Reference energy subtracted from total energy')

c        nchoos = 2
c        call upack('spec mass',sspec,j,xx,0,0,0) !I think mass==0 here
c        noerr = .false.
c        if (xx .eq. 0) then
c          i0 = partok(recrd(catbeg),'ISOTOPE=','=',xx,
c     .      'relative atomic mass (for dynamics)',-1,4,i1,
c     .      TF(17))
c          if (i0 .eq. 1) then
c            ncor = 0
c            i0 = partok(recrd(catbeg),'NCOR=','=',ncor,' ',-1,2,i1,F)
cC       ... Work out atomic mass from neutron, proton and electron masses
c            nmass = 1836.153d0/2d0
c            pmass = 1838.684d0/2d0
c            emass = 0.5d0
c            xx = z*pmass + (xx-z)*nmass + ncor*emass
c            call pack5('spec mass',sspec,j,xx,0,0,0)
c            noerr = .true.
c          endif

        nm='SPEC_ATOM_AMASS'; call gtv(trim(nm),tksw(prgn,nm),mass(j),
     .  cindx=jj,note='Nuclear mass in a.u. (for dynamics)')

C        Needs to be cleaned up
C        sw = tksw(prgn,'SPEC_ATOM_ISOTOPE')
C        if (sw/=2 .or. io_help) then
C          sw=0
C          nm='SPEC_ATOM_ISOTOPE'; call gtv(trim(nm),sw,
C     .      xx,cindx=jj, nmin=0,note=
C     .      'Relative atomic mass (for dynamics)',nout=nout, or=.true.)  !,-1,4,i1, TF(17))
C          if (nout>0) then
C            ncor = 0  !??? default???
C            nm='SPEC_ATOM_NCOR'; call gtv(trim(nm),sw,
C     .      ncor,cindx=jj ) !' ',-1,2,i1,F)
CC       ... Work out atomic mass from neutron, proton and electron masses
C            nmass = 1836.153d0/2d0
C            pmass = 1838.684d0/2d0
C            emass = 0.5d0
Cc            xx = z*pmass + (xx-z)*nmass + ncor*emass
Cc            call pack5('spec mass',sspec,j,xx,0,0,0)
C            mass(j)=z(j)*pmass + (xx-z(j))*nmass + ncor*emass
Cc            noerr = .true.
C          else
C          endif
C        endif

        nm='SPEC_ATOM_COLOUR'; call gtv(trim(nm),tksw(prgn,nm),
     .  colxbs(:,j),def_r8v=zerov,cindx=jj,note=
     .  'Colour for xbs') !,-3,4,i1,TF(37))
        xx = rmt(j); if (io_help .eq. 1) xx = NULLI
        nm='SPEC_ATOM_RADIUS'; call gtv(trim(nm),tksw(prgn,nm),
     .  radxbs(j),def_r8=xx,cindx=jj,note='Radius for xbs')
      enddo                     ! Loop over species
C ... Cleanup after looping over species data
      if (io_help==0) then
C       Maximum L-cutoff
        call imxmn(nspec,lmxa,1,i0,lmxax)
        nlmax = (lmxax+1)**2
C       Set global variables nkaph and mxorb
        nkaph = nkaph+lpzi-1
C       print *, nkaph,nlmax,nkaph*nlmax
Changenglob        xx = dglob('nkaph',dble(nkaph),1)
        globalvariables%nkaph = dble(nkaph); globalvariables%l_nkaph = globalvariables%l_nkaph +1; xx = dble(nkaph)
Changenglob        xx = dglob('nlmax',dble(nlmax),1)
        globalvariables%nlmax = dble(nlmax); globalvariables%l_nlmax = globalvariables%l_nlmax +1; xx = dble(nlmax)
Changenglob        xx = dglob('mxorb',dble(nkaph)*nlmax,1)
        globalvariables%mxorb = dble(nkaph)*nlmax; globalvariables%l_mxorb = globalvariables%l_mxorb +1; xx = dble(nkaph)*nlmax
      endif
   79 continue

C --- Site ---
      if (io_show+io_help/=0 .and. tksw(prgn,'SITE')/=2) call info0(1,1,
     .0,' --- Parameters for site data ---')

c      print *, 'goto Site in m_rdctrl'

C     SITE call is only helpful for printout
C      nm='SITE'; sw = tksw    (prgn,nm); if (sw == 2) goto 89
C      if (io_help /= 0 .or. io_show>0) write(stdo,'(1x)')
C      call gtv(trim(nm),sw,Texist=ltmp,note='Site-specific parameters')

      if (io_help .eq. 1) then
        nbas = 1
        nsite =2
      endif
      allocate(pos(3,nsite),vel(3,nsite),eula(3,nsite),vshft(nsite),
     .ips(nsite),ipl(nsite),plv(nsite),irlx(3,nsite),
     .mpole(nsite),dpole(3,nsite))
c      if (ltbe) allocate(delta(n0,nsite),ndelta(nsite))

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
      mpole = 0d0
      dpole = 0d0

      nm='SITE_FILE'; call gtv(trim(nm),tksw(prgn,nm),outs,nmin=10,
     .nout=nout, note=
     .'Name of site file containing basis and lattice information.'//
     .'%N%3fIf this token is specified, read POS from site file,'//
     .'%N%3fand VEC, EULA, VSHFT, PL, RLX, if they are present.'
     .//'%N%3fIn this case no attempt is made to read the SITE '
     .//'tokens below.')

      if (nout == 1) then
C       Get FILE lio to see what it contains
        if (procid .eq. master) then
          j = iosite(16000,3d0,0,trim(outs),i,slabl,alat,plat,nbas,
     .    nspec,xx,xx,xx,xx,xx,xx,xx)
          ii = 8000 + 140
C         File contains vel,eula,PL,rlx,vshft
          if (mod(j/32,2) .eq. 1) then
            ii = ii + 32000
          endif
c          if (ltbe) ii = ii + 64000
          j = iosite(ii,3d0,0,trim(outs),i,slabl,alat,plat,nbas,
     .    nspec,pos,vel,eula,vshft,ips,ipl,irlx)
        endif
        call mpibc1(pos,3*nsite,4,mlog,'readctrl','pos')
        call mpibc1(vel,3*nsite,4,mlog,'readctrl','vel')
        call mpibc1(eula,3*nsite,4,mlog,'readctrl','eula')
        call mpibc1(vshft,nsite,4,mlog,'readctrl','vshft')
        call mpibc1(ips,nsite,2,mlog,'readctrl','ips')
        call mpibc1(ipl,nsite,2,mlog,'readctrl','ipl')
        call mpibc1(irlx,3*nsite,2,mlog,'readctrl','irlx')
        goto 89
      endif

      if (io_help .eq. 1) then
        if (iprint() .gt. 0) write(*,383)
  383   format(/' * ',
     .  'The following tokens are input for each site. ',
     .  'Data sandwiched'/3x,'between successive occurences of ',
     .  'token ATOM apply to one site.'/3x,'Alternatively, all ',
     .  'site data can be read in via the SITE file.')
        nbas = 1
        nsite = 2
      endif

C ... Site data, one pass for each atom
      do  j = 1, nbas

        if (io_help /= 0) then
          write(stdo,'(1x)')
        elseif (io_help == 0 .and. io_show>0) then
          call info(1,0,0,' ... Site %i',j,0)
        endif

        jj=(/1,j/)
        nm='SITE_ATOM'; call gtv(trim(nm),tksw(prgn,nm),alabl,nmin=10,
     .  cindx=jj,note='Species label')
        if (io_help .ne. 1) then
          do  i = 1, nspec
c             print *, alabl,' xxx ', slabl(i)
            if (trim(alabl) == trim(slabl(i)) ) then
              ips(j) = i
              print *, alabl,' xxx ', i,ips(j)
              goto 881
            endif
          enddo
          call rxs('Category SITE referred to'//
     .    ' nonexistent species: ',alabl)
        endif
  881   continue

C       call snit
C
C   ... Site positions
        sw = tksw(prgn,'SITE_ATOM_XPOS')
        nm='SITE_ATOM_POS'; call gtv(trim(nm),tksw(prgn,nm),pos(:,j),
     .  nout=nout,cindx=jj,note='Atom coordinates, in units of alat',
     .  or=(sw.ne.2))
        if (nout == 0 .or. tksw(prgn,'SITE_ATOM_POS') == 2) then !nout=-1 if sw=2; otherwise nout=0 unless data was read
          nm='SITE_ATOM_XPOS'; call gtv(trim(nm),tksw(prgn,nm),pos(:,j),
     .    cindx=jj,note='Atom coordinates, as (fractional) '//
     .    'multiples of the lattice vectors')
          call dcopy(3,pos(1,j),1,xv,1)
          call dmpy(plat,3,1,xv,3,1,pos(1,j),3,1,3,1,3)
        endif
        nm='SITE_ATOM_DPOS'; call gtv(trim(nm),tksw(prgn,nm),
     .  xv(1:3),def_r8v=zerov,nout=nout,cindx=jj,note=
     .  'Shift in atom coordinates added to pos')
        if (nout .eq. 3) call daxpy(3,1d0,xv,1,pos(1,j),1)

        nm='SITE_ATOM_V0'; call gtv(trim(nm),tksw(prgn,nm),vel(:,j),
     .  def_r8v=zerov,cindx=jj,
     .  note='Initial velocity for molecular dynamics')

        nm='SITE_ATOM_RELAX'; call gtv(trim(nm),tksw(prgn,nm),irlx(:,j),
     .  def_i4v=(/(1,i=1,n0)/),cindx=jj,note=
     .  'relax site positions (lattice dynamics) or '//
     .  'Euler angles (spin dynamics)')

        nm='SITE_ATOM_VSHFT'; call gtv(trim(nm),tksw(prgn,nm),vshft(j),
     .  def_r8=0d0,cindx=jj,
     .  note='Constant potential shift for this site')

        nm='SITE_ATOM_ROT'; call gtv(trim(nm),tksw(prgn,nm),outs,nmin=
     .  10,cindx=jj,nout=nout,note=
     .  'Rotation of spin quantization axis at this site')
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
     .    '%15pROT= : Euler alpha = %1;6d  beta = %1;6d'//
     .    '  gamma = %1;6d',eula(1,j),eula(2,j),eula(3,j),0,0)
        endif

        nm='SITE_ATOM_PL'; call gtv(trim(nm),tksw(prgn,nm),ipl(j),
     .  def_i4=0,cindx=jj,
     .  note='Assign principal layer number to this site')

        nm='SITE_ATOM_PLV'; call gtv(trim(nm),tksw(prgn,nm),plv(j),
     .  def_i4=0,cindx=jj,
     .  note='Assign PL potential index to this site')

C$$$        nm='SITE_ATOM_DELTA'; sw = tksw(prgn,nm)
C$$$        if (sw /= 2 .and. ltbe) then
C$$$          if (io_help .eq. 1) then
C$$$            nlbji = NULLI
C$$$            nlbj = 1
C$$$          elseif (lmxb(i) == NULLI) then
C$$$            nlbji = NULLI
C$$$            nlbj = 0
C$$$          else
C$$$            i = ips(j)
C$$$            nlbj = 1+lmxb(i)
C$$$            nlbji = nlbj
C$$$          endif
C$$$          ndelta(j)=nlbj*nsp
C$$$          call gtv(trim(nm),sw,delta(1:nlbj*nsp,j),nmin=nlbji*nsp,
C$$$     .      cindx=jj,def_r8v=zerov,
C$$$     .      note='Vector of on-site energy shifts')
C$$$        endif

C        nm='SITE_ATOM_SID'; call gtv(trim(nm),tksw(prgn,nm),sid(j),
C     .    def_i4=0,cindx=jj,note=
C     .    'Site ID, used in the merging of coincident clusters')
C
      enddo

C ... Input for point multipoles
      sw = tksw(prgn,'SITE_PM')
      if (sw .lt. 2) then
        if (io_help == 1) nbasp=nbas+1
        if (nbasp .gt. nbas) then
          if (io_show>0) write
     .    (stdo,'(/'' The following input is for point multipoles:'')')
          do j = nbas+1, nbasp

            jj=(/1,j-nbas/)
C          nm='SITE_PM_POS'; call gtv(trim(nm),tksw(prgn,nm),pos(:,j),
C     .      nout=nout,cindx=jj,note='PM coordinates, in units of alat')

            nm='SITE_PM_POS'; call gtv(trim(nm),tksw(prgn,nm),pos(:,j),
     .      cindx=jj,note='PM coordinates, in units of alat')

            nm='SITE_PM_Q'; call gtv(trim(nm),tksw(prgn,nm),mpole(j),
     .      def_r8=0d0,cindx=jj,note='Point monopole moment')

            nm='SITE_PM_P'; call gtv(trim(nm),tksw(prgn,nm),dpole(:,j),
     .      def_r8v=zerov,cindx=jj,note='Point dipole moment')
          enddo
        endif
      endif

   89 continue

C --- Structure constants ---
      if (io_show+io_help/=0 .and. tksw(prgn,'STR')/=2) call info0(1,1,
     .0,' --- Parameters for structure constants ---')
      nm='STR_RMAXS'; call gtv(trim(nm),tksw(prgn,nm),str_rmax,
     .nout=nout,note='Radial cutoff for strux, in a.u.',or=T)
      if (nout == 0) then       !nout=-1 if sw=2; otherwise nout=0 unless data was read
        nm='STR_RMAX'; call gtv(trim(nm),tksw(prgn,nm),str_rmax,
     .  def_r8=0d0,note='Radial cutoff for strux, in units of avw')
        str_rmax = str_rmax*avw
      endif

C     Energies for NMTO structure constants
c      nm='STR_ENV_MODE'; call gtv(trim(nm),tksw(prgn,nm),str_mode,
c     .  def_i4=0,note='Type of envelope functions:'//
c     .  '%N%3f0 2nd generation'//
c     .  '%N%3f1 SSSW (3rd generation)'//
c     .  '%N%3f2 NMTO'//
c     .  '%N%3f3 SSSW and val-lap basis')

c      nkapsi = 1
c      if (str_mode == 2) nkapsi=2
c      if (io_help .eq. 1) nkapsi = NULLI
c      if (str_mode .eq. 0) then
c        str_nkaps = 1
c      else
c        nm='STR_ENV_NEL'; call gtv(trim(nm),tksw(prgn,nm),str_nkaps,
c     .    def_i4=nkapsi,note='Number of SSSW or NMTO energies')
c        if (str_mode .eq. 2)
c     .    call isanrg(str_nkaps,2,5,'rdctrl','NEL',.true.)
c        nkapsi = str_nkaps
c        if (io_help .eq. 1) then
c          str_nkaps = 1
c          nkapsi = NULLI
c        endif
c        nm='STR_ENV_EL';call gtv(trim(nm),tksw(prgn,nm),
c     .    str_kaps(1:str_nkaps),nmin=nkapsi,def_r8v=(/0d0,-1d0,2.3d0/),
c     .    note='SSSW or NMTO energies, atomic units')
c      endif

      nm='STR_MXNBR'; call gtv(trim(nm),tksw(prgn,nm),str_mxnbr,
     .def_i4=0,note='Max number of nbrs (for dimensioning arrays)')
c      nm='STR_SHOW'; call gtv(trim(nm),tksw(prgn,nm),str_lshow1,
c     .  def_lg=F,note='Show strux after generating them')
c      nm='STR_EQUIV'; call gtv(trim(nm),tksw(prgn,nm),str_lequiv1,
c     .  def_lg=F,note='Look for equivalent strux')

c      nm='STR_LMAXW'; call gtv(trim(nm),tksw(prgn,nm),str_lmaxw,
c     .  def_i4=-1,note=
c     .  'l-cutoff for Watson sphere, used to help localize strux')
c      nm='STR_DELRW'; call gtv(trim(nm),tksw(prgn,nm),str_drwats,
c     .  def_r8=0.1d0,note='Padding beyond cluster for Watson sphere')

c      if (io_help /= 0 .and. tksw(prgn,'STR_IINV_NIT')/=2)
c     .  call info0(1,1,0,' * IINV parameters '//
c     .  'govern iterative solutions to screened strux')
c      nm='STR_IINV_NIT'; call gtv(trim(nm),tksw(prgn,nm),iinv_nit,
c     .  def_i4=0,note='Number of iterations')
c      nm='STR_IINV_NCUT'; call gtv(trim(nm),tksw(prgn,nm),iinv_ncut,
c     .  def_i4=0,note='Number of sites for inner block')
c      nm='STR_IINV_TOL'; call gtv(trim(nm),tksw(prgn,nm),iinv_tol,
c     .  def_r8=0d0,note='Tolerance in errors')

c      nm='STR_RVL/R'; call gtv(trim(nm),tksw(prgn,nm),str_rmaxg,
c     .  def_r8=0.7d0,note='Radial cutoff for val-lap, in units of RMAX')

c      nm='STR_VLFUN'; call gtv(trim(nm),tksw(prgn,nm),str_ivl,
c     .  def_i4=0,note='Functions for val-lap basis'//
c     .  '%N%10f0 G0 + G1'//
c     .  '%N%10f1 G0 + Hsm'//
c     .  '%N%10f2 G0 + Hsm-dot')

C$$$C --- Two-center fit ---
C$$$      if (io_help .eq. 1 .and. tksw(prgn,'TCF')/=2) call info(1,1,0,
C$$$     .  ' * The following tokens are used for the 2-center fit',0,0)
C$$$      nm='TCF_NBISI'; call gtv(trim(nm),tksw(prgn,nm),tcf_nbisi,
C$$$     .  note='Mesh parameters for TCF ')
C$$$      nm='TCF_NALF'; call gtv(trim(nm),tksw(prgn,nm),tcf_nalf,
C$$$     .  note=
C$$$     .  'Polynomial order for d-dependence of TCF fit coffs ')
C$$$      nm='TCF_NCUPL'; call gtv(trim(nm),tksw(prgn,nm),tcf_ncupl,
C$$$     .  def_i4=8,note='Reduced polynomial order in TCF ')
C$$$      nm='TCF_NDUST'; call gtv(trim(nm),tksw(prgn,nm),tcf_ndust,
C$$$     .  def_i4=2,note='No. iterations to improve on TCF ')
C$$$      nm='TCF_ADEC'; call gtv(trim(nm),tksw(prgn,nm),tcf_adec,
C$$$     .  def_r8=1d0,note='Log spacing in distance for TCF ')
C$$$      nm='TCF_WZTCF'; call gtv(trim(nm),tksw(prgn,nm),tcf_wztcf,
C$$$     .  def_r8=0d0,note='Weighting of z-points in TCF ')

C --- Brillouin Zone ---
      if (io_show+io_help/=0 .and. tksw(prgn,'BZ')/=2) call info0(1,1,
     .0,' --- Parameters for Brillouin zone integration ---')

      nm='BZ_GETQP'; call gtv(trim(nm),tksw(prgn,nm),bz_lio1,
     .def_lg=F,note='Read qp from disk',or=T)
      nm='BZ_NKABC'; sw=tksw(prgn,nm); if (bz_lio1) sw = 2
      call gtv(trim(nm),sw,bz_nabc,nout=nout,
     .note='No. qp along each of 3 lattice vectors.'//
     .'%N%3fSupply one number for all vectors or a separate '//
     .'number for each vector.')
      call fill3in(nout,bz_nabc)

      nm='BZ_PUTQP'; call gtv(trim(nm),tksw(prgn,nm),bz_lio2,
     .def_lg=F,note='Write qp to disk')

      nm='BZ_BZJOB';call gtv(trim(nm),tksw(prgn,nm),bz_lshft,nout=nout,
     .def_i4v=izerv(1:1),note=
     .'0 centers BZ mesh at origin, 1 centers off origin'//
     .'%N%3fSupply one number for all vectors or a separate '//
     .'number for each vector.')
      call fill3in(nout,bz_lshft)

      nm='BZ_METAL'; call gtv(trim(nm),tksw(prgn,nm),bz_lmet,
     .def_i4=1,note=
     .'0 assume insulator; 1 save evecs on disk; '//
     .'2 use wgt from prior iter%N%3f3 always make two band passes;'//
     .' 4 BZ integration with 3-point scheme')

      nm='BZ_TETRA'; call gtv(trim(nm),tksw(prgn,nm),ctrl_lmet2,
     .def_lg=T,note='Tetrahedron integration')

      nm='BZ_N'; call gtv(trim(nm),tksw(prgn,nm),bz_n,
     .def_i4=0,note=
     .'N>0: Polynomial order for Methfessel-Paxton sampling%N%3f'//
     .'N=0: Conventional Gaussian sampling%N%3f'//
     .'N<0: Broadening by Fermi-Dirac distribution%N%3f'//
     .'To be used in conjunction with W= ; see next')
      nm='BZ_W'; call gtv(trim(nm),tksw(prgn,nm),bz_w,
     .def_r8=5d-3,note=
     .'N>=0: Line broadening for sampling integration%N%3f'//
     .'N<0 : Temperature for Fermi distribution (Ry)')
      nm='BZ_EF0'; call gtv(trim(nm),tksw(prgn,nm),bz_ef,
     .def_r8=0d0,note='Initial guess at Fermi energy')
      nm='BZ_DELEF'; call gtv(trim(nm),tksw(prgn,nm),bz_def,
     .def_r8=0.05d0,note='Initial uncertainty in Fermi energy')
      nm='BZ_ZBAK'; call gtv(trim(nm),tksw(prgn,nm),zbak(1),
     .def_r8=0d0,note='Homogeneous background charge')
      nm='BZ_SAVDOS'; call gtv(trim(nm),tksw(prgn,nm),ctrl_ldos,
     .def_i4=0,note='Choose some combination of the following:'
     .//'%N%3f1 Write DOS to directly disk (NPTS and DOS also needed)'
     .//'%N%3f2 Write weights for partial DOS'
     .//'%N%3f4 Same as (2), but weights m-resolved')
      nm='BZ_DOS'; call gtv(trim(nm),tksw(prgn,nm),bz_dosw,def_r8v=
     .(/-1d0,0d0/),note='Energy window over which DOS accumulated')
      nm='BZ_NPTS'; call gtv(trim(nm),tksw(prgn,nm),bz_ndos,def_i4=1001,
     .note='No. DOS points (sampling integration, and lmdos)')

      xx = 5d0
c      if (ltbe) xx = 5d0
      nm='BZ_EFMAX'; call gtv(trim(nm),tksw(prgn,nm),bz_efmax,
     .def_r8=xx,note='Find evecs up to efmax')
      nm='BZ_NEVMX'; call gtv(trim(nm),tksw(prgn,nm),bz_nevmx,
     .def_i4=0,note='Find at most nevmx eigenvectors'//
     .'%N%3fIf NEVMX=0, program uses internal default'//
     .'%N%3fIf NEVMX<0, no eigenvectors are generated')
      nm='BZ_ZVAL'; call gtv(trim(nm),tksw(prgn,nm),bz_zval,def_r8=0d0,
     .note='Number of electrons to accumulate in BZ integration')
      nm='BZ_NOINV'; call gtv(trim(nm),tksw(prgn,nm),ctrl_lqp1,
     .def_lg=F,note=
     .'Suppress automatic inclusion of inversion symmetry for BZ')
      nm='BZ_FSMOM'; call gtv(trim(nm),tksw(prgn,nm),bz_fsmom,
     .def_r8=0d0,note='Fixed-spin moment (fixed-spin moment method)')
c      nm='BZ_DMAT'; call gtv(trim(nm),tksw(prgn,nm),bz_lio8,
c     .  def_lg=F,note='Calculate density matrix')
      nm='BZ_INVIT'; call gtv(trim(nm),tksw(prgn,nm),ctrl_lqp2,
     .def_lg=T,note='Use inverse iteration for diagonalization')
      i0 = 0
c      if (ltbe) i0 = -1
      nm='BZ_MULL'; call gtv(trim(nm),tksw(prgn,nm),bz_lmull,
     .def_i4=i0,note='Mulliken population analysis ')
      nm='BZ_EMESH'; call gtv(trim(nm),tksw(prgn,nm),bz_semsh,nmin=10,
     .def_r8v=(/10d0,0d0,-1d0,0d0,0.01d0,0d0,0d0,0d0,1d-10,0d0/),
     .note='Mesh for energy contour integration'
     .//'%N%3fentry 1   2    3    4    5   6    7    8    9    10'
     .//'%N%9fnz mode ebot etop  ----  depend on mode ---- '
     .//'(see documentation)')

      nm='BZ_COND_MODE'; call gtv(trim(nm),tksw(prgn,nm),i0,nout=nout,
     .def_i4=0,note='Mode for type conductivity calculation'//
     .'%N%3f0: Calculate DOS D(E) instead of conductivity'//
     .'%N%3f1: ballistic conductivity <v.DIR> (must also supply DIR)'
     .//
     .'%N%3f2: <v_i . v_j> (diffusive conductivity)'//
     .'%N%6fIn this mode one or two elts of DIR should be 1, '//
     .'marking i or i,j')
      bz_lcond(1) = i0
      nm='BZ_COND_DIR'
      if (io_help /= 0 .or. i0 .gt. 0) then
        sw = tksw(prgn,nm)
      else
        sw = 2
      endif
      call gtv(trim(nm),sw,bz_lcond(2:),
     .note='Conductivity direction vector DIR, used for MODE>0')

C --- Ewald sums ---
      if (io_show+io_help/=0 .and. tksw(prgn,'EWALD')/=2)
     .call info0(1,1,0,' --- Parameters for Ewald sums ---')

      nm='EWALD_AS'; call gtv(trim(nm),tksw(prgn,nm),lat_as,
     .def_r8=2d0,note='Ewald smoothing parameter')
      nm='EWALD_TOL'; call gtv(trim(nm),tksw(prgn,nm),lat_tol,
     .def_r8=1d-8,note='Ewald tolerance')
      nm='EWALD_NKDMX'; call gtv(trim(nm),tksw(prgn,nm),lat_nkdmx,
     .def_i4=800,note='Ewald tolerance')
      nm='EWALD_RPAD'; call gtv(trim(nm),tksw(prgn,nm),lat_rpad,
     .def_r8=0d0,note='Scale rcutoff by rpad when lattice vectors'//
     .' padded in oblong geometries')

C --- Iterations (formerly MIX) ---
      vmix%b = NULLI           ! Not set
      if (tksw(prgn,'ITER')/=2) then
        if (io_show+io_help/=0)
     .  call info0(1,1,0,' --- Parameters for iterations ---')
C     Default values for smix (array has same same structure as lstra smix)
        call umix_init(vmix)
        vmix%b = 1               ! beta
        vmix%bv = 1               ! bv
C     vmix%elind =  0              ! elind
        call s8tor8('mixm',vmix%fn) ! file name
        vmix%kill =  0              ! nkill
        vmix%lxpot =  0              ! lxpot
        vmix%mmix = -1              ! mmix
        vmix%mode = 0              ! mode (0=Anderson)
C     vmix%model = 0              ! model = previous mixing mode
C     vmix%n = 0              ! n     = Number of iterations for this species
C     vmix%nitu = 1              ! nitu  = max number of LDA+U iterations
C     vmix%nmix = 1              ! nmix  = actual number mixed
        vmix%nsave = 8              ! nsave = # iter to save on disk
        call s8tor8(' ',vmix%r) ! expression for rmscst
        vmix%tolu = 0              ! tolu
        vmix%umix = 1              ! umix (mixing parm for LDA+U)
        vmix%w(1) = 1              ! w(1)
        vmix%w(2) = 1              ! w(2)
        vmix%wc = -1             ! wc
        smalit = NULLI
c      call dcopy(37,smix,1,samix,1)
        vamix = vmix
        call s8tor8('mixa',vamix%fn) ! file name

        nm='ITER_NIT';call gtv(trim(nm),tksw(prgn,nm),iter_maxit,
     .  def_i4=30,
     .  note='maximum number of iterations in self-consistency cycle')

        nm='ITER_NRMIX'; call gtv(trim(nm),tksw(prgn,nm),smalit,
     .  def_i4v=(/80,2/),note='Sphere program, (1) max iter; '//
     .  '(2) no. prior iter for Anderson mixing ')

        nm='ITER_MIX'; sw=tksw(prgn,nm); call gtv(trim(nm),sw,iter_mix,
     .  nmin=10,nout=nout,note='Mixing rules for charge mixing.  Syntax:')
        if (io_help/=0 .and. tksw(prgn,nm)/=2) print 345
  345   format(3x,
     .  'A[nmix][,b=beta][,bv=betv][,n=nit][,w=w1,w2][,nam=fn][,k=nkill]',
     .  '[;...] or'/
     .  3x,'B[nmix][,b=beta][,bv=betv][,wc=wc][,n=#][,w=w1,w2][,nam=fn]',
     .  '[,k=nkill]')
C ... Sanity check and printout

cc takao
c        mmm=adjustl(iter_mix)
c        if(mmm(1:1)=='B')
c     &      call rx(" Broyden mixing is problematic."//
c     &              " If necessary, you have to fix it (mixrho.F)")

        if (nout==1 .and. io_help==0 .and. sw/=2) then
          call r8tos8(vmix%fn,alabl)
          if (io_show==0) call pshpr(0)
          if (.not. parmxp(-1,iter_mix,len_trim(iter_mix),
     .    int(vmix%mode),int(vmix%nmix),vmix%w(1),vmix%b,vmix%elind,
     .    alabl,vmix%wc,vmix%kill,vmix%bv,xx)) then
            if (io_show==0) call poppr
            call rx( 'RDCCAT: parse in parmxp failed 1')
          endif
          if (io_show==0) call poppr
        endif

        nm='ITER_AMIX'; sw=tksw(prgn,nm); call gtv(trim(nm),sw,iter_amix,
     .  nmin=10,nout=nout,note=
     .  'Mixing rules for Euler angles.  Syntax as in MIX')
C ... Sanity check and printout
        if (nout==1 .and. io_help==0 .and. sw/=2) then
          call r8tos8(vamix%fn,alabl)
          if (io_show==0) call pshpr(0)
          if (.not. parmxp(-1,iter_amix,len_trim(iter_amix),
     .    int(vamix%mode),int(vamix%nmix),vamix%w(1),vamix%b,vamix%elind,
     .    alabl,vamix%wc,vamix%kill,vamix%bv,xx)) then
            if (io_show==0) call poppr
            call rx('READCTRL: parse in parmxp failed 2')
          endif
          if (io_show==0) call poppr
        endif

c      if (prgn /= 'LMMC') then
        nm='ITER_CONV'; call gtv(trim(nm),tksw(prgn,nm),ctrl_tol(3),
     .  def_r8=1d-4,note='Tolerance in energy change from prior'//
     .  ' iteration for self-consistency')
        nm='ITER_CONVC'; call gtv(trim(nm),tksw(prgn,nm),ctrl_tol(1),
     .  def_r8=1d-4,note='Tolerance in output-input charge for'//
     .  ' self-consistency')
c      else
cC ... CONV and CONVC are stored in the opposite order in LMMC
c        nm='ITER_CONV'; call gtv(trim(nm),tksw(prgn,nm),ctrl_tol(1),
c     .    def_r8=1d-4,note='Tolerance in energy change from prior'//
c     .    ' iteration for self-consistency')
c        nm='ITER_CONVC'; call gtv(trim(nm),tksw(prgn,nm),ctrl_tol(3),
c     .    def_r8=1d-4,note='Tolerance in output-input charge for'//
c     .    ' self-consistency')
c        ctrl_tol(2) = ctrl_tol(3)
c       nm='ITER_QTOLSP'; call gtv(trim(nm),tksw(prgn,nm),ctrl_tol(2),
c    .    def_r8=1d-4,note='Tolerance in sphere charge'//
c    .    ' for self-consistency')
c       nm='ITER_QTOLI'; call gtv(trim(nm),tksw(prgn,nm),ctrl_tol(3),
c    .    def_r8=1d-4,note='Tolerance in interstitial charge'//
c    .    ' for self-consistency')
c      endif

        i0 = 0
        nm='ITER_XIPMX'; call gtv(trim(nm),tksw(prgn,nm),i0,def_i4=0,
     .  nout=nout,note='Mix potential independently of charge:'//
     .  '%N%3fXIPMX=1: mix vin and v(qmix)'//
     .  '%N%3fXIPMX=2: mix vin and v(qout)')
        if (nout .gt. 0) then
          call isanrg(i0,0,3,'rdctrl','XIPMX',T)
          if (i0 .gt. 0) lves = 1 !it was lves= .true., but this is a bug . sctrl/lves is integer(8).
        endif
        vmix%lxpot = i0

        nm='ITER_UMIX'; call gtv(trim(nm),tksw(prgn,nm),vmix%umix,
     .  def_r8=1d0,note='Mixing parameter for densmat in LDA+U')
        nm='ITER_TOLU'; call gtv(trim(nm),tksw(prgn,nm),vmix%tolu,
     .  def_r8=0d0,note='Tolerance for densmat in LDA+U')
        nm='ITER_NITU'; call gtv(trim(nm),tksw(prgn,nm),i0,
     .  def_i4=0,note='Max number of LDA+U iterations of densmat')
        vmix%nitu = i0

      endif                     ! iterations category

C$$$C --- TB ---
C$$$      if (tksw(prgn,'TB') /= 2) then
C$$$      if (io_show+io_help/=0)
C$$$     .  call info0(1,1,0,' --- Tight binding parameters ---')
C$$$      nm='TB_OVLP'; call gtv(trim(nm),tksw(prgn,nm),ltb1,
C$$$     .  def_lg=F,note='Non orthogonal tight-binding ')
C$$$      nm='TB_CRYSF'; call gtv(trim(nm),tksw(prgn,nm),ltb2,
C$$$     .  def_lg=F,note='Crystal field terms in hamiltonian ')
C$$$      nm='TB_OVCF'; call gtv(trim(nm),tksw(prgn,nm),ltb4,
C$$$     .  def_lg=F,note='Crystal field terms in overlap ')
C$$$      nm='TB_ADDES'; call gtv(trim(nm),tksw(prgn,nm),ltb8,
C$$$     .  def_lg=F,note='Add ebarLL'' * sLL'' to hLL''')
C$$$      nm='TB_RMAXH'; call gtv(trim(nm),tksw(prgn,nm),str_rmax,
C$$$     .  def_r8=0d0,note='Hamiltonian cut-off length in units of a')
C$$$      nm='TB_FORCES'; call gtv(trim(nm),tksw(prgn,nm),ltb16,
C$$$     .  def_lg=T,note='Calculate forces ')
C$$$      nm='TB_FIJ'; call gtv(trim(nm),tksw(prgn,nm),ltb32,
C$$$     .  def_lg=F,note='To get forces by atom pair, e.g., for stresses ')
C$$$      nm='TB_3PV'; call gtv(trim(nm),tksw(prgn,nm),ltb128,
C$$$     .  def_lg=F,note='Calculate pressure')
C$$$      nm='TB_EVDISC'; call gtv(trim(nm),tksw(prgn,nm),ltb256,
C$$$     .  def_lg=T,note='Can be F '//
C$$$     .  'for insulators or to save space for metals ')
C$$$      nm='TB_PAIR'; call gtv(trim(nm),tksw(prgn,nm),ltb512,
C$$$     .  def_lg=F,note='Pair potential only')
C$$$      nm='TB_TRH'; call gtv(trim(nm),tksw(prgn,nm),ltb1024,
C$$$     .  def_lg=F,note='Calculate local projection of band energies')
C$$$      nm='TB_RHO'; call gtv(trim(nm),tksw(prgn,nm),ltb2048,
C$$$     .  def_lg=F,note='Calculate local projection of charges')
C$$$      nm='TB_U1'; call gtv(trim(nm),tksw(prgn,nm),ltb4096,
C$$$     .  def_lg=F,note='For electrostatics with L>=0 (first model)')
C$$$      nm='TB_UL'; call gtv(trim(nm),tksw(prgn,nm),ltb215,
C$$$     .  def_lg=F,note='For electrostatics with L>=0')
C$$$      nm='TB_TBU'; call gtv(trim(nm),tksw(prgn,nm),ltb213,
C$$$     .  def_lg=F,note='Tight-binding+U')
C$$$      nm='TB_GAMMA'; call gtv(trim(nm),tksw(prgn,nm),ltb217,
C$$$     .  def_lg=F,note='Do gamma-point only')
C$$$      nm='TB_MOL'; call gtv(trim(nm),tksw(prgn,nm),ltb218,
C$$$     .  def_lg=F,note='molecule: no PBCs')
C$$$      if (io_help .eq. 1) call info0(1,1,0,
C$$$     .  ' * if U1=T, UL=T or TBU=T, the next two tokens are used:')
C$$$      nm='TB_NOUAVG'; call gtv(trim(nm),tksw(prgn,nm),ltb214,
C$$$     .  def_lg=F,note='Use individual U_l from Q=; don''t average them')
C$$$      nm='TB_IODEL'; call gtv(trim(nm),tksw(prgn,nm),ltb216,
C$$$     .  def_lg=F,note='Attempt to read increments from disk')

C$$$      sw = tksw(prgn,'TB_FIT')
C$$$      if (sw /= 2) then
C$$$        if (io_help .eq. 1) call info0(1,1,0,
C$$$     .  ' * The FIT tokens are used for tight-binding fitting')
C$$$      nm='TB_FIT_MODE'; call gtv(trim(nm),tksw(prgn,nm),tb_fmode,
C$$$     .  note='Fit mode: 0 fit params, 1 fit in range of RFIT')
C$$$      if (tb_fmode>=0) call isanrg(tb_fmode,0,1,'rdctrl','FIT_MODE',T)
C$$$      if (tb_fmode>=0 .or. io_help .eq. 1) then
C$$$      nm='TB_FIT_NBFIT'; call gtv(trim(nm),tksw(prgn,nm),tb_nbfit,
C$$$     .    nmin=2,note='Min and max band number to fit')
C$$$      nm='TB_FIT_EBFIT'; call gtv(trim(nm),tksw(prgn,nm),tb_ebfit,
C$$$     .  nmin=2,note='Energy range to fit eigenvalues')
C$$$      if (io_help .eq. 1 .or. tb_fmode .eq. 1) then
C$$$      nm='TB_FIT_RFIT'; call gtv(trim(nm),tksw(prgn,nm),tb_rmfit,
C$$$     .  nmin=2,note=
C$$$     .    '(MODE=1) range for fitting Slater-Koster parameters')
C$$$      endif
C$$$      nm='TB_FIT_LAM'; call gtv(trim(nm),tksw(prgn,nm),tb_alam,
C$$$     .  def_r8=0.001d0,note=
C$$$     .  'Initial value of Levenberg-Marquardt lambda variable')
C$$$      nm='TB_FIT_SCL'; call gtv(trim(nm),tksw(prgn,nm),tb_alsc,
C$$$     .  def_r8=10d0,note=
C$$$     .  'Scale factor for Levenberg-Marquardt lambda variable')
C$$$      endif
C$$$      endif

C$$$      endif                     !TB parameters

C --- Map (do not maintain for now) ---

C$$$C --- Optics ---
C$$$      if (tksw(prgn,'OPTICS') == 2) goto 279
C$$$      if (io_show+io_help/=0)
C$$$     .  call info0(1,1,0,' --- Parameters for optics ---')

C$$$      nm='OPTICS_MODE'; call gtv(trim(nm),tksw(prgn,nm),ctrl_loptc,
C$$$     .  def_i4=0,note=
C$$$     .  '0: do nothing;'//
C$$$     .   '%N%2f-1: generate joint DOS'//
C$$$     .   '%N%3f1: generate linear eps_2'//
C$$$     .   '%N%3f2: generate linear eps_2 with on-the-fly sampling'//
C$$$     .   '%N%3f2: Add 10 for SHG')
C$$$      if (ctrl_loptc == 0) goto 279 ! doesn't skip help mode: value is NULLI
C$$$      nm='OPTICS_NPTS'; call gtv(trim(nm),tksw(prgn,nm),optic_ne,
C$$$     .  def_i4=501,note='Number of energy points')
C$$$      nm='OPTICS_WINDOW'; call gtv(trim(nm),tksw(prgn,nm),optic_window,
C$$$     . def_r8v=(/0d0,1d0/),note='Energy window over which to calc. eps')
C$$$      nm='OPTICS_FILBND'; call gtv(trim(nm),tksw(prgn,nm),optic_ocrng,
C$$$     .  def_i4v=izerv,note='Occ. energy bands from which to calc. eps')
C$$$      nm='OPTICS_EMPBND'; call gtv(trim(nm),tksw(prgn,nm),optic_unrng,
C$$$     .  def_i4v=izerv,note='Unocc energy bands from which to calc. eps')
C$$$      nm='OPTICS_PART'; call gtv(trim(nm),tksw(prgn,nm),optic_mode1,
C$$$     .  def_lg=F,note='Band-to-band decomposition of epsilon')

C$$$      optic_axes = NULLI
C$$$      nm='OPTICS_CHI2'; call gtv(trim(nm),tksw(prgn,nm),nono,
C$$$     .  Texist=ltmp,note='Parameters second harmonic generation.'//
C$$$     .    '  If present then:')
C$$$      if (ltmp) then   ! Only read if OPTICS_CHI2 present or help
C$$$        nm='OPTICS_CHI2_NCHI2'; call gtv(trim(nm),tksw(prgn,nm),
C$$$     .    optic_nchi2,def_i4=0,note=
C$$$     .    'number of direction vectors for which to calc. chi2')
C$$$        nm = 'OPTICS_CHI2_AXES'
C$$$        nn = 3*optic_nchi2
C$$$        if (io_help == 1 .and. sw /= 2) nn = NULLI
C$$$        call gtv(trim(nm),tksw(prgn,nm),optic_axes(1:nn),nmin=nn,
C$$$     .    note='a,b,c direction vectors for each of the nchi2 sets')
C$$$      endif
C$$$      nm='OPTICS_ESCISS'; call gtv(trim(nm),tksw(prgn,nm),optic_esciss,
C$$$     .  def_r8=0d0,note=
C$$$     .  'Scissors operator (energy added to unoccupied bands)')
C$$$  279 continue                  ! end of optics

C --- Dynamics ---
      if (io_show+io_help/=0 .and. tksw(prgn,'DYN')/=2) call info0(1,1,
     .0,' --- Parameters for dynamics and statics ---')

      nm='DYN_NIT'; call gtv(trim(nm),tksw(prgn,nm),nitmv,def_i4=1,
     .note='maximum number of relaxation steps (statics)'//
     .' or time steps (dynamics)')

      nm='DYN_SSTAT'; call gtv(trim(nm),tksw(prgn,nm),nono,Texist=ltmp,
     .note='Parameters for spin statics:')
      if (io_help .ne. 0 .or. ltmp) then
        nm='DYN_SSTAT_MODE'; call gtv(trim(nm),tksw(prgn,nm),sdmod,
     .  def_i4=-1,note=
     .  '%b<0:  No spin statics or dynamics'//
     .  '%N%3f0:  output Euler angles generated from density matrix'//
     .  '%N%3f1:  relax along force'//
     .  '%N%3fAdd 10   to mix angles independently of P,Q'//
     .  '%N%3fAdd 1000 to suppress self-const in P,Q')
        nm='DYN_SSTAT_SCALE'; call gtv(trim(nm),tksw(prgn,nm),sdprm(1),
     .  def_r8=1d0,note='Scale factor amplifying magnetic forces')
        if (io_help .ne. 0 .or. mod(sdmod,10) .eq. 2) then
          nm='DYN_SSTAT_TAU'; call gtv(trim(nm),tksw(prgn,nm),sdprm(2),
     .    note=
     .    'Time step (only read if Landau-Gilbert dynamics specified)')
        endif
        nm='DYN_SSTAT_MAXT'; call gtv(trim(nm),tksw(prgn,nm),sdprm(4),
     .  def_r8=0d0,note='Maximum allowed change in angle')
        nm='DYN_SSTAT_ETOL'; call gtv(trim(nm),tksw(prgn,nm),sdprm(5),
     .  def_r8=0d0,note='Set tau=0 this iter if etot-ehf>ETOL')
      endif
c      lncol16 = sdmod .ge. 0

      nm='DYN_MSTAT'; call gtv(trim(nm),tksw(prgn,nm),nono,Texist=ltmp,
     .note='Parameters for molecular statics')
      mdprm(1) = 0
      if (io_help .ne. 0 .or. ltmp) then
        nm='DYN_MSTAT_MODE'; call gtv(trim(nm),tksw(prgn,nm),i0,
     .  def_i4=0,note=
     .  '0: no relaxation  '//
     .  '4: conjugate gradients  '//
     .  '5: Fletcher-Powell  '//
     .  '6: Broyden')
        mdprm(1) = i0
        nm='DYN_MSTAT_HESS'; call gtv(trim(nm),tksw(prgn,nm),ltmp,
     .  def_lg=T,note='Read hessian matrix')
        mdprm(2) = isw(ltmp)
        nm='DYN_MSTAT_XTOL'; call gtv(trim(nm),tksw(prgn,nm),mdprm(3),
     .  def_r8=1d-3,note=
     .  'Convergence criterion in displacements'//
     .  '%N%3fXTOL>0: use length;  <0: use max val;  =0: do not use')
        nm='DYN_MSTAT_GTOL'; call gtv(trim(nm),tksw(prgn,nm),mdprm(4),
     .  def_r8=0d0,note=
     .  'Convergence criterion in gradients'//
     .  '%N%3fGTOL>0: use length;  <0: use max val;  =0: do not use')
        nm='DYN_MSTAT_STEP'; call gtv(trim(nm),tksw(prgn,nm),mdprm(5),
     .  def_r8=0.015d0,note='Initial (and maximum) step length')
        nm='DYN_MSTAT_NKILL'; call gtv(trim(nm),tksw(prgn,nm),i0,
     .  def_i4=0,note='Remove hessian after NKILL iter')
        mdprm(6) = i0
        nm='DYN_MSTAT_PDEF'; call gtv(trim(nm),tksw(prgn,nm),lat_defm,
     .  def_r8v=zerov,note='Lattice deformation modes')
      endif

C --- if DYN_MSTAT_MODE < 4 then attempt to read MD ---
      if (mdprm(1) .lt. 4) then
        nm='DYN_MD'; sw=tksw(prgn,nm)
        call gtv(trim(nm),sw,nono,Texist=ltmp,
     .  note='Parameters for molecular dynamics')
        if (io_help .ne. 0 .or. ltmp) then
          nm='DYN_MD_MODE'; call gtv(trim(nm),tksw(prgn,nm),i0,
     .    def_i4=0,note=
     .    '0: no MD '//
     .    '1: NVE '//
     .    '2: NVT '//
     .    '3: NPT')
          mdprm(1) = i0
          nm='DYN_MD_TSTEP'; call gtv(trim(nm),tksw(prgn,nm),mdprm(2),
     .    def_r8=fs,note='Time step (a.u.)')
          nm='DYN_MD_TEMP'; call gtv(trim(nm),tksw(prgn,nm),mdprm(3),
     .    def_r8=300d0*degK,note='Temperature (a.u.)')
          nm='DYN_MD_TAUP'; call gtv(trim(nm),tksw(prgn,nm),mdprm(4),
     .    def_r8=10d0*fs,note='Thermostat relaxation time (a.u.)')
          nm='DYN_MD_TIME'; call gtv(trim(nm),tksw(prgn,nm),mdprm(5),
     .    def_r8=1d6*fs,note='Total MD time (a.u.)')
          nm='DYN_MD_TAUB'; call gtv(trim(nm),tksw(prgn,nm),mdprm(6),
     .    def_r8=100d0*fs,note='Barostat relaxation time (a.u.)')
          if (tksw(prgn,'DYN_MD') .ne. 2) call info2(10,0,0,
     .    ' NB: 1 deg.K = %e a.u.; 1 fs = %d a.u.',degK,fs)
        endif
        if (sw .ne. 2 .and. mdprm(1) .gt. 3)
     .  call rx('M_RDCTRL: illegal "MODE" parameter in token "DYN_MD"')
      endif

C Content of mdprm(i)
C   i = 1: statics/dynamics switch
C          0 do neither
C          1-3 MD (1: NVE, 2: NVT; 3: NPT)
C          >3 statics (4: conjugate gradients, 5: Fletcher-Powell, 6: Broyden)
C          >100 statics with lattice deformation
C   i = 2: (stat) 1 read hessian matrix
C          (dyn)  MD timestep
C   i = 3: (stat) relaxation x-tolerance
C          (dyn)  temperature
C   i = 4: (stat) relaxation g-tolerance
C          (dyn)  thermostat relaxation time
C   i = 5: (stat) step length
C          (dyn)  total MD time
C   i = 6: (stat) Remove hessian after this many steps
C          (dyn)  barostat relaxation time

      nm='DYN_SDYN'; call gtv(trim(nm),tksw(prgn,nm),nono,Texist=lsdyn,
     .note='Parameters for spin dynamics')
C     if (io_help .ne. 0 .or. ltmp) then
      nm='DYN_SDYN_KT'; call gtv(trim(nm),tksw(prgn,nm),move_kt,
     .note='Temperature, in a.u.')
      nm='DYN_SDYN_TS'; call gtv(trim(nm),tksw(prgn,nm),move_ts,
     .note='Time step in a.u.')
      nm='DYN_SDYN_TEQU'; call gtv(trim(nm),tksw(prgn,nm),move_tsequ,
     .def_r8=0d0,note='Equilibration time, a.u.')
      nm='DYN_SDYN_TTOT'; call gtv(trim(nm),tksw(prgn,nm),move_tstot,
     .note='Duration of total simulation, in a.u.')
C     endif

      nm='DYN_MMHAM'; call gtv(trim(nm),tksw(prgn,nm),mmham,
     .note='Rules for micromagnetics hamiltonian')

C ... Global daemons
      nm='DYN_GD'; call gtv(trim(nm),tksw(prgn,nm),nono,Texist=ltmp,
     .note='Parameters for Global daemons in theromstats')
C    if (io_help .ne. 0 .or. ltmp) then
      nm='DYN_GD_NTHERM'; call gtv(trim(nm),tksw(prgn,nm),gd_nmodt,
     .def_i4=1,note='Number of GD thermostats')
      nn = gd_nmodt
      if (io_help .ne. 0 .or. tksw(prgn,nm) == 2) then
        nn = NULLI
      else
        call isanrg(nn,1,3,'rdctrl','GD_NTHERM',T)
      endif
      nm='DYN_GD_MODET'; call gtv(trim(nm),tksw(prgn,nm),
     .gd_modt(1:max(1,nn)),
     .nmin=nn,def_i4v=(/31,32,33/),note='Thermostat mode(s)')
      nm='DYN_GD_CT'; call gtv(trim(nm),tksw(prgn,nm),
     .gd_ct(1:max(1,nn)),
     .nmin=nn,def_r8v=(/1d0,1d0,1d0/),note=
     .'Corresponding thermostat coefficient(s)')
C     endif

C ... Bulirsch-Stoer integration parmaters
      nm='DYN_BSINT'; call gtv(trim(nm),tksw(prgn,nm),nono,
     .Texist=lbsprm,note='Parameters for Bulirsch-Stoer integration')
      if (io_help .ne. 0 .or. lbsprm) then
        nm='DYN_BSINT_NEW'; call gtv(trim(nm),tksw(prgn,nm),prmint_new,
     .  def_lg=T,note='Start new SD run, Bulirsch-Stoer integration')
        nm='DYN_BSINT_TOL'; call gtv(trim(nm),tksw(prgn,nm),prmint_tol,
     .  note='Tolerance in numerical integration')
        nm='DYN_BSINT_TS0'; call gtv(trim(nm),tksw(prgn,nm),prmint_ts0,
     .  def_r8=0d0,note='Minimum time step in units of TS')
        nm='DYN_BSINT_MX'; call gtv(trim(nm),tksw(prgn,nm),prmint_mx,
     .  def_i4=7,note='Maximum order of rational function extrapolation')
        nm='DYN_BSINT_MI'; call gtv(trim(nm),tksw(prgn,nm),prmint_mi,
     .  def_i4=11,note='Maximum number of midpoint rules to invoke')
        nn = prmint_mi
        if (io_help .ne. 0 .or. tksw(prgn,nm) == 2) then
          nn = NULLI
        else
          call isanrg(nn,2,11,'rdctrl','BSINT_MI',T)
        endif
        if (nout .gt. 0) nn = nout
        nm='DYN_BSINT_NSEQ'; call gtv(trim(nm),tksw(prgn,nm),
     .  prmint_nseq(1:max(1,nn)),
     .  nmin=nn,def_i4v=(/2,4,6,8,12,16,24,32,48,64,96/),
     .  note='Sequence of number of midpoint divisions')
      endif

cxxx int
C      sdyn_sdmod= sdyn(2)       !xxxxx
C      sdyn_fscal= sdyn(3)
C      sdyn_tau = sdyn(4)
C      sdyn_NoseEtot0= sdyn(5)
C      sdyn_maxtheta= sdyn(6)
C      sdyn_etol= sdyn(7)

C$$$C --- GW ---
C$$$      if (io_show+io_help/=0 .and. tksw(prgn,'GW')/=2) call info0(1,1,
C$$$     .  0,' --- Parameters for GW input (job=-1) ---')
C$$$      nm='GW_MKSIG'; call gtv(trim(nm),tksw(prgn,nm),gw_mksig,
C$$$     .  def_i4=0,note=
C$$$     .  'Make full sigma(k) for self-consistent calculations:'//
C$$$     .  '%N%3f0 do not create self-energy matrix sigma'//
C$$$     .  '%N%3f1 Make sigma_nn''(E_f)'//
C$$$     .  '%N%3f2 Make sigma_nn''(E_n/2+E_n''/2)'//
C$$$     .  '%N%3f3 Make sigma_nn''(E_n)/2 + sigma_nn''(E_n'')/2')

C$$$      nm='GW_NKABC'; call gtv(trim(nm),tksw(prgn,nm),gw_nabc,nout=nout,
C$$$     .  note='No. qp along each of 3 lattice vectors.'//
C$$$     .  '%N%3fSupply one number for all vectors or a separate '//
C$$$     .  'number for each vector.')
C$$$      call fill3in(nout,gw_nabc)
C$$$      if (gw_nabc(1) == NULLI .or. gw_nabc(1) == 0) gw_nabc = bz_nabc

C$$$      nm='GW_GCUTB'; call gtv(trim(nm),tksw(prgn,nm),gw_gcutb,
C$$$     .  def_r8=2.7d0,note='G-vector cutoff for basis functions')
C$$$      nm='GW_GCUTX'; call gtv(trim(nm),tksw(prgn,nm),gw_gcutx,
C$$$     .  def_r8=2.2d0,note='G-vector cutoff for response function')
C$$$      nm='GW_ECUTS'; call gtv(trim(nm),tksw(prgn,nm),gw_ecuts,
C$$$     .  def_r8=2.5d0,note=
C$$$     .  '(self-consistency) max energy for which to calc. sigma')
C$$$      nm='GW_NBAND'; call gtv(trim(nm),tksw(prgn,nm),gw_nband,
C$$$     .  def_i4=-NULLI,note=
C$$$     .  'No. bands to include in calc. response function')
C$$$      nm='GW_NIME'; call gtv(trim(nm),tksw(prgn,nm),gw_nime,
C$$$     .  def_i4=6,note='Number of imaginary energy points')
C$$$      nm='GW_DELRE'; call gtv(trim(nm),tksw(prgn,nm),gw_delre,
C$$$     .  def_r8v=(/.01d0,.04d0/),note=
C$$$     .  'Energy mesh spacing on real axis.'//
C$$$     .  '%N%3fSecond term is coff for nonlinear part')
C$$$      nm='GW_DELTA'; call gtv(trim(nm),tksw(prgn,nm),gw_deltax,
C$$$     .  def_r8=-1d-4,note=
C$$$     .  'Imaginary energy added to denominator in generation of Pi.'//
C$$$     .  '%N%3fNot used in Faleev''s fast integration.'//
C$$$     .  '%N%3fDelta should be negative.')
C$$$      nm='GW_DELTAW'; call gtv(trim(nm),tksw(prgn,nm),gw_deltaw,
C$$$     .  def_r8=0.02d0,note=
C$$$     .  'Width for finite difference in diff. of sigma for Z factor')
C$$$      nm='GW_GSMEAR'; call gtv(trim(nm),tksw(prgn,nm),gw_gsmear,
C$$$     .  def_r8=3d-3,note=
C$$$     .  'Broadening in pole of GF for computation of sigma')
C$$$      nm='GW_PBTOL'; call gtv(trim(nm),tksw(prgn,nm),gw_pbtol,
C$$$     .  def_r8=1d-3,note=
C$$$     .  'retain prod-basis with overlap eigenvalues>PBTOL')
C$$$      nm='GW_QOFFP'; call gtv(trim(nm),tksw(prgn,nm),gw_qoffp,
C$$$     .  def_r8=1d0,note=
C$$$     .  'k-point offset parameter for special BZ integration for xi')

cC --- Plane ---
c      if (io_show+io_help/=0 .and. tksw(prgn,'PLANE')/=2)
c     .  call info0(1,1,0,' --- Parameters defining a plane ---')
c      nm='PLANE_NORMAL'; call gtv(trim(nm),tksw(prgn,nm),
c     .  slat_plat2(7:9),def_r8v=zerov,note='Plane normal')
c      nm='PLANE_X'; call gtv(trim(nm),tksw(prgn,nm),
c     .  slat_plat2(1:3),def_r8v=zerov,note='Plane X-axis')
c      nm='PLANE_NORMAL'; call gtv(trim(nm),tksw(prgn,nm),
c     .  slat_plat2(4:6),def_r8v=zerov,note='Plane Y-axis')

      end subroutine readctrl

      subroutine readctrlpq(prgn,nclasp,nl,nsp,pnu,qnu,pp,z,ves,
     .initc,ics,clabl)
C- Read class data from ctrl file
C ----------------------------------------------------------------------
Ci Inputs
Ci   prgn
Ci   nclasp
Ci   clabl :class name
Co Outputs
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   17 Oct 07
C ----------------------------------------------------------------------
      use m_toksw
      use m_gtv
C     implicit none
C ... Passed parameters
      character*(*):: prgn
      integer:: nclasp,nl,nsp
C ... Class-related arrays
      integer:: initc(nclasp),ics(nclasp)
      real(8):: pnu(nl,nsp,nclasp),qnu(3,nl,nsp,nclasp),
     .pp(6,nl,nsp,nclasp),ves(nclasp),z(nclasp)
      character*8 clabl(nclasp)
C ... Local parameters
      character(128) :: nm
      character(256) :: outs
      logical:: ltmp,rdves,cmdopt
      integer:: stdo,lgunit,i,j,ic,lp1,jj(2),nout,nlaj,nlaji,
     .iv(100),ii,n,js
      real(8):: xv(100),xx,pi
C ... Heap
Cdwdef1       integer w(1)
Cdwdef1       common /w/ w

      stdo = lgunit(1)
      pi = 4d0*datan(1d0)

C --- Start ---
      if (tksw(prgn,'START') == 2) goto 159
C$$$      if (io_show+io_help/=0) call
C$$$     .  info0(1,1,0,' --- Read potential data (ASA) ---')

C$$$      nm='START_BEGMOM'; call gtv(trim(nm),tksw(prgn,nm),lasa3,
C$$$     .  def_i4=1,note=
C$$$     .  '0 Start from potential parameters (on file)'//
C$$$     .  '%N%3f1 Create potential and pot pars from P,Q'//
C$$$     .  '%N%3f2 Create pot pars from file potential'//
C$$$     .  '%N%3f3 Create potential from restart file')

C$$$      nm='START_FREE'; call gtv(trim(nm),tksw(prgn,nm),lasa8,
C$$$     .  def_lg=F,note=
C$$$     .  '(Free atom mode): make self-consistent density for free atoms')

      nm='START_CNTROL'; call gtv(trim(nm),tksw(prgn,nm),ltmp,def_lg=F,
     .note='if F, the remainder of this category is ignored')
      if (io_help==0 .and. .not. ltmp) goto 159

      nm='START_RDVES'; call gtv(trim(nm),tksw(prgn,nm),rdves,
     .def_lg=F,note='Read Ves(RMT) from this category')
C     Subtract current ves from pp now; add updated ves later
      if (rdves) call shftpp(nclasp,nl*nsp,pp,ves,ves,T,F)

C ... Read data for each occurence of ATOM=
C     j = index to jth occurence; ic = class index
      outs = ' '
      j = 1
      do  while (j .gt. 0)
        if (io_help /= 0) then
          write(stdo,"(/' For each class, read the following:')")
        elseif (io_show>0) then
          write(stdo,"(1x)")
C         call info(1,0,0,' ... Site %i',j,0)
        endif

        jj=(/1,j/)
        nm='START_ATOM'; call gtv(trim(nm),tksw(prgn,nm),alabl,nmin=10,
     .  cindx=jj,nout=nout,note='Class label')
        if (nout == 0 .and. io_help == 0) then
          j = 0
          cycle
        endif

        if (io_help == 0) then
          ic = 0
          do  i = 1, nclasp
            if (alabl == clabl(i)) then
              ic = i
              call awrit0('%a'//alabl//'%a,',outs,len(outs),0)
              goto 881
            endif
          enddo
          call info0(10,0,0,' rdctrl: (warning) class "'
     .    //alabl//'%a" not in class table')
          j = j+1; cycle
  881     continue
          nlaj = 1+lmxa(ic)
          nlaji = nlaj
        else                    !help mode; cycle only once
          nlaj = 1
          nlaji = NULLI
          ic = 1
          j = -1
        endif

        nm='START_ATOM_P'; call gtv(trim(nm),tksw(prgn,nm),
     .  xv(1:nlaj*nsp),nmin=nlaji*nsp,nout=nout,cindx=jj,note=
     .  'Methfessel''s Potential functions')
C       Copy to pnu, spins separately in case nlja<nl
        if (nout .gt. 0) then
          initc(ic) = 2*(initc(ic)/2) + 1
          call dpzero(pnu(1,1,ic),nl*nsp)
          call dcopy(nlaj,xv,1,pnu(1,1,ic),1)
          if (nsp .eq. 2)
     .    call dcopy(nlaj,xv(1+nlaj),1,pnu(1,2,ic),1)
C         Fill in values for l>lmxa
          if (.not. cmdopt('--keepnu',8,0,nm)) then
            call config(pnu(1,1,ic),-1,z(ic),iv,n)
            do  i = 1, nsp
              do  lp1 = 1, nl
                ii = iv(lp1)
                xx = pnu(lp1,i,ic)
C           If P=0, use free-electron value
                if (xx .eq. 0) then
                  xx = ii + 0.5d0 - datan(dble(lp1-1))/pi
C           If fractional part only, add integer part
                elseif (xx .lt. 1) then
                  xx = xx + ii
                endif
                pnu(lp1,i,ic) = xx
              enddo
            enddo
          endif
        endif

        nm='START_ATOM_Q'; call gtv(trim(nm),tksw(prgn,nm),
     .  xv(1:3*nlaj*nsp),nmin=3*nlaji*nsp,nout=nout,cindx=jj,
     .  note='Energy moments of the spherical density')
C       Copy to qnu, spins separately in case nlja<nl
        if (nout .gt. 0) then
          call dpzero(qnu(1,1,1,ic),3*nl*nsp)
          call dcopy(3*nlaj,xv,1,qnu(1,1,1,ic),1)
          if (nsp .eq. 2)
     .    call dcopy(3*nlaj,xv(1+3*nlaj),1,qnu(1,1,2,ic),1)
        endif

        nm='START_ATOM_ENU'; call gtv(trim(nm),tksw(prgn,nm),
     .  xv(1:nlaj*nsp),nmin=nlaji*nsp,nout=nout,cindx=jj,note=
     .  'Linearization energies')
        if (nout .gt. 0) then
          call dcopy(nlaj,xv,1,pp(1,1,1,ic),6)
          if (nsp .eq. 2)
     .    call dcopy(nlaj,xv(1+nlaj),1,pp(1,1,2,ic),6)
        endif

        if (rdves) then
          nm='START_ATOM_V'; call gtv(trim(nm),tksw(prgn,nm),
     .    xx,nmin=1,nout=nout,cindx=jj,note=
     .    'ASA potential at MT boundary')
          if (nout .gt. 0) ves(ic) = xx
          nm='START_ATOM_DV'; call gtv(trim(nm),tksw(prgn,nm),
     .    xx,nmin=1,nout=nout,cindx=jj,note=
     .    'Add to ASA potential at MT boundary')
          if (nout .gt. 0) ves(ic) = ves(ic) + xx
        endif

        j = j+1
      enddo
      if (rdves) call shftpp(nclasp,nl*nsp,pp,ves,ves,F,T)

C ... Printout
      if (io_show /= 0 .and. io_help == 0) then
        call awrit0('%23o%0p%N rdctrl: read P,Q for %a%b%N P,Q,ppars:',
     .  outs,80,-stdo)
        do   j = 1, nclasp
          js = ics(j)
          nlaj = 1+lmxa(j)
          if (mod(initc(j),2) .eq. 1) then
            print 335, 'class ', clabl(j), (pnu(lp1,1,j),lp1=1,nlaj)
            if (nsp .eq. 2)
     .      print 335,' spin ', '2       ',(pnu(lp1,2,j),lp1=1,nlaj)
            print 335,'      ', 'Q       ',(qnu(1,lp1,1,j),lp1=1,nlaj)
            if (nsp .eq. 2)
     .      print 335,' spin ', '2       ',(qnu(1,lp1,2,j),lp1=1,nlaj)
  335       format(4x,a,a,7f8.4)
          else
            print 335, 'class ', clabl(j) // '  missing pnu'
          endif
          if (mod(initc(j)/2,2) .eq. 1) then
            print 336,'  enu       ',      (pp(1,lp1,1,j), lp1=1,nlaj)
            if (nsp .eq. 2)
     .      print 335, ' spin ', '2       ',(pp(1,lp1,2,j),lp1=1,nlaj)
            print 336,'    c       ',      (pp(2,lp1,1,j), lp1=1,nlaj)
            if (nsp .eq. 2)
     .      print 335, ' spin ', '2       ',(pp(2,lp1,2,j),lp1=1,nlaj)
            print 336,'srdel       ',      (pp(3,lp1,1,j), lp1=1,nlaj)
            if (nsp .eq. 2)
     .      print 335, ' spin ', '2       ',(pp(3,lp1,2,j),lp1=1,nlaj)
  336       format(6x,a,7f8.4)
          else
            print 335, 'class ', clabl(j) // '  missing ppars'
          endif
          print *, ' '
        enddo
      endif

  159 continue                  ! End of START

      end subroutine readctrlpq

      end module


      subroutine fill3in(nin,res)
C- fills res(2) or res(2:3) if res(1) or res(1:2) are given
C     implicit none
      integer nin,res(3)
      if (nin==2) then
        res(3) = res(2)
      elseif (nin==1) then
        res(2:3) = res(1)
      endif
      end


      subroutine toksw_init(debug)
C- Declares program-specific tokens
      use m_toksw
C     implicit none
      character*(200):: io_sw, struc_sw
      character*(50) sscale
      character*(200):: dyn_mstat, dyn_md
      logical:: debug
C new design
c      call tadd(" IO_SHOW ",req,"LM LMF LMFGWD")
c      call tadd(" IO_SHOW ",iopt,"LM LMF LMFGWD")
C - Store data into arrays above -----
      call clear_swtok(debug)
      io_sw = ' HEADER~ CONST~ CMD~ IO_SHOW~ IO_HELP~ IO_VERBOS~ '//
     .'IO_WKP~ IO_IACTIV~ IO_TIM~ OPTIONS_Q~'
      struc_sw =
     .' STRUC STRUC_NSPEC~ STRUC_NCLASS STRUC_FILE~ STRUC_NBAS '//
     .'STRUC_PLAT STRUC_ALAT STRUC_DALAT~ STRUC_NL~ STRUC_SHEAR~ '//
     .'STRUC_ROT~ STRUC_DEFGRD~ STRUC_STRAIN~ STRUC_ALPHA~ '//
     .'STRUC_LOADV~ '
      sscale = ' SPEC_SCLWSR~ SPEC_WSRMAX~ SPEC_OMAX1~ SPEC_OMAX2~'
      dyn_mstat = ' DYN_MSTAT~ DYN_MSTAT_MODE DYN_MSTAT_HESS~ '//
     .'DYN_MSTAT_XTOL~ DYN_MSTAT_GTOL~ DYN_MSTAT_STEP~ '//
     .'DYN_MSTAT_NKILL~ DYN_MSTAT_PDEF~ '
      dyn_md = ' DYN_MD~ DYN_MD_MODE DYN_MD_TSTEP~ DYN_MD_TEMP~ '//
     .'DYN_MD_TAUP~ DYN_MD_TIME DYN_MD_TAUB~ '
C --- LMFA switches ---
      call nswadd()
      call tkadd(" LMFA::" )
      call tkadd(" VERS_LM")
C     call tkadd(" VERS_FP")
      call tkadd(trim(io_sw))
      call tkadd(trim(struc_sw))
c      call tkadd(" STRUC_NSPEC")
      call tkadd(" OPTIONS~")
      call tkadd(" OPTIONS_FRZ~")
      call tkadd(" OPTIONS_HF~")
      call tkadd(" SPEC")
      call tkadd(sscale)
C     call tkadd(" SPEC_SCLWSR~")
C     call tkadd(" SPEC_WSRMAX~")
C     call tkadd(" SPEC_OMAX1~")
C     call tkadd(" SPEC_OMAX2~")
      call tkadd(" SPEC_ATOM")
      call tkadd(" SPEC_ATOM_A~")
      call tkadd(" SPEC_ATOM_AMASS~")
      call tkadd(" SPEC_ATOM_C-HOLE~")
      call tkadd(" SPEC_ATOM_C-HQ~")
      call tkadd(" SPEC_ATOM_CSTRMX~")
      call tkadd(" SPEC_ATOM_EH~")
      call tkadd(" SPEC_ATOM_EH2~")
      call tkadd(" SPEC_ATOM_EREF~")
      call tkadd(" SPEC_ATOM_FRZWF~")
      call tkadd(" SPEC_ATOM_IDMOD~")
      call tkadd(" SPEC_ATOM_IDU~")
      call tkadd(" SPEC_ATOM_IDXDN~")
      call tkadd(" SPEC_ATOM_ISOTOPE")
      call tkadd(" SPEC_ATOM_JH~")
      call tkadd(" SPEC_ATOM_KMXA~")
      call tkadd(" SPEC_ATOM_LFOCA~")
      call tkadd(" SPEC_ATOM_LMX~")
      call tkadd(" SPEC_ATOM_LMXA~")
      call tkadd(" SPEC_ATOM_LMXL~")
      call tkadd(" SPEC_ATOM_MMOM~")
      call tkadd(" SPEC_ATOM_NR~")
      call tkadd(" SPEC_ATOM_P~")
      call tkadd(" SPEC_ATOM_PZ~")
      call tkadd(" SPEC_ATOM_Q~")
      call tkadd(" SPEC_ATOM_R")
      call tkadd(" SPEC_ATOM_R/A")
      call tkadd(" SPEC_ATOM_RCFA~")
      call tkadd(" SPEC_ATOM_RFOCA~")
      call tkadd(" SPEC_ATOM_RS3~")
      call tkadd(" SPEC_ATOM_RSMA~")
      call tkadd(" SPEC_ATOM_RSMFA~")
      call tkadd(" SPEC_ATOM_RSMG~")
      call tkadd(" SPEC_ATOM_RSMH~")
      call tkadd(" SPEC_ATOM_RSMH2~")
      call tkadd(" SPEC_ATOM_R/W")
      call tkadd(" SPEC_ATOM_UH~")
      call tkadd(" SPEC_ATOM_Z")
      call tkadd(" SITE")
      call tkadd(" SITE_FILE~")
      call tkadd(" SITE_MODE~")
      call tkadd(" SITE_ATOM")
      call tkadd(" SITE_ATOM_DPOS~")
C     call tkadd(" SITE_ATOM_PL~")
      call tkadd(" SITE_ATOM_POS")
      call tkadd(" SITE_ATOM_XPOS")
      call tkadd(" SITE_ATOM_RELAX~")
      call tkadd(" SYMGRP~")
      call tkadd(" EWALD~")
      call tkadd(" EWALD_AS~")
      call tkadd(" EWALD_NKDMX~")
      call tkadd(" EWALD_RPAD~")
      call tkadd(" EWALD_TOL~")
      call tkadd(" ITER ITER_NRMIX~")
C      call tkadd(" ITER~")
C      call tkadd(" ITER_AMIX~")
C      call tkadd(" ITER_CONV~")
C      call tkadd(" ITER_CONVC~")
C      call tkadd(" ITER_MIX~")
C      call tkadd(" ITER_XIPMX~")
C     call tkadd(" START~")
      call tkadd(" HAM~")
      call tkadd(" HAM_NSPIN~")
      call tkadd(" HAM_REL~ HAM_SO~")
      call tkadd(" HAM_GGA~")
      call tkadd(" HAM_XCFUN~")
C      call tkadd(" DYN~")
C      call tkadd(" DYN_MD~")
C      call tkadd(" DYN_MSTAT~")
C      call tkadd(" DYN_NIT~")
C      call tkadd(" DYN_PDEF~")

c This call tkadd(" SPEC_ATOM_MIX~") was missing,
c then the initialization of default values for mxcst are not supplied
c   --->this results in a strange behevior!!
c takao Sep11-tukuba 2009
      call tkadd(" SPEC_ATOM_MIX~")

C --- LMFGWD switches ---
      call nswadd()
      call tkadd(" LMFGWD::" )
      call tkadd(" VERS_LM")
      call tkadd(" VERS_FP")
      call tkadd(trim(io_sw))
      call tkadd(trim(struc_sw))
      call tkadd(" OPTIONS~")
      call tkadd(" OPTIONS_FRZ~")
      call tkadd(" OPTIONS_HF~")
      call tkadd(" HAM_NSPIN~ HAM_REL~ HAM_SO~")
      call tkadd(" SPEC")
      call tkadd(sscale)
C     call tkadd(" SPEC_SCLWSR~")
C     call tkadd(" SPEC_WSRMAX~")
C     call tkadd(" SPEC_OMAX1~")
C     call tkadd(" SPEC_OMAX2~")
      call tkadd(" SPEC_ATOM")
      call tkadd(" SPEC_ATOM_A~")
      call tkadd(" SPEC_ATOM_AMASS~")
      call tkadd(" SPEC_ATOM_C-HOLE~")
      call tkadd(" SPEC_ATOM_C-HQ~")
      call tkadd(" SPEC_ATOM_CSTRMX~")
      call tkadd(" SPEC_ATOM_EH")
      call tkadd(" SPEC_ATOM_EH2~")
      call tkadd(" SPEC_ATOM_EREF~")
      call tkadd(" SPEC_ATOM_FRZWF~")
      call tkadd(" SPEC_ATOM_IDMOD~")
      call tkadd(" SPEC_ATOM_IDU~")
      call tkadd(" SPEC_ATOM_IDXDN~")
      call tkadd(" SPEC_ATOM_ISOTOPE")
      call tkadd(" SPEC_ATOM_JH~")
      call tkadd(" SPEC_ATOM_KMXA~")
      call tkadd(" SPEC_ATOM_LFOCA~")
      call tkadd(" SPEC_ATOM_LMX~")
      call tkadd(" SPEC_ATOM_LMXA~")
      call tkadd(" SPEC_ATOM_LMXL~")
      call tkadd(" SPEC_ATOM_LMXPB~")
      call tkadd(" SPEC_ATOM_MMOM~")
      call tkadd(" SPEC_ATOM_NR~")
      call tkadd(" SPEC_ATOM_P~")
      call tkadd(" SPEC_ATOM_PBAS~")
      call tkadd(" SPEC_ATOM_PBAS2~")
      call tkadd(" SPEC_ATOM_PZ~")
      call tkadd(" SPEC_ATOM_Q~")
      call tkadd(" SPEC_ATOM_R")
      call tkadd(" SPEC_ATOM_R/A")
      call tkadd(" SPEC_ATOM_RFOCA~")
      call tkadd(" SPEC_ATOM_RS3~")
      call tkadd(" SPEC_ATOM_RSMA~")
      call tkadd(" SPEC_ATOM_RSMFA~")
      call tkadd(" SPEC_ATOM_RSMG~")
      call tkadd(" SPEC_ATOM_RSMH")
      call tkadd(" SPEC_ATOM_RSMH2~")
      call tkadd(" SPEC_ATOM_R/W")
      call tkadd(" SPEC_ATOM_UH~")
      call tkadd(" SPEC_ATOM_Z")
      call tkadd(" SITE")
      call tkadd(" SITE_FILE~")
      call tkadd(" SITE_MODE~")
      call tkadd(" SITE_ATOM")
      call tkadd(" SITE_ATOM_DPOS~")
C     call tkadd(" SITE_ATOM_PL~")
      call tkadd(" SITE_ATOM_POS")
      call tkadd(" SITE_ATOM_XPOS")
      call tkadd(" SITE_ATOM_RELAX~")
      call tkadd(" SYMGRP~")
      call tkadd(" GW~")
      call tkadd(" GW_DELRE~")
      call tkadd(" GW_DELTA~")
      call tkadd(" GW_DELTAW~")
      call tkadd(" GW_ECUTS~")
      call tkadd(" GW_GCUTB~")
      call tkadd(" GW_GCUTX~")
      call tkadd(" GW_GSMEAR~")
      call tkadd(" GW_MKSIG~")
      call tkadd(" GW_NBAND~")
      call tkadd(" GW_NIME~")
      call tkadd(" GW_NKABC~")
      call tkadd(" GW_PBTOL~")
      call tkadd(" GW_QOFFP~")
      call tkadd(" BZ")
      call tkadd(" BZ_BZJOB~")
      call tkadd(" BZ_DELEF~")
      call tkadd(" BZ_DMAT~")
      call tkadd(" BZ_DOS~")
      call tkadd(" BZ_EF0~")
      call tkadd(" BZ_EFMAX~")
      call tkadd(" BZ_FSMOM~")
      call tkadd(" BZ_GETQP~")
      call tkadd(" BZ_INVIT~")
      call tkadd(" BZ_METAL~")
      call tkadd(" BZ_MULL~")
      call tkadd(" BZ_N~")
      call tkadd(" BZ_NEVMX~")
      call tkadd(" BZ_NKABC~")
      call tkadd(" BZ_NOINV~")
      call tkadd(" BZ_NPTS~")
      call tkadd(" BZ_N.W~")
      call tkadd(" BZ_PUTQP~")
      call tkadd(" BZ_SAVDOS~")
      call tkadd(" BZ_TETRA~")
      call tkadd(" BZ_W~")
      call tkadd(" BZ_ZBAK~")
      call tkadd(" BZ_ZVAL~")
      call tkadd(" EWALD~")
      call tkadd(" EWALD_AS~")
      call tkadd(" EWALD_NKDMX~")
      call tkadd(" EWALD_RPAD~")
      call tkadd(" EWALD_TOL~")
C      call tkadd(" ITER~")
C      call tkadd(" ITER_AMIX~")
C      call tkadd(" ITER_CONV~")
C      call tkadd(" ITER_CONVC~")
C      call tkadd(" ITER_MIX~")
C      call tkadd(" ITER_XIPMX~")
C     call tkadd(" START~")
      call tkadd(" HAM")
      call tkadd(" HAM_NSPIN~")
      call tkadd(" HAM_ELIND~")
      call tkadd(" HAM_EWALD~")
      call tkadd(" HAM_FORCES~")
      call tkadd(" HAM_FRZWF~")
      call tkadd(" HAM_FTMESH")
      call tkadd(" HAM_GGA~")
      call tkadd(" HAM_GMAX")
      call tkadd(" HAM_PMAX~")
      call tkadd(" HAM_PMIN~")
      call tkadd(" HAM_PWMODE~ HAM_PWEMIN~ HAM_PWEMAX~ HAM_NPWPAD~")
      call tkadd(" HAM_REL~ HAM_SO~")
      call tkadd(" HAM_RDSIG~")
      call tkadd(" HAM_RSRNGE~")
      call tkadd(" HAM_RSSTOL~")
      call tkadd(" HAM_SIGP~")
      call tkadd(" HAM_TOL~")
      call tkadd(" HAM_XCFUN~")
C      call tkadd(" DYN~")
C      call tkadd(" DYN_MD~")
C      call tkadd(" DYN_MSTAT~")
C      call tkadd(" DYN_NIT~")
C      call tkadd(" DYN_PDEF~")
C --- LMF switches ---
      call nswadd()
      call tkadd(" LMF::" )
      call tkadd(" VERS_LM")
      call tkadd(" VERS_FP")
      call tkadd(trim(io_sw))
      call tkadd(trim(struc_sw))
      call tkadd(" OPTIONS~")
      call tkadd(" OPTIONS_FRZ~")
      call tkadd(" OPTIONS_HF~")
      call tkadd(" OPTIONS_PFLOAT~")
      call tkadd(" SPEC")
      call tkadd(sscale)
C     call tkadd(" SPEC_SCLWSR~")
C     call tkadd(" SPEC_WSRMAX~")
C     call tkadd(" SPEC_OMAX1~")
C     call tkadd(" SPEC_OMAX2~")
      call tkadd(" SPEC_ATOM")
      call tkadd(" SPEC_ATOM_A~")
      call tkadd(" SPEC_ATOM_AMASS~")
      call tkadd(" SPEC_ATOM_C-HOLE~")
      call tkadd(" SPEC_ATOM_C-HQ~")
      call tkadd(" SPEC_ATOM_CSTRMX~")
      call tkadd(" SPEC_ATOM_EH")
      call tkadd(" SPEC_ATOM_EH2~")
      call tkadd(" SPEC_ATOM_EREF~")
      call tkadd(" SPEC_ATOM_FRZWF~")
      call tkadd(" SPEC_ATOM_IDMOD~")
      call tkadd(" SPEC_ATOM_IDU~")
      call tkadd(" SPEC_ATOM_IDXDN~")
      call tkadd(" SPEC_ATOM_ISOTOPE")
      call tkadd(" SPEC_ATOM_JH~")
      call tkadd(" SPEC_ATOM_KMXA~")
      call tkadd(" SPEC_ATOM_LFOCA~")
      call tkadd(" SPEC_ATOM_LMX~")
      call tkadd(" SPEC_ATOM_LMXA~")
      call tkadd(" SPEC_ATOM_LMXL~")
      call tkadd(" SPEC_ATOM_MMOM~")
      call tkadd(" SPEC_ATOM_NR~")
      call tkadd(" SPEC_ATOM_P~")
      call tkadd(" SPEC_ATOM_PZ~")
      call tkadd(" SPEC_ATOM_Q~")
      call tkadd(" SPEC_ATOM_R")
      call tkadd(" SPEC_ATOM_R/A")
      call tkadd(" SPEC_ATOM_RFOCA~")
      call tkadd(" SPEC_ATOM_RS3~")
      call tkadd(" SPEC_ATOM_RSMA~")
      call tkadd(" SPEC_ATOM_RSMFA~")
      call tkadd(" SPEC_ATOM_RSMG~")
      call tkadd(" SPEC_ATOM_RSMH")
      call tkadd(" SPEC_ATOM_RSMH2~")
      call tkadd(" SPEC_ATOM_R/W")
      call tkadd(" SPEC_ATOM_UH~")
      call tkadd(" SPEC_ATOM_Z")
c This call tkadd(" SPEC_ATOM_MIX~") was missing,
c then the initialization of default values for mxcst are not supplied
c   --->this results in a strange behevior!!
c takao Sep11-tukuba 2009
      call tkadd(" SPEC_ATOM_MIX~")
C ... for LMFRS
      call tkadd(" SPEC_ATOM_HCR~")
      call tkadd(" SPEC_ATOM_HCR/R~")
      call tkadd(" STR~")
c      call tkadd(" STR_ENV_MODE~")
c      call tkadd(" STR_ENV_NEL~")
c      call tkadd(" STR_ENV_EL")
      call tkadd(" STR_RMAX~ STR_RMAXS~")
c      call tkadd(" STR_RVL/R~ STR_VLFUN~")

      call tkadd(" SITE")
      call tkadd(" SITE_FILE~")
      call tkadd(" SITE_MODE~")
      call tkadd(" SITE_ATOM")
      call tkadd(" SITE_ATOM_DPOS~")
C     call tkadd(" SITE_ATOM_PL~")
      call tkadd(" SITE_ATOM_POS")
      call tkadd(" SITE_ATOM_XPOS")
      call tkadd(" SITE_ATOM_RELAX~")
      call tkadd(" SYMGRP~")
      call tkadd(" BZ")
      call tkadd(" BZ_BZJOB~")
      call tkadd(" BZ_DELEF~")
      call tkadd(" BZ_DMAT~")
      call tkadd(" BZ_DOS~")
      call tkadd(" BZ_EF0~")
      call tkadd(" BZ_EFMAX~")
      call tkadd(" BZ_FSMOM~")
      call tkadd(" BZ_GETQP~")
      call tkadd(" BZ_INVIT~")
      call tkadd(" BZ_METAL~")
      call tkadd(" BZ_MULL~")
      call tkadd(" BZ_N~")
      call tkadd(" BZ_NEVMX~")
      call tkadd(" BZ_NKABC")
      call tkadd(" BZ_NOINV~")
      call tkadd(" BZ_NPTS~")
      call tkadd(" BZ_N.W~")
      call tkadd(" BZ_PUTQP~")
      call tkadd(" BZ_SAVDOS~")
      call tkadd(" BZ_TETRA~")
      call tkadd(" BZ_W~")
      call tkadd(" BZ_ZBAK~")
      call tkadd(" BZ_ZVAL~")
      call tkadd(" EWALD~")
      call tkadd(" EWALD_AS~")
      call tkadd(" EWALD_NKDMX~")
      call tkadd(" EWALD_RPAD~")
      call tkadd(" EWALD_TOL~")
      call tkadd(" ITER~")
      call tkadd(" ITER_AMIX~")
      call tkadd(" ITER_CONV~")
      call tkadd(" ITER_CONVC~")
      call tkadd(" ITER_MIX~")
      call tkadd(" ITER_NITU~ ITER_TOLU~ ITER_UMIX~")
      call tkadd(" ITER_NRMIX~")
C     call tkadd(" ITER_XIPMX~")
      call tkadd(" ITER_NIT~")
C     call tkadd(" START~")
      call tkadd(" HAM")
      call tkadd(" HAM_NSPIN~")
      call tkadd(" HAM_ELIND~")
      call tkadd(" HAM_EWALD~")
      call tkadd(" HAM_FORCES~")
      call tkadd(" HAM_FRZWF~")
      call tkadd(" HAM_FTMESH")
      call tkadd(" HAM_GGA~")
      call tkadd(" HAM_GMAX")
ctakao
      call tkadd(" HAM_OVEPS~  HAM_OVNCUT~ ")
      call tkadd(" HAM_PMAX~")
      call tkadd(" HAM_PMIN~")
      call tkadd(" HAM_PWMODE~ HAM_PWEMIN~ HAM_PWEMAX~ HAM_NPWPAD~")
      call tkadd(" HAM_REL~ HAM_SO~")
      call tkadd(" HAM_RDSIG~")
      call tkadd(" HAM_RSRNGE~")
      call tkadd(" HAM_RSSTOL~")
      call tkadd(" HAM_SIGP~")
      call tkadd(" HAM_TOL~")
      call tkadd(" HAM_XCFUN~")
      call tkadd(" DYN~")
      call tkadd(dyn_mstat)
      call tkadd(dyn_md)
      call tkadd(" DYN_NIT~")

C --- LMDOS switches ---
      call nswadd()
      call tkadd(" LMDOS::" )
      call tkadd(" VERS_LM")
      call tkadd(trim(io_sw))
      call tkadd(trim(struc_sw))
      call tkadd(" OPTIONS~")
      call tkadd(" OPTIONS_FRZ~")
      call tkadd(" OPTIONS_HF~")
C     call tkadd(" OPTIONS_SCR~")
c      call tkadd(" OPTIONS_SHARM~")
C     call tkadd(" OPTIONS_STONER~")
      call tkadd(" SPEC")
      call tkadd(sscale)
C     call tkadd(" SPEC_SCLWSR~")
C     call tkadd(" SPEC_WSRMAX~")
C     call tkadd(" SPEC_OMAX1~")
C     call tkadd(" SPEC_OMAX2~")
      call tkadd(" SPEC_ATOM")
      call tkadd(" SPEC_ATOM_A~")
      call tkadd(" SPEC_ATOM_ALPHA~")
      call tkadd(" SPEC_ATOM_CSTRMX~")
      call tkadd(" SPEC_ATOM_DV~")
      call tkadd(" SPEC_ATOM_EREF~")
      call tkadd(" SPEC_ATOM_GROUP~")
      call tkadd(" SPEC_ATOM_GRP2~")
      call tkadd(" SPEC_ATOM_HCR~")
      call tkadd(" SPEC_ATOM_HCR/R~")
      call tkadd(" SPEC_ATOM_IDMOD~")
      call tkadd(" SPEC_ATOM_IDU~")
      call tkadd(" SPEC_ATOM_IDXDN~")
      call tkadd(" SPEC_ATOM_JH~")
      call tkadd(" SPEC_ATOM_LMX~")
      call tkadd(" SPEC_ATOM_LMXA~")
      call tkadd(" SPEC_ATOM_MIX~")
      call tkadd(" SPEC_ATOM_MMOM~")
      call tkadd(" SPEC_ATOM_NR~")
      call tkadd(" SPEC_ATOM_P~")
      call tkadd(" SPEC_ATOM_Q~")
      call tkadd(" SPEC_ATOM_R")
      call tkadd(" SPEC_ATOM_R/A")
      call tkadd(" SPEC_ATOM_R/W")
      call tkadd(" SPEC_ATOM_UH~")
      call tkadd(" SPEC_ATOM_Z")
      call tkadd(" SITE")
      call tkadd(" SITE_FILE~")
      call tkadd(" SITE_MODE~")
      call tkadd(" SITE_ATOM")
      call tkadd(" SITE_ATOM_DPOS~")
      call tkadd(" SITE_ATOM_PL~")
      call tkadd(" SITE_ATOM_POS")
      call tkadd(" SITE_ATOM_XPOS")
      call tkadd(" SYMGRP~")
      call tkadd(" BZ")
      call tkadd(" BZ_BZJOB~")
      call tkadd(" BZ_COND_MODE~")
      call tkadd(" BZ_COND_DIR")
      call tkadd(" BZ_DELEF~")
      call tkadd(" BZ_DMAT~")
      call tkadd(" BZ_DOS~")
      call tkadd(" BZ_EF0~")
      call tkadd(" BZ_EFMAX~")
      call tkadd(" BZ_FSMOM~")
      call tkadd(" BZ_GETQP~")
      call tkadd(" BZ_INVIT~")
      call tkadd(" BZ_METAL~")
      call tkadd(" BZ_MULL~")
      call tkadd(" BZ_N~")
      call tkadd(" BZ_NEVMX~")
      call tkadd(" BZ_NKABC")
      call tkadd(" BZ_NOINV~")
      call tkadd(" BZ_NPTS~")
      call tkadd(" BZ_N.W~")
      call tkadd(" BZ_PUTQP~")
      call tkadd(" BZ_SAVDOS~")
      call tkadd(" BZ_TETRA~")
      call tkadd(" BZ_W~")
      call tkadd(" BZ_ZBAK~")
      call tkadd(" BZ_ZVAL~")
      call tkadd(" EWALD~")
      call tkadd(" EWALD_AS~")
      call tkadd(" EWALD_NKDMX~")
      call tkadd(" EWALD_RPAD~")
      call tkadd(" EWALD_TOL~")
      call tkadd(" HAM_NSPIN~")
C      call tkadd(" ITER~")
C      call tkadd(" ITER_AMIX~")
C      call tkadd(" ITER_CONV~")
C      call tkadd(" ITER_CONVC~")
C      call tkadd(" ITER_MIX~")
C      call tkadd(" ITER_XIPMX~")
C --- LMCHK switches ---
      call nswadd()
      call tkadd(" LMCHK::" )
      call tkadd(" VERS_LM")
      call tkadd(trim(io_sw))
      call tkadd(trim(struc_sw))
      call tkadd(" OPTIONS~")
      call tkadd(" OPTIONS_NESABC~ OPTIONS_RMAXES~ OPTIONS_RMINES~")
      call tkadd(" SPEC")
      call tkadd(sscale)
C     call tkadd(" SPEC_SCLWSR~")
C     call tkadd(" SPEC_WSRMAX~")
C     call tkadd(" SPEC_OMAX1~")
C     call tkadd(" SPEC_OMAX2~")
      call tkadd(" SPEC_ATOM")
      call tkadd(" SPEC_ATOM_A~")
      call tkadd(" SPEC_ATOM_CSTRMX~")
      call tkadd(" SPEC_ATOM_DV~")
      call tkadd(" SPEC_ATOM_LMX~")
      call tkadd(" SPEC_ATOM_MMOM~")
      call tkadd(" SPEC_ATOM_NR~")
      call tkadd(" SPEC_ATOM_P~")
      call tkadd(" SPEC_ATOM_Q~")
c      call tkadd(" SPEC_ATOM_R")
c      call tkadd(" SPEC_ATOM_R/A")
c      call tkadd(" SPEC_ATOM_R/W")
      call tkadd(" SPEC_ATOM_R~")
      call tkadd(" SPEC_ATOM_R/A~")
      call tkadd(" SPEC_ATOM_R/W~")
      call tkadd(" SPEC_ATOM_Z")
      call tkadd(" SITE")
      call tkadd(" SITE_FILE~")
      call tkadd(" SITE_MODE~")
      call tkadd(" SITE_ATOM")
      call tkadd(" SITE_ATOM_DPOS~")
      call tkadd(" SITE_ATOM_PL~")
      call tkadd(" SITE_ATOM_POS")
      call tkadd(" SITE_ATOM_XPOS")
      call tkadd(" SITE_ATOM_ROT~")
      call tkadd(" SYMGRP~")
      call tkadd(" STR~")
C     call tkadd(" STR_DELRW~")
c      call tkadd(" STR_EQUIV~")
C     call tkadd(" STR_IINV~")
C      call tkadd(" STR_LMAXW~")
C      call tkadd(" STR_MODE~")
      call tkadd(" STR_MXNBR~")
      call tkadd(" STR_RMAX~")
      call tkadd(" STR_RMAXS~")
c      call tkadd(" STR_SHOW~")
      call tkadd(" EWALD~")
      call tkadd(" EWALD_AS~")
      call tkadd(" EWALD_NKDMX~")
      call tkadd(" EWALD_RPAD~")
      call tkadd(" EWALD_TOL~")
      call tkadd(" HAM_NSPIN~")
      call tkadd(" HAM_NONCOL~ HAM_SS~")
      call tkadd(" HAM_XCFUN~")
c      call tkadd(" PGF~")
c      call tkadd(" PGF_GFOPTS~")
c      call tkadd(" PGF_MODE~")
c      call tkadd(" PGF_PLATL")
c      call tkadd(" PGF_PLATR")
c      call tkadd(" PGF_SPARSE~")

      end subroutine toksw_init
C      subroutine snit
C      return
C      end
C
C      subroutine snot
C      return
C      end

