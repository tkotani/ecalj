      module m_lmfinit
      use m_ext,only :sname               ! sname contains extension. foobar of ctrl.foobar
      use m_struc_def,only: s_spec,s_site ! spec and site structures.
      use m_MPItk,only: master_mpi
      use m_lgunit,only: stdo,stdl
      implicit none
!! All data set to run lmfp.F rdctrl2 read ctrl file and set these data
!! rdctrls made from three stages. Search the word 'Stage' in the followings.
!! At the bottom of this code, I put old document for variables such as ctrl_*, ham_* and so on (search 'old doc')
!! We can believe the old document not completely, however, it may be a help to read the code.
!! We have to move data to the modules where data is generated.
!   ! We will have to reorganize redundant data. In addition, we need to protect all variables given in this module.
!      
! memo for iand (this will be removed)
!      i=53           ! i=00110101 binary (lowest order byte)
!      j=45           ! j=00101101 binary (lowest order byte)
!      k=iand(i,j)    ! k=00100101 binary (lowest order byte)
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   ng    :number of G-vectors
Ci   gv    :list of reciprocal lattice vectors G (gvlist.f)
Ci   tol   :tolerance in wave function precision:
Ci         :Here wave function cutoff for g-vectors whose wave
Ci          functions are less than tol.
      
!!---- initial settings read from ctrl and processed -------
      integer,parameter::  noutmx=48
      logical,parameter::  T=.true., F=.false.
      integer,parameter::  NULLI=-99999,nkap0=3,mxspec=256,lstrn=10000
      integer,parameter::  n0=10,nppn=12, nab=9, nrmx=1501,nlmx=64 !bndfp
      real(8),parameter::  NULLR =-99999, fs = 20.67098d0, degK = 6.3333d-6 ! defaults for MD
!     
      integer,protected::  io_show,io_help=0,nvario=0, lat_nkqmx,nat,lxcf !,irs4
      integer(2),protected:: nono
      character(lstrn),protected:: sstrnmix,sstrnsymg
      character(256),protected:: header,symg=' ',   symgaf=' '!for Antiferro
      logical,protected :: ham_frzwf,ham_ewald
      integer,protected:: nspc,procid,nproc,master=0,nspx, !,stdo,stdl
     &     ctrl_lxcf,ctrl_nbas,ctrl_lrel,!ctrl_lves, 
     &     maxit ,ctrl_nl,ctrl_nitmv,ctrl_nspin,
     &     ctrl_nspec,ctrl_nvario,ctrl_pfloat
     &     ,ham_lxcf,gga,ftmesh(3),nmto=0,lrsigx=0
     &     ,nsp=1,lrel=1,lso=0
      real(8),protected:: pmin(n0),pmax(n0),ham_pmax(10),ham_pmin(10),
     &     ctrl_defm(6), ctrl_wsrmax,ctrl_rmaxes,
     &     ctrl_omax1(3),ctrl_omax2(3),ctrl_sclwsr,ctrl_rmines, !ctrl_zbak,,ctrl_sdprm(5)
     &     tolft,elind,vmtz, scaledsigma
     &     ,ham_elind,ham_oveps,ham_scaledsigma 
c     &     ,ham_delta_stabilize 
!! ... OPTIONS
      integer,protected :: smalit, lstonr(3)=0,nl,lpfloat 
      real(8),protected:: rmines,rmaxes
      logical,protected :: lhf,lcd4
!! ... STRUC 
      real(8),protected:: dlat,alat=NULLR,dalat=NULLR,vol,avw !lat_slat(3,3),
      integer,protected:: nbas=NULLI,nspec 
!! ... SPEC
      real(8),protected:: omax1(3),omax2(3),wsrmax,sclwsr
      character*8,protected:: slabl(mxspec)
      integer,protected:: lmxbx=-1,lmxax,nkaph,nkapi
      logical,allocatable,protected:: mxcst2(:),mxcst4(:)
      integer,allocatable,protected:: 
     .     idu(:,:),lmxb(:),lmxa(:),idmod(:,:),
     .     kmxt(:),kmxv(:),lfoca(:),lmxl(:),nxi(:),nr(:), nmcore(:),
     &  nkapii(:),nkaphh(:) 
      real(8),allocatable,protected:: rsmh1(:,:),rsmh2(:,:),eh1(:,:),eh2(:,:),
     .     rs3(:),rham(:),alpha(:,:),ehvl(:,:), uh(:,:),jh(:,:),
     .     qpol(:,:),stni(:),
     .     pnu(:,:,:),qnu(:,:,:),      pnudefault(:,:),qnudefault(:,:),qnudummy(:,:),
     .     coreq(:,:),!mass(:),
     .     rg(:),rsma(:),rfoca(:), rsmfa(:),rcfa(:,:),
     .     exi(:,:),rmt(:),pz(:,:,:),
     .     amom(:,:),spec_a(:),z(:),eref(:),rsmv(:)
      character*(8),allocatable,protected:: coreh(:)
!! ... SITE
      character(8),protected:: alabl
      real(8),allocatable,protected :: pos(:,:)!,vel(:,:)!,vshft(:)
      integer,allocatable,protected :: ips(:),irlx(:,:),ndelta(:) !,ipl(:),plv(:)
      real(8),allocatable,protected :: delta(:,:),mpole(:),dpole(:,:)
      integer,allocatable,protected ::iantiferro(:) !reserved
!! ... BZ
      integer,protected:: bz_lshft(3)=0,
     .     bz_lmet,bz_n,bz_lmull,ctrl_ldos,bz_fsmommethod
      real(8),protected:: bz_efmax,bz_zval,bz_fsmom,
     .     bz_semsh(10),zbak,bz_lcond(4),bz_range=5d0,bz_dosmax
      logical,protected:: bz_lio2,bz_tetrahedron !ctrl_lmet2 
!! ... Ewald
      real(8),protected:: lat_as,lat_tol,lat_rpad
      integer,protected:: lat_nkdmx
!! ... STR
      real(8),protected:: str_rmax=nullr 
      integer,protected:: str_mxnbr 
!! ... Iteration, MIX
      character(128),protected :: iter_mix=' ' 
      real(8),protected:: etol,qtol !ctrl_tol(3)=1d-4
      integer,protected:: iter_maxit=1
      integer,protected:: mix_kill,mix_lxpot,mix_mmix,mix_mode,mix_nmix,mix_nsave !mix_nitu,
      real(8),protected:: mix_b,mix_bv,mix_elind,mix_tolu,mix_umix,mix_w(3),mix_wc
      character(8),protected:: mix_fn
!! DYN
      real(8),protected:: mdprm(6)
      integer,protected:: nitmv
!!      
      integer,protected:: pwmode,ncutovl     !ncutovl is by takao. not in lm-7.0beta.npwpad,
      real(8),protected:: pwemax,pwemin,oveps!,delta_stabilize
      integer , allocatable,protected ::  iv_a_oips (:)   ,  lpz(:),lpzex(:),lhh(:,:),iprmb(:)
!! CG coefficient      
      real(8) , allocatable,protected  ::  rv_a_ocg (:)
      real(8) , allocatable,protected  ::  rv_a_ocy (:)
      integer , allocatable,protected  ::  iv_a_oidxcg (:)
      integer , allocatable,protected  ::  iv_a_ojcg (:) 
!!      
      integer,protected:: ham_pwmode,ham_nkaph,ham_nlibu,
     &     ctrl_noinv,nlmax,mxorb,ctrl_lfrce,bz_nevmx,
     &     ham_nbf,ham_lsig,bz_nabcin(3)=NULLI,
     &     bz_ndos              !bndfp
c     &     lat_ldist=0
c      logical,protected:: irs1x10!,irsrot !(irsrot is not tested)
      real(8),protected:: ham_seref, bz_w,!,bz_ef0,bz_def0, !iors,bndfp
     &     lat_platin(3,3),lat_alat,lat_avw,lat_tolft,lat_gmaxin
      real(8),protected::
     &     lat_gam(1:4)=[0d0,0d0,1d0,1d0], !no deformation. See lattic !lat_dist0(3,3)=0d0,
     &     ctrl_mdprm(6)
      integer,protected:: lekkl!,lcplxp !now we fix lcplxp=1 and removed from m_lmfinit 
      integer,protected:: lmaxu,nlibu
      integer,allocatable,protected::lldau(:)
      logical,protected:: lpztail=.false.
      integer,protected:: leks,lrout,plbnd,  pot_nlma,  pot_nlml,ham_nspx !,ham_ldham(16)
c      integer,protected,target::ldham(16)
      integer,protected:: nlmto !total number of MTOs
      
!! For SOC 
      real(8),protected:: socaxis(3)
     
      integer,protected:: nitrlx,natrlx,pdim
      logical,protected:: xyzfrz(3)
      integer,allocatable:: indrx_iv(:,:)
      real(8),protected:: defm(6)
      
!! sspec and ssite are unprotected (be careful)
      type(s_spec),allocatable:: v_sspec(:) !(nspec: number of species in the cell)
      type(s_site),allocatable:: v_ssite(:) !(nbas: number of atoms)
      
!! ... old memo for molecular dynamics section DYN
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
C     arg  18:        offset to bs workspace
!!ssssssssssssssssssssssssssssssssssssssssssssssssss      
      contains 
      subroutine m_lmfinit_init(prgnam)
      use m_rdfiln,only: recln,nrecs,recrd
      use m_struc_func
      use m_toksw,only:tksw
      use m_gtv
!! All the initial data are pushed into m_lmfinit
!! ----------------------------------------------------------------------
!! Inputs
!!   recrd (recln*nrecs) : preprocessed input
!!   prgnam: name of main program
!! Outputs are give in m_lmfinit
!     ! We have three stages (stage 1, stage 2 , stage 3) in this routine. Search 'stage'.
!! folloings are memo --- but not so clear. some variables are reduntant (because of historical reason).
Co   bz_*  : Brillouin Zone related
Co   ctrl_* : struct for program flow parameters; see routine uctrl
Co   ham_*  :  Hamiltonian related
Co   pot_*   : information about the potential; see routine upot
Co   lat_*   : lattice information; see routine ulat
Co   smix_*  : charge mixing parameters; see routine umix
Co  sspec : species-specific information; see routine uspec
Co  ssite : site-specific information; see routine usite
!!      
C---- note. Followings are just hints about what variables means
Cg Global variables
Cn   slabl :vector of species labels
Cg   avw   : the average Wigner-Seitz radius,
Cg   lrel  :specifies type of Schrodinger equation
Cg         :0 nonrelativistic Schrodinger equation
Cg         :1 scalar relativistic Schrodinger equation
Cg         :2 Dirac equation (reserved)
Cg   lxcf  :specifies type of XC potential.  1s digit specifies local XC:
Cg         :1 for Ceperly-Alder
Cg         :2 for Barth-Hedin (ASW fit)
Cg         :103 for PBE
Cg   mxorb :nkaph * (maximum number of lm channels in any sphere)
Cg         :Used for dimensioning the indexing arrays involved in
Cg         :assembling the hamiltonian;
Cg   nbas  :number of atoms in the basis
Cg   nkaph :The maximum number of radial functions centered at
Cg         :particular R and l channel used in the lmto basis. +1 when we have lo
Cg   nl    :1+Maximum l-cutoff for augmentation
Cg   npl   :(not set by rdctrl) number of principal layers (layer geometries)
Cg   nkaphh :The maximum number of "principal quantum" numbers centered at a particular R and l channel --- energies for one Rl
Cg         :at which augmentation (phi-phidot) functions are made. +1 is added when we have lo
Cg   nsp   :1 if not spin-polarized; otherwise 2. nspc and nspx are 
Cg   nspec :number of species
Cg   stde  :standard error file
Cg   stdl  :standard log file
Cg   stdo  :standard output file
Cr  Read input data specified by tokens (see m_toksw.F)
C ----------------------------------------------------------------------
      implicit none
c#if MPI | MPIK
      include "mpif.h"
c#endif
      character,intent(in)::  prgnam*(*)
      character strn*(recln),strn2*(recln)
      integer:: i_spec
      character fileid*64
      logical lgors,cmdopt,ltmp,ioorbp,noinv,cmdopt0,cmdopt2
      double precision dval,dglob,xx(n0*2),dgets !,ekap(6)
      integer i,is,iprint,
     .iprt,isw,ifi,ix(n0*nkap0),j,k,l,!lbas, !noinv was iand(lqp,1),irs(5)
     .lfrzw,!lgunit,
     .lrs,lstsym,ltb,nclasp,
     .nglob,scrwid,k1,k2,mpipid !nkap,lsx1, !,stde
      character*(8),allocatable::clabl(:)
      integer,allocatable:: ipc(:),initc(:),ics(:)
      real(8),allocatable:: pnuc(:,:,:),qnuc(:,:,:,:),pp(:,:,:,:),ves(:),zc(:)
      integer:: dvec1(3)=1, dvec2(3)=0
C ... basis
      double precision orbp(n0,2,nkap0)
      integer ohave,oics,opnu,opp,oqnu,osgw,osoptc,oves,owk !osordn,
      real(8):: pnux(20) ,temp33(9)
      integer:: nnn
      integer:: i_copy_size,i_spacks,iendx,inix,i_spackv
      real(8):: seref
      integer:: ib , ispec
      integer,allocatable:: wowk(:)
      logical:: isanrg,l_dummy_isanrg
      integer:: lmxcg,lmxcy,lnjcg,lnxcg,nlm
      
      integer::nout,nn,i0,ivec(10),iosite
      integer:: io_tim(2),verbos!io_iactive,
      character(256)::  a,outs
      logical::  mlog=.false.
      integer:: lmxbj,lmxaj,nlbj
      double precision vsn,vers,xv(2*n0)
      character(256*16) :: bigstr=' '
      integer:: it
c      logical :: parmxp
      logical::  debug=.false.
      integer:: lp1,lpzi
      real(8):: xxx
      real(8)::  avwsr
      integer::  izerv(n0)=(/(0,i=1,n0)/)
      real(8)::  zerov(n0)=(/(0d0,i=1,n0)/)
      integer:: ii,sw
      real(8):: dasum!,dglob
      character(128) :: nm
      real(8):: nullrv(256),d2,plat(3,3),rydberg
      integer:: nulliv(256),jj(2) !,nkapsi
      logical:: noelind
      integer:: levelinit=0
      integer:: lx,lxx
      character*256:: sss
      logical:: sexist
      integer:: ibas,ierr,lc, iqnu=0
c      logical:: ltet
      integer:: ifzbak,nn1,nn2,nnx,lmxxx,nlaj
      integer,allocatable :: iv_a_oips_bakup (:)
      integer :: ctrl_nspec_bakup,inumaf,iin,iout,ik,iprior,ibp1,indx,iposn,m,nvi,nvl,setprint0
      logical :: ipr10,fullmesh,lzz
      integer,parameter:: maxp=3
      integer,allocatable:: idxdn(:,:,:)
!! Stage 1  read ctrl file ---------------------    
      procid = mpipid(1)
      nproc  = mpipid(0)
      scrwid = 80
c      stdo = lgunit(1)
c      stdl = lgunit(2)
c      stde = stdo
      nullrv = nullr
      nulliv  =nulli
      debug = cmdopt0('--debug')
      call gtv_setrcd(recrd,nrecs,recln,stdo,stdl,stde_in=stdo) !Copy recrd to rcd in m_gtv
      call toksw_init(debug)
! io_show io_help      
      if (     master_mpi)   io_show = 1
      if (.not.master_mpi)   io_show = 0
      if (cmdopt0('--help')) io_help = 1
      if (io_help == 1) then
        write(stdo,*)' Token           Input   cast  (size,min) --------------------------'
      elseif(io_show/=0) then
        write(stdo,*)' Token           Input   cast  (size,min,read,def)     result'
      endif
      call gtv_setio(debug,io_show,io_help) ! In case io_help changed
!! IO
      i0=setprint0(30)      !initial verbose set
      nm='IO_VERBOS'; call gtv(trim(nm),tksw(prgnam,nm),verbos,
     .     note='Verbosity for printout. Set from the command-line with --pr=xxx',def_i4=30,nout=nout)
      if(io_help==1) verbos=30
      if(cmdopt2('--pr=',outs)) then !overwrite verbos via console input
         read(outs,*) verbos 
      elseif(cmdopt2('--pr',outs)) then
         read(outs,*) verbos 
      elseif(cmdopt2('-pr',outs))then
         read(outs,*) verbos
      endif
      i0 = setprint0(verbos)      !Set initial verbos
      if(.not.master_mpi) i0=setprint0(-100) !iprint() is negative except master
!! show timing      
      nm='IO_TIM';call gtv(trim(nm),tksw(prgnam,nm),io_tim,note='Turns CPU timing log. Value sets tree depth.'//
     .     '%N   Optional 2nd arg prints CPU times as routines execute.'//
     .     '%N   Args may be set through command-line: --time=#1[,#2]', def_i4v=(/1,1/),nmin=1,nout=i0)
      if(i0==1) io_tim(2)=io_tim(1)
      if (cmdopt2('--time',outs) ) then !!  Override with '--time=' commmand-line arg
         outs=trim(outs(2:))//' 999 999'
         read(outs,*)io_tim(1),io_tim(2)
         if(io_tim(1)==999) io_tim(1)=5
         if(io_tim(2)==999) io_tim(2)=2
         i0=1
      endif
      if ( i0 >=1 ) call tcinit(io_tim(2),io_tim(1),levelinit)
      call tcn('m_lmfinit')
!! CONST --- for backword compatibility. We will remove this.
      call numsyv(nvario)
      nm='CONST'; call gtv(trim(nm),tksw(prgnam,nm), bigstr, note='Constants may declared for use in expressions'//
     .     '%N   Variables may also be set from the command-line:  -vnam=#',nout=nout)
      if (nout .eq. 1) then
         i = 0
         call parsyv(bigstr,len_trim(bigstr),1999,0,i)
         if (io_show .ne. 0) call shosyv(0,0,0,stdo)
      endif
!! Struc ---
      if (tksw(prgnam,'STRUC') == 2) goto 59
      if (io_show+io_help/=0) write(stdo,*)' --- Parameters for crystal structure ---'
      nm='STRUC_ALAT';call gtv(trim(nm),tksw(prgnam,nm),alat, note= 'Scaling of lattice vectors, in a.u.')
      nm='STRUC_NBAS';call gtv(trim(nm),tksw(prgnam,nm),nbas,note='Size of basis')
      nm='STRUC_PLAT';call gtv(trim(nm),tksw(prgnam,nm),temp33,nmin=9,nout=nout,note='Primitive lattice vectors (alat unit)')
      plat= reshape(temp33,[3,3])
      avw = avwsr(plat,alat,vol,nbas)
!  'STRUC_NSPEC' removed
      if(io_help == 0) then ! .and. nspec .eq. NULLI) then
         nm='SPEC_ATOM'; sw = tksw(prgnam,nm)
         if (sw .ne. 2) then
            j = 0; nspec = 0
            do while (nspec <= 0)
               j = j+1; jj= (/1,j/)
               if (.not. debug) call pshpr(0)
               call gtv(trim(nm),0,nono,Texist=ltmp,cindx=jj)
               if (.not. debug) call poppr
               if (.not. ltmp) nspec = j-1
            enddo
            if (io_show>0) write(stdo,"(' ... found',g0,'species in SPEC category')")nspec
         endif
      endif
      nm='STRUC_DALAT'; call gtv(trim(nm),tksw(prgnam,nm),dalat,def_r8=0d0,
     &     note='added to alat after reading inputs (only affecting to SPEC_ATOM_R/A case)')
!     STRUC_NL here removed, Lattice distortion or rotatation here removed.
   59 continue
!! Options ---
      if (tksw(prgnam,'OPTIONS') == 2) goto 69
      if (io_show+io_help/=0) write(stdo,*)' --- Program Options ---'
      nm='OPTIONS_HF'; call gtv(trim(nm),tksw(prgnam,nm),lhf,def_lg=F,note='T for non-self-consistent Harris')
      nm='OPTIONS_RMINES'; call gtv(trim(nm),tksw(prgnam,nm),rmines,def_r8=1d0,note='Minimum MT radius when finding new ES')
      nm='OPTIONS_RMAXES'; call gtv(trim(nm),tksw(prgnam,nm),rmaxes,def_r8=2d0,note='Maximum MT radius when finding new ES')
      lpfloat=1
 69   continue                  ! end of OPTIONS
      
!! Hamiltonian parameters ---
      if(io_show+io_help/=0.and.tksw(prgnam,'HAM')/=2)write(stdo,*)' --- Parameters for hamiltonian ---'
      nsp = 1
      nm='HAM_NSPIN';call gtv(trim(nm),tksw(prgnam,nm),nsp,def_i4=1,note='Set to 2 for spin polarized calculations')
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
     .     '%N%3f1 for scalar relativistic Schrodinger equation'//'%N%3f2 for Dirac equation')
!  Fully relativistic => spin-orbit coupling
!  lso  =0 (no so): =1(L.S): =2(LzSz).
      if (lrel==2) lso=1
      if (nsp==2 .or. io_help/=0) then
         if (io_help /= 0) write(stdo,*)' * To read the magnetic parameters below, HAM_NSPIN must be 2'
         nm='HAM_SO'; call gtv(trim(nm),tksw(prgnam,nm), lso,def_i4=0,note=
     .        'Spin-orbit coupling (for REL=1)'//
     .        '%N%3f0 : no SO coupling'//
     .        '%N%3f1 : Add L.S to hamiltonian'//
     .        '%N%3f2 : Add Lz.Sz only to hamiltonian') !//
c     .        '%N%3f3 : Like 2, but also compute <L.S-LzSz> by perturbation')
         if (io_help==0) l_dummy_isanrg=isanrg(lso,0,2,' rdctrl:','SO',T)
      endif
      
!! SOC Spin-block matrix Aug2021 ! Taken from (A8) in Ke.Liqin2019,PhysRevB.99.054418
      nm='HAM_SOCAXIS'; call gtv(trim(nm),tksw(prgnam,nm),socaxis,nmin=3,nout=nout, def_r8v=[0d0,0d0,1d0], !default
     &     note='SOC axis! 0,0,1(default) or 1,1,0 only effective for HAM_SO=1')
! note:  We register HAM_SOCAXIS~ in toksw_init (extention ~ means optional token)
! Sanity check      
!  In general cases (except 001), we nees spin-symmetirc radial functions because we
!  use  <up|Lz|up> for the place of <up|Lz|dn>, for example.
!  If we like to take into account spin-dependent radial functions,
!     We need to calculate ohsozz%soffd, ohsopm%sdiag. See 'call gaugm' at the end of augmat.
      if(io_help==0) then
         if( sum(abs(socaxis-[0d0,0d0,1d0])) >1d-6 .and. (.not.cmdopt0('--phispinsym')))
     &     call rx('We need --phispinsym for SO=1 and HAM_SOCAXIS/=001')
      endif
!!      
      sw = tksw(prgnam,'HAM_GMAX')
      if(sw/=2) then
         nm='HAM_GMAX';call gtv(trim(nm),sw,lat_gmaxin,nmin=1,nout=nout,note='Energy cutoff for plane-wave mesh',or=T)
         if (nout .ne. 0) then
            sw = 2
         else
            lat_gmaxin = 0
         endif
         nm='HAM_FTMESH'; call gtv(trim(nm),sw,ftmesh,nout=nout, note='No. divisions for plane-wave mesh '//
     .        'along each of 3 lattice vectors.'//
     .        '%N%3fSupply one number for all vectors or a separate '//
     .        'number for each vector.')
         call fill3in(nout,ftmesh)
      endif
      nm='HAM_TOL'; call gtv(trim(nm),tksw(prgnam,nm),tolft, def_r8=1d-6, note='w.f. tolerance for FT mesh')
      nm='HAM_FRZWF'; call gtv(trim(nm),tksw(prgnam,nm),ham_frzwf,def_lg=F,
     .     note='Set to freeze augmentation wave functions for all species')
      nm='HAM_FORCES'; call gtv(trim(nm),tksw(prgnam,nm),ctrl_lfrce, def_i4=0,note=
     .     'Controls the ansatz for density shift in force calculation.'//
     .     '%N%3f-1 no force%3f0 no shift'//'%N%3f 1 free-atom shift  12 screened core+nucleus')
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
      nm='HAM_RDSIG'; call gtv(trim(nm),tksw(prgnam,nm),lrsigx,def_i4=1, note=
     .     'Controls how self-energy is added to '//
     .     'local exchange correlation functional:'//
     .     '%N%3f   0: do not read Sigma'//
     .     '%N%3f   1(or not zero): read sigm=Sigma-Vxc. Default now')
      nm='HAM_ScaledSigma'; call gtv(trim(nm),tksw(prgnam,nm),scaledsigma,
     .     def_r8=1d0, note='=\alpha_Q for QSGW-LDA hybrid. \alpha \times (\Sigma-Vxc^LDA) is added to LDA/GGA Hamiltonian.')
      nm='HAM_EWALD'; call gtv(trim(nm),tksw(prgnam,nm),ham_ewald, def_lg=.false.,note='Make strux by Ewald summation')
      nm='HAM_VMTZ'; call gtv(trim(nm),tksw(prgnam,nm),vmtz,def_r8=0d0, note='Muffin-tin zero defining wave functions')
      nm='HAM_PMIN'; call gtv(trim(nm),tksw(prgnam,nm),pmin, def_r8v=zerov,nout=nout,note=
     .     'Global minimum in fractional part of P-functions.'//
     .     '%N%3fEnter values for l=0..:'//
     .     '%N%3f0: no minimum constraint'//
     .     '%N%3f#: with #<1, floor of fractional P is #'//
     .     '%N%3f1: use free-electron value as minimum')
      nm='HAM_PMAX'; call gtv(trim(nm),tksw(prgnam,nm), pmax, def_r8v=zerov, nout=nout, note=
     .     'Global maximum in fractional part of P-functions.'//
     .     '%N%3fEnter values for l=0..:'//
     .     '%N%3f0: no maximum constraint'//
     .     '%N%3f#: with #<1, ceiling of fractional P is #')
!      We set default oveps=1d-7 16Nov2015. This was zero before the data.
      nm='HAM_OVEPS'; call gtv(trim(nm),tksw(prgnam,nm), oveps, def_r8=1d-7, nout=nout, note=
     .     'Diagonalize hamiltonian in reduced hilbert space,'//
     .     '%N%3fdiscarding part with evals of overlap < OVEPS')
      if(cmdopt0('--zmel0')) OVEPS=0d0
c      nm='HAM_STABILIZE'; call gtv(trim(nm),tksw(prgnam,nm),
c     .     delta_stabilize, def_r8=-1d0, nout=nout,note=
c     .     'Experimental. Stabilizer for Diagonalize hamiltonian (negative means unused),'//
c     .     '%N%3f "H --> H + HAM_STABILIZE*O^-1" in zhev_tk(diagonalization)')
!     ... APW basis
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
!! Symmetry group ---
      if (io_show+io_help/=0 .and. tksw(prgnam,'SYMGRP')/=2)write(stdo,*)' --- Symmetry group operations ---'
      nm='SYMGRP'; call gtv(trim(nm),tksw(prgnam,nm),symg, note='Generators for symmetry group')
!   for AF --- !june2015
      if(.not.(prgnam=='LMFA'.or.prgnam=='LMCHK')) then
        nm='SYMGRPAF'; call gtv(trim(nm),tksw(prgnam,nm),symgaf,
     .  note='One (or multiple) Extra Generator for adding anti ferro symmetry')
      endif  
!! Species
      if (io_help == 1) nspec = 1
      if (tksw(prgnam,'SPEC') == 2) goto 79
      if (io_show+io_help/=0) write(stdo,*)' --- Parameters for species data ---'
      if (io_help /= 0) write(stdo,*)' * The next four tokens apply to the automatic sphere resizer'
      nm='SPEC_SCLWSR'; call gtv(trim(nm),tksw(prgnam,nm), sclwsr, def_r8=0d0, note=
     .     'Scales sphere radii, trying to reach volume = '//
     .     'SCLWSR * cell volume'//
     .     '%N%3fSCLWSR=0 turns off this option.'//
     .     '%N%3fAdd  10  to initially scale non-ES first;'//
     .     '%N%3f or  20  to scale ES independently.')
      nm='SPEC_OMAX1'; call gtv(trim(nm),tksw(prgnam,nm),omax1, def_r8v=(/0d0,0d0,0d0/),note=
     .     'Limits max sphere overlaps when adjusting MT radii')
      nm='SPEC_OMAX2'; call gtv(trim(nm),tksw(prgnam,nm), omax2, def_r8v=(/0d0,0d0,0d0/),note=
     .     'Sphere overlap constraints of second type',nout=nout)
      nm='SPEC_WSRMAX'; call gtv(trim(nm),tksw(prgnam,nm),wsrmax, def_r8=0d0,note=
     .     'If WSRMAX is nonzero, no sphere radius may exceed its value')
      if (io_help .eq. 1) then
         write(*,382)
 382     format(/' * ',
     .        'The following tokens are input for each species. ',
     .        'Data sandwiched'/3x,'between successive occurences of ',
     .        'token ATOM apply to one species.')
         nspec = 1
      endif
!! SPEC_ATOM_*      
      if (nspec .eq. 0) goto 79
      allocate(pnu(n0,nsp,nspec),qnu(n0,nsp,nspec), pz(n0,nsp,nspec),amom(n0,nspec),idmod(n0,nspec),
     .    rsmh1(n0,nspec),eh1(n0,nspec),rsmh2(n0,nspec),eh2(n0,nspec),
     .    ehvl(n0,nspec), qpol(n0,nspec),stni(nspec),
     .    rg(nspec),rsma(nspec),rfoca(nspec),rsmfa(nspec),rcfa(2,nspec),
     .    rs3(nspec),rham(nspec),rmt(nspec),rsmv(nspec),
     .    nxi(nspec),exi(n0,nspec), spec_a(nspec),z(nspec),nr(nspec),eref(nspec),
     .    coreh(nspec),coreq(2,nspec), idxdn(n0,nkap0,nspec), idu(4,nspec),uh(4,nspec),jh(4,nspec),
     .    mxcst2(nspec),mxcst4(nspec), kmxt(nspec),kmxv(nspec),lfoca(nspec),lmxl(nspec),lmxa(nspec),lmxb(nspec),nmcore(nspec))
      allocate(lpz(nspec),lpzex(nspec))
      allocate(nkapii(nspec),nkaphh(nspec))
      lpz=0
      lpzex=0
      mxcst2=F
      mxcst4=F
      nkapii=1
      nkapi = 1
      lpzi = 0 
      qpol = NULLR
      rsmh1 = 0d0
      rsmh2 = 0d0
      eh1  = 0d0
      eh2 = 0d0
      idmod = NULLI
      ehvl = NULLR
      do 1111 j = 1, nspec
        if(debug) print *,'nspec mxcst j-loop j nspec',j,nspec
        rcfa(:,j) = NULLR; rfoca(j) = 0d0; rg(j) = 0d0
        rham(j) = NULLR; rsma(j) = 0d0; rsmfa(j) = 0d0
        spec_a(j) = NULLR; nr(j) = NULLI
        exi(:,j) = NULLR
        coreh(j) = ' '; coreq(:,j) = NULLR
        eref(j) = NULLR
        if (io_help /= 0) then
          write(stdo,'(1x)')
        elseif (io_help == 0 .and. io_show>0) then
          write(stdo,ftox)' ... Species ',j
        endif
        jj= (/1,j/)
        nm='SPEC_ATOM'; call gtv(trim(nm),tksw(prgnam,nm),slabl(j), nmin=10,cindx=jj,note='Species label')
        nm='SPEC_ATOM_Z'; call gtv(trim(nm),tksw(prgnam,nm),z(j), cindx=jj,note='Atomic number')
        sw = tksw(prgnam,'SPEC_ATOM_R')
        if (sw /= 2) then
           nout = 0
           nm='SPEC_ATOM_R';call gtv(trim(nm),sw,rmt(j),cindx=jj,nout=nout,note= 'Augmentation sphere radius rmax',or=T)
           if (nout .ne. 1) then
              nm='SPEC_ATOM_R/W';call gtv(trim(nm),sw,rmt(j),cindx=jj,nout=nout,note='rmax relative to average WS radius',or=T)
              if (nout == 1) then
                 rmt(j) =rmt(j)*avw
              else
                 nm='SPEC_ATOM_R/A';call gtv(trim(nm),sw,rmt(j),cindx=jj,nout=nout,note='rmax ratio to alat')
                 if (nout == 1) then
                    rmt(j) =rmt(j)*alat
                 else
                    rmt(j) = 0d0 !!     takao for lmchk even when R is not given. See default of LMCHK
                 endif
              endif
           endif
        endif
        if(debug) print *,'nspec bbb mxcst j-loop j nspec',j,nspec
!    Radial mesh parameters: determine default value of a
        i0 = NULLI
        xxx = NULLR
        if (io_help .ne. 1) then
           call pshpr(0)
           call rmesh(z(j),rmt(j),lrel,.false.,nrmx,xxx,i0)
           call poppr
           if (xxx .eq. .03d0) xxx = .015d0 !.025d0 jun2012 .025 to .015 as default.
        endif
        nm='SPEC_ATOM_A'; call gtv(trim(nm),tksw(prgnam,nm),spec_a(j),def_r8=xxx,cindx=jj,nout=nout,
     .       note='Radial mesh point spacing parameter')
!     Determine default NR
        if (tksw(prgnam,'SPEC_ATOM_NR') .ne. 2) then
           i0 = 0
           call pshpr(0)
           call rmesh(z(j),rmt(j),lrel,.false.,nrmx,spec_a(j),i0)
           call poppr
        endif
        nm='SPEC_ATOM_NR'; call gtv(trim(nm),tksw(prgnam,nm),nr(j),def_i4=i0,cindx=jj, note='Number of radial mesh points')
        if (nr(j) .eq. 0) nr(j) = i0
!... Basis set for lmf  1st MTO sets ----------
        nm='SPEC_ATOM_RSMH'; call gtv(trim(nm),tksw(prgnam,nm),rsmh1(1:,j),cindx=jj,nout=nout,def_r8v=zerov,
     .       note='Smoothing radii for basis. Gives l-cut max for base')!nout=n0 because of default
        nnx=nout
        do i=1,nout
           if(rsmh1(i,j)==0d0) cycle
           nnx=i
        enddo
        nm='SPEC_ATOM_EH'; call gtv(trim(nm),tksw(prgnam,nm), eh1(1:,j),nmin=nnx,cindx=jj, def_r8v=zerov,
     .       note='Kinetic energies for basis')
        nn1=nnx
!---- 2nd MTO sets       
        nm='SPEC_ATOM_RSMH2'; call gtv(trim(nm),tksw(prgnam,nm),rsmh2(1:,j),def_r8v=zerov,cindx=jj,nout=nout,
     .       Texist=ltmp,note='Basis smoothing radii, second group ')
        nnx=0
        if (ltmp) then
           nnx=nout
           do i=1,nout
              if(rsmh2(i,j)==0d0) cycle
              nnx=i
           enddo
           sw = tksw(prgnam,nm)
           if (nnx>0) then
              nkapi=2
              nkapii(j)=2
              sw = 1
           endif
           nm='SPEC_ATOM_EH2'; call gtv(trim(nm),sw,eh2(1:,j),nmin=nnx,cindx=jj, !def_r8v=zerov,
     .          note='Basis kinetic energies, second group')
        endif
        nn2=nnx
        lmxbj=max(nn1,nn2)-1
! optional cutoff for l-base. redundunt, I think.
        nm='SPEC_ATOM_LMX'; call gtv(trim(nm),tksw(prgnam,nm),lmxxx,
     .       def_i4=10,cindx=jj,nout=nout,note='optional l-cutoff for basis')
        lmxbj=min(lmxxx,lmxbj)
        lmxb(j)=lmxbj
        lmxbx=max(lmxbx,lmxbj)
!     ... Determine lmxa: floating orbitals have no lmxa
        nm='SPEC_ATOM_LMXA'; sw = tksw(prgnam,nm)
        if(debug) print *,'nspec ddd000 mxcst j rmt io_help lmxa=',j,rmt(j),io_help,lmxa(j)
        if (rmt(j) == 0.and.io_help/=1) then !takao iohelp/=1 added.
           lmxa(j) = -1
        elseif (sw == 2) then   !lmxa not read: look for subsitute
           lmxa(j) = 4
        else ! default is lmxbj
           call gtv(trim(nm),sw,lmxa(j),def_i4=lmxbj,cindx=jj,Texist=ltmp,note='l-cutoff for augmentation')
        endif
        lmxaj = lmxa(j)
!     nlaj = number of elements associated with lmxa
!     0 => no elements
!     nlaji: ditto, but used to specify number of default values
        nlaj = 1+lmxaj
        if (io_help .eq. 1) nlaj = 1
        if(debug) print *,'nspec ddd111 mxcst j-loop j nspec nlaj lmxaj=',j,nspec,nlaj,lmxaj
!     ... Parameters that depend on the existence of an augmentation sphere
!     lmxl = l-cutoff for numerical rep'sn of density in sphere
!     lfoca = mode for treating core
!     kmxt = kmax for expansion of envelope wf tails
!     kmxv = cutoff to expand smoothed potential
!     hardwired for now
        kmxv(j) = 15            !Not input
!     Cannot set default here: need set after rescaling of rmt
!     rsmv(j) = rmt(j)*.5d0   !Not input
        rsmv(j) = 0d0           !Not input
        kmxt(j) = -1            !If sought, default will be set below
        lfoca(j) = 0            !If sought, default will be reset below
        nmcore(j)=0
        lmxl(j) =  lmxaj        !Use lmxaj in case not sought (ASA:mpol)
        nxi(j) = NULLI          !If sought, default will be set below
        pnu(:,:,j)=0d0
        pz(:,:,j)=0d0
        rs3(j) = NULLI          !If sought, default will be set below
        idxdn(:,:,j) = 1
        pnu(1,1,j) = NULLI      !If sought, default will be set below
        qnu(1,1,j) = NULLI      !If sought, default will be set below
        if (nlaj .ne. 0) then
           nm='SPEC_ATOM_LMXL'; call gtv(trim(nm),tksw(prgnam,nm),lmxl(j),
     .          cindx=jj,def_i4=lmxaj,note='lmax for which to accumulate rho,V in sphere')
!     ... Set up default P,Q in absence of explicit specification
           pnu(:,:,j)=0d0
           qnu(:,:,j)=0d0
! -- takao move back default value of dev_r8v to zero june2012 --
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
!! ==== Reset default P,Q in absence of explicit specification ====
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
              if (io_show .ne. 0) call poppr
           endif
           if (nsp .eq. 2) call dcopy(n0,pnu(1,1,j),1,pnu(1,2,j),1)
           if (nsp .eq. 2 .or. io_help == 1) then
              nm='SPEC_ATOM_MMOM'; call gtv(trim(nm),tksw(prgnam,nm),
     .             qnu(1:nlaj,2,j),def_r8v=zerov,cindx=jj,note=
     .             'Starting mag. moms for each l channel.'//
     .             '%N%2f For a chanel with PZ, this is enforced to be zero.'//
     .             '%N%2f See explanation for SPEC_ATOM_Q.')
           endif
           nm='SPEC_ATOM_NMCORE'; call gtv(trim(nm),tksw(prgnam,nm),nmcore(j), !takao jun2012
     .          def_i4=0,cindx=jj,note='spin-averaged core: jun2012takao'//
     .          '%N%3f0(default): spin-polarized core'//
     .          '%N%3f1         : spin-averaged core density is from spin-averaged potential')
           nm='SPEC_ATOM_PZ'; call gtv(trim(nm),tksw(prgnam,nm),pz(1:nlaj,1,j),def_r8v=zerov,cindx=jj,note=
     .          'Starting semicore log der. parameters'//
     .          '%N%10fAdd 10 to attach Hankel tail',nout=nout) !zero default by zerov
           if(nsp==2) pz(1:n0,2,j) = pz(1:n0,1,j) !takao
          
!! lmxb corrected by pz
           nnx=nout
           do i=nout,1,-1
              if(pz(i,1,j)/=0d0) then
                 nnx=i
                 lmxb(j)=max(lmxb(j),nnx-1)
                 exit
              endif   
           enddo
!           
           if (nout>0) then
              if (dasum(nlaj,pz(1,1,j),1) /= 0) then
                 lpzi = max(lpzi,1) !,2
                 lpz(j)=1
                 if ( sum(int(pz(1:nlaj,1,j)/10))>0 ) then
                    lpzex(j)=1
                 endif   
              endif   
           endif
           if(maxval(pz(1:nout,1,j))>10d0) lpztail=.true. ! PZ +10 mode exist or not.
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
     .          def_i4=i0,cindx=jj,note='FOCA switch 0(within MT):'
     .          //'=1(frozenCore). Default: 1 for z>8;0 for z<=8') !takao Aug2010
           if(io_help==0) then
              if(lfoca(j)/=0.and.lfoca(j)/=1) then
                 call rx('LFOCA should be 0 or 1 (Aug2010): 2 is not allowed')
              endif
           endif
           nm='SPEC_ATOM_KMXA'; call gtv(trim(nm),tksw(prgnam,nm),
     .          kmxt(j),def_i4=3,cindx=jj,note='k-cutoff for projection of wave functions in sphere.')
C     Cannot set default here: need set after rescaling of rmt
C     xxx = 0.40d0*rmt(j); if (io_help .eq. 1) xxx = NULLI
           xxx = 0d0; if (io_help .eq. 1) xxx = NULLI
           nm='SPEC_ATOM_RSMA'; call gtv(trim(nm),tksw(prgnam,nm),rsma(j),
     .          def_r8=xxx,cindx=jj,note='Smoothing for projection of wave functions in sphere.'//
     .          '%N%3finput<0 => choose default * -input')
           xxx = 0d0; if (io_help .eq. 1) xxx = NULLI
           nm='SPEC_ATOM_RSMG'; call gtv(trim(nm),tksw(prgnam,nm),rg(j),
     .          def_r8=xxx,cindx=jj,note='Smoothing for projection of charge in sphere.'//
     .          '%N%3finput<0 => choose default * -input')
           xxx = 0d0; if (io_help .eq. 1) xxx = NULLI
           nm='SPEC_ATOM_RFOCA'; call gtv(trim(nm),tksw(prgnam,nm),rfoca(j),
     .          def_r8=xxx,cindx=jj,note='Smoothing for core tail.  input<0 => choose default * -input')
           xxx = 0d0; if (io_help .eq. 1) xxx = NULLI
           nm='SPEC_ATOM_RSMFA'; call gtv(trim(nm),tksw(prgnam,nm),rsmfa(j), def_r8=xxx,cindx=jj,note=
     .          'Smoothing for free atom.  input<0 => choose default * -input')
           nm='SPEC_ATOM_RCFA'; call gtv(trim(nm),tksw(prgnam,nm),rcfa(1:2,j)
     .          ,def_r8v=zerov,nmin=2,cindx=jj,note=
     .          'Cutoff radius for renormalization of free atom density'//
     .          '(WARN:takao rnatm.F is not tested).'//'%N%3fOptional 2nd argument = width'//
     .          '%N%3fRCFA<0 => renormalize potential instead of density')
C     Negative radii: convert to actual numbers
c          nm='SPEC_ATOM_IDXDN'; removed...
           if (rg(j)    .lt. 0) rg(j)    = -rg(j)*0.25d0*rmt(j)
           if (rsma(j)  .lt. 0) rsma(j)  = -rsma(j)*0.4d0*rmt(j)
           if (rfoca(j) .lt. 0) rfoca(j) = -rfoca(j)*0.4d0*rmt(j)
           if (rsmfa(j) .lt. 0) rsmfa(j) = -rsmfa(j)*0.5d0*rmt(j)
           nm='SPEC_ATOM_RS3'; call gtv(trim(nm),tksw(prgnam,nm),rs3(j),
     .          def_r8=0.5d0,cindx=jj, !takao Nov26,2010
     .          note='Minimum smoothing radius for local orbital')
          nm='SPEC_ATOM_IDMOD'; call gtv(trim(nm),tksw(prgnam,nm),
     .         idmod(1:nlaj,j),def_i4v=(/(0,i=1,n0)/), cindx=jj,note=
     .         'idmod=0 floats P to band CG, 1 freezes P, 2 freezes enu')
! mxcst2(for lmchk)   Exclude this species when auto-resizing sphere radii
! mxcst4   Freeze augmentation w.f. for this species (FP)
          nm='SPEC_ATOM_CSTRMX'; call gtv(trim(nm),tksw(prgnam,nm),
     .         mxcst2(j),cindx=jj,def_lg=F,note='Set to exclude this'//
     .         ' species when automatically resizing sphere radii (SCLWSR>0)') !for lmchk
          if (sclwsr .eq. 0) mxcst2(j) = F
          nm='SPEC_ATOM_FRZWF'; call gtv(trim(nm),tksw(prgnam,nm),mxcst4(j),cindx=jj,def_lg=F,note=
     .         'Set to freeze augmentation wave functions for this species')
        endif                    ! end of input dependent on presence of aug sphere.
       
        sw = tksw(prgnam,'SPEC_ATOM_IDU')
        if (io_help .gt. 0 .and. sw .lt. 2) write(stdo,*)' * ... The next three tokens are for LDA+U'
        nm='SPEC_ATOM_IDU'; call gtv(trim(nm),sw,idu(:,j),cindx=jj, def_i4v=(/(0,i=1,n0)/),note=
     .       'LDA+U mode:  0 nothing, 1 AMF, 2 FLL, 3 mixed; +10: nothin if sigm.* exist')
        nm='SPEC_ATOM_UH';call gtv(trim(nm),sw,uh(:,j),cindx=jj,def_r8v=zerov,note='Hubbard U for LDA+U')
        nm='SPEC_ATOM_JH';call gtv(trim(nm),sw,jh(:,j),cindx=jj,def_r8v=zerov,note='Exchange parameter J for LDA+U')
!! 2019 automatic turn off lda+u mode; Enforce uh=jh=0 when sigm exist !!
        if(procid==master .and. sum(abs(idu(:,j)))/=0) then 
           inquire(file='sigm.'//trim(sname),exist=sexist)
           if(sexist) then
              do lxx=0+1,3+1
                 if(idu(lxx,j)>10) then
                    write(stdo,"(a,2i4)")'For IDU>10 with sigm.*, we set UH=JH=0 for l,ibas=',lxx,j
                    uh(lxx,j) = 0d0
                    jh(lxx,j) = 0d0
                 endif
              enddo  
           endif
           do lxx=0+1,3+1
              idu(lxx,j) = mod(idu(lxx,j),10)
           enddo
        endif
        call mpibc1_int(idu(:,j),4,'m_lmfinit_idu')
        call mpibc1_real(uh(:,j),4,'m_lmfinit_uh') !bug fix at oct26 2021 (it was _int-->not passed )
        call mpibc1_real(jh(:,j),4,'m_lmfinit_jh')
!     Sanity checks
        if (io_help == 0 .and. lmxaj .ge. 0) then
           if(tksw(prgnam,'SPEC_ATOM_LFOCA')/= 2)l_dummy_isanrg=isanrg(lfoca(j),0,2,'rdctrl','lfoca',T)
           if(tksw(prgnam,'SPEC_ATOM_LMXL')/= 2) l_dummy_isanrg=
     .          isanrg(lmxl(j),min(0,lmxaj),max(0,lmxaj),'rdctrl','lmxl',T)
           if(tksw(prgnam,'SPEC_ATOM_KMXA')/= 2)l_dummy_isanrg=isanrg(kmxt(j),2,25,' rdctrl (warning):','kmxa',F)
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
        nkaphh(j) = nkapii(j) + lpzex(j)
 1111 continue  ! Loop over species
      if (io_help==0) then
         lmxax = maxval(lmxa) !Maximum L-cutoff
         nlmax = (lmxax+1)**2
         nkaph = nkapi + lpzi     !-1
         mxorb= nkaph*nlmax
      endif
 79   continue
!! Site ---
      if(io_show+io_help/=0.and.tksw(prgnam,'SITE')/=2)write(stdo,*)' --- Parameters for site ---'
      if(io_help .eq. 1) then
         nbas = 1
         if (iprint() .gt. 0) write(*,383)
 383     format(/' * ',
     .        'The following tokens are input for each site. ',
     .        'Data sandwiched'/3x,'between successive occurences of ',
     .        'token ATOM apply to one site.'/3x,'Alternatively, all ',
     .        'site data can be read in via the SITE file.')
      endif
      allocate(pos(3,nbas),ips(nbas),irlx(3,nbas),iantiferro(nbas))
! vshft(nbas),,ipl(nbas),plv(nbas)!mpole(nbas),dpole(3,nbas),,vel(3,nbas)
c      vel  = 0d0
      irlx = 0
      ips  = NULLI
      pos  = NULLR
!SITE_ATOM_*
      do  j = 1, nbas
         if (io_help /= 0) then
            write(stdo,'(1x)')
         elseif (io_help == 0 .and. io_show>0) then
            write(stdo,ftox)' ... Site ',j
         endif
         jj=(/1,j/)
         nm='SITE_ATOM'; call gtv(trim(nm),tksw(prgnam,nm),alabl,nmin=10, cindx=jj,note='Species label')
         if(io_help .ne. 1) then
            do  i = 1, nspec
               if (trim(alabl) == trim(slabl(i)) ) then
                  ips(j) = i
                  goto 881
               endif
            enddo
            call rxs('Category SITE referred to nonexistent species: ',alabl)
         endif
 881     continue
!  ... Site positions
         sw= tksw(prgnam,'SITE_ATOM_XPOS')
         nm='SITE_ATOM_POS'; call gtv(trim(nm),tksw(prgnam,nm),pos(:,j),
     .        nout=nout,cindx=jj,note='Atom coordinates, cartesian in alat', or=(sw.ne.2))
         if(nout==0.or.tksw(prgnam,'SITE_ATOM_POS')==2)then !nout=-1 if sw=2; otherwise nout=0 unless data was read
            nm='SITE_ATOM_XPOS'; call gtv(trim(nm),tksw(prgnam,nm),pos(:,j),
     .           cindx=jj,note='Atom coordinates, fractional in alat')
            call dcopy(3,pos(1,j),1,xv,1)
            call dmpy(plat,3,1,xv,3,1,pos(1,j),3,1,3,1,3)
         endif
c         nm='SITE_ATOM_V0'; call gtv(trim(nm),tksw(prgnam,nm),vel(:,j),def_r8v=zerov,cindx=jj,
c     .        note='Initial velocity for molecular dynamics')
         nm='SITE_ATOM_RELAX'; call gtv(trim(nm),tksw(prgnam,nm),irlx(:,j),
     .        def_i4v=(/(1,i=1,n0)/),cindx=jj,note=
     .        'relax site positions (lattice dynamics) or Euler angles (spin dynamics)')
         nm='SITE_ATOM_AF'; call gtv(trim(nm),tksw(prgnam,nm),iantiferro(j), cindx=jj,def_i4=0,
     .        note='antiferro ID:=i and -i should be af-pair, we look for space-group operation with spin-flip')
      enddo
 89   continue
!! Structure constants ---
      if(io_show+io_help/=0.and.tksw(prgnam,'STR')/=2)write(stdo,*)' --- Parameters for structure constants ---'
      nm='STR_RMAXS'; call gtv(trim(nm),tksw(prgnam,nm),str_rmax,
     .     nout=nout,note='Radial cutoff for strux, in a.u.',or=T)
      if (nout == 0) then       !nout=-1 if sw=2; otherwise nout=0 unless data was read
         nm='STR_RMAX'; call gtv(trim(nm),tksw(prgnam,nm),str_rmax,
     .        def_r8=0d0,note='Radial cutoff for strux, in units of avw')
         str_rmax = str_rmax*avw
      endif
      nm='STR_MXNBR'; call gtv(trim(nm),tksw(prgnam,nm),str_mxnbr,def_i4=0,note='Max number of nbrs (for dimensioning arrays)')
!! Brillouin Zone ---
      if(io_show+io_help/=0.and.tksw(prgnam,'BZ')/=2)write(stdo,*)' --- Parameters for Brillouin zone integration ---'
      nm='BZ_NKABC'; sw=tksw(prgnam,nm)!; if (bz_lio1) sw = 2
      call gtv(trim(nm),sw,bz_nabcin,nout=nout,
     .     note='No. qp along each of 3 lattice vectors.'//'%N%3fSupply one number for all vectors or a separate '//
     .     'number for each vector.')
      call fill3in(nout,bz_nabcin)
      nm='BZ_BZJOB';call gtv(trim(nm),tksw(prgnam,nm),bz_lshft,nout=nout,def_i4v=izerv(1:1),note=
     .     '0 centers BZ mesh at origin, 1 centers off origin'//
     .     '%N%3fSupply one number for all vectors or a separate '//
     .     'number for each vector.')
      call fill3in(nout,bz_lshft)
      nm='BZ_METAL'; call gtv(trim(nm),tksw(prgnam,nm),bz_lmet,
     .     def_i4=3,note='0 insulator only; 3 for metal (2 is for maintenance)')
      if(io_help==0.and.bz_lmet/=3.and.bz_lmet/=0.and.bz_lmet/=2) call rx('BZ_METAL error')
      nm='BZ_TETRA'; call gtv(trim(nm),tksw(prgnam,nm),bz_tetrahedron, !tetrahedron switch
     .     def_lg=T,note='Tetrahedron integration')
      if(cmdopt0('--tdos').or.cmdopt0('--pdos')) then
         write(stdo,*)' --tdos or --pdos enforces BZ_METAL=3 and BZ_TETRA=T'
         bz_lmet=3
         bz_tetrahedron=.true.
      endif
      nm='BZ_N'; call gtv(trim(nm),tksw(prgnam,nm),bz_n, def_i4=0,note=
     .     'N>0: Polynomial order for Methfessel-Paxton sampling%N%3f'//
     .     'N=0: Conventional Gaussian sampling%N%3f'//
     .     'N<0: Broadening by Fermi-Dirac distribution%N%3f'//
     .     'To be used in conjunction with W= ; see next')
      nm='BZ_W'; call gtv(trim(nm),tksw(prgnam,nm),bz_w, def_r8=5d-3,note=
     .     'If BZ_N>=0, Line broadening for sampling integration%N%3f'//
     .     'If BZ_N<0,  Temperature for Fermi distribution (Ry)')
c BZ_EF0, BZ_DELEF removed. !c!! remove writing ZBAK file here. (Write ZBAK file. sep2020)
      nm='BZ_ZBAK'; call gtv(trim(nm),tksw(prgnam,nm),zbak,def_r8=0d0,note='Homogeneous background charge')
      nm='BZ_SAVDOS'; call gtv(trim(nm),tksw(prgnam,nm),ctrl_ldos,
     .     def_i4=0,note='Choose 0 or 1: %N%3f1 Write dos.tot.* file (settings are NPTS and DOS)')
      nm='BZ_NPTS'; call gtv(trim(nm),tksw(prgnam,nm),bz_ndos,def_i4=2001,
     .     note='No. DOS points (sampling integration)')
      nm='BZ_DOSMAX'; call gtv(trim(nm),tksw(prgnam,nm),bz_dosmax,def_r8=40d0/rydberg(),
     &     note='Maximum energy to which DOS accumulated, relative to Efermi') !march 2013
      xxx = 5d0
      nm='BZ_EFMAX'; call gtv(trim(nm),tksw(prgnam,nm),bz_efmax,
     .     def_r8=xxx,note='Find evecs up to efmax')
      nm='BZ_NEVMX'; call gtv(trim(nm),tksw(prgnam,nm),bz_nevmx,
     .     def_i4=0,note='Find at most nevmx eigenvectors'//
     .     '%N%3fIf NEVMX=0, program uses internal default'//
     .     '%N%3fIf NEVMX<0, no eigenvectors are generated')
      if( cmdopt0('--tdos').or.cmdopt0('--pdos').or.cmdopt0('--zmel0')) bz_nevmx=999999
      nm='BZ_NOINV'; call gtv(trim(nm),tksw(prgnam,nm),noinv, def_lg=F,note=
     .     'Suppress automatic inclusion of inversion symmetry for BZ')
      nm='BZ_FSMOM'; call gtv(trim(nm),tksw(prgnam,nm),bz_fsmom,
     .     def_r8=NULLR,note='Fixed-spin moment (fixed-spin moment method)')
      nm='BZ_FSMOMMETHOD';call gtv(trim(nm),tksw(prgnam,nm),bz_fsmommethod,
     .     def_i4=0,note='Method of Fixed-spin moment 0:original 1:discrete')
!! Ewald sums ---
      if (io_show+io_help/=0 .and. tksw(prgnam,'EWALD')/=2) write(stdo,*)' --- Parameters for Ewald sums ---'
      nm='EWALD_AS'; call gtv(trim(nm),tksw(prgnam,nm),lat_as, def_r8=2d0,note='Ewald smoothing parameter')
      nm='EWALD_TOL'; call gtv(trim(nm),tksw(prgnam,nm),lat_tol, def_r8=1d-8,note='Ewald tolerance')
      nm='EWALD_NKDMX'; call gtv(trim(nm),tksw(prgnam,nm),lat_nkdmx, def_i4=800,note='Ewald tolerance')
!! Iterations (formerly MIX) ---
      mix_b = NULLI            ! Not set
      if (tksw(prgnam,'ITER')/=2) then
         if (io_show+io_help/=0) write(stdo,*)' --- Parameters for iterations ---'
!     Default values for smix (array has same same structure as lstra smix)
         mix_b = 1             ! beta
         mix_bv = 1            ! bv
         mix_fn='mixm'
         mix_kill =  0         ! nkill
         mix_lxpot =  0        ! lxpot
         mix_mmix = -1         ! mmix
         mix_mode = 0          ! mode (0=Anderson)
         mix_nsave = 8         ! nsave = # iter to save on disk
         mix_tolu = 0          ! tolu
         mix_umix = 1          ! umix (mixing parm for LDA+U)
         mix_w(1) = 1          ! w(1)
         mix_w(2) = 1          ! w(2)
         mix_wc = -1           ! wc
         smalit = NULLI
         nm='ITER_NIT';call gtv(trim(nm),tksw(prgnam,nm),iter_maxit,def_i4=30,
     .        note='maximum number of iterations in self-consistency cycle')
         nm='ITER_NRMIX'; call gtv(trim(nm),tksw(prgnam,nm),smalit,def_i4=80,note='Sphere program, max iter')
! we use mixrho.F ->parmxp.F. But too complicated to touch it.
         nm='ITER_MIX'; sw=tksw(prgnam,nm); call gtv(trim(nm),sw,iter_mix,
     .        nmin=10,nout=nout,note='Mixing rules for charge mixing.  Syntax:')
         if(io_help/=0 .and. tksw(prgnam,nm)/=2) print 345
 345     format(3x,'A[nmix][,b=beta][,bv=betv][,n=nit][,w=w1,w2][,nam=fn][,k=nkill]','[;...] or'/
     .          3x,'B[nmix][,b=beta][,bv=betv][,wc=wc][,n=#][,w=w1,w2][,nam=fn]','[,k=nkill]')
         nm='ITER_CONV';call gtv(trim(nm),tksw(prgnam,nm),etol,def_r8=1d-4,note=
     &    'Tolerance in energy change from prior iteration for self-consistency')
         nm='ITER_CONVC'; call gtv(trim(nm),tksw(prgnam,nm),qtol,
     .        def_r8=1d-4,note='Tolerance in output-input charge for self-consistency')
         nm='ITER_UMIX';call gtv(trim(nm),tksw(prgnam,nm),mix_umix,def_r8=.5d0,note='Mixing parameter for densmat in LDA+U') !2022mar9 default umix=0.5
         nm='ITER_TOLU';call gtv(trim(nm),tksw(prgnam,nm),mix_tolu,def_r8=0d0,note='Tolerance for densmat in LDA+U')
      endif                     ! iterations category
!! Dynamics (relaxation)---
      if(io_show+io_help/=0.and.tksw(prgnam,'DYN')/=2)write(stdo,*)' --- Parameters for dynamics and statics ---'
      nm='DYN_NIT'; call gtv(trim(nm),tksw(prgnam,nm),nitmv,def_i4=1,
     .     note='maximum number of relaxation steps (statics)'//' or time steps (dynamics)')
      nm='DYN_MSTAT'; call gtv(trim(nm),tksw(prgnam,nm),nono,Texist=ltmp,note='Parameters for molecular statics')
      mdprm(1) = 0
      if(io_help/=0 .or. ltmp) then
         if(ctrl_lfrce==0) ctrl_lfrce=1 ! takao july2012
         nm='DYN_MSTAT_MODE'; call gtv(trim(nm),tksw(prgnam,nm),i0,
     .        def_i4=0,note=
     .        '0: no relaxation  '//
     .        '4: conjugate gradients  '//
     .        '5: Fletcher-Powell  '//
     .        '6: Broyden')
         mdprm(1) = i0
         nm='DYN_MSTAT_HESS'; call gtv(trim(nm),tksw(prgnam,nm),ltmp,def_lg=T,note='Read hessian matrix')
         mdprm(2) = isw(ltmp)! T=>1 F=>0
         nm='DYN_MSTAT_XTOL'; call gtv(trim(nm),tksw(prgnam,nm),mdprm(3),def_r8=1d-3,note=
     .      'Convergence criterion in displacements'//'%N%3fXTOL>0: use length;  <0: use max val;  =0: do not use')
         nm='DYN_MSTAT_GTOL'; call gtv(trim(nm),tksw(prgnam,nm),mdprm(4),def_r8=0d0,note=
     .      'Convergence criterion in gradients'//
     .      '%N%3fGTOL>0: use length;  <0: use max val;  =0: do not use')
         nm='DYN_MSTAT_STEP'; call gtv(trim(nm),tksw(prgnam,nm),mdprm(5),def_r8=0.015d0,note=
     &      'Initial (and maximum) step length')
         nm='DYN_MSTAT_NKILL'; call gtv(trim(nm),tksw(prgnam,nm),i0,def_i4=0,note='Remove hessian after NKILL iter')
         mdprm(6) = i0
      endif
      if (io_help>0) then
         write(stdo,"(a)")'==============================================='
         call lmhelp(prgnam)
         call rx0('end of help mode')
      endif   
!!  end of read input parameter ccccccccccccccccccccccccccxxxxxxxxxxxxxxxx

!! Stage 2  Reorganize ctrl_* in module m_lmfinit ---------------------
      ctrl_noinv = isw(noinv)  ! T->1 F->0
      ctrl_lrel=lrel
      ctrl_lxcf= ham_lxcf !1 for Ceperly-Alder (VWN),  2 for Barth-Hedin (ASW fit), 103 for PBE-GGA
      maxit=iter_maxit
      ctrl_mdprm=mdprm
      nl = max(lmxbx,lmxax)+1 !max l-base l-aug +1
      nlmax=nl**2
      ctrl_nitmv=nitmv
      ctrl_nl=nl
      ctrl_nbas=nbas
      ctrl_nspec=nspec
      ctrl_nspin=nsp
      ctrl_nvario=nvario
      ctrl_omax1 = omax1
      ctrl_omax2 = omax2
      ctrl_rmaxes= rmaxes
      ctrl_rmines= rmines
      ctrl_sclwsr= sclwsr
      ctrl_wsrmax= wsrmax
      ctrl_pfloat= lpfloat
      if (dalat == NULLR) dalat=0d0
      lat_alat=alat+dalat
      lat_avw=avw
      lat_nkqmx=lat_nkdmx
      lat_platin=plat
      lat_tolft=tolft
!! set cg coefficients for lmf --- Choose dimensions for arrays
      lmxcg=8
      lmxcy=12
c      if (lmxcg .le. 6) then
c        lnjcg = 6500
c        lnxcg = 1300
c      else if (lmxcg .le. 8) then
        lnjcg = 22700
        lnxcg = 3400
c      else if (lmxcg .le. 10) then
c        lnjcg = 62200
c        lnxcg = 7400
c      else
c        call rxi('setcg: cannot handle lmxcg=',lmxcg)
c      endif
      nlm=(lmxcy+1)**2
      allocate(rv_a_ocy(abs(nlm)))
      allocate(rv_a_ocg(abs(lnjcg)))
      allocate(iv_a_ojcg(abs(lnjcg)))
      allocate(iv_a_oidxcg(abs(lnxcg)))
      call sylmnc ( rv_a_ocy , lmxcy )
      call scg ( lmxcg , rv_a_ocg , iv_a_oidxcg , iv_a_ojcg )
!!      
      ham_elind=elind
      ham_nkaph=nkaph
      ham_pmax=pmax 
      ham_pmin=pmin 
      if (procid==master) then
         inquire(file='sigm.'//trim(sname),exist=sexist)
         if (lrsigx/=0 .and. (.not.sexist) ) then
            write(stdo,*)' bndfp (warning): no sigm file found ... LDA calculation only'
            lrsigx = 0
         endif
      endif
      ham_lsig=lrsigx
      call mpibc1_int(ham_lsig,1,'bndfp_ham_lsig')
      ham_scaledsigma=scaledsigma
      ham_pwmode=pwmode
      ham_oveps=oveps
!
!! idxdn= 1 or 3  
Cr     1  Orbital is included as "active" orbitals, which means
Cr           they are included in the hamiltonian and diagonalized
Cr     3   Neglected.
cxxxxxxxxxx          
Cr     10  Orbital is a local orbital whose value and slope are
Cr         constructed to be zero at the MT boundary.
Cr         It is included in the basis.
Cr     11  Orbital is a local orbital with a smooth Hankel tail
Cr         and it is included in the basis.
      do j=1,nspec
         do  ik = 1, nkap0
            if (ik <= nkapi) then
               do lp1 = 1, lmxb(j)+1
                  if (ik==1.and.rsmh1(lp1,j)<=0) idxdn(lp1,ik,j)=3
                  if (ik==2.and.rsmh2(lp1,j)<=0) idxdn(lp1,ik,j)=3
               enddo
            endif
            if (ik == nkaph .and. sum(lpz)>0) then
               idxdn(:,ik,j)=3  !call ivset(idxdn(1,ik,j),1,n0,4)
               do  lp1  = 1, lmxb(j)+1
                  if (pz(lp1,1,j) .ne.  0) then
                     if(pz(lp1,1,j)>=10) idxdn(lp1,ik,j)=1 !11
                     if(pz(lp1,1,j)>  0) idxdn(lp1,ik,j)=1 !10
                  endif
               enddo
            endif
            if (ik > nkaph) idxdn(:,ik,j)=3 !call ivset(idxdn(1,ik,j),1,n0,4)
            idxdn(lmxb(j)+ 2:,ik,j)=3 
         enddo
      enddo
c      ham_delta_stabilize=delta_stabilize !takao sep2010
C --- Allocate and copy input to sspec ---
      allocate(v_sspec(nspec))
      do j=1,nspec
          v_sspec(j)%name=slabl(j)
          v_sspec(j)%z=z(j)
          v_sspec(j)%eh3=-0.5d0
          v_sspec(j)%vmtz=-0.5d0
          v_sspec(j)%coreh=coreh(j)
          v_sspec(j)%nmcore=nmcore(j)
          v_sspec(j)%a=spec_a(j)
          v_sspec(j)%nr=nr(j)
          v_sspec(j)%coreq=coreq(:,j)
          v_sspec(j)%nxi=nxi(j)
          v_sspec(j)%exi=exi(:,j)
          v_sspec(j)%idmod=idmod(:,j)
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
          v_sspec(j)%p = pnu(1:n0,1,j) !only for lmfa. fixed.
          v_sspec(j)%pz= pz(1:n0, 1,j) !only for lmfa. fixed.
          v_sspec(j)%q = qnu(1:n0,1:nsp,j) 
          v_sspec(j)%rcfa=rcfa(:,j)
          v_sspec(j)%rfoca=rfoca(j)
          v_sspec(j)%rg=rg(j)
          v_sspec(j)%rmt=rmt(j)
          v_sspec(j)%rs3=rs3(j)
          v_sspec(j)%rsma=rsma(j)
          v_sspec(j)%rsmfa=rsmfa(j)
          v_sspec(j)%eref=eref(j)
      enddo
      allocate(v_ssite(nbas))
      do j=1,nbas
          v_ssite(j)%spec=ips(j)
          v_ssite(j)%class=ips(j)
          v_ssite(j)%pos=pos(1:3,j)
          v_ssite(j)%relax=irlx(:,j)
          v_ssite(j)%iantiferro=iantiferro(j)
          is=v_ssite(j)%spec
          v_ssite(j)%pnu= pnu(1:n0,1:nsp,is) 
          v_ssite(j)%pz = pz(1:n0,1:nsp,is)  
      enddo
      sstrnmix=trim(iter_mix)

!! ... Suppress symmetry operations for special circumstances
!     !     Switches that automatically turn of all symops
!     ! --pdos mar2003 added. Also in lmv7.F
      lstsym = 0
      if ((mdprm(1)>=1.and.mdprm(1)<=3) .or.   
     &     cmdopt0('--cls') .or. cmdopt0('--nosym').or.cmdopt0('--pdos')) then
         symg = 'e'
         lstsym = 2             !lstsym=2: turn off symops
         ctrl_noinv=1
      endif
      sstrnsymg=trim(symg)
      if (cmdopt0('--rdbasp')) call rx('not support --rdbasp')

!!  Add dalat to alat
!!   lat_alat=(lat_alat)+dalat !this is a bug; this should be commented out. This
!!   is a bug for lm7K when I started lm7K. I had included this bug here. Fixed at 28May2010.

!! Dirac equation requires spin polarization
      if( nsp==1 .and. int(ctrl_lrel)==2 ) call rx('rdccat: Dirac equation requires NSPIN=2')
!! ... Suppress inversion when noncollinear magnetism, SX, NLO
      if(lso.ne. 0) ctrl_noinv=1 !lqp = lqp-bitand(lqp,1)+1
      
!! ... Setup for idxdn ... ctrl->lham,4 is automatic downfolding switch
c      j=3+10 ! no downfolding j=3 ; j=+10 nfp-style basis:
c      call suidx(nkaph,j,nspec,v_sspec)
!! Set default values for species data
      call defspc(v_sspec)
!!      
      lxcf = ctrl_lxcf
      nspc = 1
      if( lso==1 ) nspc = 2
     
      
      allocate( iprmb(nbas * nl**2 * maxp ) )
      iprmb=-1
      nlmto = 0
      do 110 ib = 1,nbas
         is = v_ssite(ib)%spec
         iposn = mxorb*(ib-1)
         do 1121 ik = 1, nkaph
           do  l = 0, nl-1
           do  m = -l, l
             iposn = iposn+1
             if(idxdn(l+1,ik,is)==1) then
               nlmto = nlmto+1
               iprmb(iposn) = nlmto
             endif
           enddo 
           enddo
 1121     continue
 110   continue
      nspx  = nsp
      if(lso/=0) nspx = 1
      nvi=0
      nvl=0
      do  ib = 1, nbas
        is=v_ssite(ib)%spec
        nvi = nvi + (lmxa(is)+1)**2
        nvl = nvl + (lmxl(is)+1)**2
      enddo
      pot_nlma=nvi
      pot_nlml=nvl
           
C     Make nat = number of real atoms as nbas - # sites w/ floating orbitals
      if (procid .eq. master) then
         nat = nbas
         do  i = 1, nbas
            j=v_ssite(i)%spec
            l=v_sspec(j)%lmxa
            if (l .eq. -1) nat = nat-1
         enddo
      endif
      call mpibc1_int(nat,1,'m_lmfinit_nat')
      allocate(wowk(nbas))
      wowk=0
      call pshpr(0)
      call suldau(nbas,v_sspec,v_ssite,nlibu,k,wowk)!Count LDA+U blocks (printout only)
      ham_nlibu=nlibu
      call poppr
      deallocate(wowk,pnu,qnu    ,amom,idmod,                  
     .     qpol,stni,rg,rsma,rfoca,rsmfa,rcfa,nxi,
     .     exi,coreq,rs3,rham,idxdn,
     .     rmt,idu,uh,jh,kmxv,  !kmxt,
     .     lfoca,eref,lmxl,coreh,          spec_a,z,nr,rsmv,
     .     pos,ips,irlx)  !,vshft
!! --- takao embed contents in susite here. This is only for lmf and lmfgw.
      seref = 0d0
      do ib = 1, nbas
         is = v_ssite(ib)%spec
         if (is <= 0) cycle
         seref = seref +     v_sspec(is)%eref
      enddo
      ham_seref=seref
      allocate(iv_a_oips(nbas))
      do ib=1,nbas
         iv_a_oips(ib) = v_ssite(ib)%spec
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
      if (io_help == 0 .and. io_show > 1) then
         print *, '---------- contents of sstrn ------------'
         print *, 'mix: ', trim(sstrnmix)
         print *, 'symg:', trim(sstrnsymg)
         call rx0('done show')
      endif
      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
      call MPI_BARRIER( MPI_COMM_WORLD, ierr )
      if( cmdopt0('--quit=show') ) call rx0(trim(prgnam)//' --quit=show')

      
!! Stage 3: initial settings ===================================
!!  Total energy mode setting
      if (cmdopt0('--etot')) then
        mdprm=0 !!!  No forces or dynamics
        ctrl_lfrce=0
        ctrl_mdprm=0
c        irs2=0
        maxit=1
      endif
      if(lhf) maxit= 1
!! scalar rel       
      call setcc(lrel) !lrel/=0 means scalar relativistiv c=274.074d0 in a.u.
!!
      fullmesh = cmdopt0('--fullmesh').or.cmdopt0('--fermisurface') !fullmesh stop just after do 2010 iq loop.
      lrout = 1 ! Whether to evaluate output density and/or KS energy
      leks = 1
      if (bz_nevmx .eq. -1) then !bz_nevmx ! Maximum number of eigenvalues
        lrout = 0
        leks = 0
      endif
      if(cmdopt0('--band').or.fullmesh) then !Switch to plot bands at specified qp
         ctrl_lfrce = 0
         plbnd = 1
         lrout = 0
      else
         plbnd=0
      endif
      if(lrout == 0 ) maxit = 1
      if(bz_lmet/=0 .or. bz_tetrahedron ) ctrl_ldos=1
      if(ctrl_ldos==0) bz_ndos = 1
!!
      if(lrout == 0 .and. ctrl_lfrce /= 0) then
        write(stdo,"(a)") 'Error: output density required when forces sought.\n'//
     &              '      To make output density turn off HF=t and/or NEVMX<0'
        call Rx('incompatible input. see the end of console output')
      endif
      if(lrout == 0 .and. cmdopt0('--etot')) then
        write(stdo,"(a)") 'Error: output density required with --etot switch.\n'//
     &              '      To make output density turn off HF=t and/or NEVMX<0'
        call Rx('incompatible input. see the end of console output')
      endif
      
!! LDA+U block (it was suldau.F)
      nlibu = 0
      lmaxu = 0
      allocate(lldau(nbas))
      lldau = 0
      do  ib = 1, nbas
        is = v_ssite(ib)%spec
        do  lx = 0, min(v_sspec(is)%lmxa,3)
          if (v_sspec(is)%idu(lx+1) .ne. 0) then
            if (lldau(ib) ==0) lldau(ib) = nlibu+1
            nlibu = nlibu+1
            lmaxu = max(lmaxu,lx)
          endif
        enddo
      enddo
      if (master_mpi.and.nlibu/=0) write(stdo,"(' LDA+U: block(s) lmaxu =',2i4)") nlibu,lmaxu
!! aug2012 we now fix lcplxp=1 (complex ppi integral)      
      lekkl  = ctrl_pfloat
! lhh, nkapii, nkaphh (nkaphh = nkapii(1 or 2) +1) if extented local orbital exist)
      allocate(lhh(nkap0,nspec))
      lhh=-1
      do i=1,nspec
         lmxbj = lmxb(i)
         call getiout(rsmh1(1,i), lmxbj+1,lhh(1,i))
         if(nkapii(i)==2) call getiout(rsmh2(1,i),lmxbj+1,lhh(2,i))
         if(lpz(i)==1 )   call getiout(pz(1,1,i),lmxbj+1,lhh(nkaph,i))!lh for lo
      enddo
      if(master_mpi) then
         write(stdo,*)
         write(stdo,"('mmm === MTO setting ===')")
         do i=1,nspec
          lmxbj = lmxb(i)
          write(stdo,"('mmm ispec lmxb lpzex nkapii nkaphh=',10i5)")i,lmxb(i),lpzex(i),nkapii(i),nkaphh(i)
          write(stdo,"('mmm rsmh1 ',i4,100f6.2)")i, rsmh1(1:lhh(1,i)+1,i)
          write(stdo,"('mmm   eh1 ',i4,100f6.2)")i,   eh1(1:lhh(1,i)+1,i)
          if(nkapii(i)==2) write(stdo,"('mmm rsmh2 ',i4,100f6.2)")i, rsmh2(1:lhh(2,i)+1,i)
          if(nkapii(i)==2) write(stdo,"('mmm  eh2  ',i4,100f6.2)")i,   eh2(1:lhh(2,i)+1,i)
          if(lpzex(i)==1 ) write(stdo,"('mmm pz    ',i4,100f6.2)")i,    pz(1:lhh(nkaph,i)+1,1,i)
          write(stdo,"('mmm lh    ',i4,100i3)")  lhh(1:nkaph,i)
         enddo
      endif
      
!! Atomic position Relaxation setup (MD mode)
      nitrlx= ctrl_nitmv  ! num of iteration cycle for atomic relaxiation (outer loop)
      mdprm = ctrl_mdprm  ! MD(relxation) condition setup
      defm  = ctrl_defm   !call dcopy(size(ctrl_defm),ctrl_defm,1,defm,1)
      if( nint(mdprm(1))==0) nitrlx=0 !no relaxiation. Only sc calculation for given atomic position.
      if( nint(mdprm(1))>0 .and. nint(mdprm(1))<4 ) call rx('lmf not set up for MD yet')
      if( nitrlx>0 ) then       !nitrlx >0 is for atomic position relaxiation
         allocate(indrx_iv(2,3*nbas))
         block  ! call rlxstp(natrlx, indrx_iv,  xyzfrz, pdim )
         integer:: i,j,k,ifrlx(3),lrlx 
         logical:: force,mdxx
         force = int(ctrl_lfrce) .gt. 0
         if (.not. force .or. nint(mdprm(1)) .eq. 0) goto 9299
         mdxx = nint(mdprm(1)) .le. 3
         lrlx = mod(nint(mdprm(1)),100)
         j = 0
         if (mdxx) then
            xyzfrz = .false.
            goto 9299
         elseif (force .and. mdprm(1) .ge. 100) then
            do  i = 1, 6
               if (defm(i) .eq. 1) then
                  j = j+1
                  indrx_iv(1,j) = i
               endif
            enddo
         elseif (force) then
            do  i = 1, nbas
               ifrlx=v_ssite(i)%relax
               do  k = 1, 3
                  if (ifrlx(k) .eq. 1) then
                     j = j + 1
                     indrx_iv(1,j) = k
                     indrx_iv(2,j) = i
                     xyzfrz(k) = .false.
                  endif
               enddo
            enddo
         endif
         natrlx = j
         if (natrlx .eq. 0) goto 9299
         pdim = 0
         if (.not. mdxx) then
            if (lrlx .eq. 4) pdim = natrlx*7
            if (lrlx .eq. 5) pdim = natrlx*(7+natrlx)
            if (lrlx .eq. 6) pdim = natrlx*(12+2*natrlx)
         endif
         if (iprint() .ge. 30) then
            if (lrlx .eq. 4) then
               write(stdo,*)' m_lmfinit RLXSTP: Molecular statics (conjugate gradients) ..'
            elseif (lrlx .eq. 5) then
               write(stdo,*)' m_lmfinit RLXSTP: Molecular statics (Fletcher-Powell) ..'
            else
               write(stdo,*)' m_lmfinit RLXSTP: Molecular statics (Broyden) ..'
            endif
            write(stdo,ftox)'   relaxing',natrlx,'variables,',nitrlx,'iterations'
            write(stdo,ftox)'   x-tol g-tol step=', ftom(mdprm(3:5))
         endif
 9299    continue
         end block
      endif
      call tcx('m_lmfinit')
      end subroutine m_lmfinit_init
      
!ssssssssssssssssssssss      
      subroutine getiout(a,iin,iout) !a(1:iout) can be nonzero.
      integer:: iin,iout,i
      real(8):: a(iin)
      do i=iin,1,-1
         if(a(i)>0) then
            iout=i-1 !angular mom
            exit
         endif
      enddo
      end
!ssssssssssssssssssssss      
      subroutine getiout10(a,iin,iout) !a(1:iout-1) can be >10
      real(8):: a(iin)
      integer:: iin,iout,i
      do i=iin,1,-1
         if(a(i)>10d0-1d-6) then
            iout=i-1 !angular mom
            exit
         endif
      enddo
      end

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

      end module m_lmfinit

      
!! Hereafter are list of variables. Old document, but it maybe a help.
!!mmmmmm old doc mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm      
! def      Uncertainty in Fermi level
! dosw     Energy window over which DOS accumulated
! fsmom    fixed-spin moment (fixed spin moment method)
! ef       Fermi level
! ntet     number of tetrahedra
!! --- explanation of  bz_* ---
Cr name    purpose
Cr efmax    eigenvalues above efmax are not needed
Cr lcond    conductivity, and directions
cCr lio      1  read qp
cCr          2  write qp
Cr lmet     0 insulator   1 3-point scheme  2 save evcs
Cr          3 two passes  4 wt from prior iter
Cr          (numbers, not powers of 2)
Cr lopt     2 BZMAP option
Cr lshft    shift mesh in each of three directions
Cr n        Polynomial order for Methfessel-Paxton sampling
Cr ndos     No. DOS points (sampling integration, and lmdos)
Cr nevmx    eigenvectors above nevmx are not needed
Cr nkabc    number of divisions in qp mesh
Cr nkp      Number of qp
Cr odos     offset to total dos
Cr oidtet   offset to tetrahedron corners
Cr oipq     offset to ipq array made in bzmesh
Cr opdos    offset to partial dos
Cr oqp      offset to qp array
Cr ostar    offset to array containing info on star of k
Cr owtkp    offset to array containing qp weights
Cr range    number of FWHM for sampling integration
Cr semsh    nz, modec, emin,emax, ecc,eps delta
Cr size     size of this structure
Cr stnr     parameters for Stoner model (not used)
Cr w        Line broadening for sampling integration
      
!! --- explanation of ctrl_* -----       
Cr   nclasp  (pgf) number of classes including PL -1,npl
Cr   oclabl  offset to class label array (r8tos8 format)
Cr   oics    Index to which species each class belongs
Cr   oipc    index to which class each site belongs
Cr   oipcp   index to which class each site belongs
Cr   oips    index to which species each site belongs
Cr   onrc    number of sites in each class
Cr   onrcp   number of sites in each class
Cr   ormax   sphere radius, by class
Cr   defm    Lines for deformation when relaxing plat
Cr   elin    ccor linearization energy, 2C hamiltonian
Cr   lbas    1 Hamiltonian has no screening transformation
Cr           2 Hamiltonian is nfp style
Cr           16 freeze phi,phidot for all species
Cr   lcd     1 freeze core
Cr           2 non-self-consistent Harris
Cr           4 represent full potential density on a uniform mesh
Cr           8 represent full potential density via TCF
Cr          16 unused
Cr          32 unused
Cr          64 (molecules) XC potential by FFT
C   ldos    1 make dos
cCxxxr           2 generate weights for partial dos
cCxxxr           4 generate weights for m-decompos'n of pdos
cCxxxr           8 BZ map
Cr   lfp     switches for full-potential
Cr           1 program uses full potential
Cr   lfrce   How forces are calculated
Cr           0 do not calculate forces
Cr           1 shift in FA density
Cr           2 shift in core+nucleus
Cr          10 added for shift to be screened
Cr   lmet    1 metal     2 tetrahedron
Cr   lncol   1 noncollinear magnetism
Cr           2 spin spirals
Cr           4 spin-orbit coupling
Cr           8 External magnetic field
Cr          16 mag. forces
Cr          32 spin-orbit coupling, LzSz only
Cr          64 spin-orbit coupling, LzSz + (L.S-LzSz) pert
Cr   loptc   1 generate Im(eps) (optics package)
Cr           2 generate Im(eps) w/ on-the-fly sampling integration
Cr          -1 generate joint DOS
Cr           Add 10 for SHG
Cr   lordn   1 Embedded GF  4 Vanderbilt  8 Embedded cluster
Cr   lpgf    layer GF.  First entry:
Cr           1: Generate diagonal GF and output density
Cr           2: Find Fermi left- and right-bulk Ef
Cr           3: Trace out band structure for spec'd energy mesh
Cr              In this case make Im(z)=0.
Cr           4: Ditto for right bulk.
Cr           5: Calculate current through structure
Cr              Second entry:
Cr           0  Make GF using standard approach
Cr           1  Use LU decomposition
Cr   noinv    1 do not add inversion !noinv BZ_NOINV
Cr   lrel    1 scalar relativistic
Cr           2 fully relativistic
Cr   lrs     switches concerning restart mode.
Cr           1 Read from restart file
Cr           2 Read from restart file, ascii mode
Cr           4 Read from restart file, invoke smshft
Cr           8 Write new density to restart file
Cr          16 Write new density to restart file, ascii format
Cr          32 read site positions from input file
Cr          64 read starting fermi level from input file
Cr         128 read starting pnu level from input file
Cr         256 rotate local density after reading
Cr   lscr    0 do nothing
Cr           1 Make P0(0)
Cr           2 Screen output q and ves
Cr           3 Screen output ves only
Cr           4 Use model response to screen output q
Cr             Add 1 to combine mode 1 with another mode
Cr             Add 10*k to compute intra-site contribution to
Cr             vbare each kth iteration
Cr             Add 100*k to compute response function only
Cr             each kth iteration
Cr   lstonr  second digit for graphical output
Cr   lstr    Note: no longer used; see str->lshow,str->lequiv
Cr           1 print strux   2 find equiv strux
Cr   lsx     1 Calculate screened exchange sigma
Cr          10's digit nonzero make vintra to include in sigma
Cr   ltb     1 overlap        2 crystal-field     4 ovlp+CF
Cr           8 add ebarLL    16 forces           32 fij
Cr          64 not used     128 pressure        256 evdisc
Cr         512 pair pot    1024 TrH & local E  2048 local rho
Cr        2^12 Hubbard U   2^13 No Madelung    2^14 wgt avg U
Cr        2^15 L>0 estat   2^16 disc read incr 2^17 gamma-pt
Cr   lves    1 take ves as input
Cr   lxcf    parameter defining XC functional
Cr           1s digit:
Cr           1 for Ceperly-Alder
Cr           2 for Barth-Hedin (ASW fit)
Cr           103 for PBE-GGA
Cr   maxit   max. no.  iterations in self-consistency cycle
Cr   mdprm   arg 1: 1 new dynamics  2  restart dynamics
Cr                  4 relax with conjugate gradients
Cr                  5 relax with variable metric
Cr                  6 relax with Broyden
Cr           arg 2: statics: switch
Cr                  1 read hessian matrix
Cr                  dynamics:
Cr                    number of iterations between printouts.
Cr           arg 3: (stat) relaxation x-tolerance
Cr                  (dyn)  temperature
Cr           arg 4: (stat) relaxation g-tolerance
Cr                  (dyn)  time step
Cr           arg 5: (stat) step length
Cr                  (dyn)  relaxation time
Cr           arg 6: (stat) Remove hessian after this many steps
Cr                  (dyn)  --
Cr   modep.. which dimensions are periodic
Cr   nbas    size of basis
Cr   nbasp   size of padded basis (layer geometry)
Cr   nclass  size of class
Cr   nitmv   max number of mol-dynamics iterations
Cr   nl      1 + maximum lmxa
Cr   nmap    number of maps (ASA)
Cr   npl     number of principal layers (PGF)
Cr   nbas   number of sites
Cr   nspec   number of species
Cr   nspin   number of spins
Cr   nvario  number of variables to output
Cr   omax1   sphere overlap constraints, type 1
Cr   omax2   sphere overlap constraints, type 2
Cr   pfloat  =1
Cr   rmaxes  upper limit to ES radius when finding new empty spheres
Cr   rmines  lower limit to ES radius when finding new empty spheres
Cr   sclwsr  scale wsr until reaching this fractional vol
Cr           10s digit used for assymetric treatment of ES:
Cr             0 ES and other sites are treated symmetrically
Cr             1 all sites with z>0 are resized first; then
Cr               all sites are resized.
Cr             2 all sites with z>0 are resized first; then
Cr               the ES sites only are resized
Cr   sdmod   spin dynamics mode
Crxxx   sdprm   spin dynamics parameters
Cr   sdxsi   Bulgac and Kusnezov global deamons
Cr   size    size of this structure
Cr   smalit  parameters for small iterations
Cr   tol     1 q- tolerance for self-consistency
Cr           2 e- tolerance for self-consistency
Cr   wsrmax  constraint on size of largest WS sphere
Cr   zbak    background charge and MT correction parameters
Cr                     zbak(1) = uniform background charge included
Cr                               in electrostatics but not in xc pot.
Cr                     zbak(2) = charge used to make MT correction     
Cr                               (ASA only) ==> unused!
C ----------------------------------------------------------------
!! ---------- ham_* -------------------
Cr  name    purpose
Cr  alfsi  stability factor alfsi*(S+)*S (molecules)
Cr  amgm   magnetization
Cr  bandw  Maximum size of off-diagonal, packed storage
Cr  dabc   Spacing for real-space mesh (molecules)
Cr  ehf    Harris-Foulkes energy
Cr  ehk    Hohnberg-Kohn energy
Cr  elind  Lindhard screening parameter
Cr  eterms terms making up the total energy.   For FP:
Cr          1  =  ehf
Cr          2  =  eks
Cr          3  =  utot
Cr          4  =  rhoves int (rhoval x Ves) (true-sm)
Cr          5  =  cpnves int (core+nucleus x Ves) (true-sm)
Cr          6  =  rhoeps int (rho * exc)
Cr          7  =  rhomu  int (rho * vxc)
Cr          8  =  sumec
Cr          9  =  sumtc
Cr          10 =  xcore  rhoc * total potential
Cr          11 =  valvef int (rhov * vef) (true-sm)
Cr          12 =  sumt0  sum of foca core energies
Cr          13 =  sumev  band structure energy
Cr          13 =  dq1
Cr          14 =  dq2
Cr          15 =  amom   system magnetic moment
Cr          16 =  sumev  sum of single-particle evals
Cr          17 =  rinvxt input density * external pot
Cr          18 =  rouvxt output density * external pot
Co          19 =  rhosig trace of self-energy over occ states
Cr  hord   order of polynomial approximation to ham.
Cr         In 2nd gen LMTO, relevant only in GF context.
Cr  kmto   envelope kinetic energies making up chord LMTO
Cr  lasa   not used; see ctrl lasa
Cr  ldham  vector describing hamiltonian dimensions:
Cr         1: ldim   = dimension of lmto basis
Cr         2: lidim  = ldim + size of downfolding block
Cr         3: lidhim = lidim + size of higher block
Cr         4: nspc   = number of coupled spins
Cr                   = 1 unless noncollinear magnetism
Cr         5: ldimc  = ldim * nspc
Cr         6: lidimc = lidim * nspc
Cr         7: lihdimc = lihdim * nspc
Cr         8: nspx   = number of separate spin channels
Cr                   = nsp if nspc is not 2; else 1
Cr  lmaxu   dimensioning parameter for LDA+U potential
Cr  lmxax   largest augmentation lmax in basis
Cr  lncol   1 noncollinear magnetism
Cr          2 spin spirals
Cr          4 spin-orbit coupling
Cr          8 External magnetic field
Cr         16 mag. forces
Cr         32 spin-orbit coupling, LzSz only
Cr         64 spin-orbit coupling, LzSz + (L.S-LzSz) pert
Cr         NB: ham->lncol and ctrl->lncol should be duplicates
Cr  lsig   parameters concerning self-energy
Cr          1 read sigma, Assume r.s. sigma is real
Cr          2 read sigma
Cr  ltb    not used; see ctrl->ltb
Cr  lxcf   parameter defining local XC functional.
Cr         Not used: see ctrl->lxcf
Cr  nbf    number of channels for magnetic field
Crxxxx  ndham  Largest dimension of lmto+PW basis
Cr  ndhrs  dimension of site block for r.s. h or sigma
Cr  ndofH  leading dimension to ooffH
Crxxx  neula  number of channels for euler angles
Cr         0 -> No euler angles defined
Cr  nkaph  number of repetitions of one l-quant.n. in basis
Cr  nlibu  Number of LDA+U blocks
c      
Cr  npwmin Estimate for minimum number of PWs to be calculated     <-Removable(probably)
Cr         (PW basis is q-dependent; max size not known a priori)  
Crxxx  npwpad Padding to be added to estimated max basis dimension    <-Removable(probably)
Crxxx         and subtracted from min basis dimension
Crxxx         (PW basis is q-dependent; max size not known a priori)
Cr  nqsig  Number of k-points for which self-energy known          <-Removable(probably)
Cr         (Used for interpolating between qpoints)
Crxxx  eula  Euler angles
Cr  hrs   r.s. h or sigma 
Cr  iaxs  neighbor table 
Cr  indxo orbital permutation table (pointer in w array)
Cr  lmxa  vector of augmentation l-cutoffs (pointer in w array)
Cr  ontabs table of number-of-pairs per site (pointer in w array)
Cr  oqsig  list of qp at which sigma can be computed
Cr  oveps  When diagonalizing hamiltonian, discard part of hibert space
Cr         corresponding to evals of overlap < oveps
Cr  pmax   global minimum allowed values for pnu
Cr  pmin   global minimum allowed values for pnu
Cr  pwemax High Energy cutoff for PW part of basis
Cr  pwemin Low Energy cutoff for PW part of basis
Cr  pwmode Controls PW part of basis
Cr         0 => no PW part of basis
Cr         1 => include PWs in basis
Cr         2 => include only PWs in basis
Cr  qpoff  qp offset when generating transformation of self-energy
Cr  qss    spin spiral
Cr  seref  Sum of reference energies
Cr  thrpv  3 PV for cell
Cr  udiag  diagonal-only LDA+U
c------------------------------------------------------------------------------
      
!! --- explanation of  lat_* ---
Cr name    purpose
Cr  alat    lattice parameter, in a.u.
Cr  as      dimensionless Ewald smoothing parameter
Cr  avw     average MT radius
Cr  awald   Ewald smoothing parameter
Cr  dist    deformation parameters (cf lattdf.f)
Cr  gam     lattice shear parms: gam, gx,gy,gz
Cr  gmax    cutoff gmax for Fourier transforms
Cr  ldist   switch specifying what kind of dist (cf lattdf.f)
Cr  nabc    no. divisions for F.T. mesh
Cr  ng      no. G vectors
Cr  nkd     no. direct latt. vecs. for Ewald sum
Cr  nkdmx   dimensioning for arrays holding latt. vecs
Cr  nkq     no. reciprocal latt. vecs. for Ewald sum
Cr  nkqmx   dimensioning for arrays holding latt. vecs
Cr  npgrp   Number of point symmetry group operations
Cr  nsgrp   Number of space symmetry group operations
Cr  oag     offset to symmetry group translations
Cr  obgv    phase factor sum for symmetrization of mesh rho
Cr  ocg     offset to Clebsch Gordan coeffs
Cr  ocy     offset to Ylm normalization constants
Cr  odlv    offset to direct lattice vector
Cr  ogv     offset to list of F.T. G vectors
Cr  oidxcg  offset to Clebsch Gordan indxcg
Cr  oips0   pointer to first vec in star (for symm mesh rho)
Cr  oistab  offset to site permutations table for group ops
Cr  ojcg    offset to Clebsch Gordan jcg
Cr  okv     offset to indices in list of F.T. G vectors
Cr  opos    offset to site positions
Cr  oqlv    offset to Ewald reciprocal lattice vectors
Cr  osymgr  offset to symmetry group rotation matrices
Cr  plat..  lattice vectors, units of alat
Cr  plat0.. lattice vectors before distortion
Cr  plat2.. secondary lattice vecs used in various contexts
Cr  plate.. order-N
Cr  platl.. pgf
Cr  platr.. pgf
Cr  qlat..  reciprocal lattice vectors, units 2pi/a
Cr  rpad..  truncate Ewald to rpad*rmax when lattice vector
Cr          list has to be padded in order to include at
Cr          least one lattice vector
Cr  size    size of this structure
Cr  slat    superlattice vectors
Cr  tol     Ewald tolerance
Cr  tolft   FT mesh tolerance
Cr  vol     cell volume

!! --- explanation of  umix ---
Cr   name    purpose
Cr   b       mixing beta
Cr   bl      previous mixing beta
Cr   bv      extra potential mixing
C    elind   Lindhard energy for model screening
Cr   fn      mixing file name
Cr   kill    kill the mixing file after k iterations
Cr   lxpot   decouple potential and the charge
Cr           1: mix vin and v(qmix)
Cr           2: mix vin and v(qout)
Cr   mmix    maximum number to mix
Cr   mode    1 Anderson 2 Broyden
Cr   model   previous mixing mode
Cr   n       Number of iterations for this species
Cr   nitu    max number of LDA+U itreations
Cr   nmix    actual number mixed
Cr   nsave   number of iterations to save on disk
Cr   r..     expression for rmscst
Cr   rms1    1st rms error
Cr   rms2    2nd rms error
Cr   size    total size of this struc
Cr   tj..    Anderson t's
Cr   tolu    tolerance for LDA+U
Cr   umix    mixing parameter for LDA+U
Cr   w..     Linear mixing weights
Cr   wc      Broyden weight

!! --- explanation of  move ---
Cr   off offe  name    purpose
Cr    2    1   ct      coefficients for global deamons thermostat
Cr    5    2   gyro    gyromagnetic ratio (magnetic dynamics)
Cr    6    3   kt      temperature, atomic units
Cr    7    4   modt    mode for thermostat
Cr   10    5   nmodt   number of thermostat modes
Cr   11    6   prmint  Parameters for numerical integration
Cr                     For Bulirsch-Stoer integration:
Cr                     1:   mode: (1 for BS-integration)
Cr                     2:   rst:  1 for start new integration, 0 to continue
Cr                     3:   ts0:  minimum time step size
Cr                     4:   tol:  tolerance in integration errors
Cr                     5:   mx:   order of rational function extrapolation
Cr                     6:   mi:   number of midpoint rules
Cr                     7-17 nseq: sequence of no. midpoint divisions
Cr                     18:        offset to bs workspace
Cr   32    8   tnow    duration of simulation to this point, in units of ts
Cr   33    9   ts      suggested time step size
Cr   34   10   tsequ   initial equilibration time before stats accumulated
Cr   35   11   tstot   duration of total simulation, in units of ts

!! --- explanation of  pot_* ---
Cr  bfield  global magnetic field direction and amplitude (4)
Cr  nlma    number of augmentation channels in system
Cr  nlml    (FP) total number of site density channels
Cr  nrhos   (FP) total number of site density channels
Cr  oaamom  offset to array containing local ASA mag.mom
Cr  obxc    offset to array containing local ASA Bxc dir.
Cr  odddpf  3rd derivative of potential functions (mkptfp)
Cr  oddpf   2nd derivative of potential functions (mkptfp)
Cr  oddpfr   2nd derivative of potential functions (mkptfp),
Cr          fully relativistic case
Cr  odel    (tbe) potential shifts.
Cr  odpf    energy derivative of potential functions (mkptfp)
Cr  odpfr   energy derivative of potential functions (mkptfp),
Cr          fully relativistic case
Cr  ofes    electrostatic contribution to forces
Cr  ogma    gamma-alpha in potential functions
Cr  ogmar   (gamma-alpha)*P^alpha/P^gamma fully relativistic
Cr  ogrrme  radial matrix elements of gradient.
Cr  ohab    not used now
Cr  opalp   ASA potential functions, alpha repsn (mkptfp)
Cr  opapg   p^alpha/p^gamma fully relativistic
Cr  opdel   ASA Delta, in hamiltonian order (for LDA+U)
Cr  opf     ASA potential functions (mkptfp)
Cr  opfnc   noncollinear ASA potential functions
Cr  opfr    ASA potential functions (mkptfp),
Cr          fully relativistic case
Cr  opmpol  multipole integrals of phi,phidot
Cr  opnu    P-nu (Methfessel's log derivative function)
Cr  opp     potential parameters
Cr  oppi    not used now
Cr  oppn    NMTO generation potential parameters
Cr  opprel  Dirac potential parameters
Cr  opti    inverse potential functions
Cr  oqc     Sphere core charge
Cr  oqmom   offset to multipole moments
Cr  oqnu    Energy moments of the charge density
Cr  oqpp    Multipole moments of the nonspherical charge
Cr  oqt     Sphere total charge
Cr  orhos   spin-density matrix
Cr  orhrmx  sphere electron density at rmax
Cr  osab    not used now
Cr  osgw    structure containing GW parameters
Cr  osmpot  offset to smoothed potential
Cr  osmrho  offset to smoothed density
Cr  osop    spin-orbit parameters
Cr  osoptc  structure containing optical matrix elements
Cr  osrout  offset to smoothed output density
Cr  otau    not used now
Cr  ovab    not used now
Cr  ovdif   difference ves(q,rmax) - vactual(rmax)
Cr  oves    electrostatic potential at rmax
Cr  ovintr  intra-atomic W for screening, monopole approx
Cr  ovrmax  total potential at rmax
Cr  ovshf   ASA constant potential shifts
Cr  size    size of this structure
Cr  vconst  Constant estat potential shifts, used where
Cr          the Fermi level is specified and the potential
Cr          adjusts to it.
Cr          vconst(1) = potential shift
Cr          vconst(2) = potential shift of L end region (PGF)
Cr          vconst(3) = potential shift of R end region (PGF)
!!mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm      
      

         
