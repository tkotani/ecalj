!! = Maine part of full-potential LDA/GGA/QSGW(for given sigm). Single iteration =
!!
!! Memo july2012takao
!! 1)lwtkb algolism (double path method for METAL=3, lwtkb=-1 is too complicated.)
!!   This makes this code too complicated.
!!   (second path formalism; 1st path for weight, second path for integration for density)
!!   "goto 99" mechanism should be avoided.
!!
!! 2) nevmx is given correctly ---> we need to examine a mechanism to determine nevmx in call subzi.
!!   It will be effective to give nevmx to reduce computational time-->nmx; we need to 
!!

c#define KINODEBUG
#define USEOPTIONAL 1
      subroutine bndfp(nbas,nsp,nlibu,lmaxu,lldau, ssite,sspec,slat,sctrl,sham,spot,sbz,sstrnmix,
     .  ndham,leks,lrout,lfrce,lpnu,dmxp,iter,maxit, frc,
     .  dmatu,vorb,llmfgw)
      use m_hamindex, only: napwmx,ngrp,norbmto,ibastab,ltab,ktab,offl, symops_af

c     &  ngrp, norbmto, napwmx, ltab,ktab,offl,ispec,iclasst,invgx,miat,offlrev,ibastab,
c     &  igv2,napwk,igv2rev
      use m_rdctrl, only: ncutovl,lso
      use m_globalvariables, only: globalvariables
      use m_phmbls,only :phmbls
c      use m_seneinput,only: npair,nlat,qwgt,npairmx,qsmesh, Seneinterp2 !these are given here.
c      use m_seneinput,only: Seneinterp2 !these are given here.
      use m_subzi, only: subzi
      use m_rsibl_ev, only : rsibl_ev ! to plot wavefunction in the fat band mode
      use m_w_psir, only: w_psir_cmdopt
      use m_struc_def
#if defined(__INTEL_COMPILER)
      use ifport, only: unlink
#endif
      use m_ext,only: sname !extention for file
!! ==== followings are taken from origitanl lm7beta for reference. ====
Ci   nbas  :size of basis
Ci   nsp   : number of spins
Ci   nlibu : total number of LDA+U blocks (used to dimension dmatu and vorb)
Ci   lmaxu : lmax for U used to dimension vorb and dmatu
Ci   lldau :lldau(ib)=0 => no U on this site otherwise
Ci         :U on site ib with dmat in dmats(*,lldau(ib))
Ci
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: *
Ci     Sored:    *
Ci     Passed to: dfaugm mkpot dfqkkl suclst surho sumlst hambl mullmf
Ci                makusq mkpdos addrbl mshn3p mchan vcdmel dfratm
Ci                mkrout symrho pnunew dfrce mkekin totfrc mixrho
Ci
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: *
Ci     Stored:    *
Ci     Passed to: dfaugm mkpot dfqkkl suclst surho sumlst hambl mullmf
Ci                makusq mkpdos addrbl mshn3p mchan vcdmel dfratm
Ci                mkrout symrho pnunew dfrce mkekin mixrho
Ci
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: plat qlat nabc vol nsgrp
Ci     Stored:    *
Ci     Passed to: mkpot hambl makusq addrbl vcdmel mkrout symrho dfrce
Ci                mkekin totfrc mixrho
Ci
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read: ldos,1 zbak nl
Ci     Stored:    *
Ci     Passed to: vcdmel dfrce
Ci
Ci   sham  :struct for parameters defining hamiltonian; see routine uham
Ci     Elts read: elind oindxo
Ci     Stored:    eterms ehf ehk
Ci     Passed to: mkpot hambl makusq addrbl mkrout mkehkf mkekin
Ci
Ci   spot  :struct for information about the potential; see routine upot
Ci     Elts read: oorhat osmrho osmpot osoptc nlml nlma
Co     Stored:    oorhat osmrho osmpot
Ci     Passed to: mkpot
Ci
Ci   sbz   :struct for the Brillouin Zone; see routine ubz
Ci     Elts read: nkabc nkp ntet oidtet lmet n w nevmx efmax fsmom ndos
Ci                dosw ef def oqp owtkp
Ci     Stored:    ndos dosw ef def
Ci     Passed to: *
Ci
Ci   sstrn :struct for global strings
Ci     Elts read: mix
Ci     Stored:    *
Ci     Passed to: *
Ci
Ci   ndham :dimensioning parameter, at least as large as largest
Ci         :hamiltonian dimension
Ci   leks  :>0 make the Hohnberg-Kohn-Sham energy
Ci         :>1 use the HKS forces
Ci   lrout :>0 generate output density and attendant quantities
Ci   lfrce : 0 suppress generation of forces
Ci   lpnu  : 1 make new pnu's
Ci   dmxp  :vector of mixing parameters; see mixrho.f for dmxp(1..25)
Ci         :Additionally:
Ci         :dmxp(33)  is the Lindhard parameter
Ci   iter  :current iteration number
Ci   maxit :maximum number of iterations
Ci   evl   :work array holding band eigenvalues for current spin, qp
Cio LDA+U inputs and outputs
Cio  dmatu :density matrix for LDA+U (changed upon output)
Cio  vorb  :orbital dependent LDA+U potential
Co Outputs
Co   frc   :forces.  Only calculated if lfrce>0.
Co         :If leks<2, forces are HF  forces
Co         :If leks>1, forces are HKS forces
Cl Local variables
Cl   k1,k2,k3: dimensions smrho,smpot
Cl   lpdiag:0 use standard diagonalization (zhev)
Cl         :1 use parallel diagonalization (pzhev)
Cl         :2 diagonalization done internally
Cl   lwndow:T if to make density in a specified energy window
Cl   lcplxp:0 if ppi is real; 1 if ppi is complex
Cl   jsp   :current spin index.
Cl         :In the collinear case, jsp and isp are equivalent
Cl         :In the noncollinear case, isp loops 1..2 for the
Cl         :purpose of assembling the hamiltonian.
Cl         :Once assembled, isp should not be used; and jsp=1
Cl   ispc  :2 when working on (2,2) block of noncollinear hamiltonian;
Cl         :otherwise 1
Cl   nspx: number of independent spin channels containing eigenvalues
Cl         and total DOS; nspx=nsp unless nspc=2, in which case nspx=1
Cl   onesp :do only one spin branch of isp loop (spec'd by onesp)
Cl         :also used when usual loop order (iq=1..nq, isp=1..2)
Cl         :needs to be reversed, as it does, e.g. when transforming
Cl         :sigma matrix.  Then onesp plays the role of spin index
Cl   lekkl :0 do not accumulate oeqkkl; 1 do accumulate oeqkkl
Cr Remarks
Cr   Band pass consists of:
Cr   (1) make the effective potential,
Cr   (2) generate eigenvalues (and eigenvectors if lrout)
Cr   (3) if lrout, assemble the output density by BZ integration
Cr   (4) evaluate hf (and KS, if leks) energy by BZ integration
Cr   (5) mix the output density to make a new input density.

Cu Updates before 2009.
Cu   05 Jul 08 (T. Kotani) new PW basis
Cu             Option to accumulate energy-weighted output density
Cu   27 Jun 08 Redesigned transformation of sigma to new basis
Cu   09 Jul 07 MPIK enabled to plot bands
Cu   05 Jul 07 Enable onesp to be set as switch in --band:spin1
Cu   09 Jun 07 Fixed-spin-moment, noncollinear case
Cu   16 Jan 07 First cut at I/O of sigm transformed to,from LDA basis
Cu   26 May 07 Some preparation for rotation betw/ LDA, GW basis
Cu   17 Jul 06 Some MPI changes, fixes SO case
Cu   27 Jun 06 New constraints for floating pnu
Cu   08 Jun 06 Bug fix total DOS, noncollinear case;
Cu             Mulliken works for noncollinear case
Cu   02 Jan 06 better printout of magnetic moment
Cu   27 Nov 05 LDA+U => complex potential
Cu   09 Nov 05 (wrl) Convert dmat to complex form
Cu   06 Oct 05 (A. Chantis) bug fix dos when nspc=2
Cu   25 Jul 05 bug fix partial dos combined with LDA+U
Cu   01 Jul 05 handle sites with lmxa=-1 -> no augmentation
Cu   29 Jun 05 (MvS) extended LDA+U to local orbitals
Cu   27 Apr 05 LDA+U (Lambrecht)
Cu   14 Feb 05 fixes for band plot, contour mode
Cu   03 Feb 05 (A. Chantis) implemented spin-orbit coupling by L.S
Cu   11 Jan 05 double-counting term rho*sig subtracted from ehks.
Cu   23 Dec 04 Extended to spin-coupled case
Cu   18 Nov 04 Sampling integration properly handles Fermi distribtion
Cu   25 Sep 04 (ATP) some patches for MPI parallelization
Cu    1 Sep 04 Adapted to handle complex ppi. Spin-orbit folded into ppi
Cu   12 Aug 04 First implementation of extended local orbitals
Cu   29 Jun 04 (A Chantis) First implementation of spin-orbit coupling
Cu             (Lz.Sz only)
Cu   19 Sep 03 (ATP) Modifications of CLS spectroscopy
Cu   24 May 03 New --window switch for accumulating density
Cu             in a specific energy window
Cu   24 May 03 New interpolation mode for sigma
Cu   14 Aug 02 Added file input of self-energy addition to LDA
Cu   15 Feb 02 (ATP) Added MPI parallelization
Cu   24 Oct 01 Updated mixrho, dfrce
Cu   24 Aug 01 Extended to local orbitals.
Cu   22 Apr 01 Added driver for Kotani's GW
Cu   21 Mar 01 bug fix in call to makdos
Cu   20 Mar 01 (ATP) Added Mulliken analysis, CLS
Cu   15 Feb 01 eliminated smrho, smpot from passed arguments
Cu   23 Jan 01 bug fixes connected with lrout=0
!! -------------------------------------------------------------------------
      implicit none
      logical :: llmfgw
#if MPI | MPIK
      include "mpif.h"
      integer ierr, status(MPI_STATUS_SIZE)
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*20 ext
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
      double precision starttime, endtime
c#if MPI
c      integer dims(2)
c#endif
#endif
c#if MPIK
      integer obuf,oqptr
      integer igetss,lmxh,nlmh,nelt(3),nkaph,kmax !,nlma
      integer, dimension(:),allocatable :: kpproc
      double precision sttime,entime
c#endif
c#if MPIK
      integer:: pidorigin, ndimhx_iq, nev_iq,idx
      integer:: work_mpi(10)    ! use only 1 and 2 now
c#endif
      logical mlog
      integer procid,master
      integer nbx,nbas,n0,nab,nppn,iter,maxit,ndham
c      character*(*) sstrn
      parameter ( nbx=512, n0=10, nppn=12, nab=9)
      integer k1,k2,k3,leks,lrout,lfrce,lpnu
      integer nlibu,lmaxu,lldau(nbas)
      real(8):: dmxp(33),frc(3,*)
c      type(s_sstrn):: sstrn
      character*(*):: sstrnmix
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)
      type(s_ham):: sham
      type(s_ctrl):: sctrl
      type(s_lat)::  slat
      type(s_pot)::  spot
      type(s_bz)  ::  sbz
      logical lgors,ltet,cmdopt,lwndow
      character strn*120,plbopt*120,mulopt*120,clsopt*120,strn2*120,dc*1

ckino lswtk must be 0 initially.
      integer:: ipr , iprint , ipl , isp , jsp , ispc , isum , nglob 
     .  , i , iq , ismidb ,  ldos , lmet , nspc  
     .  , lrep , lcplxp , nkabc(3) , ndos , ndos0 , nev , nevl 
     .  , nkp , nsp , numq , n1 , n2 , n3 , ntet , lwtkb , lswtk=0 , nl 
     .  , lgunit , igets , mpsord , nevmx , nvl , ngabc(3) , isw , stdo 
     .  , stdl , ifi , fxst , fopn , fopna , lfrzw , i1 , i2 , plbnd 
     .  , nfilem , iobzwt , lnoxc , lrsig ,  jobgw , nll !, lpdiag lwsig ,
     .  , iv(10) , parg , ndhamx , nspx , nk1 , nk2 , nk3 , lshft(3) 
     .  , onesp , mpipid , ldham(16) , ldim , ndimh , ndimhx !nfbn(2) , 
     .  , lekkl
      real(8) ,allocatable :: dos_rv(:,:)
      integer:: os , osoptc
      real(8),allocatable :: rv_a_owtkb(:,:,:)
      real(8),allocatable :: rv_a_oswtk(:)
      real(8) ,allocatable :: fh_rv(:)
      real(8) ,allocatable :: fes1_rv(:)
      real(8) ,allocatable :: fes2_rv(:)
      real(8) ,allocatable :: gpot0_rv(:)
      real(8) ,allocatable :: hab_rv(:)
      real(8) ,allocatable :: hbyl_rv(:,:,:)
      real(8) ,allocatable :: ppnl_rv(:,:,:,:)
      real(8) ,allocatable :: qbyl_rv(:,:,:)
      real(8) ,allocatable :: qmom_rv(:)
      real(8) ,allocatable :: sab_rv(:)
      complex(8) ,allocatable :: srout_zv(:)
      complex(8) ,allocatable :: t_zv(:,:)
      real(8) ,allocatable :: vab_rv(:)
      real(8) ,allocatable :: qm1_rv(:)
      real(8) ,allocatable :: qm2_rv(:)
      real(8) ,allocatable :: vval_rv(:)
      real(8) ,allocatable :: ww_rv(:)
      complex(8) ,allocatable :: ww_zv(:)
      integer ,allocatable :: ips_iv(:)
      real(8) ,allocatable :: orbtm_rv(:)
      integer ,allocatable :: ifbls_iv(:)
      type(s_rv1) :: sv_p_orhat1(3,nbx)
      type(s_rv1) :: sv_p_oppi(3,nbx)
      type(s_rv1) :: sv_p_oeqkkl(3,nbx)
      type(s_rv1) :: sv_p_oqkkl(3,nbx)
      type(s_rv1) :: sv_p_otau(3,nbx)
      type(s_rv1) :: sv_p_osig(3,nbx)
c      equivalence (ldim,ldham(1,1))
      logical a2bin
      integer:: ospotx , osgw
      type(s_rv1) :: sv_p_oppix(3,nbx)
      type(s_rv1) :: sv_p_otaux(3,nbx)
      type(s_rv1) :: sv_p_osigx(3,nbx)
      logical:: l_temp_rv_a_oswtk=.false., l_temp_rv_a_owtkb=.false.,
     &  l_temp_iv_a_oidtet=.false.
      double precision sigp(10)
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      double precision ef00,eferm,emax,emin,esmear,qval,qsc,sev,
     .  sev00,sev1,sumtv,alfa,vconst,qbg,ebot,def,ef0,ehar,eks,
     .  dosw(2),efmax,dgets,alat,plat(3,3),qlat(3,3),sumev(2,3),
     .  sumqv(3,2),qp(3),shftqp(3),eterms(20),elind,vol,dum,del,evtop,
     .  ecbot,vrmt(nbx),fsmom,rsrnge,xv(20),dosrng,epsovl,dval !,eomin
      equivalence (emin,dosw(1)),(emax,dosw(2))
!!    PW basis
      integer,allocatable:: igv2x(:,:)
      integer ,allocatable :: kv_iv(:)
      integer pwmode,napw
      double precision pwemin,pwemax,pwgmin,pwgmax
C     real(8),allocatable :: qpe(:,:)
      logical T, F
!!     for spin-orbit
      complex(8) ,allocatable :: auso_zv(:)
!!     for partial dos
      complex(8) ,allocatable :: ausp_zv(:)
!!     for CLS
      integer nsitmx
      parameter (nsitmx = 256)
      integer:: nlmax , i1mach , icls , isite(nsitmx) , iclsl(nsitmx) 
     .  , iclsn(nsitmx)
      complex(8) ,allocatable :: ausc_zv(:)
!!     for Mulliken
      integer:: moddos , nsites , lsites(nbx) , nchan , ng , iomoms 
     .  , nchmx , lmdim , lmxch
      logical:: mullwrite=.false.

ckino Dec.8.2011:         integer,pointer :: iv_p_oiprmb(:) =>NULL()
      real(8) ,allocatable :: doswt_rv(:)
      integer ,allocatable :: chan_iv(:),chan0_iv(:)
!!     for pzhev
      integer nblk,nprow,npcol,nmx
!!     for LDA+U
      real(8) ,allocatable :: ausu_rv(:)
      double complex vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu),
     .  dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
      integer:: i_spackv, i_spacks
      integer ::iwdummy
      logical :: debug=.false.  !debug is only for check write. However, it can cause error
!for 'make check',
!because fp/test/test.fp assumes 'something not inserted in console output'.

C For self-energy
      integer nqsig
      integer:: nbasp
      data T /.true./ F /.false./
      integer,allocatable:: nevls(:,:)
ctakao
      integer :: ifiogw,ifiese,oag,og,ib,ig,iorb,is,nsgrp,uspec_size,nspec
      logical :: nexist,sigmamode,ok !, newsigmasw
      character spid*8
      integer :: ib1,ib2,ii,ik1,ik2,ik3,j1,j2,ikx
      integer:: nhmax , nknknk , jx , nsp_dummy , ndimh_dummy
      complex(8),allocatable:: sfz(:,:,:,:,:,:),wgtq(:,:,:,:,:),sfzx(:,:,:,:,:,:),sfzmto(:,:,:,:,:,:)
      complex(8),allocatable:: hamm(:,:,:),ovlm(:,:,:),sene(:,:) ,senex(:,:),ovlmtoi(:,:),ovliovl(:,:) 
      integer:: isize_hamm(3)
      real(8):: rrr(3),rrrr
      complex(8):: pi2i
      real(8),allocatable:: range(:,:),pos(:,:)
      integer:: j1x,j2x,ix1,ix2,ix3 !ikc(3,0:1,0:1,0:1), igshift(3,0:1,0:1,0:1),
      integer,allocatable::igx(:,:,:,:)
      real(8):: wgc(0:1,0:1,0:1),qpara(3)
      real(8):: www,qmm(3),qqq(3),sumtvxxx,sevxxx
      integer:: iktm,ikt,igvecc(3),nkk(3),igvecc0(3), nnspc,osmpotbk
      logical::bzwtsf_bisection !noshorbz,
      real(8),allocatable:: qplist(:,:)
      real(8):: qpx(3)          !qpo(3),

      integer:: iqini,iqend
      complex(8),allocatable:: spotx(:,:,:,:),smrhobk(:)
      real(8),allocatable:: evlall(:,:,:),evl(:,:),dwgtall(:,:,:,:,:)
      integer:: idummy
      real(8):: epsnevec
      logical:: onlyonce=.true.
      logical:: isanrg, l_dummy_isanrg,oncewrite, savez,getz,gotosecondpass,initiqisp
      integer:: getef, i_dummy_getef,ifig, ndimsig
c      integer:: n_listwf
      integer, allocatable:: listwf(:)
      integer:: ntqxx,ispx,ix99 !,nevec unused.
      real(8)::eseavrmean
      logical:: mtosigmaonly
      real(8)::scaledsigma
      real(8):: qsfz(3,100)
      integer:: nqps,ii1,ii2,ii4,iq1 !, fxsts
      logical:: lfx

#if                             !defined(__INTEL_COMPILER)
      integer:: unlink
#endif

!! for mpi and mpik, parallel write
      character(10):: strprocid 
      integer:: nsizeohrs,nsizeoiaxs,nsizeontabs,noqsig
      integer:: ifis2,ipts,ifid
      real(8):: eee
      logical :: tdos !,pdos
      integer:: ifves,ifvesintatm,ifvesintloc,ispec,ifvbm,ifqbyl,il
      real(8):: vessm,sumvesloc,sumvesatm,vref,vesloc,evbm,ecbm
      real(8),allocatable::vesatm(:),qrmt(:,:),ql(:,:),qset(:,:),qatot(:)

      real(8):: rydberg=13.6058d0,evll,qrmtx,qsetsum,qlx(0:100),esold, vnow
      logical:: lfill=.false.,ixx,master_mpi,master_mpiset
      integer,allocatable:: ndimhx_(:),nev_(:)
      character(120):: vbmlll
      integer:: npri,ibas,ifiqb,ifiwv,lmxa,numprocs, ifimag
      real(8),parameter::    NULLR =-99999

      integer,parameter:: noutmx=48
      integer:: iout,nout,nlatout(3,noutmx),iapw
      real(8)::ppin(3)
      real(8):: rlatp(3,3),xmx2(3),vadd

      logical:: fsmode,PROCARon,allband !for --fermisurface for xcrysden.
      integer:: nkk1,nkk2,nkk3  !for fermi surface mode for xcrysden.
      integer:: ififm
      real(8):: qshort(3)
! add --fermisurface #1 #2 #3
! only implemented for lmf-MPIK now.
      integer:: ini,iend,nbxx,ival

c mar2014
      integer:: lb1,lb2,onespdummy
      integer:: iprocar,iprocar1,iprocar2,lrsigx
      complex(8),allocatable:: auspp(:,:,:,:,:)
      integer:: l,m,ilm,nlma,iband,hreal,mode,nttabs
      complex(8):: au,as,az
      real(8):: s11,s22,s33,s12,s13,s23,dwgt(100),dwgtt(100),xdat,qold(3)
      character*256:: sss,s222
      character*1000::ccc
      logical:: initbd=.true.

      integer:: nspsigm,ispsigm !nov2014
      logical:: spinoffdiag,    goto99!jan2015
      real(8),allocatable:: xdatt(:)
      logical:: fullmesh

      integer,allocatable:: ipqe(:,:,:),idtete(:,:)
      integer:: ntete,init,ient,nchanp,fopnn,nevmin,ifip,ifisyml,isyml
      real(8):: qx(3)
      real(8),allocatable:: wkd(:), pdosall(:,:,:)
      character*3::charnum3
      character*100::filenm(2)
      integer::ifbndsp(2)
      real(8):: eminp,emaxp
!! direct read syml
      integer,parameter:: nsymlmax=100
      integer:: nqp_syml(nsymlmax),nqps_syml(nsymlmax),nqpe_syml(nsymlmax)
      integer:: nqp2_syml(nsymlmax),nqp2s_syml(nsymlmax),nqp2e_syml(nsymlmax),nqp2n_syml(nsymlmax)
      logical:: masslineon(nsymlmax)

      integer:: nsyml=0,fopno,ikp,ifbnd,nsymln
      real(8):: qps_syml(3,nsymlmax), qpe_syml(3,nsymlmax),rq
      character*20 ::labeli(nsymlmax),labele(nsymlmax)
      real(8)::  totxdatt, dqsyml(nsymlmax),etolv,etolc
      integer:: ifile_handle, fxsts

      character(512),save::schar
      character(512):: aaachar
      logical:: massline , writeham=.false.
      logical,save:: siginit=.true.     ,laf
      character(10):: i2char
      character(50)::infoq
      integer::ipos,ndos_,j,ififft
      real(8)::dosmax,emaxp_,posp(3)
      integer:: jspend,jspini,jspp
      integer:: iqread,iqindex,ibread(100)=0,ibindex,nibread, job,ifih,ist,ip,ni,ix

      integer,allocatable:: npair(:,:), nlat(:,:,:,:),ib_table(:),k_table(:),l_table(:)
      integer::npairmx,ikpd,ifqplist
      real(8),allocatable:: qsmesh2(:,:,:,:)
      integer,allocatable:: nqwgt(:,:,:)
      real(8):: tpioa,pi=4d0*datan(1d0)
      complex(8):: aaaa,img=(0d0,1d0)
      character*4:: cccx

      real(8):: pnu(n0,2),pnz(n0,2),pmean
      integer:: i_copy_size
      logical:: lpz,phispinsym
!! ----------------------------
      call tcn ('bndfp')
      stdo = lgunit(1)
      stdl = lgunit(2)
      debug=cmdopt('--debugbndfp',6,0,strn)
      procid = mpipid(1) !zero for no mpi
      master = 0         !mpi master
      mlog = cmdopt('--mlog',6,0,strn) !--mlog here is taken by getarg.
      tdos = cmdopt('--tdos',6,0,strn)
      if(cmdopt('--band~col',10,0,strn)) then
        call rx(' no color mode is supported. Use PROCAR mode (--mkprocar) instead!')
      endif

#if MPI | MPIK
      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
      print *,' mpipid=',procid
#endif      

c$$$!! A trick to get the extension of ctrlfile. It is stored to trim(sname).
c$$$      if(procid==master.and.initex) then 
c$$$        ifi=fopn('temptempxxx') 
c$$$        sss=''
c$$$        INQUIRE (ifi, NAME=sss) !second run, 'sss' must contain extention, eg, '.copt'
c$$$        call fclose(ifi)
c$$$        ipos= scan(sss,'.',back=.true.) !backward search of '.' from the last.
c$$$        sname=sss(ipos+1:256)  !sname contains extention such as 'copt'.
c$$$        ifi=ifile_handle()
c$$$        open(ifi, file='temptempxxx.'//trim(sname)) 
c$$$        close(ifi,status='delete')
c$$$        initex=.false.
c$$$      endif
c$$$      if(procid==master) write(6,*) "bndfp: ctrl extension is +++"//trim(sname)//"+++"
c$$$      call mpibcc(sname,512,mlog,'bndfp','ext')

!! Write atomic positions
      if(procid==master) then
         write(stdo,"(a)")' --- bndfp: Used Atomic Positions from rst (or ctrl) ----'
         do ib=1,nbas
            is=ssite(ib)%spec
c            call r8tos8(ssite(is)%clabel,spid)
            spid=ssite(is)%clabel
            write(stdo,"('ib clabl pos =',i5,a6,3(x,f13.5))") ib, trim(spid),ssite(ib)%pos
         enddo
      endif


!     ! Write atomic positions
      print *,' procid=',procid,master
      if(procid==master) then
        write(stdo,"(a)")'== bndfp: Atomic Positions (from rst/ctrl or changing when dynamics ==)'
        do ib=1,nbas
          is=ssite(ib)%spec
          spid=ssite(is)%clabel
c          call r8tos8(ssite(is)%clabel,spid)
          write(stdo,"('  ib clabl pos =',i5,x,a6,3(x,f12.6))") ib, trim(spid),ssite(ib)%pos
        enddo
      endif

!! === New pdos mode (use --mkprocar and --fullmesh together). ===
!! We use all k points (--fullmesh), instead of using crystal symmetry.
!! pdos=cmdopt('--pdos',6,0,strn)  ! old pdos mode switch.
      PROCARon= cmdopt('--mkprocar',10,0,strn)
      fullmesh= cmdopt('--fullmesh',10,0,strn) !full mesh stop just after do 2010 iq loop.
      if(fullmesh.and.(.not.procaron)) then
         call rx(' --fullmesh is meaningful only with --mkprocar')
      endif
      if( cmdopt('--writepdos',11,0,strn) ) then
         write(6,*) '... Doing writepdos mode. Wait a while ...'
         write(6,*) '... Need "pdosinput" which is generated by --mkprocar & --fullmesh'
c         write(6,*) '... If necessary, we may make a subrouine writepdos() pararellized...'
         call writepdos(trim(sname))
         call rx0('done: end of --writepdos mode.')
      endif
   
!! === Fermi surfece mode for xcrysden. It's format is in  http://www.xcrysden.org/doc/XSF.html#2l.16
      fsmode=cmdopt('--fermisurface',14,0,strn)
      allband=cmdopt('--allband',9,0,strn)
      if(fsmode) fullmesh=.true.
      if(fullmesh) then
         nkk1=sbz%nkabc(1)
         nkk2=sbz%nkabc(2)
         nkk3=sbz%nkabc(3)
      endif   

      if(cmdopt('--vbmonly',9,0,strn)) then
         print *,'--vbmonly mode'
         vol = slat%vol
         goto 8888
      endif 
  
!! getq mode  Current version is not for spin dependent nor many restrictions.
!! spin symmetic (or nspin=1, not 2 channell binded and so on...
      if(cmdopt('--getq',6,0,strn).and.procid==master) then
         print *,'get Q from ql given by lmf'
         print *,'WARN current version is only for spin symmetric; and not more than 2(2l+1) occupancy'
         nspec=sctrl%nspec
         lmxa = maxval(sspec(1:nspec)%lmxa)
         allocate(qrmt(0:lmxa,nspec),qset(0:lmxa,nspec),qatot(nspec))
         qrmt= -1d-10
c         ifiwv  = fopna('veswavatm',-1,1)
         ifiwv= ifile_handle()
         open(ifiwv,file='veswavatm.'//trim(sname))
         do is=1,nspec
            read(ifiwv,*) qatot(is)
            print *,'qatom tot=',qatot(is)
            read(ifiwv,*)
            do 
               read(ifiwv,*)isp,il,npri,evll,qrmtx
               if(isp==9) exit  !find terminator line for each atom in veswavatm file
               if( isp > 1) cycle
               if(evll<0d0) then
                  if (qrmt(il,is) >1d-2) then
                     call rx('not yet implemented when two per l (eg. 3d and 4d) are bounded')
                  else     
                     qrmt(il,is)=qrmtx
                  endif
               endif
            enddo
            do il=0,lmxa
               if(qrmt(il,is)>0d0) print *, is,il,qrmt(il,is)
            enddo   
         enddo  
         close(ifiwv)
!! WARN. not averaged yet for given specs...
c         ifiqb  = fopna('qbyl',-1,1)
         ifiqb= ifile_handle()
         open(ifiqb,file='qbyl.'//trim(sname))
         ifiwv= ifile_handle()
         open(ifiwv,file='veswavatm.'//trim(sname))
         read(ifiqb,*)
         allocate(ql(0:lmxa,nspec))
         do ib=1,nbas
            read(ifiqb,*)ibas,isp,idummy,qlx(0:lmxa)
            is = ssite(ib)%spec
            ql(0:lmxa,is)=qlx(0:lmxa)
            write(6,"(' ibas ql=',i3,20f12.6)") ib,ql(0:lmxa,is)
         enddo
         do is=1,nspec
            qsetsum=0d0
            do il=0,lmxa
               if(qrmt(il,is)<0d0) then
                  qset(il,is)=0d0
               else
                  qset(il,is)= ql(il,is)/qrmt(il,is)
                  qsetsum = qsetsum+qset(il,is)
                  if(debug) print *,'ddddddddd ', is,il ,ql(il,is),qrmt(il,is),qset(il,is)
               endif
            enddo  
            write(*,"('not renormalized ',i3,' Q=',12f12.4)")is, qset(0:lmxa,is)
            qset = qatot(is)* qset/qsetsum
            write(*,"('    renormalized     Q=',12f12.4)") qset(0:lmxa,is)
         enddo
         close(ifiqb)
         call rx0('--getq mode done')
      endif


!! === MPI Setup continue===
      numprocs=1                !when no MPI
#if MPI | MPIK
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
      call strcop(shortname(procid),name,10,'.',i)
      namelen(procid) = i-1
#endif
      write(strprocid,'(i10)') procid
      strprocid = adjustl(strprocid) 
      ixx = master_mpiset(procid)
      savez=.false.
      getz=.false.
      if (maxit .eq. 0) return
c      siginit=.true. !commented out on mar2015. this makes reading sigma at every iteration.
      napw = 0
      ipr  = iprint() !When procid/=master, iprint=0.
      ipl  = ipr
      nsp  = globalvariables%nsp
      nspc = globalvariables%nspc
c      lso =   isw ( iand(4,int(sctrl%lncol))  .ne.0 ) 
c     &  + 2 * isw ( iand(32,int(sctrl%lncol)) .ne.0 )
c     &  + 3 * isw ( iand(64,int(sctrl%lncol)) .ne.0 )
!! lso is read by m_rdctrl.F 
!! lso  =0 (no so): =1(L.S): =2(LzSz).

      lcplxp = 0
      lekkl = int(sctrl%pfloat)
      if (lso .ne. 0) lcplxp = 1
      if (isum(nbas,lldau,1) .ne. 0) lcplxp = 1

!! Dimension maximum of Hamiltonian is ndhamx (ispx=1,nspx)
!! spinoffdiag=T case, the spin loop is too complicated --> we need to fix this.
      if(nspc==2.and.nsp==2) then !nspc==2
         spinoffdiag=.true.
         nspx=1
         ndhamx = ndham*2
      elseif(nspc==1.and.nsp==2) then
         spinoffdiag=.false.
         nspx=2
         ndhamx= ndham
      elseif(nspc==1.and.nsp==1) then
         spinoffdiag=.false.
         nspx=1
         ndhamx= ndham
      else
         call rx('bndfp: nspc==2 but nsp=1')
      endif   

!! related to the error in rseq. clean up RSEQ_ERROR at first because
!! this error is only fatal when it is finaly converged. takaoNov28_2010
      idummy=unlink('RSEQ_ERROR')

ctakao sep2010 !complicated.
c      nspc=2 for non-colinear case. 
c      allocate(evl(ndham,2))
c      call i8icopy(size(sham%ldham),sham%ldham,1,ldham,1)
      allocate(evl(ndhamx,nspx))
      ldham = sham%ldham
      ldim = ldham(1)
      if(size(sham%ldham)/=size(ldham))call rx("bndfp:size(sham%ldham)/=size(ldham)")
      pwmode = sham%pwmode
      if (procid .eq. master) then
         ifiese=ifile_handle()
         open(ifiese,file='PWMODE')
         write(ifiese,*)pwmode
         close(ifiese)
      endif 
      pwemin =sham%pwemin
      pwemax =sham%pwemax
      onesp = 0
c      nfbn(1) = 0
c      nfbn(2) = 0
      ldos = iand(1,int(sctrl%ldos))
      lrsig = int(sham%lsig)
      lfrzw = isw ( iand(16,int(sctrl%lbas)) .ne.0 )
      nspec= sctrl%nspec
      if(ldos - mod(sctrl%ldos/1, 2)/=0) call rx('ldos - mod(sctrl%ldos/1, 2)/=0')
      if(lfrzw- mod(sctrl%lbas/16,2)/=0) call rx('lfrzw- mod(sctrl%lbas/16,2)/=0')
      if(lrsig- sham%lsig/=0) call rx('lrsig- sham%lsig/=0')

      if (procid .eq. master) then
         fxsts= fxst('sigm')
         if (lrsig .ne. 0 .and. fxsts/=1) then
            write(6,*)' bndfp (warning): no sigm file found ... LDA calculation only'
            sham%lsig=0
            lrsig = 0
         endif
      endif
      call mpibc1(lrsig,1,2,mlog,'bndfp','lrsig')
      call mpibc1(sham%lsig,1,2,mlog,'bndfp','sham%lsig') 
      jobgw= -999
      epsovl = sham%oveps

      if(llmfgw) then
         i = 7
         if (cmdopt('-jobgw=',i,0,strn)) then
            if (.not. a2bin(strn,jobgw,2,0,' ',i,-1)) call
     .        rxs2('BNDFP: failed to parse "',strn(1:30),'%a"')
         endif
      endif
      if (procid .eq. master) then
        aaachar=trim(i2char(iter))//" of "//trim(i2char(maxit))
        write(6,"(a)") trim(" --- BNDFP:  begin iteration "//aaachar)
      endif
      if (nbas.gt.nbx) call rxi('bndfp: nbx exceeded, need',nbas)
      if (ndham.le.0)  call rx('bndfp: hamiltonian matrix has zero dimension')
!!
      eterms=-99d0
      eterms(19) = 0d0
      sham%eterms=eterms
      nvl  =spot%nlml
      nchan=spot%nlma
      eks = 0
      alat =slat%alat
      plat =slat%plat
      qlat =slat%qlat
      ngabc=slat%nabc
      vol = slat%vol
      k1 = n1
      k2 = n2
      k3 = n3

!! ... for BZ integration
      nkabc= sbz%nkabc
      nkp=   sbz%nkp
      ntet=  sbz%ntet
      lmet=  sbz%lmet
      mpsord=sbz%n
      esmear=sbz%w
      nevmx= sbz%nevmx
      efmax= sbz%efmax
      fsmom= sbz%fsmom
      if(debug) print *,'ttt nevmx w=',nevmx,esmear
      esmear = esmear + iabs(mpsord)
      if (mpsord .lt. 0) esmear = -esmear
      mpsord = int(dabs(esmear))
      if (esmear .lt. 0) mpsord = -mpsord
      ndos=sbz%ndos
      dosw= sbz%dosw
      ef0 = sbz%ef
      def = sbz%def
!!
      ltet = ntet .gt. 0
      lwndow = cmdopt('--window=',9,0,strn)
      if (lwndow) then
         if (lmet .eq. 0) then
            call rx(' bndfp: restart with METAL=2 for --window')
         endif
         if (.not. ltet) then
            call rx(' bndfp: restart with TETRA=T for --window')
         endif
         iq = 0
         i = parg('window=',4,strn,iq,len(strn),', ',2,2,iv,dosw)
         call info2(20,0,0,' BNDFP: generating density in energy window %2:1d',dosw,0)
         lfrce = 0
         lpnu = 0
         efmax = 1d2
         nevmx = ndham
         if (lrout .eq. 0)
     .     call rx('--window incompatible with no output density')
         call info0(20,0,0,' Delete band weights file ...')
         ifi = fopna('wkp',-1,4)
         call dfclos(ifi)
      endif
C     Always need DOS if just sampling
      if (lmet .ne. 0 .and. .not. ltet) ldos=1
      qbg = sctrl%zbak(1)
      alfa = 0

C     Switch to plot bands at specified qp
      if (cmdopt('--band',6,0,strn).or.fullmesh) then
         plbnd = 1
c         plbopt = strn(7:)
         lrout = 0
         lfrce = 0
         nkp = 0
         numq = 1
         allocate(ifbls_iv(ndhamx*2))
         ifbls_iv(:)=0
      else
         plbnd = 0
      endif
      if (lmet .eq. 1) call rx('bndfp: lmet=1 not implemented')
!! isanrg is logical function,       call isanrg(lmet,0,4,'bndfp:','lmet', .true.)
      l_dummy_isanrg=isanrg(lmet,0,4,'bndfp:','lmet', .true.)

!! == PROCAR file 
c      PROCARon= .not.cmdopt('--noprocar',10,0,strn)
      if(PROCARon) then
         iprocar1=ifile_handle()
         open(iprocar1,file='PROCAR.UP.'//trim(strprocid))
         if(nsp==2) then
            iprocar2=ifile_handle()
            open(iprocar2,file='PROCAR.DN.'//trim(strprocid))
         endif
      endif

!! spin-symmetrized pnu. aug2019. See also in pnunew and locpot
      phispinsym= cmdopt('--phispinsym',12,0,strn)
      if(phispinsym) then
         if(master_mpi()) write(6,*) 'bndfp: --phispinsym use spin-averaged potential for phi and phidot' 
         do ib = 1,nbas
            i_copy_size=size(ssite(ib)%pnu)
            call dcopy(i_copy_size,ssite(ib)%pnu,1,pnu,1)
            i_copy_size=size(ssite(ib)%pz)
            call dcopy(i_copy_size,ssite(ib)%pz,1,pnz,1)
            is = int(ssite(ib)%spec)
            lmxa = sspec(is)%lmxa
            do l=0,lmxa
              pmean = sum(pnu(l+1,1:nsp))/nsp
              if(master_mpi().and.nsp==2)   write(6,"('  ibas l=',2i3,' pnu=',2f10.5,' -->',f10.5)")
     &             ib,l,pnu(l+1,1:nsp),pmean
              pnu(l+1,1:nsp) = pmean
              lpz = pnz(l+1,1) /= 0
              if (lpz) then
                 pmean = sum(pnz(l+1,1:nsp))/nsp
                 if(master_mpi().and.nsp==2) write(6,"('  ibas l=',2i3,' pnz=',2f10.5,' -->',f10.5)")
     &                ib,l,pnz(l+1,1:nsp),pmean
                 pnz(l+1,1:nsp) = pmean
              endif
            enddo
            i_copy_size=size(ssite(ib)%pnu)
            call dcopy(i_copy_size,pnu,1,ssite(ib)%pnu,1)
            i_copy_size=size(ssite(ib)%pz)
            call dcopy(i_copy_size,pnz,1,ssite(ib)%pz,1)
         enddo   
      endif
      
      

!! jan2018kaist
      writeham= cmdopt('--writeham',10,0,strn)
      if(writeham) then
        ifih=ifile_handle()
        open(ifih,file='HamiltonianPMT.'//trim(strprocid),form='unformatted')
!!      
        if(procid == master) then ! pair-table and weight generator.  nlat,qgt are not used yet...
          nkk1=sbz%nkabc(1)
          nkk2=sbz%nkabc(2)
          nkk3=sbz%nkabc(3)
          nkp = nkk1*nkk2*nkk3
          allocate(qplist(3,nkp))
          iq=0
          do ik3=1,nkk3
          do ik2=1,nkk2
          do ik1=1,nkk1
             iq=iq+1
             qplist(:,iq) =  qlat(:,1)*dble(ik1-1)/nkk1
     &        +   qlat(:,2)*dble(ik2-1)/nkk2
     &        +   qlat(:,3)*dble(ik3-1)/nkk3
          enddo
          enddo
          enddo
          do ip=1,nkp
            write(6,"(' qplist:',i9,3f10.4)")ip, qplist(:,ip)
          enddo

c          do ix=1,3
c            write(ififft,"(3f25.15,' ! plat' )") plat(:,ix)
c          enddo
c          write(ififft,"(3i5,'  ! nk')") nkk1,nkk2,nkk3
c          write(ififft,"(i5,' !nbas')") nbas
c          do ibas=1,nbas
c            write(ififft,"(3f25.15,' ! bas' )") pos(:,ibas)
c          enddo
          allocate(npair(nbas,nbas),pos(3,nbas))
          call dcopy(nbas*3 , slat%rv_a_opos , 1 , pos , 1 )

          npairmx=nkk1*nkk2*nkk3*2 !initial size of npairmx
          do   !this loop is just to get reasonable npairmx satisfying npairmx >npair
            npairmx = npairmx + (nkk1*nkk2*nkk3+1)*.5 !+1 added Feb2014 for 1x1x1
            allocate( nlat(3,npairmx,nbas,nbas), nqwgt(npairmx,nbas,nbas) )
            call gennlat(pos,nbas,plat,nkk1,nkk2,nkk3,npairmx,ok,npair,nlat,nqwgt)
            if(ok) exit
            deallocate( nlat, nqwgt )
          enddo
          print *,'We got nlat and nqwgt'
          do ib1=1,nbas
          do ib2=1,nbas
              write(6,"(3i8,' !ib1 ib2 npair -------------')") ib1,ib2,npair(ib1,ib2)
c              if(abs(sum(qwgt(1:npair(ib1,ib2),ib1,ib2))-1d0)>1d-8) call rx('buggg:qwgt sum is not unity')
              do ni = 1,npair(ib1,ib2) 
                posp =  pos(:,ib1)-pos(:,ib2) + matmul(plat,nlat(:,ni,ib1,ib2)) ! R_i+T - R_j
                rrrr = sqrt(sum(posp**2))
                write(6,"(i6,3x,3i4,i3,x,f8.3)") ni, nlat(1:3,ni,ib1,ib2),nqwgt(ni,ib1,ib2),rrrr
              enddo
          enddo
          enddo

          ififft=ifile_handle()
          open(ififft,file='HamiltonianPMTInfo',form='unformatted')
          write(ififft) plat,nkk1,nkk2,nkk3,nbas,qlat
          write(ififft) pos,alat
          write(ififft) qplist
          write(ififft) npair,npairmx
          write(ififft) nlat,nqwgt
!! ---------------
!! delta fun check: k --> T --> k 
!!     \delta_{kk'} = \sum_{T \in T(i,j)} W_T exp( i (k-k') T)
          ikpd=1
          do ikp=1,nkp
            qp = qplist(:,ikp) - qplist(:,ikpd)
          do ib1=1,nbas
          do ib2=1,nbas
            aaaa=0d0
            do ni = 1,npair(ib1,ib2) 
                aaaa =  aaaa + 1d0/(nkp*nqwgt(ni,ib1,ib2))*exp(img*2d0*pi* sum(qp*matmul(plat,nlat(:,ni,ib1,ib2))))
            enddo
            cccx=''
            if(ikp==ikpd) cccx=' <--'
            write(6,"('\delta-fun test',3f10.4,2i3,2f23.15,a)") qplist(:,ikp),ib1,ib2,aaaa,cccx
          enddo
          enddo
          enddo
          deallocate(pos,npair)
c          stop 'xxxxxxxxxxx test end xxxxxxxxx2222'
        endif

      endif

!! --- Define local arrays used in the generation of the potential ---
      lnoxc = 0
      if (jobgw .eq. 2) lnoxc = 1
      allocate(qmom_rv(nvl))
      allocate(gpot0_rv(nvl))
      allocate(vval_rv(nchan))
      allocate(hab_rv(nab*n0*nsp*nbas))
      allocate(vab_rv(nab*n0*nsp*nbas))
      allocate(sab_rv(nab*n0*nsp*nbas))
      allocate(ppnl_rv(nppn,n0,nsp,nbas))
      allocate(fes1_rv(3*nbas))

!! dfaugm alllocates and initializes variables
      call dfaugm ( nbas , lcplxp , lso , ssite , sspec , sv_p_osig 
     .  , sv_p_otau , sv_p_oppi )

!! --- Make the potential sans XC part ---
      if(llmfgw) then ! only for lmfgw program
         if(jobgw .eq. 1 .or. jobgw .eq. -999) then
            call info(20,1,0,' Make potential without XC part ...',0,0)
            call togpr()
            i = 1 + 10*lfrzw + 100
!! dfaugm alllocates and initializes variables
            call dfaugm ( nbas , lcplxp , lso , ssite , sspec , sv_p_osigx 
     .        , sv_p_otaux , sv_p_oppix )
            allocate(spotx(k1,k2,k3,nsp))
            spotx=0d0
!! mkpot calls locpot. and locpot calls augmat.
!! augmat calculates sig,tau,ppi.
            call mkpot ( nbas , ssite , sspec , slat , sham , spot , 0 ,
     .        lcplxp , k1 , k2 , k3 , spot%zv_a_osmrho , spot%sv_a_oorhat , qbg , spotx 
     .        , qmom_rv , vconst , vrmt , sv_p_osigx , sv_p_otaux , sv_p_oppix 
     .        , ppnl_rv , hab_rv , vab_rv , sab_rv , qval , qsc , gpot0_rv 
     .        , vval_rv , fes1_rv , i , vorb , nlibu , lmaxu , lldau , nsp,0 )
            call togpr()
         endif
!! --- readin job switch
         if (jobgw .eq. -999 .and. procid==master) then
            write(stdo,*) ' lmfgw: input one of the following jobs:'
            write(stdo,*) '   0 : init mode; creates files',
     .        ' SYMOPS, LATTC, CLASS, NLAindx, ldima'
            write(stdo,*) '   1 : GW setup mode'
            write(stdo,*) '   4 : band mode '
            write(stdo,*) '   5 : eigenvalue-only mode '
            write(stdo,*) ' job?'
            read (5,*) jobgw
         endif
#if MPI|MPIK
         call mpi_barrier(MPI_comm_world,ierr)
#endif        
         call mpibc1(jobgw,1,2, mlog,'bndfp','jobgw')
      endif

C --- Make the potential and total energy terms ---
      i = 1 + 10*lfrzw + 100*lnoxc
      if (cmdopt('--wrhomt',8,0,strn)) then
         i = i + 10000
      else if (cmdopt('--wpotmt',8,0,strn)) then
         i = i + 20000
      endif
      if(llmfgw) i = i + 10000
ccccccccccccccccccccccccc
ctakao test mkekin repeat
c          print *,'vvv takao mkekin repeat init 01 ---'
c          call mkekin(nbas,ldim,ssite,sspec,slat,sham,osig,otau,oppi,
c     .      lcplxp,oqkkl,k1,k2,k3,vconst,w(osmpot),w(osmrho),sev,sumtvxxx)
c          stop 'xxxxxxxxxxxxxxxx'
ccccccccccccccccccccccccc

!! mkpot calls locpot. and locpot calls augmat.
!! augmat calculates sig,tau,ppi.
      call togpr()
      call mkpot ( nbas , ssite , sspec , slat , sham , spot , lfrce
     .  , lcplxp , k1 , k2 , k3 , spot%zv_a_osmrho , spot%sv_a_oorhat , qbg , spot%zv_a_osmpot 
     .  , qmom_rv , vconst , vrmt , sv_p_osig , sv_p_otau , sv_p_oppi 
     .  , ppnl_rv , hab_rv , vab_rv , sab_rv , qval , qsc , gpot0_rv 
     .  , vval_rv , fes1_rv , i , vorb , nlibu , lmaxu , lldau , nsp,0 )
      call togpr()
ccccccccccccccccccccccccc
ctakao test mkekin repeat
c       print *,'vvv takao mkekin repeat init 02 ---'
c       call mkekin(nbas,ldim,ssite,sspec,slat,sham,osig,otau,oppi,
c     .      lcplxp,oqkkl,k1,k2,k3,vconst,w(osmpot),w(osmrho),sevxxx,sumtvxxx)
c          stop 'xxxxxxxxxxxxxxxx  111'
ccccccccccccccccccccccccc

      elind = sham%elind
      if (elind .lt. 0d0) elind=-(3*pi**2*(qval-qsc-qbg)/vol)**.66666d0*elind
      dmxp(33) = elind
c      allocate( evlall(ndham,nsp,nkp))
      allocate( evlall(ndhamx,nspx,nkp))
      evlall = 0d0 

!! --- Setup for BZ integration ---
      if (plbnd .eq. 0) then
         call subzi(lmet,ltet,lrout.gt.0,ndham,nsp,nspc,nkp,qval-qbg,
     .     nevmx,lwtkb,eferm,rv_a_owtkb,numq,lswtk,rv_a_oswtk,ef0,def,mpsord,
     .     dabs(esmear-mpsord))
         if (lwtkb .eq. -1)
     .     call info(20,0,0,' Start first of two band passes ...',0,0)
         if (lwtkb .eq. 1) then
            if (ef0 .ne. eferm)
     .        call info(20,0,1,'%8pReplace ef0 with file ef=%;6d',eferm,0)
            ef0 = eferm
         endif
      else
         nkp = 0
         ldos = 0
         lwtkb = -1
         icls = 0
      endif

!! this section is moved from addrbl.F !nov2015
!! we have to examine lswtk mechanism. !nov2015 this is temporary fix
      if (.not.(tdos .or. PROCARon)) then  !if(ndimh*nspc /= nevl) then 
         call info0(30,0,0,' bndfp: eigenvector matrix not square ... spin weights not evaluated')
         lswtk = -2
      endif

!! --- More k-point independent local arrays ---
      if (ldos .eq. 0 .and. lmet .eq. 0) ndos = 1
      if (ndos .lt. 0) ndos = -ndos
      allocate(dos_rv(ndos,2*nsp))
      if (lrout .ne. 0) then
         call dfqkkl ( nbas , ssite , sspec , numq , sv_p_oqkkl )
         if ( lekkl .eq. 1 ) call dfqkkl ( nbas , ssite , sspec , numq 
     .     , sv_p_oeqkkl )
         allocate(srout_zv(k1*k2*k3*numq*nsp))
         allocate(fh_rv(3*nbas))
         allocate(fes2_rv(3*nbas))
      endif
      call suham2(slat,sspec,ssite,sham,spot) !,sstrn,sctrl
c      print *,'end of suham2'

!! --- Options for core level specta (CLS) ---
!!    Allocate aus for all qp
      if (cmdopt('--cls',5,0,strn)) then
         if (lmet .gt. 0 .and. (lmet .ne. 2 .and. lmet .ne. 3))
     .     call rx('For CLS restart with METAL=2 or 3')
         icls = 1
         clsopt = strn(6:)
         efmax = 1d3
         nevmx = ndhamx
         if (lrout .eq. 0) call rx('bndfp: need output density for cls')
         call suclst(nsitmx,nbas,nsp,ssite,sspec,clsopt,
     .     isite,iclsl,iclsn,nsites)
         efmax = 1d3
         nevmx = ndhamx
         nlmax = globalvariables%mxorb / globalvariables%nkaph
         if (16*3*nlmax*ndham*nbas*nsp*nkp/1000000 .gt. 24 .and. procid
     .     .eq. master) then
!! this path is go through by crn test.
c     .     call awrit6(' CLS: %iMb memory for aus: nlmax=%i ndham=%i '//
c     .     'nsistes=%i nsp=%i nkp=%i',' ',120,
c     .     i1mach(2),16*3*nlmax*ndham*nsites*nsp*nkp/1000000,
c     .     nlmax,ndham,nsites,nsp,nkp)
            aaachar=
     &        ' CLS: '//trim(i2char( 16*3*nlmax*ndham*nsites*nsp*nkp/1000000))//
     &        ' Mb memory for aus: nlmax='//trim(i2char( nlmax))//
     &        ' ndham='  //trim(i2char( ndham))//
     &        ' nsistes='//trim(i2char( nsites))//
     &        ' nsp='    //trim(i2char( nsp))//
     &        ' nkp='    //trim(i2char( nkp))
            write(6,"(a)") aaachar
         endif
         allocate(ausc_zv(3*nlmax*ndham*nsites*nsp*nkp))
         ausc_zv(:)=0.0d0
      else
         icls = 0
      endif

ctakaox june2009
c      ifiogw=9038
c      inquire(unit=ifiogw,opened=nexist)
c      if(nexist)
c     &stop 'bndfp: I assumed ifio=9038 is not used(too sloppy)'
c      open(ifiogw,file='NBANDGW')
Cki  ifiese default value is not set.  set it to 0

c$$$!! We need ESEAVR file. This is a dummy file when we use huge HAM_SIGP_EMAX such as EMAX=9999.
c$$$!! ESEAVR is called from rdsigm2-psigd.
c$$$      ifiese=0
c$$$      if(lrsig/=0) then
c$$$         ifiese=ifile_handle()
c$$$         open(ifiese,file='ESEAVR',status='old',err=127)
c$$$         goto 128
c$$$ 127     continue
c$$$         open(ifiese,file='ESEAVR')
c$$$         write(ifiese,*)"0d0 1 !dummy when no ESEAVR"
c$$$         close(ifiese)
c$$$         open(ifiese,file='ESEAVR',status='old',err=127)
c$$$      endif
c$$$ 128  continue


!! == Start loop over k points; also, re-entry for second band pass ==
!! two pass mechanism, historically. ix99=1 and ix99=2 for tetrahedron method.
!!     1st for eigenvalue (and save eigenfunctions to eigze_* files), resulting tetrahedron weight.
!!     2nd for re-read eigze_* and accumulate density.                           
!!    It will be better to remove this double path mechanism in future.
      gotosecondpass =.false.
      if(plbnd==0) then
         gotosecondpass= (lmet>=0 .and. (lmet/=4 .or. ltet) .and. lrout>0 )!this covers limited cases.
!this judges going through second pass or not.
      endif 
      ix99=0
      xdat=0d0
 99   continue
      ix99=ix99+1               !ix99=2 for second pass
      ebot = 1000d0
      call surho ( nbas , ssite , sspec , lmet , ldos , lrout , lekkl 
     .  , numq , sv_p_oqkkl , sv_p_oeqkkl , k1 , k2 , k3 , srout_zv , 
     .  ndos , dos_rv , sumev , sumqv )

      if (lfrce .gt. 0) then
         call dpzero(frc, 3*nbas*numq)
         fh_rv=0d0              ! 3 * nbas 
      endif
      if (lswtk .eq. 1) then
         efmax = 1d3
         nevmx = ndhamx
         call dpzero ( rv_a_oswtk , ndhamx * nkp )
      endif
      nl = int(sctrl%nl)

!! Mullikan mode
      mullwrite = .false.
      if(cmdopt('--mull',6,0,strn).and.procid==master) then
#if MPIK
         call rx('MULL not k-parallel yet... Need to fix this if necessary')
#endif
         nfilem = fopna('moms',-1,4)
         mullwrite = .true.
         rewind nfilem
         i = 1
C   ... If switch '--mull'; get sites, number of channels
         if (cmdopt('--mull',6,0,strn)) then
            mulopt = strn(7:)
            ng=slat%nsgrp
            efmax = 1d3
            nevmx = ndhamx
            nchmx = min(1024,nbas*nl**2)
            allocate(chan0_iv(nchmx))
            call sumlst ( 1 , nchmx , nbas , ng , ssite , sspec , mulopt 
     .        , moddos , nsites , lsites , lmxch , nchan , chan0_iv, 0 )
            if(allocated(chan_iv)) deallocate(chan_iv)
            allocate(chan_iv(nchan))
            chan_iv=chan0_iv(1:nchan)
            deallocate(chan0_iv)
            call iomomn(.true.,2,.false.,1,nspc,1,1,i)
         endif
         if (procid .eq. master) i = iomoms(-nfilem,nl,nsp,nspc,nkp,
     .     ndham,i,1,0,1,0,0,0,0,0,0d0,0d0,0d0,0d0,0d0,0d0)
      endif
      if(debug) print *,'vvvvv 22222222qqq '

C ... Case only generate bands at supplied qp: setup
      if (plbnd .ne. 0) then
C       Try and read Fermi level from file
         if (procid .eq. master) then
            ifi = fopna('wkp',-1,4)
ckino getef is integer          call getef(ifi,0,ef0)
            i_dummy_getef= getef(ifi,0,ef0) !Is this necessar???
            call fclr('wkp',ifi)
         endif
         call mpibc1(ef0,1,4,.false.,'bndfp','ef0')
      endif   

!! ... For insulator, valence band top and conduction band bottom
      evtop = -99999
      ecbot = -evtop

!! ... Setup for spin-orbit coupling
      if (lso .ne. 0) then
         nlmax = globalvariables%nlmax
         if(allocated(orbtm_rv)) deallocate(orbtm_rv)
         allocate(orbtm_rv(nl*nsp*nbas))
         orbtm_rv(:)=0.0d0
      endif

!! === Get nkp and qplist =========== feb2015
!! fatband mode ---> --mkprocar
!! plbnd/=0 ---> band plot mode.
      if(debug) write(*,*) 'bndfp: uuuuuuu1111111'
      if(allocated(qplist)) deallocate(qplist)
      if( procid==master .and. plbnd==0) then
         allocate(qplist(3,nkp))
         call dcopy(3*nkp, sbz%rv_p_oqp,1,qplist,1)
c-----------------
         ifqplist=ifile_handle()
         open(ifqplist,file='QPLIST.IBZ')
         do iq=1,nkp
            write(ifqplist,"(i5,3f23.15,3x,3f23.15)")iq, qplist(:,iq), sbz%rv_a_owtkp(iq)/2d0
         enddo
         close(ifqplist)
         if( jobgw/=0.and.jobgw/=1 )then
            write(6,*)'qplist.ibz sumcheck=',nkp,sum(sbz%rv_a_owtkp(1:nkp))
            if(abs(sum(sbz%rv_a_owtkp(1:nkp))-2d0) > 1d-12) call rx('bndfp.F: need check! wtkp sum')
         endif   
c-----------------
         !I think sbz%rv_p_oqp is not used in the following.
      elseif(procid==master) then
         if (cmdopt('--onesp',7,0,strn) .and. nspc .eq. 1) onesp = 1
         if(debug) print *,'vvvvv nkp=',nkp
ckino plotting wavefunction
c         write(*,*)'allocate listwf , size=',ndham
c         if (allocated(listwf))  deallocate(listwf)
c         allocate(listwf(ndham))
c         n_listwf=ndham; listwf=0
         if(fullmesh) then
            nkp = nkk1*nkk2*nkk3
            allocate(qplist(3,nkp))
!! Fermi surface version for xcrysden
            if(fsmode) then
               iq=0
               do ik1=1,nkk1    !ordering is differnt from procaron case...
                 do ik2=1,nkk2
                   do ik3=1,nkk3
                     iq=iq+1
                     qplist(:,iq)   =  qlat(:,1)*dble(ik1-1)/(nkk1-1)
     &                 +   qlat(:,2)*dble(ik2-1)/(nkk2-1)
     &                 +   qlat(:,3)*dble(ik3-1)/(nkk3-1)
c                     call dpscop ( qp , sbz%rv_p_oqp , 3 , 1 , 3 * iq - 2 , 1d0 )
                   enddo
                 enddo
               enddo
            endif   
!! pdos (--mkprocar and --fullmesh)
            if(procaron) then !to fit to the tetirr.F
               iq=0
               do ik3=1,nkk3
                 do ik2=1,nkk2
                   do ik1=1,nkk1
                     iq=iq+1
                     qplist(:,iq) =  qlat(:,1)*dble(ik1-1)/nkk1
     &                 +   qlat(:,2)*dble(ik2-1)/nkk2
     &                 +   qlat(:,3)*dble(ik3-1)/nkk3
c                     call dpscop ( qp , sbz%rv_p_oqp , 3 , 1 , 3 * iq - 2 , 1d0 )
                   enddo
                 enddo
               enddo
            endif
         else 
!! syml direct read for plbnd mode. See "call writeband" below. feb2015 
!!   Read efermi.lmf
            write(6,*)' --- Readin efermi.lmf --- '
            ifi=ifile_handle()
            open(unit=ifi,file='efermi.lmf',status='old',err=1012)
            read(ifi,*,err=1012) eferm
            read(ifi,*,err=1012) evtop
            read(ifi,*,err=1012) ecbot
            close(ifi) 
            goto 1013
 1012       continue
            call rx('No efermi.lmf!: copy it, or run lmf-MPIK (sc mode) to get efermi.lmf.')
 1013       continue

!! --- example of syml file ---
! ndiv qleft(1:3) qright(1:3) llabel rlabel  ndiv2 ninit2 nend2 etolv(Ry) etolc(Ry) 
! 5  0 0 0   .5 .5  .5        Gamma  L       1025  1  16     0.1      0.01 
! 5  0 0 0    1.  0  0        Gamma  X    
! 5  0 0 0   .75 .75 0        Gamma  K       1025  8  38     0.1      0.01 
!! As this shows you can add or not add line after ndiv2. These are for highly resolved calculations.
            write(6,*)' --- Readin syml file --- '
            ifisyml=ifile_handle()
            open(ifisyml,status='old',file='syml.'//trim(sname))
            nsyml=0
            nsymln=0
            nqp2_syml=0
c            massline=.false.
            do 
              if(nsyml+1>nsymlmax) call rx('bndfp: Enlarge nsymlmax')
              read(ifisyml,"(a)",end=1015) schar
              if(len(trim(schar))==0.or.schar(1:1)=='#' .or. schar(1:1)=='!' .or. schar(1:1)=='%') cycle !comment line
              read(schar,*,err=1014,end=1014) 
     &         nqp_syml(nsyml+1), qps_syml(1:3,nsyml+1), qpe_syml(1:3,nsyml+1),
     &         labeli(nsyml+1),labele(nsyml+1), 
     &         nqp2_syml(nsyml+1),nqp2s_syml(nsyml+1),nqp2e_syml(nsyml+1),etolv,etolc
              masslineon(nsyml+1)=.true.
              nqp2n_syml(nsyml+1)= nqp2e_syml(nsyml+1)-nqp2s_syml(nsyml+1)+1
              nqps_syml(nsyml+1)=1
              nqpe_syml(nsyml+1)=nqp_syml(nsyml+1)
              write(6,"(' ',i4,3f9.4,' ',3f9.4,' ',a,' ',a,'  Massl:div,init,end=',3i5)")
     &         nqp_syml(nsyml+1), qps_syml(1:3,nsyml+1), qpe_syml(1:3,nsyml+1),
     &         trim(labeli(nsyml+1)),trim(labele(nsyml+1)),
     &         nqp2_syml(nsyml+1),nqp2s_syml(nsyml+1),nqp2e_syml(nsyml+1)
              goto 1025

 1014         continue
              read(schar,*,err=1015,end=1015) 
     &         nqp_syml(nsyml+1), qps_syml(1:3,nsyml+1), qpe_syml(1:3,nsyml+1),
     &         labeli(nsyml+1),labele(nsyml+1) 
              masslineon(nsyml+1)=.false.
              nqps_syml(nsyml+1)=1
              nqpe_syml(nsyml+1)=nqp_syml(nsyml+1)
              write(6,"(' ',i4,3f9.4,' ',3f9.4,' ',a,' ',a)")
     &         nqp_syml(nsyml+1), qps_syml(1:3,nsyml+1), qpe_syml(1:3,nsyml+1),
     &         trim(labeli(nsyml+1)),trim(labele(nsyml+1))
              nsymln = nsymln+1
              nqp2n_syml(nsyml+1)= 0
 1025         continue
              if(nqp_syml(nsyml+1)==0) exit
              nsyml = nsyml + 1
            enddo
 1015       continue
            close(ifisyml)
            nkp = sum(nqp_syml(1:nsyml)+nqp2n_syml(1:nsyml))
            allocate(qplist(3,nkp))
            if(allocated(xdatt)) deallocate(xdatt)
            allocate(xdatt(nkp))
            totxdatt=0d0
            ikp=0
            do isyml=1,nsyml
              dqsyml(isyml) = dsqrt(sum((qpe_syml(1:3,isyml) -qps_syml(1:3,isyml))**2))
              do i=1,nqp_syml(isyml)+nqp2n_syml(isyml)
                ikp= ikp+1
                if(i<=nqp_syml(isyml)) then
                  rq = dble(i-1)/(nqp_syml(isyml)-1)
                else  
                  ii= i-nqp_syml(isyml)-1
                  rq = dble(nqp2s_syml(isyml)-1+ii)/(nqp2_syml(isyml)-1)
                endif
                qplist(:,ikp)= (1d0-rq)*qps_syml(1:3,isyml) +rq*qpe_syml(1:3,isyml)
                xdatt(ikp) = totxdatt + dqsyml(isyml)*rq
              enddo  
              totxdatt = totxdatt + dqsyml(isyml)
            enddo  
            write(6,"('nsyml nkp=',3i5)") nsyml,nkp
         endif
         if (nkp <= 0) call rx('bndfp: nkp<=0') ! quit if nkp==0
      endif

!! broadcase nkp and qplist
      call mpibc1( nkp,1,2,mlog,'bndfp','nkp')
      call mpibc1( onesp,1,2,mlog,'bndfp','onesp')
      if(procid/=master) allocate(qplist(3,nkp))
      call mpibc1 ( qplist, 3 * nkp , 4 , .false. , 'bndfp' , 'qp'  )
      if (procid.eq.master.and.ix99==1.and.oncewrite(11)) then
         ifqplist=ifile_handle()
         open(ifqplist,file='QPLIST')
         print *,'-------- qplist --------'
         iq=0
         do isyml=1,nsyml
         do i=1,nqp_syml(isyml) + nqp2n_syml(isyml)
            iq=iq+1
            infoq=''
            if(i==1) infoq=' <-- isyml= '//charnum3(isyml)
            if(i==nqp_syml(isyml)+1) infoq=' <-- isyml Mass= '//charnum3(isyml)
            write(6,"(i5,3f8.3,' ',a)")iq,qplist(:,iq),trim(infoq)
            write(ifqplist,"(i5,3f23.15,x,f12.6,' ',a)")iq,qplist(:,iq),xdatt(iq),trim(infoq)
         enddo
           write(ifqplist,*)
         enddo
         close(ifqplist)
      endif 

      if(procaron.and.nsyml==0.and.procid==master) then !xdatt is dummy
        if(allocated(xdatt)) deallocate(xdatt)
        allocate(xdatt(nkp))
        xdatt=0d0
      endif
      if(PROCARon) then !if(plbnd/=0) then ! march2015
        if(procid/=master) allocate(xdatt(nkp))
        call mpibc1( xdatt, nkp , 4 , .false. , 'bndfp' , 'xdatt'  )
      endif  

#if MPIK
      call mpi_barrier(MPI_COMM_WORLD,ierr) 
#endif
      if (allocated(evlall)) deallocate(evlall)
      allocate(evlall(ndhamx,nspx,nkp)) ; evlall=0d0
      nchanp=25
      if(procaron.and.fullmesh) then
         allocate(dwgtall(nchanp,nbas,ndhamx,nsp,nkp))
         dwgtall=0d0
      endif

!! --- get index for hamiltonian for m_hamindex takao june2009
!! these are used in sigm mode(QSGW).
!! memo:
!!  ib = atom index
!!  ltab= L (angular momentum index)
!!  ktab=  =1 for EH, =2 for EH2, =3 for lo
!! SYMOPS
      if(onlyonce) then 
         onlyonce=.false.
         nbasp = sctrl%nbasp
         nl    = sctrl%nl
         call gen_hamindex ( sham , sbz , ssite , sspec , slat , nl ,
     .     nbasp , sham%iv_a_oindxo , ldim , .true. , llmfgw.or.lrsig/=0,jobgw)  !mar2015
c    .     nbasp , sham%iv_a_oindxo , ldim , siginit, llmfgw.or.lrsig/=0,jobgw)  !aug2012
         if(procid==master) then
            write(6,*)
            write(6,*) ' --- Hamiltonian index ---'
c            write(6,*) ' ib l  k offl(iorb)+1  offl(iorb)+2*l+1  trim(spec)'
            allocate(ib_table(ldim),l_table(ldim),k_table(ldim))
            do iorb = 1, norbmto !Total number of MTO's
              ib   = ibastab(iorb)
              is   = ssite(ib)%spec
c              call r8tos8(sspec(is)%name,spid)
              spid=sspec(is)%name
c              write(6,"(3i3,2x,2i5,3x,a)") ib,ltab(iorb),ktab(iorb), 
c     &          offl(iorb)+1, offl(iorb)+2*ltab(iorb)+1, trim(spid)
!!
              ib_table(offl(iorb)+1: offl(iorb)+2*ltab(iorb)+1) = ib
              l_table (offl(iorb)+1: offl(iorb)+2*ltab(iorb)+1) = ltab(iorb)
              k_table (offl(iorb)+1: offl(iorb)+2*ltab(iorb)+1) = ktab(iorb)
            enddo

            write(6,*) ' i  ib l  k trim(spec)'
            do i= 1, ldim
              ib   = ib_table(i)
              is   = ssite(ib)%spec
              spid = sspec(is)%name
              write(6,"(i3,x,3i3,x,a)")i, ib_table(i),l_table(i),k_table(i),trim(spid)
            enddo
            if(writeham) then
              write(ififft)ldim,lso,nsp
              write(ififft)ib_table,l_table,k_table
              close(ififft)
            endif  
!!
         endif
      endif

!! === \Sigma-Vxc to sham%rv_a_ohrs. Neighbor table sham%iv_a_oiaxs ====
!! ndimsig is the dimension of the self-energy. 
!! ndimsig = ldim means we use only projection onto MTO spaces even when PMT. !aug2012takao
!! --- New sigma inerpolation mode. Simplified.
!! OUTPUT in this block is used in
!!   call bloch 
!!  where we use sham%iv_a_oiaxs , sham%rv_a_ohrs, sham%ndhrs and so on.
!!  In future, we should make this block to a subroutine which shows clear i/o.
      sigmamode = mod(lrsig,10) .ne. 0
      if(procid==master.and.ix99==1) print *,'sigmamode=',sigmamode
      if(sigmamode.and.siginit .and. (.not.(llmfgw.and.jobgw==0)) ) then
         ndimsig= ldim  !if(mtosigmaonly()) mode. Dimension of sigm is the size of MTOs.
         if(procid==master) then
            write(6,*)' --- goto rdsigm2 initialization section ---'
            rsrnge = sham%rsrnge
            write(6,*)'rsrnge procid=',rsrnge,procid
            ifi = fopna('sigm',-1,4)
            read(ifi,err=9995,end=9995) nspsigm,ndimh_dummy,nk1,nk2,nk3
!! takao !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        laf=allocated(symops_af) !jun2015takao
        if(laf) nspsigm=2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            allocate( qsmesh2(3,nk1,nk2,nk3) )
            if(mtosigmaonly()) then
               allocate(sfz(nk1,nk2,nk3,ndimsig,ndimsig,nspsigm))
            else
               stop 'mtosigmaonly()=T is needed in the current version sep2012'
            endif
            rewind ifi
            if (mod(pwmode,10) ==0) then
               write(6,"(a,2i5)")" --- ldim(=dim of lmto)= ",ldim
            elseif(napwmx>0) then !if undefined skip this (null in m_hamble).
               write(6,"(a,2i5)")" --- ldim(=dim of lmto), napwmx= ",ldim,napwmx
            endif   
            if (ldim .ne. ndham) then !june2009 takao
               write(6,"(a,2i5)")" --- ndham (>=max of ham dimension) =",ndham
            endif
!! sfz is the self-energy roteted in the full BZ.
!!   ndham is removed sep2012 (read from sigma file )
            call rdsigm2(nbas,nspsigm,slat,sham,sbz,ifi, 
     &        nk1,nk2,nk3,ldim,qsmesh2,sfz,  qsfz,nqps, mtosigmaonly(),ndimsig) 
            call fclose(ifi)
c$$$!! --- write sig_fbz for SYMGRP e
c$$$            if(.not.cmdopt('--wsig_fbz',10,0,strn).and.oncewrite(10)) then
c$$$               write(6,"(a)")' ---> writing sigm_fbz.* which is in full bz'
c$$$               ifis2 = fopna('sigm_fbz',-1,4)
c$$$c$$$               if(pdos) then
c$$$c$$$                  read(ifis2) nspsigm,ndimsig,nk1,nk2,nk3
c$$$c$$$               else  
c$$$               write(ifis2) nspsigm,ndimsig,nk1,nk2,nk3,nk1*nk2*nk3,0,0,0
c$$$c$$$               endif  
c$$$               do is=1,nspsigm
c$$$                  do ik1=1,nk1
c$$$                     do ik2=1,nk2
c$$$                        do ik3=1,nk3
c$$$c$$$                           if(pdos) then
c$$$c$$$                              read(ifis2) qsmesh(1:3,ik1,ik2,ik3)
c$$$c$$$                              read(ifis2) sfz(ik1,ik2,ik3,1:ndimsig,1:ndimsig,is)
c$$$c$$$                           else
c$$$                           write(ifis2) qsmesh(1:3,ik1,ik2,ik3),is
c$$$                           write(ifis2) sfz(ik1,ik2,ik3,1:ndimsig,1:ndimsig,is)
c$$$c$$$                           endif
c$$$                        enddo
c$$$                     enddo
c$$$                  enddo
c$$$               enddo
c$$$               call fclose(ifis2)
c$$$            endif

!! mesh generator.  Unused yet...->in future, we will use this.
            allocate(npair(nbas,nbas),pos(3,nbas))
            call dcopy ( nbas * 3 , slat%rv_a_opos , 1 , pos , 1 )
            npairmx=nk1*nk2*nk3*2

            if (allocated(nlat)) deallocate(nlat)
            if (allocated(nqwgt)) deallocate(nqwgt)
            do !this loop is just to determined npairmx
              npairmx= npairmx + (nk1*nk2*nk3+1)*.5 !+1 added Feb2014 for 1x1x1
              allocate( nlat(3,npairmx,nbas,nbas), nqwgt(npairmx,nbas,nbas) )
              call gennlat(pos,nbas,plat,nk1,nk2,nk3,npairmx,ok,npair,nlat,nqwgt)
              if(ok) exit
              deallocate( nlat, nqwgt )
            enddo
            deallocate(pos,npair)

!! -- We use ndimsig=ldim. This means sigm is projected on MTOs.
!!    The main input of the subroutine seneinterp is sfz. 
!!    The seneinterp stores real-space sigma to sham%rv_a_ohrs (real space representation), 
!!    which is used in the following 'call bloch' to obtain Sigm-Vxc for any q point by interpolation.
!! NOTE: sigm file is read again within seneinterp-chksgr, and compared with FFT-sfz.
c            if(pdos) then
c               ifi = fopna('sigm_fbz',-1,4)
c            else
            ifi = fopna('sigm',-1,4)
c            endif   
            call seneinterp(nbas,nspsigm,ndimsig,slat,sham,sbz,ifi,rsrnge, !aug2012takao
     &        nk1,nk2,nk3,sfz)  !NOTE: output is stored in sham%rv_a_ohrs
            call fclose(ifi)
            deallocate(sfz)
            deallocate(qsmesh2,nlat,nqwgt)
         endif             !procid==master
#ifdef MPIK
!! Broadcast date stored in call seneinterp: probably we need to make it simplified... nov2012takao 
         if(procid==master) then
            nsizeohrs  = size( sham%rv_a_ohrs )
            nsizeoiaxs = size( sham%iv_a_oiaxs )
            nsizeontabs= size(sham%iv_a_ontabs)
            noqsig = size(sham%rv_p_oqsig)
         endif
         call mpibc1(nspsigm,1,2,mlog,'bndfp','nspsigm')
         call mpibc1(sham%ndhrs,1,2,mlog,'bndfp','ndhrs')
         call mpibc1(nsizeohrs,1,2,mlog,'bndfp','nsizeohrs')
         call mpibc1(nsizeoiaxs,1,2,mlog,'bndfp','nsizeoiaxs')
         call mpibc1(nsizeontabs,1,2,mlog,'bndfp','nsizeontabs')
         call mpibc1(noqsig,1,2,mlog,'bndfp','noqsig')
!! allocate and broadcast
         if(procid/=master) then
            if(allocated(sham%rv_a_ohrs)) deallocate(sham%rv_a_ohrs)
            allocate( sham%rv_a_ohrs(nsizeohrs) )
            if(allocated(sham%iv_a_oiaxs)) deallocate(sham%iv_a_oiaxs)
            allocate( sham%iv_a_oiaxs(nsizeoiaxs))       
            if(allocated(sham%iv_a_ontabs)) deallocate(sham%iv_a_ontabs)
            allocate( sham%iv_a_ontabs(nsizeontabs))       
            if(associated(sham%rv_p_oqsig)) deallocate(sham%rv_p_oqsig)
            allocate(sham%rv_p_oqsig(noqsig))
         endif   
         call mpibc1(sham%iv_a_ontabs,nsizeontabs,2,mlog,'bndfp','ontabs')
         call mpibc1(sham%iv_a_oiaxs ,nsizeoiaxs,2,mlog,'bndfp','oiaxs')
         call mpibc1(sham%rv_a_ohrs,  nsizeohrs, 4,mlog,'bndfp','ohrs')
         call mpibc1(sham%rv_p_oqsig, noqsig,4,mlog,'bndfp','noqsig')
         call mpi_barrier(MPI_COMM_WORLD,ierr) 
#endif
         siginit=.false. !this is equivalent to broad casting.
         if(procid==master) write(6,*)' --- end of rdsigm2 initialization section ---'
      endif    ! end of sigmamode


!! == GW driver ==
      if(llmfgw) then
         nchan=spot%nlma
         call sugw(ssite , sspec , slat , sham , nbas , ndham , spot%zv_a_osmpot 
     .     , vconst , lcplxp , sv_p_osig , sv_p_otau , sv_p_oppi , ppnl_rv 
     .     , vrmt , spotx , sv_p_osigx , sv_p_otaux , sv_p_oppix , jobgw 
     .     ,nk1,nk2,nk3,sigmamode,ndimsig, numprocs,procid,master) 
! NOTE: Sigma-Vxc in real space is give by sham%rv_a_ohrs.
         deallocate(spotx)
         call tcx('bndfp')
#if MPIK
         call mpi_barrier(MPI_COMM_WORLD,ierr) 
#endif
         call rx0('sugw mode')  !exit program here normally.
      endif


!! === do loop 2010 for iq and do 2005 isp ===
      if(debug) print *,' --- goto do loop 2010 ---'
      if(allocated(nevls)) deallocate(nevls)
      allocate(nevls(nkp,nspx))
      nevls=0
C Start k and isp loop (first loop in parallel mode)
C     h,s,evecs are dimensioned ndimh in this loop (may be q dependent)
C     evl is dimensioned evl(ndham,2)
C     evlall, w(owtkb), w(oswtk) are dimensioned (ndham,nsp,nkp)
#if MPIK
      call info0(30,1,0, ' ... Start MPI k-loop')
      sttime = MPI_WTIME()
      allocate (kpproc(0:numprocs), stat=ierr)
      call dstrbp(nkp,numprocs,1,kpproc(0))
      iqini = kpproc(procid)
      iqend = kpproc(procid+1)-1
#else
      iqini=1
      iqend=nkp
#endif
      initiqisp=.true.
      if(allocated(ndimhx_)) deallocate(ndimhx_)
      if(allocated(nev_)) deallocate(nev_)
      allocate(ndimhx_(nkp),nev_(nkp))
      ndimhx_=0
      nev_   =0
      if(pwmode>0.and.pwmode<10) call shortn3_initialize(qlat,rlatp,xmx2)

      if(sham%scaledsigma/=1d0.and.sigmamode) then
        write(6,*)' Scaled Sigma method: ScaledSigma=',sham%scaledsigma
      endif   
      do 2010 iq = iqini, iqend !This is a big iq loop
         qp=qplist(:,iq)
         if(debug) print *,' do 2010 iq procid=',iq,procid
#if MPIK
         if(debug) write(*,'(a,i3,5i5)')'debug iq=',procid,iq,iqini,iqend
         if (iq .eq. kpproc(procid)) then
            if (mlog) then
               call gettime(datim)
c               call awrit4(' bndfp '//datim//' Process %i of %i on '
c     .           //shortname(procid)(1:namelen(procid))//
c     .           ' starting k-points %i to %i',' ',256,lgunit(3),
c     .           procid,numprocs,kpproc(procid),kpproc(procid+1)-1)
               aaachar=' bndfp '//datim//' Process '//
     &           trim(i2char(procid))//' of '//
     &           trim(i2char(numprocs))//' on '//
     &           shortname(procid)(1:namelen(procid))//' starting k-points '//
     &           trim(i2char( kpproc( procid)))//' to '//
     &           trim(i2char( kpproc(procid+1)-1))
               write(6,"(a)")aaachar
            endif
         endif
#endif
         if(pwmode>0.and.pwmode<10) then
            ppin=matmul(transpose(plat),qp) !basis on the qlat coordinate. qp in Cartesian.

            call shortn3(rlatp,xmx2,ppin,noutmx, nout,nlatout)
            if(debug) then
               do iout=1,nout
                  write(*,"(a,3i5,f10.4,3f8.4)")'rrrrn1 =',nlatout(:,iout),
     &              sum(matmul(qlat(:,:),ppin+nlatout(:,iout))**2),
     &              matmul(qlat(:,:),ppin+nlatout(:,iout))
               enddo
            endif
         endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cC... range check qqq is in [ikc1/nk1, (ikc1+1)/nk1],... qpara is [0,1],[0,1],[0,1]
c      print *,' --- goto qparam--- '
c      call qparam(qp,qsmesh,nk1,nk2,nk3,plat, ikc,qpara,igshift)
c      ! eight corners -->  qsmesh(1:3,ikc(1,ix1,ix2,ix3),ikc(2,ix1,ix2,ix3),ikc(3,ix1,ix2,ix3))
c      call interparam(qpara, wgc)
c      cycle
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!! ... Loop over majority, minority spins
!!     In the SO case, still a loop isp=1..2
!!     isp=1 generates just the 11 block
c$$$        qpo=qp !moved here on mar2012 only for band plot
c$$$        if(.not.noshorbz()) call shorbz(qp,qp,qlat,plat)
         do 2005 isp = 1, nsp
            if(.not.(onesp==0.or.isp==onesp)) cycle
            if(debug) print *,' bndfp: do 2005 isp loop isp pwemax pwmode=',isp,pwemax,pwmode
!! takao. Do we need to do shorbz? 
!! -->furthermore, shorbz is probably wrong.--> see new algolism of shortn3 .
c            if(.not.noshorbz()) then
c              if (isp .eq. 1 .or. isp .eq. onesp) call shorbz(qp,qp,qlat,plat)
c            endif
            if(allocated(igv2x)) deallocate(igv2x)
!! === For this qp, G vectors for PW basis and hamiltonian dimension ===
            if (pwemax .gt. 0 .and. mod(pwmode,10) .gt. 0) then
               pwgmin = dsqrt(pwemin)
               pwgmax = dsqrt(pwemax)
               call pshpr(1)
               qqq=0d0          !call dpzero(xv,3)
               if (mod(pwmode/10,10) == 1) qqq=qp !call dpcopy(qp,xv,1,3,1d0)
               call gvlst2(alat,plat,qqq,0,0,0,pwgmin,pwgmax,0,0,napw,napw,dum,dum,dum,dum)
               call poppr
               allocate(igv2x(3,napw))
               allocate(kv_iv(3*napw))
               call pshpr(iprint()-10)
               call gvlst2(alat,plat,qqq,0,0,0,pwgmin,pwgmax,0,2,napw,napw,kv_iv,dum,dum,igv2x)
               if (pwmode<10) then
                  do iapw=1,napw
                     igv2x(:,iapw)=igv2x(:,iapw)+nlatout(:,1)
                  enddo
               endif
               deallocate(kv_iv)
               call poppr
               ndimh = ldim + napw !ndimh is iq-dependent.
               if (mod(pwmode,10) .eq. 2) ndimh = napw
               if (ndimh .gt. ndham) then
                  call fexit2(-1,111,'%N Exit -1 : BNDFP: '//
     .              'ndimh=%i exceeds ndham=%i.  Set '//
     .              'HAM_NPWPAD=(ndimh-ndham) or larger.',ndimh,ndham)
               endif
            else
               ndimh = ldim
               allocate(igv2x(1,1)) !dummy
            endif
            ndimhx = ndimh*nspc !this is iq-dependent.

!!   ... Make Hamiltonian and overlap matrices
            nqsig=sham%nqsig
            ispc = min(isp,nspc) !nspc=2 for non-colinier case
            if(isp==1) then
               if(allocated(hamm)) deallocate(hamm)
               if(allocated(ovlm)) deallocate(ovlm)
               nnspc=1
               if(nspc==2) nnspc=4
               allocate(hamm(ndimh,ndimh,nnspc),ovlm(ndimh,ndimh,nnspc))
            endif
            if(ispc==1) then
               hamm=0d0
               ovlm=0d0
            endif

            i = lrsig*10
            if(debug) print *,' i lrsig=',i,lrsig

C Commented by takao lwsig is not used now.. If ncessary, need to implement it again.
C C       lwsig=-1: sigm LDA -> orbital basis
C C       lwsig= 1: sigm orbital -> LDA basis
C C       lwsig= 2: sigm orbital -> LDA basis, high energy parts replaced
C C       lwsig= 3: Write evecs of LDA hamiltonian to file
C C       lwsig= 4: Write evecs of hamiltonian to file
C C       lwsig= 5: sigm orbital -> orbital, high energy parts replaced
C C       hambls returns with evecs,evals of the LDA hamiltonian in s,h
C C       Transform sigm to LDA basis: hambls returns sigm(LDA) in h


!! === Sigma-Vxc ===
!! Generate sene(Sigma-Vxc) for given sfz.
!! Determine interpolated self-energy sene at qp from sfz. 
!! sigmat = Sigma-Vxc is generated in a basis of ndimsig (usually MTOs only) 
!!     ... Bloch transform sigm(RS)-sigm(k). :RS means realspace
!! Main input  => sham%iv_a_oiaxs,sham%rv_a_ohrs
!! Main output => sene
            if(sigmamode) then
               allocate(sene(ndimsig,ndimsig))
               hreal = 0
               if (mod(i/10,10) >= 4) hreal = 1 !this is lrsig in hambles.
               mode = 100000 + 4000 + 40*(1-hreal) + 000 + 10 !hreal=1 -> sig(T) is real
               if(debug) print *,'bndfp: goto bloch i q nl plat mxorb=',lrsig,mode,qp,ndimsig,ival (sham%iv_a_ontabs, nbas + 1)
               if(isp>nspsigm) then
                  ispsigm = nspsigm
               else
                  ispsigm=isp
               endif
!! bloch: we need hrs(means real space) and iaxs(pair table) for FFT.
               call bloch (mode,qp,globalvariables%nl,plat,globalvariables%mxorb,sham%iv_a_oindxo,1
     &           , ival ( sham%iv_a_ontabs , nbas + 1 )
     &           , sham%iv_a_oiaxs , sham%rv_a_ohrs, sham%ndhrs, ispsigm, nspsigm , ndimsig , ndimsig 
     &           , 0 , ndimsig, 0 , ndimsig , 0 , sene, iwdummy, iwdummy )
               if(debug) print *,'end of bloch sene'
! takao expected to replace bloch to be a simple version seneinterp2. But not yet...
c     call seneinterp2(qp,nbas,nk1,nk2,nk3,plat, !npair,nlat,qwgt,npairmx,qsmesh,
c     i        sfz(:,:,:,:,:,isp),ldim,ndimsig, 
c     o        sene) 
            endif


c$$$C--- This works OK for lmf mode when noshrbz=.true. because no interpolation is needed.
c$$$C    Check comment out "sene=0" and "if(j1<=ldim.and.j2<=ldim) cycle" or not ---it depends on what you want to test.
c$$$c          allocate(senex(ndimh,ndimh))
c$$$          allocate(sene(ndimh,ndimh))
c$$$c          senex=0d0
c$$$          write(6,"(a,3f13.5)")'    qp=',qp
c$$$          do ik1=1,nk1
c$$$          do ik2=1,nk2
c$$$          do ik3=1,nk3
c$$$             write(6,"(a,3f13.5)")' qsmesh=',qsmesh(1:3,ik1,ik2,ik3)
c$$$             if(sum(abs(qp-qsmesh(1:3,ik1,ik2,ik3)))<1d-6) goto 1013
c$$$          enddo
c$$$          enddo
c$$$          enddo
c$$$          stop 'xxxxx not qp=qsmesh xxxx'
c$$$ 1013    continue
c$$$          do j1=1,ndimh !ldim+napw
c$$$          do j2=1,ndimh !ldim+napw
c$$$c            if(j1<=ldim.and.j2<=ldim) cycle
c$$$c            senex(j1,j2)= sfz(ik1,ik2,ik3,j1,j2,isp)
c$$$            sene(j1,j2)= sfz(ik1,ik2,ik3,j1,j2,isp)
c$$$          enddo
c$$$          enddo
c$$$c          print *,'ldim:abs and diff sene-senex=',sum(abs(sene(1:ldim,1:ldim)))
c$$$c     &      ,sum(abs(sene(1:ldim,1:ldim)-senex(1:ldim,1g:ldim)))
c$$$c          print *,'ndimh:abs and diff sene-senex=',sum(abs(sene))
c$$$c     &      ,sum(abs(sene-senex))
c$$$c          print *, 'end of test mode of sene newsigmasw=',newsigmasw()
c$$$c          deallocate(senex)
c$$$c          stop 'zzzzzzzzzzzzzzzzzzzzzzzzzzzzzz'
c
c$$$C--- This works OK for lmf mode when noshrbz=.true. for PMT
c$$$C    Check comment out "sene=0" and "if(j1<=ldim.and.j2<=ldim) cycle" or not ---it depends on what you want to test.
c$$$c          allocate(senex(ndimh,ndimh))
c$$$          allocate(sene(ndimh,ndimh))
c$$$c          senex=0d0
c$$$          write(6,"(a,3f13.5)")'    qp=',qp
c$$$          do ik1=1,nqps
c$$$             write(6,"(a,i3,3f13.5)")' ik1 qsfz=',ik1,qsfz(1:3,ik1)
c$$$             iq1=ik1
c$$$            ii1=mod(iq1,2)+1
c$$$            ii2=mod(iq1/2,2)+1
c$$$            ii4=mod(iq1/4,2)+1
c$$$            print *,' iq1 ii4,ii2,ii1 =',iq1,ii4,ii2,ii1
c$$$             print *,' xxxxxxxxxxx this shows shorbz problem... xxxxxxxxx'
c$$$             if(sum(abs(qp-qsfz(1:3,ik1)))<1d-6) goto 1013
c$$$          enddo
c$$$          stop 'xxxxx not qp=qsmesh xxxx'
c$$$ 1013    continue
c$$$          do j1=1,ndimh !ldim+napw
c$$$          do j2=1,ndimh !ldim+napw
c$$$c            if(j1<=ldim.and.j2<=ldim) cycle
c$$$c            senex(j1,j2)= sfz(ik1,ik2,ik3,j1,j2,isp)
c$$$            sene(j1,j2)= sfz(ii4,ii2,ii1,j1,j2,isp)
c$$$          enddo
c$$$          enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

Ckino hambl calls augmbl.  augmbl calculates C_kl(=array b in the program)
Ckino finally makes F~F~=F0F0+(F1F1-F2F2), which is overlap matrix, s.
Ckino Note that F2=Hankel head at a site + Hankel tail contributions from the other site.

!! == Set up Hamiltonian. ==============
!!    Hamiltonian: hamm(1:ndimh,1:ndimh,ispc). ispc=3 means off-diagonal section when SO=1.
!!    Overlap matrix: ovlm 
!!    senex: Sigma-Vxc
!! ==========================================
            jx = 0
            if (lso == 1) jx = 10
            call hambl ( jx , nbas , ssite , sspec , slat , sham , isp ,
     .        qp , k1 , k2 , k3 , spot%zv_a_osmpot , vconst , sv_p_osig , sv_p_otau 
     .        , sv_p_oppi , lcplxp , alfa , ndimh , napw , igv2x , hamm ( 1 
     .        , 1 , ispc ) , ovlm ( 1 , 1 , ispc ) , hamm ( 1 , 1 , 3 ) )
            if(sigmamode) then
               allocate( ovlmtoi(ndimsig,ndimsig),ovliovl(ndimsig,ndimh))
               ovlmtoi = ovlm(1:ndimsig,1:ndimsig,ispc)
               call matcinv(ndimsig,ovlmtoi)
               ovliovl = matmul(ovlmtoi,ovlm(1:ndimsig,1:ndimh,ispc))
               deallocate(ovlmtoi)
               allocate(senex(ndimh,ndimh))
               senex = matmul(transpose(dconjg(ovliovl)), matmul(sene,ovliovl))
               deallocate(ovliovl,sene)
               hamm(:,:,ispc) = hamm(:,:,ispc) + sham%scaledsigma * senex
               deallocate(senex)
            endif



C   ... In the noncollinear case (spin-orbit SO=1), this marks the end of looping over
C       separate spins.
C       For the rest of the isp loop,
C       jsp=isp in the collinear case; jsp=1 in the noncollinear
C       Thus jsp should be used in place of isp
C       isp serves as a flag for the noncollinear case
            if (ispc .ne. nspc) goto 30
            jsp = isp
            if (ispc .eq. 2) jsp = 1

!! === Diagonalize and add to density ===
            if(iprint()>=30.and.ix99==1) 
     &        write(*,'(" bndfp: kpt ",i5," of ",i5, " k=",3f8.4," ndimh = nmto+napw = ",3i5,f13.5)')
     &        iq,nkp,qp,ndimh,ndimh-napw,napw
c     &   iq,nkp,qp,trim(i2char(ndimh))//' '//trim(i2char(napw))//' '//trim(i2char(ndimh-napw)),xdatt(iq)
            nmx = min(nevmx,ndimhx) !I think nmx should be given in a better manner.
!!  Calculate all eigenvalues for plotting mode, pdos, tdos modes.
c            if (nfbn(1) .gt. 0 .or.tdos) then ! or.pdos
            if (tdos) then 
               nmx = ndimhx
               efmax = 99999
            endif
            if (allocated(t_zv)) deallocate(t_zv)
            allocate(t_zv(ndimhx,nmx))
            if(debug) write(*,*)'kino, allocate t_zv' ,ndimhx,nmx,ndimhx

!! == lpdiag=0 only. ===
            if(nspc .eq. 2) then

c$$$ccccccccccccccccccccccccccccccccccccccccccccccc
c$$$c              hamm(:,:,2) = hamm(:,:,1)
c$$$c              ovlm(:,:,2) = ovlm(:,:,1)
c$$$              write(6,*)'iiix'
c$$$              do i=1,ndimh
c$$$                do j=1,ndimh
c$$$                  if(abs(ovlm(i,j,2)-ovlm(i,j,1))>1d-8) then
c$$$                  write(*,"('iiiiix2',3f8.3,2i5,2d10.2,x,2d10.2,x,2d10.2)")qp,i,j,ovlm(i,j,1),ovlm(i,j,2), ovlm(i,j,2)-ovlm(i,j,1)
c$$$                  endif
c$$$                enddo
c$$$              enddo
c$$$              write(6,*)'iiixh'
c$$$              do i=1,ndimh
c$$$                do j=1,ndimh
c$$$                  if(abs(hamm(i,j,2)-hamm(i,j,1))>1d-8) then
c$$$                  write(*,"('iiiiix2h',3f8.3,2i5,2d10.2,x,2d10.2,x,2d10.2)")qp,i,j,hamm(i,j,1),hamm(i,j,2), hamm(i,j,2)-hamm(i,j,1)
c$$$                  endif
c$$$                enddo
c$$$              enddo
c$$$cccccccccccccccccccccccccccccccccccccccccccccccc

c     call sopert ( 0 , ndimh , nspc , ww_rv , hamm , hamm )
               call sopert2 ( 0 , ndimh , nspc, hamm , hamm )
               call sopert2 ( 0 , ndimh , nspc,  ovlm , ovlm )
            endif
            if(debug) write(*,*) 'allocate ww_rv', 11*ndimhx

!! === diagonalization ===
!! It is better to clean up lwtkb-algolism in this routine. !!!!
!! savez save eigenfuncitons for second pass.(second-pass procedure should be removed in future).
            savez=.false.
            getz=.false.

            if(gotosecondpass) then
               if(ix99==1) then
                  savez=.true.  ;  getz=.false.
               else
                  savez=.false. ;  getz=.true.
               endif
               if(initiqisp) then
#if (MPI|MPIK)
!!  for parallel write, and eigze is temporaly used only in this subroutine. 
!!     (gfortran v4.7 add _{procid} at the end of the filename
!!      then how can I share the same file among the different processes? )
                  !ifig = fopna('eigze_'//trim(strprocid),-1,4)
                  !call fclose(ifig)
                  !ifig = fopna('eigze_'//trim(strprocid),-1,4)
                  inquire(file='eigze_'//trim(strprocid),opened=nexist,number=ifig)
                  if(nexist) close(ifig)
                  ifig = ifile_handle()
                  open(ifig,file='eigze_'//trim(strprocid),form='unformatted')
#else
                  !ifig = fopna('eigze',-1,4)
                  !call fclose(ifig)
                  !ifig = fopna('eigze',-1,4)
                  inquire(file='eigze',opened=nexist,number=ifig)
                  if(nexist) close(ifig)
                  ifig = ifile_handle()
                  open(ifig,file='eigze',form='unformatted')
#endif
                  initiqisp=.false.
               endif  
            endif
            
            if(PROCARon) then   !!2014mar for band weight mode
               nmx=ndimhx
               nev=ndimhx
               savez=.false.
               getz=.false.
               nmx = ndimhx
               if(allocated(t_zv)) deallocate(t_zv)
               allocate(t_zv(ndimhx,ndimhx))
            endif

!! == Diagonalize Hamiltonian ==
!! ndimhx: dimension of Hamitonian
!! hamm:Hamiltonian, ovlm: overlap matrix
!! t_zv:eigenfunciton. evl: eigenvalue.
!! nmx: input, number of requested eigenvalues(funcitons)
!! nev: out number of obtained eigenfvalues(funcitons)
!! ---
            if(getz) then   ! second pass
              read(ifig) nev,nmx
              read(ifig) evl(1:nev,jsp)
              read(ifig) t_zv(1:ndimhx,1:nmx)
            else            ! 1st pass

!! save hammovlm jan2018kaist
              if(writeham) then
                write(6,"(a,3f9.5)") "Hamiltonian: Writing hamm and ovlm for qp= ",qp
                write(ifih) qp,ndimhx,lso,epsovl,jsp
                if(lso==1) then !L.S case nspc=2 spin block
                  write(ifih) ovlm(1:ndimhx,1:ndimhx,1:4) 
                  write(ifih) hamm(1:ndimhx,1:ndimhx,1:4)
                else  !nspc=1 spin diagonal case
                  write(ifih) ovlm(1:ndimhx,1:ndimhx,1)
                  write(ifih) hamm(1:ndimhx,1:ndimhx,1)
                endif
              endif
!!
              if(epsovl<1.000001d-14) then
                call zhev_tk2( ndimhx , hamm , ovlm , nmx , nev,
     .           evl(1 , jsp ), t_zv, iq==1,savez,getz,ifig) !stock eigenfunctions z
              else              !nov2015 recover epsovl. Remove space by overlap matrix diagonalization.
                call zhev_tk3( ndimhx , hamm , ovlm , nmx , nev,
     .           evl(1 , jsp ), t_zv, iq==1,savez,getz,ifig,epsovl) !stock eigenfunctions z
              endif
              if(savez) then
                write(ifig) nev,nmx
                write(ifig) evl(1:nev,jsp)
                write(ifig) t_zv(1:ndimhx,1:nmx)
              endif  
#ifdef MPIK            
              if(writeham.and.procid==master) then
                call prtev(t_zv, ndimhx , evl(1, jsp ) , nmx , efmax , nev ) !efmax is required?
              endif
#endif
            endif

!! It is inefficient to calculate all eigenfunctions and eigenvalues for molecule cases.
!! nevl=nev: number of eigenvalues (if nmx=0, nev=n. See description in zhev_tk2.)
            nevl=nev            !ndimhx    this OK? why nevl and nev???
            evl(nevl+1:ndhamx,jsp)=1d99 !to skip these data
            nevls(iq,jsp) = nevl  !nov2014 isp and jsp is confusing...
            if(debug) print *,'iq jsp nevls=',iq,jsp,nevls(iq,jsp) !nov201
#ifndef MPIK            
            call prtev(t_zv, ndimhx , evl ( 1 , jsp ) , nmx , efmax , nev ) !efmax is required?
#else
            if (procid.eq.master.and. epsovl>=1.000001d-14.and.plbnd/=0) then
              write(6,"('      : ndimhx=',i5,' --> nev=',i5' by HAM_OVEPS ',d11.2)") ndimhx,nev,epsovl
            endif  
#endif
!! save to evall
            if(debug) print *,' nmx nev plbnd=',nmx,nev,plbnd
            evlall(1:ndhamx,jsp,iq) = evl(1:ndhamx,jsp)


!!   ... Save data for this qp to moments file
            if (plbnd .eq. 0) then
!!   --- Orbital magnetic moment (requires lso) ---
               if (lso .ne. 0 .and. lwtkb .ne. -1) then
                  if (lwtkb .eq. 0) call
     .              rx('metal weights required for orb. moment calculation')
                  allocate(auso_zv(nlmax*ndhamx*3*nsp*nbas))
                  if(debug)   write(*,*) 'allocate auso_zv',nlmax*ndhamx*3*nsp*nbas
                  auso_zv(:)=0.0d0
                  call makusq ( 0 , ssite , sspec , slat , sham , nbas , nbas , 
     .              0 , nlmax , ndham , ndimh , napw , igv2x , nev , nsp , 
     .              nspc , jsp , 1 , qp , t_zv , ppnl_rv , auso_zv )
                  call mkorbm ( ssite , sspec , jsp , nsp , nspc , nlmax , ndham
     .              , nev , rv_a_owtkb , iq , nbas , ppnl_rv , auso_zv , nl , nkp 
     .              , orbtm_rv )
                  deallocate(auso_zv)
               endif

!! Mulliken analysis
               if (lwtkb .ne. -1) then
                  if (mullwrite) then
                        call rxx(nev.ne.ndimhx,'Mulliken requires all eigenvectors')
                        call rxx(nevmx.ne.ndhamx,'Mulliken requires nevmx=nhdamx')
                        allocate(doswt_rv(nchan*nevmx*nspc))
                        doswt_rv(:)=0.0d0
                        call mullmf(nbas , ssite , sspec , sham%iv_a_oindxo , t_zv , ndimh 
     .                    , nspc , iq , jsp , moddos , nsites , lsites , lmxch , nchan 
     .                    , chan_iv , lmdim , ndham , doswt_rv )
                        call iomomn(.true.,2,.false.,1,nspc,1,1,i)
                        i1 = iomoms ( - nfilem , nl , nsp , nspc , nkp , ndimh , i , 
     .                    nspc , 1 , 1 , ndhamx , nevmx , nchan , nchan , nev , evl ( 1 
     .                    , jsp ) , 0d0 , doswt_rv , 0d0 , 0d0 , 0d0 )
                        if (allocated(doswt_rv)) deallocate(doswt_rv)
                  endif
               endif   

!! Make new density matrix dmatu for LDA+U (implementation of LDA+U is a little complicated).
               if (lwtkb .ne. -1) then
                  if (nlibu .gt. 0 .and. nev .gt. 0) then
                     if (lwtkb .eq. 0)
     .                 call rx('metal weights required for LDA+U calculation')
                     nl = int(sctrl%nl)
                     nlmax = nl*nl
                     allocate(ausu_rv(2*nlmax*ndhamx*3*nsp*nbas))
                     if(debug) write(*,*)'allocate ausu_rv',2*nlmax*ndhamx*3*nsp*nbas
                     ausu_rv(:)=0.0d0
                     call makusq ( 0 , ssite , sspec , slat , sham , nbas , nbas , 
     .                 0 , nlmax , ndham , ndimh , napw , igv2x , nev , nsp , 
     .                 nspc , jsp , 1 , qp , t_zv , ppnl_rv , ausu_rv )
                     call mkdmtu ( ssite , sspec , rv_a_owtkb , jsp , iq , nsp , nspc 
     .                 , ndham , nlmax , nbas , nev , ppnl_rv , ausu_rv , dmatu , nlibu 
     .                 , lmaxu , lldau )
                     if (allocated(ausu_rv)) deallocate(ausu_rv)
                  endif
               endif

!! Core-level spectroscopy
               if (icls .ne. 0 .and. lwtkb .ne. -1) then
                  call rxx(nspc.ne.1,'CLS not implemented in noncoll case')
#if MPI
                  call rx('CLS only k-parallel')
#endif
                  call makusq ( 0 , ssite , sspec , slat , sham , nbas , nsites 
     .              , isite , nlmax , ndham , ndimh , napw , igv2x , nev , 
     .              nsp , nspc , jsp , iq , qp , t_zv , ppnl_rv , ausc_zv   )
               elseif (lwtkb.ne.-1 .and. (lrout.ne.0 .or. lmet.eq.4)) then

!! Density mode. wtkb is the weight for iq,ib. With this setting, 
!! we can plot |\psi_{iq,ib}(\bfr)|**2 in the manner of electron density
!! lmf --density,iq=12,ib=3,4,5   (here ib=3,4,5 are superposed).
                 if(cmdopt('--density',9,0,strn)) then
                   iqread=0
                   iqindex=index(strn(10:),'iq=')+2
                   if(iqindex==2) then 
                     iqread=-999
                   else
                     read(strn(10+iqindex:),*) iqread
                   endif
                   ibindex=index(strn(10:),'ib=')+2
                   if(ibindex==2) then
                     ibread=-999
                   else
                     do nibread=1,100
                       read(strn(10+ibindex:),*,err=2019,end=2019) ibread(1:nibread)
c                       print *,'xxx ibread=',ibread(1:nibread)
                     enddo
 2019                continue
                     nibread=nibread-1
                   endif  
                   if(iqread>0) then
                     if(maxval(ibread(1:nibread))<=0) call rx('--density mode: wrong ib=foobar. Try,e.g. --density,iq=12,ib=5')
                     write(6,"('--density bandmode: psi**2 for iq=',i5,' ib=',255i5)") iqread,ibread(1:nibread)
                     rv_a_owtkb(:,:,:)=0d0
                     do ib=1,nibread
                       rv_a_owtkb(ibread(ib),:,iqread)=1d0
                     enddo
                     do ib=1,ndimh
                       if(abs(rv_a_owtkb(ib,jsp,iq))>1d-3) write(6,"('ib wtkb=',i5,2f13.6)") ib,rv_a_owtkb(ib,jsp,iq)
                     enddo
                   endif
                 endif   

!! Accumulate output density and sampling DOS
!! Even if no output rho, still call addrbl to make DOS when lmet=4
                  if(debug) print *,'goto addrbl nevl=',nevl
                  call addrbl ( ssite , sspec , slat , sham , jsp , nsp , nspc
     .              , qp , sbz%rv_a_owtkp , ndham , ndimh , napw , igv2x , lmet , lrout
     .              , lwtkb , rv_a_owtkb , lswtk , rv_a_oswtk , iq , lfrce , ldos
     .              , lekkl , k1 , k2 , k3 , spot%zv_a_osmpot , vconst , sv_p_osig , sv_p_otau
     .              , sv_p_oppi , lcplxp , numq , qval - qbg , t_zv , evl , nev
     .              , ef0 , def , esmear , emin , emax , ndos , dos_rv , srout_zv 
     .              , sumqv , sumev , sv_p_oqkkl , sv_p_oeqkkl , frc )
               endif

            endif


!!   === Save evals for this qp ===
c            elseif (plbnd .ne. 0) then
c$$$#ifndef MPIK
c$$$               if(debug) write(*,*)  'plbnd/=0 (noMPI)', ndham,jsp 
c$$$               if (nfbn(1)>0 .and. n_listwf>0) then
c$$$                  call rsibl_ev( ssite,sspec,slat,nbas,isp,qp,iq,ndimh,nspc
c$$$     .              ,  napw,igv2x,sham%iv_a_oindxo,nev,t_zv,k1,k2,k3
c$$$     .              ,  n_listwf,listwf    )
c$$$               endif
c$$$               if(debug) print *,'dddddd111'
c$$$               i = nsp
c$$$               if (onesp .ne. 0 .or. nspc .eq. 2) i = 1
c$$$c              print *,' eee1',nevl,jsp,i
c$$$c              print *,' eee2',allocated(evl)
c$$$c              print *,' eee3',evl(1:2,jsp)
c$$$c              print *,' eee3',evl(1:nevl,jsp)
c$$$               call suqlsw2(nevl,jsp,i,evl(1,jsp),qp) !takao
c$$$               if(debug) print *,'ddddd2222'
c$$$               if (nfbn(1) .ne. 0) then
c$$$                  if (ndimhx .ne. nevl)
c$$$     .              call rx('color weights not implemented when '//
c$$$     .              'nevl < hamiltonian dimension')
c$$$                  call suqlse ( ndimhx , jsp , i , ndimhx , 1 , nfbn , ifbls_iv 
c$$$     .              , ndhamx , t_zv , evl )
c$$$               endif
c$$$               if(debug) print *,'ddddd333'
c$$$               if (nfbn(2) .ne. 0) then
c$$$                  if (ndimhx .ne. nevl)
c$$$     .              call rx('color weights not implemented when '//
c$$$     .              'nevl < hamiltonian dimension')
c$$$                  call suqlse ( ndimhx , jsp , i , ndimhx , 2 , nfbn , ifbls_iv 
c$$$     .              , ndhamx , t_zv , evl )
c$$$               endif
c$$$               if(debug) print *,'ddddd4444'
c$$$#else
c               if (procid.eq.master) then
c               write(6,*)  'case plbnd.ne.0 (MPIK),procid.iq=',procid,iq
c               if(debug) write(*,*)  'nnn ', ndham,jsp 
c               endif
c              evlall(1:ndhamx,jsp,iq)= evl(1:ndhamx,jsp)
c$$$#endif
c            endif


!! Decompose to alm components. PROCAR generation.
!! ??? no nsites is used? ! Compare following codes with those of mkpdos.F (PDOS).
            if(PROCARon) then
               nlmax = nl*nl
               allocate( auspp(nlmax,ndhamx,3,nsp,nbas) )
               initbd=.false.
               qold=qp
               auspp = 0d0
               call makusq( 1 , ssite , sspec , slat , sham , nbas , nbas , 
     .           0 , nlmax , ndham , ndimh , napw , igv2x , nev , nsp , 
     .           nspc , jsp , 1 , qp , t_zv , ppnl_rv , auspp )
               if(spinoffdiag) then
                  jspini=1
                  jspend=2
               else
                  jspini=isp
                  jspend=isp
               endif
               do jspp= jspini, jspend
                  if(jspp==1) iprocar=iprocar1
                  if(jspp==2) iprocar=iprocar2
                  if(debug) print *,'iiiii =',iprocar1,iprocar2,jspp,iprocar,ef0
                  if(debug) print *,'nnnnn=',nlmax,ndham,nspc,nsp,nbas
                  ccc="ion        s       py       pz       px      dxy      dyz      dz2      dxz   dx2-y2"//
     &                 "      f-3      f-2      f-1       f0       f1       f2       f3"//
     &                 "                                                                                 tot"
                  write(iprocar,*)
                  write(iprocar,*)
                  write(iprocar,"('k-point ',i4,' :    ',3f11.8,'     weight = -------  : x =',f15.8)")iq,qp,xdatt(iq)
                  write(iprocar,*)
               do iband = 1, nev
                  write(iprocar,*)
                  write(iprocar,"('band ',i3,' # energy ',f13.8,' # occ. -----' )")iband,(evl(iband,jsp)-ef0)*rydberg
                  write(iprocar,*)
                  dwgtt=0d0
                  do ib = 1, nbas
                     is  = ssite(ib)%spec
                     ilm = 0
                     dwgt=0d0
                     do  l = 0, sspec(is)%lmxa
                        do  m = -l, l
                           ilm = ilm+1
                           au = auspp(ilm,iband,1,jspp,ib)
                           as = auspp(ilm,iband,2,jspp,ib)
                           az = auspp(ilm,iband,3,jspp,ib)
                           s11 = dconjg(au)*au*ppnl_rv(2,l+1,jspp,ib)
                           s22 = dconjg(as)*as*ppnl_rv(7,l+1,jspp,ib)
                           s33 = dconjg(az)*az*ppnl_rv(8,l+1,jspp,ib)
                           s12 = 0
                           s13 = 2*dconjg(au)*az*ppnl_rv(9,l+1,jspp,ib)
                           s23 = 2*dconjg(as)*az*ppnl_rv(10,l+1,jspp,ib)
                           dwgt(ilm)= s11+s22+s33 + s12+s13+s23
                        enddo
                     enddo
                     dwgtt = dwgtt + dwgt(1:ilm)
                     if(ib==1)  write(iprocar,"(a)") trim(ccc)
                     write(iprocar,"(i3,100(x,f8.5))")ib,(dwgt(i),i=1,nchanp),sum(dwgt)
                     if(ib==nbas) write(iprocar,"('tot',100(x,f8.5))")(dwgtt(i),i=1,nchanp),sum(dwgtt)
                     if(fullmesh) dwgtall(1:nchanp,ib,iband,jspp,iq) = dwgt(1:nchanp)
                  enddo
               enddo
               enddo
               if(debug) print *,'mmmmmmmm procid sum dwgt check=',procid,iq,jsp,sum(dwgtall),nev
               deallocate( auspp )
            endif 
            if(allocated(t_zv)) deallocate(t_zv)
!! == end loop over isp (main loop in parallel mode)==
            if(debug) print *,'bbbbbbbb before 30'
 30         continue
            deallocate(igv2x)
 2005    continue
         deallocate(hamm,ovlm)
         ndimhx_(iq)=ndimhx
         nev_(iq)   =nev
 2010 continue
         if(debug) print *, ' ---- end of do 2010 ---- ',procid
#if MPIK
      call mpibc2(ndimhx_,nkp,2,mlog,'bndfp','ndimhx_') !all reduce is not necessary.
      call mpibc2(nev_,  nkp,2,mlog,'bndfp','nev_')
#endif
      nevmin=minval(nev_(1:nkp))

      
c     this inquire is moved above (right before opening ifig files) !may2019
c     inquire(unit=ifig,opened=nexist)
c      if(nexist) close(ifig)
!! ========== end of iq loop ====================================================
      if(debug) print *, 'tttttttttt 9999999999999999999 3',procid
      if(debug.and.PROCARon) print *,'mmmm procid sum dwgt check=',procid,sum(dwgtall)



!! boltztrap data 
      if( cmdopt('--boltztrap',11,0,strn) .and. procid==master) then
        ifi=ifile_handle()
        open(unit=ifi,file='boltztrap.dat')
        write(ifi,"(3d16.8, ' # (Ry) Fermi energy given by lmf')") eferm
        write(ifi,"(d16.8, ' # alat (a.u.)')") alat
        write(ifi,"(3d16.8, ' # plat1  ')") plat(:,1)
        write(ifi,"(3d16.8, ' # plat2  ')") plat(:,2)
        write(ifi,"(3d16.8, ' # plat3  ')") plat(:,3)
        write(ifi,"(3d16.8, ' # qlat1  ')") qlat(:,1)
        write(ifi,"(3d16.8, ' # qlat2  ')") qlat(:,2)
        write(ifi,"(3d16.8, ' # qlat3  ')") qlat(:,3)
        do iq=1,nkp
        do jsp=1,nspx         !=isp.  nspx=1 for so=1
!! true q= 2pi/alat *qplist(:,iq) in cartesian
        do iband=1,ndhamx
          if(evlall(iband,jsp,iq)>1d98) cycle !evall=1d99 is dummy
            write(ifi,"('eee=',3i5,3f9.4,d13.5)") iq,jsp,iband,qplist(:,iq),evlall(iband,jsp,iq)
        enddo
        enddo
        enddo
        close(ifi)
c        call rx0(' --bolttrap mode: boltztrap.dat is generated')
      endif
      if( cmdopt('--boltztrap',11,0,strn)) then
         call rx0(' --bolttrap mode: boltztrap.dat is generated')
      endif
      
!! === pdos mode (--mkprocar and --fullmesh). ===
      if(PROCARon.and.fullmesh) then
#ifdef MPIK
         call mpi_barrier(MPI_comm_world,ierr)  
         call xmpbnd(kpproc,ndhamx,nkp,nspx,evlall) !ndhamx,nspx nov2014 bugfix
         call xmpbnd(kpproc,nbas*nchanp*ndhamx,nkp,nspx,dwgtall) 
#endif
         if(procid==master) then  
            if(debug) print *,'mmmmm tot procid sum dwgt check=',procid,sum(dwgtall)
!! tetrahedron without symmetry
            allocate(idtete(0:4,6*nkp),ipqe(nkk1,nkk2,nkk3))
            iq=0
            do ik3 = 1, nkk3 
               do ik2 = 1, nkk2
                  do ik1 = 1, nkk1
                     iq = iq+1
                     ipqe(ik1,ik2,ik3)=iq
                     qx = (ik1-1)*qlat(:,1)/nkk1 + (ik2-1)*qlat(:,2)/nkk2+ (ik3-1)*qlat(:,3)/nkk3
                     if(abs(sum(qx-qplist(:,iq)))>1d-6) then !sanity check
                        call rx("bndfp: qx/=qplist something strang")
                     endif  
                  enddo
               enddo
            enddo
            call tetirr(qlat, nkk1,nkk2,nkk3, ipqe, ntete,idtete) !tetrahedron
            write(*,*)" ntete 6*nkk1*nkk2*nkk3 nkp=",ntete,6*nkk1*nkk2*nkk3,nkp
!! doswt(1:nchanp=16,ndhamx,nsp,nkp)
!! dwgtall(1:16,nbas,ndhamx,nsp,nkp)
!!   If we will have a memory problem for dwgtall, we may store dwgtall(i,...) for a file.
!!   and repeat nchanp times.
c            eminp=-25.0/rydberg
c            emaxp=30.0/rydberg
c            ndos=5500
            write(6,"(' pdosmode: ndhamx nsp nspx =',4i7)") ndhamx, nsp, nspx,nevmin !nevmx
c            print *,' pdosmode: nchanp nbas ndos emin emax ef0=', nchanp,nbas,ndos, emin,emax,ef0
            write(6,"(' pdosmode: nchanp nbas emin emax ef0=',2i7,3f12.5)") nchanp,nbas,emin,emax,ef0 !mar2015
            ifip= fopna('pdosdata',-1,4)
            write(ifip) ndhamx,nsp,nspx,nevmin,nchanp,nbas,nkk1,nkk2,nkk3,ntete,nkp !ndos,nkp mar2015
            write(ifip) idtete,ipqe ! allocate(idtete(0:4,6*nkp),ipqe(nkk1,nkk2,nkk3))
            write(ifip) evlall  !  allocate(evlall(ndhamx,nspx,nkp))
            write(ifip) dwgtall !  allocate(dwgtall(nchanp,nbas,ndhamx,nspx,nkp))
            write(ifip) ef0 !eminp,emaxp,ef0
            call fclose(ifip)
            deallocate(idtete,ipqe) !remove this if you uncomment above lines.
         endif  
         if(writeham) call rx0('Done --writeham: --mkprocar --fullmesh are needed now. HamiltonianMTO* genereted')
         call rx0('Done pdos: --mkprocar & --fullmesh. Check by "grep k-point PROCAR.*.*"')
      endif


!! for MPI, collect data of nodes see xmpbnd.
#if MPIK
      entime = MPI_WTIME()
      call info2(30,0,0, ' ... Done MPI k-loop: %;1d seconds elapsed',entime-sttime,0)
      call info0(20,0,-1,' ... Sharing data between processes...')
      sttime = MPI_WTIME()
      call xmpbnd(kpproc,ndhamx,nkp,nspx,evlall) !ndhamx,nspx nov2014 bugfix
      if (icls .ne. 0 .and. lwtkb .ne. -1) then
         call rxx(nspc.ne.1,'CLS not implemented in noncoll case')
         call xmpbnd ( kpproc , nlmax * ndhamx * 3 * nsites , nkp , nsp, ausc_zv )
      endif
      call mpibc2(nevls,nkp*nspx,2,mlog,'bndfp','nevls')
#endif
c#ifdef KINODEBUG
c      call mpi_barrier(MPI_comm_world,ierr)
c      write(*,*)'kino afeter mpibc2 nevls'
c#endif

!!  Allreduce density-related quantities
      if(debug) print *,'goto all reduce'
      if (lrout .ne. 0) then
         call mpibc2(sumqv,6,4,mlog,'bndfp','sumqv')
         call mpibc2(sumev,6,4,mlog,'bndfp','sumev')
         call mpibc2 ( srout_zv , k1 * k2 * k3 * nsp * numq , 6 , mlog, 'bndfp' , 'smrho' )
         if (lswtk .eq. 1) then
            call mpibc2 ( rv_a_oswtk , ndhamx * nkp , 4 , mlog , 'bndfp' , 'swtk' )
         endif
C     Allreduce qkkl
         nkaph = globalvariables%nkaph
         do  ib = 1, nbas
            is = int(ssite(ib)%spec)
            lmxa=sspec(is)%lmxa
            lmxh=sspec(is)%lmxb
            kmax=sspec(is)%kmxt
            nlma = (lmxa+1)**2
            nlmh = (lmxh+1)**2
            if (lmxa .gt. -1) then
               nelt(1) = (kmax+1)*(kmax+1)*nlma*nlma
               nelt(2) = (kmax+1)*nkaph*nlma*nlmh
               nelt(3) = nkaph*nkaph*nlmh*nlmh
               do  i = 1, 3
                  call mpibc2 ( sv_p_oqkkl( i , ib )%v , nelt ( i ) * numq * nsp 
     .              * nspc , 4 , mlog , 'bndfp' , 'qkkl' )
                  if ( lekkl .eq. 1 ) call mpibc2 ( sv_p_oeqkkl( i , ib )%v , nelt 
     .              ( i ) * numq * nsp * nspc , 4 , mlog , 'bndfp' , 'eqkkl' )
               enddo
            endif
         enddo
!!    Allreduce DOS, forces, dmatu, orbtm
         if ( ndos .gt. 0 ) call mpibc2 ( dos_rv , ndos * 2 * nsp , 4 , mlog , 'bndfp' , 'dos' )
         if (lfrce .ne. 0)  call mpibc2(frc,3*nbas*numq,4,mlog,'bndfp','frc')
         if (nlibu .gt. 0)  call mpibc2(dmatu,nsp*nlibu*(lmaxu*2+1)**2,6,mlog,'bndfp','dmatu')
         if ( lso .ne. 0 .and. lwtkb .ne. - 1 ) call mpibc2 ( orbtm_rv 
     .     , nl * nsp * nbas , 4 , mlog , 'bndfp' , 'orbtm' )
      endif
!! eterms are alwasy broadcasted. comment out if conditions jan2012
      eterms=sham%eterms ! call dcopy(size(sham%eterms),sham%eterms,1,eterms,1)
      call mpibc2(eterms(19),1,4,mlog,'bndfp','rhosig')
      sham%eterms=eterms ! call dcopy(size(sham%eterms),eterms,1,sham%eterms,1)
#ifdef MPIK      
      entime = MPI_WTIME()
      call info2(20,0,0,' MPI broadcast took %;1d sec', (entime-sttime),0)
#endif

!! ===  Write bands in bands-plotting case: loop over qp getting evals from array ===
      if (plbnd .ne. 0) then
         if (procid .eq. master) then
            print *,'procid=',procid
            call info0(20,1,0,' Writing bands to bands file (MPIK) ...')
!!dec2013 --------------------------------
c            if(fsmom/=NULLR) then
c               ifimag=ifile_handle()
c               open(ifimag,file='MagField',status='unknown')
c               read(ifimag,*) vnow
c               close(ifimag)
c            endif
!! fermi surface mode
!! together with band data around Ef (Okumura, 2019/11/14)
            if(fsmode) then
               i=nsp
               if(nspc==2) i=1
               do isp=1,i
                  ifi=ifile_handle()
                  if(isp==1) open(ifi, file='fermiup.bxsf')
                  if(isp==2) open(ifi, file='fermidn.bxsf')
                  ififm=ifile_handle()
                  if(isp==1) open(ififm, file='fermiup.data') 
                  if(isp==2) open(ififm, file='fermidn.data')
                  write(ififm,*) " # Generated by bndfp.F with --fermisurface"
                  write(ififm,*) " # q-point is reduced in 1st BZ"
                  write(ifi,*) "BEGIN_INFO"
                  write(ifi,*) " # usage xcrysden --bxsf fermiup.bxsf"
                  write(ifi,*) " # http://www.xcrysden.org/doc/XSF.html#2l.16"
                  write(ifi,'(a,f9.5)') "  Fermi Energy:",ef0
                  write(ififm,'(a,f9.5)') "  # Fermi Energy [eV]:",ef0
                  write(ififm,*) " # qshort(3) band_energy[eV]"
                  write(ifi,*) "END_INFO"
                  write(ifi,*)"BEGIN_BLOCK_BANDGRID_3D"
                  write(ifi,*)"  this_is_for_xcrysden"
                  write(ifi,*)"  BEGIN_BANDGRID_3D_simple_test"
                  nbxx=0
                  do ib=1,ndhamx
                     if(allband.or.(minval(evlall(ib,isp,:))<ef0+0.5.and.
     &                 maxval(evlall(ib,isp,:))>ef0-0.5)) then
                        nbxx=nbxx+1
                     endif  
                  enddo
                  write(ifi,"(4x,i8)") nbxx
                  write(ifi,"(4x,3i8)") nkk1,nkk2,nkk3
                  write(ifi,"(4x,3(f9.5,1x))") 0d0,0d0,0d0
                  write(ifi,"(4x,3(f9.5,1x))") qlat(:,1)
                  write(ifi,"(4x,3(f9.5,1x))") qlat(:,2)
                  write(ifi,"(4x,3(f9.5,1x))") qlat(:,3)
                  call shortn3_initialize(qlat,rlatp,xmx2)
                  do ib=1,ndhamx 
                     if(allband.or.(minval(evlall(ib,isp,:))<ef0+0.5.and.
     &                 maxval(evlall(ib,isp,:))>ef0-0.5)) then
                        write(ifi,"(a,i8)")"  BAND: ",ib    
                        write(ififm,"(a,i8)")"  # BAND: ",ib    
                        write(ifi,"(3x,10(x,f9.5))")
     &                    (evlall(ib,isp,iq), iq=1,nkk1*nkk2*nkk3)
                        ! to reduce q-point to qshort
                        do iq=1,nkk1*nkk2*nkk3
                           ppin=matmul(transpose(plat),qplist(:,iq))
                           call shortn3(rlatp,xmx2,qplist(:,iq),noutmx, nout,nlatout)
                           iout=1
                           qshort(1:3)=  matmul(qlat(:,:), ppin+nlatout(:,iout))
                           write(ififm,"(4f13.5)") qshort(:),evlall(ib,isp,iq) !q-point reduced to 1stBZ
                        enddo
                     endif  
                  enddo
                  write(ifi,*)
                  write(ifi,*)"  END_BANDGRID_3D"
                  write(ifi,*)"END_BLOCK_BANDGRID_3D"
               enddo
               close(ifi)
               close(ififm)
               call rx0('done --fermisurface mode. *.bxsf for xcryden generated')
            endif     
!! write bands directly. All inputs
            if(nsyml/=0) then
              call writeband(stdo,nkp,nsyml,nspx,nsp,ix99,ndhamx, !,fsmom/=NULLR,vnow
     i        xdatt,nqp_syml,nqp2n_syml,nevls,evlall,qplist,labeli,labele,alat,eferm,evtop,ecbot,
     i        sname,nqps_syml,nqpe_syml,nsyml,dqsyml,
     i        etolv,etolc)
            endif
         endif
         if(fsmom/=NULLR) then
            write(6,"(a)") 'NOTE: Be carful! fsmom (fixed moment on). In the sc cycle, we use additional bias mag. field  '
            write(6,"(a)") 'NOTE: Use mag.field is written in a file MagField. But it is not used for --band mode! mar2015.'
         endif
         call rx0('done')
      endif

!! ===   Repeat loop for printout and goto 99 ===
!!   jsp=isp in the collinear case; jsp=1 in the noncollinear
!!     Thus jsp should be used in place of isp
!!     isp serves as a flag for the noncollinear case
      goto99=.false.
      do  iq = 1, nkp
         qp=qplist(:,iq)
cccccccccccc
c         write(stdo,"(' procid i1 qp=',2i5,3f9.4)")procid,iq,qp(:)
cccccccccccc         
         do isp = 1, nsp
            ispc = min(isp,nspc)
            if( .not. ((onesp.eq.0 .or. isp.eq.onesp) .and. (ispc.eq.nspc)) ) cycle
            jsp = isp
            if (ispc .eq. 2) jsp = 1
            if(iprint()>20.and.ix99==1) then
               write(*,'(" bndfp: kpt",i5," of",i5, " k isp=",3f8.4,i2," ndimh nev=",2i5)')
     &           iq,nkp,qp,jsp,ndimhx_(iq),nevls(iq,jsp)
            endif   
            if(iprint()>=35.and.ix99==1) then
               write(stdo,"(9f8.4)") (evlall(i,jsp,iq), i=1,nevls(iq,jsp))
            endif
            evl(1:ndhamx,jsp) = evlall(1:ndhamx,jsp,iq)
            if (mod(iq,10) .ne. 1) call pshpr(iprint()-6)
            if(debug) call info5(30,0,0,' bndfp 299loop:  kpt %i of %i, k=%3:2,5;5d', iq,nkp,qp,0,0)
            if(debug) print *,'eeeee11111111 ',procid,numprocs
            if(numprocs==1) then
               pidorigin = master
            else
               pidorigin=-999
               do idx=0,numprocs-1
                  if (kpproc(idx)<=iq .and. iq<= kpproc(idx+1)-1) then
                     pidorigin=idx
                     exit
                  endif
               enddo   
            endif
            if (procid==pidorigin) then
               work_mpi(1) =  ndimhx_(iq)
               work_mpi(2) =  nev_(iq)
            endif
            call mpibc3(work_mpi,2,2,pidorigin,mlog,'bndfp','ndimhx_iq') !broadcast
            ndimhx_iq = work_mpi(1)
            nev_iq = work_mpi(2) !1nov2012
            if (mod(iq,10) .ne. 1) call poppr
            ebot = dmin1(ebot,evl(1,jsp))
            i = max(1,nint(qval-qbg)/(3-nspc))
            evtop = max(evtop,evl(i,jsp))
            ecbot = min(ecbot,evl(i+1,jsp))
            if (lmet .eq. 0 .and. iq .eq. 1 .and. jsp .eq. 1) ef0 = evtop
            if(debug) print *,'eeeee44444444444 plbnd=',plbnd

            if (plbnd .eq. 0) then
               if (ipr.ge.10 .and. iq.eq.1 .and. ipl.gt.0)
     .           write (stdl,712) (evl(i,jsp),i=1,nev_iq)
 712           format('fp evl',8f8.4)
               if (lwtkb .ne. -1 .and. .not. lwndow) then
                  if (iq .eq. 1 .and. jsp .eq. nsp ) then !
!! comment out mar2013takao!     .and.  .not. cmdopt('--no-fixef0',11,0,strn)) then
!! Energy-Window reset. dosw=(emin,emax)
                     ef00 = ef0
                     call fixef0(qval-qbg,jsp,1,nev_iq,ndhamx,evl,dosw,ef0)
                     if(procid==master) then
                     ifid= fopn ( 'ewindow' )
                     write(ifid,"(3d23.15)") dosw(1:2),ef0
                     call fclose(ifid)
                     endif
                     if (jsp .eq. 2 .and. ef00 .ne. ef0 .and.
     .                 lwtkb .eq. 0 .and. lmet .gt. 0 .and. lrout .ne. 0) then
                        if (procid .eq. master) call info0(10,1,1,
     .                    ' ... Fermi level reset in second spin'//
     .                    ' channel ... restart band pass')
                        goto99=.true.
                        exit    !this is the case of make co test at ecalj/TestInstall/
                     endif
                  endif
C     Check for cases when nevmx is too small : i=2 => fatal error
                  i = 0
                  if (nevmx.ge.0 .and. lmet .ne. 0) then
                     dum = evl(max(nev_iq,1),jsp)
                     if (.not. ltet .and. ef0+5*dabs(esmear-mpsord) .gt. dum) i=2
                     if (lmet.eq.4 .and. ef0+def+5*dabs(esmear-mpsord).gt.dum)i=2
                  endif
                  if(i .eq. 2) then
c                     call awrit3('%N evl(nev=%i)=%;3d but '//
c     .                 'ef0=%;3d ... restart with larger efmax or nevmx',
c     .                 ' ',80,stdo,nev_iq,evl(max(nev_iq,1),jsp),ef0)
                     write(6,"(a,f13.5,f13.5)")
     &                 'evl(nev='//trim(i2char(nev_iq))//')=',
     &                 evl(max(nev_iq,1),jsp),' but ef0=',ef0
                     call rx('bndfp:... restart with larger efmax or nevmx: bndfp')
                  endif
               endif
            endif
         enddo                  ! end second loop over isp
      enddo                     !end second loop over iq
#ifdef MPIK
      if (allocated(kpproc)) deallocate(kpproc, stat=ierr)
!! takao added !I think it must be safer to add this barrier. jan2015
      call mpi_barrier(MPI_comm_world,ierr)  
#endif
      if(goto99) goto 99


!!   For now, SO by site, not by class
      if (lwtkb .eq. 1 .and. lso .ne. 0) then
         allocate(ips_iv(nbas))
         do i_spackv=1,nbas
c            call spackv_array_copy_i8_i ( 'u' , ssite ( i_spackv ) %spec 
c     .        , 1 , i_spackv + 1 - 1 , ips_iv )
            ips_iv(i_spackv)=ssite ( i_spackv ) %spec 
         enddo
!! orbital moment
         call iorbtm ( sspec , ips_iv , nl , nl , nbas , nsp , orbtm_rv  )
         if (allocated(ips_iv)) deallocate(ips_iv)
      endif
ccC ... Case generating bands: find next block of qp
cc      if (plbnd .ne. 0) goto 99
      if (ipl.gt.0) write (stdl,715) nkp,ebot,qval,qbg,esmear
 715  format('nv nkp',i5,'  ebot',f9.4,'   qval',f10.4,'  qbg',f8.4,'  smr',f8.4)


!! --- Interpolate density to Fermi energy ---
      sev = sumev(1,1)
      if (lmet .eq. 4) then
         call mshn3p ( nbas , ssite , sspec , lmet , lrout , lfrce , qval 
     .     - qbg , ef0 , def , sumqv , sumev , n1 , n2 , n3 , k1 , k2 , 
     .     k3 , srout_zv , sv_p_oqkkl , frc , lrep )
C   ... Store val q & magnetic moment in sumqv(1) and sumqv(2)
         sumqv(2,1) = sumqv(1,1) - sumqv(1,2)
         sumqv(1,1) = sumqv(1,1) + sumqv(1,2)
C   ... Eigenvalue sum including entropy term
         sev = sumev(1,1) + sumev(2,1)
C   ... Remake sev,ef linearly interpolating tabulated sampling DOS
         sev00 = sev
         ef00  = ef0
         if (ldos .ne. 0) then
            call efldos ( qval , nsp , emin , emax , ndos , dos_rv , eferm, sev1 )
            sev   = sev1
            ef0  = eferm
         endif
         if (ipr .gt. 30 .and. ldos .ne. 0)
     .     write(stdo,388) sev00,ef00,sev1,eferm,sev,ef0
 388     format(' ipol:  sev=',f12.6,'   ef=',f12.6:/
     .     ' dos:   sev=',f12.6,'   ef=',f12.6/
     .     ' use:   sev=',f12.6,'   ef=',f12.6)
         if (ipl.gt.0) write (stdl,733) ef00,eferm,ef0,sev00,sev1,sev
 733     format('nf EF:',3f9.5,'    EB:',3f12.5)
         sbz%ndos=ndos
         call dcopy(size(sbz%dosw),dosw,1,sbz%dosw,1)
         sbz%ef=ef0
         sbz%def=def
         if (lrep .eq. 1) then
            ef0 = -1
c            call awrit0('Input Fermi energy was too far off, '//
c     .        'repeat band pass',' ',80,stdo)
           write(6,*)'Input Fermi energy was too far off, repeat band pass'
#ifdef MPIK
            call mpi_barrier(MPI_comm_world,ierr)  
#endif
            goto 99
         endif
      endif

!! delete Magfield file for safe if it exist. 
#if MPIK
      call mpi_barrier(MPI_comm_world,ierr)
#endif
      if (procid .eq. master) then
         ifimag=ifile_handle()
         open(ifimag,file='MagField',status='unknown')
         close(ifimag,status='delete')
      endif  

!! --- BZ integration for fermi level, band sum and qp weights ---
      if(lmet .ge. 0 .and. (lmet .ne. 4 .or. ltet)) then
c$$$c jan2015 this is not allowed to go through the TestIntall/te
c$$$         if(ix99==2) then
c$$$           if(procid==master) write(*,*) 
c$$$           if(procid==master) write(*,*) '... Goto bzints again (output suppressed) ...'
c$$$           call pshpr(1)
c$$$c --> and call poppr before goto99
c$$$         endif
         if (lwndow) then
            allocate(ww_rv(ndham*nsp*nkp))
            eferm = min(dosw(1),dosw(2))
            call bzints ( nkabc ( 1 ) , nkabc ( 2 ) , nkabc ( 3 ) , evlall 
     .        , ww_rv , nkp , ndham , ndham , nsp , xv , xv , xv , 1 , eferm 
     .        , 2 , ntet , sbz%iv_a_oidtet , sev , dum )
            eferm = max(dosw(1),dosw(2))
            call bzints ( nkabc ( 1 ) , nkabc ( 2 ) , nkabc ( 3 ) , evlall 
     .        , rv_a_owtkb , nkp , ndham , ndham , nsp , xv , xv , xv , 1 , 
     .        eferm , 2 , ntet , sbz%iv_a_oidtet , sev , dum )
            call daxpy ( ndham * nsp * nkp , - 1d0 , ww_rv , 1 , rv_a_owtkb, 1 )
            if (allocated(ww_rv)) deallocate(ww_rv)
         else
            dosrng = 8
            if (mpsord .lt. 0) dosrng = 16
            if( debug) print *, 'tttttttttt 99999999 6 call bzwtsf'
            if(sbz%fsmommethod == 1) then !takao dec2010
!! vnow june22013 !vnow !(in Ry) contains magnetic field
!! For eigenvalus, add  -vnow/2 for isp=1, and +vnow/2 for isp=2.
               call bzwtsf2 ( ndham , ndham , nsp , nspc , nkabc ( 1 ) , nkabc 
     .           ( 2 ) , nkabc ( 3 ) , nkp , ntet , sbz%iv_a_oidtet , qval - qbg , 
     .           fsmom , lmet.ne.0 , ltet , mpsord , ndos , dabs ( esmear - mpsord 
     .           ) , dosrng , sbz%rv_a_owtkp , evlall , efmax , lswtk , rv_a_oswtk 
     .           , eferm , sev , rv_a_owtkb , sumqv ( 1 , 2 ) , lwtkb ,lfill,vnow)
            else
ckino       maybe a switch for these is lmet.ne.0
               l_temp_rv_a_oswtk=.false.
               if (.not.allocated(rv_a_oswtk)) then
                  allocate(rv_a_oswtk(1))
                  l_temp_rv_a_oswtk=.true.
               endif
               l_temp_rv_a_owtkb=.false.
               if (.not.allocated(rv_a_owtkb))then
                  allocate(rv_a_owtkb(1,1,1))
                  l_temp_rv_a_owtkb=.true.
               endif
               l_temp_iv_a_oidtet=.false.
               if (.not.allocated(sbz%iv_a_oidtet))then
                  allocate(sbz%iv_a_oidtet(1))
                  l_temp_iv_a_oidtet=.true.
               endif
!! vnow june22013
               call bzwtsf ( ndham , ndham , nsp , nspc , nkabc ( 1 ) , nkabc 
     .           ( 2 ) , nkabc ( 3 ) , nkp , ntet , sbz%iv_a_oidtet , qval - qbg , 
     .           fsmom , lmet.ne.0 , ltet , mpsord , ndos , dabs ( esmear - mpsord 
     .           ) , dosrng , sbz%rv_a_owtkp , evlall , efmax , lswtk , rv_a_oswtk 
     .           , eferm , sev , rv_a_owtkb , sumqv ( 1 , 2 ) , lwtkb ,lfill, vnow)
            endif
!! june2013 magfield is added
            if(fsmom/=NULLR.and.procid .eq. master) then
               ifimag=ifile_handle()
               open(ifimag,file='MagField',status='unknown')
               write(ifimag,"(d23.16,' !(in Ry) -vnow/2 for isp=1, +vnow/2 for isp=2')")vnow
               close(ifimag)
            endif

            if (l_temp_rv_a_oswtk) then
               deallocate(rv_a_oswtk)
               l_temp_rv_a_oswtk=.false.
            endif
            if (l_temp_rv_a_owtkb) then
               deallocate(rv_a_owtkb)
               l_temp_rv_a_owtkb=.false.
            endif
            if (l_temp_iv_a_oidtet) then
               deallocate(sbz%iv_a_oidtet);
               l_temp_iv_a_oidtet=.false.
            endif
!!         Store val charge & magnetic moment in sumqv(1..2)
            if (lmet .ne. 4) then
               sumqv(1,1) = sumqv(1,2)
               sumqv(2,1) = sumqv(2,2)
            endif
         endif
!!
         if (lmet .ne. 4) then
            ef0 = eferm
            sbz%ef=ef0
         endif
         if (lmet .gt. 0) then
            if (procid .eq. master) then
               ifi = fopna('wkp',-1,4)
               i = iobzwt ( 0 , ndhamx , nkp , nspx , eferm , rv_a_owtkb , -ifi )
               call fclr('wkp',ifi)
            endif
         endif
         goto99=.false.
         if (lwtkb .eq. -1 .and. lrout .gt. 0) then
            call info0(20,0,0,' Start second band pass ...')
            lwtkb = 1
c            if (nspc .eq. 2) lswtk = 1 !commnet out nov2015  right?
            goto99=.true.
         endif
         if (lwtkb .eq. 2 .and. lrout .gt. 0) then
            call info0(20,0,0,' New pass with constrained weights ...')
            goto99=.true.
         endif
#ifdef MPIK
         call mpi_barrier(MPI_comm_world,ierr)  
#endif
c$$$c jan2015 did not work for Testinstall/te
c$$$         call poppr
         if(goto99) goto 99
      endif
      if(allocated(sfz)) deallocate(sfz)
!! == 99 loop over =============================


!! ... Save Fermi level, nonmetal or sampling integration
      if (lmet .eq. 0 .or. .not. (lmet .ne. 4 .or. ltet)) then
         if (lmet .eq. 0) then
            ef0 = (evtop+ecbot)/2
c            if (ipr .ge. 20) call awrit3(' Highest occ. level = %,5;5d '//
c     .        ' Lowest unocc. = %,5;5d  average = %,5;5d',' ',80,
c     .        stdo,evtop,ecbot,ef0)
            write(6,"(' Highest occ. level = ',f13.5,
     &        ' Lowest unocc. = ',f13.5, ' average =',f13.5)") evtop,ecbot,ef0
         endif
         if (procid .eq. master) then
            ifi = fopna('wkp',-1,4)
            i = iobzwt ( 1 , ndham , nkp , nsp , ef0 , rv_a_owtkb , -ifi )
            call fclr('wkp',ifi)
         endif
      endif

!! ... Cleanup asa-style moments file, print table of DOS channels
!! now only for --mull mullikan analysis mode.
      if (mullwrite) then
         i = iomoms(-nfilem,nl,nsp,nspc,nkp,ndimh,11,1,nkp*nsp+1,1,
     .     ndham,ndimh,nchan,nchan,ndimh,0d0,0d0,0d0,0d0,eferm,0d0)
         call fclose(nfilem)
         call mchan ( lmdim , ssite , sspec , nsp , nsites , lsites , 
     .     0 , 0 , 0 , 0 , chan_iv )
      endif

!! ... Average forces so net force on system is zero (APW case)
      if (lfrce .ne. 0 .and. napw .ne. 0) then
         call dpzero(xv,3)
         do  i1 = 1, nbas
            do  i = 1, 3
               xv(i) = xv(i) + frc(i,i1)/nbas
            enddo
         enddo
         do  i1 = 1, nbas
            do  i = 1, 3
               frc(i,i1) = frc(i,i1) - xv(i)
            enddo
         enddo
      endif


!! === Generate DOS on disk ===
      if(.not. lwndow .and. procid==master) then
         ndos0 = sbz%ndos !ndos0 is only used as a swith wether DOS or integrated-DOS.
         if ( iand(1,int(sctrl%ldos)) .ne.0 ) then
c  if(iprint()>10) call awrit1('%x%N ... Generating %?#n<0#integrated#total# DOS',' ',80,lgunit(1),ndos0)
            if(ndos0<0) write(6,"(a)")'... Generating integrated DOS' !, ndos0=,ndos0
            if(ndos0>0) write(6,"(a)")'... Generating total DOS'      !  ndos0=,ndos0

!! Default and external option. 
            if(tdos) then
              dosmax = sbz%dosmax
c            if (cmdopt('--emint=',7,strn)) then    !not yer for dosw(1)
c              j=10
c              if (a2bin(strn,eminp_,4,0,' ',j,len(strn))) eminp=eminp_
c            endif
              if (cmdopt('-emax=',6,0,strn)) then
                j = 6
                if (a2bin(strn,emaxp_,4,0,' ',j,len(strn))) dosmax=emaxp_/rydberg
              endif
              if(cmdopt('-ndos=',6,0,strn))  then
                j = 6
                if (a2bin(strn,ndos_,2,0,' ',j,len(strn))) ndos=ndos_
                deallocate(dos_rv)
                allocate(dos_rv(ndos,2*nsp))
              endif
            endif
            if(cmdopt('--tdostetf',10,0,strn)) ltet=.false. !oct2015

            if (ltet) then
!! determine eferm
               call bzints ( nkabc ( 1 ) , nkabc ( 2 ) , nkabc ( 3 ) , evlall 
     .           , dum , nkp , ndhamx , ndhamx , nspx , dosw ( 1 ) , dosw (
     .           2 ) , dos_rv , ndos , eferm , 1 , ntet , sbz%iv_a_oidtet , dum , dum )
!! repeat again for given range emax=dosw(2). mar2012
               if(tdos) dosw(2) = dosmax + eferm !.or.pdos
               call bzints ( nkabc ( 1 ) , nkabc ( 2 ) , nkabc ( 3 ) , evlall 
     .           , dum , nkp , ndhamx , ndhamx , nspx , dosw ( 1 ) , dosw (
     .           2 ) , dos_rv , ndos , eferm , 1 , ntet , sbz%iv_a_oidtet , dum , dum )!job=1 integrated dos
               if ( ndos0 > 0 ) call xxxdif ( dosw ( 1 ) , dosw ( 2 ) , ndos 
     .           , nspx , 0 , dos_rv ) !derivative 'integrated dos -> DOS'
               del = 0d0
            else
               if(tdos) dosw(2) = dosmax + eferm  !.or.pdos
               if ( mpsord .ge. 100) mpsord = mod(mpsord,100)
               if ( ndos0 > 0 ) call makdos ( nkp , ndhamx , ndhamx , nspx
     .           , sbz%rv_a_owtkp , evlall , mpsord , dabs ( esmear - mpsord ) 
     .           , - 6d0 , dosw ( 1 ) , dosw ( 2 ) , ndos , dos_rv )
               if ( ndos0 < 0 ) call maknos ( nkp , ndhamx , ndhamx , nspx
     .           , sbz%rv_a_owtkp , evlall , mpsord , dabs ( esmear - mpsord ) 
     .           , - 6d0 , dosw ( 1 ) , dosw ( 2 ) , ndos , dos_rv )
               del = mpsord+dabs(esmear-mpsord)
            endif
            if ( nspc .eq. 2 ) call dscal ( ndos , .5d0 , dos_rv , 1 )
!! ... dos file io
            call dosio ( dos_rv , ndos , nspx , ndos , 1 , dosw(1) 
     .        , dosw(2) , nspx , eferm , del , 1 , - fopn ( 'dos' ) )
            call fclose(fopn('dos'))
!! ... write total dos file. mar2013
c            ifi=fopn('dos.tot')
            ifi= ifile_handle()     
            open(unit=ifi,file='dos.tot.'//trim(sname))
c            write(ifi,"('#lm ordering. See the end of lmdos. relative to efermi=',d14.5)") eferm
            do ipts=1,ndos
               eee=dosw(1)+ (ipts-1d0)*(dosw(2)-dosw(1))/(ndos-1d0)-eferm
               write(ifi,"(255(f13.5,x))")eee,(dos_rv(ipts,isp),isp=1,nsp)
            enddo
c            call fclose(ifi)
            close(ifi)
            if (allocated(orbtm_rv)) deallocate(orbtm_rv)
            if (allocated(dos_rv)) deallocate(dos_rv)
         endif
      endif
      if(tdos) call rx0('Done tdos mode:')

!! efermi.lmf 
      if(procid==master.and.plbnd==0) then
         ifi= ifile_handle()     
         open(unit=ifi,file='efermi.lmf')
         write(ifi,"(d16.8, ' # (Ry) Fermi energy given by lmf')") eferm
         write(ifi,"(d16.8, ' # (Ry) Top of Valence')") evtop
         write(ifi,"(d16.8, ' # (Ry) Bottom of conduction')") ecbot
         write(ifi,"(d16.8, ' # number of electrons qval-qbg=')") qval-qbg
         write(ifi,"(d16.8, ' # mag. mom (nup-ndown)=')") sumqv(2,1) 
         write(ifi,
     &     "('# NOTE: These valuse are overwritten by lmf(-MPIK) in SC loop')")
         write(ifi,"('# NOTE: Not by band plot mode')")
         close(ifi)             
      endif

!! --- Core-level spectroscopy ---
      if (icls .ne. 0) then
         if (procid .eq. master) then
            eferm=sbz%ef
            call vcdmel ( sctrl%nl , ssite , sspec , slat , nlmax , ndham , 
     .        ndimh , nkp , nsp , nspc , eferm , evlall , ausc_zv , 
     .        nsites , isite , iclsl , iclsn )
            call rx0('done generating core level spectra')
         else
            call fexit(0,0,' ',0)
         endif
      endif
!! --- Assemble output density ---
      call dfratm ( nbas , ssite , sspec , sv_p_orhat1 )
      allocate(qbyl_rv(n0,nsp,nbas))
      allocate(hbyl_rv(n0,nsp,nbas))
      if(debug) write(*,*) 'kino: allocate hbyl_rv,',n0*nbas*nsp

!!     --window: Put output density into orhat and smrho, and exit
      if (lwndow) then
         call mkrout ( ssite , sspec , slat , sham , nbas , nsp , ldim
     .     , lekkl , sv_p_oqkkl , sv_p_oeqkkl , spot%sv_a_oorhat , hab_rv , sab_rv 
     .     , qbyl_rv , hbyl_rv , lrout )
         call zcopy( k1 * k2 * k3 * nsp , srout_zv , 1 , spot%zv_a_osmrho 
     .     , 1 )
         call symrho(ssite , sspec , slat , lfrce , spot%zv_a_osmrho , spot%sv_a_oorhat 
     .     , qbyl_rv , hbyl_rv , frc )
         goto 999
      endif

      call mkrout ( ssite , sspec , slat , sham , nbas , nsp , ldim 
     .  , lekkl , sv_p_oqkkl , sv_p_oeqkkl , sv_p_orhat1 , hab_rv , sab_rv 
     .  , qbyl_rv , hbyl_rv , lrout )
C
      if (lrout .ne. 0) then
C   ... Symmetrize output density and forces
         call symrho ( ssite , sspec , slat , lfrce , srout_zv , sv_p_orhat1 
     .     , qbyl_rv , hbyl_rv , frc )
C   ... New boundary conditions pnu for phi and phidot
C       call pshpr(iprint()-12)
         if (lpnu .gt. 0) then
            call dcopy(size(sham%pmin),sham%pmin,1,xv,1)
            call dcopy(size(sham%pmax),sham%pmax,1,xv(11),1)
            call pnunew ( nbas , nsp , ssite , sspec , xv , xv ( 11 ) , lfrzw 
     .        , hab_rv , sab_rv , qbyl_rv , hbyl_rv, eferm )
         endif
C       call poppr
      endif

!! --- print out qbyl ---
      if(procid==master) then
         ifqbyl = ifile_handle() !fopna('qbyl',-1,2)
         open(ifqbyl,file='qbyl.'//trim(sname))
         write(ifqbyl,"(i5)")nbas
         do  ib = 1, nbas
           ispec=ssite(ib)%spec
           write(ifqbyl,"(i4, i3,i2,10f12.6)")ib, ispec, sspec(ispec)%lmxa, (sum(qbyl_rv(il,1:nsp,ib)),il=1,sspec(ispec)%lmxa+1)
         enddo  
         close(ifqbyl)
      endif

!! --- Evaluate Harris energy ---
      call mkehkf(1,sham,sev,sumqv(2,1),sumtv,ehar)
!! --- Evaluate KS total energy, correction to Harris force ---
      if (lrout .ne. 0) then
!!   ... Correction to harris force
        job=int(sctrl%lfrce)
        if(job>0 ) then
          call dfrce (job, ssite , sspec , slat , k1 , k2 , k3 , nvl !sctrl , 
     .     , spot%sv_a_oorhat , sv_p_orhat1 , elind , qmom_rv , spot%zv_a_osmrho , 
     .     srout_zv , fh_rv )
        endif
!!   ... Evaluate KS total energy and output magnetic moment
         eks = 0d0
         if (leks .ge. 1) then
            call togpr()
            call mkekin ( nbas , ldim , ssite , sspec , slat , sham , sv_p_osig 
     .        , sv_p_otau , sv_p_oppi , lcplxp , sv_p_oqkkl , k1 , k2 , k3
     .        , vconst , spot%zv_a_osmpot , srout_zv , sev , sumtv )

ccccccccccccccccccccccccc
ctakao test mkekin repeat
c         print *,'vvv takao mkekin repeat 000 --- oldpot x oldrho'
c         call mkekin(nbas,ldim,ssite,sspec,slat,sham,osig,otau,oppi,
c     .      lcplxp,oqkkl,k1,k2,k3,vconst,w(osmpot),w(osmrho),sevxxx,sumtvxxx)
c       stop 'xxxxxxxxxxxxxxxx  111aaa'
ccccccccccccccccccccccccc
ccccccccccccccccccccc
c          print *,'vvv: takao test osmpotbk 000 '
c          call defdc(osmpotbk,k1*k2*k3*2)
c          call dcopy(k1*k2*k3*4,w(osmpot),1,w(osmpotbk),1)
ccccccccccccccccccccc
c          call pshpr(ipr-4)

!! mkpot calls locpot. and locpot calls augmat.
!! augmat calculates sig,tau,ppi.
c          print *,'goto mkpot procid=',procid
            call mkpot ( nbas , ssite , sspec , slat , sham , spot , lfrce
     .        , lcplxp, k1, k2, k3, srout_zv , sv_p_orhat1 , qbg , spot%zv_a_osmpot 
     .        , qmom_rv , vconst , vrmt , sv_p_osig , sv_p_otau , sv_p_oppi 
     .        , ppnl_rv , hab_rv , vab_rv , sab_rv , qval , qsc , gpot0_rv 
     .        , vval_rv , fes2_rv , 0 , vorb , nlibu , lmaxu , lldau , nsp,0 )

           if(cmdopt('--density',9,0,strn)) then
#ifdef MPIK
             call mpi_barrier(MPI_comm_world,ierr)  
#endif
             call rx0('end of --density mode')
           endif  
             
ccccccccccccccccccccccccccccccccccc
ctakao test mkekin repeat
c          print *,'vvv takao mkekin repeat 001x ---newpot x newrho'
c          call mkekin(nbas,ldim,ssite,sspec,slat,sham,osig,otau,oppi,
c     .      lcplxp,oqkkl,k1,k2,k3,vconst,w(osmpot),srout_zv,sev,sumtvxxx)
c         print *,'vvv takao mkekin repeat 001 --- newpot x oldrho'
c         call mkekin(nbas,ldim,ssite,sspec,slat,sham,osig,otau,oppi,
c     .      lcplxp,oqkkl,k1,k2,k3,vconst,w(osmpot),w(osmrho),sev,sumtvxxx)
c          print *,'vvv takao mkekin repeat 002 --- oldpot x oldrho'
c          call mkekin(nbas,ldim,ssite,sspec,slat,sham,osig,otau,oppi,
c     .      lcplxp,oqkkl,k1,k2,k3,vconst,w(osmpotbk),w(osmrho),sev,sumtvxxx)
c          call rlse(osmpotbk)
c          stop 'xxxxxxxxxxxxxxxx'
cccccccccccccccccccccccccccccccc
            call mkehkf(2,sham,sev,sumqv(2,1),sumtv,eks)
            call togpr()
         endif
         if(debug) print *, 'tttttttttt 9999999999999999999 9'

!! --- Add together force terms ---
!! fes1_rv: contribution to HF forces from estat + xc potential   This is for input  density  !=3rd term in (B.5) in JPSJ.84.034705
!! fes2_rv: contribution to KS forces from estat + xc potential   This is for output density
!! fh_rv      : 1st term  in (B.5)
!! input frc  : 2nd term in (B.5). Puley term
!! output frc  : totlal
         if(lfrce .gt. 0) call totfrc( nbas , ssite , slat , leks 
     .     , fes1_rv , fes2_rv , fh_rv , frc )


!!   --- Mix input and output densities ---
         allocate(qm1_rv(nvl))
         allocate(qm2_rv(nvl))

cccccccccccccccccccccc
cc          print *,'vvv takao mkekin copy test xxx1---'
cc          allocate(smrhobk(k1*k2*k3*2))
cc1          call dcopy(k1*k2*k3*4,w(osmrho),1,smrhobk,1)
cc2          call dcopy(k1*k2*k3*2,srout_zv,1,smrhobk,1)
ccccccccccccccccccccccc

c         call upacks('strn mix',i1,i2)
         call mixrho ( ssite , sspec , slat , nsp , iter , sstrnmix
c         call mixrho ( ssite , sspec , slat , nsp , iter , sstrn ( i1:i2)
     .     , qval - qbg , elind , sv_p_orhat1 , spot%sv_a_oorhat , k1 , k2
     .     , k3 , dmxp , srout_zv , spot%zv_a_osmrho )

ccccccccccccccccccccccccc
ctakao test mkekin repeat
cc          print *,'vvv takao mkekin copy test xxx2---'
cc          call dcopy(k1*k2*k3*2,smrhobk,1,w(osmrho),1)
cc          deallocate(smrhobk)
c          print *,'vvv takao mkekin repeat 2 ---'
c          call dcopy(k1*k2*k3*4,w(osmpotbk),1,w(osmpot),1)
c          call mkekin(nbas,ldim,ssite,sspec,slat,sham,osig,otau,oppi,
c     .      lcplxp,oqkkl,k1,k2,k3,vconst,w(osmpot),w(osmrho),sev,sumtvxxx)
c          stop 'xxxxxxxxxxxxxxxx'
ccccccccccccccccccccccccc

      else
         eks = 0d0
      endif
      sham%ehf=ehar
      sham%ehk=eks
 999  continue

!! === VBM CBM section. Mar2013takao === just output
 8888 continue
      if (procid .eq. master) then
         ifvesintatm=fopna('vesintatm',-1,0) !generated by lmfa. electrostatic potential integrals.
         rewind ifvesintatm
         allocate(vesatm(sctrl%nspec))
         do ispec=1,sctrl%nspec
            if(sspec(ispec)%z==0d0 .and. sspec(ispec)%rmt==0d0 ) cycle
            read(ifvesintatm,*,end=9898) vesatm(ispec)
         enddo
         ifves=fopna('vessm',-1,1)
         rewind ifves
         read(ifves,*) vessm !smooth part of electrostatic pot. given in smves.F. (es on MT average).
c   ifvesintloc=fopna('vesintloc',-1,1) !generated in mkpot-locpot.F electrostatic pot. integral.
         ifvesintloc = ifile_handle()
         open(ifvesintloc,file='vesintloc',form='formatted',status='old')
c         rewind ifvesintloc
         sumvesloc=0d0
         sumvesatm=0d0
         do ib=1,nbas
            ispec=ssite(ib)%spec
            if(sspec(ispec)%z==0d0 .and. sspec(ispec)%rmt==0d0 )  cycle
            read(ifvesintloc,*) vesloc
            sumvesloc = sumvesloc + vesloc
            sumvesatm = sumvesatm + vesatm(ssite(ib)%spec)
         enddo   
         deallocate(vesatm)
         print *,'vessm=',vessm
         vref= sumvesatm/vol - vessm -sumvesloc/vol
         write(*,"('### Add vref to eigenvalues to get relative to vaccum=',f12.6)") vref
         write(*,"('### mean estatic pot by superposition of atoms(eV)=',f12.6)")sumvesatm/vol*rydberg
         if(cmdopt('--vbmonly',9,0,strn)) then
            ifvbm = fopna('vbmcbm',-1,1)
            read(ifvbm,"(a)") vbmlll
            print *,'readin vbmbm--> ',vbmlll
            read(vbmlll(9:),*)evbm
            read(vbmlll(30:),*)ecbm
            read(vbmlll(59:),*)esold
            write(*,*)evbm,ecbm,esold
            evbm= evbm +sumvesatm/vol*rydberg-esold
            ecbm= ecbm +sumvesatm/vol*rydberg-esold
         elseif(lfill) then
            evbm=(evtop+vref)*rydberg
            ecbm=(ecbot+vref)*rydberg
         else
            evbm=(eferm+vref)*rydberg
            ecbm=(eferm+vref)*rydberg
         endif   
         ifvbm = fopna('vbmcbm',-1,2)
         write(*,    "('### VBM= ',f12.6,' eV CBM= ',f12.6,' eV estatic_pot=',f12.6,' eV')")  
     &     evbm, ecbm, sumvesatm/vol*rydberg
         write(ifvbm,"('### VBM= ',f12.6,' eV CBM= ',f12.6,' eV estatic_pot=',f12.6,' eV')")  
     &     evbm, ecbm, sumvesatm/vol*rydberg
         call fclose(ifvbm)
         call fclose(ifves)
         call fclose(ifvesintatm)
c        call fclose(ifvesintloc)
         close(ifvesintloc)
 9898    continue
         print *,' not readin vesintatm.* ---> not show VBM. or you need to repeat new lmfa again.'
         if(cmdopt('--vbmonly',9,0,strn)) call rx0('end of --vbmonly mode vrefolny')
      endif
!!
      if (allocated(qm2_rv)) deallocate(qm2_rv)
      if (allocated(qm1_rv)) deallocate(qm1_rv)
      if (allocated(hbyl_rv)) deallocate(hbyl_rv)
      if (allocated(qbyl_rv)) deallocate(qbyl_rv)
      if (allocated(fes1_rv)) deallocate(fes1_rv)
      if (allocated(ppnl_rv)) deallocate(ppnl_rv)
      if (allocated(sab_rv)) deallocate(sab_rv)
      if (allocated(vab_rv)) deallocate(vab_rv)
      if (allocated(hab_rv)) deallocate(hab_rv)
      if (allocated(vval_rv)) deallocate(vval_rv)
      if (allocated(gpot0_rv)) deallocate(gpot0_rv)
      if (allocated(qmom_rv)) deallocate(qmom_rv)
      call tcx('bndfp')
c      if ( ifiese>0 ) close(ifiese)
      if (allocated(ifbls_iv)) deallocate(ifbls_iv)
      deallocate(evlall)
      return
 9995 continue
c      call rx('bndfp: readin error on sigm or sigm_fbz(--pdos) file')
      end subroutine bndfp


#if MPIK
      subroutine xmpbnd(kpproc,ndham,nkp,nsp,eb)
C- Collect eb from various processors (MPI)
C ----------------------------------------------------------------------
Ci Inputs
Ci   kpproc
Ci   ndham :leading dimension of eb
Ci   nkp   :number of irreducible k-points (bzmesh.f)
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   eb    :energy bands; alias eband
Co Outputs
C ----------------------------------------------------------------------
      implicit none
      include "mpif.h"
      integer kpproc(0:*),ndham,nkp,nsp
      double precision eb(ndham,*)
      integer i,ista,iend
      integer procid,master
      integer numprocs, ierr
      integer, dimension(:),allocatable :: offset,length
      real(8) ,allocatable :: buf_rv(:)
      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      allocate (offset(0:numprocs), stat=ierr)
      allocate (length(0:numprocs), stat=ierr)
      offset(0) = 0
      do  i = 0, numprocs-1
        ista = kpproc(i)
        iend = kpproc(i+1)-1
        length(i) = (iend - ista + 1)*nsp*ndham
        offset(i+1) = offset(i) + length(i)
      enddo
      ista = kpproc(procid)
      allocate(buf_rv(ndham*nkp*nsp))
#ifdef KINODEBUG
      write(*,*)'kino allocate bug_rv',ndham*nkp*nsp
#endif
      call mpi_allgatherv ( eb ( 1 , 1 + nsp * ( ista - 1 ) ) , length 
     .( procid ) , mpi_double_precision , buf_rv , length , offset 
     ., mpi_double_precision , mpi_comm_world , ierr )
      call dcopy ( ndham * nsp * nkp , buf_rv , 1 , eb , 1 )
      if (allocated(buf_rv)) deallocate(buf_rv)
      deallocate(offset, stat=ierr)
      deallocate(length, stat=ierr)
      end subroutine xmpbnd
#endif

c      subroutine xxxbfp
c
c      end subroutine xxxbfp

c      subroutine getnevec(ndham,nsp,nkp,nev,wtkb, nevec)
c      implicit none
c      integer:: ndham,nsp,nkp,nevec,i,nev
c      real(8):: wtkb(ndham,nsp,nkp)
c      print *,'nev sum(wtkb)=',nev,sum(wtkb)
c      do  i = nev, 1, -1 !this is the same as nvec search in addrbl.
c        nevec = i
c        if (maxval(wtkb(i,:,:)) .gt. 1d-6) exit
c      enddo
c      end

      logical function master_mpi()
      implicit none
      logical,save:: mm=.true.
      logical::  master_mpiset
      integer:: procid
      master_mpi=mm
      return
      entry master_mpiset(procid)
      mm=.false.
      if(procid==0) mm=.true.
      end

!! ----------------------------------
      subroutine writepdos(ext)
!! == readin pdosinput, and print out pdos files. ==
      implicit none
      integer:: ifip,ndhamx,nsp,nspx,nevmin,nchanp,nbas,nkk1,nkk2,nkk3,ntete,ndos,nkp
     & ,ibas,jsp,ifi,init,iend,ipts,fopna,fopnn,j,ndos_,ichan,isp,itet,ksp,i,ib
      integer,allocatable::idtete(:,:),ipqe(:,:,:)
      real(8),allocatable:: evlall(:,:,:),dwgtall(:,:,:,:,:),pdosp(:,:),pdosalla(:,:,:,:)
      real(8)::eminp,emaxp,ef0,eee,eminp_,emaxp_
      character*3::charnum3
c      character*100::filenm(2)
      real(8):: rydberg=13.6058d0, bin,eigen(4),vvv,wt,bin2
      character strn*120
      character*(*)::ext
      logical::a2bin,cmdopt,mlog

      integer, dimension(:),allocatable :: kpproc
      integer::numprocs,procid,ierr,itete,iteti,mpipid,ifile_handle
#if MPIK
      include "mpif.h"
#endif

cfopna('Pdosdata',-1,4)
c      ifip= ifile_handle()
      print *,' pdosdata file=','pdosdata.'//trim(ext)
      open(newunit=ifip,form='unformatted',file='pdosdata.'//trim(ext))
      read(ifip) ndhamx,nsp,nspx,nevmin,nchanp,nbas,nkk1,nkk2,nkk3,ntete,nkp !ndos,nkp mar2015
      allocate(idtete(0:4,6*nkp),ipqe(nkk1,nkk2,nkk3))
      allocate(evlall(ndhamx,nspx,nkp))             
      allocate(dwgtall(nchanp,nbas,ndhamx,nspx,nkp)) 
      read(ifip) idtete !,ipqe 
      read(ifip) evlall         
      read(ifip) dwgtall        
      read(ifip) ef0 !may2015
      close(ifip)
c      call fclose(ifip)

#if MPIK
      mlog = cmdopt('--mlog',6,0,strn) !--mlog here is taken by getarg.
      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      allocate (kpproc(0:numprocs), stat=ierr)
      call dstrbp(ntete,numprocs,1,kpproc(0))
      iteti = kpproc(procid)
      itete = kpproc(procid+1)-1
c      print *,'ppppp',numprocs,procid,iteti,itete
#else
      iteti = 1
      itete = ntete
#endif

!! Default and external option.
      eminp=-25.0/rydberg
      emaxp=30.0/rydberg
      ndos=5500
      if (cmdopt('-emin=',6,0,strn)) then
        j=6
        if (a2bin(strn,eminp_,4,0,' ',j,len(strn))) eminp=eminp_/rydberg
      endif
      if (cmdopt('-emax=',6,0,strn)) then
        j =6
        if (a2bin(strn,emaxp_,4,0,' ',j,len(strn))) emaxp=emaxp_/rydberg
      endif
      if(cmdopt('-ndos=',6,0,strn))  then
        j = 6
        if (a2bin(strn,ndos_,2,0,' ',j,len(strn))) ndos=ndos_
      endif
      
!! dostet ---
c      call dostet(ndhamx,nsp,nspx,nevmin,nchanp*nbas,nkk1,nkk2,nkk3,ntete,idtete,evlall,
c     &  dwgtall, ndos, eminp+ef0, emaxp+ef0,.false.,wkd,pdosall)
      bin = (emaxp - eminp) / (ndos - 1)
      vvv = ( 3d0  -  nsp ) / ( nkk1 * nkk2 * nkk3 * 6d0 )/4d0
      allocate(pdosalla(ndos,nsp,nchanp,nbas))
!! --- Loop over tetrahedra ---
      do itet = iteti, itete
      do isp = 1, nspx
      do ib = 1, nevmin
        eigen(1:4) = evlall(ib,isp,idtete(1:4,itet))
        if( minval(eigen) > emaxp+ef0 ) cycle
        do ibas = 1,nbas
        do ichan = 1, nchanp
          wt = sum(dwgtall(ichan,ibas,ib,isp,idtete(1:4,itet))) * idtete(0,itet) * vvv
          call slinz(wt,eigen,eminp+ef0,emaxp+ef0,pdosalla(1,isp,ichan,ibas),ndos)
        enddo
        enddo
      enddo
      enddo
      enddo
#if MPIK
      call mpibc2( pdosalla, ndos*nsp*nchanp*nbas, 4 , mlog , 'writepdos' , 'pdosalla' )
#endif
      if(procid==0) then !master only
      allocate(pdosp (ndos,nchanp))
      bin2 = 2d0 * bin
      do isp =1,nspx
      do ibas=1,nbas
        ifi = fopnn(trim(  'dos.isp'//char(48+isp)//'.site'//charnum3(ibas)  ))
        write(ifi,"('#lm ordering. See the end of lmdos. relative to efermi')")
!!      DOS from finite difference of NOS 
        do ichan=1,nchanp
          do i = 2, ndos - 1
            pdosp(i,ichan)=(pdosalla(i+1,isp,ichan,ibas) - pdosalla(i-1,isp,ichan,ibas))/bin2
          enddo
          pdosp(1,ichan)    = pdosp(2,ichan)
          pdosp(ndos,ichan) = pdosp(ndos-1,ichan)
        enddo
        do ipts=1,ndos
          eee = eminp+ (ipts-1d0)*(emaxp-eminp)/(ndos-1d0)
          write(ifi,"(255(f13.5,x))")eee,pdosp(ipts,1:nchanp)
        enddo
        call fclose(ifi)
      enddo
      enddo   
      deallocate(pdosp,pdosalla)
      endif
      end
!!
!! ----------------------------------
      subroutine writeband(stdo,nkp,nsyml,nspx,nsp,ix99,ndhamx,
     & xdatt,nqp_syml,nqp2n_syml,nevls,evlall,qplist,labeli,labele,alat,
     & eferm,evtop,ecbot,sname,nqps_syml,nqpe_syml,nsymln,dqsyml,
     & etolv,etolc)
!! == write band file ==
      implicit none
      integer,intent(in):: nkp,nsyml,nspx,ndhamx,nsp,  nevls(nkp,nspx),ix99,stdo,nsymln
      real(8),intent(in):: evlall(ndhamx,nspx,nkp),qplist(3,nkp),xdatt(nkp),alat  !vnow, 
      integer,intent(in):: nqp_syml(nsyml),nqp2n_syml(nsyml),nqps_syml(nsyml),nqpe_syml(nsyml)
      real(8),intent(in):: eferm,evtop,ecbot ! evtop is max of n-th band. !evbot is bottom of bands upper than n+1 
                                             ! if evtop <evbot ---> insulator 
      real(8),intent(in):: dqsyml(nsyml),etolv,etolc
      character*20,intent(in)::labeli(nsyml),labele(nsyml)
      character(*),intent(in)::sname

      integer:: ifbndo,ikp,isyml,jsp,fopnn
      character*300::filenm(2),bchar
      real(8):: rydberg=13.6058d0 !,vadd
      character*3::charnum3
      integer::ifbndsp(2),iprint,iq,i,ifile_handle,ifglt,ifbndsp_nearef(2),ifmass(2),ifmglt,ifmglt2
      character*300::aaa,addx
      real(8):: emin=-18.0,emax=16.0 !eV for default plotting.
      integer:: ikps,ne,ifi,fopn,ix,nee,ibb,ilt,imin,imax
      real(8),allocatable::diffeb(:),diff2eb(:)
      real(8)::polinta,tpiba,eee,eqq,qqq,dEdkatef
      character*100::fname,massfile,fname2
      integer,allocatable::ikpoff(:)
      real(8),allocatable::disoff(:)
      character(100),allocatable:: fnameb(:,:),fnamem(:,:,:)
      logical:: initii
      logical:: semiconband,metalband
      integer:: idat,ichangesign,ifglts(2),iqplist
      character(100)::acrossef
      character(13)::massd,mass2d
      logical:: scd
      real(8)::kef
      real(8),allocatable:: kabs(:)

      scd= evtop <ecbot
      tpiba = 2d0*4d0*datan(1d0)/alat

!! ikpoffset 
      allocate(ikpoff(nsyml+1),disoff(nsyml+1))
      ikpoff(1)=0
      disoff(1)=0
      do isyml = 2,nsyml+1
        ikpoff(isyml) = ikpoff(isyml-1)+ nqp_syml(isyml-1)+nqp2n_syml(isyml-1)
        disoff(isyml) = disoff(isyml-1)+dqsyml(isyml-1)
        write(stdo,"('ikpoff=',2i5)") isyml,ikpoff(isyml)
      enddo

!! write bandplot.glt for gnuplot
      allocate(fnameb(nsyml,nspx))
      do jsp = 1, nspx
        ifglt =ifile_handle()
        ifglts(jsp)=ifglt
        fname='bandplot.isp'//char(48+jsp)
        open(unit=ifglt, file=trim(fname)//'.glt')
        write(ifglt,'(a)')'#set terminal postscript enhanced color eps'
        write(ifglt,'(a)')'#set output "'//trim(fname)//'.band.eps"'
        write(ifglt,'(a)')'set xzeroaxis'
        write(ifglt,'(a)')'set grid'
        write(ifglt,'(a)')'set ylabel "Energy-Efermi(eV)"'
        write(ifglt,'(a)')'# This is given written in subroutine writeband in lm7K/fp/bndfp.F'
        if(nspx==1) addx='"'
        if(nspx==2) addx=' isp='//char(48+jsp)//'"'
        write(ifglt,'(a)')'set title "Band '//trim(sname)//trim(addx)
        write(ifglt,'(a,F12.5,a,F12.5,a)') 'set yrange [',emin,':',emax,']'
        write(ifglt,'(a,F12.5,a)') 'set xrange [0.0:',disoff(nsyml+1),']'
!!
        write(ifglt,'(a,$)') 'set xtics ('
        do isyml = 1,nsyml   !symmetry line start
          write(ifglt,'(a,F15.10,",\")')"'"//trim(labeli(isyml))//"'",disoff(isyml)
        enddo
        write(ifglt,'(a,F15.10,")")') "'"//trim(labele(nsyml))//"'",disoff(nsyml+1)
!!
        write(ifglt,'(a,$)') 'set x2tics ('
        do isyml = 2,nsyml   !symmetry line end
          write(ifglt,'(a,F15.10,",\")')"'"//trim(labele(isyml-1))//"'",disoff(isyml)
        enddo
        write(ifglt,'(a,F15.10,")")')"'"//trim(labele(nsyml))//"'",disoff(nsyml+1)

        write(ifglt,'(a,F12.5)') 'tpia=2*3.1415926/',alat
        write(ifglt,"('plot \')")
        do isyml = 1,nsyml
          fnameb(isyml,jsp)='bnd'//charnum3(isyml)//'.'//'spin'//char(48+jsp)
          if(isyml/=1) write(ifglt,'(",\")') 
          write(ifglt,"(a,$)")  '"'//trim(fnameb(isyml,jsp))//'" u ($2):($3) lt 1 pt 1 w lp' 
        enddo  
      enddo

!! Write bnd**.isp*.glt, open qplist.dat 
      iqplist = ifile_handle()
      open(iqplist,file='qplist.dat')
      do 4111 isyml = 1,nsyml  
        ne = nqp_syml(isyml) !+nqp2n_syml(isyml)
        do iq = 1,ne
          ikp = ikpoff(isyml)+iq
          do jsp = 1, nspx
            if (iprint() > 20.and.ix99==1) then
               write(stdo,'(" bndfp: kpt",i5," of",i5, " k jsp=",3f9.5,i2," nev=",i5)')
     &         ikp,nkp,qplist(:,ikp),jsp,nevls(ikp,jsp)
            endif
            if(iprint() >=35.and.ix99==1) then
               write(stdo,"(9f8.4)") (evlall(i,jsp,ikp), i=1,nevls(ikp,jsp))
            endif
          enddo
        enddo
        do 4113 jsp = 1, nspx   ! ispx index.
          ifbndsp(jsp) =ifile_handle()
          open(unit=ifbndsp(jsp),file=fnameb(isyml,jsp))
          write(ifbndsp(jsp),"('#  ',i5,f10.5,15x,'QPE(ev)',11x,'1st-deri')") nkp,eferm
          ibb=0
          do 5113 i=1,minval(nevls(:,:)) !band index
!! take derivatives
            allocate(diffeb(ne))
            diffeb(1)  = 0d0 
            diffeb(ne) = 0d0 
            do iq = 2,ne-1
              ikp = ikpoff(isyml)+iq
              diffeb(iq)= (evlall(i,jsp,ikp+1)-evlall(i,jsp,ikp-1))/(xdatt(ikp+1)-xdatt(ikp-1))/tpiba
            enddo    
            do iq = 1,ne
              ikp = ikpoff(isyml)+iq
              write(ifbndsp(jsp),"(i5, (1x,f12.8), 2(1x,f16.10),x,3f9.5)") 
     &         i,xdatt(ikp),(evlall(i,jsp,ikp)-eferm)*rydberg,diffeb(iq), qplist(:,ikp)

!! write qplist.dat (xaxis, q, efermi) used for band plot
              if(i==1.and.jsp==1) then
                if(iq==1.and.isyml==1) write(iqplist,"(f16.10,' ! ef')") eferm
                ikp = ikpoff(isyml)+iq
                write(iqplist,"(f16.10,x,3f9.5,' !x q')") 
     &          xdatt(ikp), qplist(:,ikp)
              endif
            
           enddo

            deallocate(diffeb)
            write(ifbndsp(jsp),*)
 5113     continue
          close(ifbndsp(jsp))
 4113   continue
 4111 continue
      close(iqplist)

c$$$
c$$$!! write massplotband.glt for gnuplot
c$$$!! write massplot.glt for gnuplot (curvature, not the slope for metal).
c$$$      allocate(fnamem(minval(nevls(:,:)),nsymln+1:nsyml,nspx))
c      do jsp = 1, nspx
c$$$c        ifmglt=ifile_handle()
c$$$c        fname='massplot.isp'//char(48+jsp)
c$$$c        open(unit=ifmglt,file=trim(fname)//'.glt')
c$$$c        ifmglt2=ifile_handle()
c$$$        fname2='massplotband.isp'//char(48+jsp)
c$$$c        open(unit=ifmglt2,file=trim(fname2)//'.glt')
c$$$c   set terminal postscript enhanced color eps
c$$$c   set output 'temp.eps'
c$$$c        write(ifmglt,'(a)')'set terminal postscript enhanced color eps'
c$$$c        write(ifmglt,'(a)')'set output "'//trim(fname)//'.eps"'
c$$$c        write(ifmglt,'(a,i1)')'set multiplot layout 1,',nsyml-nsymln
c$$$c        write(ifmglt,'(a)')'set xzeroaxis'
c$$$c        write(ifmglt,'(a)')'set grid'
c$$$c        write(ifmglt,'(a)')'set ylabel "mass/mass(electron)"'
c$$$c        write(ifmglt,'(a)')'set xlabel " |q|/(2pi/alat) "'
c$$$c        write(ifmglt,'(a)')'# If necessary, see subroutine writeband given in lm7K/fp/bndfp.F for your purpose!'
c$$$
c$$$c$$$        write(ifglt,'(a)')'set terminal postscript enhanced color eps'
c$$$c$$$        write(ifglt,'(a)')'set output "'//trim(fname2)//'.eps"'
c$$$c$$$        write(ifglt,'(a,i1)')'set multiplot layout 1,',nsyml-nsymln
c$$$c$$$        write(ifglt,'(a)')'set xzeroaxis'
c$$$c$$$        write(ifglt,'(a)')'set grid'
c$$$c$$$        write(ifglt,'(a)')'set ylabel " e - E(fermi) (eV)"'
c$$$c$$$        write(ifglt,'(a)')'set xlabel " |q| " '
c$$$c$$$        write(ifglt,'(a)')'# If necessary, see subroutine writeband given in lm7K/fp/bndfp.F for your purpose!'
c$$$c$$$        if(nspx==1) addx=''
c$$$c$$$        if(nspx==2) addx=' isp='//char(48+jsp)
c$$$c$$$c        write(ifmglt,*)
c        do isyml = nsymln+1,nsyml
c$$$          ne = nqpe_syml(isyml)-nqps_syml(isyml)+1
c$$$c$$$c          if(scd)         aaa='# insul u ($5):($8)'
c$$$c$$$c          if(.not.scd)    aaa='# metal u ($5):($9)'
c$$$c$$$c          if(.not.scd)    aaa='# metal u ($5):($8)'
c$$$c$$$c          write(ifmglt,'(a)')'set title "Mass '//trim(sname)//trim(addx)//
c$$$c$$$c     &    ' '//trim(labeli(isyml))//'--'//trim(labele(isyml))//' '//trim(aaa)//'"'
c$$$c$$$c          write(ifmglt,'(a,F12.5,a)') 'set xrange [0.0:',xdatt(ikpoff(isyml+1))-disoff(isyml),']'
c$$$c$$$c          write(ifmglt,"('plot \')")
c$$$c$$$          write(ifglt,'(a)')'set title "Band '//trim(sname)//trim(addx)//
c$$$c$$$     &    ' '//trim(labeli(isyml))//'--'//trim(labele(isyml))//'"'
c$$$c$$$          write(ifglt,'(a,F12.5,a)') 'set xrange [0.0:',xdatt(ikpoff(isyml+1))-disoff(isyml),']'
c$$$          write(ifglt,'(a,F12.5)') 'tpia=2*3.1415926/',alat
c$$$          write(ifglt,"('plot \')")
c$$$          ibb=0
c$$$          ikps= ikpoff(isyml)+1
c$$$          initii=.true.
c          do i=1,minval(nevls(:,:))
c$$$            eee = evlall(i,jsp,ikps)
c$$$            semiconband = scd .and. evtop-etolv<eee .and. eee<ecbot+etolv .and. isyml>nsymln !check i-th band is neare VCM and CBM
c$$$            idat= ichangesign(evlall(i,jsp,ikps:ikps+ne-1)-eferm,ne) !for metal crosspoint point across eferm
c$$$            metalband  = idat>-1 .and. isyml>nsymln !check i-th band is near VCM and CBM
c$$$            if(semiconband.or.metalband) then
c               ibb = ibb+1
c               fnamem(ibb,isyml,jsp)='Band'//charnum3(ibb)//'Syml'//charnum3(isyml)//'Spin'//char(48+jsp)//'.mass'
c$$$c               if(.not.initii) write(ifmglt,'(",\")') 
c$$$c               if(scd)         aaa='" u ($5):($8) lt'
c$$$c               if(.not.scd)    aaa='" u ($5):($9) lt'
c$$$c               write(bchar,"(a,i2,a,i2,a,i3,a)") '"'//trim(fnamem(ibb,isyml,jsp))//trim(aaa),
c$$$c     &         ibb,' pt ',ibb,' w lp ti "band=',i,'"'
c$$$c               write(ifmglt,"(a,$)")trim(bchar)
c$$$               if(.not.initii) write(ifglt,'(",\")') 
c$$$               aaa='" u (tpia*$5):($6) lt'
c$$$               write(bchar,"(a,i2,a,i2,a,i3,a)") '"'//trim(fnamem(ibb,isyml,jsp))//trim(aaa),
c$$$     &         ibb,' pt ',ibb,' w lp ti "band=',i,'"'
c$$$               write(ifglt,"(a,$)")trim(bchar)
c$$$               initii=.false.
c$$$            endif
c          enddo
c$$$c          write(ifmglt,*)
c$$$c          write(ifmglt,'(a)')'unset ylabel'
c$$$          write(ifglt,*)
c$$$          write(ifglt,'(a)')'unset ylabel'
c$$$        enddo  
c$$$c        close(ifmglt)
c$$$c        close(ifmglt2)
c      enddo


!! Write Band*Syml*Spin*.dat : These are between [VBM(left)-etolv(Ry), CBM(left) + etolc(Ry)]
      allocate(fnamem(minval(nevls(:,:)),nsyml,nspx))
      do jsp = 1, nspx
      do isyml = 1,nsyml
        ibb=0
        do i=1,minval(nevls(:,:))
         ibb = ibb+1
         fnamem(ibb,isyml,jsp)=
     &    'Band'//charnum3(ibb)//'Syml'//charnum3(isyml)//'Spin'//char(48+jsp)//'.dat'
        enddo
      enddo
      enddo
      do 2111 isyml = 1,nsyml  
        ne = nqp2n_syml(isyml)
        if(ne==0) cycle
        allocate(diffeb(ne),diff2eb(ne))
!! diffeb and diff2eb are 1st and 2nd derivatives, 
!! dEi(k)/dk d^2Ei(k)/d^2k along a symmetry line (isyml index).
!! But be careful to determine effective mass, because it is not analytic at Gamma point.
!! E.g, see Kittel's book.
        do 2112 jsp = 1, nspx   ! ispx index.
          ibb=0
          do 2113 i = 1,minval(nevls(:,:)) !band index
!! take derivatives
            diffeb(1)  = 0d0 
            diffeb(ne) = 0d0 
            diff2eb(1) = 0d0 
            diff2eb(ne)= 0d0 
            do iq = 2,ne-1
              ikp = ikpoff(isyml)+ iq + nqp_syml(isyml)
              diffeb(iq)= (evlall(i,jsp,ikp+1)-evlall(i,jsp,ikp-1))/(xdatt(ikp+1)-xdatt(ikp-1))/tpiba
            enddo    
!! NOTE: effective mass. 2d0 is because emass in Ry unit is 1/2
            ikps= ikpoff(isyml)+1 + nqp_syml(isyml)
            eee = evlall(i,jsp,ikps)
c            print *,'jjjjjjjj',jsp,isyml,ikps,i,eee
            semiconband = scd.and.evtop-etolv<eee .and. eee<ecbot+etolv 
              !check i-th band is neare VCM and CBM
            metalband = ichangesign(evlall(i,jsp,ikps:ikps+ne-1)-eferm,ne) >-1 
              !for metal crosspoint point across eferm
              !check i-th band is near VCM and CBM
            if(semiconband.or.metalband) then
               ifmass(jsp) = ifile_handle()
               ibb = ibb+1
               if(metalband) then
                 imin=max(2,idat-1)  !diffeb is given from 2 to ne-1
                 imax=min(ne-1,idat+2)
                 allocate( kabs(imin:imax))
                 kabs(imin:imax) = (xdatt(ikps+imin-1:ikps+imax-1)-disoff(isyml))*tpiba ! tpiba is 2pi/alat. |k|. left-end is zero.
                 kef = polinta(0d0, evlall(i,jsp,ikps+imin-1:ikps+imax-1)-eferm, kabs, imax-imin+1) !we have k at Ef |k|=qef
                 dEdkatef = polinta(kef, kabs, diffeb(imin:imax), imax-imin+1)
                 deallocate(kabs)
c                 write(stdo,"('Effective Mass for metal: ',a,
c     &            ' iq isyml |k|_ef/(2pi/alat) dEdk 2*k/dEdk=', 2i4,' ',f8.3,' ',f8.3,' ',f8.3)") 
c     &            trim(fnamem(ibb,isyml,jsp)),idat,isyml, kef/tpiba, dEdkatef, 2d0*kef/dEdkatef
               endif
c               print *,'ffff ', ibb,isyml,jsp,trim(fnamem(ibb,isyml,jsp))
               open(unit=ifmass(jsp),file=trim(fnamem(ibb,isyml,jsp)))
               write(ifmass(jsp),"(a)") '# mass1 is for metal, mass2 is for semiconductor'
               write(ifmass(jsp),"(a)") '# isyml,ib,iq,isp, |k|/(2pi*alat), QPE-EF, QPE-QPE(start),'
     &              //' mass2=2*(2*(QPE-QPE(start))/|k|**2), mass1=2*|k|/(dE/dk)'
               do ix=1,ne
                 acrossef=''
                 if(ix==idat)   then
                   write(bchar,"(f13.8,' ',f8.3)") kef/tpiba,2d0*kef/dEdkatef
                   acrossef=' <--- Ef. kef/(2pi/alat), 2*|k|/(dE/dk)_kef= '//trim(bchar)
                 endif  
                 if(ix==idat+1) acrossef=' <--- Ef.'
                 qqq= (xdatt(ikps+ix-1)-disoff(isyml))*tpiba ! tpiba is 2pi/alat. |k|. left-end is zero.
                 eqq= (evlall(i,jsp,ikps+ix-1)-evlall(i,jsp,ikps)) 
                              ! in rydberg !eqq is relative to the eval at ikps(left end point)
                 if(ix==1.or.ix==ne) then
                   massd= '   --------  '
                   mass2d='   --------  '
                 else
                   write(massd,"(f13.8)")  2d0*qqq/diffeb(ix)
                   write(mass2d,"(f13.8)") 2d0/(2d0*eqq/qqq**2)
                 endif
                 write(ifmass(jsp),"(3i4,i2,' ',f13.8,' ',f13.8,' ',f13.8,' ',a,' ',a,a)")
     &           isyml,i,ix,jsp,(xdatt(ikps+ix-1)-disoff(isyml))*tpiba, 
     &           (evlall(i,jsp,ikps+ix-1)-eferm)*rydberg,
     &           (evlall(i,jsp,ikps+ix-1)-evlall(i,jsp,ikps))*rydberg,
     &           mass2d,  massd,  trim(acrossef)
c    &         isyml,i,jsp,(evlall(i,jsp,ikps)- eferm)*rydberg, (2d0/(2d0*eqq/qqq**2) xxxxxxxxxxxxxx
               enddo  
               close(ifmass(jsp))
               write(aaa,"(a,f13.5,a)") '" u ($5/tpia+',disoff(isyml),'):($6) lt'
               write(bchar,"(a,i2,a,i2,a,i3,a)") '"'//trim(fnamem(ibb,isyml,jsp))//trim(aaa),
     &          ibb+1,' pt ',ibb+1,' w lp'
               write(ifglts(jsp),'(",\")') 
               write(ifglts(jsp),"(a,$)")trim(bchar)
            endif
 2113     continue
 2112   continue
        deallocate(diffeb,diff2eb)
 2111 continue

      do jsp=1,nspx
      write(ifglts(jsp),"(/,a)")'pause -1'
      write(ifglts(jsp),*)
      close(ifglts(jsp))
      enddo

!! bnds.* is only for backward compatibility.
      ifbndo = fopnn('bnds')
      write(ifbndo,"(i5,f10.5,i6)") sum(nqp_syml(1:nsyml)),eferm,0
      do 3111 isyml = 1,nsyml  
        ne = nqp_syml(isyml)
        write(ifbndo,"(i5)") ne
        do iq = 1,ne
          ikp = ikpoff(isyml) + iq
          do jsp = 1, nspx
            write(ifbndo,"(3f10.5,i6)") qplist(:,ikp),nevls(ikp,jsp)
            write(ifbndo,"(10f8.4)")(evlall(i,jsp,ikp),i=1,nevls(ikp,jsp))
          enddo
        enddo
 3111 continue
      write(ifbndo,"(i6)") 0
      call fclose(ifbndo)
      end

