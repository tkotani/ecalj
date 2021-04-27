!! = Main part of full-potential LDA/GGA/QSGW(for given sigm). Single iteration =
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
      subroutine bndfp(nlibu,lmaxu,lldau, ssite,sspec, !sstrnmix, !,sbz,sctrl,sham,slat spot,
     .     leks,lrout,lfrce,lpnu,dmxp,iter,maxit, frc, !ndham,
     .     dmatu,vorb,llmfgw)
      
      use m_supot,only:  zv_a_osmrho,zv_a_osmpot,sv_a_oorhat
      use m_suham,only:iv_a_oindxo, pot_nlma, pot_nlml, ndham=>ham_ndham
c      use m_seneinterp,only:,seneinterp,rv_a_ohrs,iv_a_oiaxs,iv_a_ontabs,rv_p_oqsig
      use m_hamindex, only: napwmx,ngrp,norbmto,ibastab,ltab,ktab,offl, symops_af,symops
      use m_lmfinit, only: ncutovl,lso,bz_nabc,bz_ndos,bz_def,bz_ef,bz_efmax,bz_fsmom,
     &     bz_nkp,bz_dosmax,bz_ntet,bz_w,bz_dosw,bz_nevmx,bz_lmet,bz_fsmommethod,bz_n,
     &     ctrl_nspec,ctrl_pfloat,ctrl_ldos,ctrl_zbak,ctrl_lfrce,ctrl_nbasp, !ctrl_nl,
     &     ham_ldham,ham_pwmode,ham_pwemin,ham_pwemax,ham_lsig,ham_oveps,ham_eterms,ham_elind,ham_rsrnge,
     &     ham_ndhrs,ham_scaledsigma,ham_pmin,ham_pmax,ham_ehf,ham_ehk, !ham_nqsig,
     &     lat_vol,lat_alat,lat_plat,lat_qlat,lat_nabc,rv_a_opos ,stdo,stdl,procid,master,
     &     nsp,nspc,nkaph,nlmax,nl,nbas,nsp, ham_frzwf,
     &     sstrnmix
c      use m_lmfinit,only: globalvariables
      use m_subzi, only: subzi
      use m_rsibl_ev, only : rsibl_ev ! to plot wavefunction in the fat band mode
      use m_w_psir, only: w_psir_cmdopt
      use m_struc_def
      use m_ext,only: sname     !file extension. Open a file like file='ctrl.'//trim(sname)
      use m_mkqp,only: iv_a_ostar,rv_a_owtkp,rv_p_oqp,iv_a_oipq,iv_a_oidtet
      use m_seneinterp,only: Senebroadcast,Rdsigm2,ndimsig,nspsigm,Sigmainit,Getsene,sene,Dsene
!! ==== followings are just notations in this code (in cases wrong).
Ci   nbas  : size of basis
Ci   nsp   : number of spins
Ci   nlibu : total number of LDA+U blocks (used to dimension dmatu and vorb)
Ci   lmaxu : lmax for U used to dimension vorb and dmatu. lmaxu=2 if d is, but lmaxu=3 if f is included.
Ci   lldau :lldau(ib)=0 => no U on this site otherwise 
Ci         :U on site ib with dmat in dmats(*,lldau(ib))
Ci   Following are old. Moved to lmfinit
Ci     ssite :struct for site-specific information; see routine usite
Ci     sspec :struct for species-specific information; see routine uspec
Ci     slat  :struct for lattice information; see routine ulat
Ci     sctrl :struct for program flow parameters; see routine uctrl
Ci     sham  :struct for parameters defining hamiltonian; see routine uham
Ci     spot  :struct for information about the potential; see routine upot
Ci     sbz   :struct for the Brillouin Zone; see routine ubz
c      
Ci   ndham :dimensioning parameter, at least as large as largest
Ci         :hamiltonian dimension
Ci   leks  :>0 make the Hohnberg-Kohn-Sham energy
Ci         :>1 use the HKS forces
Ci   lrout :>0 generate output density and attendant quantities
Ci   lfrce : 0 suppress generation of forces
Ci   lpnu  : 1 make new pnu's
Ci   sstrnmix:mixing
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
ccCl   lwndow:T if to make density in a specified energy window
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
      
Cr Remarks ---- now improving...
Cr   Band pass consists of:
Cr   (1) make the effective potential,
Cr   (2) generate eigenvalues (and eigenvectors if lrout)
Cr   (3) if lrout, assemble the output density by BZ integration
Cr   (4) evaluate hf (and KS, if leks) energy by BZ integration
Cr   (5) mix the output density to make a new input density.

      
Cu Updates before github 2009. See history github ecalj after 2009.
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
Cuxxx   24 May 03 New --window switch for accumulating density
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
#endif
c#if MPIK
      integer obuf,oqptr
      integer igetss,lmxh,nlmh,nelt(3),kmax !,nlma
      integer, dimension(:),allocatable :: kpproc
      double precision sttime,entime
c#endif
c#if MPIK
      integer:: pidorigin, ndimhx_iq, nev_iq,idx
      integer:: work_mpi(10)    ! use only 1 and 2 now
c#endif
      logical mlog
c      integer procid,master
      integer nbx,n0,nab,nppn,iter,maxit!,ndham !nbas
c      character*(*) sstrn
      parameter ( nbx=512, n0=10, nppn=12, nab=9)
      integer k1,k2,k3,leks,lrout,lfrce,lpnu
      integer nlibu,lmaxu,lldau(nbas)
      real(8):: dmxp(33),frc(3,*)
c      type(s_sstrn):: sstrn
c      character*(*):: sstrnmix
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)
      logical lgors,ltet,cmdopt!,lwndow
      character strn*120,plbopt*120,mulopt*120,clsopt*120,strn2*120,dc*1
ckino lswtk must be 0 initially.
      integer:: ipr , iprint , ipl , isp , jsp , ispc , isum , nglob 
     .  , i , iq , ismidb ,  ldos , lmet !, nspc  , nsp
     .  , lrep , lcplxp , nkabc(3) , ndos , ndos0 , nev , nevl 
     .  , nkp  , numq , n1 , n2 , n3 , ntet , lwtkb , lswtk=0 !, nl 
     .  , lgunit , igets , mpsord , nevmx , nvl , ngabc(3) , isw  !stdo , stdl 
     .  , ifi , fxst , fopn , fopna , lfrzw , i1 , i2 , plbnd 
     .  , nfilem , iobzwt , lrsig ,  jobgw , nll, lnoxc !, lpdiag lwsig  , 
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
c      logical:: l_temp_rv_a_oswtk=.false., l_temp_rv_a_owtkb=.false.
c     &  l_temp_iv_a_oidtet=.false.
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
      integer::  i1mach , icls , isite(nsitmx) , iclsl(nsitmx) 
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
c      integer nqsig
c      integer:: nbasp
      data T /.true./ F /.false./
      integer,allocatable:: nevls(:,:)
ctakao
      integer :: ifiogw,ifiese,og,ib,ig,iorb,is,nsgrp,uspec_size,nspec
      logical :: nexist,sigmamode,ok !, newsigmasw
      character spid*8
      integer :: ib1,ib2,ii,ik1,ik2,ik3,j1,j2,ikx
      integer:: nhmax , nknknk , jx , nsp_dummy , ndimh_dummy
      complex(8),allocatable:: hamm(:,:,:),ovlm(:,:,:) ,senex(:,:),ovlmtoi(:,:),ovliovl(:,:) 
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
      integer:: getef, i_dummy_getef,ifig!, ndimsig
c      integer:: n_listwf
      integer, allocatable:: listwf(:)
      integer:: ntqxx,ispx,ix99 !,nevec unused.
      real(8)::eseavrmean
c      logical:: mtosigmaonly
      real(8)::scaledsigma
c      real(8):: qsfz(3,100)
      integer:: nqps,ii1,ii2,ii4,iq1 !, fxsts
      logical:: lfx

      integer:: unlink
!! for mpi and mpik, parallel write
      character(10):: strprocid 
      integer:: nsizeohrs,nsizeoiaxs,nsizeontabs,noqsig
      integer:: ifis2,ipts,ifid
      real(8):: eee
      logical :: tdos !,pdos
      integer:: ifves,ifvesintatm,ifvesintloc,ispec,ifvbm,ifqbyl,il

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

      integer:: ispsigm !nov2014
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
      logical,save:: siginit=.true.     !,laf
      character(10):: i2char
      character(50)::infoq
      integer::ipos,ndos_,j,ififft
      real(8)::dosmax,emaxp_,posp(3)
      integer:: jspend,jspini,jspp
      integer:: iqread,iqindex,ibread(100)=0,ibindex,nibread, job,ifih,ist,ip,ni,ix

      integer,allocatable:: npair(:,:), nlat(:,:,:,:),ib_table(:),k_table(:),l_table(:)
      integer::npairmx,ikpd,ifqplist
      integer,allocatable:: nqwgt(:,:,:)
      real(8):: tpioa,pi=4d0*datan(1d0)
      complex(8):: aaaa,img=(0d0,1d0)
      character*4:: cccx

      real(8):: pnu(n0,2),pnz(n0,2),pmean,symxx(3,3),qvec(3),dee,dosi(2)
      integer:: i_copy_size, ndummy1,ndummy2,ndummy3,ncount,nbandx,ifii
      logical:: lpz,phispinsym,ljobgw
      
!! ----------------------------
      call tcn ('bndfp')
      debug    = cmdopt('--debugbndfp',6,0,strn)
      mlog     = cmdopt('--mlog',6,0,strn) !--mlog here is taken by getarg.
      tdos     = cmdopt('--tdos',6,0,strn)
      PROCARon = cmdopt('--mkprocar',10,0,strn) !write PROCAR(vasp format).
      fullmesh = cmdopt('--fullmesh',10,0,strn) !full mesh stop just after do 2010 iq loop.
      fsmode   = cmdopt('--fermisurface',14,0,strn)
      allband  = cmdopt('--allband',9,0,strn)
      if(fullmesh.and.(.not.procaron)) call rx(' --fullmesh is meaningful only with --mkprocar')
#if MPI | MPIK
      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
#endif      
!! Write atomic positions
      if(procid == master) then
         write(stdo,"(a)")' --- bndfp: Atomic Positions from rst/ctrl. '
         write(stdo,"(a)")'   Note: rst and ctrl may contain different atomic positions when DYN mode'
         do ib=1,nbas
            is  = ssite(ib)%spec   !species index
            spid= ssite(is)%clabel !atom label
            write(stdo,"('ib clabl pos =',i5,a6,3(x,f13.5))") ib, trim(spid),ssite(ib)%pos
         enddo
      endif
!! Fermi surfece mode for xcrysden. It's format is in  http://www.xcrysden.org/doc/XSF.html#2l.16
      if(fsmode) fullmesh=.true.
      if(fullmesh) then
         nkk1=bz_nabc(1)
         nkk2=bz_nabc(2)
         nkk3=bz_nabc(3)
      endif   
!! Calculate VBM and CBM relative to vaccum (but in a simple approximaiton. need fixing.).
      if(cmdopt('--vbmonly',9,0,strn)) goto 8888
!! getq mode.  Current version is not for spin dependent nor many restrictions.
!! spin symmetic (or nspin=1, not 2 channell binded and so on...
      if(cmdopt('--getq',6,0,strn)) then
         if(procid==master) call getqmode()
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
      getz =.false.
      if (maxit .eq. 0) return
c      siginit=.true. !commented out on mar2015. this makes reading sigma at every iteration.
      napw = 0
      ipr  = iprint() !When procid/=master, iprint=0.
      ipl  = ipr
!! lso  =0 (no so): =1(L.S): =2(LzSz).
      lcplxp = 0
      lekkl  = ctrl_pfloat
      if (lso .ne. 0)  lcplxp = 1
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
!! related to the error in rseq. Remove RSEQ_ERROR at first because
!! this error is only fatal when it is finaly converged. takaoNov28_2010
      idummy = unlink('RSEQ_ERROR') 
ctakao sep2010 !complicated.
c      nspc=2 for non-colinear case. 
c      allocate(evl(ndham,2))
c      call i8icopy(size(ham_ldham),ham_ldham,1,ldham,1)
      allocate(evl(ndhamx,nspx))
      ldham = ham_ldham
      ldim  = ldham(1)
      if(size(ham_ldham)/=size(ldham))call rx("bndfp:size(ham_ldham)/=size(ldham)")
      pwmode = ham_pwmode
      if (procid .eq. master) then
         ifiese=ifile_handle()
         open(ifiese,file='PWMODE')
         write(ifiese,*)pwmode
         close(ifiese)
      endif 
      pwemin = ham_pwemin
      pwemax = ham_pwemax
      onesp = 0
      ldos  = ctrl_ldos !0 or 1
      lrsig = ham_lsig  !12 for read sigma mode
      if(ham_frzwf)then
         lfrzw = 1  !freeze all augmentation wave
      else
         lfrzw = 0
      endif
      nspec = ctrl_nspec
      if (procid .eq. master) then
         inquire(file='sigm.'//trim(sname),exist=nexist)
         if (lrsig .ne. 0 .and. (.not.nexist)) then
            write(6,*)' bndfp (warning): no sigm file found ... LDA calculation only'
            ham_lsig=0
            lrsig = 0
         endif
      endif
      call mpibc1(lrsig,1,2,mlog,'bndfp','lrsig')
      call mpibc1(ham_lsig,1,2,mlog,'bndfp','ham_lsig') 
      epsovl = ham_oveps
      if (procid .eq. master) then
        aaachar=trim(i2char(iter))//" of "//trim(i2char(maxit))
        write(6,"(a)") trim(" --- BNDFP:  begin iteration "//aaachar)
      endif
      if (nbas.gt.nbx) call rxi('bndfp: nbx exceeded, need',nbas)
      if (ndham.le.0)  call rx('bndfp: hamiltonian matrix has zero dimension')
      eterms=-99d0
      eterms(19) = 0d0
      ham_eterms=eterms
      nvl  =pot_nlml
      nchan=pot_nlma
      eks = 0
      alat =lat_alat
      plat =lat_plat
      qlat =lat_qlat
      ngabc=lat_nabc
      vol = lat_vol
      k1 = n1
      k2 = n2
      k3 = n3
!! ... for BZ integration
      nkabc= bz_nabc
      nkp=   bz_nkp
      ntet=  bz_ntet
      lmet=  bz_lmet 
      mpsord=bz_n
      esmear=bz_w
      nevmx= bz_nevmx
      efmax= bz_efmax
      fsmom= bz_fsmom
      if(debug) print *,'ttt nevmx w=',nevmx,esmear
      esmear = esmear + iabs(mpsord)
      if (mpsord .lt. 0) esmear = -esmear
      mpsord = int(dabs(esmear))
      if (esmear .lt. 0) mpsord = -mpsord
      ndos= bz_ndos
      dosw = bz_dosw
      ef0  = bz_ef
      def  = bz_def
      ltet = ntet .gt. 0
      if (lmet /= 0 .and. .not. ltet) ldos=1 ! ldot=1(T) when Metal=T and Tetra=F
                                             ! We need DOS if just sampling
      qbg = ctrl_zbak(1)
      alfa = 0
C     Switch to plot bands at specified qp
      if (cmdopt('--band',6,0,strn).or.fullmesh) then
         plbnd = 1
         lrout = 0
         lfrce = 0
         nkp = 0
         numq = 1
         allocate(ifbls_iv(ndhamx*2))
         ifbls_iv(:)=0
      else
         plbnd = 0
      endif
!! == PROCAR file 
      if(PROCARon) then
         iprocar1=ifile_handle()
         open(iprocar1,file='PROCAR.UP.'//trim(strprocid))
         if(nsp==2) then
            iprocar2=ifile_handle()
            open(iprocar2,file='PROCAR.DN.'//trim(strprocid))
         endif
      endif
!! spin-symmetrized pnu. aug2019. See also in pnunew and locpot
      phispinsym = cmdopt('--phispinsym',12,0,strn)
      if(phispinsym) then
         if(master_mpi())write(6,*)' --phispinsym use spin-averaged potential for phi and phidot' 
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
!! Write out Hamiltonian stated at jan2018kaist
      writeham= cmdopt('--writeham',10,0,strn)
      if(writeham) then
        open(newunit=ifih,file='HamiltonianPMT.'//trim(strprocid),form='unformatted')
        if(procid == master) then ! pair-table and weight generator.  nlat,qgt are not used yet...
          nkk1=bz_nabc(1)
          nkk2=bz_nabc(2)
          nkk3=bz_nabc(3)
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
          allocate(npair(nbas,nbas),pos(3,nbas))
          call dcopy(nbas*3 , rv_a_opos , 1 , pos , 1 )
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
c             if(abs(sum(qwgt(1:npair(ib1,ib2),ib1,ib2))-1d0)>1d-8) call rx('buggg:qwgt sum is not unity')
              do ni = 1,npair(ib1,ib2) 
                posp =  pos(:,ib1)-pos(:,ib2) + matmul(plat,nlat(:,ni,ib1,ib2)) ! R_i+T - R_j
                rrrr = sqrt(sum(posp**2))
                write(6,"(i6,3x,3i4,i3,x,f8.3)") ni, nlat(1:3,ni,ib1,ib2),nqwgt(ni,ib1,ib2),rrrr
              enddo
          enddo
          enddo
          open(newunit=ififft,file='HamiltonianPMTInfo',form='unformatted')
          write(ififft) plat,nkk1,nkk2,nkk3,nbas,qlat
          write(ififft) pos,alat
          write(ififft) qplist
          write(ififft) npair,npairmx
          write(ififft) nlat,nqwgt
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
        endif
      endif
!! --- Define local arrays used in the generation of the potential ---
      allocate(qmom_rv(nvl))
      allocate(gpot0_rv(nvl))
      allocate(vval_rv(nchan))
      allocate(hab_rv(nab*n0*nsp*nbas))
      allocate(vab_rv(nab*n0*nsp*nbas))
      allocate(sab_rv(nab*n0*nsp*nbas))
      allocate(ppnl_rv(nppn,n0,nsp,nbas))
      allocate(fes1_rv(3*nbas))
!! dfaugm alllocates and initializes variables
      call dfaugm ( nbas , lcplxp , lso , ssite , sspec , sv_p_osig, sv_p_otau , sv_p_oppi )
!! --- Make the potential sans XC part ---
      jobgw=-999
      if(llmfgw) then ! only for lmfgw program
         if(procid==master) then
            write(stdo,*) ' lmfgw: input one of the following jobs:'
            write(stdo,*) '   0 : init mode; creates files SYMOPS, LATTC, CLASS, NLAindx, ldima'
            write(stdo,*) '   1 : GW setup mode'
            write(stdo,*) ' job 0 or 1?'
            read (5,*) jobgw
            if(.not.( jobgw==0.or.jobgw==1)) call rx(' lmf-MPIK is only for jobgw=0,1 now')
         endif
!! --- generate potential without XC. This mkep to obtain
         if(jobgw/=0) then
            write(stdo,"(a)")' Make potential without XC part ...'
            lnoxc=1
            i = 1 + 10*lfrzw + 100*lnoxc ! Adding 100 means excluding XC(LDA) part. nolxc=T
              ! We obtain sv_p_osigx , sv_p_otaux , sv_p_oppix smpotx without XC by the 'call mkpot'.
              ! Other output are overwritten by the 2nd call of mkpot 
            call dfaugm(nbas, lcplxp , lso, ssite, sspec, sv_p_osigx, sv_p_otaux, sv_p_oppix)
            allocate(spotx(k1,k2,k3,nsp))
            spotx=0d0
            call mkpot ( nbas , ssite , sspec ,   0 , 
     .           lcplxp , k1 , k2 , k3 , zv_a_osmrho , sv_a_oorhat , qbg , spotx 
     .           , qmom_rv , vconst , vrmt , sv_p_osigx , sv_p_otaux , sv_p_oppix 
     .           , ppnl_rv , hab_rv , vab_rv , sab_rv , qval , qsc , gpot0_rv 
     .           , vval_rv , fes1_rv , i , vorb , nlibu , lmaxu , lldau , nsp,0 )
         endif
#if MPI|MPIK
         call mpi_barrier(MPI_comm_world,ierr)
         call mpibc1(jobgw,1,2, mlog,'bndfp','jobgw')
#endif        
      endif
C --- Make the potential and total energy terms ---
      i = 1 + 10*lfrzw !+ 100*(lnoxc=0)
      if (cmdopt('--wrhomt',8,0,strn)) then
         i = i + 10000 !write rhomt mode
      else if (cmdopt('--wpotmt',8,0,strn)) then
         i = i + 20000 !write potmt mode
      endif
      if(llmfgw) i = i + 10000
!! mkpot calls locpot. and locpot calls augmat. augmat calculates sig,tau,ppi.
      call togpr()
      call mkpot ( nbas , ssite , sspec ,   lfrce!sham ,slat ,spot , 
     .  , lcplxp , k1 , k2 , k3 , zv_a_osmrho , sv_a_oorhat , qbg , zv_a_osmpot 
     .  , qmom_rv , vconst , vrmt , sv_p_osig , sv_p_otau , sv_p_oppi 
     .  , ppnl_rv , hab_rv , vab_rv , sab_rv , qval , qsc , gpot0_rv 
     .  , vval_rv , fes1_rv , i , vorb , nlibu , lmaxu , lldau , nsp,0 )
      call togpr()
      elind = ham_elind
      if (elind .lt. 0d0) elind=-(3*pi**2*(qval-qsc-qbg)/vol)**.66666d0*elind
      dmxp(33) = elind
      allocate( evlall(ndhamx,nspx,nkp))
      evlall = 0d0 

!! --- Setup for BZ integration ---
      if (plbnd .eq. 0) then
         call subzi(lmet,ltet,lrout.gt.0,ndham,nsp,nspc,nkp,qval-qbg,
     .     nevmx,lwtkb,eferm,rv_a_owtkb,numq,lswtk,rv_a_oswtk,ef0,def,mpsord,
     .     dabs(esmear-mpsord))
         write(stdo,*)' Start first of two band passes ...'
      else
         nkp = 0
         ldos = 0
         lwtkb = -1
         icls = 0
      endif
c      if(lwtkb /= -1) call rx('bndfp: we need lwtkb=-1')
c      if(numq/=1)     call rx('bndfp: we need numq=1')
      if(.not.(tdos .or. PROCARon)) then !if(ndimh*nspc /= nevl) then 
         write(stdo,*)' bndfp: lswtk=-2 ... spin weights not evaluated'
         lswtk = -2
      endif
!! --- More k-point independent local arrays ---
      if (ldos .eq. 0 .and. lmet .eq. 0) ndos = 1
      if (ndos .lt. 0) ndos = -ndos
      allocate(dos_rv(ndos,2*nsp))
      if (lrout .ne. 0) then
         call dfqkkl ( nbas , ssite , sspec , numq , sv_p_oqkkl )
         if ( lekkl .eq. 1 ) call dfqkkl ( nbas , ssite , sspec , numq , sv_p_oeqkkl )
         allocate(srout_zv(k1*k2*k3*numq*nsp))
         allocate(fh_rv(3*nbas))
         allocate(fes2_rv(3*nbas))
      endif
c      if(allocated(frc)) deallocate(frc)
      call suham2(sspec,ssite)
!! --- Options for core level specta (CLS) ---
      if (cmdopt('--cls',5,0,strn)) then
         if (lmet /=2) call rx('For CLS restart with METAL=2')
         icls = 1
         clsopt = strn(6:)
         efmax = 1d3
         nevmx = ndhamx
         if (lrout .eq. 0) call rx('bndfp: need output density for cls')
         call suclst(nsitmx,nbas,nsp,ssite,sspec,clsopt,
     .     isite,iclsl,iclsn,nsites)
         efmax = 1d3
         nevmx = ndhamx
         if (16*3*nlmax*ndham*nbas*nsp*nkp/1000000 .gt. 24 .and. procid.eq. master) then
!! this path is go through by crn test.
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

!! --- get index for hamiltonian for m_hamindex takao june2009
!! these are used in sigm mode(QSGW).
!! memo:
!!  ib = atom index
!!  ltab= L (angular momentum index)
!!  ktab=  =1 for EH, =2 for EH2, =3 for lo
!! SYMOPS
      if(onlyonce) then 
         onlyonce=.false.
         if(procid==master) then
            call gen_hamindex(ssite, sspec, ldim, .true., jobgw==0)
c           write(6,*) ' ib l  k offl(iorb)+1  offl(iorb)+2*l+1  trim(spec)'
            allocate(ib_table(ldim),l_table(ldim),k_table(ldim))
            do iorb = 1, norbmto !Total number of MTO's
               ib   = ibastab(iorb)
               is   = ssite(ib)%spec
               spid=sspec(is)%name
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
         endif
      endif
!!
C ... Case only generate bands at supplied qp: setup
      if (plbnd .ne. 0) then
         if (procid .eq. master) then !Try and read Fermi level from file
            inquire(file='wkp.'//trim(sname),exist=nexist)
            if(nexist) then
               open(newunit=ifi,file='wkp.'//trim(sname),form='unformatted')
               read(ifi) ndummy1,ndummy2,ndummy3,ef0
               close(ifi)
            endif
         endif   
         call mpibc1(ef0,1,4,.false.,'bndfp','ef0')
      endif   
!! ... Setup for spin-orbit coupling
      if (lso /= 0) then
c         if(allocated(orbtm_rv)) deallocate(orbtm_rv)
         allocate(orbtm_rv(nl*nsp*nbas))
      endif
!     !
!! === Get nkp and qplist =========== feb2015
!! fatband mode ---> --mkprocar
!! plbnd/=0 ---> band plot mode.
      if(debug) write(*,*) 'bndfp: uuuuuuu1111111'
      if(allocated(qplist)) deallocate(qplist)
      if( procid==master .and. plbnd==0) then
         allocate(qplist(3,nkp))
         call dcopy(3*nkp, rv_p_oqp,1,qplist,1)
c-----------------
         open(newunit=ifqplist,file='QPLIST.IBZ')
         do iq=1,nkp
            write(ifqplist,"(i5,3f23.15,3x,3f23.15)")iq, qplist(:,iq), rv_a_owtkp(iq)/2d0
         enddo !!I think sbz%rv_p_oqp is not used in the following.
         close(ifqplist)
      elseif(procid==master) then
         if (cmdopt('--onesp',7,0,strn) .and. nspc .eq. 1) onesp = 1
         if(debug) print *,'vvvvv nkp=',nkp
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
      if (procid.eq.master .and.oncewrite(11)) then
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

      if(PROCARon) then         !if(plbnd/=0) then ! march2015
        if(procid/=master) allocate(xdatt(nkp))
        call mpibc1( xdatt, nkp , 4 , .false. , 'bndfp' , 'xdatt'  )
      endif  
#if MPIK
      call mpi_barrier(MPI_COMM_WORLD,ierr) 
#endif
      if (allocated(evlall)) deallocate(evlall)
      allocate(evlall(ndhamx,nspx,nkp)) ; evlall=0d0
      nchanp=25 !total of s,p,d,f
      if(procaron.and.fullmesh) then
         allocate(dwgtall(nchanp,nbas,ndhamx,nsp,nkp))
         dwgtall=0d0
      endif
      
!! === \Sigma-Vxc to ham_rv_a_ohrs. Neighbor table ham_iv_a_oiaxs ====
!! ndimsig is the dimension of the self-energy. 
!! ndimsig = ldim means we use only projection onto MTO spaces even when PMT. !aug2012takao
!! --- New sigma inerpolation mode. Simplified.
!! OUTPUT in this block is used in
!!   call bloch 
!!  where we use ham_iv_a_oiaxs , ham_rv_a_ohrs, ham_ndhrs and so on.
!!  In future, we should make this block to a subroutine which shows clear i/o.
      sigmamode = mod(lrsig,10) .ne. 0
      if(procid==master.and.ix99==1) print *,'sigmamode=',sigmamode
      if( llmfgw.and.jobgw==0 ) then
         continue
      elseif( sigmamode.and.siginit ) then
         call sigmainit(ldim,procid,master)
         siginit=.false. !this is equivalent to broad casting.
      endif                     ! end of sigmamode
!! == GW driver ==
      if(llmfgw) then
         nchan = pot_nlma
         call sugw(ssite , sspec ,   nbas , ndham , zv_a_osmpot 
     .     , vconst , lcplxp , sv_p_osig , sv_p_otau , sv_p_oppi , ppnl_rv 
     .     , vrmt , spotx    , sv_p_osigx, sv_p_otaux, sv_p_oppix, jobgw 
     .     , nk1,nk2,nk3,sigmamode,ndimsig, numprocs,procid,master,lso) 
! NOTE: Sigma-Vxc in real space is give by rv_a_ohrs.
         deallocate(spotx)
         call tcx('bndfp')
#if MPIK
         call mpi_barrier(MPI_COMM_WORLD,ierr) 
#endif
         call rx0('sugw mode')  !exit program here normally.
      endif
      if(ham_scaledsigma/=1d0.and.sigmamode) then
        write(6,*)' Scaled Sigma method: ScaledSigma=',ham_scaledsigma
      endif
      
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
      if(pwmode>0.and.pwmode<10) call shortn3_initialize(qlat,rlatp,xmx2)
      allocate(ndimhx_(nkp),nev_(nkp),nevls(nkp,nspx))
      ndimhx_=0
      nev_   =0                
      nevls  =0
      
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
      ix99 = 0
      xdat = 0d0
 99   continue
      ix99 = ix99+1               !ix99=2 for second pass
      ebot = 1000d0
      call surho ( nbas , ssite , sspec , lmet , ldos , lrout , lekkl 
     .  , numq , sv_p_oqkkl , sv_p_oeqkkl , k1 , k2 , k3 , srout_zv , 
     .  ndos , dos_rv , sumev , sumqv )
      if (lfrce .gt. 0) then
         frc(1:3,1:nbas*numq) = 0d0
         fh_rv= 0d0              ! 3 * nbas 
      endif
      if (lswtk .eq. 1) then
         efmax = 1d3
         nevmx = ndhamx
         call dpzero ( rv_a_oswtk , ndhamx * nkp )
      endif
!! --> Mullikan mode was here
!! ... For insulator, valence band top and conduction band bottom
      evtop = -99999
      ecbot = -evtop
      if(lso/=0) orbtm_rv=0d0
      initiqisp=.true.
C Start k and isp loop (first loop in parallel mode)
C     h,s,evecs are dimensioned ndimh in this loop (may be q dependent)
C     evl is dimensioned evl(ndham,2)
C     evlall, w(owtkb), w(oswtk) are dimensioned (ndham,nsp,nkp)
!! === do loop 2010 for iq and do 2005 isp ===
      if(debug) print *,' --- goto do loop 2010 ---'
      do 2010 iq = iqini, iqend !This is a big iq loop
         qp=qplist(:,iq)
         if(debug) print *,' do 2010 iq procid=',iq,procid
#if MPIK
         if(debug) write(*,'(a,i3,5i5)')'debug iq=',procid,iq,iqini,iqend
         if (iq .eq. kpproc(procid)) then
            if(mlog) then
               call gettime(datim)
               aaachar=' bndfp '//datim//' Process '//
     &              trim(i2char(procid))//' of '//trim(i2char(numprocs))//' on '//
     &              shortname(procid)(1:namelen(procid))//' starting k-points '//
     &              trim(i2char( kpproc( procid)))//' to '//trim(i2char( kpproc(procid+1)-1))
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
         if(allocated(igv2x)) deallocate(igv2x)
         if (pwemax .gt. 0 .and. mod(pwmode,10) .gt. 0) then
            pwgmin = dsqrt(pwemin)
            pwgmax = dsqrt(pwemax)
            qqq = 0d0 
            if (mod(pwmode/10,10) == 1) qqq = qp !pwmode 1 in 10th digit means q-dependent nw
            call pshpr(1)
            call gvlst2(alat,plat,qqq,0,0,0,pwgmin,pwgmax,0,0,napw,napw,dum,dum,dum,dum) !get napw
            allocate(igv2x(3,napw), kv_iv(3*napw))
            call gvlst2(alat,plat,qqq,0,0,0,pwgmin,pwgmax,0,2,napw,napw,kv_iv,dum,dum,igv2x)
            call poppr
            if (pwmode<10) then
               do iapw=1,napw
                  igv2x(:,iapw)=igv2x(:,iapw)+nlatout(:,1)
               enddo
            endif
            deallocate(kv_iv)   !we only keep 
            ndimh = ldim + napw !ndimh is iq-dependent for pwmode >=10
            if (mod(pwmode,10) .eq. 2) ndimh = napw !APW-only mode
            if (ndimh .gt. ndham) then
               call rx('bndfp:something wrong ndimh>ndham')
            endif
         else
            ndimh = ldim
            allocate(igv2x(1,1)) !dummy
         endif
         ndimhx = ndimh*nspc    !this is iq-dependent.
         nmx    = min(nevmx,ndimhx) 
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
!! takao removed shorbz(shortest qp) here. !I think not necessary
!! === For this qp, G vectors for PW basis and hamiltonian dimension ===
            ispc = min(isp,nspc) !nspc=2 for non-colinier case
            jsp = isp
            if (ispc .eq. 2) jsp = 1
ccccccccccccccccccccccccccccccc
            if(ix99==2.and.(nspc==2.and.isp==1)) cycle
            if(ix99==2) goto 1919
ccccccccccccccccccccccccccccccc            
!!   ... Make Hamiltonian and overlap matrices
c            nqsig=ham_nqsig
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
!! === Sigma-Vxc ===
!! Generate sene(Sigma-Vxc) for given sfz.
!! Determine interpolated self-energy sene at qp from sfz. 
!! sigmat = Sigma-Vxc is generated in a basis of ndimsig (usually MTOs only) 
!!     ... Bloch transform sigm(RS)-sigm(k). :RS means realspace
!! Main input  => ham_iv_a_oiaxs,ham_rv_a_ohrs
!! Main output => sene
            if(sigmamode) then
               call getsene(lrsig,qp,isp) !Get sene = Sigma-Vxc 
            endif
!!note kino
!! hambl calls augmbl.  augmbl calculates C_kl(=array b in the program)
!! finally makes F~F~=F0F0+(F1F1-F2F2), which is overlap matrix, s.
!! Note that F2=Hankel head at a site + Hankel tail contributions from the other site.

!! == Set up Hamiltonian. ==============
!!    Hamiltonian: hamm(1:ndimh,1:ndimh,ispc). ispc=3 means off-diagonal section when SO=1.
!!    Overlap matrix: ovlm 
!!    senex: Sigma-Vxc
!! ==========================================
            jx = 0
            if (lso == 1) jx = 10
            call hambl ( jx , nbas , ssite , sspec ,  isp ,!sham ,slat , 
     .        qp , k1 , k2 , k3 , zv_a_osmpot , vconst , sv_p_osig , sv_p_otau 
     .        , sv_p_oppi , lcplxp , alfa , ndimh , napw , igv2x , hamm ( 1 
     .           , 1 , ispc ) , ovlm ( 1 , 1 , ispc ) , hamm ( 1 , 1 , 3 ) )
c            print *,'hamble sumcheck =',sum(abs(ovlm)),sum(abs(hamm))
            if(sigmamode) then
               allocate( ovlmtoi(ndimsig,ndimsig),ovliovl(ndimsig,ndimh))
               ovlmtoi = ovlm(1:ndimsig,1:ndimsig,ispc)
               call matcinv(ndimsig,ovlmtoi)
               ovliovl = matmul(ovlmtoi,ovlm(1:ndimsig,1:ndimh,ispc))
               deallocate(ovlmtoi)
               allocate(senex(ndimh,ndimh))
               senex = matmul(transpose(dconjg(ovliovl)), matmul(sene,ovliovl))
               call dsene()
               deallocate(ovliovl)!,sene)
               hamm(:,:,ispc) = hamm(:,:,ispc) + ham_scaledsigma * senex
               deallocate(senex)
            endif
C   ... In the noncollinear case (spin-orbit SO=1), this marks the end of looping over
C       separate spins.
C       For the rest of the isp loop,
C       jsp=isp in the collinear case; jsp=1 in the noncollinear
C       Thus jsp should be used in place of isp
C       isp serves as a flag for the noncollinear case
            if (ispc .ne. nspc) cycle !do2005 goto 30
c            jsp = isp
c            if (ispc .eq. 2) jsp = 1
            
!! === Diagonalize and add to density ===
            if(iprint()>=30.and.ix99==1) 
     &        write(*,'(" bndfp: kpt ",i5," of ",i5, " k=",3f8.4," ndimh = nmto+napw = ",3i5,f13.5)')
     &        iq,nkp,qp,ndimh,ndimh-napw,napw
            nmx = min(nevmx,ndimhx) !lager is safer, but I think nmx can be given in a better manner.
c!!  Calculate all eigenvalues for plotting mode, pdos, tdos modes.
c            if (tdos) then 
c               nmx = ndimhx
c               efmax = 99999
c            endif
            if(debug) write(*,*)'kino, allocate t_zv' ,ndimhx,nmx,ndimhx,ix99
!! == lpdiag=0 only. ===
            if(nspc .eq. 2) then
               call sopert2 ( 0 , ndimh , nspc, hamm , hamm )
               call sopert2 ( 0 , ndimh , nspc,  ovlm , ovlm )
            endif
            if(debug) write(*,*) 'allocate ww_rv', 11*ndimhx
            
!! === diagonalization ===
!! It is better to clean up lwtkb-algolism in this routine. !!!!
!     ! savez save eigenfuncitons for second pass.(second-pass procedure should be removed in future).
ccccccccccccccccccccccccccccccc            
 1919       continue
ccccccccccccccccccccccccccccccc            
!!  Calculate all eigenvalues for plotting mode, pdos, tdos modes.
            if (tdos) then 
               nmx = ndimhx
               efmax = 99999
            endif
            if (allocated(t_zv)) deallocate(t_zv)
            allocate(t_zv(ndimhx,nmx))
            savez=.false.
            getz =.false.
            if(gotosecondpass) then
               if(ix99==1) then
                  savez=.true.  ;  getz=.false.
               else
                  savez=.false. ;  getz=.true.
               endif
               if(initiqisp) then
#if (MPI|MPIK)
           !!  for parallel write, and eigze is temporaly used only in this subroutine. 
                  inquire(file='eigze_'//trim(strprocid),opened=nexist,number=ifig)
                  if(nexist) close(ifig)
                  open(newunit=ifig,file='eigze_'//trim(strprocid),form='unformatted')
#else
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
              if(epsovl<1.000001d-14) then
                call zhev_tk2( ndimhx , hamm , ovlm , nmx , nev,
     .           evl(1 , jsp ), t_zv, iq==1,savez,getz,ifig) !stock eigenfunctions z
              else    !nov2015 recover epsovl. Remove space by overlap matrix diagonalization.
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
            if( procid.eq.master.and. epsovl>=1.000001d-14.and.plbnd/=0) then
              write(6,"('      : ndimhx=',i5,' --> nev=',i5' by HAM_OVEPS ',d11.2)") ndimhx,nev,epsovl
            endif  
#endif
!  ! save to evall
            if(debug) print *,' nmx nev plbnd=',nmx,nev,plbnd
            evlall(1:ndhamx,jsp,iq) = evl(1:ndhamx,jsp)
            
!!   ... Save data for this qp to moments file
            if( plbnd .eq. 0) then
!!   --- Orbital magnetic moment (requires lso) ---
               if (lso .ne. 0 .and. lwtkb .ne. -1) then
                  if (lwtkb .eq. 0) call rx('metal weights required for orb. moment calculation')
                  allocate(auso_zv(nlmax*ndhamx*3*nsp*nbas))
                  if(debug)   write(*,*) 'allocate auso_zv',nlmax*ndhamx*3*nsp*nbas
                  auso_zv(:)=0.0d0
                  call makusq ( 0 , ssite , sspec ,  nbas , nbas , 
     .              0 , nlmax , ndham , ndimh , napw , igv2x , nev , nsp , 
     .              nspc , jsp , 1 , qp , t_zv , ppnl_rv , auso_zv )
                  call mkorbm ( ssite , sspec , jsp , nsp , nspc , nlmax , ndham
     .              , nev , rv_a_owtkb , iq , nbas , ppnl_rv , auso_zv , nl , nkp 
     .              , orbtm_rv )
                  deallocate(auso_zv)
               endif

c$$$!! Mulliken analysis was here
c$$$               if (lwtkb .ne. -1) then
c$$$                  if (mullwrite) then
c$$$                        call rxx(nev.ne.ndimhx,'Mulliken requires all eigenvectors')
c$$$                        call rxx(nevmx.ne.ndhamx,'Mulliken requires nevmx=nhdamx')
c$$$                        allocate(doswt_rv(nchan*nevmx*nspc))
c$$$                        doswt_rv(:)=0.0d0
c$$$                        call mullmf(nbas , ssite , sspec , ham_iv_a_oindxo , t_zv , ndimh 
c$$$     .                    , nspc , iq , jsp , moddos , nsites , lsites , lmxch , nchan 
c$$$     .                    , chan_iv , lmdim , ndham , doswt_rv )
c$$$                        call iomomn(.true.,2,.false.,1,nspc,1,1,i)
c$$$                        i1 = iomoms ( - nfilem , nl , nsp , nspc , nkp , ndimh , i , 
c$$$     .                    nspc , 1 , 1 , ndhamx , nevmx , nchan , nchan , nev , evl ( 1 
c$$$     .                    , jsp ) , 0d0 , doswt_rv , 0d0 , 0d0 , 0d0 )
c$$$                        if (allocated(doswt_rv)) deallocate(doswt_rv)
c$$$                  endif
c$$$               endif   

!! Make new density matrix dmatu for LDA+U (implementation of LDA+U is a little complicated).
               if (lwtkb .ne. -1) then
                  if (nlibu .gt. 0 .and. nev .gt. 0) then
                     if (lwtkb .eq. 0)
     .                 call rx('metal weights required for LDA+U calculation')
c                     nl = int(ctrl_nl)
                     nlmax = nl*nl
                     allocate(ausu_rv(2*nlmax*ndhamx*3*nsp*nbas))
                     if(debug) write(*,*)'allocate ausu_rv',2*nlmax*ndhamx*3*nsp*nbas
                     ausu_rv(:)=0.0d0
                     call makusq ( 0 , ssite , sspec ,   nbas , nbas , !sham ,slat ,
     .                 0 , nlmax , ndham , ndimh , napw , igv2x , nev , nsp , 
     .                 nspc , jsp , 1 , qp , t_zv , ppnl_rv , ausu_rv )
                     call mkdmtu ( ssite , sspec , rv_a_owtkb , jsp , iq , nsp , nspc 
     .                 , ndham , nlmax , nbas , nev , ppnl_rv , ausu_rv , dmatu , nlibu 
     .                 , lmaxu , lldau )
                     if (allocated(ausu_rv)) deallocate(ausu_rv)
                  endif
               endif
!! Core-level spectroscopy
               if(icls .ne. 0 .and. lwtkb .ne. -1) then
                  call rxx(nspc.ne.1,'CLS not implemented in noncoll case')
#if MPI
                  call rx('CLS only k-parallel')
#endif
                  call makusq ( 0 , ssite , sspec ,  nbas , nsites !sham ,slat , 
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
                 if(debug) print *,'goto addrbl nevl=',nevl,lwtkb
                 if(lwtkb>=0) then
                   call addrbl ( ssite , sspec ,  jsp , nsp , nspc !sham ,slat , 
     .              , qp , rv_a_owtkp , ndham , ndimh , napw , igv2x , lmet , lrout
     .              , lwtkb , rv_a_owtkb , lswtk , rv_a_oswtk , iq , lfrce , ldos
     .              , lekkl , k1 , k2 , k3 , zv_a_osmpot , vconst , sv_p_osig , sv_p_otau
     .              , sv_p_oppi , lcplxp , numq , qval - qbg , t_zv , evl , nev
     .              , ef0 , def , esmear , emin , emax , ndos , dos_rv , srout_zv 
     .              , sumqv , sumev , sv_p_oqkkl , sv_p_oeqkkl , frc )
                 endif     
               endif
            endif
            

!!   === Save evals for this qp ===
c            elseif (plbnd .ne. 0) then
c$$$#ifndef MPIK
c$$$               if(debug) write(*,*)  'plbnd/=0 (noMPI)', ndham,jsp 
c$$$               if (nfbn(1)>0 .and. n_listwf>0) then
c$$$                  call rsibl_ev( ssite,sspec,slat,nbas,isp,qp,iq,ndimh,nspc
c$$$     .              ,  napw,igv2x,ham_iv_a_oindxo,nev,t_zv,k1,k2,k3
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
               call makusq( 1 , ssite , sspec ,   nbas , nbas , !sham ,slat ,
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
c 30         continue
c            deallocate(igv2x)
 2005    continue
         if(allocated(hamm)) deallocate(hamm,ovlm)
         ndimhx_(iq)=ndimhx
         nev_(iq)   =nev
 2010 continue
      if(debug) write(stdo,"(' ---- end of do 2010 ---- ',2i5)") procid,ix99
!!----------------------
!! fixed at sep2020. nev_ and ndimhx_ are not correctly passed via MPI      
      if(numprocs/=1.and.ix99==1) then
         do iq = 1, nkp
            pidorigin=-999
            do idx=0,numprocs-1
               if (kpproc(idx)<=iq .and. iq<= kpproc(idx+1)-1) then
                  pidorigin=idx
                  exit
               endif
            enddo   
            call MPI_BCAST(ndimhx_(iq),1,MPI_INTEGER, pidorigin,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(nev_(iq),1,MPI_INTEGER,    pidorigin,MPI_COMM_WORLD,ierr)
         enddo   
      endif
      nevmin=minval(nev_(1:nkp))
!! ========== end of iq loop ====================================================

      
!! === pdos mode (--mkprocar and --fullmesh). ===
      if(debug.and.PROCARon) print *,'mmmm procid sum dwgt check=',procid,sum(dwgtall)
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
            open(newunit=ifip,file='pdosdata.'//trim(sname),form='unformatted')
            write(ifip) ndhamx,nsp,nspx,nevmin,nchanp,nbas,nkk1,nkk2,nkk3,ntete,nkp !ndos,nkp mar2015
            write(ifip) idtete,ipqe ! allocate(idtete(0:4,6*nkp),ipqe(nkk1,nkk2,nkk3))
            write(ifip) evlall  !  allocate(evlall(ndhamx,nspx,nkp))
            write(ifip) dwgtall !  allocate(dwgtall(nchanp,nbas,ndhamx,nspx,nkp))
            write(ifip) ef0 !eminp,emaxp,ef0
            close(ifip)
            if(nkp/=nkk1*nkk2*nkk3) call rx('pdosmode but nkp/=nkk1*nkk2*nkk3')
            if( cmdopt('--tetraw',8,0,strn)) then
               open(newunit=ifip,file='tetradata.dat',form='unformatted')
               write(ifip) ndhamx,nkp,ntete !ndos,nkp mar2015
               write(ifip) idtete
               close(ifip)
               open(newunit=ifip,file='qplistf.dat')
               do iq=1,nkp
                  write(ifip,"(3f15.8)") qplist(:,iq)
               enddo   
               close(ifip)
            endif               !! we need evlall(ndhamx,nkp) for qplist
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
c         call xmpbnd ( kpproc , nlmax * ndhamx * 3 * nsites , nkp , nsp, ausc_zv )
      endif
      call mpibc2(nevls,nkp*nspx,2,mlog,'bndfp','nevls')
#endif
!!
!! boltztrap data 
      if( cmdopt('--boltztrap',11,0,strn) .and. procid==master) then
         iqread=0
         iqindex = index(strn(12:),'nb=')+2
         nbandx = ndhamx
         if(iqindex/=2) then 
            read(strn(12+iqindex:),*) nbandx
            nbandx=min(nbandx,ndhamx)
         endif
         open(newunit=ifi,file='efermi.lmf') !readin fermi energy from efermi.lmf
         read(ifi,*)  eferm
         close(ifi)
        open(newunit=ifi,file=trim(sname)//'.intrans_template.boltztrap')
        write(ifi,"(a)")'GENE          # format '
        write(ifi,"(a)")'0 0 0 0       # iskip (not presently used) idebug setgap shiftgap'
        write(ifi,"(f20.16,a,f9.4,a)") eferm,' 0.0005 0.4 ', qval-qbg,
     &  '       # efermi.lmf (Ry), energygrid, energy span around Fermilevel, number of electrons'
        write(ifi,"(a)")'CALC          # CALC (calculate expansion coeff), NOCALC read from file'
        write(ifi,"(a)")'5             # lpfac, number of latt-points per k-point'
        write(ifi,"(a)")'BOLTZ         # run mode (only BOLTZ is supported)'
        write(ifi,"(a)")'0.15          # (efcut) energy range of chemical potential'
        write(ifi,"(a)")'800.0 50.0  # Tmax, temperature grid'
        write(ifi,"(a)")'-1.0  # energyrange of bands given individual DOS output sig_xxx and dos_xxx (xxx is band number)'
        write(ifi,"(a)")'TETRA'
        close(ifi)
        open(newunit=ifi,file=trim(sname)//'.struct.boltztrap')
        write(ifi,"(a)") trim(sname)
        write(ifi,"(3d16.8, ' # plat1  ')") plat(:,1)*alat
        write(ifi,"(3d16.8, ' # plat2  ')") plat(:,2)*alat
        write(ifi,"(3d16.8, ' # plat3  ')") plat(:,3)*alat
        write(ifi,*) ngrp
        do ig=1,ngrp
           symxx = matmul(transpose(qlat), matmul(symops(:,:,ig),plat))
           if(abs(sum( nint(symxx(:,:))-symxx(:,:) )) >1d-6) call rx('bndfp: sym not integers --- probably strange')
           write(ifi,"(9i3)") ((nint(symxx(i,j)),j=1,3),i=1,3)
        enddo
        close(ifi)   
        do jsp=1,nspx         !=isp.  nspx=1 for so=1
           open(newunit=ifi,file=trim(sname)//'.energy.isp'//trim(i2char(jsp))//'.boltztrap')
           write(ifi,"(a)") trim(sname)
           write(ifi,"(i10)") nkp
        do iq=1,nkp
!! true q= 2pi/alat *qplist(:,iq) in cartesian
           qvec= matmul(transpose(plat),qplist(:,iq)) !qvec in qlat unit
           ncount=0
           do iband=1,nbandx
              if(evlall(iband,jsp,iq)>1d98) cycle
              ncount=ncount+1
           enddo
           write(ifi,"(3f15.8,i10)") qvec,ncount
           do iband=1,nbandx
              if(evlall(iband,jsp,iq)>1d98) cycle !evall=1d99 is dummy
              write(ifi,"(d23.16)") evlall(iband,jsp,iq)
           enddo   
        enddo
           close(ifi)
        enddo
      endif
      if( cmdopt('--boltztrap',11,0,strn)) then
         call rx0(' end of --boltztrap mode: boltztrap.* are generated')
      endif
      
      
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
c         nkaph = globalvariables%nkaph
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
      eterms=ham_eterms ! call dcopy(size(ham_eterms),ham_eterms,1,eterms,1)
      call mpibc2(eterms(19),1,4,mlog,'bndfp','rhosig')
      ham_eterms=eterms ! call dcopy(size(ham_eterms),eterms,1,ham_eterms,1)
#ifdef MPIK      
      entime = MPI_WTIME()
      call info2(20,0,0,' MPI broadcast took %;1d sec', (entime-sttime),0)
#endif


!! ===  Write bands in bands-plotting case: loop over qp getting evals from array ===
      if(plbnd .ne. 0) then
         if (procid .eq. master) then
            print *,'procid=',procid
!! fermi surface mode together with band data around Ef (Okumura, 2019/11/14)
            if(fsmode) call writefs()
!! write bands directly. All inputs
            write(stdo,*)' Writing bands to bands file ...'
            if(nsyml/=0) then
              call writeband(stdo,nkp,nsyml,nspx,nsp,ix99,ndhamx, !,fsmom/=NULLR,vnow
     i        xdatt,nqp_syml,nqp2n_syml,nevls,evlall,qplist,labeli,labele,alat,eferm,evtop,ecbot,
     i        sname,nqps_syml,nqpe_syml,nsyml,dqsyml,
     i        etolv,etolc)
            endif
         endif
         if(fsmom/=NULLR) then
         write(6,"(a)")'NOTE: Be carful! fsmom (fixed moment on). In sc cycle, we use additional bias mag. field  '
         write(6,"(a)")'NOTE: Use mag.field is written in MagField. But it is not used for --band mode! mar2015.'
         endif
         if(fsmode) call rx0('done --fermisurface mode. *.bxsf for xcryden generated')
         call rx0('plot band mode done')
      endif
      
!! ===   Repeat loop for printout and goto 99 ===
!!   jsp=isp in the collinear case; jsp=1 in the noncollinear
!!     Thus jsp should be used in place of isp
!!     isp serves as a flag for the noncollinear case
      goto99=.false.
      do  iq = 1, nkp
         qp=qplist(:,iq)
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
            nev_iq    = nev_(iq) 
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
               if (lwtkb .ne. -1) then ! .and. .not. lwndow) then
                  if (iq .eq. 1 .and. jsp .eq. nsp ) then !
!! Energy-Window reset. dosw=(emin,emax)
                     ef00 = ef0
                     call fixef0(qval-qbg,jsp,1,nev_iq,ndhamx,evl,dosw,ef0)
                     if(procid==master) then
                        ifid= fopn ( 'ewindow' )
                        write(ifid,"(3d23.15)") dosw(1:2),ef0
                        call fclose(ifid)
                     endif
                     if (jsp .eq. 2 .and. ef00 .ne. ef0 .and.
     .                    lwtkb .eq. 0 .and. lmet .gt. 0 .and. lrout .ne. 0) then
                        if (procid .eq. master) call info0(10,1,1,
     .                       ' ... Fermi level reset in second spin'//
     .                       ' channel ... restart band pass')
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
c      if (allocated(kpproc)) deallocate(kpproc, stat=ierr)
      call mpi_barrier(MPI_comm_world,ierr)  !!! takao added I think it must be safer to add this barrier. jan2015
#endif
      print *,'goto99=',goto99
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
         bz_ndos=ndos
         bz_dosw = dosw
         bz_ef=ef0
         bz_def=def
         if (lrep .eq. 1) then
            ef0 = -1
           write(stdo,*)'Input Fermi energy was too far off, repeat band pass'
#ifdef MPIK
            call mpi_barrier(MPI_comm_world,ierr)  
#endif
            goto 99
         endif
      endif
#if MPIK
      call mpi_barrier(MPI_comm_world,ierr)
#endif
      if (procid .eq. master) then
         ierr=unlink('MagField') !delete
      endif  
!! --- BZ integration for fermi level, band sum and qp weights ---
      if(lmet .ge. 0 .and. (lmet .ne. 4 .or. ltet)) then
c$$$         if (lwndow) then
c$$$            allocate(ww_rv(ndham*nsp*nkp))
c$$$            eferm = min(dosw(1),dosw(2))
c$$$            call bzints ( nkabc ( 1 ) , nkabc ( 2 ) , nkabc ( 3 ) , evlall 
c$$$     .        , ww_rv , nkp , ndham , ndham , nsp , xv , xv , xv , 1 , eferm 
c$$$     .        , 2 , ntet , iv_a_oidtet , sev , dum )
c$$$            eferm = max(dosw(1),dosw(2))
c$$$            call bzints ( nkabc ( 1 ) , nkabc ( 2 ) , nkabc ( 3 ) , evlall 
c$$$     .        , rv_a_owtkb , nkp , ndham , ndham , nsp , xv , xv , xv , 1 , 
c$$$     .        eferm , 2 , ntet , iv_a_oidtet , sev , dum )
c$$$            call daxpy ( ndham * nsp * nkp , - 1d0 , ww_rv , 1 , rv_a_owtkb, 1 )
c$$$            if (allocated(ww_rv)) deallocate(ww_rv)
c$$$         else
         dosrng = 8
         if (mpsord .lt. 0) dosrng = 16
         if( debug) print *, 'tttttttttt 99999999 6 call bzwtsf'
         if(bz_fsmommethod == 1) then !takao dec2010
!     ! vnow june22013 !vnow !(in Ry) contains magnetic field
!     ! For eigenvalus, add  -vnow/2 for isp=1, and +vnow/2 for isp=2.
            call bzwtsf2 ( ndham , ndham , nsp , nspc , nkabc ( 1 ) , nkabc 
     .           ( 2 ) , nkabc ( 3 ) , nkp , ntet , iv_a_oidtet , qval - qbg , 
     .           fsmom , lmet.ne.0 , ltet , mpsord , ndos , dabs ( esmear - mpsord 
     .           ) , dosrng , rv_a_owtkp , evlall , efmax , lswtk , rv_a_oswtk 
     .           , eferm , sev , rv_a_owtkb , sumqv ( 1 , 2 ) , lwtkb ,lfill,vnow)
         else
!     ! vnow june22013
            call bzwtsf ( ndham , ndham , nsp , nspc , nkabc ( 1 ) , nkabc 
     .           ( 2 ) , nkabc ( 3 ) , nkp , ntet , iv_a_oidtet , qval - qbg , 
     .           fsmom , lmet.ne.0 , ltet , mpsord , ndos , dabs ( esmear - mpsord 
     .           ) , dosrng , rv_a_owtkp , evlall , efmax , lswtk , rv_a_oswtk 
     .           , eferm , sev , rv_a_owtkb , sumqv ( 1 , 2 ) , lwtkb ,lfill, vnow)
         endif
!     ! june2013 magfield is added
         if(fsmom/=NULLR.and.procid .eq. master) then
            ifimag=ifile_handle()
            open(ifimag,file='MagField',status='unknown')
            write(ifimag,"(d23.16,' !(in Ry) -vnow/2 for isp=1, +vnow/2 for isp=2')")vnow
            close(ifimag)
         endif
!     !         Store val charge & magnetic moment in sumqv(1..2)
         if (lmet .ne. 4) then
            sumqv(1,1) = sumqv(1,2)
            sumqv(2,1) = sumqv(2,2)
         endif
c$$$         endif
!!
         if (lmet .ne. 4) then
            ef0 = eferm
            bz_ef=ef0
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
c$$$  call poppr
         if(goto99) goto 99
      endif
c      if(allocated(sfz)) deallocate(sfz)
!! == ix99 loop over =============================

!! ... Save Fermi level, nonmetal or sampling integration
      if (lmet .eq. 0 .or. .not. (lmet .ne. 4 .or. ltet)) then
         if (lmet .eq. 0) then
            ef0 = (evtop+ecbot)/2
            write(6,"(' Highest occ. level = ',f13.5,
     &        ' Lowest unocc. = ',f13.5, ' average =',f13.5)") evtop,ecbot,ef0
         endif
         if (procid .eq. master) then
            ifi = fopna('wkp',-1,4)
            i = iobzwt ( 1 , ndham , nkp , nsp , ef0 , rv_a_owtkb , -ifi )
            call fclr('wkp',ifi)
         endif
      endif
c$$$ -->!! Here ... Cleanup asa-style moments file, print table of DOS channels

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
c      if(.not. lwndow .and. procid==master) then
      if( procid==master) then
         ndos0 = bz_ndos !ndos0 is only used as a swith wether DOS or integrated-DOS.
         if ( ldos/=0 ) then
            if(ndos0<0) write(6,"(a)")'... Generating integrated DOS' !, ndos0=,ndos0
            if(ndos0>0) write(6,"(a)")'... Generating total DOS'      !  ndos0=,ndos0
!! Default and external option. 
            if(tdos) then
              dosmax = bz_dosmax
            endif
            if(cmdopt('--tdostetf',10,0,strn)) ltet=.false. !oct2015
            if (ltet) then
!! determine eferm
               call bzints ( nkabc ( 1 ) , nkabc ( 2 ) , nkabc ( 3 ) , evlall 
     .           , dum , nkp , ndhamx , ndhamx , nspx , dosw ( 1 ) , dosw (
     .           2 ) , dos_rv , ndos , eferm , 1 , ntet , iv_a_oidtet , dum , dum )
!! repeat again for given range emax=dosw(2). mar2012
               if(tdos) dosw(2) = dosmax + eferm !.or.pdos (dos between dosw(1) and dosw(2)
               call bzints ( nkabc ( 1 ) , nkabc ( 2 ) , nkabc ( 3 ) , evlall 
     .           , dum , nkp , ndhamx , ndhamx , nspx , dosw ( 1 ) , dosw (
     .           2 ) , dos_rv , ndos , eferm , 1 , ntet , iv_a_oidtet , dum , dum )!job=1 integrated dos
               if ( ndos0 > 0 ) call xxxdif ( dosw ( 1 ) , dosw ( 2 ) , ndos 
     .           , nspx , 0 , dos_rv ) !derivative 'integrated dos -> DOS'
               del = 0d0
            else
               if(tdos) dosw(2) = dosmax + eferm  !.or.pdos
               if ( mpsord .ge. 100) mpsord = mod(mpsord,100)
               if ( ndos0 > 0 ) call makdos ( nkp , ndhamx , ndhamx , nspx
     .           , rv_a_owtkp , evlall , mpsord , dabs ( esmear - mpsord ) 
     .           , - 6d0 , dosw ( 1 ) , dosw ( 2 ) , ndos , dos_rv )
               if ( ndos0 < 0 ) call maknos ( nkp , ndhamx , ndhamx , nspx
     .           , rv_a_owtkp , evlall , mpsord , dabs ( esmear - mpsord ) 
     .           , - 6d0 , dosw ( 1 ) , dosw ( 2 ) , ndos , dos_rv )
               del = mpsord+dabs(esmear-mpsord)
            endif
            if ( nspc .eq. 2 ) call dscal ( ndos , .5d0 , dos_rv , 1 )
            open(newunit=ifi,file='dos.tot.'//trim(sname))
            open(newunit=ifii,file='dosi.tot.'//trim(sname))
            dee=(dosw(2)-dosw(1))/(ndos-1d0)
            do ipts=1,ndos
               eee=dosw(1)+ (ipts-1d0)*(dosw(2)-dosw(1))/(ndos-1d0)-eferm
               dosi(1:nsp)= dosi(1:nsp) + dos_rv(ipts,isp)*dee
               write(ifi,"(255(f13.5,x))") eee,  (dos_rv(ipts,isp),isp=1,nsp)
               write(ifii,"(255(f13.5,x))") eee, (dosi(isp),isp=1,nsp)
            enddo
            close(ifi)
            close(ifii)
c            if (allocated(orbtm_rv)) deallocate(orbtm_rv)
c            if (allocated(dos_rv)) deallocate(dos_rv)
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
            eferm=bz_ef
            call vcdmel ( nl , ssite , sspec ,  nlmax , ndham , !slat ,
     .        ndimh , nkp , nsp , nspc , eferm , evlall , ausc_zv , 
     .        nsites , isite , iclsl , iclsn ,dosw) !sbz,
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
c$$$!!     --window: Put output density into orhat and smrho, and exit
c$$$      if (lwndow) then
c$$$         call mkrout ( ssite , sspec ,  nbas , nsp , ldim  !sham , slat ,
c$$$     .     , lekkl , sv_p_oqkkl , sv_p_oeqkkl , sv_a_oorhat , hab_rv , sab_rv 
c$$$     .     , qbyl_rv , hbyl_rv , lrout )
c$$$         call zcopy( k1 * k2 * k3 * nsp , srout_zv , 1 , zv_a_osmrho 
c$$$     .     , 1 )
c$$$         call symrho(ssite , sspec ,  lfrce , zv_a_osmrho , sv_a_oorhat !slat ,
c$$$     .     , qbyl_rv , hbyl_rv , frc )
c$$$         goto 999
c$$$      endif

      call mkrout ( ssite , sspec ,  nbas , nsp , ldim  !sham ,slat , 
     .  , lekkl , sv_p_oqkkl , sv_p_oeqkkl , sv_p_orhat1 , hab_rv , sab_rv 
     .  , qbyl_rv , hbyl_rv , lrout )
C
      if (lrout .ne. 0) then
C   ... Symmetrize output density and forces
         call symrho ( ssite , sspec ,  lfrce , srout_zv , sv_p_orhat1 !slat ,
     .     , qbyl_rv , hbyl_rv , frc )
C   ... New boundary conditions pnu for phi and phidot
C       call pshpr(iprint()-12)
         if (lpnu .gt. 0) then
            call dcopy(size(ham_pmin),ham_pmin,1,xv,1)
            call dcopy(size(ham_pmax),ham_pmax,1,xv(11),1)
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
      call mkehkf(1,sev,sumqv(2,1),sumtv,ehar) !sham,
!! --- Evaluate KS total energy, correction to Harris force ---
      if (lrout .ne. 0) then
!!   ... Correction to harris force
        job=int(ctrl_lfrce)
        if(job>0 ) then
          call dfrce (job, ssite , sspec ,  k1 , k2 , k3 , nvl !sctrl , slat ,
     .     , sv_a_oorhat , sv_p_orhat1 , elind , qmom_rv , zv_a_osmrho , 
     .     srout_zv , fh_rv )
        endif
!!   ... Evaluate KS total energy and output magnetic moment
         eks = 0d0
         if (leks .ge. 1) then
            call togpr()
            call mkekin ( nbas , ldim , ssite , sspec ,  sv_p_osig  !sham ,slat , 
     .        , sv_p_otau , sv_p_oppi , lcplxp , sv_p_oqkkl , k1 , k2 , k3
     .        , vconst , zv_a_osmpot , srout_zv , sev , sumtv )
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
            call mkpot ( nbas , ssite , sspec ,   lfrce  !sham ,slat ,spot , 
     .        , lcplxp, k1, k2, k3, srout_zv , sv_p_orhat1 , qbg , zv_a_osmpot 
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
            call mkehkf(2,sev,sumqv(2,1),sumtv,eks) !sham,
            call togpr()
         endif
         
!! --- Add together force terms ---
!! fes1_rv: contribution to HF forces from estat + xc potential   This is for input  density  !=3rd term in (B.5) in JPSJ.84.034705
!! fes2_rv: contribution to KS forces from estat + xc potential   This is for output density
!! fh_rv      : 1st term  in (B.5)
!! input frc  : 2nd term in (B.5). Puley term
!! output frc  : totlal
         if(lfrce .gt. 0) call totfrc( nbas , ssite ,  leks !slat ,
     .     , fes1_rv , fes2_rv , fh_rv , frc )
!!   --- Mix input and output densities ---
         allocate(qm1_rv(nvl))
         allocate(qm2_rv(nvl))
         call mixrho ( ssite , sspec ,  nsp , iter , sstrnmix 
     .     , qval - qbg , elind , sv_p_orhat1 , sv_a_oorhat , k1 , k2
     .     , k3 , dmxp , srout_zv , zv_a_osmrho )
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
      ham_ehf=ehar
      ham_ehk=eks
 999  continue
 8888 continue
      if(procid==master) call vbmmode()
      if(cmdopt('--vbmonly',9,0,strn)) call rx0('end of --vbmonly mode vrefolny')      
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
      if (allocated(ifbls_iv)) deallocate(ifbls_iv)
      deallocate(evlall)
 9995 continue
      call tcx('bndfp')
      return
      
!sssssssssssssssssssssssssssssssssssssssssssss
      contains
      subroutine vbmmode()!no output
      implicit none
      integer:: ifvesatm,ispec,ifves,ifvesintloc
      real(8):: vol,vessm,sumvesloc,sumvesatm,vref,vesloc,evbm,ecbm
      real(8),allocatable::vesatm(:)
!     ! === VBM CBM section. Mar2013takao === just output
      print *,'--vbmonly mode'
      vol = lat_vol
      open(newunit=ifvesintatm,file='vesintatm.'//trim(sname),status='old',err=9898)
!! vesintatm is given by lmfa. electrostatic potential integrals.
      allocate(vesatm(ctrl_nspec))
      do ispec=1,ctrl_nspec
         if(sspec(ispec)%z==0d0 .and. sspec(ispec)%rmt==0d0 ) cycle
         read(ifvesintatm,*,err=9898,end=9898) vesatm(ispec)
      enddo
      close(ifvesintatm)
      open(newunit=ifves,file='vessm.'//trim(sname),status='old')
      read(ifves,*) vessm       !smooth part of electrostatic pot. given in smves.F. (es on MT average).
      close(ifves)
c     ifvesintloc=fopna('vesintloc',-1,1) !generated in mkpot-locpot.F electrostatic pot. integral.
      open(newunit=ifvesintloc,file='vesintloc.'//trim(sname),status='old',err=9898)
c     rewind ifvesintloc
      sumvesloc=0d0
      sumvesatm=0d0
      do ib=1,nbas
         ispec=ssite(ib)%spec
         if(sspec(ispec)%z==0d0 .and. sspec(ispec)%rmt==0d0 )  cycle
         read(ifvesintloc,*) vesloc
         sumvesloc = sumvesloc + vesloc
         sumvesatm = sumvesatm + vesatm(ssite(ib)%spec)
      enddo   
      close(ifvesintloc)
      deallocate(vesatm)
      write(*,*)
c     print *,'vessm=',vessm
      vref= sumvesatm/vol - vessm -sumvesloc/vol
      write(*,"('### VBM: Add vref to eigval to estimate eV relative to vaccum. vref(eV)=',f12.6)") vref*Rydberg
      write(*,"('### VBM: Mean estatic pot by superposition of atoms(eV)=',f12.6)")sumvesatm/vol*rydberg
      if(cmdopt('--vbmonly',9,0,strn)) then
         open(newunit=ifvbm,file='vbmcbm.'//trim(sname))
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
      open(newunit=ifvbm,file='vbmcbm.'//trim(sname))
      write(*,    "('### VBM: VBM= ',f12.6,' eV CBM= ',f12.6,' eV estatic_pot=',f12.6,' eV')")  
     &     evbm, ecbm, sumvesatm/vol*rydberg
      write(ifvbm,"('### VBM: VBM= ',f12.6,' eV CBM= ',f12.6,' eV estatic_pot=',f12.6,' eV')")  
     &     evbm, ecbm, sumvesatm/vol*rydberg
      close(ifvbm)
 9898 continue
      print *,' if vesintatm.'//trim(sname)//' ---> VBM: is not shown. Need to repeat new lmfa for VBM:'
      end subroutine vbmmode
!sssssssssssssssssssssssssssssssssssssssssssss
      subroutine writefs() !no output fermi surface mode (eigenvalues in full BZ). Not return variables.
      implicit none
      integer:: i,isp,ififm,ifile_handle,nbxx
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
     &           maxval(evlall(ib,isp,:))>ef0-0.5)) then
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
     &           maxval(evlall(ib,isp,:))>ef0-0.5)) then
               write(ifi,"(a,i8)")"  BAND: ",ib    
               write(ififm,"(a,i8)")"  # BAND: ",ib    
               write(ifi,"(3x,10(x,f9.5))")
     &              (evlall(ib,isp,iq), iq=1,nkk1*nkk2*nkk3)
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
      end subroutine writefs
!ssssssssssssssssssssssssssssssssssssssssssssssssssss      
      subroutine getqmode()  !no output. getq mode just output. Not return variables.
      implicit none
      integer:: lmxa,nspec,is,ifiwv
      real(8),allocatable::qrmt(:,:),ql(:,:),qset(:,:),qatot(:)
         print *,'get Q from ql given by lmf'
         print *,'WARN current version is only for spin symmetric; and not more than 2(2l+1) occupancy'
         nspec=ctrl_nspec
         lmxa = maxval(sspec(1:nspec)%lmxa)
         allocate(qrmt(0:lmxa,nspec),qset(0:lmxa,nspec),qatot(nspec))
         qrmt= -1d-10
         open(newunit=ifiwv,file='veswavatm.'//trim(sname))
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
         open(newunit=ifiqb,file='qbyl.'//trim(sname))
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
c      endif
      end subroutine getqmode
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
