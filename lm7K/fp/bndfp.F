#define USEOPTIONAL 1
C oevl -> evlall
      subroutine bndfp(nbas,nsp,nlibu,lmaxu,lldau,
     .ssite,sspec,slat,sctrl,sham,spot,sbz,sstrn,
     .ndham,leks,lrout,lfrce,lpnu,dmxp,iter,maxit, !evl,
     . frc,
     .dmatu,vorb,llmfgw)
      use m_hamindex, only: 
     &  ngrp, norbt, napwmx, ltab,ktab,offl,ispec,iclasst,invgx,miat,offlrev,ibastab,
     &  igv2,napwk,igv2rev
      use m_rdctrl, only: ncutovl
      use m_struc_def
      use m_globalvariables
      use m_phmbls
      use m_seneinput,only: npair,nlat,qwgt,npairmx,qsmesh !these are given here.
#if defined(__INTEL_COMPILER)
      use ifport, only: unlink
#endif
#if MPI | MPIK
      use m_events !Event numbers: include "events.ins"
#endif
      use m_subzi 
      use m_rsibl_ev, only : rsibl_ev  ! to plot wavefunction in the fat band mode 
      use m_w_psir , only: w_psir_cmdopt
C- One band pass, full-potential hamiltonian
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   nsp   : number of spins
Ci   nlibu : total number of LDA+U blocks (used to dimension dmatu and vorb)
Ci   lmaxu : lmax for U used to dimension vorb and dmatu
Ci   lldau :lldau(ib)=0 => no U on this site otherwise
Ci         :U on site ib with dmat in dmats(*,lldau(ib))
Ci
Ci   ssite :struct for site-specific information; see routine usite
Ci     Elts read: *
Ci     Stored:    *
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
Cl   lwsig :special modes to handling reading/writing of sigma or evecs
Cl         :1  Rotates sigm to LDA basis; saves in file sigm2
Cl         :-1 reads sigm (assumed to be written in the LDA basis),
Cl         :   rotates it to orbital basis, stores the result in file sigm2
Cl         :   lwsig=-1 is the inverse operation of lwsig=1.
Cl         :2  Similar to lwsig=1, except
Cl         :   low- and high- energy blocks replaced by diagonal parts
Cl         :3  Writes evals,evecs of LDA hamiltonian to file,
Cl         :   saves in file 'evec'
Cl         :4  Writes evals,evecs of hamiltonian to file,
Cl         :   saves in file 'evec'
Cl         :5  Writes sigm as in lwsig=2 with low- and high- energy blocks
Cl         :   replaced by diagonal parts, but writes sigm in orbital bas.
Cl   nsmidb:smallest value of nmax encountered in truncating sigma
Cl         :   (only used for printout)
Cr Remarks
Cr   Band pass consists of:
Cr   (1) make the effective potential,
Cr   (2) generate eigenvalues (and eigenvectors if lrout)
Cr   (3) if lrout, assemble the output density by BZ integration
Cr   (4) evaluate hf (and KS, if leks) energy by BZ integration
Cr   (5) mix the output density to make a new input density.
Cu Updates
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
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      logical :: llmfgw
#if MPI | MPIK
      include "mpif.h"
#if MPE
      include "mpef.h"
#endif
      integer numprocs, ierr, status(MPI_STATUS_SIZE)
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*20 ext
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
      double precision starttime, endtime
#if MPI
      integer dims(2)
#endif
#if MPIK
ckino 05Feb09
      integer:: pidorigin, ndimhx_iq, nev_iq
      integer:: work_mpi(10)  ! use only 1 and 2 now
#endif
c#if MPE
cC Event numbers:
c      include "events.ins"
c#endif
#endif
      logical mlog
      integer procid,master
      integer nbx,nbas,n0,nab,nppn,iter,maxit,ndham
      character*(*) sstrn
      parameter ( nbx=512, n0=10, nppn=12, nab=9)
      integer k1,k2,k3,leks,lrout,lfrce,lpnu
      integer nlibu,lmaxu,lldau(nbas)
      real(8):: dmxp(33),frc(3,nbas)
ctakao
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)

c     .  spot(1),sbz(1),sham(1),slat(1),sctrl(1)

c      double precision ssite(1),sspec(1)
      type(s_ham):: sham
      type(s_ctrl):: sctrl
      type(s_lat)::  slat
      type(s_pot)::  spot
      type(s_bz)  ::  sbz
c      type(s_site):: ssitev(nbas)
CKI      type(s_spec),allocatable:: sspecv(:)


C ... Local parameters
      logical lgors,ltet,cmdopt,lwndow
      character strn*120,plbopt*120,mulopt*120,clsopt*120,strn2*120,dc*1
c      real(8) ,allocatable :: eomin_rv(:)

ckino lswtk must be 0 initially.
      integer:: ipr , iprint , ipl , isp , jsp , ispc , isum , nglob 
     ., i , iopq , iq , isqp , ismidb , nsmidb , ldos , lmet , nspc 
     ., lrep , lcplxp , lso , nkabc(3) , ndos , ndos0 , nev , nevl 
     ., nkp , nsp , numq , n1 , n2 , n3 , ntet , lwtkb , lswtk=0 , nl 
     ., lgunit , igets , mpsord , nevmx , nvl , ngabc(3) , isw , stdo 
     ., stdl , ifi , fxst , fopn , fopna , lfrzw , i1 , i2 , plbnd 
     ., nfilem , iobzwt , lnoxc , lrsig , lwsig , jobgw , nll , lpdiag 
     ., iv(10) , parg , ndhamx , nspx , nk1 , nk2 , nk3 , lshft(3) 
     ., onesp , nfbn(2) , mpipid , ldham(8,2) , ldim , ndimh , ndimhx 
     ., lekkl 
      real(8) ,allocatable :: dos_rv(:)

              integer:: os , osoptc 
ckino change rv_p_owtkb allocatable
ckino              real(8),pointer :: rv_p_owtkb(:,:,:)=>NULL()
ckino Dec.8.2011:                real(8),allocatable :: rv_p_owtkb(:,:,:)
              real(8),allocatable :: rv_a_owtkb(:,:,:)

ckino change rv_p_oswtk allocatable             real(8),pointer :: rv_p_oswtk(:)=>NULL()
ckino Dec.8.2011:               real(8),allocatable :: rv_p_oswtk(:)
             real(8),allocatable :: rv_a_oswtk(:)

ckino Dec.16.2011:              integer,pointer :: iv_p_oidtet(:)=>NULL()

ckino Dec.8.2011:             type(s_rv1),pointer :: sv_p_oorhat(:)=>NULL()

ckino Dec.8.2011:            complex(8),pointer :: zv_p_osmrho(:)=>NULL()

ckino Dec.8.2011:           complex(8),pointer :: zv_p_osmpot(:)=>NULL()

ckino Dec.8.2011:          real(8),pointer :: rv_p_owtkp(:)=>NULL()

ckino Dec.8.2011:         real(8),pointer :: rv_p_oqp(:)=>NULL()

      real(8) ,allocatable :: fh_rv(:)
      real(8) ,allocatable :: fes1_rv(:)
      real(8) ,allocatable :: fes2_rv(:)
      real(8) ,allocatable :: gpot0_rv(:)
      real(8) ,allocatable :: hab_rv(:)
      real(8) ,allocatable :: hbyl_rv(:)
      real(8) ,allocatable :: ppnl_rv(:)
      real(8) ,allocatable :: qbyl_rv(:)
      real(8) ,allocatable :: qmom_rv(:)
      real(8) ,allocatable :: sab_rv(:)
      complex(8) ,allocatable :: srout_zv(:)
      complex(8) ,allocatable :: t_zv(:)
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

      equivalence (ldim,ldham(1,1))

      logical a2bin
         integer:: ospotx , osgw 
         type(s_rv1) :: sv_p_oppix(3,nbx)

        type(s_rv1) :: sv_p_otaux(3,nbx)

       type(s_rv1) :: sv_p_osigx(3,nbx)


ckino, used for dummy allocations 
ckino Dec.8.2011:         logical:: l_temp_rv_p_oswtk=.false., l_temp_rv_p_owtkb=.false.,
ckino Dec.8.2011:         logical:: l_temp_rv_a_oswtk=.false., l_temp_rv_p_owtkb=.false.,
       logical:: l_temp_rv_a_oswtk=.false., l_temp_rv_a_owtkb=.false.,
ckino Dec.28.2011:       .    l_temp_iv_p_oidtet=.false.
     .    l_temp_iv_a_oidtet=.false.

ckino end 
      double precision sigp(10)

      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      !equivalence (nk1,nkabc(1)), (nk2,nkabc(2)), (nk3,nkabc(3))
      double precision ef00,eferm,emax,emin,esmear,qval,qsc,sev,
     .sev00,sev1,sumtv,alfa,vconst,qbg,ebot,pi,def,ef0,ehar,eks,
     .dosw(2),efmax,dgets,alat,plat(3,3),qlat(3,3),sumev(2,3),
     .sumqv(3,2),qp(3),shftqp(3),eterms(20),elind,vol,dum,del,evtop,
     .ecbot,vrmt(nbx),fsmom,rsrnge,xv(20),dosrng,epsovl,dval !,eomin
      equivalence (emin,dosw(1)),(emax,dosw(2))

C     For now: PW basis
      integer,allocatable:: igv2x(:,:)
      integer ,allocatable :: kv_iv(:)

      integer pwmode,napw
      double precision pwemin,pwemax,pwgmin,pwgmax

C     real(8),allocatable :: qpe(:,:)

      logical T, F
C     for spin-orbit
      complex(8) ,allocatable :: auso_zv(:)

C     for partial dos
      complex(8) ,allocatable :: ausp_zv(:)

C     for CLS
      integer nsitmx
      parameter (nsitmx = 256)
      integer:: nlmax , i1mach , icls , isite(nsitmx) , iclsl(nsitmx) 
     ., iclsn(nsitmx) 
      complex(8) ,allocatable :: ausc_zv(:)

C     for Mulliken
       integer:: moddos , nsites , lsites(nbx) , nchan , ng , iomoms 
     . , nchmx , lmdim , lmxch 
ckino Dec.8.2011:         integer,pointer :: iv_p_oiprmb(:) =>NULL()

      real(8) ,allocatable :: doswt_rv(:)
      integer ,allocatable :: chan_iv(:),chan0_iv(:)

C     for pzhev
      integer nblk,nprow,npcol,nmx
C     for LDA+U
      real(8) ,allocatable :: ausu_rv(:)

      double complex vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu),
     .dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)

      integer:: i_copy_size, i_spackv, i_spacks
      integer ::iwdummy 


      logical :: debug=.true. !debug is only for check write. However, it can cause error
      !for 'make check',
      !because fp/test/test.fp assumes 'something not inserted in console output'.
#if MPIK
      integer obuf,oqptr
      integer igetss,lmxa,lmxh,nlma,nlmh,nelt(3),nkaph,kmax
      integer, dimension(:),allocatable :: kpproc
      double precision sttime,entime
#endif

C      for debugging and testing
C      integer ocg,ojcg,oidxcg,ocy
C      call upack('lat ocg ojcg oidxcg ocy',slat,ocg,ojcg,oidxcg,ocy,0)
C      call thxpbl(slat,w(ocg),w(oidxcg),w(ojcg),w(ocy))
C #if TESTRWF
C       integer ooffH,oidxsh,oistab,oag,og
C       integer nk1,nk2,nk3,ogstar,is(3),ifac(3),nsgrp
C       double precision rb(9),qb(9)
C       logical llshft(3)
C       integer oipq
C #endif

C For self-energy
      integer nqsig
ckino Dec.8.2011:         real(8),pointer :: rv_p_oqsig(:) =>NULL()


       integer:: nbasp 
ckino Dec.8.2011:         integer,pointer :: iv_p_oidxsh(:) =>NULL()


C     Data statements
      data T /.true./ F /.false./

      integer,allocatable:: nevls(:,:)

ctakao
      integer :: ifiogw,ifiese,oag,og,ib,ig,iorb,is,nsgrp,uspec_size,nspec
      logical :: nexist,sigmamode,ok, newsigmasw
      character spid*8
      integer :: ib1,ib2,ii,ik1,ik2,ik3,j1,j2
       integer:: nhmax , nknknk , jx , nsp_dummy , ndimh_dummy 
ckino Dec.8.2011:         real(8),pointer :: rv_p_opos(:) =>NULL()

      complex(8),allocatable:: sfz(:,:,:,:,:,:),wgtq(:,:,:,:,:),sfzx(:,:,:,:,:,:)
      complex(8),allocatable:: hamm(:,:,:),ovlm(:,:,:),sene(:,:) !,hamsom(:,:)
      integer:: isize_hamm(3)
      real(8):: rrr(3)
      complex(8):: pi2i
      real(8),allocatable:: range(:,:),pos(:,:)


      integer:: j1x,j2x,ix1,ix2,ix3 !ikc(3,0:1,0:1,0:1), igshift(3,0:1,0:1,0:1),
      integer,allocatable::igx(:,:,:,:)
      real(8):: wgc(0:1,0:1,0:1),qpara(3)

c      complex(8)::sss
      real(8):: www,qmm(3),qqq(3),sumtvxxx,sevxxx
      integer:: iktm,ikt,igvecc(3),nkk(3),igvecc0(3), nnspc,osmpotbk
      logical:: noshorbz,siginit,bzwtsf_bisection
      real(8),allocatable:: qplist(:,:)
      real(8):: qpo(3),qpx(3)

      integer:: iqini,iqend
      complex(8),allocatable:: spotx(:,:,:,:),smrhobk(:)
      real(8),allocatable:: evlall(:,:,:),evl(:,:)
      
      integer:: nevec,idummy
      real(8):: epsnevec
      logical:: onlyonce=.true.

      logical:: isanrg, l_dummy_isanrg,oncewrite

      integer:: getef, i_dummy_getef


ckino
      integer:: n_listwf
      integer, allocatable:: listwf(:)

      integer(4):: ntqxx,ispx
      real(8)::eseavrmean

#if !defined(__INTEL_COMPILER)
      integer:: unlink 
#endif


ckino begin

      nevec=0
C     call wkfast(.false.)
      if (maxit .eq. 0) return

C ... MPI setup
      procid = mpipid(1)
      master = 0
      mlog = cmdopt('--mlog',6,0,strn)
#if MPI | MPIK
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
      call strcop(shortname(procid),name,10,'.',i)
      namelen(procid) = i-1
#endif



!!== Setup ==
      call tcn ('bndfp')
      siginit=.true.
      napw = 0
      ipr  = iprint()
      ipl  = ipr
Changenglob      nsp  = nglob('nsp')
      nsp  = globalvariables%nsp
Changenglob      nspc = nglob('nspc')
      nspc = globalvariables%nspc
C     npqn = nglob('npqn')
      lso = isw ( iand(4,int(sctrl%lncol)) .ne.0 ) + 2 * isw ( iand(32,int(sctrl%lncol)) .ne.0 
     .) + 3 * isw ( iand(64,int(sctrl%lncol)) .ne.0 )

      lcplxp = 0
      lekkl = int(sctrl%pfloat)

      if (lso .ne. 0) lcplxp = 1
      if (isum(nbas,lldau,1) .ne. 0) lcplxp = 1
      ndhamx = ndham*nspc
      nspx = nsp / nspc

!! related to the error in rseq. clean up RSEQ_ERROR at first because
!! this error is only fatal when it is finaly converged. takaoNov28_2010
      idummy=unlink('RSEQ_ERROR')

ctakao sep2010 !complicated.
c      I think nspc=2 for non-colinear case.
c      allocate(evl(ndham,2))
      allocate(evl(ndhamx,nspx))

      i_copy_size=size(sham%ldham)
      call i8icopy(i_copy_size,sham%ldham,1,ldham,1)
      pwmode=sham%pwmode
      pwemin=sham%pwemin
      pwemax=sham%pwemax

      onesp = 0
      nfbn(1) = 0
      nfbn(2) = 0
      stdo = lgunit(1)
      stdl = lgunit(2)

c ipr and stdl are processor dependent in mpi modes.
c      print *,'zzzzzzzzzz =',ipl,stdl,lgunit(1),lgunit(2),lgunit(3)
      ldos = iand(1,int(sctrl%ldos))
      lrsig = int(sham%lsig)
      lfrzw = isw ( iand(16,int(sctrl%lbas)) .ne.0 )
c      print *,'xxx' ,lfrzw,sctrl%lbas
c      print *,'ttt1: lrsig,ldos,lfrzw=',lrsig,ldos,lfrzw
c      print *,'ttt: lrsig,ldos,lfrzw=',lrsig,ldos,lfrzw
c      print *,'ttt: sham%lsig       =',sham%lsig
      if(ldos - mod(sctrl%ldos/1, 2)/=0) stop 'ttt xxx111'
      if(lfrzw- mod(sctrl%lbas/16,2)/=0) stop 'ttt xxx222'
      if(lrsig- sham%lsig/=0) stop 'ttt xxx333'
      nspec= sctrl%nspec

CKI!      allocate(sspecv(nspec+1)) ! I don't know why but nspec+1 is necessary to avoid run time error.
CKI      allocate(sspecv(nspec+1)) !+1 makes it safer. (see below)
CKI      do i=1,nspec
CKIc Be careful for this procedure; the size uspec_size() is probably larger than the size of s_sspec.
CKIc Thus we needed to allocate it as sspecv(nspec+1).
CKIc Without this, deallocate(sspecv) caused a problem---this may mean memeory destruction.
CKI         call dcopy(uspec_size(),sspec(1+(i-1)*uspec_size()),1,sspecv(i),1)
CKI         if(debug) print *,'ttt: sspecv%idmod=',i,sspecv(i)%nr,sspecv(i)%a
CKI      enddo
CKI      deallocate(sspecv)
c      stop 'xxxxxxxxxxxxxxx aaaaaaaaaaa'
c      print *,'ttt fxstsigm=',fxst('sigm')



      if (procid .eq. master) then
        if (lrsig .ne. 0 .and. fxst('sigm') .ne. 1) then
          call info0(1,1,0,' bndfp (warning): '//
     .    'no sigm file found ... LDA calculation only')

          sham%lsig=0

          lrsig = 0
        endif
      endif

      call mpibc1(lrsig,1,2,mlog,'rdsigm','lrsig')
      jobgw= -999
      epsovl = sham%oveps

      lpdiag = isw(cmdopt('--pdiag',7,0,strn))

#if MPIK
      lpdiag = 0
#endif
      if(llmfgw) then
        i = 7
        if (cmdopt('-jobgw=',i,0,strn)) then
          if (.not. a2bin(strn,jobgw,2,0,' ',i,-1)) call
     .    rxs2('BNDFP: failed to parse "',strn(1:30),'%a"')
        endif
      endif


ckino Dec.8.2011:         sv_p_oorhat => spot%sv_p_oorhat 

ckino Dec.8.2011:         zv_p_osmrho => spot%zv_p_osmrho 

ckino Dec.8.2011:         zv_p_osmpot => spot%zv_p_osmpot 

c       osoptc=spot%osoptc

C     icls = igets('optic cls',w(osoptc))
      pi = 4d0*datan(1d0)
      if (iprint() .ge. 20) call awrit2('%N --- BNDFP:  '//
     .'begin iteration %i of %i ---',' ',80,stdo,iter,maxit)
      if (nbas.gt.nbx) call rxi('bndfp: nbx exceeded, need',nbas)
      if (ndham.le.0)
     .call rx('bndfp: hamiltonian matrix has zero dimension')
c

#if MPI
C MPI Process configuration
      if (lpdiag .eq. 1) then
        nblk = 16
        dims(1) = 0
        dims(2) = 0
        call MPI_DIMS_CREATE(numprocs,2,dims,ierr)
        npcol = dims(1)
        nprow = dims(2)
        if (iprint() .ge. 30) then
          call awrit3(
     .    ' MPI creating process configuration .. nprow=%i npcol=%i,'//
     .    ' blocking factor %i',
     .    ' ',256,lgunit(1),nprow,npcol,nblk)
        endif
      endif
#endif

      call dvset(eterms,1,20,-99d0)
      eterms(19) = 0d0
      i_copy_size=size(sham%eterms)
      call dcopy(i_copy_size,eterms,1,sham%eterms,1)
      nvl=spot%nlml
      nchan=spot%nlma
      eks = 0

c      i_copy_size=size(slat%plat)
c      call dcopy(i_copy_size,slat%plat,1,plat,1)
c      i_copy_size=size(slat%qlat)
c      call dcopy(i_copy_size,slat%qlat,1,qlat,1)
      alat=slat%alat
      plat=slat%plat
      qlat=slat%qlat
      i_copy_size=size(slat%nabc)
      call i8icopy(i_copy_size,slat%nabc,1,ngabc,1)
      vol = slat%vol

c      call fftz30(n1,n2,n3,k1,k2,k3)
      k1 = n1
      k2 = n2
      k3 = n3


C ... for BZ integration
      i_copy_size=size(sbz%nkabc)
      call i8icopy(i_copy_size,sbz%nkabc,1,nkabc,1)
      nkp=sbz%nkp
      ntet=sbz%ntet
ckino Dec.16.2011:         iv_p_oidtet => sbz%iv_p_oidtet 

      lmet=sbz%lmet
      mpsord=sbz%n
      esmear=sbz%w
      nevmx=sbz%nevmx
      efmax=sbz%efmax
      fsmom=sbz%fsmom

      print *,'ttt nevmx w=',nevmx,esmear
      esmear = esmear + iabs(mpsord)
      if (mpsord .lt. 0) esmear = -esmear
      mpsord = int(dabs(esmear))
      if (esmear .lt. 0) mpsord = -mpsord


      ndos=sbz%ndos
      i_copy_size=size(sbz%dosw)
      call dcopy(i_copy_size,sbz%dosw,1,dosw,1)
      ef0=sbz%ef
      def=sbz%def


ckino Dec.8.2011:         rv_p_oqp => sbz%rv_p_oqp 

ckino Dec.8.2011:         rv_p_owtkp => sbz%rv_p_owtkp 


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
        call info2(20,0,0,
     .  ' BNDFP: generating density in energy window %2:1d',dosw,0)
        lfrce = 0
        lpnu = 0
        efmax = 1d2
        nevmx = ndham
        if (lrout .eq. 0)
     .  call rx('--window incompatible with no output density')
        call info0(20,0,0,' Delete band weights file ...')
        ifi = fopna('wkp',-1,4)
        call dfclos(ifi)
      endif
C     Always need DOS if just sampling
      if (lmet .ne. 0 .and. .not. ltet) ldos=1
      qbg = sctrl%zbak(1)

      alfa = 0

C     Switch to plot bands at specified qp
      if (cmdopt('--band',6,0,strn)) then
        plbnd = 1
        plbopt = strn(7:)
        lrout = 0
        lfrce = 0
        nkp = 0
        numq = 1
        allocate(ifbls_iv(ndham*nspc*2))
        ifbls_iv(:)=0

      else
        plbnd = 0
      endif

C     Sanity checks
      if (lmet .eq. 1) call rx('bndfp: lmet=1 not implemented')
Ckino isanrg is logical function,       call isanrg(lmet,0,4,'bndfp:','lmet', .true.)
      l_dummy_isanrg=isanrg(lmet,0,4,'bndfp:','lmet', .true.)

C --- Define local arrays used in the generation of the potential ---
C  10 continue
      lnoxc = 0
      if (jobgw .eq. 2) lnoxc = 1
      allocate(qmom_rv(nvl))
      allocate(gpot0_rv(nvl))
      allocate(vval_rv(nchan))
      allocate(hab_rv(nab*n0*nsp*nbas))
      allocate(vab_rv(nab*n0*nsp*nbas))
      allocate(sab_rv(nab*n0*nsp*nbas))
      allocate(ppnl_rv(nppn*n0*nsp*nbas))
      allocate(fes1_rv(3*nbas))

Ckino dfaugm alllocates and initializes variables
         call dfaugm ( nbas , lcplxp , lso , ssite , sspec , sv_p_osig 
     .   , sv_p_otau , sv_p_oppi ) 




C --- Make the potential sans XC part ---
      if(llmfgw) then
        if (jobgw .eq. 1 .or. jobgw .eq. -999) then
          call info(20,1,0,' Make potential without XC part ...',0,0)
          call togpr()
          i = 1 + 10*lfrzw + 100
Ckino dfaugm alllocates and initializes variables
             call dfaugm ( nbas , lcplxp , lso , ssite , sspec , sv_p_osigx 
     .       , sv_p_otaux , sv_p_oppix ) 



c          call defcc(ospotx, -k1*k2*k3*nsp)
          allocate(spotx(k1,k2,k3,nsp))
          spotx=0d0
c          call defcc(osmrox, -k1*k2*k3*nsp)
c          call dcopy(k1*k2*k3*nsp*2,w(osmrho),1,w(osmrox),1)
Ckino mkpot calls locpot. and locpot calls augmat. 
Ckino augmat calculates sig,tau,ppi.
               call mkpot ( nbas , ssite , sspec , slat , sham , spot , 0 , 
ckino Dec.8.2011:       .         lcplxp , k1 , k2 , k3 , zv_p_osmrho , sv_p_oorhat , qbg , spotx 
ckino Dec.8.2011:       .         lcplxp , k1 , k2 , k3 , zv_p_osmrho , spot%sv_p_oorhat , qbg , spotx 
ckino Dec.28.2011:       .         lcplxp , k1 , k2 , k3 , spot%zv_p_osmrho , spot%sv_p_oorhat , qbg , spotx 
ckino Dec.28.2011:       .         lcplxp , k1 , k2 , k3 , spot%zv_a_osmrho , spot%sv_p_oorhat , qbg , spotx 
     .         lcplxp , k1 , k2 , k3 , spot%zv_a_osmrho , spot%sv_a_oorhat , qbg , spotx 
     .         , qmom_rv , vconst , vrmt , sv_p_osigx , sv_p_otaux , sv_p_oppix 
     .         , ppnl_rv , hab_rv , vab_rv , sab_rv , qval , qsc , gpot0_rv 
     .         , vval_rv , fes1_rv , i , vorb , nlibu , lmaxu , lldau , nsp,0 ) 






          call togpr()
c          call rlse(osmrox)
        endif
C --- readin job switch
        if (jobgw .eq. -999 ) then
          write(stdo,*) ' lmfgw: input one of the following jobs:'
c        write(stdo,*) '  -1 : creates files',
c     .    ' GWinput, QIBZ, KPTin1BZ'
          write(stdo,*) '   0 : init mode; creates files',
     .    ' SYMOPS, LATTC, CLASS, NLAindx, ldima'
          write(stdo,*) '   1 : GW setup mode'
          write(stdo,*) '   4 : band mode '
          write(stdo,*) '   5 : eigenvalue-only mode '
          write(stdo,*) ' job?'
          read (5,*) jobgw
        endif
C ... Mode -1 ... only after mkpot because qval made in mkpot
c      if (jobgw .eq. -1) then
c        call sugwin(ssite,sspec,slat,sham,w(osgw),sbz,nbas,qval)
c        call rx0('lmfgw, job -1')
c      endif
      endif



c xxxxxxxxxxx new sigma block














C --- Make the potential and total energy terms ---
      i = 1 + 10*lfrzw + 100*lnoxc
      if (cmdopt('--wrhomt',8,0,strn)) then
        i = i + 10000
      else if (cmdopt('--wpotmt',8,0,strn)) then
        i = i + 20000
      endif

      if(llmfgw) i = i + 10000

ctakao sep2009
c      print *,'vvv'

ccccccccccccccccccccccccc
ctakao test mkekin repeat
c          print *,'vvv takao mkekin repeat init 01 ---'
c          call mkekin(nbas,ldim,ssite,sspec,slat,sham,osig,otau,oppi,
c     .      lcplxp,oqkkl,k1,k2,k3,vconst,w(osmpot),w(osmrho),sev,sumtvxxx)
c          stop 'xxxxxxxxxxxxxxxx'
ccccccccccccccccccccccccc

      call togpr()
Ckino mkpot calls locpot. and locpot calls augmat. 
Ckino augmat calculates sig,tau,ppi.
            call mkpot ( nbas , ssite , sspec , slat , sham , spot , lfrce 
ckino Dec.8.2011:       .      , lcplxp , k1 , k2 , k3 , zv_p_osmrho , sv_p_oorhat , qbg , zv_p_osmpot 
ckino Dec.8.2011:       .      , lcplxp , k1 , k2 , k3 , zv_p_osmrho , spot%sv_p_oorhat , qbg , zv_p_osmpot 
ckino Dec.8.2011:       .      , lcplxp , k1 , k2 , k3 , spot%zv_p_osmrho , spot%sv_p_oorhat , qbg , zv_p_osmpot 
ckino Dec.28.2011:       .      , lcplxp , k1 , k2 , k3 , spot%zv_p_osmrho , spot%sv_p_oorhat , qbg , spot%zv_p_osmpot 
ckino Dec.28.2011:       .      , lcplxp , k1 , k2 , k3 , spot%zv_p_osmrho , spot%sv_p_oorhat , qbg , spot%zv_a_osmpot 
ckino Dec.28.2011:       .      , lcplxp , k1 , k2 , k3 , spot%zv_a_osmrho , spot%sv_p_oorhat , qbg , spot%zv_a_osmpot 
     .      , lcplxp , k1 , k2 , k3 , spot%zv_a_osmrho , spot%sv_a_oorhat , qbg , spot%zv_a_osmpot 
     .      , qmom_rv , vconst , vrmt , sv_p_osig , sv_p_otau , sv_p_oppi 
     .      , ppnl_rv , hab_rv , vab_rv , sab_rv , qval , qsc , gpot0_rv 
     .      , vval_rv , fes1_rv , i , vorb , nlibu , lmaxu , lldau , nsp,0 ) 







      call togpr()

ccccccccccccccccccccccccc
ctakao test mkekin repeat
c       print *,'vvv takao mkekin repeat init 02 ---'
c       call mkekin(nbas,ldim,ssite,sspec,slat,sham,osig,otau,oppi,
c     .      lcplxp,oqkkl,k1,k2,k3,vconst,w(osmpot),w(osmrho),sevxxx,sumtvxxx)
c          stop 'xxxxxxxxxxxxxxxx 111'
ccccccccccccccccccccccccc


C Commented by takao lwsig should be renewed. July2009
C C     lwsig0 = isign(1,lwsig) * mod(iabs(lwsig),10)
C C     lwsig1 = iabs(lwsig/10)
C       if (cmdopt('--evec',6,0,strn)) then
C         if (lwsig .ne. 0) call rxi('bandfp: --evec not allowed with '
C      .      //'--wsig:trans=',i)
C         strn2 = strn(7:)
C         dc = strn2(1:1)
C         call evcflg(dc,strn2,i,shftqp)
C         if (i .gt. 10) lwsig = 2 + mod(i,10)
C       endif
C
C
C C     Transformation modes: read qp and jump to start
C       if (lwsig .ne. 0) then
Ckino isanrg is logical function, C         call isanrg(lwsig,-1,5,'BNDFP:','lwsig',.true.)
C         l_dummy_isanrg=isanrg(lwsig,-1,5,'BNDFP:','lwsig',.true.)
C         onesp = 0
C         lwtkb = 0
C
C         ifi = fopna('qpts',-1,0)
C         call getqp(0,ifi,nkp,nkabc,lshft,i,w,w,w)
C         call defrr(oqp,3*nkp)
C         call getqp(1,ifi,nkp,nkabc,lshft,i,w(oqp),w,w)
C
C         goto 50
C       endif

      elind = sham%elind

      if (elind .lt. 0d0)
     .elind=-(3*pi**2*(qval-qsc-qbg)/vol)**.66666d0*elind
      dmxp(33) = elind

c      call defrr(oevlall, ndham*nsp*nkp)
      allocate( evlall(ndham,nsp,nkp))

C --- Setup for BZ integration ---
      if (plbnd .eq. 0) then

cccccccccccccccccccccc
c         print *,'nsp,nspc,nkp,nkp=',nsp,nspc,nkp,nkp
cccccccccccccccccccccc


        call subzi(lmet,ltet,lrout.gt.0,ndham,nsp,nspc,nkp,qval-qbg,
c     .    nevmx,lwtkb,eferm,oevlall,owtkb,numq,lswtk,oswtk,ef0,def,mpsord,
ckino Dec.8.2011:       .  nevmx,lwtkb,eferm,rv_p_owtkb,numq,lswtk,rv_p_oswtk,ef0,def,mpsord,
ckino Dec.8.2011:       .  nevmx,lwtkb,eferm,rv_p_owtkb,numq,lswtk,rv_a_oswtk,ef0,def,mpsord,
     .  nevmx,lwtkb,eferm,rv_a_owtkb,numq,lswtk,rv_a_oswtk,ef0,def,mpsord,
     .  dabs(esmear-mpsord))

        if (lwtkb .eq. -1)
     .  call info(20,0,0,' Start first of two band passes ...',0,0)
        if (lwtkb .eq. 1) then
          if (ef0 .ne. eferm)
     .    call info(20,0,1,'%8pReplace ef0 with file ef=%;6d',eferm,0)
          ef0 = eferm
        endif
c#if MPIK
c        if (epsovl .ne. 0) then
c          allocate(eomin_rv(nsp*nkp))
c
c        endif
c#endif
      else
        nkp = 0
        ldos = 0
        lwtkb = -1
        icls = 0
      endif
c      print *,'ttt after subzi nevmx=',nevmx

c$$$C ... Set switch to write file sigii
c$$$#ifndef LMFGW
c$$$      if (lrsig .ne. 0 .and. plbnd .eq. 0 .and. procid .eq. master) then
c$$$        i = fopna('sigii',-1,0)
c$$$        rewind i
c$$$        call pack1('ham sigp:9.9',sham,1d0)
c$$$      endif
c$$$#endif

C --- More k-point independent local arrays ---
      if (ldos .eq. 0 .and. lmet .eq. 0) ndos = 1
      if (ndos .lt. 0) ndos = -ndos
      allocate(dos_rv(ndos*2*nsp))
      if (lrout .ne. 0) then
         call dfqkkl ( nbas , ssite , sspec , numq , sv_p_oqkkl ) 

         if ( lekkl .eq. 1 ) call dfqkkl ( nbas , ssite , sspec , numq 
     .   , sv_p_oeqkkl ) 

        allocate(srout_zv(k1*k2*k3*numq*nsp))
        allocate(fh_rv(3*nbas))
        allocate(fes2_rv(3*nbas))
      endif
      call suham2(sctrl,slat,sspec,ssite,sham,spot,sstrn)
      print *,'end of suham2'

C --- Options for core level specta (CLS) ---
C     Allocate aus for all qp
      if (cmdopt('--cls',5,0,strn)) then
        if (lmet .gt. 0 .and. (lmet .ne. 2 .and. lmet .ne. 3))
     .  call rx('For CLS restart with METAL=2 or 3')
        icls = 1
        clsopt = strn(6:)
        efmax = 1d3
        nevmx = ndhamx
        if (lrout .eq. 0) call rx('bndfp: need output density for cls')
        call suclst(nsitmx,nbas,nsp,ssite,sspec,clsopt,
     .  isite,iclsl,iclsn,nsites)
        efmax = 1d3
        nevmx = ndhamx
Changenglob        nlmax = nglob('mxorb') / nglob('nkaph')
        nlmax = globalvariables%mxorb / globalvariables%nkaph
        if (16*3*nlmax*ndham*nbas*nsp*nkp/1000000 .gt. 24 .and. procid
     .  .eq. master)
     .  call awrit6(' CLS: %iMb memory for aus: nlmax=%i ndham=%i '//
     .  'nsistes=%i nsp=%i nkp=%i',' ',120,
     .  i1mach(2),16*3*nlmax*ndham*nsites*nsp*nkp/1000000,
     .  nlmax,ndham,nsites,nsp,nkp)
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
      ifiese=0
      if(lrsig/=0) then
        ifiese=9039
        inquire(unit=ifiese,opened=nexist)
        if(nexist) 
     &  stop 'bndfp: I assumed ifio=9039 is not used(too sloppy)'
        open(ifiese,file='ESEAVR',status='old')
      endif




C --- Start loop over k points; also, re-entry for second band pass ---
   99 continue
      ebot = 1000d0
        call surho ( nbas , ssite , sspec , lmet , ldos , lrout , lekkl 
     .  , numq , sv_p_oqkkl , sv_p_oeqkkl , k1 , k2 , k3 , srout_zv , 
     .  ndos , dos_rv , sumev , sumqv ) 



      if (lfrce .gt. 0) then
        call dpzero(frc, 3*nbas*numq)
        call dpzero ( fh_rv , 3 * nbas ) 

      endif


      if (lswtk .eq. 1) then
        efmax = 1d3
        nevmx = ndhamx
ckino Dec.8.2011:           call dpzero ( rv_p_oswtk , ndhamx * nkp ) 
         call dpzero ( rv_a_oswtk , ndhamx * nkp ) 

      endif

C --- Setup moments file : write header ---
      nl = int(sctrl%nl)

      nfilem = 0
      if (procid .eq. master) then
        nfilem = fopna('moms',-1,4)
      endif
C     if (cmdopt('--mull',6,0,strn)) nfilem = fopna('moms',-1,4)
C     if (cmdopt('--pdos',6,0,strn)) nfilem = fopna('moms',-1,4)

      if(debug) print *,'vvvvv 000aaa 22222222qqq '
      if (nfilem .gt. 0) then
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
     .    , moddos , nsites , lsites , lmxch , nchan , chan0_iv, 0 ) 

          if(allocated(chan_iv)) deallocate(chan_iv)
          allocate(chan_iv(nchan))
          chan_iv=chan0_iv(1:nchan)
          deallocate(chan0_iv)
          call iomomn(.true.,2,.false.,1,nspc,1,1,i)
        endif
        if(debug) print *,'vvvvv 111aaa 22222222qqq '

C   ... If switch '--pdos'; get sites, number of channels
        if (cmdopt('--pdos',6,0,strn)) then
          mulopt = strn(7:)

          ng=slat%nsgrp

          nchmx = min(1024,nbas*nl**2)
          allocate(chan_iv(nchmx))

          nll = 0
           call sumlst ( 0 , nchmx , nbas , ng , ssite , sspec , mulopt 
     .     , moddos , nsites , lsites , lmxch , nchan , iwdummy , nll ) 

          nlmax = (lmxch+1)**2
          i = nlmax*ndham*3*nsp*nbas
          if (16*i*nkp/1000000 .gt. 24 .and. procid .eq. master)
     .    call info(20,0,0,' PDOS: %iMb memory for aus: nlmax=%i',
     .    i/1000000,nlmax)
c           allocate(ausp_zv(i))
c           if (i<0) ausp_zv(:)=0.0d0

c           if (allocated(ausp_zv)) deallocate(ausp_zv)

          if (cmdopt('--mull',6,0,strn))
     .    call rx('--pdos and --mull not allowed in conjunction')
          if (cmdopt('--cls',5,0,strn))
     .    call rx('--pdos and --cls not allowed in conjunction')
          call iomomn(.true.,2,.false.,1,nspc,1,1,i)
        endif
        if(debug) print *,'vvvvv 22222222qqq111 '
        if (procid .eq. master) i = iomoms(-nfilem,nl,nsp,nspc,nkp,
     .  ndham,i,1,0,1,0,0,0,0,0,0d0,0d0,0d0,0d0,0d0,0d0)
      endif

      if(debug) print *,'vvvvv 22222222qqq '

C ... Case only generate bands at supplied qp: setup
      if (plbnd .ne. 0) then
C       Try and read Fermi level from file
        if (procid .eq. master) then
          ifi = fopna('wkp',-1,4)
ckino getef is integer          call getef(ifi,0,ef0)
          i_dummy_getef= getef(ifi,0,ef0)
          call fclr('wkp',ifi)
        endif

        call mpibc1(ef0,1,4,.false.,'bndfp','ef0')
        iopq = 0
C       suqlst in MPIK mode; returns cumulative number of k-points
#if MPIK
        iopq = 2
#endif
        if (cmdopt('--onesp',7,0,strn) .and. nspc .eq. 1) onesp = 1
        i = nsp
        if (onesp .ne. 0 .or. nspc .eq. 2) i = 1
C       In parallel mode, suqlst call only serves to generate nkp

cccccccccccccccccccccccccccccccccccccccccccccccccccc
        if(debug) print *,'vvvvv nkp=',nkp
cccccccccccccccccccccccccccccccccccccccccccccccccccccc

ckino plotting wavefunction
        write(*,*)'allocate listwf , size=',ndham
        if (allocated(listwf))  deallocate(listwf) 
        allocate(listwf(ndham))
        n_listwf=ndham; listwf=0

        if (procid .eq. master) then
ckino  
           write(*,*) 'calling suqlst'
           call suqlst ( plbopt , iopq , ndhamx , ef0 , i , iwdummy , nfbn 
     .     , ifbls_iv , nkp , qp , onesp ) 

ckino       fat band and write psi
           n_listwf=0
ckino           write(*,*)'checking --writewf option'
ckino 
ckino  Usage: also use fat band option
ckino    --band~col=1:8~fn=syml --writewf:band=1:10
ckino         if you want to write wavefunction (F0) from the 1st to 10th band
ckino
ckino         q-list must be set in the syml file
ckino         For example, if you want to write the wavefunction at the gamma point
ckino            ' 2   0 0 0  1 0 0 '
ckino            ' 0   0 0 0  0 0 0 '
ckino           is necessary. It also write the wavefunction at k=(1 0 0). 
ckino
           strn=' '
           if (cmdopt('--writewf',9,0,strn)) then
              call w_psir_cmdopt(strn(10:),ndham,n_listwf,listwf)
              write(*,*)'--writewf option, ndham,n_listwf=',ndham,n_listwf
              write(*,'(8I5)')  listwf(1:n_listwf)
           endif

        endif
cccccccccccccccccccccccccccccccccccccccccccccccccccc
        if(debug) print *,'vvvvv nkp222=',nkp,iopq
c        stop 'xxxxxxxxxxxxxxxxxxxxxx'
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
        call mpibc1(nkp,1,2,mlog,'bndfp','nkp')
        if (nkp .le. 0) call rx0('bndfp')
        call mpibc1(nfbn,2,2,mlog,'bndfp','nfbn')
        call mpibc1(onesp,1,2,mlog,'bndfp','onesp')

ckino plotting wavefunction
        call mpibc1(n_listwf,1,2, mlog,'bndfp','n_listwf')
        if (n_listwf>0) then
          call mpibc1(listwf,n_listwf,2, mlog,'bndfp','listwf')
        endif

C MPIK: Setup to assemble all k-points into single list with qp table
#if MPIK
        if (nfbn(1) .gt. 0 .or. nfbn(2) .gt. 0) then
          call rx('Cannot use color weights with MPIK')
        endif

        call suqlsm(i)
C       Re-allocate qp and evl arrays
ckino Dec.8.2011:           allocate(rv_p_oqp(abs(3*nkp)))
ckino Jan.4.2012
         if (associated(sbz%rv_p_oqp)) deallocate(sbz%rv_p_oqp)
         allocate(sbz%rv_p_oqp(abs(3*nkp)))
ckino Dec.8.2011:           if (3*nkp<0) rv_p_oqp(:)=0.0d0
         if (3*nkp<0) sbz%rv_p_oqp(:)=0.0d0

c        call defrr(oevlall, ndham*nsp*nkp)

c        if (epsovl .ne. 0) then
c          allocate(eomin_rv(nsp*nkp))
c        endif
        call info2(20,1,1,
     .  ' bndfp:  MPIK band plotting mode %i:  %i q-points',i,nkp)
C       Loop through all qp; accumulate vector of qp.
C       Use i2 in place of nkp to preserve nkp
        if (procid .eq. master) then
C       iq = running index to big qp list, i1 = index to current line
          iq = 0
          i2 = 0
  199     continue
          i = 1
          call pshpr(0)
           call suqlst ( plbopt , 1 , ndhamx , ef0 , i , iwdummy , nfbn 
     .     , ifbls_iv , i2 , qp , onesp ) 


          call poppr
          if (i2 .gt. 0) then
            call pshpr(0)
            do  i1 = 1, i2
              iq = iq+1
               call suqlst ( plbopt , 1 , ndhamx , ef0 , i , iwdummy , nfbn 
     .         , ifbls_iv , i2 , qp , onesp ) 


ckino Dec.8.2011:                 call dpscop ( qp , rv_p_oqp , 3 , 1 , 3 * iq - 2 , 1d0 ) 
               call dpscop ( qp , sbz%rv_p_oqp , 3 , 1 , 3 * iq - 2 , 1d0 ) 

            enddo
            call poppr
            call suqlsm(i)
            if (i .ne. 3) goto 199
          endif
        endif
ckino Dec.8.2011:           call mpibc1 ( rv_p_oqp , 3 * nkp , 4 , .false. , 'bndfp' , 'qp' 
         call mpibc1 ( sbz%rv_p_oqp , 3 * nkp , 4 , .false. , 'bndfp' , 'qp' 
     .   ) 

#endif
      endif

C ... For insulator, valence band top and conduction band bottom
      evtop = -99999
      ecbot = -evtop

C ... Setup for spin-orbit coupling
      if (lso .ne. 0) then
C        call sumlst(10,nchmx,nbas,ng,ssite,sspec,mulopt,moddos,
C     .    nsites,lsites,lmxch,nchan,w,nll)
Changenglob        nlmax = nglob('nlmax')
        nlmax = globalvariables%nlmax
        if(allocated(orbtm_rv)) deallocate(orbtm_rv)
        allocate(orbtm_rv(nl*nsp*nbas))
        orbtm_rv(:)=0.0d0
      endif

C ... Setup for case sigma or evecs written to disk.
C     No integrated quantities accumulated.  qp read from file
C     Also sigma must be written in (iq,isp) order, opposite to the
C     order in which they would be generated here.
C     Requires second loop over (iq,isp) and filtering of isp
C     First pass should have onesp=0 and lwtkb=0
C     Next line is re-entry point for 2nd spin when writing sigma
C
C Commented by takao july2009
C   50 continue
C       if (lwsig .ne. 0) then
C         onesp = onesp + 1
C
C         call info5(20,1,0,' BNDFP:  '//
C      .    '%?#(n==4)#Write evals,evecs to file##%-1j'//
C      .    '%?#(n==3)#Write LDA evals,evecs to file##%-1j'//
C      .    '%?#(n==1|n==2)#Write sigm(LDA)##%-1j'//
C      .    '%?#(n==5)#Write sigm(orb)##%-1j'//
C      .    '%?#(n==-1)#Rotate sigm back to orbital basis##%-1j'//
C      .    '%j for %i qp%?#n==2#, spin 2##',
C      .    lwsig,nkp,onesp,0,0)
C
C C       endif
C       endif

C --- Get qplist. This is only for no MPI mode yet...
      if(allocated(qplist)) deallocate(qplist)
      allocate(qplist(3,nkp))
      do iq = 1, nkp
        if (plbnd .ne. 0) then
          i = nsp
          if (onesp .ne. 0 .or. nspc .eq. 2) i = 1
           call suqlst ( plbopt , 0 , ndhamx , ef0 , i , iwdummy , nfbn 
     .     , ifbls_iv , nkp , qp , onesp ) 
        else
ckino Dec.8.2011:             call dpscop ( rv_p_oqp , qp , 3 , 3 * iq - 2 , 1 , 1d0 ) 
           call dpscop ( sbz%rv_p_oqp , qp , 3 , 3 * iq - 2 , 1 , 1d0 ) 
        endif
        qplist(:,iq)=qp
      enddo
      print *,'-------- qplist --------'
      do iq=1,nkp
        write(6,"(i5,3f8.3)")iq,qplist(:,iq)
      enddo


!! --- New sigma inerpolation mode. Simplified.
      if(debug) print *,'vvv: new sigma=',newsigmasw()
      sigmamode = mod(lrsig,10) .ne. 0
      print *,'sigmamode=',sigmamode
c --- get index for hamiltonian for m_hamindex takao june2009
c this must be used for the case of QSGW.
c      if(sigmamode.and.iter==1) then

c      if(iter==1.and.onlyonce.and.sham%pwmode.ge.10) then
c      if(iter==1.and.(onlyonce.or.plbnd/=0)) then
      if(iter==1.and.onlyonce) then
        onlyonce=.false.
c takao to get all qp. I think suqlst looks too complicated...
c this block is just to get list of qp. suqlst should be cleaned up in future.
        nbasp=sctrl%nbasp
        nl=sctrl%nl
ckino Dec.8.2011:          iv_p_oidxsh => sham%iv_p_oindxo 
        print *,'goto gen_hamindex nkp=',nkp
        call gen_hamindex ( sham , sbz , ssite , sspec , slat , nl , 
ckino Dec.8.2011:       .   nbasp , iv_p_oidxsh , ldim , qplist , nkp , siginit, .not.sigmamode) 
ckino Dec.27.2011:      .   nbasp , sham%iv_p_oindxo , ldim , qplist , nkp , siginit, .not.sigmamode) 
ctakao feb2012          .   nbasp , sham%iv_a_oindxo , ldim , qplist , nkp , siginit, .not.sigmamode) 
c     .   nbasp , sham%iv_a_oindxo , ldim , qplist , nkp , siginit)
     .   nbasp , sham%iv_a_oindxo , ldim , siginit)

#if 1
Ckino ib=atom index
Ckino ltab=L 
Ckino ktab=  =1 for EH, =2 for EH2, =3 for  cMTO
        write(6,*)
        write(6,*) ' --- Hamiltonian index ---'
        write(6,*) ' ib l  k offl(iorb)+1  offl(iorb)+2*l+1  trim(spec)'
        do  iorb = 1, norbt
          ib = ibastab(iorb)

          is=ssite(ib)%spec

          do i_spacks=is,is
            call spacks_copy('u',sspec(i_spacks)%name,is,is,spid,i_spacks)
          enddo

          write(6,"(3i3,2x,2i5,3x,a)") 
     &    ib,ltab(iorb),ktab(iorb), 
     &    offl(iorb)+1, offl(iorb)+2*ltab(iorb)+1, trim(spid)
        enddo
#endif

        write(6,*) ' ngrp=', ngrp
        !do ig = 1, ngrp
        !write(6,*) ig
        !do i=1,3
        !  write(6,"(3f13.5,3x,f13.5)") symops(i,1:3,ig),ag(i,ig)
        !enddo
        !enddo
      endif

c --- generate \Sigma for q points on all the mesh points. takao
c      nhmax = ldim + napwmx !napwmx is in m_hamindex
c      lwsig=0
      if(sigmamode.and.siginit) then
        ifi = 1
        if (procid .eq. master) ifi = fopna('sigm',-1,4)
        read(ifi,err=99,end=99) nsp_dummy,ndimh_dummy,nk1,nk2,nk3
        allocate( sfz(nk1,nk2,nk3,ndham,ndham,nsp),qsmesh(3,nk1,nk2,nk3) )
        rewind ifi
        !rsrnge = dgets('ham rsrnge',sham) !Real-space range
        !print *,'rsrnge=',rsrnge
        !print *,'vvv',nk1,nk2,nk3,ndham,ndham,nsp
        !stop 'xxx'
        rsrnge = sham%rsrnge
        print *,'rsrnge=',rsrnge
        write(6,"(a,2i5)")" --- ldim(=dim of lmto), napwmx= ",ldim,napwmx
        if (ldim .ne. ndham) then !june2009 takao
          write(6,"(a,2i5)")" --- ndham (>=maxinum of ham dimension) =",ndham
        endif
c        if(ndham/=ndham) stop "we now assume ndham = ndmam"
        write(6,*)"Now we go into rdsigm2 for PMT mode...lrsig=" ,lrsig

        call rdsigm2(nbas,nsp,ndham,slat,sham,sbz,ifi,
     &  nk1,nk2,nk3,ldim,qsmesh,sfz) !sfz is the self-energy roteted in the full BZ.
        ! rdsigm2 use m_hamindex (napwk
c        deallocate(sfz)
        call fclose(ifi)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c$$$cc This is original version for readin sigm--->ohrs or so are used in hambls
c$$$cc rdsigm2+seneinterp is equivalent with rdsigm.
c$$$c --- original
c$$$C ... Read and store self-energy addition to LDA potential
c$$$C     Read lwsig: special modes that write sig, Z to disk.
c$$$         lwsig = 0
c$$$         if (mod(lrsig,10) .ne. 0) then
c$$$C       Real-space range
c$$$         rsrnge = dgets('ham rsrnge',sham)
c$$$         ifi = 1
c$$$         if (procid .eq. master) ifi = fopna('sigm',-1,4)
c$$$         if (ldim .ne. ndham) then !june2009 takao
c$$$           write(6,*)"ldim(dim of lmto)  ndham(max of ham dimension)=",ldim,ndham
c$$$           write(6,*)"go to rdsigm for PMT mode..."
c$$$           call rx('rdsigm not ready for PW basis')
c$$$         endif
c$$$         call rdsigm(lrsig,nbas,nsp,ndham,slat,sham,sbz,ifi,rsrnge,lwsig)
c$$$         if (lwsig .ne. 0) then
c$$$          lrout = 0
c$$$          lfrce = 0
c$$$          lwtkb = 0
c$$$          call rxx(plbnd.ne.0,'incompatible options, lwsig and plbnd')
c$$$         endif
c$$$         call fclose(ifi)
c$$$         call phmbl3(1,0,0,0,0,0,0,0)
c$$$         endif
c$$$         print *,'lwsig=',lwsig
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c --- pair-table and weight generator
c      nknknk = nk1*nk2*nk3
c      allocate(npair(nbas,nbas),range(nbas,nbas),pos(3,nbas))
c      call upack2('lat opos',slat,opos,0)
c      call dcopy(nbas*3,w(opos),1,pos,1)
c      call pairtakao(1,pos,nbas,plat,nknknk,npair, npairmx,range) !,nlat,qwgt)
c      allocate( nlat(3,npairmx,nbas,nbas), qwgt(npairmx,nbas,nbas))
c      call pairtakao(2,pos,nbas,plat,nknknk,npair, npairmx,range, nlat,qwgt)
c      deallocate(range)
c
        allocate(npair(nbas,nbas),pos(3,nbas))

ckino Dec.8.2011:           rv_p_opos => slat%rv_p_opos 


ckino Dec.8.2011:           call dcopy ( nbas * 3 , rv_p_opos , 1 , pos , 1 ) 
ckino Jan.04.2012:           call dcopy ( nbas * 3 , slat%rv_p_opos , 1 , pos , 1 ) 
        call dcopy ( nbas * 3 , slat%rv_a_opos , 1 , pos , 1 ) 
        npairmx=nk1*nk2*nk3*2
        do
          npairmx= npairmx + nk1*nk2*nk3*.5
          allocate( nlat(3,npairmx,nbas,nbas), qwgt(npairmx,nbas,nbas) )
          call gennlat(pos,nbas,plat,nk1,nk2,nk3,npairmx,ok,npair,nlat,qwgt)
          if(ok) exit
          deallocate( nlat, qwgt )
        enddo
C... This is needed only for hambls.F case (store ohrs)
        if(.not.newsigmasw()) then
          allocate( sfzx(nk1,nk2,nk3,ndham,ndham,nsp) )
          sfzx=sfz
          ifi = 1
          if (procid .eq. master) ifi = fopna('sigm',-1,4)
          call seneinterp(nbas,nsp,ndham,slat,sham,sbz,ifi,rsrnge,
     &    nk1,nk2,nk3,sfzx) !,qp,isp,sene)
          call fclose(ifi)
          deallocate(sfzx)
        endif
        siginit=.false.
      endif !end of sigmamode if


C--- GW driver must be after rdsigm2 for sigm mode
      if(llmfgw) then
C --- GW driver ---
c       osgw=spot%osgw
        nchan=spot%nlma
C ... Modes 0..4
ckino Dec.8.2011:                 call sugw ( ssite , sspec , slat , sham , nbas , ndham , zv_p_osmpot 
ckino Dec.28.2011:                 call sugw ( ssite , sspec , slat , sham , nbas , ndham , spot%zv_p_osmpot 
               call sugw ( ssite , sspec , slat , sham , nbas , ndham , spot%zv_a_osmpot 
     .         , vconst , lcplxp , sv_p_osig , sv_p_otau , sv_p_oppi , ppnl_rv 
     .         , vrmt , spotx , sv_p_osigx , sv_p_otaux , sv_p_oppix , jobgw 
     .         ,nk1,nk2,nk3,sfz,sigmamode) 
        deallocate(spotx)
c        call rlse(oqmom)
C      if (jobgw .eq. 1) then
C        jobgw = 2
C        call info(20,1,0,' ... Remake potential without vxc',0,0)
C        goto 10
C      endif
        call tcx('bndfp')
        call rx0('bndfp')
      endif



C --- do loop for iq and isp ---------------
      if(debug) print *,' --- goto do loop 2010 ---'
ctakao
      if(allocated(nevls)) deallocate(nevls)
      allocate(nevls(nkp,nsp))
      nevls=0
C Start k and isp loop (first loop in parallel mode)
C     h,s,evecs are dimensioned ndimh in this loop (may be q dependent)
C     evl is dimensioned evl(ndham,2)
C     evlall, w(owtkb), w(oswtk) are dimensioned (ndham,nsp,nkp)
#if MPE & MPIK
      ierr = MPE_LOG_EVENT(EVENT_START_KLOOP,procid,"k-loop")
#endif

#if MPIK
      call info0(30,1,0, ' ... Start MPI k-loop')
      sttime = MPI_WTIME()
      allocate (kpproc(0:numprocs), stat=ierr)
      call dstrbp(nkp,numprocs,1,kpproc(0))
      iqini = kpproc(procid)
      iqend = kpproc(procid+1)-1
c      do 2010 iq = kpproc(procid), kpproc(procid+1)-1
c        if (iq .eq. kpproc(procid)) then
c          if (mlog) then
c            call gettime(datim)
c            call awrit4(' bndfp '//datim//' Process %i of %i on '
c     .        //shortname(procid)(1:namelen(procid))//
c     .        ' starting k-points %i to %i',' ',256,lgunit(3),
c     .        procid,numprocs,kpproc(procid),kpproc(procid+1)-1)
c          endif
c        endif
#else
c      nsmidb = ndham
      iqini=1
      iqend=nkp
#endif

      do 2010 iq = iqini, iqend
#if MPIK
#ifdef KINODEBUG
        write(*,'(a,i3,5i5)')'debug iq=',procid,iq,iqini,iqend
#endif
        if (iq .eq. kpproc(procid)) then
          if (mlog) then
            call gettime(datim)
            call awrit4(' bndfp '//datim//' Process %i of %i on '
     .      //shortname(procid)(1:namelen(procid))//
     .      ' starting k-points %i to %i',' ',256,lgunit(3),
     .      procid,numprocs,kpproc(procid),kpproc(procid+1)-1)
          endif
        endif
#endif
        isqp = nsp*(iq-1)
C       Get qp either from qp list or read from suqlst

#ifndef MPIK
c        if (plbnd .ne. 0) then
c          i = nsp
c          if (onesp .ne. 0 .or. nspc .eq. 2) i = 1
c          call suqlst(plbopt,0,ndhamx,ef0,i,w,nfbn,w(oifbls),nkp,qp,
c     .      onesp)
c        else
c          call dpscop(w(oqp),qp,3,3*iq-2,1,1d0)
c        endif
        qp=qplist(:,iq)
#else
ckino Dec.8.2011:           call dpscop ( rv_p_oqp , qp , 3 , 3 * iq - 2 , 1 , 1d0 ) 
         call dpscop ( sbz%rv_p_oqp , qp , 3 , 3 * iq - 2 , 1 , 1d0 ) 

#endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cC... range check qqq is in [ikc1/nk1, (ikc1+1)/nk1],... qpara is [0,1],[0,1],[0,1]
c      print *,' --- goto qparam--- '
c      call qparam(qp,qsmesh,nk1,nk2,nk3,plat, ikc,qpara,igshift)
c      ! eight corners -->  qsmesh(1:3,ikc(1,ix1,ix2,ix3),ikc(2,ix1,ix2,ix3),ikc(3,ix1,ix2,ix3))
c      call interparam(qpara, wgc)
c      cycle
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C   ... Loop over majority, minority spins
C       In the SO case, still a loop isp=1..2
C       isp=1 generates just the 11 block
        if(debug) print *,'goto spin loop nsp=',nsp
        do  isp = 1, nsp
          if(debug)  print *,'vvvvvv isp loop isp=',isp
          if (onesp .eq. 0 .or. isp .eq. onesp) then
ctakao do we need to do shorbz? --> furthermore, shorbz is probably wrong.---> see new algolism shortn3.
            qpo=qp
            if(.not.noshorbz()) then
              if (isp .eq. 1 .or. isp .eq. onesp) call shorbz(qp,qp,qlat,plat)
            endif
            if(debug)  print *,'vvvvvv xxx cccccccccc1 pwemax pwmode',pwemax,pwmode
            if(allocated(igv2x)) deallocate(igv2x)
C   ... For this qp, G vectors for PW basis and hamiltonian dimension
            if (pwemax .gt. 0 .and. mod(pwmode,10) .gt. 0) then
              pwgmin = dsqrt(pwemin)
              pwgmax = dsqrt(pwemax)
              call pshpr(1)
              call dpzero(xv,3)
              if (mod(pwmode/10,10) .eq. 1) call dpcopy(qp,xv,1,3,1d0)
              !print *,'xxx xxx 1111'
              call gvlst2(alat,plat,xv,0,0,0,pwgmin,pwgmax,0,
     .        0,napw,napw,dum,dum,dum,dum)
              call poppr
c              call defi(oigv2,3*napw)
              allocate(igv2x(3,napw))
              allocate(kv_iv(3*napw))

              call pshpr(iprint()-10)
              !print *,'xxx xxx 2222'
              call gvlst2 ( alat , plat , xv , 0 , 0 , 0 , pwgmin , pwgmax 
     .        , 0 , 2 , napw , napw , kv_iv , dum , dum , igv2x ) 

              !print *,'xxx xxx 3333'
              deallocate(kv_iv)

              call poppr
              ndimh = ldim + napw  !ndimh is iq-dependent.
              if (mod(pwmode,10) .eq. 2) ndimh = napw
              if (ndimh .gt. ndham) then
                call fexit2(-1,111,'%N Exit -1 : BNDFP: '//
     .          'ndimh=%i exceeds ndham=%i.  Try increasing '//
     .          'input NPWPAD',ndimh,ndham)
              endif
            else
              ndimh = ldim
              allocate(igv2x(1,1)) !dummy
            endif
            ndimhx = ndimh*nspc !this is iq-dependent.
            if(debug)  print *,'vvvvvv xxx cccccccccc2222'

C   ... Make Hamiltonian and overlap matrices
            nqsig=sham%nqsig
ckino Dec.8.2011:               rv_p_oqsig => sham%rv_p_oqsig 
cki            if (oqsig .eq. 0) oqsig = 1
            ispc = min(isp,nspc) !nspc=2 for non-colinier case
            !nspc=1
            if(debug) print *,'zzzzzzzzzzzzz11111111'
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
C Commented by takao lwsig should be renewed.
C C       lwsig=-1: sigm LDA -> orbital basis
C C       lwsig= 1: sigm orbital -> LDA basis
C C       lwsig= 2: sigm orbital -> LDA basis, high energy parts replaced
C C       lwsig= 3: Write evecs of LDA hamiltonian to file
C C       lwsig= 4: Write evecs of hamiltonian to file
C C       lwsig= 5: sigm orbital -> orbital, high energy parts replaced
C C       hambls returns with evecs,evals of the LDA hamiltonian in s,h
C C       Transform sigm to LDA basis: hambls returns sigm(LDA) in h
C         if (lwsig .eq. 1) then
C           i = i + 3000
C         elseif (lwsig .eq. 2) then
C           i = i + 4000
C         elseif (lwsig .eq. -1 .or. lwsig .eq. 3) then
C           i = i + 1000
C           if (nspc .eq. 2)
C      .      call rx('transform sigm not implemented in noncoll case')
C         elseif (lwsig .eq. -1 .or. lwsig .eq. 4) then
C           i = i + 6000
C         elseif (lwsig .eq. 5) then
C           i = i + 5000
C         endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccc
ctakaox
c        if(iq==1) then
c              ifio = fopna('ovlmat',-1,4)
c        endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccc

! Determine interpolated self-energy sene at qp from sfz
            if(sigmamode) then
              if(debug) print *,' goto seneinterp2: qp=',qp
              allocate(sene(ndimh,ndimh))
              call seneinterp2(qp,nbas,nk1,nk2,nk3,plat, !npair,nlat,qwgt,npairmx,qsmesh,
     i        sfz(:,:,:,:,:,isp),ldim,ndimh, 
     o        sene) 
            endif

C cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C C--- This works OK for lmf mode when noshrbz=.true. because no interpolation is needed.
C C    Check comment out "sene=0" and "if(j1<=ldim.and.j2<=ldim) cycle" or not ---it depends on what you want to test.
C           ! sene=0d0
C           write(6,"(a,3f13.5)")'    qp=',qp
C           do ik1=1,nk1
C           do ik2=1,nk2
C           do ik3=1,nk3
C              write(6,"(a,3f13.5)")' qsmesh=',qsmesh(1:3,ik1,ik2,ik3)
C              if(sum(abs(qp-qsmesh(1:3,ik1,ik2,ik3)))<1d-6) goto 1013
C           enddo
C           enddo
C           enddo
C           stop 'xxxxx not qp=qsmesh xxxx'
C  1013    continue
C           do j1=1,ndimh !ldim+napw
C           do j2=1,ndimh !ldim+napw
C             if(j1<=ldim.and.j2<=ldim) cycle
C             sene(j1,j2)= sfz(ik1,ik2,ik3,j1,j2,isp)
C           enddo
C           enddo
C C cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

Ckino hambl calls augmbl.  augmbl calculates C_kl(=array b in the program)
Ckino finally makes F~F~=F0F0+(F1F1-F2F2), which is overlap matrix, s. 
Ckino Note that F2=Hankel head at a site + Hankel tail contributions from the other site.

C --- Set up hamiltonian.
ccccc     if (mod(mode,10) .eq. 0) then
            if(debug) print *,' mode ispc=',i,ispc
            if( newsigmasw() ) then
              print *,' --- hambl --- sigmamode=',sigmamode
              jx = 0
              if (lso == 1) jx = 10
                  call hambl ( jx , nbas , ssite , sspec , slat , sham , isp , 
ckino Dec.8.2011:       .            qp , k1 , k2 , k3 , zv_p_osmpot , vconst , sv_p_osig , sv_p_otau 
ckino Dec.28.2011:      .            qp , k1 , k2 , k3 , spot%zv_p_osmpot , vconst , sv_p_osig , sv_p_otau 
     .            qp , k1 , k2 , k3 , spot%zv_a_osmpot , vconst , sv_p_osig , sv_p_otau 
     .            , sv_p_oppi , lcplxp , alfa , ndimh , napw , igv2x , hamm ( 1 
     .            , 1 , ispc ) , ovlm ( 1 , 1 , ispc ) , hamm ( 1 , 1 , 3 ) ) 
              if(sigmamode) hamm(:,:,ispc) = hamm(:,:,ispc)+sene
              if(sigmamode) deallocate(sene)
ccccc     .    h(1,1,ispc),s(1,1,ispc),h(1,1,3))
c      endif
c      lrsig = mod(mode/10,10)
c      lwsig = mod(mode/1000,10)
c      ismidb = 0
Ckino isanrg is logical function, c      call isanrg(lwsig,0,5,'hambls:','lwsig',T)
c      l_dummy_isanrg=isanrg(lwsig,0,5,'hambls:','lwsig',T)

C ... Store LDA eigenvectors, eigenvalues in s and h and exit
c      if (lwsig .eq. 1) then
c        allocate(zqi(ndimh,ndimh))
c        allocate(z0(ndimh,ndimh))
c        allocate(e(ndimh))
c        call zhevx(ndimh,ndimh,h(1,1,ispc),s(1,1,ispc),1,T,ndimh,9d9,i,
c     .    z0,F,e,ndimh,zqi)
c        call dcopy(ndimh**2*2,zqi,1,s(1,1,ispc),1)
c        call dcopy(ndimh,e,1,h,1)
c        deallocate(z0,zqi,e)
cc        goto 999
c      endif
c      if (lrsig .eq. 0) goto 999
            else
C ... Interpolation of sene= \Sigm-Vxc at given q.
ckino                    isize_hamm(3)=size(hamm)
                   isize_hamm(3)=size(hamm,dim=3)
                   call hambls ( i , nbas , ssite , sspec , slat , sham , isp , 
ckino Dec.8.2011:       .             ispc , qp , k1 , k2 , k3 , rv_p_oqsig , nqsig , zv_p_osmpot , 
ckino Dec.8.2011:       .             ispc , qp , k1 , k2 , k3 , rv_p_oqsig , nqsig , spot%zv_p_osmpot , 
ckino Dec.28.2011:       .             ispc , qp , k1 , k2 , k3 , sham%rv_p_oqsig , nqsig , spot%zv_p_osmpot , 
     .             ispc , qp , k1 , k2 , k3 , sham%rv_p_oqsig , nqsig , spot%zv_a_osmpot , 
     .             vconst , sv_p_osig , sv_p_otau , sv_p_oppi , lcplxp , lso , alfa 
     .             , ndimh , napw , igv2x , hamm , ovlm ,isize_hamm, ismidb ) 





c        nsmidb = min(nsmidb,ismidb) !????
              if(sigmamode) deallocate(sene)
              if(oncewrite(3)) print *,' end of hambls mode=',i
            endif

ccccccccccccccccccccccccccccccccccc
c        print *,'xxxxxxxxxxxxxxx AAAAAAAAAAQQQQQQQQQ'
c        goto 2011
cccccccccccccccccccccccccccccccccccc


C       call zprm('evec, LDA',2,w(os),ndimh,ndimh,ndimh)

cccccccccccccccccc
ctakaox test
c        if(iq==1) then
c              ifio = fopna('ovlmat',-1,4)
c        endif
c        call dpdump(w(os),ndimh**2*2,-ifio)
ccccccccccccccccc

C Commented by takao lwsig should be renewed july 2009 takao
C C   ... Mode to write sigma or evecs to file
C         if (lwsig .ne. 0) then
C C         Write header information
C           call upack('bz lshft',sbz,lshft,0,0,0,0)
C           if (iq .eq. 1 .and. isp .eq. 1) then
C             if (lwsig .eq. 1 .or. lwsig .eq. 2 .or. lwsig .eq. 5) then
C               ifi = fopna('sigm2',-1,4)
C             elseif (lwsig .eq. 3 .or. lwsig .eq. 4) then
C               ifi = fopna('evec',-1,4)
C             endif
C             rewind ifi
C             call iosigh(0,lwsig,nsp,ndimh,nk1,nk2,nk3,nkp,
C      .        lshft(1),lshft(2),lshft(3),-ifi)
C           endif
C
C C         lwsig=1,2 : dump sigma(LDA basis) into file sigm2
C C         lwsig=3,4 : dump LDA eigenvectors into file evec
C C         lwsig=5   : dump sigma(orb basis) into file sigm2
C           if (lwsig .ge. 1 .and. lwsig .le. 5) then
C             if (lwsig .eq. 1 .or. lwsig .eq. 2 .or. lwsig .eq. 5) then
C               ifi = fopna('sigm2',-1,4)
C             else
C               ifi = fopna('evec',-1,4)
C             endif
C             write(ifi) qp
C             call dpdump(w(oh),ndimh,-ifi)
C             call dpdump(w(os),ndimh**2*2,-ifi)
C
C           elseif (lwsig .ne. 0) then
C             call rxi('BNDFP: band mode lwsig=',lwsig)
C           endif
C           goto 30
C         endif

C   ... In the noncollinear case, this marks the end of looping over
C       separate spins.
C       For the rest of the isp loop,
C       jsp=isp in the collinear case; jsp=1 in the noncollinear
C       Thus jsp should be used in place of isp
C       isp serves as a flag for the noncollinear case
            if (ispc .ne. nspc) goto 30
            jsp = isp
            if (ispc .eq. 2) jsp = 1

            if (i .eq. -1) lpdiag = 2

            if(debug) print *,'bndfp: lpdiag=',lpdiag

Ccomment out this test block takao
C C   ... debugging ... test of trotwf
C #if TESTRWF
C         call upack('ham ldham ooffH oindxo',sham,ldham,ooffH,oidxsh,0,0)
C         call upack('lat plat oistab oag osymgr nsgrp',slat,plat,oistab,
C      .    oag,og,nsgrp)
C         call upack('bz oqp owtkp oipq ostar',sbz,oqp,owtkp,oipq,ogstar,
C      .    0)
C         call upack('bz nkp nkabc lshft',sbz,nkp,nkabc,lshft,0,0)
C C   ... Make is,ifac,qb,qlat
C         do  8  i = 1, 3
C     8   llshft(i) = lshft(i) .ne. 0
C         call pshpr(0)
C         call bzmsh0(plat,llshft,0,nk1,nk2,nk3,is,ifac,rb,qb)
C         call poppr
C
C         rsrnge = dgets('ham rsrnge',sham)
C         call trothf(nl,nbas,jsp,nsp,w(ooffH),w(oidxsh),w(oistab),nk1,
C      .    nk2,nk3,w(oipq),nkp,w(oqp),w(owtkp),rsrnge,ndimh,0,slat,ssite,
C      .    sspec,sham,w(osmpot),osig,otau,oppi,vconst,w(og),w(oag),nsgrp,
C      .    w(ogstar),lshft,ifac,qb)
C C        call trotwf(ssite,sspec,slat,sham,w(osmpot),vconst,osig,otau,
C C     .    oppi,nl,nbas,jsp,w(ooffH),w(oidxsh),w(oistab),nk1,nk2,nk3,
C C     .    w(oipq),nkp,w(oqp),ndimh,plat,w(og),w(oag),w(ogstar),
C C     .    ifac,qb)
C         call rx('done')
C #endif


C   ... Diagonalize and add to density
            if (lpdiag==0.or.lpdiag==2) then
              if (allocated(t_zv)) deallocate(t_zv)
              allocate(t_zv(ndimhx**2))
            endif
#ifdef KINODEBUG
             write(*,*)'kino, allocate t_zv' ,ndimhx**2
#endif

            if (mod(iq,10) .ne. 1) call pshpr(iprint()-6)
#ifndef MPIK
            call info5(30,0,0,' bndfp:  kpt %i of %i, k=%3:2,5;5d'//
     .      '%?#n#   ndimh = %i##',
     .      iq,nkp,qp,mod(pwmode/10,10),ndimh)
#endif
            nmx = min(nevmx,ndimhx)
            nmx = max(nmx,nevec) !nevec is given in the case of double path formalism
                                 !Nov25 2010.
            if (lwtkb .eq. -1) nmx = -1
C       Need all eigenvalues if 'fat bands' plotting mode
            if (nfbn(1) .gt. 0) then
              nmx = ndimhx
              efmax = 99999
            endif

            if(napw/=0.and.ncutovl/=0) then
              if(lpdiag ==1 .or. lpdiag==2) then
                stop 'bndfp: lpdiag-PW mode with OVNCUT/=0 is not yet'
              endif
            endif
            if(debug) print *,'lpdiag=',lpdiag
            if(lpdiag .eq. 1) then
              call rx('parallel diag is not yet.'//
     &                'Need to implement. See bndfp.F.')
C$$$          call rxx(nspc.ne.1,'parallel diag not implemented for noncol')
C$$$#if MPE & MPI
C$$$          ierr = MPE_LOG_EVENT(EVENT_START_PZHEV,procid,"pzhev")
C$$$#endif
C$$$          if (iq .eq. 1) then
C$$$            call info0(20,0,0,' bndfp:  diagonalise with SCALALPACK ..')
C$$$          endif
C$$$c         call pzhev(T,ndimh,oh,os,nblk,nprow,npcol,efmax,nmx,nev,
C$$$          call pzhev(T,ndimh,hamm,ovlm,nblk,nprow,npcol,efmax,nmx,nev,
C$$$     .               evl(1,jsp),ot)
C$$$
C$$$          nevl = ndimh
C$$$#if MPE & MPI
C$$$          ierr = MPE_LOG_EVENT(EVENT_END_PZHEV,procid,"pzhev")
C$$$#endif
            elseif (lpdiag .eq. 2) then
              call dcopy ( ndimhx * * 2 * 2 , ovlm , 1 , t_zv , 1 ) 
              nev = ndimhx
              nmx = ndimhx
              nevl = ndimhx
#ifdef USEOPTIONAL
               call phmbls ( mode=2 , ndimh=ndimhx , eval=evl ( 1 , jsp ) 
     .         , zhz= hamm ) 
#else
               call phmbls ( 2 , ndimhx , evl ( 1 , jsp ) , iwdummy , iwdummy 
     .         , iwdummy , iwdummy , iwdummy , hamm ) 
#endif
              if(debug) print *,'vvv lwtkb=',lwtkb
              if (lwtkb .ne. -1) then
                allocate(ww_zv(ndimhx**2*2))
#ifdef KINODEBUG
                write(*,*)'kino allocate ww_zv', ndimhx**2*2
#endif
                if (ldim .ne. ndham) call
     .          rx('blsig needs separate dim for wtkp,evec')
                  call blsig ( 1 + lrsig * 10 , nbas , sham , isp , nsp , nspc 
ckino Dec.8.2011:       .            , plat , qp , lwtkb , qval - qbg , iq , rv_p_owtkp , rv_p_owtkb 
ckino Dec.8.2011:       .            , plat , qp , lwtkb , qval - qbg , iq , sbz%rv_p_owtkp , rv_p_owtkb 
ckino Dec.27.2011:       .            , plat , qp , lwtkb , qval - qbg , iq , sbz%rv_p_owtkp , rv_a_owtkb 
     .            , plat , qp , lwtkb , qval - qbg , iq , sbz%rv_a_owtkp , rv_a_owtkb 
     .            , ndimh , t_zv , ww_zv ) 
                deallocate(ww_zv)
              endif
              if(debug) print *,'vvv2 lwtkb=',lwtkb
            else
              if(nspc .eq. 2) then
c                call sopert ( 0 , ndimh , nspc , ww_rv , hamm , hamm ) 
                call sopert2 ( 0 , ndimh , nspc, hamm , hamm ) 
                call sopert2 ( 0 , ndimh , nspc,  ovlm , ovlm ) 
              endif
              if(debug)  print *,' vvv zzzzz 555'
              allocate(ww_rv(11*ndimhx))
#ifdef KINODEBUG
              write(*,*) 'allocate ww_rv', 11*ndimhx
#endif
              if(napw == 0) then
c               call zhev(ndimhx,w(oh),w(os),T,T,nmx,efmax,nev,w(oww),F,-1,
c                print *,'goto zhev ndimh ndimhx=',ndimh,ndimhx
                call zhev ( ndimhx , hamm , ovlm , t , t , nmx , efmax , nev 
     .          , ww_rv , f , - 1 , evl ( 1 , jsp ) , t_zv ) 
c                print *,'end of zhev'
                nevl = ndimhx
              else
                nevl = -1
                ww_rv(1)=99999d0
c                call dvset ( ww_rv , 1 , 1 , 99999d0 ) 
c               call zhevo(ndimhx,ndimhx,w(oh),w(os),T,nmx,efmax,epsovl,
c     .         nevl,nev,evl(1,jsp),w(oww),ndimhx,w(ot))
!! takao's version of diagonization.
                call zhev_tk ( ndimhx - napw , ndimhx , hamm , ovlm , nmx , epsovl 
     .          , ncutovl , nevl , nev , evl ( 1 , jsp ) , ww_rv , ndimhx 
     .          , t_zv, iq==1, sham%delta_stabilize) 
c                print *,'end of zhev_tk'
c                eomin = dval ( ww_rv , 1 ) 
              endif
              deallocate(ww_rv)
            endif

!!       Pad evals between ndimh and ndham with a large positive number
!!       to avoid mixing up integration routines
            if (ndhamx .gt. nevl .and. nspc .eq. 2) then
              call dvset(evl,1+nevl,ndhamx,99999d0)
            elseif (ndham*nspc .gt. nevl) then
              call dvset(evl(1,jsp),1+nevl,ndham,99999d0)
            endif
            nevls(iq,isp)=nevl
            if(debug) print *,'iq isp nevls=',iq,isp,nevls(iq,isp)

#ifndef MPIK
c            if (epsovl .ne. 0) then
c              call info5(30,0,0,' Overlap''s smallest eigenvalue: %;3g.  '//
c     .        '%?#(n>0)#H dim reduced from %i to %i#H dim not reduced#',
c     .        eomin,ndimhx-nevl,ndimhx,nevl,0)
c            endif
            if(debug) print *,' goto print nevl nev nmx=',nevl,nev,nmx
            call prtev ( t_zv , nevl , evl ( 1 , jsp ) , nmx , efmax , 
     .      nev ) 

            if(debug) print *,' end of print ev'
            if (iprint() .ge. 110) then
              call yprm('evals',1,evl(1,jsp),1,ndhamx,nevl,1)
              if ( nev .gt. 0 ) call zprm ( 'evecs' , 2 , t_zv , ndimhx 
     .        , ndimhx , nev ) 

            endif
#endif
            if (mod(iq,10) .ne. 1) call poppr

C ...  In k-parallel mode, defer this section until all qp available
#ifndef MPIK
            ebot = dmin1(ebot,evl(1,jsp))
            i = max(1,nint(qval-qbg)/(3-nspc))
            evtop = max(evtop,evl(i,jsp))
            ecbot = min(ecbot,evl(i+1,jsp))
            if (lmet .eq. 0 .and. iq .eq. 1 .and. jsp .eq. 1) ef0 = evtop
            if (plbnd .eq. 0) then
#if ! (MPI | MPIK | MPE)
              if (ipr.ge.10 .and. iq.eq.1 .and. ipl.gt.0)
     .        write (stdl,712) (evl(i,jsp),i=1,nev)
  712         format('fp evl',8f8.4)
#endif
              if (lwtkb .ne. -1 .and. .not.lwndow) then
                if (iq .eq. 1 .and. jsp .eq. nsp .and.
     .          .not. cmdopt('--no-fixef0',11,0,strn)) then
                  ef00 = ef0
                  call fixef0(qval-qbg,jsp,1,nevl,ndham,evl,dosw,ef0)
                  if (jsp .eq. 2 .and. ef00 .ne. ef0 .and.
     .            lwtkb .eq. 0 .and. lmet .gt. 0 .and. lrout .ne. 0) then
                    if (procid .eq. master)
     .              call info0(10,1,1,
     .              ' ... Fermi level reset in second spin'//
     .              ' channel ... restart band pass')
                    goto 99
                  endif
                endif
!!takao comment out because of nevec algolism. Really OK? Nov25 2010
!!we may assume double-path method (weight first, then accumulate density(second band path)).
cC         Check for cases when nevmx is too small : i=2 => fatal error
c                i = 0
c                if (nevmx.ge.0 .and. lmet .ne. 0) then
c                  dum = evl(max(nev,1),jsp)
cC           if (ef0 .ge. dum) i = 2
c                  if (.not. ltet .and. ef0+5*dabs(esmear-mpsord) .gt. dum) i=2
c                  if (lmet.eq.4 .and. ef0+def+5*dabs(esmear-mpsord).gt.dum)i=2
c                endif
c                print *,' ttt: nev=',nev,ndimhx
c                if (i .eq. 2) then
c                  if (iprint() .gt. 0)
c     .            call awrit3('%N evl(nev=%i)=%;3d but '//
c     .            'ef0=%;3d ... restart with larger efmax or nevmx',
c     .            ' ',80,stdo,nev,evl(max(nev,1),jsp),ef0)
c                  call rx('bndfp')
c                endif
              endif
            endif
#endif
c##### !MPIK

            if( debug) print *, 'tttttttttt 11111 plbnd=',plbnd
C   ... Save data for this qp to moments file
            if (plbnd .eq. 0) then
C   ... Copy eigenvalues into array containing ev for all qp
              call dpscop(evl(1,jsp),evlall,ndhamx,
     .        1,1+ndham*(jsp-1+isqp),1d0)
c#if MPIK
c              if (epsovl .ne. 0) then
c                call dpscop ( eomin , eomin_rv , 1 , 1 , 1 + ( jsp - 1 + isqp 
c     .          ) , 1d0 ) 
c
c              endif
c#endif

              if( debug) print *, 'tttttttttt 22222'

C        call prmx('ev',evlall,ndham,ndham,isqp+2)
C        call prmx('ev',evlall,ndham,ndham,nkp*nsp)





C   --- Orbital magnetic moment (requires lso) ---
              if (lso .ne. 0 .and. lwtkb .ne. -1) then
                if (lwtkb .eq. 0) call
     .          rx('metal weights required for orb. moment calculation')
                allocate(auso_zv(nlmax*ndhamx*3*nsp*nbas))
#ifdef KINODEBUG
                write(*,*) 'allocate auso_zv',nlmax*ndhamx*3*nsp*nbas
#endif
                auso_zv(:)=0.0d0

                call makusq ( 0 , ssite , sspec , slat , sham , nbas , nbas , 
     .          0 , nlmax , ndham , ndimh , napw , igv2x , nev , nsp , 
     .          nspc , jsp , 1 , qp , t_zv , ppnl_rv , auso_zv ) 

                 call mkorbm ( ssite , sspec , jsp , nsp , nspc , nlmax , ndham 
ckino Dec.8.2011:       .           , nev , rv_p_owtkb , iq , nbas , ppnl_rv , auso_zv , nl , nkp 
     .           , nev , rv_a_owtkb , iq , nbas , ppnl_rv , auso_zv , nl , nkp 
     .           , orbtm_rv ) 


                deallocate(auso_zv)

              endif

              if( debug) print *, 'tttttttttt 333'

C   --- Mulliken analysis and partial DOS ---
              if (lwtkb .ne. -1) then
                if (procid .eq. master) then
                  if (cmdopt('--mull',6,0,strn)) then
#if MPIK
                    call rx('MULL not k-parallel.')
#endif
                    call rxx(nev.ne.ndimhx,'Mulliken requires all eigenvectors')
                    call rxx(nevmx.ne.ndhamx,'Mulliken requires nevmx=nhdamx')
                    allocate(doswt_rv(nchan*nevmx*nspc))
                    doswt_rv(:)=0.0d0


ckino Dec.8.2011:                       iv_p_oiprmb => sham%iv_p_oindxo 


ckino Dec.8.2011:                       call mullmf ( nbas , ssite , sspec , iv_p_oiprmb , t_zv , ndimh 
ckino Dec.27.2011:                       call mullmf ( nbas , ssite , sspec , sham%iv_p_oindxo , t_zv , ndimh 
                     call mullmf ( nbas , ssite , sspec , sham%iv_a_oindxo , t_zv , ndimh 
     .               , nspc , iq , jsp , moddos , nsites , lsites , lmxch , nchan 
     .               , chan_iv , lmdim , ndham , doswt_rv ) 


                    call iomomn(.true.,2,.false.,1,nspc,1,1,i)
                    i1 = iomoms ( - nfilem , nl , nsp , nspc , nkp , ndimh , i , 
     .              nspc , 1 , 1 , ndhamx , nevmx , nchan , nchan , nev , evl ( 1 
     .              , jsp ) , 0d0 , doswt_rv , 0d0 , 0d0 , 0d0 ) 

                    if (allocated(doswt_rv)) deallocate(doswt_rv)

                  elseif (cmdopt('--pdos',6,0,strn)) then
#if MPIK
                    call rx('PDOS not k-parallel.')
#endif
                    call rxx(nspc.ne.1,
     &              'pdos not implemented in noncoll case')
                    i = nlmax*ndham*3*nsp*nbas
                    allocate(ausp_zv(i))
                    ausp_zv(:)=0.0d0

                    call makusq ( 1 , ssite , sspec , slat , sham , nbas , nbas , 
     .              0 , nlmax , ndham , ndimh , napw , igv2x , nev , nsp , 
     .              nspc , jsp , 1 , qp , t_zv , ppnl_rv , ausp_zv ) 

                    allocate(doswt_rv(nchan*ndimhx*nspc))
                    doswt_rv(:)=0.0d0

                    call mkpdos ( moddos , ssite , sspec , jsp , nsp , nlmax , ndham 
     .              , nev , nchan , 1 , lsites , nsites , ppnl_rv , ausp_zv 
     .              , doswt_rv ) 

                    call iomomn(.true.,2,.false.,1,nspc,1,1,i)
#if MPIK
                    i1 = iomoms ( - nfilem , nl , nsp , nspc , 2 , ndimh , 1000 + 
     .              i , 1 , iq , isp , ndham , ndimh , nchan , nchan , nev , evl 
     .              ( 1 , jsp ) , 0d0 , doswt_rv , 0d0 , 0d0 , 0d0 ) 

#else
                    i1 = iomoms ( - nfilem , nl , nsp , nspc , 2 , ndimh , i , 1 
     .              , 1 , 1 , ndham , ndimh , nchan , nchan , nev , evl ( 1 , jsp 
     .              ) , 0d0 , doswt_rv , 0d0 , 0d0 , 0d0 ) 

#endif
                    if (allocated(doswt_rv)) deallocate(doswt_rv)

                    if (allocated(ausp_zv)) deallocate(ausp_zv)

                  elseif (nfilem .gt. 0) then
                    write (nfilem) 0, ndimhx
                    call dpdump(evl(1,jsp),ndimhx,-nfilem)
                  endif
                endif

                if( debug) print *, 'tttttttttt 5555'

C   ... Make new density matrix dmatu for LDA+U
                if (nlibu .gt. 0 .and. nev .gt. 0) then
                  if (lwtkb .eq. 0)
     .            call rx('metal weights required for LDA+U calculation')
                  nl = int(sctrl%nl)

                  nlmax = nl*nl
                  allocate(ausu_rv(2*nlmax*ndhamx*3*nsp*nbas))
#ifdef KINODEBUG
                  write(*,*)'allocate ausu_rv',2*nlmax*ndhamx*3*nsp*nbas
#endif
                  ausu_rv(:)=0.0d0

                  call makusq ( 0 , ssite , sspec , slat , sham , nbas , nbas , 
     .            0 , nlmax , ndham , ndimh , napw , igv2x , nev , nsp , 
     .            nspc , jsp , 1 , qp , t_zv , ppnl_rv , ausu_rv ) 

ckino Dec.8.2011:                     call mkdmtu ( ssite , sspec , rv_p_owtkb , jsp , iq , nsp , nspc 
                   call mkdmtu ( ssite , sspec , rv_a_owtkb , jsp , iq , nsp , nspc 
     .             , ndham , nlmax , nbas , nev , ppnl_rv , ausu_rv , dmatu , nlibu 
     .             , lmaxu , lldau ) 


                  if (allocated(ausu_rv)) deallocate(ausu_rv)

                endif
              endif


              if( debug) print *, 'tttttttttt 666666'

C   ... for core-level spectroscopy
              if (icls .ne. 0 .and. lwtkb .ne. -1) then
                call rxx(nspc.ne.1,'CLS not implemented in noncoll case')
#if MPI
                call rx('CLS only k-parallel')
#endif
                call makusq ( 0 , ssite , sspec , slat , sham , nbas , nsites 
     .          , isite , nlmax , ndham , ndimh , napw , igv2x , nev , 
     .          nsp , nspc , jsp , iq , qp , t_zv , ppnl_rv , ausc_zv   ) 

C       Accumulate output density and sampling DOS
C       Even if no output rho, still call addrbl to make DOS when lmet=4
              elseif (lwtkb.ne.-1 .and. (lrout.ne.0 .or. lmet.eq.4)) then
C         dum=frc(3,1)
                if(debug) print *,'goto addrbl nevl=',nevl
                call addrbl ( ssite , sspec , slat , sham , jsp , nsp , nspc 
ckino Dec.8.2011:       .          , qp , rv_p_owtkp , ndham , ndimh , napw , igv2x , lmet , lrout 
ckino Dec.27.2011:       .          , qp , sbz%rv_p_owtkp , ndham , ndimh , napw , igv2x , lmet , lrout 
     .          , qp , sbz%rv_a_owtkp , ndham , ndimh , napw , igv2x , lmet , lrout 
ckino Dec.8.2011:       .          , lwtkb , rv_p_owtkb , lswtk , rv_p_oswtk , iq , lfrce , ldos 
ckino Dec.8.2011:       .          , lwtkb , rv_p_owtkb , lswtk , rv_a_oswtk , iq , lfrce , ldos 
     .          , lwtkb , rv_a_owtkb , lswtk , rv_a_oswtk , iq , lfrce , ldos 
ckino Dec.8.2011:       .          , lekkl , k1 , k2 , k3 , zv_p_osmpot , vconst , sv_p_osig , sv_p_otau 
ckino Dec.28.2011:       .          , lekkl , k1 , k2 , k3 , spot%zv_p_osmpot , vconst , sv_p_osig , sv_p_otau 
     .          , lekkl , k1 , k2 , k3 , spot%zv_a_osmpot , vconst , sv_p_osig , sv_p_otau 
c    .          , sv_p_oppi , lcplxp , numq , qval - qbg , t_zv , evl , nevl 
     .          , sv_p_oppi , lcplxp , numq , qval - qbg , t_zv , evl , nev 
                     !Nov25 2010 nev should contain calculated number of eigenvectors.
     .          , ef0 , def , esmear , emin , emax , ndos , dos_rv , srout_zv 
     .          , sumqv , sumev , sv_p_oqkkl , sv_p_oeqkkl , frc ) 
C          print 399,iq,frc(3,1),frc(3,1)-dum
C  399     format(' after addrbl: frc(3,1)=',i4,2f12.6)
              endif

#ifndef MPIK
            elseif (plbnd .ne. 0) then

ckino
              write(*,*)  'case plbnd.ne.0 (no MPIK)'

ckino             nfbn(1)>0 means fat band mode => nmx=ndimh, all the evec are calculated 
ckino             n_listwf>0 means listwf is set by --writewf
              if (nfbn(1)>0 .and. n_listwf>0) then

                call rsibl_ev( ssite,sspec,slat,nbas,isp,qp,iq,ndimh,nspc
ckino Dec.27.2011:       .          ,  napw,igv2x,sham%iv_p_oindxo,nev,t_zv,k1,k2,k3
     .          ,  napw,igv2x,sham%iv_a_oindxo,nev,t_zv,k1,k2,k3
     .          ,  n_listwf,listwf    )

              endif

              i = nsp
              if (onesp .ne. 0 .or. nspc .eq. 2) i = 1
c          call suqlsw(nevl,jsp,i,evl(1,jsp))
              call suqlsw2(nevl,jsp,i,evl(1,jsp),qpo) !takao
              if (nfbn(1) .ne. 0) then
                if (ndimhx .ne. nevl)
     .          call rx('color weights not implemented when '//
     .          'nevl < hamiltonian dimension')
                call suqlse ( ndimhx , jsp , i , ndimhx , 1 , nfbn , ifbls_iv 
     .          , ndhamx , t_zv , evl ) 

              endif
              if (nfbn(2) .ne. 0) then
                if (ndimhx .ne. nevl)
     .          call rx('color weights not implemented when '//
     .          'nevl < hamiltonian dimension')
                call suqlse ( ndimhx , jsp , i , ndimhx , 2 , nfbn , ifbls_iv 
     .          , ndhamx , t_zv , evl ) 

              endif
#else
C   ... Save evals for this qp
            elseif (plbnd .ne. 0) then
ckino
              write(*,*)  'case plbnd.ne.0 (no MPIK)'

              call dpscop(evl(1,jsp),evlall,ndhamx,
     .        1,1+ndham*(jsp-1+isqp),1d0)
c              if (epsovl .ne. 0) then
c                call dpscop ( eomin , eomin_rv , 1 , 1 , 1 + ( jsp - 1 + isqp 
c     .          ) , 1d0 ) 
c
c              endif
#endif

            endif

c        if (ot .ne. os) call rlse(ot)
            if ( lpdiag = = 0.or.lpdiag = = 2 ) then 
              if (allocated(t_zv)) deallocate(t_zv)
            endif 

            if( debug) print *, 'tttttttttt 9999999999999999999',procid,iq

C end loop over isp (main loop in parallel mode)
   30       continue
            deallocate(igv2x)
          endif
        enddo
        if( debug) print *, 'tttttttttt 9999999999999999999 2',procid,iq

c        call rlse(oh)
C end loop over iq (main loop in parallel mode)
c      enddo
        deallocate(hamm,ovlm)!,hamsom)
 2010 continue

      if( debug) print *, 'tttttttttt 9999999999999999999 3',procid,iq

ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 2011 continue
c      call mpi_finalize(ierr)
c      stop 'xxxxxxxxxxxxxx BBBBB'
ccccccccccccccccccccccccccccccccccccccccccccccccccccc


c takao
cccccccccccccccccccccccccc
ctakaox
c      call fclose(ifio)
cccccccccccccccccccccccccc



C Commented by takao lwsig should be renewed july 2009
C C Second pass over iq for second spin when writing sigma
C       if (lwsig .ne. 0) then
C         if (onesp .lt. nsp) goto 50
C         if (lwsig .eq. 1 .or. lwsig .eq. 2 .or. lwsig .eq. 5) then
C           ifi = fopna('sigm2',-1,4)
C           call fclose(ifi)
C           if (lwsig .eq. 5)
C      .    call rx0('BNDFP:  sigm(orb basis) saved in file sigm2')
C           call rx0('BNDFP:  sigm(LDA basis) saved in file sigm2')
C         else if (lwsig .eq. 3 .or. lwsig .eq. 4) then
C           ifi = fopna('evec',-1,4)
C           call fclose(ifi)
C           call rx0('BNDFP:  LDA evecs saved in file evec')
C         endif
C C        if (lwsig .eq. -1) then
C C        ifi = fopna('sigm3',-1,4)
C C        call fclose(ifi)
C C        call rx0('BNDFP:  sigm(orbital basis) saved in file sigm3')
C C        endif
C C        if (lwsig .eq. 3) then
C C        call rx0('BNDFP:  U(LDA-QP) saved in file sigm2')
C C        endif
C C        if (lwsig .eq. 4) then
C C        call rx0('BNDFP:  H(QP) saved in file sigm2')
C C        endif
C       endif

c$$$C ... Set switch to write file sigii
c$$$#ifndef LMFGW
c$$$      if (lrsig .ne. 0 .and. plbnd .eq. 0 .and. procid .eq. master) then
c$$$        call upack1('ham sigp',sham,sigp)
c$$$        if (mpipid(0) .le. 1) call phmbl3(12,0,nsmidb,0,0,sigp,qp,qp(2))
c$$$        call fclose(fopna('sigii',-1,0))
c$$$      endif
c$$$#endif


C --- Second loop over qp, needed to make k-parallelisation possible: ---
C     You can't do this until you have all the evals.
#if MPIK
#if MPE
      ierr = MPE_LOG_EVENT(EVENT_END_KLOOP,procid,"k-loop")
#endif
      entime = MPI_WTIME()
      call info2(30,0,0, ' ... Done MPI k-loop: %;1d seconds elapsed',
     .entime-sttime,0)
      call info0(20,0,-1,' ... Sharing data between processes...')
      sttime = MPI_WTIME()

C     Share bands among all processors
      call xmpbnd(kpproc,ndham,nkp,nsp,evlall)
c      if (epsovl .ne. 0) then
c        call xmpbnd ( kpproc , 1 , nkp , nsp , eomin_rv ) 
c      endif
      if (icls .ne. 0 .and. lwtkb .ne. -1) then
        call rxx(nspc.ne.1,'CLS not implemented in noncoll case')
        call xmpbnd ( kpproc , nlmax * ndham * 3 * nsites , nkp , nsp 
     .  , ausc_zv ) 

      endif

ctakao
      call mpibc2(nevls,nkp*nsp,2,mlog,'bndfp','nevls')

#ifdef KINODEBUG
      call mpi_barrier(MPI_comm_world,ierr)
      write(*,*)'kino afeter mpibc2 nevls'
#endif

ckino 05Feb09, use kpproc later, move it.
c      deallocate(kpproc, stat=ierr)

C     Allreduce density-related quantities
      if (lrout .ne. 0) then
        call mpibc2(sumqv,6,4,mlog,'bndfp','sumqv')
        call mpibc2(sumev,6,4,mlog,'bndfp','sumev')
        call mpibc2 ( srout_zv , k1 * k2 * k3 * nsp * numq , 6 , mlog 
     .  , 'bndfp' , 'smrho' ) 

        if (lswtk .eq. 1) then
ckino Dec.8.2011:             call mpibc2 ( rv_p_oswtk , ndhamx * nkp , 4 , mlog , 'bndfp' 
           call mpibc2 ( rv_a_oswtk , ndhamx * nkp , 4 , mlog , 'bndfp' 
     .     , 'swtk' ) 

        endif
C     Allreduce qkkl
Changenglob        nkaph = nglob('nkaph')
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
     .         * nspc , 4 , mlog , 'bndfp' , 'qkkl' ) 

ctakao March1-2009
               if ( lekkl .eq. 1 ) call mpibc2 ( sv_p_oeqkkl( i , ib )%v , nelt 
     .         ( i ) * numq * nsp * nspc , 4 , mlog , 'bndfp' , 'eqkkl' ) 

            enddo
          endif
        enddo
C     Allreduce DOS, forces, dmatu, orbtm
        if ( ndos .gt. 0 ) call mpibc2 ( dos_rv , ndos * 2 * nsp , 
     .  4 , mlog , 'bndfp' , 'dos' ) 

        if (lfrce .ne. 0)
     .  call mpibc2(frc,3*nbas*numq,4,mlog,'bndfp','frc')
        if (nlibu .gt. 0)
     .  call mpibc2(dmatu,nsp*nlibu*(lmaxu*2+1)**2,6,mlog,'bndfp','dmatu')
        if ( lso .ne. 0 .and. lwtkb .ne. - 1 ) call mpibc2 ( orbtm_rv 
     .  , nl * nsp * nbas , 4 , mlog , 'bndfp' , 'orbtm' ) 

      endif
      if (lpdiag .eq. 2 .and. plbnd .eq. 0) then

        i_copy_size=size(sham%eterms)
        call dcopy(i_copy_size,sham%eterms,1,eterms,1)

        call mpibc2(eterms(19),1,4,mlog,'bndfp','rhosig')

        i_copy_size=size(sham%eterms)
        call dcopy(i_copy_size,eterms,1,sham%eterms,1)

      endif

      entime = MPI_WTIME()
      call info2(20,0,0,' MPI broadcast took %;1d sec',
     .(entime-sttime),0)

C     Write bands in bands-plotting case: loop over qp getting evals from array
      if (plbnd .ne. 0) then
        call info0(20,1,0,' Writing bands to bands file (MPIK) ...')
        if (procid .eq. master) then
C     iq = running index to big qp list, i1 = index to current line
          iq = 0
          i2 = 0
  299     continue
          i = nsp
          if (onesp .ne. 0 .or. nspc .eq. 2) i = 1
!      call suqlst(plbopt,0,ndhamx,ef0,i,w,nfbn,w(oifbls),i2,qp,onesp)
          if (i2 .le. 0) call rx0('bndfp')
          do  i1 = 1, i2
            iq = iq+1
            isqp = nsp*(iq-1)
!        call suqlst(plbopt,1,ndhamx,ef0,i,w,nfbn,w(oifbls),i2,qp,onesp)
            qp=qplist(:,iq)
!        call dpscop(qp,w(oqp),3,1,3*iq-2,1d0)
            do  isp = 1, nsp
              ispc = min(isp,nspc)
C         jsp=isp in the collinear case; jsp=1 in the noncollinear
C         Thus jsp should be used in place of isp
C         isp serves as a flag for the noncollinear case
              if ((onesp.eq.0 .or. isp.eq.onesp) .and. (ispc.eq.nspc)) then
                jsp = isp
                if (ispc .eq. 2) jsp = 1
                call dpscop(evlall,evl(1,jsp),ndhamx,
     .          1+ndham*(jsp-1+isqp),1,1d0)
                if (mod(i1,10) .ne. 1) call pshpr(iprint()-6)
                call info5(30,0,0,' bndfp:  kpt %i of %i, k=%3:2,5;5d',
     .          i1,nkp,qp,0,0)
                if (mod(i1,10) .ne. 1) call poppr
                i = nsp
                if (onesp .ne. 0 .or. nspc .eq. 2) i = 1
ctakao
c            call suqlsw(ndimhx,jsp,i,evl(1,jsp))
                print *,'nevls(iq,jsp)=',iq,jsp,nevls(iq,jsp)
                call suqlsw2(nevls(iq,jsp),jsp,i,evl(1,jsp),qp)
              endif
            enddo
          enddo
          if (i .ne. 3) goto 299
        endif
        call rx0('done')
      endif

C     Repeat loop for printout.  Put evals back into local array
      do  iq = 1, nkp
        isqp = nsp*(iq-1)
!        call dpscop(w(oqp),qp,3,3*iq-2,1,1d0)
        qp=qplist(:,iq)

        if(.not.noshorbz()) then
          call shorbz(qp,qp,qlat,plat)
        endif
        do  isp = 1, nsp
          ispc = min(isp,nspc)
C       jsp=isp in the collinear case; jsp=1 in the noncollinear
C       Thus jsp should be used in place of isp
C       isp serves as a flag for the noncollinear case
          if ((onesp.eq.0 .or. isp.eq.onesp) .and. (ispc.eq.nspc)) then
            jsp = isp
            if (ispc .eq. 2) jsp = 1
            call dpscop(evlall,evl(1,jsp),ndhamx,
     .      1+ndham*(jsp-1+isqp),1,1d0)
c            if (epsovl .ne. 0) then
c              call dpscop ( eomin_rv , eomin , 1 , 1 + ( jsp - 1 + isqp 
c     .        ) , 1 , 1d0 ) 
c            endif
            if (mod(iq,10) .ne. 1) call pshpr(iprint()-6)
            call info5(30,0,0,' bndfp:  kpt %i of %i, k=%3:2,5;5d',
     .      iq,nkp,qp,0,0)
c            if (epsovl .ne. 0) then
c              call info5(30,0,0,' Overlap''s smallest eigenvalue: %;3g',
c     .        eomin,0,0,0,0)
c            endif
ckino start, set origin  for bcast, 05Feb09
c        --- use ndimhx_iq and nev_iq below
c        --- I am afraid that this part can't execute correctly in the case of MPIK.
c        --- I fix wrong values of ndimhx and nev for iq.
            pidorigin=0
            if (kpproc(procid)<=iq .and. iq<= kpproc(procid+1)-1) then
              pidorigin=procid
            endif
c        ---  allreduce
            call mpibc2(pidorigin,1,2,mlog,'bndfp','pidorigin')
c        ---  pidorigin node has data
            if (procid.eq.pidorigin) then
              work_mpi(1) =  ndimhx
              work_mpi(2) =  nev
            endif
c        ---   bcast
            call mpibc3(work_mpi,2,2,pidorigin,mlog,'bndfp','ndimhx_iq')
            ndimhx_iq = work_mpi(1)
            nev_iq = work_mpi(1)

cki this causes SEGV, delete Jul08, 2010
cki            if (procid.eq.master) then
cki              call prtev(w(1),ndimhx_iq,evl(1,jsp),nevmx,efmax,nev_iq)
cki            endif
cki end
ckino       call prtev(w(1),ndimhx,evl(1,jsp),nevmx,efmax,nev)
ckino end

            if (mod(iq,10) .ne. 1) call poppr
            ebot = dmin1(ebot,evl(1,jsp))
            i = max(1,nint(qval-qbg)/(3-nspc))
            evtop = max(evtop,evl(i,jsp))
            ecbot = min(ecbot,evl(i+1,jsp))
            if (lmet .eq. 0 .and. iq .eq. 1 .and. jsp .eq. 1) ef0 = evtop
            if (plbnd .eq. 0) then
              if (ipr.ge.10 .and. iq.eq.1 .and. ipl.gt.0)
     .        write (stdl,712) (evl(i,jsp),i=1,nev_iq)
ckino 05Feb09     .    write (stdl,712) (evl(i,jsp),i=1,nev)
  712         format('fp evl',8f8.4)
              if (lwtkb .ne. -1 .and. .not. lwndow) then
                if (iq .eq. 1 .and. jsp .eq. nsp .and.
     .          .not. cmdopt('--no-fixef0',11,0,strn)) then
                  ef00 = ef0
                  call fixef0(qval-qbg,jsp,1,ndimh,ndham,evl,dosw,ef0)
                  if (jsp .eq. 2 .and. ef00 .ne. ef0 .and.
     .            lwtkb .eq. 0 .and. lmet .gt. 0 .and. lrout .ne. 0) then
                    if (procid .eq. master) call info0(10,1,1,
     .              ' ... Fermi level reset in second spin'//
     .              ' channel ... restart band pass')
                    goto 99
                  endif
                endif
C         Check for cases when nevmx is too small : i=2 => fatal error
                i = 0
                if (nevmx.ge.0 .and. lmet .ne. 0) then
                  dum = evl(max(nev_iq,1),jsp)
ckino 05Feb09            dum = evl(max(nev,1),jsp)
C           if (ef0 .ge. dum) i = 2
                  if (.not. ltet .and. ef0+5*dabs(esmear-mpsord) .gt. dum) i=2
                  if (lmet.eq.4 .and. ef0+def+5*dabs(esmear-mpsord).gt.dum)i=2
                endif
                if (i .eq. 2) then
                  call awrit3('%N evl(nev=%i)=%;3d but '//
     .            'ef0=%;3d ... restart with larger efmax or nevmx',
     .            ' ',80,stdo,nev_iq,evl(max(nev_iq,1),jsp),ef0)
ckino 05Feb09     .        ' ',80,stdo,nev,evl(max(nev,1),jsp),ef0)
                  call rx('bndfp')
                endif
              endif
            endif
C end second loop over isp
          endif
        enddo
C end second loop over iq
      enddo
ckino 05Feb09
      if (allocated(kpproc)) deallocate(kpproc, stat=ierr)
#endif
ctakao
      deallocate(nevls)
c###### second loop over qp (parallel k-points mode)



      if( debug) print *, 'tttttttttt 9999999999999999999 4'



C     For now, SO by site, not by class
      if (lwtkb .eq. 1 .and. lso .ne. 0) then
        allocate(ips_iv(nbas))

        i_copy_size=1;
        do i_spackv=1,nbas
          call spackv_array_copy_i8_i ( 'u' , ssite ( i_spackv ) %spec 
     .    , i_copy_size , i_spackv + 1 - 1 , ips_iv ) 

        enddo

        call iorbtm ( sspec , ips_iv , nl , nl , nbas , nsp , orbtm_rv 
     .  ) 

        if (allocated(ips_iv)) deallocate(ips_iv)

      endif

      if( debug) print *, 'tttttttttt 9999999999999999999 5'

C ... Case generating bands: find next block of qp
      if (plbnd .ne. 0) goto 99

      if (ipl.gt.0) write (stdl,715) nkp,ebot,qval,qbg,esmear
  715 format('nv nkp',i5,'  ebot',f9.4,'   qval',f10.4,'  qbg',f8.4,
     .'  smr',f8.4)
C ... End of k point loop

C     call zprm3('smrho after k-point loop',0,w(osrout),k1,k2,k3)

C --- Interpolate density to Fermi energy ---
      sev = sumev(1,1)
      if (lmet .eq. 4) then
         call mshn3p ( nbas , ssite , sspec , lmet , lrout , lfrce , qval 
     .   - qbg , ef0 , def , sumqv , sumev , n1 , n2 , n3 , k1 , k2 , 
     .   k3 , srout_zv , sv_p_oqkkl , frc , lrep ) 


C   ... Store val q & magnetic moment in sumqv(1) and sumqv(2)
        sumqv(2,1) = sumqv(1,1) - sumqv(1,2)
        sumqv(1,1) = sumqv(1,1) + sumqv(1,2)
C   ... Eigenvalue sum including entropy term
        sev = sumev(1,1) + sumev(2,1)
C   ... Remake sev,ef linearly interpolating tabulated sampling DOS
        sev00 = sev
        ef00  = ef0
        if (ldos .ne. 0) then
          call efldos ( qval , nsp , emin , emax , ndos , dos_rv , eferm 
     .    , sev1 ) 

          sev   = sev1
          ef0  = eferm
        endif
        if (ipr .gt. 30 .and. ldos .ne. 0)
     .  write(stdo,388) sev00,ef00,sev1,eferm,sev,ef0
  388   format(' ipol:  sev=',f12.6,'   ef=',f12.6:/
     .  ' dos:   sev=',f12.6,'   ef=',f12.6/
     .  ' use:   sev=',f12.6,'   ef=',f12.6)
        if (ipl.gt.0) write (stdl,733) ef00,eferm,ef0,sev00,sev1,sev
  733   format('nf EF:',3f9.5,'    EB:',3f12.5)

        sbz%ndos=ndos
        i_copy_size=size(sbz%dosw)
        call dcopy(i_copy_size,dosw,1,sbz%dosw,1)
        sbz%ef=ef0
        sbz%def=def

        if (lrep .eq. 1) then
          ef0 = -1
          call awrit0('Input Fermi energy was too far off, '//
     .    'repeat band pass',' ',80,stdo)
          goto 99
        endif
      endif

      if( debug) print *, 'tttttttttt 9999999999999999999 6',lmet,ltet,lwndow
C --- BZ integration for fermi level, band sum and qp weights ---
      if (lmet .ge. 0 .and. (lmet .ne. 4 .or. ltet)) then
        if (lwndow) then
          allocate(ww_rv(ndham*nsp*nkp))

          eferm = min(dosw(1),dosw(2))
           call bzints ( nkabc ( 1 ) , nkabc ( 2 ) , nkabc ( 3 ) , evlall 
     .     , ww_rv , nkp , ndham , ndham , nsp , xv , xv , xv , 1 , eferm 
ckino Dec.16.2011:       .     , 2 , ntet , iv_p_oidtet , sev , dum ) 
ckino Dec.28.2011:       .     , 2 , ntet , sbz%iv_p_oidtet , sev , dum ) 
     .     , 2 , ntet , sbz%iv_a_oidtet , sev , dum ) 


C         call prmx('w(min)',w(oww),ndham,ndham,nkp*nsp)
          eferm = max(dosw(1),dosw(2))
            call bzints ( nkabc ( 1 ) , nkabc ( 2 ) , nkabc ( 3 ) , evlall 
ckino Dec.8.2011:       .      , rv_p_owtkb , nkp , ndham , ndham , nsp , xv , xv , xv , 1 , 
     .      , rv_a_owtkb , nkp , ndham , ndham , nsp , xv , xv , xv , 1 , 
ckino Dec.16.2011:       .      eferm , 2 , ntet , iv_p_oidtet , sev , dum ) 
ckino Dec.28.2011:       .      eferm , 2 , ntet , sbz%iv_p_oidtet , sev , dum ) 
     .      eferm , 2 , ntet , sbz%iv_a_oidtet , sev , dum ) 

C         call prmx('w(max)',w(owtkb),ndham,ndham,nkp*nsp)
ckino Dec.8.2011:             call daxpy ( ndham * nsp * nkp , - 1d0 , ww_rv , 1 , rv_p_owtkb, 1 ) 
           call daxpy ( ndham * nsp * nkp , - 1d0 , ww_rv , 1 , rv_a_owtkb, 1 ) 

C         call prmx('w',w(owtkb),ndham,ndham,nkp*nsp)
          if (allocated(ww_rv)) deallocate(ww_rv)

        else
          dosrng = 8
          if (mpsord .lt. 0) dosrng = 16
              if( debug) print *, 'tttttttttt 99999999 6 call bzwtsf'

              if(sbz%fsmommethod == 1) then !takao dec2010
               call bzwtsf2 ( ndham , ndham , nsp , nspc , nkabc ( 1 ) , nkabc 
ckino Dec.16.2011:       .         ( 2 ) , nkabc ( 3 ) , nkp , ntet , iv_p_oidtet , qval - qbg , 
ckino Dec.28.2011:       .         ( 2 ) , nkabc ( 3 ) , nkp , ntet , sbz%iv_p_oidtet , qval - qbg , 
     .         ( 2 ) , nkabc ( 3 ) , nkp , ntet , sbz%iv_a_oidtet , qval - qbg , 
     .         fsmom , lmet.ne.0 , ltet , mpsord , ndos , dabs ( esmear - mpsord 
ckino Dec.8.2011:       .         ) , dosrng , rv_p_owtkp , evlall , efmax , lswtk , rv_p_oswtk 
ckino Dec.8.2011:       .         ) , dosrng , sbz%rv_p_owtkp , evlall , efmax , lswtk , rv_p_oswtk 
ckino Dec.27.2011:       .         ) , dosrng , sbz%rv_p_owtkp , evlall , efmax , lswtk , rv_a_oswtk 
     .         ) , dosrng , sbz%rv_a_owtkp , evlall , efmax , lswtk , rv_a_oswtk 
ckino Dec.8.2011:       .         , eferm , sev , rv_p_owtkb , sumqv ( 1 , 2 ) , lwtkb ) 
     .         , eferm , sev , rv_a_owtkb , sumqv ( 1 , 2 ) , lwtkb ) 
              else
ckino       maybe a switch for these is lmet.ne.0 
ckino Dec.8.2011:                 l_temp_rv_p_oswtk=.false.
               l_temp_rv_a_oswtk=.false.
ckino change rv_p_oswtk allocatable               if (.not.associated(rv_p_oswtk)) then
ckino Dec.8.2011:                 if (.not.allocated(rv_p_oswtk)) then
               if (.not.allocated(rv_a_oswtk)) then
ckino Dec.8.2011:                     allocate(rv_p_oswtk(1))
                   allocate(rv_a_oswtk(1))
ckino Dec.8.2011:                     l_temp_rv_p_oswtk=.true.
                   l_temp_rv_a_oswtk=.true.
               endif
ckino Dec.8.2011:                 l_temp_rv_p_owtkb=.false.
               l_temp_rv_a_owtkb=.false.
ckino bug fix oswtk -> owtkb
ckino change rv_p_owtkb allocatable
ckino               if (.not.associated(rv_p_oswtk))then
ckino                   allocate(rv_p_oswtk(1))
ckino Dec.8.2011:                 if (.not.allocated(rv_p_owtkb))then
               if (.not.allocated(rv_a_owtkb))then
ckino Dec.8.2011:                     allocate(rv_p_owtkb(1,1,1))
                   allocate(rv_a_owtkb(1,1,1))
ckino Dec.8.2011:                     l_temp_rv_p_owtkb=.true.
                   l_temp_rv_a_owtkb=.true.
               endif
ckino Dec.28.2011:                 l_temp_iv_p_oidtet=.false.
               l_temp_iv_a_oidtet=.false.
ckino Dec.16.2011:                 if (.not.associated(iv_p_oidtet))then
ckino Dec.28.2011:                 if (.not.associated(sbz%iv_p_oidtet))then
               if (.not.allocated(sbz%iv_a_oidtet))then
ckino Dec.16.2011:                     allocate(iv_p_oidtet(1))
ckino Dec.28.2011:                     allocate(sbz%iv_p_oidtet(1))
                   allocate(sbz%iv_a_oidtet(1))
ckino Dec.28.2011:                     l_temp_iv_p_oidtet=.true.
                   l_temp_iv_a_oidtet=.true.
               endif
ckino end
               call bzwtsf ( ndham , ndham , nsp , nspc , nkabc ( 1 ) , nkabc 
ckino Dec.16.2011:       .         ( 2 ) , nkabc ( 3 ) , nkp , ntet , iv_p_oidtet , qval - qbg , 
ckino Dec.28.2011:       .         ( 2 ) , nkabc ( 3 ) , nkp , ntet , sbz%iv_p_oidtet , qval - qbg , 
     .         ( 2 ) , nkabc ( 3 ) , nkp , ntet , sbz%iv_a_oidtet , qval - qbg , 
     .         fsmom , lmet.ne.0 , ltet , mpsord , ndos , dabs ( esmear - mpsord 
ckino Dec.8.2011:       .         ) , dosrng , rv_p_owtkp , evlall , efmax , lswtk , rv_p_oswtk 
ckino Dec.8.2011:       .         ) , dosrng , sbz%rv_p_owtkp , evlall , efmax , lswtk , rv_p_oswtk 
ckino Dec.27.2011:       .         ) , dosrng , sbz%rv_p_owtkp , evlall , efmax , lswtk , rv_a_oswtk 
     .         ) , dosrng , sbz%rv_a_owtkp , evlall , efmax , lswtk , rv_a_oswtk 
ckino Dec.8.2011:       .         , eferm , sev , rv_p_owtkb , sumqv ( 1 , 2 ) , lwtkb ) 
     .         , eferm , sev , rv_a_owtkb , sumqv ( 1 , 2 ) , lwtkb ) 
              endif
ckino
ckino Dec.8.2011:                if (l_temp_rv_p_oswtk) then
              if (l_temp_rv_a_oswtk) then
ckino change rv_p_oswtk allocatable                 deallocate(rv_p_oswtk); rv_p_oswtk=>NULL()
ckino Dec.8.2011:                   deallocate(rv_p_oswtk)
                 deallocate(rv_a_oswtk)
ckino add .false
ckino Dec.8.2011:                   l_temp_rv_p_oswtk=.false.
                 l_temp_rv_a_oswtk=.false.
              endif
ckino Dec.8.2011:                if (l_temp_rv_p_owtkb) then
              if (l_temp_rv_a_owtkb) then
ckino change rv_p_owtkb allocatable                 deallocate(rv_p_owtkb); rv_p_owtkb=>NULL()
ckino Dec.8.2011:                   deallocate(rv_p_owtkb)
                 deallocate(rv_a_owtkb)
ckino add .false
ckino Dec.8.2011:                   l_temp_rv_p_owtkb=.false.
                 l_temp_rv_a_owtkb=.false.
              endif
ckino Dec.28.2011:                if (l_temp_iv_p_oidtet) then
              if (l_temp_iv_a_oidtet) then
ckino Dec.16.2011:                   deallocate(iv_p_oidtet); iv_p_oidtet=>NULL()
ckino Dec.28.2011:                   deallocate(sbz%iv_p_oidtet); sbz%iv_p_oidtet=>NULL()
                 deallocate(sbz%iv_a_oidtet); 
ckino add .false
ckino Dec.28.2011:                   l_temp_iv_p_oidtet=.false.
                 l_temp_iv_a_oidtet=.false.
              endif
ckinoend
              if( debug) print *, 'tttttttttt 99999999 6 call bzwtsf'

C         Store val charge & magnetic moment in sumqv(1..2)
          if (lmet .ne. 4) then
            sumqv(1,1) = sumqv(1,2)
            sumqv(2,1) = sumqv(2,2)
          endif
!!=== takao automatic search for nevmx (this can be just outside of endif below? I failed why?) ===
!!Nov25 2010
ckino change rv_p_owtkb_allocatable         if(associated(rv_p_owtkb)) then
ckino Dec.8.2011:            if(allocated(rv_p_owtkb)) then
          if(allocated(rv_a_owtkb)) then
            do  i = ndham, 1, -1 !this is the same as nvec search in addrbl.
            nevec = i
ckino Dec.8.2011:              if (maxval(rv_p_owtkb(i,:,:)) .gt. epsnevec()) exit
            if (maxval(rv_a_owtkb(i,:,:)) .gt. epsnevec()) exit
            enddo
c          call getnevec(ndham,nsp,nkp,nev,rv_p_owtkb, nevec)
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
ckino Dec.8.2011:               i = iobzwt ( 0 , ndhamx , nkp , nspx , eferm , rv_p_owtkb , -ifi ) 
             i = iobzwt ( 0 , ndhamx , nkp , nspx , eferm , rv_a_owtkb , -ifi ) 
            call fclr('wkp',ifi)
          endif
        endif
        if (lwtkb .eq. -1 .and. lrout .gt. 0) then
          call info0(20,0,0,' Start second band pass ...')
          lwtkb = 1
          if (nspc .eq. 2) lswtk = 1
          goto 99
        endif
        if (lwtkb .eq. 2 .and. lrout .gt. 0) then
          call info0(20,0,0,' New pass with constrained weights ...')
          goto 99
        endif
      endif

C ... release sigma-related things.
      if(sigmamode) then
        deallocate(sfz,qsmesh,npair,pos,nlat,qwgt)
      endif

      if( debug) print *, 'tttttttttt 9999999999999999999 7'

C ... Save Fermi level, nonmetal or sampling integration
      if (lmet .eq. 0 .or. .not. (lmet .ne. 4 .or. ltet)) then
        if (lmet .eq. 0) then
          ef0 = (evtop+ecbot)/2
          if (ipr .ge. 20) call awrit3(' Highest occ. level = %,5;5d '//
     .    ' Lowest unocc. = %,5;5d  average = %,5;5d',' ',80,
     .    stdo,evtop,ecbot,ef0)
        endif
        if (procid .eq. master) then
          ifi = fopna('wkp',-1,4)
ckino Dec.8.2011:             i = iobzwt ( 1 , ndham , nkp , nsp , ef0 , rv_p_owtkb , -ifi ) 
           i = iobzwt ( 1 , ndham , nkp , nsp , ef0 , rv_a_owtkb , -ifi ) 
          call fclr('wkp',ifi)
        endif
      endif

C ... Cleanup asa-style moments file, print table of DOS channels
      if (nfilem .gt. 0) then
        i = iomoms(-nfilem,nl,nsp,nspc,nkp,ndimh,11,1,nkp*nsp+1,1,
     .  ndham,ndimh,nchan,nchan,ndimh,0d0,0d0,0d0,0d0,eferm,0d0)
        call fclose(nfilem)
        if (iprint() .ge. 10 .and. cmdopt('--mull',6,0,strn)) then
          call mchan ( lmdim , ssite , sspec , nsp , nsites , lsites , 
     .    0 , 0 , 0 , 0 , chan_iv ) 

        endif
      endif
      if (allocated(chan_iv)) deallocate(chan_iv)

C ... Average forces so net force on system is zero (APW case)
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

C --- Generate DOS on disk ---
      if (.not. lwndow) then

        ndos0=sbz%ndos

        if ( iand(1,int(sctrl%ldos)) .ne.0 ) then

c          allocate(dos_rv(3*ndos))

          if (iprint()>=30) call awrit1('%x%N ... Generating %?#n<0#'
     .    //'integrated#total# DOS',' ',80,lgunit(1),ndos0)
          if (ltet) then
             call bzints ( nkabc ( 1 ) , nkabc ( 2 ) , nkabc ( 3 ) , evlall 
     .       , dum , nkp , ndhamx , ndhamx , nspx , dosw ( 1 ) , dosw ( 
ckino Dec.16.2011:       .       2 ) , dos_rv , ndos , eferm , 1 , ntet , iv_p_oidtet , dum , 
ckino Dec.28.2011:       .       2 ) , dos_rv , ndos , eferm , 1 , ntet , sbz%iv_p_oidtet , dum , 
     .       2 ) , dos_rv , ndos , eferm , 1 , ntet , sbz%iv_a_oidtet , dum , 
     .       dum ) 


            if ( ndos0 .gt. 0 ) call xxxdif ( dosw ( 1 ) , dosw ( 2 ) , ndos 
     .      , nspx , 0 , dos_rv ) 

            del = 0d0
          else
            if (mpsord .ge. 100) mpsord = mod(mpsord,100)
             if ( ndos0 .gt. 0 ) call makdos ( nkp , ndhamx , ndhamx , nspx 
ckino Dec.8.2011:       .       , rv_p_owtkp , evlall , mpsord , dabs ( esmear - mpsord ) 
ckino Dec.27.2011:       .       , sbz%rv_p_owtkp , evlall , mpsord , dabs ( esmear - mpsord ) 
     .       , sbz%rv_a_owtkp , evlall , mpsord , dabs ( esmear - mpsord ) 
     .       , - 6d0 , dosw ( 1 ) , dosw ( 2 ) , ndos , dos_rv ) 


             if ( ndos0 .lt. 0 ) call maknos ( nkp , ndhamx , ndhamx , nspx 
ckino Dec.8.2011:       .       , rv_p_owtkp , evlall , mpsord , dabs ( esmear - mpsord ) 
ckino Dec.27.2011:       .       , sbz%rv_p_owtkp , evlall , mpsord , dabs ( esmear - mpsord ) 
     .       , sbz%rv_a_owtkp , evlall , mpsord , dabs ( esmear - mpsord ) 
     .       , - 6d0 , dosw ( 1 ) , dosw ( 2 ) , ndos , dos_rv ) 


            del = mpsord+dabs(esmear-mpsord)
          endif
          if ( nspc .eq. 2 ) call dscal ( ndos , .5d0 , dos_rv , 1 ) 

          call dosio ( dos_rv , ndos , nspx , ndos , 1 , dosw ( 1 ) 
     .    , dosw ( 2 ) , nspx , eferm , del , 1 , - fopn ( 'DOS' ) ) 

          call fclose(fopn('DOS'))
          if (allocated(orbtm_rv)) deallocate(orbtm_rv)
c           if (allocated(ausc_zv)) deallocate(ausc_zv)
c           if (allocated(fes2_rv)) deallocate(fes2_rv)
c           if (allocated(fh_rv)) deallocate(fh_rv)
c           if (allocated(srout_zv)) deallocate(srout_zv)
          if (allocated(dos_rv)) deallocate(dos_rv)
c           call rlse(oqp) --->wrong rlse. Break some arrays

        endif
      endif


      if( debug) print *, 'tttttttttt 9999999999999999999 7'


C --- Core-level spectroscopy ---
      if (icls .ne. 0) then
        if (procid .eq. master) then

          eferm=sbz%ef

          call vcdmel ( sctrl , ssite , sspec , slat , nlmax , ndham , 
     .    ndimh , nkp , nsp , nspc , eferm , evlall , ausc_zv , 
     .    nsites , isite , iclsl , iclsn ) 

          call rx0('done generating core level spectra')
        else
          call fexit(0,0,' ',0)
        endif
      endif

C --- Assemble output density ---
       call dfratm ( nbas , ssite , sspec , sv_p_orhat1 ) 

      allocate(qbyl_rv(n0*nbas*nsp))
      allocate(hbyl_rv(n0*nbas*nsp))
#ifdef KINODEBUG
      write(*,*) 'kino: allocate hbyl_rv,',n0*nbas*nsp
#endif
C     --window: Put output density into orhat and smrho, and exit
      if (lwndow) then
           call mkrout ( ssite , sspec , slat , sham , nbas , nsp , ldim 
ckino Dec.8.2011:       .     , lekkl , sv_p_oqkkl , sv_p_oeqkkl , sv_p_oorhat , hab_rv , sab_rv 
ckino Dec.28.2011:       .     , lekkl , sv_p_oqkkl , sv_p_oeqkkl , spot%sv_p_oorhat , hab_rv , sab_rv 
     .     , lekkl , sv_p_oqkkl , sv_p_oeqkkl , spot%sv_a_oorhat , hab_rv , sab_rv 
     .     , qbyl_rv , hbyl_rv , lrout ) 




ckino Dec.8.2011:           call zcopy ( k1 * k2 * k3 * nsp , srout_zv , 1 , zv_p_osmrho 
ckino Dec.28.2011:           call zcopy ( k1 * k2 * k3 * nsp , srout_zv , 1 , spot%zv_p_osmrho 
         call zcopy ( k1 * k2 * k3 * nsp , srout_zv , 1 , spot%zv_a_osmrho 
     .   , 1 ) 


ckino Dec.8.2011:            call symrho ( ssite , sspec , slat , lfrce , zv_p_osmrho , sv_p_oorhat 
ckino Dec.8.2011:            call symrho ( ssite , sspec , slat , lfrce , zv_p_osmrho , spot%sv_p_oorhat 
ckino Dec.28.2011:            call symrho ( ssite , sspec , slat , lfrce , spot%zv_p_osmrho , spot%sv_p_oorhat 
ckino Dec.28.2011:            call symrho ( ssite , sspec , slat , lfrce , spot%zv_a_osmrho , spot%sv_p_oorhat 
          call symrho ( ssite , sspec , slat , lfrce , spot%zv_a_osmrho , spot%sv_a_oorhat 
     .    , qbyl_rv , hbyl_rv , frc ) 



        goto 999
      endif
      call mkrout ( ssite , sspec , slat , sham , nbas , nsp , ldim 
     ., lekkl , sv_p_oqkkl , sv_p_oeqkkl , sv_p_orhat1 , hab_rv , sab_rv 
     ., qbyl_rv , hbyl_rv , lrout ) 


      if( debug) print *, 'tttttttttt 9999999999999999999 7'


C
      if (lrout .ne. 0) then
C   ... Symmetrize output density and forces
         call symrho ( ssite , sspec , slat , lfrce , srout_zv , sv_p_orhat1 
     .   , qbyl_rv , hbyl_rv , frc ) 



C   ... New boundary conditions pnu for phi and phidot
C       call pshpr(iprint()-12)
        if (lpnu .gt. 0) then

          i_copy_size=size(sham%pmin)
          call dcopy(i_copy_size,sham%pmin,1,xv,1)


          i_copy_size=size(sham%pmax)
          call dcopy(i_copy_size,sham%pmax,1,xv(11),1)

          call pnunew ( nbas , nsp , ssite , sspec , xv , xv ( 11 ) , lfrzw 
     .    , hab_rv , sab_rv , qbyl_rv , hbyl_rv, eferm ) 

        endif
C       call poppr
      endif


      if( debug) print *, 'tttttttttt 9999999999999999999 8'

C --- Evaluate Harris energy ---
c      call pshpr(ipr-10) !this was needed to go through test.fp (just for check write.) Why? --> but no problem... sep2009
c      print *,'vvv takao goto mkehfk mode1'
      call mkehkf(1,sham,sev,sumqv(2,1),sumtv,ehar)
c      call poppr()

C --- Evaluate KS total energy, correction to Harris force ---
      if (lrout .ne. 0) then

C   ... Correction to harris force
           call dfrce ( ssite , sspec , slat , sctrl , k1 , k2 , k3 , nvl 
ckino Dec.8.2011:       .     , sv_p_oorhat , sv_p_orhat1 , elind , qmom_rv , zv_p_osmrho , 
ckino Dec.8.2011:       .     , spot%sv_p_oorhat , sv_p_orhat1 , elind , qmom_rv , zv_p_osmrho , 
ckino Dec.28.2011:       .     , spot%sv_p_oorhat , sv_p_orhat1 , elind , qmom_rv , spot%zv_p_osmrho , 
ckino Dec.28.2011:       .     , spot%sv_p_oorhat , sv_p_orhat1 , elind , qmom_rv , spot%zv_a_osmrho , 
     .     , spot%sv_a_oorhat , sv_p_orhat1 , elind , qmom_rv , spot%zv_a_osmrho , 
     .     srout_zv , fh_rv ) 





C   ... Evaluate KS total energy and output magnetic moment
        eks = 0d0
        if (leks .ge. 1) then
          call togpr()
c          print *,'vvv takao mkekin 000 --- oldpot x newrho'
               call mkekin ( nbas , ldim , ssite , sspec , slat , sham , sv_p_osig 
     .         , sv_p_otau , sv_p_oppi , lcplxp , sv_p_oqkkl , k1 , k2 , k3 
ckino Dec.8.2011:       .         , vconst , zv_p_osmpot , srout_zv , sev , sumtv ) 
ckino Dec.28.2011:       .         , vconst , spot%zv_p_osmpot , srout_zv , sev , sumtv ) 
     .         , vconst , spot%zv_a_osmpot , srout_zv , sev , sumtv ) 







ccccccccccccccccccccccccc
ctakao test mkekin repeat
c         print *,'vvv takao mkekin repeat 000 --- oldpot x oldrho'
c         call mkekin(nbas,ldim,ssite,sspec,slat,sham,osig,otau,oppi,
c     .      lcplxp,oqkkl,k1,k2,k3,vconst,w(osmpot),w(osmrho),sevxxx,sumtvxxx)
c       stop 'xxxxxxxxxxxxxxxx 111aaa'
ccccccccccccccccccccccccc

ccccccccccccccccccccc
c          print *,'vvv: takao test osmpotbk 000 '
c          call defdc(osmpotbk,k1*k2*k3*2)
c          call dcopy(k1*k2*k3*4,w(osmpot),1,w(osmpotbk),1)
ccccccccccccccccccccc

ctakao sep2009
c          call pshpr(ipr-4)
Ckino mkpot calls locpot. and locpot calls augmat. 
Ckino augmat calculates sig,tau,ppi.
               call mkpot ( nbas , ssite , sspec , slat , sham , spot , lfrce 
ckino Dec.8.2011:       .         , lcplxp , k1 , k2 , k3 , srout_zv , sv_p_orhat1 , qbg , zv_p_osmpot 
ckino Dec.28.2011:       .         , lcplxp , k1 , k2 , k3 , srout_zv , sv_p_orhat1 , qbg , spot%zv_p_osmpot 
     .         , lcplxp , k1 , k2 , k3 , srout_zv , sv_p_orhat1 , qbg , spot%zv_a_osmpot 
     .         , qmom_rv , vconst , vrmt , sv_p_osig , sv_p_otau , sv_p_oppi 
     .         , ppnl_rv , hab_rv , vab_rv , sab_rv , qval , qsc , gpot0_rv 
     .         , vval_rv , fes2_rv , 0 , vorb , nlibu , lmaxu , lldau , nsp,0 ) 






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
c          call poppr()
c          call pshpr(ipr-20)
c          print *,'vvv takao goto mkehfk mode2'
          call mkehkf(2,sham,sev,sumqv(2,1),sumtv,eks)
          call togpr()
        endif

      if( debug) print *, 'tttttttttt 9999999999999999999 9'

C   --- Add together force terms ---
        if ( lfrce .gt. 0 ) call totfrc ( nbas , ssite , slat , leks 
     .  , fes1_rv , fes2_rv , fh_rv , frc ) 


C   --- Mix input and output densities ---
        allocate(qm1_rv(nvl))

        allocate(qm2_rv(nvl))
#ifdef KINODEBUG
         write(*,*) 'kino allocate qm2_rv',nvl
#endif

        call upacks('strn mix',i1,i2)
#if MPE
        ierr = MPE_LOG_EVENT(EVENT_START_MIXRHO,procid,"mixrho")
#endif
cccccccccccccccccccccc
cc          print *,'vvv takao mkekin copy test xxx1---'
cc          allocate(smrhobk(k1*k2*k3*2))
cc1          call dcopy(k1*k2*k3*4,w(osmrho),1,smrhobk,1)
cc2          call dcopy(k1*k2*k3*2,srout_zv,1,smrhobk,1)
ccccccccccccccccccccccc
           call mixrho ( ssite , sspec , slat , nsp , iter , sstrn ( i1:i2 
ckino Dec.8.2011:       .     ) , qval - qbg , elind , sv_p_orhat1 , sv_p_oorhat , k1 , k2 
ckino Dec.28.2011:       .     ) , qval - qbg , elind , sv_p_orhat1 , spot%sv_p_oorhat , k1 , k2 
     .     ) , qval - qbg , elind , sv_p_orhat1 , spot%sv_a_oorhat , k1 , k2 
ckino Dec.8.2011:       .     , k3 , dmxp , srout_zv , zv_p_osmrho ) 
ckino Dec.28.2011:       .     , k3 , dmxp , srout_zv , spot%zv_p_osmrho ) 
     .     , k3 , dmxp , srout_zv , spot%zv_a_osmrho ) 
ccccccccccccccc


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

#if MPE
        ierr = MPE_LOG_EVENT(EVENT_END_MIXRHO,procid,"mixrho")
#endif
      else
        eks = 0
      endif

      sham%ehf=ehar
      sham%ehk=eks


  999 continue
      if (allocated(qm2_rv)) deallocate(qm2_rv)
      if (allocated(qm1_rv)) deallocate(qm1_rv)
      if (allocated(hbyl_rv)) deallocate(hbyl_rv)
      if (allocated(qbyl_rv)) deallocate(qbyl_rv)
c      if (allocated(eomin_rv)) deallocate(eomin_rv)
      if (allocated(fes1_rv)) deallocate(fes1_rv)
      if (allocated(ppnl_rv)) deallocate(ppnl_rv)
      if (allocated(sab_rv)) deallocate(sab_rv)
      if (allocated(vab_rv)) deallocate(vab_rv)
      if (allocated(hab_rv)) deallocate(hab_rv)
      if (allocated(vval_rv)) deallocate(vval_rv)
      if (allocated(gpot0_rv)) deallocate(gpot0_rv)
      if (allocated(qmom_rv)) deallocate(qmom_rv)
c       call rlse(oqp) !release all the w(oxxx) array allocated within this routine.
      call tcx('bndfp')
ctakaox june2009
c      close(ifiogw)
Cki
      if ( ifiese>0 ) close(ifiese)
c      call xxxbfp ! Patch to avoid PGI compiler bug on AMD processor
      if (allocated(ifbls_iv)) deallocate(ifbls_iv)
#if MPIK
cki       if (associated(rv_p_oqp)) deallocate(rv_p_oqp)

#endif
c      call rlse(oevlall)
      deallocate(evlall)
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
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   10 Jul 06
C ----------------------------------------------------------------------
      implicit none
      include "mpif.h"
C ... Passed parameters
      integer kpproc(0:*),ndham,nkp,nsp
      double precision eb(ndham,*)
C ... Local parameters
      integer i,ista,iend
      integer procid,master
      integer numprocs, ierr
      integer, dimension(:),allocatable :: offset,length
      real(8) ,allocatable :: buf_rv(:)

C ... Heap

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

C      if (procid .eq. master) then
C        print *, procid,eb(1,1),eb(1,nsp*nkp)
C      endif
C      call rx('done')
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

