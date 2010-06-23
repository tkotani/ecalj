Cchp1 cmdvar= integer,iv_p_,oiprmb
Cchp1 cmdvar= integer,iv_p_,oidxsh
Cchp1 cmdvar= integer,iv_p_,oindxo
Cchp1 cmdvar= real(8),rv_p_,owtkp
Cchp1 cmdvar= real(8),rv_p_,oqsig
Cchp1 cmdvar= real(8),rv_p_,oqp
      subroutine bndfp(nbas,nsp,nlibu,lmaxu,lldau,
     .ssite,sspec,slat,sctrl,sham,spot,sbz,sstrn,
     .ndham,leks,lrout,lfrce,lpnu,dmxp,iter,maxit,evl,frc,
     .dmatu,vorb,llmfgw)
      use m_hamindex, only: 
     &ngrp, norbt, napwmx, ltab,ktab,offl,ispec,iclasst,invgx,miat,offlrev,ibastab,
     &igv2,napwk,igv2rev,ibasindex,getikt
      use m_rdctrl, only: ncutovl
      use m_struc_def

#if MPI | MPIK
      use m_events !Event numbers: include "events.ins"
#endif
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
      real(8):: dmxp(33),evl(ndham,2),frc(3,nbas)
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
      real(8) ,allocatable :: eomin_rv(:)

      integer:: ipr , iprint , ipl , isp , jsp , ispc , isum , nglob 
     ., i , iopq , iq , isqp , ismidb , nsmidb , ldos , lmet , nspc 
     ., lrep , lcplxp , lso , nkabc(3) , ndos , ndos0 , nev , nevl 
     ., nkp , nsp , numq , n1 , n2 , n3 , ntet , lwtkb , lswtk , nl 
     ., lgunit , igets , mpsord , nevmx , nvl , ngabc(3) , isw , stdo 
     ., stdl , ifi , fxst , fopn , fopna , lfrzw , i1 , i2 , plbnd 
     ., nfilem , iobzwt , lnoxc , lrsig , lwsig , jobgw , nll , lpdiag 
     ., iv(10) , parg , ndhamx , nspx , nk1 , nk2 , nk3 , lshft(3) 
     ., onesp , nfbn(2) , mpipid , ldham(8,2) , ldim , ndimh , ndimhx 
     ., lekkl 
      real(8) ,allocatable :: dos_rv(:)

Cchp1       integer:: oevl , os , oidtet , owtkb , oswtk , owtkp , oorhat 
Cchp1      ., osoptc , osmpot , osmrho , oqp 
Cchp1        integer:: oevl , os , oidtet , owtkb , oswtk , owtkp , oorhat 
Cchp1      . , osoptc , osmpot , osmrho 
        integer:: oevl , os , oidtet , owtkb , oswtk , oorhat , osoptc 
     .  , osmpot , osmrho 
        real(8),pointer :: rv_p_owtkp(:)

       real(8),pointer :: rv_p_oqp(:)

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

      integer osig(3,nbx),otau(3,nbx),oppi(3,nbx),oqkkl(3,nbx),
     .oeqkkl(3,nbx),orhat1(3,nbx)
      equivalence (ldim,ldham(1,1))

      logical a2bin
      integer osigx(3,nbx),otaux(3,nbx),oppix(3,nbx),ospotx,osgw !osmrox,

      double precision sigp(10)

      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      !equivalence (nk1,nkabc(1)), (nk2,nkabc(2)), (nk3,nkabc(3))
      double precision ef00,eferm,emax,emin,esmear,qval,qsc,sev,
     .sev00,sev1,sumtv,alfa,vconst,qbg,ebot,pi,def,ef0,ehar,eks,
     .dosw(2),efmax,dgets,alat,plat(3,3),qlat(3,3),sumev(2,3),
     .sumqv(3,2),qp(3),shftqp(3),eterms(20),elind,vol,dum,del,evtop,
     .ecbot,vrmt(nbx),fsmom,rsrnge,xv(20),dosrng,epsovl,eomin,dval
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
Cchp1       integer:: oiprmb , moddos , nsites , lsites(nbx) , nchan , ng 
Cchp1      ., iomoms , nchmx , lmdim , lmxch 
       integer:: moddos , nsites , lsites(nbx) , nchan , ng , iomoms 
     . , nchmx , lmdim , lmxch 
       integer,pointer :: iv_p_oiprmb(:)

      real(8) ,allocatable :: doswt_rv(:)
      integer ,allocatable :: chan_iv(:),chan0_iv(:)

C     for pzhev
      integer nblk,nprow,npcol,nmx
C     for LDA+U
      real(8) ,allocatable :: ausu_rv(:)

      double complex vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu),
     .dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)

      integer:: i_copy_size, i_spackv, i_spacks
C ... Heap
      integer w(1)
      common /w/ w

      logical :: debug=.false. !debug is only for check write. However, it can cause error
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
Cchp1       integer oqsig
       real(8),pointer :: rv_p_oqsig(:)


Cchp1       integer oidxsh,nbasp
       integer:: nbasp 
       integer,pointer :: iv_p_oidxsh(:)


C     Data statements
      data T /.true./ F /.false./

      integer,allocatable:: nevls(:,:)

ctakao
      integer :: ifiogw,ifiese,oag,og,ib,ig,iorb,is,nsgrp,uspec_size,nspec
      logical :: nexist,sigmamode,ok, newsigmasw
      character spid*8
      integer :: ib1,ib2,ii,ik1,ik2,ik3,j1,j2
      integer :: nhmax,nknknk,npairmx,opos,jx,nsp_dummy,ndimh_dummy
      complex(8),allocatable:: sfz(:,:,:,:,:,:),wgtq(:,:,:,:,:),sfzx(:,:,:,:,:,:)
      complex(8),allocatable:: hamm(:,:,:),ovlm(:,:,:),sene(:,:) !,hamsom(:,:)
      real(8):: rrr(3)
      complex(8):: pi2i
      integer,allocatable:: npair(:,:), nlat(:,:,:,:)
      real(8),allocatable:: qwgt(:,:,:),range(:,:),qsmesh(:,:,:,:),pos(:,:)


      integer:: ikc(3,0:1,0:1,0:1), igshift(3,0:1,0:1,0:1),j1x,j2x,ix1,ix2,ix3
      integer,allocatable::igx(:,:,:,:)
      real(8):: wgc(0:1,0:1,0:1),qpara(3)

      complex(8)::sss
      real(8):: www,qmm(3),qqq(3),sumtvxxx,sevxxx
      integer:: iktm,ikt,igvecc(3),nkk(3),igvecc0(3), nnspc,osmpotbk
      logical:: noshorbz,siginit
      real(8),allocatable:: qplist(:,:)
      real(8):: qpo(3),qpx(3)

      integer:: iqini,iqend
      complex(8),allocatable:: spotx(:,:,:,:)

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



C --- Setup ---
      call tcn ('bndfp')
      siginit=.true.
      napw = 0
      ipr  = iprint()
      ipl  = ipr
      nsp  = nglob('nsp')
      nspc = nglob('nspc')
C     npqn = nglob('npqn')
      lso = isw ( iand(4,int(sctrl%lncol)) .ne.0 ) + 2 * isw ( iand(32,int(sctrl%lncol)) .ne.0 
     .) + 3 * isw ( iand(64,int(sctrl%lncol)) .ne.0 )

      lcplxp = 0
      lekkl = int(sctrl%pfloat)

      if (lso .ne. 0) lcplxp = 1
      if (isum(nbas,lldau,1) .ne. 0) lcplxp = 1
      ndhamx = ndham*nspc
      nspx = nsp / nspc

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

ctakao's test for replacement.
      ldos = iand(1,int(sctrl%ldos))

      lrsig = int(sham%lsig)

      lfrzw = isw ( iand(16,int(sctrl%lbas)) .ne.0 )

c      print *,'xxx' ,lfrzw,sctrl%lbas
c      print *,'ttt1: lrsig,ldos,lfrzw=',lrsig,ldos,lfrzw
      print *,'ttt: lrsig,ldos,lfrzw=',lrsig,ldos,lfrzw
      print *,'ttt: sham%lsig       =',sham%lsig
      if(ldos - mod(sctrl%ldos/1, 2)/=0) stop 'ttt xxx111'
      if(lfrzw- mod(sctrl%lbas/16,2)/=0) stop 'ttt xxx222'
      if(lrsig- sham%lsig/=0) stop 'ttt xxx333'
c$$$ctakao test
      nspec= sctrl%nspec
ctakao
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


      oorhat=spot%oorhat
      osmrho=spot%osmrho
      osmpot=spot%osmpot
c       osoptc=spot%osoptc

C     icls = igets('optic cls',w(osoptc))
      pi = 4d0*datan(1d0)
      if (iprint() .ge. 20) call awrit2('%N --- BNDFP:  '//
     .'begin iteration %i of %i ---',' ',80,stdo,iter,maxit)
      if (nbas.gt.nbx) call rxi('bndfp: nbx exceeded, need',nbas)
      if (ndham.le.0)
     .call rx('bndfp: hamiltonian matrix has zero dimension')


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

      alat=slat%alat
      i_copy_size=size(slat%plat)
      call dcopy(i_copy_size,slat%plat,1,plat,1)
      i_copy_size=size(slat%qlat)
      call dcopy(i_copy_size,slat%qlat,1,qlat,1)
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
      oidtet=sbz%oidtet
      lmet=sbz%lmet


      mpsord=sbz%n
      esmear=sbz%w
      nevmx=sbz%nevmx
      efmax=sbz%efmax
      fsmom=sbz%fsmom

c      print *,' w=',esmear
      esmear = esmear + iabs(mpsord)
      if (mpsord .lt. 0) esmear = -esmear
      mpsord = int(dabs(esmear))
      if (esmear .lt. 0) mpsord = -mpsord


      ndos=sbz%ndos
      i_copy_size=size(sbz%dosw)
      call dcopy(i_copy_size,sbz%dosw,1,dosw,1)
      ef0=sbz%ef
      def=sbz%def


Cchp1       oqp=sbz%oqp
       rv_p_oqp => sbz%rv_p_oqp 

Cchp1       owtkp=sbz%owtkp
       rv_p_owtkp => sbz%rv_p_owtkp 


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
      call isanrg(lmet,0,4,'bndfp:','lmet', .true.)

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

      call dfaugm(nbas,lcplxp,lso,ssite,sspec,osig,otau,oppi)

C --- Make the potential sans XC part ---
      if(llmfgw) then
        if (jobgw .eq. 1 .or. jobgw .eq. -999) then
          call info(20,1,0,' Make potential without XC part ...',0,0)
          call togpr()
          i = 1 + 10*lfrzw + 100
          call dfaugm(nbas,lcplxp,lso,ssite,sspec,osigx,otaux,oppix)
c          call defcc(ospotx, -k1*k2*k3*nsp)
          allocate(spotx(k1,k2,k3,nsp))
          spotx=0d0
c          call defcc(osmrox, -k1*k2*k3*nsp)
c          call dcopy(k1*k2*k3*nsp*2,w(osmrho),1,w(osmrox),1)
          call mkpot ( nbas , ssite , sspec , slat , sham , spot , 0 , 
     .    lcplxp , k1 , k2 , k3 , w ( osmrho ) , w ( oorhat ) , qbg , spotx 
     .    , qmom_rv , vconst , vrmt , osigx , otaux , oppix , ppnl_rv 
     .    , hab_rv , vab_rv , sab_rv , qval , qsc , gpot0_rv 
     .    , vval_rv , fes1_rv , i , vorb , nlibu , lmaxu , lldau 
     .    , 0 ) 

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
      call mkpot ( nbas , ssite , sspec , slat , sham , spot , lfrce 
     ., lcplxp , k1 , k2 , k3 , w ( osmrho ) , w ( oorhat ) , qbg , 
     .w ( osmpot ) , qmom_rv , vconst , vrmt , osig , otau , oppi 
     ., ppnl_rv , hab_rv , vab_rv , sab_rv , qval , qsc 
     ., gpot0_rv , vval_rv , fes1_rv , i , vorb , nlibu , 
     .lmaxu , lldau , 0 ) 

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
C         call isanrg(lwsig,-1,5,'BNDFP:','lwsig',.true.)
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

      call defrr(oevl, ndham*nsp*nkp)
C --- Setup for BZ integration ---
      if (plbnd .eq. 0) then

cccccccccccccccccccccc
c         print *,'nsp,nspc,nkp,nkp=',nsp,nspc,nkp,nkp
cccccccccccccccccccccc


        call subzi(lmet,ltet,lrout.gt.0,ndham,nsp,nspc,nkp,qval-qbg,
c     .    nevmx,lwtkb,eferm,oevl,owtkb,numq,lswtk,oswtk,ef0,def,mpsord,
     .  nevmx,lwtkb,eferm,owtkb,numq,lswtk,oswtk,ef0,def,mpsord,
     .  dabs(esmear-mpsord))

        if (lwtkb .eq. -1)
     .  call info(20,0,0,' Start first of two band passes ...',0,0)
        if (lwtkb .eq. 1) then
          if (ef0 .ne. eferm)
     .    call info(20,0,1,'%8pReplace ef0 with file ef=%;6d',eferm,0)
          ef0 = eferm
        endif
#if MPIK
        if (epsovl .ne. 0) then
          allocate(eomin_rv(nsp*nkp))

        endif
#endif
      else
        nkp = 0
        ldos = 0
        lwtkb = -1
        icls = 0
      endif

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
        call dfqkkl(nbas,ssite,sspec,numq,oqkkl)
        if (lekkl .eq. 1)
     .  call dfqkkl(nbas,ssite,sspec,numq,oeqkkl)
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
        nlmax = nglob('mxorb') / nglob('nkaph')
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
      ifiogw=9038
      inquire(unit=ifiogw,opened=nexist)
      if(nexist) 
     &stop 'bndfp: I assumed ifio=9038 is not used(too sloppy)'
      open(ifiogw,file='NBANDGW')
Cki  ifiese default value is not set.  set it to 0
      ifiese=0
      if(lrsig/=0) then
        ifiese=9039
        inquire(unit=ifiese,opened=nexist)
        if(nexist) 
     &  stop 'bndfp: I assumed ifio=9039 is not used(too sloppy)'
        open(ifiese,file='ESEAVR',status='old')
      endif


      if(llmfgw) then
C --- GW driver ---

c       osgw=spot%osgw
        nchan=spot%nlma

C ... Modes 0..4
        call sugw ( ssite , sspec , slat , sham , nbas , ndham , w ( 
     .  osmpot ) , vconst , lcplxp , osig , otau , oppi , ppnl_rv 
     .  , vrmt , spotx , osigx , otaux , oppix , jobgw ) 

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



C --- Start loop over k points; also, re-entry for second band pass ---
   99 continue
      ebot = 1000d0
      call surho ( nbas , ssite , sspec , lmet , ldos , lrout , lekkl 
     ., numq , oqkkl , oeqkkl , k1 , k2 , k3 , srout_zv , ndos , 
     .dos_rv , sumev , sumqv ) 

      if (lfrce .gt. 0) then
        call dpzero(frc, 3*nbas*numq)
        call dpzero ( fh_rv , 3 * nbas ) 

      endif

      if (lswtk .eq. 1) then
        efmax = 1d3
        nevmx = ndhamx
        call dpzero(w(oswtk),ndhamx*nkp)
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
          call sumlst(0,nchmx,nbas,ng,ssite,sspec,mulopt,moddos,
     .    nsites,lsites,lmxch,nchan,w,nll)
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
          call getef(ifi,0,ef0)
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
        if (procid .eq. master) then
          call suqlst ( plbopt , iopq , ndhamx , ef0 , i , w , nfbn , ifbls_iv 
     .    , nkp , qp , onesp ) 

        endif
cccccccccccccccccccccccccccccccccccccccccccccccccccc
        if(debug) print *,'vvvvv nkp222=',nkp,iopq
c        stop 'xxxxxxxxxxxxxxxxxxxxxx'
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
        call mpibc1(nkp,1,2,mlog,'bndfp','nkp')
        if (nkp .le. 0) call rx0('bndfp')
        call mpibc1(nfbn,2,2,mlog,'bndfp','nfbn')
        call mpibc1(onesp,1,2,mlog,'bndfp','onesp')

C MPIK: Setup to assemble all k-points into single list with qp table
#if MPIK
        if (nfbn(1) .gt. 0 .or. nfbn(2) .gt. 0) then
          call rx('Cannot use color weights with MPIK')
        endif

        call suqlsm(i)
C       Re-allocate qp and evl arrays
Cchp1         call defrr(oqp,  3*nkp)
         allocate(rv_p_oqp(abs(3*nkp)))
         if (3*nkp<0) rv_p_oqp(:)=0.0d0

c        call defrr(oevl, ndham*nsp*nkp)
        if (epsovl .ne. 0) then
          allocate(eomin_rv(nsp*nkp))
        endif
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
          call suqlst ( plbopt , 1 , ndhamx , ef0 , i , w , nfbn , ifbls_iv 
     .    , i2 , qp , onesp ) 

          call poppr
          if (i2 .gt. 0) then
            call pshpr(0)
            do  i1 = 1, i2
              iq = iq+1
              call suqlst ( plbopt , 1 , ndhamx , ef0 , i , w , nfbn , ifbls_iv 
     .        , i2 , qp , onesp ) 

Cchp1               call dpscop(qp,w(oqp),3,1,3*iq-2,1d0)
               call dpscop ( qp , rv_p_oqp , 3 , 1 , 3 * iq - 2 , 1d0 ) 

            enddo
            call poppr
            call suqlsm(i)
            if (i .ne. 3) goto 199
          endif
        endif
Cchp1         call mpibc1(w(oqp),3*nkp,4,.false.,'bndfp','qp')
         call mpibc1 ( rv_p_oqp , 3 * nkp , 4 , .false. , 'bndfp' , 'qp' 
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
        nlmax = nglob('nlmax')
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
      do  iq = 1, nkp
        if (plbnd .ne. 0) then
          i = nsp
          if (onesp .ne. 0 .or. nspc .eq. 2) i = 1
          call suqlst ( plbopt , 0 , ndhamx , ef0 , i , w , nfbn , ifbls_iv 
     .    , nkp , qp , onesp ) 

        else
Cchp1           call dpscop(w(oqp),qp,3,3*iq-2,1,1d0)
           call dpscop ( rv_p_oqp , qp , 3 , 3 * iq - 2 , 1 , 1d0 ) 

        endif
        qplist(:,iq)=qp
      enddo
      print *,'-------- qplist --------'
      do iq=1,nkp
        write(6,"(i5,3f8.3)")iq,qplist(:,iq)
      enddo

C --- New sigma inerpolation mode. Simplified.
      if(debug) print *,'vvv: new sigma=',newsigmasw()
      sigmamode = mod(lrsig,10) .ne. 0
c --- get index for hamiltonian for m_hamindex takao june2009
c this must be used for the case of QSGW.
      if(sigmamode.and.iter==1) then
c takao to get all qp. I think suqlst looks too complicated...
c this block is just to get list of qp. suqlst should be cleaned up in future.

        nbasp=sctrl%nbasp
        nl=sctrl%nl


Cchp1         oidxsh=sham%oindxo
         iv_p_oidxsh => sham%iv_p_oindxo 


        print *,'goto gen_hamindex nkp=',nkp
Cchp1         call gen_hamindex(sham,sbz,ssite,sspec,slat,nl,nbasp,w(oidxsh),ldim,qplist,nkp,siginit)
         call gen_hamindex ( sham , sbz , ssite , sspec , slat , nl , 
     .   nbasp , iv_p_oidxsh , ldim , qplist , nkp , siginit ) 

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
        if (ldim .ne. ndham) then !june2009 takao
          write(6,"(a,2i5)")" --- ldim(dim of lmto)   ; napwmx= ",ldim,napwmx
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

        opos=slat%opos

        call dcopy(nbas*3,w(opos),1,pos,1)
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



C --- do loop for iq and isp ---------------
      if(debug) print *,' --- goto do loop 2010 ---'
ctakao
      if(allocated(nevls)) deallocate(nevls)
      allocate(nevls(nkp,nsp))
      nevls=0
C Start k and isp loop (first loop in parallel mode)
C     h,s,evecs are dimensioned ndimh in this loop (may be q dependent)
C     evl is dimensioned evl(ndham,2)
C     w(oevl), w(owtkb), w(oswtk) are dimensioned (ndham,nsp,nkp)
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
Cchp1         call dpscop(w(oqp),qp,3,3*iq-2,1,1d0) !this should be modified as in non-MPIK mode.
         call dpscop ( rv_p_oqp , qp , 3 , 3 * iq - 2 , 1 , 1d0 ) 

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

c        call defrr(oh,-ndhamx**2*2)
c        call defrr(os,-ndhamx**2*2)

        if(debug) print *,'goto spin loop nsp=',nsp
        do  isp = 1, nsp
          if(debug)  print *,'vvvvvv isp loop isp=',isp
          if (onesp .eq. 0 .or. isp .eq. onesp) then

cccccccccccccc
ctakao do we need to do shorbz? --> furthermore, shorbz is probably wrong.---> see new algolism shortn3.
            qpo=qp
            if(.not.noshorbz()) then
              if (isp .eq. 1 .or. isp .eq. onesp) call shorbz(qp,qp,qlat,plat)
            endif


ccccccccccccccccccccccccc
c        qp=(/-1d0/3d0,1d0/3d0,1d0/3d0/)
c        qp=(/.5d0,.5d0,-.5d0/)
c        qp=(/-1d0,0d0,0d0/)
ccccccccccccccccccccccccc


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
     .        0,0,napw,dum,dum,dum,dum)
              call poppr
c              call defi(oigv2,3*napw)
              allocate(igv2x(3,napw))
              allocate(kv_iv(3*napw))

              call pshpr(iprint()-10)
              !print *,'xxx xxx 2222'
              call gvlst2 ( alat , plat , xv , 0 , 0 , 0 , pwgmin , pwgmax 
     .        , 0 , 2 , 0 , napw , kv_iv , dum , dum , igv2x ) 

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
Cchp1             oqsig=sham%oqsig
             rv_p_oqsig => sham%rv_p_oqsig 

#ifdef KINODEBUG
             write(*,'(a,i3,i3,i12)')'oqsig=',procid,iq,%loc(rv_p_oqsig)
#endif
cki            if (oqsig .eq. 0) oqsig = 1

            ispc = min(isp,nspc) !nspc=2 for non-colinier case
            !nspc=1

            if(debug) print *,'zzzzzzzzzzzzz11111111'
c ---
            if(isp==1) then
              if(allocated(hamm)) deallocate(hamm)
              if(allocated(ovlm)) deallocate(ovlm)
c      if(allocated(hamsom)) deallocate(hamsom)
              nnspc=1
              if(nspc==2) nnspc=4
              allocate(hamm(ndimh,ndimh,nnspc),ovlm(ndimh,ndimh,nnspc))
c      if(lso/=0) then
c        allocate(hamsom(ndimh,ndimh))
c      else
c        allocate(hamsom(1,1)) !dummy
c      endif
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
c      print *,' mmm: qp=',qp


            if(debug) print *,' mmm1: qp=',qp
            if(sigmamode) then
c --- Generate interpolation weight wgtq.
              allocate( wgtq(nbas,nbas,nk1,nk2,nk3) )
              pi2i = 2d0*4d0*datan(1d0)*(0d0,1d0)
              wgtq = 0d0
              !print *,' mmm2: qp=',qp
              do ib1=1,nbas
                do ib2=1,nbas
                  do ii = 1,npair(ib1,ib2)
                    rrr = matmul(plat, nlat(:,ii,ib1,ib2))
ccccccccccccccccccccccccccccccccccc
c          write(6,"(a, 3i5,2x,3f8.3)")'vvv1: qwgt=', ib1,ib2,ii
ccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccc
c          qpx=(/-0.50000D+00, -0.50000D+00, -0.50000D+00/)
ccccccccccccccccccccccccccccccccccccccccccc

c wgtq is interpolation weight for MTO.
                    do ik1=1,nk1
                      do ik2=1,nk2
                        do ik3=1,nk3
                          wgtq(ib1,ib2,ik1,ik2,ik3) = wgtq(ib1,ib2,ik1,ik2,ik3) 
     &                    + qwgt(ii,ib1,ib2)*exp( -pi2i*sum( (qp-qsmesh(:,ik1,ik2,ik3))*rrr ) )
c     &      + qwgt(ii,ib1,ib2)*exp( -pi2i*sum( (qpx-qsmesh(:,ik1,ik2,ik3))*rrr ) )
                        enddo
                      enddo
                    enddo
                  enddo
                enddo
              enddo
              !deallocate( nlat, qwgt )
cccccccccccccccccccc
c        do ib1=1,nbas
c        do ib2=1,nbas
c        do ik1=1,nk1
c        do ik2=1,nk2
c        do ik3=1,nk3
c          write(6,"(a, 2i3,3i8,2x,2f8.3)")'vvv: wgtq=', ib1,ib2,ik1,ik2,ik3, wgtq(ib1,ib2,ik1,ik2,ik3)
c       enddo
c       enddo
c        enddo
c        enddo
c        enddo
cccccccccccccccccccc

c ---  MTO part of interpolated sigm(i,j)
              allocate(sene(ndimh,ndimh))
              sene=0d0
c       open(3012,file="sfz")
              do j1=1,ldim
                do j2=1,ldim
                  ib1= ibasindex(j1)
                  ib2= ibasindex(j2)
                  do ik1=1,nk1
                    do ik2=1,nk2
                      do ik3=1,nk3
c           write(3012,"(6i4,2d13.5)")ik1,ik2,ik3,j1,j2,isp,sfz(ik1,ik2,ik3,j1,j2,isp)
c           if(j1==1.and.j2==1) then
c          if( sum(abs(qp-qsmesh(:,ik1,ik2,ik3)))<1d-8) then
c            write(6,"(a, 3f8.3,2x,3f8.3,3i6)")'vvv: qp=',qp,qsmesh(:,ik1,ik2,ik3),ik1,ik2,ik3
c          endif
c            endif
                        sene(j1,j2)= sene(j1,j2) + wgtq(ib1,ib2,ik1,ik2,ik3)*sfz(ik1,ik2,ik3,j1,j2,isp)
                      enddo
                    enddo
                  enddo
                enddo
              enddo
c       close(3012)
              deallocate(wgtq)

ccccccccccccccccccccccccccccc
              ikt=getikt(qp) ! index for qp
              if(napwk(ikt)==0) then
                print *,' napw=0 skip napw part '
                goto 1019
              endif
ccccccccccccccccccccccccccccc

C --- APW part of interpolarion sigm(i,j)
c... range check qqq is in [ikc1/nk1, (ikc1+1)/nk1],... qpara is [0,1],[0,1],[0,1]
              print *,' --- goto qparam--- '
              call qparam(qp,qsmesh,nk1,nk2,nk3,plat, ikc,qpara,igshift)
              ! eight corners of a cube -->
              !  do ix1=0,1
              !  do ix2=0,1
              !  do ix3=0,1
              !    qsmesh(1:3,ikc(1,ix1,ix2,ix3),ikc(2,ix1,ix2,ix3),ikc(3,ix1,ix2,ix3))
              !  enddo
              !  enddo
              !  enddo
c... interpolation weight at eight corners of a cube
              call interparam(qpara, wgc)
              print *,' --- end of interparam ---'
cccccccccccccccccccccccccc
c        wgc=0d0
c        wgc(0,0,0)=1d0
cccccccccccccccccccccccccc

c... obtain G index for eight corners of a cube. qp is in the cube.
              allocate( igx(ndimh,0:1,0:1,0:1) )
              ikt=getikt(qp) ! index for qp
              print *,'ikt qp=',ikt,qp
              nkk= (/nk1,nk2,nk3/)

              do ig=1,napwk(ikt) !qp+G
                do ix1=0,1
                  do ix2=0,1
                    do ix3=0,1
                      igvecc0 = igv2(:,ig,ikt) + igshift(:,ix1,ix2,ix3)
                      qmm = qsmesh(1:3,ikc(1,ix1,ix2,ix3),ikc(2,ix1,ix2,ix3),ikc(3,ix1,ix2,ix3))
                      iktm = getikt(qmm) !index for qmm
                      qqq = matmul(qp-qmm,plat)

                      igvecc = igvecc0 + qqq + 10000 + 1d-8 ! +10000 is to make igvecc+qqq positve. and then subtract -10000
                      igvecc = igvecc - 10000              !

                      igx(ig,ix1,ix2,ix3) = igv2rev(igvecc(1),igvecc(2),igvecc(3),iktm)
                      !igx at eight corners. If missing, 9999999 !see suham.F
                      ! if(ix1==0.and.ix2==0.and.ix3==0) then
                      !   write(6,"(a,3f8.3)")'=', igvecc0 + qqq- igvecc
                      !   print *,'iktm qmm=',iktm,qmm
                      !   write(6,"(a,20i8)")'ix1 ix2 ix3 igvecc=', igvecc,igshift(:,ix1,ix2,ix3),igv2(:,ig,ikt)
                      !   write(6,"(a,10i8)")'ix1 ix2 ix3 ig igx=', ix1,ix2,ix3, ig, igx(ig,ix1,ix2,ix3)
                      ! endif
                    enddo
                  enddo
                enddo
              enddo
c
              print *,'goto sss section'
              do j1=1,ndimh
                do j2=1,ndimh
                  if(j1<=ldim.and.j2<=ldim) cycle
ccccccccccccccc
c         if(j1>ldim.and.j2>ldim) cycle
ccccccccccccccc
                  sss=0d0
                  www=0d0
                  do ix1=0,1
                    do ix2=0,1
                      do ix3=0,1
                        !if(wgc(ix1,ix2,ix3)==0d0) cycle
                        if(j1<=ldim) then
                          j1x= j1
                        else
                          j1x= igx(j1-ldim,ix1,ix2,ix3)+ldim  !G1
                        endif
                        if(j2<=ldim) then
                          j2x= j2
                        else
                          j2x= igx(j2-ldim,ix1,ix2,ix3)+ldim  !G2
                        endif
                        !print *,'j1x j2x=',j1-ldim,j2-ldim,j1x,j2x,ik1,ik2,ik3
                        !if(j1x == 999999+ldim) cycle
                        !if(j2x == 999999+ldim) cycle
                        if(j1x == 999999+ldim) goto 888
                        if(j2x == 999999+ldim) goto 888
                        ik1= ikc(1,ix1,ix2,ix3)
                        ik2= ikc(2,ix1,ix2,ix3)
                        ik3= ikc(3,ix1,ix2,ix3)
c          if(wgc(ix1,ix2,ix3)/=0d0) then
c            write(6,"(a,3i5,a,3i5,d13.6)")
c     &      ' ik1 ik2 ik3=',ik1,ik2,ik3,' ix1 ix2 ix3=',ix1,ix2,ix3,wgc(ix1,ix2,ix3)
c          endif
                        !print *,' ix=',ix1,ix2,ix3
                        sss = sss+ wgc(ix1,ix2,ix3)* sfz(ik1,ik2,ik3,j1x,j2x,isp)
                        www = www+ wgc(ix1,ix2,ix3)
                        !print *,' xxx'
                      enddo                                         ! sigma at eight corners
                    enddo
                  enddo
                  !print *,'j1 j2 www=',j1,j2,www
                  sene(j1,j2)= sene(j1,j2) + sss/www

 888              continue
                enddo
              enddo
              deallocate(igx)
 1019         continue
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


C --- Set up hamiltonian.
ccccc     if (mod(mode,10) .eq. 0) then
            if(debug) print *,' mode ispc=',i,ispc
            if( newsigmasw() ) then
              print *,' --- hambl --- sigmamode=',sigmamode
              jx = 0
              if (lso == 1) jx = 10
              call hambl(jx,nbas,ssite,sspec,slat,sham,isp,qp,k1,k2,k3,
     .        w(osmpot),vconst,osig,otau,oppi,lcplxp,alfa,ndimh,napw,igv2x,
     .        hamm(1,1,ispc),ovlm(1,1,ispc),hamm(1,1,3))
              if(sigmamode) hamm(:,:,ispc) = hamm(:,:,ispc)+sene
              if(sigmamode) deallocate(sene)
ccccc     .    h(1,1,ispc),s(1,1,ispc),h(1,1,3))
c      endif
c      lrsig = mod(mode/10,10)
c      lwsig = mod(mode/1000,10)
c      ismidb = 0
c      call isanrg(lwsig,0,5,'hambls:','lwsig',T)

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
Cchp1               call hambls(i,nbas,ssite,sspec,slat,sham,isp,ispc,qp,k1,k2,k3,
Cchp1      .        w(oqsig),nqsig,w(osmpot),vconst,osig,otau,oppi,lcplxp,lso,
Cchp1      .        alfa,ndimh,napw,igv2x,hamm,ovlm,ismidb)
               call hambls ( i , nbas , ssite , sspec , slat , sham , isp , 
     .         ispc , qp , k1 , k2 , k3 , rv_p_oqsig , nqsig , w ( osmpot ) 
     .         , vconst , osig , otau , oppi , lcplxp , lso , alfa , ndimh , 
     .         napw , igv2x , hamm , ovlm , ismidb ) 

c        nsmidb = min(nsmidb,ismidb) !????
              if(sigmamode) deallocate(sene)
              print *,' end of hambls mode=',i
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

            if (mod(iq,10) .ne. 1) call pshpr(iprint()-6)
#ifndef MPIK
            call info5(30,0,0,' bndfp:  kpt %i of %i, k=%3:2,5;5d'//
     .      '%?#n#   ndimh = %i##',
     .      iq,nkp,qp,mod(pwmode/10,10),ndimh)
#endif
            nmx = min(nevmx,ndimhx)
            if (lwtkb .eq. -1) nmx = -1
C       Need all eigenvalues if 'fat bands' plotting mode
            if (nfbn(1) .gt. 0) then
              nmx = ndimhx
              efmax = 99999
            endif

ctakao
            if(napw/=0.and.ncutovl/=0) then
              if(lpdiag ==1 .or. lpdiag==2) then
                stop 'bndfp: lpdiag-PW mode with OVNCUT/=0 is not yet'
              endif
            endif
c
            if(debug) print *,'lpdiag=',lpdiag



            if(lpdiag .eq. 1) then
              call rx('parallel diag is not yet. Need to recover followings.')
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
              call phmbls(2,ndimhx,evl(1,jsp),w,w,w,w,w,hamm)!w(oh))
              if(debug) print *,'vvv lwtkb=',lwtkb
              if (lwtkb .ne. -1) then
                allocate(ww_zv(ndimhx**2*2))

                if (ldim .ne. ndham) call
     .          rx('blsig needs separate dim for wtkp,evec')
Cchp1                 call blsig ( 1 + lrsig * 10 , nbas , sham , isp , nsp , nspc 
Cchp1      .          , plat , qp , lwtkb , qval - qbg , iq , w ( owtkp ) , w ( owtkb 
Cchp1      .          ) , ndimh , t_zv , ww_zv ) 
                 call blsig ( 1 + lrsig * 10 , nbas , sham , isp , nsp , nspc 
     .           , plat , qp , lwtkb , qval - qbg , iq , rv_p_owtkp , w ( owtkb 
     .           ) , ndimh , t_zv , ww_zv ) 


                deallocate(ww_zv)

              endif
              if(debug) print *,'vvv2 lwtkb=',lwtkb
            else
              if(debug)  print *,' vvv zzzzz vvvv lpdiag else'

              if (nspc .eq. 2) then
c                 allocate(ww_rv(ndimhx**2))
                if(debug)  print *,' vvv zzzzz 111'
c                 call sopert ( 0 , ndimh , nspc , ww_rv , hamm , hamm ) 
                call sopert2 ( 0 , ndimh , nspc, hamm , hamm ) 
                if(debug)  print *,' vvv zzzzz 22'
c                 call sopert ( 0 , ndimh , nspc , ww_rv , hamm , hamm ) 
                call sopert2 ( 0 , ndimh , nspc,  ovlm , ovlm ) 
                if(debug)  print *,' vvv zzzzz 333'

c                 deallocate(ww_rv)
                if(debug)  print *,' vvv zzzzz 444'

              endif
              if(debug)  print *,' vvv zzzzz 555'

              allocate(ww_rv(11*ndimhx))
c takao
c          if (epsovl .eq. 0) then
c

              if(debug)  print *,' vvv zzzzzzzzzzzzzz bbbbbbbb'
c          if(.false.) then
              if(napw == 0) then
c            call zhev(ndimhx,w(oh),w(os),T,T,nmx,efmax,nev,w(oww),F,-1,
                !print *,'goto zhev ndimh ndimhx=',ndimh,ndimhx
                call zhev ( ndimhx , hamm , ovlm , t , t , nmx , efmax , nev 
     .          , ww_rv , f , - 1 , evl ( 1 , jsp ) , t_zv ) 

                !print *,'end of zhev'
                nevl = ndimhx
              else
                nevl = -1
                call dvset ( ww_rv , 1 , 1 , 99999d0 ) 

c            call zhevo(ndimhx,ndimhx,w(oh),w(os),T,nmx,efmax,epsovl,
c     .        nevl,nev,evl(1,jsp),w(oww),ndimhx,w(ot))
c takao version. Not in lm-7.0beta
c            call zhev_tk(ndimhx-napw,ndimhx,w(oh),w(os),nmx,epsovl,
                call zhev_tk ( ndimhx - napw , ndimhx , hamm , ovlm , nmx , epsovl 
     .          , ncutovl , nevl , nev , evl ( 1 , jsp ) , ww_rv , ndimhx 
     .          , t_zv , iq = = 1 ) 

                print *,'end of zhev_tk'
                eomin = dval ( ww_rv , 1 ) 

              endif
              deallocate(ww_rv)
            endif


C       Pad evals between ndimh and ndham with a large positive number
C       to avoid mixing up integration routines
            if (ndhamx .gt. nevl .and. nspc .eq. 2) then
              call dvset(evl,1+nevl,ndhamx,99999d0)
            elseif (ndham*nspc .gt. nevl) then
              call dvset(evl(1,jsp),1+nevl,ndham,99999d0)
            endif
ctakao
            nevls(iq,isp)=nevl

            if(debug) print *,'iq isp nevls=',iq,isp,nevls(iq,isp)

#ifndef MPIK
            if (epsovl .ne. 0) then
              call info5(30,0,0,' Overlap''s smallest eigenvalue: %;3g.  '//
     .        '%?#(n>0)#H dim reduced from %i to %i#H dim not reduced#',
     .        eomin,ndimhx-nevl,ndimhx,nevl,0)
            endif
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
C         Check for cases when nevmx is too small : i=2 => fatal error
                i = 0
                if (nevmx.ge.0 .and. lmet .ne. 0) then
                  dum = evl(max(nev,1),jsp)
C           if (ef0 .ge. dum) i = 2
                  if (.not. ltet .and. ef0+5*dabs(esmear-mpsord) .gt. dum) i=2
                  if (lmet.eq.4 .and. ef0+def+5*dabs(esmear-mpsord).gt.dum)i=2
                endif
                if (i .eq. 2) then
                  if (iprint() .gt. 0)
     .            call awrit3('%N evl(nev=%i)=%;3d but '//
     .            'ef0=%;3d ... restart with larger efmax or nevmx',
     .            ' ',80,stdo,nev,evl(max(nev,1),jsp),ef0)
                  call rx('bndfp')
                endif
              endif
            endif
#endif
c##### !MPIK

            if( debug) print *, 'tttttttttt 11111 plbnd=',plbnd
C   ... Save data for this qp to moments file
            if (plbnd .eq. 0) then
C   ... Copy eigenvalues into array containing ev for all qp
              call dpscop(evl(1,jsp),w(oevl),ndhamx,
     .        1,1+ndham*(jsp-1+isqp),1d0)
#if MPIK
              if (epsovl .ne. 0) then
                call dpscop ( eomin , eomin_rv , 1 , 1 , 1 + ( jsp - 1 + isqp 
     .          ) , 1d0 ) 

              endif
#endif

              if( debug) print *, 'tttttttttt 22222'

C        call prmx('ev',w(oevl),ndham,ndham,isqp+2)
C        call prmx('ev',w(oevl),ndham,ndham,nkp*nsp)





C   --- Orbital magnetic moment (requires lso) ---
              if (lso .ne. 0 .and. lwtkb .ne. -1) then
                if (lwtkb .eq. 0) call
     .          rx('metal weights required for orb. moment calculation')
                allocate(auso_zv(nlmax*ndhamx*3*nsp*nbas))
                auso_zv(:)=0.0d0

                call makusq ( 0 , ssite , sspec , slat , sham , nbas , nbas , 
     .          0 , nlmax , ndham , ndimh , napw , igv2x , nev , nsp , 
     .          nspc , jsp , 1 , qp , t_zv , ppnl_rv , auso_zv ) 

                call mkorbm ( ssite , sspec , jsp , nsp , nspc , nlmax , ndham 
     .          , nev , w ( owtkb ) , iq , nbas , ppnl_rv , auso_zv , nl 
     .          , nkp , orbtm_rv ) 

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


Cchp1                     oiprmb=sham%oindxo
                     iv_p_oiprmb => sham%iv_p_oindxo 


Cchp1                     call mullmf ( nbas , ssite , sspec , w ( oiprmb ) , t_zv , 
Cchp1      .              ndimh , nspc , iq , jsp , moddos , nsites , lsites , lmxch , 
Cchp1      .              nchan , chan_iv , lmdim , ndham , doswt_rv ) 
                     call mullmf ( nbas , ssite , sspec , iv_p_oiprmb , t_zv , ndimh 
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
                    call rxx(nspc.ne.1,'pdos not implemented in noncoll case')
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
                  ausu_rv(:)=0.0d0

                  call makusq ( 0 , ssite , sspec , slat , sham , nbas , nbas , 
     .            0 , nlmax , ndham , ndimh , napw , igv2x , nev , nsp , 
     .            nspc , jsp , 1 , qp , t_zv , ppnl_rv , ausu_rv ) 

                  call mkdmtu ( ssite , sspec , w ( owtkb ) , jsp , iq , nsp , 
     .            nspc , ndham , nlmax , nbas , nev , ppnl_rv , ausu_rv , 
     .            dmatu , nlibu , lmaxu , lldau ) 

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
     .          nsp , nspc , jsp , iq , qp , t_zv , ppnl_rv , ausc_zv 
     .          ) 


C       Accumulate output density and sampling DOS
C       Even if no output rho, still call addrbl to make DOS when lmet=4
              elseif (lwtkb.ne.-1 .and. (lrout.ne.0 .or. lmet.eq.4)) then
C         dum=frc(3,1)
Cchp1                 call addrbl ( ssite , sspec , slat , sham , jsp , nsp , nspc 
Cchp1      .          , qp , w ( owtkp ) , ndham , ndimh , napw , igv2x , lmet 
Cchp1      .          , lrout , lwtkb , w ( owtkb ) , lswtk , w ( oswtk ) , iq , lfrce 
Cchp1      .          , ldos , lekkl , k1 , k2 , k3 , w ( osmpot ) , vconst , osig 
Cchp1      .          , otau , oppi , lcplxp , numq , qval - qbg , t_zv , evl , 
Cchp1      .          nevl , ef0 , def , esmear , emin , emax , ndos , dos_rv , 
Cchp1      .          srout_zv , sumqv , sumev , oqkkl , oeqkkl , frc ) 
                 call addrbl ( ssite , sspec , slat , sham , jsp , nsp , nspc 
     .           , qp , rv_p_owtkp , ndham , ndimh , napw , igv2x , lmet , lrout 
     .           , lwtkb , w ( owtkb ) , lswtk , w ( oswtk ) , iq , lfrce , ldos 
     .           , lekkl , k1 , k2 , k3 , w ( osmpot ) , vconst , osig , otau 
     .           , oppi , lcplxp , numq , qval - qbg , t_zv , evl , nevl , ef0 
     .           , def , esmear , emin , emax , ndos , dos_rv , srout_zv , sumqv 
     .           , sumev , oqkkl , oeqkkl , frc ) 


C          print 399,iq,frc(3,1),frc(3,1)-dum
C  399     format(' after addrbl: frc(3,1)=',i4,2f12.6)
              endif
#ifndef MPIK
            elseif (plbnd .ne. 0) then
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
              call dpscop(evl(1,jsp),w(oevl),ndhamx,
     .        1,1+ndham*(jsp-1+isqp),1d0)
              if (epsovl .ne. 0) then
                call dpscop ( eomin , eomin_rv , 1 , 1 , 1 + ( jsp - 1 + isqp 
     .          ) , 1d0 ) 

              endif
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
      call xmpbnd(kpproc,ndham,nkp,nsp,w(oevl))
      if (epsovl .ne. 0) then
        call xmpbnd ( kpproc , 1 , nkp , nsp , eomin_rv ) 

      endif
      if (icls .ne. 0 .and. lwtkb .ne. -1) then
        call rxx(nspc.ne.1,'CLS not implemented in noncoll case')
        call xmpbnd ( kpproc , nlmax * ndham * 3 * nsites , nkp , nsp 
     .  , ausc_zv ) 

      endif

ctakao
      call mpibc2(nevls,nkp*nsp,2,mlog,'bndfp','nevls')

ckino 05Feb09, use kpproc later, move it.
c      deallocate(kpproc, stat=ierr)

C     Allreduce density-related quantities
      if (lrout .ne. 0) then
        call mpibc2(sumqv,6,4,mlog,'bndfp','sumqv')
        call mpibc2(sumev,6,4,mlog,'bndfp','sumev')
        call mpibc2 ( srout_zv , k1 * k2 * k3 * nsp * numq , 6 , mlog 
     .  , 'bndfp' , 'smrho' ) 

        if (lswtk .eq. 1) then
          call mpibc2(w(oswtk),ndhamx*nkp,4,mlog,'bndfp','swtk')
        endif
C     Allreduce qkkl
        nkaph = nglob('nkaph')
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
              call mpibc2(w(oqkkl(i,ib)),nelt(i)*numq*nsp*nspc,4,
     .        mlog,'bndfp','qkkl')
ctakao March1-2009
              if (lekkl .eq. 1) 
     .        call mpibc2(w(oeqkkl(i,ib)),nelt(i)*numq*nsp*nspc,4,
     .        mlog,'bndfp','eqkkl')
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
                call dpscop(w(oevl),evl(1,jsp),ndhamx,
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
            call dpscop(w(oevl),evl(1,jsp),ndhamx,
     .      1+ndham*(jsp-1+isqp),1,1d0)
            if (epsovl .ne. 0) then
              call dpscop ( eomin_rv , eomin , 1 , 1 + ( jsp - 1 + isqp 
     .        ) , 1 , 1d0 ) 

            endif
            if (mod(iq,10) .ne. 1) call pshpr(iprint()-6)
            call info5(30,0,0,' bndfp:  kpt %i of %i, k=%3:2,5;5d',
     .      iq,nkp,qp,0,0)
            if (epsovl .ne. 0) then
              call info5(30,0,0,' Overlap''s smallest eigenvalue: %;3g',
     .        eomin,0,0,0,0)
            endif
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

            if (procid.eq.master) then
              call prtev(w(1),ndimhx_iq,evl(1,jsp),nevmx,efmax,nev_iq)
            endif
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
     .  - qbg , ef0 , def , sumqv , sumev , n1 , n2 , n3 , k1 , k2 , 
     .  k3 , srout_zv , oqkkl , frc , lrep ) 

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

C --- BZ integration for fermi level, band sum and qp weights ---
      if (lmet .ge. 0 .and. (lmet .ne. 4 .or. ltet)) then
        if (lwndow) then
          allocate(ww_rv(ndham*nsp*nkp))

          eferm = min(dosw(1),dosw(2))
          call bzints ( nkabc ( 1 ) , nkabc ( 2 ) , nkabc ( 3 ) , w ( oevl 
     .    ) , ww_rv , nkp , ndham , ndham , nsp , xv , xv , xv , 1 , 
     .    eferm , 2 , ntet , w ( oidtet ) , sev , dum ) 

C         call prmx('w(min)',w(oww),ndham,ndham,nkp*nsp)
          eferm = max(dosw(1),dosw(2))
          call bzints(nkabc(1),nkabc(2),nkabc(3),w(oevl),w(owtkb),nkp,
     .    ndham,ndham,nsp,xv,xv,xv,1,eferm,2,ntet,w(oidtet),sev,dum)
C         call prmx('w(max)',w(owtkb),ndham,ndham,nkp*nsp)
          call daxpy ( ndham * nsp * nkp , - 1d0 , ww_rv , 1 , w ( owtkb 
     .    ) , 1 ) 

C         call prmx('w',w(owtkb),ndham,ndham,nkp*nsp)
          if (allocated(ww_rv)) deallocate(ww_rv)

        else
          dosrng = 8
          if (mpsord .lt. 0) dosrng = 16
Cchp1           call bzwtsf(ndham,ndham,nsp,nspc,nkabc(1),nkabc(2),nkabc(3),
Cchp1      .    nkp,ntet,w(oidtet),qval-qbg,fsmom,lmet.ne.0,ltet,mpsord,
Cchp1      .    ndos,dabs(esmear-mpsord),dosrng,w(owtkp),w(oevl),efmax,
Cchp1      .    lswtk,w(oswtk),eferm,sev,w(owtkb),sumqv(1,2),lwtkb)
           call bzwtsf ( ndham , ndham , nsp , nspc , nkabc ( 1 ) , nkabc 
     .     ( 2 ) , nkabc ( 3 ) , nkp , ntet , w ( oidtet ) , qval - qbg 
     .     , fsmom , lmet.ne.0 , ltet , mpsord , ndos , dabs ( esmear - 
     .     mpsord ) , dosrng , rv_p_owtkp , w ( oevl ) , efmax , lswtk , 
     .     w ( oswtk ) , eferm , sev , w ( owtkb ) , sumqv ( 1 , 2 ) , lwtkb 
     .     ) 


C         Store val charge & magnetic moment in sumqv(1..2)
          if (lmet .ne. 4) then
            sumqv(1,1) = sumqv(1,2)
            sumqv(2,1) = sumqv(2,2)
          endif
        endif
        if (lmet .ne. 4) ef0 = eferm
        if ( lmet .ne. 4 ) then

          sbz%ef=ef0
        endif

        if (lmet .gt. 0) then
          if (procid .eq. master) then
            ifi = fopna('wkp',-1,4)
            i = iobzwt(0,ndhamx,nkp,nspx,eferm,w(owtkb),-ifi)
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
        deallocate(sfz,qsmesh,npair,pos)
      endif

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
          i = iobzwt(1,ndham,nkp,nsp,ef0,w(owtkb),-ifi)
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
            call bzints ( nkabc ( 1 ) , nkabc ( 2 ) , nkabc ( 3 ) , w ( oevl 
     .      ) , dum , nkp , ndhamx , ndhamx , nspx , dosw ( 1 ) , dosw ( 
     .      2 ) , dos_rv , ndos , eferm , 1 , ntet , w ( oidtet ) , dum 
     .      , dum ) 

            if ( ndos0 .gt. 0 ) call xxxdif ( dosw ( 1 ) , dosw ( 2 ) , ndos 
     .      , nspx , 0 , dos_rv ) 

            del = 0d0
          else
            if (mpsord .ge. 100) mpsord = mod(mpsord,100)
Cchp1             if ( ndos0 .gt. 0 ) call makdos ( nkp , ndhamx , ndhamx , nspx 
Cchp1      .      , w ( owtkp ) , w ( oevl ) , mpsord , dabs ( esmear - mpsord 
Cchp1      .      ) , - 6d0 , dosw ( 1 ) , dosw ( 2 ) , ndos , dos_rv ) 
             if ( ndos0 .gt. 0 ) call makdos ( nkp , ndhamx , ndhamx , nspx 
     .       , rv_p_owtkp , w ( oevl ) , mpsord , dabs ( esmear - mpsord ) 
     .       , - 6d0 , dosw ( 1 ) , dosw ( 2 ) , ndos , dos_rv ) 


Cchp1             if ( ndos0 .lt. 0 ) call maknos ( nkp , ndhamx , ndhamx , nspx 
Cchp1      .      , w ( owtkp ) , w ( oevl ) , mpsord , dabs ( esmear - mpsord 
Cchp1      .      ) , - 6d0 , dosw ( 1 ) , dosw ( 2 ) , ndos , dos_rv ) 
             if ( ndos0 .lt. 0 ) call maknos ( nkp , ndhamx , ndhamx , nspx 
     .       , rv_p_owtkp , w ( oevl ) , mpsord , dabs ( esmear - mpsord ) 
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



C --- Core-level spectroscopy ---
      if (icls .ne. 0) then
        if (procid .eq. master) then

          eferm=sbz%ef

          call vcdmel ( sctrl , ssite , sspec , slat , nlmax , ndham , 
     .    ndimh , nkp , nsp , nspc , eferm , w ( oevl ) , ausc_zv , 
     .    nsites , isite , iclsl , iclsn ) 

          call rx0('done generating core level spectra')
        else
          call fexit(0,0,' ',0)
        endif
      endif

C --- Assemble output density ---
      call dfratm(nbas,ssite,sspec,orhat1)
      allocate(qbyl_rv(n0*nbas*nsp))
      allocate(hbyl_rv(n0*nbas*nsp))
C     --window: Put output density into orhat and smrho, and exit
      if (lwndow) then
        call mkrout ( ssite , sspec , slat , sham , nbas , nsp , ldim 
     .  , lekkl , oqkkl , oeqkkl , w ( oorhat ) , hab_rv , sab_rv 
     .  , qbyl_rv , hbyl_rv , lrout ) 

        call zcopy ( k1 * k2 * k3 * nsp , srout_zv , 1 , w ( osmrho 
     .  ) , 1 ) 

        call symrho ( ssite , sspec , slat , lfrce , w ( osmrho ) , w 
     .  ( oorhat ) , qbyl_rv , hbyl_rv , frc ) 

        goto 999
      endif
      call mkrout ( ssite , sspec , slat , sham , nbas , nsp , ldim 
     ., lekkl , oqkkl , oeqkkl , orhat1 , hab_rv , sab_rv , qbyl_rv 
     ., hbyl_rv , lrout ) 

C
      if (lrout .ne. 0) then
C   ... Symmetrize output density and forces
        call symrho ( ssite , sspec , slat , lfrce , srout_zv , orhat1 
     .  , qbyl_rv , hbyl_rv , frc ) 


C   ... New boundary conditions pnu for phi and phidot
C       call pshpr(iprint()-12)
        if (lpnu .gt. 0) then

          i_copy_size=size(sham%pmin)
          call dcopy(i_copy_size,sham%pmin,1,xv,1)


          i_copy_size=size(sham%pmax)
          call dcopy(i_copy_size,sham%pmax,1,xv(11),1)

          call pnunew ( nbas , nsp , ssite , sspec , xv , xv ( 11 ) , lfrzw 
     .    , hab_rv , sab_rv , qbyl_rv , hbyl_rv ) 

        endif
C       call poppr
      endif


C --- Evaluate Harris energy ---
c      call pshpr(ipr-10) !this was needed to go through test.fp (just for check write.) Why? --> but no problem... sep2009
c      print *,'vvv takao goto mkehfk mode1'
      call mkehkf(1,sham,sev,sumqv(2,1),sumtv,ehar)
c      call poppr()

C --- Evaluate KS total energy, correction to Harris force ---
      if (lrout .ne. 0) then

C   ... Correction to harris force
        call dfrce ( ssite , sspec , slat , sctrl , k1 , k2 , k3 , nvl 
     .  , w ( oorhat ) , orhat1 , elind , qmom_rv , w ( osmrho ) , 
     .  srout_zv , fh_rv ) 


C   ... Evaluate KS total energy and output magnetic moment
        eks = 0d0
        if (leks .ge. 1) then
          call togpr()
c          print *,'vvv takao mkekin 000 --- oldpot x newrho'
          call mkekin ( nbas , ldim , ssite , sspec , slat , sham , osig 
     .    , otau , oppi , lcplxp , oqkkl , k1 , k2 , k3 , vconst , w ( 
     .    osmpot ) , srout_zv , sev , sumtv ) 


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
          call mkpot ( nbas , ssite , sspec , slat , sham , spot , lfrce 
     .    , lcplxp , k1 , k2 , k3 , srout_zv , orhat1 , qbg , w ( osmpot 
     .    ) , qmom_rv , vconst , vrmt , osig , otau , oppi , ppnl_rv 
     .    , hab_rv , vab_rv , sab_rv , qval , qsc , gpot0_rv 
     .    , vval_rv , fes2_rv , 0 , vorb , nlibu , lmaxu , lldau 
     .    , 0 ) 

ccccccccccccccccccccccccccccccccccc
ctakao test mkekin repeat
c          print *,'vvv takao mkekin repeat 001x ---newpot x newrho'
c          call mkekin(nbas,ldim,ssite,sspec,slat,sham,osig,otau,oppi,
c     .      lcplxp,oqkkl,k1,k2,k3,vconst,w(osmpot),w(osrout),sev,sumtvxxx)
c          print *,'vvv takao mkekin repeat 001 --- newpot x oldrho'
c          call mkekin(nbas,ldim,ssite,sspec,slat,sham,osig,otau,oppi,
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

C   --- Add together force terms ---
        if ( lfrce .gt. 0 ) call totfrc ( nbas , ssite , slat , leks 
     .  , fes1_rv , fes2_rv , fh_rv , frc ) 


C   --- Mix input and output densities ---
        allocate(qm1_rv(nvl))

        allocate(qm2_rv(nvl))

        call upacks('strn mix',i1,i2)
#if MPE
        ierr = MPE_LOG_EVENT(EVENT_START_MIXRHO,procid,"mixrho")
#endif
        call mixrho ( ssite , sspec , slat , nsp , iter , sstrn ( i1:i2 
     .  ) , qval - qbg , elind , orhat1 , w ( oorhat ) , k1 , k2 , k3 
     .  , dmxp , srout_zv , w ( osmrho ) ) 


ccccccccccccccccccccccccc
ctakao test mkekin repeat
c          print *,'vvv takao mkekin repeat 2 ---'
c          call dcopy(k1*k2*k3*2,w(osmpotbk),1,w(osmpot),1)
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
      if (allocated(eomin_rv)) deallocate(eomin_rv)
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
      close(ifiogw)
Cki
      if ( ifiese>0 ) close(ifiese)
c      call xxxbfp ! Patch to avoid PGI compiler bug on AMD processor
      if (allocated(ifbls_iv)) deallocate(ifbls_iv)
#if MPIK
Cchp1       call rlse(oqp)
cki       if (associated(rv_p_oqp)) deallocate(rv_p_oqp)

#endif
      call rlse(oevl)
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
      integer w(1)
      common /w/ w

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


