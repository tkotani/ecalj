      subroutine bndfp(nbas,nsp,nlibu,lmaxu,lldau,
     .  ssite,sspec,slat,sctrl,sham,spot,sbz,sstrn,
     .  ndham,leks,lrout,lfrce,lpnu,dmxp,iter,maxit,evl,frc,
     .  dmatu,vorb)
      use m_rdctrl, only: ncutovl
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
C     implicit none
C ... Passed parameters
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
      integer oeomin
ckino 05Feb09
      integer:: pidorigin, ndimhx_iq, nev_iq
      integer:: work_mpi(10)  ! use only 1 and 2 now
#endif
#if MPE
C Event numbers:
      include "events.ins"
#endif
#endif
      logical mlog
      integer procid,master
      integer nbx,nbas,n0,nab,nppn,iter,maxit,ndham
      character*(*) sstrn
      parameter ( nbx=512, n0=10, nppn=12, nab=9)
      integer k1,k2,k3,leks,lrout,lfrce,lpnu
      integer nlibu,lmaxu,lldau(nbas)
      double precision ssite(1),sspec(1),slat(1),sctrl(1),sham(1),
     .  spot(1),sbz(1),dmxp(33),evl(ndham,2),frc(3,nbas)
C ... Local parameters
      logical lgors,ltet,cmdopt,lwndow
      character strn*120,plbopt*120,mulopt*120,clsopt*120,strn2*120,dc*1
      integer ipr,iprint,ipl,isp,jsp,ispc,isum,nglob,i,iopq,iq,isqp,
     .  ismidb,nsmidb,ldos,lmet,nspc,lrep,lcplxp,lso,nkabc(3),ndos,
     .  ndos0,nev,nevl,nkp,nsp,numq,n1,n2,n3,ntet,lwtkb,lswtk,nl,lgunit,
     .  igets,mpsord,nevmx,nvl,odos,ngabc(3),isw,stdo,stdl,ifi,fxst,
     .  fopn,fopna,lfrzw,i1,i2,plbnd,nfilem,iobzwt,lnoxc,lrsig,lwsig,
     .  jobgw,nll,lpdiag,iv(10),parg,ndhamx,nspx,nk1,nk2,nk3,lshft(3),
     .  onesp,nfbn(2),mpipid,ldham(8,2),ldim,ndimh,ndimhx,lekkl
      integer oevl,ofh,ofes1,ofes2,ogpot0,oh,ohab,ohbyl,oppnl,
     .  oqbyl,oqmom,os,osab,osrout,ot,ovab,oidtet,owtkb,oswtk,owtkp,
     .  oqm1,oqm2,oorhat,osoptc,osmpot,osmrho,ovval,oww,oqp,oips,oorbtm,
     .  oifbls
      integer osig(3,nbx),otau(3,nbx),oppi(3,nbx),oqkkl(3,nbx),
     .  oeqkkl(3,nbx),orhat1(3,nbx)
      equivalence (ldim,ldham(1,1))
#if LMFGW
      logical a2bin
      integer osigx(3,nbx),otaux(3,nbx),oppix(3,nbx),ospotx,osmrox,osgw
#else
      double precision sigp(10)
#endif
      equivalence (n1,ngabc(1)),(n2,ngabc(2)),(n3,ngabc(3))
      equivalence (nk1,nkabc(1)), (nk2,nkabc(2)), (nk3,nkabc(3))
      double precision ef00,eferm,emax,emin,esmear,qval,qsc,sev,
     .  sev00,sev1,sumtv,alfa,vconst,qbg,ebot,pi,def,ef0,ehar,eks,
     .  dosw(2),efmax,dgets,alat,plat(3,3),qlat(3,3),sumev(2,3),
     .  sumqv(3,2),qp(3),shftqp(3),eterms(20),elind,vol,dum,del,evtop,
     .  ecbot,vrmt(nbx),fsmom,rsrnge,xv(20),dosrng,epsovl,eomin,dval
      equivalence (emin,dosw(1)),(emax,dosw(2))

C     For now: PW basis
      integer oigv2,okv
      integer pwmode,napw
      double precision pwemin,pwemax,pwgmin,pwgmax

C     real(8),allocatable :: qpe(:,:)

      logical T, F
C     for spin-orbit
      integer oauso
C     for partial dos
      integer oausp
C     for CLS
      integer nsitmx
      parameter (nsitmx = 256)
      integer nlmax,oausc,i1mach,icls,isite(nsitmx),iclsl(nsitmx),
     .        iclsn(nsitmx)
C     for Mulliken
      integer oiprmb,odoswt,moddos,nsites,lsites(nbx),nchan,ng,iomoms,
     .        nchmx,lmdim,ochan,lmxch
C     for pzhev
      integer nblk,nprow,npcol,nmx
C     for LDA+U
      integer oausu
      double complex vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu),
     .              dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
C ... Heap
      integer w(1)
      common /w/ w

#if MPIK
      integer obuf,oqptr
      integer ib,is,igetss,lmxa,lmxh,nlma,nlmh,nelt(3),nkaph,kmax
      integer, dimension(:),allocatable :: kpproc
      double precision sttime,entime
#endif
C      for debugging and testing
C      integer ocg,ojcg,oidxcg,ocy
C      call upack('lat ocg ojcg oidxcg ocy',slat,ocg,ojcg,oidxcg,ocy,0)
C      call thxpbl(slat,w(ocg),w(oidxcg),w(ojcg),w(ocy))
#if TESTRWF
      integer ooffH,oidxsh,oistab,oag,og
      integer nk1,nk2,nk3,ogstar,is(3),ifac(3),nsgrp
      double precision rb(9),qb(9)
      logical llshft(3)
      integer oipq
#endif
C For self-energy
      integer nqsig
      integer oqsig

C     Data statements
      data T /.true./ F /.false./

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
      napw = 0
      ipr  = iprint()
      ipl  = ipr
      nsp  = nglob('nsp')
      nspc = nglob('nspc')
C     npqn = nglob('npqn')
      lso =   isw(lgors('ctrl lncol,4',sctrl))
     .    + 2*isw(lgors('ctrl lncol,32',sctrl))
     .    + 3*isw(lgors('ctrl lncol,64',sctrl))
      lcplxp = 0
      lekkl = igets('ctrl pfloat',sctrl)
      if (lso .ne. 0) lcplxp = 1
      if (isum(nbas,lldau,1) .ne. 0) lcplxp = 1
      ndhamx = ndham*nspc
      nspx = nsp / nspc
      call upack4('ham ldham pwmode pwemin pwemax',sham,ldham,pwmode,
     .  pwemin,pwemax)
      onesp = 0
      lwsig = 0
      nfbn(1) = 0
      nfbn(2) = 0
      stdo = lgunit(1)
      stdl = lgunit(2)
      ldos = igets('ctrl ldos,1',sctrl)
      lfrzw= isw(lgors('ctrl lbas,16',sctrl))
      lrsig= igets('ham lsig',sham)
      if (procid .eq. master) then
      if (lrsig .ne. 0 .and. fxst('sigm') .ne. 1) then
        call info0(1,1,0,' bndfp (warning): '//
     .    'no sigm file found ... LDA calculation only')
        call pack1('ham lsig',sham,0)
        lrsig = 0
      endif
      endif
      call mpibc1(lrsig,1,2,mlog,'rdsigm','lrsig')
      jobgw= -999
      epsovl = dgets('ham oveps',sham)
      lpdiag = isw(cmdopt('--pdiag',7,0,strn))
#if MPIK
      lpdiag = 0
#endif
#if LMFGW
      i = 7
      if (cmdopt('-jobgw=',i,0,strn)) then
        if (.not. a2bin(strn,jobgw,2,0,' ',i,-1)) call
     .    rxs2('BNDFP: failed to parse "',strn(1:30),'%a"')
      endif
#endif
      call upack('pot oorhat osmrho osmpot osoptc',spot,oorhat,osmrho,
     .  osmpot,osoptc,0)
C     icls = igets('optic cls',w(osoptc))
      pi = 4d0*datan(1d0)
      if (iprint() .ge. 20) call awrit2('%N --- BNDFP:  '//
     .  'begin iteration %i of %i ---',' ',80,stdo,iter,maxit)
      if (nbas.gt.nbx) call rxi('bndfp: nbx exceeded, need',nbas)
      if (ndham.le.0)
     .  call rx('bndfp: hamiltonian matrix has zero dimension')

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
     .     ' MPI creating process configuration .. nprow=%i npcol=%i,'//
     .     ' blocking factor %i',
     .     ' ',256,lgunit(1),nprow,npcol,nblk)
        endif
      endif
#endif

      call dvset(eterms,1,20,-99d0)
      eterms(19) = 0d0
      call pack1('ham eterms',sham,eterms)
      call upack2('pot nlml nlma',spot,nvl,nchan)
      eks = 0
      call upack('lat alat plat qlat nabc',slat,alat,plat,qlat,ngabc,0)
      vol = dgets('lat vol',slat)
      call fftz30(n1,n2,n3,k1,k2,k3)

C ... for BZ integration
      call upack('bz nkabc nkp ntet oidtet lmet',sbz,nkabc,nkp,ntet,
     .  oidtet,lmet)
      call upack('bz n w nevmx efmax fsmom',sbz,mpsord,esmear,nevmx,
     .  efmax,fsmom)
      esmear = esmear + iabs(mpsord)
      if (mpsord .lt. 0) esmear = -esmear
      mpsord = int(dabs(esmear))
      if (esmear .lt. 0) mpsord = -mpsord

      call upack('bz ndos dosw ef def',sbz,ndos,dosw,ef0,def,0)
      call upack('bz oqp owtkp',sbz,oqp,owtkp,0,0,0)
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
     .    ' BNDFP: generating density in energy window %2:1d',dosw,0)
        lfrce = 0
        lpnu = 0
        efmax = 1d2
        nevmx = ndham
        if (lrout .eq. 0)
     .    call rx('--window incompatible with no output density')
        call info0(20,0,0,' Delete band weights file ...')
        ifi = fopna('wkp',-1,4)
        call dfclos(ifi)
      endif
C     Always need DOS if just sampling
      if (lmet .ne. 0 .and. .not. ltet) ldos=1
      qbg = dgets('ctrl zbak',sctrl)
      alfa = 0

C     Switch to plot bands at specified qp
      if (cmdopt('--band',6,0,strn)) then
        plbnd = 1
        plbopt = strn(7:)
        lrout = 0
        lfrce = 0
        nkp = 0
        numq = 1
        call defi(oifbls, -ndham*nspc*2)
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
      call defrr (oqmom,  nvl)
      call defrr (ogpot0, nvl)
      call defrr (ovval,  nchan)
      call defrr (ohab,   nab*n0*nsp*nbas)
      call defrr (ovab,   nab*n0*nsp*nbas)
      call defrr (osab,   nab*n0*nsp*nbas)
      call defrr (oppnl,  nppn*n0*nsp*nbas)
      call defrr (ofes1,  3*nbas)
      call dfaugm(nbas,lcplxp,lso,ssite,sspec,osig,otau,oppi)

C --- Make the potential sans XC part ---
#if LMFGW
      if (jobgw .eq. 1 .or. jobgw .eq. -999) then
        call info(20,1,0,' Make potential without XC part ...',0,0)
        call togpr()
        i = 1 + 10*lfrzw + 100
        call dfaugm(nbas,lcplxp,lso,ssite,sspec,osigx,otaux,oppix)
        call defcc(ospotx, -k1*k2*k3*nsp)
        call defcc(osmrox, -k1*k2*k3*nsp)
        call dcopy(k1*k2*k3*nsp*2,w(osmrho),1,w(osmrox),1)
        call mkpot(nbas,ssite,sspec,slat,sham,spot,0,lcplxp,k1,k2,k3,
     .    w(osmrho),w(oorhat),qbg,w(ospotx),w(oqmom),vconst,vrmt,osigx,
     .    otaux,oppix,w(oppnl),w(ohab),w(ovab),w(osab),qval,qsc,
     .    w(ogpot0),w(ovval),w(ofes1),i,vorb,nlibu,lmaxu,lldau)
        call togpr()
        call rlse(osmrox)
      endif
#endif

C --- Make the potential and total energy terms ---
      i = 1 + 10*lfrzw + 100*lnoxc
      if (cmdopt('--wrhomt',8,0,strn)) then
        i = i + 10000
      else if (cmdopt('--wpotmt',8,0,strn)) then
        i = i + 20000
      endif
#if LMFGW
      i = i + 10000
#endif
      call togpr()
      call mkpot(nbas,ssite,sspec,slat,sham,spot,lfrce,lcplxp,k1,k2,k3,
     .  w(osmrho),w(oorhat),qbg,w(osmpot),w(oqmom),vconst,vrmt,osig,
     .  otau,oppi,w(oppnl),w(ohab),w(ovab),w(osab),qval,qsc,
     .  w(ogpot0),w(ovval),w(ofes1),i,vorb,nlibu,lmaxu,lldau)
      call togpr()

C ... Read and store self-energy addition to LDA potential
C     Read lwsig: special modes that write sig, Z to disk.
      lwsig = 0
      if (mod(lrsig,10) .ne. 0) then
C       Real-space range
        rsrnge = dgets('ham rsrnge',sham)
        ifi = 1
        if (procid .eq. master) ifi = fopna('sigm',-1,4)
        if (ldim .ne. ndham) call rx('rdsigm not ready for PW basis')
        call rdsigm(lrsig,nbas,nsp,ndham,slat,sham,sbz,ifi,rsrnge,lwsig)
        if (lwsig .ne. 0) then
          lrout = 0
          lfrce = 0
          lwtkb = 0
          call rxx(plbnd.ne.0,'incompatible options, lwsig and plbnd')
        endif
        call fclose(ifi)
        call phmbl3(1,0,0,0,0,0,0,0)
      endif
C     lwsig0 = isign(1,lwsig) * mod(iabs(lwsig),10)
C     lwsig1 = iabs(lwsig/10)
      if (cmdopt('--evec',6,0,strn)) then
        if (lwsig .ne. 0) call rxi('bandfp: --evec not allowed with '
     .      //'--wsig:trans=',i)
        strn2 = strn(7:)
        dc = strn2(1:1)
        call evcflg(dc,strn2,i,shftqp)
        if (i .gt. 10) lwsig = 2 + mod(i,10)
      endif

C     Transformation modes: read qp and jump to start
      if (lwsig .ne. 0) then
        call isanrg(lwsig,-1,5,'BNDFP:','lwsig',.true.)
        onesp = 0
        lwtkb = 0

        ifi = fopna('qpts',-1,0)
        call getqp(0,ifi,nkp,nkabc,lshft,i,w,w,w)
        call defrr(oqp,3*nkp)
        call getqp(1,ifi,nkp,nkabc,lshft,i,w(oqp),w,w)

        goto 50
      endif

      elind = dgets('ham elind',sham)
      if (elind .lt. 0d0)
     .  elind=-(3*pi**2*(qval-qsc-qbg)/vol)**.66666d0*elind
      dmxp(33) = elind

C --- Setup for BZ integration ---
      if (plbnd .eq. 0) then
        call subzi(lmet,ltet,lrout.gt.0,ndham,nsp,nspc,nkp,qval-qbg,
     .    nevmx,lwtkb,eferm,oevl,owtkb,numq,lswtk,oswtk,ef0,def,mpsord,
     .    dabs(esmear-mpsord))
        if (lwtkb .eq. -1)
     .    call info(20,0,0,' Start first of two band passes ...',0,0)
        if (lwtkb .eq. 1) then
          if (ef0 .ne. eferm)
     .    call info(20,0,1,'%8pReplace ef0 with file ef=%;6d',eferm,0)
          ef0 = eferm
        endif
#if MPIK
        if (epsovl .ne. 0) then
          call defrr(oeomin,nsp*nkp)
        endif
#endif 
      else
        nkp = 0
        ldos = 0
        lwtkb = -1
        icls = 0
      endif

C ... Set switch to write file sigii
#ifndef LMFGW
      if (lrsig .ne. 0 .and. plbnd .eq. 0 .and. procid .eq. master) then
        i = fopna('sigii',-1,0)
        rewind i
        call pack1('ham sigp:9.9',sham,1d0)
      endif
#endif

C --- More k-point independent local arrays ---
      if (ldos .eq. 0 .and. lmet .eq. 0) ndos = 1
      if (ndos .lt. 0) ndos = -ndos
      call defrr (odos, ndos*2*nsp)
      if (lrout .ne. 0) then
        call dfqkkl(nbas,ssite,sspec,numq,oqkkl)
        if (lekkl .eq. 1)
     .  call dfqkkl(nbas,ssite,sspec,numq,oeqkkl)
        call defcc(osrout,  k1*k2*k3*numq*nsp)
        call defrr (ofh,    3*nbas)
        call defrr (ofes2,  3*nbas)
      endif
      call suham2(sctrl,slat,sspec,ssite,sham,spot,sstrn)

C --- Options for core level specta (CLS) ---
C     Allocate aus for all qp
      if (cmdopt('--cls',5,0,strn)) then
        if (lmet .gt. 0 .and. (lmet .ne. 2 .and. lmet .ne. 3))
     .    call rx('For CLS restart with METAL=2 or 3')
        icls = 1
        clsopt = strn(6:)
        efmax = 1d3
        nevmx = ndhamx
        if (lrout .eq. 0) call rx('bndfp: need output density for cls')
        call suclst(nsitmx,nbas,nsp,ssite,sspec,clsopt,
     .              isite,iclsl,iclsn,nsites)
        efmax = 1d3
        nevmx = ndhamx
        nlmax = nglob('mxorb') / nglob('nkaph')
        if (16*3*nlmax*ndham*nbas*nsp*nkp/1000000 .gt. 24 .and. procid
     .                                                     .eq. master)
     .  call awrit6(' CLS: %iMb memory for aus: nlmax=%i ndham=%i '//
     .              'nsistes=%i nsp=%i nkp=%i',' ',120,
     .             i1mach(2),16*3*nlmax*ndham*nsites*nsp*nkp/1000000,
     .             nlmax,ndham,nsites,nsp,nkp)
        call defcc(oausc,-3*nlmax*ndham*nsites*nsp*nkp)
      else
        icls = 0
      endif

#if LMFGW
C --- GW driver ---
      call upack2('pot osgw nlma',spot,osgw,nchan)
      if (jobgw .eq. -999 ) then
        write(stdo,*) ' lmfgw: input one of the following jobs:'
        write(stdo,*) '  -1 : creates files',
     .    ' GWinput, QIBZ, KPTin1BZ'
        write(stdo,*) '   0 : init mode; creates files',
     .                ' SYMOPS, LATTC, CLASS, NLAindx, ldima'
        write(stdo,*) '   1 : GW setup mode'
        write(stdo,*) '   4 : band mode '
        write(stdo,*) '   5 : eigenvalue-only mode '
        write(stdo,*) ' job?'
        read (5,*) jobgw
      endif

C ... Mode -1 ... only after mkpot because qval made in mkpot
      if (jobgw .eq. -1) then
        call sugwin(ssite,sspec,slat,sham,w(osgw),sbz,nbas,qval)
        call rx0('lmfgw, job -1')
      endif

C ... Modes 0..4
      call sugw(ssite,sspec,slat,sham,nbas,ndham,
     .  w(osmpot),vconst,lcplxp,osig,otau,oppi,w(oppnl),vrmt,
     .  w(ospotx),osigx,otaux,oppix,jobgw)
      call rlse(oqmom)
C      if (jobgw .eq. 1) then
C        jobgw = 2
C        call info(20,1,0,' ... Remake potential without vxc',0,0)
C        goto 10
C      endif
      call tcx('bndfp')
      call rx0('bndfp')
#endif

C --- Start loop over k points; also, re-entry for second band pass ---
   99 continue
      ebot = 1000d0
      call surho(nbas,ssite,sspec,lmet,ldos,lrout,lekkl,numq,
     .  oqkkl,oeqkkl,k1,k2,k3,w(osrout),ndos,w(odos),sumev,sumqv)
      if (lfrce .gt. 0) then
        call dpzero(frc, 3*nbas*numq)
        call dpzero(w(ofh),3*nbas)
      endif

      if (lswtk .eq. 1) then
        efmax = 1d3
        nevmx = ndhamx
        call dpzero(w(oswtk),ndhamx*nkp)
      endif

C --- Setup moments file : write header ---
      nl = igets('ctrl nl',sctrl)
      nfilem = 0
      if (procid .eq. master) then
        nfilem = fopna('moms',-1,4)
      endif
C     if (cmdopt('--mull',6,0,strn)) nfilem = fopna('moms',-1,4)
C     if (cmdopt('--pdos',6,0,strn)) nfilem = fopna('moms',-1,4)
      if (nfilem .gt. 0) then
        rewind nfilem
        i = 1
C   ... If switch '--mull'; get sites, number of channels
        if (cmdopt('--mull',6,0,strn)) then
          mulopt = strn(7:)
          call upack('lat nsgrp',slat,ng,0,0,0,0)
          efmax = 1d3
          nevmx = ndhamx
          nchmx = min(1024,nbas*nl**2)
          call defi(ochan,nchmx)
          call sumlst(1,nchmx,nbas,ng,ssite,sspec,mulopt,moddos,
     .      nsites,lsites,lmxch,nchan,w(ochan),0)
          call redfi(ochan,nchan)
          call iomomn(.true.,2,.false.,1,nspc,1,1,i)
        endif
C   ... If switch '--pdos'; get sites, number of channels
        if (cmdopt('--pdos',6,0,strn)) then
          mulopt = strn(7:)
          call upack('lat nsgrp',slat,ng,0,0,0,0)
          nchmx = min(1024,nbas*nl**2)
          call defi(ochan,nchmx)
          nll = 0
          call sumlst(0,nchmx,nbas,ng,ssite,sspec,mulopt,moddos,
     .      nsites,lsites,lmxch,nchan,w,nll)
          nlmax = (lmxch+1)**2
          i = nlmax*ndham*3*nsp*nbas
          if (16*i*nkp/1000000 .gt. 24 .and. procid .eq. master)
     .    call info(20,0,0,' PDOS: %iMb memory for aus: nlmax=%i',
     .      i/1000000,nlmax)
          call defcc(oausp,i)
          call rlse(oausp)
          if (cmdopt('--mull',6,0,strn))
     .      call rx('--pdos and --mull not allowed in conjunction')
          if (cmdopt('--cls',5,0,strn))
     .      call rx('--pdos and --cls not allowed in conjunction')
          call iomomn(.true.,2,.false.,1,nspc,1,1,i)
        endif
        if (procid .eq. master) i = iomoms(-nfilem,nl,nsp,nspc,nkp,
     .    ndham,i,1,0,1,0,0,0,0,0,0d0,0d0,0d0,0d0,0d0,0d0)
      endif

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
        if (procid .eq. master) then
          call suqlst(plbopt,iopq,ndhamx,ef0,i,w,nfbn,w(oifbls),nkp,qp,
     .      onesp)
        endif
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
        call defrr(oqp,  3*nkp)
        call defrr(oevl, ndham*nsp*nkp)
        if (epsovl .ne. 0) then
          call defrr(oeomin,nsp*nkp)
        endif
        call info2(20,1,1,
     .    ' bndfp:  MPIK band plotting mode %i:  %i q-points',i,nkp)
C       Loop through all qp; accumulate vector of qp.
C       Use i2 in place of nkp to preserve nkp
        if (procid .eq. master) then
C       iq = running index to big qp list, i1 = index to current line
          iq = 0
          i2 = 0
  199     continue
          i = 1
          call pshpr(0)
          call suqlst(plbopt,1,ndhamx,ef0,i,w,nfbn,w(oifbls),i2,qp,
     .      onesp)
          call poppr
          if (i2 .gt. 0) then
            call pshpr(0)
            do  i1 = 1, i2
              iq = iq+1
              call suqlst(plbopt,1,ndhamx,ef0,i,w,nfbn,w(oifbls),i2,qp,
     .          onesp)
              call dpscop(qp,w(oqp),3,1,3*iq-2,1d0)
            enddo
            call poppr
            call suqlsm(i)
            if (i .ne. 3) goto 199
          endif
        endif
        call mpibc1(w(oqp),3*nkp,4,.false.,'bndfp','qp')
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
        call defrr(oorbtm, -nl*nsp*nbas)
      endif

C ... Setup for case sigma or evecs written to disk.
C     No integrated quantities accumulated.  qp read from file
C     Also sigma must be written in (iq,isp) order, opposite to the
C     order in which they would be generated here.
C     Requires second loop over (iq,isp) and filtering of isp
C     First pass should have onesp=0 and lwtkb=0
C     Next line is re-entry point for 2nd spin when writing sigma
   50 continue
      if (lwsig .ne. 0) then
        onesp = onesp + 1

        call info5(20,1,0,' BNDFP:  '//
     .    '%?#(n==4)#Write evals,evecs to file##%-1j'//
     .    '%?#(n==3)#Write LDA evals,evecs to file##%-1j'//
     .    '%?#(n==1|n==2)#Write sigm(LDA)##%-1j'//
     .    '%?#(n==5)#Write sigm(orb)##%-1j'//
     .    '%?#(n==-1)#Rotate sigm back to orbital basis##%-1j'//
     .    '%j for %i qp%?#n==2#, spin 2##',
     .    lwsig,nkp,onesp,0,0)

C       endif
      endif

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
      do iq = kpproc(procid), kpproc(procid+1)-1
        if (iq .eq. kpproc(procid)) then
          if (mlog) then
            call gettime(datim)
            call awrit4(' bndfp '//datim//' Process %i of %i on '
     .        //shortname(procid)(1:namelen(procid))//
     .        ' starting k-points %i to %i',' ',256,lgunit(3),
     .        procid,numprocs,kpproc(procid),kpproc(procid+1)-1)
          endif
        endif
#else
      nsmidb = ndham
      do  iq = 1, nkp
#endif
        isqp = nsp*(iq-1)
C       Get qp either from qp list or read from suqlst
#ifndef MPIK
        if (plbnd .ne. 0) then
          i = nsp
          if (onesp .ne. 0 .or. nspc .eq. 2) i = 1
          call suqlst(plbopt,0,ndhamx,ef0,i,w,nfbn,w(oifbls),nkp,qp,
     .      onesp)
        else
          call dpscop(w(oqp),qp,3,3*iq-2,1,1d0)
        endif
#else
        call dpscop(w(oqp),qp,3,3*iq-2,1,1d0)
#endif

C   ... Loop over majority, minority spins
C       In the SO case, still a loop isp=1..2
C       isp=1 generates just the 11 block
        call defrr(oh,-ndhamx**2*2)
        call defrr(os,-ndhamx**2*2)
        do  isp = 1, nsp
        if (onesp .eq. 0 .or. isp .eq. onesp) then
        if (isp .eq. 1 .or. isp .eq. onesp) call shorbz(qp,qp,qlat,plat)

C   ... For this qp, G vectors for PW basis and hamiltonian dimension
        if (pwemax .gt. 0 .and. mod(pwmode,10) .gt. 0) then
          pwgmin = dsqrt(pwemin)
          pwgmax = dsqrt(pwemax)
          call pshpr(1)
          call dpzero(xv,3)
          if (mod(pwmode/10,10) .eq. 1) call dpcopy(qp,xv,1,3,1d0)
          call gvlst2(alat,plat,xv,0,0,0,pwgmin,pwgmax,0,
     .      0,0,napw,dum,dum,dum,dum)
          call poppr
          call defi(oigv2,3*napw)
          call defi(okv,3*napw)
          call pshpr(iprint()-10)
          call gvlst2(alat,plat,xv,0,0,0,pwgmin,pwgmax,0,
     .      2,0,napw,w(okv),dum,dum,w(oigv2))
          call rlse(okv)
          call poppr
          ndimh = ldim + napw
          if (mod(pwmode,10) .eq. 2) ndimh = napw
          if (ndimh .gt. ndham) then
            call fexit2(-1,111,'%N Exit -1 : BNDFP: '//
     .        'ndimh=%i exceeds ndham=%i.  Try increasing '//
     .        'input NPWPAD',ndimh,ndham)
          endif
        else
          ndimh = ldim
          oigv2 = 1
        endif
        ndimhx = ndimh*nspc

C   ... Make Hamiltonian and overlap matrices
        call upack2('ham nqsig oqsig',sham,nqsig,oqsig)
        if (oqsig .eq. 0) oqsig = 1

        ispc = min(isp,nspc)
        i = lrsig*10
C       lwsig=-1: sigm LDA -> orbital basis
C       lwsig= 1: sigm orbital -> LDA basis
C       lwsig= 2: sigm orbital -> LDA basis, high energy parts replaced
C       lwsig= 3: Write evecs of LDA hamiltonian to file
C       lwsig= 4: Write evecs of hamiltonian to file
C       lwsig= 5: sigm orbital -> orbital, high energy parts replaced
C       hambls returns with evecs,evals of the LDA hamiltonian in s,h
C       Transform sigm to LDA basis: hambls returns sigm(LDA) in h
        if (lwsig .eq. 1) then
          i = i + 3000
        elseif (lwsig .eq. 2) then
          i = i + 4000
        elseif (lwsig .eq. -1 .or. lwsig .eq. 3) then
          i = i + 1000
          if (nspc .eq. 2)
     .      call rx('transform sigm not implemented in noncoll case')
        elseif (lwsig .eq. -1 .or. lwsig .eq. 4) then
          i = i + 6000
        elseif (lwsig .eq. 5) then
          i = i + 5000
        endif

        call hambls(i,nbas,ssite,sspec,slat,sham,isp,ispc,qp,k1,k2,k3,
     .    w(oqsig),nqsig,w(osmpot),vconst,osig,otau,oppi,lcplxp,lso,
     .    alfa,ndimh,napw,w(oigv2),w(oh),w(os),ismidb)
        nsmidb = min(nsmidb,ismidb)

C       call zprm('evec, LDA',2,w(os),ndimh,ndimh,ndimh)

C   ... Mode to write sigma or evecs to file
        if (lwsig .ne. 0) then
C         Write header information
          call upack('bz lshft',sbz,lshft,0,0,0,0)
          if (iq .eq. 1 .and. isp .eq. 1) then
            if (lwsig .eq. 1 .or. lwsig .eq. 2 .or. lwsig .eq. 5) then
              ifi = fopna('sigm2',-1,4)
            elseif (lwsig .eq. 3 .or. lwsig .eq. 4) then
              ifi = fopna('evec',-1,4)
            endif
            rewind ifi
            call iosigh(0,lwsig,nsp,ndimh,nk1,nk2,nk3,nkp,
     .        lshft(1),lshft(2),lshft(3),-ifi)
          endif

C         lwsig=1,2 : dump sigma(LDA basis) into file sigm2
C         lwsig=3,4 : dump LDA eigenvectors into file evec
C         lwsig=5   : dump sigma(orb basis) into file sigm2
          if (lwsig .ge. 1 .and. lwsig .le. 5) then
            if (lwsig .eq. 1 .or. lwsig .eq. 2 .or. lwsig .eq. 5) then
              ifi = fopna('sigm2',-1,4)
            else
              ifi = fopna('evec',-1,4)
            endif
            write(ifi) qp
            call dpdump(w(oh),ndimh,-ifi)
            call dpdump(w(os),ndimh**2*2,-ifi)

          elseif (lwsig .ne. 0) then
            call rxi('BNDFP: band mode lwsig=',lwsig)
          endif
          goto 30
        endif

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

C   ... debugging ... test of trotwf
#if TESTRWF
        call upack('ham ldham ooffH oindxo',sham,ldham,ooffH,oidxsh,0,0)
        call upack('lat plat oistab oag osymgr nsgrp',slat,plat,oistab,
     .    oag,og,nsgrp)
        call upack('bz oqp owtkp oipq ostar',sbz,oqp,owtkp,oipq,ogstar,
     .    0)
        call upack('bz nkp nkabc lshft',sbz,nkp,nkabc,lshft,0,0)
C   ... Make is,ifac,qb,qlat
        do  8  i = 1, 3
    8   llshft(i) = lshft(i) .ne. 0
        call pshpr(0)
        call bzmsh0(plat,llshft,0,nk1,nk2,nk3,is,ifac,rb,qb)
        call poppr

        rsrnge = dgets('ham rsrnge',sham)
        call trothf(nl,nbas,jsp,nsp,w(ooffH),w(oidxsh),w(oistab),nk1,
     .    nk2,nk3,w(oipq),nkp,w(oqp),w(owtkp),rsrnge,ndimh,0,slat,ssite,
     .    sspec,sham,w(osmpot),osig,otau,oppi,vconst,w(og),w(oag),nsgrp,
     .    w(ogstar),lshft,ifac,qb)
C        call trotwf(ssite,sspec,slat,sham,w(osmpot),vconst,osig,otau,
C     .    oppi,nl,nbas,jsp,w(ooffH),w(oidxsh),w(oistab),nk1,nk2,nk3,
C     .    w(oipq),nkp,w(oqp),ndimh,plat,w(og),w(oag),w(ogstar),
C     .    ifac,qb)
        call rx('done')
#endif

C   ... Diagonalize and add to density
        if (lpdiag .eq. 0) then
          call defcc (ot,ndimhx**2)
        elseif (lpdiag .eq. 2) then
          ot  = os
        endif
        if (mod(iq,10) .ne. 1) call pshpr(iprint()-6)
#ifndef MPIK
        call info5(30,0,0,' bndfp:  kpt %i of %i, k=%3:2,5;5d'//
     .    '%?#n#   ndimh = %i##',
     .    iq,nkp,qp,mod(pwmode/10,10),ndimh)
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
        if (lpdiag .eq. 1) then
          call rxx(nspc.ne.1,'parallel diag not implemented for noncol')
#if MPE & MPI
          ierr = MPE_LOG_EVENT(EVENT_START_PZHEV,procid,"pzhev")
#endif
          if (iq .eq. 1) then
            call info0(20,0,0,' bndfp:  diagonalise with SCALALPACK ..')
          endif
          call pzhev(T,ndimh,oh,os,nblk,nprow,npcol,efmax,nmx,nev,
     .               evl(1,jsp),ot)
          nevl = ndimh
#if MPE & MPI
          ierr = MPE_LOG_EVENT(EVENT_END_PZHEV,procid,"pzhev")
#endif
        elseif (lpdiag .eq. 2) then
          nev = ndimhx
          nmx = ndimhx
          nevl = ndimhx
          call phmbls(2,ndimhx,evl(1,jsp),w,w,w,w,w,w(oh))
          if (lwtkb .ne. -1) then
            call defrr(oww,  ndimhx**2*2)
            if (ldim .ne. ndham) call
     .        rx('blsig needs separate dim for wtkp,evec')
            call blsig(1+lrsig*10,nbas,sham,isp,nsp,nspc,plat,qp,
     .        lwtkb,qval-qbg,iq,w(owtkp),w(owtkb),ndimh,w(ot),w(oww))
            call rlse(oww)
          endif

        else
          if (nspc .eq. 2) then
            call defrr(oww,  ndimhx**2)
            call sopert(0,ndimh,nspc,w(oww),w(oh),w(oh))
            call sopert(0,ndimh,nspc,w(oww),w(os),w(os))
            call rlse(oww)
          endif
          call defrr(oww,  11*ndimhx)
c takao
c          if (epsovl .eq. 0) then
          if(napw == 0) then
            call zhev(ndimhx,w(oh),w(os),T,T,nmx,efmax,nev,w(oww),F,-1,
     .        evl(1,jsp),w(ot))
            nevl = ndimhx
          else
            nevl = -1
            call dvset(w(oww),1,1,99999d0)
c            call zhevo(ndimhx,ndimhx,w(oh),w(os),T,nmx,efmax,epsovl,
c     .        nevl,nev,evl(1,jsp),w(oww),ndimhx,w(ot))
c takao version. Not in lm-7.0beta
            call zhev_tk(ndimhx-napw,ndimhx,w(oh),w(os),nmx,epsovl,
     .        ncutovl,nevl,nev,evl(1,jsp),w(oww),ndimhx,w(ot),iq==1)

            eomin = dval(w(oww),1)
          endif
          call rlse(oww)
        endif
C       Pad evals between ndimh and ndham with a large positive number
C       to avoid mixing up integration routines
        if (ndhamx .gt. nevl .and. nspc .eq. 2) then
          call dvset(evl,1+nevl,ndhamx,99999d0)
        elseif (ndham*nspc .gt. nevl) then
          call dvset(evl(1,jsp),1+nevl,ndham,99999d0)
        endif

#ifndef MPIK
        if (epsovl .ne. 0) then
          call info5(30,0,0,' Overlap''s smallest eigenvalue: %;3g.  '//
     .    '%?#(n>0)#H dim reduced from %i to %i#H dim not reduced#',
     .      eomin,ndimhx-nevl,ndimhx,nevl,0)
        endif
        call prtev(w(ot),nevl,evl(1,jsp),nmx,efmax,nev)
        if (iprint() .ge. 110) then
          call yprm('evals',1,evl(1,jsp),1,ndhamx,nevl,1)
          if (nev .gt. 0) call zprm('evecs',2,w(ot),ndimhx,ndimhx,nev)
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
     .  write (stdl,712) (evl(i,jsp),i=1,nev)
  712   format('fp evl',8f8.4)
#endif
        if (lwtkb .ne. -1 .and. .not.lwndow) then
          if (iq .eq. 1 .and. jsp .eq. nsp .and.
     .      .not. cmdopt('--no-fixef0',11,0,strn)) then
            ef00 = ef0
            call fixef0(qval-qbg,jsp,1,nevl,ndham,evl,dosw,ef0)
            if (jsp .eq. 2 .and. ef00 .ne. ef0 .and.
     .        lwtkb .eq. 0 .and. lmet .gt. 0 .and. lrout .ne. 0) then
              if (procid .eq. master)
     .          call info0(10,1,1,
     .          ' ... Fermi level reset in second spin'//
     .          ' channel ... restart band pass')
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
     .        call awrit3('%N evl(nev=%i)=%;3d but '//
     .        'ef0=%;3d ... restart with larger efmax or nevmx',
     .        ' ',80,stdo,nev,evl(max(nev,1),jsp),ef0)
            call rx('bndfp')
          endif
        endif
        endif
#endif
c##### !MPIK

C   ... Save data for this qp to moments file
        if (plbnd .eq. 0) then
C   ... Copy eigenvalues into array containing ev for all qp
        call dpscop(evl(1,jsp),w(oevl),ndhamx,
     .              1,1+ndham*(jsp-1+isqp),1d0)
#if MPIK
        if (epsovl .ne. 0) then
          call dpscop(eomin,w(oeomin),1,1,1+(jsp-1+isqp),1d0)
        endif
#endif

C        call prmx('ev',w(oevl),ndham,ndham,isqp+2)
C        call prmx('ev',w(oevl),ndham,ndham,nkp*nsp)

C   --- Orbital magnetic moment (requires lso) ---
        if (lso .ne. 0 .and. lwtkb .ne. -1) then
          if (lwtkb .eq. 0) call
     .      rx('metal weights required for orb. moment calculation')
          call defcc(oauso,-nlmax*ndhamx*3*nsp*nbas)
          call makusq(0,ssite,sspec,slat,sham,nbas,nbas,0,nlmax,ndham,
     .      ndimh,napw,w(oigv2),nev,nsp,nspc,jsp,1,qp,w(ot),w(oppnl),
     .      w(oauso))
          call mkorbm(ssite,sspec,jsp,nsp,nspc,nlmax,ndham,nev,w(owtkb),
     .      iq,nbas,w(oppnl),w(oauso),nl,nkp,w(oorbtm))
          call rlse(oauso)
        endif

C   --- Mulliken analysis and partial DOS ---
        if (lwtkb .ne. -1) then
        if (procid .eq. master) then
        if (cmdopt('--mull',6,0,strn)) then
#if MPIK
          call rx('MULL not k-parallel.')
#endif
          call rxx(nev.ne.ndimhx,'Mulliken requires all eigenvectors')
          call rxx(nevmx.ne.ndhamx,'Mulliken requires nevmx=nhdamx')
          call defrr(odoswt,-nchan*nevmx*nspc)
          call upack('ham oindxo',sham,oiprmb,0,0,0,0)
          call mullmf(nbas,ssite,sspec,w(oiprmb),w(ot),ndimh,nspc,
     .                iq,jsp,moddos,nsites,lsites,lmxch,nchan,w(ochan),
     .                lmdim,ndham,w(odoswt))
          call iomomn(.true.,2,.false.,1,nspc,1,1,i)
          i1 = iomoms(-nfilem,nl,nsp,nspc,nkp,ndimh,i,nspc,1,1,ndhamx,
     .      nevmx,nchan,nchan,nev,evl(1,jsp),0d0,w(odoswt),0d0,0d0,0d0)
          call rlse(odoswt)
        elseif (cmdopt('--pdos',6,0,strn)) then
#if MPIK
          call rx('PDOS not k-parallel.')
#endif
          call rxx(nspc.ne.1,'pdos not implemented in noncoll case')
          i = nlmax*ndham*3*nsp*nbas
          call defcc(oausp,-i)
          call makusq(1,ssite,sspec,slat,sham,nbas,nbas,0,nlmax,ndham,
     .      ndimh,napw,w(oigv2),nev,nsp,nspc,jsp,1,qp,w(ot),w(oppnl),
     .      w(oausp))
          call defrr(odoswt,-nchan*ndimhx*nspc)
          call mkpdos(moddos,ssite,sspec,jsp,nsp,nlmax,ndham,nev,
     .      nchan,1,lsites,nsites,w(oppnl),w(oausp),w(odoswt))
          call iomomn(.true.,2,.false.,1,nspc,1,1,i)
#if MPIK
          i1 = iomoms(-nfilem,nl,nsp,nspc,2,ndimh,1000+i,1,iq,isp,
     .      ndham,ndimh,nchan,nchan,nev,evl(1,jsp),0d0,w(odoswt),0d0,
     .      0d0,0d0)
#else
          i1 = iomoms(-nfilem,nl,nsp,nspc,2,ndimh,i,1,1,1,ndham,ndimh,
     .      nchan,nchan,nev,evl(1,jsp),0d0,w(odoswt),0d0,0d0,0d0)
#endif
          call rlse(odoswt)
          call rlse(oausp)
        elseif (nfilem .gt. 0) then
          write (nfilem) 0, ndimhx
          call dpdump(evl(1,jsp),ndimhx,-nfilem)
        endif
        endif

C   ... Make new density matrix dmatu for LDA+U
        if (nlibu .gt. 0 .and. nev .gt. 0) then
          if (lwtkb .eq. 0)
     .      call rx('metal weights required for LDA+U calculation')
          nl = igets('ctrl nl',sctrl)
          nlmax = nl*nl
          call defrr(oausu,-2*nlmax*ndhamx*3*nsp*nbas)
          call makusq(0,ssite,sspec,slat,sham,nbas,nbas,0,nlmax,ndham,
     .      ndimh,napw,w(oigv2),nev,nsp,nspc,jsp,1,qp,w(ot),w(oppnl),
     .      w(oausu))
          call mkdmtu(ssite,sspec,w(owtkb),jsp,iq,nsp,nspc,ndham,nlmax,
     .      nbas,nev,w(oppnl),w(oausu),dmatu,nlibu,lmaxu,lldau)
          call rlse(oausu)
        endif
        endif

C   ... for core-level spectroscopy
        if (icls .ne. 0 .and. lwtkb .ne. -1) then
          call rxx(nspc.ne.1,'CLS not implemented in noncoll case')
#if MPI
            call rx('CLS only k-parallel')
#endif
          call makusq(0,ssite,sspec,slat,sham,nbas,nsites,isite,nlmax,
     .      ndham,ndimh,napw,w(oigv2),nev,nsp,nspc,jsp,iq,qp,w(ot),
     .      w(oppnl),w(oausc))

C       Accumulate output density and sampling DOS
C       Even if no output rho, still call addrbl to make DOS when lmet=4
        elseif (lwtkb.ne.-1 .and. (lrout.ne.0 .or. lmet.eq.4)) then
C         dum=frc(3,1)
          call addrbl(ssite,sspec,slat,sham,jsp,nsp,nspc,qp,w(owtkp),
     .      ndham,ndimh,napw,w(oigv2),lmet,lrout,lwtkb,w(owtkb),lswtk,
     .      w(oswtk),iq,lfrce,ldos,lekkl,k1,k2,k3,w(osmpot),vconst,osig,
     .      otau,oppi,lcplxp,numq,qval-qbg,w(ot),evl,nevl,ef0,def,
     .      esmear,emin,emax,ndos,w(odos),w(osrout),sumqv,sumev,oqkkl,
     .      oeqkkl,frc)
C          print 399,iq,frc(3,1),frc(3,1)-dum
C  399     format(' after addrbl: frc(3,1)=',i4,2f12.6)
        endif
#ifndef MPIK
        elseif (plbnd .ne. 0) then
          i = nsp
          if (onesp .ne. 0 .or. nspc .eq. 2) i = 1
          call suqlsw(nevl,jsp,i,evl(1,jsp))
          if (nfbn(1) .ne. 0) then
            if (ndimhx .ne. nevl)
     .        call rx('color weights not implemented when '//
     .        'nevl < hamiltonian dimension')
            call suqlse(ndimhx,jsp,i,ndimhx,1,nfbn,w(oifbls),ndhamx,
     .        w(ot),evl)
          endif
          if (nfbn(2) .ne. 0) then
            if (ndimhx .ne. nevl)
     .        call rx('color weights not implemented when '//
     .        'nevl < hamiltonian dimension')
            call suqlse(ndimhx,jsp,i,ndimhx,2,nfbn,w(oifbls),ndhamx,
     .        w(ot),evl)
          endif
#else
C   ... Save evals for this qp
        elseif (plbnd .ne. 0) then
          call dpscop(evl(1,jsp),w(oevl),ndhamx,
     .                1,1+ndham*(jsp-1+isqp),1d0)
          if (epsovl .ne. 0) then
            call dpscop(eomin,w(oeomin),1,1,1+(jsp-1+isqp),1d0)
          endif
#endif

        endif

        if (ot .ne. os) call rlse(ot)

C end loop over isp (main loop in parallel mode)
   30   continue
        endif
        enddo
        call rlse(oh)
C end loop over iq (main loop in parallel mode)
      enddo

C Second pass over iq for second spin when writing sigma
      if (lwsig .ne. 0) then
        if (onesp .lt. nsp) goto 50
        if (lwsig .eq. 1 .or. lwsig .eq. 2 .or. lwsig .eq. 5) then
          ifi = fopna('sigm2',-1,4)
          call fclose(ifi)
          if (lwsig .eq. 5)
     .    call rx0('BNDFP:  sigm(orb basis) saved in file sigm2')
          call rx0('BNDFP:  sigm(LDA basis) saved in file sigm2')
        else if (lwsig .eq. 3 .or. lwsig .eq. 4) then
          ifi = fopna('evec',-1,4)
          call fclose(ifi)
          call rx0('BNDFP:  LDA evecs saved in file evec')
        endif
C        if (lwsig .eq. -1) then
C        ifi = fopna('sigm3',-1,4)
C        call fclose(ifi)
C        call rx0('BNDFP:  sigm(orbital basis) saved in file sigm3')
C        endif
C        if (lwsig .eq. 3) then
C        call rx0('BNDFP:  U(LDA-QP) saved in file sigm2')
C        endif
C        if (lwsig .eq. 4) then
C        call rx0('BNDFP:  H(QP) saved in file sigm2')
C        endif
      endif

C ... Set switch to write file sigii
#ifndef LMFGW
      if (lrsig .ne. 0 .and. plbnd .eq. 0 .and. procid .eq. master) then
        call upack1('ham sigp',sham,sigp)
        if (mpipid(0) .le. 1) call phmbl3(12,0,nsmidb,0,0,sigp,qp,qp(2))
        call fclose(fopna('sigii',-1,0))
      endif
#endif

C --- Second loop over qp, needed to make k-parallelisation possible: ---
C     You can't do this until you have all the evals.
#if MPIK
#if MPE
      ierr = MPE_LOG_EVENT(EVENT_END_KLOOP,procid,"k-loop")
#endif
      entime = MPI_WTIME()
      call info2(30,0,0, ' ... Done MPI k-loop: %;1d seconds elapsed',
     .  entime-sttime,0)
      call info0(20,0,-1,' ... Sharing data between processes...')
      sttime = MPI_WTIME()

C     Share bands among all processors
      call xmpbnd(kpproc,ndham,nkp,nsp,w(oevl))
      if (epsovl .ne. 0) then
        call xmpbnd(kpproc,1,nkp,nsp,w(oeomin))
      endif
      if (icls .ne. 0 .and. lwtkb .ne. -1) then
        call rxx(nspc.ne.1,'CLS not implemented in noncoll case')
        call xmpbnd(kpproc,nlmax*ndham*3*nsites,nkp,nsp,w(oausc))
      endif

ckino 05Feb09, use kpproc later, move it.
c      deallocate(kpproc, stat=ierr)

C     Allreduce density-related quantities
      if (lrout .ne. 0) then
      call mpibc2(sumqv,6,4,mlog,'bndfp','sumqv')
      call mpibc2(sumev,6,4,mlog,'bndfp','sumev')
      call mpibc2(w(osrout),k1*k2*k3*nsp*numq,6,mlog,'bndfp','smrho')
      if (lswtk .eq. 1) then
        call mpibc2(w(oswtk),ndhamx*nkp,4,mlog,'bndfp','swtk')
      endif
C     Allreduce qkkl
      nkaph = nglob('nkaph')
      do  ib = 1, nbas
        is = igetss('site spec',ib,ssite)
        call upack('spec lmxa lmxb kmxt',sspec,is,lmxa,lmxh,kmax,0)
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
     .      call mpibc2(w(oeqkkl(i,ib)),nelt(i)*numq*nsp*nspc,4,
     .        mlog,'bndfp','eqkkl')
          enddo
        endif
      enddo
C     Allreduce DOS, forces, dmatu, orbtm
      if (ndos .gt. 0)
     .  call mpibc2(w(odos),ndos*2*nsp,4,mlog,'bndfp','dos')
      if (lfrce .ne. 0)
     .  call mpibc2(frc,3*nbas*numq,4,mlog,'bndfp','frc')
      if (nlibu .gt. 0)
     .call mpibc2(dmatu,nsp*nlibu*(lmaxu*2+1)**2,6,mlog,'bndfp','dmatu')
      if (lso .ne. 0 .and. lwtkb .ne. -1)
     .  call mpibc2(w(oorbtm),nl*nsp*nbas,4,mlog,'bndfp','orbtm')
      endif
      if (lpdiag .eq. 2 .and. plbnd .eq. 0) then
        call upack1('ham eterms',sham,eterms)
        call mpibc2(eterms(19),1,4,mlog,'bndfp','rhosig')
        call pack1('ham eterms',sham,eterms)
      endif

      entime = MPI_WTIME()
      call info2(20,0,0,' MPI broadcast took %;1d sec',
     .  (entime-sttime),0)

C     Write bands in bands-plotting case: loop over qp getting evals from array
      if (plbnd .ne. 0) then
      call info0(20,1,0,' Writing bands to bands file (MPIK) ...')
      if (procid .eq. master) then
C     iq = running index to big qp list, i1 = index to current line
      iq = 0
      i2 = 0
  299 continue
      i = nsp
      if (onesp .ne. 0 .or. nspc .eq. 2) i = 1
      call suqlst(plbopt,0,ndhamx,ef0,i,w,nfbn,w(oifbls),i2,qp,onesp)
      if (i2 .le. 0) call rx0('bndfp')
      do  i1 = 1, i2
        iq = iq+1
        isqp = nsp*(iq-1)
        call suqlst(plbopt,1,ndhamx,ef0,i,w,nfbn,w(oifbls),i2,qp,onesp)
        call dpscop(qp,w(oqp),3,1,3*iq-2,1d0)
        do  isp = 1, nsp
          ispc = min(isp,nspc)
C         jsp=isp in the collinear case; jsp=1 in the noncollinear
C         Thus jsp should be used in place of isp
C         isp serves as a flag for the noncollinear case
          if ((onesp.eq.0 .or. isp.eq.onesp) .and. (ispc.eq.nspc)) then
            jsp = isp
            if (ispc .eq. 2) jsp = 1
            call dpscop(w(oevl),evl(1,jsp),ndhamx,
     .        1+ndham*(jsp-1+isqp),1,1d0)
            if (mod(i1,10) .ne. 1) call pshpr(iprint()-6)
            call info5(30,0,0,' bndfp:  kpt %i of %i, k=%3:2,5;5d',
     .        i1,nkp,qp,0,0)
            if (mod(i1,10) .ne. 1) call poppr
            i = nsp
            if (onesp .ne. 0 .or. nspc .eq. 2) i = 1
            call suqlsw(ndimhx,jsp,i,evl(1,jsp))
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
        call dpscop(w(oqp),qp,3,3*iq-2,1,1d0)
        call shorbz(qp,qp,qlat,plat)
        do  isp = 1, nsp
        ispc = min(isp,nspc)
C       jsp=isp in the collinear case; jsp=1 in the noncollinear
C       Thus jsp should be used in place of isp
C       isp serves as a flag for the noncollinear case
        if ((onesp.eq.0 .or. isp.eq.onesp) .and. (ispc.eq.nspc)) then
        jsp = isp
        if (ispc .eq. 2) jsp = 1
        call dpscop(w(oevl),evl(1,jsp),ndhamx,
     .    1+ndham*(jsp-1+isqp),1,1d0)
        if (epsovl .ne. 0) then
          call dpscop(w(oeomin),eomin,1,1+(jsp-1+isqp),1,1d0)
        endif
        if (mod(iq,10) .ne. 1) call pshpr(iprint()-6)
        call info5(30,0,0,' bndfp:  kpt %i of %i, k=%3:2,5;5d',
     .    iq,nkp,qp,0,0)
        if (epsovl .ne. 0) then
          call info5(30,0,0,' Overlap''s smallest eigenvalue: %;3g',
     .      eomin,0,0,0,0)
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
     .    write (stdl,712) (evl(i,jsp),i=1,nev_iq)
ckino 05Feb09     .    write (stdl,712) (evl(i,jsp),i=1,nev)
  712   format('fp evl',8f8.4)
        if (lwtkb .ne. -1 .and. .not. lwndow) then
          if (iq .eq. 1 .and. jsp .eq. nsp .and.
     .      .not. cmdopt('--no-fixef0',11,0,strn)) then
            ef00 = ef0
            call fixef0(qval-qbg,jsp,1,ndimh,ndham,evl,dosw,ef0)
            if (jsp .eq. 2 .and. ef00 .ne. ef0 .and.
     .        lwtkb .eq. 0 .and. lmet .gt. 0 .and. lrout .ne. 0) then
              if (procid .eq. master) call info0(10,1,1,
     .          ' ... Fermi level reset in second spin'//
     .          ' channel ... restart band pass')
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
     .        'ef0=%;3d ... restart with larger efmax or nevmx',
     .        ' ',80,stdo,nev_iq,evl(max(nev_iq,1),jsp),ef0)
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
c###### second loop over qp (parallel k-points mode)

C     For now, SO by site, not by class
      if (lwtkb .eq. 1 .and. lso .ne. 0) then
        call defi(oips,nbas)
        call spackv(10,'site spec',ssite,1,nbas,w(oips))
        call iorbtm(sspec,w(oips),nl,nl,nbas,nsp,w(oorbtm))
        call rlse(oips)
      endif

C ... Case generating bands: find next block of qp
      if (plbnd .ne. 0) goto 99

      if (ipl.gt.0) write (stdl,715) nkp,ebot,qval,qbg,esmear
  715 format('nv nkp',i5,'  ebot',f9.4,'   qval',f10.4,'  qbg',f8.4,
     .   '  smr',f8.4)
C ... End of k point loop

C     call zprm3('smrho after k-point loop',0,w(osrout),k1,k2,k3)

C --- Interpolate density to Fermi energy ---
      sev = sumev(1,1)
      if (lmet .eq. 4) then
        call mshn3p(nbas,ssite,sspec,lmet,lrout,lfrce,qval-qbg,ef0,def,
     .    sumqv,sumev,n1,n2,n3,k1,k2,k3,w(osrout),oqkkl,frc,lrep)
C   ... Store val q & magnetic moment in sumqv(1) and sumqv(2)
        sumqv(2,1) = sumqv(1,1) - sumqv(1,2)
        sumqv(1,1) = sumqv(1,1) + sumqv(1,2)
C   ... Eigenvalue sum including entropy term
        sev = sumev(1,1) + sumev(2,1)
C   ... Remake sev,ef linearly interpolating tabulated sampling DOS
        sev00 = sev
        ef00  = ef0
        if (ldos .ne. 0) then
          call efldos(qval,nsp,emin,emax,ndos,w(odos),eferm,sev1)
          sev   = sev1
          ef0  = eferm
        endif
        if (ipr .gt. 30 .and. ldos .ne. 0)
     .  write(stdo,388) sev00,ef00,sev1,eferm,sev,ef0
  388   format(' ipol:  sev=',f12.6,'   ef=',f12.6:/
     .         ' dos:   sev=',f12.6,'   ef=',f12.6/
     .         ' use:   sev=',f12.6,'   ef=',f12.6)
        if (ipl.gt.0) write (stdl,733) ef00,eferm,ef0,sev00,sev1,sev
  733   format('nf EF:',3f9.5,'    EB:',3f12.5)
        call pack5('bz ndos dosw ef def',sbz,ndos,dosw,ef0,def,0)
        if (lrep .eq. 1) then
          ef0 = -1
          call awrit0('Input Fermi energy was too far off, '//
     .      'repeat band pass',' ',80,stdo)
          goto 99
        endif
      endif

C --- BZ integration for fermi level, band sum and qp weights ---
      if (lmet .ge. 0 .and. (lmet .ne. 4 .or. ltet)) then
        if (lwndow) then
          call defrr(oww,ndham*nsp*nkp)
          eferm = min(dosw(1),dosw(2))
          call bzints(nkabc(1),nkabc(2),nkabc(3),w(oevl),w(oww),nkp,
     .      ndham,ndham,nsp,xv,xv,xv,1,eferm,2,ntet,w(oidtet),sev,dum)
C         call prmx('w(min)',w(oww),ndham,ndham,nkp*nsp)
          eferm = max(dosw(1),dosw(2))
          call bzints(nkabc(1),nkabc(2),nkabc(3),w(oevl),w(owtkb),nkp,
     .      ndham,ndham,nsp,xv,xv,xv,1,eferm,2,ntet,w(oidtet),sev,dum)
C         call prmx('w(max)',w(owtkb),ndham,ndham,nkp*nsp)
          call daxpy(ndham*nsp*nkp,-1d0,w(oww),1,w(owtkb),1)
C         call prmx('w',w(owtkb),ndham,ndham,nkp*nsp)
          call rlse(oww)
        else
          dosrng = 8
          if (mpsord .lt. 0) dosrng = 16
          call bzwtsf(ndham,ndham,nsp,nspc,nkabc(1),nkabc(2),nkabc(3),
     .      nkp,ntet,w(oidtet),qval-qbg,fsmom,lmet.ne.0,ltet,mpsord,
     .      ndos,dabs(esmear-mpsord),dosrng,w(owtkp),w(oevl),efmax,
     .      lswtk,w(oswtk),eferm,sev,w(owtkb),sumqv(1,2),lwtkb)
C         Store val charge & magnetic moment in sumqv(1..2)
          if (lmet .ne. 4) then
            sumqv(1,1) = sumqv(1,2)
            sumqv(2,1) = sumqv(2,2)
          endif
        endif
        if (lmet .ne. 4) ef0 = eferm
        if (lmet .ne. 4) call pack1('bz ef',sbz,ef0)
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

C ... Save Fermi level, nonmetal or sampling integration
      if (lmet .eq. 0 .or. .not. (lmet .ne. 4 .or. ltet)) then
        if (lmet .eq. 0) then
          ef0 = (evtop+ecbot)/2
          if (ipr .ge. 20) call awrit3(' Highest occ. level = %,5;5d '//
     .      ' Lowest unocc. = %,5;5d  average = %,5;5d',' ',80,
     .      stdo,evtop,ecbot,ef0)
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
     .    ndham,ndimh,nchan,nchan,ndimh,0d0,0d0,0d0,0d0,eferm,0d0)
        call fclose(nfilem)
        if (iprint() .ge. 10 .and. cmdopt('--mull',6,0,strn)) then
          call mchan(lmdim,ssite,sspec,nsp,nsites,lsites,0,0,0,0,
     .               w(ochan))
        endif
      endif

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
      call upack('bz ndos',sbz,ndos0,0,0,0,0)
      if (lgors('ctrl ldos,1',sctrl)) then
        call defdr(odos,3*ndos)
        if (iprint() .ge. 30) call awrit1('%x%N ... Generating %?#n<0#'
     .    //'integrated#total# DOS',' ',80,lgunit(1),ndos0)
        if (ltet) then
          call bzints(nkabc(1),nkabc(2),nkabc(3),w(oevl),dum,nkp,ndhamx,
     .      ndhamx,nspx,dosw(1),dosw(2),w(odos),ndos,eferm,1,ntet,
     .      w(oidtet),dum,dum)
          if (ndos0 .gt. 0)
     .      call xxxdif(dosw(1),dosw(2),ndos,nspx,0,w(odos))
          del = 0d0
        else
          if (mpsord .ge. 100) mpsord = mod(mpsord,100)
          if (ndos0 .gt. 0)
     .      call makdos(nkp,ndhamx,ndhamx,nspx,w(owtkp),w(oevl),mpsord,
     .      dabs(esmear-mpsord),-6d0,dosw(1),dosw(2),ndos,w(odos))
          if (ndos0 .lt. 0)
     .      call maknos(nkp,ndhamx,ndhamx,nspx,w(owtkp),w(oevl),mpsord,
     .      dabs(esmear-mpsord),-6d0,dosw(1),dosw(2),ndos,w(odos))
          del = mpsord+dabs(esmear-mpsord)
        endif
        if (nspc .eq. 2) call dscal(ndos,.5d0,w(odos),1)
        call dosio(w(odos),ndos,nspx,ndos,1,dosw(1),dosw(2),nspx,
     .    eferm,del,1,-fopn('DOS'))
        call fclose(fopn('DOS'))
        call rlse(odos)
      endif
      endif

C --- Core-level spectroscopy ---
      if (icls .ne. 0) then
        if (procid .eq. master) then
          call upack1('bz ef',sbz,eferm)
          call vcdmel(sctrl,ssite,sspec,slat,nlmax,ndham,ndimh,nkp,nsp,
     .      nspc,eferm,w(oevl),w(oausc),nsites,isite,iclsl,iclsn)
          call rx0('done generating core level spectra')
        else
          call fexit(0,0,' ',0)
        endif
      endif

C --- Assemble output density ---
      call dfratm(nbas,ssite,sspec,orhat1)
      call defrr(oqbyl,  n0*nbas*nsp)
      call defrr(ohbyl,  n0*nbas*nsp)
C     --window: Put output density into orhat and smrho, and exit
      if (lwndow) then
        call mkrout(ssite,sspec,slat,sham,nbas,nsp,ldim,lekkl,oqkkl,
     .    oeqkkl,w(oorhat),w(ohab),w(osab),w(oqbyl),w(ohbyl),lrout)
        call zcopy(k1*k2*k3*nsp,w(osrout),1,w(osmrho),1)
        call symrho(ssite,sspec,slat,lfrce,w(osmrho),w(oorhat),w(oqbyl),
     .    w(ohbyl),frc)
        goto 999
      endif
      call mkrout(ssite,sspec,slat,sham,nbas,nsp,ldim,lekkl,oqkkl,
     .  oeqkkl,orhat1,w(ohab),w(osab),w(oqbyl),w(ohbyl),lrout)
C
      if (lrout .ne. 0) then
C   ... Symmetrize output density and forces
        call symrho(ssite,sspec,slat,lfrce,w(osrout),orhat1,w(oqbyl),
     .    w(ohbyl),frc)

C   ... New boundary conditions pnu for phi and phidot
C       call pshpr(iprint()-12)
        if (lpnu .gt. 0) then
          call upack1('ham pmin',sham,xv)
          call upack1('ham pmax',sham,xv(11))
          call pnunew(nbas,nsp,ssite,sspec,xv,xv(11),lfrzw,w(ohab),
     .      w(osab),w(oqbyl),w(ohbyl))
        endif
C       call poppr
      endif

C --- Evaluate Harris energy ---
      call mkehkf(1,sham,sev,sumqv(2,1),sumtv,ehar)

C --- Evaluate KS total energy, correction to Harris force ---
      if (lrout .ne. 0) then

C   ... Correction to harris force
        call dfrce(ssite,sspec,slat,sctrl,k1,k2,k3,nvl,w(oorhat),orhat1,
     .    elind,w(oqmom),w(osmrho),w(osrout),w(ofh))

C   ... Evaluate KS total energy and output magnetic moment
        eks = 0d0
        if (leks .ge. 1) then
          call togpr()
          call mkekin(nbas,ldim,ssite,sspec,slat,sham,osig,otau,oppi,
     .      lcplxp,oqkkl,k1,k2,k3,vconst,w(osmpot),w(osrout),sev,sumtv)
          call pshpr(ipr-20)
          call mkpot(nbas,ssite,sspec,slat,sham,spot,lfrce,lcplxp,
     .      k1,k2,k3,w(osrout),orhat1,qbg,w(osmpot),w(oqmom),vconst,
     .      vrmt,osig,otau,oppi,w(oppnl),w(ohab),w(ovab),w(osab),
     .      qval,qsc,w(ogpot0),w(ovval),w(ofes2),0,vorb,nlibu,lmaxu,
     .      lldau)
          call poppr()
          call mkehkf(2,sham,sev,sumqv(2,1),sumtv,eks)
          call togpr()
        endif

C   --- Add together force terms ---
        if (lfrce .gt. 0)
     .    call totfrc(nbas,ssite,slat,leks,w(ofes1),w(ofes2),w(ofh),frc)

C   --- Mix input and output densities ---
        call defrr(oqm1, nvl)
        call defrr(oqm2, nvl)
        call upacks('strn mix',i1,i2)
#if MPE
        ierr = MPE_LOG_EVENT(EVENT_START_MIXRHO,procid,"mixrho")
#endif
        call mixrho(ssite,sspec,slat,nsp,iter,sstrn(i1:i2),qval-qbg,
     .    elind,orhat1,w(oorhat),k1,k2,k3,dmxp,w(osrout),w(osmrho))
#if MPE
        ierr = MPE_LOG_EVENT(EVENT_END_MIXRHO,procid,"mixrho")
#endif
      else
        eks = 0
      endif
      call pack2('ham ehf ehk',sham,ehar,eks)

  999 continue
      call rlse(oqmom)
      call tcx('bndfp')

      call xxxbfp ! Patch to avoid PGI compiler bug on AMD processor
      end
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
      integer obuf
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
      call defrr(obuf,ndham*nkp*nsp)
      call MPI_ALLGATHERV(eb(1,1+nsp*(ista-1)),
     ,  length(procid),MPI_DOUBLE_PRECISION,w(obuf),length,
     ,  offset,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
      call dcopy(ndham*nsp*nkp,w(obuf),1,eb,1)
      call rlse(obuf)

      deallocate(offset, stat=ierr)
      deallocate(length, stat=ierr)

C      if (procid .eq. master) then
C        print *, procid,eb(1,1),eb(1,nsp*nkp)
C      endif
C      call rx('done')
      end
#endif
      subroutine xxxbfp
      end

