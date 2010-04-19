      subroutine bndasa(sctrl,slat,sspec,ssite,sarray,sbz,sham,spot,
     .  sstrn,vorb,dmatu,efermi,eband,nbmax,nevmx,qnu,
     .  sumev,rhos,amag,aamom)
C- k-integration of bands over BZ for moments.
C ----------------------------------------------------------------------
Ci Inputs
Ci   sctrl :struct containing parameters governing program flow
Ci     Elts read: nbas nclass nl nspec nspin lncol lsx lham lasa lgen3
Ci                loptc ldos
Ci     Stored:    lham
Ci     Passed to: nmpot secmat asaddq asaopm optint
Ci
Ci   slat  :struct containing information about the lattice
Ci     Elts read: avw alat plat ocg ojcg oidxcg nsgrp oistab oag osymgr
Ci                vol
Ci     Passed to: secmat asaopm
Ci
Ci   sspec :struct containing species-specific information
Ci     Passed to: nmpot secmat iorbtm
Ci
Ci   ssite :struct containing site-specific information
Ci     Passed to: nmpot
Ci
Ci   sarray:structure containing offsets to various arrays
Ci     Elts read: oics oipc onrc opos ormax oclabl
Ci     Passed to: secmat asaddq asaopm
Ci
Ci   sbz   :struct containing parameters for the BZ
Ci     Elts read: nkabc nkp ntet oidtet lmet n range w nevmx efmax lio
Ci                ndos dosw oqp owtkp zval
Ci     Passed to: optint
Ci
Ci   sham  :struct containing information about the hamiltonian
Ci     Elts read: nmto qss oeula neula ldham oindxo
Ci     Passed to: nmpot secmat asaddq
Ci
Ci   spot  :struct containing information about the potential
Ci     Elts read: opp osop vmtz opmpol oqpp osoptc
Ci     Stored:    oppn
Ci     Passed to: nmpot secmat asaddq asaopm
Ci
Ci   sstrn :struct containing global strings
Ci
Ci   nbmax :leading dimension of eband
Ci   vorb  :LDA+U potential
Co Outputs
Co   efermi:Fermi energy
Co   eband :energy bands; alias eb (sec*.f)
Co   nevmx :largest number of eigenvectors found
Co   qnu   :energy-weighted moments of the sphere charges
Co   sumev :sum of eigenvalues
Co   rhos  :spin density-matrix
Co   amag  :system magnetic moment
Co   aamom :local magnetic moments
Co   dmatu :site density-matrix corresponding to LDA+U blocks
Cl Local variables and switches
Cl   nspc: number of coupled spins (1 unless noncoll and/or S-O coupling).
Cl   nspx: number of independent spin channels containing eigenvalues
Cl         and total DOS; nspx=nsp unless nspc=2, in which case nspx=1
Cl   nev:  number of eigenvalues (and eigenvectors) calculated
Cl   nchan:number of l (or lm) channels summed over all classes,
Cl         e.g. nlo*nclass.
Cl   liomom:0 suppress writing dos weights to moments file
Cl         :1 write dos weights to moments file
Cl   lwtkb :0 weights are neither required a priori nor available
Cl         :1 weights are required a priori, and are available
Cl         :  (having read from disk or from a prior band pass)
Cl         :-1 weights are required a priori, but were not read
Cl   lnsph=1 make coffs for nonspherical density
Cl   ldos =a combination of the following integers
Cl         1 make dos
Cl         2 generate weights for partial dos
Cl         4 generate weights for m-decompos'n of pdos
Cl   lrho=T make the on-site spin-density matrix, needed for mag. forces.
Cl   ldens=T makes the density matrix (under construction)
Cl   lsx1   if 1, a response function will be computed
Cl          from eigenvectors created by bndasa; bndasa must generate
Cl          eigenvectors appropriate for this purpose.
Cl   nfilem is the logical units for moments file.
Cr Remarks
Cm  Memory:
Cm   Without exploiting hermiticity, and without folding down, the
Cm   minimum storage requirement is three complex matrices (H, O and Z)
Cm   dimension ndim^2. The last of these is not needed until the work-
Cm   space for the real-space structure constants has been released.
Cm   If ndim is separated into lower and intermediate sets
Cm   (ndim = ldim + idim) then if the i-waves remain in the basis the
Cm   requirement is 3ldim^2 + 6(ldim*idim) + 3idim^2. If the i-waves are
Cm   folded down the requirement becomes 3ldim^2 + 2(ldim*idim) + idim^2
Cm   but these must be declared simultaneously with the work-space for
Cm   the real-space structure constants.
Cf  Files:
Cf    For each qp, bndasa writes eband, and then accmom (and also accsm
Cf    if spin-pol bands coupled), provided metal=T and lwtkb=F.
Cf    If metal=F, writes doswt instead.
Cf    See iomoms for structure of the moments file.
Cb Bugs
Cb   lwtkb should be patterned after bndfp.
Cb   ldos should be decoupled from lwtkb (see makwts).  Requires caller to
Cb   allocate appropriate memory for accwt.
Cu Updates
Cu   26 Oct 08 (ATP) first attempt at inverse Bloch transform for band fitting
Cu   04 Jun 08 Orbital moments now generated when lmet=2 or lmet=3
Cu   08 Nov 07 (J. Xu) New LDA+U
Cu   08 Jun 06 optional weights for band plotting
Cu   13 Jul 05 (wrl) Nonlinear optics re-instated
Cu   27 Apr 04 Optics branch works with MPI
Cu   25 Apr 04 Reads efermi when lmet=2 and qp weights available.
Cu   21 Apr 04 Additions for an m-dependent spin-density matrix
Cu   20 Feb 04 (U.J.Rho) Optics code can to on-the-fly sampling integration
Cu   23 Sep 03 SX patterned after GW.  sigm(rs) obtained
Cu             by call to rdsigh.
Cu   06 Sep 03 (jek) extended optics code to metals
Cu   10 Mar 03 Poke double-counting from appl. field into eterms
Cu    1 Mar 03 Revised calculation of output magnetization
Cu             in noncollinear case; new double counting
Cu   27 Jan 03 In noncollinear mode, always accumulates density
Cu             matrix and prints output magnetization from it.
Cu    9 Jan 03 (K. Belashchenko) parallelize bndasa over k-points
Cu   30 May 02 Allow reading SX sigma in alpha repsn
Cu   18 Jan 02 Redesigned to handle accumulation of m-resolved weights
Cu             New command-line argument --pdos; see sumlst for syntax
Cu   03 Feb 01 patched code for m-resolved moments.
Cu   28 Apr 98 optical matrix elements and dielectric eps added
Cu    4 May 98 calls atfold so high not automatically folded to interm.
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nbmax,mpsord
      double precision sarray(1),sbz(1),sctrl(1),slat(1),ssite(1),
     .  sspec(1),sham(1),spot(1)
      character*(*) sstrn
      double precision qnu(1),eband(nbmax,*),efermi,vmtz(2),sumev,
     .  rhos(2,*),amag(3),aamom(*)
C ... For LDA+U; dimensions packed in sham
C     integer nlibu,lmaxu
C     double complex vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,2,nlibu)
C     double complex dmatu(-lmaxu:lmaxu,-lmaxu:lmaxu,2,nlibu)
      double complex vorb(*),dmatu(*)
C ... Local parameters
      character*120 strn,plbopt
      logical bittst,lgors,metal,ltet,onesp,cmdopt,lrhos,
     .  ldens,savvec,T,F,a2bin
      integer bitand,fopn,fopna,i,iscr,j,iobzwt,ifi,
     .  igets,ikp,kkp,iomoms,iprint,iprt,iskp,isp,isw,lbzio,
     .  ldham(16),ldim,ldimx,ldos,lgen3,lgunit,lham,lasa,lasas,lhams,
     .  lidim,lihdim,lmet,lncol,lnsph,loptic,lrout,lsx,lsx1,
     .  lrsig,lwtkb,lswtk,n1,n2,n3,nbas,nbmaxx,nclass,ndos,ndos0,nbf,
     .  neul,nrhos,nev,nevmx,nevsav,nevx,nfilem,nfstg,nschan,nkabc(3),
     .  nkp,nl,nlo,nlspc,nmto,nppn,nqpp,nsp,nspc,nspec,nspx,ntet,numq,
     .  plbnd,stdo,stdl,izval,moddos,nll,liomom,nfbn(2),nmom
      integer procid, master
C     Used for site-resolved partial dos (--pdos)
      integer nsitmx,lmxch,nchan
      parameter (nsitmx=1024)
      integer lsite(nsitmx)
      parameter (nppn=12)
      integer o,oaccsm,oaccwt,odos,oics,oidtet,oidxsh,oipc,onrc,oo,opos,
     .  ormax,oqp,owk,owtkb,oswtk,owtkp,ozll,osoptc,oeula,obdots,omagf,
     .  opp,osop,obxc,oclabl,oppn,oevl,oifbls
      double precision dum(20),avw,del,alat,plat(9),qspirl(4),fsmom,
     .  efloc(2),efmax,swidth,srnge,dosw(2),qp(3),dgets,ckbas,cksumf,
     .  xx(10),zval,ef0,swtkb,evtop,ecbot,eterms(20)
      parameter(T=.true., F=.false.)
      equivalence (n1,nkabc(1)),(n2,nkabc(2)),(n3,nkabc(3))
      equivalence (ldim,ldham(1)),(lidim,ldham(2)),(lihdim,ldham(3))
      equivalence (ldimx,ldham(5))
      equivalence (nspx,ldham(8))
      character*80 outs
C ... For optics
      integer nfilo,nfiup,nemlo,nemup,ocrng(2),unrng(2),nfilm,nempm,
     .  npopt,oeps,ooptmt,ovelmt,odnpp,oorbtm,owk1,owk3
      equivalence (nfilo,ocrng(1)),(nfiup,ocrng(2))
      equivalence (nemlo,unrng(1)),(nemup,unrng(2))
C ... For LDA+U
      integer nlibu,lmaxu,ludiag
      integer oidu
C ... For inverse Bloch transform
      integer mxorb,nglob
      integer ostrx,ohrs,oors,oh1rs,oh2rs,oo2rs,os,odon,
     .        ostrxs,ohrss,oorss

#if STONER
      integer l,nchan,odoswt,oeband
#endif
C ... For the density matrix
      integer nsite,oalph,ontab,oiax,odmat,odmats
C ... For nonspherical output density
      integer ocg,ojcg,oidxcg,opmpol,oqpp
C ... For symmetry operations
      integer nsgrp,oistab,oag,osymgr
C     integer nstar,oqstar,ormat,oistab,igstar,ig
C     double precision qlat(3,3)
C ... Heap
      integer w(1)
      common /w/ w
C     bittst(n,bit) = (mod(n,bit+bit) - mod(n,bit) .eq. bit)
C ... External calls
      external amagnc,asaddq,asaopm,awrit0,awrit1,awrit2,awrit3,awrit7,
     .         bdotsr,bzints,bzwtsf,dcopy,defdc,defdr,defrr,dosio,
     .         dpdump,dpscop,dpzero,dscal,dvset,fclose,fclr,fexit2,info,
     .         info0,iomomn,iomomq,iorbtm,iostr,isanrg,lsets,makdos,
     .         maknos,moment,nmpot,optint,pack1,pack2,poppr,pptrns,
     .         pshpr,query,rlse,rsmsym,rx,rx0,rx1,secmat,setpr,shostr,
     .         subzi,sumlst,suqlst,suqlsw,togpr,upack,upack1,upack2,
     .         xxxdif
#if MPI | MPIK
      include "mpif.h"
      integer numprocs, ierr, status(MPI_STATUS_SIZE)
      integer MAX_PROCS
      parameter (MAX_PROCS = 100)
      integer resultlen
      character*(MPI_MAX_PROCESSOR_NAME) name
      character*10 shortname(0:MAX_PROCS-1)
      character*20 ext
      character*26 datim
      integer namelen(0:MAX_PROCS-1)
      integer ista, iend
      double precision buf
C      integer length(0:MAX_PROCS-1)
C      integer offset(0:MAX_PROCS-1)
      double precision sttime,entime
      logical mlog
#if MPIK
      integer , dimension(:),allocatable :: kpproc,offset,length
      integer obuf,oqptr
#endif

      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
      call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
      call strcop(shortname(procid),name,10,'.',i)
      namelen(procid) = i-1
      master = 0
#else
      procid = 0
      master = 0
#endif

C --- Unpack parameters and parameter setup ---
      call tcn('bndasa')
      stdo = lgunit(1)
C     stde = lgunit(1)
      stdl = lgunit(2)
      liomom = 1
#if MPI | MPIK
      liomom = 0
      mlog = cmdopt('--mlog',6,0,strn)
#endif
      call upack('bz nkabc nkp ntet oidtet lmet',sbz,nkabc,nkp,ntet,
     .  oidtet,lmet)
      call upack('bz n range w nevmx efmax',sbz,mpsord,srnge,
     .  swidth,nevmx,efmax)
      call upack('bz lio ndos dosw oqp owtkp',sbz,lbzio,ndos0,dosw,oqp,
     .  owtkp)
      zval = dgets('bz zval',sbz)
      call upack2('bz zval fsmom',sbz,zval,fsmom)
      call upack('lat avw alat plat',slat,avw,alat,plat,0,0)
      call upack('lat ocg ojcg oidxcg',slat,ocg,ojcg,oidxcg,0,0)
      call upack('lat nsgrp oistab oag osymgr',slat,nsgrp,oistab,oag,
     .  osymgr,0)
      call upack('ctrl nbas nclass nl nspec nspin',sctrl,nbas,nclass,nl,
     .  nspec,nsp)
      call upack('ham qss oeula neula obdots eterms',sham,qspirl,oeula,
     .  neul,obdots,eterms)
      call upack('ham omagf nbf nmto',sham,omagf,nbf,nmto,0,0)
      call upack('ham nlibu lmaxu udiag',sham,nlibu,lmaxu,ludiag,0,0)
      call dvset(eterms,17,18,0d0)
      call upack('pot opp osop obxc nrhos',spot,opp,osop,obxc,nrhos,0)
      call upack('ctrl lncol lsx',sctrl,lncol,lsx,0,0,0)
      lrsig= igets('ham lsig',sham)
      iscr = mod(igets('ctrl lscr',sctrl),10)
      lsx1 = mod(lsx,2)
      if (mod(iscr,2) .eq. 1) lsx1 = 1
C     Hang onto lham, lasa to restore after exiting bndasa
      call upack('ctrl lham lasa lgen3 loptc',sctrl,lhams,lasas,lgen3,
     .  loptic,0)
      lnsph = isw(lgors('ctrl lasa,32',sctrl))
      call upack('array oics oipc onrc opos ormax',sarray,oics,oipc,
     .  onrc,opos,ormax)
      call upack('pot vmtz opmpol oqpp',spot,vmtz,opmpol,oqpp,0,0)
      call upack('ham ldham oindxo',sham,ldham,oidxsh,0,0,0)
      call defrr(o,1)
      if (lnsph .ne. 0 .or. nlibu .gt. 0) then
        i = nl**2
        nqpp = (i*(i+1))/2
        call zinit(w(oqpp),nqpp*4*nsp*nbas)
      else
        nqpp = 0
        oqpp = 1
      endif
      if (lmet .gt. 3) call rx1('BNDASA: METAL=%i not implemented',lmet)
      ef0 = 0
C     Flag indicating Fermi level is not known
      efermi = -99

C ... Get nldau for ldau
      if (nlibu .ne. 0) then
        call defi(oidu,4*nbas)
        call getidu(nbas,sspec,ssite,w(oidu))
      else
        call defi(oidu,1)
      endif

C ... Switch to plot bands at specified qp
      if (cmdopt('--band',6,0,strn)) then
        plbnd = 1
        plbopt = strn(7:)
        nevmx = -1
        loptic = 0
        numq = 1
        call defi(oifbls, -ldimx*2*2)
        lmet = 0
        lsx1 = 0
      else
        plbnd = 0
      endif

      nspc = 1
      if (bitand(lncol,1+2+4+8) .ne. 0) nspc = 2
C     lrhos  = bittst(lncol,16) .and. nspc .eq. 2 .and. nevmx .ge. 0
      lrhos  = nspc .eq. 2 .and. nevmx .ge. 0
      metal  = lmet .ne. 0
      ltet  = ntet .gt. 0
      ldos   = igets('ctrl ldos,7',sctrl)
      nsite = 0
      if (cmdopt('--pdos',6,0,strn)) then
        nll = nl
        call sumlst(0,0,nbas,nsgrp,ssite,sspec,strn(7:),moddos,
     .    nsite,lsite,lmxch,nchan,w,nll)
        if (moddos.eq.2 .or. moddos.eq.5 .and. ldos.lt.4) then
          call lsets('ctrl ldos',sctrl,.true.,4)
          ldos = igets('ctrl ldos,7',sctrl)
        endif
      endif
      nlo = nl
      if (ldos .ge. 4) nlo = nl*nl
      if (liomom .eq. 0 .and. lmet .eq. 1) call rx('BNDASA: METAL=1 '
     .  //'incompatible with suppressing I/O to moms file')
      if (nsite .ne. 0 .and. lmet .eq. 1 .and. iprint().ge.10) then
        call info(10,1,0,' ***WARNING*** bndasa: partial DOS mode not'//
     .    ' compatible with metal=1:',0,0)
        call info(10,0,0,' Program will not accumulate output moments.',
     .    0,0)
      endif
      if (nsgrp.ne.1 .and. nlo.ne.nl .and. iprint().ge.10) then
        call info(10,1,0,' ***WARNING*** bndasa: ASA does not '//
     .    'symmetrize m-resolved moments.',0,0)
        call info(10,0,0,' Suppress symops for m-resolved DOS.',0,0)
      endif

      ldens  = bittst(lbzio,8)
      onesp  = cmdopt('--onesp',7,0,outs)
      ndos   = iabs(ndos0)
      nlspc = nl * nsp * nclass
      oaccwt = 1
      oaccsm = 1
      nevsav = -1
      lrout  = 1
      if (nevmx .lt. 0) lrout = 0
      nbmaxx = nbmax
      if (nspc .eq. 2) nbmaxx = nbmax*nsp

C --- Setup pertaining to screened exchange ---
      if (lrsig+lsx1 .ne. 0) then
C   ... Turn on Gamma rep, ortho z; KKR rho; no downfolding
        if (lsx1 .eq. 1) then
          call lsets('ctrl lham',sctrl,T,128)
          call lsets('ctrl lham',sctrl,T,32)
          call lsets('ctrl lham',sctrl,F,4)
          call lsets('ctrl lasa',sctrl,F,256)
          nevmx = nbmax
          efmax = 999
        endif
C   ... SX -> Gamma rep
        if (lrsig .ne. 0) then
          call lsets('ctrl lham',sctrl,T,128)
        endif
      endif

C     Unpack (possibly altered) lham, lasa
      call upack2('ctrl lham lasa',sctrl,lham,lasa)
      call pack2('ham lham',sham,lham,0)

C --- Printout at entry ---
      if (iprint() .ge. 30) then
        call awrit8('%x%N --- BNDASA : band pass ('//
     .    '%?#n#sigma, ##'//
     .    '%?#n#two-c-H, ##'//
     .    '%?#n#ortho-z, ##'//
     .    '%?#n#H-gamma, ##'//
     .    '%?#n#KKR-qout, ##'//
     .    '%?#n#NMTO=%-1j%i, ##'//
     .    '%?#n#noncol, ##'//
     .    '%b%b%?#p>24#)## ---',outs,80,-stdo,
     .    lrsig,bitand(lham,1),bitand(lham,32),bitand(lham,128),
     .    256-bitand(lasa,256),nmto,lncol,0)
      endif

C ... Starting Euler angles
      if (bittst(lncol,1)) then
        if (iprt(2) .ne. iprt(1)) then
          call pshpr(iprt(2))
        else
          call pshpr(iprint()-0)
        endif
C       call setpr(50)
        call amagnc(nbas,nl,w(oipc),rhos,nrhos,qnu,w(oeula),neul,1,amag,
     .    aamom)
        if (bittst(lncol,8)) then
          call bdotsr(nbas,nl,w(oipc),rhos,nrhos,qnu,
     .      w(obdots),lihdim,w(oidxsh),1,dum)
C         ehterm(5) = dum(1)
          eterms(17) = dum(1)
        endif
        if (neul .gt. 1 .or. bittst(lncol,8)) then
          if (.not. bittst(lncol,8)) nbf = 99
          call bsrhos(nbas,nl,w(oipc),rhos,nrhos,qnu,w(opp),
     .      w(osop),w(oeula),neul,w(obxc),w(omagf),nbf,lihdim,w(oidxsh),
     .      1,dum)
          eterms(17) = dum(2)
        endif

        call poppr
      endif

C ... Sanity checks
      call isanrg(nevmx,loptic-1,9999,'optics sought and','nevmx',T)

C ... Allocate space for density matrix.  For now, same range as s
      if (ldens) then
        ckbas = cksumf(w(opos),3*nbas)
        call iostr(8,'STR',nl,nbas,1,dum,0,ckbas,-1,nsite,oalph,
     .    oiax,ontab,odmat)
        call rlse(odmat)
        call defdc(odmat,-w(ontab+nbas)*nl**4)
      endif
C      if (.true.) then
C        call gtstr0(nkp,igstar,w(owtkp),nstar)
C        call defrr(oqstar,3*nstar)
C        call defrr(ormat,nl**4*nstar)
C        call defi(oistab,nbas*nstar)
C        qp(1) = -3.75d-01
C        qp(2) =  1.25d-01
C        qp(3) = -3.75d-01
C      endif

C --- Memory allocation and defaults for optics calculations ---
      if (loptic .ne. 0) then
        call upack1('pot osoptc',spot,osoptc)
        call upack3('optic ocrng unrng ne',w(osoptc),ocrng,unrng,npopt)
        if (metal) then
          if (nfilo .eq. 0) nfilo = 1
          if (nfiup .eq. 0) nfiup = ldimx
          if (nemlo .eq. 0) nemlo = 1
          if (nemup .eq. 0) nemup = ldimx
        else
          if (nfilo .eq. 0) nfilo = 1
          if (nfiup .eq. 0) nfiup = (nspc*nint(zval)+1)/2
          if (nemlo .eq. 0) nemlo = nfiup+1
          if (nemup .eq. 0) nemup = ldimx
        endif
        if (nfiup .gt. nemup) then
          call info2(30,0,0,' (warning, optics) last '//
     .      'filled state (=%i) > last empty state (=%i)',nfiup,nemup)
          nfiup = nemup
        endif
        if (nemlo .lt. nfilo) then
          call info2(30,0,0,' (warning, optics) first '//
     .      'empty state (=%i) < first filled state (=%i)',nemlo,nfilo)
          nemlo = nfilo
        endif
        if (npopt .eq. 0) npopt = 400
        call pack3('optic ocrng unrng ne',w(osoptc),ocrng,unrng,npopt)
        nfilm = nfiup-nfilo+1
        nempm = nemup-nemlo+1
        ooptmt = 1
        ovelmt = 1
        if (loptic .gt. 1) then
          call defrr(ooptmt,-3*nfilm*nempm*nspx)
          call defrr(ovelmt,-3*nfilm*nspx)
          call defdr(oeps,-(1+nspx*3)*npopt)
        elseif (loptic .gt. 0) then
          call defrr(ooptmt,-3*nfilm*nempm*nkp*nspx)
          call defrr(ovelmt,-3*nfilm*nkp*nspx)
        endif
C        nevmx = nbmaxx
C        efmax = 999
      endif

C --- Setup for BZ integration ---
      if (plbnd .eq. 0) then
C       call pshpr(0)
        call subzi(lmet,ntet.gt.0,lrout.gt.0,ldim,nsp,nspc,nkp,zval,
     .    nevmx,lwtkb,efermi,oevl,owtkb,numq,lswtk,oswtk,ef0,0d0,mpsord,
     .    swidth)
C       call poppr
C       bndasa takes eband as a passed argument
        call rlse(oevl)
C       Chuck wtkb if it won't be needed until after bands
        if (lmet .eq. 1) call rlse(owtkb)
      else
        nkp = 0
        ldos = 0
        loptic = 0
        lwtkb = 0
        nevmx = -1
        nfstg = 1
      endif

C --- NMTO setup ---
      if (lgen3 .ne. 0) then
        if (nevmx .ge. 0 .and. iprint() .ge. 30) then
          call awrit0(' (warning) no output density yet in 3rd gen',
     .      outs,80,stdo)
        endif
        nevmx = -1
        nfstg = 1

C   ... NMTO Potential parameters
        call upack1('array oclabl',sarray,oclabl)
        call defdr(oppn,nppn*lihdim*nmto*nsp)
        call togpr
        call nmpot(40,sctrl,sspec,ssite,sham,spot,lidim,
     .    lihdim,w(oidxsh),w(oclabl),w(oppn))
        call togpr
        call pack1('pot oppn',spot,oppn)
      endif

C --- More k-point independent local arrays and initialization ---
C     if (lrhos) call defdc(oaccsm,nlo*nclass*ldimx*nspc*3)
      if (lrhos) call dpzero(rhos,2*3*nrhos*2*2*nclass)
      call dvset(eband,1,nbmax*nsp*nkp,9999d0)
      if (loptic .ne. 0 .and. nevmx .lt. nemup) then
        if (iprint() .ge. 10) call awrit2(' bndasa (warning): '//
     .    'nevmx reset from %i to %i',' ',80,stdo,nevmx,nemup)
        nevmx = nemup
      endif
      if (lrout .ne. 0) call dpzero(qnu,3*nlspc)
      oorbtm = 1
      if (nspc .eq. 2) call defrr(oorbtm,-nlo*2*nclass)
      savvec = lsx1.eq.1 .or. bittst(lham,64) .and. nevmx.ge.0

C --- Start loop over k points; also, re-entry for second band pass ---
      nevsav = nevmx
   99 continue
      call defdc(ozll,ldimx**2)
      nevmx = nevsav
      if (lwtkb .eq. -1) nevmx = -1
C     Information to write to moments file
      nmom = nspc
      if (lwtkb .eq. 0 .and. metal) nmom = nspc+2
      call iomomn(metal,ldos,lrhos,nevmx,nspc,nmom,lwtkb,nfstg)

      if (lwtkb .eq. -1)
     .  call info(20,0,0,' Start first of two band passes ...',0,0)

C      if (lwtkb .eq. 1) then
C        call setpr(110)
C      endif

C ... Case bands at usual qp list
      if (plbnd .eq. 0) then

C       Moments file setup
        if (liomom .ne. 0) then
          nfilem = fopna('MOMS',-1,4)
          rewind nfilem
          write (nfilem) nl, nsp, nspc, nkp, ldim, nfstg
        else
          nfilem = 0
        endif
C       Eigenvector file setup
        if (savvec) then
          i = fopna('EVEC',-1,4)
          rewind (i)
          write(i) nkp,nsp
        endif

C ... Case only bands at supplied qp generated: setup
      else

C   ... Try and read Fermi level from file
        if (procid .eq. master) then
        ifi = fopna('wkp',-1,4)
        i = iobzwt(1,i,i,i,ef0,w,ifi)
        call fclr('wkp',ifi)
        endif
        call mpibc1(ef0,1,4,.false.,'bndasa','ef0')
C   ... Override Fermi level with command-line value
        i = 4
        if (cmdopt('-ef=',i,0,strn)) then
          if (.not. a2bin(strn,ef0,4,0,' ',i,-1)) call
     .      rxs2('BNDASA: failed to parse "',strn(1:30),'%a"')
        endif
        i = nspx
        call suqlst(plbopt,0,ldimx,ef0,i,w,nfbn,w(oifbls),nkp,qp,onesp)
        if (nkp .le. 0) call rx0('bndasa: finished bands')

      endif

      call pshpr(iprint())
      nevx = 0
      izval = zval/2
      evtop = -9999
      ecbot = -evtop

C ... Loop over k-points

#if MPE & MPIK
      ierr = MPE_LOG_EVENT(EVENT_START_KLOOP,procid,"k-loop")
#endif
#if MPIK
      call info0(30,1,0, ' ... Start MPI k-loop')
      allocate (kpproc(0:numprocs), stat=ierr)
      call dstrbp(nkp,numprocs,1,kpproc(0))
      sttime = MPI_WTIME()
      do  10  ikp = kpproc(procid), kpproc(procid+1)-1
C       print *, 'start ikp',procid,ikp
        if (ikp .eq. kpproc(procid)) then
          if (mlog) then
            call gettime(datim)
            call awrit4(' bndasa '//datim//' Process %i of %i on '
     .        //shortname(procid)(1:namelen(procid))//
     .        ' starting k-points %i to %i',' ',256,lgunit(3),
     .        procid,numprocs,kpproc(procid),kpproc(procid+1)-1)
          endif
        endif
#else
      ostrx  = 1
      ohrs   = 1
      oors   = 1
      oh1rs  = 1
      oh2rs  = 1
      oo2rs  = 1
      ostrxs = 1
      ohrss  = 1
      oorss  = 1
      if (cmdopt('--invbl',7,0,outs)) then
        ckbas = cksumf(w(opos),3*nbas)
        call iostr(8,'STR',nl,nbas,1,dum,0,ckbas,-1,nsite,oalph,
     .             oiax,ontab,os)
        nsite = w(ontab+nbas)
        call rlse(os)
        call defdr(ostrx,-nl**4*nsite)
        call defdr(ohrs,-nl**4*nsite)
        call defdr(oors,-nl**4*nsite)
        call defdr(oh1rs,-nl**4*nsite)
        call defdr(oh2rs,-nl**4*nsite)
        call defdr(oo2rs,-nl**4*nsite)

        call defdr(ostrxs,-nl**4*nsite)
        call defdr(ohrss,-nl**4*nsite)
        call defdr(oorss,-nl**4*nsite)
      endif
      do  10  ikp = 1, nkp
#endif

      iskp = nsp*(ikp-1)
      call setpr(iprt(1))
#if MPIK
      if (ikp .gt. 1 .and. mod(ikp,100) .ne. 1) call setpr(iprt(1)-26)
#else
      if (ikp .gt. 1 .and. mod(ikp,100) .ne. 1) call setpr(iprt(1)-6)
#endif

C ... Get qp either from qp list or read from suqlst
      if (plbnd .ne. 0) then
        i = nspx
        call suqlst(plbopt,0,ldimx,ef0,i,w,nfbn,w(oifbls),nkp,qp,onesp)
        iskp = 1
        if (nfbn(1) .gt. 0) then
          nevmx = ldimx
          efmax = 9999
        endif
      else
        call dpscop(w(oqp),qp,3,3*ikp-2,1,1d0)
      endif

C ... debugging to retrieve info about star of q
C      call gtstr(111,ikp,qp,w(osymgr),w(oag),igstar,w(owtkp),nl,plat,
C     .  nbas,w(opos),nstar,w(oqstar),w(ormat),w(oistab))
C      goto 10

      do  11  isp = 1, nsp
C   ... No second spin for coupled-spin case
        if (nspc .eq. 2 .and. isp .eq. 2) goto 11
C   ... Skip over second spin in collinear antiferromagnetic case
        if (onesp .and. isp .eq. 2) then
          call dcopy(ldim,eband(1,1+iskp),1,eband(1,2+iskp),1)
          goto 11
        endif

C   --- Eigenvalues and vectors of ASA hamiltonian ---
        call secmat(sctrl,sspec,sham,spot,sarray,slat,ikp,nkp,qp,
     .    w(owtkp),isp,lrsig,nevmx,efmax,nlibu,lmaxu,w(oidu),ludiag,
     .    vorb,nev,w(ozll),eband(1,isp+iskp),w(ostrx),w(ohrs),w(oors))
        evtop = max(evtop,eband(izval,isp+iskp))
        ecbot = min(ecbot,eband(izval+1,isp+iskp))

        if (cmdopt('--invbl',7,0,outs)) goto 10
        if (plbnd .eq. 0) then

C          if (check .ne. 0) then
C          call fixef0(zval,isp,nspc,nevx,nbmax,eband,dosw,ef0)
C          if (ef0+.1d0.gt.eband(max(nev,1),isp) .and. nevmx.ge.0 .and.
C     .        lmet.ne.0) then
C            call awrit3('%N evl(nev=%i)=%;3d but '//
C     .        'ef0=%;3d ... restart with larger efmax or nevmx',
C     .        ' ',80,stdo,nev,eband(max(nev,1),isp),ef0)
C            call rx('bndfp')
C          endif
C          endif
C
C          if (ef0+.1d0.gt.evl(max(nev,1),isp) .and. nevmx.ge.0 .and.
C     .        lmet .ne. 0) then
C            call awrit3('%N evl(nev=%i)=%;3d but '//
C     .        'ef0=%;3d ... restart with larger efmax or nevmx',
C     .        ' ',80,stdo,nev,evl(max(nev,1),isp),ef0)
C            call rx('bndfp')
C          endif

        else
          call suqlsw(ldimx,isp,nspx,eband(1,isp+iskp))
          if (nfbn(1) .ne. 0) then
            if (nmto .eq. 0)
     .        call ztoyy(w(ozll),ldimx,ldimx,ldimx,ldimx,0,1)
            call suqlse(ldimx,isp,nspx,ldimx,1,nfbn,w(oifbls),ldimx,
     .        w(ozll),eband(1,isp+iskp))
            if (nfbn(2) .ne. 0)
     .      call suqlse(ldimx,isp,nspx,ldimx,2,nfbn,w(oifbls),ldimx,
     .        w(ozll),eband(1,isp+iskp))
          endif
          goto 11
        endif

C       Sanity checks
        if (loptic .ne. 0 .and. lwtkb .ge. 0 .and. nev .lt. nemup) then
          call rx('bndasa: too few eigenvectors generated'//
     .      '...increase NEVMX or EFMAX in BZ')
        endif

C   ... Write energy bands to moms file and jump to end of k-loop if no evecs
        nevx = max(nevx,nev)
        if (nevmx .le. 0 .and. liomom .ne. 0) then
          write (nfilem) 0, ldimx
          call dpdump(eband(1,isp+iskp),ldimx,-nfilem)
          goto 11

        endif

C   --- Contribution to output density and related quantities from qp --
C   ... Save eigenvectors
        if (savvec) then
          i = fopna('EVEC',-1,4)
C         Record: eb header.  Read 1st eb with mc -r:br:s=1
          write(i) ldimx,1,1
C         Record: eband
          call dpdump(eband(1,isp+iskp),ldimx,-i)
C         Record: evec header.  Read 1st evec with mc -r:br:s=3
          write(i) ldimx,nev,2
C         Record: evecs
          call dpdump(w(ozll),ldimx**2*2,-i)
        endif

        if (lrout .ne. 0) then
        odnpp = 1
        if (loptic .gt. 0) call defrr(odnpp,nl**2*nbas*nspc*4*nev)
        call asaddq(sctrl,slat,sham,spot,sarray,ikp,nkp,isp,nev,nbmax,
     .    eband(1,isp+iskp),w(ozll),
     .    nlibu,lmaxu,w(oidu),vorb,
     .    nfilem,lwtkb,w(owtkb),w(owtkp),zval,
     .    metal,loptic.gt.0,nqpp,ldos,lrhos,nsite,lsite,0,qnu,w(oqpp),
     .    w(odnpp),rhos,w(oorbtm))
        endif

C   ... Optical matrix elements
        if (loptic .gt. 0 .and. lwtkb .ne. -1) then
          if (nemup .gt. nev) call fexit2(-1,111,' Exit -1 BNDASA: '
     .      //'optics needs %i evecs but calculated only %i',nemup,nev)
          kkp = ikp
          if (mod(loptic,10) .gt. 1) then
            kkp = 1
            if (isp .eq. 1) then
              call dpzero(w(ooptmt),3*nfilm*nempm*nspx)
              call dpzero(w(ovelmt),3*nfilm*nspx)
            endif
          endif
          call asaopm(sctrl,sarray,spot,slat,w(osoptc),nl,isp,nsp,
     .      nspc,nclass,nbas,ldimx,nev,qp,kkp,nkp,eband(1,isp+iskp),
     .      nbmaxx,w(odnpp),nfilm,nempm,w(ooptmt),w(ovelmt))
          call rlse(odnpp)

C     ... On-the-fly sampling integration of optics
          if (mod(loptic, 10).gt. 1 .and. nevmx .gt. 0) then
            call defdr(owk1,           npopt)
            call defdr(owk3,         3*npopt)
            call optint(sctrl,w(osoptc),sbz,eband,nbmaxx,isp,nspx,
     .        nspc,ikp,efermi,w(oidtet),w(owtkp),dgets('lat vol',slat),
     .        nfilm,nempm,npopt,w(owk1),w(owk3),w(ooptmt),w(oeps))
            call rlse(owk1)
          endif
        endif

   11 continue
C ... End of k-point loop
   10 continue
      if (cmdopt('--invbl',7,0,outs)) then
        call upack1('array oclabl',sarray,oclabl)
        ckbas = cksumf(w(opos),3*nbas)
        call iostr(8,'STR',nl,nbas,1,dum,0,ckbas,-1,nsite,oalph,
     .             oiax,ontab,os)
        nsite = w(ontab+nbas)
        mxorb = nglob('mxorb')
C ... symmetrise real space sqrdel S sqrdel
        call rsmsym(0,plat,mxorb,w(oidxsh),ldim,nbas,w(opos),nl,nsp,
     .              1,nsite,w(ontab),w(oiax),w(osymgr),w(oistab),
     .              nsgrp,nl*nl,w(ostrx),w(ostrxs))
C ... symmetrise real space H
        call rsmsym(0,plat,mxorb,w(oidxsh),ldim,nbas,w(opos),nl,nsp,
     .              1,nsite,w(ontab),w(oiax),w(osymgr),w(oistab),
     .              nsgrp,nl*nl,w(ohrs),w(ohrss))
C ... symmetrise real space O
        call rsmsym(0,plat,mxorb,w(oidxsh),ldim,nbas,w(opos),nl,nsp,
     .              1,nsite,w(ontab),w(oiax),w(osymgr),w(oistab),
     .              nsgrp,nl*nl,w(oors),w(oorss))
        call defdr(odon, nbas*nl**4)
        call wrirsh(T,T,bittst(lham,128),nl,nsp,nclass,nbas,plat,
     .              w(opos),alat,w(oclabl),w(oipc),w(opp),w(oalph),
     .              nsite,w(oiax),w(ostrxs),w(ohrss),w(oorss),w(oh1rs),
     .              w(oh2rs),w(oo2rs),w(odon))
        call iostr(1,'HAM',nl,nbas,1,dum,0,ckbas,-1,nsite,oalph,
     .             oiax,ontab,ohrss)
        call iostr(1,'HAM1',nl,nbas,1,dum,0,ckbas,-1,nsite,oalph,
     .             oiax,ontab,oh1rs)
        call iostr(1,'HAM2',nl,nbas,1,dum,0,ckbas,-1,nsite,oalph,
     .             oiax,ontab,oh2rs)
        call iostr(1,'OVL',nl,nbas,1,dum,0,ckbas,-1,nsite,oalph,
     .             oiax,ontab,oorss)
        call iostr(1,'OVL2',nl,nbas,1,dum,0,ckbas,-1,nsite,oalph,
     .             oiax,ontab,oo2rs)
        call fexit(0,9,'Real space hamiltonian written to disc',0)
      endif
      call poppr
      call rlse(ozll)

C --- Extract dmatu ---
      if (nlibu .gt. 0) then
        call asadmu(nbas,nsp,lmaxu,nl,w(oidu),nqpp,w(oqpp),ldim,
     .    w(oidxsh),dmatu)
      endif

C ... Gathering the results from all processes
#if MPIK
      entime = MPI_WTIME()
      call info2(30,0,0, ' ... Done MPI k-loop: %;1d seconds elapsed',
     .  (entime-sttime),0)

      allocate (offset(0:numprocs), stat=ierr)
      allocate (length(0:numprocs), stat=ierr)
      offset(0) = 0

      if (procid .eq. master) then
        call info0(20,0,-1,' Sharing data between processes...')
        sttime = MPI_WTIME()
      endif
      call MPI_ALLREDUCE(evtop,buf,1,MPI_DOUBLE_PRECISION,MPI_MAX,
     ,  MPI_COMM_WORLD,ierr)
      evtop = buf
      call MPI_ALLREDUCE(ecbot,buf,1,MPI_DOUBLE_PRECISION,MPI_MIN,
     ,  MPI_COMM_WORLD,ierr)
      ecbot = buf

C     Collect eband from various processors
      do  i = 0, numprocs-1
        ista = kpproc(i)
        iend = kpproc(i+1)-1
        length(i) = (iend - ista + 1)*nsp*nbmax
        offset(i+1) = offset(i) + length(i)
C       print *, 'ista=',i,ista,iend,offset(i),length(i)
      enddo
C     print *, offset(0:numprocs-1)
C     print *, length(0:numprocs-1)
      ista = kpproc(procid)
      call defrr(obuf,nbmax*nkp*nsp)
C     print 974, procid, eband(1,1:nkp)
      call MPI_ALLGATHERV(eband(1,1+nsp*(ista-1)),
     ,  length(procid),MPI_DOUBLE_PRECISION,w(obuf),length,
     ,  offset,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
      call dcopy(nbmax*nsp*nkp,w(obuf),1,eband,1)
      call rlse(obuf)
C     print 974, procid, eband(1,1:nkp)
C 974 format(i4/(6f12.6))

C     Add the qnu contributed by various processors
      if (lwtkb .ne. -1 .and. lrout .gt. 0) then
        call defrr(obuf,3*nlspc)
        call MPI_ALLREDUCE(qnu,w(obuf),3*nlspc,MPI_DOUBLE_PRECISION,
     ,    MPI_SUM,MPI_COMM_WORLD,ierr)
        call dcopy(3*nlspc,w(obuf),1,qnu,1)
        call rlse(obuf)
      endif

C     Add the rhos contributed by various processors
      if (lwtkb .ne. -1 .and. lrout .gt. 0 .and. lrhos) then
        i = 2*3*nrhos*2*2*nclass
        call defrr(obuf,i)
        call MPI_ALLREDUCE(rhos,w(obuf),i,MPI_DOUBLE_PRECISION,
     ,    MPI_SUM,MPI_COMM_WORLD,ierr)
        call dcopy(i,w(obuf),1,rhos,1)
        call rlse(obuf)
      endif

C     Add the optical matrix elements contributed by various processors
      if (lwtkb.ne.-1 .and. lrout.gt.0 .and. mod(loptic,10).eq.1) then
        i = 3*nfilm*nempm*nkp*nspx
        if (loptic .gt. 1) i = 3*nfilm*nempm*nspx
        call defrr(obuf,i)
        call MPI_ALLREDUCE(w(ooptmt),w(obuf),i,MPI_DOUBLE_PRECISION,
     ,    MPI_SUM,MPI_COMM_WORLD,ierr)
        call dcopy(i,w(obuf),1,w(ooptmt),1)
        call rlse(obuf)
      endif

C     Add the response function contributed by various processors
      if (lwtkb.ne.-1 .and. lrout.gt.0 .and. mod(loptic,10).eq.2) then
        i = (1+nspx*3)*npopt
        call defrr(obuf,i)
        call MPI_ALLREDUCE(w(oeps),w(obuf),i,MPI_DOUBLE_PRECISION,
     ,    MPI_SUM,MPI_COMM_WORLD,ierr)
        call dcopy(i,w(obuf),1,w(oeps),1)
        call rlse(obuf)
      endif

      deallocate(kpproc, stat=ierr)
      deallocate(offset, stat=ierr)
      deallocate(length, stat=ierr)

      if (procid .eq. master) then
        entime = MPI_WTIME()
        call info2(20,0,0,' MPI broadcast took %;1d sec',
     .    (entime-sttime),0)
      endif

C     Printout of eigenvalues
      call info0(30,0,0,' Printout of eigenvalues ...')
      call pshpr(iprint())
      do  ikp = 1, nkp
        call dpscop(w(oqp),qp,3,3*ikp-2,1,1d0)
        if (ikp .gt. 1 .and. mod(ikp,100) .ne. 1) call setpr(iprt(1)-6)
        iskp = nsp*(ikp-1)
C       call setpr(iprt(1))
        do  isp = 1, nsp
C     ... No second spin for coupled-spin case
          if (nspc .ne. 2 .or. isp .ne. 2) then
            j = min(9,ldimx)
            if (iprint() .ge. 35) j = ldimx
            if (iprint() .ge. 30) then
              call awrit3(' kpt %i of %i, k=%3:2,5;5d',
     .          ' ',80,stdo,ikp,nkp,qp)
              write(stdo,'(255(9f8.4:/))') (eband(i,isp+iskp), i=1,j)
            endif
            if (procid .eq. master) then
              call awrit3(' kpt %i of %i, k=%3:2,5;5d',
     .          ' ',80,stdl,ikp,nkp,qp)
              write(stdl,'(255(9f8.4:/))') (eband(i,isp+iskp), i=1,j)
            endif
          endif
        enddo
      enddo
      call poppr

#endif

C ... Finish on-the-fly sampling integration of optics 
      if (loptic.eq.2 .and. lwtkb.ge.0 .and. procid.eq.master) then
        call optin2(sctrl,w(osoptc),sbz,nspx,nspc,
     .    dgets('lat vol',slat),npopt,w(oeps))
      endif

C ... Case generating bands: find next block of qp
      if (plbnd .ne. 0) goto 99

C --- Interpolate density to Fermi energy (lmet=4) ---
C     sev = sumev(1,1)
      if (lmet .eq. 4) then
        call rx('bndasa not ready for lmet=4')
      endif

C --- BZ integration for fermi level, band sum and qp weights ---
      if (lmet .ne. 4 .or. ntet .gt. 0) then
        if (lmet .le. 1) call defdr(owtkb,ldim*nsp*nkp)
        call bzwtsf(nbmax,ldim,nsp,nspc,n1,n2,n3,nkp,ntet,w(oidtet),
     .    zval,fsmom,metal,ltet,mpsord,ndos,swidth,srnge,w(owtkp),
     .    eband,efmax,-2,w,efloc,sumev,w(owtkb),xx,lwtkb)
        efermi = efloc(1)

        if (lmet .ne. 4) then
          ef0 = efermi
          call pack1('bz ef',sbz,ef0)
        endif
        if (lwtkb .eq. -1 .and. lrout .gt. 0) then
          call info(20,0,0,' Start second band pass ...',0,0)
          lwtkb = 1
          goto 99
        endif
      endif

C ... Save k-point weights and Fermi level in wkp file
      if (lmet .eq. 0) then
C        if (iprint() .ge. 20)
C     .  call awrit3(' Highest occ. level = %,5;5d '//
C     .  ' Lowest unocc. = %,5;5d  diff =  %,5;5d',' ',80,stdo,
C     .  evtop,ecbot,ecbot-evtop)
      endif
      if (procid .eq. master) then
        ifi = fopna('wkp',-1,4)
        j = 0
        if (lmet .eq. 0) j = 1
        i = iobzwt(j,ldimx,nkp,nspx,efermi,w(owtkb),-ifi)
        call fclr('wkp',ifi)
      endif

C --- Generate density matrix ---
      if (ldens) then
        call defdc(odmats,-nsite*nl**4)
        call rx('bndasa: fix call to rsmsym')
C        call rsmsym(1,plat,nbas,w(opos),nl,nsp,1,nsite,
C     .    w(ontab),w(oiax),w(osymgr),w(oag),nsgrp,1,w(odmat),w(odmats))
        call query('V>=40 to display rs s',-1,0)
        call shostr(nl,nsite,nbas,1,plat,w(opos),11000,w(oalph),
     .    w(oiax),w(ontab),w(odmats),1,dum,1,1d0)
        call rx0('done showing dmat')
      endif

C --- e'vector decomp'n and moments from supplied weights ---
      if (nevmx .le. 0 .or. .not. metal .or. lwtkb .ne. 0) goto 1002
      if (efermi .gt. efmax) call fexit2(-1,111,' Exit -1 BNDASA: '//
     .  'E_f = %1;6d exceeds efmax = %1;6d',efermi,efmax)
C ... Reread header in moments file, resetting internal variables
      if (liomom .eq. 0)
     .  call rx('BNDASA: no moments file ... can''t make output moms')
      call iomomq(nfilem,2,nl,nsp,nspc,nkp,ldim,0,i,nbmax,nlo*nclass,
     .  nrhos*nclass,nev,w,w,w,w,efermi,xx)
      i = nlo*nclass*nevx*mod(nfstg/10,10)
C     this should never happen
      if (mod(nfstg/10,10) .ne. (nspc+2)*nspc)
     .  call rx('BNDASA: moments file mismatch ... bug in bndasa')
      call defdr(oaccwt,-i)
      i = nrhos*nclass*nevx*mod(nfstg/10,10)
      if (lrhos) call defdc(oaccsm,-i)
      if (lrhos) call dpzero(rhos,2*3*nrhos*2*2*nclass)
      call defdr(owk,nevx)
      swtkb = 0
      nschan = mod(nfstg/10,10)
      do  40  ikp = 1, nkp
      do  40  isp = 1, nsp
        if (onesp .and. isp .eq. 2) goto 40
        if (nspc .eq. 2 .and. isp .eq. 2) goto 40
        i = iomoms(nfilem,dum,nsp,nspc,1,ldim,nfstg,nschan,1,1,nevx,
     .    nevx,nlo*nclass,nrhos*nclass,nev,w(owk),w(oaccwt),w(oaccwt),
     .    w(oaccsm),dum,dum)
        i = isw(lrhos)
        call moment(i,nl,nlo,nrhos,isp,nsp,nspc,nkp,ldim,nevx,nev,ikp,
     .    w(owtkb),nclass,w(oaccwt),w(oaccsm),qnu,rhos,w(oorbtm),
     .    swtkb)
   40 continue

C     If qnu only l-dependent, approx patch moment in rhos
      if (lrhos) call amagn2(0,nl,nlo,nclass,nbas,w(oipc),w(oeula),neul,
     .  w(obxc),qnu,nrhos,rhos)

      if (dabs(zval-swtkb/nspc) .gt. 1d-6) then
        call awrit3(' (warning) sum-of-weights (%;6d) does not'//
     .    ' match valence charge (%;6d)'//
     .    '%?#n<0#%N Consider increasing NEVMX##',
     .    ' ',160,stdo,swtkb,zval,nevmx-ldimx)
      endif
      if (onesp) call dscal(3*nl*nsp*nclass,2d0,qnu,1)
 1002 continue
      call rlse(owtkb)
      if (liomom .ne. 0) then
        write (nfilem) efermi, vmtz(1)
        call fclose(nfilem)
      endif
      if (nevmx.le.0 .or. .not. metal .or. nspc.ne.2) then
      else
        call mpibc2(w(oorbtm),nlo*2*nclass,4,.false.,'bndasa','orbtm')
        call iorbtm(sspec,w(oics),nl,nlo,nclass,nsp,w(oorbtm))
        call rlse(oorbtm)
      endif

C --- Restore potential parameters to orthogonal rep'n ---
      call defdr(oo,nlspc)
      call pptrns(1,nl,w(oipc),nclass,nsp,dum,1,w(opp),w(oo))
      call rlse(oo)

      if (nevmx .gt. 0 .and. mod(loptic,10) .eq. 1
     .    .and. procid .eq. master) then
C        if (.not. ltet)
C     .    call rx('optical integrals require tetrahedron integration')
C       if (bittst(loptic,1) .or. loptic .lt. 0) then
        call defdr(owk1,            npopt)
        call defdr(oeps,-(1+nspx*3)*npopt)
        call defdr(owk3,          3*npopt)
        call optint(sctrl,w(osoptc),sbz,eband,nbmaxx,1,nspx,nspc,
     .    -1,efermi,w(oidtet),w(owtkp),dgets('lat vol',slat),nfilm,
     .    nempm,npopt,w(owk1),w(owk3),w(ooptmt),w(oeps))
        call optin2(sctrl,w(osoptc),sbz,nspx,nspc,
     .    dgets('lat vol',slat),npopt,w(oeps))
        call rlse(owk1)
C       call veltet(n1,n2,n3,eband,nkp,nbmax,nsp,efermi,ntet,w(oidtet),
C    .    nfilo,nfiup,nfilm,w(ovelmt))
      endif

      if (bitand(loptic/10,2) .ne. 0) then
        call info0(30,0,0,' BNDASA: nonlinear optics')
        call intdrv(sctrl,w(osoptc),sbz,eband,nbmax,nsp,nspc,efermi,
     .    w(oidtet),dgets('lat vol',slat),nfilm,nempm,w(ooptmt))
      endif
C     call veltet(n1,n2,n3,eband,nkp,nbmax,nsp,efermi,ntet,w(oidtet),
C    .  nfilo,nfiup,nfilm,w(ovelmt))

C --- Make local DOS for Stoner magnetism ---
#if STONER
       if (lgors('ctrl lstonr,1',sctrl)) then
         stop 'upack zos,stoner index'
        if (nlo .ne. nl)
     .    call rx('BNDASA: (class,l,m) moments not allowed with STONER')
        call rxx (nl .lt. 3,'BNDASA: nl must be > 2 for Stoner')
        if (iprint() .ge. 30) call awrit3(
     .   ' BNDASA: Make number-of-states functions for Stoner model'//
     .   ', %i bins in (%d,%d)',' ',100,stdo,ndos-1,dosw(1),
     .   dosw(2))
        i = 0
        do  14  ic = 0, nclass-1
        do  14  l = 0, nl-1
        if (idxdn(l+nl*ic) .ne. 3) i = i+1
   14   if (l .eq. 2) index(ic) = i-1
C       write (*,157) (ic, index(ic), ic = 0, nclass-1)
C 157   format ('ic=',i3,' index=',i3)
        if (liomom .ne. 0) call rx('BNDASA: no moments file')
        nfilem = fopna('MOMS',-1,4)
        call defdr(owk,ndos)
        nchan = nlo*nclass
        if (ltet) then
          call defdr(oeband,ldim*nsp*nkp)
          call defdr(odoswt,-nchan*ldim*nsp*nkp)
          nfstg = 11
          call iomomq(nfilem,22,nl,nsp,nspc,nkp,ldim,nfstg,i,ldim,
     .      nchan,nchan,nev,w(oeband),w,w(odoswt),w,efermi,xx)
          call rxx(i.ne.nkp,'BNDASA:  moments file missing qpts')
          call dostet(ldim,nsp,nsp,nev,nchan,n1,n2,n3,ntet,
     .      w(oidtet),w(oeband),w(odoswt),ndos,dosw(1),
     .      dosw(2),.true.,w(owk),zos)
        else
          call defdr(oeband,ldim)
          call defdr(odoswt,nchan*ldim)
          nfstg = 0
          call iomomq(nfilem,12,nl,nsp,nspc,nkp,ldim,nfstg,i,ldim,
     .      nchan,nchan,nev,w,w,w,w,efermi,xx)
          call rxx(i.ne.nkp,'BNDASA:  moments file missing qpts')
          call dosspl(nfilem,ldim,nsp,1,nchan,mpsord,swidth,nkp,
     .      w(owtkp),w(oeband),w(odoswt),ndos,dosw(1),dosw(2),.true.,
     .      w(owk),zos)
        endif
        call rlse(owk)
C   ... Scale zos to get per spin
        call dscal(ndos*nchan,0.5d0,zos,1)
      endif
#else
      if (lgors('ctrl lstonr,1',sctrl))
     .  call rx('BNDASA: recompile with ccomp -dSTONER ...')
#endif

C --- Generate DOS on disk ---
      if (lgors('ctrl ldos,1',sctrl)) then
        call defdr(odos,3*ndos)
        if (iprint() .ge. 30) call awrit1('%x%N ... Generating %?#n<0#'
     .    //'integrated#total# DOS',' ',80,stdo,ndos0)
        if (ltet) then
          call bzints(n1,n2,n3,eband,dum,nkp,ldimx,nbmaxx,nspx,
     .      dosw(1),dosw(2),w(odos),ndos,efermi,1,ntet,w(oidtet),
     .      dum,xx)
          if (ndos0 .gt. 0)
     .      call xxxdif(dosw(1),dosw(2),ndos,nspx,0,w(odos))
          del = 0d0
        else
          if (mpsord .ge. 100) mpsord = mod(mpsord,100)
          if (ndos0 .gt. 0)
     .      call makdos(nkp,ldimx,nbmaxx,nspx,w(owtkp),eband,mpsord,
     .      swidth,-srnge,dosw(1),dosw(2),ndos,w(odos))
          if (ndos0 .lt. 0)
     .      call maknos(nkp,ldimx,nbmaxx,nspx,w(owtkp),eband,mpsord,
     .      swidth,-srnge,dosw(1),dosw(2),ndos,w(odos))
          del = mpsord+swidth
        endif
        if (nspc .eq. 2) call dscal(ndos,.5d0,w(odos),1)
        call dosio(w(odos),ndos,nspx,ndos,1,dosw(1),dosw(2),nspx,
     .    efermi,del,1,-fopn('DOS'))
        call fclose(fopn('DOS'))
        call rlse(odos)
      endif

C --- Magnetic moments corresponding to (output) density matrix ---
      if (lrhos) then
        if (iprt(2) .ne. iprt(1)) then
          call pshpr(iprt(2))
        else
          call pshpr(iprint()-0)
        endif
C       call setpr(50)
        call amagnc(nbas,nl,w(oipc),rhos,nrhos,qnu,w(oeula),neul,0,amag,
     .    aamom)
        if (bittst(lncol,8)) then
          call bdotsr(nbas,nl,w(oipc),rhos,nrhos,qnu,
     .      w(obdots),lihdim,w(oidxsh),0,dum)
C         ehterm(6) = dum(1)
          eterms(18) = dum(1)
        endif
        if (neul .gt. 1 .or. bittst(lncol,8)) then
          if (.not. bittst(lncol,8)) nbf = 99
          call bsrhos(nbas,nl,w(oipc),rhos,nrhos,qnu,w(opp),
     .      w(osop),w(oeula),neul,w(obxc),w(omagf),nbf,lihdim,w(oidxsh),
     .      0,dum)
          eterms(18) = dum(2)
        endif
       call poppr
      endif

C --- Make Magnetization of r
C      call asvsph2(sctrl,slat,sarray,sspec,sham,spot,rhos,
C     .  max(neul,nl),1,0,1) 
C ... To check quantities integrated throughout the sphere:
c      call asvsphint(sctrl,slat,sarray,sspec,sham,spot,rhos,
c     .  max(neul,nl),1,0,1) 

C     call info0(30,0,0,' ')
      call rlse(o)

C ... Repack d.c. terms from applied field
      call pack1('ham eterms',sham,eterms)

C ... Restore struc variables we changed
      call pack2('ctrl lham lasa',sctrl,lhams,lasas)
      call pack2('ham lham',sham,lhams,0)

      call tcx('bndasa')
      end

C      subroutine gethkp(nsp,isp,nkp,ikp,ldim,hk,hkp)
C      
C      integer nsp,isp,nkp,ikp,i,j
C      double precision hk(ldim,ldim*2)
C      double precision hkp(ldim*2,ldim)
C      
C      call dpzero(hkp,ldim*ldim*2)
C      do i = 1, ldim
C        do j = 1, ldim
C          
C          hkp(i*2-1,j) = hk(i,j)
C          hkp(i*2,j) = hk(i,j+ldim)
C        end do
C      end do
C
C      end 

