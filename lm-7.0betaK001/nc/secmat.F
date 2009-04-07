#define NC 1
#define BLAS3 1
      subroutine secmat(sctrl,sspec,sham,spot,sarray,slat,ikp,nkp,
     .  qp,wtkp,isp,lrsig,nevmx,efmax,nlibu,lmaxu,idu,ludiag,vorb,
     .  nev,z,eb,strx,rhrs,rors)
C- Set up Hamiltonian and Overlap and diagonalize secular matrix.
C ----------------------------------------------------------------------
Ci Inputs
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read: nbas nl nspin nclass lgen3 lham lncol loptc
Ci     Stored:    *
Ci     Passed to: secmtn lgors secm2c secmt2
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: *
Ci     Stored:    *
Ci     Passed to: secmtn
Ci   sham  :struct for parameters defining hamiltonian; see routine uham
Ci     Elts read: oeula neula qss obdots oindxo ldham ndhrs ontabs
Ci                oiaxs ohrs
Ci     Stored:    *
Ci     Passed to: secmtn
Ci   spot  :struct for information about the potential; see routine upot
Ci     Elts read: oppn osop opp opti vmtz
Ci     Stored:    *
Ci     Passed to: secmtn
Ci   sarray:structure containing offsets to various arrays
Ci     Elts read: oips oipc ormax
Ci     Stored:    *
Ci     Passed to: *
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: plat opos
Ci     Stored:    *
Ci     Passed to: secmtn
Ci   ikp   :k-point label
Ci   nkp   :number of irreducible k-points (bzmesh.f)
Ci          Needed mainly for nonlocal exchange potential
Ci   qp    :k-point
Ci   wtkp  :k-point weights, including spin degeneracy (bzmesh.f)
Ci         :(only used with --invbl to make rsrs,rors)
Ci   isp   :current spin channel (1 or 2)
Ci   lrsig :0 no nonlocal (screened-exchange) potential
Ci         :1 or 2 add screened-exchange potential sigm
Ci         :       as bloch sum of sham->hrs
Ci         :-1     add screened-exchange potential sigm
Ci         :       by reading sigma from disk file 'sigm'
Ci   nevmx :largest number of eigenvectors to find
Ci   efmax :largest eigenvalue for which to find eigenvectors
Co Outputs
Co   nev   :actual number of eigenvectors generated
Co   z     :eigenvectors
Co   eb    :energy bands; alias eband
Co   strx,rsrs,rors : inverse Bloch transformed S, H and O (case --invbl)
Cr Remarks
Cu Updates
Cu   13 Nov 08 (ATP) Inverse Bloch transform
Cu   15 Nov 07 Improved LDA+U; works in noncollinear case
Cu   08 Nov 07 (J. Xu) LDA+U, first cut
Cu   18 Jun 04 (A Chantis) fully relativistic hamiltonian
Cu   10 Oct 03 SX sigma can be spin polarized
Cu   23 Sep 03 SX patterned after GW.  sigm(rs) stored in sham->hrs
Cu   14 Feb 03 Added applied magnetic field
Cu   03 Oct 01 bug fix, 2nd gen case when gam-rep & neglected orbitals
Cu   30 Aug 00 Added NMTO ASA Hamiltonian
Cu   14 Sep 99 rewrite argument list using structures
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer isp,nevmx,nev,ikp,nkp,lrsig
      double precision sctrl(1),sspec(1),sham(1),spot(1),sarray(1),
     .  slat(1)
      double precision eb(*),qp(3),wtkp(nkp),z(*),efmax
C     For LDA+U
      integer nlibu,lmaxu,idu(4,*),ludiag
      double complex vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,2,nlibu)

C ... Local parameters
C     logical lgors
      integer lgen3,ldham(16),nbas,nclass,nl,nsp,neul,ldim,
     .  lidim,lihdim,fopnT,nfilet,ndim,nglob
      equivalence (ldim,ldham(1)),(lidim,ldham(2)),(lihdim,ldham(3))
      integer oips,oipc,oclabl,oeula,opp,opti,opph,oppn,oidxsh,osop,
     .  obas,ormax,obs
      double precision qss(4),vmtz,cksumf,ckbas,plat(3,3),avw,dglob
C     For self-energy sigma
      integer ndhrs
      integer ontabs,oiaxs,ohrs
C     For inverse Bloch transform
      double precision strx(1),rhrs(1),rors(1)

C heap:
      integer w(1)
      common /w/ w

      call tcn('secmat')

C --- Setup ---
      call upack('ctrl nbas nl nspin nclass lgen3',sctrl,nbas,nl,nsp,
     .  nclass,lgen3)
      call upack('array oips oipc ormax oclabl',sarray,oips,oipc,ormax,
     .  oclabl,0)
      call upack('ham oeula neula qss obdots',sham,oeula,neul,qss,obs,0)
      call upack('ham oindxo ldham',sham,oidxsh,ldham,0,0,0)
      call upack('lat plat opos',slat,plat,obas,0,0,0)
      ckbas = cksumf(w(obas),3*nbas)
      avw = dglob('avw',0d0,0)
      ndim = nbas * nglob('mxorb')

      if (lgen3 .ne. 0) then
        call upack('pot oppn osop',spot,oppn,osop,0,0,0)
        call secmtn(sctrl,sspec,slat,sham,spot,nl,nsp,nbas,w(oips),
     .    w(oidxsh),qss,w(oeula),neul,ikp,nkp,ldim,lidim,lihdim,qp,
     .    w(oppn),w(osop),isp,nevmx,efmax,nev,z,eb)
      else
C       nfilet = fopn('TMP')
        nfilet = fopnT('TMP',-1,4,0)
        rewind nfilet
        call upack('pot opp osop opti vmtz',spot,opp,osop,opti,vmtz,0)
        call defdr(opph,5*lihdim*nsp)
        call makpph(nl,nsp,nbas,lihdim,w(oipc),w(oidxsh),w(opp),w(opph))
        call u2pph(0,nbas,lmaxu,nsp,nlibu,idu,vorb,w(oidxsh),
     .    ldim,lihdim,1d0,w(opph))
        call upack('ham ndhrs ontabs oiaxs ohrs',sham,ndhrs,ontabs,
     .    oiaxs,ohrs,0)
C       Fast version works, but comment out for generality; use secmt2
C        if (lgors('ctrl lham,3',sctrl)) then
C          if (lrsig .ne. 0) call rx('update lrsig for call to 2C')
C          call secm2c(sctrl,ckbas,plat,nl,nsp,nbas,nclass,w(oipc),
C     .      w(oidxsh),qss,w(oeula),neul,ikp,nkp,ldim,lidim,lihdim,ndim,
C     .      qp,w(opp),w(osop),vmtz,w(ormax),avw,isp,w(opti),w(opph),
C     .      nfilet,nevmx,efmax,nev,lrsig,z,eb)
C        else
        call secmt2(sctrl,spot,ckbas,plat,nl,nsp,nbas,nclass,w(oipc),
     .    nlibu,lmaxu,idu,ludiag,vorb,w(oclabl),
     .    w(oidxsh),qss,w(oeula),neul,w(obs),ikp,wtkp,nkp,ldim,lidim,
     .    lihdim,ndim,qp,w(opp),w(osop),vmtz,w(ormax),avw,isp,w(opti),
     .    w(opph),nfilet,nevmx,efmax,nev,
     .    lrsig,ndhrs,w(ontabs),w(oiaxs),w(ohrs),
     .    z,eb,strx,rhrs,rors)
C        endif
        call rlse(opph)
      endif
      call tcx('secmat')
      end
      subroutine secmt2(sctrl,spot,ckbas,plat,nl,nsp,nbas,nclass,ipc,
     .  nlibu,lmaxu,idu,ludiag,vorb,dclabl,indxsh,qss,
     .  eula,neul,bdots,ikp,wtkp,nkp,ldim,lidim,lihdim,ndim,qp,pp,
     .  sop,vmtz,wsr,avw,isp,pti,pph,nfilet,nevmx,efmax,nev,lrsig,ndsig,
     .  ntabs,iaxs,sigrs,z,eb,strx,rhrs,rors)
C- Hamiltonian and Overlap, 2nd generation lmto
C ----------------------------------------------------------------
Ci Inputs:
Ci   sctrl :struct containing parameters governing program flow
Ci     Elts read: lham lncol loptc
Ci   ckbas :checksum for basis vectors;  make with cksum(bas,3*nbas)
Ci   plat  :primitive lattice vectors, in units of alat
Ci   nl    :(global maximum l) + 1
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nbas  :size of basis
Ci   nclass:number of inequivalent classes
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   indxsh:permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   qss   :Parameters defining spin-spiral
Ci   eula  :Euler angles for noncollinear spins
Ci   neul  :1 if Euler angles are l-independent, nl otherwise
Ci   bdots :(magnetic field . Pauli spin matrices), downfolding order
Ci   ikp   :k-point label
Ci   wtkp  :k-point weights, including spin degeneracy (bzmesh.f)
Ci   nkp   :number of irreducible k-points (bzmesh.f)
Ci          Needed mainly for nonlocal exchange potential
Ci   ldim  :dimension of hamiltonian matrix (makidx.f)
Ci   lidim :number of lower+intermediate orbitals
Ci   lihdim:number of lower+intermediate+higher orbitals
Ci   ndim  :nl*nl*nbas
Ci   qp    :k-point
Ci   pp    :potential parameters (atomsr.f)
Ci   sop   :spin-orbit parameters (atomsr.f)
Ci   vmtz  :muffin-tin zero (asamad.f)
Ci   wsr   :Wigner-Seitz radius, in a.u. (input; alias rmax)
Ci   avw   :length scale, usu. average Wigner-Seitz sphere radius
Ci   isp   :current spin channel (1 or 2)
Ci   pti   :inverse potential functions needed for downfolding
Ci   pph   :potential parameters in downfolding order (makpph.f), alpha rep'n
Ci   nfilet:logical units for temporary file.  File has structure:
Ci         sil (idim>0 and nevmx>0)
Ci         i-wave contr. to h (idim>0 and nspc=2)
Ci         sll, beta repsn (nevmx>0)
Ci   nevmx :max. no. evecs to generate.
Ci         -1 suppresses generation of z
Ci         -2 Do not diagonalize, but return overlap in z,
Ci            allocate oc for hamiltonian and place there
Ci   efmax :largest eigenvalue for which to find eigenvectors
Co Outputs:
Ci   ccd:  diagonal matrices for 1- 2- & 3-centre CCOR integrals
Co   eigenvalues and eigenvectors are returned in eb, z
Co   nev:  number of evecs generated.
Co   eb    :energy bands
Co   z     :eigenvectors.
Co         :z is used as a work array, whether or not evecs generated
Cr Remarks
Cr  *Orbitals may be folded down as described in the Kanpur notes or
Cr   they may be discarded before diagonalising. Such orbitals are
Cr   indicated by IDXDN(i) = 2 or 3 respectively. If the automatic
Cr   downfolding switch is set, orbitals are assigned to a set using
Cr   a set of rules determined in ATFOLD.
Cr
Cr   Folding down is about the inverse, unscreened potential function
Cr   at e_nu, and this is the only choice allowed (the choice uncouples
Cr   the moments of the charge density of the lower and intermediate
Cr   sets, and is needed in order to ensure a representation-independent
Cr   construction of the i-wave eigenvectors, 3-centre integrals
Cr   and beta^dot).
Cr
Cr   Downfolding automatically turns on the combined correction.
Cr
Cr   bittst(lham,8) can be used to transform
Cr   structure constants to an arbitrary representation.
Cr
Cr   Hybridisation is turned off when bittst(lham,16) is set
Cr
Cr   Dimensions of pph,eb,z are doubled when spins are coupled.
Cu Updates
Cu   14 Feb 03 Added applied magnetic field
Cu   01 Jun 01 bug fix, iwave avg (noncol case) when lihdim<ndim
Cu   14 Sep 99 renamed from original secmat.
Cu   28 Apr 98 uses spherical harmonics when lsph set.
C ----------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nl,nsp,nbas,nclass,neul,isp,nfilet,ldim,lidim,lihdim,ndim,
     .  nevmx,nev,ipc(nbas),indxsh(ndim),ikp,nkp,lrsig,ludiag
      double precision ckbas,plat(3,3),eb(ldim*2),pp(6,nl,nsp,*),vmtz,
     .  avw,qp(3),z(ldim,ldim*2),pti(ndim,nsp),eula(nbas,neul,3),wsr(*),
     .  pph(5,lihdim,nsp),efmax,qss(4),sop(0:nl-1,nsp,nsp,9,*),sctrl(1),
     .  spot(1),dclabl(1),wtkp(nkp)
C     For LDA+U
      integer nlibu,lmaxu,idu(4,nbas)
      double complex vorb(-lmaxu:lmaxu,-lmaxu:lmaxu,nsp,nlibu)
C     For inverse Bloch transform      
      double precision wtqp,strx(1),rhrs(1),rors(1)
C     double complex bdots(2,2,lihdim)
      double precision bdots(2,2,2,lihdim)
C     For SX
      integer ndsig,ntabs(nbas+1),iaxs(*)
      double precision sigrs(2,ndsig,ndsig,nsp,*)
C ... Local parameters
      logical iostr,ndnfld,ccor,lss,lnc,lso,lbf,lx,F,T,bittst,lgors,ltmp
      parameter (F=.false., T=.true.)
      integer bit,i,ifi,i1mach,i2,idim,ii,iprint,j,l2,fopna,
     .  lbloch,ld2,ldimx,lgunit,lham,li,linv,lncol,loptic,lov,n,
     .  nl2,nlspc,nsite,ncsw
      integer oa,oacopy,oadot,oalph,occd,odiawk,ogma,oh,oiax,ontab,oo,
     .  os,osii,osil,osll,osod,osoph,owk,owk2,owk3,oxsidt,opprel,oohp1
#if IREP
      integer oxsi,onewa
#endif
C     For inverse Bloch transform
      integer oiaxb,nsmax
      parameter (nsmax=5000)
      logical invb,cmdopt
      character*20 outs
      double precision kap2(20),qpq(3)
C     For self-energy sigma
      integer hreal,nttabs,mxorb,nglob,lrel,igets

C external calls:
      external dpzero,defi,defdr,rlse,dscal,diagno,igets,
     .   iostr,bloch,lgunit,trS,trSdot,invbl
C heap:
      integer w(1)
      common /w/ w
      bittst(n,bit) = (mod(n,bit+bit) - mod(n,bit) .eq. bit)

C     print *, '!!'
C     qp(1) = .1d0
C     qp(2) = .2d0
C     qp(3) = .3d0

C --- Setup ---
      invb = cmdopt('--invbl',7,0,outs)

      call upack('ctrl lham lncol loptc',sctrl,lham,lncol,loptic,0,0)
      nl2 = nl**2
      idim = lidim - ldim
C     hdim = lihdim - lidim
      nlspc = nl * nsp * nclass

C     Noncollinear switches
      lnc = lncol .ne. 0
      lss = bittst(lncol,2)
      lso = bittst(lncol,4)
      lbf = bittst(lncol,8)
C     In noncollinear case use average of up+down pot. functions.
C     print *, '!!';  if (lnc .or. nsp .eq. 2) then
      if (lnc) then
        call saiwav(ldim,lidim,1,ndim,pti)
        call saiwav(ldim,lidim,5,lihdim,pph)
      endif

C     Possibly rotate to spherical harmonics when making Bloch sum
      lbloch = 0
      if (bittst(lham,256)) lbloch = 1000

C     Combined correction required if downfolding
      ndnfld = idim .eq. 0
      ccor = .not. ndnfld .or. lgors('ctrl lasa,4',sctrl)
      if (.not. ccor  .and.  iprint()  .ge. 30  .and.  ikp .eq. 1)
     . print *, 'SECMAT : Combined Correction switched off'

C     Diagonalize by inverse iteration, or not
      linv = 0
      if (nevmx .gt. 0 .and. lgors('ctrl lqp,2',sctrl)) linv = 1

C ... Sanity checks
      if (.not. ndnfld .and. bittst(lham,128))
     .  call rx('SECMAT: no downfolding in gamma rep')

C     Some dimensioning parameters and memory allocation
      ldimx = ldim
      ld2 = ldim**2
      if (lnc) ldimx = 2*ldim
      l2 = ldimx**2
      i2 = idim**2
      li = ldim * idim
      call zinit(z,l2)
      if (ccor) call defdr(occd,-3*lihdim)
      if (.not. ccor) call defdr(occd,1)
      call defdc(osll,l2)
      call defdr(owk,2*ndim)
      if ( .not.  ndnfld ) call defdc(oh,-l2)
      if (lss) then
        i = 2
      else
        i = 1
      endif
C     call wkprnt(1)
      call defdc(osii,i2*i)
      call defdc(osil,li*i)
      call defdc(oa,li*i)

C --- Get screened strux from disc and Bloch-transform them ---
C     w(oalph) is alpha(1..nl,nbas)
      if (invb) call defi(oiaxb,10*nsmax)
      ltmp = iostr(8,'STR',nl,nbas,1,kap2,0,ckbas,-1,nsite,oalph,
     .  oiax,ontab,os)
      nsite = w(ontab+nbas)
      call dcopy(3,qp,1,qpq,1)
C ... make copy of iax for invbl
      if (invb) then
        call rxx(nsite.gt.nsmax,' Increase nsmax in SECMT2')
        call icopy(10*nsite,w(oiax),1,w(oiaxb),1)
      endif
      if (lss) then
        qpq(1) = qp(1) + qss(1)/2
        qpq(2) = qp(2) + qss(2)/2
        qpq(3) = qp(3) + qss(3)/2
        call bloch(lbloch,qpq,nl,plat,nl**2,indxsh,1,nsite,w(oiax),
     .    w(os),nl2,1,1,ldim,ldim,idim,ldim,idim,ldim,0,w(osll),w(osil),
     .    w(osii))
C       call yprm('Sll(2)',12,w(osll),ldim*ldim,ldim,ldim,ldim)
        call pvsec2(2*ld2,w(osll))
        call pvsec2(2*li, w(osil))
        call pvsec2(2*i2, w(osii))
        qpq(1) = qp(1) - qss(1)/2
        qpq(2) = qp(2) - qss(2)/2
        qpq(3) = qp(3) - qss(3)/2
      endif

      call bloch(lbloch,qpq,nl,plat,nl**2,indxsh,1,nsite,w(oiax),w(os),
     .  nl2,1,1,ldim,ldim,idim,ldim,idim,ldim,0,w(osll),w(osil),w(osii))

C      if (hdim .ne. 0) then
C        call defcc(oshl, hdim*ldim)
C        call bloch(lbloch+50000,qpq,nl,plat,nl**2,indxsh,1,nsite,w(oiax)
C     .    ,w(os),nl2,1,1,lidim,ldim,hdim,1,hdim,ldim,0,w,w(oshl),w)
C      endif

      call rlse(oiax)

C ... debugging ... put sk into z and exit
C      print *, 'debugging ... put zk into z'
C      call dcopy(ld2*2,w(osll),1,z,1)
C      call rlse(osll)
C      return

C --- Transform potential parameters to gamma representation ---
      if (bittst(lham,128)) then
        call defdr(ogma,ndim)
C   ... Put pp in gamma; make vector of gamma's
        call defdr(oo,nlspc)
        i = -3
        if (isp .eq. 2) i = -4
        if (lnc) i = -5
        if (lgors('ctrl lasa,512',sctrl)) i = -5
        call pptrns(i,nl,ipc,nclass,nsp,w(ogma),nbas,pp,w(oo))
C       call prmx('gamma',w(ogma),ndim,ndim,1)
C       call prmx('alpha',w(oalph),ndim,ndim,1)
        call rlse(oo)
        call makpph(nl,nsp,nbas,lihdim,ipc,indxsh,pp,pph)
        call u2pph(0,nbas,lmaxu,nsp,nlibu,idu,vorb,indxsh,ldim,
     .    lihdim,1d0,pph)
C   ... put gamma-alpha in gam, gamma in alpha
        call daxpy(ndim,-1d0,w(oalph),1,w(ogma),1)
        call daxpy(ndim, 1d0,w(ogma),1,w(oalph),1)
C       call prmx('gamma-alpha',w(ogma),ndim,ndim,1)
C --- Or else transform potential parameters to alpha rep'n ---
      else
        call defdr(oo,nlspc)
        call pptrns(0,nl,ipc,nclass,nsp,w(oalph),1,pp,w(oo))
        call rlse(oo)
      endif

C --- Change representation interactively ---
#if IREP
      if (bittst(lham,8)) then
        call rxx(bittst(lham,128),'already changed to gamma rep')
        call rxx(.not. ndnfld,'NEW_REP: unset downfolding switches')
        if (ikp.eq.1) then
            print *, 'SECMAT : Transforming to new representation'
        endif
        call defdr(onewa,ldim)
        call defdr(oo,nlspc)
        call newalp(ldim,nl,nbas,ipc,nsp,isp,pp,w(onewa))
        call pptrns(0,nl,ipc,nclass,nsp,w(onewa),1,pp,w(oo))
        call rlse(oo)
      endif
#else
      if (bittst(lham,8)) call rx('uncomment IREP in secmat')
#endif

C --- S^alpha -> S^beta ---
      if (.not. ndnfld) then

C   ... Make a copy of screening alphas, transform i-waves to beta=pti
        call defdr(oacopy,ndim)
        call dcopy(ndim,w(oalph),1,w(oacopy),1)
        call makbet(ndim,ldim,idim,indxsh,pti(1,isp),w(owk),w(oalph))

C   ... Transform S^alpha to S^beta for downfolding the i-waves
        call shffle(.true.,ndim,indxsh,w(oacopy),w(owk))
        if (lss) then
C         Transform s(2)=s(q+qss/2) first
          call pvsec2(2*ld2,w(osll))
          call pvsec2(2*li, w(osil))
          call pvsec2(2*i2, w(osii))
          call pvsec2(2*li, w(oa))
          call defdr(owk2,ndim)
          call dcopy(ndim,w(oacopy),1,w(owk2),1)
          call trS(ldim,idim,w(owk2),pti(1,isp),
     .      w(osll),w(osil),w(osii),w(oa),w(owk))
C         call yprm('Sil2 (beta)',02,w(osil),idim*ldim,idim,idim,ldim)
          call rlse(owk2)
C         Restore s(1)=sdot(q-qss/2) to first array, and transform
          call pvsec2(2*ld2,w(osll))
          call pvsec2(2*li, w(osil))
          call pvsec2(2*i2, w(osii))
          call pvsec2(2*li, w(oa))
        endif
        call trS(ldim,idim,w(oacopy),pti(1,isp),
     .    w(osll),w(osil),w(osii),w(oa),w(owk))
        if (iprint() .ge. 110) then
          call yprm('Sll (beta)',02,w(osll),ldim*ldim,ldim,ldim,ldim)
          call yprm('Sil (beta)',02,w(osil),idim*ldim,idim,idim,ldim)
        endif
      endif


C --- Save S_il in temporary file, if needed ---
      if (nevmx .gt. 0) then
        i = ldim
        if (lss) i = ldim*2
        if (idim .gt. 0) call dpdump(w(osil),idim*i*2,-nfilet)
      endif

C --- Make 3-centre hamiltonian integrals over the i-waves ---
      if (idim .gt. 0) then
       call i3cntr(ldim,idim,lnc,w(osil),pph(1,1,isp),vmtz,w(oh),w(owk))
       if (lnc) then
C        SS: 3-center terms for s(q+qss)
         if (lss) then
           call pvsec2(2*li, w(osil))
           call pvsec2(2*ld2, w(oh))
           call i3cntr(ldim,idim,lnc,w(osil),pph(1,1,isp),vmtz,w(oh),
     .       w(owk))
           call pvsec2(2*li, w(osil))
           call pvsec2(2*ld2, w(oh))
         endif
         call defrr(owk2,ldim*2)
         call dcopy(ldim,pph(3,1,2),5,w(owk2),1)
         call dpscop(w(owk2),w(owk2),ldim,1,1+ldim,1d0)
         call dcopy(ldim,pph(3,1,1),5,w(owk2),1)
         ncsw = 2000
         if (lss) ncsw = ncsw + 20000
         call rotspn(ncsw,nbas,nl,indxsh,eula,neul,qss(4),
     .     w(owk2),w(owk2),ldim,ldim,ldim,ldim,ldim,w(oh),z)
         call dpdump(z,l2*2,-nfilet)
         call rlse(owk2)
         if (iprint() .ge. 110)
     .     call yprm('iwaves H',12,z,ldimx**2,ldimx,ldimx,ldimx)
        endif
      endif

C --- Make S-dot(alpha) for the combined correction; store in z ---
      if (ccor) then
        ltmp = iostr(8,'SDOT',nl,nbas,1,kap2,0,ckbas,-1,nsite,
     .    oadot,oiax,ontab,os)
        nsite = w(ontab+nbas)
        if (lss) then
          qpq(1) = qp(1) + qss(1)/2
          qpq(2) = qp(2) + qss(2)/2
          qpq(3) = qp(3) + qss(3)/2
          call bloch(lbloch,qpq,nl,plat,nl**2,indxsh,1,nsite,w(oiax),
     .      w(os),nl2,1,1,ldim,ldim,idim,ldim,idim,ldim,0,z(1+2*ld2,1),
     .      w(osil),w(osii))
          call pvsec2(2*li, w(osil))
          call pvsec2(2*i2, w(osii))
C         call yprm('Sd(2)',12,z(1+2*ld2,1),ldim*ldim,ldim,ldim,ldim)
          qpq(1) = qp(1) - qss(1)/2
          qpq(2) = qp(2) - qss(2)/2
          qpq(3) = qp(3) - qss(3)/2
          call bloch(lbloch,qpq,nl,plat,nl**2,indxsh,1,nsite,w(oiax),
     .      w(os),nl2,1,1,ldim,ldim,idim,ldim,idim,ldim,0,z,w(osil),
     .      w(osii))
        else
          if (iprint() .ge. 100) print *, '... bloch sum sdot'
          call bloch(lbloch,qpq,nl,plat,nl**2,indxsh,1,nsite,w(oiax),
     .      w(os),nl2,1,1,ldim,ldim,idim,ldim,idim,ldim,0,z,w(osil),
     .      w(osii))
        endif
        call rlse(oiax)
      endif

C --- Transform sll, sdot to gamma representation ---
      if (bittst(lham,128)) then
        call defdr(owk3,2*ld2)
        call defdr(oxsidt,-ndim)

C   ... SS case: s(q+qss/2) also rotated to gamma
        if (lss) then
C         Use copies of gamma-alpha, xsidot since mksbet overwrites them
          call dpcopy(w(ogma),w(owk),1,ndim,1d0)
          call defdr(owk2,ndim)
          call dpcopy(w(oxsidt),w(owk2),1,ndim,1d0)
C         swap 1st and 2nd B.T. sll, to avoid more workspace
          call pvsec2(2*ld2,w(osll))
          if (ndim .ne. ldim) call rx('secmat: need shuffle wk')
          call mksbet(ccor,ldim,w(owk),w(owk2),w(osll),z(1+2*ld2,1),
     .      w(owk3))
          if (iprint() .ge. 110) then
            call yprm('Sll(2) (gamma)',02,w(osll),ldim*ldim,ldim,ldim,
     .        ldim)
            if (ccor) call yprm('Sdt(2) (gamma)',02,z(1+2*ld2,1),
     .        ldim*ldim,ldim,ldim,ldim)
          endif
          call pvsec2(2*ld2,w(osll))
        endif

        call shffle(.true.,ndim,indxsh,w(ogma),w(owk))
        call shffle(.true.,ndim,indxsh,w(oxsidt),w(owk))
        call mksbet(ccor,ldim,w(ogma),w(oxsidt),w(osll),z,w(owk3))
        call shffle(.false.,ndim,indxsh,w(ogma),w(owk))
        call shffle(.false.,ndim,indxsh,w(oxsidt),w(owk))

        if (iprint() .ge. 110) then
          call yprm('Sll (gamma)',02,w(osll),ldim*ldim,ldim,ldim,ldim)
          if (ccor) call yprm('Sdt (gam)',02,z,ldim*ldim,ldim,ldim,ldim)
        endif
        call rlse(owk3)
      endif

C --- Change representation interactively ---
#if IREP
      if (bittst(lham,8)) then
        call defdr(oxsi,ldim)
        call dmadd(w(onewa),1,1,1d0,w(oalph),1,1,-1d0,w(oxsi),1,1,
     .             ldim,1)
        call defdr(oxsidt,-ldim)
        call defdr(owk3,2*ld2)
        call mksbet(ccor,ldim,w(oxsi),w(oxsidt),w(osll),z,w(owk3))
        call rlse(owk3)
        if (iprint() .ge. 110 .or. F) then
          call yprm('Sll (new rep)',02,w(osll),ldim*ldim,ldim,ldim,ldim)
          if (ccor) call yprm('Sdt',02,z,ldim*ldim,ldim,ldim,ldim)
        endif
        call rlse(oxsi)
        call dcopy(ldim,w(onewa),1,w(oalph),1)
      endif
#endif

C ... Save S_ll in temporary file
      if (nevmx .gt. 0) call dpdump(w(osll),ldim*i*2,-nfilet)

C --- Rotate S-dot(alpha) to S-dot(beta) for downfolding ---
      if (ccor) then
C   ... Diagonal matrices for 1- 2- and 3-centre CCOR integrals
        call makdia(nl,nbas,lihdim,indxsh,lidim,ipc,wsr,avw,w(oalph),
     .    w(oadot),w(occd))
        call dscal(3*lihdim,avw**2,w(occd),1)

        if (.not. ndnfld) then
C     ... w^-2(beta^dot - alpha^dot) for the i-waves
          call shffle(.true.,ndim,indxsh,w(oadot),w(owk))
          call mkbdot(ldim,idim,lihdim,lnc,w(occd),pph(1,1,isp),
     .      avw,w(oadot))

C     ... Transform w^-2 S^alpha^dot to the beta representation
          if (lss) then
C           Transform sdot(2)=s(q+qss/2) first
            call pvsec2(2*li, w(osil))
            call pvsec2(2*i2, w(osii))
            call pvsec2(2*li, w(oa))
            call defdr(owk2,lidim)
            call dpcopy(w(oadot),w(owk2),1,lidim,1d0)
            call trSdot(ldim,idim,w(oacopy),w(owk2),w(oa),
     .        z(1+2*ld2,1),w(osil),w(osii))
            call rlse(owk2)
C           Restore sdot(1)=sdot(q-qss/2) to first array
            call pvsec2(2*li, w(osil))
            call pvsec2(2*i2, w(osii))
            call pvsec2(2*li, w(oa))
          endif
          call trSdot(ldim,idim,w(oacopy),w(oadot),w(oa),
     .      z,w(osil),w(osii))
          call rlse(osii)
        endif

C   ... Remove w^-2 scaling from Sdot
        call dscal(2*ld2,-avw**2,z,1)
        if (lss) call dscal(2*ld2,-avw**2,z(1+2*ld2,1),1)
      endif

      if (ndnfld .or. .not. ccor) call defdc(oh,-l2)

C --- Noncollinear and/or S-O three-center hamiltonian ---
      if (nsp .eq. 2 .and. lnc) then
        if (lso .or. lbf) then
          call defrr(osoph,-32*ndim)
          call mksoph(nl,nsp,nbas,lihdim,ipc,indxsh,sop,w(osoph))
          call defdr(osod,ldim*12)
          call mksod(ldim,lihdim,ndim,indxsh,pph,w(osoph),w(osod))
        endif

#if NC
        lrel = igets('ctrl lrel',sctrl)
        if (lrel .eq. 2) then 
          if ( .not. bittst(lham,128)) 
     .      call rx('for fully relativistic set gamma=T') 
          call upack('pot opprel',spot,opprel,0,0,0,0)
          if (lgors('ctrl lham,3',sctrl)) then
            call hmfr2c(nbas,nl,indxsh,qss,eula,neul,ldim,lihdim,
     .        ipc,nsp,w(opprel),w(osll),w(oh),z)
          else
            call hmfr3c(nbas,nl,indxsh,qss,eula,neul,ldim,lihdim,
     .        ipc,nsp,w(opprel),pp,sop,w(osll),w(oh),z)
          endif

        else
#endif
          if (lss.and.lso) call rx('S-O coupling incompatible with SS')
          call hmltnc(0,nbas,nl,indxsh,qss,eula,neul,pph,w(osod),T,
     .      ccor,lss,lnc,lso,w(occd),vmtz,ldim,lihdim,w(osll),w(oh),z,
     .      w(owk))
C     ... Add (1+oh)+ [V(LDA+U)-diagonal part] (1+oh) to the LDA hamiltonian
          if (nlibu .gt. 0 .and. ludiag .eq. 0) then
          call defdc(oohp1,-l2)
          call hmltnc(4,nbas,nl,indxsh,qss,eula,neul,pph,w(osod),T,
     .      ccor,lss,lnc,lso,w(occd),vmtz,ldim,lihdim,w(osll),
     .      w(oohp1),z,w(owk))
          call asldau(nbas,ldim,lmaxu,isp,nsp,2,nlibu,idu,
     .      indxsh,w(oohp1),vorb,w(oh))
          endif
#if NC
        endif
#endif

C        call yprm('H after hmltnc',12,w(oh),ldimx*ldimx,ldimx,ldimx,
C     .    ldimx)

C   ... External magnetic field
        if (lbf) then
          call ztoyy(w(osll),ldimx,ldimx,ldimx,ldimx,0,1)
          call ztoyy(w(oh),ldimx,ldimx,ldimx,ldimx,0,1)
C         call zprm('h',2,w(oh),ldimx,ldimx,ldimx)
          call hmladb(ldim,lihdim,ndim,indxsh,pph,bdots,w(osoph),
     .      w(osll),w(oh))
          call ztoyy(w(oh),ldimx,ldimx,ldimx,ldimx,1,0)
        endif

C        call yprm('H after hmltnc',12,w(oh),ldimx*ldimx,ldimx,ldimx,
C     .    ldimx)
        if (lso) call rlse(osoph)

C   ... Add iwaves to hamiltonian
        if (lnc .and. idim .gt. 0) then
           rewind nfilet
           if (nevmx .gt. 0) call dpdump(w(osll),1,nfilet)
           call dpdump(w(osll),l2*2,nfilet)
           call daxpy(l2*2,1d0,w(osll),1,w(oh),1)
        endif

C --- Collinear three-center hamiltonian ---
      else
        if (nlibu .gt. 0) then
        endif
        call hmltns(0,ccor,w(occd),vmtz,ldim,lihdim,pph(1,1,isp),
     .    w(osll),w(oh),z,w(owk))
        if (bittst(lham,16)) call remhyb(ldim,nl,nbas,w(oh),z)
C   ... Add LDA+U Hamiltonian
        if (nlibu .gt. 0 .and. ludiag .eq. 0) then
c         call yprm('test H',12,w(oh),ldimx*ldimx,ldimx,ldimx,ldimx)
C         This one merely adds V(LDA+U)  to H
C          call asaddu(3,nbas,ldim,lmaxu,isp,nsp,1,nlibu,idu,
C     .      indxsh,vorb,w(oh),w)
C         This one adds (1+oh)+ V(LDA+U) (1+oh) to the LDA hamiltonian
          call defdc(oohp1,-l2)
          call hmltns(4,ccor,w(occd),vmtz,ldim,lihdim,pph(1,1,isp),
     .      w(osll),w(oohp1),w,w(owk))
          call asldau(nbas,ldim,lmaxu,isp,nsp,1,nlibu,idu,
     .      indxsh,w(oohp1),vorb,w(oh))
c         call ztoyy(w(oh),ldim,ldim,ldim,ldim,1,0)
        endif
C       call yprm('after LDA+U',12,w(oh),ldimx*ldimx,ldimx,ldimx,ldimx)
      endif

      if (iprint() .ge. 110 .or. F) then
        call yprm('3-C H',12,w(oh),ldimx*ldimx,ldimx,ldimx,ldimx)
        call yprm('3-C O',12,z,ldimx*ldimx,ldimx,ldimx,ldimx)
      endif

C --- Add sigma potential ---
C ... For now, assume sigm not spin-pol; use j instead of ldimx
      j = ldimx
      j = ldim
      if (lrsig .ge. 1) then
C       Call to bloch : 104110=perm orb, transpose, no add, c*16
        hreal = 0
        if (lrsig .ge. 4) hreal = 1
        i = 100000 + 0000 + 40*(1-hreal)
        nttabs = ntabs(nbas+1)
        mxorb = nglob('mxorb')
        call bloch(lbloch+i,qp,nl,plat,mxorb,indxsh,1,nttabs,
     .    iaxs,sigrs,ndsig,isp,nsp,j,j,0,j,0,j,0,w(osll),w,w)
      elseif (lrsig .eq. -1) then
        ifi = fopna('SIGM',-1,4)
        read(ifi)
        call dpdump(w(osll),ld2*2,ifi)
      endif
      if (lrsig .ne. 0) then
C       To guarantee Hermitian
        call dosymm(2,w(osll),j,j)
        if (iprint() .ge. 110) then
          call yprm('SX sigma',2,w(osll),j*j,j,j,j)
        endif
C   ... For now, since sigma not sp
        if (lso) then
          call pvsec1(ldim,w(osll),w(oh))
        else
          call daxpy(l2*2,-1d0,w(osll),1,w(oh),1)
        endif
        if (iprint() .ge. 110)
     .    call yprm('3-C H + sx',2,w(oh),ldimx*ldimx,ldimx,ldimx,ldimx)
      endif

C --- Return with hamiltonian in w(oc) if nevmx is -2 ---
      if (nevmx .eq. -2) then
        call rx('not ready here for nevmx=-2')
C        oc = osll
C        call dcopy(l2*2,w(oh),1,w(oc),1)
C        call rlse(owk)
        return
      endif
      if (invb) then
        wtqp = wtkp(ikp)
        mxorb = nglob('mxorb')
        call ibloch(0,qp,wtqp,plat,mxorb,indxsh,1,nsite,w(oiaxb),strx,
     .              strx,strx,nl*nl,isp,nsp,ldim,ldim,ldim,ldim,w(osll))
        call ibloch(0,qp,wtqp,plat,mxorb,indxsh,1,nsite,w(oiaxb),rhrs,
     .              rhrs,rhrs,nl*nl,isp,nsp,ldim,ldim,ldim,ldim,w(oh))
        call ibloch(0,qp,wtqp,plat,mxorb,indxsh,1,nsite,w(oiaxb),rors,
     .              rors,rors,nl*nl,isp,nsp,ldim,ldim,ldim,ldim,z)
C        call invbl(T,T,nbas,nl,nsp,ldim,plat,nsite,w(oiaxb),
C     .             indxsh,ikp,nkp,qp,w(osll),w(oh),z,strx,rhrs,rors)
        call rlse(oiaxb)
        goto 999
      endif

C --- Eigenvalues and eigenvectors of H ---
#if BLAS3
      lx = .true.
#else
      lx = .false.
#endif
      call defdr(odiawk,ldimx*5)
      if (linv .ne. 0) call redfrr(odiawk,ldimx*11)
      call dcopy(2*l2,z,1,w(osll),1)
      lov = 1
      if (bittst(lham,32)) lov = 2
#if SIGMO
      if (lrsig .ne. 0) lov = 3
#endif
      call diagno(ldimx,w(oh),w(osll),w(odiawk),lx,lov,linv,nevmx,efmax,
     .  nev,z,eb)

C --- Printout ---
      if (iprint() .ge. 30) then
        j = min(9,ldimx)
        if (iprint() .ge. 35) j = ldimx
#if LINUX_PGI
        do  18  ii = 1, 1
#else
        do  18  ii = 1, 2
#endif
        call awrit3(' SECMAT:  kpt %i of %i, k=%3:2,5;5d',
     .    ' ',80,lgunit(ii),ikp,nkp,qp)
   18   write(lgunit(ii),'(255(9f8.4:/))') (eb(i), i=1,j)
        if (iprint() .ge. 36 .and. nev .gt. 0) call awrit5(
     .    ' nev, nevmx, ldim=  %i  %i  %i  ev(nev) = %1;5d  efmax '//
     .    '= %1;5d',' ',80,i1mach(2),nev,nevmx,ldimx,eb(nev),efmax)
        call ftflsh(lgunit(1))
      endif
      if (iprint() .ge. 110 .and. lov .ne. 3) then
        outs = 'evec'
        if (lov .eq. 2) outs = 'ortho evec'
        call yprm(outs,2,z,ldimx*ldimx,ldimx,ldimx,nev)
        call yprm('eval',1,eb,ldimx*1,ldimx,nev,1)
        call query('V<110 to skip matrix printing',-1,0)
      endif

  999 continue 
      call rlse(occd)

      end
      subroutine pvsec1(ldim,sigm,h)
C     implicit none
      integer ldim,i,j
      double precision sigm(ldim,ldim,2),h(ldim,2,ldim,2,2)

      do  10  j = 1, ldim
      do  10  i = 1, ldim
        h(i,1,j,1,1) = h(i,1,j,1,1) - sigm(i,j,1)
        h(i,2,j,2,1) = h(i,2,j,2,1) - sigm(i,j,1)
        h(i,1,j,1,2) = h(i,1,j,1,2) - sigm(i,j,2)
        h(i,2,j,2,2) = h(i,2,j,2,2) - sigm(i,j,2)
   10 continue

      end
      subroutine pvsec2(ld,s)
C     implicit none
      integer ld
      double precision s(ld,2)

      if (ld .eq. 0) return
      call dswap(ld,s,1,s(1,2),1)
      end

