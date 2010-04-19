C- TBE Brillouin zone integral for total energy, forces
      subroutine tbzint(prgnam,sbz,sctrl,sham,spot,slat,smix,
     .  sspec,ssite,sstr,sarray,stb,sstrn,ltbe)
C- TBE self-consistency loop
Cu      Aug 06 (ATP) TB+U: nsp is already used for S-O (JK) so
C                    use nspc for this as these are coupled spins;
C                    then nsp=2 will denote either TB+U or TB-L
C                    with spin polarisation. We will also have
C                    nsp1=2 if nsp==2 | nspc==2 for dimensioning
Cu Updates
Cu   04 Jun 08 (ATP) Handles molecules as well as solids
Cu    8 Jun 07 (MvS) Merged Klepeis's additions to TB package
Cu   15 Feb 02 (ATP) Added MPI parallelization
C switches in ltb:
C bit  decimal  token       switch
C  0      1      OVLP       non orthogonal TB
C  1      2      CRYSF      crystal field (empirical)
C  2      4      OVCF        - ditto - with overlap
C  3      8      ADDES      add e*S to H
C  4     16      FORCES     calculate force
C  5     32      FIJ        force on atom i due to atom j
C  6     64      DONLY      dos only
C  7    128      3PV        calculate pressure
C  8    256      EVDISC     keep e'vecs on disc for BZ integration
C  9    512      PAIR       pair potential only
C 10   1024      TRH        calculate Tr[rho][H]
C 11   2048      RHO        calculate local charges
C 12   4096      not used
C 13   8192      TBU        TB+U
C 14  16384      NOUAVG     Old TB-L don't average U
C 15  32768      UL         New TB-L using delta-n and average U (PRL)
C 16  65536      IODEL      Read delta-H or delta-rho from disc
C 17 131072      GAMMA      gamma-point only
C 18 262144      MOL        do a molecule, or cluster

C     implicit none

C#ifdefC MPI
C      include "mpif.h"
C      integer numprocs, status(MPI_STATUS_SIZE)
C      integer MAX_PROCS
C      parameter (MAX_PROCS = 100)
C      integer resultlen
C      character*(MPI_MAX_PROCESSOR_NAME) name
C      character*10 shortname(0:MAX_PROCS-1)
C      character*20 ext
C      character*26 datim
C      integer namelen(0:MAX_PROCS-1)
C#endif

      character*(*) sstrn, prgnam*8
      double precision sarray(1),sbz(1),sctrl(1),slat(1),smix(1),
     .  sspec(1),ssite(1),sstr(1),sham(1),spot(1),stb(1)
      integer ltbe
C ... Heap
      integer w(1)
      common /w/ w
C ... Local variables
      character*120 outs, ch*1
      character*20 zvfrom(3)
      logical bittst,cmdopt,lgors,T,F
      parameter (T=.true., F=.false.)
      integer oeband,oics,oidmod,oinitc,oipc,olmx,onrc,onrcp,opnu,oqc,
     .  oqnu,oqt,ormax,ovdif,oves,owk,oz,oclabl,odq,oidxdn,oistab,
     .  osymgr,oiwk,odlv,oqlv,odlv2,oqlv2,ostni
      integer bitand,fopn,i,j,i1,i2,ifi,igets,ipr,iprint,i1mach,lasa,
     .  lpgf,lgunit,lncol,nkd,nkq,lstnr(3),maxit,nband,nbas,nclasp,
     .  nclspp,nclass,nsgrp,nl,nspin,nlspc,nspc,nsp,isp,nsp1,
     .  nspec,sdmod,fopna,fopnn,fxst,inxsh(4),izvfr,nkp,nwin
C     MPI
      integer procid, master, mpipid
C ... For the Gaunt coefficients
      integer indxcg(1300),jcg(6500)
      double precision cg(6500),cy(289),gaunt(9,9,25),ak(9,9,9,9,3)
      double precision dgets,d1mach
C ... Total energy and potential.  For ASA Harris functional:
      double precision thrpv,amgm,efermi(2),emad,zval
C ... Spin dynamics
      double precision sdprm(6)
C ... for PGF
      integer npadl,npadr,nbasp
      integer opgfsl
C     integer nzp,npl
C     integer ozp,owz,ogll,olgii,oogll,opgplp
C     double precision semsh(10)
C ... tight-binding specific.  Equivalence order needed for isave
      integer ltb,memx,mxnbr,oeref,oerep,oawk,oawk1,oawk2,odelta,otj,
     .        otj1,otj2,ovel,pdim,is,orhoit,orhlit,oqit,orho0,orhl0,
     .        orhlm0,ommom0,ommom,omit,oidu,ouh,ojh
      parameter (memx=500)
      double precision alat,rmaxh,emg,energy(7),sumev,entrpy,etot,efree,
     .                 alpha,erep,vol,avw,ecorr,rmax,ppdip(3),fmax,tpvq,
     .                 tpv
      equivalence (amgm,energy(1)), (sumev,energy(2)), (erep,energy(3)),
     .            (etot,energy(4)), (emad,energy(5)),
     .            (thrpv,energy(6)),(entrpy,energy(7))
C ... for relaxation or dynamics
      logical md,xyzfrz(3),defer
      integer mdsw,nf,icom,nitrlx,natrlx,nvar,itrlx0,itrlx,ierr,
     .  oindrx,ow,op,opos,oips
      character*3 ensbl
      double precision mdprm(6),tstep,temp,pext,accrlx,
     .  taup,taub,ekin,tkin,cons,time0,time,mdtime,zeta,logs,eps,
     .  veps,zacc,zsqacc,alat0
      equivalence (tstep,mdprm(2)), (temp,mdprm(3)),
     .  (taup,mdprm(4)), (mdtime,mdprm(5)), (taub,mdprm(6))
C ... For strux
      integer nsites,oalpha,oiax,ortab,orham,oh,oh0,mxcsiz,nttab,ontab,
     .  os2
      integer lmxst,nkdmx,nkqmx
      double precision ekap(20),plat(3,3),plat0(3,3),qlat(3,3),as,awald,
     .  ewtol,rpad,gam(4),gx,gy,gz,gt
      equivalence (gam(1), gx), (gam(2), gy), (gam(3), gz), (gam(4), gt)
C ... For O-N
      integer mordrn
C ... for charge mixing
      logical parmxp,uconv,lc
      double precision beta,betav(2),tjmax,ddot,cnvg,rmsdel,dum(3),dum2
      integer nelts,neltsm,neltst,mmix,nitmax,it,it0,nkill,broy
C ... for bzmaps,dos
      logical bzmp
      integer mull,npln,nwmx,nqmx,onq,onw,ovx,ovy,oxx,oyy,ozz,ozos
      double precision wtbzm,ckbas,cksumf
C ... for TBBND
      integer nfilem,ndim,ldim,idim,nlb,nlbmx,ifib,nq,iq,nevmx,nev
      parameter (nlbmx=2000)
      double precision q(3),q1(3),q2(3),evl(nlbmx),xx,efmax
      integer oidxsh,oe
      integer fopno
C ... for the tight-binding hamiltonian
      logical swtmp, rl, mixrho, mixh, mixQ, mixQ2
      logical ovlp,cryf,ocryf,addsll,force,pv,pair,rdham,
     .  MOL,UL,TBU,sc,lso,tbe,tbbnd,tbfit,iodel,ltmp,iostr
      integer n,m,nterm,memode,nset,npair,getef,nlme,nlmesp
      integer otabme,ov0,odecay,oitab,oidec,nfilin,otabcf,otabov,otbocf,
     .  odeccf,odecov,odcocf,oitbcf,oitbov,oitocf,oidcf,oidov,oidocf,
     .  odh,opot0,odhcf,odov,odovcf,of,ofe,ofnou,ofnom,oesite,ohso,
     .  orho,oindex,oemag,onpm,oiam,odel,odelL,oamass,onpr,ovso,
     .  odrhos,odrosl,orhon,orhoc,oqpol,oqmpol,ovm,ovu,ovj,odelh,
     .  ofwk,oew,oov,orholm,oqp,oiq1,oivso,oifit

      data oitab /1/ oidec /1/ odeccf /1/ otabcf /1/ otabov /1/
     .  odecov /1/ oitbcf /1/ oidcf /1/ oidocf /1/ oidov /1/ oew /1/
     .  onq /1/ onw /1/ ovx /1/ ovy /1/ oxx /1/ oyy /1/ ozz /1/
     .  ozos /1/ oindex /1/ oemag /1/
 
      data xyzfrz /F,F,F/ npln /0/ nwmx /0/ nqmx /0/
C     data autime /0.048377d0/
      data zvfrom /'ZVAL token','START category','Atomic number'/
C#ifdefC MPI
C      call MPI_COMM_RANK( MPI_COMM_WORLD, procid, ierr )
C      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
C      call MPI_GET_PROCESSOR_NAME(name, resultlen, ierr)
C      call strcop(shortname(procid),name,10,'.',i)
C      namelen(procid) = i-1
C      master = 0
C      call MPI_BARRIER( MPI_COMM_WORLD, ierr )
C#endif
      procid = mpipid(1)
      master = 0

      tbe   = ltbe .eq. 1
      tbbnd = ltbe .eq. 2
      tbfit = ltbe .eq. 4

      call upack('ctrl nbas nclass nl nspec nspin',sctrl,nbas,nclass,nl,
     .  nspec,nspin)
      call upack('array ohave oics oclabl onrc onrcp',sarray,oinitc,
     .  oics,oclabl,onrc,onrcp)
      call upack('array npadl npadr nclasp oips oipc',sarray,npadl,
     .  npadr,nclasp,oips,oipc)
      call upack('pot opnu oqnu oqc oqt oves',spot,opnu,oqnu,oqc,oqt,
     .  oves)
      call upack('ctrl maxit lasa lncol',sctrl,maxit,lasa,lncol,0,0)
      lpgf = igets('ctrl lpgf',sctrl)
      call upack('ctrl sdmod sdprm lstonr',sctrl,sdmod,sdprm,lstnr,0,0)
      call upack('ctrl tol',sctrl,cnvg,0,0,0,0)
      call upack('bz lmull',sbz,mull,0,0,0,0)

      call sp2cls('spec idmod',sspec,w(oics),nl,1,nclasp,oidmod)
      call sp2cls('spec stni',sspec,w(oics),1,1,nclasp,ostni)
      call upack('lat alat plat avw',slat,alat,plat,avw,0,0)
      call upack('lat vol oistab osymgr nsgrp opos',slat,vol,oistab,
     .           osymgr,nsgrp,opos)
      call upack('ctrl nitmv mdprm ltb',sctrl,nitrlx,mdprm,ltb,0,0)
      call upack('str mxnbr rmax',sstr,mxnbr,rmaxh,0,0,0)
      call upack('lat awald nkd nkq',slat,awald,nkd,nkq,0,0)
      call upack('lat odlv oqlv gam',slat,odlv,oqlv,gam,0,0)
      call upack('lat as tol rpad nkdmx nkqmx',slat,as,ewtol,rpad,nkdmx,
     .           nkqmx)
      call dpzero(betav,2)

C --- make larger arrays for lattice vectors for the case of NPT ---
      call defdr(odlv2,  -3*nkdmx)
      call defdr(oqlv2,  -3*nkqmx)
      call dcopy(3*nkd,w(odlv),1,w(odlv2),1)
      call dcopy(3*nkq,w(oqlv),1,w(oqlv2),1)
      
C ... switches
      ovlp   = bittst(ltb,1)
      cryf   = bittst(ltb,2)
      ocryf  = bittst(ltb,4)
      addsll = bittst(ltb,8) .and. ovlp
      force  = bittst(ltb,16)
C     stress = bittst(ltb,32)
      pv     = bittst(ltb,128)
      pair   = bittst(ltb,512)
C     trh    = bittst(ltb,2**10)
      TBU    = bittst(ltb,2**13)
      UL     = bittst(ltb,2**15)
      MOL    = bittst(ltb,2**18)
      iodel  = bittst(ltb,2**16)
      sc = UL .or. TBU
      lso    = bittst(lncol,4)
      bzmp = lgors('ctrl ldos,8',sctrl)
C --- defer force calculation until self consistent ---
      defer = T
      if (cmdopt('--nodefer',9,0,outs)) then
        defer = F
      endif

C --- read hamiltonian and overlap as strux from disc ---
      rdham = F
      if (cmdopt('--rdham',7,0,outs)) then
        rdham = T
        sc = F
      endif
      
C --- retain nspin for dimensioning certain arrays, eg qnu ---
      nsp1 = nspin
      if (TBU .or. UL) then
        if (nspin .eq. 1 .and. TBU)
     .    call rx('TBZINT: for TB+U restart with NSPIN=2')
        nsp = nspin
        nspc = 1
      else
        nspc = nspin
        nsp = 1
      endif
C ... Modes of mixing
C 1. Mix density matrix
      mixrho = F
C 2. Mix multipole moments and magnetic moments separately: beta, betav
      mixQ2 = cmdopt('--mxq2',6,0,outs)
C 3. Mix multipole moments and magnetic moments togther: beta
      mixQ = cmdopt('--mxq',5,0,outs) .or. mixQ2
C 4. Mix hamiltonian (default); cannot use with overlap
      mixh = T
      if (cmdopt('--mxr',5,0,outs)) then
        mixrho = T
        mixh = F
      endif
      lc = F
      if (pair .or. ovlp .or. mixQ .or. mixQ2) then
        call rxx(bittst(ltb,2**14),' TBZINT: need ave. U for Q mixing')
        lc = cmdopt('--lc',4,0,outs)
        mixh = F
      endif

      ckbas = cksumf(w(opos),3*nbas)
      mdsw = nint(mdprm(1))
      md = mdsw .gt. 0 .and. mdsw .lt. 4
      zval = dgets('bz zval',sbz)
      izvfr = 1
      nlspc = nl*nspc*nclass
      nbasp = nbas + npadl + npadr
      nclspp = 2*nclasp-nclass
C ... Hardwire no order-N  for now
      mordrn = 0
      call defrr(ovdif,-nclspp)
C ... eref, lmx,z by class
      call sp2cls('spec eref',sspec,w(oics),1,1,nclasp,oeref)
      call redfrr(oeref,2*nclasp)
      call dvset(w(oeref),1+nclasp,2*nclasp,0d0)
      call sp2cls('spec z',sspec,w(oics),1,1,nclasp,oz)
      call sp2cls('spec lmxa',sspec,w(oics),1,1,nclasp,olmx)
      call upack1('array ormax',sarray,ormax)
      call sp2cls('spec idxdn',sspec,w(oics),nl,1,nclasp,oidxdn)
C ... Spin-orbit coupling parameters stashed in alpha
C     call sp2cls('spec alpha',sspec,w(oics),nl,1,nclasp,ovso)
      call defdr(odel,-nl*nsp1*nbasp)
      call spackv(10,'site delta',ssite,1,nbasp,w(odel))

C ... TB polarisation parameters, potentials and Hubbard U's and J's
      if (sc) then
        call sp2cls('spec qpol',sspec,w(oics),10,1,nclasp,oqpol)
        call defdr(oqmpol,9*nbas)
        call defdr(ovm,25*nbas)
        call defdr(ovu,nl*nbas)
        call defdr(ovj,-2*nbas)
        call defdr(odelh,-nl**4*nbas*nsp)
        call defdr(odelL,-nl**4*nbas*nsp)
        if (ovlp) call defdr(opot0,-nbas)
        call defdr(ouh,-4*nbas)
        call defdr(ojh,-4*nbas)
        call defi(oidu,-4*nbas)
        if (nsp .eq. 2) then
          call getujh(sspec,nl,nbas,w(oipc),w(oclabl),
     .                w(oidu),w(ouh),w(ojh))
        endif
      endif

      if (.not. (tbe .or. tbbnd)) goto 499

C ... tbe-specific initialization
C ... for alpha, see Phys Rev B, 53, 15381 (1996)
      alpha = 0d0
      entrpy = 0d0
      tpvq = 0
      if (bzmp) then
        call rxx(bzmp,'TBZINT: not ready for bzmaps')
        call rx('bzmp branch not checked')
        ifi = fopno('BZPL')
        rewind ifi
        nkp = 0
   93   read (ifi,*,err=96,end=96) nwin,npln
        do  94  i = 1, nwin
          read (ifi,*,err=96,end=96) xx,xx
   94   continue
        do  95  i = 1, npln
          read (ifi,*,err=96,end=96) xx,xx,xx,xx,xx,n,
     .                               xx,xx,xx,xx,xx,m,xx
          nkp = nkp + n*m
   95   continue
        goto 93
   96   call rxx(nkp .eq. 0,prgnam//': empty or badly formed BZPL file')
      else
        nkp = igets('bz nkp',sbz)
      endif
      nband = nbas*nl**2
      call defdr(oeband,nband*nsp1*nkp)
      call subasi(sctrl,sspec,sham)

C ... Read positions from file
      if (cmdopt('--rpos=',7,0,outs))
     .    call iopos(F,0,outs(8:),nbasp,w(opos))

C ... Write positions to file, nit=0
      nitrlx = igets('ctrl nitmv',sctrl)

      if (nitrlx .eq. 0) then
C ... Write positions to file, nit=0
        if (cmdopt('--wpos=',7,0,outs))
     .    call iopos(T,0,outs(8:),nbasp,w(opos))
        call rx(prgnam//': no iterations')
      endif

C ... Make only DOS weights; not total energy, E_F, etc.
      if (lgors('ctrl ldos,64',sctrl)) goto 499

C --- Get zval. Priority is: 1. ZVAL from CTRL; 2. from moms; 3. from Z
      if (zval .eq. 0) then
        izvfr = 2
        call getmom(nsp1,nl,w(oqnu),nclass,w(onrc),w(oidxdn),zval)
        if (zval .eq. 0) then
          izvfr = 3
          call pshpr(0)
          call defdr(odq,nclasp)
          call getq(nsp1,nl,w(olmx),nclasp,w(oz),w(opnu),w(oqnu),
     .      w(oics),sspec,w(oqc),w(oqt),w(odq))
C          call getq(nsp1,nl,w(olmx),nclasp,w(oz),w(opnu),w(oqnu),0,0,
C     .      w(oqc),w(oqt),w(odq),0d0,0d0)
          call getzv(nclasp,w(onrcp),w(oz),w(oqc),sbz,sctrl,zval)
          call rlse(odq)
          call poppr
        endif
      endif
      if (iprint() .ge. 30) then
        call awrit1('%N TBZINT: zval=%d electrons from '//
     .              zvfrom(izvfr),' ',128,i1mach(2),zval)
      endif

C --- Set up relaxation parameters ---
      call defi(oindrx,6*nbas)
      call rlxstp(sctrl,ssite,natrlx,nvar,w(oindrx),xyzfrz,pdim)
      call rlse(oindrx)
      icom = 0
      if (nvar .ne. 0) then
        call defi(oindrx,2*natrlx)
        call defdr(ow,nvar*nvar)
        call defdr(op,pdim)
      endif

C --- Set up force/charge symmetrization ---
  499 continue

C --- Read tight-binding matrix elements ---
      if (procid .eq. master) then
      nfilin = fopna('CTRL',-1,1)
      endif
C     An upper bound to the number of coefficients to be read in
      nterm = 9
C     Note: ME are spin dependent in noncollinear case: 
C     Then groups of three blocks (++, --, +-)
      nlmesp = nlme(nl)*(2*nspc-1)
      call defdr(otabme,-nterm*nlmesp*3*memx)
      call defdr(oV0,-9)
      call defdr(odecay,-nlmesp*3*nclass**2)
C --- First pass to get dimensions ---
      call rdtbh(nfilin,nclass,w(oclabl),w(oz),nlmesp,ltb,nterm,1,F,0,
     .  0,w(otabme),w(otabcf),w(otabov),w(otbocf),w(odecay),w(odeccf),
     .  w(odecov),w(odcocf),w(oitab),w(oitbcf),w(oitbov),w(oitocf),
     .  w(oidec),w(oidcf),w(oidov),w(oidocf),w(oV0),nset,memode)
      call rlse(otabme)
      nterm = 1
      if (mod(memode,10) .eq. 4) nterm = 9
      if (mod(memode,10) .eq. 5) nterm = 5
      if (mod(memode,10) .eq. 6) nterm = 2
      if (mod(memode,10) .eq. 7) nterm = 4
      if (mod(memode,10) .eq. 6 .and. nl .ne. 2)
     .  call rx(prgnam//': need nl=2 for memode 6')
      if (tbfit) then
        if (memode .ge. 10)
     .    call rx(prgnam//': canonical ham not allowed for fit')
        if (mordrn .ne. 0)
     .    call rx(prgnam//': order N not implemented for fit')
        call rxx(lpgf .ne. 0,prgnam//': PGF not implemented for fit')
        call rxx(UL,prgnam//': UL=T not implemented for fit')
        call rxx(nitrlx.gt.1,prgnam//': cannot do relaxation with fit')
      endif
      call defdr(otabme,-nterm*nlmesp*nset)
      call defdr(odecay,-nlmesp*nset)
      if (tbfit) then
        call defi(oitab,-nterm*nlmesp*nset)
        call defi(oidec,-nlmesp*nset)
      endif
      if (cryf) then
        call defdr(otabcf,-nterm*nlmesp*nset)
        call defdr(odeccf,-nlmesp*nset)
        if (tbfit) then
          call defi(oitbcf,-nterm*nlmesp*nset)
          call defi(oidcf,-nlmesp*nset)
        endif
      endif
      if (ovlp) then
        call defdr(otabov,-nterm*nlmesp*nset)
        call defdr(odecov,-nlmesp*nset)
        if (tbfit) then
          call defi(oitbov,-nterm*nlmesp*nset)
          call defi(oidov,-nlmesp*nset)
        endif
      endif
      if (ocryf) then
        call defdr(otbocf,-nterm*nlmesp*nset)
        call defdr(odcocf,-nlmesp*nset)
        if (tbfit) then
          call defi(oitocf,-nterm*nlmesp*nset)
          call defi(oidocf,-nlmesp*nset)
        endif
      endif
      call defdr(oV0,-9*nset)
      call defi(onpm,2*nclass)
      call defi(oiam,4*3*nclass**2)
C --- Second pass to get remaining parameters and pointer arrays ---
      npair = nset
      call rdtbh(nfilin,nclass,w(oclabl),w(oz),nlmesp,ltb,nterm,2,tbfit,
     .  w(oiam),w(onpm),w(otabme),w(otabcf),w(otabov),w(otbocf),
     .  w(odecay),w(odeccf),w(odecov),w(odcocf),w(oitab),w(oitbcf),
     .  w(oitbov),w(oitocf),w(oidec),w(oidcf),w(oidov),w(oidocf),w(oV0),
     .  npair,memode)
      call rlse(oiam)
      call defi(oiam,3*npair)
      call fclose(nfilin)
C --- Printout ---
      if (iprint() .gt. 30) then
        if (tbe .or. tbbnd) then
          i = 1
          if (memode .eq. 6) i = 0
        else
          i = 0
        endif
        call shotbm(w(otabme),w(odecay),w(oV0),nset,npair,nterm,nl,nspc,
     .    nlmesp,memode,nclass,w(oclabl),w(oiam),i,'Hamiltonian')
        if (cryf .and. memode .ne. 1) call shotbm(w(otabcf),
     .    w(odeccf),w(oV0),nset,npair,nterm,nl,nspc,nlmesp,memode,
     .    nclass,w(oclabl),w(oiam),0,'Hamiltonian crystal field')
        if (ovlp .and. memode .ne. 1) call shotbm(w(otabov),
     .    w(odecov),w(oV0),nset,npair,nterm,nl,nspc,nlmesp,memode,
     .    nclass,w(oclabl),w(oiam),0,'overlap')
        if (ocryf .and. memode .ne. 1) call shotbm(w(otbocf),
     .    w(odcocf),w(oV0),nset,npair,nterm,nl,nspc,nlmesp,memode,
     .    nclass,w(oclabl),w(oiam),0,'overlap crystal field')
      endif

C --- Set up for spin-orbit matrix elements ---
      if (lso) then
C       call sp2cls('spec alpha',sspec,w(oics),nl,1,nclasp,ovso)
        call sp2cls('spec vso',sspec,w(oics),nl,1,nclasp,ovso)
C       call shstru('spec',sspec,1,nspec)
C       call prmx('vso',w(ovso),nl,nl,nclasp)
        call defdc(ohso, -4*nl**4)
        call skhso(nl,w(ohso))
      else
        call defrr(ovso,1)
      endif

C --- Make mxnbr ---
      if (rmaxh .gt. d1mach(3)) then
        if (mxnbr .eq. 0) then
          mxnbr = 4*nint((rmaxh*alat/avw)**3*nbasp)
        else
          mxnbr = mxnbr*nbas
        endif
      endif

C --- Get LL' on-site hamiltonian melts from disk
      if (iodel .and. sc) then
C#ifdefC MPI
C        if (procid .eq. master) then
C          call iodelL(ltb,1,nl,nbas,nsp,w(odel),w(odelL))
C        endif
C        call MPI_BCAST(w(odelL),nl**4*nsp*nbas,MPI_DOUBLE_PRECISION,
C     .                 master,MPI_COMM_WORLD,ierr)
C#else
        call iodelL(ltb,1,nl,nbas,nsp,w(odel),w(odelL))
C#endif
      endif

      nelts = nl*nspc
      if (nspc .eq. 2 .and. nsp .eq. 2)
     .    call rx(' UL=T not set up for spin orbit coupling yet')

C --- Set up for tight-binding MD ---
      if (nvar .eq. 0 .and. .not. md) nitrlx = 1
      itrlx0 = 0
      if (md) then
C#ifdefC MPI
C          if (procid .eq. master) ifi = fopn('STRT')
C#else
        ifi = fopn('STRT')
C#endif
C --- choose MD ensemble ---
        ensbl = 'NVE'
        if (mdsw .eq. 2) then
          ensbl = 'NVT'
        endif
        if (mdsw .eq. 3) then
          call rxx(.not.pv,'TBZINT set 3PV=T for NPT')
          ensbl = 'NPT'
          call dcopy(9,plat,1,plat0,1)
          alat0 = alat
        endif
        nitrlx = mdtime / tstep + 1
        time = 0d0
C --- hardwire external pressure to zero for now ---
        pext = 0d0
        eps = 0d0
        call defdr(ovel,3*nbas)
        call defdr(oamass,-nspec)
        call spackv(10,'spec mass',sspec,1,nspec,w(oamass))
        if (ddot(nspec,w(oamass),1,w(oamass),1) .eq. 0)
     .    call rx('TBZINT: masses must be all nonzero')
        if (.not. cmdopt('--st',4,0,outs)) then
C#ifdefC MPI
C            if (procid .eq. master) then
C              call iostrt(nbas,nf,itrlx0,w(opos),w(ovel),zeta,time0,
C     .                    ifi,ierr)
C            endif
C            if (ierr .eq. 0) then
C              call MPI_BCAST(itrlx0,1,MPI_INTEGER,master,
C     .                       MPI_COMM_WORLD,ierr)
C              call MPI_BCAST(zeta,1,MPI_DOUBLE_PRECISION,master,
C     .                       MPI_COMM_WORLD,ierr)
C              call MPI_BCAST(w(opos),3*nbas,MPI_DOUBLE_PRECISION,master,
C     .                       MPI_COMM_WORLD,ierr)
C              call MPI_BCAST(w(ovel),3*nbas,MPI_DOUBLE_PRECISION,master,
C     .                       MPI_COMM_WORLD,ierr)
C            endif
C#else
          call iostrt(nbas,nf,itrlx0,w(opos),w(ovel),eps,zeta,zacc,
     .                zsqacc,veps,time0,ifi,ierr)
C#endif
          time = time0
          if (ensbl .eq. 'NPT') then
C --- remake vectors in scaled alat ---
            alat = alat0 * exp(eps)
            lmxst = 6
            call pshpr(0)
            call lattc(as,ewtol,rpad,alat,alat,plat,gx,gy,gz,gt,
     .                plat,qlat,lmxst,vol,awald,w(odlv2),nkd,w(oqlv2),
     .                nkq,nkdmx,nkqmx)
            call poppr 
          endif
          call zercmv(nbas,w(ovel),w(oamass),w(oips))
        endif
        if (cmdopt('--st',4,0,outs) .or. ierr .eq. 1) then
          if (ierr .eq. 1 .and. iprint() .ge. 10) then
            print *,' TBZINT *warning* could not read strt file '//
     .              'starting new MD'
          endif
          call initv(ensbl,nbas,nf,tstep,w(oips),w(oamass),temp,pext,
     .               taup,taub,nspec,1,' ',w(oclabl),zeta,zacc,zsqacc,
     .               veps,w(ovel))
        endif
        logs = zacc * tstep
        mdtime = time0 + mdtime
      endif
C --- Make Gaunt (CG) coefficients ---
      if (sc) then
        call sylmnc(cy,16)
        call scg(6,cg,indxcg,jcg)
        call makcg9(indxcg,jcg,cg,gaunt,ak)
      endif

C --- set up arrays for forces, density matrix ---
      call defi(oiwk,nbas)
      call defdr(oerep,2*nbas)
      call defdr(oesite,-nbas*nsp1)
      call defdr(ofnou,-4*nbas)
      call defdr(ofnom,-4*nbas)
      call defdr(of,-3*nbas)
      call defdr(ofe,-3*nbas)
      call defdr(orho,-6*nbas)
      call defdr(ommom,-nbas)
      call defdr(orholm,-nl**2*2*nbas)
      if (sc) then
        call defdr(orhoc,-nl**4*nbas)
        if (ovlp) then
          call defi(oidxsh,nbas*nl**2)
          call pshprt(0)
          call makidx(nl,1,1,nbas,0,sspec,w(oips),-1,w(oidxsh),inxsh)
          call popprt
          ldim = inxsh(1)
          call rlse(oidxsh)
          call defdr(odrosl,-3*ldim**2)
        endif
      endif
      if (TBU) call defdr(orhon,-nl**2*nl**2*nsp*nbas)
C .. set up work arrays for mixing
      if (sc) then
        call upacks('strn mix',i1,i2)
        tjmax = 10d0
        broy = 0
        dum(1) = 1
        dum(2) = 1
        dum2   = 0
        rmsdel = 1d0
        nkill = 0
        call pshpr(80)
        if (.not. parmxp(1,sstrn(i1:i2),i2-i1+1,broy,mmix,dum,beta,
     .      dum2,outs,dum2,nkill,betav,rmsdel))
     .      call rx(prgnam//': parse in parmxp failed')
        call poppr
C .. input Mulliken charges (l and lm resolved) ..
        call defdr(ommom0,-nbas)
        call defdr(orhl0,-nl*2*nbas)
        call defdr(orhlm0,-nl**2*2*nbas)
        neltst = nl**4*nsp*nbas
        if (mixrho) then
          if (ovlp .or. nsp .eq. 2) then
            call defdr(orhlit,-nl*nsp*nbas*(mmix+2)*2)
            neltst = neltst + nl*nsp*nbas
          endif
          call defdr(orho0,-nl**4*nsp*nbas)
          call defdr(orhoit,-nl**4*nsp*nbas*(mmix+2)*2)
          call mkrho0(1,nl,nbas,nclass,w(oipc),w(oclabl),nsp,
     .                w(oqnu),w(orho0),w(orhl0),w(orhlm0),w(ommom0))
        elseif (mixh) then
          call defdr(odelta,-nl**4*nsp*nbas*(mmix+2)*2)
          if (ovlp) then
            neltst = neltst + nbas*(nbas-1)/2
          endif
        elseif (mixQ) then
          if (pair) then
            call dcopy(6*nbas,0d0,0,w(orho),1)
          else
            if (mixQ2) then
              nelts  = 9 * nbas
              neltsm = nbas
              neltst = (9 + (nsp-1))*nbas
              call defdr(oqit,-nelts*(mmix+2)*2)
              omit = 1
              if (nsp .eq. 2) then
                call defdr(omit,-neltsm*(mmix+2)*2)
              endif
              call defdr(oawk1,-nelts*(mmix+2)*2)
              call defdr(oawk2,-neltsm*(mmix+2)*2)
              call defdr(otj1,mmix)
              call defdr(otj2,mmix)
            else
              neltst = (9 + (nsp-1))*nbas
              call defdr(oqit,-neltst*(mmix+2)*2)
              call defdr(omit,-nbas*(mmix+2)*2)
              call defdr(oawk,-neltst*(mmix+2)*2)
              call defdr(otj,mmix)
            endif
            call mkrho0(0,nl,nbas,nclass,w(oipc),w(oclabl),nsp,
     .                  w(oqnu),w,w(orho),w(orhlm0),w(ommom0))
            if (lc) nelts = nbas
          endif
        else
          call rx(' TBZINT: No mixing option set')
        endif
        if (.not. mixQ) then
          call defdr(oawk,-neltst*(mmix+2)*2)
          call defdr(otj,mmix)
        endif
      endif

C ----------------- Relaxation or MD loop ---------------------
      do  2  itrlx = itrlx0+1, itrlx0+nitrlx

        if (.not. iodel .and. sc) then
          call dcopy(nl**4*nbas*nsp,0d0,0,w(odelL),1)
        endif

C --- Make tight-binding Hamiltonian ---
        call defi(onpr,2*nbas)
C --- Embedded cluster scheme (order N) ---
        if (mordrn .eq. 1) then
          call rx('no order-N scheme implemented')
        endif

C   --- Get neighbor table iax for each atom in the cluster ---
        if (lpgf .ne. 0) then
          i = 2
          j = 1
        else
          i = 3
          j = -1
        endif
C   ... Make nttab,ontab,oiax
        call defdr(opgfsl, nbas)

        if (rdham) then
          ckbas = cksumf(w(opos),3*nbas)
          call pshpr(100)
          ltmp = iostr(8,'HAM',nl,nbas,1,ekap,0,ckbas,-1,nttab,oalpha,
     .                 oiax,ontab,oh)
          oh0 = oh
          if (ovlp) then
            ltmp = iostr(8,'OVL',nl,nbas,1,ekap,0,ckbas,-1,nttab,
     .                   oalpha,oiax,ontab,oov)
          endif
          call poppr
        else
          if (rmaxh .gt. d1mach(3)) then
            mxcsiz = mxnbr
            if (MOL) then
              call defdr(ontab, nbas+1)
              call defdr(oiax,  300*nbas*10)
              call defdr(ortab, 300*nbas*3)
              call defdr(orham, nspec)
              call dcopy(nspec,alat*rmaxh/2d0,0,w(orham),1)
              call hpair(nbas,w(oips),alat,plat,w(opos),w(orham),nttab,
     .                   w(ontab),w(oiax),w(ortab),mxcsiz,1)
              call rlse(oiax)
              call defrr(oiax, nttab*10)
            else
              call pairs(nbas,nbasp,1d0,plat,rmaxh/2d0,w(opos),
     .                   -1,i,j,w(opgfsl),nttab,ontab,oiax,mxcsiz)
            endif
          else
            if (MOL) then
              call rx(' set RMAXH > 0 in CTRL')
            endif
            mxcsiz = 0
            call defrr(ormax, nspec)
            do  is = 1, nspec
              call upack('spec rham',sspec,is,rmax,0,0,0)
              call dvset(w(ormax),is,is,rmax/2d0)
            enddo
            call pairs(nbas,nbasp,1d0,plat,w(ormax),w(opos),
     .                  w(oips),i,j,w(opgfsl),nttab,ontab,oiax,mxcsiz)
          endif
        endif
C   ... Patch iax(6) for the padded basis
        if (nbasp .gt. nbas)
     .    call pairp6(nbas,npadl,npadr,w(oiax),w(ontab))
        if (mordrn .eq. 1) then
          call rx('pairec commented out')
C          mxnbr = 2*rmaxs**3*nbasp
C          call redfi(oiax, niax*mxnbr)
C          call pairec(nbas,nbase,0,w(oipc),alat/alat,plate,w(obas),
C     .      w(orham),nttab,w(ontab),w(oiax),w(owk),mxcsiz)
C          call redfi(oiax, niax*nttab)
        endif
        nsites = nttab
        call mkiaxd(nsites,w(olmx),w(oips),w(oiax))

C   ... tb programs not converted to ntab from npr yet ...
        call npr2tb(1,nbasp,w(onpr),w(ontab))

C   --- Slater-Koster Hamiltonian ---
        if (.not. rdham) then
          call tbham(nsp,nspc,sctrl,sarray,slat,F,0,w,w,nlmesp,
     .         mod(memode,10),w(odecay),w(odeccf),w(odecov),w(odcocf),
     .         w(oiam),w(onpm),nterm,nset,w(otabme),w(otabcf),w(otabov),
     .         w(otbocf),nsites,w(onpr),oiax,oh,oh0,odh,oov,odov,odhcf,
     .         odovcf)
C#ifdefC TBPG
C        call swapad(nl**2,nbas,npadl,npadr,w(ontab),w(oiax),w(oh))
C#endif
          call tbdiag(ltb,nbasp,nl,nsp,nspc,nsp1,w(oipc),nsites,
     .         w(onpr),w(oinitc),w(oqnu),nl*nsp1,w(odel),w(odelL),w(oh),
     .         w(oh0),w(oov))
        endif

C   --- Add to Hamiltonian MEs: hLL' -> hLL' + ebarLL' * sLL' ---
        if (addsll) then
          call defdr(owk,nl*nbas)
          if (tbe) then
            call pshpr(iprint())
          else
            call pshpr(0)
          endif
          call addes(nsites,nl,nspc,nbas,ltb,w(onpr),w(oiax),w(oov),
     .      w(odov),w(owk),w(oh),w(odh))
          call poppr
          call rlse(owk)
        endif

C   --- Write Hamiltonian to STR file ---
        if (lpgf .ne. 0 .or. cmdopt('-dumph',6,0,outs)) then
          call npr2tb(0,nbasp,w(onpr),w(ontab))
          call defdr(oalpha,-nl**2*nbasp)
          ltmp = iostr(1,'STR',nl,nbasp,1,ekap,0,ckbas,-1,nsites,
     .      oalpha,oiax,onpr,oh)
          nsites = w(onpr+nbasp)
          call fclose(fopn('STR'))
          call npr2tb(1,nbasp,w(onpr),w(ontab))
          call rlse(oalpha)
        endif

C   --- Printout ---
        if (iprint() .ge. 50 .and. .not. pair) then
          call npr2tb(0,nbasp,w(onpr),w(ontab))
          call defdr(oalpha,-nbasp*nl**2)
C     ... Printout Hamiltonian
          print *
          print *,'Real-space TB hamiltonian :'

          if (TBU) then
            do j = 1, 2
              print *
              if (j .eq. 1) print *,'TB+U spin up'
              if (j .eq. 2) print *,'TB+U spin down'
              os2 = oh + (j-1)*nl**4*nsites*i1mach(18)
              call shostr(nl**2,nsites,nbasp,plat,w(opos),0,
     .          w(oalpha),w(oiax),w(ontab),w(os2),1,1,0d0,0,1d0)
C     ...  Printout Overlap
              if (ovlp .and. j .eq. 1) then
                print *
                print *
                print *, 'Real-space Overlap matrix'
                os2 = oov + (j-1)*nl**4*nsites*i1mach(18)
                call shostr(nl**2,nsites,nbasp,plat,w(opos),0,
     .            w(oalpha),w(oiax),w(ontab),w(os2),1,1,0d0,0,1d0)
              endif
            enddo
          else
            do  j = 1, nspc**2
              if (nspc .eq. 2) print *
              if (j .eq. 1 .and. nspc .eq. 2)
     .          print *,'Spin: Up-Up, w/o spin-orbit'
              if (j .eq. 2) print *,'Spin: Up-Down, w/o spin-orbit'
              if (j .eq. 3) print *,'Spin: Down-Up, w/o spin-orbit'
              if (j .eq. 4) print *,'Spin: Down-Down, w/o spin-orbit'
              os2 = oh + (j-1)*nl**4*nsites*i1mach(18)
              call shostr(nl**2,nsites,nbasp,plat,w(opos),0,
     .          w(oalpha),w(oiax),w(ontab),w(os2),1,1,0d0,0,1d0)
C     ...  Printout Overlap
              if (ovlp) then
                print *, 'overlap matrix'
                os2 = oov + (j-1)*nl**4*nsites*i1mach(18)
                call shostr(nl**2,nsites,nbasp,plat,w(opos),0,
     .            w(oalpha),w(oiax),w(ontab),w(os2),1,1,0d0,0,1d0)
              endif
            enddo
          endif
          if (iprint() .gt. 50 .and. bitand(ltb,16+128) .ne. 0) then
C       ... Printout Hamiltonian derivatives
            print *
            print *, 'Real-space gradH :'
            do  102  i = 1, 4
            print *
            if (i .eq. 1) print *, 'x-component:'
            if (i .eq. 2) print *, 'y-component:'
            if (i .eq. 3) print *, 'z-component:'
            if (i .eq. 4) print *, 'radial component:'
            do  102  j = 1, nspc**2
              if (nspc .eq. 2) print *
              if (j .eq. 1 .and. nspc .eq. 2) print *,'Spin: Up-Up'
              if (j .eq. 2) print *,'Spin: Up-Down'
              if (j .eq. 3) print *,'Spin: Down-Up'
              if (j .eq. 4) print *,'Spin: Down-Down'
              os2 = odh
     .          + (nspc**2*(i-1)+(j-1))*nl**4*nsites*i1mach(18)
              call shostr(nl**2,nsites,nbasp,plat,w(opos),0,
     .          w(oalpha),w(oiax),w(ontab),w(os2),1,1,0d0,0,1d0)
C         ... Printout Crystal Field derivatives
              if (cryf) then
                print *,'Crystal Field gradH'
                os2 = odhcf
     .            + (nspc**2*(i-1)+(j-1))*nl**4*nsites*i1mach(18)
                call shostr(nl**2,nsites,nbasp,plat,w(opos),0,
     .            w(oalpha),w(oiax),w(ontab),w(os2),1,1,0d0,0,1d0)
              endif
C         ... Printout Overlap derivatives
              if (ovlp) then
                print *,'gradO'
                os2 = odov
     .            + (nspc**2*(i-1)+(j-1))*nl**4*nsites*i1mach(18)
                call shostr(nl**2,nsites,nbasp,plat,w(opos),0,
     .            w(oalpha),w(oiax),w(ontab),w(os2),1,1,0d0,0,1d0)
C           ... Printout Overlap Crystal Field derivatives
                if (ocryf) then
                  print *,'Crystal Field gradO'
                  os2 = odovcf
     .              + (nspc**2*(i-1)+(j-1))*nl**4*nsites*i1mach(18)
                  call shostr(nl**2,nsites,nbas,plat,w(opos),0,
     .              w(oalpha),w(oiax),w(ontab),w(os2),1,1,0d0,0,1d0)
                endif
              endif
  102       continue
          endif
          call rlse(oalpha)
          call npr2tb(1,nbasp,w(onpr),w(ontab))
        endif
        if (.not. sc .and. (tbbnd .or. tbfit)) goto 3

C#ifdefC TBPG
C          if (lpgf .ne. 0) then
CC       ... Setup energy mesh ...
C            nzp = nint(semsh(1))
C            call defdc(ozp,nzp)
C            call defdc(owz,nzp)
C            call emesh(semsh,w(ozp),w(owz))
CC     ... ASA planar GF ...
C            call defi(oogll,npl+2)
C            call defi(olgii,npl+2)
C             call defdr(opp,-6*nlspc)
C          call rx('patch call to pgfasa')
C             call pgfasa(plat,w(opos),alat,nl,nbas,nspc,nclass,w(oiclsp)
C     .         ,w(onrc),clabl,wsr,rmaxs,w(olmx),nband,switch,w(opp)
C     .         ,w(ovshfp),w(oqnu),totmom,qspirl,w(oeula),nkxyz(1)
C     .         ,nkxyz(2),nl,norder,width,range,drange,npts,nkp,w(oqp)
C     .         ,w(owght),semsh,w(ozp),w(owz),w(oeband),w(oidxdn),efermi
C     .         ,vmtz,elin,w(oz),w(oqc),ntet,w(ontet),sumev,w(opgslp)
C     .         ,w(opgplp),gfopts,.true.,w(olgii),w(oogll))
C            call rx0(prgnam//' done')
C          else
C            call rx('tbpg: pgf not set')
C          endif
C#endif
        if (lgors('ctrl quit,2',sctrl)) call rx0(prgnam//
     .    ': Q=ATOM encountered')
C#ifdefC STONER
C        if (lstnr(1).ne.0) then
C          if (nsp .eq. 2) call rx(prgnam//': for STONER=T set NSPIN=1')
C          call defdr(ozos,iabs(npts)*nlspc)
C          call defi(oindex,nclass)
C          call defdr(onbar,mnpts)
C          call defdr(oewk,mnpts)
C          call defdr(omwk,mnpts)
C          call defdr(ommom,-nbas)
C          call defdr(oemag,-nbas)
C        endif
C#endif
        sumev = 0d0
        thrpv = 0d0
        emad = 0d0
        emg = 0d0
C --- Vanderbilt order N ---
        if (mordrn .eq. 2) then
          if (bitand(ltb,1+2+2**7+2**13).ne. 0 .or. nspc.eq.2)
     .      call fexit2(-1,111,' Exit -1 : VDB not implemented for'//
     .      ' ltb=%i, nspc=%i',ltb,nspc)
          if (lstnr(1).ne.0 .or. bitand(ltb,1024+2048+4096).ne.0)
     .      call fexit2(-1,111,' Exit -1 : VDB incompatible with'//
     .      ' ltb=%i, stoner=%l',ltb,lstnr)
C#ifdefC VDB
C          rmx(3) = rmaxr + rmaxh
C          call defi(onp3,-3*nbas)
C          call defi(otnp3,-3*nbas)
C          call defi(oinorb,-nbas)
C          stop 'convert iclass to ipc, vdbptr'
C          call vdbptr(nbas,nspc,w(oipc),nclass,nl,nhs,nrs,
C     .      w(oidxdn),w(oinorb),alat,plat,w(opos),
C     .      rmx,w(onpr),w(onp3),w(otnp3),tnpr,
C     .      oiax3,oihx,oirx,oirhx,ohs,odhs,odms,odm)
C          stop 'convert iclass to ipc, shorth,shordh,shortr'
C          call shorth(nbas,tnpr,w(otnp3),w(onp3),
C     .      w(oiax3),nl,w(oipc),w(oinorb),
C     .      w(oh),oh,ohs)
C          call shordh(nbas,tnpr,w(otnp3),w(onp3),
C     .      w(oiax3),nl,w(oipc),w(oinorb),nhs,
C     .      w(odh),odh,odhs)
C          call shortr(nbas,tnpr,w(otnp3),w(onp3),
C     .      w(oiax3),nl,w(oipc),w(oinorb),
C     .      w(odm),odm,odms)
CC ... Can't release oh here anymore ...
CC          call rlse(oh)
C          call cgmin(nbas,nl,tnpr,nrs,
C     .             oinorb,onp3,otnp3,
C     .             oiax3,oihx,oirx,oirhx,ohs,odms,cgmtol,sumev)
C          call hfvdb(nbas,nl,tnpr,nhs,
C     .      oinorb,onp3,otnp3,
C     .      oiax3,oihx,oirx,odms,odhs,of)
C#else
          call rxx(.true.,prgnam//': Vanderbilt not implemented')
C#endif VDB
        endif
C --- Self-consistency loop ---
        nitmax = 1
        it0 = 1
        uconv = .not. sc
        if (sc) then
          nitmax = maxit
          if (pair) nitmax = 1
          if (nitmax .eq. 1) uconv = T
          if (mixQ) it0 = 0
        endif
  700   do  1  it = it0, nitmax
C --- device to supress force calculation until self consistent ---
          if (defer .and. sc .and. .not. uconv) then
            if (iprint() .ge. 20) print *,
     .          'TBZINT: defer force calculation until self consistent'
            pv = F
            force = F
          endif
          if ((.not. sc .and. pair) .or. mordrn .gt. 0) then
            entrpy = 0d0
            call dcopy(3*nbas,0d0,0,w(of),1)
            goto 500
          endif
          if (it .gt. 0) then
          call bndtb(sctrl,sarray,sbz,slat,sspec,nbas,nl,nspc,nsp,nsp1,
     .      w(olmx),w(oidxdn),nclass,w(oips),w(oipc),w(onrc),pv,force,
     .      nband,zval,mull,npln,nwmx,nqmx,w(onw),w(oew),w(onq),w(ovx),
     .      w(ovy),w(oxx),w(oyy),w(ozz),wtbzm,nsites,w(oiax),w(onpr),
     .      xyzfrz,w(oh),w(oh0),w(odh),w(odhcf),w(ovso),w(ohso),w(oov),
     .      w(odov),w(odovcf),w(opot0),w(oiwk),w(oeband),efermi,sumev,
     .      entrpy,w(of),thrpv,w(oesite),w(orho),w(orholm),w(orhoc),
     .      w(orhon),w(ozos),w(oindex),w(odrosl),ldim)
          endif
          if (bittst(ltb,2**6)) goto 3
C     ... Electrostatic terms 
          if (sc) then
            if (it .gt. 1) then
              call upacks('strn mix',i1,i2)
              broy = 0
              call pshpr(80)
              if (.not. parmxp(it,sstrn(i1:i2),i2-i1+1,broy,mmix,dum,
     .          beta,dum,outs,dum,nkill,betav,rmsdel))
     .          call rx(prgnam//': parse in parmxp failed')
              call poppr
            endif
            if (nkill .eq. -999) then
              rmsdel = 1d-16
              goto 500
            endif
            if (nitmax .gt. 1 .and. mixrho) then
C         ... mix density matrix (and Mulliken charges) now
              call rhomix(ovlp,nl,nsp,nbas,w(oclabl),w(oipc),it,nitmax,
     .                    cnvg,mmix,neltst,beta,tjmax,w(otj),w(orhl0),
     .                    w(orho),w(orho0),w(orhoc),w(orhlit),w(orhoit),
     .                    w(oawk),rmsdel)
            endif
            call tbmpol(nbas,nsp,nl,w(oqnu),w(oipc),ltb,gaunt,
     .                  w(oqpol),w(orho),w(orhoc),w(oqmpol),w(ommom))
            if (nitmax .gt. 1 .and. mixQ) then
C         ... mix multipole moments now
              if (mixQ2) then
                call qmix2(lc,nbas,nsp,ltb,w(oclabl),w(oipc),it,nitmax,
     .                    cnvg,mmix,nkill,nelts,beta,betav,tjmax,
     .                    w(otj1),w(otj2),w(oqmpol),w(oqit),w(ommom),
     .                    w(ommom0),w(omit),w(oawk1),w(oawk2),rmsdel)
              else
                call qmix(lc,nbas,nsp,ltb,w(oclabl),w(oipc),it,nitmax,
     .                    cnvg,mmix,nkill,neltst,beta,tjmax,w(otj),
     .                    w(oqmpol),w(oqit),w(ommom),w(ommom0),w(omit),
     .                    w(oawk),rmsdel)
              endif
            endif
            call defdr(ofwk,3*nbas)
            call defdr(odrhos,nbas*nbas*3)
            call tbesel(mixQ,ssite,ltb,nbas,nl,nsp,nclass,w(oipc),
     .        w(oclabl),w(oidxdn),w(oz),w(osymgr),nsgrp,w(oistab),
     .        w(opos),awald,alat,vol,w(odlv2),nkd,w(oqlv2),nkq,
     .        w(oqpol),indxcg,jcg,cg,cy,gaunt,ak,w(orho),w(ommom),
     .        w(odrosl),ldim,w(odrhos),w(orhoc),w(orhon),w(oqnu),
     .        w(ostni),w(oidu),w(ouh),w(ojh),w(oqmpol),w(ovm),w(ovu),
     .        w(ovj),w(odelh),w(opot0),ecorr,w(ofwk),w(ofe),w(ofnou),
     .        w(ofnom),ppdip,tpvq)
            call rlse(ofwk)
            if (pair) then
              entrpy = 0d0
              call dcopy(3*nbas,0d0,0,w(of),1)
              goto 500
            endif
            if (nitmax .gt. 1) then
              if (mixrho .or. mixQ) then
                call tbadh2(nl,nbas,w(opos),plat,nsp,nsites,
     .                      w(oipc),w(oiax),w(onpr),w(odelh),w(ovm),
     .                      w(ovu),w(oh0),w(oh),w(oov))
              else
C           ... mix hamiltonian matrix elements now
                call dhmix(neltst,nl,nsp,nbas,w(oidxdn),w(oclabl),
     .                    w(oipc),it,nitmax,cnvg,mmix,nkill,beta,tjmax,
     .                    w(otj),w(odelh),w(odelL),w(odelta),w(oawk),
     .                    rmsdel)
                call tbadh1(nl,nbas,nsp,nsites,w(oiax),it,mmix,
     .                      w(onpr),w(odelta),w(oh0),w(oh))
              endif
            endif
          endif
C --- Generalised Stoner model ---
C#ifdefC STONER
C          if ((lgors('ctrl lstonr,1',sctrl))) then
C            call dpzero(w(ommom),nbas)
C            call dpzero(w(oemag),nbas)
c             stop 'iclass in stoner'
C            call stoner(nl,nbas,nclass,clabl,w(oipc),w(oidxdn),efermi,
C     .        drange,iabs(npts),w(ozos),w(oindex),w(ostni),w(oammx),
C     .        mnpts,switch(44),w(onbar),w(oewk),w(omwk),w(ommom),
C     .        w(oemag))
C          endif
C#else
      if (lgors('ctrl lstonr,1',sctrl))
     .      call rx('LM: recompile with ccomp -dSTONER ...')
C#endif

  500     continue
C --- TB total energy ---
          if (it .gt. 0 .or. pair) then
          call tbtote(sctrl,slat,sarray,ssite,nl,nsp1,w(oidxdn),nclass,
     .      w(onrc),w(oclabl),pv,force,nbas,w(oiax),w(onpr),w(oiam),
     .      w(onpm),w(oV0),w(oqnu),sumev,alpha,entrpy,emad,ecorr,
     .      w(orho),w(ostni),w(ouh),w(ojh),w(ommom),w(oemag),w(oesite),
     .      w(oerep),ppdip,thrpv,tpvq,w(of),w(ofe),fmax,w(ofnou),
     .      w(ofnom),w(oeref),erep,etot,efree,emg,amgm)
          endif
          if (pair) goto 600
          if (it .gt. 1 .and. sc) then
            if (rmsdel .lt. cnvg .or. it .eq. maxit) then
C#ifdefC MPI
C            if (procid .eq. master) then
C#endif
              call iodelL(ltb,0,nl,nbas,nsp,w(odel),w(odelL))
C#ifdefC MPI
C            endif
C#endif
              if (rmsdel .gt. cnvg .and. .not. uconv) then
                if (iprint() .gt. 10) call awrit3(
     .            ' (warning) cnvg=%g, rms=%g, after %i iterations.',' '
     .            ,80,i1mach(2),cnvg,rmsdel,it)
              endif
              uconv = .true.
C --- device to supress force calculation until self consistent ---
              if (defer) then
                if ((bittst(ltb,2**7) .and. .not. pv) .or.
     .              (bittst(ltb,2**4) .and. .not. force)) then
                  it0 = 1
                  nitmax = 1
                  pv     = bittst(ltb,2**7)
                  force  = bittst(ltb,2**4)
                  goto 700
                else
                  goto 600
                endif
              else
                goto 600
              endif
            endif
          endif
C   --- end self consistency loop ---
    1   continue
  600   continue
        if (tbbnd) goto 3

C   --- Molecular statics ---
        if (nvar .ne. 0) then
          call query('step ',4,tstep)
          call query('accrlx ',4,accrlx)
          call relax(prgnam,sctrl,ssite,sspec,itrlx,w(oindrx),natrlx,
     .               nvar,w(of),w(op),w(ow),nelts,w(odel),w(opos),icom)
          ifi = fopn('MV')
          if (cmdopt('--st',4,0,outs)
     .            .and. (md .and. itrlx .eq. itrlx0+1) ) then
            if (iprint() .gt. 10) print*, 'Starting new mv file ..'
          else
            if (iprint() .gt. 10) print*, 'Appending to mv file ..'
            call poseof(ifi)
          endif
          call awrit2('frame etot= %;4d Ry fmax=%;4da.u.',' ',128,ifi,
     .                etot,fmax)
          call iomv(nbas,w(opos),alat,ifi)
          call fclose(ifi)
C   --- Molecular dynamics ---
        elseif (md) then
          call verlet(ensbl,itrlx,nbas,nf,w(oamass),tstep,w(of),w(oips),
     .                w(oclabl),temp,efree,thrpv,vol,pext,sc,w(oqmpol),
     .                taup,taub,mdtime,time,alat,w(opos),w(ovel),zeta,
     .                eps,veps,logs,zacc,zsqacc,ekin,tkin,cons)
          if (ensbl .eq. 'NPT') then
C            call dscal(9,exp(eps),plat0,1)
C            call dcopy(9,plat0,1,plat,1)
            alat = alat0 * exp(eps)
            lmxst = 6
            call pshpr(0)
            call lattc(as,ewtol,rpad,alat,alat,plat,gx,gy,gz,gt,
     .                 plat,qlat,lmxst,vol,awald,w(odlv2),nkd,w(oqlv2),
     .                 nkq,nkdmx,nkqmx)
            call poppr
            if (iprint() .ge. 10) then
              call awrit3(' TBZINT scaling, eps=%d, alat=%d vol=%d',' ',
     .                    128,i1mach(2),eps,alat,vol)
            endif
            call pack3('lat alat plat vol',slat,alat,plat,vol)
          endif
C#ifdefC MPI
C            if (procid .eq. master) then
C#endif
          ifi = fopn('STRT')
          call iostrt(nbas,nf,itrlx,w(opos),w(ovel),eps,zeta,zacc,
     .                zsqacc,veps,time,-ifi,ierr)
          call fclose(ifi)
C#ifdefC MPI
C            endif
C#endif
        endif

C   ... Write positions to file
        if (cmdopt('--wpos=',7,0,outs))
     .    call iopos(T,0,outs(8:),nbasp,w(opos))

C --- Write variables to SAVE file ---
C#ifdefC MPI
C        if (procid .eq. master) then
C#endif
        ifi = fopn('SAVE')
        print *
        outs = 'mmom,sev,erep,etot,emad,3pv,TS'
        ch = 't'
        if (sc .and. uconv)
     .                   ch = 'c'
        if (icom .eq. 1) ch = 'r'
        if (md)          ch = 'm'
        call iosave(ch,outs,energy,-ifi,igets('ctrl nvario',sctrl))
        call iosave(ch,outs,energy,-lgunit(1),
     .    igets('ctrl nvario',sctrl))
C        call iosave(ch,outs,energy(2),-ifi,4)
C        call iosave(ch,outs,energy(2),-lgunit(1),4)
        call fclose(ifi)
C#ifdefC MPI
C        endif
C#endif
        if (icom .eq. 1) then
            call info0(20,1,1,' '//prgnam//'%a:  Relaxation complete')
          call rlse(onpr)
          goto 3
        endif
        call rlse(onpr)
C --- End of relaxation/MD loop ---
    2 continue
    3 continue
      if (tbe) then
        return
      endif
C --- End of TBE branch ---

C --- TBBND and TBFIT branch ---
      if (.not. (tbbnd .or. tbfit)) goto 5
C --- Get efmax, nevmx
      if (tbbnd) then
        nevmx = 0
        call upack1('bz efmax',sbz,efmax)
C ---   Get Fermi energy ---
        nfilem = fopn('BAND')
        efermi(1) = 999
        i = getef(nfilem,111,efermi(1))
        call info2(21,0,0,'%N Found fermi energy : %d',efermi(1),0)
        call fclose(nfilem)
      endif

C --- Set up permutation matrix for Bloch ---
      ndim = nbas*nl**2
      call defi(oidxsh,ndim)
      if (iprint() .ge. 30) print *
      call makidx(nl,1,1,nbas,0,sspec,w(oips),-1,w(oidxsh),inxsh)
      ldim = inxsh(1)*nspc
      idim = inxsh(2)*nspc - ldim
      call rxx(idim .ne. 0,prgnam//': bad orbital switches')
      nlb  = ldim
      call rxx(nlb .gt. nlbmx,'increase nlbmx in main')
      call defdc(oe,ldim**2)

      if (tbbnd) then
        if (fxst('SYML') .ne. 1) then
          print *, ' '
          print *,' TBZINT can''t find symmetry line file!'
          print *, ' '
          goto 5
        endif
        ifi =  fopno('SYML')
        ifib = fopnn('BNDS')
        write(ifib,112) nlb,efermi(1),0
  112   format(i5,f10.5,i5)
      endif
      if (tbfit) then
C       Tony used nsp=1 for nspc=2.  Won't work until this is fixed.
        if (nspc .eq. 2  .or. nsp .eq. 2)
     .    call rx('need to restore nsp=2 to SO coupled case')
        ifib = fopno('BNDS')
        read (ifib,*) nlb
        call pshpr(0)
C       call suqlsr(1,ifib,nlb,.false.,.false.,nq,w,w)
        call suqlsr(1,ifib,max(nspc,nsp),nlb,ldim,1,ldim,.false.,
     .    .false.,nq,w,w)
        call poppr
        call defdr(oqp, 3*nq)
        call defdr(oeband, -ldim*max(nspc,nsp)*nq)
        rewind ifib
        read (ifib,*) nlb
        call suqlsr(6,ifib,max(nspc,nsp),nlb,ldim,1,ldim,.false.,.false.
     .    ,nq,w(oqp),w(oeband))
        call sp2cls('spec iq1',sspec,w(oics),nl*nsp,1,nclasp,oiq1)
        call sp2cls('spec ivso',sspec,w(oics),nl,1,nclasp,oivso)
        call defi(oifit, 2*nq)
        call fitmrq(sctrl,sarray,slat,stb,nterm,nlmesp,nset,nl,nsp,nspc,
     .    nclass,npair,nq,ldim,nelts,mod(memode,10),mxnbr,idim,w(oclabl)
     .    ,w(oitab),w(oidec),w(oitbcf),w(oidcf),w(oitbov),w(oidov),
     .    w(oitocf),w(oidocf),w(oiq1),w(oivso),w(oidxdn),w(oiam),w(onpm)
     .    ,w(olmx),w(oidxsh),rmaxh,w(oeband),w(odel),w(oqp),w(ohso),
     .    w(oV0),w(opnu),w(oifit),w(otabme),w(odecay),w(otabcf),
     .    w(odeccf),w(otabov),w(odecov),w(otbocf),w(odcocf),w(oqnu),
     .    w(ovso),w(opot0),rl)
        call rlse(oifit)
        return
      endif

C --- Loop through all lines of k-points ---
      rewind ifi
   91 continue
      read(ifi,*,end=92,err=92) nq,q1,q2
      if (iprint() .ge. 30) write(*,113) nq,q1,q2
  113 format(/' nq=',i3,'   q1=',3f7.4,'   q2=',3f7.4)
      write(ifib,114) nq*nsp
  114 format(i5)
      do  iq = 1, nq
      do  isp = 1, nsp
        xx = dble(iq-1)/dble(nq-1)
        q(1) = xx*q2(1) + (1-xx)*q1(1)
        q(2) = xx*q2(2) + (1-xx)*q1(2)
        q(3) = xx*q2(3) + (1-xx)*q1(3)
        if (iprint() .ge. 40) write(*,115) q
  115   format('   Q=',3F9.6)

C --- Generate the eigenvalues into evl ---
        ipr = 0
        if (iprint() .gt. 40) ipr = 3
        call secmtb(sctrl,plat,nbas,nl,nspc,nsp,isp,w(olmx),w(oipc),
     .    w(oidxsh),ldim,nevmx,efmax,iq,nq,q,nsites,w(oiax),
     .    w(onpr),w(oh),w(ovso),w(ohso),w(oov),w(opot0),sc,ipr,F,nev,
     .    w(oe),evl,w(1),w(1))
        write(ifib,116) q, (evl(i),i=1,ldim)
  116   format(3f10.5/(10f8.4))
      enddo
      enddo
      goto 91
   92 continue
      if (iprint() .ge. 30) print *
      write(ifib,114) 0

    5 continue
      end

