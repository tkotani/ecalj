Cgetarg...info...           structure ['atscpp', 'ctrl', 'sctrl']
Cgetarg...info...           structure ['atscpp', 'spec', 'sspec']
Cgetarg...info...           use_to_add ['atscpp', 'm_struc_def']
#define ONE_ATOM_UNIT 1
      subroutine atscpp(sctrl,sspec,is,ic,dclabl,rmax,imake,nl,nsp,
     .  initc,rhozbk,avw,pnu,qnu,ves,dv,eula,neul,bxc,ekap,pp,pprel,
     .  etot,sumev,qtot,amgm,amom,rhrmx,vrmax,thrpv,sop,gradrm,pmpol,
     .  vintra,clabl,eterms)
       
       use m_struc_def  !Cgetarg

C- Make one atom self-consistent and generate new potential parameters
C ----------------------------------------------------------------
Ci Inputs
Ci   sctrl :struct containing parameters governing program flow
Ci   sspec :struct containing species-specific information
Ci   is    :species index
Ci   ic    :class index
Ci   dclabl:class names, packed as real numbers
Ci   rmax  :augmentation radius, in a.u.
Ci   imake :0  Read double counting terms from atom files.
Ci         :   No potential or pot pars calculated; no atom file written
Ci         :   This is a 'cheap' version of imake=4, where
Ci         :   instead of computing the d.c. terms from the potential,
Ci         :   they are copied from atom files.
Ci         :1  Make self-consistent potential from given pnu,qnu
Ci         :   No potential parameters calculated.
Ci         :2  Like imake=0, but make ppars from disk potential
Ci         :   Passed pnu,qnu written to disk.
Ci         :   NB: potential, ppars need not be related to pnu,qnu
Ci         :3  Make self-consistent potential from given pnu,qnu
Ci         :   and potential parameters from resulting potential
Ci         :4  Make sphere double-counting terms from disk potential
Ci         :   and supplied moments P,Q.  No internal self-consistency
Ci         :   in the potential, nor are potential parameters generated;
Ci         :   no atom file written.  Using this mode to make terms
Ci         :   in Kohn-Sham energy, e.g. v_in, rho_out.
Ci   nl    :(global maximum l) + 1, for dimensioning pp
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   initc : record of what parameters are available.
Ci           1 P,Q   2 pp   4 sop   8 vintra  16 pmpol  32 gradrm
Ci   rhozbk: constant nuclear background density (jellium)
Ci   avw   :length scale, usu. average Wigner-Seitz sphere radius
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
Ci          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Ci   qnu   :energy-weighted moments of the sphere charges
Ci   ves   :l=0 electrostatic potential at rmax (see Remarks)
Ci   dv    :constant potential shift, passed through to vxc0sp
Ci   ekap  :LMTO energy
Ci   rhozbk:constant background charge, excluded in electrostatics
Ci Inputs/Outputs
Cio  initc : record of what parameters are available.
Cio        : 1 P,Q   2 pp   4 sop   8 vintra  16 pmpol  32 gradrm
Cio        : may be modiifed on output
Co Outputs
Co   etot,sumev: total and band structure energy of atom, which
Co         :can be made within KKR viewpoint if moments P,Q
Co         :generate sphere potential
Co   qtot  :total charge within sphere
Co   amgm  :difference between spin up and spin down charge (mag mom)
Co   amom  :(noncollinear case) projection of moment along bxc
Co   rhrmx :density at rmax (not generated for imake=4)
Co   vrmax :total l=0 potential at rmax boundary, by class
Co         :(not generated for imake=4)
Co   sop   :matrix elements of spin orbit coupling
Co   gradrm:radial matrix elements of gradient operator; see rgrme
Co   pmpol :integral (phi-or-phidot * phi-or-phidot * r**l) :
Co         :matrix elements of w.f. * wf * r**l for multipole moments
Co   pp    :(depending on imake) potential parameters
Co   bxc   :(noncollinear) orientation of XC field
Co         :NB: if imake=4, bxc is an input, not an output
Co   pprel :(depending on imake) relativistic potential parameters
Co   eterms:integrals for the total E are accumulated for this sphere
Co         :(1)  ehar   --- not touched here
Co         :(2)  eks    --- not touched here
Co         :(3)  utot   = total electrostatic energy
Co         :(4)  valves --- not used by ASA
Co         :(5)  cpnves --- not used by ASA
Co         :(6)  rhoexc = rho * exc
Co         :(7)  rhovxc = rho * vxc (not needed for total energy)
Co         :(8)  sumec  = sum-of-core eigenvalues
Co         :(9)  sumtc  = sum-of-core K.E (not needed for total energy)
Co         :(10) xcore  = rhoc * total potential
Co         :(11) valvef = rhov * total potential
Co         :(12) sumt0  --- not used by ASA
Co         :(13) dq1    --- not used by ASA
Co         :(14) dq2    --- not used by ASA
Co         :(15) amgm   = system magnetic moment
Co         :(16) sumev  = sphere sum-of-eigenvalues
Co         :(17) rinvxt --- not touched here
Co         :(18) rouvxt --- not touched here
Co         :(19) bmval  = M<B> : magnetic contribution to valvef (asadc)
Cl Local variables
Cl   lfree :true for free atom (program overrides passed rmax);
Cl   lfrz  :true for frozen core
Cl   lso   :true if to calc spin-orbit coupling parameters.
Cl   lintra:true to calculate intra-atomic dC_i/dq_j
Cl   nrmix :1, maximum number of iterations towards self-consistency
Cl          2, number of previous iterations Anderson mixing for charge
Cl             nrmix>10: no amix; set beta to nrmix/100
Cl    nr   :number of points on the radial mesh
Cl   havedc:T if double-counting terms returned; otherwise F
Cr Remarks
Cr   Charges qtot and amgm are always calculated.  Sign
Cr   convention of charge is + for electrons, - for nucleii.
Cr
Cr   Potential parameters and total energy are calculated from atomsc
Cr   subject to boundary condition that the potential at rmax
Cr   is zero, regardless of the total net charge inside the sphere.
Cr   See subroutine madpot for discussion of how this choice affects
Cr   the Madelung energy.
Cr
Cr   The potential parameters, however may be shifted by the
Cr   electrostatic energy ves depending on the value of idmod.
Cr   In the default mode of operation (idmod 0 or 1), potential
Cr   parameters pp are calculated for the enu corresponding to the
Cr   pnu; the pp's enu and c are shifted however by the electrostatic
Cr   potential ves so that the pp's as returned to caller correspond
Cr   to electrostatic potential = ves at rmax.  If idmod is 2, the
Cr   potential parameters are generated around enu given by the input
Cr   enu (again assuming that enu corresponds to electrostatic
Cr   potential = ves at rmax).  There is no shifting of enu and c
Cr   unless the parameters are generated internally from potpar.
Cr
Cr   Regardless of the value of idmod, self-consistency is achieved
Cr   keeping pnu and moments fixed.  This means that if any idmod is 2
Cr   the potential and moments are self-consistent with respect
Cr   to the potential parameters.
Cr
Cb Bugs
Cb   Inconsistency in lmx(ic) always assumed to be nl-1 (see sop)
Cu Updates
Cu   26 Oct 08 Do not implement LDAU, IDU modes 4,5 here (move to suldau)
Cu   21 Dec 05 (wrl) allow potential shifts to mimic LDA+U
Cu   29 Sep 04 Reads/writes relativistic ppar's
Cu   18 Jun 04 (A Chantis) working fully relativistic code
Cu   21 Apr 04 Changes for l- and m-dependent XC fields
Cu    4 Apr 04 Additions for magnetic field in the presence of
Cu             orbital-dependent XC field.
Cu   19 Sep 03 (ATP) Enabled partial core occupation (core holes)
Cu   18 Mar 03 (A Chantis) relativistic potential parameters.
Cu             Altered argument list.
Cu   28 Feb 03 imake=4 implemented.  New argument eterms
Cu   22 Feb 03 Make and Printout <v_xc> inside sphere
Cu   15 Feb 03 SO parameters now include matrix element for ASA Bfield
Cu   22 Dec 01 Adjustments to accomodate changes in phidx
Cu   28 Apr 98 potpar can make gradrm
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nrmx
      parameter (nrmx=1501)
      integer is,ic,nl,nsp,nr,imake,initc,neul
Cgetarg       double precision dclabl(ic),sspec(1),sctrl(1),rmax(ic)
       real(8):: dclabl(ic) , rmax(ic) 
       type(s_spec)::sspec(*)
       type(s_ctrl)::sctrl

      double precision pnu(ic),qnu(ic),pp(6,nl,nsp,ic),bxc(3,ic),
     .  pprel(4,nl,2*nl,2,2,ic),thrpv,rhozbk,avw
      double precision etot,qtot(ic),amgm,amom,eterms(20),
     .  rhrmx(ic),vrmax(2,ic),ves(ic),dv(ic),eula(neul,3),
     .  pmpol(nl,nl,2*nl-1,3,nsp,ic),vintra(nl*nl,nsp*nsp,ic),
     .  sop(0:nl-1,nsp,nsp,9,ic),gradrm(4,2,2,nl,nsp,ic)
      character clabl*8
C Heap allocation
      integer w(1)
      common /w/ w
C Local parameters
      logical lrell,lfree,lgrad,lfrz,lso,lbf,lintra,loptc,lgors,lmpol,
     .  lgdd,lneedv,havedc
      integer ifi,nitmax,isw,ncmx,nvmx,ipr,lmx,igets,fopna,i,j,nlspic,
     .  lgunit,nmix,nn,iwk(10),nrmix(2),n0,nglob,lrel,kcor,lcor
      integer ov,orofi,orhoi,orho,orhoc,og,ogp,owgt
      parameter (n0=10, ncmx=50, nvmx=20)
      integer idmod(n0),idu(4)
      double precision a,sumec,sumtc,sumev,ekin,utot,rhoeps,rhomu,rhov,
     .  rmx,ekap,exc(2),thrpvl(10),xx(2),z,qc,ec(ncmx),ev(nvmx),avvxc,
     .  xcore,qval,pnuloc(nl*nsp),qnuloc(3*nl*nsp),bhat(3),bmval,
     .  ddot,qcor(2)
      double precision uh(4),jh(4)
C     double precision pz(n0,2),qz(3,n0,3)
C     integer idmoz(n0)
      character job*3, lbl*8, outs*80, outs1*80
      logical lscf
      logical aiogen,aiomom,aiopar,aiopot,aiosop,aiorme,aiocor,aiova,
     .  aiomp,sw,swc,havev,haveso,havegr,havemp,haveva,cmdopt
      external lgunit,aiogen,aiomom,aiopar,aiopot,aiocor,aiosop,aiova,
     .  aiorme,asprjq

C --- Setup ---
C     lscf: T, make sphere self-consistent given moments
C           F, only calculate charges qc, qtot and amgm
      call tcn('atscpp')
      call getpr(ipr)
      call r8tos8(dclabl(ic),clabl)
      lrel   = nglob('lrel')
      lscf   = imake.eq.1 .or. imake.eq.3
      lneedv = imake.eq.2 .or. imake.eq.4

C$$$Cgetarg       lfree  = lgors('ctrl lasa,8',sctrl)
C$$$       lfree = iand(8,int(sctrl%lasa)) .ne.0 
C$$$Cgetarg       lgdd   = .not. lgors('ctrl lasa,128',sctrl)
C$$$       lgdd = .not. iand(128,int(sctrl%lasa)) .ne.0 
C$$$Cgetarg       lmpol  = lgors('ctrl lasa,32',sctrl)
C$$$       lmpol = iand(32,int(sctrl%lasa)) .ne.0 
ctakao removed asa
       lfree = .false.
       lgdd  = .true.
       lmpol = .false.



Cgetarg       lfrz   = lgors('ctrl lcd,1',sctrl)
       lfrz = iand(1,int(sctrl%lcd)) .ne.0 

Cgetarg       lso    = lgors('ctrl lncol,4',sctrl)
       lso = iand(4,int(sctrl%lncol)) .ne.0 

Cgetarg       lbf    = lgors('ctrl lncol,8',sctrl)
       lbf = iand(8,int(sctrl%lncol)) .ne.0 

Cgetarg       lrell  = lgors('ctrl lrel,-1',sctrl)
       lrell = iand(-1,int(sctrl%lrel)) .ne.0 

Cgetarg       lintra = mod(igets('ctrl lscr',sctrl),100) .ge. 10 .or.
Cgetarg      .         mod(igets('ctrl lsx',sctrl)/10,2) .ne. 0
       lintra = mod ( int(sctrl%lscr) , 100 ) .ge. 10 .or. mod ( int(sctrl%lsx) 
     . /10 , 2 ) .ne. 0 

Cgetarg       loptc  = igets('ctrl loptc',sctrl) .gt. 0
       loptc = int(sctrl%loptc) .gt. 0 

      call dpzero(bhat,3)
Cgetarg...info...           integer :: nrmix(2) , integer(8) :: sctrl%smalit 2
Cgetarg       call upack('ctrl smalit',sctrl,nrmix,0,0,0,0)
       
       i_copy_size=size(sctrl%smalit) 
       call i8icopy(i_copy_size,sctrl%smalit,1,nrmix,1) 

Cgetarg...info...           integer :: idmod(n0) , integer(8) :: sspec%idmod 10
Cgetarg...info...           integer :: lmx , integer(8) :: sspec%lmxa 1
Cgetarg...info...           real(8) :: z , real(8) :: sspec%z 1
Cgetarg       call upack('spec idmod lmxa z',sspec,is,idmod,lmx,z,0)
       
       i_copy_size=size(sspec(is)%idmod) 
       call i8icopy(i_copy_size,sspec(is)%idmod,1,idmod,1) 
       lmx=sspec(is)%lmxa
       z=sspec(is)%z

C     This is now handled by LDA+U
C     call upack('spec idu uh jh',sspec,is,idu,uh,jh,0)
      call iinit(idu,4)
Cgetarg...info...           real(8) :: a , real(8) :: sspec%a 1
Cgetarg...info...           integer :: nr , integer(8) :: sspec%nr 1
Cgetarg       call upack('spec a nr',sspec,is,a,nr,0,0)
       
       a=sspec(is)%a
       nr=sspec(is)%nr

      etot = 0
      amgm = 0
      sumev = 0
      call dpzero(exc,2)
      havedc = .false.
      nitmax = nrmix(1)
      nmix = nrmix(2)
      thrpv = 0
      rmx = rmax(ic)
      if (lfree) then
        rmx = 50d0
        if (z .lt. 10) rmx = 25
        if (z .le.  6) rmx = 20
        if (nsp .eq. 2 .and. nmix .le. 10) nmix = 0
      endif
      nlspic = nl*nsp*(ic-1)
      lgrad = nglob('lxcf')/10 .ne. 0
      call rmesh(z,rmx,isw(lrell),isw(lgrad),nrmx,a,nr)
C     if (ipr .gt. 20) print *
#if ONE_ATOM_UNIT
      ifi = fopna(clabl,30,0)
#else
      ifi = fopna(clabl,30+ic,0)
#endif
      call defdr(ov,   nr*nsp)
      call defdr(orofi,nr*2)
      call defdr(owgt, nr*2)
      call defdr(orho, nr*nsp)
      call defdr(orhoi,nr*nsp)
      call defdr(orhoc,nr*nsp)
      call defdr(og,   nr*2)
      call defdr(ogp,  nr*2*4)
      haveso = mod(initc/4,2) .eq. 1
      haveva = mod(initc/8,2) .eq. 1
      havemp = mod(initc/16,2) .eq. 1
      havegr = mod(initc/32,2) .eq. 1
      job = 'gue'
      havev = aiopot(nr,nsp,a,rmx,bhat,w(ov),ifi)
C     if no bhat read in, flag that it is not there
      if (ddot(3,bhat,1,bhat,1) .eq. 0) bhat(1) = -99d0
      if (havev) job = 'pot'
      call awrit0('%x     ... available from input:',outs,80,0)
      if (mod(initc,2) .eq. 1) call awrit0('%a  p,q',outs,80,0)
      if (mod(initc/2,2) .eq. 1) call awrit0('%a  ppar',outs,80,0)
      if (haveso) call awrit0('%a  so-par',outs,80,0)
      if (havegr) call awrit0('%a  rgrad-me',outs,80,0)
      if (haveva) call awrit0('%a  vintra',outs,80,0)
      if (havemp) call awrit0('%a  mp-par',outs,80,0)
      if (havev) call awrit0('%a  V',outs,80,0)
      swc = aiocor(nr,nsp,a,rmx,w(orhoc),sumec,sumtc,ifi)
      if (swc) call awrit0('%a  core',outs,80,0)
C     sw=T flags that some required information is missing, for printout
      sw = (lneedv .and. .not. havev) .or.
     .     (lscf .and. mod(initc,2) .ne. 1)
      i = -lgunit(1)
      call awrit0('%x ATSCPP: class '//clabl,outs1,80,0)
      if (imake .eq. 0) call awrit0('%a   reading '//
     .  'double-counting terms from GEN',outs1,80,0)
      if (lscf) call awrit0('%a  making sphere potential from P,Q',
     .  outs1,80,0)
      if (lneedv) call awrit0(
     .  '%a  reading sphere potential from file',outs1,80,0)
      if (ipr .ge. 50 .or. sw) call awrit0('%a',outs1,-80,i)
      if (ipr .ge. 50 .or. sw) call awrit0('%a',outs,-80,i)
      if (lscf .and. mod(initc,2) .ne. 1)
     .  call rx('atscpp: missing P,Q for class '//clabl)
      if (imake .eq. 2 .and. .not. havev)
     .  call rx('atscpp: missing potential for class '//clabl)
      if (lfrz .and. .not. swc)
     .  call rx('atscpp: missing core for frozen core '//clabl)
      if (lfree .and. imake .ne. 1)
     .  call rx('atscpp: lfree=T compatible only with imake=1')
      if (imake .eq. 4 .and. .not. havev) return

C ... Scaling of moments by Euler angle rotation
      j = 23
C     Case bhat is input
      if (imake .eq. 4) j = 21
C     Rigid spin: no need to project moments onto Bxc
      if (neul .eq. 1) j = j-10
C     Collinear : only copy qnu to qnuloc
      if (neul .le. 0) j = 10
      if (cmdopt('--keepsignm',11,0,outs)) j = j+100
      call asprjq(j,clabl,nl,nsp,eula,neul,pnu(1+nlspic),
     .  qnu(1+3*nlspic),pnuloc,qnuloc,bxc(1,ic),amom)

C     Extract core hole parameters
      call gtpcor(sspec,is,kcor,lcor,qcor)

C --- Create potential from P,Q --
      if (lscf) then
        call getqvc(nsp,nl,lmx,z,pnuloc,qnuloc,0,0,kcor,lcor,qcor,qc,
     .    qtot(ic),amgm,0d0,0d0)
        ec(1) = 0
        call awrit3('%xATOM='//clabl//'%10pZ=%d  %?;n==1;frz;Qc=%d;',
     .    outs,80,0,z,lfrz,qc)
        call awrit6('%a  R=%1,6;6d  Qv=%1;6d  %??n==2?mom=%1;5d?%0d?'//
     .    '%a  a=%d  nr=%i',outs,80,0,rmx,qtot(ic),nsp,amgm,a,nr)
        do  33  j = 1, 2
   33   if (ipr .ge. 20 .or. j .eq. 2)
     .      call awrit0('%a',outs,-80,-lgunit(j))
        call atomsc(lgdd,nl,nsp,lmx,z,rhozbk,kcor,lcor,qcor,rmx,a,nr,
     .    w(orofi),ec,ev,pnuloc,qnuloc,idmod,w(ov),dv(ic),w(orhoi),
     .    w(orho),w(orhoc),nmix,qc,sumec,sumtc,sumev,ekin,utot,rhoeps,
     .    etot,amgm,rhrmx(ic),vrmax(1,ic),qtot(ic),exc,job,nitmax,lfrz)
        call radmwt(0,rmx,a,nr,w(orofi),w(owgt))
C        print *, '!!', ic,vrmax(1,ic)
C        call pshpr(99)
        call asadc(ic,nr,nsp,z,w(orofi),a,w(owgt),w(ov),w(orhoc),
     .    w(orho),rhov,bmval,xcore,rhoeps,rhomu,
     .    utot,qval)
C       call poppr
        eterms(3)  =  utot
        eterms(6)  =  rhoeps
        eterms(7)  =  rhomu
C       Suppress qc * ves(ic) to xcore and sumec; else add to both
        eterms(8)  =  sumec + ves(ic)*qc*1
        eterms(10) =  xcore + ves(ic)*qc*1
        eterms(9)  =  sumtc
C       Add qval * ves(ic) to rhov, since not incl. in asadc
        eterms(11) =  rhov + ves(ic)*qval
C       NB: in noncollinear case, this is moment || B
        eterms(15) =  amgm
C       Add qval * ves(ic) to sumev, since not incl. in atomsc
        eterms(16) =  sumev + ves(ic)*qval
C       Magnetic contribution to rhov
        eterms(19) =  bmval
        havedc = .true.
        if (.not. lfrz) swc = .true.

C --- Compute double-counting terms for HK total energy ---
      else if (imake .eq. 4) then
        call pshpr(0)
        call getqvc(nsp,nl,lmx,z,pnuloc,qnuloc,0,0,kcor,lcor,qcor,qc,
     .    qtot(ic),amgm,0d0,0d0)
        ec(1) = 0
        call atomsc(lgdd,nl,nsp,lmx,z,rhozbk,kcor,lcor,qcor,rmx,a,nr,
     .    w(orofi),ec,ev,pnuloc,qnuloc,idmod,w(ov),dv(ic),w(orhoi),
     .    w(orho),w(orhoc),nmix,qc,sumec,sumtc,sumev,ekin,utot,rhoeps,
     .    etot,amgm,xx(1),xx(2),qtot(ic),exc,job,0,lfrz)
        call radmwt(0,rmx,a,nr,w(orofi),w(owgt))
C        print *, '!!'
C        call setpr(99)
        call asadc(ic,nr,nsp,z,w(orofi),a,w(owgt),w(ov),w(orhoc),
     .    w(orho),rhov,bmval,xcore,rhoeps,rhomu,
     .    utot,qval)
        call poppr
        eterms(3)  =  utot
        eterms(6)  =  rhoeps
        eterms(7)  =  rhomu
C       Suppress qc * ves(ic) to xcore and sumec; else add to both
        eterms(8)  =  sumec + ves(ic)*qc*1
        eterms(10) =  xcore + ves(ic)*qc*1
        eterms(9)  =  sumtc
C       Add qval * ves(ic) to rhov, since not incl. in asadc
        eterms(11) =  rhov + ves(ic)*qval
        eterms(15) =  amgm
C       Magnetic contribution to rhov
        eterms(19) =  bmval
        havedc = .true.
        goto 99
C --- Else, try and read atomic parameters by aiogen ---
      else
        call getqvc(nsp,nl,lmx,z,pnuloc,qnuloc,0,0,kcor,lcor,qcor,qc,
     .    qtot(ic),amgm,0d0,0d0)
        havedc = aiogen(lbl,xx,xx,nn,nn,nn,nn,xx,xx,xx,vrmax(1,ic),
     .    sumec,sumev,thrpv,ekin,utot,rhoeps,etot,ifi)
        if (havedc) then
          qval = z - qc + qtot(ic)
          eterms(3)  =  utot
          eterms(6)  =  rhoeps
          eterms(8)  =  sumec + ves(ic)*qc*1
C         Add qval * ves(ic) to sumev, since not incl. in atomsc
          eterms(16) =  sumev + ves(ic)*qval
C         Not enough info to resolve xcore, rhov; stuff into rhov
          eterms(10) =  0 + ves(ic)*qc*1
          eterms(11) =  eterms(16) + sumec - ekin
        endif

C       vrmax from v(r) takes precedence, if it's available
        if (havev) then
          call dpscop(w(ov),xx,1,nr,1,1d0)
          if (nsp .eq. 2) then
            call dpscop(w(ov),xx,1,2*nr,2,1d0)
            vrmax(1,ic) = (xx(1)+xx(2))/2 - 2*Z/rmx
            vrmax(2,ic) =  xx(1)-xx(2)
          else
            vrmax(1,ic) = xx(1) - 2*Z/rmx
            vrmax(2,ic) = 0
          endif
        endif
      endif

C --- Make sphere potential parameters ---
      avvxc = 0
      if (imake .eq. 2 .or. imake .eq. 3) then
        call icopy(nl,idmod,1,iwk,1)

C        if (cmdopt('-elin=',6,0,outs)) then
C          i = 6
C          call rxx(.not. a2bin(outs,elin,4,0,' ',i,len(outs)),
C     .      'atomsr: failed to parse'//outs)
C          do  15  i = 0, nl*nsp-1
C            if (mod(idmod(1+mod(i,nl)),10) .eq. 0) then
C              pp(1,i+1,1,ic) = elin
C              iwk(mod(i,nl)+1) = 2
C            endif
C   15     continue
C        endif

C   --- Average exchange-correlation field ---
        if (nsp .eq. 2) then
          call radmwt(0,rmax(ic),a,nr,w(orofi),w(owgt))
          call savvxc(nr,w(orho),w(orhoc),w(ov),w(owgt),avvxc)
        endif

C   --- Second-generation potential parameters ---
        call radmsh(rmx,a,nr,w(orofi))
        call setcc(lrell)
C       print *, '!!'
C       call potpar(nl,nsp,lmx,z,rmx,avw,ekap,lso.or.lbf.or..true.,
        call potpar(nl,nsp,lmx,z,rmx,avw,ekap,lso.or.lbf.or.neul.gt.1,
     .    loptc,lmpol,a,
     .    nr,w(orofi),w(ov),pnuloc,iwk,ves(ic),exc,
     .    qnuloc,idu,uh,jh,thrpv,thrpvl,w(og),w(ogp),pp(1,1,1,ic),
     .    pprel(1,1,1,1,1,ic),sop(0,1,1,1,ic),pmpol(1,1,1,1,1,ic),
     .    gradrm(1,1,1,1,1,ic))
        haveso = haveso .or. (lso.or.lbf.or.neul.gt.1)
        havegr = havegr .or. loptc
        havemp = havemp .or. lmpol
C   ... Shift enu and c by crystal electrostatic potential
        call daxpy(nl*nsp,1d0,ves(ic),0,pp(1,1,1,ic),6)
        call daxpy(nl*nsp,1d0,ves(ic),0,pp(2,1,1,ic),6)
C   ... Second-order hamiltonian: set p^gamma = pph(4) = 0
Cgetarg         if (lgors('ctrl lham,3',sctrl)) then
         if ( iand(3,int(sctrl%lham)) .ne.0 ) then 

          do  18  i = 0, nl*nsp-1
   18     pp(4,i+1,1,ic) = 0
        endif

C   --- Intraatomic Coulomb d^2E/dq_i dq_j ---
        if (lintra) then
          call v0intr(nl,nsp,lmx,z,rhozbk,rmx,a,nr,w(orofi),
     .      pnuloc,qnuloc,w(ov),w(orhoi),w(orho),
     .      w(orhoc),w(og),w(ogp),nmix,nitmax,qc,lfrz,avw,ekap,2,
     .      vintra(1,1,ic))
          haveva = .true.
        endif
      endif

C --- Reopen atomic file and write atomic data ---
      if (imake .eq. 1) then
        initc = 1
      elseif (imake .ne. 0) then
      call dfclos(ifi)
#if ONE_ATOM_UNIT
      ifi = fopna(clabl,30,0)
#else
      ifi = fopna(clabl,30+ic,0)
#endif
      if (havedc)
     .  sw = aiogen(clabl,z,rmax(ic),lmx,nsp,lrell,nr,a,qc,qtot(ic),
     .  vrmax(1,ic),sumec,sumev,thrpv,ekin,utot,rhoeps,etot,-ifi)
      if (mod(initc,2) .eq. 1)
     .  sw = aiomom(clabl,pnu(1+nlspic),qnu(1+3*nlspic),idmod,
     .  nl,lmx,nsp,rhrmx(ic),vrmax(1,ic),-ifi)
      sw = aiopar(clabl,lrel,pp(1,1,1,ic),pprel(1,1,1,1,1,ic),ves(ic),
     .  nl,lmx,nsp,-ifi)
      if (haveso) sw = aiosop(clabl,sop(0,1,1,1,ic),nl,lmx,nsp,-ifi)
      if (havegr) sw = aiorme(clabl,gradrm(1,1,1,1,1,ic),nl,nsp,-ifi)
      if (haveva) sw = aiova(clabl,vintra(1,1,ic),nl,lmx,nsp,-ifi)
      if (havemp) sw=aiomp(clabl,pmpol(1,1,1,1,1,ic),nl,2*nl-2,nsp,-ifi)
      sw = aiopot(nr,nsp,a,rmax(ic),bhat,w(ov),-ifi)
      if (swc) sw = aiocor(nr,nsp,a,rmax(ic),w(orhoc),sumec,sumtc,-ifi)
      call fclose(ifi)
C ... Now both pp's and moments are available
      initc = 3
      if (haveso) initc = initc+4
      if (haveva) initc = initc+8
      if (havemp) initc = initc+16
      if (havegr) initc = initc+32
      endif

C --- Printout of atomic parameters ---
C  20 continue
      if (ipr .gt. 30 .and. imake .ne. 1) then
        write(lgunit(1),'(1x)')
        do  25  j = 1, 2
         if (havedc) then
         call awrit4(' v_rmax= %,6d%23petot= %,6d%?!n!%4f<v_xc>= %,6d',
     .     ' ',80,lgunit(j),vrmax(1,ic),etot,isw(avvxc.ne.0d0),avvxc)
         if (ipr .gt. 31 .and. lscf) then
           call awrit3(' thrpv=  %,6d%23pby l:%n:1,6d',
     .       ' ',80,lgunit(j),thrpv,lmx+1,thrpvl(1))
         endif
         endif
         write(lgunit(j),'(1x)')
         if (.not. lfree) sw = aiopar(clabl,lrel,pp(1,1,1,ic),
     .     pprel(1,1,1,1,1,ic),ves(ic),nl,lmx,nsp,-lgunit(j))
         if (haveso) sw =
     .     aiosop(clabl,sop(0,1,1,1,ic),nl,lmx,nsp,-lgunit(j))
         if (haveva) sw =
     .     aiova(clabl,vintra(1,1,ic),nl,lmx,nsp,-lgunit(j))
         if (havemp) sw =
     .     aiomp(clabl,pmpol(1,1,1,1,1,ic),nl,2*nl-2,nsp,-lgunit(j))
   25  continue
      endif

C --- Cleanup ---
   99 continue
      call fclr(clabl,ifi)
      call rlse(ov)
      call tcx('atscpp')

Cgetarg       end
       end subroutine atscpp 

      subroutine atomsc(lgdd,nl,nsp,lmax,z,rhozbk,kcor,lcor,qcor,rmax,a,
     .  nr,rofi,ec,ev,pnu,qnu,idmod,v,dv,rhoin,rho,rhoc,nrmix,qc,sumec,
     .  sumtc,sumev,ekin,utot,rhoeps,etot,amgm,rhrmx,vrmax,qtot,exrmax,
     .  job,niter,lfrz)

C- Makes an atomic sphere self-consistent and get atomic charges
C ----------------------------------------------------------------
Ci Inputs
Ci   lgdd  :T  add q2 phi phidd into the density
Ci         :F  add q2 <phi phidd> phi phi into the density
Ci         :NB: both produce the same integrated density.
Ci         :Which one should be used depends on the context.
Ci   nl    :leading dimension of pnu,qnu.
Ci         :Also, total charge cutoff; see Remarks
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   lmax  :maximum l for this site (but see Remarks)
Ci   z     :nuclear charge
Ci   rhozbk:constant nuclear background density (jellium) added to z
Ci   kcor  :(partial core occupation) p.q.n for occupation
Ci   lcor  :(partial core occupation) l quantum for occupation
Ci   qcor  :(partial core occupation) core charge and moment
Ci   rmax  :potential, density calculated to rmax
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :                 -//-
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,,
Ci         :pnu = .5 - atan(dnu)/pi + (princ.quant.number).
Ci   qnu   :energy moments of charge (see Remarks)
Ci   idmod :0,1 or 2, specifing how the enu is set for an l-channel
CCi   pz   :pnu for second panel (npan=2)
CCi   qz   :qnu for second panel (npan=2)
CCi   idmoz:idmod for second panel; specifies how semicore state
CCi        :is treated
CCi   npan :>1 for two-panel calculations
Ci   v     :spherical potential (job='pot', else input v not used)
Ci   dv    :constant potential shift added to sphere.
Ci         :Not used now.
Ci   rhoin:input density (job='rho', else input rhoin not used)
Ci         :used internally as a work array
Ci   nrmix :number of prior densities to mix to accelerate convergence
Ci         :to self-consistency (Anderson mixing); see Remarks.
Ci         :nrmix<0, linear mixing, with mixing beta |nrmix/100|
Ci   job   :(see Remarks)
Ci         :job='pot': start with potential
Ci         :job='rho': start with rhoin.
Ci         :job='gue': start with guessed charge density.
Ci   niter :number of iterations to attempt convergence to
Ci          self-consistency (see Remarks)
Ci   lfrz  :0 for soft core, >0 for frozen core
Cio Input/Outputs
Cio  ec    :core eigenvalues.  On input, these are guessed values.
Cio        :if input ec(1) = 0, atomsc makes an internal choice fo ec.
Cio  ev    :valence eigenvalues  On input, these are guessed values.
Cio        :if input ec(1) = 0, atomsc makes an internal choice for ev.
Co Outputs
Co   rofi  :dimensioned (nr,2).
Co         :(*,1) radial mesh points
Co         :(*,2) weights
Co   rho   :spherical charge density = 4 pi r**2 rhotrue
Co   rhoc  :core charge density (unchanged if lfrz=1)
Co   qc:   :core electronic charge
Co   sumec :core single particle energy
Co   sumtc :core kinetic energy (unchanged if lfrz=1)
Co   ekin  :total kinetic energy
Co   utot  :electrostatic energy
Co   rhoeps:exchange-correlation energy
Co   etot  :sphere total energy
Co   amgm  :difference between spin up and spin down charge
Co   rhrmx :true density at rmax
Co   vrmax :true potential at rmax
Co   qtot  :net charge in sphere
Co   exrmax:exchange-correlation energy at rmax
Cr Remarks
Cr   Boundary conditions pnu, moments qnu, and the electrostatic
Cr   potential at rmax are all that is required to uniquely determine
Cr   a self-consistent spherical charge and potential.  'job'
Cr   determines how the starting potential is made, but the final
Cr   potential  and density should be independent of the initial choice.
Cr
Cr   atomsc uses the boundary condition that the potential at rmax
Cr   is zero, regardless of the total net charge inside the sphere.
Cr   See subroutine madpot for discussion of how this choice affects
Cr   the Madelung energy.
Cr
Cr   Sphere total energy is sum of K.E., Hartree energy, XC energy:
Cr      etot = ekin + utot + rhoeps
Cr   The kinetic energy is computed via double-counting terms
Cr     ekin = sumev + sumec + dsumec - rhovh - rhomu
Cb Bugs
Cb   Total energy terms need to be cleaned up and simplified.
C ----------------------------------------------------------------
C     implicit none
C Passed parameters
      logical lfrz,lgdd
      character job*3
      integer nr,nsp,nl,nrmix,niter,kcor,lcor,ncmx,nvmx,lmax,
     .  idmod(0:nl-1)
      parameter (ncmx=50, nvmx=20)
      double precision ec(ncmx),ev(nvmx),rofi(nr,2),v(nr,nsp),
     .  rho(nr,nsp),rhoc(nr,nsp),rhoin(nr,nsp),pnu(nl,nsp),
     .  qnu(3,nl,nsp),z,rmax,a,qc,vrmax(2),exrmax(2),dv,rhrmx,
     .  rhozbk,qcor(2),amgm,ekin,etot,qtot,rhoeps,sumec,sumev,sumtc,
     .  utot
C     double precision z(nl,2),qz(3,nl,3)
C     integer idmoz(nl)
C  Local parameters
      logical last,cmdopt,ltmp
      integer iprint,ir,isp,iter,jmix,lgunit,lrel,nglob,ipr1,l,ii,stdo
      character strn*10
      double precision b,ddot,decay,dl,dold,drdi,drho,dsumec,ea,fac,pi,
     .  rho0t,rhomu,rhovh,rmsdel,ro,sum,tolrsq,vhrmax,vnucl,vrhoc,vsum,
     .  zvnucl,rvh(2),rho0(2),reps(2),rmu(2),sec(2),stc(2),sev(2),dasum
C for Anderson mixing:
      integer amix,nmix,orho
      double precision norm(10,10),awk(10,2),beta,beta1
      parameter (beta = 0.3d0)
c Heap allocation
      integer w(1)
      common /w/ w
C tolch is tolerance for change in the charge density, tolv for rseq
      double precision tolch,tl
      parameter (tolch=5d-5,tolrsq=1d-12)
      external lgunit

      if (lmax .ge. nl) call rx('atomsc:  lmax too large')
      lrel = nglob('lrel')
      stdo = lgunit(1)
      pi = 4d0*datan(1d0)
      b = rmax/(dexp(a*(nr-1)) - 1)
      vhrmax = 0d0
      nmix = min(max(nrmix,0),10)
      sec(2) = 0d0
      sec(1) = 0d0
      stc(2) = 0d0
      stc(1) = 0d0
      call defdr(orho,nr*nsp*2*(nmix+2))

C --- Core charge, radial mesh points and weights ---
      if (kcor .ne. 0) then
      if (qcor(1) .ne. 0 .or. qcor(2) .ne. 0) then
        call info5(30,0,0,'%9fAdd core hole:  kcor=%i  lcor=%i'//
     .    '  qcor=%d  amom=%d',kcor,lcor,qcor,qcor(2),0)
      endif
      endif
      call getqvc(nsp,nl,lmax,z,pnu,qnu,0,0,kcor,lcor,qcor,qc,qtot,amgm,
     .          ec,ev)
C --- Guesses for core and valence eigenvalues ---
      if (ec(1) .eq. 0)
     .call getqvc(nsp,nl,lmax,z,pnu,qnu,ncmx,nvmx,kcor,lcor,qcor,qc,
     .            qtot,amgm,ec,ev)
      call radmsh(rmax,a,nr,rofi)
      call radwgt(rmax,a,nr,rofi(1,2))

C --- Moments printout ---
      if (iprint() .ge. 20) then
        ltmp=dasum(lmax+1,qnu(2,1,1),3)+dasum(lmax+1,qnu(3,1,1),3).eq.0
        if (ltmp) then
          call awrit5('  Pl= %n:-1,1;5#8d%?!n==2! spn 2  %n:-1,1;5#8d',
     .      ' ',180,stdo,lmax+1,pnu,nsp,lmax+1,pnu(1,2))
          call dcopy(lmax+1,qnu,3,awk,1)
          call dcopy(lmax+1,qnu(1,1,nsp),3,awk(1,2),1)
          call awrit5('  Ql= %n:-1,1;5#8d%?!n==2! spn 2  %n:-1,1;5#8d',
     .      ' ',180,stdo,lmax+1,awk,nsp,lmax+1,awk(1,2))
        elseif (iprint() .ge. 30) then
          write(stdo,891)
  891     format(
     .      '   l',8x,'pl',11x,'q0',11x,'q1',11x,'q2',5x,' id ',6x,'dl')
          do  11  isp = 1, nsp
            do  11  l = 0, lmax
              dl = dtan(pi*(.5d0 - pnu(l+1,isp)))
              if (dabs(dl) .gt. 9999) dl = 0
              write(stdo,100) l,pnu(l+1,isp),(qnu(ii,l+1,isp),ii=1,3),
     .          idmod(l),dl
  100         format(i4,4f13.7,i4,f13.7)
   11     continue
        endif
      endif

C --- Initial charge density ---
      if (job .eq. 'pot') then
        call newrho(z,lrel,lgdd,nl,1,lmax,a,b,nr,rofi,v,rhoin,rhoc,kcor,
     .    lcor,qcor,pnu,qnu,sec,stc,sev,ec,ev,tolrsq,nsp,lfrz,000)
       if (niter .eq. 0) then
          call dcopy(nr*nsp,rhoin,1,rho,1)
          return
        endif
      else if (job .eq. 'gue') then
        decay = 1d0+z/10d0
        decay = dmin1(decay,5d0)
        decay = 5
C         write(6,335) decay
C  335   format(/' initialize exponential density with decay',f7.3)
        sum = 0d0
        do  4  ir = 1, nr
        ro = dexp(-decay*rofi(ir,1))*rofi(ir,1)**2
        rhoin(ir,1) = ro
    4   sum = sum+ro*a*(rofi(ir,1)+b)
        fac = z/(sum*nsp)
        do  5  ir = 1, nr
    5   rhoin(ir,1) = rhoin(ir,1)*fac
        if (nsp .eq. 2) then
          do  6  ir = 1, nr
    6     rhoin(ir,2) = rhoin(ir,1)
        endif
      else if (job .ne. 'rho') then
        call rx('atomsc: job not pot|rho|gue')
      endif
      call dpscop(rhoin,w(orho),nr*nsp,1,1+nr*nsp*(nmix+2),1d0)

C --- Start self-consistency loop ---
      drho = 100d0
      last = .false.
      if (iprint() .ge. 30) write(stdo,341)
      jmix = 0
      dold = 1
      beta1 = beta
      if (nrmix .lt. 0) beta1 = dble(-nrmix)/100
      do  35  iter = 1, niter
        tl = tolrsq
C       Hartree potential
        call addzbk(rofi,nr,nsp,rhoin,rhozbk,-1d0)
        call poiss0(z,a,b,rofi,rhoin,nr,vhrmax,v,rvh,vsum,nsp)
C       call prrmsh('vh atom ',rofi,v,nr,nr,nsp)

        call addzbk(rofi,nr,nsp,rhoin,rhozbk,1d0)
        vnucl = v(1,1)
C       Exchange-correlation potential
        if (last .and. iprint() .ge. 50) call pshpr(80)
        call vxc0sp(a,b,rofi,rhoin,nr,v,rho0,reps,rmu,nsp,exrmax)
        if (last .and. iprint() .ge. 50) call poppr
C       Get rhrmx, exrmax
        fac = 4*pi*rofi(nr,1)**2
        rhrmx = rhoin(nr,1)/fac
        if (nsp .eq. 2) then
          exrmax(2) = exrmax(1)
          rhrmx = rhrmx + rhoin(nr,2)/fac
        endif
        ipr1 = 0
        if (last .and. iprint().ge.40) ipr1 = 1
        if (last .and. iprint().gt.40) ipr1 = 2
        call newrho(z,lrel,lgdd,nl,1,lmax,a,b,nr,rofi,v,rho,rhoc,
     .    kcor,lcor,qcor,pnu,qnu,sec,stc,sev,ec,ev,tl,nsp,lfrz,ipr1)
        drho = 0
        sum = 0
        vrhoc = 0
        rho0t = 0
        do  40  isp = 1, nsp
          rho0t = rho0t + rho0(isp)
          sum = sum + ddot(nr,rofi(1,2),1,rho(1,isp),1)
          do  42  ir = 1, nr
          drdi = a*(rofi(ir,1) + b)
   42     drho = drho + rofi(ir,2)/drdi*dabs(rho(ir,isp)-rhoin(ir,isp))
          do  41  ir = 2, nr
          ea = (v(ir,isp)-2*z/rofi(ir,1))
   41     vrhoc = vrhoc + rofi(ir,2)*ea*rhoc(ir,isp)
   40   continue
C       call prrmsh('rho for atom ',rofi,rho,nr,nr,nsp)
        call dcopy(nr*nsp,rho,1,w(orho),1)
        jmix = amix(nr*nsp,min(jmix,nmix),nmix,0,beta1,iprint()-70,
     .              .9D0,norm,awk(1,2),w(orho),awk,rmsdel)
        call dpscop(w(orho),rhoin,nr*nsp,1+nr*nsp*(nmix+2),1,1d0)
C
        if (last) goto 90
        if (iprint() .ge. 41 .or. iprint() .ge. 30 .and.
     .     (drho .lt. tolch .or. iter .eq. niter-1 .or. iter .eq. 1))
     .   write(stdo,340) iter,sum,drho,vnucl,rho0t,vsum,beta1
  340   format(i5,f12.6,1p,e12.3,0p,f14.4,e14.4,f14.4,f7.2)
  341   format(/'  iter     qint',9x,'drho',10x,'vh0',10x,'rho0',
     .  10x,'vsum',5x,'beta')
        last = (drho .lt. tolch .or. iter .eq. niter-1)
        if (iprint() .gt. 100)  call query(' ',-1,0)
        jmix = jmix+1
C       Beta for next iteration
        beta1 = min(max((1-drho/dold)/beta1,beta1-.2d0,beta),
     .              1d0,beta1+.2d0)
        if (nmix .gt. 0 .and. drho .lt. 1) beta1 = 1
        if (nrmix .lt. 0) beta1 = dble(-nrmix)/100
        dold = drho
   35 continue
   90 continue
      call rlse(orho)
C     End of iteration loop
      if (iprint() .ge. 30) write(stdo,'(1x)')

C --- Collect terms for total energy ---
      rhoeps = 0
      rhomu  = 0
      sumev  = 0
      if (.not. lfrz) then
        sumec = 0
        sumtc = 0
      endif
      rhovh  = 0
      do  80  isp = 1, nsp
        if (nsp .eq. 2 .and. iprint() .gt. 30)
     .    write(stdo,230) isp,v(nr,isp)-2*z/rmax,
     .   sev(isp),sec(isp),rvh(isp),reps(isp),rmu(isp)
  230   format(' Spin',i2,':',
     .  /' vrmax=',f12.5,'    sumev= ',f12.5,'    sumec=',f12.5,
     .  /' rhovh=',f12.5,'    rhoeps=',f12.5,'    rhomu=',f12.5)
        rhoeps = rhoeps + reps(isp)
        rhomu = rhomu + rmu(isp)
        sumev = sumev + sev(isp)
        if (.not. lfrz) then
          sumec = sumec + sec(isp)
          sumtc = sumtc + stc(isp)
        endif
        rhovh = rhovh + rvh(isp)
  80  continue
      zvnucl = -z*vnucl
      utot = .5d0*(rhovh + zvnucl)
C ... Correction to core eigenvalues if sumec not obtained from this V
      dsumec = vrhoc - (sumec-sumtc)
      ekin = sumev + sumec + dsumec - rhovh - rhomu
      etot = ekin + utot + rhoeps
      if (iprint() .ge. 40) write(stdo,139) sumev,sumec,vnucl,rhovh,
     .  zvnucl,utot,rhomu,rhoeps,dsumec,ekin,sumtc,etot
  139 format(/' sumev=',f13.6,'    sumec =',f13.6,'   vnucl =',f13.6
     .       /' rhovh=',f13.6,'    zvnucl=',f13.6,'   utot  =',f13.6
     .       /' rhomu=',f13.6,'    rhoeps=',f13.6,'   dsumec=',f13.6
     .       /' ekin= ',f13.6,'    tcore =',f13.6,'   etot  =',f13.6)
      vrmax(1) = -2*z/rmax
      do  55  isp = 1, nsp
   55 vrmax(1) = vrmax(1) + v(nr,isp)/nsp
      vrmax(2) = 0
      if (nsp .eq. 2) vrmax(2) = v(nr,1)-v(nr,2)

C --- write out rho if requested ---
      if (cmdopt('--dumprho',8,0,strn)) then
        call prrmsh('rho for atom ',rofi,rho,nr,nr,nsp)
        call prrmsh('pot for atom ',rofi,v,nr,nr,nsp)
        call prrmsh('rhoc for atom ',rofi,rhoc,nr,nr,nsp)
        call defrr(orho,-nr*nl*nsp)
        call newrho(z,lrel,lgdd,nl,nl,lmax,a,b,nr,rofi,v,w(orho),rhoc,
     .    kcor,lcor,qcor,pnu,qnu,sec,stc,sev,ec,ev,tl,nsp,lfrz,000)
        call prrmsh('rhol for atom ',rofi,w(orho),nr,nr,nl*nsp)
        call rlse(orho)
      endif

Cgetarg       end
       end subroutine atomsc 

      subroutine addzbk(rofi,nr,nsp,rho,rhozbk,scale)

C     implicit none
      integer nr,nsp
      double precision rofi(1),rho(nr,1),rhozbk,scale
      integer ir,isp
      double precision s

      if (rhozbk .eq. 0) return
      s = 16*datan(1d0)*scale*rhozbk
      do  10  isp = 1, nsp
      do  10  ir = 2, nr
   10 rho(ir,isp) = rho(ir,isp) + s*rofi(ir)*rofi(ir)
Cgetarg       end
       end subroutine addzbk 

      subroutine poiss0(z,a,b,rofi,rho,nr,vhrmax,v,rhovh,vsum,nsp)

C- Hartree potential for spherical rho.
C  ---------------------------------------------------------------------
Ci Inputs:
Ci   z     :nuclear charge
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   rho   :spherical charge density times 4*pi*r*r
Ci   nr    :number of mesh points
Ci   vhrmax:value of v(nr)
Ci   nsp   :=1 spin degenerate, =2 non-degenerate
Co Outputs:
Co   rofi  :radial mesh points
Co   v     :spherical Hartree potential
Co   vsum  :integral over that potential which is zero at rmax.
Co   rhovh :integral of rho*vh, where vh is the electrostatic
Co         :potential from both electronic and nuclear charge
Cr Remarks:
Cr    Solves Poisson's equation for given spherical charge density
Cr    and a specified value vhrmax at rofi(nr).
Cr    rho =  4*pi*r*r*rhotrue  but  v = vtrue
C ----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      integer nr,nsp
      double precision rofi(nr),v(nr,nsp),rho(nr,nsp),rhovh(nsp)
      double precision pi,ea,rpb,a,b,z
      double precision rmax,r2,r3,r4,f2,f3,f4,x23,x34,cc,bb,dd,df,drdi,
     .        r,srdrdi,g,f,y2,y3,y4,ro,vnow,wgt,a2b4,vhrmax,vsum,vhat0
      integer ir,isp

      pi = 4d0*datan(1d0)
      ea = dexp(a)
      rpb = b
      do  5  ir = 1, nr
        rofi(ir) = rpb - b
        rpb = rpb*ea
    5 continue
      rmax = rofi(nr)

C --- Approximate rho/r**2 by cc + bb*r + dd*r*r near zero  ---
      r2 = rofi(2)
      r3 = rofi(3)
      r4 = rofi(4)
      f2 = 0d0
      f3 = 0d0
      f4 = 0d0
      do  75  isp = 1, nsp
        f2 = f2 + rho(2,isp)/r2**2
        f3 = f3 + rho(3,isp)/r3**2
        f4 = f4 + rho(4,isp)/r4**2
   75 continue
      x23 = (r3*r3*f2 - r2*r2*f3)/(r3 - r2)
      x34 = (r4*r4*f3 - r3*r3*f4)/(r4 - r3)
      cc = (r2*x34 - r4*x23) / (r3*(r2 - r4))
      bb = ((r2+r3)*x34 - (r3+r4)*x23) / (r3*r3*(r4-r2))
      dd = (f2 - bb*r2 - cc)/r2**2

C --- Numerov for inhomogeneous solution ---
      a2b4 = a*a/4d0
      v(1,1) = 1d0
      df = 0d0
      do  12  ir = 2, 3
        r = rofi(ir)
        drdi = a*(r + b)
        srdrdi = dsqrt(drdi)
        v(ir,1) = v(1,1) - r*r*(cc/3d0 + r*bb/6d0 + r*r*dd/10d0)
        g = v(ir,1)*r/srdrdi
        f = g*(1d0 - a2b4/12d0)
        if (ir .eq. 2) y2 = -2d0*f2*r2*drdi*srdrdi
        if (ir .eq. 3) y3 = -2d0*f3*r3*drdi*srdrdi
        df = f - df
   12 continue
      do 13  ir = 4, nr
        r = rofi(ir)
        drdi = a*(r + b)
        srdrdi = dsqrt(drdi)
        ro = 0d0
        do  76  isp = 1, nsp
          ro = ro + rho(ir,isp)
   76   continue
        y4 = -2d0*drdi*srdrdi*ro/r
        df = df + g*a2b4 + (y4 + 10d0*y3 + y2)/12d0
        f = f + df
        g = f/(1d0 - a2b4/12d0)
        v(ir,1) = g*srdrdi/r
        y2 = y3
        y3 = y4
   13 continue

C --- Add constant to get v(nr) = vhrmax ---
      vnow = v(nr,1) - 2d0*z/rmax
      do  27  ir = 1, nr
        v(ir,1) = v(ir,1) + (vhrmax - vnow)
   27 continue

C --- Integral vh and  rho * vh ---
      rhovh(1) = 0d0
      rhovh(nsp) = 0d0
      vsum = 0d0
      vhat0 = 0d0
      do  30  ir = 2, nr
        r = rofi(ir)
        drdi = a*(r + b)
        wgt = 2*(mod(ir+1,2) + 1)/3d0
        if (ir .eq. nr) wgt = 1d0/3d0
        ro = 0d0
        do  31  isp = 1, nsp
          rhovh(isp) =  rhovh(isp)
     .                + wgt*drdi*rho(ir,isp)*(v(ir,1) - 2d0*z/r)
          ro = ro + rho(ir,isp)
   31   continue
        vhat0 = vhat0 + wgt*drdi*ro*(1d0/r - 1d0/rmax)
        vsum = vsum + wgt*drdi*r*r*(v(ir,1) - vhrmax)
   30 continue
      vsum = 4d0*pi*(vsum - z*rmax*rmax)
      vhat0 = 2d0*vhat0 + 2d0*z/rmax + vhrmax
      v(1,1) = vhat0

C --- Copy to second spin channel if spin polarized ---
      if (nsp .eq. 1) return
      do  84  ir = 1, nr
        v(ir,2) = v(ir,1)
   84 continue
Cgetarg       end
       end subroutine poiss0 


C      subroutine fctp0(l,rofi,v,z,nr,nctp0,xrim,xmin,nsave)
CC- Initialize things for FCTP, which finds classical turning point
CC  ---------------------------------------------------
CCi Inputs
CCi   l     :angular quantum number
CCi   rofi  :radial mesh points
CCi   v     :spherical potential (atomsr.f)
CCi   z     :nuclear charge
CCi   nr    :number of radial mesh points
CCo Outputs:
CCo   nctp0:
CCo   xrim:
CCo   xmin:
CCo   nsave: Estimate for nctp
CC  ---------------------------------------------------
C      implicit none
C      integer l,nctp0,nr,nsave
C      double precision xmin,xrim,z,v(nr),rofi(nr)
CC Local
C      integer ir
C      double precision fllp1,r,x,xlast
C      fllp1 = l*(l + 1)
C      r = rofi(10)
C      x = fllp1/r/r - 2d0*z/r + v(10)
C      ir = 10
C   80 ir = ir + 1
C      xlast = x
C      r = rofi(ir)
C      x = fllp1/r/r - 2d0*z/r + v(ir)
C      if (x .le. xlast .and. ir .lt. nr) goto 80
C      nctp0 = ir - 1
C      xmin = xlast
C      r = rofi(nr)
C      xrim = fllp1/r/r - 2d0*z/r + v(nr)
C      if (xmin .ge. xrim-3d0) nctp0 = nr
C      if (xmin .ge. xrim-3d0) xmin = xrim
C      nsave = (nctp0 + nr)/2
C      end
C      subroutine fctp(e,nctp,nctp0,xrim,xmin,nsave,l,rofi,v,z,nr,a,b)
CC- Finds classical turning point for wave function
CC ----------------------------------------------------------------------
CCi Inputs
CCi   e     :energy
CCi   nctp0 :
CCi   xrim  :
CCi   xmin  :
CCi   nsave: Estimate for nctp
CCi   l     :l quantum number
CCi   rofi  :radial mesh points
CCi   v     :spherical potential
CCi   z     :nuclear charge
CCi   nr    :number of radial mesh points
CCi   a     :mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
CCi   b     :mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
CCo Outputs:
CCo   nctp:  rofi(nctp) is classical turning point
CCo   nsave: New estimate for NCTP, next iteration
CCr Remarks
CC  ---------------------------------------------------
C      implicit none
C      integer l,nctp,nctp0,nr,nsave
C      double precision a,b,e,xmin,xrim,z,v(nr),rofi(nr)
CC Local
C      integer irep,k,n1,n2,nlast,ntry
C      double precision dfdr,dvdr,fllp1,fntry,fofr,r,rtry,vme
C
C      fllp1 = l*(l+1)
C      if (nctp0 .eq. nr) then
C        nctp = nr
C        return
C      endif
C      if (e .gt. xrim) then
C        nctp = nr
C        return
C        endif
C      if (e .lt. xmin) then
C        nctp = 2
C        return
C      endif
C      n1 = nctp0
C      n2 = nr
C      nctp = nsave
C      nlast = -10
C      do  20  irep = 1, 20
C      if (nctp.gt.n2 .or. nctp.lt.n1) nctp=(n1+n2+1)/2
C        r = rofi(nctp)
C        vme = (v(nctp)-e)
C        if (nctp .gt. nr) call rx('fctp0: nctp gt nr')
C        k = min(nctp,nr-1)
C        dvdr = (v(k+1)-v(k-1))/(2.d0*a*(r+b))
CC       dvdr = (v(nctp+1)-v(nctp-1))/(2.d0*a*(r+b))
C        fofr = fllp1/r/r - 2.d0*z/r + vme
C        dfdr = -2.d0*fllp1/r/r/r+2.d0*z/r/r + dvdr
C        rtry = r - fofr/dfdr
C        rtry = dmax1(rtry,rofi(2))
C        fntry = dlog(rtry/b+1.d0)/a + 1.d0
C        ntry = fntry + .5d0
CC|      write(6,810) irep,n1,nctp,n2,fntry,ntry
CC|810   format(i6,'   n1,nctp,n2=',3i5,'   ntry=',f8.3,i6)
C        if (nlast .eq. nctp) goto 98
C        if (fofr .gt. 0.d0) n2=nctp
C        if (fofr .lt. 0.d0) n1=nctp
C        nlast = nctp
C        nctp = ntry
C   20 continue
C   98 if (nctp .eq. nctp0+1) nctp = 2
C      nsave = nctp
C      end

      subroutine fctp0(l,nr,rofi,v,z,nctp0)

C- Initialize things for FCTP, which finds classical turning point
C ----------------------------------------------------------------------
Ci Inputs:
Ci   l     :angular momentum
Ci   rofi  :radial mesh points
Ci   v     :spherical potential
Ci   z     :nuclear charge
Ci   nr    :number of mesh points
Co Outputs:
Co   nctp0 :minimum of effective potential
Co          or nr if v(nr) > vmin + 3 Ry
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer l,nr,nctp0
      double precision z,v(nr),rofi(nr)
C Local variables:
      integer ir,irmin
      double precision veff,zz,fllp1,vi,vim1
      parameter (irmin=11)
C Statement functions:
      veff(ir)=fllp1/(rofi(ir)*rofi(ir))-zz/rofi(ir)+v(ir)

      zz = z+z
      fllp1 = l*(l+1)

      ir = irmin
      vi = veff(ir)
      vim1 = veff(ir-1)
   10 if (vi.le.vim1 .and. ir.lt.nr) then
        ir = ir+1
        vim1 = vi
        vi = veff(ir)
        goto 10
      endif

      nctp0 = ir-1
      if (veff(nctp0).ge.veff(nr)-3.d0) nctp0 = nr

Cgetarg       end
       end subroutine fctp0 

      subroutine fctp(a,b,e,l,nctp0,nr,rofi,v,z,nctp)

C- Finds classical turning point
C  ---------------------------------------------------------------------
Ci Inputs:
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :                 -//-
Ci   e     :energy
Ci   l     :angular momentum
Ci   nctp0 :minimum of effective potential (fctp0)
Ci   nr    :number of mesh points
Ci   rofi  :radial mesh points
Ci   v     :spherical potential
Ci   z     :nuclear charge
Ci Inputs/Outputs:
Cio  nctp  :rofi(nctp) is classical turning point
Cu Updates
Cu   13 Jun 00 Added safety checks to guard against jumps in potential
C  ---------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer nctp,nctp0,l,nr
      double precision e,rofi(nr),v(nr),z,a,b
C Local variables:
      integer ir,irep,n1,n2,nlast,ntry
      double precision r,veff,fllp1,fofr,dfdr,rtry,zz
C Intrinsic functions:
      intrinsic dlog,dmax1,min0
C Statement functions:
      veff(ir)=fllp1/(rofi(ir)*rofi(ir))-zz/rofi(ir)+v(ir)

      zz=z+z
      fllp1 = l*(l+1)

      if (nctp0.eq.nr .or. e.gt.veff(nr)) then
        nctp = nr
      elseif (e.lt.veff(nctp0)) then
        nctp = 2
      else
        n1 = nctp0
        n2 = nr-1
        nlast = -10
        do  10  irep = 1, 20
          if (nctp.gt.n2 .or. nctp.lt.n1) nctp = (n1 + n2)/2
          r = rofi(nctp)
          fofr =  veff(nctp)-e
          dfdr =-(fllp1+fllp1)/r/r/r + zz/r/r +
     .           (v(nctp+1) - v(nctp-1))/(2.d0*a*(r+b))
          rtry = dmax1(r-fofr/dfdr,rofi(2))
          ntry = dlog(rtry/b+1.d0)/a + 1.5d0

C         If there was a large change, check for safety
          if (nlast .eq. nctp .and. iabs(ntry-nctp) .gt. 10) then
            if (fofr .gt. 0) then
              do  ntry = nctp-1, nctp0+1, -1
                fofr =  veff(ntry)-e
                if (fofr .lt. 0) goto 20
                nctp = ntry
                n2 = ntry
              enddo
   20         continue
            else
              do  ntry = nctp+1, nr-1
                fofr =  veff(ntry)-e
                if (fofr .gt. 0) goto 21
                nctp = ntry
                n1 = ntry
              enddo
   21         continue
            endif
          endif

C         Exit point
          if (nlast .eq. nctp) then
            if (nctp .eq. nctp0+1) nctp = 2
            return
          endif
          if (fofr .gt. 0.d0)  then
            n2 = nctp
          else
            n1 = nctp
          endif
          nlast = nctp
          nctp = min0(ntry,nr-1)
   10   continue
        if (nctp .eq. nctp0+1) nctp = 2
      endif
Cgetarg       end
       end subroutine fctp 


      subroutine newrho(z,lrel,lgdd,nl,nlr,lmax,a,b,nr,rofi,v,rho,rhoc,
     .  kcor,lcor,qcor,pnu,qnu,sumec,sumtc,sumev,ec,ev,tol,nsp,lfrz,ipr)

C- Makes spherical charge density for a spherical potential.
C  ---------------------------------------------------
Ci Inputs:
Ci   z     :nuclear charge
Ci   lrel  :0 for nonrelativistic, 1 for relativistic
Ci   lgdd  :T q2 is coefficient to phidot**2 + phi*phidotdot
Ci         :F q2 is coefficient to phidot**2 - p phi**2
Ci         :Both produce the same integrated density; see Remarks.
Ci         :lgdd=F follows Stuttgart conventions.
Ci   nl    :(global maximum l) + 1
Ci   nlr   :second dimension of rho:
Cr         :1 if spherical part of rho is to be generated
Ci         :nl if generated rho is to be decomposed by l
Ci         :In the latter case, the core is not calculated
Ci   lmax  :maximum l for a given site
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nr    :number of radial mesh points
Ci   rofi  :radial mesh points
Ci   v     :spherical potential (atomsr.f)
Ci   kcor  :(partial core occupation) p.q.n for occupation
Ci   lcor  :(partial core occupation) l quantum for occupation
Ci   qcor  :(partial core occupation) core charge and moment
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,,
Ci         :pnu = .5 - atan(dnu)/pi + (princ.quant.number).
Ci   qnu   :energy moments of charge (see Remarks)
Ci   tol   :precision to which wave functions are integrated.
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   g     :normalized wave function times r
Ci   gp    :energy derivatives of g
Ci   lfrz  :T, do not make core rho
Ci   ipr   :0 no printout
Ci         :1 summary printout of core
Ci         :2 detailed printout of core
Co Outputs:
Co   rho   :spherical charge density times 4*pi*r*r
Co   rhoc  :core density times 4*pi*r*r
Co   sumec :sum of core eigenvalues
Co   sumtc :core kinetic energy = sumec - v*rhoc
Co   sumev :sum of valence eigenvalues
Co   ec    :core eigenvalues
Co   ev    :valence eigenvalues
Cr Remarks:
Cr   rho is determined by boundary conditions pnu and moments qnu.
Cr   For rmax>10 sets phidot, phidotdot are made zero.
Cr
Cr   Recall that val,slo correspond to u = r*phi, so
Cr     rmax * slo / val = D_u = D_phi + 1.
Cr     For val=rmax  slo = D + 1
Cr
Cr   Switch lgdd concerns the contribution of <phidot|phidot> to the
Cr   sphere charge.  The zeroth moment may be defined as the amount of
Cr   charge in channel l, i.e.
Cr     q^0 = (amount of <phi|phi>) +  p* (amount of <phidot|phidot>)
Cr   where
Cr     p = <phidot|phidot>
Cr   Then q^0 is not the amount of phi*phi inside the sphere, but rather
Cr   but instead, it is
Cr     (q_0 - p q_2)
Cr   since q2 is the amount of <phidot|phidot>.  The charge density is
Cr     rho= (q_0 - p q_2) phi*phi + 2 q_1 phi*phidot + q_2 phidot*phidot
Cr   This is the Stuttgart convention (lgdd=F)
Cr
Cr   Methfessel convention:  to avoid explicit dependence of rho on
Cr   p, he approximated p*phi*phi with -phi*phidotdot (they have the
Cr   same integrated charge).  Then
Cr     rho= q_0 phi*phi + 2 q_1 phi*phidot +
Cr          q_2 (phidot*phidot + phi*phidotdot)
C  ---------------------------------------------------
C     implicit none
C Passed parameters
      logical lgdd,lfrz
      integer nl,nlr,lmax,nr,nsp,lrel,kcor,lcor,ipr
      double precision z,a,b,tol,sumev(nsp),sumec(nsp),sumtc(nsp),ec(1),
     .  ev(1),qcor(2),v(nr,nsp),rofi(nr),rho(nr,nlr,nsp),rhoc(nr,nsp),
     .  qnu(3,nl,nsp),pnu(nl,nsp)
C Local parameters
      logical free
      integer konfig(0:10),l,isp,ir,ival,nn,nre,jr,lmaxc,lr,k,nrmx
      parameter (nrmx=1501)
      double precision rocrit,pi,eb1,eb2,q0,q1,q2,rmax,eval,dl,val(5),
     .  slo(5),sum,ro,phi,dphi,phip,dphip,p,fllp1,r,tmc,gfac,c,q00
      double precision g(2*nrmx),gp(2*nrmx*4)
C     Speed of light, or infinity in nonrelativistic case
      common /cc/ c

      call setcc(lrel)

      if (nr .gt. nrmx) call rxi(' newrho: increase nrx, need',nr)
      lr = 1
      rocrit = 0.002d0/4
      pi = 4d0*datan(1d0)
      eb1 = -50d0
      eb2 =  50d0
      rmax = rofi(nr)
      free = (rmax .gt. 9.99d0)
      call config(pnu,lmax,z,konfig,lmaxc)
      if (kcor .gt. 0) then
        lmaxc = max(lmaxc,lcor)
        konfig(lcor) = max(konfig(lcor),kcor+1)
      endif

C --- Calculate core density ---
      if (nlr .eq. 1) then
        if (.not. lfrz) then
          call dpzero(rhoc,nr*nsp)
          call rhocor(0,z,lmaxc,nsp,konfig,a,b,nr,rofi,v,g,
     .                kcor,lcor,qcor,tol,ec,sumec,sumtc,rhoc,0d0,ipr)
        endif
        call dcopy(nr*nsp,rhoc,1,rho,1)
      endif

C --- Loop over valence states ---
      ival = 0
      do  20  isp = 1, nsp
        sumev(isp) = 0d0
        do  20  l = 0, lmax

        if (nlr .eq. nl) lr = l+1
        q0 = max(qnu(1,l+1,isp),0d0)
        q1 = qnu(2,l+1,isp)
        q2 = qnu(3,l+1,isp)
        if (q0 .lt. 1d-6) goto 20
        nn = int(pnu(l+1,isp)) - l - 1
        ival = ival+1
        eval = ev(ival)

        val(1) = rmax
        dl = dtan(pi*(0.5d0 - pnu(l+1,isp)))
        slo(1) = dl + 1

        if (free) val(1) = 1d-30
        if (free) slo(1) = -val(1)
        call rseq(eb1,eb2,eval,tol,z,l,nn,val,slo,v(1,isp),g,
     .            sum,a,b,rofi,nr,nre)
        ev(ival) = eval
        sumev(isp) = sumev(isp) + eval*q0 + q1
        ro = g(nr)**2
        if (.not.free .and. ro .lt. rocrit) write(*,766) l,nn,nre,ro
  766   format(' NEWRHO (warning): PHP,PHPP set to zero,l,nn,nre,rho=',
     .    3i5,2f8.4)
        if (free .or. ro .lt. rocrit) then
          call dpzero(gp,8*nr)
          p = 0
        else
          val(1) = val(1)/dsqrt(sum)
          slo(1) = slo(1)/dsqrt(sum)
C          call phidot(z,l,v(1,isp),eval,a,b,rofi,nr,g,val,slo,tol,
C     .                nn,gp,phi,dphi,phip,dphip,p)
          call phidx(1,z,l,v(1,isp),0d0,0d0,rofi,nr,2,tol,eval,val,
     .      slo,nn,g,gp,phi,dphi,phip,dphip,p,0d0,0d0,0d0,0d0)
        endif
        fllp1 = l*(l+1)
C  ...  Case add q2 phi phidd rho
        if (lgdd) then
        k = 2*nr
        do  21  ir = 2, nre
          jr = ir + nr
          r = rofi(ir)
          tmc = c - (v(ir,isp) - 2d0*z/r - eval)/c
          gfac = 1d0 + fllp1/(tmc*r)**2
          rho(ir,lr,isp) =  rho(ir,lr,isp) +
     .      q0*(gfac*g(ir)**2 + g(jr)**2) +
     .    2*q1*(gfac*g(ir)*gp(ir) + g(jr)*gp(jr)) +
     .      q2*(gfac*(gp(ir)**2 + g(ir)*gp(ir+k)) +
     .          gp(jr)**2 + g(jr)*gp(jr+k))
  21    continue
C  ...  Case add -p q2 phi phi into rho
        else
          q00 = q0-p*q2
          do  22  ir = 2, nre
            jr = ir + nr
            r = rofi(ir)
            tmc = c - (v(ir,isp) - 2d0*z/r - eval)/c
            gfac = 1d0 + fllp1/(tmc*r)**2
            rho(ir,lr,isp) = rho(ir,lr,isp) +
     .        q00*(gfac*g(ir)**2 + g(jr)**2) +
     .      2*q1*(gfac*g(ir)*gp(ir) + g(jr)*gp(jr)) +
     .        q2*(gfac*gp(ir)**2 + gp(jr)**2)
   22     continue
        endif
  20  continue
Cgetarg       end
       end subroutine newrho 

      subroutine rhocor(isw,z,lmax,nsp,konfig,a,b,nr,rofi,v,g,kcor,lcor,
     .                  qcor,tol,ec,sumec,sumtc,rho,gcore,ipr)

C- Generates the (spherical) charge density from the core states
C ----------------------------------------------------------------------
Ci Inputs:
Ci   isw   :1 return the core wave functions in gcore
Ci   z     :nuclear charge
Ci   lmax  :maximum l
Ci   nsp   :=1 for non-spin-polarized =2 for spin-polarized calculations
Ci   konfig:core configuration. Core orbitals are specified by:
Ci         :  1, 2, ..., konf(0)-1 for s            
Ci         :  2, 3, ..., konf(1)-1 for p            
Ci         :  3, 4, ..., konf(2)-1 for d, and so on.
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nr    :number of radial mesh points
Ci   rofi  :radial mesh points
Ci   v     :spherical potential (electronic contribution)
Ci   g     :work array holding normalized wave function times r
Ci   kcor  :(partial core occupation) p.q.n for occupation
Ci   lcor  :(partial core occupation) l quantum for occupation
Ci   qcor  :(partial core occupation) core charge and moment
Ci   tol   :wave function tolerance
Ci   ipr   :0 no printout
Ci         :1 summary printout
Ci         :2 detailed printout
Cio Inputs/Outputs:
Cio  ec    :guessed core eigenvalues (input)
Cio         core eigenvalues (output)
Co Outputs:
Co   rho   :spherical charge density times 4*pi*r*r
Co         :the spherical charge density of the core states is added
Co   sumec :sum of core eigenvalues
Co   sumtc :core kinetic energy
Co   gcore :(isw=0) not used
Co         :(isw=1) core wave functions
Cl Local variables
Cl   deg   :orbital occupation number
Cr Remarks:
Cu Updates
Cu   19 Apr 01 core wave functions may be saved in gcore
C  ---------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer isw,ipr,lmax,konfig(0:lmax),nr,nsp,kcor,lcor
      double precision a,b,z,g(nr,2),rho(nr,nsp),rofi(nr),tol,qcor(2),
     .  sumec(nsp),sumtc(nsp),v(nr,nsp),gcore(nr,2,*),ec(*)
C ... Local parameters
      character*1 pqn(9),ang(6)
      integer icore,isp,konf,l,ll,nodes,nre,stdo,lgunit
      double precision deg,dlml,e1,e2,ecor0,ecore,qcore,rhorim,rmax,
     .  rorim,slo,sum,tcore,val,vrho,kappa2
C ... Heap
      data pqn /'1','2','3','4','5','6','7','8','9'/
      data ang /'s','p','d','f','g','h'/

      stdo = lgunit(1)
      rmax  = rofi(nr)
      b     = rmax/(dexp(a*nr-a)-1.d0)
      e1    = -2.5d0*z*z - 5d0
      e2    = 20d0
      icore = 0
      do  80  isp = 1, nsp
        sumec(isp) = 0d0
        sumtc(isp) = 0d0
        qcore = 0d0
        rhorim = 0d0
        if (ipr .ge. 2) write(stdo,757)
        if(ipr.ge.2 .and. nsp.eq.2)
     .    write(stdo,'('' spin'',i2,'':'')') isp
  757 format(/' state  chg          ecor0',10x,'ecore',10x,'tcore',
     .     4x,'nre',2x,'rho(rmax)')
        do  10  l = 0, lmax
        do  10  konf = l, konfig(l)-2
          deg = (2*(2*l+1))/nsp
          if (konf+1.eq.kcor .and. l.eq.lcor) then
            deg = deg + (qcor(1)+(3-2*isp)*qcor(2))/nsp
          endif
          icore = icore+1
          nodes = konf - l
          ecor0 = ec(icore)
          val = 1d-30
          slo = -val
          call rseq(e1,e2,ecor0,tol,z,l,nodes,val,slo,v(1,isp),
     .      g,sum,a,b,rofi,nr,nre)
          ecore = ecor0
C     ... Correct core energy by using hankel bc's
          kappa2 = ecor0 - v(nr,isp) + 2*z/rmax
          if (nre .eq. nr .and. kappa2 .lt. 0d0) then
            dlml = -1d0 - dsqrt(-kappa2)*rmax
            do  31  ll = 1, l
   31       dlml = -kappa2*rmax*rmax/dlml - (2*ll+1)
            slo = val*(dlml+l+1)/rmax
            call rseq(e1,e2,ecore,tol,z,l,nodes,val,slo,v(1,isp),
     .        g,sum,a,b,rofi,nr,nre)
          endif
          ec(icore) = ecore
          if (isw .eq. 1) call dcopy(2*nr,g,1,gcore(1,1,icore),1)

C     --- Add to rho, make integral v*rho ---
          call xyrhsr(ecore,l,z,a,b,nr,nre,g,rofi,v(1,isp),rho(1,isp),
     .      deg,vrho,rorim)
          rhorim = rhorim + rorim
          tcore = ecore - vrho
          qcore = qcore + deg
          sumec(isp) = sumec(isp) + deg*ecore
          sumtc(isp) = sumtc(isp) + deg*tcore
          if (ipr.ge.2) write(stdo,758)
     .      pqn(konf+1),ang(l+1),deg,ecor0,ecore,tcore,nre,rorim
  758     format(1x,2a1,f8.2,3f15.6,i7,f9.5)
   10   continue
        if (ipr.gt.0) write(stdo,230) qcore,sumec(isp),sumtc(isp),rhorim
  230   format(' sum q=',f5.2,'  sum ec=',f12.5,'  sum tc=',f12.5,
     .    '  rho(rmax)',f8.5)
   80 continue
Cgetarg       end
       end subroutine rhocor 


      subroutine xyrhsr(ecore,l,z,a,b,nr,nre,g,rofi,v,rho,deg,vrho,
     .  rhormx)

C- Make density and integrate potential*density for one core state
C ----------------------------------------------------------------------
Ci Inputs
Ci   ecore :core eigenvalue
Ci   l     :l quantum number
Ci   z     :nuclear charge
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nr    :number of radial mesh points
Ci   nre   :Make density to the smaller of nr and nre
Ci   g     :normalized wave function times r
Ci   rofi  :radial mesh points
Ci   v     :electronic part of spherical potential
Ci   deg   :occupation number
Co Outputs
Co   rho   :contribution core density from this state from 1..min(nr,nre)
Co   vrho  :integral
Co   rhormx:contribution to true density at nr from this state
Cr Remarks
Cr   xyrhsr makes the density at points 1..min(nr,nre), and the integral
Cr   of rho*density from 1..nre.  Thus:
Cr     if nre .eq. nr, the density and integral are the same
Cr     if nre .lt. nr  g**2 is negligible at nr (deep core state)
Cr     if nre .gt. nr  (large sphere) contribution to rho*v is included
Cu Updates
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nr,nre,l
      double precision a,b,z,ecore,g(nr,2),rofi(1),v(nr),rho(nr),rhormx
C ... Local parameters
      integer nrmx,ir
      double precision fpi,fllp1,vrho,r,wgt,tmc,c,gfac,rhoir,deg,rmax
C     Speed of light, or infinity in nonrelativistic case
      common /cc/ c

      fpi = 16d0*datan(1d0)
      fllp1 = l*(l+1)
      vrho = 0
      nrmx = min0(nr,nre)
C ... Make rho, and integrate vrho for points 1..nrmx
      do  11  ir = 2, nrmx
        r = rofi(ir)
        wgt = 2*(mod(ir+1,2)+1)
        if (ir .eq. nre) wgt = 1
        tmc = c - (v(ir) - 2*z/r - ecore)/c
        gfac = 1 + fllp1/(tmc*r)**2
        rhoir = gfac*g(ir,1)**2 + g(ir,2)**2
        vrho = vrho + wgt*rhoir * (v(ir)-2*z/r) * (r+b)
        rho(ir) = rho(ir) + deg*rhoir
   11 continue
      rmax = rofi(nr)
      rhormx = 0
C     nfp has the following line:
C     if (nre .ge. nr) rhormx = deg*g(nr,1)**2/(rmax*rmax*fpi)
      if (nre .ge. nr) rhormx = deg*rhoir/(fpi*rmax**2)

C ... Integrate rho*v from nrmx+1 .. nre
      do  12  ir = nrmx+1, nre
        r = rofi(ir)
        wgt = 2*(mod(ir+1,2)+1)
        if (ir .eq. nre) wgt = 1
        tmc = c - (v(ir) - 2*z/r - ecore)/c
        gfac = 1 + fllp1/(tmc*r)**2
        rhoir = gfac*g(ir,1)**2 + g(ir,2)**2
        vrho = vrho + wgt*rhoir * (v(ir)-2*z/r) * (r+b)
   12 continue
      vrho = vrho*a/3

Cgetarg       end
       end subroutine xyrhsr 


      subroutine potpar(nl,nsp,lmx,z,rmax,avw,ekap,lso,loptc,lmpol,a,nr,
     .  rofi,v,pnu,idmod,ves,exc,qnu,idu,uh,jh,thrpv,thrpvl,g,gp,pp,
     .  pprel,sop,pmpol,gradrm)

C- Generates potential parameters for given potential
C  ---------------------------------------------------
Ci Inputs:
Ci   nl    :(global maximum l) + 1
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   lmx   :lmx(j) = maximum l for atom j
Ci   z     :nuclear charge
Ci   rmax  :augmentation radius, in a.u.
Ci   avw   :length scale, usu. average Wigner-Seitz sphere radius
Ci   ekap  :muffin-tin zero plus kap2
Ci   lso   :if true, make spin-orbit coupling parms sop
Ci   loptc :if true, make matrix elements <phi grad phi> etc
Ci   lmpol :if true, make multipole moments of phi,phidot
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nr    :number of radial mesh points
Ci   rofi  :radial mesh points
Ci   v     :spherical potential
Ci   pnu   :Determines b.c. for particular l channel; see remarks
Ci   idmod :0 or 1 : linearization energy enu for l channel determined
Ci         :         by input log derivative spec'd by pnu
Ci         :2      : linearization energy enu for l channel determined
Ci         :         from pp(1)-ves; pnu is output; see remarks
Ci         :See also remarks
Ci   ves   :used only when idmod=2 (see remarks)
Ci   exc   :XC energy used to calc. thrpv
Ci   qnu   :energy-weighted moments of the sphere charges
Ci         :used only in the calculation of 3pV
Ci   idu   :idu(l+1)=0 => this l has no nonlocal U matrix
Ci                   1..3 this l has no nonlocal U matrix
Ci                   4    majority and minority C,E shifted by U; see remarks
Ci                   5    majority and minority C,E shifted by U and J; see remarks
Ci         :Note: this routine doesn't implement LDA+U; see remarks
Ci   uh    :U in LDA+U
Ci   jh    :J in LDA+U
Co Outputs:
Co   pp:   :enu (see remarks), c, srdel, p, gamma (aka q), alpha=gamma
Co         :c = center of band = enu + omega-
Co         :sqrt(delta) has proper sign (that of phi-).
Co         :pp(1) enu
Co         :pp(2) calpha
Co         :pp(3) srdel = sqrt(delta) with proper sign (that of phi-).
Co         :pp(4) palpha
Co         :pp(5) gamma, or Q in Varenna notes
Co         :pp(6) alpha, or qbar in Varenna notes
Co   pprel :potential parameters from the Dirac equation
Co         :pprel(1,l,mu,isp,jsp) cgamma(l,mu,isp,jsp) ?
Co         :pprel(2,l,mu,isp,jsp) gamma(l,mu,isp,jsp)
Co         :pprel(3,l,mu,isp,jsp) delta(l,mu,isp,jsp)
Co         :pprel(4,l,mu,isp,jsp) pgamma(l,mu,isp,jsp)
Co   thrpv : 3 PV for sphere
Co   thrpvl: 3 PV for sphere, decomposed by l-channel and spin
Co   g     : wave function; see rseq
Co   gp    : phidot, except for ...
Co   sop   : spin-orbit coupling parameters
Co   pmpol :integral (phi-or-phidot * phi-or-phidot * r**l) :
Co         :matrix elements of w.f. * wf * r**l for multipole moments
Co   gradrm:
Cr Remarks:
Cr   pnu = .5 - atan(dnu)/pi + (princ.quant.number).
Cr   pp's are generated for a specified potential.
Cr   The enu about which the pp's are generated is calculated from
Cr   the input pnu.  Alternatively, by setting idmod=2 for a specified
Cr   l channel the enu is not calculated but taken from pp(1)-ves, in
Cr   which case pnu is not used.
Cr
Cr   Potential parameters generated from (ekap=0):
Cr     omega(-) = -(phi/phidot) (-l-1-dnu)/(-l-1-dnudot)
Cr     omega(+) = -(phi/phidot) (l-dnu)/(l-dnudot)
Cr     phi(+) = phi + omega(+) phidot
Cr     phi(-) = phi + omega(-) phidot
Cr     C  = e + omega(-)
Cr     Delta = phi(-)**2 sdivw rmax/2
Cr     Gamma = (phi(-)/phi(+)) sdivw / (2 (2l+1))
Cr     where sdivw = (rmax/avw)**(l+l+1)
Cr   These conventions differ by those of MSM by the scale factor
Cr   sdivw for delta and gamma.
Cr   For general ekap (see Kanpur notes)
Cr     C = e - W{K,phi} / W{K,phidot}
Cr       = e - phi/phidot * (D{phi}-D{K})/(D{phidot}-D{K})
Cr     sqrt(delta) = - sqrt(avw/2) / W{K,phidot}
Cr                 = - sqrt(avw/2) /wsr / K_l / phidot /(D{phidot}-D{K})
Cr     gamma = W{J,phidot} / W{K,phidot}
Cr           = J_l/K_l * (D{phidot}-D{J})/(D{phidot}-D{K})
Cr   Uses W{a,b} = r*r*(ab' - a'b) = r*a*b*[D{b}-D{a}]
Cr
Cr   LDA+U parameters:  No LDA+U is implemented in this code.
Cr   For now, parameters are just a device to shift pot pars by a constant.
Cr   You can shift (enu,C) of both spins by U (idu=4)
Cr        or shift (enu,C) of spin1 by U, spin2 by J (idu=5)
Cl Local variables
Cl   lpzi  :flags how local orbitals is to be treated in current channel
Cl         :0 no local orbital gz (none implemented for now)
Cr Updates
Cu  21 Dec 05 (wrl) potential shifts to mimic LDA+U
Cu  11 Jul 05 Updated soprm call to render compatible with fp
Cu  18 Jun 04 (A Chantis) relativistic potential parameters
Cu   4 Apr 04 New matrix elements of <phi Bxc phi>
Cu  05 Jan 04 bug fix making SO matrix elements when nl ne lmx+1
Cu  07 Feb 03 When calc. SO matrix elements, also make <phi|phi>
Cu            for inclusion of external magnetic field.
Cu  20 Apr 98 Optical matrix elements adapted from Sergey Rashkeev
Cu            MvS spin-orbit matrix elements adapted from V. Antropov
C  ---------------------------------------------------
C     implicit none
C Passed parameters
      logical lso,lmpol,loptc
      integer nl,lmx,nr,nsp,idmod(0:nl-1),idu(4)
      double precision uh(4),jh(4)
      double precision z,ekap,rmax,avw,a,ves,thrpv,thrpvl(0:nl-1,nsp)
      double precision rofi(1),v(nr,nsp),g(1),gp(nr,4),pp(6,0:nl-1,nsp),
     .                 pnu(0:nl-1,nsp),qnu(3,0:nl-1,nsp),exc(2),
     .                 sop(0:nl-1,nsp,nsp,9),pmpol(nl,nl,2*nl-1,3,nsp),
     .                 gradrm(4,2,nl,2,nsp),pprel(4,0:nl-1,2*nl,2,2)
C Local parameters
      integer i,j,l,nn,nre,n0,stdo
      parameter (n0=10)
      integer lpzi(0:n0)
      double precision eb1,eb2,pi,b,e,val(5),slo(5),sum,sdivw,
     .  phi,dphi,phip,dphip,dlphi,dlphip,phplus,phmins,
     .  el(2),ql(2),e1me2,e1pe2,e1e2,d,dlmins,dlplus,
     .  fi(0:9),gi(0:9),ptmp1(2,2),ptmp2(2,2),ptmp3(2,2),ptmp4(2,2),
     .  ptmp5(2,2),ptmp6(2,2),ptmp7(2,2),ptmp8(2,2)
C     double precision omegam,omegap
      double precision wk,wkdot,wjdot,scl,enu(0:8,2),fac1,fac2
      integer ophi,odphi,owk,odv,owgt

C Heap allocation
      integer w(1)
      common /w/ w
      double precision tol
      parameter (tol=1d-12)
      external rsq1,rseq,phidx

C ... For Dirac equation
      integer mumax,imu,lrel,nglob
      double precision measur,phisq(0:lmx,nsp),
     .  gmt(2,2),fmt(2,2),gmtde(2,2),fmtde(2,2),gsmt(2,2)
      double precision srd(0:lmx,2), srddot(0:lmx,2)
      double precision frd(0:lmx,2), frddot(0:lmx,2)
      double precision srg(0:lmx,2), srgdot(0:lmx,2)
      double precision frg(0:lmx,2), frgdot(0:lmx,2)

      stdo = nglob('stdo')
      lrel = nglob('lrel')
      do  l = 0, min(lmx,3)
C       LDA+U-like shifts in C
        if (idu(l+1) .le. 3) then
        elseif (idu(l+1) .le. 4) then
          call info2(20,0,0,'%10fpotpar l=%i: enu and C shifted by %d',
     .              l,uh(l+1))
C          if (idu(l+1) .le. 3)
C     .      call rx('POTPAR: LDA+U not implemented in ASA')
        elseif (idu(l+1) .eq. 5) then
          call info5(20,0,0,'%10fpotpar l=%i: enu+ and C+ shifted by'//
     .      ' %d;  enu- and C- shifted by %d',l,uh(l+1),jh(l+1),0,0)
        endif
      enddo

      if (lso .or. lmpol .or. loptc) then
        call defdr(ophi,nr*nl*nsp)
        call defdr(odphi,nr*nl*nsp)
        call defdr(owk,nr*4)
        call defdr(odv,nr)
        call defdr(owgt,nr)
        call radwgt(rmax,a,nr,w(owgt))
        call iinit(lpzi,n0+1)
      endif

      thrpv = 0
      eb1 = -20d0
      eb2 =  20d0
      pi = 4d0*datan(1d0)
      b = rmax/(dexp(a*nr - a) - 1d0)
      call bessl2(ekap*rmax**2,0,lmx+1,fi,gi)

      do  80  i = 1, nsp
        do  10  l = 0, lmx
          if (mod(idmod(l),10) .eq. 2) then
            e = pp(1,l,i) - ves
            call rsq1(0,e,l,z,v(1,i),nr,g,val,slo,nn,a,b,rofi,nr)
          else
            e = -0.5d0
            nn = int(pnu(l,i))-l-1
            val(1) = rmax
            slo(1) = 1 + dtan(pi*(0.5d0 - pnu(l,i)))
          endif
          call rseq(eb1,eb2,e,tol,z,l,nn,val,slo,v(1,i),
     .              g,sum,a,b,rofi,nr,nre)
          val(1) = val(1)/dsqrt(sum)
          slo(1) = slo(1)/dsqrt(sum)
          call phidx(1,z,l,v(1,i),0d0,0d0,rofi,nr,2,tol,e,val,slo,nn,
     .      g,gp,phi,dphi,phip,dphip,pp(4,l,i),0d0,0d0,0d0,0d0)
          phisq(l,i) = phi*phi
C     ... Keep local copies of phi and phidot for SO coupling
          if (lso .or. lmpol .or. loptc) then
            call dpscop(g,w(ophi),nr,1,1+nr*(l+nl*(i-1)),1d0)
            call dpscop(gp,w(odphi),nr,1,1+nr*(l+nl*(i-1)),1d0)
            enu(l,i) = e
          endif
          dlphi =  rmax*dphi/phi
          dlphip = rmax*dphip/phip
          sdivw = (rmax/avw)**(l+l+1)
C     ... Following valid for ekap=0 only:
C         omegam = -(phi/phip)*(-l-1-dlphi)/(-l-1-dlphip)
C         omegap = -(phi/phip)*(l-dlphi)/(l-dlphip)
C         phplus = phi + omegap*phip
C         phmins = phi + omegam*phip
C         pp(1,l,i) = e
C         pp(2,l,i) = e + omegam
C         pp(3,l,i) = phmins*dsqrt(sdivw*rmax/2)
C         pp(5,l,i) = phmins/(2*(2*l+1)*phplus)*sdivw
C         pp(6,l,i) = pp(5,l,i)
C     ... The following for general ekap:
C         scl should be 1.  Scale srdel to make compatible w/ ekap=0
          wk     = (dlphi-l)*gi(l)  + (l+l+1)*gi(l+1)
          wkdot  = (dlphip-l)*gi(l) + (l+l+1)*gi(l+1)
          wjdot  = (dlphip-l)*fi(l) + ekap*rmax**2*fi(l+1)/(l+l+1)
C         omegam = -(phi/phip)*wk/wkdot
          scl    = phi*phip*(dlphi-dlphip)*rmax
C     ... Force enu to stay really fixed for idmod=2
          if (mod(idmod(l),10) .ne. 2) then
            pp(1,l,i) = e
          else
            pp(1,l,i) = pp(1,l,i) - ves
          endif
          pp(2,l,i) = e - (phi/phip)*wk/wkdot
          pp(3,l,i) = -dsqrt(sdivw/2/rmax)/phip/wkdot*scl
          pp(5,l,i) = sdivw*wjdot/wkdot
          pp(6,l,i) = pp(5,l,i)
C     ... Constant shift for potential parameters, mimicking LDA+U
          if (idu(l+1) .eq. 4) then
             pp(1,l,i) = pp(1,l,i) + uh(l+1)
             pp(2,l,i) = pp(2,l,i) + uh(l+1)
          endif
C     ... Spin-dependent constant shift for enu,C, mimicking LDA+U
          if (idu(l+1) .eq. 5) then
            if (i .eq. 1) then
              pp(1,l,i) = pp(1,l,i) + uh(l+1)
              pp(2,l,i) = pp(2,l,i) + uh(l+1)
            elseif (i .eq. 2) then
              pp(1,l,i) = pp(1,l,i) + jh(l+1)
              pp(2,l,i) = pp(2,l,i) + jh(l+1)
            endif
          endif
C --- Calculate 3PV ala MSM J. Phys. F 12, 141, Eqn 2.24 ---
C NB: formula valid only for ekap=0
          thrpvl(l,i) = 0
          if (qnu(1,l,i).le.0 .or. qnu(3,l,i).le.0 .or. exc(i) .eq. 0)
     .      goto 20
          d     = qnu(2,l,i)**2 - qnu(1,l,i)*qnu(3,l,i)
          e1pe2 = qnu(2,l,i)*qnu(3,l,i)/d
          e1e2  = qnu(3,l,i)**2/d
C The following two apply if have also third moment ...
C         e1pe2 = (qnu(2,l,i)*qnu(3,l,i) - qnu(1,l,i)*qnu(4,l,i))/d
C         e1e2  = (qnu(3,l,i)**2 - qnu(2,l,i)*qnu(4,l,i))/d
          e1me2 = dsqrt(max(e1pe2**2 - 4*e1e2,0d0))
          if (e1me2 .ne. 0) then
            el(1) = (e1pe2 + e1me2)/2
            el(2) = (e1pe2 - e1me2)/2
            ql(1) = ( qnu(2,l,i) - el(2)*qnu(1,l,i))/e1me2
            ql(2) = (-qnu(2,l,i) + el(1)*qnu(1,l,i))/e1me2
          else
            el(1) = e1pe2
            el(2) = e1pe2
            ql(1) = 0
            ql(2) = 0
          endif

C First order approximation to phi, dphi
          phplus = phi + el(1)*phip
          phmins = phi + el(2)*phip
          dlplus = dlphi - el(1)/phi**2/rmax
          dlmins = dlphi - el(2)/phi**2/rmax

#ifndef FIRST_ORDER_PRESSURE
C Get phi, dphi by reintegrating wave function
          call rsq1(0,el(1)+e,l,z,v(1,i),nr,g,val,slo,nn,a,b,rofi,nr)
          dlplus = rmax*(slo(1)-val(1)/rmax)/val(1)
          call gintsr(g,g,a,b,nr,z,el(1)+e,l,v(1,i),rofi,sum)
          phplus = val(1)/rmax/dsqrt(sum)
          call rsq1(0,el(2)+e,l,z,v(1,i),nr,g,val,slo,nn,a,b,rofi,nr)
          dlmins = rmax*(slo(1)-val(1)/rmax)/val(1)
          call gintsr(g,g,a,b,nr,z,el(2)+e,l,v(1,i),rofi,sum)
          phmins = val(1)/rmax/dsqrt(sum)
#endif

          thrpvl(l,i) =
     .      rmax*ql(1)*phplus**2*
     .      (dlplus*(dlplus+1) - l*(l+1) + rmax**2*(el(1)+e-exc(i))) +
     .      rmax*ql(2)*phmins**2*
     .      (dlmins*(dlmins+1) - l*(l+1) + rmax**2*(el(2)+e-exc(i)))
          thrpv = thrpv + thrpvl(l,i)

   20     continue


   10   continue
   80 continue

C --- Make spin-orbit matrix elements of phi, phidot ---
      if (lso) then
C       call pshpr(51)

C   ... Gradient of average v
        call dpcopy(v,w(owk),1,nr,1d0/nsp)
        if (nsp .eq. 2) call dpadd(w(owk),v(1,2),1,nr,.5d0)
        call radgra(a,b,nr,rofi,w(owk),w(odv))

C   ... Calculate s-o parameters
        call soprm(5,lpzi,w(ophi),w(odphi),w,nr,nsp,nl-1,lmx,v,w(odv),
     .    enu,w,z,rofi,w(owgt),w(owk),sop,w)

C --- Fully relativistic potential parameters ---
      if (lrel .eq. 2) then
        do  l = 0, lmx
          mumax= 2*(l+1)
          do  imu = 1, mumax
C           mu = dble(imu-l) - 1.5d0
            measur = rmax*0.5d0*(phisq(l,1) + phisq(l,2))
            call rdeq(pp(1,l,1),pp(1,l,2),sop,z,v,rofi,nr,2,a,b,l,lmx,
     .        imu,measur,gmt,fmt,gmtde,fmtde,gsmt,pprel)

           call fdpp(pp(1,l,1),pp(1,l,2),sop,ves,gmt,fmt,gmtde,fmtde,z,
     .        rmax,avw,l,lmx,imu,srd(l,1),srd(l,2),frd(l,1),frd(l,2),
     .        srddot(l,1),srddot(l,2),frddot(l,1),frddot(l,2),
     .        srg(l,1),srg(l,2),frg(l,1),frg(l,2),srgdot(l,1),
     .        srgdot(l,2),frgdot(l,1),frgdot(l,2),pprel,gsmt)

C      ... Test: Give scalar Relativistic small parameter
           if(z .ne. 0) then
             do i = 1, 2
               do j = 1, 2
                 ptmp4(i,j) = 0d0
               enddo
               ptmp4(i,i) = pp(4,l,i)
             enddo
             call clebsh(l,imu,ptmp4,ptmp8)
             do i = 1, 2
               do j = 1, 2
                 pprel(4,l,imu,i,j) = ptmp8(i,j)
               enddo
             enddo
           endif

           if (z .eq. 0) then
             do i = 1, 2
               do j = 1, 2
                 ptmp1(i,j) = 0d0
                 ptmp2(i,j) = 0d0
                 ptmp3(i,j) = 0d0
                 ptmp4(i,j) = 0d0
               enddo
               ptmp1(i,i) = pp(2,l,i)+ves
               ptmp2(i,i) = pp(5,l,i)
               ptmp3(i,i) = pp(3,l,i)
               ptmp4(i,i) = pp(4,l,i)
             enddo
             call clebsh(l,imu,ptmp1,ptmp5)
             call clebsh(l,imu,ptmp2,ptmp6)
             call clebsh(l,imu,ptmp3,ptmp7)
             call clebsh(l,imu,ptmp4,ptmp8)
             do i = 1, 2
               do j = 1, 2
                 pprel(1,l,imu,i,j) = ptmp5(i,j)
                 pprel(2,l,imu,i,j) = ptmp6(i,j)
                 pprel(3,l,imu,i,j) = ptmp7(i,j)
                 pprel(4,l,imu,i,j) = ptmp8(i,j)
C                p2 for the empty spheres
C                pprel(4,l,imu,i,j) = 0.d0
               enddo
             enddo
           endif

          enddo
        enddo
      endif


C   ... Matrix elements for constant magnetic field
        call dvset(w(odv),1,nr,1d0)
        call soprm(2,lpzi,w(ophi),w(odphi),w,nr,nsp,nl-1,lmx,v,w(odv),
     .    enu,w,z,rofi,w(owgt),w(owk),sop(0,1,1,4),w)
C       Correct <dot||dot> term from (already known) diagonal part
C       <phi|phi> should be unity; <phi|dot> should be zero
        do  110  l = 0, lmx
          fac1 = pp(4,l,1)/sop(l,1,1,6)
          fac2 = pp(4,l,2)/sop(l,2,2,6)
          sop(l,1,1,6) = sop(l,1,1,6)*sqrt(fac1*fac1)
          sop(l,1,2,6) = sop(l,1,2,6)*sqrt(fac1*fac2)
          sop(l,2,1,6) = sop(l,2,1,6)*sqrt(fac2*fac1)
          sop(l,2,2,6) = sop(l,2,2,6)*sqrt(fac2*fac2)
  110   continue

C   ... Matrix elements of XC field
        call dpcopy(v,w(odv),1,nr,0.5d0)
        call dpadd(w(odv),v(1,nsp),1,nr,-0.5d0)
C       call prrmsh('vxc ',rofi,w(odv),nr,nr,1)
        call soprm(2,lpzi,w(ophi),w(odphi),w,nr,nsp,nl-1,lmx,v,w(odv),
     .    enu,w,z,rofi,w(owgt),w(owk),sop(0,1,1,7),w)

C       debugging: generate M<B>
C        print *, '!!'
C        if (nsp .eq. 2) then
C          wk = 0
C          do  l = 0, lmx
C            ql(1) = qnu(1,l,1) - pp(4,l,1)*qnu(3,l,1)
C            ql(2) = qnu(1,l,2) - pp(4,l,2)*qnu(3,l,2)
C            fac1 =
C     .        ql(1)*sop(l,1,1,7) +
C     .        qnu(2,l,1)*sop(l,1,1,8) +
C     .        qnu(3,l,1)*sop(l,1,1,9)
C            fac2 =
C     .        ql(2)*sop(l,2,2,7) +
C     .        qnu(2,l,2)*sop(l,2,2,8) +
C     .        qnu(3,l,2)*sop(l,2,2,9)
C            wk = wk + fac1-fac2
C            write(stdo,333) l, fac1,fac2,fac1-fac2
C  333       format(' l=',i1,'   q+<B>=',f12.6,'   q-<B>=',f12.6,
C     .             '   M<B>=',f12.6)
C          enddo
C          write(stdo,'('' <Bxc*M>'',f12.6)') wk
C        endif
      endif

C --- Matrix elements of wave function and its gradient ---
      if (loptc) then
        call rgrme(nsp,nl-1,lmx,nr,rofi,w(owgt),w(ophi),w(odphi),gradrm)
      endif

C --- Multipole moments of phi, phidot ---
      if (lmpol) then
C   ... Orthonormalize phi,phidot
        call soprm(4,lpzi,w(ophi),w(odphi),w,nr,nsp,nl-1,lmx,v,w,enu,w,
     .    z,rofi,w(owgt),w(owk),sop,w)
        call mpint(w(ophi),w(odphi),nl-1,lmx,2*lmx,nr,nsp,
     .    rofi,rofi(1+nr),pmpol)
      endif

      if (lso .or. lmpol) call rlse(ophi)
Cgetarg       end
       end subroutine potpar 

C      subroutine phidot(z,l,v,e,a,b,rofi,nr,g,val,slo,tol,nn,
C     .  gp,phi,dphi,phip,dphip,p)
CC- Generate Phidot,Phidotdot for a prescribed energy
CC ----------------------------------------------------------------
CCi Inputs:
CCi   z:     nuclear charge
CCi   l:     l quantum number for this g
CCi   v:     potential
CCi   a,b:   defines shifted logarithmic mesh (rmesh.f)
CCi   rofi:  list of points; must be consistent with a,b
CCi   nr:    number of mesh points
CCi   e:     Energy
CCi   val:   Value of r * wave function at sphere radius (rseq)
CCi   slo:   Derivative of r * wave function at sphere radius (rseq)
CCi   g:     Wave function times r normalized so that int (g*g) dr = 1
CCi   tol:   precision to which wave function is integrated
CCi   nn:    number of nodes
CCo Outputs:
CCo   gp:    first four energy derivatives to G
CCo   phi:   wave function at rmax, i.e. g/rmax
CCo   dphi:  slope of wave function at rmax, i.e. d(g/rmax)/dr
CCo   phip:  energy derivative of wave function at rmax
CCo   dphip: energy derivative of slope of wave function at rmax
CCo   p:     <gp**2> (potential parameter)
CCr Remarks:
CCr   This version makes energy derivatives by numerical differentiation
CCr   of wave function phi, and has the same calling sequence as the
CCr   analytic version phidot.  The only difference is that this routine
CCr   returns four derivatives of phi, whereas the analytic version only
CCr   returns two and is applicable only to the nonrelativistic case.
CC ----------------------------------------------------------------
C      implicit none
CC passed parameters
C      integer l,nr,nn
C      double precision z,e,a,b,val,slo,phi,dphi,phip,dphip,p,tol
C      double precision v(1),rofi(1)
C      double precision g(2*nr),gp(2*nr,4)
CC local variables
C      integer nre,i,iprint,nptdif
C      double precision rmax,eb1,eb2,dele,ddde,sum1,
C     .                 vali(5),sloi(5),ei(4),de1,de2,del1,del2
C      parameter (nptdif = 2)
C
C      rmax = rofi(nr)
C      eb1 = -50d0
C      eb2 = 20d0
C
Cc      dele = tol**.2D0
C      dele = .002d0
C      if (tol .gt. 1d-9 .and. iprint() .ge. 0)
C     .  print *, 'phidot: tol too high for reliable num. diff'
C
C      ddde = -rmax/g(nr)**2
C      ei(1) = 1
C      ei(2) = -1
C      ei(3) = 1.5d0
C      ei(4) = -1.5d0
C      do  10  i = 1, nptdif
C        sloi(i) = slo + dele*ei(i)*ddde*val/rmax
C        ei(i) = e + dele*ei(i)
C        call rseq(eb1,eb2,ei(i),tol,z,l,nn,val,sloi(i),v,gp(1,i),
C     .            sum1,a,b,rofi,nr,nre)
C        vali(i) = val/dsqrt(sum1)
C        sloi(i) = sloi(i)/dsqrt(sum1)
C   10 continue
C      de1  = (ei(1) - ei(2))/2
C      del1 = (ei(1) + ei(2))/2 - e
C      de2  = (ei(3) - ei(4))/2
C      del2 = (ei(3) + ei(4))/2 - e
CC     Energy derivatives of value and slope
C      call dfphi(de1,del1,de2,del2,1,val,vali,nptdif.eq.4)
C      call dfphi(de1,del1,de2,del2,1,slo,sloi,nptdif.eq.4)
C
C      call dfphi(de1,del1,de2,del2,2*nr,g,gp,nptdif.eq.4)
C      call gintsr(gp,gp,a,b,nr,z,e,l,v,rofi,p)
C
CC ... Get phi,dphi from val,slo = (r*phi),(r*phi)' at rmax
C      phi = val/rmax
C      dphi = (slo - phi)/rmax
C      phip = vali(1)/rmax
C      dphip = (sloi(1) - phip)/rmax
C
C      if (iprint() .ge. 111) print 749, phi,dphi,phip,dphip
C  749 format(' PHIDOT:  phi,phip,phip,dphip=',4f12.6)
C      end
      subroutine gintsr(g1,g2,a,b,nr,z,e,l,v,rofi,sum)

C- Integrate inner product of two wave equations
C ----------------------------------------------------------------
Ci   g1,g2 :First and second radial wave functions
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nr    :number of radial mesh points
Ci   z     :nuclear charge
Ci   e     :energy
Ci   l     :l quantum number of g1,g2
Ci   v     :spherical potential
Ci   rofi  :radial mesh points
Co Outputs:
Co   sum   :inner product
Cr Remarks:
Cr   Uses Simpson's rule
C ----------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nr,l
      double precision a,b,z,e,sum,g1(nr,2),g2(nr,2),v(nr),rofi(nr)
C ... Local parameters
      integer i,ir
      double precision fllp1,c,r
      double precision tmc,fi
C     Speed of light, or infinity in nonrelativistic case
      common /cc/ c
      tmc(i,r) = c - (v(i) - 2d0*z/r - e)/c
      fi(i,r) = (r+b)*(g1(i,1)*g2(i,1)*(1 + fllp1/(tmc(i,r)*r)**2)
     .  + g1(i,2)*g2(i,2))

      fllp1 = l*(l+1)
      sum = 0d0
      do  10  ir = 2, nr-1, 2
   10 sum = sum + fi(ir,rofi(ir))
      sum = 2*sum
      do  11  ir = 3, nr-2, 2
   11 sum = sum + fi(ir,rofi(ir))
      sum = (2*sum + fi(nr,rofi(nr)))*a/3

Cgetarg       end
       end subroutine gintsr 

      subroutine gintsl(g1,g2,a,b,nr,rofi,sum)

C- Integrate inner product of two wave equations, large component only
C ----------------------------------------------------------------------
Ci Inputs
Ci   g1    :first wave function
Ci   g2    :second wave function
Ci   a     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   b     :the mesh points are given by rofi(i) = b [e^(a(i-1)) -1]
Ci   nr    :number of radial mesh points
Ci   rofi  :radial mesh points
Co Outputs
Co   sum:   inner product
Cr Remarks
Cr   Like gintsr, but uses large component only (corresponds to c->infty)
Cu Updates
Cu   20 Apr 01
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nr
      double precision a,b,g1(nr),g2(nr),rofi(nr)
C ... Local parameters
      integer ir
      double precision sum,r

      sum = 0
      do  10  ir = 2, nr-1,2
      r = rofi(ir)
  10  sum = sum+(r+b)*(g1(ir)*g2(ir))
      sum = sum+sum
      do  11  ir = 3, nr-2, 2
      r = rofi(ir)
  11  sum = sum+(r+b)*(g1(ir)*g2(ir))
      sum = sum+sum
      r = rofi(nr)
      sum = sum+(r+b)*(g1(nr)*g2(nr))
      sum = sum*a/3
Cgetarg       end
       end subroutine gintsl 

      subroutine asprjq(mode,clabl,nl,nsp,eula,neul,pnu,qnu,
     .  pnuloc,qnuloc,bhat,amom)

C- Find average magnetization axis and project ASA moments onto it
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :1s digit affects bhat
Ci         :0 use bhat as input
Ci         :1 Print out mag. dir. from eula and qnu; make amom
Ci         :2 make bhat from eula and qnu; make amom
Ci         :3 combination of 1+2
Ci         :10s digit affects qnuloc,pnuloc
Ci         :0 do nothing to qnuloc,pnuloc
Ci         :1 copy pnu->pnuloc and qnu->qnuloc
Ci         :2 copy and rotate to bhat coordinates
Ci         :10s digit affects B,qnuloc when magnetization < 0
Ci         :1 scale B by -1 when
Ci   mode  :0 do nothing; just return
Ci   clabl :class name (for printout only)
Ci   nl    :(global maximum l) + 1
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   eula  :Euler angles for noncollinear spins
Ci   neul  :1, nl, or nl**2 if Euler angles are: l-independent,
Ci         :l-dependent, or lm-dependent, respectively
Ci   pnu   :boundary conditions.  If Dl = log. deriv. at rmax,
Ci          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Ci   qnu   :energy-weighted moments of the sphere charges
Cio Inputs/Outputs
Cio  bhat  :direction vector for B-field; see 1s digit mode
Co Outputs
Co   amom  :projection of magnetization onto bhat
Co   pnuloc:projection of pnu onto bhat
Co   qnuloc:projection of qnu onto bhat
Cr Remarks
Cr   If asprjq computes the B-field, it is taken to be parallel to
Cr   the average magnetization, which is computed by summing the
Cr   (vector) magnetization over all orbitals.
Cu Updates
Cu   06 Apr 04 First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,nl,nsp,neul
      double precision pnu(nl,nsp),qnu(3,nl,nsp),eula(neul,3),
     .  amom,bhat(3),pnuloc(nl,nsp),qnuloc(3,nl,nsp)
      character clabl*8
C ... Local parameters
      logical lwrite
      integer i,k,l,m,ilm,stdo,lgunit,ipr,PRT1,PRT2,mode0,mode1,mode2
      double precision ql,al,alpha,beta,gamma,hatm(3),ploc,qloc(3),
     .  dploc,aloc(3),dotp,ddot,dsqrt,rotm(3,3),sal
      character strn*5
      parameter (PRT1=40,PRT2=50)

      stdo = lgunit(1)
      call getpr(ipr)

C --- Construct the average magnetization (direction, amplitude) ---
      mode0 = mod(mode,10)
      mode2 = mod(mode/100,10)
      if (mode0 .ne. 0 .and. nsp .eq. 2) then
      lwrite = .true.
      if (mod(mode0,2) .eq. 0 .or. ipr .lt. PRT2) lwrite = .false.
      if (neul .le. 1) lwrite = .false.
      strn = '    l'
      if (neul .eq. nl*nl) strn = '  ilm'
      if (lwrite .and. mode0 .ge. 2) write(stdo,345) strn
      if (lwrite .and. mode0 .lt. 2) write(stdo,345) strn,'mhat.bxc'
  345 format(a,'     ql       mom',9x,'alpha     beta      gamma',
     .  19x,'mhat':17x,a)

      ilm = 0
      call dpzero(aloc,3)
      sal = 0
      do  l = 0, nl-1
        do   m = -l, l
          ilm = ilm+1
          if (neul .eq. nl) then
            alpha = eula(l+1,1)
            beta  = eula(l+1,2)
            gamma = eula(l+1,3)
          elseif (neul .eq. nl*nl) then
            alpha = eula(ilm,1)
            beta  = eula(ilm,2)
            gamma = eula(ilm,3)
          elseif (neul .eq. 1) then
            alpha = eula(1,1)
            beta  = eula(1,2)
            gamma = eula(1,3)
          else
            call rxi('atscpp: bad value neul=',neul)
          endif

C         Charge, magnetic moments, quantization axis for these angles
          ql = qnu(1,l+1,1)+qnu(1,l+1,2)
          al = qnu(1,l+1,1)-qnu(1,l+1,2)
          sal = sal + al/(2*l+1)
C         The local magnetization points along V = zhat(loc).
C         In global coordinates V = rotm^-1 zhat(loc) because
C         rotm*V = zhat(loc) when V points along M as
C         described in eua2rm.  V is then
C         hatm(1) = dcos(alpha)*dsin(beta)
C         hatm(2) = dsin(alpha)*dsin(beta)
C         hatm(3) = dcos(beta)
          call eua2rm(alpha,beta,gamma,rotm)
          hatm(1) = rotm(3,1)
          hatm(2) = rotm(3,2)
          hatm(3) = rotm(3,3)

C         Add to total magnetization
          call daxpy(3,al/(2*l+1),hatm,1,aloc,1)

C         Printout
          if (neul .eq. nl**2 .and. lwrite .and. mode0 .ge. 2) then
            write(stdo,'(i5,2f10.6,3x,3f10.6,3x,3f10.6,3x,f10.6)')
     .        ilm, ql/(2*l+1), al/(2*l+1), alpha, beta, gamma, hatm
          elseif (neul .eq. nl**2 .and. lwrite .and. mode0 .lt. 2) then
            write(stdo,'(i5,2f10.6,3x,3f10.6,3x,3f10.6,3x,f10.6)')
     .        ilm, ql/(2*l+1), al/(2*l+1), alpha, beta, gamma, hatm,
     .        ddot(3,hatm,1,bhat,1)
          elseif (neul .eq. nl .and. lwrite .and. mode0 .ge. 2) then
            write(stdo,'(i5,2f10.6,3x,3f10.6,3x,3f10.6,3x,f10.6)')
     .        l, ql, al, alpha, beta, gamma, hatm
            lwrite = .false.
          elseif (neul .eq. nl .and. lwrite .and. mode0 .lt. 2) then
            write(stdo,'(i5,2f10.6,3x,3f10.6,3x,3f10.6,3x,f10.6)')
     .        l, ql, al, alpha, beta, gamma, hatm, ddot(3,hatm,1,bhat,1)
            lwrite = .false.
          endif

        enddo
        lwrite = .true.
        if (mod(mode0,2) .eq. 0 .or. ipr .lt. PRT2) lwrite = .false.
      enddo
      amom = dsqrt(ddot(3,aloc,1,aloc,1))

      if (amom .ne. 0) then
C       Assign bhat to point along magnetization direction
        if (mode0 .ge. 2) call dpcopy(aloc,bhat,1,3,1/amom)

C       If sum moments < 0, optionally reverse B
        if (sal .lt. 0 .and. mode2 .eq. 1 .and. mode0 .ge. 2) then
          if (mode0 .ge. 2) call dpcopy(aloc,bhat,1,3,-1/amom)
        endif

C       Printout of bhat, moment, angle
        call info5(PRT1,0,0,' ATOM='//clabl//
     .    '%a  bhat=%3;10,6D  |M|=%;6,6d  bhat.M/|M|=%;6d',bhat,
     .    amom,ddot(3,aloc,1,bhat,1)/amom,0,0)

      else
        bhat(1) = 0
        bhat(2) = 0
        bhat(3) = 1
      endif

C     End of block to contruct bhat
      endif

C --- Rotate the qnu along projection bhat ---
      mode1 = mod(mode/10,10)
      if (mode1 .eq. 0) return

      if (nsp .eq. 1 .or. mode1 .eq. 1 .or.
     .    ddot(3,bhat,1,bhat,1) .eq. 0d0) then
        call dcopy(3*nl*nsp,qnu,1,qnuloc,1)
        call dcopy(1*nl*nsp,pnu,1,pnuloc,1)
        return
      endif

      call dpzero(qnuloc,3*nl*nsp)
      ilm = 0
      do  l = 0, nl-1
        call dpzero(aloc,3)
        dploc = 0
        do   m = -l, l
          ilm = ilm+1
          if (neul .eq. nl) then
            alpha = eula(l+1,1)
            beta  = eula(l+1,2)
C           gamma = eula(l+1,3)
          elseif (neul .eq. nl*nl) then
            alpha = eula(ilm,1)
            beta  = eula(ilm,2)
C           gamma = eula(ilm,3)
          elseif (neul .eq. 1) then
            alpha = eula(1,1)
            beta  = eula(1,2)
C           gamma = eula(1,3)
          else
            call rxi('atscpp: bad value neul=',neul)
          endif

C         Charge, magnetic moments, quantization axis for these angles
          hatm(1) = dcos(alpha)*dsin(beta)
          hatm(2) = dsin(alpha)*dsin(beta)
          hatm(3) = dcos(beta)
          dotp = hatm(1)*bhat(1) + hatm(2)*bhat(2) + hatm(3)*bhat(3)

          do  i = 1, 3
            aloc(i) = aloc(i) + (qnu(i,l+1,1)-qnu(i,l+1,2))/(2*l+1)*dotp
          enddo
          dploc = dploc + (pnu(l+1,1)-pnu(l+1,2))/(2*l+1)*dotp

        enddo

        do  i = 1, 3
          qloc(i) = qnu(i,l+1,1) + qnu(i,l+1,2)
C         aloc(i) = qnu(i,l+1,1) - qnu(i,l+1,2)
          qnuloc(i,l+1,1) = (qloc(i) + aloc(i))/2
          qnuloc(i,l+1,2) = (qloc(i) - aloc(i))/2
        enddo
        ploc = pnu(l+1,1) + pnu(l+1,2)
        pnuloc(l+1,1) = (ploc + dploc)/2
        pnuloc(l+1,2) = (ploc - dploc)/2

      enddo

      if (ipr .ge. PRT2) then
        write(stdo,'(''  l isp'',19x,''qnu'',37x,''qloc'')')
        do  k = 1, 2
        do  l = 0, nl-1
          write(stdo,'(2i3,1x,3f13.7:1x,3f13.7)')
     .      l,k,(qnu(i,l+1,k),i=1,3),(qnuloc(i,l+1,k),i=1,3)
        enddo
        enddo
      endif

Cgetarg       end
       end subroutine asprjq 


