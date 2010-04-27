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
