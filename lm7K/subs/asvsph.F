      subroutine asvsph(sctrl,slat,sarray,sspec,sham,spot,vrl,imake,
     .  ehterm,lhave)
C- Make ASA sphere potential, potential parameters
C ----------------------------------------------------------------------
Ci Inputs
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read: nbas nbasp nclass nspec zbak nl nspin lncol loptc
Ci                lpgf lves lrel
Ci     Stored:
Ci     Passed to: lgors atscpp asamad
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: vol avw
Ci     Stored:
Ci     Passed to: asamad
Ci   sarray:structure containing offsets to various arrays
Ci     Elts read: nclasp ormax ohave oics oclabl onrcp oipc
Ci     Stored:
Ci     Passed to: asamad
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: eref
Ci     Stored:
Ci     Passed to: atscpp asamad
Ci   sham  :struct for parameters defining hamiltonian; see routine uham
Ci     Elts read: nmto kmto eterms
Ci     Stored:    eterms ehk thrpv seref amgm
Ci     Passed to:
Ci   spot  :struct for information about the potential; see routine upot
Ci     Elts read: orhrmx oves ovdif ovrmax ovintr opnu oqnu oqt opprel
Ci                opp osop ogrrme opmpol
Ci     Stored:    vmtz
Ci     Passed to: asamad
Ci   vrl   :(pgf) difference in potential between left- and right- leads.
Ci         :      See lmasa-gf.f
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
Ci         :   in the potential, no potential parameters generated;
Ci         :   ves = estat potential for supplied Q is not retained.
Ci         :   no atom file written.  Using this mode to make terms
Ci         :   in Kohn-Sham energy, e.g. v_in, rho_out.
Co Outputs
Co   ehterm:(1) ehkk (no Mad)
Co          (2) sphere band sum (Ves=0 at RMT)
Co          (3) sum Q_i V(R)
Co          (4) emad
Co          NOTE: ehterm is obsolete and will be phased out in favor of
Co          sham->eterms
Co   sham->eterms integrals for the total energy are accumulated
Co         :(1)  ehar   --- not touched here
Co         :(2)  eks    --- not touched here
Co         :(3)  utot   = total electrostatic energy
Co         :(4)  valves --- not used by ASA
Co         :(5)  cpnves --- not used by ASA
Co         :(6)  rhoexc = rho * exc
Co         :(7)  rhovxc = rho * vxc
Co         :(8)  sumec  = sum-of-core eigenvalues
Co         :(9)  sumtc  --- not used by ASA
Co         :(10) xcore  = rhoc * total potential
Co         :(11) valvef = rhov * total potential
Co         :(12) sumt0  --- not touched here
Co         :(13) dq1    --- not touched here
Co         :(14) dq2    --- not touched here
Co         :(15) amom   =  total magnetic moment
Co         :(16) sumev  =  sum-of-eigenvalues from moments
Co         :(17) rinvxt --- not touched here
Co         :(18) rouvxt --- not touched here
Co         :(19) bmval  = M<B> : magnetic contribution to valvef (asadc)
Co   lhave :1s digit
Co          0 insufficient input available to generate potentials
Co            sought
Co          1 All potentials sought were calculated or read from disk
Co         10s digit
Co          0 insufficient input available to generate ehkk,ehterm
Co          1 ehkk, ehterm have been generated
Co   See Remarks
Cl Local variables
Cl    ehkk  :ASA total energy as sum-of-sphere energies
Cl          :This is the HK total energy in versions 6.12 and earlier
Cr Remarks
Cu Updates
Cu   21 Jul 07 (pgf) vne->vrl (for inequivalent left- and right- end layeers)
Cu   10 Feb 04 (S.Faleev) vrl added to argument list; passed 
Cu              to asamad for non-equilibrium mode
Cu    1 Apr 04 Sphere program for orbital-dependent XC field
Cu   26 Apr 03 Added MPI calls
Cu   18 Mar 03 (A Chantis) relativistic potential parameters.
Cu   10 Mar 03 First cut at making proper Kohn-Sham energy.
Cu   17 Feb 03 Added double-counting terms to ehk for applied field
Cu   15 Feb 03 B-field requires SO parameters
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer lhave
      double precision sctrl(1),slat(1),sarray(1),sspec(1),sham(1),
     .  spot(1),ehterm(5),vrl
C ... Local parameters
      logical bittst,lgors,lehk,havedc
      integer oclabl,ogrrme,oics,oinitc,oeula,oipc,onrcp,opmpol,opnu,
     .  opp,oqnu,oqt,orhrmx,ormax,osop,ovdif,oves,ovwk,ovintr,ovrmax,
     .  owk,opprel,oeulat,obxc
      integer ic,nclass,nclasp,nclspp,nbas,nl,nsp,igets,i,iclbsj,
     .  nbasp,lpgf,imake,is,lncol,loptc,lnsph,isw,lgunit,lves,k,ival,
     .  nrclas,nspec,neul
      integer mpipid,procid
      double precision thrpv,amgm,ehkk,emad,trumad,vmtz(2),zbak(2),
     .  dgetss,wk(5),dval,amag(3),dsqrt,ddot
C ... For atscpp
      character clabl*8
      double precision amgmat,amom,ehkat,sevat,thrpva,vol,avw,seref
C ... For third-generation LMTO
      integer nmto
      double precision kmto(20)
C     double precision qat,zat,dval
C ... Heap
      integer w(1)
      common /w/ w

C     This ordering must match sham->eterms; see uham
      double precision eterms(20),eterma(20)

      call upack('ctrl nbas nbasp nclass nspec zbak',sctrl,
     .  nbas,nbasp,nclass,nspec,zbak)
      call upack('ctrl nl nspin lncol loptc',sctrl,nl,nsp,lncol,loptc,0)
      call upack('array nclasp ormax',sarray,nclasp,ormax,0,0,0)
      call upack('array ohave oics oclabl onrcp oipc',sarray,oinitc,
     .  oics,oclabl,onrcp,oipc)
      call upack('pot orhrmx oves ovdif ovrmax ovintr',spot,orhrmx,
     .  oves,ovdif,ovrmax,ovintr)
      call upack('pot opnu oqnu oqt opprel',spot,opnu,oqnu,oqt,opprel,0)
      call upack('pot opp osop obxc ogrrme opmpol',spot,
     .  opp,osop,obxc,ogrrme,opmpol)
      call upack('ham nmto kmto',sham,nmto,kmto,0,0,0)
      call upack('ham oeula neula',sham,oeula,neul,0,0,0)
      call upack('lat vol avw',slat,vol,avw,0,0,0)
      lpgf = igets('ctrl lpgf',sctrl)
      lnsph = isw(lgors('ctrl lasa,32',sctrl))
      lves = 2*igets('ctrl lves,1',sctrl)
      if (lpgf .ne. 0) then
        lves = lves+100
      endif

      ehkk = 0
      amgm = 0
      thrpv = 0
      seref = 0
      nclspp = 2*nclasp-nclass
      call dpzero(ehterm,4)
      call dpzero(amag,3)
      call upack1('ham eterms',sham,eterms)
C     asvsph only modifies eterms(3:11)
      call dvset(eterms,3,11,0d0)
      call dvset(eterms,15,16,0d0)

C ... MPI: only master does sphere program
      procid = mpipid(1)
      if (procid .eq. 0) then

      call togprt()
      call defrr(owk,-nclspp)
      lehk = .true.
      lhave = 1
      havedc = .true.
      do  40  ic = 1, nclasp
        nrclas = ival(w(onrcp),ic)
C   ... Classes to skip
        i = iclbsj(ic,w(oipc),-nbasp,1)
        if (i .lt. 1) goto 40
C   ... What to make
        if (lpgf .eq. 2 .and. i .le. nbas) goto 40
C       is = w(oics+ic-1)
        is = ival(w(oics),ic)

C       call pshpr(80)
        call dpzero(eterma,20)
        if (neul .gt. 0) then
          call defrr(oeulat,neul*3)
          call dmscop(w(oeulat),1,w(oeula),nbas,i,i,1,neul*3,1,1,1d0)
        endif
        call atscpp(sctrl,sspec,is,ic,w(oclabl),w(ormax),imake,
     .    nl,nsp,w(oinitc+ic-1),zbak(1)/vol,avw,w(opnu),w(oqnu),w(oves),
     .    w(owk),w(oeulat),neul,w(obxc),0d0,w(opp),w(opprel),ehkat,
     .    sevat,w(oqt),amgmat,amom,w(orhrmx),w(ovrmax),thrpva,w(osop),
     .    w(ogrrme),w(opmpol),w(ovintr),clabl,eterma)
        if (neul .gt. 0) then
          call rlse(oeulat)
          call dpsadd(amag,w(obxc),3,1,3*ic-2,amom*nrclas)
        endif
C       If estat energy changed, double-counting terms are available
        havedc = havedc .and. eterma(3) .ne. 0
        if (havedc) then
          eterma(4) = dval(w(oves),ic) * dval(w(oqt),ic) / 2
          call dpadd(eterms,eterma,3,11,dble(nrclas))
          call dpadd(eterms,eterma,15,16,dble(nrclas))
          call dpadd(eterms,eterma,19,19,dble(nrclas))
        else
          call dvset(eterms,3,11,0d0)
        endif
C       call poppr
        seref = seref + dgetss('spec eref',is,sspec)*nrclas
        if (imake .eq. 4) goto 40

C   ... Check that potential and other pars available if to continue
        k = w(oinitc+ic-1)
        if (.not. lgors('ctrl lgen3',sctrl)) then
        if (mod(k/2,2) .ne. 1
     .    .or.(bittst(lncol,4).or.bittst(lncol,8)).and.(mod(k/4,2)).ne.1
     .    .or. loptc.gt.0 .and. (mod(k/32,2)) .ne. 1
     .    .or. lnsph .ne. 0 .and. (mod(k/16,2)) .ne. 1) then
          call awrit0(' LM (warning) class '//clabl//
     .     '%a missing pot, s-orbit or optical parms',' ',-80,lgunit(1))
          lhave = 0
          goto 40
        endif
        endif

C   ... Only accumulate into ehkk if ehkat is set
C       qat = dval(w(oqt),ic)
C       zat = dgetss('spec z',is,sspec)
        if (havedc .and. lehk) then
          ehkat = ehkat - dgetss('spec eref',is,sspec)
          ehkk = ehkk + ehkat*nrclas
        else
          lehk = .false.
        endif
        call atsev(ic,w(onrcp),w(oves),nl,nsp,w(oqnu),sevat,ehterm(2))
        amgm = amgm + amgmat*nrclas
        thrpv = thrpv + thrpva*nrclas
   40 continue
      call rlse(owk)
      call togprt()

C --- Make ves, emad, vmtz for new sphere Q, potential ---
C     imake=4: don't update pot->ves; use temporary array
      if (imake .eq. 4) then
        call defrr(ovwk,-nclspp)
      else
        ovwk = oves
      endif
      call asamad(sctrl,spot,sarray,slat,sspec,lves,
     .  w(opnu),w(oqnu),vrl,w(ovwk),emad,trumad,vmtz)
C ... V[q] needed for total energy if a potential shift
      call dpzero(w(ovdif),nclspp)
      if (mod(lves,100) .eq. 2) then
        call pshpr(0)
        call asamad(sctrl,spot,sarray,slat,sspec,0,w(opnu),w(oqnu),
     .    vrl,w(ovdif),wk(1),wk(2),wk(3))
        call daxpy(nclspp,-1d0,w(ovwk),1,w(ovdif),1)
        call poppr
      endif
      if (imake .eq. 4) call rlse(ovwk)

C ... Add Madelung contribution to total electrostatic energy
      eterms(3) = eterms(3) + emad
C ... Total magnetization, noncollinear case
      if (neul .gt. 0) then
        amgm = dsqrt(ddot(3,amag,1,amag,1))
        eterms(15) = amgm
      endif
C ... Repack double-counting terms
      call pack1('ham eterms',sham,eterms)

C --- Make NMTO potential parameters ---
C      if (lgors('ctrl lgen3',sctrl)) then
C        call ppmesh(sctrl,sspec,slat,nl,nmto,kmto,w(oclabl),
C     .    vmtz,w(oics),w(oppn))
C      endif

C --- Save energy for HF (to recover if information lost) ---
C terms:
C (1) ehkk (no Mad) (2) sumeV(VH=0)
C (3) sum Q_i V(R) (4) emad
C (5) d.c. terms from applied fields.
      ehterm(4) = emad
      if (lehk) then
        ehterm(1) = ehkk
        ehkk = ehkk + emad - ehterm(5)
      endif

      thrpv = thrpv + trumad
      call pack1('pot vmtz',spot,vmtz)
      call pack5('ham ehk thrpv seref amgm',sham,ehkk,thrpv,seref,amgm,
     .  0)
      if (lehk) lhave = lhave+10

      endif

C --- MPI broadcast everything passed out of asvsph ---
C     sarray and contents of sarray pointers
      call mpibc1(w(oinitc),nclasp,2,.false.,'asvsph','initc')
C     sham
      call mpibc1(sham,nint(sham(1))*nbas,4,.false.,'asvsph','sham')
C     spot and and contents of spot pointers
      call mpibc1(w(orhrmx),nclasp,4,.false.,'asvsph','rhrmx')
      call mpibc1(w(oves),nclspp,4,.false.,'asvsph','ves')
      call mpibc1(w(ovdif),nclspp,4,.false.,'asvsph','vdif')
      call mpibc1(w(ovrmax),2*nclspp,4,.false.,'asvsph','vrmax')
      if (ovintr .gt. 0) then
        i = nclasp*(nl*nsp)**2
        call mpibc1(w(ovintr),i,4,.false.,'asvsph','vintr')
      endif
      call mpibc1(w(opnu),nl*nsp*nclasp,4,.false.,'asvsph','pnu')
      call mpibc1(w(oqnu),3*nl*nsp*nclasp,4,.false.,'asvsph','qnu')
      call mpibc1(w(oqt),nclasp,4,.false.,'asvsph','qt')
      call mpibc1(w(opp),6*nl*nsp*nclasp,4,.false.,'asvsph','pp')
C     lrel = igets('ctrl lrel',sctrl)
C     if (lrel .eq. 2) then
      if (opprel .gt. 0) then
        call mpibc1(w(opprel),32*nl*nl*max(nclspp,nspec),4,.false.,
     .    'asvsph','pprel')
      endif
C     if (igets('ctrl lncol',sctrl) .gt. 0) then
      if (osop .gt. 0) then
       call mpibc1(w(osop),9*nl*nsp*nsp*nclasp,4,.false.,'asvsph','sop')
      endif
C     if (igets('ctrl loptc',sctrl) .gt. 0) then
      if (ogrrme .gt. 0) then
       call mpibc1(w(ogrrme),8*nl*nsp*nclasp,4,.false.,'asvsph','grrme')
      endif
C     if (lgors('ctrl lasa,32',sctrl)) then
      if (opmpol .gt. 0) then
       i = (2*nl-1)*nl**2*3*nsp*nclasp
       call mpibc1(w(opmpol),i,4,.false.,'asvsph','pmpol')
      endif
      call mpibc1(spot,nint(spot(1)),4,.false.,'asvsph','spot')
C     Remaining arguments passed to asvsph
      call mpibc1(ehterm,5,4,.false.,'asvsph','ehterm')
      call mpibc1(lhave,1,2,.false.,'asvsph','lhave')

      end

