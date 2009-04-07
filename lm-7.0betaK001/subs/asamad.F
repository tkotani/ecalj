#define NSPH 1
      subroutine asamad(sctrl,spot,sarray,slat,sspec,lq,p,q,
     .  vrl,ves,emad,trumad,vmtz)
C- Calculate ASA Madelung potential, Madelung energy, vmtz
C ----------------------------------------------------------------
Ci Inputs
Ci   sctrl :struct for program flow parameters; see routine uctrl
Ci     Elts read: nbas nl nspec nspin lpgf lmet nclass zbak
Ci     Stored:    *
Ci     Passed to: lgors asavqm
Ci   spot  :struct for information about the potential; see routine upot
Ci     Elts read: omad ovrmax ovshf oqpp vmtz0
Ci     Stored:    *
Ci     Passed to: asavqm
Ci   sarray:structure containing offsets to various arrays
Ci     Elts read: oics oclabl nclasp onrcp oipcp npadl npadr ormax
Ci     Stored:    *
Ci     Passed to: asavqm
Ci   slat  :struct for lattice information; see routine ulat
Ci     Elts read: vol alat as tol nkdmx nkqmx opos platl platr
Ci     Duplicate: vol
Ci     Stored:    *
Ci     Passed to: asavqm
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: lmxa z
Ci     Stored:    *
Ci     Passed to: asavqm
Ci   lq   :  0 input p,q are pnu,qnu; make qt from pnu,qnu
Ci           1 input q is qt; p is not used
Ci           2 Take ves as input; evaluate vmtz only
Ci          10s digit nonzero substitutes nonspherical ves for 0,1
Ci             iff qpp available
Ci             NB: qpp should be passed as argument, and lq more
Ci                 sensibly arranged.
Ci          100s digit applies to the layer branch (lpgf>0)
Ci               0 No special treatment: estat potential of all layers
Ci                 computed from sphere charges
Ci               1 Estat potential of sites in embedding region
Ci                 computed from Madelung, but the potential of
Ci                 sites in L- and R- regions are computed as though
Ci                 in bulk with periodic boundary conditions.
Ci               2 Display the deviation in the potential in L- and R-
Ci                 PL as computed by the two approaches above.
Ci               3 Add to vconst(1) a best average of this deviation.
Ci   p     :boundary conditions.  If Dl = log. deriv. at rmax,
Ci          pnu = .5 - atan(Dl)/pi + (princ.quant.number).
Ci          Used only to compute nuclear contribution to sphere charge
Ci   q     :energy-weighted moments of the sphere charges.
Ci          Used to compute electronic contribution to sphere charge
Cl   vrl   :Difference in estat potential in the left- and right- electrodes (PGF)
Ci         :vrl = bias across device: + [vconst(3)-vconst(2)]
Ci         :First term = ef(R)-efermi, where ef(R) is fermi energy of right lead
Ci         :Second term is needed because of how Poisson's eqn is solved (see lmasa-gf.f)
Ci  sarray: requires arrays oics,onrcp,oipcp,ovrmax,oz
Co Outputs (see Remarks for qualifications)
Co   ves:  Hartree potential on sphere surfaces
Co   emad: Madelung energy using convention as discussed in remarks
Co   trumad:True interatomic Madelung energy; see madpot
Co   vmtz: muffin-tin zero (vmtz(2) = vmtz(spin-up)-vmtz(spin-dn))
Cl Local variables
Cl   These variables only apply to the layer case
Cl   vavo :electrostatic potential of the layers computed from
Cl        :multilayer
Cl        :(1) = avg ves in L
Cl        :(2) = avg ves in R
Cl   vavn :electrostatic potential computed from charges for L- and R-
Cl        :endpoints individually , under the condition each endpoint
Cl        :is subject to its own periodic boundary conditions as
Cl        :L- and R- layers were separate (periodic) materials.
Cl        :(1) = avg ves in L
Cl        :(2) = avg ves in R
Cr Remarks
Cr   asamad generates the ASA electrostatic potential by Ewald
Cr   summation (see madpot.f, below).  There is an additional
Cr   potential that may be added; this has no fundamental signifance
Cr   but has a practical significance in various contexts.  First,
Cr   double counting of the potential shifts must be eliminated in
Cr   the total energy.  Second, in some cases the electrostatic
Cr   potential is not sought, for models or to facilitate convergence
Cr   to self-consistency.  Last, in cases when the energy mesh or
Cr   Fermi level is explicit, as in Green's function calculations,
Cr   the mesh is affected by potential shifts.
Cr
Cr  *Case 1s digit lq eq 2:
Cr   The input ves is used.
Cr
Cr  *The USUAL band structure case (lq ne 2):
Cr   Electrostatic potential calculated directly by Ewald summation.
Cr
Cr  *The LAYER Green's function case (lq eq 0,1):
Cr   Case separate bulk calculations (lpgf eq 2):
Cr   The electrostatic potentials are calculated for the left-
Cr   and right- bulk PL independently.
Cr
Cr   Case interface calculation (lpgf eq 1):
Cr   The electrostatic potentials initially calculated as in the usual
Cr   band-structure case.  Were the left and right bulk-like PL
Cr   really bulk-like, this would be sufficient, apart from a possibly
Cr   constant shift for insulating L or R.  (The shift is zero for
Cr   a metallic end PL).  But there may arise unwanted deviations owing
Cr   to the nonlocal contribution from other layers.  The potentials
Cr   in these layers are obtained from the electrostatic potentials of
Cr   the bulk, plus a constant shift, calculated as:
Cr       the average electrostatic potential of the interface -
Cr       the average electrostatic potential of the bulk PL
Cr   The caller can force either L or R, or both shifts be zero.
Cr   Establish this condition by ctrl->lmet, bit 2 (L) and bit 3 (R).
Cu Updates
Cu   21 Jul 07 (pgf) vne->vrl (for inequivalent left- and right- end layeers)
Cu   10 Feb 04 (S.Faleev) sbz struc added to argument list
Cu              to handle non-equilibrium mode
Cu   22 Mar 03 Revised electrostatics for layer case;
Cu             estat potential is be shifted by pot->vconst
Cu   19 Feb 02 Some changes adapt to updates in layer code
C ----------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer lq
      double precision sctrl(1),spot(1),sarray(1),slat(1),
     .  sspec(1),ves(20),emad,trumad,vmtz(2),p(1),q(1),vrl
C ... Local parameters
      logical lves,lgors,metL,metR,bittst
      integer nbas,nl,nsp,nspec,nclasp,nclspp,npadl,npadr,nkdmx,nkqmx,
     .  lpgf,iprint,nbaspp,lmxst,nkd,nkq,ic,ic1,ic2,jc1,jc2,ib,nglob,
     .  nclass,lmet,lnsph,isw,igets,nlmf,ival,stdo
      integer oics,oclabl,onrcp,oipcp,omad,ovrmax,ormax,odq,olmx,oz,oqc,
     .  oqt,odlv,oqlv,opos,opos2,omad2,ovshft,oves,oqpp,ovval,oqmp,ovh
      double precision rhrmx(2),dgets,platl(3,3),platr(3,3),alat,awald0,
     .  tol,vol,plat(9),qlat(9),awald,xx,vavo(5),vavn(5),dval,dum(3),
     .  vglob,vconst(3),rmsv(3)
C ... Parameters for non-equilibrium mode (S.F.)
      double precision plat0(3,3)
C     double precision vadd,fac,vshfl,vshfr,vshfc
C     character*8 clabl
C ... Heap
      integer w(1)
      common /w/ w

C      call pshpr(30)
C      print *, '!! asamad',lq

C ... Setup
      call upack('ctrl nbas nl nspec nspin',sctrl,nbas,nl,nspec,nsp,0)
      lpgf = igets('ctrl lpgf',sctrl)
      call upack('ctrl lmet nclass zbak',sctrl,lmet,nclass,rhrmx,0,0)
      call dscal(2,1/dgets('lat vol',slat),rhrmx,1)
      call upack('array oics oclabl nclasp onrcp oipcp',sarray,oics,
     .  oclabl,nclasp,onrcp,oipcp)
      call upack('array npadl npadr',sarray,npadl,npadr,0,0,0)
      call upack('pot vconst omad ovrmax ovshf',spot,vconst,omad,ovrmax,
     .  ovshft,0)
      call upack('lat alat plat opos',slat,alat,plat0,opos,0,0)
      stdo = nglob('stdo')

      nclspp = 2*nclasp-nclass
      nbaspp = 2*(npadl+npadr) + nbas
      call upack1('array ormax',sarray,ormax)
      call defdr(oqt,nclspp)
      if (mod(lq,100) .eq. 0) then
        call defdr(odq,nclspp)
        call defdr(oqc,nclspp)
        call sp2cls('spec lmxa',sspec,w(oics),1,1,nclspp,olmx)
        call sp2cls('spec z',sspec,w(oics),1,1,nclspp,oz)
        call getq(nsp,nl,w(olmx),nclspp,w(oz),p,q,w(oics),
     .    sspec,w(oqc),w(oqt),w(odq))
C       call snot('qt',nbaspp,w(oqt),w(oipcp))
        call rlse(odq)
      elseif (mod(lq,100) .eq. 1) then
        call dcopy(nclspp,q,1,w(oqt),1)
      else
        call dpzero(w(oqt),nclspp)
      endif
      lves = mod(lq,100) .eq. 2

C ... Copy charges and potential to doubly padded layer
C      if (nclasp .ne. nclass) then
C        call dpscop(w(oqt),w(oqt),nclasp-nclass,nclass+1,nclasp+1,1d0)
C        call dpscop(ves,ves,nclasp-nclass,nclass+1,nclasp+1,1d0)
C        call dpscop(w(ovrmax),w(ovrmax),2*(nclasp-nclass),
C     .    2*nclass+1,2*nclasp+1,1d0)
CC        call prmx('vrmax',w(ovrmax),2,2,nclspp)
CC        call prmx('qt',w(oqt),nclspp,nclspp,1)
C      endif

C      call setpr(50)
C      call dvset(w(oqt),nclass+1,nclass+1,.2d0)
*      call dvset(w(oqt),nclass+2,nclass+2,-.2d0)

#if NSPH
C --- Make potential from multipole moments ---
      lnsph = isw(lgors('ctrl lasa,32',sctrl))
      if (lnsph .ne. 0 .and. .not. lves) then
      oqpp = igets('pot oqpp',spot)
      if (dval(w(oqpp),1) .ge. 0d0) then
        nlmf = (2*nl-1)**2
        call defrr(ovval, -nlmf*nbas)
        call defrr(oqmp,   nlmf*nbas)
        call defrr(ovh,  nbas)
        call asavqm(0,sctrl,spot,slat,sarray,sspec,nlmf,w(ovh),w(oqmp),
     .    w(ovval))
        call dpzero(ves,nclasp)
        do  60  ib = 1, nbas+(npadl+npadr)
          ic = w(oipcp+ib-1)
          xx = 1/dble(ival(w(onrcp),ic))
          ves(ic) = ves(ic) + dval(w(ovh),ib)/ival(w(onrcp),ic)
          xx = dval(w(ovh),ib)
   60   continue
        lves = .true.
        call rlse(ovval)
      endif
      endif
#endif

C --- Electrostatic potential for charges as given ---
      if (lpgf .ne. 2) then
        if (lpgf .ne. 0) call info(20,1,0,
     .  ' Electrostatics for embedded L-C-R system:',0,0)
        vmtz(1) = dgets('pot vmtz0',spot)

C       call snot('rmax',nbaspp,w(ormax),w(oipcp))
C       Kirill's fix (March 2006) superseded by update in lmasa-gf
C        call madpot(nbaspp,1,nclspp,w(onrcp),w(oipcp),w(oclabl),
C     .    w(oqt),vconst,rhrmx,w(ormax),w(omad),w(opos),plat0,vrl+
C     .    (vconst(3)-vconst(2)),lves,
        call madpot(nbaspp,1,nclspp,w(onrcp),w(oipcp),w(oclabl),
     .    w(oqt),vconst,rhrmx,w(ormax),w(omad),w(opos),plat0,vrl,lves,
     .    w(ovrmax),ves,emad,trumad,vmtz)
        if (abs(vrl) .ge. 1d-5) call dplmom(nbaspp,w(opos),alat,plat0,
     .                           w(oipcp),w(oqt))
C       call snot('ves[qt]',nbaspp,ves,w(oipcp))
C   ... Nothing further if potentials are input
        if (lves) goto 99
      endif

C --- Layer GF: es pot for end layers and pot shift for ch. neutrality
      if (lpgf .ne. 0 .and. mod(lq/100,100) .gt. 0) then
        call upack('lat as tol vol',slat,awald0,tol,vol,0,0)      
        call upack('lat nkdmx nkqmx',slat,nkdmx,nkqmx,0,0,0)      
        call upack2('lat platl platr',slat,platl,platr)
C       nbasp = nbas + npadl + npadr
C   ... Save original ves(Mad) for printout
        call defrr(oves, nclasp)
        call dpcopy(ves,w(oves),1,nclasp,1d0)

C   ... Average potentials for interface L,R, to calculate shift
        ic1 = nclasp
        ic2 = 1
C       vshfl = dval(w(ovshft),2)
        do  10  ib = nbas, nbas+npadl-1
          ic1 = min(ic1,w(oipcp+ib))
          ic2 = max(ic2,w(oipcp+ib))
   10   continue
        jc1 = nclasp
        jc2 = 1
C       vshfr = dval(w(ovshft),4)
        do  12  ib = nbas+npadl, nbas+npadl+npadr-1
          jc1 = min(jc1,w(oipcp+ib))
          jc2 = max(jc2,w(oipcp+ib))
   12   continue
C       vshfc = dval(w(ovshft),6)
        call pvmadp(ic1,ic2,w(onrcp),0,w,w(ormax),ves,vavo(1),rmsv)
        call pvmadp(jc1,jc2,w(onrcp),0,w,w(ormax),ves,vavo(2),rmsv(2))

C        print *, ic1,ic2
C        print *, jc1,jc2
C       stop

C        call dvset(w(oqt),2+ic1,2+ic1,.1d0)
C        call dvset(w(oqt),2+ic1+1,2+ic1+1,-.1d0)

C   ... Electrostatic potential for bulk L
        call defrr(odlv,  3*nkdmx)
        call defrr(oqlv,  3*nkqmx)
        lmxst = 6
        call pshpr(0)
        call lattc(awald0,tol,0d0,alat,alat,platl,0d0,0d0,1d0,1d0,plat,
     .    qlat,lmxst,vol,awald,w(odlv),nkd,w(oqlv),nkq,nkdmx,nkqmx)
        call defrr(opos2, max(npadl,npadr)*3)
        call dpscop(w(opos),w(opos2),3*npadl,3*nbas+1,1,1d0)
c        call prmx('pos',w(opos2),3,3,npadl)
        call defdr(omad2,max(npadl,npadr)**2)
        call madmat(npadl,w(opos2),awald,alat,vol,w(odlv),nkd,
     .    w(oqlv),nkq,w(omad2))
        call poppr
c       call prmx('mad',w(omad2),npadl,npadl,npadl)
        vavn(1) = 0
        call info(20,1,0,
     .    ' ES potential for infinitely repeating L princ layer:',0,0)
        call madpot(npadl,ic1,ic2,w(onrcp),w(oipcp+nbas),w(oclabl),
     .    w(oqt),vconst(2),rhrmx,w(ormax),w(omad2),w,w,0d0,lves,
     .    w(ovrmax),ves,emad,xx,vavn(1))
C       Keep vavn(1) for sanity check (see if unchanged after bulk R)
        call pvmadp(ic1,ic2,w(onrcp),0,w,w(ormax),ves,dum(1),rmsv)

C   ... Electrostatic potential for bulk R
        call pshpr(0)
        call lattc(awald0,tol,0d0,alat,alat,platr,0d0,0d0,1d0,1d0,plat,
     .    qlat,lmxst,vol,awald,w(odlv),nkd,w(oqlv),nkq,nkdmx,nkqmx)
        call dpscop(w(opos),w(opos2),3*npadr,3*(nbas+npadl)+1,1,1d0)
        call madmat(npadr,w(opos2),awald,alat,vol,w(odlv),nkd,
     .    w(oqlv),nkq,w(omad2))
        call poppr
        vavn(2) = 0
        call info(20,1,0,
     .    ' ES potential for infinitely repeating R princ layer:',0,0)
        call madpot(npadr,jc1,jc2,w(onrcp),w(oipcp+nbas+npadl),
     .    w(oclabl),w(oqt),vconst(3),rhrmx,w(ormax),w(omad2),w,w,0d0,
     .    lves,w(ovrmax),ves,emad,trumad,vavn(2))

        call pvmadp(ic1,ic2,w(onrcp),0,w,w(ormax),ves,vavn(1),rmsv)
        call pvmadp(jc1,jc2,w(onrcp),0,w,w(ormax),ves,vavn(2),rmsv(2))

C       Sanity check: did <ves(L)> change after calc ves(R)?
        if (abs(dum(1)-vavn(1)) .gt. 1d-6 .and. iprint() .ge. 10) then
          call info2(10,1,0,'  asamad (warning): adjusting vconst(L):'//
     .      ' average ves(L) altered by %,6;6d%N%8f after ves(R) '//
     .      'calculated (possible misuse of equivalent classes)',
     .      vavn(1)-dum(1),0)
          vconst(2) = vconst(2)+vavn(1)-dum(1)
          call pack1('pot vconst',spot,vconst)
        endif

C   ... RMS deviation in the two ways to compute ves[padded]
        if (lpgf .ne. 2) then
          call dpadd(ves,w(oves),1,nclasp,-1d0)
          call pvmadp(ic1,ic2,w(onrcp),0,w,w(ormax),ves,dum(1),rmsv)
          call pvmadp(jc1,jc2,w(onrcp),0,w,w(ormax),ves,dum(2),rmsv(2))
          call dpadd(ves,w(oves),1,nclasp, 1d0)

          if (iprint() .ge. 30) then
            metL = bittst(lmet,8)
            metR = bittst(lmet,4)
            write(stdo,334)
  334       format(/' Deviations in end potentials:'/' region met',
     .        '  <ves>Bulk','   <ves>layer     Diff',6x,'RMS diff')
            write(stdo,335) 'L', metL,
     .        vavn(1), vavo(1), vavo(1)-vavn(1), dsqrt(rmsv(1))
            write(stdo,335) 'R', metR,
     .        vavn(2), vavo(2), vavo(2)-vavn(2), dsqrt(rmsv(2))
  335       format(2x,a,L7,5f12.6)
          endif
          rmsv(3) =
     .          dsqrt((rmsv(1)**2*npadl+rmsv(2)**2*npadr)/(npadl+npadr))
          call info5(20,1,0,' RMS pot difference in '//
     .      'L PL = %,6;6d  in R PL = %,6;6d  total = %,6;6d',
     .      rmsv,rmsv(2),rmsv(3),0,0)
          vglob = vconst(1) + (vavn(1) - vavo(1) - vavo(2) + vavn(2))/2

          call info0(40,1,0,' vconst should be fixed by charge '//
     .      'neutrality requirements.'//
     .      '%N%8fIt can found iteratively in a GF pass using mode 1.')

          if (iprint() .lt. 30)
     .      call info2(20,0,0,' vconst is now %,6;6d  vconst '//
     .      'minimizing RMS diff to end layers = %,6;6d',vconst,vglob)
          if (iprint() .ge. 30) then
          call info2(20,0,0,' vconst that minimizes potential '//
     .      'mismatch to end layers =%;10,6D',vglob,0)
          call info2(20,0,0,' vconst is now (estimate to'//
     .      ' satisfy charge neutrality)  =%;10,6D',vconst,0)
          call info2(20,0,0,' difference%45f=%;10,6D',vglob-vconst(1),0)
          endif

          if (iprint() .ge. 40) call query('vconst=',4,vconst)

          call pack1('pot vconst',spot,vconst)
        else
          call info2(20,1,0,' shifts added to bulk ES potential: '//
     .      ' vconst(L)=%,6;6d   vconst(R)=%,6;6d',vconst(2),vconst(3))
        endif

C   OLD
CC   --- Global potential shift, layer case ---
CC       We are permitted to shift all layers by a single constant
CC       An end layer that is metallic should not have any potential
CC       shifts to satisfy charge neutrality.  However, this requirement
CC       cannot be exactly satisfied if both end layers are metallic,
CC       because there is only one free parameter.  In that case we pick
CC       the average of the two layers.  If only one is metallic, it
CC       determines the global potential shift.  If neither is metallic,
CC       choose the shift from the central layer shift.
C        if (lpgf .ne. 2 .and. mod(lq/100,10) .ge. 2) then
C
CC         fac=0 if show shifts, but don't add to vshft; fac=1 to add
C          fac = 0
C          if (mod(lq/100,10) .ge. 3) fac = 1
C
CC     ... Compute global potential shift according to comments above
CC         and assign vavo(3..5) to shifts we will actually add
CC         to existing values.
CC         L layer: vglob would be vavn(1) - vavo(1)
CC         R layer: vglob would be vavn(2) - vavo(2)
CC         C layer: vglob would be dval(w(ovshft),6) ??
C          metL = bittst(lmet,8)
C          metR = bittst(lmet,4)
C
C          if (metL .and. metR) then
C            vglob   = (vavn(1) - vavo(1) - vavo(2) + vavn(2))/2
C            vavo(3) = 0*vshfl
C            vavo(4) = 0*vshfr
C            vavo(5) = vglob
C          elseif (metL) then
C            vglob   = (vavn(1) - vavo(1))
C            vavo(3) = 0*vshfl
C            vavo(4) = vavo(4) - (vavo(1) - vavn(1))
C            vavo(5) = vglob
C          elseif (metR) then
C            vglob   = (vavn(2) - vavo(2))
C            vavo(3) = vavo(3) + (vavn(2) - vavo(2))
C            vavo(4) = 0*vshfr
C            vavo(5) = vglob
C          else
C            vglob = vavo(5)
C            stop 'not ready'
C            vavo(3) = vavo(3) + vglob
C            vavo(4) = vavo(3) + vglob
C            vavo(5) = 0
C          endif
C
C          if (iprint() .ge. 20) then
C            write(stdo,334) metL,metR,vglob
C  334       format(/' Global potential shift to best satisfy charge',
C     .        ' neutrality:'/' met(L) = ',L1,'  met(R) = ',L1,
C     .        '  global vshft =',f11.6//
C     .        ' region <ves>Bulk  <ves>layer     Diff',
C     .        4x,'<vshft_in>',2x,'<vshft_out>')
C            write(stdo,335) 'L', vavn(1), vavo(1), vavo(1)-vavn(1),
C     .        vshfl, vshfl+vavo(3)*fac
C            write(stdo,335) 'R', vavn(2), vavo(2), vavo(2)-vavn(2),
C     .        vshfr, vshfr+vavo(4)*fac
C            write(stdo,336) 'C', vshfc,vshfc*(1-fac)+vavo(5)*fac
C  335       format(2x,a,1x,5f12.6)
C  336       format(2x,a,37x,2f12.6)
C          endif
C          vshfl = vshfl + vavo(3)*fac
C          vshfr = vshfr + vavo(4)*fac
CC         vshfc = vshfc + vavo(5)*fac
C          vshfc = vshfc*(1-fac)+vavo(5)*fac
C          call dvset(w(ovshft),2,2,vshfl)
C          call dvset(w(ovshft),4,4,vshfr)
C          call dvset(w(ovshft),6,6,vshfc)
C          if (fac .ne. 0) call dvset(w(ovshft),8,8,1d0)
C
C          if (fac .ne. 0) then
C            write(stdo,500)
C            write(lgunit(2),500)
C            do  20  ic = 1, jc2
C              vadd = vshfc
C              if (ic .gt. nclass) vadd = vshfl
C              if (ic .gt. ic2) vadd = vshfr
C              call r8tos8(dval(w(oclabl),ic),clabl)
Cc              ves(ic) = ves(ic) + vadd
C              if (iprint() .ge. 30) then
C                do  21  i = 1, 2
C   21           write(lgunit(i),501) clabl,dval(w(oqt),ic),
C     .              dval(w(oves),ic),ves(ic),ves(ic)+vadd,ves(ic)+vadd+
C     .              dval(w(ovrmax),2*ic-1)
C  501           format(1x,a,5f11.6)
C  500           format(/' Class',8x,'Qtot    ',
C     .            'Vh(Iface)   Vh(Bulk)   Vh+Shift    V(Rmax)')
C              endif
C   20       continue
C          endif
C        endif

      endif

   99 call rlse(oqt)

C      if (iprint() .ge. 20) pause
C      call poppr

      end

      subroutine madpot(nbas,ic1,ic2,nrc,ipc,dclabl,qt,vconst,rhrmx,
     .  rmax,dmad,bas,plat,vrl,noves,vrmax,ves,emad,trumad,vmtz)
C- Calculate Madelung potential, Madelung energy, vmtz
C ----------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   ic1,ic2: calculate Madelung potential for classes ic1 .. ic2
Ci   nrc   :number of memembers of each class
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   dclabl:class name, packed as a real number
Ci   qt    :electronic charge, less nuclear charge
Ci         :(excluding homogeneous nuclear charge in rhrmx(2))
Ci   vconst:constant potential shift added to estat potentials
Ci   rhrmx :constant background density (for Muffin-Tin correction)
Ci          rhrmx(1) = constant background density
Ci          rhrmx(2) = density to be incorporated in MT correction
Ci   rmax  :radius at which potential is evaluated, in a.u.
Ci   dmad  :Madelung matrix
Ci   noves :if true, evaluate vmtz only (using input ves)
Ci   vrmax :non-estat (XC) potential at wsr, for calc. vmtz
Co Outputs
Co   ves:  Hartree potential on sphere surfaces
Co   emad: Madelung energy using convention as discussed in remarks
Co   trumad:True interatomic Madelung energy
Co   vmtz: muffin-tin zero (vmtz(2) = vmtz(up)-vmtz(dn))
Cr Remarks
Cr   Potential at site i, due to point charges from all points in the
Cr   lattice except i itself is given by (see madmat):
Cr
Cr   (1)  V_i(0) = 2 \sum_j M_ij DQ_j
Cr   DQ_j are the total charges and j runs over all sites in the
Cr   basis.  Rydberg units, with electronic charge positive,
Cr   nuclear charge negative.  The average Hartree potential at
Cr   the sphere radius rmax is
Cr
Cr   (2) V_i(rmax) = V_i(0)  +  2 DQ_i / rmax
Cr
Cr   where V_intra,i = 2 DQ_i / rmax is the intra-atomic contribution
Cr   to the potential.  The atomic code calculates potentials and total
Cr   energies subject to the boundary condition that the Hartree
Cr   potential at rmax is zero, regardless of the total charge inside
Cr   the sphere.  Then there is an onsite contribution to the
Cr   sphere energy,
Cr
Cr   (3) E_intra = 1/2 \int rho V_intra d^3r = 1/2 DQ_i V_intra,i
Cr
Cr   and the usual interatomic terms, \sum_i V_i(0) DQ_i, so that
Cr   the net Madelung contribution to the total energy is
Cr
Cr   (4) E_mad = 1/2 \sum_i V_i(rmax) DQ_i
Cr
Cr   Alternatively, the charges outside sphere i shift the potential
Cr   inside sphere i by a constant V_i(rmax).  This contributes
Cr   an extra term to the total energy Eqn 4, and also the
Cr   potential parameters for sphere i are shifted by V_i(rmax).
Cr
Cr   ---- Addition of background charge to sphere charges ---
Cr
Cr   ? The sphere charges consist of the electronic-nuclear charge qt,
Cr   less the smeared-out nuclear charge rhrmx(1)*vol
Cr
Cr   ---- MUFFIN-TIN Correction and homogeneous background ---
Cr
Cr   If the net system charge is a sum of nuclear charges, sphere
Cr   charges and a constant background charge of density n, the
Cr   charge DQ remains the net charge inside each sphere.  The
Cr   Madelung potential is calculated with charges Ztwid_i =
Cr   -DQ_i - Qtwid, ie the total charge in each sphere excluding
Cr   the charge Qtwid in sphere i from the homogeneous
Cr   background.  Eqn (2) still describes the potential at rmax.
Cr   Then for each R the Hartree energy of (Ztwid + homogeneous
Cr   background) is subtracted, since we use the Hartree energy
Cr   of the true density there. This term is
Cr
Cr   (5)  -3 Ztwid Qtwid / R  +  6 Qtwid^2 / 5 R
Cr
Cu Updates
Cu   21 Mar 03 New vconst
C ----------------------------------------------------------------
C     implicit none
      integer nbas,ic1,ic2
      logical noves
      double precision vconst,vrl,emad,qt(ic2),vrmax(2,ic2),rmax(ic2),
     .  bas(3,nbas),plat(3,3),dmad(nbas,nbas),ves(ic2),vmtz(2),
     .  trumad,dclabl(ic2),rhrmx(2)
      integer ipc(nbas),nrc(ic2)
C Local variables
      character*8 clabl
      integer i,ibas,jbas,ic,jc,iprint,lgunit
      double precision vadd,sumdq,fpi
      integer ix,scrwid,isw,awrite
      double precision dq,qih,qtw,ebak,ezv,sezv,sebak,sqb,avves,rmsv(2)
      parameter (scrwid=80)
C     Parameters for non-equilibrium mode
      double precision qlat(3,3),vol,b3(3),Lz,wk(nbas),xx

C ... Net, inhomogeneous, homo, background electronic & nuclear charges
      dq(ix)  = qt(ix) - fpi/3*rmax(ix)**3*rhrmx(2)
      qih(ix) = qt(ix) - fpi/3*rmax(ix)**3*rhrmx(1)
      qtw(ix) = fpi/3*rmax(ix)**3*(rhrmx(1)-rhrmx(2))
C     qbk(ix) = fpi/3*rmax(ix)**3*rhrmx(1)
C     zbk(ix) = fpi/3*rmax(ix)**3*rhrmx(2)

      fpi = 16*datan(1d0)
C ... Setup for linear potential drop across plat(3)
      if (vrl .ne. 0) then
        call dinv33(plat,1,qlat,vol)
        do  i = 1, 3
          b3(i) = qlat(i,3)/sqrt(qlat(1,3)**2+qlat(2,3)**2+qlat(3,3)**2)
        enddo
        Lz = b3(1)*plat(1,3) + b3(2)*plat(2,3) + b3(3)*plat(3,3)
      endif

      if (noves) then
        sumdq = 0
        if (iprint() .ge. 30) then
          do  131  i = 1, 2
  131     write(lgunit(i),502)
        endif
        do  130  ic = ic1, ic2
        sumdq = sumdq + nrc(ic)*dq(ic)
        if (iprint() .ge. 30) call r8tos8(dclabl(ic),clabl)
        do  130  i = 1, 2
          if (iprint() .ge. 30)
     .      write(lgunit(i),501) clabl,dq(ic),qtw(ic),ves(ic),
     .      ves(ic)+vrmax(1,ic)
  130     continue
  502     format(/' MADPOT: potential not calculated.  Input Q,V:'/
     .      ' Class        Qtot       Qbak       V0        V0+Vxc')
          emad = 0
          trumad = 0
          goto 40
      endif

      avves = 0
      do  10  ibas = 1, nbas
        ic = ipc(ibas)
        if (ic .lt. ic1 .or. ic .gt. ic2) goto 10
        ves(ic) = vconst
        if (vrl .ne. 0) then
          xx =( bas(1,ibas)*b3(1)+bas(2,ibas)*b3(2)+bas(3,ibas)*b3(3))
     .      /Lz*vrl
          ves(ic) = ves(ic) + xx
          wk(ic)  = xx
        endif
        do  11  jbas = 1, nbas
          jc = ipc(jbas)
          if (jc .lt. ic1 .or. jc .gt. ic2) goto 11
          ves(ic) = ves(ic) + 2*qih(jc)*dmad(ibas,jbas)
   11   continue
        avves = avves + (ves(ic) + 2*dq(ic)/rmax(ic))/nbas
   10 continue
      if (iprint() .ge. 30) then
        do  12  i = 1, 2
          if (vrl .eq. 0) then
            write(lgunit(i),500)
          else
            write(lgunit(i),600)
          endif
   12   continue
      endif
      emad   = 0
      trumad = 0
      sumdq  = 0
      sebak  = 0
      sezv   = 0
      sqb    = 0
      do  30  ic = ic1, ic2
        if (iprint() .ge. 30) call r8tos8(dclabl(ic),clabl)
        vadd = 2*dq(ic)/rmax(ic)
        ebak = (3*qih(ic)*qtw(ic) + 6*qtw(ic)**2/5)/rmax(ic)
        ezv = 0.5d0*qih(ic)*ves(ic)
        if (iprint() .ge. 30 .and. ebak .eq. 0) then
          do  15  i = 1, 2
            if (vrl .eq. 0) then
              write(lgunit(i),501) clabl,dq(ic),qtw(ic),ves(ic),
     .          ves(ic)+vadd,ves(ic)+vadd+vrmax(1,ic)
            else
              write(lgunit(i),601) clabl,dq(ic),qtw(ic),wk(ic),
     .          ves(ic),ves(ic)+vadd,ves(ic)+vadd+vrmax(1,ic)
            endif
   15     continue
        elseif (iprint() .ge. 30 .and. ebak .ne. 0) then
          do  16  i = 1, 2
            if (vrl .eq. 0) then
              write(lgunit(i),501) clabl,dq(ic),qtw(ic),ves(ic),
     .          ves(ic)+vadd,ves(ic)+vadd+vrmax(1,ic),ezv-ebak
            else
              write(lgunit(i),601) clabl,dq(ic),qtw(ic),wk(ic),ves(ic),
     .          ves(ic)+vadd,ves(ic)+vadd+vrmax(1,ic),ezv-ebak
            endif
   16     continue
        endif
  501   format(1x,a,5f11.6:'  ezv-ebak=',f9.6)
  601   format(1x,a,6f11.6:'  ezv-ebak=',f9.6)
  500   format(/' Class',8x,
     .    'Qtot       Qbak       Vmad     Vh(Rmax)    V(Rmax)')
  600   format(/' Class',8x,
     .  'Qtot       Qbak      vrl*z/L     Vmad     Vh(Rmax)    V(Rmax)')

        ebak = nrc(ic)*ebak
        ezv =  nrc(ic)*ezv
        sumdq = sumdq + nrc(ic)*dq(ic)
        trumad = trumad + ezv - ebak
        emad  = emad + ezv - ebak + .5d0*dq(ic)*vadd*nrc(ic)
        sezv  = sezv + ezv
        sebak = sebak + ebak
        sqb   = sqb + qtw(ic)*nrc(ic)
        ves(ic) = ves(ic) + vadd
        if (rhrmx(1) .ne. 0) ves(ic) = ves(ic) - avves
   30 continue

C --- Reset Madelung potential if total system charge .gt. .5 ---
      if (dabs(sumdq) .gt. .01d0 .and. iprint() .ge. 10) print
     .    '('' (warning)  system has net charge: q='',f10.6)', sumdq
C      if (dabs(sumdq) .gt. .5d0) then
C        emad = 0
C        if (iprint() .ge. 10) print
C     .    '('' madpot:  net charge too big; reset ves to zero'')'
C        call dpzero(ves(ic1),ic2-ic1+1)
C      endif

C --- Make vmtz by averaging over surfaces of spheres ---
   40 continue
      vmtz(2) = 0
      if (vmtz(1) .eq. 0)
     .  call pvmadp(ic1,ic2,nrc,1,vrmax,rmax,ves,vmtz,rmsv)

C --- Printout ---
      if (iprint() .ge. 20 .and. .not. noves) then
        do  22  i = 1, 2
C   22   call awrit6(' Sum Q=%,6d  Emad=%,6d(%,6d)  '//
C     .      'Vmtz=%,6d%?;n; (up-dn=%,6d);;',' ',scrwid,lgunit(i),
C     .      sumdq, emad, trumad,vmtz,isw(vmtz(2).ne.0),vmtz(2))
        ix = awrite(' Sum Q=%,6d  Emad=%,6d(%,6d)'//
     .    '%?#n# Vconst=%,6;6d#%j#  Vmtz=%,6d%?;n; (up-dn=%,6d);;',' ',
     .    scrwid,lgunit(i),sumdq, emad, trumad,
     .    isw(vconst.ne.0),vconst,
     .    vmtz,isw(vmtz(2).ne.0),vmtz(2))
   22   continue

C  440   format(' Sum Q=',f10.6,'  Emad=',f10.6,'(',f10.6,')'/
C     .         ' Vmtz=',f9.6,' (up-dn)=',f9.6)
        if (rhrmx(1) .ne. 0) then
          print 441, sqb, sezv, sebak, sezv-sebak
  441     format(' SUMQB=',f9.6,'  EMAD=',f10.6,'(RB)',
     .      f10.6,'(SPHR)','   DIFF=',f9.6)
          print 442, avves
  442     format(' Shifting average VES by',f10.6)
        endif

      endif
      if (iprint() .ge. 20 .and. noves) then
        do  122  i = 1, 2
  122   call awrit4(' Sum Q=%,6d  Vmtz=%,6d%?;n; (up-dn=%,6d);;',' ',
     .      scrwid,lgunit(i),sumdq, vmtz,isw(vmtz(2).ne.0),vmtz(2))
C        do  122  i = 1, 2
C  122   write(lgunit(i),540) sumdq, vmtz
C  540   format(' SUM Q=',f10.6,'  VMTZ=',f9.6,'   VMTZ(UP-DN)=',f9.6)
      endif

      end
      subroutine pvmadp(ic1,ic2,nrc,lvrmx,vrmax,rmax,ves,vmtz,rmsv)
C- Make vmtz
C ----------------------------------------------------------------------
Ci Inputs
Ci   ic1   :starting class index
Ci   ic2   :ending class index
Ci   nrc   :number of members of this calss
Ci   lvrmx :1 if to add vrmax; 0 if to just average ves
Ci   rmax  :augmentation radius, in a.u.,
Ci   ves   :estat potential, by class
Co Outputs
Co   vmtz  :muffin-tin zero
Co   rmsv  :sum ves(ic)**2 wgt(ic)
Cu Updates
Cu   22 Mar 03 New argument rmsv
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer lvrmx,ic1,ic2,nrc(ic2)
      double precision vrmax(2,ic2),rmax(ic2),ves(ic2),vmtz(2),rmsv(2)
C ... Local parameters
      integer ic
      double precision swgt,wgt

      vmtz(1) = 0d0
      vmtz(2) = 0d0
      swgt = 0d0
      rmsv(1) = 0
      rmsv(2) = 0
      do  50  ic = ic1, ic2
        wgt = nrc(ic)*rmax(ic)**3
        if (lvrmx .eq. 1) then
          vmtz(1) = vmtz(1) + wgt*(vrmax(1,ic) + ves(ic))
          vmtz(2) = vmtz(2) + wgt*vrmax(2,ic)
          rmsv(1) = rmsv(1) + wgt*ves(ic)**2
          rmsv(2) = rmsv(2) + wgt*(vrmax(1,ic) + ves(ic))**2
        else
          vmtz(1) = vmtz(1) + wgt*ves(ic)
          rmsv(1) = rmsv(1) + wgt*ves(ic)**2
        endif
        swgt = swgt + wgt
   50 continue
      vmtz(1) = vmtz(1)/swgt
      vmtz(2) = vmtz(2)/swgt
      rmsv(1) = rmsv(1)/swgt
      rmsv(2) = rmsv(2)/swgt

      end

      subroutine dplmom(nbas,bas,alat,plat,ipc,qt)
C- Prints out contribution to dipole moment n plat(3) from excess charge
C ----------------------------------------------------------------------
Ci Inputs
Ci   nbas  :size of basis
Ci   bas   :basis vectors, in units of alat
Ci   alat  :length scale of lattice and basis vectors, a.u.
Ci   plat  :primitive lattice vectors, in units of alat
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   qt    :electronic charge, less nuclear charge
Ci         :(excluding homogeneous nuclear charge in rhrmx(2))
Co Outputs
Co   Nothing is output: routine is informational only
Cl Local variables
Cl         :
Cr Remarks
Cr   
Cu Updates
Cu   15 Sep 04 
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nbas,ipc(nbas)
      double precision bas(3,nbas),alat,plat(3,3),qt(*)
C ... Local parameters
      integer ibas,i,ic,stdo,nglob,iprint
      double precision qlat(3,3),vol,b3(3),pi,z,Lz,dplm,delv

      if (iprint() .lt. 20) return

      stdo = nglob('stdo')
      pi = 4*datan(1d0)
      call dinv33(plat,1,qlat,vol)
      do  i = 1, 3
       b3(i) = qlat(i,3)/sqrt(qlat(1,3)**2+qlat(2,3)**2+qlat(3,3)**2)
      enddo
      Lz = b3(1)*plat(1,3) + b3(2)*plat(2,3) + b3(3)*plat(3,3)

      dplm = 0d0
      do  ibas = 1, nbas
        ic = ipc(ibas)
        z = b3(1)*bas(1,ibas) + b3(2)*bas(2,ibas) + b3(3)*bas(3,ibas)
        dplm = dplm + 2*z*qt(ic)
      enddo
      delv = 4*pi*dplm*Lz/vol/alat
      write(stdo,333) delv
  333 format(' Dipole in plat(3) from charge distribution =',
     .  f10.6)
      end

C      subroutine uefld(nbas,bas,alat,plat,dmad)
CC- Add contribution of uniform electric field E=4*pi/vol*sum_i Q_i*z_i
CC- to Madelung constants due dipole moment of supercell
CC ----------------------------------------------------------------------
CCi Inputs
CCi Inputs
CCi   nbas  :size of basis
CCi   bas   :basis vectors, in units of alat
CCi   alat  :length scale of lattice and basis vectors, a.u.
CCi   plat  :primitive lattice vectors, in units of alat
CCo Input/Outputs
CCo   dmad  :constant E-field contribution added to Madelung matrix
CCu Updates
CCu   9 Feb 04 (S.Faleev) First created
CC ----------------------------------------------------------------------
C      implicit none
CC ... Passed parameters
C      integer nbas
C      double precision alat,plat(3,3),bas(3,nbas),dmad(nbas,nbas)
CC ... Local parameters
C      integer i,j,ibas,jbas
C      double precision pi,b3(3),s(3),qlat(3,3),vol,z1,z2
C
C      pi = 4*datan(1d0)
C      call dinv33(plat,1,qlat,vol)
C      vol = dabs(vol)
C      do i = 1,3
C        b3(i) = qlat(i,3)/sqrt(qlat(1,3)**2+qlat(2,3)**2+qlat(3,3)**2)
C      enddo
C
CC ... Checks
C      call dpzero(s,3)
C      do j = 1,3
C        do i = 1,3
C          s(j) = s(j) + b3(i)*plat(i,j)
C        enddo
C      enddo
C      if ( (abs(s(1))+abs(s(2))) .gt. 1d-10 .or. s(3) .le. 0d0 )
C     .  call rx('uefld: wrong qlat(*,3)')
C
C      do  jbas = 1,nbas
C        do  ibas = 1, nbas
C          z1 = 0d0
C          z2 = 0d0
C          do  i = 1, 3
C            z1 = z1 + bas(i,ibas)*b3(i)
C            z2 = z2 + bas(i,jbas)*b3(i)
C          enddo
C          dmad(ibas,jbas) = dmad(ibas,jbas) + 4d0*pi*z1*z2/(vol*alat)
C        enddo
C      enddo
C
CCCCCCCCCCCCCCCC TEMP CCCCCCCCCCCCC
Cc      write(*,'(a,i4,4f12.6)')'uefld:nbas,b3,vol=',nbas,b3,vol
Cc      write(*,'(a,9f8.4)')    'uefld:plat=',plat
Cc      do i = 1,nbas
Cc      write(*,'(a,i4,3f12.6)')'uefld: ibas,bas=',i,bas(1,i),bas(2,i),
Cc     .        bas(3,i)
Cc      enddo
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      end

C      subroutine snot(strn,nbas,fc,ipc)
C      implicit none
C      character *(*) strn
C      integer nbas,ipc(nbas)
C      double precision fc(*)
C      integer ib
C      double precision fb(nbas)
C
C      do  ib = 1, nbas
C        fb(ib) = fc(ipc(ib))
C      enddo
C      call prmx(strn,fb,nbas,nbas,1)
C      end

