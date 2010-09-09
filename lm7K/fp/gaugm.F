CSFCPP#define F90 1
      subroutine gaugm(nr,nsp,lso,lcplxp,rofi,rwgt,lmxa,lmxl,nlml,vsm,
     .gpotb,gpot0,hab,vab,sab,sodb,qum,vum,cg,jcg,indxcg,
     .nf1,nf1s,lmx1,lx1,f1,x1,v1,d1,
     .nf2,nf2s,lmx2,lx2,f2,x2,v2,d2,
     .lmux,sig,tau,nlx1,nlx2,ppi,
     .lmaxu,vumm,lldau,idu)
C- Generic routine to make augmentation matrices
C ----------------------------------------------------------------------
Ci Inputs
Ci   nr    :number of radial mesh points
Ci   rofi  :radial mesh points
Ci   rwgt  :radial mesh weights
Ci   lmxa  :augmentation L-cutoff
Ci   lmxl  :l-cutoff for density, potential on the radial mesh
Ci   nlml  :(lmxl+1)*(lmxl+1)
Ci   vsm   :smooth local potential
Ci   gpotb :integrals of compensating gaussians * local smooth estat
Ci         :pot calculated from the compensated smoothed local density
Ci   gpot0 :integrals of local gaussians * phi0~ (smves.f)
Ci         :phi0~ is the estatic potential of the interstitial
Ci         :smooth density including compensating local charges.
Ci   hab   :matrix elements of the hamiltonian for the true wave
Ci         :functions and true (spherical part of) potential,
Ci         :where w.f. are expressed in forms of (ul,sl,gz) as
Ci         :described in potpus.f
Ci   vab   :corresponding matrix elements of the (spherical part of) V
Ci   sab   :corresponding matrix elements of the overlap
Ci   sodb  :corresponding matrix elements of SO = 2/(c^2) dV/dr*(1/r)
Ci   qum   :moments (u,s) * (u,s) * r**l (momusl.f)
Ci   vum   :integrals ((u,s,gz) * (true potential) * (u,s,gz))
Ci         :for the full (nonspherical) true potential (momusl.f)
Ci   cg    :Clebsch Gordon coefficients, stored in condensed form (scg)
Ci   jcg   :L q.n. for the C.G. coefficients stored in condensed form
Ci   indxcg:index for Clebsch Gordon coefficients
Ci   ...    Two sets of functions f1 and f2 are input, given as:
Ci   nf1   :number of 'bra' function types for each l
Ci   nf1s  :number of 'bra' function types for each l
Ci         :which are formed by linear combination of (u,s)
Ci   lmx1  :dimensions f1
Ci   lx1   :l-cutoffs for each of the nf1 functions
Ci   f1    :'bra' unaugmented envelope fn (radial part . r).
Ci         :Must be zero in channels that have no envelope functions.
Ci         :See Remarks
Ci   x1    :radial part of grad(f1/r) times r
Ci   v1    :values of f1 at rofi(nr) (not multiplied by r)
Ci   d1    :slopes of f1 at rofi(nr) (not multiplied by r)
Ci   nf2   :number of 'ket' function types for each l
Ci   nf2s  :number of 'ket' function types for each l
Ci         :which are formed by linear combination of (u,s)
Ci   lmx2  :dimensions f2
Ci   lx2   :l-cutoffs for each of the nf2 functions
Ci   f2    :'ket' unaugmented envelope fn (radial part . r).
Ci         :Must be zero in channels that have no envelope functions
Ci         :See Remarks.
Ci   x2    :radial part of grad(f2/r) times r
Ci   v2    :values of f2 at rofi(nr) (not multiplied by r)
Ci   d2    :slopes of f2 at rofi(nr) (not multiplied by r)
Ci   lmux  :l-cutoff for sigma,tau, and spherical part of ppi
Ci   nlx1  :dimensions ppi
Ci   nlx2  :dimensions ppi
Ci   lso   :if 2 calculate LzSz matrix elements, if 1 calculate LS
Ci  lcplxp :0 if ppi is real; 1 if ppi is complex
Ci   ...   The following are associated with LDA+U
Ci   lmaxu :dimensioning parameter for U matrix
Ci   vumm  : orbital dependent pot matrix in (u,s) representation
Ci   lldau :lldau=0 => no U on this site otherwise
Ci         :U on this site
Ci   idu   :idu(l+1)=1 => this l has a nonlocal U matrix
Co Outputs
Co   sig   :augmentation overlap integrals; see Remarks.
Co   tau   :augmentation kinetic energy integrals
Co   ppi   :augmentation potential integrals
Co         :NB: tau is added to pi, so ppi = (kinetic energy + potential)
Co         :ppi is returned real if lcplxp=0, complex if lcplxp=1
Co         :In the noncollinear case:
Co         :ppi(:,:,:,:,isp,1,1) = spin-diagonal part of potential
Co         :ppi(:,:,:,:,1,2,1) = up-down   (12) block of potential
Co         :ppi(:,:,:,:,2,2,1) = down-down (21) block of potential
Co         :In the SO case: the SO hamiltonian is added to ppi
Co         :and also stored in ppi(:,:,:,:,:,:,2)
Cl Local variables
Cl   ppi0  :contribution to ppi from spherical part of potential
Cl   qm    :multipole moments; see Remarks
Cl   hso   :matrix elements of SO hamiltonian inside sphere; see Remarks
Cl         :hso(:,:,:,:,isp,1) = spin-diagonal part of potential
Cl         :hso(:,:,:,:,1,2)   = up-down   (12) block of potential
Cl         :hso(:,:,:,:,2,2)   = down-down (21) block of potential
Cl         :On exit, hso is added to ppi; see ppi
Cr Remarks
Cr   This subroutine assembles the various terms in the computation of
Cr   sigma, tau, pi that comprise the local (augmented) part of the
Cr   hamiltonian.  See also Remarks in routine augmat.f.  Equations
Cr   references are found in methods paper; see book by Springer:
Cr      M. Methfessel, M. van Schilfgaarde, and R. A. Casali,
Cr      Lecture Notes in Physics, {\bf 535}. H. Dreysse,
Cr      ed. (Springer-Verlag, Berlin) 2000.
Cr
Cr   Augmentation integrals sig, tau, ppi are constructed as a
Cr   difference between augmented function products and unaugmented ones.
Cr
Cr     sig = integral f1~*f2~ - f1^*f2^.
Cr     tau = integral f1~ (-nabla) f2~ - f1^ (-nabla) f2^.
Cr         = integral -1* (grad f1~ grad f2~ - grad f1^ grad f2^)
Cr           since  f1~ grad f2~ = f1^ grad f2^ at the MT boundary.
Cr           NB: Integration by parts also works for true local orbitals
Cr           (without smooth parts ^) because val, slo = 0 at rmt.
Cr     ppi = integral f1~ vtrue f2~ - f1^ vsm f2^.
Cr
Cr   f1^, f2^ are one-center expansions to the envelope function.
Cr
Cr   The augmented part of sig and tau and ppi (spherical potential only)
Cr   is obtained from linear combinations of sab, hab, and vab.
Cr   The smoothed part is calculated directly by numerical integration.
Cr
Cr   hab,vab,sab are matrix elements of h, spherical-v, and 1 of (u,s)
Cr   where u and s are linear combinations of and phi,phidot defined as:
Cr   u has val=1, slo=1 at rmax, s has val=0, slo=1.  Matrix elements
Cr   of the potential are divided into a spherical part and the
Cr   remainder: the (dominant) spherical term is computed with the
Cr   relativistic wave function including the small component, while
Cr   the residual nonspherical part does not.
Cr
Cr   There two additional contributions to ppi.  One involves integrals
Cr
Cr   QkkLL'M integral g_M (V0~-V2~), where the integrals are passed as
Cr   gpot0 and gpotb.   QkkLL'M (cf the Springer book chapter Eq 27)
Cr   is evaluated in pvagm3 (QkkLL'M is called qm here).
Cr      QkkLL'M gpot0 corresponds to Eq. 28;
Cr      QkkLL'M gpotb corresponds to the last term in Eq. 29.
Cr
Cr   QkkLL'M  = integrals ((u,s) * 1 * (u,s)) for valence functions,
Cr   and corresponding integrals for local orbitals, is obtained from
Cr   input qum.
Cr
Cr   The second term involves matrix elements of v-v(l=0) which are
Cr   obtained in the augmented term from linear combinations of input
Cr      vum = integrals ((u,s) * vtrue * (u,s))
Cr
Cr   Notes on local orbitals:
Cr   They occur either as true local orbitals (val,slo=0 at rmt) or
Cr   as extended orbitals, where a tail is attached at the MT sphere.
Cr   In the latter case, a smooth part may be subtracted, to
Cr   be included later as part of matrix elements of the interstitial
Cr   potential.  Caller must ensure that f1,f2 for all orbitals with
Cr   no smooth part are zero; and for extended local orbitals where
Cr   a smooth part is to be subtracted, that they match continuously
Cr   and differentiably to the smooth functions f1,f2 at rmt.
Cr
Cr   Notes on Spin-Orbit Coupling
Cr
Cr   The (spin-diagonal) LzSz inside augmentation sphere
Cr   has matrix elements in real harmonics:
Cr    <l,m|Lz|l,m'> = +i|m| delta_m,-m if m>0
Cr                    -i|m| delta_m,-m if m<0
Cr   where delta the Kronecker delta. The Hamiltonian is ordered
Cr   -m,...,0,...,m.  The p block for example looks like:
Cr          m:   -1       0     1
Cr         -1     0       0   i<|so|>
Cr          0     0       0      0
Cr          1 -i<|so|>    0      0
Cr
Cr   The (spin-off-diagonal) LxSx+LySy inside augmentation sphere.
Cr   has these matrix elements:
Cr     Spin up-down block is 1/2<l,m|L-|l,m'>; L- = Lx - i Ly
Cr     Spin down-up block is 1/2<l,m|L+|l,m'>; L+ = Lx + i Ly
Cr     The p and blocks have the following form: <p,m|L-|p,m'> =
Cr               -1                    0               1
Cr      -1        0           -i*a1*<|so|>/sqrt(2)     0
Cr       0 i*a2*<|so|>/sqrt(2)         0       -a2*<|so|>/sqrt(2)
Cr       1        0              a1*<|so|>/sqrt(2)     0
Cr     where a1 = sqrt{(l-m)*(l+m+1)} and a2 = sqrt{(l+m)*(l-m+1)}
Cr     <p,m|L+|p,m'> is the hermitian conjugate of <p,m|L-|p,m'>.
Cr
Cr   Notes on LDA+U
Cr   Call to paugnl adds nonlocal potential part to ppi matrix
Cr
Cu Updates
Cu   27 Nov 05 LDA+U => complex potential
Cu   09 Nov 05 (wrl) Convert dmat to complex form
Cu   26 Sep 05 (A. Chantis) Bug fix: local orbitals in conjunction w/ SO
Cu   29 Jun 05 (MvS) redesign to keep hso separated
Cu   27 Apr 05 LDA+U (Lambrecht)
Cu   03 Feb 05 (A. Chantis) calculate matrix elements of L.S
Cu             inside augmentation sphere
Cu    1 Sep 04 Adapted to handle complex ppi; so folded into ppi
Cu   12 Aug 04 Added treatment for extended local orbitals.
Cu             Envelopes f1,f2 must be zero for all channels that
Cu             have no smooth counterparts to subtract.
Cu   20 Jul 04 bug fix in pvagm3
Cu   29 Jun 04 (A. Chantis) calculate matrix elements of LzSz
Cu             inside augmentation sphere
Cu   14 Sep 01 Extended to local orbitals.  Altered argument list.
Cu   17 May 00 Adapted from nfp gaugm.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      use m_globalvariables
      integer nr,nsp,lmxa,lmxl,nlml,lmux,lso,lcplxp,
     .nf1,nf1s,lmx1,nlx1,lx1(nf1),jcg(1),indxcg(1),n0,nab,
     .nf2,nf2s,lmx2,nlx2,lx2(nf2)
      integer lmaxu,lldau,idu(4)
      parameter (n0=10,nab=9)
      double complex vumm(-lmaxu:lmaxu,-lmaxu:lmaxu,nab,2,0:lmaxu)
      double precision rofi(nr),rwgt(nr),vsm(nr,nlml,nsp),cg(1),
     .qum(0:lmxa,0:lmxa,0:lmxl,6,nsp),gpotb(1),gpot0(1),
     .hab(nab,n0,nsp),vab(nab,n0,nsp),sab(nab,n0,nsp),
     .f1(nr,0:lmx1,nf1s),v1(0:lmx1,nf1s),d1(0:lmx1,nf1s),
     .f2(nr,0:lmx2,nf2s),v2(0:lmx2,nf2s),d2(0:lmx2,nf2s),
     .x1(nr,0:lmx1,nf1s),x2(nr,0:lmx2,nf2s),
     .sig(nf1,nf2,0:lmux,nsp),tau(nf1,nf2,0:lmux,nsp),
     .ppi(nf1,nf2,nlx1,nlx2,nsp,1),vum(0:lmxa,0:lmxa,nlml,6,nsp)
C     Spin-Orbit related
      double precision sodb(nab,n0,nsp,2)
C ... Local parameters
      integer ilm1,ilm2,nspc,nglob
      integer i1,i2,ilm,l,ll,nlm,nlm1,nlm2,i,iprint,stdo,lgunit
      double precision pi,y0
      double precision vsms(nr),ppi0(nf1*nf2*(lmux+1)),
     .qm(nf1*nf2*(lmx1+1)*(lmx2+1)*(lmxl+1)),
     .sum((lmx1+1)*(lmx2+1)*nlml)
CSFCPP#if F90
      complex(8),allocatable:: hso(:,:,:,:,:,:),ppiz(:,:,:,:,:,:)
CSFCPP#else
CSFCPP      double complex hso(nf1,nf2,nlx1,nlx2,2,2)
CSFCPP      double complex ppiz(nf1,nf2,nlx1,nlx2,2,2)
CSFCPP#endif
#if DEBUG
      integer fopna
      integer kount
      save kount
      data kount /0/
      kount = kount + 1
#endif

      stdo = lgunit(1)
      pi = 4d0*datan(1d0)
      y0 = 1d0/dsqrt(4d0*pi)
Changenglob      nspc = nglob('nspc')
      nspc = globalvariables%nspc

CSFCPP#if F90
      allocate(hso(nf1,nf2,nlx1,nlx2,2,2),ppiz(nf1,nf2,nlx1,nlx2,2,2))
CSFCPP#endif
C     for now, ppi is complex only if S.O coupling
      if ((lso .ne. 0 .or. lldau .ne. 0) .and. lcplxp .eq. 0)
     .call rx('gaugm: incompatible lso,lcplxp')
      if (lso .ne. 0) then
        call dpzero(hso,nf1*nf2*nlx1*nlx2*2*2)
      endif
      if (lldau .ne. 0) then
        call dpzero(ppiz,2*nf1*nf2*nlx1*nlx2*2*2)
      endif
      do  i = 1, nsp

C     Spherical part of the smooth potential
        call dpcopy(vsm(1,1,i),vsms,1,nr,y0)

C --- Make sig, tau, ppi0 = spherical part of ppi ---
C     pvagm2: smooth f1^*f2^ part of sig and corresponding tau, ppi0
        call pvagm2(nf1,lmx1,lx1,f1,x1,nf2,lmx2,lx2,f2,x2,
     .  nr,rofi,rwgt,vsms,lmux,sig(1,1,0,i),tau(1,1,0,i),ppi0)
C     print *, '!! zero ppi0 (sm), sig(sm)'
C     ppi0 = 0
C     sig(1:nf1,1:nf2,0:lmux,i) = 0

C     pvagm1: augm. f1~f2~ part of sig and corresponding tau,ppi0
C     for orbitals constructed out of (u,s)
        call pvagm1(nf1,nf1s,lmx1,lx1,nlx1,v1,d1,i,
     .  nf2,nf2s,lmx2,lx2,nlx2,v2,d2,lso,
     .  hab(1,1,i),vab(1,1,i),sab(1,1,i),
     .  sodb(1,1,i,1),sodb(1,1,i,2),lmux,
     .  sig(1,1,0,i),tau(1,1,0,i),ppi0,
     .  hso(1,1,1,1,i,1),hso(1,1,1,1,i,2))
C     pvagm1: augm. f1~f2~ part of sig and corresponding tau,ppi0
C     for local orbitals
        call pvaglc(nf1,nf1s,lmx1,lx1,nlx1,v1,d1,i,
     .  nf2,nf2s,lmx2,lx2,nlx2,v2,d2,lso,
     .  hab(1,1,i),vab(1,1,i),sab(1,1,i),
     .  sodb(1,1,i,1),sodb(1,1,i,2),
     .  lmux,sig(1,1,0,i),tau(1,1,0,i),ppi0,
     .  hso(1,1,1,1,i,1),hso(1,1,1,1,i,2))

C     Symmetrize tau,sig,ppi0 (debugging)
C     call pvagm4(nf1,lx1,nf2,lx2,tau(1,1,0,i))
C     call pvagm4(nf1,lx1,nf2,lx2,sig(1,1,0,i))
C     call pvagm4(nf1,lx1,nf2,lx2,ppi0)

C --- Contribution to ppi from non-spherical potential ---
        call dpzero(ppi(1,1,1,1,i,1),nf1*nf2*nlx1*nlx2)
C     paug2: smooth integral f1^ (-vsm) f2^ for nonspherical part of vsm
C     Also include local orbitals; require f=0 if no smooth part.
        call paug2(nr,nlml,vsm(1,1,i),rwgt,cg,jcg,indxcg,nf1,nf1,lmx1,
     .  lx1,f1,nf2,nf2,lmx2,lx2,f2,sum,nlx1,nlx2,ppi(1,1,1,1,i,1))
C     print *, '!! zero ppi (sm)'
C     ppi(:,:,:,:,i) = 0
C     paug1: augm integral f1~ vtrue f2~ for nonspherical part of vtrue
        call paug1(nf1,nf1s,lmx1,lx1,v1,d1,nf2,nf2s,lmx2,lx2,v2,d2,lmxa,
     .  nlml,cg,jcg,indxcg,vum(0,0,1,1,i),nlx1,nlx2,ppi(1,1,1,1,i,1))

C --- Moments qm = (f1~*f2~ - f1^*f2^) r^m Y_m ---
C     Needed for the term qm * (gpot0-gpotb) added to ppi
        call pvagm3(nf1,nf1s,lmx1,lx1,f1,v1,d1,nf2,nf2s,lmx2,lx2,f2,v2,d2,
     .  nr,rofi,rwgt,lmxa,qum(0,0,0,1,i),lmxl,qm)
C     print *, '!! qm=0'
C     qm=0

C --- Add any LDA+U contribution that exists
        if (lldau .ne. 0) then
          call paugnl(nf1,nf1s,lmx1,lx1,v1,d1,nf2,nf2s,lmx2,lx2,v2,d2,
     .    lmaxu,vumm,nlx1,nlx2,ppiz(1,1,1,1,i,1),i,idu)
        endif

C --- Assemble ppi from ppi0, non-spherical part and multipole contr ---
        call paug3(nf1,lmx1,lx1,nf2,lmx2,lx2,lmxl,nlml,cg,jcg,
     .  indxcg,qm,gpotb,gpot0,lmux,ppi0,nlx1,nlx2,ppi(1,1,1,1,i,1))

C     if (kount .eq. 16) call snut(nf1,nlx1,ppi)

C --- Add tau into ppi ---
        do  i1 = 1, nf1
          do  i2 = 1, nf2
            nlm1 = (lx1(i1)+1)**2
            nlm2 = (lx2(i2)+1)**2
            nlm = min0(nlm1,nlm2)
            do  ilm = 1, nlm
              l = ll(ilm)
              ppi(i1,i2,ilm,ilm,i,1) = ppi(i1,i2,ilm,ilm,i,1)
     .        + tau(i1,i2,l,i)
            enddo
          enddo
        enddo

C     call pshpr(70)
C     stdo = 88

#if DEBUG
        call info2(0,0,0,
     .  ' GAUGM dumping ppi: %ix%i radial functions',nf1,nf2)
        i1 = fopna('ppi',-1,4)
        write(i1) ppi
#endif

C --- Print diagonal augmentation matrices ---
        if (iprint() .ge. 70) then
          write(stdo,501) i
          do  i1 = 1, nf1
            do  i2 = 1, nf2
              nlm1 = (lx1(i1)+1)**2
              nlm2 = (lx2(i2)+1)**2
              nlm = min0(nlm1,nlm2)
              write(stdo,*) ' '
              do  ilm = 1, nlm
                l=ll(ilm)
                write(stdo,500) i1,i2,ilm,sig(i1,i2,l,i),tau(i1,i2,l,i),
     .          ppi(i1,i2,ilm,ilm,i,1)
  500           format(2i4,i6,3f14.8)
  501           format(/'# diagonal aug matrices spin',i2,':'
     .          /'# i1  i2   ilm',8x,'sig',11x,'tau',11x,'ham')
              enddo
            enddo
          enddo
        endif

C ... Print LL' potential matrix (which includes tau)
        if (iprint() .ge. 80) then
          write(stdo,'(''#pi matrix:'')')
          write(stdo,'(''#i1 i2  ilm1 ilm2     ham'')')
          do  i1 = 1, nf1
            do  i2 = 1, nf2
              nlm1 = (lx1(i1)+1)**2
              nlm2 = (lx2(i2)+1)**2
              do  ilm1 = 1, nlm1
                do  ilm2 = 1, nlm2
                  if (dabs(ppi(i1,i2,ilm1,ilm2,i,1)) .ge. 1d-8)
     .            write(stdo,320) i1,i2,ilm1,ilm2,ppi(i1,i2,ilm1,ilm2,i,1)
C 320         format(2i3,2i5,f15.10)
  320             format(2i3,2i5,f14.8)

                enddo
              enddo
            enddo
          enddo
        endif

      enddo

C --- Copy to complex potential matrix when needed ---
      if (lso .ne. 0) then
        if (lso .eq. 3) call rx('fix dim of ppi to add so')
C       1. Put (true) hso into (complex version of) ppi(:,:,:,:,:,:,2)
C       2. convert ppiz(1,1,1,1,1:2,1,1) into complex form
C       3. Add hso=ppiz(:,:,:,:,:,:,2) to ppiz(:,:,:,:,:,:,1)
        call ppi2z(12,nsp,nspc,nf1,nf2,nlx1,nlx2,hso,ppi,ppi)
C       i1 = nf1*nf2*nlx1*nlx2
C       call zprm('ppiz',2,ppi,i1,i1,4)
      else if (lcplxp .ne. 0) then
        call ppi2z(1,nsp,nspc,nf1,nf2,nlx1,nlx2,ppiz,ppi,ppi)
C       call zprm('ppiz',2,ppi,nf1*nf2,nf1*nf2,nlx1*nlx2)
C     else
C       call prmx('ppi',ppi,nf1*nf2,nf1*nf2,nlx1*nlx2*nsp*nspc)
      endif

      if (lldau .ne. 0) then
        call ppi2z(7,nsp,nspc,nf1,nf2,nlx1,nlx2,ppiz,ppi,ppi)
C       call zprm('ppiz',2,ppi,nf1*nf2,nf1*nf2,nlx1*nlx2*nsp*nspc)
      endif

C     call poppr
C     pause

CSFCPP#if F90
      deallocate(hso,ppiz)
CSFCPP#endif

      end

      subroutine pvagm1(nf1,nf1s,lmx1,lx1,nlx1,v1,d1,isp,
     .nf2,nf2s,lmx2,lx2,nlx2,v2,d2,lso,
     .hab,vab,sab,sodb,sondb,lmux,sig,tau,ppi,
     .hsozz,hsopm)
C- Add augmentation part of sig and tau and ppi (spherical pot)
C ----------------------------------------------------------------------
Ci Inputs
Ci   nf1   :number of 'bra' function types for each l
Ci   nf1s  :number of 'bra' function types for each l
Ci         :which are formed by linear combination of (u,s)
Ci   lmx1  :dimensions f1
Ci   lx1   :l-cutoffs for each of the nf1 functions
Ci   f1    :'bra' unaugmented envelope fn (radial part . r).
Ci         :Must be zero in channels that have no envelope functions.
Ci         :See Remarks
Ci   v1    :values of f1 at rofi(nr) (not multiplied by r)
Ci   d1    :slopes of f1 at rofi(nr) (not multiplied by r)
Ci   nf2   :number of 'ket' function types for each l
Ci   nf2s  :number of 'ket' function types for each l
Ci         :which are formed by linear combination of (u,s)
Ci   lmx2  :dimensions f2
Ci   lx2   :l-cutoffs for each of the nf2 functions
Ci   f2    :'ket' unaugmented envelope fn (radial part . r).
Ci         :Must be zero in channels that have no envelope functions
Ci         :See Remarks.
Ci   v2    :values of f2 at rofi(nr) (not multiplied by r)
Ci   d2    :slopes of f2 at rofi(nr) (not multiplied by r)
Ci   isp   :current spin index; used for SO coupling
Ci   lso   :if nonzero calculate hsozz = LzSz matrix elements
Ci         :if 1 calculate also  hsopm
Ci   hab   :hamiltonian matrix elements of (u,s); see Remarks
Ci   vab   :potential matrix elements of (u,s)
Ci   sab   :overlap matrix elements of (u,s)
Ci   sodb  :corresponding matrix elements of SO = 2/(c^2) dV/dr*(1/r)
Ci         :see comments for routine gaugm
Ci   sondb :corresponding matrix elements of SO = 2/(c^2) dV/dr*(1/r)
Ci         :for off-diagonal part; see comments for routine gaugm
Ci   lmux  :l-cutoff for sig,tau,ppi
Ci         :Usually min (largest lmx1, largest lmx2)
Cio Inputs/Outputs
Cio  sig   :On input, sig is smoothed counterpart of overlap
Cio        :  matrix <f1^ f2^> (see pvagm2)
Cio        :On output, sig is overwritten by
Cio        :  <f1~ f2~> - <f1^ f2^>
Cio        :  where the first term is calculated from sab
Cio  tau   :On input, tau is smoothed counterpart of k.e.
Cio        :  matrix <f1^ -nabla f2^> (see pvagm2)
Cio        :On output, tau is overwritten by
Cio        :  <f1~ -nabla f2~> - <f1^ -nabla f2^>
Cio        :  where the first term is calculated from hab-vab
Cio  ppi   :On input, ppi is smoothed counterpart of pot
Cio        :  matrix <f1^ v2~ f2^> (see pvagm2)
Cio        :On output, ppi is overwritten by
Cio        :  <f1~ v1(l=0) f2~> - <f1^ v2~(l=0) f2^>
Cio        :  where the first term is calculated from vab
Co Outputs
Co  hsozz  :spin-orbit contribution to ppi, spin diagonal block
Co         :(LzSz part)
Co         :isp=1: makes up-up block
Co         :isp=2: makes down-down block
Co  hsopm  :spin-orbit contribution to ppi, spin off-diagonal block
Co         :(LxSx + LySy part)
Co         :isp=1: makes up-down block
Co         :isp=2: makes down-up block
Cu Updates
Cu   03 Feb 05 (A. Chantis) calculate matrix elements of L.S
Cr Remarks
Cr   See Remarks for routine gaugm, above.
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer lmux,lmx1,lmx2,nf1,nf2,nf1s,nf2s,lx1(nf1),lx2(nf2),n0,nab,
     .nlx1,nlx2,isp
      parameter (n0=10,nab=9)
      double precision hab(nab,0:n0-1),sab(nab,0:n0-1),vab(nab,0:n0-1),
     .v1(0:lmx1,nf1),d1(0:lmx1,nf1),
     .v2(0:lmx2,nf2),d2(0:lmx2,nf2),
     .sig(nf1,nf2,0:lmux),tau(nf1,nf2,0:lmux),ppi(nf1,nf2,0:lmux)
C     Spin-Orbit related
      double precision sodb(nab,0:n0-1),sondb(nab,0:n0-1)
      double precision tmp(nf1,nf2,nlx1,nlx2),a1,a2
      double complex hsozz(nf1,nf2,nlx1,nlx2),hsopm(nf1,nf2,nlx1,nlx2)
C ... Local parameters
      integer i1,i2,lmax1,lmax2,lmax,l
C     Spin-Orbit related
      integer m1,m2,l1,l2,lso
      double complex img

      if (lso .ne. 0) then
        img = dcmplx(0d0,1d0)
        call dpzero(tmp,   nf1*nf2*nlx1*nlx2)
        call dpzero(hsopm, nf1*nf2*nlx1*nlx2*2)
        call dpzero(hsozz, nf1*nf2*nlx1*nlx2*2)
      endif

      do  i1 = 1, nf1s
        do  i2 = 1, nf2s
          lmax1 = lx1(i1)
          lmax2 = lx2(i2)
          lmax = min0(lmax1,lmax2)
          do  l = 0, lmax

            sig(i1,i2,l) = -sig(i1,i2,l)
     .      + v1(l,i1)*sab(1,l)*v2(l,i2)
     .      + v1(l,i1)*sab(2,l)*d2(l,i2)
     .      + d1(l,i1)*sab(3,l)*v2(l,i2)
     .      + d1(l,i1)*sab(4,l)*d2(l,i2)

            tau(i1,i2,l) = -tau(i1,i2,l)
     .      + v1(l,i1)*(hab(1,l)-vab(1,l))*v2(l,i2)
     .      + v1(l,i1)*(hab(2,l)-vab(2,l))*d2(l,i2)
     .      + d1(l,i1)*(hab(3,l)-vab(3,l))*v2(l,i2)
     .      + d1(l,i1)*(hab(4,l)-vab(4,l))*d2(l,i2)

            ppi(i1,i2,l) = -ppi(i1,i2,l)
     .      + v1(l,i1)*vab(1,l)*v2(l,i2)
     .      + v1(l,i1)*vab(2,l)*d2(l,i2)
     .      + d1(l,i1)*vab(3,l)*v2(l,i2)
     .      + d1(l,i1)*vab(4,l)*d2(l,i2)

          enddo
        enddo

      enddo

C ... Spin-Orbit matrix elements in real harmonics.
C     This is LxSx+LySy part
      if (lso .eq. 1) then
        do  i1 = 1, nf1s
          do  i2 = 1, nf2s
            lmax1 = lx1(i1)
            lmax2 = lx2(i2)
            lmax = min0(lmax1,lmax2)
            l1 = 0
            l2 = 0
            do  l = 0, lmax
              do  m1 = -l, l
                l1 = l1 + 1
                if (m1 .ge. (-l+1)) l2 = l2 - (2*l + 1)
                do  m2 = -l, l

                  l2 = l2 + 1
                  a1 = dsqrt(dble((l-abs(m2))*(l+abs(m2)+1)))
                  a2 = dsqrt(dble((l+abs(m2))*(l-abs(m2)+1)))
                  tmp(i1,i2,l1,l2) =
     .            (v1(l,i1)*sondb(1,l)*v2(l,i2)
     .            + v1(l,i1)*sondb(2,l)*d2(l,i2)
     .            + d1(l,i1)*sondb(3,l)*v2(l,i2)
     .            + d1(l,i1)*sondb(4,l)*d2(l,i2))

C         ... Spin up-down block <l,m|L-|l,m'>
                  if (isp .eq. 1) then

C             Case A
                    if (abs(m2) .gt. 1 .and. (abs(m2)+1) .le. l) then

                      if (m2 .gt. 0) then
                        if ((abs(m2)-1) .eq. m1) hsopm(i1,i2,l1,l2) =
     .                  (-1)**(2*m2-1)*a2*0.5d0*tmp(i1,i2,l1,l2)

                        if ((abs(m2)-1) .eq. -m1) hsopm(i1,i2,l1,l2) =
     .                  img*(-1)**(2*m2-1)*a2*0.5d0*tmp(i1,i2,l1,l2)

                        if ((abs(m2)+1) .eq. m1) hsopm(i1,i2,l1,l2) =
     .                  a1*0.5d0*tmp(i1,i2,l1,l2)

                        if ((abs(m2)+1) .eq. -m1) hsopm(i1,i2,l1,l2) =
     .                  -img*a1*0.5d0*tmp(i1,i2,l1,l2)

                      else

                        if ((abs(m2)-1) .eq. m1) hsopm(i1,i2,l1,l2) =
     .                  -img*(-1)**(2*m2-1)*a2*0.5d0*tmp(i1,i2,l1,l2)

                        if ((abs(m2)-1) .eq. -m1) hsopm(i1,i2,l1,l2) =
     .                  (-1)**(2*m2-1)*a2*0.5d0*tmp(i1,i2,l1,l2)

                        if ((abs(m2)+1) .eq. m1) hsopm(i1,i2,l1,l2) =
     .                  img*a1*0.5d0*tmp(i1,i2,l1,l2)

                        if ((abs(m2)+1) .eq. -m1) hsopm(i1,i2,l1,l2) =
     .                  a1*0.5d0*tmp(i1,i2,l1,l2)
                      endif

                    endif

C             Case B
                    if (abs(m2) .gt. 1 .and. (abs(m2)+1) .gt. l) then

                      if (m2 .gt. 0) then
                        if ((abs(m2)-1) .eq. m1) hsopm(i1,i2,l1,l2) =
     .                  (-1)**(2*m2+1)*a2*0.5d0*tmp(i1,i2,l1,l2)
                        if ((abs(m2)-1) .eq. -m1) hsopm(i1,i2,l1,l2) =
     .                  img*(-1)**(2*m2+1)*a2*0.5d0*tmp(i1,i2,l1,l2)
                      endif

                      if (m2 .lt. 0) then
                        if ((abs(m2)-1) .eq. m1) hsopm(i1,i2,l1,l2) =
     .                  -img*(-1)**(2*m2-1)*a2*0.5d0*tmp(i1,i2,l1,l2)
                        if ((abs(m2)-1) .eq. -m1) hsopm(i1,i2,l1,l2) =
     .                  (-1)**(2*m2-1)*a2*0.5d0*tmp(i1,i2,l1,l2)
                      endif

                    endif

C             Case C
                    if (abs(m2) .eq. 1 .and. (abs(m2)+1) .le. l) then
                      if (m2 .gt. 0) then
                        if (m1 .eq. 0) hsopm(i1,i2,l1,l2) =
     .                  (-1)**m2*a2*dsqrt(0.5d0)*tmp(i1,i2,l1,l2)
                        if ((abs(m2)+1) .eq. m1) hsopm(i1,i2,l1,l2) =
     .                  a1*0.5d0*tmp(i1,i2,l1,l2)
                        if ((abs(m2)+1) .eq. -m1) hsopm(i1,i2,l1,l2) =
     .                  -img*a1*0.5d0*tmp(i1,i2,l1,l2)
                      endif

                      if (m2 .lt. 0) then
                        if (m1 .eq. 0) hsopm(i1,i2,l1,l2) =
     .                  -img*(-1)**m2*a2*dsqrt(0.5d0)*tmp(i1,i2,l1,l2)
                        if ((abs(m2)+1) .eq. m1) hsopm(i1,i2,l1,l2) =
     .                  img*a1*0.5d0*tmp(i1,i2,l1,l2)
                        if ((abs(m2)+1) .eq. -m1) hsopm(i1,i2,l1,l2) =
     .                  a1*0.5d0*tmp(i1,i2,l1,l2)
                      endif

                    endif

C             Case D
                    if (abs(m2) .eq. 1 .and. (abs(m2)+1) .gt. l) then
                      if (m2 .gt. 0) then
                        if (m1 .eq. 0) hsopm(i1,i2,l1,l2) =
     .                  (-1)**m2*a2*dsqrt(0.5d0)*tmp(i1,i2,l1,l2)
                      endif

                      if (m2 .lt. 0) then
                        if (m1 .eq. 0) hsopm(i1,i2,l1,l2) =
     .                  -img*(-1)**m2*a2*dsqrt(0.5d0)*tmp(i1,i2,l1,l2)
                      endif
                    endif

C             Case m=0
                    if (abs(m2) .eq. 0) then
                      if (m1 .eq. 1) hsopm(i1,i2,l1,l2) =
     .                a1*dsqrt(0.5d0)*tmp(i1,i2,l1,l2)
                      if (m1 .eq. -1) hsopm(i1,i2,l1,l2) =
     .                -img*a1*dsqrt(0.5d0)*tmp(i1,i2,l1,l2)
                    endif

C         ... Spin down-up block <l,m|L+|l,m'>
                  else
C               Case A
                    if (abs(m2) .gt. 1 .and. (abs(m2)+1) .le. l) then

                      if (m2 .gt. 0) then
                        if ((abs(m2)-1) .eq. m1) hsopm(i1,i2,l1,l2) =
     .                  a2*0.5d0*tmp(i1,i2,l1,l2)

                        if ((abs(m2)-1) .eq. -m1) hsopm(i1,i2,l1,l2) =
     .                  -img*a2*0.5d0*tmp(i1,i2,l1,l2)

                        if ((abs(m2)+1) .eq. m1) hsopm(i1,i2,l1,l2) =
     .                  (-1)**(2*m2+1)*a1*0.5d0*tmp(i1,i2,l1,l2)

                        if ((abs(m2)+1) .eq. -m1) hsopm(i1,i2,l1,l2) =
     .                  img*(-1)**(2*m2+1)*a1*0.5d0*tmp(i1,i2,l1,l2)

                      else

                        if ((abs(m2)-1) .eq. m1) hsopm(i1,i2,l1,l2) =
     .                  img*a2*0.5d0*tmp(i1,i2,l1,l2)

                        if ((abs(m2)-1) .eq. -m1) hsopm(i1,i2,l1,l2) =
     .                  a2*0.5d0*tmp(i1,i2,l1,l2)

                        if ((abs(m2)+1) .eq. m1) hsopm(i1,i2,l1,l2) =
     .                  -img*a1*(-1)**(2*m2+1)*0.5d0*tmp(i1,i2,l1,l2)

                        if ((abs(m2)+1) .eq. -m1) hsopm(i1,i2,l1,l2) =
     .                  (-1)**(2*m2+1)*a1*0.5d0*tmp(i1,i2,l1,l2)
                      endif

                    endif

C               Case B
                    if (abs(m2) .gt. 1 .and. (abs(m2)+1) .gt. l) then

                      if (m2 .gt. 0) then
                        if ((abs(m2)-1) .eq. m1) hsopm(i1,i2,l1,l2) =
     .                  a2*0.5d0*tmp(i1,i2,l1,l2)
                        if ((abs(m2)-1) .eq. -m1) hsopm(i1,i2,l1,l2) =
     .                  -img*a2*0.5d0*tmp(i1,i2,l1,l2)
                      endif

                      if (m2 .lt. 0) then
                        if ((abs(m2)-1) .eq. m1) hsopm(i1,i2,l1,l2) =
     .                  img*a2*0.5d0*tmp(i1,i2,l1,l2)
                        if ((abs(m2)-1) .eq. -m1) hsopm(i1,i2,l1,l2) =
     .                  a2*0.5d0*tmp(i1,i2,l1,l2)
                      endif

                    endif

C               Case C
                    if (abs(m2) .eq. 1 .and. (abs(m2)+1) .le. l) then
                      if (m2 .gt. 0) then
                        if (m1 .eq. 0) hsopm(i1,i2,l1,l2) =
     .                  a2*dsqrt(0.5d0)*tmp(i1,i2,l1,l2)
                        if ((abs(m2)+1) .eq. m1) hsopm(i1,i2,l1,l2) =
     .                  (-1)**(2*m2+1)*a1*0.5d0*tmp(i1,i2,l1,l2)
                        if ((abs(m2)+1) .eq. -m1) hsopm(i1,i2,l1,l2) =
     .                  img*(-1)**(2*m2+1)*a1*0.5d0*tmp(i1,i2,l1,l2)
                      endif

                      if (m2 .lt. 0) then
                        if (m1 .eq. 0) hsopm(i1,i2,l1,l2) =
     .                  img*a2*dsqrt(0.5d0)*tmp(i1,i2,l1,l2)
                        if ((abs(m2)+1) .eq. m1) hsopm(i1,i2,l1,l2) =
     .                  -img*(-1)**(2*m2+1)*a1*0.5d0*tmp(i1,i2,l1,l2)
                        if ((abs(m2)+1) .eq. -m1) hsopm(i1,i2,l1,l2) =
     .                  (-1)**(2*m2+1)*a1*0.5d0*tmp(i1,i2,l1,l2)
                      endif

                    endif

C               Case D
                    if (abs(m2) .eq. 1 .and. (abs(m2)+1) .gt. l) then
                      if (m2 .gt. 0) then
                        if (m1 .eq. 0) hsopm(i1,i2,l1,l2) =
     .                  a2*dsqrt(0.5d0)*tmp(i1,i2,l1,l2)
                      endif

                      if (m2 .lt. 0) then
                        if (m1 .eq. 0) hsopm(i1,i2,l1,l2) =
     .                  img*a2*dsqrt(0.5d0)*tmp(i1,i2,l1,l2)
                      endif
                    endif

C               Case m=0
                    if (abs(m2) .eq. 0) then
                      if (m1 .eq. 1) hsopm(i1,i2,l1,l2) =
     .                -a1*dsqrt(0.5d0)*tmp(i1,i2,l1,l2)
                      if (m1 .eq. -1) hsopm(i1,i2,l1,l2) =
     .                -img*a1*dsqrt(0.5d0)*tmp(i1,i2,l1,l2)
                    endif

                  endif


                enddo
C         Ends loop over L
              enddo
            enddo
C       Ends loop over orbital pairs
          enddo
        enddo

      endif

C ... LzSz part
      if (lso .ne. 0) then

        do  i1 = 1, nf1s
          do  i2 = 1, nf2s
            lmax1 = lx1(i1)
            lmax2 = lx2(i2)
            lmax = min0(lmax1,lmax2)
            l1 = 0
            l2 = 0
            do  l = 0, lmax
              do  m1 = -l, l
                l1 = l1 + 1
                if (m1 .ge. (-l+1)) l2 = l2 - (2*l + 1)
                do  m2 = -l, l
                  l2 = l2 + 1

                  if (m1 .ne. 0 .and. m2 .ne. 0) then
                    if (l1 .lt. l2 .and. m1 .eq. -m2) then
                      hsozz(i1,i2,l1,l2) = abs(m1)*
     .                (v1(l,i1)*dcmplx(0d0,sodb(1,l))*v2(l,i2)
     .                + v1(l,i1)*dcmplx(0d0,sodb(2,l))*d2(l,i2)
     .                + d1(l,i1)*dcmplx(0d0,sodb(3,l))*v2(l,i2)
     .                + d1(l,i1)*dcmplx(0d0,sodb(4,l))*d2(l,i2))
                    endif
                    if (l1 .gt. l2 .and. m1 .eq. -m2) then
                      hsozz(i1,i2,l1,l2) = -abs(m1)*
     .                (v1(l,i1)*dcmplx(0d0,sodb(1,l))*v2(l,i2)
     .                + v1(l,i1)*dcmplx(0d0,sodb(2,l))*d2(l,i2)
     .                + d1(l,i1)*dcmplx(0d0,sodb(3,l))*v2(l,i2)
     .                + d1(l,i1)*dcmplx(0d0,sodb(4,l))*d2(l,i2))
                    endif
                  endif

                enddo
              enddo
            enddo
          enddo
        enddo

      endif

C      if (lso .ne. 0) then
C        print *, '!! zero out hsopm'
C        call dpzero(hsopm, nf1*nf2*nlx1*nlx2*2)
CC       call dpzero(hsozz, nf1*nf2*nlx1*nlx2*2)
C      endif


C      print *, 'pvaug1'
C      do  i1 = 1, nf1s
C        do  i2 = 1, nf2s
C          lmax1 = lx1(i1)
C          lmax2 = lx2(i2)
C          lmax = min0(lmax1,lmax2)
C          write(*,*) ' '
C          do  l = 0, lmax
C            write(*,500) i1,i2,l,sig(i1,i2,l),tau(i1,i2,l)
C  500       format(2i4,i6,3f14.8)
C          enddo
C        enddo
C      enddo
C      stop

      end

      subroutine pvaglc(nf1,nf1s,lmx1,lx1,nlx1,v1,d1,isp,
     .nf2,nf2s,lmx2,lx2,nlx2,v2,d2,lso,
     .hab,vab,sab,sodb,sondb,
     .lmux,sig,tau,ppi,hsozz,hsopm)
C- Augmentation part of sig and tau and ppi, loc. orbitals (spher. pot)
C ----------------------------------------------------------------------
Ci Inputs
Ci   nf1   :number of 'bra' function types for each l
Ci   nf1s  :number of 'bra' function types for each l
Ci         :which are formed by linear combination of (u,s)
Ci   lmx1  :dimensions f1
Ci   lx1   :l-cutoffs for each of the nf1 functions
Ci   nlx1  :dimensions hsozz,hsopm
Ci   v1    :values of f1 at rofi(nr) (not multiplied by r)
Ci   d1    :slopes of f1 at rofi(nr) (not multiplied by r)
Ci   isp   :(used only when adding to SO hamiltonian)
Ci         :isp=1 => make <l,m|L-|l,m'> of hsopm
Ci         :isp=2 => make <l,m|L+|l,m'> of hsopm
Ci   nf2   :number of 'ket' function types for each l
Ci   nf2s  :number of 'ket' function types for each l
Ci         :which are formed by linear combination of (u,s)
Ci   lmx2  :dimensions f2
Ci   lx2   :l-cutoffs for each of the nf2 functions
Ci   nlx2  :dimensions hsozz,hsopm
Ci   v2    :values of f2 at rofi(nr) (not multiplied by r)
Ci   d2    :slopes of f2 at rofi(nr) (not multiplied by r)
Ci   lso   :1 include L.S coupling; 2 include LzSz part only
Ci   hab   :hamiltonian matrix elements of (u,s); see Remarks
Ci         :Leading dimension of hab corresponds to these entries:
Ci         :(1)=<ul|h|ul>   (5)=<ul|h|gz>   (9)=<gz|h|sl>
Ci         :(2)=<ul|h|sl>   (6)=<sl|h|gz>
Ci         :(3)=<sl|h|ul>   (7)=<gz|h|gz>
Ci         :(4)=<sl|h|sl>   (8)=<gz|h|ul>
Ci   vab   :potential matrix elements of (u,s)
Ci   sab   :overlap matrix elements of (u,s)
Ci   sodb
Ci   sondb
Ci   lmux  :l-cutoff for sig,tau,ppi
Ci         :Usually min (largest lmx1, largest lmx2)
Cio Inputs/Outputs
Cio  sig   :On input, sig is smoothed counterpart of overlap
Cio        :  matrix <f1^ f2^> (see pvagm2)
Cio        :On output, sig is overwritten by
Cio        :  <f1~ f2~> - <f1^ f2^>
Cio        :  where the first term is calculated from sab
Cio  tau   :On input, tau is smoothed counterpart of k.e.
Cio        :  matrix <f1^ -nabla f2^> (see pvagm2)
Cio        :On output, tau is overwritten by
Cio        :  <f1~ -nabla f2~> - <f1^ -nabla f2^>
Cio        :  where the first term is calculated from hab-vab
Cio  ppi   :On input, ppi is smoothed counterpart of pot
Cio        :  matrix <f1^ v2~ f2^> (see pvagm2)
Cio        :On output, ppi is overwritten by
Cio        :  <f1~ v1(l=0) f2~> - <f1^ v2~(l=0) f2^>
Cio        :  where the first term is calculated from vab
Cio  hsozz :
Cio  hsopm :
Cr Remarks
Cr   See Remarks for routine gaugm, above.
Cr
Cr   Symmetrization of kinetic energy.  The kinetic energy matrix
Cr   should always be symmetric.  In general the
Cr   difference T_fg - T_gf between orbitals f and g is
Cr    T_fg - T_gf = -[r^2 (fg' - f'g)]_rmax = -W{f,g;rmax}
Cr   which is nonzero.  However, the kinetic energy matrix always
Cr   takes of two forms:
Cr      1) it consists of a difference T_fg = <f T g> - <f^ T g^>
Cr         between true and smoothed functions where f matches f^
Cr         and g matches g^ at rmax.  In this case, T_fg - T_gf = 0
Cr         because W{f,g;rmax} = W{f^,g^;rmax}
Cr      2) either f or g is a local orbital with value,slope=0.
Cr         Then -W{f,g;rmax} = 0.
Cr
Cr   See also Remarks for routine gaugm, above.
Cu Updates
Cu   26 Sep 05 (A. Chantis) added matrix elements for SO coupling
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer lmux,lmx1,lmx2,nf1,nf2,nf1s,nf2s,lx1(nf1),lx2(nf2),n0,nab,
     .isp,nlx1,nlx2
      parameter (n0=10,nab=9)
      double precision hab(nab,0:n0-1),sab(nab,0:n0-1),vab(nab,0:n0-1),
     .v1(0:lmx1,nf1),d1(0:lmx1,nf1),
     .v2(0:lmx2,nf2),d2(0:lmx2,nf2),
     .sig(nf1,nf2,0:lmux),tau(nf1,nf2,0:lmux),ppi(nf1,nf2,0:lmux)
C     Spin-Orbit related
      double precision sodb(nab,0:n0-1),sondb(nab,0:n0-1)
      double complex hsozz(nf1,nf2,nlx1,nlx2),hsopm(nf1,nf2,nlx1,nlx2)
C ... Local parameters
      integer i1,i2,lmax1,lmax2,lmax,l
C     Spin-Orbit related
      integer m1,m2,l1,l2,lso
      double precision tmp(nf1,nf2,nlx1,nlx2),tmp1(nf1,nf2,nlx1,nlx2),
     .a1,a2
      double complex img

      if (nf1s .ge. nf1 .and. nf2s .ge. nf2) return

      if (lso .ne. 0) then
        img = dcmplx(0d0,1d0)
        call dpzero(tmp,   nf1*nf2*nlx1*nlx2)
        call dpzero(tmp1,  nf1*nf2*nlx1*nlx2)
      endif

      do  i1 = 1, nf1
        do  i2 = 1, nf2
          if (i1 .gt. nf1s .or. i2 .gt. nf2s) then
            lmax1 = lx1(i1)
            lmax2 = lx2(i2)
            lmax = min0(lmax1,lmax2)
            l1 = 0
            l2 = 0
            do  l = 0, lmax

C            print *, '!! zero'
C            v2(l,i2) = 0
C            d2(l,i2) = 1
C            v1(l,i1) = 0
C            d1(l,i1) = 1

C       ... hzz,vzz,szz
              if (i1 .gt. nf1s .and. i2 .gt. nf2s) then
                sig(i1,i2,l) = -sig(i1,i2,l) + sab(7,l)
                tau(i1,i2,l) = -tau(i1,i2,l) + (hab(7,l)-vab(7,l))
                ppi(i1,i2,l) = -ppi(i1,i2,l) + vab(7,l)

C       ... hzu,vzu,szs
              elseif (i1 .gt. nf1s) then
                sig(i1,i2,l) = -sig(i1,i2,l)
     .          + sab(8,l)*v2(l,i2)
     .          + sab(9,l)*d2(l,i2)
                tau(i1,i2,l) = -tau(i1,i2,l)
     .          + (hab(8,l)-vab(8,l))*v2(l,i2)
     .          + (hab(9,l)-vab(9,l))*d2(l,i2)
                ppi(i1,i2,l) = -ppi(i1,i2,l)
     .          + vab(8,l)*v2(l,i2)
     .          + vab(9,l)*d2(l,i2)

C       ... huz,vuz,ssz
              elseif (i2 .gt. nf2s) then
                sig(i1,i2,l) = -sig(i1,i2,l)
     .          + v1(l,i1)*sab(5,l)
     .          + d1(l,i1)*sab(6,l)
                tau(i1,i2,l) = -tau(i1,i2,l)
     .          + v1(l,i1)*(hab(5,l)-vab(5,l))
     .          + d1(l,i1)*(hab(6,l)-vab(6,l))
                ppi(i1,i2,l) = -ppi(i1,i2,l)
     .          + v1(l,i1)*vab(5,l)
     .          + d1(l,i1)*vab(6,l)
              endif

C       ... Spin-Orbit matrix elements in real harmonics.
              if (lso .ne. 0) then
                do  m1 = -l, l
                  l1 = l1 + 1
                  if (m1 .ge. (-l+1)) l2 = l2 - (2*l + 1)
                  do  m2 = -l, l
                    l2 = l2 + 1
C             ... hso_zz
                    if (i1 .gt. nf1s .and. i2 .gt. nf2s) then
                      tmp1(i1,i2,l1,l2) = sodb(7,l)
                    elseif (i1 .gt. nf1s) then
C             ... hso_zu
                      tmp1(i1,i2,l1,l2) = sodb(8,l)*v2(l,i2)
     .                + sodb(9,l)*d2(l,i2)
                    elseif (i2 .gt. nf2s) then
C             ... hso_uz
                      tmp1(i1,i2,l1,l2) = v1(l,i1)*sodb(5,l)
     .                + d1(l,i1)*sodb(6,l)
                    endif
                    if (m1 .ne. 0 .and. m2 .ne. 0) then
                      if (l1 .lt. l2 .and. m1 .eq. -m2) then

                        hsozz(i1,i2,l1,l2) =  abs(m1)
     .                  *dcmplx(0d0,tmp1(i1,i2,l1,l2))
                      endif
                      if (l1 .gt. l2 .and. m1 .eq. -m2) then
                        hsozz(i1,i2,l1,l2) = -abs(m1)
     .                  *dcmplx(0d0,tmp1(i1,i2,l1,l2))
                      endif
                    endif

C       ... This is LxSx+LySy part
                    if (lso .eq. 1) then

                      a1 = dsqrt(dble((l-abs(m2))*(l+abs(m2)+1)))
                      a2 = dsqrt(dble((l+abs(m2))*(l-abs(m2)+1)))

C         ... hso_zz
                      if (i1 .gt. nf1s .and. i2 .gt. nf2s) then
                        tmp(i1,i2,l1,l2) = sondb(7,l)
                      elseif (i1 .gt. nf1s) then
C         ... hso_zu
                        tmp(i1,i2,l1,l2) = sondb(8,l)*v2(l,i2)
     .                  + sondb(9,l)*d2(l,i2)
                      elseif (i2 .gt. nf2s) then
C         ... hso_uz
                        tmp(i1,i2,l1,l2) = v1(l,i1)*sondb(5,l)
     .                  + d1(l,i1)*sondb(6,l)
                      endif

C         ... Spin up-down block <l,m|L-|l,m'>
                      if (isp .eq. 1) then

C               Case A
                        if (abs(m2) .gt. 1 .and. (abs(m2)+1) .le. l) then

                          if (m2 .gt. 0) then
                            if ((abs(m2)-1) .eq. m1) hsopm(i1,i2,l1,l2) =
     .                      (-1)**(2*m2-1)*a2*0.5d0*tmp(i1,i2,l1,l2)

                            if ((abs(m2)-1) .eq. -m1) hsopm(i1,i2,l1,l2) =
     .                      img*(-1)**(2*m2-1)*a2*0.5d0*tmp(i1,i2,l1,l2)

                            if ((abs(m2)+1) .eq. m1) hsopm(i1,i2,l1,l2) =
     .                      a1*0.5d0*tmp(i1,i2,l1,l2)

                            if ((abs(m2)+1) .eq. -m1) hsopm(i1,i2,l1,l2) =
     .                      -img*a1*0.5d0*tmp(i1,i2,l1,l2)

                          else

                            if ((abs(m2)-1) .eq. m1) hsopm(i1,i2,l1,l2) =
     .                      -img*(-1)**(2*m2-1)*a2*0.5d0*tmp(i1,i2,l1,l2)

                            if ((abs(m2)-1) .eq. -m1) hsopm(i1,i2,l1,l2) =
     .                      (-1)**(2*m2-1)*a2*0.5d0*tmp(i1,i2,l1,l2)

                            if ((abs(m2)+1) .eq. m1) hsopm(i1,i2,l1,l2) =
     .                      img*a1*0.5d0*tmp(i1,i2,l1,l2)

                            if ((abs(m2)+1) .eq. -m1) hsopm(i1,i2,l1,l2) =
     .                      a1*0.5d0*tmp(i1,i2,l1,l2)
                          endif

                        endif

C               Case B
                        if (abs(m2) .gt. 1 .and. (abs(m2)+1) .gt. l) then

                          if (m2 .gt. 0) then
                            if ((abs(m2)-1) .eq. m1) hsopm(i1,i2,l1,l2) =
     .                      (-1)**(2*m2+1)*a2*0.5d0*tmp(i1,i2,l1,l2)
                            if ((abs(m2)-1) .eq. -m1) hsopm(i1,i2,l1,l2) =
     .                      img*(-1)**(2*m2+1)*a2*0.5d0*tmp(i1,i2,l1,l2)
                          endif

                          if (m2 .lt. 0) then
                            if ((abs(m2)-1) .eq. m1) hsopm(i1,i2,l1,l2) =
     .                      -img*(-1)**(2*m2-1)*a2*0.5d0*tmp(i1,i2,l1,l2)
                            if ((abs(m2)-1) .eq. -m1) hsopm(i1,i2,l1,l2) =
     .                      (-1)**(2*m2-1)*a2*0.5d0*tmp(i1,i2,l1,l2)
                          endif

                        endif

C               Case C
                        if (abs(m2) .eq. 1 .and. (abs(m2)+1) .le. l) then
                          if (m2 .gt. 0) then
                            if (m1 .eq. 0) hsopm(i1,i2,l1,l2) =
     .                      (-1)**m2*a2*dsqrt(0.5d0)*tmp(i1,i2,l1,l2)
                            if ((abs(m2)+1) .eq. m1) hsopm(i1,i2,l1,l2) =
     .                      a1*0.5d0*tmp(i1,i2,l1,l2)
                            if ((abs(m2)+1) .eq. -m1) hsopm(i1,i2,l1,l2) =
     .                      -img*a1*0.5d0*tmp(i1,i2,l1,l2)
                          endif

                          if (m2 .lt. 0) then
                            if (m1 .eq. 0) hsopm(i1,i2,l1,l2) =
     .                      -img*(-1)**m2*a2*dsqrt(0.5d0)*tmp(i1,i2,l1,l2)
                            if ((abs(m2)+1) .eq. m1) hsopm(i1,i2,l1,l2) =
     .                      img*a1*0.5d0*tmp(i1,i2,l1,l2)
                            if ((abs(m2)+1) .eq. -m1) hsopm(i1,i2,l1,l2) =
     .                      a1*0.5d0*tmp(i1,i2,l1,l2)
                          endif

                        endif

C               Case D
                        if (abs(m2) .eq. 1 .and. (abs(m2)+1) .gt. l) then
                          if (m2 .gt. 0) then
                            if (m1 .eq. 0) hsopm(i1,i2,l1,l2) =
     .                      (-1)**m2*a2*dsqrt(0.5d0)*tmp(i1,i2,l1,l2)
                          endif

                          if (m2 .lt. 0) then
                            if (m1 .eq. 0) hsopm(i1,i2,l1,l2) =
     .                      -img*(-1)**m2*a2*dsqrt(0.5d0)*tmp(i1,i2,l1,l2)
                          endif
                        endif

C               Case m=0
                        if (abs(m2) .eq. 0) then
                          if (m1 .eq. 1) hsopm(i1,i2,l1,l2) =
     .                    a1*dsqrt(0.5d0)*tmp(i1,i2,l1,l2)
                          if (m1 .eq. -1) hsopm(i1,i2,l1,l2) =
     .                    -img*a1*dsqrt(0.5d0)*tmp(i1,i2,l1,l2)
                        endif

C         ... Spin down-up block <l,m|L+|l,m'>
                      else
C               Case A
                        if (abs(m2) .gt. 1 .and. (abs(m2)+1) .le. l) then

                          if (m2 .gt. 0) then
                            if ((abs(m2)-1) .eq. m1) hsopm(i1,i2,l1,l2) =
     .                      a2*0.5d0*tmp(i1,i2,l1,l2)

                            if ((abs(m2)-1) .eq. -m1) hsopm(i1,i2,l1,l2) =
     .                      -img*a2*0.5d0*tmp(i1,i2,l1,l2)

                            if ((abs(m2)+1) .eq. m1) hsopm(i1,i2,l1,l2) =
     .                      (-1)**(2*m2+1)*a1*0.5d0*tmp(i1,i2,l1,l2)

                            if ((abs(m2)+1) .eq. -m1) hsopm(i1,i2,l1,l2) =
     .                      img*(-1)**(2*m2+1)*a1*0.5d0*tmp(i1,i2,l1,l2)

                          else

                            if ((abs(m2)-1) .eq. m1) hsopm(i1,i2,l1,l2) =
     .                      img*a2*0.5d0*tmp(i1,i2,l1,l2)

                            if ((abs(m2)-1) .eq. -m1) hsopm(i1,i2,l1,l2) =
     .                      a2*0.5d0*tmp(i1,i2,l1,l2)

                            if ((abs(m2)+1) .eq. m1) hsopm(i1,i2,l1,l2) =
     .                      -img*a1*(-1)**(2*m2+1)*0.5d0*tmp(i1,i2,l1,l2)

                            if ((abs(m2)+1) .eq. -m1) hsopm(i1,i2,l1,l2) =
     .                      (-1)**(2*m2+1)*a1*0.5d0*tmp(i1,i2,l1,l2)
                          endif

                        endif

C               Case B
                        if (abs(m2) .gt. 1 .and. (abs(m2)+1) .gt. l) then

                          if (m2 .gt. 0) then
                            if ((abs(m2)-1) .eq. m1) hsopm(i1,i2,l1,l2) =
     .                      a2*0.5d0*tmp(i1,i2,l1,l2)
                            if ((abs(m2)-1) .eq. -m1) hsopm(i1,i2,l1,l2) =
     .                      -img*a2*0.5d0*tmp(i1,i2,l1,l2)
                          endif

                          if (m2 .lt. 0) then
                            if ((abs(m2)-1) .eq. m1) hsopm(i1,i2,l1,l2) =
     .                      img*a2*0.5d0*tmp(i1,i2,l1,l2)
                            if ((abs(m2)-1) .eq. -m1) hsopm(i1,i2,l1,l2) =
     .                      a2*0.5d0*tmp(i1,i2,l1,l2)
                          endif

                        endif

C               Case C
                        if (abs(m2) .eq. 1 .and. (abs(m2)+1) .le. l) then
                          if (m2 .gt. 0) then
                            if (m1 .eq. 0) hsopm(i1,i2,l1,l2) =
     .                      a2*dsqrt(0.5d0)*tmp(i1,i2,l1,l2)
                            if ((abs(m2)+1) .eq. m1) hsopm(i1,i2,l1,l2) =
     .                      (-1)**(2*m2+1)*a1*0.5d0*tmp(i1,i2,l1,l2)
                            if ((abs(m2)+1) .eq. -m1) hsopm(i1,i2,l1,l2) =
     .                      img*(-1)**(2*m2+1)*a1*0.5d0*tmp(i1,i2,l1,l2)
                          endif

                          if (m2 .lt. 0) then
                            if (m1 .eq. 0) hsopm(i1,i2,l1,l2) =
     .                      img*a2*dsqrt(0.5d0)*tmp(i1,i2,l1,l2)
                            if ((abs(m2)+1) .eq. m1) hsopm(i1,i2,l1,l2) =
     .                      -img*(-1)**(2*m2+1)*a1*0.5d0*tmp(i1,i2,l1,l2)
                            if ((abs(m2)+1) .eq. -m1) hsopm(i1,i2,l1,l2) =
     .                      (-1)**(2*m2+1)*a1*0.5d0*tmp(i1,i2,l1,l2)
                          endif

                        endif

C               Case D
                        if (abs(m2) .eq. 1 .and. (abs(m2)+1) .gt. l) then
                          if (m2 .gt. 0) then
                            if (m1 .eq. 0) hsopm(i1,i2,l1,l2) =
     .                      a2*dsqrt(0.5d0)*tmp(i1,i2,l1,l2)
                          endif

                          if (m2 .lt. 0) then
                            if (m1 .eq. 0) hsopm(i1,i2,l1,l2) =
     .                      img*a2*dsqrt(0.5d0)*tmp(i1,i2,l1,l2)
                          endif
                        endif

C               Case m=0
                        if (abs(m2) .eq. 0) then
                          if (m1 .eq. 1) hsopm(i1,i2,l1,l2) =
     .                    -a1*dsqrt(0.5d0)*tmp(i1,i2,l1,l2)
                          if (m1 .eq. -1) hsopm(i1,i2,l1,l2) =
     .                    -img*a1*dsqrt(0.5d0)*tmp(i1,i2,l1,l2)
                        endif

                      endif
C           End of block for (lso .eq. 1)
                    endif
cc          hsopm(i1,i2,l1,l2)=dcmplx(-dble(hsopm(i1,i2,l1,l2)),
cc     .                        dimag(hsopm(i1,i2,l1,l2)))
C               End loop over (m1,m2)
                  enddo
                enddo
C           End of block for lso
              endif

            enddo
          endif
        enddo
      enddo

C      print *, 'pvaglc'
C      print *, ' i1  i2     l      sig           tau           ppi0'
C      do  i1 = 1, nf1
C        do  i2 = 1, nf2
C          if (i1 .gt. nf1s .or. i2 .gt. nf2s) then
C          lmax1 = lx1(i1)
C          lmax2 = lx2(i2)
C          lmax = min0(lmax1,lmax2)
C          if (lmax .ge. 0) write(*,*) ' '
C          do  l = 0, lmax
C            write(*,500) i1,i2,l,sig(i1,i2,l),
C     .        tau(i1,i2,l)+0*ppi(i1,i2,l),ppi(i1,i2,l)
C  500       format(2i4,i6,3f14.8)
C          enddo
C        endif
C        enddo
C      enddo
C      stop
      end

      subroutine pvagm2(nf1,lmx1,lx1,f1,x1,nf2,lmx2,lx2,f2,x2,
     .nr,rofi,rwgt,vsms,lmux,sig,tau,ppi)
C- Smooth part of sig, tau, ppi (spherical part of local smooth pot)
C ----------------------------------------------------------------------
Ci Inputs
Ci   nf1   :number of 'bra' function types for each l
Ci   lmx1  :dimensions f1
Ci   lx1   :l-cutoffs for each of the nf1 functions
Ci   f1    :'bra' unaugmented envelope fn (radial part * r).
Ci         :Must be zero in channels that have no envelope functions.
Ci         :See Remarks
Ci   x1    :(radial derivative of 'bra' functions) * r = r*d(f1/r)/dr
Ci   nf2   :number of 'ket' function types for each l
Ci   lmx2  :dimensions f2
Ci   lx2   :l-cutoffs for each of the nf2 functions
Ci   f2    :'ket' unaugmented envelope fn (radial part . r).
Ci         :Must be zero in channels that have no envelope functions
Ci         :See Remarks.
Ci   x2    :(radial derivative of 'ket' functions) * r = r*d(f2/r)/dr
Ci   nr    :number of radial mesh points
Ci   rofi  :radial mesh points
Ci   rwgt  :radial mesh weights
Ci   vsms  :spherical smooth potential V2~
Ci   lmux  :maximum l for which to calculate sig,tau,ppi.
Ci         :Usually min (largest lmx1, largest lmx2)
Co Outputs
Co   sig   :overlap matrix <f1 f2>, with f1,f2 on a radial mesh
Co   tau   :kinetic energy matrix <f1 -nabla f2>
Co         := <grad f1 grad f2> + l(l+1)< f1 r^-2 f2> + surface term
Co   ppi   :potential matrix integral <f1 vsms f2>, spherical part of vsms
Cr Remarks
Cr   This routine computes the matrix elements of smoothed functions
Cr     <f1^ f2^>,  -<grad f1^ grad f2^>,   <f1^ (V2~_l=0) f2^>
Cr   which correspond to the second term in Eq. 21 for overlap,
C    and the second half of the first term Eq. 29 of the
Cr   Springer book chapter.  (pvagm2 computes only ppi matrix element
Cr   for the spherical part of V2~).  Note that there are three
Cr   flavors of sig,tau,ppi as described in the  Remarks in augmat.f:
Cr        P op P     H op P      H op H
Cr   with op = one of (1, -nabla, or V2~) and V2~ is the
Cr   one-center repsn'f of the smooth potential.
Cr   This routine makes one of these three; which one depends on
Cr   the functions f1^,f2^ passed to pvagm2.
Cr
Cr   sig_kL,k'L' for k=1..nf1 and k'=1..nf2 is diagonal in LL' and
Cr   depends only on l.  Only sig(nf1,nf2,0..l) is stored.  Ditto for
Cr   tau and ppi (for spherical part of potential) treated here.
Cr
Cr   Formula for kinetic energy.  If f1=r*phi1, x1=r*phi1'  and
Cr   f2=r*phi2,x2=r*phi2', the kinetic energy in the sphere to radius
Cr   R for channel l, excluding angular part, is
Cr     T_fg =  int_0^R (phi1) (-1/r d^2/dr^2) (r*phi2) r^2 dr
Cr          = -int_0^R (r*phi1) (d^2/dr^2) (r*phi2) dr
Cr          = -[(r phi1) d(r*phi2)/dr]_R
Cr            +int_0^R d(r*phi1)/dr * d(r*phi2)/dr
Cr          = -[r^2 phi1 dphi2/dr]_R + int_0^R r*(dphi1/dr)*r*(dphi2/dr)
Cr          = -[f1*x2]_R + int_0^R r*x1 r*x2 dr
Cr     The fourth step follows after some simple algebra.
Cu Updates
Cu   20 Jul 04 Added treatment for extended local orbitals.
Cu             Envelopes f1,f2 must be zero for all channels that
Cu             have no smooth counterparts to subtract.
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer lmux,lmx1,lmx2,nf1,nf2,nr,lx1(nf1),lx2(nf2)
      double precision rofi(nr),rwgt(nr),vsms(nr),
     .f1(nr,0:lmx1,nf1),x1(nr,0:lmx1,nf1),
     .f2(nr,0:lmx2,nf2),x2(nr,0:lmx2,nf2),
     .ppi(nf1,nf2,0:lmux),sig(nf1,nf2,0:lmux),tau(nf1,nf2,0:lmux)
C ... Local parameters
      integer i1,i2,lmax1,lmax2,lmax,l,i
      double precision sum,sim,tum,vum,xbc

      call dpzero(sig, nf1*nf2*(lmux+1))
      call dpzero(tau, nf1*nf2*(lmux+1))
      call dpzero(ppi, nf1*nf2*(lmux+1))
C     print *, '!! skip pvagm2'
C     return

      do  i1 = 1, nf1
        do  i2 = 1, nf2
          lmax1 = lx1(i1)
          lmax2 = lx2(i2)
          lmax = min0(lmax1,lmax2)
          do  l = 0, lmax
            sum = 0d0
            sim = 0d0
            tum = 0d0
            vum = 0d0
            do  i = 2, nr
              sum = sum + rwgt(i)*f1(i,l,i1)*f2(i,l,i2)
              sim = sim + rwgt(i)*f1(i,l,i1)*f2(i,l,i2)/rofi(i)**2
              tum = tum + rwgt(i)*x1(i,l,i1)*x2(i,l,i2)
              vum = vum + rwgt(i)*f1(i,l,i1)*f2(i,l,i2)*vsms(i)
            enddo
            sig(i1,i2,l) = sum
C           See Remarks for surface term in kinetic energy
            xbc = f1(nr,l,i1) * x2(nr,l,i2)
            tau(i1,i2,l) = tum + l*(l+1)*sim - xbc
            ppi(i1,i2,l) = vum
          enddo
        enddo
      enddo

C      call dvset(sig,1,nf1*nf2*(lmux+1),-99d0)

C      print *, 'pvagm2'
C      print *, ' i1  i2     l      sig           tau           ppi0'
C      do  i1 = 1, nf1
C        do  i2 = 1, nf2
C          lmax1 = lx1(i1)
C          lmax2 = lx2(i2)
C          lmax = min0(lmax1,lmax2)
C          if (lmax .ge. 0) write(*,*) ' '
C          do  l = 0, lmax
C            write(*,500) i1,i2,l,sig(i1,i2,l),tau(i1,i2,l),ppi(i1,i2,l)
C  500       format(2i4,i6,3f14.8)
C          enddo
C        enddo
C      enddo
C      stop

C     call prrmsh('f1',rofi,f1(1,0,1),nr,nr,lx1(1)+1)
C     call prrmsh('f2',rofi,f1(1,0,2),nr,nr,lx1(2)+1)

      end

      subroutine pvagm3(nf1,nf1s,lmx1,lx1,f1,v1,d1,nf2,nf2s,lmx2,lx2,f2,
     .v2,d2,nr,rofi,rwgt,lmxa,qum,lmxl,qm)
C- Moments of f1~*f2~ - f1*f2
C ----------------------------------------------------------------------
Ci Inputs
Ci   nf1   :number of 'bra' function types for each l
Ci   nf1s  :number of 'bra' function types for each l
Ci         :which are formed by linear combination of (u,s)
Ci   lmx1  :dimensions f1
Ci   lx1   :l-cutoffs for each of the nf1 functions
Ci   f1    :'bra' unaugmented envelope fn (radial part . r).
Ci         :Must be zero in channels that have no envelope functions.
Ci         :See Remarks
Ci   v1    :values of f1 at rofi(nr) (not multiplied by r), for
Ci         :functions 1..nf1s
Ci   d1    :slopes of f1 at rofi(nr) (not multiplied by r), for
Ci         :functions 1..nf1s
Ci   nf2   :number of 'ket' function types for each l
Ci   nf2s  :number of 'ket' function types for each l
Ci         :which are formed by linear combination of (u,s)
Ci   lmx2  :dimensions f2
Ci   lx2   :l-cutoffs for each of the nf2 functions
Ci   f2    :'ket' unaugmented envelope fn (radial part . r).
Ci         :Must be zero in channels that have no envelope functions
Ci         :See Remarks.
Ci   v2    :values of f2 at rofi(nr) (not multiplied by r), for
Ci         :functions 1..nf2s
Ci   d2    :slopes of f2 at rofi(nr) (not multiplied by r), for
Ci         :functions 1..nf2s
Ci   nr    :number of radial mesh points
Ci   rofi  :radial mesh points
Ci   rwgt  :radial mesh weights
Ci   lmxa  :augmentation l-cutoff
Ci   qum   :integrals (u,s) * (u,s) * r**l (see momusl)
Ci         :qum(l1,l2,l,1) = int (u_l1 u_l2) * r**l
Ci         :qum(l1,l2,l,2) = int (u_l1 s_l2) * r**l
Ci         :qum(l1,l2,l,3) = int (s_l1 s_l2) * r**l
Ci         :qum(l1,l2,l,4) = int (u_l1 g_l2) * r**l
Ci         :qum(l1,l2,l,5) = int (s_l1 g_l2) * r**l
Ci         :qum(l1,l2,l,6) = int (g_l1 g_l2) * r**l
Ci         :u and s are augmented functions defined as:
Ci         :u has val=1, slo=1 at rmax, s has val=0, slo=1
Ci   lmxl  :l-cutoff for density, potential on the radial mesh
Co Outputs
Co   qm    :integrals (f1~*f2~ - f1*f2) r**l
Cr Remarks
Cr   The qm correspond to Qkk'LL'M of Eq. 27 in
Cr      M. Methfessel, M. van Schilfgaarde, and R. A. Casali,
Cr      Lecture Notes in Physics, {\bf 535}. H. Dreysse,
Cr      ed. (Springer-Verlag, Berlin) 2000.
Cr   However, see comments in augmat.f about indices k,k'.
Cr
Cr   f1~ and f2~ are not supplied explicitly, but their contribution
Cr   to qm is supplied through qum.
Cr
Cr   If f1~ (f2~) is a valence orbital, it is the particular linear
Cr   combination of u,s that has value v1 (v2) and slope d1 (d2)
Cr   at rofi(nr); see qum(*,1..3), i.e. f1~ = V1*u1 + D1*s1
Cr
Cr   If f1~ (f2~) is a local orbital, the function is given directly;
Cr   no linear combination of functions is needed.
Cr
Cr   Define the following:
Cr   Q12uu = qum(l1,l2,lm,1)
Cr   Q12ss = qum(l1,l2,lm,3)
Cr   Q12zz = qum(l1,l2,lm,6)
Cr   Q12us = qum(l1,l2,lm,2) = Q21su;   Q12su = qum(l2,l1,lm,2)
Cr   Q12uz = qum(l1,l2,lm,4) = Q21zu;   Q12zu = qum(l2,l1,lm,4)
Cr   Q12sz = qum(l1,l2,lm,5) = Q21zs;   Q12zs = qum(l2,l1,lm,5)
Cr
Cr   Then
Cr   f1~ f2~ = V1*V2*Q12uu + V1*D2*Q12us + D1*V2*Q12su + D1*D2*Q12ss
Cr             if f1~ and f2~ are both valence orbitals;
Cr   f1~ f2~ = V1*Q12uz + D1*Q12sz
Cr             if f1~ = valence, f2~ = local;
Cr   f1~ f2~ = V2*Q12zu + D2*Q12zs = V2*Q21uz + D2*Q21sz
Cr             if f1~ = local, f2~ = valence;
Cr   f1~ f2~ = Q12zz
Cr             if f1~ = local, f2~ = local.
Cr
Cr   In principle, we must distinguish types of extended local orbitals.
Cr   Some types have no corresponding smooth part, e.g. ones with no
Cr   extension outside the MT region, while others do.  This distinction
Cr   is automatically incorporated if f1 and f2 are zero in the channels
Cr   where there no smooth part is to be subtracted.
Cr
Cr   The smooth contribution is computed by numerical integration.
Cu Updates
Cu   20 Jul 04 Added treatment for extended local orbitals.
Cu             Envelopes f1,f2 must be zero for all channels that
Cu             have no smooth counterparts to subtract.
Cu   20 Jul 04 bug fix: improper indices in qum for local orbitals
Cu   14 Sep 01 Added treatment for local orbitals.  Altered argument list.
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer nf1,nf2,nf1s,nf2s,lmx1,lmx2,lmxa,lmxl,nr
      integer lx1(nf1),lx2(nf2)
      double precision rofi(nr),rwgt(nr),
     .qm(nf1,nf2,0:lmx1,0:lmx2,0:lmxl),
     .f1(nr,0:lmx1,nf1),v1(0:lmx1,nf1s),d1(0:lmx1,nf1s),
     .f2(nr,0:lmx2,nf2),v2(0:lmx2,nf2s),d2(0:lmx2,nf2s),
     .qum(0:lmxa,0:lmxa,0:lmxl,6)
C ... Local parameters
      integer i1,i2,l1,l2,lm,i
      double precision sum,sam

      call dpzero(qm,nf1*nf2*(lmx1+1)*(lmx2+1)*(lmxl+1))

      do  i1 = 1, nf1s
        do  i2 = 1, nf2s
          do  l1 = 0, lx1(i1)
            do  l2 = 0, lx2(i2)
              do  lm = 0, lmxl

C           ... Part from smooth functions
                sum = 0d0
                do  i = 2, nr
                  sum = sum+rwgt(i)*f1(i,l1,i1)*f2(i,l2,i2)*rofi(i)**lm
                enddo

C           ... Part from augmented functions
                sam = v1(l1,i1) * v2(l2,i2) * qum(l1,l2,lm,1)
     .          + v1(l1,i1) * d2(l2,i2) * qum(l1,l2,lm,2)
     .          + d1(l1,i1) * v2(l2,i2) * qum(l2,l1,lm,2)
     .          + d1(l1,i1) * d2(l2,i2) * qum(l1,l2,lm,3)

                qm(i1,i2,l1,l2,lm) = sam-sum

              enddo
            enddo
          enddo
        enddo
      enddo

C --- Moments involving local orbitals ---
      if (nf1s .ge. nf1 .and. nf2s .ge. nf2) return
      do  i1 = 1, nf1
        do  i2 = 1, nf2
          if (i1 .gt. nf1s .or. i2 .gt. nf2s) then
            do  l1 = 0, lx1(i1)
              do  l2 = 0, lx2(i2)
                do  lm = 0, lmxl

C           ... Part from smooth functions (zero for true loc. orbitals)
                  sum = 0d0
                  do  i = 2, nr
                    sum = sum+rwgt(i)*f1(i,l1,i1)*f2(i,l2,i2)*rofi(i)**lm
                  enddo

C               Both f1~ and f2~ are local orbitals
                  if (i1 .gt. nf1s .and. i2 .gt. nf2s) then
                    sam = qum(l1,l2,lm,6)
C               f1~ is local, f2~ is linear combination of (u,s)
                  elseif (i1 .gt. nf1s) then
                    sam = v2(l2,i2) * qum(l2,l1,lm,4)
     .              + d2(l2,i2) * qum(l2,l1,lm,5)
C               f1~ is linear combination of (u,s), f2~ is local
                  elseif (i2 .gt. nf2s) then
                    sam = v1(l1,i1) * qum(l1,l2,lm,4)
     .              + d1(l1,i1) * qum(l1,l2,lm,5)
                  endif

                  if (sam .eq. 0 .and. sum .ne. 0) then
                    call rx(
     .              'gaugm: inconsistent treatment of local orbitals')
                  endif

C               if (sam .ne. 0) print 43, i1,i2,l1,l2,lm,sam,sum,sum/sam
C  43           format(5i5,3f12.6)

                  qm(i1,i2,l1,l2,lm) = sam-sum

                enddo
              enddo
            enddo
          endif
        enddo
      enddo

      end

      subroutine ppi2z(job,nsp,nspc,nf1,nf2,nlx1,nlx2,hso,ppir,ppiz)
C- Copies pi matrix to/from complex form, possibly with SO part
C ----------------------------------------------------------------------
Ci Inputs
Ci   job   : 0  do nothing
Ci         : 10's digit for copying into ppiz(2).  This step done first.
Ci         : 1  ppiz(:,:,:,:,isp,1,2) <- (3/2-isp)*hso(:,:,:,:,isp,1)
Ci         :    for isp=1..nsp, and if nspc=2:
Ci         :    ppiz(:,:,:,:,isp,2,2) <- hso(:,:,:,:,isp,2)
Ci         : 2  hso(:,:,:,:,isp,jsp) <- ppiz(:,:,:,:,isp,jsp,2)
Ci         :    for isp=1..nsp and jsp=1..nspc (reverse of 1)
Ci         :    (never checked)
Ci
Ci         : 1s digit for copying (or addition) into ppiz(1)
Ci         : 1  ppiz(:,:,:,:,isp,1,1) <- ppir(:,:,:,:,isp)
Ci         :    for isp=1..nsp
Ci         : 2  ppiz(:,:,:,:,isp,jsp,1) <-
Ci         :    ppir(:,:,:,:,isp) + ppiz(:,:,:,:,isp,jsp,2)
Ci         :    for isp=1..nsp and jsp=1..nspc
Ci         :    (Note: this copy is made after 10's digit copy)
Ci         : 3  ppir(:,:,:,:,isp) <- ppiz(:,:,:,:,isp,1,1)
Ci         :    for isp=1..nsp  (reverse of 1)
Ci         : 4  ppir(:,:,:,:,isp) <-
Ci         :    ppiz(:,:,:,:,isp,1,1) - ppiz(:,:,:,:,isp,1,2)
Ci         :    for isp=1..nsp and jsp=1..nspc (reverse of 2)
Ci         : 5  ppiz(:,:,:,:,isp,jsp,1) += ppiz(:,:,:,:,isp,jsp,2)
Ci         :    for isp=1..nsp and jsp=1..nspc
Ci         : 6  ppiz(:,:,:,:,isp,jsp,1) -= ppiz(:,:,:,:,isp,jsp,2)
Ci         :    for isp=1..nsp and jsp=1..nspc (reverse of 5)
Ci         : 7  ppiz(:,:,:,:,isp,jsp,1) += hso(:,:,:,:,isp,jsp)
Ci         :    for isp=1..nsp and jsp=1..nspc
Ci         : 8  ppiz(:,:,:,:,isp,jsp,1) -= hso(:,:,:,:,isp,jsp)
Ci         :    for isp=1..nsp and jsp=1..nspc (reverse of 7)
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   nf1   :dimensions ppi
Ci   nf2   :dimensions ppi
Ci   nlx1  :dimensions ppi
Ci   nlx2  :dimensions ppi
Cio Inputs/Outputs
Cio  hso   :spin-orbit hamiltonian, or some complex potential
Cio        :Used only for 10's digit mode nonzero
Cio  ppir  :real potential matrix
Cio  ppiz  :complex potential matrix + SO hamiltonian
Cio        :ppiz(:,:,:,:,isp,1,1) = spin-diagonal part of potential
Cio        :ppiz(:,:,:,:,1,2,1) = up-down   (12) block of potential
Cio        :ppiz(:,:,:,:,2,2,1) = down-down (21) block of potential
Cio        :ppiz(:,:,:,:,isp,jsp,2) : like ppiz(:,:,:,:,isp,jsp,1)
Cio                                 : but holds just hso
Cl Local variables
Cl         :
Cr Remarks
Cr   ppi and ppiz may occupy the same address space
Cr
Cu Updates
Cu   29 Jun 05  First created
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer job,nsp,nspc,nf1,nf2,nlx1,nlx2
      double precision ppir(nf1,nf2,nlx1,nlx2,nsp)
      double complex   ppiz(nf1,nf2,nlx1,nlx2,nsp,nspc,2)
      double complex    hso(nf1,nf2,nlx1,nlx2,nsp,nspc)
C ... Local parameters
      integer i1,i
      double precision sz
      double complex tmp(nf1,nf2,nlx1,nlx2,4)

      i1 = nf1*nf2*nlx1*nlx2

C ... Copy SO hamiltonian to/from ppiz(:,:,:,:,:,:,2)
      if (mod(job/10,10) .eq. 1) then
C       call zprm('hso',2,hso,i1,i1,nsp*nspc)
        call dcopy(2*i1*nsp*nspc,hso,1,ppiz(1,1,1,1,1,1,2),1)
C       if (mod(job/10,10) .eq. 1) then
        do  i = 1, nsp
          sz = 1.5d0 - i
          call dpcopy(hso(1,1,1,1,i,1),ppiz(1,1,1,1,i,1,2),1,i1*2,sz)
        enddo
C       endif
C       call zprm('ppiz(2)',2,ppiz(1,1,1,1,1,1,2),i1,i1,nsp*nspc)

      elseif (mod(job/10,10) .eq. 2) then
C       call zprm('ppiz(2)',2,ppiz(1,1,1,1,1,1,2),i1,i1,nsp*nspc)
        call dcopy(2*i1*nsp*nspc,ppiz(1,1,1,1,1,1,2),1,hso,1)
        do  i = 1, nsp
          sz = 1.5d0 - i
          call dpcopy(ppiz(1,1,1,1,i,1,2),hso(1,1,1,1,i,1),1,i1*2,1/sz)
        enddo
C       call zprm('hso',2,hso,i1,i1,nsp*nspc)
      endif

C ... ppiz <- ppir  or ppir + ppiz(:,:,:,:,:,:,2)
      if (mod(job,10) .eq. 1 .or. mod(job,10) .eq. 2) then
C       call prmx('ppir',ppir,i1,i1,nsp)
        call dpzero(tmp,2*i1*nsp)
C       Copy in two steps to avoid possible overwrite of ppi
        call dcopy(i1*nsp,ppir,1,tmp,2)
        call dpzero(ppiz,2*i1*nsp*nspc)
        call dcopy(2*i1*nsp,tmp,1,ppiz,1)
        if (mod(job,10) .eq. 2) then
          call daxpy(2*i1*nsp*nspc,1d0,ppiz(1,1,1,1,1,1,2),1,ppiz,1)
C         call zprm('ppiz+so',2,ppiz,i1,i1,nsp*nspc)
        endif

C        call zprm('ppiz(1)',2,ppiz,nf1*nf2,nf1*nf2,nlx1*nlx2)
C        call zprm('ppiz(2)',2,ppiz(1,1,1,1,2,1,1),
C     .    nf1*nf2,nf1*nf2,nlx1*nlx2)

C ... ppir <- ppiz or ppir - ppiz(:,:,:,:,:,:,2)
      elseif (mod(job,10) .eq. 3 .or. mod(job,10) .eq. 4) then

C       call zprm('ppiz',2,ppiz,i1,i1,nsp*nspc)
C       Copy in two steps to avoid possible overwrite of ppi
        call dcopy(2*i1*nsp,ppiz,1,tmp,1)
        call dcopy(i1*nsp,tmp,2,ppir,1)

C       subtract Re(hso), spin-diagonal blocks
        if (mod(job,10) .eq. 4) then
          call daxpy(i1*nsp,-1d0,ppiz(1,1,1,1,1,1,2),2,ppir,1)
        endif
C       call prmx('ppir',ppir,i1,i1,nsp)

      elseif (mod(job,10) .eq. 5 .or. mod(job,10) .eq. 6) then
C       call zprm('starting ppiz',2,ppiz,i1,i1,nsp*nspc)
        sz = 1
        if (mod(job,10) .eq. 6) sz = -1
        call daxpy(2*i1*nsp*nspc,sz,ppiz(1,1,1,1,1,1,2),1,ppiz,1)
C       call zprm('ending ppiz',2,ppiz,i1,i1,nsp*nspc)

      elseif (mod(job,10) .eq. 7 .or. mod(job,10) .eq. 8) then
C       call zprm('starting ppiz',2,ppiz,i1,i1,nsp*nspc)
        sz = 1
        if (mod(job,10) .eq. 6) sz = -1
        call daxpy(2*i1*nsp*nspc,sz,hso,1,ppiz,1)

C        call zprm('ppiz(1)',2,ppiz,nf1*nf2,nf1*nf2,nlx1*nlx2)
C        call zprm('ppiz(2)',2,ppiz(1,1,1,1,2,1,1),
C     .    nf1*nf2,nf1*nf2,nlx1*nlx2)

      endif

      end

C      subroutine pvagm4(nf1,lx1,nf2,lx2,tau)
CC- Symmetrize kinetic energy matrix
CC ----------------------------------------------------------------------
CCi Inputs
CCi   nf1   :number of 'bra' function types for each l
CCi   lx1   :l-cutoffs for each of the nf1 functions
CCi   f1    :'bra' unaugmented envelope fn (radial part . r).
CCi         :Must be zero in channels that have no envelope functions.
CCi         :See Remarks
CCi   nf2   :number of 'ket' function types for each l
CCi   lx2   :l-cutoffs for each of the nf2 functions
CCi   f2    :'ket' unaugmented envelope fn (radial part . r).
CCi         :Must be zero in channels that have no envelope functions
CCi         :See Remarks.
CCio Inputs/Outputs
CCio  tau   :Kinetic energy matrix is symmetrized.  See Remarks.
CCr Remarks
CCr   Symmetrization of kinetic energy.  The kinetic energy matrix
CCr   should always be symmetric.  In general the difference
CCr   T_fg - T_gf  between orbitals f and g is
CCr     T_fg - T_gf = -[r^2 (fg' - f'g)]_rmax = -W{f,g;rmax}
CCr   which is nonzero.  However, the kinetic energy matrix always
CCr   takes of two forms:
CCr      1) it consists of a difference T_fg-T^_fg = <f T g> - <f^ T g^>
CCr         between true and smoothed functions where f matches f^
CCr         and g matches g^ at rmax.  In this case,
CCr            T_fg - T^_fg - (T_gf - T^_fg) = 0
CCr         because W{f,g;rmax} = W{f^,g^;rmax}
CCr      2) either f or g is a local orbital with value,slope=0.
CCr         Then -W{f,g;rmax} = 0.
CCr
CCr   If the matrix is properly constructed it will be symmetric and
CCr   this routine is unnecessary.  However, when f or g is an extended
CCr   local orbital, only the average (T_fg+T_gf)/2 is kept
CCr   (see potpus.f, hab(5) and hab(6)), even while the smooth
CCr   contribution is calculated properly. Thus in such a case
CCr     tau(f,g) = (T_fg+T_gf)/2 - T^_fg = T_fg - dT - T^_fg
CCr     tau(g,f) = (T_fg+T_gf)/2 - T^_gf = T_gf + dT - T^_gf
CCr   where dT = (T_fg-T_gf)/2
CCr   Therefore the correct (and symmetric) matrix elements are
CCr   simply the average of input values (tau(f,g)+tau(g,f))/2 since
CCr     T_fg-T^_fg = T_gf-T^_gf
CCr                = tau(f,g)+dT = tau(g,f)-dT = (tau(f,g)+tau(g,f))/2
CCr
CCr  NB: this routine can only be applied wheh the functions f,g
CCr  are taken from the same group.
CC ----------------------------------------------------------------------
C      implicit none
CC ... Passed parameters
C      integer nf1,nf2,lx1(nf1),lx2(nf2)
C      double precision tau(nf1,nf2,0:*)
CC ... Local parameters
C      integer i1,i2,lmax1,lmax2,lmax,l
C      double precision tbar
C
C      do  i1 = 1, nf1
C        do  i2 = 1, nf2
C          lmax1 = lx1(i1)
C          lmax2 = lx2(i2)
C          lmax = min0(lmax1,lmax2)
C          do  l = 0, lmax
C            if (abs(tau(i1,i2,l)-tau(i2,i1,l)).gt. 1d-6) then
C              print *, 'hi!'
C            write(*,500) i1,i2,l,tau(i1,i2,l),tau(i1,i2,l)-tau(i2,i1,l)
C
C            endif
C            tbar = (tau(i1,i2,l) + tau(i2,i1,l))/2
C            tau(i1,i2,l) = tbar
C            tau(i2,i1,l) = tbar
C          enddo
C        enddo
C      enddo
C
C      print *, 'pvagm4'
C      do  i1 = 1, nf1
C        do  i2 = 1, nf2
C          lmax1 = lx1(i1)
C          lmax2 = lx2(i2)
C          lmax = min0(lmax1,lmax2)
C          write(*,*) ' '
C          do  l = 0, lmax
C            write(*,500) i1,i2,l,tau(i1,i2,l),tau(i1,i2,l)-tau(i2,i1,l)
C  500       format(2i4,i6,3f14.8)
C          enddo
C        enddo
C      enddo
Cc      pause
C
C      end

