      subroutine makwts(metal,lwtkb,lrhos,ldos,lkkrm,lorth,nl,nlo,nrhos,
     .  isp,nsp,nspc,nbas,ikp,nkp,nev,wtkp,zval,pph,ldim,idim,
     .  lihdim,indxsh,ldnpp,dnpp,nqpp,qpp,zll,zil,c,ipc,nclass,nrclas,
     .  wk1,wk2,wk3,eband,wtkb,qnu,doswt,dosw2,rhos)
C- Calculate weights l- or L- decomposed into spheres by class
C ----------------------------------------------------------------------
Ci Inputs
Ci   metal :T => metal, F => nonmetal
Ci   lwtkb :energy weights wtkb are known.
Ci         :If lwtkb is set, makwts does the following:
Ci         :  * accumulates moments (and rhos) directly using wtkb.
Ci         :  * makes only the 0th energy moment of doswt; see Remarks
Ci         :Otherwise, makwts makes energy moments 0,1,2 of doswt
Ci   lrhos :T, accumulate dosw2 for site diagonal spin-density matrix
Ci   ldos  :NO LONGER USED
Ci         :a combination of the following integers
Ci         : 1 make dos
Ci         : 2 generate weights for partial dos
Ci         : 4 generate weights for m-decompos'n of pdos
Ci   lkkrm :T accumulate moments 'KKR' style; see Remarks
Ci         :F accumulate moments 'basis' style; see Remarks
Ci   lorth :assume c=0,zil,pgam=0 (for two-center hamiltonian)
Ci   nl    :(global maximum l) + 1
Ci   nlo   :either nl or nl*nl, depending on whether integrated
Ci         :quantities are m-resolved
Ci   isp   :current spin channel (1 or 2)
Ci   nsp   :2 for spin-polarized case, otherwise 1
Ci   nspc  :2 if spin-up and spin-down channels are coupled; else 1.
Ci   nbas  :size of basis
Ci   ikp   :k-point label
Ci   nkp   :number of irreducible k-points (bzmesh.f)
Ci   nev   :actual number of eigenvectors generated
Ci   wtkp  :weight of k-point, including spin degeneracy (bzmesh.f)
Ci         :NB: only abs(wtkp) is used; see bzmesh.f
Ci   zval  :valence charge
Ci   pph   :potential parameters in downfolding order (makpph.f)
Ci         : pph(1) = enu
Ci         : pph(2) = calpha
Ci         : pph(3) = sqrdel
Ci         : pph(4) = palpha
Ci         : pph(5) = oalp
Ci   ldim  :dimension of hamiltonian matrix (makidx.f)
Ci   idim  :dimension of block of orbitals to be folded down (makidx.f)
Ci   lidim :dimensions pph
Ci   indxsh:permutations ordering orbitals into l+i+h blocks (makidx.f)
Ci   ldnpp :T, return phi,phid decomposition of the norm
Ci   nqpp  :nonzero, accumulate qpp
Ci   zll   :ll block of eigenvectors in alpha representation
Ci   zil   :il block of eigenvectors in alpha representation
Ci   c     :(C-enu + srdel*S*srdel) * zll = h * zll = amount of phidot^a
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci         :ipc(1)=0 -> makwts resolves by site, not by class
Ci         :In that case, nclass should be nbas
Ci   nclass:number of inequivalent classes
Ci   nrclas:nrclas(i) = number of atoms in class i
Ci         :ipc(1)=0 -> nrclas is not used
Ci   wk1   :work array of dimension 2(ldim+idim)*nspc
Ci   wk2   :work array of dimension 2(ldim+idim)*nspc
Ci   wk3   :work array of dimension 2(ldim+idim)*nspc
Ci   eband :energy bands; alias eb (sec*.f)
Ci   wtkb  :(ltkb=t) sampling weights associated with band index for
Ci         :this qp; see bzwts.f, which enable direct accumulation
Ci         :of qnu,rhos
Ci         :(lwtkb=f) not used.
Co Outputs  -- must all be initialised to zero on first entry --
Co   doswt :weights for decomposition of partial dos into sphere parts.
Co         :(lwtkb=f and metal=t) doswt(1..3) are made storing energy
Co         :           moments 0..2, and if nspc=2, doswt(4) is made
Co         :           storing the orbital moment.
Co         :(lwtkb=t or metal=f) doswt(1) is made storing energy moment
Co         :           0 only; and if nspc=2, doswt(2) is made storing
Co         :           the orbital moment.
Co         :doswt is used for one of two independent purposes:
Co         :1. Weights for the ASA energy moments.  In the standard ASA
Co         :   self-consistency cycle, specification of the 0th, 1st,
Co         :   2nd energy moments of partial DOS integrated over the
Co         :   sphere, qnu, is sufficient to completely specify a sphere
Co         :   density; see routine atomsc.  If the weights for each
Co         :   eigenvector are known in advance (lwtkb=T), makwts
Co         :   accumulates the contribution from each eigenvector at
Co         :   this qp and spin to qnu directly, and doswt isn't used
Co         :   for making qnu.  If lwtkb=F, the partial dos weights
Co         :   are be saved for all three energy moments, for subseqent
Co         :   calculation of qnu (see moment.f).  Note that for the ASA
Co         :   the qnu need only be l-resolved, and the m index in
Co         :   doswt can be contracted for this purpose.
Co         :
Co         :2. Decomposition of partial DOS inside MT sphere.  To
Co         :   generate this or related quantities, the zeroth moment of
Co         :   doswt at each q.p. together with the eigenvalues; see
Co         :   routines asados and contet for generation of partial dos.
Co         :   The weights may be resolved into Rl or Rlm.  If weights
Co         :   are contracted in m, nlo=nl; if kept m-resolved,
Co         :   nlo=nl*nl.  In the latter case makwts doesn't know about
Co         :   symmetry operations and it is the caller's responsibility
Co         :   to appropriately allocate arrays qnu,doswt,dosw2,rhos
Co         :   by site, and the indexing arrays ipc and nrclas should
Co         :   contain one class per site.
Co         :
Co         :   While it is possible to both resolve doswt by m and keep
Co         :   moments 0..2, it is expensive in memory requirements.
Co         :   It is advised that if context 2 is used with m resolution
Co         :   that makwts is called with lwtkb=T, so that doswt is
Co         :   used only for the second context.
Co         :
Co   qnu   :energy moments of charges in RL channels (see description of
Co         :doswt above).  Accumulated if metal=F or if lwtkb=T.
Co         :Otherwise, qnu may be generated from doswt (moment.f)
Co   dosw2 :(lrhos=T) weights of off-diagonal (in spin space) of doswt
Co   dnpp  :(ldnpp=T) phi,phidot decomposition of the norm for each evec
Co         :(used by optics branch)
Co   qpp   :(nqpp>0) coefficients to all wave function products, needed for
Co         :density matrix inside spheres or nonspherical density.
Co         :qpp(*,1,ib) phi-phi   products for (L<=L')
Co         :qpp(*,2,ib) phi-dphi  products for all (L,L')
Co         :qpp(*,4,ib) dphi-dphi products for (L<=L')
Cr Remarks
Cr   Moments are over the occupied partial densities of states; i.e.,
Cr   they are l-decomposed moments of the charge density in each sphere.
Cr   The decomposition of the norm is given in KANPUR NOTES (eq. 2.66);
Cr   i.e. inside the sphere, the wave function is
Cr      psi = zll (phi_l + h phd_l^alpha) + zil (phi_i)
Cr          = zll (phi_l + h (phd_l + o phi_l)) + zil (phi_i)
Cr
Cr   Thus, decomposition of wave function products is as follows:
Cr     qpp(L,L') = (1+oh)z+_RL (1+oh)z_RL'  phi-phi
Cr                 (1+oh)z+_RL (hz)_RL'     phi-dphi
Cr                    (hz)+_RL (hz)_RL'     dphi-dphi
Cr
Cr   ASA moments: Charge from the high waves and combined correction is
Cr   renormalized evenly into the lower channels; see variable rnorm.
Cr   ASA moments are accumulated in one of two ways.
Cr   *Fixed-basis style (corresponds to Stuttgart ASA code)
Cr     Moments are coefficients to (phi**2, phi*phidot, phidot**2)
Cr     as defined in KANPUR eq. 2.67, but renormalized for ccor.
Cr     Define :  m0 = z+z,   m1 = z+ c,   m2 = c+ c
Cr
Cr     In representation alpha (^a), the charge density is
Cr      (1)  n(r) = m0^a\phi^2 + 2 m1^a\phi\phidot^a + m2^a\phidot^a^2
Cr
Cr     Moments qnu(0..2) are related to m(0..2) except that:
Cr      (a) they correspond to the orthogonal representation.
Cr      (b) qnu(0) is taken by convention to be the electronic charge
Cr          q0 = m0 + p^gamma m2,  rather than m0.
Cr
Cr     The relation of \phidot^a to \phidot (in gamma) is
Cr      (2) \phidot^a = \phidot + o^a\phi
Cr
Cr     which leads to
Cr      (3) n(r) = (m0^a + 2 m1^a o^a + m2^a o^a^2) \phi^2 +
Cr                 (2 m1^a + 2 m2^a o^a) \phi\phidot + m2^a\phidot^2
Cr
Cr     and allows the identification of the moments in gamma as
Cr      (4) m0 = m0^a  +  2 m1^a o^a  +  m2^a o^a^2
Cr          m1 = m1^a  +  m2^a o^a
Cr          m2 = m2^a
Cr
Cr     makwts accumulates q(0..2) directly.  For the lower waves:
Cr          q0 = sum ((1+oh)z)+ (1+oh)z  +  (hz)+ p^gamma (hz)
Cr          q1 = sum ((1+oh)z)+ (hz)
Cr          q2 = sum (hz)+ (hz)
Cr
Cr   *KKR style (corresponds to Methfessel's ASA)
Cr     m(0..2) are true power moments accumulated over a renormalised
Cr     DOS.  This mode is consistent with Methfessel's"compressed atom"
Cr     philosophy.  Here, the energy-dependent partial wave is
Cr      (5) phi(e) - enu = phi(enu) + (e-enu) phidot(enu)
Cr
Cr    m(0..2) are now the 0th, 1st and second moments of the DOS, and
Cr    are later identified with the coefficients to the density, Eq. 1.
Cr    This approach is not exactly equivalent to the first approach
Cr    when the combined correction is included
Cr          m0 = sum ((1+oh)z)+ (1+oh)z
Cr          m1 = sum (e-enu) m0
Cr          m2 = sum (e-enu)^2 m0
Cr    In practice the code computes
Cr          q0 = sum ((1+oh)z)+ (1+oh)z  +  (hz)+ p^gamma (hz)
Cr          q1 = sum q0 (e-enu)
Cr          q2 = sum q0 (e-enu)^2
Cr    which erroneously includes the p^gamma term in q1,q2
Cr
Cb Bugs
Cb   KKR style accumulation of q1 and q2 is slightly erroneous.
Cu Updates
Cu   08 Nov 07 (J. Xu) qpp is complex; diag qpp is made in noncoll case
Cu   21 Apr 04 rhos,dosw2 can be m-dependent.  New argument list
Cu   04 Jan 04 bug fix for noncollinear nonmetallic case
Cu   04 May 03 bug fix: rhos accumulates correctly when lwtbk=T
Cu   18 Jan 02 some functionality removed, for better modularity
Cu             Altered argument list.
Cu   17 Dec 01 Average downfolded moments channels, noncollinear case
Cu   16 May 01 Bug fix for magnetic insulator case
Cu   17 Apr 01 Added switch lkkrm
Cu   03 Feb 01 logical dos changed into ldos
Cu   19 Jan 01 works for downfolded orbitals, noncollinear case
Cu   20 Dec 00 dnpp may be accumulated for the noncollinear case
Cu   28 Apr 98 makwts can make dnpp and/or qpp
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      logical metal,lwtkb,lrhos,ldnpp,lkkrm,lorth
      integer nl,nlo,nsp,nspc,nbas,nrhos,isp,ikp,nkp,ldim,idim,lihdim,
     .  nev,nqpp,ipc(nbas),nclass,nrclas(nclass),indxsh(ldim),ldos
      double precision qnu(0:2,0:nl-1,nsp,nclass),zval,
     .  pph(5,lihdim,nspc),eband(nev),
     .  wk1(2,ldim+idim,nspc),zll(ldim,nspc,ldim*nspc,2),wtkb(ldim),
     .  wk2(2,ldim+idim,nspc),zil(idim,nspc,ldim*nspc,2),wtkp(nkp),
     .  wk3(2,ldim+idim,nspc),c(ldim,nspc,ldim*nspc,2),
     .  doswt(0:nlo-1,nclass,nev,nspc,2+nspc),
     .  dnpp(nl**2,nbas,nspc,2,2,nev),
     .  dosw2(2,nrhos,nclass,nev,0:2),rhos(2,0:2,nrhos,2,2,nclass)
      double complex qpp(nqpp,4,nsp,nbas)
C Local variables
      logical lall,noclas
      integer i,ib,ic,id,idq,ilm,ilm1,ilm2,iprint,iqpd,iqpp,j,j1,
     .  j2,jsp,k,k2,ksp,l,l1,l2,lm,m,m1,m2,nstate,lidim,nrc
      double precision rnorm,wt,e,w0,w1,w2,xx,xx2,bync,pgam

C --- Setup ---
      lall = nlo .eq. nl**2
      lidim = ldim+idim
      if (.not. lall .and. nlo .ne. nl) call rx('makwts: bad input nl')
      nstate = (int(zval) + 1)/2
      if (.not. metal) then
        if (nstate*2 .ne. int(zval) .and. iprint() .ge. 10 ) print *,
     .  'MAKWTS (warning) : odd no. of electrons in non-metal'
C       if (lrhos) call rx('makwts: rhos not allowed with insulator')
      endif
      if (ldnpp) call dpzero(dnpp,nl**2*nspc*nbas*4*nev)
      noclas = ipc(1) .eq. 0
      if (noclas .and. nclass .ne. nbas)
     .  call rx('makwts: nclass ne nbas')

C --- For each eigenvector, do ---
      do  5  ib = 1, nev

C   ... Sampling weight for this evec, used in cases when it is known
        if (lwtkb) then
          wt = wtkb(ib)
        else
          wt = 0
          if (.not. metal .and. ib .le. nstate*nspc)
     .      wt = abs(wtkp(ikp))/nsp
        endif

C   --- wk1 = (1+oh) z;  wk2 = sqrt(p^gam) h z;  wk3 = h z  ---
        rnorm = 0
        do  jsp = 1, nspc
C         Three-center hamiltonian
          if (.not. lorth) then
          do  j = 1, ldim
            wk1(1,j,jsp) = zll(j,jsp,ib,1) + pph(5,j,jsp)*c(j,jsp,ib,1)
            wk1(2,j,jsp) = zll(j,jsp,ib,2) + pph(5,j,jsp)*c(j,jsp,ib,2)
C           Add a tiny number in case p^gamma is zero
            pgam = dsqrt(max(pph(4,j,jsp)-pph(5,j,jsp)**2,1d-15))
            wk2(1,j,jsp) = pgam*c(j,jsp,ib,1)
            wk2(2,j,jsp) = pgam*c(j,jsp,ib,2)
            wk3(1,j,jsp) = c(j,jsp,ib,1)
            wk3(2,j,jsp) = c(j,jsp,ib,2)
          enddo
          do  j = ldim+1, lidim
            wk1(1,j,jsp) = zil(j-ldim,jsp,ib,1)
            wk1(2,j,jsp) = zil(j-ldim,jsp,ib,2)
            wk2(1,j,jsp) = 0
            wk2(2,j,jsp) = 0
            wk3(1,j,jsp) = 0
            wk3(2,j,jsp) = 0
          enddo
C         Two-center hamiltonian: p^gam==0, h z = 0
          else
          call rxx(.not.lkkrm,'makwts lorth must be used with kkr mode')
          call dpzero(wk1(1,1,jsp),2*lidim)
          call dpzero(wk2(1,1,jsp),2*lidim)
          call dpzero(wk3(1,1,jsp),2*lidim)
          do  j = 1, ldim
            wk1(1,j,jsp) = zll(j,jsp,ib,1)
            wk1(2,j,jsp) = zll(j,jsp,ib,2)
          enddo
          endif
C         call dpzero(wk1,2*(lidim))
C         call zprm('cofp=(1+oh)z',2,wk1(1,1,jsp),lidim,lidim,1)
C         call zprm('cofd=hz',2,wk2,lidim,lidim,1)

C     ... Normalization not exactly 1 because of combined correction
          call yydotc(lidim,wk1(1,1,jsp),wk1(2,1,jsp),2,
     .                wk1(1,1,jsp),wk1(2,1,jsp),2,xx,xx2)
          rnorm = rnorm + xx
          call yydotc(lidim,wk2(1,1,jsp),wk2(2,1,jsp),2,
     .                wk2(1,1,jsp),wk2(2,1,jsp),2,xx,xx2)
          rnorm = rnorm + xx
        enddo

C   ... ASA renormalization owing to combined correction
        rnorm = 1d0 / rnorm

C   --- ASA moments, or weights for moments and DOS for this band ---
        id = 0
        do  i = 1, nbas
          if (noclas) then
            ic = i
            nrc = 1
          else
C           Workaround for intel compiler bug
C           ic = ipc(i)
C           nrc = nrclas(ic)
            call pvmakw(i,ipc,nrclas,ic,nrc)
          endif
          ilm = 0
          do  l = 0, nl-1
          do  m = -l, l
C           ilm = compound l+m index; lm = ilm or l, depending on lall
            ilm = ilm + 1
            lm = ilm-1
            if (.not. lall) lm = l
            k2 = l+1
            if (nrhos .eq. nl*nl) k2 = ilm
            id = id + 1
            j = indxsh(id)

            if (j .le. lidim) then

C       ... Site-diagonal, spin density matrix rho_12
            if (lrhos) then
              e = eband(ib) - pph(1,j,1)/2 - pph(1,j,2)/2
C             Re (rho_12) = Re (z+_1 z_2) for o and p^gam = 0
              w0 = wk1(1,j,1)*wk1(1,j,2) + wk1(2,j,1)*wk1(2,j,2) +
     .             wk2(1,j,1)*wk2(1,j,2) + wk2(2,j,1)*wk2(2,j,2)
              w0 = w0 * rnorm / nrc
              dosw2(1,k2,ic,ib,0) = dosw2(1,k2,ic,ib,0) + w0
              dosw2(1,k2,ic,ib,1) = dosw2(1,k2,ic,ib,1) + w0*e
              dosw2(1,k2,ic,ib,2) = dosw2(1,k2,ic,ib,2) + w0*e*e
C             Im (rho_12) = Im (z+_1 z_2) for o and p^gam = 0
C             The sign of rhos(12) depends on the definition of the
C             spin-density matrix, e.g. whether rhos=z z+  or  rhos=z+ z
C             Heuristically this convention is consistent with the
C             standard definitions.
              w0 = wk1(2,j,1)*wk1(1,j,2) - wk1(1,j,1)*wk1(2,j,2) +
     .             wk2(2,j,1)*wk2(1,j,2) - wk2(1,j,1)*wk2(2,j,2)
C             Alternatively, change the sign of w0
C             w0 = wk1(1,j,1)*wk1(2,j,2) - wk1(2,j,1)*wk1(1,j,2) +
C    .             wk2(1,j,1)*wk2(2,j,2) - wk2(2,j,1)*wk2(1,j,2)
              w0 = w0 * rnorm / nrc
              dosw2(2,k2,ic,ib,0) = dosw2(2,k2,ic,ib,0) + w0
              dosw2(2,k2,ic,ib,1) = dosw2(2,k2,ic,ib,1) + w0*e
              dosw2(2,k2,ic,ib,2) = dosw2(2,k2,ic,ib,2) + w0*e*e
            endif

            bync = rnorm / nrc
C           jsp = spin index for evecs; ksp=spin index for qnu
            do  jsp = 1, nspc
              ksp = max(jsp,isp)

C             w0 = ((1+oh) z)+ ((1+oh) z)  +  p^gam (h z)+ (h z)
              w0 = wk1(1,j,jsp)**2 + wk1(2,j,jsp)**2 +
     .             wk2(1,j,jsp)**2 + wk2(2,j,jsp)**2
C             w1 = ((1+oh) z)+ (h z)
              w1 = wk1(1,j,jsp)*wk3(1,j,jsp) + wk1(2,j,jsp)*wk3(2,j,jsp)
C             w2 = ((1+oh) z)+ ((1+oh) z)
              w2 = wk3(1,j,jsp)**2 + wk3(2,j,jsp)**2

C             Spin-average up-, down- downfolded moments, noncol case
              if (nspc .eq. 2 .and. j .gt. ldim) then
                w0 = wk1(1,j,3-jsp)**2 + wk1(2,j,3-jsp)**2 +
     .               wk2(1,j,3-jsp)**2 + wk2(2,j,3-jsp)**2 + w0
                w0 = w0/2
              endif

              w0 = w0 * bync
              w1 = w1 * bync
              w2 = w2 * bync

              if (lwtkb .or. .not. metal) then

C             Band weight known a priori : add to doswt,qnu KKR style
              if (lkkrm .or. j .gt. ldim) then
                doswt(lm,ic,ib,jsp,1) = doswt(lm,ic,ib,jsp,1) + w0
                if (nspc .eq. 2)
     .          doswt(lm,ic,ib,jsp,2) = doswt(lm,ic,ib,jsp,2) + w0*m
                e = eband(ib) - pph(1,j,jsp)

C                if (l.eq.3 .and. j.ge.10 .and. j.le.16 .and.
C     .            w0 .gt. 0.1d0 .and. abs(wt) .gt. 0.001) then
C                  print 543, ikp, isp, ib, j,
C     .              wt, e, w0, eband(ib), pph(1,j,jsp)
C  543             format(4i3,6f12.6)
C                endif

                qnu(0,l,ksp,ic) = qnu(0,l,ksp,ic) + wt*w0
                qnu(1,l,ksp,ic) = qnu(1,l,ksp,ic) + wt*w0*e
                qnu(2,l,ksp,ic) = qnu(2,l,ksp,ic) + wt*w0*e*e

C             Band weight known a priori : add to doswt,qnu LMTO style
              else
                doswt(lm,ic,ib,jsp,1) = doswt(lm,ic,ib,jsp,1) + w0
                if (nspc .eq. 2)
     .          doswt(lm,ic,ib,jsp,2) = doswt(lm,ic,ib,jsp,2) + w0*m
                qnu(0,l,ksp,ic) = qnu(0,l,ksp,ic) + wt*w0
                qnu(1,l,ksp,ic) = qnu(1,l,ksp,ic) + wt*w1
                qnu(2,l,ksp,ic) = qnu(2,l,ksp,ic) + wt*w2
              endif

              if (lrhos) then
                do  k = 0, 2
                rhos(1,k,k2,ksp,ksp,ic)=rhos(1,k,k2,ksp,ksp,ic)+wt*w0
                if (l .eq. m .or. lall) then
                rhos(1,k,k2,ksp,3-ksp,ic) = rhos(1,k,k2,ksp,3-ksp,ic)
     .                                    + dosw2(1,k2,ic,ib,k)*wt
                rhos(2,k,k2,ksp,3-ksp,ic) = rhos(2,k,k2,ksp,3-ksp,ic)
     .                               + dosw2(2,k2,ic,ib,k)*wt*(3-2*ksp)
                endif
                enddo
              endif

C             Band weight not known a priori : add to doswt KKR style
C             Metal : add to doswt KKR style
              elseif (metal .and. (lkkrm .or. j .gt. ldim)) then
                e = eband(ib) - pph(1,j,jsp)
                doswt(lm,ic,ib,jsp,1) = doswt(lm,ic,ib,jsp,1) + w0
                doswt(lm,ic,ib,jsp,2) = doswt(lm,ic,ib,jsp,2) + w0*e
                doswt(lm,ic,ib,jsp,3) = doswt(lm,ic,ib,jsp,3) + w0*e*e
                if (nspc .eq. 2)
     .          doswt(lm,ic,ib,jsp,4) = doswt(lm,ic,ib,jsp,4) + w0*m

C             Metal : add to doswt LMTO style
              elseif (metal) then
                doswt(lm,ic,ib,jsp,1) = doswt(lm,ic,ib,jsp,1) + w0
                doswt(lm,ic,ib,jsp,2) = doswt(lm,ic,ib,jsp,2) + w1
                doswt(lm,ic,ib,jsp,3) = doswt(lm,ic,ib,jsp,3) + w2
                if (nspc .eq. 2)
     .          doswt(lm,ic,ib,jsp,4) = doswt(lm,ic,ib,jsp,4) + w0*m
              endif
            enddo
            endif
          enddo
          enddo
        enddo

C   --- Accumulate nonspherical moments using passed wtkb ---
        if (nqpp .gt. 0 .and. (metal .or. ib .le. nstate)) then
C       if (nqpp .gt. 0 .and. (metal .or. ib .le. 1)) then

          if (.not. lwtkb .and. metal) call
     .      rx('MAKWTS: missing weights for nonspherical rho ')
C         if (nspc .eq. 2) call rx('qnsp not noncollinear')

C         jsp = spin index for evecs; ksp=spin index for qpp
          do  jsp = 1, nspc
          ksp = max(jsp,isp)
          id = 0
          idq = 0
          do  i = 1, nbas

          iqpp = 0
          iqpd = 0
          ilm1 = 0
          do  l1 = 0, nl-1
          do  m1 = -l1, l1
          ilm1 = ilm1+1
          j1 = indxsh(id+ilm1)
          if (j1 .gt. lidim) goto 36

          ilm2 = 0
          do   l2 = 0, nl-1
          do   m2 = -l2, l2
          ilm2 = ilm2+1
          j2 = indxsh(id+ilm2)
          if (j2 .gt. lidim) goto 37

          if (ilm2 .le. ilm1) then
            iqpp = iqpp+1
c           if (i .eq. 1 .and. iqpp .eq. 2) print *, qpp(iqpp,1,i)
            qpp(iqpp,1,ksp,i) = qpp(iqpp,1,ksp,i) + wt*rnorm*dcmplx(
     .        wk1(1,j1,jsp)*wk1(1,j2,jsp)+wk1(2,j1,jsp)*wk1(2,j2,jsp),
     .       -wk1(1,j1,jsp)*wk1(2,j2,jsp)+wk1(2,j1,jsp)*wk1(1,j2,jsp))
            qpp(iqpp,4,ksp,i) = qpp(iqpp,4,ksp,i) + wt*rnorm*dcmplx(
     .        wk3(1,j1,jsp)*wk3(1,j2,jsp)+wk3(2,j1,jsp)*wk3(2,j2,jsp),
     .       -wk3(1,j1,jsp)*wk3(2,j2,jsp)+wk3(2,j1,jsp)*wk3(1,j2,jsp))
C           print*, qpp(iqpp,1,ksp,i)
          endif

          iqpd = iqpd+1
          qpp(iqpd,2,ksp,i) = qpp(iqpd,2,ksp,i) + wt*rnorm*dcmplx(
     .      wk1(1,j1,jsp)*wk3(1,j2,jsp)+wk1(2,j1,jsp)*wk3(2,j2,jsp),
     .     -wk1(1,j1,jsp)*wk3(2,j2,jsp)+wk1(2,j1,jsp)*wk3(1,j2,jsp))

   37     continue
          enddo
          enddo
   36     continue
          enddo
          enddo
          id = id + nl**2
          idq = idq + nl**4
          enddo
          enddo
        endif

C   --- Projection of phi, phidot for this eigenvector ---
        if (ldnpp) then
          id = 0
          do  40  i = 1, nbas
          lm = 0
          do  41  l = 0, nl-1
          do  42  m = -l, l
            lm = lm + 1
            id = id + 1
            j = indxsh(id)
            if (j .gt. lidim) goto 42
            do  43  jsp = 1, nspc
              dnpp(lm,i,jsp,1,1,ib) = wk1(1,j,jsp)
              dnpp(lm,i,jsp,1,2,ib) = wk1(2,j,jsp)
              dnpp(lm,i,jsp,2,1,ib) = wk3(1,j,jsp)
              dnpp(lm,i,jsp,2,2,ib) = wk3(2,j,jsp)
   43       continue
   42     continue
   41     continue
   40     continue
        endif

C ... End of loop over evals
    5 continue

      end
      subroutine pvmakw(i,ipc,nrclas,ic,nrc)
C- Kernel called by makwts as workaround to intel compiler
C     implicit none
      integer i,ic,nrc,ipc(i),nrclas(ic)
      ic = ipc(i)
      nrc = nrclas(ic)
      end

