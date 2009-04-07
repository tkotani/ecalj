      subroutine hmltnc(mode,nbas,nl,indxsh,qspirl,eula,neul,pph,sod,
     .  lasa,ccor,lss,lnc,lso,ccd,vmtz,ldim,lihdim,sk,hk,ok,wk)
C- Three-center ASA noncollinear hamiltonian and overlap
C ---------------------------------------------------------------------
Ci Inputs
Ci   mode  :0 hamiltonian-> hk and overlap -> ok (sk overwritten)
Ci         :1 small h -> hk
Ci         :2 1+oh -> hk
Ci         :3 small h -> hk, input sk is noncol srdel * S * srdel
Ci         :4 1+oh -> hk, input sk is noncol srdel * S * srdel
Ci         :  modes 1-4: 
Ci   nbas  :size of basis
Ci   nl    :(global maximum l) + 1
Ci   indxsh:permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   qspirl:Parameters defining spin-spiral
Ci   eula  :Euler angles for noncollinear spins
Ci   neul  :1 if Euler angles are l-independent, nl otherwise
Ci   pph   :potential parameters in downfolding order (makpph.f)
Ci         :pph(1..5,i,is): parms for ith RL and is(th) spin channel.
Ci         :pph(1) : enu
Ci         :pph(2) : calpha
Ci         :pph(3) : sqrdel
Ci         :pph(4) : palpha
Ci         :pph(5) : oalp
Ci   sod   :diagonal parms for S/O ham., LL block.
Ci         :sod(*,1..3,isp): 1-, 2-, 3- center terms, ++,-- blocks
Ci         :sod(*,1,   3):   1- center terms, L- S+ block
Ci         :sod(*,2..3,3):   2- center terms, L- S+ block
Ci         :sod(*,4,   3):   3- center terms, L- S+ block
Ci   lasa  :switch telling whether to add ASA H
Ci   ccor  :switch telling whether to add combined correction
Ci   lss   :switch for spin spiral
Ci   lnc   :switch turning on noncollinear hamiltonian
Ci   lso   :switch turning on S/O coupling
Ci   ccd:  :diagonal matrices for 1-, 2-, 3-center terms for the
Ci         :combined corr;they are terms in parentheses in Kanpur notes
Ci         :eq. 3.87 multiplied by w^2. (see makdia.f)
Ci   vmtz  :muffin-tin zero for combined correction (asamad.f)
Ci   ldim  :dimension of the hamiltonian
Ci   lihdim:dimensions ccd and pph
Ci   sk    :structure constants, s^beta
Ci   ok    :S^beta-dot
Ci   hk    :?i-wave 3-centre integrals (see i3cntre)
Ci   wk    :work array of length ldim*2
Co Outputs
Co   hk,ok,:(mode 0) hamiltonian-> hk and overlap -> ok
Co         :(mode 1) small h -> hk
Co         :(mode 2) 1+oh -> hk
Co   sk    :(mode 0) sk is changed to noncol sqrdel*sk*sqrdel
Co         :(mode>0) sk is changed to noncol (C-enu)+sqrdel*sk*sqrdel
Cr Remarks
Cr
Cr   Inside the the sphere, a basis function is a linear combination
Cr   of phi's and dot's (including downfolded orbitals):
Cr     | psi> = | phi_L> + | phidot_L> h_LL + | phi_I> (h)_IL
Cr            = | phi_L>(1+oh_LL) + |dot_L> h_LL + | phi_I> (h)_IL
Cr   The first form uses phidot = phidot^alpha; the second form uses
Cr     phidot^alpha = phidot^gamma + o phi; and calls 'dot'=phidot^gamma
Cr   Note that <phi|phi>=1 <phidot|phi>=o <phidot|phidot>=p
Cr   Considering the LL block only, the ASA part of the overlap is:
Cr    <psi|psi>_ASA = <phi|phi> + h<phidot|phi> + h.c.
Cr                    + h <phidot|phidot> h
Cr                  = 1 + ho + oh + hph
Cr
Cr   To work directly with  D = srdel S srdel, rather
Cr   that h = C-enu + D, the diagonal parts connecting C-enu
Cr   in the one-, two-, three-center terms are reshuffled.  Thus
Cr    <psi|psi>_ASA = 1 + ho + oh + hph
Cr                  = 1 + (C-e+D)o + o(C-e+D) + (C-e+D) p (C-e+D)
Cr                  = 1 + 2(C-e)o + (C-e)^2 p      (one center)
Cr                  + D(o+p(C-e)) + (o+p(C-e))D    (two center)
Cr                  + D p D                        (three center)
Cr
Cr   The hamiltonian <psi|H|psi>_ASA has a corresponding structure
Cr   with similar 1-, 2- and 3- center terms; but the diagonal parts
Cr   are calculated from <phi|H|phi>, <phidot|H|phi>, <phidot|H|phidot>.
Cr
Cr
Cr   In the noncollinear case the same formulas apply, with the
Cr   following extensions:
Cr     1.  potential parameters have a ++ and a -- part
Cr     2.  D is a 2x2 matrix in spin space
Cr     3.  There may be a spin-orbit hamiltonian which has a structure
Cr         similar to the noncollinear hamiltonian, i.e. terms one-,
Cr         two- and three-center in D, but the diagonal bits couple
Cr         (l,m) to (l,m+/-1).
Cr     4.  Contribution from the applied field.  This is done separately
Cr         in routine hmladb.f, which see.
Cr
Cr   ?Routine untested for ldim ne lihdim (downfolding not implemented)
Cr   <k|k> = <kappa|kappa>
Cu Updates
Cu   15 Nov 07 New mode, to generate noncollinear h or 1+oh
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nbas,nl,indxsh(*),lihdim,ldim,neul,mode
      logical lasa,ccor,lss,lnc,lso
      double precision ccd(lihdim,0:2),pph(5,lihdim,2),wk(ldim,2),
     .  eula(nbas,neul,3),sk(ldim,2,ldim,2*2),hk(ldim,2,ldim,2*2),
     .  ok(ldim,2,ldim,2*2),sod(ldim,4,3),vmtz,qspirl(4)
C Local parameters
      integer i,j,l2,i1,j1,nasa,iprint,ncsw
      double precision xx

      call rx('hmltnc: noncollinear not installed')
      end

