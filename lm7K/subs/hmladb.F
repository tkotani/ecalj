#define F90 1
      subroutine hmladb(ldim,lihdim,ndim,indxsh,pph,bdots,soph,sk,hk)
C- Add ASA magnetic hamiltonian
C ---------------------------------------------------------------------
Ci Inputs
Ci   ldim  :lower dimension of hamiltonian matrix (makidx.f)
Ci   lihdim:dimensions pph,soph,bdots
Ci   ndim  :dimensions of indxsh; total no. of orbitals.
Ci   indxsh:permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci
Ci   pph   :potential parameters in downfolding order (makpph.f)
Ci         :pph(1..5,i,is): parms for ith RL and is(th) spin channel.
Ci         :pph(1) : enu
Ci         :pph(2) : calpha
Ci         :pph(3) : sqrdel
Ci         :pph(4) : palpha
Ci         :pph(5) : oalp
Ci
Ci   bdots :(magnetic field . Pauli matrix), in downfolding order
Co Outputs
Co   bsd   :bsd(i,1..3,is1,is2): diagonal arrays for 1-, 2-, 3- center
Co         :terms for the matrix elements of (b.s) in the ith RL
Co         :channel; see Remarks.
Cr Remarks
Cr *Note: this routine uses complex*16 arithmetic.
Cr
Cr *It makes diagonal matrices belonging to the 1-, 2- and 3-center
Cr  terms in the matrix elements of the magnetic field.  Its form is
Cr  (B.s)_(i,j) = <psi_i| B . Sigma | psi_j>.
Cr
Cr  Inside the the sphere, a basis function is a linear combination
Cr  of phi's and dot's (including downfolded orbitals):
Cr    | psi> = | phi_L> + | phidot_L> h_LL + | phi_I> (h)_IL
Cr           = | phi_L>(1+oh_LL) + |dot_L> h_LL + | phi_I> (h)_IL
Cr  The first form uses phidot = phidot^alpha; the second form uses
Cr    phidot^alpha = phidot^gamma + o phi; and calls 'dot'=phidot^gamma
Cr  This subroutine combines radial matrix elements between phi's and
Cr  dot's (phidot^gamma) and the B.sigma operator into diagonal matrices
Cr  associated with the 1-, 2- and 3-center terms comprising matrix
Cr  elements (B.s)_(i,j).  The field is assumed constant within each
Cr  sphere, so the matrix elements (B.s)_(i,j) are very similar to the
Cr  matrix elements of the unit operator (overlap matrix), which is:
Cr   <psi|psi>_ASA = <phi|phi> + h<phidot|phi> + h.c.
Cr                   + h <phidot|phidot> h
Cr                 = 1 + ho + oh + hph
Cr
Cr  To work directly with  D = srdel S srdel, rather
Cr  that h = C-enu + D, the diagonal parts of both the overlap and B.s
Cr  are reshuffled.  Thus we get
Cr   <psi|psi>_ASA = 1 + ho + oh + hph
Cr                 = 1 + (C-e+D)o + o(C-e+D) + (C-e+D) p (C-e+D)
Cr                 = 1 + 2(C-e)o + (C-e)^2 p      (one center)
Cr                 + D(o+p(C-e)) + (o+p(C-e))D    (two center)
Cr                 + D p D                        (three center)
Cr
Cr  Matrix elements of (B-sigma) in the ASA:
Cr
Cr  Wave function is  phi + (C-enu + D) phidot, where
Cr    phidot = phidot^alpha = phidot^gamma + o phi = dot + o phi
Cr    D = U sqrt(delta) S sqrt(delta) U+,
Cr
Cr  and U being the rotation matrix mixing spin components of the
Cr  tails.  Writing phi= phi^0 and phidot=phi^1, the ASA matrix
Cr  elements are
Cr
Cr  ME =
Cr  <(phi^0+(C-enu+D)phi^1)_s1 |(B.s)| (phi^0+(C-enu+D) phi^1)_s2 >
Cr
Cr  Write
Cr  <B.s>^ij_s1,s2 = <phi^i_s1|(B.s)|phi^j_s2>
Cr                 = <B.s>_s1,s2 <phi^i_s1|phi^j_s2>
Cr                 = <B.s>_s1,s2 o^ij_s1,s2
Cr
Cr  Then
Cr  ME = <B.s>^00_s1,s2
Cr     + (C-enu)_s1 <B.s>^10_s1,s2 + <B.s>^01_s1,s2 (C-enu)_s2
Cr     + (C-enu)_s1 <B.s>^11_s1,s2 (C-enu)_s2
Cr     + sum_s3 D_s1,s3 (<B.s>^10_s3,s2 + <B.s>^11_s3,s2 (C-enu)_s2)
Cr     + sum_s3 (<B.s>^01_s1,s3 + (C-enu)_s1 <B.s>^11_s3,s2) D_s3,s2
Cr     + sum_s3,s4 D_s1,s3 <B.s>^11_s3,s4 D_s4,s1
Cr
Cr  Suppressing spin indices (2x2 matrices in spin space) we can write
Cr  ME = B0 + B1 D + D B1(+) + D B2 D  where
Cr
Cr  Here (+) means complex conjugate, and exchange all indices:  Note
Cr  <B.s>^ij_s1,s2 = <B.s>_s1,s2 <phi^i_s1|phi^j_s2>
Cr                 = (<B.s>_s2,s1)* <phi^j_s2|phi^i_s1>*
Cr                 = <B.s>^ji_s2,s1 = [<B.s>^ij_s1,s2](+)
Cr
Cr  B0 = <B.s>^00 + (C-enu) <B.s>^10 + <B.s>^01 (C-enu)
Cr     + (C-enu) <B.s>^11 (C-enu)
Cr     = <B.s>^00 + (C-enu) <B.s>^10 + h.c. + (C-enu) <B.s>^11 (C-enu)
Cr     because <B.s>^10_s1,s2 = (<B.s>^10_s1,s2)+ = <B.s>^01_s2,s1
Cr     since operator (B.s) is hermitian.
Cr
Cr  B1 = <B.s>^01 + (C-enu) <B.s>^11
Cr
Cr  B1(+) = <B.s>^10 + <B.s>^11 (C-enu)
Cr
Cr  B2 = <B.s>^11
Cr
Cr  If (B.s)_s1,s2 is replaced by the unit operator, delta_s1,s2,
Cr  the above expression reduces to the ASA overlap matrix:
Cr    <B.s>^00_s1,s2 -> delta_s1,s2
Cr    <B.s>^10_s1,s2 -> o_s1,s2 delta_s1,s2 = o
Cr    <B.s>^01_s1,s2 -> o_s2,s1 delta_s1,s2 = o
Cr    <B.s>^11_s1,s2 -> p_s1,s2 delta_s1,s2 = p
Cr  So B0, B1, B2 become diagonal in spin space with
Cr    B0 -> 1 + (C-enu) o + o (C-enu) + (C-enu) p (C-enu)
Cr    B1 -> o + (C-enu) p = (B1)+
Cr    B2 -> p
Cr
Cr *For checking: create these files:
Cr     sij = <phi^i|phi^j> : i=0..1, j=0..1; each an ldim x 4 matrix
Cr     o   = overlap, pph(5)
Cr     p   = <phidot|phidot> delta_s1,s2 (small parameter, pph(4))
Cr     cme = (pph(2)-pph(1)) delta_s1,s2
Cr     O1c = diagonal contr. to overlap matrix
Cr     O2c = diagonal array mult. D for 2C overlap matrix
Cr     (this files will change with choice of magnetic field)
Cr     B0  = 1-center contr. to <B.S>
Cr  Checks to make for space part of matrix elements:
Cr  1.  show that s11(diagonal) = p
Cr      mc s11 -real p -real -- -coll 1,4 -px
Cr  2.  show that s11(off-diagonal) ~ geometric mean of diagonal
Cr      mc s11 -real -e2 'sqrt(x1*x4)' x2 -e1 x1/x2
Cr  3.  show that s01(diagonal) = o
Cr      mc o s01 -- -coll 1,4 -px
Cr  4.  show that s01(off-diagonal) ~ geometric mean of diagonal
Cr      mc s01 -real -e2 'x1*x4' 'x2*x3' -e1 x1/x2
Cr  5.  show that s01 = s10(+)
Cr      mc s01 s10 -coll 1,3,2,4 -- -px
Cr  ... in these checks, make bfield=0,0,1 so that result should
Cr      be diagonal in sigma, and (11)=o(+)/2, (22)=o(-)/2
Cr  6.  mc -f10f15.10 ../out.mnpt -s2 -real -e2 x1 -x4 O$n{}c --
Cr  ... in these checks, make bfield=1,0,0; now (12,21) almost (11,-22)
Cr  6.  mc -f10f15.10 B0 -s2 -coll 2,3 O$n{}c -- -real -e3 x1 x2 x1+x2
Cr  ... in these checks, make bfield=0,1,0; now (12,21) almost (11,-22)
Cr  6.  mc -f10f15.10 B0 -s2 -coll 2,3 -s0,-1 -real -e2 x1 -x2 O$n{}c -real --
Cu Updates
Cu   14 Feb 03 First created
C ---------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer ldim,lihdim,ndim,indxsh(ndim)
      double precision pph(5,lihdim,2),soph(7,2,2,lihdim)
      double complex bdots(2,2,lihdim)
      double complex sk(ldim,2,ldim,2),hk(ldim,2,ldim,2)
C ... Local parameters
      integer i,j,k,ldimx,is1,is2,i1,j1,k1,iprint
      double precision o(2),CmE(2),s00,s01,s10,s11
      double complex b00,b01,b10,b11
      double complex B0(ldim,2,2),B1(ldim,2,2),Bd(ldim,2,2),B2(ldim,2,2)
#if F90
      complex(8),allocatable:: wk(:,:,:,:)
#else
      double complex wk(ldim,2,ldim,2)
#endif

      call rx('Applied B-field not installed')

      end

