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
Ci   soph  :matrix elements coupling up- and down- spatial wave functions
Ci   sk    :scaled structure constants, sqrdel*s^beta*sqrdel
Ci   bdots :(magnetic field . Pauli matrix), in downfolding order
Co Outputs
Co   hk    :Magnetic field term B.sigma is added to hk
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

C     print *, '!!h'
C     hk = 0

      ldimx = ldim*2

C --- Make the diagonal matrices B0, B1, B2 ---
      do  k = 1, ndim

        i = indxsh(k)
        if (i .gt. ldim) goto 16
        do  is1 = 1, 2
          o(is1)   = pph(5,i,is1)
          CmE(is1) = pph(2,i,is1)-pph(1,i,is1)
        enddo

        do  is1 = 1, 2
        do  is2 = 1, 2

C     ... Spatial parts of the matrix elements
          s00 = soph(4,is1,is2,i)
C         <phi|dot + o phi> = <phi|dot> + <phi|phi> o
          s01 = soph(5,is1,is2,i) + soph(4,is1,is2,i)*o(is2)
C         <dot + o phi|phi> = <dot|phi> + o <phi|phi>
C                           = <phi|dot>(+) + (<phi|phi> o)(+)
C                           = <phi|dot>(+) + o <phi|phi>
          s10 = soph(5,is2,is1,i) + o(is1)*soph(4,is1,is2,i)
C         <dot + o phi|dot + o phi> =
C         <dot|dot> + o <phi|dot> + <dot|phi> o + o <phi|phi> o =
C         <dot|dot> + o <phi|dot> + [o <phi|dot>](+) + o <phi|phi> o
          s11 = soph(6,is1,is2,i)
     .        + o(is1)*soph(5,is1,is2,i)
     .        + o(is2)*soph(5,is2,is1,i)
     .        + o(is1)*soph(4,is1,is2,i)*o(is2)

C     ... Spatial parts * spin parts
          b00 = bdots(is1,is2,i)*s00
          b01 = bdots(is1,is2,i)*s01
          b10 = bdots(is1,is2,i)*s10
          b11 = bdots(is1,is2,i)*s11

          B0(i,is1,is2) = b00 + CmE(is1)*b10 + b01*CmE(is2)
     .                  + CmE(is1)*b11*CmE(is2)

          B1(i,is1,is2) = b01 + CmE(is1)*b11
C         Bd should be the hermitian conjugate of B1.
          Bd(i,is1,is2) = b10 + b11*CmE(is2)

          B2(i,is1,is2) = b11


C          B0(i,is1,is2) = 0
C          B0(i,1,1) = CmE(1)
C          B0(i,2,2) = CmE(2)
C          B0(i,is1,is2) = s11

        enddo
        enddo

   16   continue
      enddo

C     call zprm('B0',2,B0,ldim,ldim,4)
C     call zprm('B1',2,B1,ldim,ldim,4)
C     call zprm('Bd',2,Bd,ldim,ldim,4)
C     call zprm('B2',2,B2,ldim,ldim,4)

#if F90
      allocate(wk(ldim,2,ldim,2))
#endif

C --- H += D B2 D ---
      call dpzero(wk,ldim**2*8)
      do  i1 = 1, 2
      do  j1 = 1, 2
      do  k1 = 1, 2
        do  j = 1, ldim
        do  i = 1, ldim
          wk(i,i1,j,j1) = wk(i,i1,j,j1) + B2(i,i1,k1)*sk(i,k1,j,j1)
        enddo
        enddo
      enddo
      enddo
      enddo
C     call zprm('wk in hmladb',2,wk,ldimx,ldimx,ldimx)
C     call zprm('s in hmladb',2,sk,ldimx,ldimx,ldimx)
      call zgemm('N','N',ldimx,ldimx,ldimx,(1d0,0d0),sk,ldimx,wk,ldimx,
     .           (1d0,0d0),hk,ldimx)
C     call zprm('<B.s>(3C)',2,hk,ldimx,ldimx,ldimx)

C --- H += B1 D + D B1(+) ---
      do  i1 = 1, 2
      do  j1 = 1, 2
      do  k1 = 1, 2
        do  j = 1, ldim
        do  i = 1, ldim
          hk(i,i1,j,j1) = hk(i,i1,j,j1)
     .                  + B1(i,i1,k1)*sk(i,k1,j,j1)
     .                  + sk(i,i1,j,k1)*Bd(j,k1,j1)
        enddo
        enddo
      enddo
      enddo
      enddo
C     call zprm('<B.s>(2C)',2,hk,ldimx,ldimx,ldimx)

C --- H += B0 ---
      do  i1 = 1, 2
      do  j1 = 1, 2
        do  i = 1, ldim
          hk(i,i1,i,j1) = + hk(i,i1,i,j1) + B0(i,i1,j1)
        enddo
      enddo
      enddo

#if F90
      deallocate(wk)
#endif

      if (iprint() .gt. 110) then
        call zprm('H + <B.s>',2,hk,ldimx,ldimx,ldimx)
      endif

      end

