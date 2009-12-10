      subroutine mkbdia(ldim,lihdim,ndim,indxsh,pph,soph,sod)
C- Diagonal arrays for one-, two- and three- center terms applied B
C ---------------------------------------------------------------------
Ci Inputs
Ci   ldim  :lower dimension of hamiltonian matrix (makidx.f)
Ci   lihdim:dimensions pph
Ci   ndim  :dimensions of indxsh and soph; total no. of orbitals.
Ci   indxsh:permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci
Ci   pph   :potential parameters in downfolding order (makpph.f)
Ci         :pph(1..5,i,s): parameters for ith RL and sth spin channel.
Ci         :pph(1) : enu    
Ci         :pph(2) : calpha 
Ci         :pph(3) : sqrdel 
Ci         :pph(4) : palpha 
Ci         :pph(5) : oalp   
Ci
Ci   soph  :Spin-orbit parameters, in downfolding order (mksoph.f)
Ci         :soph(*,*,*,i) are parameters for ith RL channel.
Ci         :soph(1..3,s1,s2,i) are radial integrals for the spin-orbit
Ci         :              perturbation (computed in soprm.f). s1 and s2
Ci         :              are spin indices; 1st index is for cases:
Ci         :              <phi Hso phi>, <phi Hso dot>, <dot Hso dot>
Ci         :soph(4,1,1,i) l quantum number
Ci         :soph(4,2,1,i) m quantum number
Co Outputs
Co   sod(*,1..3,isp): arrays for 1-, 2-, 3- center terms, ++,-- blocks
Co   sod(*,1,3):    for 1- center terms, L- S+ block
Co   sod(*,2..3,3): for 2- center terms, L- S+ block
Co   sod(*,4,3):    for 3- center terms, L- S+ block
Cr Remarks
Cr   Makes diagonal matrices belonging to the 1-, 2- and 3-center
Cr   terms in the matrix elements of the SO hamiltonian.  Its form is
Cr   (H_so)_(i,j) = <psi_i| I L . S | psi_j>.
Cr   where L.S = Lz.Sz + .5(L+.S- + L-.S+) and I = (2/(c^2) dV/dr*(1/r).
Cr   Inside the the sphere, a basis function is a linear combination
Cr   of phi's and dot's (including downfolded orbitals):
Cr     | psi> = | phi_L> + | phidot_L> h_LL + | phi_I> (h)_IL
Cr            = | phi_L>(1+oh_LL) + |dot_L> h_LL + | phi_I> (h)_IL
Cr   The first form uses phidot = phidot^alpha; the second form uses
Cr     phidot^alpha = phidot^gamma + o phi; and calls 'dot'=phidot^gamma
Cr   This subroutine combines radial matrix elements I between phi's
Cr   and dot's (phidot^gamma) and the L.S operator into diagonal
Cr   matrices associated with the 1-, 2- and 3-center terms comprising
Cr   matrix elements (H_so)_(i,j).
Cr   For the (L-S+) block, the SO coupling uses:
Cr     L+ |lm> = sqrt((l-m)*(l+m+1)) |l,m+1>
Cr     L- |lm> = sqrt((l+m)*(l-m+1)) |l,m-1>
Cr     <phi1_l,m | L-.S+ | phi2_l,m+1> = F/2 < phi1_l | I | phi2_l>
Cr   where F = sqrt((l+m+1)*(l-m))
Cr   L-.S+ couples l,m only to l,m+1 (2l elements for each l).  The
Cr   following evaluates (2l+1) elements foreach l, exploiting the
Cr   fact that <lm | L-S+ | lm+1> is zero for m=l because F=0, and
Cr   similarly <lm-1 | L-S+ | lm> is zero for m=-l.
Cr
Cr   To work directly with  D = srdel S srdel, rather
Cr   that h = C-enu + D, the diagonal parts wi connecting C-enu
Cr   in the one-, two-, three-center terms are reshuffled.  Thus
Cr     w1 + w2 (C-e+D) + h.c. + (C-e+D) w3 (C-e+D)  becomes
Cr     w1+2*w2*(C-e)+w3(C-e)**2 + (w2+(C-e)w3) D + h.c. + D w3 D.
C ---------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer ldim,lihdim,ndim,indxsh(ndim)
      double precision pph(5,lihdim,2),soph(4,2,2,ndim),sod(ldim,4,3)
C ... Local parameters
      double precision sp,sz,flmp,flmm
      parameter (sp=.5d0)
      integer i,j,k,isp
      flmp(i) = dsqrt((soph(4,1,1,i)+soph(4,2,1,i)+1)*
     .                (soph(4,1,1,i)-soph(4,2,1,i)))
      flmm(i) = dsqrt((soph(4,1,1,i)+soph(4,2,1,i))*
     .                (soph(4,1,1,i)-soph(4,2,1,i)+1))

C ------ Make up-up, down-down blocks of S-O hamiltonian ------
      do  10  isp = 1, 2
        sz = 1.5d0 - isp

C   ... For 3C  ho <phi H phi> oh + 2 ho <phi H dot> h + h <dot H dot> h
        do  16  k = 1, ndim
          i = indxsh(k)
          if (i .gt. ldim) goto 16
          sod(i,3,isp) =
     .      soph(4,2,1,i)*sz*soph(1,isp,isp,i)*pph(5,i,isp)**2 +
     .      soph(4,2,1,i)*sz*soph(2,isp,isp,i)*pph(5,i,isp)*2 +
     .      soph(4,2,1,i)*sz*soph(3,isp,isp,i)
   16   continue

C   ... For 2C <phi H_so phi> oh + <phi H_so dot> h + h.c.
        do  14  k = 1, ndim
          i = indxsh(k)
          if (i .gt. ldim) goto 14
          sod(i,2,isp) =
     .      soph(4,2,1,i)*sz*soph(1,isp,isp,i)*pph(5,i,isp) +
     .      soph(4,2,1,i)*sz*soph(2,isp,isp,i)
   14   continue

C   ... For 1C  <phi H_so phi>
        do  12  k = 1, ndim
          i = indxsh(k)
          if (i .gt. ldim) goto 12
          sod(i,1,isp) = soph(1,isp,isp,i)*soph(4,2,1,i)*sz
   12   continue

C ...   Add diagonal (C-enu) subtracted from 2- and 3-center

        do  18  k = 1, ndim
          i = indxsh(k)
          if (i .gt. ldim) goto 18
          sod(i,1,isp) = sod(i,1,isp) +
     .      sod(i,3,isp)*(pph(2,i,isp)-pph(1,i,isp))**2 +
     .      2*sod(i,2,isp)*(pph(2,i,isp)-pph(1,i,isp))
   18   continue
        do  19  k = 1, ndim
          i = indxsh(k)
          if (i .gt. ldim) goto 19
          sod(i,2,isp) = sod(i,2,isp) +
     .      sod(i,3,isp)*(pph(2,i,isp)-pph(1,i,isp))
   19   continue
   10   continue

C ... For 3C contribution to L-S+ block
      do  22  k = 1, ndim-1
        i = indxsh(k)
        if (i .gt. ldim) goto 22
        sod(i,4,3) =
     .    (pph(5,i,1)*(soph(1,1,2,i)*pph(5,i+1,2) + soph(2,1,2,i))
     .    + soph(2,2,1,i)*pph(5,i+1,2) + soph(3,1,2,i))*sp*flmp(i)
   22 continue
      sod(ndim,4,3) = 0

C ... For 2C contribution to L-S+ block
      do  24  k = 1, ndim-1
      i = indxsh(k)
      if (i .gt. ldim) goto 24
      sod(i,2,3) = sp*flmp(i)*(soph(1,1,2,i)*pph(5,i,1) + soph(2,2,1,i))
   24 continue
      do  25  k = 2, ndim
      i = indxsh(k)
      if (i .gt. ldim) goto 25
      sod(i,3,3) = sp*flmm(i)*(soph(1,1,2,i)*pph(5,i,2) + soph(2,1,2,i))
   25 continue
      sod(ndim,2,3) = 0
      sod(1,3,3) = 0

C ... For 1C contribution to L-S+ block
      do  26  k = 1, ndim-1
      i = indxsh(k)
      if (i .gt. ldim) goto 26
      sod(i,1,3) = soph(1,1,2,i)*flmp(i)*sp
   26 continue
      sod(ndim,1,3) = 0

C ... Add diagonal (C-enu) subtracted from 2- and 3-center
      if (.true.) then
        do  28  k = 1, ndim-1
          i = indxsh(k)
          if (i .gt. ldim) goto 28
          sod(i,1,3) = sod(i,1,3) +
     .    (pph(2,i+1,1)-pph(1,i+1,1))*sod(i,2,3) +
     .    (pph(2,i+1,2)-pph(1,i+1,2))*sod(i+1,3,3) +
     .    (pph(2,i,1)-pph(1,i,1))*(pph(2,i+1,2)-pph(1,i+1,2))*sod(i,4,3)
   28   continue
        do  29  k = 1, ndim-1
        i = indxsh(k)
        if (i .gt. ldim) goto 29
        sod(i,2,3) = sod(i,2,3) + (pph(2,i+1,2)-pph(1,i+1,2))*sod(i,4,3)
   29   continue
        do  30  k = 2, ndim
        i = indxsh(k)
        if (i .gt. ldim) goto 30
        j = i-1
        sod(i,3,3) = sod(i,3,3) + (pph(2,j,1)-pph(1,j,1))*sod(j,4,3)
   30   continue
      endif

C     call yprm('sod',2,sod,0,ldim,ldim,12)

      end


