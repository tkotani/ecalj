      subroutine mksod(ldim,lihdim,ndim,indxsh,pph,soph,sod)
C- Diagonal arrays for one-, two- and three- center terms in the SO ham
C ---------------------------------------------------------------------
Ci Inputs
Ci   ldim  :lower dimension of hamiltonian matrix (makidx.f)
Ci   lihdim:dimensions pph
Ci   ndim  :dimensions of indxsh and soph; total no. of orbitals.
Ci   indxsh:permutations ordering orbitals in l+i+h blocks (makidx.f)
Ci   pph   :potential parameters in downfolding order (makpph.f)
Ci   soph (see remarks)
Co Outputs
Co   sod(*,1..3,isp): arrays for 1-, 2-, 3- center terms, ++,-- blocks
Co   sod(*,1,3): for 1- center terms,    L- S+ block
Co   sod(*,2..3,3): for 2- center terms, L- S+ block
Co   sod(*,4,3):    for 3- center terms, L- S+ block
Cr Remarks
Cr   Makes diagonal matrices belonging to the 1-, 2- and 3-center
Cr   terms in the matrix elements of the SO hamiltonian.  Its form is
Cr   (H_so)_(i,j) = <psi_i| I L . S | psi_j>.
Cr   where L.S = Lz.Sz + .5(L+.S- + L-.S+) and I = (2/(c^2) dV/dr*(1/r).
Cr   Inside the the sphere, a basis function is a linear combination
Cr   of phi's and dot's (including downfolded orbitals):
Cr     | psi> = | phi_L> + | phidot_L>^alpha h_LL + | phi_I> (h)_IL
Cr            = | phi_L>(1+oh_LL) + |dot_L> h_LL + | phi_I> (h)_IL
Cr   This subroutine combines radial matrix elements I between phi's
Cr   and dot's (phidot^gamma) and the L.S operator into diagonal
Cr   matrices associated with the 1-, 2- and 3-center terms comprising
Cr   matrix elements (H_so)_(i,j).
Cr   For the (L-S+) block, the SO coupling uses:
Cr   L+ |lm> = sqrt((l-m)*(l+m+1)) |l,m+1>
Cr   L- |lm> = sqrt((l+m)*(l-m+1)) |l,m-1>
Cr   <phi1_l,m | L-.S+ | phi2_l,m+1> = F/2 < phi1_l | I | phi2_l>
Cr   where F = sqrt((l+m+1)*(l-m))
Cr   L-.S+ couples l,m only to l,m+1 (2l elements for each l).  The
Cr   following evaluates (2l+1) elements foreach l, exploiting the
Cr   fact that <lm | L-S+ | lm+1> is zero for m=l because F=0, and
Cr   similarly <lm-1 | L-S+ | lm> is zero for m=-l.
Cr
Cr   To work directly with  D = srdel S srdel, rather
Cr   that h = C-enu + D, the diagonal parts connecting C-enu
Cr   in the one-, two-, three-center terms are reshuffled.  Thus
Cr   w1 + w2 (C-e+D) + h.c. + (C-e+D) w3 (C-e+D)  becomes
Cr   w1+2*w2*(C-e)+w3(C-e)**2 + (w2+(C-e)w3) D + h.c. + D w3 D.
C ---------------------------------------------------------------------
C     implicit none
      integer ldim,lihdim,ndim,indxsh(ndim)
      double precision pph(5,lihdim,2),soph(4,2,2,ndim),sod(ldim,4,3)
      double precision sp,sz,flmp,flmm
      parameter (sp=.5d0)
      integer i,j,k,isp

      call rx('SO coupling not installed')

      end


