      subroutine hmltso(lasa,ccor,lso,ldim,idim,pph,
     .  sod,diags,ndim,vmtz,sk,hk,ok,wk)
C- ASA + spin-orbit hamiltonian, and overlap
C ---------------------------------------------------------------------
Ci Inputs
Ci   lasa,ccor,lso: switches telling whether to add (1) asa H,
Ci                  (2) combined correction and (3) spin-orbit H
Ci   ldim,idim, dimension of lower and intermediate waves
Ci   pph: vector of pot pars, make by makpph
Ci   sod: diagonal arrays for SO coupling, made by mksod:
Co   sod(*,1..3,isp): arrays for 1-, 2-, 3- center terms, ++,-- blocks
Co   sod(*,1,3): for 1- center terms,    L- S+ block
Co   sod(*,2..3,3): for 2- center terms, L- S+ block
Co   sod(*,4,3):    for 3- center terms, L- S+ block
Cr Remarks
Cr   (H_so)_(i,j) = <psi_i| I L . S | psi_j>.
Cr   where L.S = Lz.Sz + .5(L+.S- + L-.S+) and I = (2/(c^2) dV/dr*(1/r).
Cr   Inside the the sphere, a basis function is a linear combination
Cr   of phi's and dot's (including downfolded orbitals):
Cr     | psi> = | phi_L> + | phidot_L> h_LL + | phi_I> (h)_IL
Cr            = | phi_L>(1+oh_LL) + |dot_L> h_LL + | phi_I> (h)_IL
Cr   Subroutine mksod combines radial matrix elements I between phi's
Cr   and dot's (phidot^gamma) and the L.S operator into diagonal
Cr   matrices associated with the 1-, 2- and 3-center terms comprising
Cr   matrix elements (H_so)_(i,j).
Cr
Cr   The up-up and down-down (Lz.Sz) terms calculated first
Cr   the up-down (L-S+) calculated next; (L+S-) obtained from (L-S+).
Cr
Cr   The ASA hamiltonian works directly with D = srdel S srdel, rather
Cr   that h = C-e+D.  To combine SO with the ASA H, the diagonal parts
Cr   of the one-, two-, three-center terms are reshuffled.  Thus
Cr   w1 + w2 (C-e+D) + h.c. + (C-e+D) w3 (C-e+D)  becomes
Cr   w1+2*w2*(C-e)+w3(C-e)**2 + (w2+(C-e)w3) D + h.c. + D w3 D.
C ---------------------------------------------------------------------
C     implicit none
C Passed parameters
      logical lasa,ccor,lso
      integer ldim,idim,ndim
      double precision vmtz,
     .  sk(ldim,2,ldim,2),pph(5,ldim,2),diags(ndim,0:2),
     .  ok(ldim,2,ldim,2),hk(ldim,2,ldim,2),sod(ldim,4,3),wk(ldim,2)
C Local parameters
      integer i,j,i1,j1,l2,iprint,nasa

      call rx('hmltso not installed')
      end


