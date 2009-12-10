      subroutine makbet(ndim,ldim,idim,indxsh,pti,wk,alpha)
C- make screening constants for downfolding
C ----------------------------------------------------------------------
Ci Inputs
Ci   ndim  :dimension of alpha and pti
Ci   ldim  :dimension of hamiltonian matrix (makidx.f)
Ci   idim  :dimension of block of orbitals to be folded down (makidx.f)
Ci   indxsh:permutation indices ordering orbitals in downfolding order 
Ci          (makidx.f)
Ci   pti   :inverse potential functions (suham.f,makpti.f)
Ci   wk    :work array of dimension ndim
Ci   alpha :tight-binding screening parameters
Co Outputs
Co   alpha over the i-waves is changed to [P^0(e_nu)]^-1
Cr Remarks
Cu Updates
C-----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer ndim,ldim,idim,indxsh(ndim)
      double precision pti(ndim),alpha(ndim),wk(ndim)

      call shffle(.true.,ndim,indxsh,alpha,wk)
      call dcopy(idim,pti(1+ldim),1,alpha(1+ldim),1)
      call shffle(.false.,ndim,indxsh,alpha,wk)
      end

