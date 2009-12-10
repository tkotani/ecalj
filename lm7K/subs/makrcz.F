      subroutine makrcz(ldim,idim,nev,pph,sil,zll,zil)
C- Make rectangular eigenvectors for intermediate waves in beta rep'n
C-----------------------------------------------------------------------
Ci Inputs
Ci   ldim  :dimension of hamiltonian matrix (makidx.f)
Ci   idim  :dimension of block of orbitals to be folded down (makidx.f)
Ci   nev   :actual number of eigenvectors generated
Ci   pph   :potential parameters in downfolding order (makpph.f)
Ci         :for the current spin.  Not used in noncollinear case.
Ci   sil   :structure constants (beta rep'n)
Ci   zll   :(ll-block of eigenvectors in beta rep'n).  In noncollinear
Ci         :case, zll should be U d_l zll; see Remarks.
Co Outputs
Co   zil   :i-l block of eigenvectors in beta rep'n.  In noncollinear
Co         :case, zil is returned as s_il U d_l z_il; see Remarks.
Co   sil   :is decorated with potential parameters; see Remarks
Co   pph   :iwaves are overwritten (not touched in noncollinear case)
Cr Remarks
Cr   Intermediate wave eigenvectors z_i are obtained from the lower
Cr   wave ones as follows (NB the 1st term RHS is rep'n independent)
Cr      H_il = d_i S_il d_l
Cr      z_il = H_il  z_ll
Cr   where
Cr      d_i = sqrdel_i/(C - e_nu)_i
Cr      d_l - sqrdel_l
Cr   In the noncollinear case, must substitute S_il -> U+ S_il U
Cr   but in practice this is done not by calling makrcz but by 
Cr   operation right-to-left matrix operations:
Cr     z_il <-  d_i U_i+ S_il U_l d_l z_ll
Cr   See addhnq.f
C-----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer ldim,idim,nev
      double precision pph(5,*),sil(idim,ldim,2),
     .                 zll(ldim,ldim,2),zil(idim,ldim,2)

C Local parameters
      integer i,l2,li

      call tcn('makrcz')

      l2 = ldim**2
      li = ldim * idim
      do  1  i = ldim+1, ldim+idim
        pph(3,i) = pph(3,i) / ( pph(2,i) - pph(1,i) )
    1 continue
      call makdsd(0,idim,ldim,idim,ldim,0,ldim,pph,sil,sil)
      call zmpy(sil,idim,1,li,zll,ldim,1,l2,zil,idim,1,li,
     .          idim,nev,ldim)
      call tcx('makrcz')
      end

