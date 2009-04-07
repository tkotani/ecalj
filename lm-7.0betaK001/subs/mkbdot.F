      subroutine mkbdot(ldim,idim,lihdim,lavg,ccd,pph,avw,adot)
C- make w^-2(beta^dot - alpha^dot) for i-waves to fold down
C ----------------------------------------------------------------------
Ci Inputs
Ci   ldim  :dimension of hamiltonian matrix (makidx.f)
Ci   idim  :dimension of block of orbitals to be folded down (makidx.f)
Ci   lihdim:dimensions ccd and pph
Ci   lavg  :F use pph(1)
Ci         :T use (pph(1)+pph(2))/2
Ci   ccd   :ccd(*,2) = bilinear terms in eq. 3.87 Kanpur notes
Ci   pph   :potential parameters in downfolding order (makpph.f)
Ci   avw   :length scale, usu. average Wigner-Seitz sphere radius
Ci   adot  :(kappa*avw)^2-derivative of tight-binding screening const.
Co Output:
Co   adot: (intermediate waves) changed to w^-2(beta^dot-alpha^dot)
Cr Remarks:
Cr   The condition of downfolding that 3-centre overlap integrals over
Cr   the i-waves should vanish is
Cr      w^-2(beta^dot - alpha^dot) + D(2) + w^-2 DEL / (C - e_nu)^2 = 0
Cr   where D(2) is the diagonal matrix in parentheses in the fourth
Cr   term in eq. 3.87 Kanpur notes.
Cr   NB the last term is representation - independent
Cu Updates
Cu   14 Sep 99 pass entire ccd array through; changed arg list.
Cu   17 Jan 00 add lavg argument
C-----------------------------------------------------------------------
C     implicit none
      logical lavg
      integer ldim,idim,lihdim
      double precision ccd(lihdim,0:2),pph(5,lihdim,2),adot(lihdim),avw
      integer i
      double precision xx

      do  1  i = ldim+1, ldim+idim
        if (lavg) then
          xx = ((pph(3,i,1)/(pph(2,i,1)-pph(1,i,1)))**2
     .       +  (pph(3,i,2)/(pph(2,i,2)-pph(1,i,2)))**2)/2
        else
          xx = (pph(3,i,1)/(pph(2,i,1)-pph(1,i,1)))**2
        endif
        adot(i) = - ( ccd(i,2) + xx ) / avw**2
    1 continue
      end

