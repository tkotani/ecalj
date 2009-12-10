      subroutine hml2nc(nbas,nl,indxsh,qspirl,eula,neul,
     .  pph,ccor,lss,lnc,ccd,wk,vmtz,elin,ldim,sk,sd,hk)
C- Generate two-center ASA noncollinear hamiltonian
C ---------------------------------------------------------------------
Ci Inputs
Ci   ccor,lss,lnc switches for comb. corr, spin spiral and noncollinear
Ci   qspirl(4): rotation angle of spin spiral
Ci   eula:  euler angles of noncollinear hamiltonian
Ci   pph:  vector of potential parameters (see makpph)
Ci   ccd, diagonal matrices constant, linear and
Ci   bilinear in S^alpha: they are the terms in parentheses in eq.3.87
Ci   in the Kanpur notes multiplied by w^2; vmtz, muffin-tin zero;
Ci   ldim: dimension of the hamiltonian
Ci   sk,sd: structure constants (sdot used only if ccor is true).
Ci   wk: work array of length ldim*2
Co Outputs
Co   hk
Cr Remarks
Cr   downfolding not implemented.
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer nbas,neul,nl,indxsh(*),ldim
      logical ccor,lss,lnc
      double precision ccd(ldim,0:2),eula(nbas,neul,3),pph(5,ldim,2),
     .  sk(ldim,ldim,2,2),sd(ldim,ldim,2,2),hk(ldim,2,ldim,2,2),
     .  vmtz,elin,qspirl(4),wk(ldim,2)
C Local parameters
      integer i,j,ncsw,i1,owk
      double precision xx
C heap:
      integer w(1)
      common /w/ w

      call rx('hml2nc not installed')

      end

