      subroutine hmlt2c(ccor,diags,vmtz,eln,ldim,ndim,pph,sk,sdot,wk)
C- Generate ASA and CCOR part of Hamiltonian, neglecting 3rd order terms
C ---------------------------------------------------------------------
Ci Inputs
Ci   ccor, logical to switch of combined correction
Ci   diags, diagonal matrices constant and linear in S for CC;
Ci   see Kanpur notes 3.87 (scaled by w^-2); see makdia
Ci   pph:  vectors of potential parameters (see makpph)
Ci   ndim: 2nd dimension of pph = nbas * nl^2
Ci   ldim: dimension of the hamiltonian
Ci   sdot: S-dot * ; sk: S
Ci   wk:work array of length ldim
Co Outputs
Co   hk; and sk is changed to sqrdel*s*sqrdel
Cr Remarks
Co   Combined correction is to first order in VMTZ-eln
C ----------------------------------------------------------------------
C     implicit none
C Passed parameters
      integer ndim,ldim
      logical ccor
      double precision pph(5,ndim),wk(ldim),diags(ndim,0:2),vmtz,eln
      double precision sk(ldim,ldim*2),sdot(ldim,ldim*2)
C Local parameters
      integer i,j,l2

      call tcn('hmlt2c')
      l2 = ldim**2

C --- Make d*S*d and d*Sdot*d with d = sqrt(delta) ---
      call makdsd(0,ldim,ldim,ldim,ldim,0,0,pph,sk,sk)
      if (ccor) then
        call makdsd(0,ldim,ldim,ldim,ldim,0,0,pph,sdot,sdot)
C       call prmx('sdot',sdot,ldim,ldim,ldim)
        call dscal(2*l2,vmtz-eln,sdot,1)
      endif

C --- H += d*S*d* + (vmtz-eln)*<k|k>_linear ---
      do  6  i = 1, ldim
    6 wk(i) = .5d0
      if (ccor) call daxpy(ldim,vmtz-eln,diags(1,1),1,wk,1)
      do  8  j = 1, ldim
      do  8  i = 1, ldim
        sk(i,j) = sdot(i,j) + sk(i,j)*(wk(i)+wk(j))
        sk(l2+i,j) = sdot(l2+i,j) + sk(l2+i,j)*(wk(i)+wk(j))
    8 continue

C --- H += c + (vmtz-eln)*<k|k>_constant ---
      call dcopy(ldim,pph(2,1),5,wk,1)
      if (ccor) then
        do  10  i = 1, ldim
   10   wk(i) = wk(i) + (vmtz-eln)*diags(i,0)*pph(3,i)**2
      endif
      do  11  i = 1, ldim
   11 sk(i,i) = sk(i,i) + wk(i)

C     call prmx('h',sk,ldim,ldim,ldim)

      call tcx('hmlt2c')

      end

