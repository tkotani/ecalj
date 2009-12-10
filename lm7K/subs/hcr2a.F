      subroutine hcr2a(sspec)
C- Set hcr according to alpha, E=0
C ----------------------------------------------------------------------
Ci Inputs
Ci   lalp  :T, generate hcr corresponding tp species alpha
Ci   loka  :Use Andersen's definitions for strux
Ci   nkap  :number of energies for which strux are calculated
Ci   kap2  :interstitial kinetic energy
Ci   kap2  :kinetic energies for which strux are calculated
Ci   sspec :struct for species-specific information; see routine uspec
Ci     Elts read: lmxa alpha hcr
Ci     Stored:    hcr
Ci     Passed to:
Co Outputs
Cl Local variables
Cl         :
Cr Remarks
Cr
Cu Updates
Cu   08 Aug 06  Redesigned 
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      double precision sspec(1)
C Local
      logical lalp
      integer lmxx,nspec,nglob,is,lmxa,igetss,l,stdo,nl,ipr
      parameter(lmxx=20)
      double precision alpha(lmxx),adot(lmxx),hcr(lmxx),kap2
      double precision fi(-1:lmxx),gi(-1:lmxx),avw,dglob

      stdo = nglob('stdo')
      call getpr(ipr)
      nspec = nglob('nspec')
      nl = nglob('nl')
      avw = dglob('avw',0d0,0)
      lalp = .true.
      kap2 = 0

      if (ipr .ge. 50) then
        call info0(50,1,0,
     .    ' Hard core radii corresponding to 2nd generation alpha''s')
        call awrit1(
     .  ' spec   energy    hcr ...%npalpha/adot',
     .  ' ',120,stdo,19+10*nl)
      endif
      do  10  is = 1, nspec

        lmxa = igetss('spec lmxa',is,sspec)
        if (lmxa .ge. lmxx) call rx('hcr2a: increase lmxx')
        if (lalp) then
          call upack2('spec alpha',sspec,is,alpha)
          call bessl2(0d0,-1,lmxa+1,fi(-1),gi(-1))
          do  12  l = 0, lmxa
   12     hcr(l+1) = avw*(alpha(l+1)/fi(l)*gi(l))**(1d0/dble(2*l+1))
          call pack2('spec hcr',sspec,is,hcr)
        endif
        call upack2('spec hcr',sspec,is,hcr)

C      This will make alphas for any energy at this hcr(l) (not used)
        do  22  l = 0, lmxa
          call bessl2(kap2*hcr(l+1)**2,-1,l+1,fi(-1),gi(-1))
          alpha(l+1) = fi(l)/gi(l)*(hcr(l+1)/avw)**(l+l+1)
          adot(l+1)  = -0.5d0*(hcr(l+1)/avw)**(l+l+3)*
     .          (fi(l+1)*gi(l)/(l+l+1)+fi(l)*gi(l-1)/(l+l-1))/
     .          (gi(l)*gi(l))
   22   continue

C       ik = 1
c       call dscal(lmxa+1,1/avw,hcr,1)
        call info8(50,0,0,
     .    '%,4i %5p%;10,6D '//
     .    '%n;10,6D'//
     .    '%47p%n;10,6D',
     .    is,kap2,lmxa+1,hcr,lmxa+1,alpha,0,0)
        call info2(50,0,0,'%47p%n;10,6D',lmxa+1,adot)

   10 continue
      end

