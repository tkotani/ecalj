      subroutine alwats(ndimL,kap2,ldot,lmaxw,rwats,alphv,adotv)
C- Append to alpha vector alpha for Watson sphere
C ----------------------------------------------------------------------
Ci Inputs:
Ci   ndimL :offset to beginning of Watson sphere alpha
Ci   clabl :name of the different inequivalent atoms
Ci   kap2  :(kappa*average WS-radius)**2
Ci   ldot  :=0 do not generate adotv
Ci         :=1 return in adotv energy derivative of alphv
Ci         :>1 return in adotv energy derivative of alphv^-1
Ci   lmaxw :maximum l for Watson-sphere
Co Outputs:
Co   alphv :screening constants for Watson-sphere orbitals
Co   adotv :(kappa*avw)^2-derivative of alphv
Co   rwats :Watson sphere radius
Cr Remarks:
Cr   First, the distance from central atom to outermost atom is found.
Cr   Then the screening constants of the Watson orbitals are calculated.
Cr   N^al_omega(rwats)=N^0_omega(rwats)-alphv*J^0_omega(rwats)=!0
Cu Updates
Cu   20 Sep 00 return in adotv e-derivative of alphv^-1 for ldot=2
C ----------------------------------------------------------------------
C     implicit none
C Passed variables:
      integer lmaxw,ndimL,ldot
      double precision kap2,rwats,alphv(*),adotv(*)
C Local variables:
      integer iprint,l,lgunit,nlmax,ilm,m,ilmd,ipr,stdo
      parameter (nlmax=20)
      double precision fi(-1:nlmax),gi(-1:nlmax),drr2,fac,r2,rfac,
     .                 alphom,adoth,alpl(0:nlmax),adotl(0:nlmax)
C External calls:
      external  bessl2,drr2,iprint,lgunit
C Intrinsic functions:
      intrinsic dsqrt

      call getpr(ipr)
      stdo = lgunit(1)

      if (ipr .gt. 40) then
        call awrit4(' alwats:  kap = %d  rwats = %d  lmax = %i'//
     .    '%?#n>=50#%N    l    alph(hom)    adot(hom)##'
     .    ,' ',100,stdo,kap2,rwats,lmaxw,ipr)
      endif

      r2 = rwats**2
      rfac = 1d0/rwats

C --- Make Watson alpha's, and possibly energy derivatives ---
      call bessl2(kap2*r2,-1,lmaxw+1,fi(-1),gi(-1))
      ilm = 0
      do  20  l = 0, lmaxw
        rfac = rfac*r2
        alphom = fi(l)/gi(l)*rfac
        alphom = 1/alphom
        ilmd = ilm
        alpl(l) = alphom
        do  22  m = -l, l
        ilm = ilm+1
   22   alphv(ilm+ndimL) = alphom
        if (ldot .ne. 0) then
          adoth = -0.5d0*rfac*r2*(fi(l+1)*gi(l)/(l+l+1)
     .      +fi(l)*gi(l-1)/(l+l-1))/(gi(l)*gi(l))
          adoth = -adoth*alphom*alphom
          adotl(l) = adoth
          fac = 1
          do  24  m = -l, l
          ilmd = ilmd+1
          if (ldot .gt. 1) then
            fac = 1/alphv(ilmd+ndimL)**2
          endif
   24     adotv(ilmd+ndimL) = adoth*fac
          if (ipr .ge. 50)
     .    call awrit3('%,5i   %;10F   %;10F',' ',80,stdo,l,alphom,adoth)
        elseif (ipr .ge. 50) then
          call awrit2('%,5i   %;10F',' ',80,stdo,l,alphom)
        endif
   20 continue

C ... Printout
      if (ipr .gt. 40 .and. ipr .lt. 50) then
      call awrit2(' ahom =%n:1;7F',' ',80,stdo,lmaxw,alpl)
      if(ldot.ne.0)call awrit2(' adot =%n:1;7F',' ',80,stdo,lmaxw,adotl)
      endif

      end

