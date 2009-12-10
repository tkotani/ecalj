      subroutine setofl(mode,ssite,sspec,nbas,nvl,iiv0)
C- Sets vector of offsets for lmxl
C ----------------------------------------------------------------------
Ci Inputs
Ci   mode  :not used now.
Ci   ssite :struct containing site-specific information
Ci   sspec :struct containing species-specific information
Ci   nbas  :size of basis
Co Outputs
Co    nvl  :dimension of array containing charge channels
Co    iiv0 :offsets to start of channels for sites 1..nbas
Cr Remarks
Cr   Useful for, e.g. multiple threaded code.
Cu Updates
Cu   1 May 00 Adapted from nfp setofl.f
C ----------------------------------------------------------------------
C     implicit none
C ... Passed parameters
      integer mode,nbas,nvl,iiv0(nbas)
      double precision ssite(1),sspec(1)
C ... Local parameters
      integer ib,is,iv0,lmxl,nlm,igetss

      if (mode .ne. 0) call rx('setofl: bad mode')
      iv0 = 0
      do  12  ib = 1, nbas
        is = igetss('site spec',ib,ssite)
        lmxl = igetss('spec lmxl',is,sspec)
        nlm = (lmxl+1)**2
        iiv0(ib) = iv0
        iv0 = iv0+nlm
   12 continue
      nvl = iv0
      end

