Cgetarg...info...           structure ['setofl', 'site', 'ssite']
Cgetarg...info...           structure ['setofl', 'spec', 'sspec']
Cgetarg...info...           use_to_add ['setofl', 'm_struc_def']
      subroutine setofl(mode,ssite,sspec,nbas,nvl,iiv0)

      use m_struc_def  !Cgetarg

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
Cgetarg       double precision ssite(1),sspec(1)
      type(s_site)::ssite(*)
      type(s_spec)::sspec(*)

C ... Local parameters
      integer ib,is,iv0,lmxl,nlm,igetss

      if (mode .ne. 0) call rx('setofl: bad mode')
      iv0 = 0
      do  12  ib = 1, nbas
Cgetarg         is = igetss('site spec',ib,ssite)
        is = int(ssite(ib)%spec)

Cgetarg         lmxl = igetss('spec lmxl',is,sspec)
        lmxl = int(sspec(is)%lmxl)

        nlm = (lmxl+1)**2
        iiv0(ib) = iv0
        iv0 = iv0+nlm
   12 continue
      nvl = iv0
Cgetarg       end
      end subroutine setofl


