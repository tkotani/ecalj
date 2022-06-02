subroutine setofl(mode,ssite,sspec,nbas,nvl,iiv0)

  use m_struc_def  !Cgetarg

  !- Sets vector of offsets for lmxl
  ! ----------------------------------------------------------------------
  !i Inputs
  !i   mode  :not used now.
  !i   ssite :struct containing site-specific information
  !i   sspec :struct containing species-specific information
  !i   nbas  :size of basis
  !o Outputs
  !o    nvl  :dimension of array containing charge channels
  !o    iiv0 :offsets to start of channels for sites 1..nbas
  !r Remarks
  !r   Useful for, e.g. multiple threaded code.
  !u Updates
  !u   1 May 00 Adapted from nfp setofl.f
  ! ----------------------------------------------------------------------
  !     implicit none
  ! ... Passed parameters
  integer :: mode,nbas,nvl,iiv0(nbas)
  type(s_site)::ssite(*)
  type(s_spec)::sspec(*)

  ! ... Local parameters
  integer :: ib,is,iv0,lmxl,nlm,igetss

  if (mode /= 0) call rx('setofl: bad mode')
  iv0 = 0
  do  12  ib = 1, nbas
     is = int(ssite(ib)%spec)

     lmxl = int(sspec(is)%lmxl)

     nlm = (lmxl+1)**2
     iiv0(ib) = iv0
     iv0 = iv0+nlm
12 enddo
  nvl = iv0
end subroutine setofl


