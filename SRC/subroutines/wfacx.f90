!> Get weight in window [el,eh] for a Gaussian smearing
!! with wfacx(x) = \\int_el^eh 1/sqrt(2)/esmr exp( -(x-ek/esmr)**2))
pure real(8) function wfacx(el,eh, ek,esmr)
  implicit none
  real(8), intent(in) :: el,eh,ek,esmr
  real(8) ::ekp,ekm,derfcx,wfacx_old
  if(esmr==0d0) then
     wfacx=0d0
     if(el <= ek .AND. ek <eh ) wfacx=1d0
     return
  endif
  wfacx = .5d0*erfc(-(eh-ek)/sqrt(2d0)/esmr) - .5d0*erfc(-(el-ek)/sqrt(2d0)/esmr)
end function wfacx
module  m_wfac
  !! Get weight in window [el,eh] for a Gaussian smearing
  !! with wfacx(x) = \\int_el^eh 1/sqrt(2)/esmr exp( -(x-ek/esmr)**2))
  public wfacx2,weavx2
contains
  pure real(8) function wfacx2(e1,e2, ek,esmr)
    implicit none
    real(8), intent(in) :: e1, e2, ek, esmr
    real(8) ::ekp,ekm,el,eh,derfcx,wfacx_old
    real(8),parameter :: ewidthcut=1d-6
    el=min(e1,e2)  !May2006
    eh=max(e1,e2)
    if(eh-el< ewidthcut) then !July2006
       wfacx2=0d0
       return
    endif
    if(esmr==0d0) then
       wfacx2=0d0
       if(el <= ek .AND. ek <eh ) wfacx2=1d0
       return
    endif
    wfacx2 = .5d0*erfc(-(eh-ek)/sqrt(2d0)/esmr) - .5d0*erfc(-(el-ek)/sqrt(2d0)/esmr)
  END function wfacx2
  !! Avaraged energy in window[el, eh] for a gaussian smearing.
  !! This now stops (sanity check) if ek wtt is less than 1d-10
  pure real(8) function weavx2(e1,e2, ek,esmr)
    implicit none
    real(8), intent(in) :: e1, e2, ek, esmr
    real(8) ::ekp,ekm,el,eh,derfcx,wtt,sig2
    real(8),parameter:: pi=3.1415926535897932d0
    el=min(e1,e2)  
    eh=max(e1,e2)
    wtt=    0.5d0*erfc(-(eh-ek)/sqrt(2d0)/esmr) -0.5d0*erfc(-(el-ek)/sqrt(2d0)/esmr)
    sig2= 2d0*esmr**2
    weavx2 = ek+ esmr/sqrt(2d0*pi) *( -exp(-(eh-ek)**2/sig2) + exp(-(el-ek)**2/sig2) )/wtt
  END function weavx2
end module m_wfac
