!> Get weight in window [el,eh] for a Gaussian smearing
!! with wfacx(x) = \\int_el^eh 1/sqrt(2)/esmr exp( -(x-ek/esmr)**2))
real(8) function wfacx(el,eh, ek,esmr)
  implicit none
  real(8), intent(in) :: el,eh,ek,esmr
  real(8) ::ekp,ekm,derfcx,wfacx_old
  if(esmr==0d0) then
     wfacx=0d0
     if(el <= ek .AND. ek <eh ) wfacx=1d0
     return
  endif
  wfacx = 0.5d0*erfc(-(eh-ek)/sqrt(2d0)/esmr) &
       -  0.5d0*erfc(-(el-ek)/sqrt(2d0)/esmr)
END function wfacx

!> Now the ordering of e1,e2 does not matter.May2006
!! Get weight in window [el,eh] for a Gaussian smearing
!! with wfacx(x) = \\int_el^eh 1/sqrt(2)/esmr exp( -(x-ek/esmr)**2))
real(8) function wfacx2(e1,e2, ek,esmr)
  implicit none
  real(8), intent(in) :: e1, e2, ek, esmr
  real(8) ::ekp,ekm,el,eh,derfcx,wfacx_old
  real(8):: ewidthcut=1d-6
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
  wfacx2 = 0.5d0*erfc(-(eh-ek)/sqrt(2d0)/esmr) &
       -   0.5d0*erfc(-(el-ek)/sqrt(2d0)/esmr)
END function wfacx2

!> Now the ordering of e1,e2 does not matter.May2006
!! Avaraged energy in window[el, eh] for a gaussian smearing.
!! This now stops (sanity check) if ek wtt is less than 1d-10
real(8) function weavx2(e1,e2, ek,esmr)
  implicit none
  real(8), intent(in) :: e1, e2, ek, esmr
  real(8) ::ekp,ekm,el,eh,derfcx,wtt,sig2
  real(8),parameter:: pi=3.1415926535897932d0
  el=min(e1,e2)  !May2006
  eh=max(e1,e2)
  if(esmr==0d0) then
     weavx2=1d99
     !! until 2021, we have no experience to go though this sanity check.
     !        if(el <= ek .and. ek <eh ) then
     !          weavx2=0.5d0*(el+eh)
     !          ! At feb2006, Takao found that
     !          !   "weavx=0.5d0*(el+eh)" should be (maybe) replaced by weavx=ek
     !          ! However, Takao think it is no way to do this line---so no problem.
     !          print *, " weavx2: .5*(el+eh) is wrong!"//
     !     &     " this is supposed not to go. Ask to takao:known at feb2006"
     !          call rx( " weavx2: this is not to go through. Ask to takao")
     !        endif
     return
  endif
  wtt=    0.5d0*erfc(-(eh-ek)/sqrt(2d0)/esmr) &
       -  0.5d0*erfc(-(el-ek)/sqrt(2d0)/esmr)
  sig2= 2d0*esmr**2
  if(wtt<1d-10) then
     weavx2=9999999999d99
     return
     call rx( "weavx2: wtt<1d-10")
  endif
  weavx2 = ek+ esmr/sqrt(2d0*pi) &
       *( -exp(-(eh-ek)**2/sig2) + exp(-(el-ek)**2/sig2) )/wtt
  if(weavx2<el) then
     write(6,*)'weavx2:err',el,eh, ek,weavx2
     call rx( " weavx2:  weavx < el")
  elseif(eh<weavx2) then
     write(6,*)'weavx2:err',el,eh, ek,weavx2
     call rx( " weavx2:  eh < weavx")
  endif
END function weavx2


!$$$!> Avaraged energy in window[el, eh] for a gaussian smearing.
!$$$!! This now stops (sanity check) if ek wtt is less than 1d-10
!$$$      real(8) function weavx(el,eh, ek,esmr)
!$$$      implicit none
!$$$      real(8) ::ekp,ekm,el,eh,ek,esmr,derfcx,wtt,sig2
!$$$      real(8),parameter:: pi=3.1415926535897932d0
!$$$      if(esmr==0d0) then
!$$$        weavx=1d99
!$$$        if(el <= ek .and. ek <eh ) then
!$$$          weavx=0.5d0*(el+eh)
!$$$          ! At feb2006, Takao found that
!$$$          !   "weavx=0.5d0*(el+eh)" should be (maybe) replaced by weavx=ek
!$$$          ! However, Takao think it is no way to do this line---so no problem.
!$$$          print *," weavx: .5*(el+eh) is wrong!"//
!$$$     &     " this is supposed not to go. Ask to takao:known at feb2006"
!$$$          call rx( " weavx: this is not to go through. Ask to takao")
!$$$        endif
!$$$        return
!$$$      endif
!$$$      wtt=    0.5d0*derfcx(-(eh-ek)/sqrt(2d0)/esmr)
!$$$     &       -0.5d0*derfcx(-(el-ek)/sqrt(2d0)/esmr)
!$$$      sig2= 2d0*esmr**2
!$$$      if(wtt<1d-10) then
!$$$        weavx=9999999999d99
!$$$        return
!$$$        call rx( "weavx: wtt<1d-10")
!$$$      endif
!$$$      weavx = ek+ esmr/sqrt(2d0*pi)
!$$$     &  *( -exp(-(eh-ek)**2/sig2) + exp(-(el-ek)**2/sig2) )/wtt
!$$$      if(weavx<el) then
!$$$        write(6,*)'weavx:err',el,eh, ek,weavx
!$$$        call rx( " weavx:  weavx < el")
!$$$      elseif(eh<weavx) then
!$$$        write(6,*)'weavx:err',el,eh, ek,weavx
!$$$        call rx( " weavx:  eh < weavx")
!$$$      endif
!$$$      end

!$$$      real(8) function wfacx_old(el,eh, ek,esmr)
!$$$      implicit none
!$$$      real(8) ::ekp,ekm,el,eh,ek,esmr,wfacx
!$$$      wfacx = 0d0
!$$$      if(esmr==0d0) then
!$$$        if(el<ek .and.ek<eh) then
!$$$          wfacx =1d0
!$$$        endif
!$$$        return
!$$$      endif
!$$$      ekp = ek+0.5d0*esmr
!$$$      ekm = ek-0.5d0*esmr
!$$$      if(ekm<el) then
!$$$        if(el<ekp.and.ekp<eh) then
!$$$          wfacx=(ekp-el)/esmr
!$$$        elseif(eh<=ekp) then
!$$$          wfacx=(eh-el)/esmr
!$$$        endif
!$$$      elseif(ekm < eh) then
!$$$        if(ekp<eh) then
!$$$          wfacx=1d0
!$$$        elseif(eh<=ekp) then
!$$$          wfacx=(eh-ekm)/esmr
!$$$        endif
!$$$      endif
!$$$      wfacx_old=wfacx
!$$$      end
