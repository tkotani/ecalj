!> Get weight in window [el,eh] for a Gaussian smearing
!! with wfacx(x) = \\int_el^eh 1/sqrt(2)/esmr exp( -(x-ek/esmr)**2))
      real(8) function wfacx(el,eh, ek,esmr)
      implicit none
      real(8), intent(in) :: el,eh,ek,esmr
      real(8) ::ekp,ekm,derfcx,wfacx_old
      if(esmr==0d0) then
        wfacx=0d0
        if(el <= ek .and. ek <eh ) wfacx=1d0
        return
      endif
      wfacx = 0.5d0*derfcx(-(eh-ek)/sqrt(2d0)/esmr)
     &       -0.5d0*derfcx(-(el-ek)/sqrt(2d0)/esmr)
      end

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
        if(el <= ek .and. ek <eh ) wfacx2=1d0
        return
      endif
      wfacx2 = 0.5d0*derfcx(-(eh-ek)/sqrt(2d0)/esmr)
     &       -0.5d0*derfcx(-(el-ek)/sqrt(2d0)/esmr)
      end
      
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
c        if(el <= ek .and. ek <eh ) then
c          weavx2=0.5d0*(el+eh)
c          ! At feb2006, Takao found that
c          !   "weavx=0.5d0*(el+eh)" should be (maybe) replaced by weavx=ek
c          ! However, Takao think it is no way to do this line---so no problem.
c          print *, " weavx2: .5*(el+eh) is wrong!"//
c     &     " this is supposed not to go. Ask to takao:known at feb2006"
c          call rx( " weavx2: this is not to go through. Ask to takao")
c        endif
        return
      endif
      wtt=    0.5d0*derfcx(-(eh-ek)/sqrt(2d0)/esmr)
     &       -0.5d0*derfcx(-(el-ek)/sqrt(2d0)/esmr)
      sig2= 2d0*esmr**2
      if(wtt<1d-10) then
        weavx2=9999999999d99
        return
        call rx( "weavx2: wtt<1d-10")
      endif
      weavx2 = ek+ esmr/sqrt(2d0*pi)
     &  *( -exp(-(eh-ek)**2/sig2) + exp(-(el-ek)**2/sig2) )/wtt
      if(weavx2<el) then
        write(6,*)'weavx2:err',el,eh, ek,weavx2
        call rx( " weavx2:  weavx < el")
      elseif(eh<weavx2) then
        write(6,*)'weavx2:err',el,eh, ek,weavx2
        call rx( " weavx2:  eh < weavx")
      endif
      end


c$$$!> Avaraged energy in window[el, eh] for a gaussian smearing.
c$$$!! This now stops (sanity check) if ek wtt is less than 1d-10
c$$$      real(8) function weavx(el,eh, ek,esmr)
c$$$      implicit none
c$$$      real(8) ::ekp,ekm,el,eh,ek,esmr,derfcx,wtt,sig2
c$$$      real(8),parameter:: pi=3.1415926535897932d0
c$$$      if(esmr==0d0) then
c$$$        weavx=1d99
c$$$        if(el <= ek .and. ek <eh ) then
c$$$          weavx=0.5d0*(el+eh)
c$$$          ! At feb2006, Takao found that
c$$$          !   "weavx=0.5d0*(el+eh)" should be (maybe) replaced by weavx=ek
c$$$          ! However, Takao think it is no way to do this line---so no problem.
c$$$          print *," weavx: .5*(el+eh) is wrong!"//
c$$$     &     " this is supposed not to go. Ask to takao:known at feb2006"
c$$$          call rx( " weavx: this is not to go through. Ask to takao")
c$$$        endif
c$$$        return
c$$$      endif
c$$$      wtt=    0.5d0*derfcx(-(eh-ek)/sqrt(2d0)/esmr)
c$$$     &       -0.5d0*derfcx(-(el-ek)/sqrt(2d0)/esmr)
c$$$      sig2= 2d0*esmr**2
c$$$      if(wtt<1d-10) then
c$$$        weavx=9999999999d99
c$$$        return
c$$$        call rx( "weavx: wtt<1d-10")
c$$$      endif
c$$$      weavx = ek+ esmr/sqrt(2d0*pi)
c$$$     &  *( -exp(-(eh-ek)**2/sig2) + exp(-(el-ek)**2/sig2) )/wtt
c$$$      if(weavx<el) then
c$$$        write(6,*)'weavx:err',el,eh, ek,weavx
c$$$        call rx( " weavx:  weavx < el")
c$$$      elseif(eh<weavx) then
c$$$        write(6,*)'weavx:err',el,eh, ek,weavx
c$$$        call rx( " weavx:  eh < weavx")
c$$$      endif
c$$$      end

c$$$      real(8) function wfacx_old(el,eh, ek,esmr)
c$$$      implicit none
c$$$      real(8) ::ekp,ekm,el,eh,ek,esmr,wfacx
c$$$      wfacx = 0d0
c$$$      if(esmr==0d0) then
c$$$        if(el<ek .and.ek<eh) then
c$$$          wfacx =1d0
c$$$        endif
c$$$        return
c$$$      endif
c$$$      ekp = ek+0.5d0*esmr
c$$$      ekm = ek-0.5d0*esmr
c$$$      if(ekm<el) then
c$$$        if(el<ekp.and.ekp<eh) then
c$$$          wfacx=(ekp-el)/esmr
c$$$        elseif(eh<=ekp) then
c$$$          wfacx=(eh-el)/esmr
c$$$        endif
c$$$      elseif(ekm < eh) then
c$$$        if(ekp<eh) then
c$$$          wfacx=1d0
c$$$        elseif(eh<=ekp) then
c$$$          wfacx=(eh-ekm)/esmr
c$$$        endif
c$$$      endif
c$$$      wfacx_old=wfacx
c$$$      end
