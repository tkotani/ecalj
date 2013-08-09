      real(8) function wfacx2(e1,e2, ek,esmr)
c Now the ordering of e1,e2 does not matter.May2006
c Get weight in window [el,eh] for a Gaussian smearing
c with wfacx(x) = \int_el^eh 1/sqrt(2)/esmr exp( -(x-ek/esmr)**2))
      implicit none
      real(8) ::ekp,ekm,el,eh,ek,esmr,derfcx,wfacx_old,e1,e2
      logical:: GaussSmear !external
      real(8):: ewidthcut=1d-6
      el=min(e1,e2)  !May2006
      eh=max(e1,e2)
      if(eh-el< ewidthcut) then !July2006
        wfacx2=0d0
        return
      endif
      if(.not.GaussSmear()) then
        wfacx2 = wfacx_old(el,eh, ek,esmr)
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

      real(8) function weavx2(e1,e2, ek,esmr)
c Now the ordering of e1,e2 does not matter.May2006
c Avaraged energy in window[el, eh] for a gaussian smearing.
c This now stops (sanity check) if ek wtt is less than 1d-10
      implicit none
      real(8) ::ekp,ekm,el,eh,ek,esmr,derfcx,wtt,sig2,e1,e2
      real(8),parameter:: pi=3.1415926535897932d0
      el=min(e1,e2)  !May2006
      eh=max(e1,e2)
      if(esmr==0d0) then
        weavx2=1d99
        if(el <= ek .and. ek <eh ) then
          weavx2=0.5d0*(el+eh)
          ! At feb2006, Takao found that
          !   "weavx=0.5d0*(el+eh)" should be (maybe) replaced by weavx=ek
          ! However, Takao think it is no way to do this line---so no problem.
          print *, " weavx2: .5*(el+eh) is wrong!"//
     &     " this is supposed not to go. Ask to takao:known at feb2006"
Cstop2rx 2013.08.09 kino          stop " weavx2: this is not to go through. Ask to takao"
          call rx( " weavx2: this is not to go through. Ask to takao")
        endif
        return
      endif
      wtt=    0.5d0*derfcx(-(eh-ek)/sqrt(2d0)/esmr)
     &       -0.5d0*derfcx(-(el-ek)/sqrt(2d0)/esmr)
      sig2= 2d0*esmr**2
      if(wtt<1d-10) then
        weavx2=9999999999d99
        return
Cstop2rx 2013.08.09 kino        stop "weavx2: wtt<1d-10" !feb2006 I changed---but this is just for sanity check.
        call rx( "weavx2: wtt<1d-10")
      endif
      weavx2 = ek+ esmr/sqrt(2d0*pi)
     &  *( -exp(-(eh-ek)**2/sig2) + exp(-(el-ek)**2/sig2) )/wtt
      if(weavx2<el) then
        write(6,*)'weavx2:err',el,eh, ek,weavx2
Cstop2rx 2013.08.09 kino        stop " weavx2:  weavx < el"
        call rx( " weavx2:  weavx < el")
      elseif(eh<weavx2) then
        write(6,*)'weavx2:err',el,eh, ek,weavx2
Cstop2rx 2013.08.09 kino        stop " weavx2:  eh < weavx"
        call rx( " weavx2:  eh < weavx")
      endif
      end



      real(8) function wfacx(el,eh, ek,esmr)
c Get weight in window [el,eh] for a Gaussian smearing
c with wfacx(x) = \int_el^eh 1/sqrt(2)/esmr exp( -(x-ek/esmr)**2))
      implicit none
      real(8) ::ekp,ekm,el,eh,ek,esmr,derfcx,wfacx_old
      logical:: GaussSmear !external

      if(.not.GaussSmear()) then
        wfacx = wfacx_old(el,eh, ek,esmr)
        return
      endif
      if(esmr==0d0) then
        wfacx=0d0
        if(el <= ek .and. ek <eh ) wfacx=1d0
        return
      endif
      wfacx = 0.5d0*derfcx(-(eh-ek)/sqrt(2d0)/esmr)
     &       -0.5d0*derfcx(-(el-ek)/sqrt(2d0)/esmr)
      end

      real(8) function weavx(el,eh, ek,esmr)
c Avaraged energy in window[el, eh] for a gaussian smearing.
c This now stops (sanity check) if ek wtt is less than 1d-10
      implicit none
      real(8) ::ekp,ekm,el,eh,ek,esmr,derfcx,wtt,sig2
      real(8),parameter:: pi=3.1415926535897932d0
      if(esmr==0d0) then
        weavx=1d99
        if(el <= ek .and. ek <eh ) then
          weavx=0.5d0*(el+eh)
          ! At feb2006, Takao found that
          !   "weavx=0.5d0*(el+eh)" should be (maybe) replaced by weavx=ek
          ! However, Takao think it is no way to do this line---so no problem.
          print *," weavx: .5*(el+eh) is wrong!"//
     &     " this is supposed not to go. Ask to takao:known at feb2006"
Cstop2rx 2013.08.09 kino          stop " weavx: this is not to go through. Ask to takao"
          call rx( " weavx: this is not to go through. Ask to takao")
        endif
        return
      endif
      wtt=    0.5d0*derfcx(-(eh-ek)/sqrt(2d0)/esmr)
     &       -0.5d0*derfcx(-(el-ek)/sqrt(2d0)/esmr)
      sig2= 2d0*esmr**2
      if(wtt<1d-10) then
        weavx=9999999999d99
        return
Cstop2rx 2013.08.09 kino        stop "weavx: wtt<1d-10" !feb2006 I changed---but this is just for sanity check.
        call rx( "weavx: wtt<1d-10")
      endif
      weavx = ek+ esmr/sqrt(2d0*pi)
     &  *( -exp(-(eh-ek)**2/sig2) + exp(-(el-ek)**2/sig2) )/wtt
      if(weavx<el) then
        write(6,*)'weavx:err',el,eh, ek,weavx
Cstop2rx 2013.08.09 kino        stop " weavx:  weavx < el"
        call rx( " weavx:  weavx < el")
      elseif(eh<weavx) then
        write(6,*)'weavx:err',el,eh, ek,weavx
Cstop2rx 2013.08.09 kino        stop " weavx:  eh < weavx"
        call rx( " weavx:  eh < weavx")
      endif

ccccccccccccccccccccccccccccc
c      write (6,"('el eh ek es w eav=',4f10.5,2d14.6)")
c     &  el,eh,ek,esmr,wtt,weavx
c      stop 'weavx test end'
ccccccccccccccccccccccccccccc
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real(8) function wfacx_old(el,eh, ek,esmr)
      implicit none
      real(8) ::ekp,ekm,el,eh,ek,esmr,wfacx
      wfacx = 0d0
c
      if(esmr==0d0) then
        if(el<ek .and.ek<eh) then
          wfacx =1d0
        endif
        return
      endif
c
      ekp = ek+0.5d0*esmr
      ekm = ek-0.5d0*esmr
      if(ekm<el) then
        if(el<ekp.and.ekp<eh) then
          wfacx=(ekp-el)/esmr
        elseif(eh<=ekp) then
          wfacx=(eh-el)/esmr
        endif
      elseif(ekm < eh) then
        if(ekp<eh) then
          wfacx=1d0
        elseif(eh<=ekp) then
          wfacx=(eh-ekm)/esmr
        endif
      endif
      wfacx_old=wfacx
      end
