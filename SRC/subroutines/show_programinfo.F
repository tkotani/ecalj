      subroutine show_programinfo(io)
      implicit none
      integer:: io
      character(6),parameter:: info='INFO: '
      character(256):: FFL=___FFLAGS___
      character(256):: UA=___UNAME_A___ 
      character(256):: LL=___LIBLOC___ 
      character(512):: GC1=___GIT_COMMIT1___
      character(512):: GC2=___GIT_COMMIT2___
      character(512):: GC3=___GIT_COMMIT3___
      write(io,'(a,a,a)') info,'system:',trim(UA)//' '//___ETC_ISSUE___
      write(io,'(a,a,a)') info,'compilar: ',___PLATFORM___//' '//___FC_VERSION___
      write(io,'(a,a,a)') info,'FFLAGS: ',trim(FFL)
      write(io,'(a,a,a)') info,'LIBLOC: ',trim(LL)
      write(io,'(a,a,a)') info,'git: ',trim(GC1)
      write(io,'(a,a,a)') info,'   : ',trim(GC2)
      write(io,'(a,a,a)') info,'   : ',trim(GC3)
      write(io,'(a,a,a)') info,'linked at: ',___LINK_TIME___
      end subroutine show_programinfo

c        program test
c       call show_programinfo(6)
c       end
