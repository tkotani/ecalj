subroutine show_programinfo(io)
  implicit none
  integer:: io
  character(6),parameter:: info='INFO: '
  character(256):: FFL=___FFLAGS___
  character(256):: UA=___UNAME_A___
  character(256):: EI=___ETC_ISSUE___
  character(256):: FV=___FC_VERSION___
  character(256):: PF=___PLATFORM___
  character(256):: LL=___LIBLOC___
  character(512):: GC1=___GIT_COMMIT1___
  character(512):: GC2=___GIT_COMMIT2___
  character(512):: GC3=___GIT_COMMIT3___
  write(io,'(a,a,a,a)') info,trim(EI)
  write(io,'(a,a,a,a)') info,trim(FV)
  write(io,'(a,a,a)') info, trim(FFL)
  write(io,'(a,a,a)') info,'MATH: ', trim(LL)
  write(io,'(a,a,a)') info,'git: ',trim(GC1)
  write(io,'(a,a,a)') info,'   : ',trim(GC2)
  write(io,'(a,a,a)') info,'   : ',trim(GC3)
  write(io,'(a,a,a)') info,'linked at ',trim(___LINK_TIME___)
end subroutine show_programinfo
!        program test
!       call show_programinfo(6)
!       end
