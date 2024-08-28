# 1 "../subroutines/show_programinfo.template"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 1 "<command-line>" 2
# 1 "../subroutines/show_programinfo.template"
subroutine show_programinfo(io)
  implicit none
  integer:: io
  character(6),parameter:: info='INFO: '
  character(256):: FFL=""
  character(256):: UA="kugui1 4.18.0-513.24.1.el8_9.x86_64 x86_64 GNU/Linux"
  character(256):: EI="\S"
  character(256):: FV=""
  character(256):: PF=""
  character(256):: LL=""
  character(512):: GC1="commit d33ca4ffa96a8ca87643b5f87fe68a5954c4cc69"
  character(512):: GC2="Author: Masao Obata <obata.mso@gmail.com>"
  character(512):: GC3="Date:   Tue Aug 27 20:56:12 2024 +0900"
  write(io,'(a,a,a,a)') info,trim(EI)
  write(io,'(a,a,a,a)') info,trim(FV)
  write(io,'(a,a,a)') info, trim(FFL)
  write(io,'(a,a,a)') info,'MATH: ', trim(LL)
  write(io,'(a,a,a)') info,'git: ',trim(GC1)
  write(io,'(a,a,a)') info,'   : ',trim(GC2)
  write(io,'(a,a,a)') info,'   : ',trim(GC3)
  write(io,'(a,a,a)') info,'linked at ',trim("Wed Aug 28 12:09:35 JST 2024")
end subroutine show_programinfo
! program test
! call show_programinfo(6)
! end
