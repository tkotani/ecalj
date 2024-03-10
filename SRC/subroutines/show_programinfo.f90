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
  character(256):: UA="tt480s 5.4.0-172-generic x86_64 GNU/Linux"
  character(256):: EI="Ubuntu 20.04.6 LTS \n \l"
  character(256):: FV="GNU Fortran (Ubuntu 9.4.0-1ubuntu1~20.04.2) 9.4.0"
  character(256):: PF=""
  character(256):: LL=""
  character(512):: GC1="commit eb3d8d88ce22e2218738fc428cd4a21858e5bd81"
  character(512):: GC2="Author: Takao Kotani <takaokotani@gmail.com>"
  character(512):: GC3="Date:   Sun Mar 10 21:43:59 2024 +0900"
  write(io,'(a,a,a,a)') info,trim(EI)
  write(io,'(a,a,a,a)') info,trim(FV)
  write(io,'(a,a,a)') info, trim(FFL)
  write(io,'(a,a,a)') info,'MATH: ', trim(LL)
  write(io,'(a,a,a)') info,'git: ',trim(GC1)
  write(io,'(a,a,a)') info,'   : ',trim(GC2)
  write(io,'(a,a,a)') info,'   : ',trim(GC3)
  write(io,'(a,a,a)') info,'linked at ',trim("Sun Mar 10 22:53:36 JST 2024")
end subroutine show_programinfo
! program test
! call show_programinfo(6)
! end
