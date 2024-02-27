# 0 "../subroutines/show_programinfo.template"
# 0 "<built-in>"
# 0 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 0 "<command-line>" 2
# 1 "../subroutines/show_programinfo.template"
subroutine show_programinfo(io)
  implicit none
  integer:: io
  character(6),parameter:: info='INFO: '
  character(256):: FFL=""
  character(256):: UA="kr1 6.2.0-37-generic x86_64 GNU/Linux"
  character(256):: EI="Ubuntu 22.04.3 LTS \n \l"
  character(256):: FV="GNU Fortran (Ubuntu 11.4.0-1ubuntu1~22.04) 11.4.0"
  character(256):: PF=""
  character(256):: LL=""
  character(512):: GC1="commit 3ce8043ade7c346cf1a72a412dc87f7e1f2cd3f7"
  character(512):: GC2="Author: Takao Kotani <takaokotani@gmail.com>"
  character(512):: GC3="Date:   Sat Feb 24 22:00:30 2024 +0900"
  write(io,'(a,a,a,a)') info,trim(EI)
  write(io,'(a,a,a,a)') info,trim(FV)
  write(io,'(a,a,a)') info, trim(FFL)
  write(io,'(a,a,a)') info,'MATH: ', trim(LL)
  write(io,'(a,a,a)') info,'git: ',trim(GC1)
  write(io,'(a,a,a)') info,'   : ',trim(GC2)
  write(io,'(a,a,a)') info,'   : ',trim(GC3)
  write(io,'(a,a,a)') info,'linked at ',trim("Tue Feb 27 20:10:46 JST 2024")
end subroutine show_programinfo
! program test
! call show_programinfo(6)
! end
